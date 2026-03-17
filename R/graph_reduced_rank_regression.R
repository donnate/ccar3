library(magrittr)
library(tidyr)
library(pracma) # Required for pinv
library(caret) # Required for createFolds
library(foreach)

.resolve_preprocess_mode <- function(standardize = NULL, preprocess = NULL) {
  if (!is.null(preprocess)) {
    return(match.arg(preprocess, c("scale", "center", "none")))
  }

  if (is.character(standardize)) {
    return(match.arg(standardize, c("scale", "center", "none")))
  }

  if (is.logical(standardize) && length(standardize) == 1L) {
    return(if (isTRUE(standardize)) "scale" else "center")
  }

  if (is.null(standardize)) {
    return("center")
  }

  stop("`standardize` must be logical or one of: 'scale', 'center', 'none'.", call. = FALSE)
}

.preprocess_matrix <- function(M, mode) {
  M <- as.matrix(M)
  if (mode == "none") {
    return(M)
  }

  M <- scale(M, center = TRUE, scale = identical(mode, "scale"))
  M <- as.matrix(M)
  M[!is.finite(M)] <- 0
  M
}

.as_double_matrix <- function(M) {
  M <- as.matrix(M)
  storage.mode(M) <- "double"
  M
}

.make_graph_invSx_apply <- function(A, rho, ridge = 1e-8, verbose = FALSE) {
  A <- .as_double_matrix(A)
  n <- nrow(A)
  p <- ncol(A)

  if (p <= n) {
    gram <- crossprod(A) + rho * diag(p)
    R <- tryCatch(chol(gram + ridge * diag(p)), error = function(e) NULL)

    solve_gram <- if (!is.null(R)) {
      function(B) backsolve(R, forwardsolve(t(R), B))
    } else {
      if (verbose) message("Cholesky failed in graph ADMM solve; using solve().")
      function(B) solve(gram, B)
    }

    return(list(
      apply = function(rhs) solve_gram(.as_double_matrix(rhs)),
      uses_woodbury = FALSE
    ))
  }

  K <- diag(n) + tcrossprod(A) / rho
  R <- tryCatch(chol(K + ridge * diag(n)), error = function(e) NULL)

  solve_K <- if (!is.null(R)) {
    function(B) backsolve(R, forwardsolve(t(R), B))
  } else {
    if (verbose) message("Cholesky failed in graph Woodbury solve; using solve().")
    function(B) solve(K, B)
  }

  list(
    apply = function(rhs) {
      rhs <- .as_double_matrix(rhs)
      middle <- solve_K(A %*% rhs)
      rhs / rho - crossprod(A, middle) / (rho^2)
    },
    uses_woodbury = TRUE
  )
}

.prepare_graph_cache <- function(Gamma, gamma_ridge, Gamma_dagger = NULL) {
  if (is.null(Gamma_dagger)) {
    Gamma <- if (inherits(Gamma, "dgCMatrix")) {
      Gamma
    } else {
      methods::as(Matrix::Matrix(Gamma, sparse = TRUE), "dgCMatrix")
    }

    L_eps <- Matrix::crossprod(Gamma) + gamma_ridge * Matrix::Diagonal(n = ncol(Gamma))
    L_fac <- Matrix::Cholesky(L_eps, LDL = FALSE)

    solve_L <- function(rhs) {
      as.matrix(Matrix::solve(L_fac, .as_double_matrix(rhs), system = "A"))
    }

    return(list(
      mode = "regularized",
      Gamma = Gamma,
      gamma_ridge = gamma_ridge,
      solve_L = solve_L
    ))
  }

  list(
    mode = "explicit",
    Gamma = Gamma,
    Gamma_dagger = Gamma_dagger,
    gamma_ridge = gamma_ridge
  )
}

.transform_graph_design <- function(graph_cache, X = NULL, solved_tX = NULL) {
  if (graph_cache$mode == "regularized") {
    if (is.null(solved_tX)) {
      if (is.null(X)) stop("Either `X` or `solved_tX` must be provided.", call. = FALSE)
      solved_tX <- graph_cache$solve_L(t(X))
    }

    return(list(
      XGamma_dagger = as.matrix(Matrix::t(graph_cache$Gamma %*% solved_tX)),
      XPi = graph_cache$gamma_ridge * t(solved_tX),
      solved_tX = solved_tX
    ))
  }

  if (is.null(X)) {
    stop("`X` must be provided when using an explicit Gamma_dagger.", call. = FALSE)
  }

  XGamma_dagger <- as.matrix(X %*% graph_cache$Gamma_dagger)
  list(
    XGamma_dagger = XGamma_dagger,
    XPi = as.matrix(X - XGamma_dagger %*% graph_cache$Gamma)
  )
}

.reconstruct_graph_coefficients <- function(graph_cache, B, Projection) {
  if (graph_cache$mode == "regularized") {
    return(list(
      gamma_dagger_B = graph_cache$solve_L(Matrix::crossprod(graph_cache$Gamma, .as_double_matrix(B))),
      Pi_projection = graph_cache$gamma_ridge * graph_cache$solve_L(Projection)
    ))
  }

  list(
    gamma_dagger_B = graph_cache$Gamma_dagger %*% B,
    Pi_projection = Projection - as.matrix(graph_cache$Gamma_dagger %*% (graph_cache$Gamma %*% Projection))
  )
}

.project_graph_rrr_coefficients <- function(Bhat, X, r, Sx = NULL) {
  r_eff <- min(r, nrow(Bhat), ncol(Bhat))
  if (r_eff < 1 || r_eff >= ncol(Bhat)) {
    return(Bhat)
  }

  if (is.null(Sx)) {
    sx_svd <- svd(X / sqrt(nrow(X)), nu = 0, nv = min(nrow(X), ncol(X)))
    sqrt_Sx_B <- sx_svd$v %*% (sx_svd$d * crossprod(sx_svd$v, Bhat))
  } else {
    svd_Sx <- svd(as.matrix(Sx), nu = 0, nv = ncol(Sx))
    sqrt_Sx_B <- svd_Sx$v %*% (sqrt(pmax(svd_Sx$d, 0)) * crossprod(svd_Sx$v, Bhat))
  }

  weighted_svd <- .rrr_safe_rank_svd(sqrt_Sx_B, r_eff, prefer_sparse = FALSE)
  if (is.null(weighted_svd)) {
    return(Bhat)
  }

  Vr <- weighted_svd$v[, seq_len(r_eff), drop = FALSE]
  Bhat %*% Vr %*% t(Vr)
}

.fit_graph_rrr_preprocessed <- function(X, Y,
                                        Sx = NULL, Sy = NULL,
                                        lambda = 0, r,
                                        LW_Sy = TRUE, rho = 10,
                                        niter = 1e4, thresh = 1e-4,
                                        thresh_0 = 1e-6,
                                        verbose = FALSE,
                                        graph_cache,
                                        graph_design = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)

  if (q > n) {
    stop("The X/Y swapping heuristic (for q > n) is not compatible with the graph constraint Gamma. Swap X and Y in your arguments.")
  }

  if (n < min(p, q)) {
    warning("Both X and Y are high dimensional; method may be unstable.")
  }

  fro_norm <- function(M) sqrt(sum(M^2))

  if (is.null(Sy)) {
    Sy <- crossprod(Y) / n
    if (LW_Sy) Sy <- as.matrix(corpcor::cov.shrink(Y))
  }

  svd_Sy <- svd(as.matrix(Sy))
  sqrt_inv_Sy <- svd_Sy$u %*% diag(ifelse(svd_Sy$d > 1e-4, 1 / sqrt(svd_Sy$d), 0)) %*% t(svd_Sy$v)
  tilde_Y <- Y %*% sqrt_inv_Sy

  if (is.null(graph_design)) {
    graph_design <- .transform_graph_design(graph_cache, X = X)
  }

  XGamma_dagger <- graph_design$XGamma_dagger
  XPi <- graph_design$XPi

  Projection <- tryCatch(
    qr.solve(XPi, tilde_Y),
    error = function(e) {
      XtX <- .as_double_matrix(crossprod(XPi) + graph_cache$gamma_ridge * diag(ncol(XPi)))
      rhs_proj <- .as_double_matrix(crossprod(XPi, tilde_Y))
      solve(XtX, rhs_proj)
    }
  )
  new_Ytilde <- tilde_Y - XPi %*% Projection

  if (verbose) {
    cat("Norm of Projection:", fro_norm(Projection), "\n")
  }

  new_p <- ncol(XGamma_dagger)
  prod_xy <- crossprod(XGamma_dagger, new_Ytilde) / n
  invSx_info <- .make_graph_invSx_apply(XGamma_dagger / sqrt(n), rho, verbose = verbose)
  solve_invSx <- invSx_info$apply

  U_dual <- matrix(0, new_p, q)
  Z <- matrix(0, new_p, q)
  B <- matrix(0, new_p, q)

  for (i in seq_len(niter)) {
    Z_old <- Z

    B <- solve_invSx(prod_xy + rho * (Z - U_dual))
    Z <- B + U_dual

    norms <- sqrt(rowSums(Z^2))
    shrink_factors <- pmax(0, 1 - lambda / (rho * norms))
    shrink_factors[is.nan(shrink_factors)] <- 0
    Z <- Z * shrink_factors

    U_dual <- U_dual + B - Z

    primal <- fro_norm(Z - B)
    dual <- fro_norm(Z_old - Z)
    if (verbose) {
      cat("ADMM iter:", i, "Primal:", primal, "Dual:", dual, "\n")
    }
    if (max(primal / sqrt(new_p), dual / sqrt(new_p)) < thresh) break
  }

  reconstruction <- .reconstruct_graph_coefficients(graph_cache, B, Projection)
  gamma_dagger_B <- reconstruction$gamma_dagger_B
  Pi_projection <- reconstruction$Pi_projection

  B_opt <- gamma_dagger_B + Pi_projection
  B_opt[abs(B_opt) < thresh_0] <- 0

  if (verbose) {
    cat("Norm of Pi Projection:", fro_norm(Pi_projection), "\n")
    cat("Norm of Gamma_dagger %*% B:", fro_norm(gamma_dagger_B), "\n")
    cat("Norm of XPi:", fro_norm(XPi), "\n")
    cat("rank of XPi:", qr(XPi)$rank, "\n")
  }

  B_post <- .project_graph_rrr_coefficients(B_opt, X, r, Sx = Sx)
  post <- .postprocess_rrr_fit(
    B_post, X, Y, sqrt_inv_Sy, r,
    ridge_whiten = max(thresh_0, 1e-8),
    thresh_0 = thresh_0,
    verbose = verbose
  )

  list(
    U = post$U,
    V = post$V,
    loss = post$loss,
    cor = post$cor,
    Lambda = post$Lambda,
    B_opt = B_post
  )
}


cca_graph_rrr_cv_folds <- function(X, Y, Gamma,
                                Sx = NULL, Sy = NULL, kfolds = 5,
                                lambda = 0.01, r = 2, 
                                standardize = FALSE,
                                preprocess = NULL,
                                LW_Sy = TRUE, rho = 10,
                                niter = 1e4, thresh = 1e-4,
                                thresh_0 = 1e-6,
                                Gamma_dagger = NULL,
                                graph_cache = NULL,
                                solved_tX = NULL,
                                folds = NULL,
                                return_fold_values = FALSE) {
  preprocess_mode <- .resolve_preprocess_mode(standardize = standardize, preprocess = preprocess)
  if (is.null(folds)) {
    folds <- caret::createFolds(1:nrow(Y), k = kfolds, list = TRUE)
  }

  rmse <- vapply(seq_along(folds), function(i) {
    n_full <- nrow(X)
    X_train <- X[-folds[[i]], ]
    Y_train <- Y[-folds[[i]], ]
    X_val <- X[folds[[i]], ]
    Y_val <- Y[folds[[i]], ]

   # We don't even need to create X_train and Y_train explicitly for the covariance calculation
    n_train <- n_full - nrow(X_val)

    # Reuse downdated empirical covariances only when shrinkage is disabled.
    # For LW_Sy = TRUE, each fold recomputes the shrinkage estimate from Y_train.
    Sx_train <- if (is.null(Sx)) NULL else (n_full * Sx - crossprod(X_val)) / n_train
    Sy_train <- if (LW_Sy || is.null(Sy)) NULL else (n_full * Sy - crossprod(Y_val)) / n_train


    
    tryCatch({
      if (!is.null(graph_cache)) {
        graph_design <- if (!is.null(solved_tX) && graph_cache$mode == "regularized") {
          .transform_graph_design(
            graph_cache,
            solved_tX = solved_tX[, -folds[[i]], drop = FALSE]
          )
        } else {
          NULL
        }

        fit <- .fit_graph_rrr_preprocessed(
          X_train, Y_train,
          Sx = Sx_train, Sy = Sy_train,
          lambda = lambda, r = r,
          LW_Sy = LW_Sy, rho = rho,
          niter = niter, thresh = thresh,
          thresh_0 = thresh_0,
          verbose = FALSE,
          graph_cache = graph_cache,
          graph_design = graph_design
        )
      } else {
        fit <- cca_graph_rrr(X_train, Y_train, Gamma,
                             Sx = Sx_train, Sy = Sy_train,
                             lambda = lambda, r = r,
                             standardize = standardize, preprocess = preprocess_mode,
                             LW_Sy = LW_Sy, rho = rho, niter = niter, thresh = thresh,
                             thresh_0 = thresh_0,
                             Gamma_dagger = Gamma_dagger)
      }

      sqrt(mean((X_val %*% fit$U - Y_val %*% fit$V)^2))
    }, error = function(e) {
      message("Error in fold ", i, ": ", conditionMessage(e))
      NA_real_
    })
  }, numeric(1))

  if (return_fold_values) {
    return(rmse)
  }

  if (all(is.na(rmse))) return(1e8)
  mean(rmse, na.rm = TRUE)
}

#' Graph-regularized Reduced-Rank Regression for Canonical Correlation Analysis with cross validation
#'
#' Solves a sparse canonical correlation problem using a graph-constrained reduced-rank regression formulation.
#' The problem is solved via an ADMM approach.
#'
#' @param X Matrix of predictors (n x p)
#' @param Y Matrix of responses (n x q)
#' @param Gamma Graph constraint matrix (g x p)
#' @param kfolds Number of folds for cross-validation
#' @param parallelize Whether to parallelize cross-validation
#' @param lambdas Grid of regularization parameters to test for sparsity
#' @param r Target rank
#' @param standardize Backward-compatible preprocessing flag: TRUE = `"scale"`, FALSE = `"center"`.
#' @param preprocess Preprocessing mode. One of `"scale"` (center + scale), `"center"` (center only), or `"none"`.
#' @param LW_Sy Whether to apply Ledoit-Wolf shrinkage to Sy
#' @param rho ADMM penalty parameter
#' @param niter Maximum number of ADMM iterations
#' @param thresh Convergence threshold for ADMM
#' @param thresh_0 Threshold for small values in the coefficient matrix (default 1e-6)
#' @param verbose Whether to print diagnostic output
#' @param Gamma_dagger Optional pseudoinverse of Gamma (computed if NULL)
#' @param nb_cores Number of cores to use for parallelization. Defaults to min(kfolds, available cores minus 1).
#'
#' @return A list with elements:
#' \describe{
#'   \item{U}{Canonical direction matrix for X (p x r)}
#'   \item{V}{Canonical direction matrix for Y (q x r)}
#'   \item{lambda}{Optimal regularisation parameter lambda chosen by CV}
#'   \item{rmse}{Mean squared error of prediction (as computed in the CV)}
#'   \item{cor}{Canonical covariances}
#'   \item{lambda_x}{Alias of the selected `lambda`}
#'   \item{lambda_x_se}{Foldwise standard error at the selected `lambda`}
#'   \item{lambda_y}{Placeholder for symmetry with two-penalty interfaces}
#'   \item{lambda_y_se}{Placeholder for symmetry with two-penalty interfaces}
#'   \item{resultsx}{Backward-compatible alias of `cv_summary`}
#'   \item{cv_summary}{Data frame with one row per lambda containing mean RMSE and its foldwise standard error}
#'   \item{cv_folds}{Data frame with fold-level RMSE values for each lambda}
#'   \item{Lambda}{Canonical correlations from the final fit}
#'   \item{B}{Estimated reduced-rank coefficient matrix from the final fit}
#'   \item{fit}{Final fit at the selected lambda}
#' }
#' @importFrom foreach foreach %dopar%

#' @export
cca_graph_rrr_cv <- function(X, Y, Gamma,
                             r = 2, 
                             lambdas = 10^seq(-3, 1.5, length.out = 10),
                             kfolds = 5, parallelize = FALSE,
                             standardize = TRUE, LW_Sy = TRUE,
                             preprocess = NULL,
                             rho = 10, niter = 1e4, thresh = 1e-4,
                             thresh_0 = 1e-6, verbose = FALSE,
                             Gamma_dagger = NULL, nb_cores= NULL) {
  preprocess_mode <- .resolve_preprocess_mode(standardize = standardize, preprocess = preprocess)

  if (nrow(X) < min(ncol(X), ncol(Y))) {
    warning("Both X and Y are high dimensional; method may fail.")
  }

  # Apply preprocessing once before CV. Each fold then reuses these transformed matrices.
  X <- .preprocess_matrix(X, preprocess_mode)
  Y <- .preprocess_matrix(Y, preprocess_mode)

  n <- nrow(X)
  Sx <- NULL
  Sy <- if (LW_Sy) NULL else crossprod(Y) / n
  graph_cache <- .prepare_graph_cache(Gamma, max(thresh_0, 1e-6), Gamma_dagger = Gamma_dagger)
  solved_tX <- if (graph_cache$mode == "regularized") graph_cache$solve_L(t(X)) else NULL
  full_graph_design <- if (!is.null(solved_tX)) .transform_graph_design(graph_cache, solved_tX = solved_tX) else NULL

  folds <- caret::createFolds(1:nrow(Y), k = kfolds, list = TRUE)
  .rrr_chol_jitter <- .rrr_chol_jitter
  .rrr_invsqrt_psd <- .rrr_invsqrt_psd
  .rrr_whiten_factor <- .rrr_whiten_factor
  .rrr_safe_rank_svd <- .rrr_safe_rank_svd
  .postprocess_rrr_fit <- .postprocess_rrr_fit
  .project_graph_rrr_coefficients <- .project_graph_rrr_coefficients
  .resolve_preprocess_mode <- .resolve_preprocess_mode
  .preprocess_matrix <- .preprocess_matrix
  .make_graph_invSx_apply <- .make_graph_invSx_apply
  .prepare_graph_cache <- .prepare_graph_cache
  .transform_graph_design <- .transform_graph_design
  .reconstruct_graph_coefficients <- .reconstruct_graph_coefficients
  .fit_graph_rrr_preprocessed <- .fit_graph_rrr_preprocessed
  cca_graph_rrr <- cca_graph_rrr

  environment(.postprocess_rrr_fit) <- environment()
  environment(.project_graph_rrr_coefficients) <- environment()
  environment(.fit_graph_rrr_preprocessed) <- environment()
  environment(cca_graph_rrr) <- environment()

  graph_rrr_cv_folds <- cca_graph_rrr_cv_folds
  environment(graph_rrr_cv_folds) <- environment()
  safe_mean_na <- .safe_mean_na
  safe_sd_na <- .safe_sd_na

  run_cv <- function(lambda) {
    fold_values <- graph_rrr_cv_folds(
      X, Y, Gamma,
      Sx = Sx, Sy = Sy,
      kfolds = kfolds,
      lambda = lambda, r = r,
      preprocess = "none",
      LW_Sy = LW_Sy, rho = rho,
      niter = niter, thresh = thresh,
      thresh_0 = thresh_0,
      Gamma_dagger = Gamma_dagger,
      graph_cache = graph_cache,
      solved_tX = solved_tX,
      folds = folds,
      return_fold_values = TRUE
    )

    list(
      summary = data.frame(
        lambda = lambda,
        rmse = safe_mean_na(fold_values),
        se = safe_sd_na(fold_values) / sqrt(sum(!is.na(fold_values))),
        stringsAsFactors = FALSE
      ),
      fold = data.frame(
        lambda = lambda,
        fold = seq_along(fold_values),
        rmse = fold_values,
        stringsAsFactors = FALSE
      )
    )
  }

    if (parallelize) {
      if (!requireNamespace("doParallel", quietly = TRUE)) {
      stop("Package 'doParallel' must be installed to use the parallelization option.",
          call. = FALSE)
      }

      if (!requireNamespace("crayon", quietly = TRUE)) {
      stop("Package 'crayon' must be installed to use the parallelization option.",
          call. = FALSE)
      }
    # --- GRACEFUL PARALLEL SETUP ---
      available_cores_str <- Sys.getenv("SLURM_CPUS_PER_TASK")
      available_cores <- if (available_cores_str == "") {
        parallel::detectCores()
      } else {
        as.integer(available_cores_str)
      }
      available_cores <- max(1L, available_cores)
      if (is.null(nb_cores)) {
        nb_cores_effective <- min(kfolds, length(lambdas), available_cores - 1L)
      } else {
        nb_cores_effective <- min(as.integer(nb_cores), kfolds, length(lambdas))
      }
      nb_cores_effective <- max(1L, nb_cores_effective)
      cl <- setup_parallel_backend(nb_cores_effective)
      
      if (!is.null(cl)) {
        # If the cluster was created successfully, register it and plan to stop it
        if (verbose)  { cat(crayon::green("Parallel backend successfully registered.\n")) }

        initialize_parallel_workers(cl, verbose = verbose)
        doParallel::registerDoParallel(cl)
        on.exit(cleanup_parallel_backend(cl), add = TRUE)
      } else {
        # If setup_parallel_backend returned NULL, print a warning and proceed serially
        warning("All parallel setup attempts failed. Proceeding in serial mode.", immediate. = TRUE)
        parallelize <- FALSE # Ensure %dopar% runs serially
    }
   }

  graph_rrr_parallel_exports <- c(
    ".rrr_chol_jitter", ".rrr_invsqrt_psd",
    ".rrr_whiten_factor", ".rrr_safe_rank_svd",
    ".postprocess_rrr_fit", ".project_graph_rrr_coefficients",
    ".resolve_preprocess_mode", ".preprocess_matrix",
    ".make_graph_invSx_apply", ".prepare_graph_cache",
    ".transform_graph_design", ".reconstruct_graph_coefficients",
    ".fit_graph_rrr_preprocessed"
  )

   if (parallelize) {
    results_list <- foreach(
      lambda = lambdas,
      .packages = c("Matrix", "caret"),
      .export = graph_rrr_parallel_exports
    ) %dopar% run_cv(lambda)
  } else {
    results_list <- lapply(lambdas, run_cv)
  }

  cv_summary <- dplyr::bind_rows(lapply(results_list, `[[`, "summary"))
  cv_folds <- dplyr::bind_rows(lapply(results_list, `[[`, "fold"))

  cv_summary$rmse[is.na(cv_summary$rmse) | cv_summary$rmse == 0] <- 1e8
  cv_summary <- cv_summary %>% dplyr::filter(rmse > 1e-5)

  opt_lambda <- cv_summary$lambda[which.min(cv_summary$rmse)]
  opt_lambda <- ifelse(is.na(opt_lambda), 0.1, opt_lambda)

  final <- .fit_graph_rrr_preprocessed(
    X, Y,
    Sx = Sx, Sy = Sy,
    lambda = opt_lambda, r = r,
    LW_Sy = LW_Sy, rho = rho, niter = niter,
    thresh = thresh, thresh_0 = thresh_0,
    verbose = verbose,
    graph_cache = graph_cache,
    graph_design = full_graph_design
  )

  list(
    U = final$U,
    V = final$V,
    lambda = opt_lambda,
    rmse = cv_summary$rmse,
    cor = final$cor,
    lambda_x = opt_lambda,
    lambda_x_se = cv_summary$se[match(opt_lambda, cv_summary$lambda)],
    lambda_y = NA_real_,
    lambda_y_se = NA_real_,
    resultsx = cv_summary,
    cv_summary = cv_summary,
    cv_folds = cv_folds,
    Lambda = final$Lambda,
    B = final$B_opt,
    fit = final
  )
}



#' Graph-regularized Reduced-Rank Regression for Canonical Correlation Analysis
#'
#' Solves a sparse canonical correlation problem using a graph-constrained reduced-rank regression formulation.
#' The problem is solved via an ADMM approach.
#'
#' @param X Matrix of predictors (n x p)
#' @param Y Matrix of responses (n x q)
#' @param Gamma Graph constraint matrix (g x p)
#' @param Sx Optional covariance matrix for X. Kept for backward compatibility; the graph fit now postprocesses directly from `X` and does not need to form `Sx`.
#' @param Sy Optional covariance matrix for Y. If NULL, computed similarly; optionally shrunk via Ledoit-Wolf
#' @param Sxy Optional cross-covariance matrix (not currently used)
#' @param lambda Regularization parameter for sparsity
#' @param r Target rank
#' @param standardize Backward-compatible preprocessing flag: TRUE = `"scale"`, FALSE = `"center"`.
#' @param preprocess Preprocessing mode. One of `"scale"` (center + scale), `"center"` (center only), or `"none"`.
#' @param LW_Sy Whether to apply Ledoit-Wolf shrinkage to Sy
#' @param rho ADMM penalty parameter
#' @param niter Maximum number of ADMM iterations
#' @param thresh Convergence threshold for ADMM
#' @param verbose Whether to print diagnostic output
#' @param thresh_0 Threshold for small values in the coefficient matrix (default 1e-6)
#' @param Gamma_dagger Optional pseudoinverse of Gamma (computed if NULL)
#'
#' @return A list with elements:
#' \describe{
#'   \item{U}{Canonical direction matrix for X (p x r)}
#'   \item{V}{Canonical direction matrix for Y (q x r)}
#'   \item{cor}{Canonical covariances}
#'   \item{loss}{The prediction error 1/n * \| XU - YV\|^2}
#'   \item{Lambda}{Canonical correlations}
#'   \item{B_opt}{Estimated reduced-rank coefficient matrix}
#' }
#' @export
cca_graph_rrr <- function(X, Y, Gamma,
                          Sx = NULL, Sy = NULL, Sxy = NULL,
                          lambda = 0, r,
                          standardize = FALSE,
                          preprocess = NULL,
                          LW_Sy = TRUE, rho = 10,
                          niter = 1e4, thresh = 1e-4,
                          thresh_0 = 1e-6,
                          verbose = FALSE, Gamma_dagger = NULL) {
  preprocess_mode <- .resolve_preprocess_mode(standardize = standardize, preprocess = preprocess)
  X <- .preprocess_matrix(X, preprocess_mode)
  Y <- .preprocess_matrix(Y, preprocess_mode)
  graph_cache <- .prepare_graph_cache(Gamma, max(thresh_0, 1e-6), Gamma_dagger = Gamma_dagger)
  .fit_graph_rrr_preprocessed(
    X, Y,
    Sx = Sx, Sy = Sy,
    lambda = lambda, r = r,
    LW_Sy = LW_Sy, rho = rho,
    niter = niter, thresh = thresh,
    thresh_0 = thresh_0,
    verbose = verbose,
    graph_cache = graph_cache
  )
}
