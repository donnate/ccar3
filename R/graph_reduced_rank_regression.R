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


cca_graph_rrr_cv_folds <- function(X, Y, Gamma,
                                Sx = NULL, Sy = NULL, kfolds = 5,
                                lambda = 0.01, r = 2, 
                                standardize = FALSE,
                                preprocess = NULL,
                                LW_Sy = TRUE, rho = 10,
                                niter = 1e4, thresh = 1e-4,
                                thresh_0 = 1e-6,
                                Gamma_dagger = NULL,
                                folds = NULL,
                                return_fold_values = FALSE) {
  preprocess_mode <- .resolve_preprocess_mode(standardize = standardize, preprocess = preprocess)
  if (is.null(folds)) {
    folds <- caret::createFolds(1:nrow(Y), k = kfolds, list = TRUE)
  }

  rmse <- foreach(i = seq_along(folds), .combine = c, .packages = c('Matrix')) %do% {
    n_full <- nrow(X)
    X_train <- X[-folds[[i]], ]
    Y_train <- Y[-folds[[i]], ]
    X_val <- X[folds[[i]], ]
    Y_val <- Y[folds[[i]], ]

   # We don't even need to create X_train and Y_train explicitly for the covariance calculation
    n_train <- n_full - nrow(X_val)

    # --- DOWNDATE TRICK ---
    #Sx_train <- if (is.null(Sx)) NULL else (n_full * Sx - crossprod(X_val)) / n_train
    #Sy_train <- if (is.null(Sy)) NULL else (n_full * Sy - crossprod(Y_val)) / n_train


    
    tryCatch({
      fit <- cca_graph_rrr(X_train, Y_train, Gamma,
                           Sx = NULL, Sy = NULL,
                           lambda = lambda, r = r,
                           standardize = standardize, preprocess = preprocess_mode,
                           LW_Sy = LW_Sy, rho = rho, niter = niter, thresh = thresh,
                           thresh_0 = thresh_0,
                           Gamma_dagger = Gamma_dagger)

      print(paste("Fold completed", i, ". Evaluating on validation set..."))
      print(colMeans(X_val %*% fit$U))
      print(colMeans(Y_val %*% fit$V))

      sqrt(mean((X_val %*% fit$U - Y_val %*% fit$V)^2))
    }, error = function(e) {
      message("Error in fold ", i, ": ", conditionMessage(e))
      return(NA)
    })
  }

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

  folds <- caret::createFolds(1:nrow(Y), k = kfolds, list = TRUE)

  safe_mean <- function(x) {
    out <- mean(x, na.rm = TRUE)
    if (is.nan(out)) NA_real_ else out
  }

  safe_sd <- function(x) {
    out <- stats::sd(x, na.rm = TRUE)
    if (is.nan(out)) NA_real_ else out
  }

  run_cv <- function(lambda) {
    fold_values <- cca_graph_rrr_cv_folds(
      X, Y, Gamma,
      Sx = NULL, Sy = NULL,
      kfolds = kfolds,
      lambda = lambda, r = r,
      preprocess = "none",
      LW_Sy = LW_Sy, rho = rho,
      niter = niter, thresh = thresh,
      thresh_0 = thresh_0,
      Gamma_dagger = Gamma_dagger,
      folds = folds,
      return_fold_values = TRUE
    )

    list(
      summary = data.frame(
        lambda = lambda,
        rmse = safe_mean(fold_values),
        se = safe_sd(fold_values) / sqrt(sum(!is.na(fold_values))),
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

        doParallel::registerDoParallel(cl)
        on.exit(parallel::stopCluster(cl), add = TRUE)
      } else {
        # If setup_parallel_backend returned NULL, print a warning and proceed serially
        warning("All parallel setup attempts failed. Proceeding in serial mode.", immediate. = TRUE)
        parallelize <- FALSE # Ensure %dopar% runs serially
    }
   }

   if (parallelize) {
    results_list <- foreach(lambda = lambdas, .packages = c("Matrix", "caret")) %dopar% run_cv(lambda)
  } else {
    results_list <- lapply(lambdas, run_cv)
  }

  cv_summary <- dplyr::bind_rows(lapply(results_list, `[[`, "summary"))
  cv_folds <- dplyr::bind_rows(lapply(results_list, `[[`, "fold"))

  cv_summary$rmse[is.na(cv_summary$rmse) | cv_summary$rmse == 0] <- 1e8
  cv_summary <- cv_summary %>% dplyr::filter(rmse > 1e-5)

  opt_lambda <- cv_summary$lambda[which.min(cv_summary$rmse)]
  opt_lambda <- ifelse(is.na(opt_lambda), 0.1, opt_lambda)

  final <- cca_graph_rrr(X, Y, Gamma, Sx = NULL, Sy = NULL,
                         lambda = opt_lambda, r = r,
                         preprocess = "none",
                         LW_Sy = LW_Sy, rho = rho, niter = niter,
                         thresh = thresh, Gamma_dagger = Gamma_dagger, thresh_0=thresh_0)

  list(
    U = final$U,
    V = final$V,
    lambda = opt_lambda,
    rmse = cv_summary$rmse,
    cor = sapply(1:r, function(i) cov(X %*% final$U[, i], Y %*% final$V[, i])),
    lambda_x = opt_lambda,
    lambda_y = NA_real_,
    lambda_x_se = cv_summary$se[match(opt_lambda, cv_summary$lambda)],
    lambda_y_se = NA_real_,
    cv_summary = cv_summary,
    cv_folds = cv_folds,
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
#' @param Sx Optional covariance matrix for X. If NULL, computed as t(X) %*% X / n
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

  sx_svd <- NULL
  if (is.null(Sx)) {
    if (p > n) {
      # For p > n, avoid explicitly forming p x p covariance.
      sx_svd <- svd(X / sqrt(n))
    } else {
      Sx <- crossprod(X) / n
    }
  }

  if (is.null(Sy)) {
    Sy <- crossprod(Y) / n
    if (LW_Sy) Sy <- as.matrix(corpcor::cov.shrink(Y))
  }

  # Inverse square root of Sy
  svd_Sy <- svd(as.matrix(Sy))
  sqrt_inv_Sy <- svd_Sy$u %*% diag(ifelse(svd_Sy$d > 1e-4, 1 / sqrt(svd_Sy$d), 0)) %*% t(svd_Sy$v)
  tilde_Y <- Y %*% sqrt_inv_Sy

  gamma_ridge <- max(thresh_0, 1e-6)

  if (is.null(Gamma_dagger)) {
    if (!inherits(Gamma, "sparseMatrix")) {
      Gamma <- Matrix::Matrix(Gamma, sparse = TRUE)
    }

    L_eps <- Matrix::crossprod(Gamma) + gamma_ridge * Matrix::Diagonal(n = ncol(Gamma))
    L_fac <- Matrix::Cholesky(L_eps, LDL = FALSE)

    solve_L <- function(rhs) {
      as.matrix(Matrix::solve(L_fac, rhs, system = "A"))
    }

    # Apply the regularized pseudoinverse without explicitly forming it.
    XGamma_dagger <- as.matrix(Matrix::t(Gamma %*% solve_L(t(X))))

    # XPi = X (I - Gamma^dagger Gamma) = gamma_ridge * X (Gamma'Gamma + gamma_ridge I)^(-1)
    XPi <- gamma_ridge * t(solve_L(t(X)))
  } else {
    XGamma_dagger <- as.matrix(X %*% Gamma_dagger)
    XPi <- as.matrix(X - XGamma_dagger %*% Gamma)
    solve_L <- NULL
  }

  # Remove projection on Pi
  Projection <- tryCatch(
    qr.solve(XPi, tilde_Y),
    error = function(e) {
      XtX <- crossprod(XPi) + gamma_ridge * diag(ncol(XPi))
      solve(XtX, crossprod(XPi, tilde_Y))
    }
  )
  new_Ytilde <- tilde_Y - XPi %*% Projection
  if (verbose) {
    cat("Norm of Projection:", fro_norm(Projection), "
")
  }

  # ADMM initialization
  new_p <- ncol(XGamma_dagger)
  prod_xy <- crossprod(XGamma_dagger, new_Ytilde) / n

  if (is.null(Gamma_dagger)) {
    A <- XGamma_dagger / sqrt(n)
    K <- diag(nrow(A)) + tcrossprod(A) / rho
    solve_invSx <- function(rhs) {
      rhs <- as.matrix(rhs)
      tmp <- A %*% rhs
      middle <- solve(K, tmp)
      rhs / rho - crossprod(A, middle) / (rho^2)
    }
  } else {
    invSx <- solve(crossprod(XGamma_dagger) / n + rho * diag(new_p))
    solve_invSx <- function(rhs) invSx %*% rhs
  }

  U_dual <- matrix(0, new_p, q)
  Z <- matrix(0, new_p, q)
  B <- matrix(0, new_p, q)

  for (i in seq_len(niter)) {
    Z_old <- Z

    B <- solve_invSx(prod_xy + rho * (Z - U_dual))
    Z <- B + U_dual

    # Group soft-thresholding
    norms <- sqrt(rowSums(Z^2))
    shrink_factors <- pmax(0, 1 - lambda / (rho * norms))
    shrink_factors[is.nan(shrink_factors)] <- 0
    Z <- Z * shrink_factors

    U_dual <- U_dual + B - Z

    primal <- fro_norm(Z - B)
    dual <- fro_norm(Z_old - Z)
    if (verbose) {
      cat("ADMM iter:", i, "Primal:", primal, "Dual:", dual, "
")
    }
    if (max(primal / sqrt(new_p), dual / sqrt(new_p)) < thresh) break
  }

  # Reconstruct full coefficient matrix
  if (is.null(Gamma_dagger)) {
    gamma_dagger_B <- solve_L(crossprod(Gamma, B))
    Pi_projection <- gamma_ridge * solve_L(Projection)
  } else {
    gamma_dagger_B <- Gamma_dagger %*% B
    Pi_projection <- Projection - as.matrix(Gamma_dagger %*% (Gamma %*% Projection))
  }

  B_opt <- gamma_dagger_B + Pi_projection
  B_opt[abs(B_opt) < thresh_0] <- 0

  if (verbose) {
    cat("Norm of Pi Projection:", fro_norm(Pi_projection), "
")
    cat("Norm of Gamma_dagger %*% B:", fro_norm(gamma_dagger_B), "
")
    cat("Norm of XPi:", fro_norm(XPi), "
")
    cat("rank of XPi:", qr(XPi)$rank, "
")
  }

  # Final CCA step: apply sqrt(Sx) without explicitly building Sx when p > n.
  if (!is.null(sx_svd)) {
    sqrt_Sx_B <- sx_svd$v %*% (sx_svd$d * crossprod(sx_svd$v, B_opt))
  } else {
    svd_Sx <- svd(as.matrix(Sx))
    sqrt_Sx_B <- svd_Sx$u %*% (sqrt(pmax(svd_Sx$d, 0)) * crossprod(svd_Sx$u, B_opt))
  }

  sol <- svd(sqrt_Sx_B)
  if (verbose) {
    cat("Singular values of sqrt(Sx) %*% B:
")
    print(sol$d)
  }
  V <- sqrt_inv_Sy %*% sol$v[, 1:r]
  inv_D <- diag(ifelse(sol$d[1:r] > 1e-4, 1 / sol$d[1:r], 0))
  U <- B_opt %*% sol$v[, 1:r] %*% inv_D

  list(
    U = U,
    V = V,
    loss = mean((Y %*% V - X %*% U)^2),
    cor = sapply(1:r, function(i) cov(X %*% U[, i], Y %*% V[, i]))
  )
}
