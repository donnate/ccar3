.resolve_group_preprocess_mode <- function(standardize = NULL, preprocess = NULL) {
  if (is.logical(preprocess) && length(preprocess) == 1L) {
    if (isTRUE(preprocess)) {
      preprocess <- NULL
    } else {
      return("none")
    }
  }

  if (!is.null(preprocess)) {
    return(match.arg(preprocess, c("scale", "center", "none")))
  }

  if (is.character(standardize)) {
    return(match.arg(standardize, c("scale", "center", "none")))
  }

  if (is.logical(standardize) && length(standardize) == 1L) {
    return(if (isTRUE(standardize)) "scale" else "center")
  }

  if (is.numeric(standardize) && length(standardize) == 1L && is.finite(standardize)) {
    return(if (isTRUE(as.logical(standardize))) "scale" else "center")
  }

  if (is.null(standardize)) {
    return("center")
  }

  stop("`standardize` must be logical or one of: 'scale', 'center', 'none'.", call. = FALSE)
}

.preprocess_group_matrix <- function(M, mode) {
  M <- as.matrix(M)
  if (mode == "none") {
    return(M)
  }

  M <- scale(M, center = TRUE, scale = identical(mode, "scale"))
  M <- as.matrix(M)
  M[!is.finite(M)] <- 0
  M
}


# Helper: ADMM-based group sparse solver
solve_group_rrr_admm <- function(X, tilde_Y, groups, lambda, Sx = NULL,
                                 rho = 1, niter = 1e4,
                                 thresh = 0.00001, thresh_0 = 1e-6, verbose = FALSE,
                                 invSx_apply = NULL,
                                 matrix_free_threshold = 4000L,
                                 cg_tol = 1e-6,
                                 cg_maxiter = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(tilde_Y)
  prod_xy <- crossprod(X, tilde_Y) / n

  if (is.null(invSx_apply)) {
    invSx_info <- make_invSx_apply(
      X, Sx, rho, verbose = verbose,
      matrix_free_threshold = matrix_free_threshold,
      cg_tol = cg_tol,
      cg_maxiter = cg_maxiter
    )
    invSx_apply <- invSx_info$apply
    if (verbose) message("ADMM linear solver strategy: ", invSx_info$strategy)
  }
  
  U <- matrix(0, p, q)
  Z <- matrix(0, p, q)
  
  for (i in seq_len(niter)) {
    U_old <- U; Z_old <- Z
    B <- invSx_apply(prod_xy + rho * (Z - U))
    Z <- B + U
    
    for (g in seq_along(groups)) {
      idx <- groups[[g]]
      norm_grp <- sqrt(sum(Z[idx, ]^2))
      if (norm_grp < lambda * sqrt(length(idx)) / rho) {
        Z[idx, ] <- 0
      } else {
        Z[idx, ] <- (1 - (lambda * sqrt(length(idx)) / rho) / norm_grp) * Z[idx, ]
      }
    }
    
    U <- U + B - Z
    
    if (verbose) cat("ADMM iter", i, "Primal: ", norm(Z - B), "Dual:", norm(Z_old - Z), "\n")
    if (max(norm(Z - B), norm(Z_old - Z)) < thresh) break
  }
  
  B_opt <- B
  B_opt[abs(B_opt) < thresh_0] <- 0
  return(B_opt)
}


# Helper: CVXR-based group sparse solver
solve_group_rrr_cvxr <- function(X, tilde_Y, groups, lambda, thresh_0 = 1e-6) {

  if (!nzchar(system.file(package = "CVXR"))) {
    stop("Package 'CVXR' must be installed to use the CVX solver.",
         call. = FALSE)
  }

  p <- ncol(X); q <- ncol(tilde_Y)
  n <- nrow(X)
  cvxr_variable <- getExportedValue("CVXR", "Variable")
  cvxr_norm <- getExportedValue("CVXR", "norm")
  cvxr_minimize <- getExportedValue("CVXR", "Minimize")
  cvxr_sum_squares <- getExportedValue("CVXR", "sum_squares")
  cvxr_problem <- getExportedValue("CVXR", "Problem")
  cvxr_psolve <- getExportedValue("CVXR", "psolve")
  cvxr_value <- getExportedValue("CVXR", "value")

  B <- cvxr_variable(c(p, q))
  
  penalty_exprs <- lapply(groups, function(g) {
    cvxr_norm(B[g, , drop = FALSE], "F")
  })
  penalty <- Reduce(`+`, penalty_exprs)  # or do.call("+", penalty_exprs)
  
  objective <- cvxr_minimize(
    1/n * cvxr_sum_squares(tilde_Y - X %*% B) + lambda * penalty
  )
  
  prob <- cvxr_problem(objective)
  res <- cvxr_psolve(prob)

  B_opt <- cvxr_value(B)
  B_opt[abs(B_opt) < thresh_0] <- 0
  return(B_opt)
}



#' Group-Sparse Canonical Correlation via Reduced-Rank Regression
#'
#' Performs group-sparse reduced-rank regression for CCA using either ADMM or CVXR solvers.
#'
#' @param X Predictor matrix (n x p)
#' @param Y Response matrix (n x q)
#' @param groups List of index vectors defining groups of predictors
#' @param Sx Optional covariance matrix for X; if NULL computed internally
#' @param Sy Optional covariance matrix for Y; if NULL computed internally
#' @param Sxy Optional cross covariance matrix for X and Y; if NULL computed internally
#' @param lambda Regularization parameter
#' @param r Target rank
#' @param standardize Backward-compatible preprocessing flag: TRUE = `"scale"`, FALSE = `"center"`.
#' @param preprocess Preprocessing mode. One of `"scale"` (center + scale), `"center"` (center only), or `"none"`.
#' @param LW_Sy Whether to apply Ledoit-Wolf shrinkage to Sy (default TRUE)
#' @param solver Either "ADMM" or "CVXR". The `"CVXR"` backend requires the
#'   optional package `CVXR`.
#' @param rho ADMM parameter
#' @param niter Maximum number of ADMM iterations
#' @param thresh Convergence threshold for ADMM
#' @param thresh_0 tolerance for declaring entries non-zero
#' @param matrix_free_threshold For ADMM: when both `n` and `p` are at least this value, use a matrix-free conjugate-gradient solve instead of forming a dense linear system.
#' @param cg_tol Relative tolerance for the matrix-free conjugate-gradient solve used in ADMM.
#' @param cg_maxiter Maximum iterations for the matrix-free conjugate-gradient solve. Defaults to `min(p, 1000)`.
#' @param verbose Print diagnostics
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
cca_group_rrr <- function(X, Y, groups, 
                          Sx = NULL, Sy = NULL, Sxy = NULL, 
                          lambda = 0,  r,
                          standardize = FALSE,
                          preprocess = NULL,
                          LW_Sy = TRUE, solver = "ADMM",
                          rho = 1, niter = 1e4, thresh = 1e-4, thresh_0=1e-6,
                          matrix_free_threshold = 4000L,
                          cg_tol = 1e-6,
                          cg_maxiter = NULL,
                          verbose = FALSE) {
  preprocess_mode <- .resolve_group_preprocess_mode(standardize = standardize, preprocess = preprocess)
  X <- .preprocess_group_matrix(X, preprocess_mode)
  Y <- .preprocess_group_matrix(Y, preprocess_mode)

  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  if (n < min(p, q)) warning("Both X and Y are high-dimensional; method may be unstable.")

  if (is.null(Sy)) {
    Sy <- crossprod(Y) / n
    if (LW_Sy) Sy <- as.matrix(corpcor::cov.shrink(Y, verbose=verbose))
  }
  
  svd_Sy <- svd(as.matrix(Sy))
  sqrt_inv_Sy <- svd_Sy$u %*% diag(ifelse(svd_Sy$d > 1e-4, 1 / sqrt(svd_Sy$d), 0)) %*% t(svd_Sy$v)
  tilde_Y <- Y %*% sqrt_inv_Sy
  
  B_opt <- switch(solver,
                  "ADMM" = solve_group_rrr_admm(X, tilde_Y, groups=groups, lambda=lambda,
                                                Sx = Sx,
                                                rho=rho, niter=niter, 
                                                thresh=thresh, thresh_0=thresh_0, verbose=verbose,
                                                matrix_free_threshold = matrix_free_threshold,
                                                cg_tol = cg_tol,
                                                cg_maxiter = cg_maxiter),
                  "CVXR" = solve_group_rrr_cvxr(X, tilde_Y, groups, lambda, thresh_0=thresh_0),
                  stop("Unsupported solver: choose either 'ADMM' or 'CVXR'")
  )

  post <- .postprocess_rrr_fit(
    B_opt, X, Y, sqrt_inv_Sy, r,
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
    B_opt = B_opt
  )
}



# Cross-validated loss for group-penalized CCA
cca_group_rrr_cv_folds <- function(X, Y, groups, Sx = NULL, Sy = NULL, kfolds = 5, 
                                   lambda = 0.01, r = 2, 
                                   preprocess = NULL,
                                   LW_Sy = FALSE, solver = "ADMM", 
                                   rho = 1, niter = 1e4, thresh = 1e-4,
                                   thresh_0 = 1e-6,
                                   matrix_free_threshold = 4000L,
                                   cg_tol = 1e-6,
                                   cg_maxiter = NULL,
                                   verbose = FALSE,
                                   folds = NULL,
                                   return_fold_values = FALSE) {
  preprocess_mode <- .resolve_group_preprocess_mode(standardize = FALSE, preprocess = preprocess)
  X <- .preprocess_group_matrix(X, preprocess_mode)
  Y <- .preprocess_group_matrix(Y, preprocess_mode)
  if (is.null(folds)) {
    folds <- .create_cv_folds(nrow(Y), kfolds)
  }

  rmse <- vapply(seq_along(folds), function(i) {
    n <- nrow(X)
    X_train <- X[-folds[[i]], ]; Y_train <- Y[-folds[[i]], ]
    X_val <- X[folds[[i]], ]; Y_val <- Y[folds[[i]], ]
    n_train <- n - nrow(X_val)
    
    if (!is.null(Sx)) {
      Sx_train <- (n * Sx - crossprod(X_val)) / n_train
    } else {
      Sx_train <- NULL
    }

    if (!is.null(Sy)) {
      Sy_train <- (n * Sy - crossprod(Y_val)) / n_train
    } else {
      Sy_train <- NULL
    }

    tryCatch({
      fit <- cca_group_rrr(X_train, Y_train, groups, Sx = Sx_train, Sy = Sy_train,
                           lambda = lambda, r = r, 
                           standardize = FALSE, preprocess = "none",
                           LW_Sy = LW_Sy, solver = solver,
                           rho = rho, niter = niter, thresh = thresh, thresh_0 = thresh_0,
                           matrix_free_threshold = matrix_free_threshold,
                           cg_tol = cg_tol,
                           cg_maxiter = cg_maxiter,
                           verbose = FALSE)
      mean((X_val %*% fit$U - Y_val %*% fit$V)^2)
    }, error = function(e) {
      message("Error in fold ", i, ": ", conditionMessage(e))
      NA_real_
    })
  }, numeric(1))
  
  if (return_fold_values) return(rmse)
  if (all(is.na(rmse))) return(1e8)
  mean(rmse, na.rm = TRUE)
}


#' Group-Sparse Canonical Correlation via Reduced-Rank Regression with CV
#'
#' Performs group-sparse reduced-rank regression for CCA using either ADMM or CVXR solvers.
#'
#' @param X Predictor matrix (n x p)
#' @param Y Response matrix (n x q)
#' @param groups List of index vectors defining groups of predictors
#' @param lambdas Grid of regularization parameters to try out
#' @param r Target rank
#' @param kfolds Nb of folds for the CV procedure
#' @param parallelize Whether to use parallel processing (default is FALSE)
#' @param standardize Backward-compatible preprocessing flag: TRUE = `"scale"`, FALSE = `"center"`.
#' @param preprocess Preprocessing mode. One of `"scale"` (center + scale), `"center"` (center only), or `"none"`.
#' @param LW_Sy Whether to apply Ledoit-Wolf shrinkage to Sy (default TRUE)
#' @param solver Either "ADMM" or "CVXR". The `"CVXR"` backend requires the
#'   optional package `CVXR`.
#' @param rho ADMM parameter
#' @param niter Maximum number of ADMM iterations
#' @param thresh Convergence threshold for ADMM
#' @param verbose Print diagnostics
#' @param thresh_0 tolerance for declaring entries non-zero
#' @param matrix_free_threshold For ADMM: when both `n` and `p` are at least this value, use a matrix-free conjugate-gradient solve instead of forming a dense linear system.
#' @param cg_tol Relative tolerance for the matrix-free conjugate-gradient solve used in ADMM.
#' @param cg_maxiter Maximum iterations for the matrix-free conjugate-gradient solve. Defaults to `min(p, 1000)`.
#' @param nb_cores Number of cores to use for parallelization (default is all available cores minus 1)
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
cca_group_rrr_cv <- function(X, Y, groups, r = 2, 
                             lambdas = 10^seq(-3, 1.5, length.out = 10),
                             kfolds = 5, parallelize = FALSE, standardize = FALSE,
                             preprocess = NULL,
                             LW_Sy = TRUE, solver = "ADMM", rho = 1,
                             thresh_0 = 0,
                             niter = 1e4, thresh = 1e-4,
                             matrix_free_threshold = 4000L,
                             cg_tol = 1e-6,
                             cg_maxiter = NULL,
                             verbose = FALSE,
                             nb_cores = NULL) {
  preprocess_mode <- .resolve_group_preprocess_mode(standardize = standardize, preprocess = preprocess)
  
  if (nrow(X) < min(ncol(X), ncol(Y))) {
    warning("Both X and Y are high-dimensional; method may be unstable.")
  }
  
  # Apply preprocessing once before CV. Each fold reuses transformed matrices.
  X <- .preprocess_group_matrix(X, preprocess_mode)
  Y <- .preprocess_group_matrix(Y, preprocess_mode)
  
  n <- nrow(X)
  Sx <- NULL
  Sy <- if (LW_Sy) NULL else crossprod(Y) / n

  # Rebind the current source definitions into a local environment so PSOCK
  # workers use the updated functions rather than any stale installed copy.
  .rrr_cg_solve <- .rrr_cg_solve
  .rrr_chol_jitter <- .rrr_chol_jitter
  .rrr_invsqrt_psd <- .rrr_invsqrt_psd
  .rrr_whiten_factor <- .rrr_whiten_factor
  .rrr_safe_rank_svd <- .rrr_safe_rank_svd
  make_invSx_apply <- make_invSx_apply
  .postprocess_rrr_fit <- .postprocess_rrr_fit
  .resolve_group_preprocess_mode <- .resolve_group_preprocess_mode
  .preprocess_group_matrix <- .preprocess_group_matrix
  .create_cv_folds <- .create_cv_folds
  solve_group_rrr_admm <- solve_group_rrr_admm
  solve_group_rrr_cvxr <- solve_group_rrr_cvxr
  cca_group_rrr <- cca_group_rrr
  cca_group_rrr_cv_folds <- cca_group_rrr_cv_folds
  safe_mean_na <- .safe_mean_na
  safe_sd_na <- .safe_sd_na
  environment(make_invSx_apply) <- environment()
  environment(.postprocess_rrr_fit) <- environment()
  environment(solve_group_rrr_admm) <- environment()
  environment(solve_group_rrr_cvxr) <- environment()
  environment(cca_group_rrr) <- environment()
  environment(cca_group_rrr_cv_folds) <- environment()
  
  run_cv <- function(lambda) {
    fold_values <- cca_group_rrr_cv_folds(
      X, Y, groups, Sx = Sx, Sy = Sy, kfolds = kfolds,
      lambda = lambda, r = r,
      preprocess = "none",
      LW_Sy = LW_Sy, solver = solver,
      rho = rho, niter = niter, thresh = thresh, thresh_0 = thresh_0,
      matrix_free_threshold = matrix_free_threshold,
      cg_tol = cg_tol,
      cg_maxiter = cg_maxiter,
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

  folds <- .create_cv_folds(n, kfolds)
  
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
      cl <- setup_parallel_backend(nb_cores)
      
      if (!is.null(cl)) {
        # If the cluster was created successfully, register it and plan to stop it
        if (verbose){
          cat(crayon::green("Parallel backend successfully registered.\n"))
        }

        initialize_parallel_workers(cl, verbose = verbose)
        doParallel::registerDoParallel(cl)
        on.exit(cleanup_parallel_backend(cl), add = TRUE)
      } else {
        # If setup_parallel_backend returned NULL, print a warning and proceed serially
        warning("All parallel setup attempts failed. Proceeding in serial mode.", immediate. = TRUE)
        parallelize <- FALSE # Ensure %dopar% runs serially
    }
   }

   if (parallelize) {


    parallel_packages <- if (solver == "CVXR") c("Matrix", "CVXR") else c("Matrix")
    results_list <- foreach(
      lambda = lambdas,
      .packages = parallel_packages,
      .export = c(
        "cca_group_rrr_cv_folds", "cca_group_rrr",
        "solve_group_rrr_admm", "solve_group_rrr_cvxr",
        "make_invSx_apply", ".rrr_cg_solve", ".rrr_chol_jitter",
        ".rrr_invsqrt_psd", ".rrr_whiten_factor", ".rrr_safe_rank_svd",
        ".postprocess_rrr_fit", ".resolve_group_preprocess_mode",
        ".preprocess_group_matrix", ".create_cv_folds",
        ".safe_mean_na", ".safe_sd_na"
      )
    ) %dopar% run_cv(lambda)
  } else {
    results_list <- lapply(lambdas, run_cv)
  }
  
  cv_summary <- dplyr::bind_rows(lapply(results_list, `[[`, "summary"))
  cv_folds <- dplyr::bind_rows(lapply(results_list, `[[`, "fold"))

  cv_summary$rmse[is.na(cv_summary$rmse) | cv_summary$rmse == 0] <- 1e8
  cv_summary <- cv_summary %>% dplyr::filter(rmse > 1e-5)
  
  opt_lambda <- cv_summary$lambda[which.min(cv_summary$rmse)]
  if (is.na(opt_lambda)) opt_lambda <- 0.1
  
  final <- cca_group_rrr(X, Y, groups, Sx = Sx, Sy = Sy, lambda = opt_lambda,
                         r = r, preprocess = "none", thresh = thresh, thresh_0 = thresh_0,
                         LW_Sy = LW_Sy, solver = solver, rho = rho, niter = niter,
                         matrix_free_threshold = matrix_free_threshold,
                         cg_tol = cg_tol,
                         cg_maxiter = cg_maxiter,
                         verbose=verbose)
  
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
