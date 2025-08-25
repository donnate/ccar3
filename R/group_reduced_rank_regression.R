#Required libraries
library(dplyr)
library(magrittr)
library(foreach)


# Helper: ADMM-based group sparse solver
solve_group_rrr_admm <- function(X, tilde_Y, Sx, groups, lambda, rho=1, niter=1e4, 
                                 thresh=0.00001, thresh_0=1e-6, verbose = FALSE) {
  p <- ncol(X); q <- ncol(tilde_Y)
  Sx_tot <- Sx
  prod_xy <- matmul(t(X), tilde_Y)/ nrow(X)
  invSx <- solve(Sx_tot + rho * diag(p))
  
  U <- matrix(0, p, q)
  Z <- matrix(0, p, q)
  
  for (i in seq_len(niter)) {
    U_old <- U; Z_old <- Z
    B <- invSx %*% (prod_xy + rho * (Z - U))
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

  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop("Package 'CVXR' must be installed to use the CVX solver.",
         call. = FALSE)
  }

  p <- ncol(X); q <- ncol(tilde_Y)
  n <- nrow(X)
  B <- CVXR::Variable(p, q)
  
  penalty_exprs <- lapply(groups, function(g) {
    CVXR::norm(B[g, , drop = FALSE], "F")
  })
  penalty <- Reduce(`+`, penalty_exprs)  # or do.call("+", penalty_exprs)
  
  objective <- CVXR::Minimize(
    1/n *  CVXR::sum_squares(tilde_Y - X %*% B) + lambda * penalty
  )
  
  prob <- CVXR::Problem(objective)
  res <- CVXR::solve(prob)
  
  B_opt <- res$getValue(B)
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
#' @param standardize Whether to scale variables
#' @param LW_Sy Whether to apply Ledoit-Wolf shrinkage to Sy (default TRUE)
#' @param solver Either "ADMM" or "CVXR"
#' @param rho ADMM parameter
#' @param niter Maximum number of ADMM iterations
#' @param thresh Convergence threshold for ADMM
#' @param thresh_0 tolerance for declaring entries non-zero
#' @param verbose Print diagnostics
#'
#' @return A list with elements:
#' \describe{
#'   \item{U}{Canonical direction matrix for X (p x r)}
#'   \item{V}{Canonical direction matrix for Y (q x r)}
#'   \item{cor}{Canonical covariances}
#'   \item{loss}{The prediction error 1/n * \| XU - YV\|^2}
#' }
#' @export
cca_group_rrr <- function(X, Y, groups, 
                          Sx = NULL, Sy = NULL, Sxy = NULL, 
                          lambda = 0,  r,
                          standardize = FALSE, 
                          LW_Sy = TRUE, solver = "ADMM",
                          rho = 1, niter = 1e4, thresh = 1e-4, thresh_0=1e-6, verbose = FALSE) {
  
  n <- nrow(X); p <- ncol(X); q <- ncol(Y)
  
  if (standardize) {
    X <- scale(X); Y <- scale(Y)
  }
  #  else {
  #   X <- scale(X, scale = FALSE); Y <- scale(Y, scale = FALSE)
  # }
  
  if (n < min(p, q)) warning("Both X and Y are high-dimensional; method may be unstable.")
  if (is.null(Sx)) Sx <- t(X) %*% X / n
  if (is.null(Sy)) {
    Sy <- t(Y) %*% Y / n
    if (LW_Sy) Sy <- as.matrix(corpcor::cov.shrink(Y, verbose=verbose))
  }
  
  svd_Sy <- svd(Sy)
  sqrt_inv_Sy <- svd_Sy$u %*% diag(ifelse(svd_Sy$d > 1e-4, 1 / sqrt(svd_Sy$d), 0)) %*% t(svd_Sy$v)
  tilde_Y <- Y %*% sqrt_inv_Sy
  
  B_opt <- switch(solver,
                  "ADMM" = solve_group_rrr_admm(X, tilde_Y, Sx, groups=groups, lambda=lambda, rho=rho, niter=niter, 
                                                thresh=thresh, thresh_0=thresh_0, verbose=verbose),
                  "CVXR" = solve_group_rrr_cvxr(X, tilde_Y, groups, lambda, thresh_0=thresh_0),
                  stop("Unsupported solver: choose either 'ADMM' or 'CVXR'")
  )
  
  svd_Sx <- svd(Sx)
  sqrt_Sx <- svd_Sx$u %*% diag(ifelse(svd_Sx$d > 1e-4, sqrt(svd_Sx$d), 0)) %*% t(svd_Sx$u)
  sqrt_inv_Sx <- svd_Sx$u %*% diag(ifelse(svd_Sx$d > 1e-4, 1 / sqrt(svd_Sx$d), 0)) %*% t(svd_Sx$u)
  
  sol <- svd(sqrt_Sx %*% B_opt)
  V <- sqrt_inv_Sy %*% sol$v[, 1:r]
  inv_D <- diag(ifelse(sol$d[1:r] > 1e-4, 1 / sol$d[1:r], 0))
  U <- B_opt %*% sol$v[, 1:r] %*% inv_D
  
  loss <- mean((Y %*% V - X %*% U)^2)
  
  list(
    U = U,
    V = V,
    loss = loss,
    cor = sapply(1:r, function(i) cov(X %*% U[, i], Y %*% V[, i]))
  )
}



# Cross-validated loss for group-penalized CCA
cca_group_rrr_cv_folds <- function(X, Y, groups, Sx = NULL, Sy = NULL, kfolds = 5, 
                                   lambda = 0.01, r = 2, 
                                   standardize = FALSE, LW_Sy = FALSE, solver = "ADMM", 
                                   rho = 1, niter = 1e4, thresh = 1e-4,
                                   thresh_0=1e-6, 
                                   verbose = FALSE) {
  
  folds <- caret::createFolds(1:nrow(Y), k = kfolds, list = TRUE)
  if (solver == "CVXR"){
    if (!requireNamespace("CVXR", quietly = TRUE)) {
       stop("Package 'CVXR' must be installed to use the CVXR solver.", call. = FALSE)
    }
    packages_list <- c("CVXR", "Matrix")
  } else {
    packages_list <- c()
  }
  rmse <- foreach(i = seq_along(folds), .combine = c, .packages = packages_list) %do% {
    n <- nrow(X)
    X_train <- X[-folds[[i]], ]; Y_train <- Y[-folds[[i]], ]
    X_val <- X[folds[[i]], ]; Y_val <- Y[folds[[i]], ]
    n_train <- n - nrow(X_val)
    
    if (is.null(Sx) == FALSE) {
      Sx_train <- (n * Sx - crossprod(X_val)) / n_train
    } else {
      Sx_train <- NULL
    }
    
    
    tryCatch({
      fit <- cca_group_rrr(X_train, Y_train, groups, Sx = Sx_train, Sy = NULL,
                           lambda = lambda, r = r, 
                           standardize = FALSE, LW_Sy = LW_Sy, solver = solver,
                           rho = rho, niter = niter, thresh = thresh, thresh_0=thresh_0,verbose = FALSE)
      mean((X_val %*% fit$U - Y_val %*% fit$V)^2)
    }, error = function(e) {
      message("Error in fold ", i, ": ", conditionMessage(e))
      return(NA)
    })
  }
  
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
#' @param standardize Whether to scale variables
#' @param LW_Sy Whether to apply Ledoit-Wolf shrinkage to Sy (default TRUE)
#' @param solver Either "ADMM" or "CVXR"
#' @param rho ADMM parameter
#' @param niter Maximum number of ADMM iterations
#' @param thresh Convergence threshold for ADMM
#' @param verbose Print diagnostics
#' @param thresh_0 tolerance for declaring entries non-zero
#' @param nb_cores Number of cores to use for parallelization (default is all available cores minus 1)
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
cca_group_rrr_cv <- function(X, Y, groups, r = 2, 
                             lambdas = 10^seq(-3, 1.5, length.out = 10),
                             kfolds = 5, parallelize = FALSE, standardize = FALSE,
                             LW_Sy = TRUE, solver = "ADMM", rho = 1,
                             thresh_0 = 1e-6,
                             niter = 1e4, thresh = 1e-4, verbose = FALSE,
                             nb_cores = NULL) {
  
  if (nrow(X) < min(ncol(X), ncol(Y))) {
    warning("Both X and Y are high-dimensional; method may be unstable.")
  }
  
  X <- if (standardize) scale(X) else X #scale(X, scale = FALSE)
  Y <- if (standardize) scale(Y) else Y #scale(Y, scale = FALSE)
  
  Sx = matmul(t(X), X) / nrow(X)
  #Sy <- if (LW_Sy) as.matrix(corpcor::cov.shrink(Y, verbose=verbose )) else t(Y) %*% Y / nrow(Y)
  
  run_cv <- function(lambda) {
    rmse <- cca_group_rrr_cv_folds(X, Y, groups, Sx = Sx, Sy = NULL, kfolds = kfolds,
                                   lambda = lambda, r = r, 
                                   standardize = FALSE, LW_Sy = LW_Sy, solver = solver,
                                   rho = rho, niter = niter, thresh = thresh, thresh_0 = thresh_0)
    data.frame(lambda = lambda, rmse = rmse)
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
      cl <- setup_parallel_backend(nb_cores)
      
      if (!is.null(cl)) {
        # If the cluster was created successfully, register it and plan to stop it
        cat(crayon::green("Parallel backend successfully registered.\n"))
        doParallel::registerDoParallel(cl)
        on.exit(parallel::stopCluster(cl), add = TRUE)
      } else {
        # If setup_parallel_backend returned NULL, print a warning and proceed serially
        warning("All parallel setup attempts failed. Proceeding in serial mode.", immediate. = TRUE)
        parallelize <- FALSE # Ensure %dopar% runs serially
    }
   }

   if (parallelize) {


    results <- foreach(lambda = lambdas, .combine = rbind, .packages = c( "Matrix")) %dopar% run_cv(lambda)
  } else {
    results <- purrr::map_dfr(lambdas, run_cv)
  }
  
  results$rmse[is.na(results$rmse) | results$rmse == 0] <- 1e8
  results <- results %>% dplyr::filter(rmse > 1e-5)
  
  opt_lambda <- results$lambda[which.min(results$rmse)]
  if (is.na(opt_lambda)) opt_lambda <- 0.1
  
  final <- cca_group_rrr(X, Y, groups, Sx = Sx, Sy = NULL, lambda = opt_lambda,
                         r = r, standardize = FALSE, thresh = thresh, thresh_0 = thresh_0,
                         LW_Sy = LW_Sy, solver = solver, rho = rho, niter = niter, verbose=verbose)
  
  list(
    U = final$U,
    V = final$V,
    lambda = opt_lambda,
    #resultsx = results,
    rmse = results$rmse,
    cor = sapply(1:r, function(i) cov(X %*% final$U[, i], Y %*% final$V[, i]))
  )
}


