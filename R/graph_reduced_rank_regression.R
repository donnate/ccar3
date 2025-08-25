library(magrittr)
library(tidyr)
library(SMUT) # Required for eigenMapMatMult
library(pracma) # Required for pinv
library(caret) # Required for createFolds


matmul <- function(A, B) {
  if (requireNamespace("SMUT", quietly = TRUE)) {
    # Use the fast C++ multiplication from SMUT
    SMUT::eigenMapMatMult(A, B)
  } else {
    # Fallback to base R multiplication
    A %*% B
  }
}



cca_graph_rrr_cv_folds <- function(X, Y, Gamma,
                                Sx = NULL, Sy = NULL, kfolds = 5,
                                lambda = 0.01, r = 2, 
                                standardize = FALSE, 
                                LW_Sy = FALSE, rho = 10,
                                niter = 1e4, thresh = 1e-4,
                                thresh_0 = 1e-6,
                                Gamma_dagger = NULL) {
  
  folds <- caret::createFolds(1:nrow(Y), k = kfolds, list = TRUE)

  rmse <- foreach(i = seq_along(folds), .combine = c, .packages = c('CVXR', 'Matrix')) %do% {
    n_full <- nrow(X)
    X_train <- X[-folds[[i]], ]
    Y_train <- Y[-folds[[i]], ]
    X_val <- X[folds[[i]], ]
    Y_val <- Y[folds[[i]], ]

   # We don't even need to create X_train and Y_train explicitly for the covariance calculation
    n_train <- n_full - nrow(X_val)

    # --- DOWNDATE TRICK ---
    Sx_train <- (n_full * Sx - crossprod(X_val)) / n_train
    Sy_train <- (n_full * Sy - crossprod(Y_val)) / n_train


    
    tryCatch({
      fit <- cca_graph_rrr(X_train, Y_train, Gamma,
                           Sx_train, Sy_train, lambda, r=r,
                           standardize=standardize, 
                           LW_Sy=LW_Sy, rho=rho, niter=niter, thresh=thresh,
                           thresh_0 = thresh_0,
                           Gamma_dagger = Gamma_dagger)

      mean((X_val %*% fit$U - Y_val %*% fit$V)^2)
    }, error = function(e) {
      message("Error in fold ", i, ": ", conditionMessage(e))
      return(NA)
    })
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
#' @param standardize Whether to center and scale X and Y (default FALSE = center only)
#' @param LW_Sy Whether to apply Ledoit-Wolf shrinkage to Sy
#' @param rho ADMM penalty parameter
#' @param niter Maximum number of ADMM iterations
#' @param thresh Convergence threshold for ADMM
#' @param thresh_0 Threshold for small values in the coefficient matrix (default 1e-6)
#' @param verbose Whether to print diagnostic output
#' @param Gamma_dagger Optional pseudoinverse of Gamma (computed if NULL)
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
cca_graph_rrr_cv <- function(X, Y, Gamma,
                             r = 2, 
                             lambdas = 10^seq(-3, 1.5, length.out = 10),
                             kfolds = 5, parallelize = FALSE,
                             standardize = TRUE, LW_Sy = FALSE,
                             rho = 10, niter = 1e4, thresh = 1e-4,
                             thresh_0 = 1e-6, verbose = FALSE,
                             Gamma_dagger = NULL, nb_cores= NULL) {

  if (nrow(X) < min(ncol(X), ncol(Y))) {
    warning("Both X and Y are high dimensional; method may fail.")
  }

  if (standardize) {
    X <- scale(X)
    Y <- scale(Y)
  } 

  n <- nrow(X) # Define n
  Sx <- crossprod(X) / n # Pre-compute full Sx
  Sy <- crossprod(Y) / n # Pre-compute full Sy



  run_cv <- function(lambda) {
    rmse <- cca_graph_rrr_cv_folds(X, Y, Gamma,
                                 Sx = NULL, Sy = NULL,
                                 kfolds = kfolds,
                                 lambda = lambda, r = r,
                                 standardize = FALSE,
                                 LW_Sy = LW_Sy, rho = rho,
                                 niter = niter, thresh = thresh,
                                 thresh_0 = thresh_0,
                                 Gamma_dagger = Gamma_dagger)
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
    results <-  foreach(lambda = lambdas, .combine = rbind, .packages = c('CVXR', 'Matrix')) %dopar% run_cv(lambda)
  } else {
    results <- purrr::map_dfr(lambdas, run_cv)
  }

  results$rmse[is.na(results$rmse) | results$rmse == 0] <- 1e8
  results <- results %>% dplyr::filter(rmse > 1e-5)

  opt_lambda <- results$lambda[which.min(results$rmse)]

  final <- cca_graph_rrr(X, Y, Gamma, Sx = NULL, Sy = NULL,
                         lambda = opt_lambda, r = r,
                         standardize = FALSE,
                         LW_Sy = LW_Sy, rho = rho, niter = niter,
                         thresh = thresh, Gamma_dagger = Gamma_dagger, thresh_0=thresh_0)

  list(
    U = final$U,
    V = final$V,
    lambda = opt_lambda,
    #resultsx = results,
    rmse = results$rmse,
    cor = sapply(1:r, function(i) cov(X %*% final$U[, i], Y %*% final$V[, i]))
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
#' @param standardize Whether to center and scale X and Y (default FALSE = center only)
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
                          LW_Sy = TRUE, rho = 10,
                          niter = 1e4, thresh = 1e-4,
                          thresh_0 = 1e-6,
                          verbose = FALSE, Gamma_dagger = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)

if (q > n) {
  stop("The X/Y swapping heuristic (for q > n) is not compatible with the graph constraint Gamma. Swap X and Y in your arguments.")
}

  if (standardize) {
    X <- scale(X)
    Y <- scale(Y)
  } 
  # else {
  #   X <- scale(X, scale = FALSE)
  #   Y <- scale(Y, scale = FALSE)
  # }

  if (n < min(p, q)) {
    warning("Both X and Y are high dimensional; method may be unstable.")
  }

  # Covariance matrices
  if (is.null(Sx)) Sx <- t(X) %*% X / n
  if (is.null(Sy)) {
    Sy <- t(Y) %*% Y / n
    if (LW_Sy) Sy <- as.matrix(corpcor::cov.shrink(Y))
  }

  # Inverse square root of Sy
  svd_Sy <- svd(Sy)
  sqrt_inv_Sy <- svd_Sy$u %*% diag(ifelse(svd_Sy$d > 1e-4, 1 / sqrt(svd_Sy$d), 0)) %*% t(svd_Sy$v)
  tilde_Y <- Y %*% sqrt_inv_Sy

  # Graph penalty
  if (is.null(Gamma_dagger)) Gamma_dagger <- pracma::pinv(Gamma)
  Pi <- diag(p) - Gamma_dagger %*% Gamma
  XPi <- X %*% Pi
  XGamma_dagger <- X %*% Gamma_dagger

  # Remove projection on Pi
  Projection <- pracma::pinv(t(XPi) %*% XPi) %*% t(XPi) %*% tilde_Y
  new_Ytilde <- tilde_Y - XPi %*% Projection

  # Add structure penalty
  Sx_tot <-  Sx

  # ADMM initialization
  new_p <- ncol(XGamma_dagger)
  prod_xy <- t(XGamma_dagger) %*% new_Ytilde / n
  invSx <- solve(t(XGamma_dagger) %*% XGamma_dagger / n + rho * diag(new_p))
  U <- matrix(0, new_p, q)
  Z <- matrix(0, new_p, q)

  for (i in seq_len(niter)) {
    U_old <- U; 
    Z_old <- Z

    B <- invSx %*% (prod_xy + rho * (Z - U))
    Z <- B + U

    # Group soft-thresholding
    norms <- sqrt(rowSums(Z^2))
    shrink_factors <- pmax(0, 1 - lambda / (rho * norms))
    shrink_factors[is.nan(shrink_factors)] <- 0
    Z <- Z * shrink_factors

    U <- U + B - Z

    # Convergence check
    if (verbose) {
      cat("ADMM iter:", i, "Primal:", norm(Z - B), "Dual:", norm(Z_old - Z), "\n")
    }
    if (max(norm(Z - B, "F") / sqrt(new_p), norm(Z_old - Z, "F") / sqrt(new_p)) < thresh) break
  }

  # Reconstruct full coefficient matrix
  B_opt <- Gamma_dagger %*% B + Pi %*% Projection
  B_opt[abs(B_opt) < thresh_0] <- 0

  # Final CCA step: SVD of transformed coefficient matrix
  svd_Sx <- svd(Sx_tot)
  sqrt_inv_Sx <- svd_Sx$u %*% diag(ifelse(svd_Sx$d > 1e-4, 1 / sqrt(svd_Sx$d), 0)) %*% t(svd_Sx$u)
  sqrt_Sx <- svd_Sx$u %*% diag(ifelse(svd_Sx$d > 1e-4, sqrt(svd_Sx$d), 0)) %*% t(svd_Sx$u)

  sol <- svd(sqrt_Sx %*% B_opt)
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
