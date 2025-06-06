library(dplyr)
library(tidyr)
library(Matrix)
library(glmnet)
library(gglasso)
library(rrpack)
library(foreach)
library(doParallel)
library(CVXR)
library(caret) # Required for createFolds
library(parallel) # Required for detectCores


cca_graph_rrr_cv_folds <- function(X, Y, Gamma,
                                Sx = NULL, Sy = NULL, kfolds = 5,
                                lambda = 0.01, r = 2, Kx = NULL,
                                do.scale = FALSE, lambda_Kx = 0,
                                LW_Sy = FALSE, rho = 10,
                                niter = 1e4, thresh = 1e-4,
                                Gamma_dagger = NULL) {
  
  folds <- createFolds(1:nrow(Y), k = kfolds, list = TRUE)
  no_cores <- detectCores() - 2
  registerDoParallel(cores = no_cores)

  rmse <- foreach(i = seq_along(folds), .combine = c, .packages = c('CVXR', 'Matrix')) %dopar% {
    X_train <- X[-folds[[i]], ]
    Y_train <- Y[-folds[[i]], ]
    X_val <- X[folds[[i]], ]
    Y_val <- Y[folds[[i]], ]
    
    tryCatch({
      fit <- cca_graph_rrr(X_train, Y_train, Gamma,
                           Sx, Sy, lambda, Kx, r,
                           do.scale, lambda_Kx,
                           LW_Sy, rho, niter, thresh,
                           Gamma_dagger)

      mean((X_val %*% fit$U - Y_val %*% fit$V)^2)
    }, error = function(e) {
      message("An error occurred in fold ", i)
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
#' @param Sx Optional covariance matrix for X. If NULL, computed as t(X) %*% X / n
#' @param Sy Optional covariance matrix for Y. If NULL, computed similarly; optionally shrunk via Ledoit-Wolf
#' @param kfolds Number of folds for cross-validation
#' @param parallelize Whether to parallelize cross-validation
#' @param Sxy Optional cross-covariance matrix (not currently used)
#' @param param_lambda Grid of regularization parameters to test for sparsity
#' @param Kx Optional penalty matrix to be added to Sx
#' @param r Target rank
#' @param do.scale Whether to center and scale X and Y (default FALSE = center only)
#' @param lambda_Kx Regularization strength for Kx
#' @param LW_Sy Whether to apply Ledoit-Wolf shrinkage to Sy
#' @param rho ADMM penalty parameter
#' @param niter Maximum number of ADMM iterations
#' @param thresh Convergence threshold for ADMM
#' @param verbose Whether to print diagnostic output
#' @param Gamma_dagger Optional pseudoinverse of Gamma (computed if NULL)
#'
#' @return A list with elements:
#' \describe{
#'   \item{U}{Canonical direction matrix for X (p x r)}
#'   \item{V}{Canonical direction matrix for Y (q x r)}
#'   \item{lambda}{Optimal regularisation parameter lambda chosen by CV}
#'   \item{rmse}{Mean squared error of prediction (as computed in the CV)}
#'   \item{cor}{Canonical covariances}
#' }
#' @export
cca_graph_rrr_cv <- function(X, Y, Gamma,
                             r = 2, Kx = NULL, lambda_Kx = 0,
                             param_lambda = 10^seq(-3, 1.5, length.out = 10),
                             kfolds = 5, parallelize = FALSE,
                             do.scale = FALSE, LW_Sy = FALSE,
                             rho = 10, niter = 1e4, thresh = 1e-4,
                             Gamma_dagger = NULL) {

  if (nrow(X) < min(ncol(X), ncol(Y))) {
    warning("Both X and Y are high dimensional; method may fail.")
  }

  run_cv <- function(lambda) {
    rmse <- cca_graph_rrr_cv_folds(X, Y, Gamma,
                                 Sx = NULL, Sy = NULL,
                                 kfolds = kfolds,
                                 lambda = lambda, r = r,
                                 Kx = Kx, lambda_Kx = lambda_Kx,
                                 do.scale = do.scale,
                                 LW_Sy = LW_Sy, rho = rho,
                                 niter = niter, thresh = thresh,
                                 Gamma_dagger = Gamma_dagger)
    data.frame(lambda = lambda, rmse = rmse)
  }

  results <- if (parallelize) {
    no_cores <- detectCores() - 5
    registerDoParallel(cores = no_cores)
    foreach(lambda = param_lambda, .combine = rbind, .packages = c('CVXR', 'Matrix')) %dopar% run_cv(lambda)
  } else {
    map_dfr(param_lambda, run_cv)
  }

  results$rmse[is.na(results$rmse) | results$rmse == 0] <- 1e8
  results <- results %>% filter(rmse > 1e-5)

  opt_lambda <- results$lambda[which.min(results$rmse)]

  final <- cca_graph_rrr(X, Y, Gamma, Sx = NULL, Sy = NULL,
                         lambda = opt_lambda, Kx = Kx, r = r,
                         lambda_Kx = lambda_Kx, do.scale = do.scale,
                         LW_Sy = LW_Sy, rho = rho, niter = niter,
                         thresh = thresh, Gamma_dagger = Gamma_dagger)

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
#' @param Kx Optional penalty matrix to be added to Sx
#' @param r Target rank
#' @param do.scale Whether to center and scale X and Y (default FALSE = center only)
#' @param lambda_Kx Regularization strength for Kx
#' @param LW_Sy Whether to apply Ledoit-Wolf shrinkage to Sy
#' @param rho ADMM penalty parameter
#' @param niter Maximum number of ADMM iterations
#' @param thresh Convergence threshold for ADMM
#' @param verbose Whether to print diagnostic output
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
                          lambda = 0, Kx = NULL, r,
                          do.scale = FALSE, lambda_Kx = 0,
                          LW_Sy = FALSE, rho = 10,
                          niter = 1e4, thresh = 1e-4,
                          verbose = FALSE, Gamma_dagger = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)

  if (q > n) {
    tmp <- X; X <- Y; Y <- tmp
  }

  if (do.scale) {
    X <- scale(X)
    Y <- scale(Y)
  } else {
    X <- scale(X, scale = FALSE)
    Y <- scale(Y, scale = FALSE)
  }

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
  if (is.null(Gamma_dagger)) Gamma_dagger <- pinv(Gamma)
  Pi <- diag(p) - Gamma_dagger %*% Gamma
  XPi <- X %*% Pi
  XGamma_dagger <- X %*% Gamma_dagger

  # Remove projection on Pi
  Projection <- pinv(t(XPi) %*% XPi) %*% t(XPi) %*% tilde_Y
  new_Ytilde <- tilde_Y - XPi %*% Projection

  # Add structure penalty
  Sx_tot <- if (!is.null(Kx)) Sx + lambda_Kx * as.matrix(Kx) else Sx

  # ADMM initialization
  new_p <- ncol(XGamma_dagger)
  prod_xy <- t(XGamma_dagger) %*% new_Ytilde / n
  invSx <- solve(t(XGamma_dagger) %*% XGamma_dagger / n + rho * diag(new_p))
  U <- matrix(0, new_p, q)
  Z <- matrix(0, new_p, q)

  for (i in seq_len(niter)) {
    U_old <- U; Z_old <- Z

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
  B_opt[abs(B_opt) < 1e-5] <- 0

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
