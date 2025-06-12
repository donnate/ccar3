# Required libraries
library(tidyverse)
library(Matrix)
library(glmnet)
library(gglasso)
library(rrpack)
library(foreach)
library(doParallel)

# Helper functions
compute_sqrt_inv <- function(S, threshold = 1e-4) {
  svd_S <- svd(S)
  svd_S$u %*% diag(sapply(svd_S$d, function(x) ifelse(x > threshold, 1 / sqrt(x), 0))) %*% t(svd_S$v)
}

compute_sqrt <- function(S, threshold = 1e-4) {
  svd_S <- svd(S)
  svd_S$u %*% diag(sapply(svd_S$d, function(x) ifelse(x > threshold, sqrt(x), 0))) %*% t(svd_S$u)
}

# Helper: ADMM-based group sparse solver
solve_rrr_admm <- function(X, tilde_Y, Sx, lambda, rho=1, niter=10, thresh, verbose = FALSE) {
  p <- ncol(X); q <- ncol(tilde_Y)
  n <- nrow(X)
  Sx_tot <- Sx
  
  prod_xy <- t(X) %*% tilde_Y / nrow(X)
  invSx <- solve(Sx_tot + rho * diag(p))
  
  U <- Z <- matrix(0, p, q)

  prod_xy <- crossprod(X, tilde_Y) / n
  invSx <- solve(Sx_tot + rho * diag(p))
  for (i in seq_len(niter)) {
    U_old <- U; Z_old <- Z
    B <- invSx %*% (prod_xy + rho * (Z - U))
    Z <- B + U
    norm_col <- sqrt(rowSums(Z^2))
    shrinkage <- pmax(0, 1 - (lambda / rho) / norm_col)
    shrinkage[is.nan(shrinkage)] <- 0
    Z <- sweep(Z, 1, shrinkage, '*')
    U <- U + B - Z
    if (verbose) cat("ADMM iter", i, "Primal: ", norm(Z - B), "Dual:", norm(Z_old - Z), "\n")
    if (max(c(norm(Z - B) / sqrt(p), norm(Z_old - Z) / sqrt(p))) < thresh) break
  }
  B_opt <- B
    
  B_opt[abs(B_opt) < 1e-5] <- 0
  return(B_opt)
}


# Helper: CVXR-based group sparse solver
solve_rrr_cvxr <- function(X, tilde_Y, lambda) {
  p <- ncol(X); q <- ncol(tilde_Y)
  n <- nrow(X)
  B <- CVXR::Variable(p, q)

  objective <- CVXR::Minimize(1 / n * CVXR::sum_squares(tilde_Y - X %*% B) + lambda * sum(CVXR::norm2(B, axis = 1)))
  result <- CVXR::solve(CVXR::Problem(objective))
  B_opt <- result$getValue(B)
  
  B_opt[abs(B_opt) < 1e-5] <- 0
  return(B_opt)
}


#' Canonical Correlation Analysis via Reduced Rank Regression (RRR)
#'
#' Estimates canonical directions using various RRR solvers and penalties.
#'
#' @param X Matrix of predictors.
#' @param Y Matrix of responses.
#' @param Sx Optional X covariance matrix.
#' @param Sy Optional Y covariance matrix.
#' @param lambda Regularization parameter.
#' @param r Rank of the solution.
#' @param highdim Boolean for high-dimensional regime.
#' @param solver Solver type: "rrr", "CVX", or "ADMM".
#' @param LW_Sy Whether to use Ledoit-Wolf shrinkage for Sy.
#' @param standardize Logical; should X and Y be scaled.
#' @param rho ADMM parameter.
#' @param niter Maximum number of iterations for ADMM.
#' @param thresh Convergence threshold.
#' @param verbose Logical for verbose output.
#'
#' @return A list with elements:
#' \describe{
#'   \item{U}{Canonical direction matrix for X (p x r)}
#'   \item{V}{Canonical direction matrix for Y (q x r)}
#'   \item{cor}{Canonical covariances}
#'   \item{loss}{The prediction error 1/n * \| XU - YV\|^2}
#' }
#' @export
cca_rrr <- function(X, Y, Sx=NULL, Sy=NULL,
                    lambda = 0, 
                    r, highdim=TRUE, 
                    solver="ADMM",
                    LW_Sy = FALSE,
                    standardize = TRUE,
                    rho=1,
                    niter=1e4,
                    thresh = 1e-4, verbose=FALSE) {
  
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)

  if (q > n) {
    tmp <- X
    X <- Y
    Y <- tmp
    tmp <- p
    p <- q
    q <- tmp
  }

  X <- if (standardize) scale(X) else scale(X, scale = FALSE)
  Y <- if (standardize) scale(Y) else scale(Y, scale = FALSE)

  if (is.null(Sx)) Sx <- crossprod(X) / n
  if (is.null(Sy)) {
    Sy <- crossprod(Y) / n
    if (LW_Sy) Sy <- as.matrix(corpcor::cov.shrink(Y, verbose=verbose))
  }

  sqrt_inv_Sy <- compute_sqrt_inv(Sy)
  tilde_Y <- Y %*% sqrt_inv_Sy
  Sx_tot <- Sx
  Sxy <- crossprod(X, tilde_Y) / n

  if (!highdim) {
    if(verbose){print("Not using highdim")}
    B_OLS <- solve(Sx_tot, Sxy)
    sqrt_Sx <- compute_sqrt(Sx)
    sqrt_inv_Sx <- compute_sqrt_inv(Sx)
    sol <- svd(sqrt_Sx %*% B_OLS)
    V <- sqrt_inv_Sy %*% sol$v[, 1:r]
    U <- sqrt_inv_Sx %*% sol$u[, 1:r]
  } else {
    if (solver == "CVX") {
      if(verbose){ print("Using CVX solver")}
      B_opt <- solve_rrr_cvxr(X, tilde_Y, lambda) 
    } else if (solver == "ADMM") {
      if(verbose){print("Using ADMM solver")}
      B_opt <- solve_rrr_admm(X, tilde_Y, Sx, lambda=lambda, rho=rho, niter=niter, thresh=thresh, verbose = FALSE)
    } else {
      if(verbose){print("Using gglasso solver")}
      fit <- rrpack::cv.srrr(tilde_Y, X, nrank = r, method = "glasso", nfold = 2,
                     modstr = list("lamA" = rep(lambda, 10), "nlam" = 10))
      B_opt <- fit$coef
    }

    B_opt[abs(B_opt) < 1e-5] <- 0
    active_rows <- which(rowSums(B_opt^2) > 0)
    if (length(active_rows) > r - 1) {
      sqrt_Sx <- compute_sqrt(Sx[active_rows, active_rows])
      sol <- svd(sqrt_Sx %*% B_opt[active_rows, ])
      V <- sqrt_inv_Sy %*% sol$v[, 1:r]
      inv_D <- diag(sapply(sol$d[1:r], function(d) ifelse(d < 1e-4, 0, 1 / d)))
      U <- B_opt %*% sol$v[, 1:r] %*% inv_D
    } else {
      U <- matrix(0, p, r)
      V <- matrix(0, q, r)
    }
  }

  loss <- mean((Y %*% V - X %*% U)^2)
  canon_corr <- sapply(seq_len(r), function(i) cov(X %*% U[, i], Y %*% V[, i]))

  list(U = U, V = V, loss = loss, cor = canon_corr)
}




#' Cross-validated Canonical Correlation Analysis via RRR
#'
#' Performs cross-validation to select optimal lambda, fits CCA_rrr.
#' Canonical Correlation Analysis via Reduced Rank Regression (RRR)
#' @param X Matrix of predictors.
#' @param Y Matrix of responses.
#' @param Sx Optional X covariance matrix.
#' @param Sy Optional Y covariance matrix.
#' @param r Rank of the solution.
#' @param kfolds Number of folds for cross-validation.
#' @param lambdas Sequence of lambda values for cross-validation.
#' @param parallelize Logical; should cross-validation be parallelized?
#' @param highdim Boolean for high-dimensional regime.
#' @param solver Solver type: "rrr", "CVX", or "ADMM".
#' @param LW_Sy Whether to use Ledoit-Wolf shrinkage for Sy.
#' @param standardize Logical; should X and Y be scaled.
#' @param rho ADMM parameter.
#' @param niter Maximum number of iterations for ADMM.
#' @param thresh Convergence threshold.
#' @param verbose Logical for verbose output.
#'
#' @return A list with elements:
#' \describe{
#'   \item{U}{Canonical direction matrix for X (p x r)}
#'   \item{V}{Canonical direction matrix for Y (q x r)}
#'   \item{lambda}{Optimal regularisation parameter lambda chosen by CV}
#'   \item{rmse}{Mean squared error of prediction (as computed in the CV)}
#'   \item{cor}{Canonical correlations}
#' }
#' @importFrom foreach foreach %dopar%
#' @export
cca_rrr_cv <- function(X, Y, 
                       r=2, 
                       lambdas=10^seq(-3, 1.5, length.out = 100),
                       kfolds=14,
                       solver="ADMM",
                       parallelize = FALSE,
                       LW_Sy = FALSE,
                       standardize=TRUE,
                       rho=1,
                       niter=1e4,
                       thresh = 1e-4, verbose=FALSE) {

  X <- if (standardize) scale(X) else scale(X, scale = FALSE)
  Y <- if (standardize) scale(Y) else scale(Y, scale = FALSE)
  n <- nrow(X)
  Sx <- crossprod(X) / n
  Sy <- if (LW_Sy) as.matrix(corpcor::cov.shrink(Y, verbose = FALSE)) else crossprod(Y) / n

  cv_function <- function(lambda) {
    cca_rrr_cv_folds(X, Y, Sx=NULL, Sy=NULL, kfolds=kfolds, 
                  lambda=lambda, r=r, solver=solver, 
                  standardize=standardize, rho=rho, niter=niter, thresh=thresh)
  }

  if (parallelize && solver %in% c("CVX", "CVXR", "ADMM")) {
    no_cores <- parallel::detectCores() - 2
    doParallel::registerDoParallel(cores=no_cores)
    resultsx <- foreach(lambda=lambdas, .combine=rbind, .packages=c('CVXR','Matrix')) %dopar% {
      data.frame(lambda=lambda, rmse=cv_function(lambda))
    }
  } else {
    resultsx <- tibble(lambda = lambdas) %>% 
      mutate(rmse = purrr::map_dbl(lambda, cv_function))
  }

  resultsx <- resultsx %>% 
    mutate(rmse = ifelse(is.na(rmse) | rmse == 0, 1e8, rmse)) %>%
    filter(rmse > 1e-5)

  opt_lambda <- resultsx$lambda[which.min(resultsx$rmse)]
  opt_lambda <- ifelse(is.na(opt_lambda), 0.1, opt_lambda)

  final <- cca_rrr(X, Y, Sx=NULL, Sy=NULL, lambda=opt_lambda, r=r,
                   highdim=TRUE, solver=solver,
                   standardize=standardize, LW_Sy=LW_Sy, rho=rho, niter=niter, 
                   thresh=thresh, verbose=verbose)

  sqrt_inv_Sy <- compute_sqrt_inv(Sy)
  list(U = final$U, 
       V = sqrt_inv_Sy %*% final$V,
       lambda = opt_lambda,
       #resultsx = resultsx,
       rmse = resultsx$rmse,
       cor = sapply(1:r, function(i) cov(X %*% final$U[,i], Y %*% final$V[,i])))
}

#' Cross-validation Fold Evaluation for CCA_rrr
#'
#' Evaluates a single value of lambda using k-fold CV.
#'
#' @return Average RMSE across folds.
cca_rrr_cv_folds <- function(X, Y, Sx, Sy, kfolds=5,
                          lambda=0.01,
                          r=2,
                          standardize=TRUE,
                          solver = "ADMM",
                          rho=1,
                          LW_Sy = TRUE,
                          niter=1e4,
                          thresh = 1e-4) {
  folds <- caret::createFolds(1:nrow(Y), k = kfolds, list = TRUE)

  no_cores <- parallel::detectCores() - 2
  doParallel::registerDoParallel(cores=no_cores) 

  rmse <- foreach(i = seq_along(folds), .combine = c, .packages = c('CVXR', 'Matrix')) %dopar% {
    X_train <- X[-folds[[i]], ]
    Y_train <- Y[-folds[[i]], ]
    X_val <- X[folds[[i]], ]
    Y_val <- Y[folds[[i]], ]

    tryCatch({
      final <- cca_rrr(X_train, Y_train, Sx=NULL, Sy=NULL, highdim=TRUE,
                       lambda=lambda, r=r, solver=solver,
                       LW_Sy=LW_Sy, standardize=standardize, rho=rho, niter=niter, thresh=thresh,
                       verbose=FALSE)
      mean((X_val %*% final$U - Y_val %*% final$V)^2)
    }, error = function(e) {
      message("Error in fold ", i, ": ", conditionMessage(e))
      NA
    })
  }

  if (mean(is.na(rmse)) == 1) return(1e8)
  mean(rmse, na.rm = TRUE)
}

