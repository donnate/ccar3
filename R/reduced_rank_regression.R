# Required libraries
library(dplyr)
library(magrittr)
library(foreach)

.rrr_chol_jitter <- function(G, jitter = 1e-8, max_tries = 6) {
  p <- nrow(G)
  for (k in 0:max_tries) {
    out <- tryCatch(
      chol(G + (10^k) * jitter * diag(p)),
      error = function(e) NULL
    )
    if (!is.null(out)) return(out)
  }
  stop("chol() failed even with jitter.")
}

.rrr_invsqrt_psd <- function(G, eps = 1e-8) {
  G <- (G + t(G)) / 2
  ed <- eigen(G, symmetric = TRUE)
  vals <- pmax(ed$values, eps)
  ed$vectors %*% diag(1 / sqrt(vals), nrow = length(vals)) %*% t(ed$vectors)
}

.rrr_whiten_factor <- function(G, ridge = 1e-8, verbose = FALSE) {
  if (any(!is.finite(G))) {
    if (verbose) cat("\nNon-finite entries in whitening Gram matrix; sanitizing.\n")
    G[!is.finite(G)] <- 0
    G <- (G + t(G)) / 2 + ridge * diag(nrow(G))
  }

  tryCatch({
    R <- .rrr_chol_jitter(G, jitter = ridge)
    backsolve(R, diag(nrow(G)))
  }, error = function(e) {
    if (verbose) cat("\nchol_jitter failed; using eigen whitening fallback.\n")
    .rrr_invsqrt_psd(G, eps = ridge)
  })
}

.rrr_safe_rank_svd <- function(M, k, prefer_sparse = FALSE) {
  if (k < 1) return(NULL)
  nr <- nrow(M)
  nc <- ncol(M)
  if (nr < 1 || nc < 1) return(NULL)
  k <- min(k, nr, nc)

  out <- NULL

  if (requireNamespace("RSpectra", quietly = TRUE) && k < min(nr, nc)) {
    out <- tryCatch({
      if (prefer_sparse && requireNamespace("Matrix", quietly = TRUE)) {
        RSpectra::svds(Matrix::Matrix(M, sparse = TRUE), k)
      } else {
        RSpectra::svds(M, k)
      }
    }, error = function(e) NULL)

    if (!is.null(out)) {
      ok <- all(is.finite(out$u)) && all(is.finite(out$v)) && all(is.finite(out$d))
      if (!ok) out <- NULL
    }
  }

  if (is.null(out)) {
    out <- tryCatch(svd(as.matrix(M), nu = k, nv = k), error = function(e) NULL)
    if (is.null(out)) return(NULL)
    out$u <- out$u[, seq_len(k), drop = FALSE]
    out$v <- out$v[, seq_len(k), drop = FALSE]
    out$d <- out$d[seq_len(k)]
  }

  out
}

.postprocess_rrr_fit <- function(Bhat, X, Y, sqrt_inv_Sy, r,
                                 ridge_whiten = 1e-8, thresh_0 = 1e-6,
                                 verbose = FALSE) {
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)

  if (all(abs(Bhat) < thresh_0)) {
    XU <- matrix(0, n, r)
    YV <- matrix(0, n, r)
    return(list(
      U = matrix(0, p, r),
      V = matrix(0, q, r),
      Lambda = rep(0, r),
      cor = rep(0, r),
      loss = mean((YV - XU)^2)
    ))
  }

  r_eff <- min(r, nrow(Bhat), ncol(Bhat))
  SVDs <- .rrr_safe_rank_svd(Bhat, r_eff, prefer_sparse = TRUE)
  if (is.null(SVDs)) {
    stop("SVD failed during reduced-rank postprocessing.", call. = FALSE)
  }

  U0 <- SVDs$u
  V0 <- sqrt_inv_Sy %*% SVDs$v

  XU0 <- X %*% U0
  YV0 <- Y %*% V0

  GX <- crossprod(XU0) / n + ridge_whiten * diag(r_eff)
  GY <- crossprod(YV0) / n + ridge_whiten * diag(r_eff)

  U <- U0 %*% .rrr_whiten_factor(GX, ridge = ridge_whiten, verbose = verbose)
  V <- V0 %*% .rrr_whiten_factor(GY, ridge = ridge_whiten, verbose = verbose)

  XU <- X %*% U
  YV <- Y %*% V
  cor <- diag(crossprod(XU, YV) / n)

  neg <- cor < 0
  if (any(neg)) {
    V[, neg] <- -V[, neg, drop = FALSE]
    YV[, neg] <- -YV[, neg, drop = FALSE]
    cor[neg] <- -cor[neg]
  }

  ord <- order(cor, decreasing = TRUE)
  cor <- cor[ord]
  U <- U[, ord, drop = FALSE]
  V <- V[, ord, drop = FALSE]
  XU <- XU[, ord, drop = FALSE]
  YV <- YV[, ord, drop = FALSE]

  if (r_eff < r) {
    U <- cbind(U, matrix(0, p, r - r_eff))
    V <- cbind(V, matrix(0, q, r - r_eff))
    cor <- c(cor, rep(0, r - r_eff))
  }

  list(
    U = U,
    V = V,
    Lambda = cor,
    cor = cor,
    loss = mean((YV - XU)^2)
  )
}


# Helper: efficient (Sx + rho I)^{-1} * W without forming large Sx when p >> n
make_invSx_apply <- function(X, Sx, rho, ridge = 1e-8, verbose = FALSE) {
  p <- ncol(X)
  if (!is.null(Sx)) {
    invSx <- solve(Sx + rho * diag(p))
    return(list(
      apply = function(W) invSx %*% W,
      uses_woodbury = FALSE
    ))
  }

  n <- nrow(X)
  if (n < 1 || p < 1) stop("X must have positive dimensions.")

  # Woodbury: (rho I + X'X/n)^{-1} = (1/rho)I - (1/rho^2) X' (I + XX'/(n rho))^{-1} X / n
  M <- diag(n) + tcrossprod(X) / (n * rho)
  R <- tryCatch(chol(M + ridge * diag(n)), error = function(e) NULL)

  solve_M <- if (!is.null(R)) {
    function(B) backsolve(R, forwardsolve(t(R), B))
  } else {
    if (verbose) message("Cholesky failed in Woodbury solve; using solve().")
    function(B) solve(M, B)
  }

  list(
    apply = function(W) {
      XW <- (X %*% W) / n
      tmp <- solve_M(XW)
      (1 / rho) * W - (1 / rho^2) * crossprod(X, tmp)
    },
    uses_woodbury = TRUE
  )
}



# Helper: ADMM-based group sparse solver
solve_rrr_admm <- function(X, tilde_Y, Sx, lambda, rho=1, niter=10, thresh, verbose = FALSE, thresh_0 = 1e-6,
                           invSx_apply = NULL) {
  p <- ncol(X); q <- ncol(tilde_Y)
  n <- nrow(X)

  if (is.null(invSx_apply)) {
    invSx_info <- make_invSx_apply(X, Sx, rho, verbose = verbose)
    invSx_apply <- invSx_info$apply
  }

  U <- Z <- matrix(0, p, q)
  
  prod_xy <- crossprod(X, tilde_Y) / n
  for (i in seq_len(niter)) {
    U_old <- U; Z_old <- Z
    B <- invSx_apply(prod_xy + rho * (Z - U))
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
  
  B_opt[abs(B_opt) < thresh_0] <- 0
  return(B_opt)
}


# Helper: CVXR-based group sparse solver
solve_rrr_cvxr <- function(X, tilde_Y, lambda, thresh_0=1e-6) {
  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop("Package 'CVXR' must be installed to use the CVX solver.",
         call. = FALSE)
  }

  p <- ncol(X); q <- ncol(tilde_Y)
  n <- nrow(X)
  B <- CVXR::Variable(p, q)
  
  objective <- CVXR::Minimize(1 / n * CVXR::sum_squares(tilde_Y - X %*% B) + lambda * sum(CVXR::norm2(B, axis = 1)))
  result <- CVXR::solve(CVXR::Problem(objective))
  B_opt <- result$getValue(B)
  
  B_opt[abs(B_opt) < thresh_0] <- 0
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
#' @param thresh_0 For the ADMM solver: Set entries whose absolute value is below this to 0 (default 1e-6).
#' @param verbose Logical for verbose output.
#'
#' @return A list with elements:
#' \itemize{
#'   \item U: Canonical direction matrix for X (p x r)
#'   \item V: Canonical direction matrix for Y (q x r)
#'   \item cor: Canonical covariances
#'   \item loss: The prediction error 1/n * || XU - YV ||^2
#' }

#' @export
cca_rrr <- function(X, Y, Sx=NULL, Sy=NULL,
                    lambda = 0, 
                    r, highdim=TRUE, 
                    solver="ADMM",
                    LW_Sy = TRUE,
                    standardize = TRUE,
                    rho=1,
                    niter=1e4,
                    thresh = 1e-4, thresh_0 = 1e-6,
                    verbose=FALSE) {
  
  n <- nrow(X)
  p <- ncol(X)
  
  X <- if (standardize) scale(X) else X #scale(X, scale = FALSE)
  Y <- if (standardize) scale(Y) else Y # scale(Y, scale = FALSE)
  
  use_highdim <- isTRUE(highdim) || p > n
  skip_Sx <- use_highdim && is.null(Sx) && p > n
  if (is.null(Sx) && !skip_Sx) Sx <- crossprod(X) / n
  if (is.null(Sy)) {
    Sy <- crossprod(Y) / n
    if (LW_Sy) Sy <- as.matrix(corpcor::cov.shrink(Y, verbose=verbose))
  }
  
  sqrt_inv_Sy <- compute_sqrt_inv(Sy)
  tilde_Y <- Y %*% sqrt_inv_Sy
  
  if (!use_highdim) {
    if(verbose){print("Not using highdim")}
    B_fit <- solve(Sx, crossprod(X, tilde_Y) / n)
  } else {
    if (solver == "CVX") {
      if(verbose){ print("Using CVX solver")}
      B_fit <- solve_rrr_cvxr(X, tilde_Y, lambda, thresh_0=thresh_0) 
    } else if (solver == "ADMM") {
      if(verbose){print("Using ADMM solver")}
      B_fit <- solve_rrr_admm(X, tilde_Y, Sx, lambda=lambda, rho=rho, niter=niter, thresh=thresh, thresh_0=thresh_0,
                              verbose = FALSE)
    } else {
      if(verbose){print("Using gglasso solver")}
        if (!requireNamespace("rrpack", quietly = TRUE)) {
          stop("Package 'rrpack' must be installed to use the rrpack solver.",
              call. = FALSE)
        }
      fit <- rrpack::cv.srrr(tilde_Y, X, nrank = r, method = "glasso", nfold = 2,
                             modstr = list("lamA" = rep(lambda, 10), "nlam" = 10))
      B_fit <- fit$coef
    }
  }
  B_fit[abs(B_fit) < thresh_0] <- 0
  post <- .postprocess_rrr_fit(
    B_fit, X, Y, sqrt_inv_Sy, r,
    ridge_whiten = max(thresh_0, 1e-8),
    thresh_0 = thresh_0,
    verbose = verbose
  )
  
  list(U = post$U, V = post$V, Lambda = post$Lambda, loss = post$loss, cor = post$cor, B_opt = B_fit)
}




#' Cross-validated Canonical Correlation Analysis via RRR
#'
#' Performs cross-validation to select optimal lambda, fits CCA_rrr.
#' Canonical Correlation Analysis via Reduced Rank Regression (RRR)
#' @param X Matrix of predictors.
#' @param Y Matrix of responses.
#' @param r Rank of the solution.
#' @param kfolds Number of folds for cross-validation.
#' @param lambdas Sequence of lambda values for cross-validation.
#' @param parallelize Logical; should cross-validation be parallelized?
#' @param solver Solver type: "rrr", "CVX", or "ADMM".
#' @param LW_Sy Whether to use Ledoit-Wolf shrinkage for Sy.
#' @param standardize Logical; should X and Y be scaled.
#' @param rho ADMM parameter.
#' @param niter Maximum number of iterations for ADMM.
#' @param thresh Convergence threshold.
#' @param verbose Logical for verbose output.
#' @param thresh_0 tolerance for declaring entries non-zero
#' @param nb_cores Number of cores to use for parallelization. Defaults to min(kfolds, available cores minus 1).
#'
#' @return A list with elements:
#' \itemize{
#'   \item U: Canonical direction matrix for X (p x r)
#'   \item V: Canonical direction matrix for Y (q x r)
#'   \item lambda: Optimal regularisation parameter lambda chosen by CV
#'   \item rmse: Mean squared error of prediction (as computed in the CV)
#'   \item cor: Canonical correlations
#' }
#' @importFrom foreach foreach %dopar%
#' @importFrom foreach foreach %do%
#' @export
cca_rrr_cv <- function(X, Y, 
                       r=2, 
                       lambdas=10^seq(-3, 1.5, length.out = 100),
                       kfolds=14,
                       solver="ADMM",
                       parallelize = FALSE,
                       LW_Sy = TRUE,
                       standardize=TRUE,
                       rho=1,
                       thresh_0=1e-6,
                       niter=1e4,
                       thresh = 1e-4, verbose=FALSE,
                       nb_cores = NULL) {
  
  X <- if (standardize) scale(X) else X #scale(X, scale = FALSE)
  Y <- if (standardize) scale(Y) else Y #scale(Y, scale = FALSE)
  n <- nrow(X)
  p <- ncol(X)
  Sx <- if (p > n) NULL else matmul(t(X), X) / n
  Sy <- if (LW_Sy) NULL else crossprod(Y) / n
  folds <- caret::createFolds(1:nrow(Y), k = kfolds, list = TRUE)
  cv_function <- function(lambda) {
    #print(lambda)
    cca_rrr_cv_folds(X, Y, Sx=Sx, Sy=Sy, kfolds=kfolds, 
                     LW_Sy = LW_Sy,
                     lambda=lambda, r=r, solver=solver, 
                     standardize=FALSE, rho=rho, niter=niter, thresh=thresh,
                     thresh_0=thresh_0, folds = folds)
  }
  
  if (parallelize && solver %in% c("CVX", "CVXR", "ADMM")) {

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
        if (verbose){
          cat(crayon::green("Parallel backend successfully registered.\n"))
        }
        doParallel::registerDoParallel(cl)
        on.exit(parallel::stopCluster(cl), add = TRUE)
      } else {
        # If setup_parallel_backend returned NULL, print a warning and proceed serially
        warning("All parallel setup attempts failed. Proceeding in serial mode.", immediate. = TRUE)
        parallelize <- FALSE # Ensure %dopar% runs serially
    }
   }

  if (parallelize && solver %in% c("CVX", "CVXR", "ADMM")) {
    if (solver  %in% c("CVX", "CVXR" )){
        if (!requireNamespace("CVXR", quietly = TRUE)) {
         stop("Package 'CVXR' must be installed to use the CVXR/CVX solver.",
         call. = FALSE)
          }

      resultsx <- foreach(lambda=lambdas, .combine=rbind, .packages=c('CVXR','Matrix')) %dopar% {
      data.frame(lambda=lambda, rmse=cv_function(lambda))
    }

    }else{
      resultsx <- foreach(lambda=lambdas, .combine=rbind) %dopar% {
      data.frame(lambda=lambda, rmse=cv_function(lambda))
    }
    }
    
  } else {
    resultsx <- data.frame(lambda = lambdas)
    resultsx$rmse <- sapply(lambdas, cv_function)
    
  }
  
  resultsx <- resultsx %>% 
    dplyr::mutate(rmse = ifelse(is.na(rmse) | rmse == 0, 1e8, rmse)) %>%
    dplyr::filter(rmse > 1e-5)
  
  opt_lambda <- resultsx$lambda[which.min(resultsx$rmse)]
  opt_lambda <- ifelse(is.na(opt_lambda), 0.1, opt_lambda)

  final <- cca_rrr(X, Y, Sx=NULL, Sy=NULL, lambda=opt_lambda, r=r,
                   highdim=TRUE, solver=solver,
                   standardize=FALSE, LW_Sy=LW_Sy, rho=rho, niter=niter, 
                   thresh=thresh, thresh_0=thresh_0,verbose=verbose)


  return(list(U = final$U, 
       V = final$V,
       lambda = opt_lambda,
       resultsx = resultsx,
       rmse = resultsx$rmse,
       cor = final$cor,
       Lambda = final$Lambda,
       B = final$B_opt
       ))
}




cca_rrr_cv_folds <- function(X, Y, Sx, Sy, kfolds=5,
                             lambda=0.01,
                             r=2,
                             standardize=FALSE,
                             solver = "ADMM",
                             rho=1,
                             LW_Sy = TRUE,
                             niter=1e4,
                             thresh_0=1e-6,
                             thresh = 1e-4,
                             folds = NULL) {
  if (is.null(folds)) {
    folds <- caret::createFolds(1:nrow(Y), k = kfolds, list = TRUE)
  }
  
  rmse <- foreach(i = seq_along(folds), .combine = c) %do% { 
    
    n <- nrow(X)
    X_train <- X[-folds[[i]], ]; Y_train <- Y[-folds[[i]], ]
    X_val <- X[folds[[i]], ]; Y_val <- Y[folds[[i]], ]
    n_train <- n - nrow(X_val)
    
    if (is.null(Sx) == FALSE) {
      Sx_train <- (n * Sx - crossprod(X_val)) / n_train
    } else {
      Sx_train <- NULL
    }

    if (is.null(Sy) == FALSE) {
      Sy_train <- (n * Sy - crossprod(Y_val)) / n_train
    } else {
      Sy_train <- NULL
    }
    
    tryCatch({
      final <- cca_rrr(X_train, Y_train, Sx=Sx_train, Sy=Sy_train, highdim=TRUE,
                       lambda=lambda, r=r, solver=solver,
                       LW_Sy=LW_Sy, standardize=FALSE, rho=rho, niter=niter, 
                       thresh=thresh, thresh_0=thresh_0,
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
