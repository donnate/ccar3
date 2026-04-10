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

.resolve_preprocess_mode <- function(standardize = NULL, preprocess = NULL) {
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

.rrr_match_mode <- function(mode) {
  if (missing(mode) || is.null(mode) || length(mode) == 0L || all(is.na(mode))) {
    return("sqrtm_norm")
  }

  mode <- tolower(as.character(mode[[1]]))
  if (identical(mode, "new")) {
    mode <- "sqrtm_norm"
  } else if (identical(mode, "old")) {
    mode <- "product_norm"
  }

  match.arg(mode, c("sqrtm_norm", "product_norm"))
}

.rrr_match_cv_metric <- function(cv_metric) {
  if (missing(cv_metric) || is.null(cv_metric) || length(cv_metric) == 0L || all(is.na(cv_metric))) {
    metric <- "mse"
  } else {
    metric <- tolower(as.character(cv_metric[[1]]))
  }

  if (metric %in% c("mse", "rmse", "prediction")) {
    return("mse")
  }

  if (metric %in% c("correlation", "cor", "corr", "association")) {
    return("correlation")
  }

  stop(
    "Unsupported `cv_metric`. Use one of: 'mse' or 'correlation'.",
    call. = FALSE
  )
}

.rrr_correlation_score <- function(X_scores, Y_scores) {
  X_scores <- as.matrix(X_scores)
  Y_scores <- as.matrix(Y_scores)

  if (nrow(X_scores) == 0L || nrow(Y_scores) == 0L ||
      ncol(X_scores) == 0L || ncol(Y_scores) == 0L) {
    return(0)
  }

  x_sd <- apply(X_scores, 2, stats::sd)
  y_sd <- apply(Y_scores, 2, stats::sd)

  keep_x <- which(is.finite(x_sd) & x_sd > 0)
  keep_y <- which(is.finite(y_sd) & y_sd > 0)

  if (length(keep_x) == 0L || length(keep_y) == 0L) {
    return(0)
  }

  C <- stats::cor(
    X_scores[, keep_x, drop = FALSE],
    Y_scores[, keep_y, drop = FALSE]
  )
  C[!is.finite(C)] <- 0

  sum(svd(C, nu = 0, nv = 0)$d)
}


.postprocess_rrr_fit <- function(Bhat, X, Y, sqrt_inv_Sy, r,
                                 ridge_whiten = 0, thresh_0 = 1e-6,
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


.postprocess_rrr_fit_old <- function(Bhat, X, Y, sqrt_inv_Sy, r,
                                 ridge_whiten = 0, thresh_0 = 1e-6,
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
  active_rows <- which(rowSums(Bhat^2) > 0)
  Sx_active = crossprod(X[, active_rows, drop = FALSE]) / n
  SVDs <- .rrr_safe_rank_svd(Bhat, r_eff, prefer_sparse = TRUE)
  if (is.null(SVDs)) {
    stop("SVD failed during reduced-rank postprocessing.", call. = FALSE)
  }
  if (length(active_rows) > r - 1) {
      sqrt_Sx <- compute_sqrt(Sx_active)
      sol <- svd(sqrt_Sx %*% Bhat[active_rows, ])
      V <- sqrt_inv_Sy %*% sol$v[, 1:r]
      inv_D <- diag(sapply(sol$d[1:r], function(d) ifelse(d < 1e-4, 0, 1 / d)))
      U <- Bhat %*% sol$v[, 1:r] %*% inv_D
    } else {
      U <- matrix(0, p, r)
      V <- matrix(0, q, r)
    }

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


.rrr_cg_solve <- function(A_mul, b, M_inv = NULL, tol = 1e-6,
                          maxiter = NULL, x0 = NULL) {
  p <- length(b)
  if (is.null(maxiter)) {
    maxiter <- min(p, 1000L)
  }

  x <- if (is.null(x0)) numeric(p) else as.numeric(x0)
  r <- as.numeric(b) - as.numeric(A_mul(x))
  bnorm <- sqrt(sum(b^2))
  if (!is.finite(bnorm) || bnorm < .Machine$double.eps) {
    return(x)
  }

  z <- if (is.null(M_inv)) r else M_inv(r)
  rz_old <- drop(crossprod(r, z))
  if (!is.finite(rz_old) || rz_old <= 0) {
    return(x)
  }

  p_dir <- z
  tol_abs <- tol * max(1, bnorm)

  for (iter in seq_len(maxiter)) {
    Ap <- as.numeric(A_mul(p_dir))
    denom <- drop(crossprod(p_dir, Ap))
    if (!is.finite(denom) || abs(denom) < .Machine$double.eps) {
      break
    }

    alpha <- rz_old / denom
    x <- x + alpha * p_dir
    r <- r - alpha * Ap
    if (sqrt(sum(r^2)) <= tol_abs) {
      break
    }

    z <- if (is.null(M_inv)) r else M_inv(r)
    rz_new <- drop(crossprod(r, z))
    if (!is.finite(rz_new) || rz_new <= 0) {
      break
    }

    beta <- rz_new / rz_old
    p_dir <- z + beta * p_dir
    rz_old <- rz_new
  }

  x
}

# Helper: efficient (Sx + rho I)^{-1} * W without forcing an explicit Sx build.
# For very large n and p, use a matrix-free CG solve that only applies X and t(X).
make_invSx_apply <- function(X, Sx = NULL, rho, ridge = 1e-8, verbose = FALSE,
                             matrix_free_threshold = 4000L,
                             cg_tol = 1e-6,
                             cg_maxiter = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  if (n < 1 || p < 1) stop("X must have positive dimensions.")

  solve_from_factor <- function(G) {
    R <- .rrr_chol_jitter(G, jitter = ridge)
    function(B) backsolve(R, forwardsolve(t(R), as.matrix(B)))
  }

  if (!is.null(Sx)) {
    if (verbose) message("Using direct Cholesky solve with supplied Sx.")
    solve_gram <- solve_from_factor(Sx + rho * diag(p))
    return(list(
      apply = solve_gram,
      strategy = "direct",
      uses_woodbury = FALSE
    ))
  }

  if (min(n, p) >= matrix_free_threshold) {
    if (verbose) message("Using matrix-free CG solve for ADMM update.")

    diag_precond <- rho + colSums(X^2) / n
    diag_precond[!is.finite(diag_precond) | diag_precond < ridge] <- ridge

    A_mul <- function(v) {
      as.numeric(rho * v + crossprod(X, X %*% v) / n)
    }
    M_inv <- function(v) {
      as.numeric(v / diag_precond)
    }

    last_sol <- NULL
    return(list(
      apply = function(W) {
        W <- as.matrix(W)
        out <- matrix(0, p, ncol(W))
        if (!is.null(last_sol) && identical(dim(last_sol), dim(out))) {
          out[,] <- last_sol
        }

        for (j in seq_len(ncol(W))) {
          out[, j] <- .rrr_cg_solve(
            A_mul, W[, j],
            M_inv = M_inv,
            tol = cg_tol,
            maxiter = cg_maxiter,
            x0 = out[, j]
          )
        }

        last_sol <<- out
        out
      },
      strategy = "matrix_free_cg",
      uses_woodbury = FALSE
    ))
  }

  if (p <= n) {
    if (verbose) message("Using direct Cholesky solve from X without precomputing Sx upstream.")
    solve_gram <- solve_from_factor(crossprod(X) / n + rho * diag(p))
    return(list(
      apply = solve_gram,
      strategy = "direct",
      uses_woodbury = FALSE
    ))
  }

  if (verbose) message("Using Woodbury solve for ADMM update.")
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
    strategy = "woodbury",
    uses_woodbury = TRUE
  )
}



# Helper: ADMM-based group sparse solver
solve_rrr_admm <- function(X, tilde_Y, Sx = NULL, lambda, rho=1, niter=10, thresh, verbose = FALSE, thresh_0 = 0,
                           invSx_apply = NULL,
                           matrix_free_threshold = 4000L,
                           cg_tol = 1e-6,
                           cg_maxiter = NULL) {
  p <- ncol(X); 
  q <- ncol(tilde_Y)
  n <- nrow(X)

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

  U <- Z <- matrix(0, p, q)
  
  prod_xy <- crossprod(X, tilde_Y) / n ### size is p x q, so okay since q <<< n
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
solve_rrr_cvxr <- function(X, tilde_Y, lambda, thresh_0=1e-6,
                           rho = 1, niter = 1e4, thresh = 1e-4,
                           verbose = FALSE,
                           matrix_free_threshold = 4000L,
                           cg_tol = 1e-6,
                           cg_maxiter = NULL) {
  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop("Package 'CVXR' must be installed to use the CVX solver.",
         call. = FALSE)
  }

  p <- ncol(X); q <- ncol(tilde_Y)
  n <- nrow(X)

  fallback_admm <- function(reason = NULL) {
    if (verbose && !is.null(reason)) {
      message("Falling back to ADMM from CVXR: ", reason)
    }
    solve_rrr_admm(
      X, tilde_Y, Sx = NULL, lambda = lambda, rho = rho,
      niter = niter, thresh = thresh, thresh_0 = thresh_0,
      verbose = verbose,
      matrix_free_threshold = matrix_free_threshold,
      cg_tol = cg_tol,
      cg_maxiter = cg_maxiter
    )
  }

  # CVXR becomes impractical quickly because canonicalization densifies the
  # quadratic objective. For larger problems, the ADMM path is both faster and
  # dramatically more memory efficient.
  if (n * p * q > 5e6) {
    return(fallback_admm("problem size is too large for CVXR canonicalization"))
  }

  B <- CVXR::Variable(c(p, q))
  
  row_penalty <- Reduce(
    `+`,
    lapply(seq_len(p), function(j) CVXR::norm(B[j, , drop = FALSE], 2))
  )
  objective <- CVXR::Minimize(
    1 / n * CVXR::sum_squares(tilde_Y - X %*% B) + lambda * row_penalty
  )
  B_opt <- tryCatch({
    problem <- CVXR::Problem(objective)
    result <-  CVXR::psolve(problem)
    result$getValue(B)
  }, error = function(e) {
    msg <- conditionMessage(e)
    if (grepl("memory|vector memory|sparse->dense|cannot allocate", msg, ignore.case = TRUE)) {
      return(fallback_admm(msg))
    }
    stop(e)
  })
  
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
#' @param mode Mode for postprocessing the RRR solution. One of `"sqrtm_norm"`
#'   (default) or `"product_norm"`. Legacy aliases `"new"` and `"old"` are
#'   also accepted. The former whitens the canonical variates to have identity
#'   covariance, while the latter does not whiten and instead returns the raw
#'   SVD factors of the RRR solution. The `"product_norm"` mode may be more
#'   interpretable in some cases but can yield canonical variates with very
#'   different scales and is not guaranteed to be numerically stable when the
#'   RRR solution is very low-rank or nearly low-rank.
#' @param highdim Boolean for high-dimensional regime.
#' @param solver Solver type: "rrr", "CVX", or "ADMM".
#' @param LW_Sy Whether to use Ledoit-Wolf shrinkage for Sy.
#' @param standardize Backward-compatible preprocessing flag: TRUE = `"scale"`, FALSE = `"center"`.
#' @param preprocess Preprocessing mode. One of `"scale"` (center + scale), `"center"` (center only), or `"none"`.
#'   Logical values are still accepted for backward compatibility: `TRUE` uses `standardize`, `FALSE` skips preprocessing.
#' @param rho ADMM parameter.
#' @param niter Maximum number of iterations for ADMM.
#' @param thresh Convergence threshold.
#' @param thresh_0 For the ADMM solver: Set entries whose absolute value is below this to 0 (default 1e-6).
#' @param matrix_free_threshold For ADMM: when both `n` and `p` are at least this value, use a matrix-free conjugate-gradient solve instead of forming a dense linear system.
#' @param cg_tol Relative tolerance for the matrix-free conjugate-gradient solve used in ADMM.
#' @param cg_maxiter Maximum iterations for the matrix-free conjugate-gradient solve. Defaults to `min(p, 1000)`.
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
                    mode = "sqrtm_norm",
                    standardize = FALSE,
                    preprocess = NULL,
                    rho=1,
                    niter=1e4,
                    thresh = 1e-4, thresh_0 = 0,
                    matrix_free_threshold = 4000L,
                    cg_tol = 1e-6,
                    cg_maxiter = NULL,
                    verbose=FALSE) {
  
  preprocess_mode <- .resolve_preprocess_mode(standardize = standardize, preprocess = preprocess)
  mode <- .rrr_match_mode(mode)
  X <- .preprocess_matrix(X, preprocess_mode)
  Y <- .preprocess_matrix(Y, preprocess_mode)

  n <- nrow(X)
  p <- ncol(X)
  
  if (is.null(Sx) && !highdim) Sx <- crossprod(X) / n
  if (is.null(Sy)) {
    Sy <- crossprod(Y) / n
    if (LW_Sy) Sy <- as.matrix(corpcor::cov.shrink(Y, verbose=verbose))
  }
  
  sqrt_inv_Sy <- compute_sqrt_inv(Sy)
  tilde_Y <- Y %*% sqrt_inv_Sy
  
  if (!highdim) {
    if (verbose) print("Not using highdim")
    B_fit <- solve(Sx, crossprod(X, tilde_Y) / n)
  } else {
    if (solver == "CVX") {
      if(verbose){ print("Using CVX solver")}
      B_fit <- solve_rrr_cvxr(
        X, tilde_Y, lambda, thresh_0 = thresh_0,
        rho = rho, niter = niter, thresh = thresh, verbose = verbose,
        matrix_free_threshold = matrix_free_threshold,
        cg_tol = cg_tol,
        cg_maxiter = cg_maxiter
      )
    } else if (solver == "ADMM") {
      if(verbose){print("Using ADMM solver")}
      B_fit <- solve_rrr_admm(X, tilde_Y, Sx = Sx,
                              lambda=lambda, rho=rho, niter=niter, 
                              thresh=thresh, thresh_0=thresh_0,
                              verbose = verbose,
                              matrix_free_threshold = matrix_free_threshold,
                              cg_tol = cg_tol,
                              cg_maxiter = cg_maxiter)
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
  if(mode =="product_norm"){
    post <- .postprocess_rrr_fit_old(
    B_fit, X, Y, sqrt_inv_Sy, r,
    ridge_whiten = max(thresh_0, 1e-18),
    thresh_0 = thresh_0,
    verbose = verbose
  )

  }else{
    post <- .postprocess_rrr_fit(
    B_fit, X, Y, sqrt_inv_Sy, r,
    ridge_whiten = max(thresh_0, 1e-18),
    thresh_0 = thresh_0,
    verbose = verbose
  )

  }


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
#' @param mode Mode for postprocessing the RRR solution. One of `"sqrtm_norm"`
#'   (default) or `"product_norm"`. Legacy aliases `"new"` and `"old"` are
#'   also accepted. The former whitens the canonical variates to have identity
#'   covariance, while the latter does not whiten and instead returns the raw
#'   SVD factors of the RRR solution. The `"product_norm"` mode may be more
#'   interpretable in some cases but can yield canonical variates with very
#'   different scales and is not guaranteed to be numerically stable when the
#'   RRR solution is very low-rank or nearly low-rank.
#' @param solver Solver type: "rrr", "CVX", or "ADMM".
#' @param LW_Sy Whether to use Ledoit-Wolf shrinkage for Sy.
#' @param standardize Backward-compatible preprocessing flag: TRUE = `"scale"`, FALSE = `"center"`.
#' @param preprocess Preprocessing mode. One of `"scale"` (center + scale), `"center"` (center only), or `"none"`.
#'   Logical values are still accepted for backward compatibility: `TRUE` uses `standardize`, `FALSE` skips preprocessing.
#' Data are preprocessed once up front, and fold fits reuse those transformed matrices without re-centering or re-scaling inside each fold.
#' @param rho ADMM parameter.
#' @param niter Maximum number of iterations for ADMM.
#' @param thresh Convergence threshold.
#' @param cv_metric Cross-validation metric. Use `"mse"` to minimize held-out prediction
#'   error or `"correlation"` to maximize held-out association between `X %*% U`
#'   and `Y %*% V`.
#' @param verbose Logical for verbose output.
#' @param thresh_0 tolerance for declaring entries non-zero
#' @param matrix_free_threshold For ADMM: when both `n` and `p` are at least this value, use a matrix-free conjugate-gradient solve instead of forming a dense linear system.
#' @param cg_tol Relative tolerance for the matrix-free conjugate-gradient solve used in ADMM.
#' @param cg_maxiter Maximum iterations for the matrix-free conjugate-gradient solve. Defaults to `min(p, 1000)`.
#' @param nb_cores Number of cores to use for parallelization. Defaults to min(kfolds, available cores minus 1).
#'
#' @return A list with elements:
#' \describe{
#'   \item{U}{Canonical direction matrix for X (p x r)}
#'   \item{V}{Canonical direction matrix for Y (q x r)}
#'   \item{lambda}{Optimal regularisation parameter lambda chosen by CV}
#'   \item{rmse}{Backward-compatible optimization objective. For `cv_metric = "mse"`
#'   this is the held-out mean squared error; for `cv_metric = "correlation"` it is
#'   `-cv_score` so that smaller still means better.}
#'   \item{cv_score}{Raw held-out cross-validation score averaged across folds.}
#'   \item{cv_metric}{The metric used to score lambdas during cross-validation.}
#'   \item{cor}{Canonical correlations at the selected lambda}
#'   \item{lambda_x}{Alias of the selected `lambda`}
#'   \item{lambda_x_se}{Foldwise standard error at the selected `lambda`}
#'   \item{lambda_y}{Placeholder for symmetry with two-penalty interfaces}
#'   \item{lambda_y_se}{Placeholder for symmetry with two-penalty interfaces}
#'   \item{resultsx}{Backward-compatible alias of `cv_summary`}
#'   \item{cv_summary}{Data frame with one row per lambda containing the mean CV score
#'   and its foldwise standard error.}
#'   \item{cv_folds}{Data frame with fold-level CV scores for each lambda.}
#'   \item{Lambda}{Canonical correlations from the final fit}
#'   \item{B}{Estimated reduced-rank coefficient matrix from the final fit}
#'   \item{fit}{Final fit at the selected lambda}
#' }
#' @importFrom foreach foreach %dopar%
#' @importFrom foreach foreach %do%
#' @export
cca_rrr_cv <- function(X, Y, 
                       r=2, 
                       lambdas=10^seq(-3, 1.5, length.out = 100),
                       kfolds=10,
                       solver="ADMM",
                       mode = "sqrtm_norm",
                       parallelize = FALSE,
                       LW_Sy = TRUE,
                       standardize = FALSE,
                       preprocess = NULL,
                       cv_metric = "mse",
                       rho=1,
                       thresh_0=0,
                       niter=1e4,
                       matrix_free_threshold = 4000L,
                       cg_tol = 1e-6,
                       cg_maxiter = NULL,
                       thresh = 1e-4, verbose=FALSE,
                       nb_cores = NULL) {
  preprocess_mode <- .resolve_preprocess_mode(standardize = standardize, preprocess = preprocess)
  mode <- .rrr_match_mode(mode)
  cv_metric <- .rrr_match_cv_metric(cv_metric)
  X <- .preprocess_matrix(X, preprocess_mode)
  Y <- .preprocess_matrix(Y, preprocess_mode)

  n <- nrow(X)
  Sx <- NULL
  Sy <- if (LW_Sy) NULL else crossprod(Y) / n
  folds <- .create_cv_folds(n, kfolds)
  safe_mean_na <- .safe_mean_na
  safe_sd_na <- .safe_sd_na
  cv_function <- function(lambda) {
    fold_values <- cca_rrr_cv_folds(
      X, Y, Sx = Sx, Sy = Sy, kfolds = kfolds,
      LW_Sy = LW_Sy,
      preprocess = "none",
      cv_metric = cv_metric,
      mode = mode,
      lambda = lambda, r = r, solver = solver,
      rho = rho, niter = niter, thresh = thresh,
      matrix_free_threshold = matrix_free_threshold,
      cg_tol = cg_tol,
      cg_maxiter = cg_maxiter,
      thresh_0 = thresh_0, folds = folds,
      return_fold_values = TRUE
    )

    fold_objective <- if (cv_metric == "mse") fold_values else -fold_values

    list(
      summary = data.frame(
        lambda = lambda,
        rmse = safe_mean_na(fold_objective),
        cv_score = safe_mean_na(fold_values),
        cv_metric = cv_metric,
        se = safe_sd_na(fold_values) / sqrt(sum(!is.na(fold_values))),
        stringsAsFactors = FALSE
      ),
      fold = data.frame(
        lambda = lambda,
        fold = seq_along(fold_values),
        rmse = fold_objective,
        cv_score = fold_values,
        cv_metric = cv_metric,
        stringsAsFactors = FALSE
      )
    )
  }

  foreach <- foreach::foreach
  `%dopar%` <- foreach::`%dopar%`
  
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
        initialize_parallel_workers(cl, verbose = verbose)
        doParallel::registerDoParallel(cl)
        on.exit(cleanup_parallel_backend(cl), add = TRUE)
      } else {
        # If setup_parallel_backend returned NULL, print a warning and proceed serially
        warning("All parallel setup attempts failed. Proceeding in serial mode.", immediate. = TRUE)
        parallelize <- FALSE # Ensure %dopar% runs serially
    }
   }

  if (parallelize && solver %in% c("CVX", "CVXR", "ADMM")) {
    rrr_parallel_exports <- c(
      "cca_rrr_cv_folds", "cca_rrr", "solve_rrr_admm", "solve_rrr_cvxr",
      "make_invSx_apply", ".rrr_cg_solve", ".rrr_chol_jitter", ".resolve_preprocess_mode",
      ".preprocess_matrix", ".create_cv_folds",
      ".rrr_invsqrt_psd", ".rrr_whiten_factor", ".rrr_safe_rank_svd",
      ".postprocess_rrr_fit", "compute_sqrt_inv", ".rrr_match_cv_metric",
      ".rrr_correlation_score",
      ".safe_mean_na", ".safe_sd_na"
    )
    if (solver  %in% c("CVX", "CVXR" )){
        if (!requireNamespace("CVXR", quietly = TRUE)) {
         stop("Package 'CVXR' must be installed to use the CVXR/CVX solver.",
         call. = FALSE)
          }

      results_list <- foreach(
        lambda = lambdas,
        .packages = c("CVXR", "Matrix"),
        .export = rrr_parallel_exports
      ) %dopar% {
        cv_function(lambda)
      }

    }else{
      results_list <- foreach(
        lambda = lambdas,
        .export = rrr_parallel_exports
      ) %dopar% {
        cv_function(lambda)
      }
    }
    
  } else {
    results_list <- lapply(lambdas, cv_function)
  }
  
  cv_summary <- dplyr::bind_rows(lapply(results_list, `[[`, "summary"))
  cv_folds <- dplyr::bind_rows(lapply(results_list, `[[`, "fold"))

  if (cv_metric == "mse") {
    cv_summary <- cv_summary %>% 
      dplyr::mutate(rmse = ifelse(is.na(rmse) | rmse == 0, 1e8, rmse)) %>%
      dplyr::filter(rmse > 1e-5)
  } else {
    cv_summary <- cv_summary[is.finite(cv_summary$cv_score), , drop = FALSE]
  }
  
  resultsx <- cv_summary

  if (nrow(cv_summary) == 0L) {
    stop("No valid lambda values were found during cross-validation.", call. = FALSE)
  }

  opt_lambda <- if (cv_metric == "mse") {
    cv_summary$lambda[which.min(cv_summary$rmse)]
  } else {
    cv_summary$lambda[which.max(cv_summary$cv_score)]
  }
  opt_lambda <- ifelse(is.na(opt_lambda), 0.1, opt_lambda)

  final <- cca_rrr(X, Y, Sx=NULL, Sy=NULL, 
                   lambda=opt_lambda, r=r,
                   highdim=TRUE, solver=solver,
                   mode =mode,
                   standardize=FALSE, preprocess="none",
                   LW_Sy=LW_Sy, rho=rho, niter=niter, 
                   matrix_free_threshold = matrix_free_threshold,
                   cg_tol = cg_tol,
                   cg_maxiter = cg_maxiter,
                   thresh=thresh, thresh_0=thresh_0,
                   verbose=verbose)


  return(list(U = final$U, 
       V = final$V,
       lambda = opt_lambda,
       rmse = cv_summary$rmse,
       cv_score = cv_summary$cv_score,
       cv_metric = cv_metric,
       cor = final$cor,
       lambda_x = opt_lambda,
       lambda_x_se = cv_summary$se[match(opt_lambda, cv_summary$lambda)],
       lambda_y = NA_real_,
       lambda_y_se = NA_real_,
       resultsx = resultsx,
       cv_summary = cv_summary,
       cv_folds = cv_folds,
       Lambda = final$Lambda,
       B = final$B_opt,
       fit = final
       ))
}




cca_rrr_cv_folds <- function(X, Y, Sx, Sy, kfolds=5,
                             lambda=0.01,
                             r=2,
                             solver = "ADMM",
                             mode = "sqrtm_norm",
                             standardize = FALSE,
                             preprocess = "none",
                             cv_metric = "mse",
                             rho=1,
                             LW_Sy = TRUE,
                             niter=1e4,
                             thresh_0=0,
                             matrix_free_threshold = 4000L,
                             cg_tol = 1e-6,
                             cg_maxiter = NULL,
                             thresh = 1e-4,
                             folds = NULL,
                             return_fold_values = FALSE) {
  if (is.null(folds)) {
    folds <- .create_cv_folds(nrow(Y), kfolds)
  }

  preprocess_mode <- .resolve_preprocess_mode(standardize = FALSE, preprocess = preprocess)
  mode <- .rrr_match_mode(mode)
  cv_metric <- .rrr_match_cv_metric(cv_metric)
  X <- .preprocess_matrix(X, preprocess_mode)
  Y <- .preprocess_matrix(Y, preprocess_mode)

  rmse <- vapply(seq_along(folds), function(i) {
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
      final <- cca_rrr(X_train, Y_train, 
                       Sx=Sx_train, Sy=Sy_train, highdim=TRUE,
                       mode =mode,
                       lambda=lambda, r=r, solver=solver,
                       LW_Sy=LW_Sy, standardize=FALSE, preprocess="none",
                       rho=rho, niter=niter, 
                       matrix_free_threshold = matrix_free_threshold,
                       cg_tol = cg_tol,
                       cg_maxiter = cg_maxiter,
                       thresh=thresh, thresh_0=thresh_0,
                       verbose=FALSE)
      if (cv_metric == "mse") {
        mean((X_val %*% final$U - Y_val %*% final$V)^2)
      } else {
        .rrr_correlation_score(X_val %*% final$U, Y_val %*% final$V)
      }
    }, error = function(e) {
      message("Error in fold ", i, ": ", conditionMessage(e))
      NA_real_
    })
  }, numeric(1))
  
  if (return_fold_values) return(rmse)
  if (all(is.na(rmse))) {
    return(if (cv_metric == "mse") 1e8 else NA_real_)
  }
  mean(rmse, na.rm = TRUE)
}
