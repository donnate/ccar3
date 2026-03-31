# -----------------------------
# Helper: robust chol with jitter
# -----------------------------
chol_jitter <- function(G, jitter = 1e-8, max_tries = 6) {
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

# -----------------------------
# Helper: inverse square root for PSD matrices (eigen fallback)
# -----------------------------
invsqrt_psd <- function(G, eps = 1e-8) {
  G <- (G + t(G)) / 2
  ed <- eigen(G, symmetric = TRUE)
  vals <- pmax(ed$values, eps)
  ed$vectors %*% diag(1 / sqrt(vals), nrow = length(vals)) %*% t(ed$vectors)
}

# -----------------------------
# Helper: whitening factor G^{-1/2}, with chol then eigen fallback
# -----------------------------
whiten_factor <- function(G, ridge = 1e-8, verbose = FALSE) {
  if (any(!is.finite(G))) {
    if (verbose) cat("\nNon-finite entries in whitening Gram matrix; sanitizing.\n")
    G[!is.finite(G)] <- 0
    G <- (G + t(G)) / 2 + ridge * diag(nrow(G))
  }
  tryCatch({
    R <- chol_jitter(G, jitter = ridge)
    backsolve(R, diag(nrow(G)))
  }, error = function(e) {
    if (verbose) cat("\nchol_jitter failed; using eigen whitening fallback.\n")
    invsqrt_psd(G, eps = ridge)
  })
}

# -----------------------------
# Helper: stable rank-k SVD with fallback
# -----------------------------
safe_rank_svd <- function(M, k, prefer_sparse = FALSE) {
  if (k < 1) return(NULL)
  nr <- nrow(M); nc <- ncol(M)
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

#' Efficient CCA for Two High-Dimensional Views
#'
#' Fits sparse canonical directions with an ADMM-based reduced-rank regression
#' formulation tailored to the setting where both views are high-dimensional.
#'
#' @param X Predictor matrix (n x p).
#' @param Y Response matrix (n x q).
#' @param lambda Regularization parameter.
#' @param groups Optional group structure for blockwise sparsity.
#' @param r Target rank.
#' @param standardize Whether to scale variables after centering.
#' @param rho ADMM penalty parameter.
#' @param B0 Optional warm start for the coefficient matrix.
#' @param eps Convergence tolerance for ADMM.
#' @param maxiter Maximum number of ADMM iterations.
#' @param verbose Whether to print diagnostics.
#' @param epsilon_sv Numerical threshold used to discard near-zero singular values.
#' @param ridge_whiten Ridge added when whitening Gram matrices.
#'
#' @return A list containing the estimated canonical directions, canonical
#'   correlations, the fitted coefficient matrix, preprocessing metadata, and
#'   convergence information.
#' @export
ecca <- function(
  X, Y,
  lambda = 0,
  groups = NULL,
  r = 2,
  standardize = FALSE,
  rho = 1,
  B0 = NULL,
  eps = 1e-4,
  maxiter = 500,
  verbose = TRUE,
  epsilon_sv = 1e-8,
  ridge_whiten = 1e-8
) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (nrow(X) != nrow(Y)) stop("X and Y must have the same number of rows.")

  n  <- nrow(X)
  p0 <- ncol(X)
  q0 <- ncol(Y)
  x_names <- colnames(X)
  y_names <- colnames(Y)

  if (verbose) cat("eps =", eps, " maxiter =", maxiter, "\n")

  # ---- preprocess: center, drop 0-var, optional scale
  x_center <- colMeans(X)
  y_center <- colMeans(Y)
  X <- sweep(X, 2, x_center, "-")
  Y <- sweep(Y, 2, y_center, "-")

  sx <- matrixStats::colSds(X)
  sy <- matrixStats::colSds(Y)

  keepX <- which(is.finite(sx) & sx > 0)
  keepY <- which(is.finite(sy) & sy > 0)

  dropX <- setdiff(seq_len(p0), keepX)
  dropY <- setdiff(seq_len(q0), keepY)

  X <- X[, keepX, drop = FALSE]
  Y <- Y[, keepY, drop = FALSE]

  x_scale <- rep(1, p0)
  y_scale <- rep(1, q0)
  if (standardize) {
    X <- sweep(X, 2, sx[keepX], "/")
    Y <- sweep(Y, 2, sy[keepY], "/")
    x_scale[keepX] <- sx[keepX]
    y_scale[keepY] <- sy[keepY]
  }

  p <- ncol(X)
  q <- ncol(Y)

  if (p < 1 || q < 1) {
    return(list(
      U = matrix(NA_real_, p0, r),
      V = matrix(NA_real_, q0, r),
      cor = rep(0, r),
      loss = Inf,
      Bhat = matrix(0, p0, q0),
      keepX = keepX, keepY = keepY, dropX = dropX, dropY = dropY,
      center = list(X = x_center, Y = y_center),
      scale  = list(X = x_scale,  Y = y_scale),
      converged = FALSE
    ))
  }

  # ---- map coord groups from ORIGINAL coords to filtered coords
  if (!is.null(groups)) {
    mapX <- integer(p0); mapX[keepX] <- seq_along(keepX)
    mapY <- integer(q0); mapY[keepY] <- seq_along(keepY)

    groups <- lapply(groups, function(idx) {
      if (is.matrix(idx) && ncol(idx) == 2) {
        ix <- mapX[idx[, 1]]
        iy <- mapY[idx[, 2]]
        ok <- (ix > 0) & (iy > 0)
        if (!any(ok)) return(matrix(integer(0), ncol = 2))
        cbind(ix[ok], iy[ok])
      } else {
        idx
      }
    })
  }

  # ---- linearize groups once + precompute sqrt(|g|)
  group_sqrt <- NULL
  if (!is.null(groups)) {
    groups_lin <- vector("list", length(groups))
    group_sqrt <- numeric(length(groups))
    for (g in seq_along(groups)) {
      idx <- groups[[g]]
      if (is.matrix(idx) && ncol(idx) == 2) {
        groups_lin[[g]] <- as.integer(idx[, 1] + (idx[, 2] - 1L) * p)
        group_sqrt[g] <- sqrt(nrow(idx))
      } else {
        groups_lin[[g]] <- idx
        group_sqrt[g] <- sqrt(length(idx))
      }
    }
    groups <- groups_lin
  }

  # ---- reduced bases via SVD
  EDx <- svd(X, nu = 0, nv = min(n, p))
  EDy <- svd(Y, nu = 0, nv = min(n, q))

  Ux_all <- EDx$v
  Uy_all <- EDy$v
  Lx_all <- (EDx$d^2) / n
  Ly_all <- (EDy$d^2) / n

  idx_x <- which(Lx_all > epsilon_sv)
  idx_y <- which(Ly_all > epsilon_sv)
  if (length(idx_x) < 1 || length(idx_y) < 1) {
    U_full <- matrix(0, p0, r)
    V_full <- matrix(0, q0, r)
    if (!is.null(x_names)) rownames(U_full) <- x_names
    if (!is.null(y_names)) rownames(V_full) <- y_names
    return(list(
      U = U_full, V = V_full, cor = rep(0, r), loss = Inf,
      Bhat = matrix(0, p0, q0),
      keepX = keepX, keepY = keepY, dropX = dropX, dropY = dropY,
      center = list(X = x_center, Y = y_center),
      scale  = list(X = x_scale,  Y = y_scale),
      converged = FALSE
    ))
  }

  Ux <- Ux_all[, idx_x, drop = FALSE]; Lx <- Lx_all[idx_x]
  Uy <- Uy_all[, idx_y, drop = FALSE]; Ly <- Ly_all[idx_y]
  UyT <- t(Uy)

  kx <- length(Lx)
  ky <- length(Ly)

  b <- outer(Lx, Ly) + rho

  # B1 = (X Ux)^T (Y Uy) / n
  XUx <- X %*% Ux
  YUy <- Y %*% Uy
  B1  <- crossprod(XUx, YUy) / n
  rm(XUx, YUy)

  # choose multiplication order once
  use_left_Low    <- (p * ky <= kx * q)
  use_left_Ztilde <- (kx * q <= p * ky)

  # ---- ADMM (low-memory)
  if (!is.null(B0)) {
    Z <- B0[keepX, keepY, drop = FALSE]
  } else {
    Z <- matrix(0, p, q)
  }

  # Ztilde = Ux^T Z Uy with order choice
  if (use_left_Ztilde) {
    Ztilde <- crossprod(Ux, Z) %*% Uy
  } else {
    Ztilde <- crossprod(Ux, Z %*% Uy)
  }
  Htilde <- matrix(0, kx, ky)

  converged <- FALSE
  tau_base <- lambda / rho   

  for (iter in 1:maxiter) {

    proj   <- Ztilde - Htilde
    Btilde <- (B1 + rho * proj) / b
    M      <- Btilde - proj

    # Low = Ux M Uy^T
    if (use_left_Low) {
      tmp <- Ux %*% M
      Low <- tmp %*% UyT
    } else {
      tmp <- M %*% UyT
      Low <- Ux %*% tmp
    }

    # Z_in = Low + Z   (no Z_old copy)
    Z_in <- Low + Z

    if (is.null(groups)) {
      Z <- soft_thresh(Z_in, tau_base)
    } else {
      Z <- Z_in
      for (g in seq_along(groups)) {
        idx <- groups[[g]]
        if (length(idx) == 0) next
        Z[idx] <- soft_thresh2(Z[idx], group_sqrt[g] * tau_base)
      }
    }

    # update Ztilde
    if (use_left_Ztilde) {
      Ztilde_new <- crossprod(Ux, Z) %*% Uy
    } else {
      Ztilde_new <- crossprod(Ux, Z %*% Uy)
    }

    Htilde <- Htilde + (Btilde - Ztilde_new)

    r_norm <- sqrt(sum((Btilde - Ztilde_new)^2))
    s_norm <- rho * sqrt(sum((Ztilde_new - Ztilde)^2))

    eps_pri  <- eps * (1 + max(sqrt(sum(Btilde^2)), sqrt(sum(Ztilde_new^2))))
    eps_dual <- eps * (1 + sqrt(sum(Htilde^2)))

    if (verbose && (iter %% 50 == 0)) {
      cat("\niter:", iter, " r:", signif(r_norm, 3), " s:", signif(s_norm, 3))
    }

    Ztilde <- Ztilde_new
    if (r_norm < eps_pri && s_norm < eps_dual) {
      converged <- TRUE
      break
    }
  }

  if (verbose) {
    if (!converged) cat("\nADMM did not converge (hit maxiter).\n")
    else cat("\nADMM converged.\n")
  }

  # ---- postprocessing (small SVD)
  r_eff <- min(r, nrow(Ztilde), ncol(Ztilde))
  if (r_eff < 1) {
    U_full <- matrix(0, p0, r)
    V_full <- matrix(0, q0, r)
    if (!is.null(x_names)) rownames(U_full) <- x_names
    if (!is.null(y_names)) rownames(V_full) <- y_names
    Bhat_full <- matrix(0, p0, q0); Bhat_full[keepX, keepY] <- Z
    return(list(
      U = U_full, V = V_full, cor = rep(0, r), loss = Inf,
      Bhat = Bhat_full,
      keepX = keepX, keepY = keepY, dropX = dropX, dropY = dropY,
      center = list(X = x_center, Y = y_center),
      scale  = list(X = x_scale,  Y = y_scale),
      converged = converged
    ))
  }

  # Csmall <- Ztilde
  # Csmall <- sweep(Csmall, 1, sqrt(Lx), "*")
  # Csmall <- sweep(Csmall, 2, sqrt(Ly), "*")

  # --- SVD on sparse Z (preserves row sparsity) ---
  SVDs <- safe_rank_svd(Z, r_eff, prefer_sparse = TRUE)
  if (is.null(SVDs)) {
    if (verbose) cat("\nSVD failed in ecca(). Returning NA directions.\n")
    U_full <- matrix(NA_real_, p0, r)
    V_full <- matrix(NA_real_, q0, r)
    if (!is.null(x_names)) rownames(U_full) <- x_names
    if (!is.null(y_names)) rownames(V_full) <- y_names
    Bhat_full <- matrix(0, p0, q0)
    Bhat_full[keepX, keepY] <- Z
    return(list(
      U = U_full,
      V = V_full,
      cor = rep(0, r),
      loss = Inf,
      Bhat = Bhat_full,
      keepX = keepX, keepY = keepY, dropX = dropX, dropY = dropY,
      center = list(X = x_center, Y = y_center),
      scale  = list(X = x_scale,  Y = y_scale),
      converged = FALSE
    ))
  }

  U0 <- SVDs$u
  V0 <- SVDs$v

  # --- whitening ---
  XU0 <- X %*% U0
  YV0 <- Y %*% V0

  GX <- crossprod(XU0) / n + ridge_whiten * diag(r_eff)
  GY <- crossprod(YV0) / n + ridge_whiten * diag(r_eff)

  WGX <- tryCatch(
    whiten_factor(GX, ridge = ridge_whiten, verbose = verbose),
    error = function(e) NULL
  )
  WGY <- tryCatch(
    whiten_factor(GY, ridge = ridge_whiten, verbose = verbose),
    error = function(e) NULL
  )
  if (is.null(WGX) || is.null(WGY)) {
    if (verbose) cat("\nWhitening failed in ecca(). Returning NA directions.\n")
    U_full <- matrix(NA_real_, p0, r)
    V_full <- matrix(NA_real_, q0, r)
    if (!is.null(x_names)) rownames(U_full) <- x_names
    if (!is.null(y_names)) rownames(V_full) <- y_names
    Bhat_full <- matrix(0, p0, q0)
    Bhat_full[keepX, keepY] <- Z
    return(list(
      U = U_full,
      V = V_full,
      cor = rep(0, r),
      loss = Inf,
      Bhat = Bhat_full,
      keepX = keepX, keepY = keepY, dropX = dropX, dropY = dropY,
      center = list(X = x_center, Y = y_center),
      scale  = list(X = x_scale,  Y = y_scale),
      converged = FALSE
    ))
  }

  U <- U0 %*% WGX
  V <- V0 %*% WGY

  XU <- X %*% U
  YV <- Y %*% V

  cor <- diag(crossprod(XU, YV) / n)

  # if (requireNamespace("RSpectra", quietly = TRUE) && r_eff < min(dim(Csmall))) {
  #   SVDs <- RSpectra::svds(Csmall, r_eff)
  # } else {
  #   SVDs <- svd(Csmall, nu = r_eff, nv = r_eff)
  # }

  # U0 <- Ux %*% SVDs$u
  # V0 <- Uy %*% SVDs$v

  # GX <- crossprod(X %*% U0) / n + ridge_whiten * diag(r_eff)
  # GY <- crossprod(Y %*% V0) / n + ridge_whiten * diag(r_eff)

  # Rx <- chol_jitter(GX, jitter = ridge_whiten)
  # Ry <- chol_jitter(GY, jitter = ridge_whiten)

  # U <- U0 %*% backsolve(Rx, diag(r_eff))
  # V <- V0 %*% backsolve(Ry, diag(r_eff))

  # XU <- X %*% U
  # YV <- Y %*% V
  # cor <- diag(crossprod(XU, YV) / n)

  neg <- cor < 0
  if (any(neg)) {
    V[, neg] <- -V[, neg, drop = FALSE]
    cor[neg] <- -cor[neg]
    YV[, neg] <- -YV[, neg, drop = FALSE]
  }

  ord <- order(cor, decreasing = TRUE)
  cor <- cor[ord]
  U <- U[, ord, drop = FALSE]
  V <- V[, ord, drop = FALSE]
  XU <- XU[, ord, drop = FALSE]
  YV <- YV[, ord, drop = FALSE]

  loss <- sum((XU - YV)^2) / n

  # ---- inflate to original coords
  U_full <- matrix(0, p0, r_eff)
  V_full <- matrix(0, q0, r_eff)
  U_full[keepX, ] <- U
  V_full[keepY, ] <- V

  if (!is.null(x_names)) rownames(U_full) <- x_names
  if (!is.null(y_names)) rownames(V_full) <- y_names
  colnames(U_full) <- paste0("comp", seq_len(r_eff))
  colnames(V_full) <- paste0("comp", seq_len(r_eff))

  Bhat_full <- matrix(0, p0, q0)
  Bhat_full[keepX, keepY] <- Z

  list(
    U = U_full,
    V = V_full,
    Ux = Ux, Uy = Uy, Lx = Lx, Ly = Ly,
    cor = cor,
    loss = loss,
    Bhat = Bhat_full,
    keepX = keepX, keepY = keepY, dropX = dropX, dropY = dropY,
    center = list(X = x_center, Y = y_center),
    scale  = list(X = x_scale,  Y = y_scale),
    converged = converged
  )
}

# -----------------------------
# ecca_across_lambdas (lambda path, CV-friendly)
# -----------------------------
ecca_across_lambdas <- function(
  X, Y, lambdas = 0, groups = NULL, r = 2,
  Sx = NULL, Sy = NULL, Sxy = NULL,                 # kept for compatibility; ignored
  standardize = TRUE,
  rho = 1, B0 = NULL, eps = 1e-4, maxiter = 500, verbose = TRUE,
  dense = TRUE, optimized = FALSE,                  # kept for compatibility; ignored
  epsilon_sv = 1e-8, ridge_whiten = 1e-8,
  warm_start = TRUE,
  lambda_order = c("decreasing", "increasing", "given"),
  log_prefix = "",
  admm_print_every = 50L,
  lambda_print = TRUE,
  collect_logs = FALSE,
  preprocess = TRUE,
  return_uv = TRUE,
  X_val = NULL, Y_val = NULL,
  scoring_method = c("mse", "trace")
){
  scoring_method <- match.arg(scoring_method)
  lambda_order <- match.arg(lambda_order)

  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (nrow(X) != nrow(Y)) stop("X and Y must have same nrow().")
  logs <- character(0)

  lambdas <- as.numeric(lambdas)
  L <- length(lambdas)
  if (L < 1) stop("lambdas must have length >= 1")

  ord <- switch(
    lambda_order,
    decreasing = order(lambdas, decreasing = TRUE),
    increasing = order(lambdas, decreasing = FALSE),
    given      = seq_along(lambdas)
  )
  inv_ord <- order(ord)
  lambdas_run <- lambdas[ord]

  p0 <- ncol(X); q0 <- ncol(Y)
  x_names <- colnames(X); y_names <- colnames(Y)

  keepX <- seq_len(p0); keepY <- seq_len(q0)
  centerX <- rep(0, p0); centerY <- rep(0, q0)
  scaleX  <- rep(1, p0); scaleY  <- rep(1, q0)

  if (preprocess) {
    centerX <- colMeans(X); centerY <- colMeans(Y)
    X <- sweep(X, 2, centerX, "-")
    Y <- sweep(Y, 2, centerY, "-")

    sx <- matrixStats::colSds(X)
    sy <- matrixStats::colSds(Y)

    keepX <- which(is.finite(sx) & sx > 0)
    keepY <- which(is.finite(sy) & sy > 0)

    X <- X[, keepX, drop = FALSE]
    Y <- Y[, keepY, drop = FALSE]

    if (standardize) {
      X <- sweep(X, 2, sx[keepX], "/")
      Y <- sweep(Y, 2, sy[keepY], "/")
      scaleX[keepX] <- sx[keepX]
      scaleY[keepY] <- sy[keepY]
    }

    # map coord groups (original -> filtered)
    if (!is.null(groups)) {
      mapX <- integer(p0); mapX[keepX] <- seq_along(keepX)
      mapY <- integer(q0); mapY[keepY] <- seq_along(keepY)

      groups <- lapply(groups, function(idx) {
        if (is.matrix(idx) && ncol(idx) == 2) {
          ix <- mapX[idx[, 1]]
          iy <- mapY[idx[, 2]]
          ok <- (ix > 0) & (iy > 0)
          if (!any(ok)) return(matrix(integer(0), ncol = 2))
          cbind(ix[ok], iy[ok])
        } else {
          idx
        }
      })
    }

    # preprocess validation too
    if (!is.null(X_val) && !is.null(Y_val)) {
      if (!is.matrix(X_val)) X_val <- as.matrix(X_val)
      if (!is.matrix(Y_val)) Y_val <- as.matrix(Y_val)
      X_val <- sweep(X_val, 2, centerX, "-")
      Y_val <- sweep(Y_val, 2, centerY, "-")
      X_val <- X_val[, keepX, drop = FALSE]
      Y_val <- Y_val[, keepY, drop = FALSE]
      if (standardize) {
        X_val <- sweep(X_val, 2, scaleX[keepX], "/")
        Y_val <- sweep(Y_val, 2, scaleY[keepY], "/")
      }
    }
  } else {
    if (!is.null(X_val) && !is.null(Y_val)) {
      if (!is.matrix(X_val)) X_val <- as.matrix(X_val)
      if (!is.matrix(Y_val)) Y_val <- as.matrix(Y_val)
    }
  }

  n <- nrow(X)
  p <- ncol(X); q <- ncol(Y)

  # linearize groups once + precompute sqrt(|g|)
  group_sqrt <- NULL
  if (!is.null(groups)) {
    groups_lin <- vector("list", length(groups))
    group_sqrt <- numeric(length(groups))
    for (g in seq_along(groups)) {
      idx <- groups[[g]]
      if (is.matrix(idx) && ncol(idx) == 2) {
        groups_lin[[g]] <- as.integer(idx[, 1] + (idx[, 2] - 1L) * p)
        group_sqrt[g] <- sqrt(nrow(idx))
      } else {
        groups_lin[[g]] <- idx
        group_sqrt[g] <- sqrt(length(idx))
      }
    }
    groups <- groups_lin
  }

  if (p < 1 || q < 1) {
    out_scores <- rep(Inf, L)
    if (!return_uv) return(list(scores = out_scores[inv_ord], lambdas = lambdas))
    Uout <- lapply(seq_len(L), function(.) matrix(NA_real_, p0, r))
    Vout <- lapply(seq_len(L), function(.) matrix(NA_real_, q0, r))
    return(list(U = Uout, V = Vout, scores = out_scores[inv_ord], lambdas = lambdas))
  }

  # reduced bases
  EDx <- svd(X, nu = 0, nv = min(n, p))
  EDy <- svd(Y, nu = 0, nv = min(n, q))

  Ux_all <- EDx$v
  Uy_all <- EDy$v
  Lx_all <- (EDx$d^2) / n
  Ly_all <- (EDy$d^2) / n

  idx_x <- which(Lx_all > epsilon_sv)
  idx_y <- which(Ly_all > epsilon_sv)
  if (length(idx_x) < 1 || length(idx_y) < 1) {
    out_scores <- rep(Inf, L)
    if (!return_uv) return(list(scores = out_scores[inv_ord], lambdas = lambdas))
    Uout <- lapply(seq_len(L), function(.) matrix(NA_real_, p0, r))
    Vout <- lapply(seq_len(L), function(.) matrix(NA_real_, q0, r))
    return(list(U = Uout, V = Vout, scores = out_scores[inv_ord], lambdas = lambdas))
  }

  Ux <- Ux_all[, idx_x, drop = FALSE]; Lx <- Lx_all[idx_x]
  Uy <- Uy_all[, idx_y, drop = FALSE]; Ly <- Ly_all[idx_y]
  UyT <- t(Uy)

  kx <- length(Lx); ky <- length(Ly)
  b  <- outer(Lx, Ly) + rho

  # B1 once
  XUx <- X %*% Ux
  YUy <- Y %*% Uy
  B1  <- crossprod(XUx, YUy) / n
  rm(XUx, YUy)

  use_left_Low    <- (p * ky <= kx * q)
  use_left_Ztilde <- (kx * q <= p * ky)

  compute_Ztilde <- function(Zmat) {
    if (use_left_Ztilde) crossprod(Ux, Zmat) %*% Uy
    else crossprod(Ux, Zmat %*% Uy)
  }
  compute_Low <- function(M) {
    if (use_left_Low) (Ux %*% M) %*% UyT
    else Ux %*% (M %*% UyT)
  }

  out_scores <- rep(NA_real_, L)
  out_cor    <- vector("list", L)
  out_loss   <- rep(NA_real_, L)
  out_U      <- if (return_uv) vector("list", L) else NULL
  out_V      <- if (return_uv) vector("list", L) else NULL

  if (!is.null(B0)) {
    Z <- if (preprocess) B0[keepX, keepY, drop = FALSE] else B0
  } else {
    Z <- matrix(0, p, q)
  }
  admm_print_every <- as.integer(admm_print_every)
  if (!is.finite(admm_print_every) || admm_print_every < 1L) admm_print_every <- 50L

  for (t in seq_len(L)) {
    lam <- lambdas_run[t]
    t_lambda0 <- Sys.time()
    if (verbose && lambda_print) {
      cat(sprintf("\nStarting lambda %d/%d: %g", t, L, lam))
      cat(sprintf("\n%s lambda %d/%d = %g", log_prefix, t, L, lam))
      utils::flush.console()
      if (collect_logs) {
        logs <- c(logs,
                  sprintf("Starting lambda %d/%d: %g", t, L, lam),
                  sprintf("%s lambda %d/%d = %g", log_prefix, t, L, lam))
      }
    }


    if (!warm_start && t > 1) Z <- matrix(0, p, q)

    Ztilde <- compute_Ztilde(Z)
    Htilde <- matrix(0, kx, ky)

    tau_base <- lam / rho
    converged <- FALSE

    for (iter in seq_len(maxiter)) {
      proj   <- Ztilde - Htilde
      Btilde <- (B1 + rho * proj) / b
      M      <- Btilde - proj

      Low <- compute_Low(M)

      Z_in <- Low + Z

      if (is.null(groups)) {
        Z <- soft_thresh(Z_in, tau_base)
      } else {
        Z <- Z_in
        for (g in seq_along(groups)) {
          idx <- groups[[g]]
          if (length(idx) == 0) next
          Z[idx] <- soft_thresh2(Z[idx], group_sqrt[g] * tau_base)
        }
      }

      Ztilde_new <- compute_Ztilde(Z)
      Htilde <- Htilde + (Btilde - Ztilde_new)

      r_norm <- sqrt(sum((Btilde - Ztilde_new)^2))
      s_norm <- rho * sqrt(sum((Ztilde_new - Ztilde)^2))

      eps_pri  <- eps * (1 + max(sqrt(sum(Btilde^2)), sqrt(sum(Ztilde_new^2))))
      eps_dual <- eps * (1 + sqrt(sum(Htilde^2)))

      if (verbose && iter %% 100 == 0) {
        cat("\nlambda", signif(lam, 3), "iter", iter,
            "r", signif(r_norm, 3), "s", signif(s_norm, 3))
        if (collect_logs) {
          logs <- c(logs, sprintf("lambda %s iter %d r %s s %s",
                                  signif(lam, 3), iter, signif(r_norm, 3), signif(s_norm, 3)))
        }
      }
          # inside ADMM iter loop
      if (verbose && (iter %% admm_print_every == 0L)) {
        cat(sprintf("\n%s   iter %d  r=%g  s=%g", log_prefix, iter, r_norm, s_norm))
        utils::flush.console()
        if (collect_logs) {
          logs <- c(logs, sprintf("%s   iter %d  r=%g  s=%g", log_prefix, iter, r_norm, s_norm))
        }
      }

      Ztilde <- Ztilde_new
      if (r_norm < eps_pri && s_norm < eps_dual) {
        converged <- TRUE
        break
      }
    }

    # postprocessing small SVD
    Csmall <- Ztilde
    Csmall <- sweep(Csmall, 1, sqrt(Lx), "*")
    Csmall <- sweep(Csmall, 2, sqrt(Ly), "*")

    r_eff <- min(r, nrow(Csmall), ncol(Csmall))
    if (r_eff < 1) {
      out_scores[t] <- Inf
      out_loss[t]   <- Inf
      out_cor[[t]]  <- rep(0, r)
      if (return_uv) {
        out_U[[t]] <- matrix(0, p0, r)
        out_V[[t]] <- matrix(0, q0, r)
      }
      next
    }

    SVDs <- safe_rank_svd(Z, r_eff, prefer_sparse = TRUE)
    if (is.null(SVDs)) {
      out_scores[t] <- Inf
      out_loss[t]   <- Inf
      out_cor[[t]]  <- rep(0, r)
      if (return_uv) {
        out_U[[t]] <- matrix(0, p0, r)
        out_V[[t]] <- matrix(0, q0, r)
      }
      if (collect_logs) {
        logs <- c(logs, sprintf("%s SVD failed at lambda=%g; assigned Inf score", log_prefix, lam))
      }
      next
    }

    Uhat <- SVDs$u
    Vhat <- SVDs$v

    # --- whitening ---
    XU0 <- X %*% Uhat
    YV0 <- Y %*% Vhat

    GX <- crossprod(XU0) / n + ridge_whiten * diag(r_eff)
    GY <- crossprod(YV0) / n + ridge_whiten * diag(r_eff)

    WGX <- tryCatch(
      whiten_factor(GX, ridge = ridge_whiten, verbose = FALSE),
      error = function(e) NULL
    )
    WGY <- tryCatch(
      whiten_factor(GY, ridge = ridge_whiten, verbose = FALSE),
      error = function(e) NULL
    )
    if (is.null(WGX) || is.null(WGY)) {
      out_scores[t] <- Inf
      out_loss[t]   <- Inf
      out_cor[[t]]  <- rep(0, r)
      if (return_uv) {
        out_U[[t]] <- matrix(0, p0, r)
        out_V[[t]] <- matrix(0, q0, r)
      }
      if (collect_logs) {
        logs <- c(logs, sprintf("%s whitening failed at lambda=%g; assigned Inf score", log_prefix, lam))
      }
      next
    }

    Uhat <- Uhat %*% WGX
    Vhat <- Vhat %*% WGY

    XU <- X %*% Uhat
    YV <- Y %*% Vhat

    cor <- diag(crossprod(XU, YV) / n)

    neg <- cor < 0
    if (any(neg)) {
      Vhat[, neg] <- -Vhat[, neg, drop = FALSE]
      cor[neg] <- -cor[neg]
      YV[, neg] <- -YV[, neg, drop = FALSE]
    }

    ord_cor <- order(cor, decreasing = TRUE)
    cor <- cor[ord_cor]
    Uhat <- Uhat[, ord_cor, drop = FALSE]
    Vhat <- Vhat[, ord_cor, drop = FALSE]
    XU <- XU[, ord_cor, drop = FALSE]
    YV <- YV[, ord_cor, drop = FALSE]

    if (!is.null(X_val) && !is.null(Y_val)) {
      XU_val <- X_val %*% Uhat
      YV_val <- Y_val %*% Vhat
      if (scoring_method == "mse") {
        out_scores[t] <- mean((XU_val - YV_val)^2)
      } else {
        # Faster than cor(): compute only paired (diagonal) correlations.
        XUc <- scale(XU_val, center = TRUE, scale = FALSE)
        YVc <- scale(YV_val, center = TRUE, scale = FALSE)
        num <- colSums(XUc * YVc)
        den <- sqrt(colSums(XUc^2) * colSums(YVc^2))
        den[den < 1e-12 | !is.finite(den)] <- 1e-12
        corr_diag <- num / den
        out_scores[t] <- -mean(corr_diag)
      }
    } else {
      out_scores[t] <- NA_real_
    }

    out_loss[t]  <- sum((XU - YV)^2) / n
    out_cor[[t]] <- cor

    if (return_uv) {
      U_full <- matrix(0, p0, r_eff)
      V_full <- matrix(0, q0, r_eff)
      U_full[keepX, ] <- Uhat
      V_full[keepY, ] <- Vhat
      if (!is.null(x_names)) rownames(U_full) <- x_names
      if (!is.null(y_names)) rownames(V_full) <- y_names
      out_U[[t]] <- U_full
      out_V[[t]] <- V_full
    }
    if (verbose && lambda_print) {
      cat(sprintf("\nFinished lambda %d/%d in %.1fs",
                  t, L, as.numeric(difftime(Sys.time(), t_lambda0, units="secs"))))
      if (collect_logs) {
        logs <- c(logs, sprintf("Finished lambda %d/%d in %.1fs",
                                t, L, as.numeric(difftime(Sys.time(), t_lambda0, units="secs"))))
      }
    }
  }

  # restore original lambda order
  out_scores <- out_scores[inv_ord]
  out_loss   <- out_loss[inv_ord]
  out_cor    <- out_cor[inv_ord]
  if (return_uv) {
    out_U <- out_U[inv_ord]
    out_V <- out_V[inv_ord]
  }

  if (!return_uv) {
    out <- list(scores = out_scores, lambdas = lambdas)
    if (collect_logs) out$logs <- logs
    return(out)
  }

  if (length(lambdas) == 1) {
    out <- list(U = out_U[[1]], V = out_V[[1]],
                cor = out_cor[[1]], loss = out_loss[1],
                scores = out_scores, lambdas = lambdas,
                center = list(X = centerX, Y = centerY),
                scale  = list(X = scaleX,  Y = scaleY),
                keepX = keepX, keepY = keepY)
    if (collect_logs) out$logs <- logs
    return(out)
  }

  out <- list(U = out_U, V = out_V,
              cor = out_cor, loss = out_loss,
              scores = out_scores, lambdas = lambdas,
              center = list(X = centerX, Y = centerY),
              scale  = list(X = scaleX,  Y = scaleY),
              keepX = keepX, keepY = keepY,
              Ux = Ux, Uy = Uy, Lx = Lx, Ly = Ly)
  if (collect_logs) out$logs <- logs
  out
}

# -----------------------------
# ecca.eval (CV) - with seamless group mapping + faster parallel behavior
# -----------------------------
ecca.eval <- function(
  X, Y, lambdas = 0, groups = NULL, r = 2,
  standardize = TRUE,
  rho = 1, B0 = NULL, nfold = 5,
  eps = 1e-4, maxiter = 500, verbose = TRUE,
  parallel = TRUE, nb_cores = NULL, set_seed_cv = NULL,
  scoring_method = c("mse", "trace"),
  cv_use_median = FALSE,
  dense = TRUE, optimized = FALSE,
  epsilon_sv = 1e-8, ridge_whiten = 1e-8
){
  scoring_method <- match.arg(scoring_method)
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (nrow(X) != nrow(Y)) stop("X and Y must have same nrow().")

  n <- nrow(X)
  if (n < 2 * nfold) {
    if (verbose) cat("\nWarning: n too small for CV; using first lambda.\n")
    scores <- data.frame(lambda = lambdas, mse = NA_real_, se = NA_real_)
    return(list(scores = scores, lambda.min = lambdas[1], lambda.1se = lambdas[1]))
  }

  # global preprocessing
  p0 <- ncol(X); q0 <- ncol(Y)

  X <- scale(X, center = TRUE, scale = FALSE)
  Y <- scale(Y, center = TRUE, scale = FALSE)

  sx <- matrixStats::colSds(X)
  sy <- matrixStats::colSds(Y)
  keepX <- which(is.finite(sx) & sx > 0)
  keepY <- which(is.finite(sy) & sy > 0)

  X <- X[, keepX, drop = FALSE]
  Y <- Y[, keepY, drop = FALSE]

  if (standardize) {
    X <- sweep(X, 2, sx[keepX], "/")
    Y <- sweep(Y, 2, sy[keepY], "/")
  }

  n <- nrow(X)

  # IMPORTANT: map groups to the filtered coordinates ONCE (since preprocess=FALSE in folds)
  if (!is.null(groups)) {
    mapX <- integer(p0); mapX[keepX] <- seq_along(keepX)
    mapY <- integer(q0); mapY[keepY] <- seq_along(keepY)

    groups <- lapply(groups, function(idx) {
      if (is.matrix(idx) && ncol(idx) == 2) {
        ix <- mapX[idx[, 1]]
        iy <- mapY[idx[, 2]]
        ok <- (ix > 0) & (iy > 0)
        if (!any(ok)) return(matrix(integer(0), ncol = 2))
        cbind(ix[ok], iy[ok])
      } else {
        idx
      }
    })
  }

  if (!is.null(set_seed_cv)) set.seed(set_seed_cv)
  folds <- .create_cv_folds(n, nfold)

  # parallel setup
  if (parallel) {
    cores <- if (is.null(nb_cores)) max(1L, parallel::detectCores() - 1L) else as.integer(nb_cores)
    cores <- min(cores, length(folds))
    Sys.setenv(
      OMP_NUM_THREADS = 1,
      OPENBLAS_NUM_THREADS = 1,
      MKL_NUM_THREADS = 1,
      VECLIB_MAXIMUM_THREADS = 1
    )

    cl <- setup_parallel_backend(cores, verbose = verbose)
    if (is.null(cl)) {
      warning("All parallel setup attempts failed. Proceeding in serial mode.", immediate. = TRUE)
      parallel <- FALSE
    }
  }

  if (parallel) {
    initialize_parallel_workers(cl, verbose = verbose)
    on.exit(cleanup_parallel_backend(cl), add = TRUE)
    if (verbose) {
      cat("\nRunning CV in parallel on", cores, "cores\n")
    }
    ecca_across_lambdas_current <- ecca_across_lambdas
    safe_rank_svd_current <- safe_rank_svd
    whiten_factor_current <- whiten_factor
    chol_jitter_current <- chol_jitter
    invsqrt_psd_current <- invsqrt_psd
    soft_thresh_current <- soft_thresh
    soft_thresh2_current <- soft_thresh2
    safe_rank_svd <- safe_rank_svd_current
    whiten_factor <- whiten_factor_current
    chol_jitter <- chol_jitter_current
    invsqrt_psd <- invsqrt_psd_current
    soft_thresh <- soft_thresh_current
    soft_thresh2 <- soft_thresh2_current
    collect_logs_parallel <- isTRUE(verbose)
    ecca_parallel_fold_worker <- function(fold_id) {
      tryCatch({
        fold <- folds[[fold_id]]
        pid <- Sys.getpid()
        prefix <- sprintf("[pid=%d fold=%d]", pid, fold_id)

        tr <- setdiff(seq_len(n), fold)
        Xtr <- X[tr, , drop = FALSE]
        Ytr <- Y[tr, , drop = FALSE]
        Xva <- X[fold, , drop = FALSE]
        Yva <- Y[fold, , drop = FALSE]

        fit <- ecca_across_lambdas_current(
          Xtr, Ytr, lambdas = lambdas, groups = groups, r = r,
          standardize = FALSE, rho = rho, B0 = B0,
          eps = eps, maxiter = maxiter,
          verbose = collect_logs_parallel,
          log_prefix = prefix,
          admm_print_every = 50L,
          lambda_print = collect_logs_parallel,
          epsilon_sv = epsilon_sv, ridge_whiten = ridge_whiten,
          preprocess = FALSE,
          return_uv = FALSE,
          collect_logs = collect_logs_parallel,
          X_val = Xva, Y_val = Yva,
          scoring_method = scoring_method
        )

        list(
          fold_id = fold_id,
          scores = fit$scores,
          logs = fit$logs
        )
      }, error = function(e) e)
    }

    parallel::clusterExport(
      cl,
      varlist = c(
        "X", "Y", "folds", "n",
        "lambdas", "groups", "r",
        "rho", "B0", "eps", "maxiter",
        "epsilon_sv", "ridge_whiten", "scoring_method",
        "collect_logs_parallel",
        "ecca_across_lambdas_current",
        "safe_rank_svd",
        "whiten_factor",
        "chol_jitter",
        "invsqrt_psd",
        "soft_thresh",
        "soft_thresh2",
        "ecca_parallel_fold_worker"
      ),
      envir = environment()
    )

    cv_list <- parallel::parLapplyLB(cl, seq_along(folds), ecca_parallel_fold_worker)

    if (verbose) {
      ok_folds <- cv_list[!vapply(cv_list, inherits, logical(1), "error")]
      if (length(ok_folds) > 0) {
        for (res in ok_folds) {
          if (!is.null(res$logs) && length(res$logs) > 0) {
            cat("\n")
            cat(paste(res$logs, collapse = "\n"))
            cat("\n")
          }
        }
        utils::flush.console()
      }
    }

    is_err <- vapply(cv_list, inherits, logical(1), "error")

    if (all(is_err)) {
      msgs <- unique(vapply(cv_list, function(e) {
        msg <- tryCatch(conditionMessage(e), error = function(...) "")
        if (!nzchar(msg)) msg <- paste(class(e), collapse = "/")
        msg
      }, character(1)))
      stop("All CV folds failed. Example error(s):\n- ",
           paste(utils::head(msgs, 5), collapse = "\n- "))
    }

    if (any(is_err)) warning(sum(is_err), " folds failed.")
    cv_list <- lapply(cv_list[!is_err], function(x) x$scores)
    if (length(cv_list) == 0) stop("All CV folds failed.")
    scores.cv <- do.call(cbind, cv_list)
    n_success <- ncol(scores.cv)

    

  } else {
    cv_list <- lapply(folds, function(fold) {
      tryCatch({
        tr <- setdiff(seq_len(n), fold)
        Xtr <- X[tr, , drop = FALSE]
        Ytr <- Y[tr, , drop = FALSE]
        Xva <- X[fold, , drop = FALSE]
        Yva <- Y[fold, , drop = FALSE]

        fit <- ecca_across_lambdas(
          Xtr, Ytr, lambdas = lambdas, groups = groups, r = r,
          standardize = FALSE, rho = rho, B0 = B0,
          eps = eps, maxiter = maxiter, verbose = TRUE,
          epsilon_sv = epsilon_sv, ridge_whiten = ridge_whiten,
          preprocess = FALSE,
          return_uv = FALSE,
          X_val = Xva, Y_val = Yva,
          scoring_method = scoring_method
        )
        fit$scores
      }, error = function(e) e)
    })

    is_err <- vapply(cv_list, inherits, logical(1), "error")
    if (any(is_err)) warning(sum(is_err), " folds failed.")
    cv_list <- cv_list[!is_err]
    if (length(cv_list) == 0) stop("All CV folds failed.")
    scores.cv <- do.call(cbind, cv_list)
    n_success <- ncol(scores.cv)
  }

  if (!cv_use_median) {
    mse <- rowMeans(scores.cv, na.rm = TRUE)
  } else {
    mse <- apply(scores.cv, 1, median, na.rm = TRUE)
  }
  se <- matrixStats::rowSds(scores.cv, na.rm = TRUE) / sqrt(n_success)
  se[!is.finite(se)] <- NA_real_
  scores <- data.frame(lambda = lambdas, mse = mse, se = se)

  lambda.min <- scores$lambda[which.min(scores$mse)]
  upper <- scores$mse[which.min(scores$mse)] + scores$se[which.min(scores$mse)]
  lambda.1se <- max(scores$lambda[scores$mse <= upper])

  list(scores = scores, lambda.min = lambda.min, lambda.1se = lambda.1se)
}

#' Cross-Validated Efficient CCA
#'
#' Selects a regularization parameter for [ecca()] by cross-validation and
#' refits the final model at the selected value.
#'
#' @param X Predictor matrix (n x p).
#' @param Y Response matrix (n x q).
#' @param lambdas Candidate regularization values.
#' @param groups Optional group structure for blockwise sparsity.
#' @param r Target rank.
#' @param standardize Whether to scale variables after centering.
#' @param rho ADMM penalty parameter.
#' @param B0 Optional warm start for the coefficient matrix.
#' @param nfold Number of cross-validation folds.
#' @param select Selection rule for the final lambda. One of `"lambda.min"` or `"lambda.1se"`.
#' @param eps Convergence tolerance for the final ADMM refit.
#' @param maxiter Maximum iterations for the final ADMM refit.
#' @param verbose Whether to print diagnostics.
#' @param maxiter_cv Maximum iterations used inside the cross-validation fits.
#' @param parallel Whether to parallelize cross-validation.
#' @param nb_cores Number of worker processes to use when `parallel = TRUE`.
#' @param set_seed_cv Optional random seed for fold generation.
#' @param scoring_method Cross-validation score to optimize. One of `"mse"` or `"trace"`.
#' @param cv_use_median Whether to aggregate fold scores with the median instead of the mean.
#' @param dense Retained for backward compatibility.
#' @param optimized Retained for backward compatibility.
#' @param epsilon_sv Numerical threshold used to discard near-zero singular values.
#' @param ridge_whiten Ridge added when whitening Gram matrices.
#'
#' @return A list with the final fit, selected lambda, and cross-validation
#'   scores when more than one lambda is supplied.
#' @export
ecca.cv <- function(
  X, Y, lambdas = 0, groups = NULL, r = 2, standardize = FALSE,
  rho = 1, B0 = NULL, nfold = 5, select = "lambda.min",
  eps = 1e-3, maxiter = 1000, verbose = FALSE, maxiter_cv = 300,
  parallel = FALSE, nb_cores = NULL, set_seed_cv = NULL,
  scoring_method = c("mse", "trace"), cv_use_median = FALSE,
  dense = TRUE, optimized = FALSE,
  epsilon_sv = 1e-8, ridge_whiten = 1e-8
){
  scoring_method <- match.arg(scoring_method)

  eval <- NULL
  if (length(lambdas) > 1) {
    eval <- ecca.eval(
      X, Y, lambdas = lambdas, groups = groups, r = r,
      standardize = standardize,
      rho = rho, B0 = B0, nfold = nfold,
      eps = eps, maxiter = maxiter_cv, verbose = verbose,
      parallel = parallel, nb_cores = nb_cores, set_seed_cv = set_seed_cv,
      scoring_method = scoring_method, cv_use_median = cv_use_median,
      epsilon_sv = epsilon_sv, ridge_whiten = ridge_whiten
    )
    lambda.opt <- if (select == "lambda.1se") eval$lambda.1se else eval$lambda.min
  } else {
    lambda.opt <- as.numeric(lambdas)
  }

  if (verbose) cat("\nselected lambda:", lambda.opt, "\n")

  final_verbose <- TRUE
  cat("\nRefitting at selected lambda with ADMM logs...\n")
  utils::flush.console()

  fit <- ecca_across_lambdas(
    X, Y, lambdas = lambda.opt, groups = groups, r = r,
    standardize = standardize,
    rho = rho, B0 = B0, eps = eps, maxiter = maxiter, verbose = final_verbose,
    admm_print_every = 1L,
    lambda_print = TRUE,
    epsilon_sv = epsilon_sv, ridge_whiten = ridge_whiten,
    preprocess = TRUE,
    return_uv = TRUE
  )
  cat("\nFinal ADMM refit completed.\n")
  utils::flush.console()

  out <- list(
    U = fit$U,
    V = fit$V,
    cor = fit$cor,
    loss = fit$loss,
    lambda.opt = lambda.opt
  )
  if (!is.null(eval)) out$cv.scores <- eval$scores
  out
}
