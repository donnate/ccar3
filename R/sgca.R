# Optional fast matmul if SMUT is installed
matmul <- function(A, B) {
  if (requireNamespace("SMUT", quietly = TRUE)) {
    SMUT::eigenMapMatMult(A, B)
  } else {
    A %*% B
  }
}

soft_threshold <- function(X, T) {
  # Elementwise soft-thresholding allowing scalar or matrix thresholds
  sign(X) * pmax(abs(X) - T, 0)
}

# Utility: symmetric PD inverse square-root (for normalization)
sym_inv_sqrt <- function(S, eps = 1e-10) {
  ee <- eigen((S + t(S))/2, symmetric = TRUE)
  vals <- pmax(ee$values, eps)
  V <- ee$vectors
  matmul(matmul(V, diag(1/sqrt(vals), nrow = length(vals))), t(V))
}

# eigs for top-r symmetric (RSpectra if available)
top_eigs_sym <- function(A, r) {
  if (requireNamespace("RSpectra", quietly = TRUE) && r < nrow(A)) {
    out <- RSpectra::eigs_sym(A, k = r, which = "LM")
    list(values = Re(out$values), vectors = Re(out$vectors))
  } else {
    ev <- eigen((A + t(A))/2, symmetric = TRUE)
    list(values = ev$values[seq_len(r)], vectors = ev$vectors[, seq_len(r), drop = FALSE])
  }
}

# ---- ADMM solver for Step 1 (C) + Steps 2-3 (U) ----
admm_sgca <- function(Sigma, Sigma0, lambda, r,
                      rho = 1,
                      max_iter = 4000,
                      abs_tol = 1e-4, rel_tol = 1e-3,
                      sparsity_threshold = 1e-4,
                      penalize = c("all", "offdiag"),
                      weight = NULL,         # optional matrix of weights for weighted L1
                      warm_start = NULL,     # list(C=..., Z=..., U=...)
                      adapt_rho = TRUE, mu = 10, tau_incr = 2, tau_decr = 2,
                      verbose = FALSE) {
  
  penalize <- match.arg(penalize)
  p <- nrow(Sigma0)
  stopifnot(ncol(Sigma0) == p, nrow(Sigma) == p, ncol(Sigma) == p, r >= 1, r <= p)
  
  # EVD of Sigma0: Sigma0 = U0 diag(lam2) U0^T, where lam2 = (lambda_i)^2 in your notation.
  evd0 <- eigen((Sigma0 + t(Sigma0))/2, symmetric = TRUE)
  U0 <- evd0$vectors
  lam2 <- pmax(evd0$values, 0)                    # lam2 = λ_i^2 (nonnegative)
  Lam2 <- diag(lam2, nrow = p)
  
  # Precompute Sigma in the U0 basis: ~ constant term in RHS
  Sigma_tilde <- matmul(t(U0), matmul(Sigma, U0))
  const_rhs <- matmul(matmul(Lam2, Sigma_tilde), Lam2)
  
  # Denominator per (i,j): rho + lam2_i * lam2_j  (NOT lam2^2)
  denom <- rho + outer(lam2, lam2, `*`)
  
  # Initialize
  C <- if (!is.null(warm_start$C)) warm_start$C else matrix(0, p, p)
  Z <- if (!is.null(warm_start$Z)) warm_start$Z else C
  U <- if (!is.null(warm_start$U)) warm_start$U else matrix(0, p, p) # scaled dual
  
  # Penalization mask
  mask <- matrix(1, p, p)
  if (penalize == "offdiag") diag(mask) <- 0
  
  # Weights
  if (is.null(weight)) {
    W <- matrix(1, p, p)
  } else {
    stopifnot(all(dim(weight) == c(p, p)))
    W <- weight
  }
  W <- W * mask  # only applied where penalized
  
  # Residual histories
  primal_res <- numeric(max_iter)
  dual_res   <- numeric(max_iter)
  eps_pri_v  <- numeric(max_iter)
  eps_dual_v <- numeric(max_iter)
  
  for (iter in seq_len(max_iter)) {
    Z_prev <- Z
    
    # ----- C update (closed form in U0-basis) -----
    # RHS = rho * U0^T (Z - U) U0 + Lam2 * (U0^T Sigma U0) * Lam2
    rhs_tilde <- rho * matmul(t(U0), matmul(Z - U, U0)) + const_rhs
    C_tilde <- rhs_tilde / denom          # elementwise division
    C <- matmul(U0, matmul(C_tilde, t(U0)))
    
    # ----- Z update (soft-thresholding) -----
    Z_tilde <- C + U
    if (all(W == 1) && all(mask == 1)) {
      # fastest path: shrink all
      Z <- soft_threshold(Z_tilde, lambda / rho)
    } else {
      Z <- Z_tilde
      shr_idx <- (mask == 1)
      Z[shr_idx] <- soft_threshold(Z_tilde[shr_idx], (lambda / rho) * W[shr_idx])
    }
    
    # ----- dual update (scaled) -----
    U <- U + (C - Z)
    
    # ----- diagnostics & stopping -----
    r_norm <- norm(C - Z, type = "F")
    s_norm <- rho * norm(Z - Z_prev, type = "F")
    # epsilons per Boyd et al. (scaled form)
    eps_pri  <- sqrt(p * p) * abs_tol + rel_tol * max(norm(C, "F"), norm(Z, "F"))
    eps_dual <- sqrt(p * p) * abs_tol + rel_tol * rho * norm(U, "F")
    
    primal_res[iter] <- r_norm
    dual_res[iter]   <- s_norm
    eps_pri_v[iter]  <- eps_pri
    eps_dual_v[iter] <- eps_dual
    
    if (verbose && (iter %% 50 == 0L)) {
      cat(sprintf("iter %5d  r=%.3e  s=%.3e  eps_pri=%.3e  eps_dual=%.3e  rho=%.3g\n",
                  iter, r_norm, s_norm, eps_pri, eps_dual, rho))
    }
    
    if (r_norm <= eps_pri && s_norm <= eps_dual) {
      converged <- TRUE
      break
    }
    
    # ----- adaptive rho (optional) -----
    if (adapt_rho) {
      if (r_norm > mu * s_norm) {
        rho <- rho * tau_incr
        U   <- U / tau_incr
        denom <- rho + outer(lam2, lam2, `*`)
      } else if (s_norm > mu * r_norm) {
        rho <- rho / tau_decr
        U   <- U * tau_decr
        denom <- rho + outer(lam2, lam2, `*`)
      }
    }
    
    converged <- FALSE
  }
  
  # Step 2: eigenvectors of Sigma0^{1/2} C Sigma0^{1/2}
  Sigma0_sqrt <- matmul(U0, matmul(diag(sqrt(lam2), nrow = p), t(U0)))
  target <- matmul(Sigma0_sqrt, matmul(C, Sigma0_sqrt))
  eigU <- top_eigs_sym(target, r)
  U_svd <- eigU$vectors
  
  # Step 3: normalization  U <- U_svd (U_svd^T Sigma0 U_svd)^{-1/2}
  B <- matmul(t(U_svd), matmul(Sigma0, U_svd))
  B_inv_sqrt <- sym_inv_sqrt(B)
  U_canon <- matmul(U_svd, B_inv_sqrt)
  
  # Sparsities (fraction below threshold)
  C_sparsity <- mean(abs(C) < sparsity_threshold)
  U_sparsity <- mean(abs(U_canon) < sparsity_threshold)
  
  list(
    C = C,
    Z = Z,
    U_dual = U,
    U = U_canon,                # canonical directions (columns)
    C_sparsity = C_sparsity,
    U_sparsity = U_sparsity,
    iter = iter,
    converged = converged,
    primal = primal_res[seq_len(iter)],
    dual   = dual_res[seq_len(iter)],
    eps_pri = eps_pri_v[seq_len(iter)],
    eps_dual = eps_dual_v[seq_len(iter)],
    rho = rho, lambda = lambda
  )
}
