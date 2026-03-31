# environment to store whether we've already warned
.ccar3_internal_env <- new.env(parent = emptyenv())
.ccar3_internal_env$smut_warned <- FALSE


matmul <- function(A, B) {
  if (requireNamespace("SMUT", quietly = TRUE)) {
    # Use the fast C++ multiplication from SMUT
    SMUT::eigenMapMatMult(A, B)
  } else {
      # If not installed, warn once per session
      if (!.ccar3_internal_env$smut_warned) {
        packageStartupMessage("Using base R %*% for matrix multiplication. 
          Install 'SMUT' for faster performance.")
        .ccar3_internal_env$smut_warned <- TRUE
      }
    # Fallback to base R multiplication
    A %*% B
  }
}

compute_sqrt_inv <- function(S, threshold = 1e-4) {
  svd_S <- svd(S)
  svd_S$u %*% diag(sapply(svd_S$d, function(x) ifelse(x > threshold, 1 / sqrt(x), 0))) %*% t(svd_S$v)
}

compute_sqrt <- function(S, threshold = 1e-4) {
  svd_S <- svd(S)
  svd_S$u %*% diag(sapply(svd_S$d, function(x) ifelse(x > threshold, sqrt(x), 0))) %*% t(svd_S$u)
}

soft_thresh = function(A, lambda){
  #A * pmax(1 - lambda / abs(A), 0)
  sign(A) * pmax(abs(A) - lambda, 0)
}


fnorm = function(A){
  sqrt(sum(A^2))
}



soft_thresh_group = function(A, lambda){
  norm_A <- sqrt(sum(A^2))
  if (norm_A == 0) return(A)
  A * pmax(0, 1 - lambda / norm_A)
}



soft_thresh2 <- function(A, lambda){
  norm_A <- sqrt(sum(A^2))
  if(norm_A == 0){
    return(A)
  }
  result = A * pmax(1 - lambda/(norm_A), 0)
  return(result)
}


rmat = function(n, p){
  matrix(rnorm(n * p), n, p)
}


