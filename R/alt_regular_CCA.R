

#' Function to perform regular (low dimensional) canonical correlation analysis (CCA
#' @param X Matrix of predictors (n x p)
#' @param Y Matrix of responses (n x q)
#' @param rank Number of canonical components to extract
#' @return A list with elements:
#' \describe{
#'   \item{U}{Canonical direction matrix for X (p x r)}
#'   \item{V}{Canonical direction matrix for Y (q x r)}
#'   \item{cor}{Canonical covariances}
#' }
#' @export
regular_cca <- function(X, Y, rank){
  Sigma_x = t(X)%*% X
  Sigma_y = t(Y)%*% Y
  Sigma_xy = t(X)%*% Y
  svd_x = svd(Sigma_x)
  inv_sqrt_Sigma_x = svd_x$u %*% diag(sapply(svd_x$d, function(x){ifelse(x<1e-5, 0, 1/sqrt(x))})) %*% t(svd_x$v)
  svd_y = svd(Sigma_y)
  inv_sqrt_Sigma_y = svd_y$u %*% diag(sapply(svd_y$d, function(x){ifelse(x<1e-5, 0, 1/sqrt(x))})) %*% t(svd_y$v)
  svd_for_cca = svd(inv_sqrt_Sigma_x %*% Sigma_xy %*% inv_sqrt_Sigma_y)
  xcoef = inv_sqrt_Sigma_x %*% svd_for_cca$u[, 1:rank]
  ycoef = inv_sqrt_Sigma_y %*% svd_for_cca$v[, 1:rank]
  cancor <- svd_for_cca$d
  return(list(xcoef=xcoef, ycoef=ycoef, cancor = cancor))
}
