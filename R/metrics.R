#########################################
# AUXILIARLY FUNCTIONS SIMULATION STUDY #
#########################################

#' @title Metrics for subspaces
#' @description Calculate principal angles between subspace spanned by the columns of a and the subspace spanned by the columns of b
#' @param a A matrix whose columns span a subspace.
#' @param b A matrix whose columns span a subspace.
#' @return  a vector of principal angles (in radians)
#' @examples
#' a <- matrix(rnorm(9), 3, 3)
#' b <- matrix(rnorm(9), 3, 3)
#' principal_angles(a, b)
#' @export
principal_angles <- function(a, b){
  ### Calculate principal angles between subspace spanned by the columns of a and the subspace spanned by the columns of b

  angles=matrix(0, ncol=ncol(a), nrow=1)
  qa= qr.Q(qr(a))
  qb= qr.Q(qr(b))
  C=svd(t(qa)%*%qb)$d
  rkA=qr(a)$rank;rkB=qr(b)$rank
  if (rkA<=rkB){
    B = qb - qa%*%(t(qa)%*%qb);
  } else {B = qa - qb%*%(t(qb)%*%qa)}
  S=svd(B)$d
  S=sort(S)

  for (i in 1:min(rkA, rkB)){
    if (C[i]^2 < 0.5) {angles[1, i]=acos(C[i])}
    else if (S[i]^2 <=0.5) {angles[1, i]=asin(S[i])}
  }
  angles=t(angles)

  ##OUTPPUT
  out=list(angles=angles)
  return(out)
}




#' @title True Positive Rate (TPR)
#' 
#' @description This is a function that compares the structure of two matrices A and B.  It outputs the number of entries that A and B have in common that are different from zero. A and B need to have the same number of rows and columns
#' @param A A matrix (the estimator).
#' @param B A matrix (assumed to be the ground truth).
#' @param tol tolerance for declaring the entries non zero.
#' @return True Positive Rate (nb of values that are non zero in both A and B / (nb of values that are non zero in A))
#' @examples
#' A <- matrix(c(1, 0, 0, 1, 1, 0), nrow = 2)
#' B <- matrix(c(1, 0, 1, 1, 0, 0), nrow = 2)
#' TPR(A, B)  
#' @export
TPR  <-  function(A, B, tol=1e-4){

  A[which(abs(A) <tol)] =0
  B[which(abs(B) <tol)] =0
  out  <-  sum((A!=0)*(B!=0))/max(1,sum(A!=0))
  return(out)
}



#' @title False Positive Rate (TPR)
#' @description This is a function that compares the structure of two matrices A and B. It outputs the number of entries where A  is not zero but Bis. A and B need to have the same number of rows and columns
#' @param A A matrix.
#' @param B A matrix (assumed to be the ground truth).
#' @param tol tolerance for declaring the entries non zero.
#' @return False Positive Rate (nb of values that are non zero in A and  zero in B / (nb of values that are non zero in A))
#' @examples
#' A <- matrix(c(1, 0, 0, 1, 1, 0), nrow = 2)
#' B <- matrix(c(1, 0, 1, 1, 0, 0), nrow = 2)
#' FPR(A, B)  
#' @export
FPR  <-  function(A, B, tol=1e-4){
  A[which(abs(A) <tol)] =0
  B[which(abs(B) <tol)] =0
  out  <-  sum((A!=0)*(B==0))/max(1,sum(A!=0))
  return(out)
}



#' True Negative Rate (TNR)
#' 
#' @description This is a function that compares the structure of two matrices A and B. It outputs the number of entries where A and B are both 0. A and B need to have the same number of rows and columns
#' @param A A matrix.
#' @param B A matrix (assumed to be the ground truth)..
#' @param tol tolerance for declaring the entries non zero.
#' @return True Negative Rate (nb of values that are zero in A and  zero in B / (nb of values that are  zero in A))
#' @export
TNR  <-  function(A, B, tol=1e-4){
  # This is a function that compares the structure of two matrices A and B
  # It outputs the number of entries that A and B have in common that are zero #
  # A and B need to have the same number of rows and columns
  A[which(abs(A) <tol)] =0
  B[which(abs(B) <tol)] =0
  out  <-  sum((A==0)*(B==0))/max(1,sum(A==0))
  return(out)
}



#' @title Subdistance between subspaces
#' @description Calculate subdistance between subspace spanned by the columns of a and the subspace spanned by the columns of b
#' @param A A matrix whose columns span a subspace.
#' @param B A matrix whose columns span a subspace.
#' @return subdistance between the two subspaces spanned by the matrices A and B, defined as min(O orthogonal) ||AO-B||_F
#' @export
subdistance <- function(A, B){
  svdresult = svd(t(A) %*% B);
  U = svdresult$u
  V = svdresult$v
  O = U %*% t(V);
  l = norm(A %*% O-B, type = c('F'));
  return(l)
}



#' @title SinTheta distance between subspaces
#' @description Calculate the distance spanned by the columns of A and the subspace spanned by the columns of B, defined as ||UU^T - VV^T||_F / sqrt(2)
#' @param U A matrix whose columns span a subspace.
#' @param V A matrix whose columns span a subspace.
#' @return sinTheta distance between the two subspaces spanned by the matrices A and B, defined as ||UU^T - VV^T||_F / sqrt(2)
#' @export
sinTheta<- function(U, V){
  l = 1/sqrt(2) * norm(U %*% t(U)-V %*% t(V), type = c('F'));
  return(l)
}