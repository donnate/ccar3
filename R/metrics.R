#########################################
# AUXILIARLY FUNCTIONS SIMULATION STUDY #
#########################################


#' Calculate principal angles between subspace spanned by the columns of a and the subspace spanned by the columns of b
#' @param a A matrix whose columns span a subspace.
#' @param b A matrix whose columns span a subspace.
#' 
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
}


#' True Positive Rate (TPR)
#'  This is a function that compares the structure of two matrices A and B.
#' It outputs the number of entries that A and B have in common that are different from zero. A and B need to have the same number of rows and columns
#' @param a A matrix.
#' @param b A matrix.
#' 
#' @export
TPR  <-  function(A, B, tol=1e-4){

  A[which(abs(A) <tol)] =0
  B[which(abs(B) <tol)] =0
  out  <-  sum((A!=0)*(B!=0))/max(1,sum(A!=0))
  return(out)
}


#' False Positive Rate (TPR)
#' This is a function that compares the structure of two matrices A and B.
#' It outputs the number of entries where A  is not zero but Bis.
#' A and B need to have the same number of rows and columns
#' @param a A matrix.
#' @param b A matrix.
#' 
#' @export
FPR  <-  function(A, B, tol=1e-4){
  A[which(abs(A) <tol)] =0
  B[which(abs(B) <tol)] =0
  out  <-  sum((A!=0)*(B==0))/max(1,sum(A!=0))
  return(out)
}


#' True Negative Rate (TPR)
#' This is a function that compares the structure of two matrices A and B.
#' It outputs the number of entries where A and B are both 0.
#' A and B need to have the same number of rows and columns
#' @param a A matrix.
#' @param b A matrix.
#' 
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


#' Calculate subdistance between subspace spanned by the columns of a and the subspace spanned by the columns of b
#' @param a A matrix whose columns span a subspace.
#' @param b A matrix whose columns span a subspace.
#' 
#' @export
subdistance <- function(A, B){
  svdresult = svd(t(A) %*% B);
  U = svdresult$u
  V = svdresult$v
  O = U %*% t(V);
  l = norm(A %*% O-B, type = c('F'));
  return(l)
}


#' sinTheta spanned by the columns of a and the subspace spanned by the columns of b
#' @param a A matrix whose columns span a subspace.
#' @param b A matrix whose columns span a subspace.
#' 
#' @export
sinTheta<- function(U, V){
  l = 1/sqrt(2) * norm(U %*% t(U)-V %*% t(V), type = c('F'));
  return(l)
}