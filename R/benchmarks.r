
# source("R/alt_Parkhomenko.R")
# source("R/alt_SAR.R")
# source("R/alt_Waaijenborg.R")
# source("R/alt_Witten_CrossValidation.R")
# source("R/utils.R")

library(PMA)
#' @title Additional Benchmarks for Sparse CCA Methods
#' @param X_train Matrix of predictors (n x p)
#' @param Y_train  Matrix of responses (n x q)
#' @param S Optional covariance matrix (default is NULL, which computes it from X_train and Y_train)
#' @param rank Target rank for the CCA (default is 2)
#' @param method.type Type of method to use for Sparse CCA (default is "FIT_SAR_CV"). Choices include "FIT_SAR_BIC", "FIT_SAR_CV", "Witten_Perm", "Witten.CV", "Waaijenborg-Author", "Waaijenborg-CV", and "SCCA_Parkhomenko".
#' @param kfolds Number of cross-validation folds (default is 5)
#' @param lambdax Vector of sparsity parameters for X (default is a sequence from 0 to 1 with step 0.1)
#' @param lambday Vector of sparsity parameters for Y (default is a sequence from 0 to 1 with step 0.1)
#' 
#' 
#'
#' @return A matrix (p+q)x r containing the canonical directions for X and Y.
#' @export
sparse_CCA_benchmarks <- function(X_train, Y_train, S=NULL, 
                              rank=2, kfolds=5, method.type="FIT_SAR_CV",
                              lambdax = 10^seq(from=-3,to=2,length=10),
                              lambday = c(0, 1e-7, 1e-6, 1e-5),
                              standardize = TRUE){

  X_train = as.matrix(data.frame(X_train) %>% mutate_all(~replace_na(., mean(., na.rm = TRUE))))
  Y_train = as.matrix(data.frame(Y_train) %>% mutate_all(~replace_na(., mean(., na.rm = TRUE))))
  p1 <- dim(X_train)[2]
  p2 <- dim(Y_train)[2]
  p <- p1 + p2;
  n <- nrow(X_train)
  pp <- c(p1,p2);
  if(is.null(S)){
    S = cov(cbind(X_train, Y_train))
  }
  
  if (method.type=="FIT_SAR_BIC"){
    method<-SparseCCA(X=X_train,Y=Y_train,rank=rank,
                           lambdaAseq=lambdax,
                           lambdaBseq=lambday,
                           max.iter=100, conv=10^-2,
                           selection.criterion=1, n.cv=kfolds,
                           standardize=standardize)
    a_estimate = rbind(method$uhat, method$vhat)
    
  }
  if(method.type=="FIT_SAR_CV"){
    method<-SparseCCA(X=X_train,Y=Y_train,rank=rank,
                          lambdaAseq=lambdax,
                          lambdaBseq=lambday,
                          max.iter=100,conv=10^-2, selection.criterion=2, 
                          n.cv=kfolds,
                          standardize=standardize)
    a_estimate = rbind(method$uhat, method$vhat)
    
  }
  if (method.type=="Witten_Perm"){
    Witten_Perm <- PMA::CCA.permute(x=X_train,z=Y_train,
                               typex="standard",typez="standard", 
                               penaltyxs =lambdax[which(lambdax < 1)],
                               penaltyzs = lambday[which(lambday < 1)],
                               standardize = standardize,
                               nperms=50)
    
    method<-PMA::CCA(x=X_train, z=Y_train, typex="standard",typez="standard",K=rank,
                         penaltyx=Witten_Perm$bestpenaltyx,penaltyz=Witten_Perm$bestpenaltyz,trace=FALSE,
                         standardize = standardize)
    a_estimate = rbind(method$u, method$v)
  }
  if(method.type=="Witten.CV"){
    Witten_CV<-Witten.CV(X=X_train,Y=Y_train, n.cv=5,lambdax=lambdax[which(lambdax < 1)],
                         lambday=lambdax[which(lambdax < 1)])
    
    method <-PMA::CCA(x=X_train,z=Y_train,typex="standard",typez="standard",
                 K=rank,penaltyx=Witten_CV$lambdax.opt,
                 penaltyz=Witten_CV$lambday.opt,trace=FALSE,
                 standardize = standardize)
    a_estimate = rbind(method$u, method$v)
    
  }
  if(method.type=="Waaijenborg-Author"){
    method<-Waaijenborg(X=X_train,Y=Y_train,
                        lambdaxseq=lambdax,
                        lambdayseq=lambday,
                        rank=rank, selection.criterion=1,
                        standardize = standardize)
    a_estimate = rbind(method$vhat, method$uhat)
    
  }
  if(method.type=="Waaijenborg-CV"){
    method<-Waaijenborg(X=X_train,
                        Y=Y_train,lambdaxseq=lambdax,
                        lambdayseq=lambday,
                        rank=rank, selection.criterion=2,
                        standardize = standardize)
    a_estimate = rbind(method$vhat, method$uhat)
    
  }
  if(method.type=="SCCA_Parkhomenko"){
    method<- SCCA_Parkhomenko(x.data=X_train, y.data=Y_train, Krank=rank,
                              lambda.v.seq = lambdax[which(lambdax < 2)],
                              lambda.u.seq = lambday[which(lambday < 2)],
                              standardize = standardize)
    a_estimate = rbind(method$uhat, method$vhat)
    
  }
  a_estimate <- gca_to_cca(a_estimate, S, pp)
  return(a_estimate)
}



gca_to_cca <-
  function(a_estimate, S, pp){
    p1 = pp[1];
    p2 = pp[2];
    p = p1 + p2;
    nnz_indices = which(apply(a_estimate, 1, function(x) sqrt(sum(x^2))) >0)
    nnz_indices_x = nnz_indices[which(nnz_indices<(p1+1))]
    nnz_indices_y = nnz_indices[which(nnz_indices>(p1))]
    ### Make sure things are normalized
    if (length(which(nnz_indices<(p1+1)))>0){
      sigmaxhat = S[nnz_indices_x,nnz_indices_x];
      a_estimate[nnz_indices_x,] = a_estimate[nnz_indices_x,] %*% pracma::sqrtm(t(a_estimate[nnz_indices_x,]) %*% sigmaxhat %*% a_estimate[nnz_indices_x,])$Binv;
    }
    if (length(nnz_indices_y)>0){
      sigmayhat = S[nnz_indices_y,nnz_indices_y];
      a_estimate[nnz_indices_y,] = a_estimate[nnz_indices_y,] %*% pracma::sqrtm(t(a_estimate[nnz_indices_y,]) %*% sigmayhat %*% a_estimate[nnz_indices_y,])$Binv;
    }
    
    u_estimate = a_estimate[1:p1,]
    v_estimate = a_estimate[(p1+1):p,] 
    l = list("U" = u_estimate, "V" = v_estimate)
    return(l)
  }
