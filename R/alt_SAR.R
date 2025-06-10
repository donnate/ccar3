# LIBRARIES
library(glmnet)
library(MASS)
#library(pls)
library(CCA)


#' Function to perform Sparse CCA based on Wilms and Croux (2018)
#'  REFERENCE Wilms, I., & Croux, C. (2018). Sparse canonical correlation analysis using alternating regressions. Journal of Computational and Graphical Statistics, 27(1), 1-10.
#' @param X Matrix of predictors (n x p)
#' @param Y Matrix of responses (n x q)
#' @param lambdaAseq Vector of sparsity parameters for X (default is a sequence from 0 to 1 with step 0.1)
#' @param lambdaBseq Vector of sparsity parameters for Y (default is a sequence from 0 to 1 with step 0.1)
#' @param rank Number of canonical components to extract
#' @param selection.criterion Criterion for selecting the optimal tuning parameter (1 for minimizing difference between test and training sample correlation, 2 for maximizing test sample correlation)
#' @param n.cv Number of cross-validation folds (default is 5)
#' @param A.initial Initial value for the canonical vector A (default is NULL, which uses a canonical ridge solution)
#' @param B.initial Initial value for the canonical vector B (default is NULL, which uses a canonical ridge solution)
#' @param max.iter Maximum number of iterations for convergence (default is 20)
#' @param conv Convergence threshold (default is 1e-2)
#' @return A list with elements:
#' \describe{
#'   \item{U}{Canonical direction matrix for X (p x r)}
#'   \item{V}{Canonical direction matrix for Y (q x r)}
#'   \item{loss}{Mean squared error of prediction}
#'   \item{cor}{Canonical covariances}
#' }
#' @export
SparseCCA <- function(X, Y, lambdaAseq=seq(from=1, to=0.01, by=-0.01),
                    lambdaBseq=seq(from=1, to=0.01, by=-0.01),
                    rank, selection.criterion=1, n.cv=5, A.initial=NULL,
                    B.initial=NULL, max.iter=20, conv=10^-2, standardize=TRUE){
  ### Function to perform Sparse Canonical Correlation Analysis using alternating regressions

  ### INPUT
  # X                   : (nxp) data matrix
  # Y                   : (nxq) data matrix
  # lambdaAseq          : grid of sparsity parameters for lambdaA
  # lambdaBseq          : grid of sparsity parameters for lambdaB
  # rank                : number of canonical vector pairs to extract
  # selection.criterion : 1 for BIC,  2 for Cross-validation to select sparsity parameters
  # n.cv                : n.cv-fold cross-validation
  # A.initial           : starting value for the canonical vector A
  # B.initial           : starting value for the canonical vector B
  # max.iter            : maximum number of iterations
  # conv                : tolerance value for convergence

  ### OUTPUT
  # ALPHA               : (pxr) estimated canonical vectors correpsonding to the first data set
  # BETA                : (qxr) estimated canonical vectors correpsonding to the second data set
  # cancors             : r estimated canonical correlations
  # U_ALL               : (nxr) estimated canonical variates corresponding to the first data set
  # V_ALL               : (nxr) estimated canonical variates corresponding to the second data set
  # lamdbaA             : value of the sparsity parameter lambdaA
  # lamdbaB             : value of the sparsity parameter lambdaB
  # it                  : number of iterations


  ### STORE RESULTS
  ALPHA_ALL <- matrix(NA, ncol=rank, nrow=ncol(X))
  BETA_ALL <- matrix(NA, ncol=rank, nrow=ncol(Y))
  U_ALL <- matrix(NA, ncol=rank, nrow=nrow(X))
  V_ALL <- matrix(NA, ncol=rank, nrow=nrow(Y))
  cancors <- matrix(NA, ncol=rank, nrow=1)
  lambdaA_ALL <- matrix(NA, ncol=rank, nrow=max.iter)
  lambdaB_ALL <- matrix(NA, ncol=rank, nrow=max.iter)
  iterations <- matrix(NA, ncol=rank, nrow=1)
  obj.it <- matrix(NA, ncol=rank, nrow=max.iter+1)



  ### START CODE
  if (standardize){
    X = scale(X, center = TRUE, scale=TRUE)
    Y = scale(Y, center = TRUE, scale=TRUE)
  }

  # Starting Values: Canonical Ridge Solution
  if(is.null(A.initial)){
    cancor_regpar <- estim.regul_crossvalidation(X, Y, n.cv=n.cv,
                                                 lambda1grid = lambdaAseq,
                                                 lambda2grid = lambdaBseq)
    cancor_regul <- CCA::rcc(X, Y, cancor_regpar$lambda1.optim, cancor_regpar$lambda2.optim)
    A.initial_ridge <- matrix(cancor_regul$xcoef[, 1:rank], ncol=rank, nrow=ncol(X))
    B.initial_ridge <- matrix(cancor_regul$ycoef[, 1:rank], ncol=rank, nrow=ncol(Y))
    A.initial <- apply(A.initial_ridge, 2, NORMALIZATION_UNIT)
    B.initial <- apply(B.initial_ridge, 2, NORMALIZATION_UNIT)
  }

  # Sequentially extract the r canonical vector pairs
  for (i.r in 1:rank){

    if (i.r==1){ # for r=1: start from original data sets
      X_data <- X
      Y_data <- Y
    }

    # STARTING VALUES: canonical ridge
    A.STARTING <- matrix(A.initial[, i.r], ncol=1)
    B.STARTING <- matrix(B.initial[, i.r], ncol=1)

    # CONVERGENCE CRITERION
    obj.initial <- mean((X_data%*%matrix(A.STARTING, ncol=1)-Y_data%*%matrix(B.STARTING, ncol=1))^2)
    obj.it[1, i.r] <- obj.initial

    # INITIALIZE CONVERGENCE PARAMETERS
    it <- 1
    diff.obj <- conv*10

    # FROM i until convergence: canonical vectors
    while( (it<max.iter) & (diff.obj>conv) ){

      # Estimating A conditional on B
      FIT.A <- alternating.regression(Xreg=X_data, Yreg=Y_data%*%B.STARTING,
                                    lambdaseq=lambdaAseq,
                                    selection.criterion=selection.criterion,
                                    n.cv=n.cv)
      AHAT_FINAL <- FIT.A$COEF_FINAL
      lambdaA_ALL[it, i.r] <- FIT.A$LAMBDA_FINAL

      # Estimating B conditional on A
      FIT.B <- alternating.regression(Xreg=Y_data,
                                    Yreg=X_data%*%AHAT_FINAL,
                                    lambdaseq=lambdaBseq,
                                    selection.criterion=selection.criterion,
                                    n.cv=n.cv)
      BHAT_FINAL <- FIT.B$COEF_FINAL
      lambdaB_ALL[it, i.r] <- FIT.B$LAMBDA_FINAL

      # Check convergence
      obj.new <- mean((X_data%*%matrix(AHAT_FINAL, ncol=1)-Y_data%*%matrix(BHAT_FINAL, ncol=1))^2)
      obj.it[it+1, i.r] <- obj.new
      diff.obj <- abs(obj.new-obj.initial)/max(obj.initial, 1e-4)

      # Updated starting values
      B.STARTING <- BHAT_FINAL
      A.STARTING <- AHAT_FINAL
      obj.initial <- obj.new
      it <- it+1
    } # end while-loop

    # Number of ITERATIONS
    iterations[1, i.r] <- it

    # CANONICAL VARIATES after convergence
    Uhat <- X_data%*%AHAT_FINAL
    Vhat <- Y_data%*%BHAT_FINAL


    # Express canonical vectors in terms of ORIGINAL DATA MATRICES
    if (i.r==1){# FIRST DIMENSION

      # Final estimates of canonical vectors,  variates and canonical correlation
      ALPHA_ALL[, i.r] <- AHAT_FINAL
      BETA_ALL[, i.r] <- BHAT_FINAL
      U_ALL[, i.r] <- Uhat
      V_ALL[, i.r] <- Vhat
      cancors[1, i.r] <- cor(Uhat, Vhat)

      # Deflated data matrices
      X_data <-  round(X_data  - Uhat%*%solve(t(Uhat)%*%Uhat)%*%t(Uhat)%*%X_data, 10)
      Y_data <-  round(Y_data - Vhat%*%solve(t(Vhat)%*%Vhat)%*%t(Vhat)%*%Y_data, 10)

      # Sparsity parameters
      lambdaA_FINAL <- FIT.A$LAMBDA_FINAL
      lambdaB_FINAL <- FIT.B$LAMBDA_FINAL


    } else {# HIGHER ORDER DIMENSIONS

      # A expressed in terms of original data set X
      FIT.Aorig <- alternating.regression(Yreg=Uhat, Xreg=X, lambdaseq=lambdaAseq,
                                        selection.criterion=selection.criterion,
                                        n.cv=n.cv)
      ALPHAhat <- FIT.Aorig$COEF_FINAL
      lambdaA_FINAL <- FIT.Aorig$LAMBDA_FINAL

      # B expressed in terms of original data set Y
      FIT.Borig <- alternating.regression(Yreg=Vhat, Xreg=Y,
                                        lambdaseq=lambdaBseq,
                                        selection.criterion=selection.criterion,
                                        n.cv=n.cv)
      BETAhat <- FIT.Borig$COEF_FINAL
      lambdaB_FINAL <- FIT.Borig$LAMBDA_FINAL


      # Final estimates of canonical vectors,  variates and canonical correlation
      ALPHA_ALL[, i.r] <- ALPHAhat
      BETA_ALL[, i.r] <- BETAhat
      Uhat <- X%*%ALPHAhat
      Vhat <- Y%*%BETAhat
      U_ALL[, i.r] <- Uhat
      V_ALL[, i.r] <- Vhat
      cancors[1, i.r] <- cor(Uhat, Vhat)

      # Deflated data matrices: regress original data sets on all previously found canonical variates
      X_data <-  X  - U_ALL[, 1:i.r]%*%solve(t(U_ALL[, 1:i.r])%*%U_ALL[, 1:i.r])%*%t(U_ALL[, 1:i.r])%*%X
      Y_data <-  Y -  V_ALL[, 1:i.r]%*%solve(t(V_ALL[, 1:i.r])%*%V_ALL[, 1:i.r])%*%t(V_ALL[, 1:i.r])%*%Y
    }
  } # END FOR-LOOP

  ##OUTPUT
  out <- list(uhat=ALPHA_ALL, vhat=BETA_ALL, cancors=cancors, U_ALL=U_ALL, V_ALL=V_ALL, lambdaA=lambdaA_FINAL, lambdaB=lambdaB_FINAL, it=iterations)

}



alternating.regression <- function(Xreg, Yreg,
                                 lambdaseq=seq(from=1, to=0.01, by=-0.01),
                                 selection.criterion=1, n.cv=5){
  ### Function to perform sparse alternating regression

  ### INPUT
  #Xreg               : design matrix
  #Yreg               : response
  #lambdaseq          : sequence of sparsity parameters
  #selection.criterion: 1 for BIC,  2 for Cross-validation to select sparsity parameter
  #n.cv               : n.cv-fold cross-validation

  ### OUTPUT
  #COEF_FINAL         : estimated coefficients
  #LAMBDA_FINAL       : optimal sparsity parameter


  ##Standardize
  Xreg = scale(Xreg, center = TRUE, scale=TRUE)
  Xreg_st <- matrix(Xreg, ncol=ncol(Xreg))
  for (i.variable in 1:ncol(Xreg)){
    if (is.na(apply(Xreg_st, 2, sum)[i.variable])==T) {
      Xreg_st[, i.variable] <- 0}
  }

  ##LASSO FIT
  LASSOFIT <- glmnet(y=Yreg, x=Xreg_st, family="gaussian", lambda=lambdaseq, intercept=T)
  if (is.integer(which(LASSOFIT$df!=0)) && length(which(LASSOFIT$df!=0)) == 0L) {
    # Smaller lambda sequence necessary
    LASSOFIT <- glmnet(y=Yreg, x=Xreg_st, family="gaussian", intercept=T)
    COEFhat <- matrix(LASSOFIT$beta[, which(LASSOFIT$df!=0)[1:length(lambdaseq)]], nrow=ncol(Xreg_st)) # estimated coefficients
    LAMBDA <- LASSOFIT$lambda[which(LASSOFIT$df!=0)[1:length(lambdaseq)]] # lambda values
  } else {
    COEFhat <- matrix(LASSOFIT$beta[, which(LASSOFIT$df!=0)], nrow=ncol(Xreg_st)) # estimated coefficients
    LAMBDA <- LASSOFIT$lambda[which(LASSOFIT$df!=0)] # lambda values
  }

  ##Selection of sparsity parameter
  if (selection.criterion==1){ #BIC
    BICvalues <- apply(COEFhat, 2, BIC, Y.data=Yreg, X.data=Xreg_st) # BIC
    COEF_FINAL <- matrix(COEFhat[, which.min(BICvalues)], ncol=1) # Final coefficient estimates
    COEF_FINAL[which(apply(Xreg, 2, sd)!=0), ] <- COEF_FINAL[which(apply(Xreg, 2, sd)!=0), ]/apply(Xreg, 2, sd)[which(apply(Xreg, 2, sd)!=0)]
    COEF_FINAL <- apply(COEF_FINAL, 2, NORMALIZATION_UNIT)
    LAMBDA_FINAL <- LAMBDA[which.min(BICvalues)]

  } else {
    if (selection.criterion==2){ #Cross-validation
      n = nrow(Xreg_st)
      n.cv.sample <- trunc(n/n.cv)
      whole.sample <- seq(1, n)
      X.data <- Xreg_st
      Y.data <- Yreg
      cvscore <- matrix(NA, ncol=length(LAMBDA), nrow=n.cv)
      for (i in 1:n.cv){
        testing.sample <- whole.sample[((i-1)*n.cv.sample+1):(i*n.cv.sample)]
        training.sample <- whole.sample[-(((i-1)*n.cv.sample+1):(i*n.cv.sample))]
        Xcv = X.data[training.sample,  ]
        Ycv = Y.data[training.sample,  ]
        LASSOFIT_cv <- glmnet(y=Ycv, x=Xcv, family="gaussian", lambda=LAMBDA)
        COEF_cv <- LASSOFIT_cv$beta
        cvscore[i, ] <- apply(COEF_cv, 2, testsample.correlation, Xdata=X.data[testing.sample,  ], yscore = Y.data[testing.sample,  ] )
      }

      CVscore.mean <- apply(cvscore, 2, mean) # cv score
      COEF_FINAL <- matrix(COEFhat[, which.max(CVscore.mean)], ncol=1) # Final coefficient estimates
      COEF_FINAL[which(apply(Xreg, 2, sd)!=0), ] <- COEF_FINAL[which(apply(Xreg, 2, sd)!=0), ]/apply(Xreg, 2, sd)[which(apply(Xreg, 2, sd)!=0)]
      COEF_FINAL <- apply(COEF_FINAL, 2, NORMALIZATION_UNIT)
      LAMBDA_FINAL <- LAMBDA[which.max(CVscore.mean)]
    } else {
      stop("selection.criterion needs to be equal to 1 (BIC),  2 or 3 (Cross-validation)")
    }
  }

  ##OUTPUT
  out <- list(COEF_FINAL=COEF_FINAL, LAMBDA_FINAL=LAMBDA_FINAL)
}

NORMALIZATION_UNIT <- function(U){
  ### Normalize a vector U to have norm one
  if (is.null(dim(U))) U <- matrix(U, ncol = 1)


  length.U <- as.numeric(sqrt(t(U)%*%U))
  if(length.U==0){length.U <- 1}
  Uunit <- U/length.U
}

BIC <- function(U, Y.data, X.data){
  ### Calculate value of Bayesian information Criterion

  NEGLOGLIK <- sum(diag((1/nrow(Y.data))*(t(Y.data-X.data%*%U)%*%(Y.data-X.data%*%U))))
  BIC <- 2*NEGLOGLIK+ length(which(U!=0))*log(nrow(Y.data))
  return(BIC)
}

testsample.correlation <- function(U, Xdata, yscore){
  ### Calculate correlation in test sample
  if (is.null(dim(U))) U <- matrix(U, ncol = 1)

  xscore = Xdata%*%U
  if (all(U==0)) {
    return(0)
  } else {
    return(abs(cor(xscore, yscore)))
  }
}


############################################################
# CANONICAL RIDGE (used as starting value in SAR algorithm #
############################################################

estim.regul_crossvalidation <- function (X,  Y,  lambda1grid = NULL,
                                         lambda2grid = NULL, n.cv=5){

  
  
  if (is.null(lambda1grid)) {
    lambda1grid <- seq(0.001, 1, length = 5)
  } else {
    lambda1grid <- as.vector(lambda1grid)
  }
  
  if (is.null(lambda2grid)) {
    lambda2grid <- seq(0.001, 1, length = 5)
  } else {
    lambda2grid <- as.vector(lambda2grid)
  }
  
  lambda.grid <- expand.grid(lambda1=lambda1grid, lambda2=lambda2grid)
  
  cvscores <- mapply(function(l1, l2) {
    RCC_crossvalidation(X=X, Y=Y, lambda1=l1, lambda2=l2, n.cv=n.cv)$cv
  }, lambda.grid$lambda1, lambda.grid$lambda2)
  
  best_idx <- which.max(cvscores)
  lambda1.optim <- lambda.grid$lambda1[best_idx]
  lambda2.optim <- lambda.grid$lambda2[best_idx]
  cv.optim <- cvscores[best_idx]
  
  out <- list(lambda1.optim=lambda1.optim, lambda2.optim=lambda2.optim, cv.optim=cv.optim)
  return(out)
}


RCC_crossvalidation <- function (X,  Y,  lambda1,  lambda2, n.cv) {
# AUXILIARY FUNCTION CANONICAL RIDGE: n.cv-fold cross-validation
  ### INPUT
  # X                   : (nxp) data matrix
  # Y                   : (nxq) data matrix
  # lambda1             : tuning parameter for lambda1
  # lambda2             : tuning parameter for lambda2
  # n.cv                : n.cv-fold cross-validation

  ### OUTPUT
  # cv                  : value test sample canonical correlation

  n = nrow(X)
  ncvsample <- trunc(n/n.cv)
  ncvindex <- matrix(rep(1:n.cv), nrow=1)

  CVfit=apply(ncvindex, MARGIN=2, FUN=cvfunction, n=n,
              Xmatrix=X, Ymatrix=Y, lambda1=lambda1,
              lambda2=lambda2, ncvsample=ncvsample)
  cv <- sum(CVfit)/n.cv
  out=list(cv=cv)
}

cvfunction <- function(U, n, Xmatrix, Ymatrix, lambda1, lambda2, ncvsample){
 # AUXILIARY FUNCTION CANONICAL RIDGE: Canonical ridge fit
  cv <- 0
  whole.sample <- seq(1, n)
  testing.sample <- whole.sample[((U-1)*ncvsample+1):(U*ncvsample)]
  training.sample <- whole.sample[-(((U-1)*ncvsample+1):(U*ncvsample))]
  Xcv = Xmatrix[training.sample,  ]
  Ycv = Ymatrix[training.sample,  ]
  res = CCA::rcc(Xcv,  Ycv,  lambda1,  lambda2)

  xscore = Xmatrix[testing.sample,  ] %*% res$xcoef[,  1]
  yscore = Ymatrix[testing.sample,  ] %*% res$ycoef[,  1]
  cv <-  cv + abs(cor(xscore, yscore, use="pairwise"))
  return(cv)
}

