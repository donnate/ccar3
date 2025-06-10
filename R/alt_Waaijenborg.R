# LIBRARIES
library(glmnet)
library(MASS)
library(pls)


#' Function to perform Sparse CCA based on Waaijenborg et al. (2008)
#'  REFERENCE Waaijenborg et al. (2008), "Quantifying the Association between Gene Expressions and DNA-Markers by Penalized Canonical Correlation Analysis" in Statistical Applications in Genetics and Molecular Biology, Volume 7, Issue 1, Article 3
#'
#' @param X Matrix of predictors (n x p)
#' @param Y Matrix of responses (n x q)
#' @param n.cv Number of cross-validation folds (default is 5)
#' @param lambdaxseq Vector of sparsity parameters for X (default is a sequence from 0 to 1 with step 0.1)
#' @param lambdayseq Vector of sparsity parameters for Y (default is a sequence from 0 to 1 with step 0.1)
#' @param rank Number of canonical components to extract
#' @param selection.criterion Criterion for selecting the optimal tuning parameter (1 for minimizing difference between test and training sample correlation, 2 for maximizing test sample correlation)
#' @return A list with elements:
#' \describe{
#'   \item{U}{Canonical direction matrix for X (p x r)}
#'   \item{V}{Canonical direction matrix for Y (q x r)}
#'   \item{loss}{Mean squared error of prediction}
#'   \item{cor}{Canonical covariances}
#' }
#' @export
Waaijenborg<-function(X, Y, lambdaxseq=matrix(seq(from=1,to=0,by=-0.01),nrow=1),
                      lambdayseq=matrix(seq(from=1,to=0,by=-0.01),nrow=1), rank,
                      selection.criterion=1,n.cv=5,max.iter=20,conv=10^-3, standardize=TRUE){
  ### Function to perform Sparse CCA based on Waaijenborg et al. (2008)
  # REFERENCE Waaijenborg et al. (2008), "Quantifying the Association between Gene Expressions and DNA-Markers by Penalized Canonical Correlation Analysis" in Statistical Applications in Genetics and Molecular Biology, Volume 7, Issue 1, Article 3
  if (standardize){
    X <- scale(X, center = TRUE, scale=TRUE)
    Y <- scale(Y, center = TRUE, scale=TRUE)
  }
  ### STORE RESULTS
  lambdaxseq=matrix(lambdaxseq,nrow=1)
  lambdayseq=matrix(lambdayseq,nrow=1)
  
  u_ALL<-matrix(NA,ncol=rank,nrow=ncol(Y))
  v_ALL<-matrix(NA,ncol=rank,nrow=ncol(X))
  ksi_ALL<-matrix(NA,ncol=rank,nrow=nrow(X))
  omega_ALL<-matrix(NA,ncol=rank,nrow=nrow(Y))
  cancors<-matrix(NA,ncol=rank,nrow=1)
  # Store sparsity parameters
  lambdax_ALL<-matrix(NA,ncol=rank,nrow=max.iter)
  lambday_ALL<-matrix(NA,ncol=rank,nrow=max.iter)
  iterations<-matrix(NA,ncol=rank,nrow=1)
  
  
  ### START CODE
  
  # Starting Values
  p<-ncol(Y)
  u.starting<-rep(1/p,p)
  q<-ncol(X)
  v.starting<-rep(1/q,q)
  

  # Sequentially extract the r canonical vector pairs
  for (i.r in 1:rank){
    
    if (i.r==1){ # for r=1: start from original data sets
      X_data<-X 
      Y_data<-Y
    } 
    
    # STARTING VALUES
    u.initial<-u.starting
    v.initial<-v.starting
    
    # STANDARDIZE DATA (Data is already standardized as an option)
    X_data_st<- X
    for (i.variable in 1:ncol(X_data)){
      if (is.na(apply(X_data_st,2,sum)[i.variable])==T) {
        X_data_st[,i.variable]<-0}    
    }
    
    Y_data_st<-Y
    for (i.variable in 1:ncol(Y_data)){
      if (is.na(apply(Y_data_st,2,sum)[i.variable])==T) {
        Y_data_st[,i.variable]<-0}    
    }
    
    
    # CANONICAL VARIATES
    ksi<-X_data_st%*%v.initial
    omega<-Y_data_st%*%u.initial
    
    # INITIALIZE CONVERGENCE PARAMETERS
    it<-1
    diff.u<-conv*10
    diff.v<-conv*10
    
    
    while( (it<max.iter) & ((diff.u>conv) & (diff.v>conv))) {
      
      # Estimate canonical vector u
      U.hat<-apply(lambdayseq,2,UST,a=t(ksi)%*%Y_data_st)
      U.hat.reduced<- U.hat[,which(apply(U.hat,2,sum)!=0), drop = FALSE]
      lambday.reduced<-matrix(lambdayseq[which(apply(U.hat,2,sum)!=0)],nrow=1)
      
      # Select tuning parameter
      if (selection.criterion==1){ #Minimize difference between test and training sample correlation 
            n = nrow(Y_data_st)
            n.cv.sample<-trunc(n/n.cv)
            whole.sample<-seq(1,n)
            Y.data<-Y_data_st
            ksi.data<-ksi
            cvscore<-matrix(NA,ncol=length(lambday.reduced),nrow=n.cv)
            for (i in 1:n.cv){
              testing.sample<-whole.sample[((i-1)*n.cv.sample+1):(i*n.cv.sample)]
              training.sample<-whole.sample[-(((i-1)*n.cv.sample+1):(i*n.cv.sample))]
              Ycv = Y.data[training.sample, ]
              ksicv = ksi[training.sample, ]
              U.hat<-apply(lambday.reduced,2,UST,a=t(ksicv)%*%Ycv)
              
              cvscore[i,]<-apply(U.hat,2,delta.correlation,Xtest=Y.data[testing.sample, ],Ytest=ksi[testing.sample, ],Xtrain=Y.data[training.sample, ],Ytrain=ksi[training.sample, ])  
            
            } 
            CVscore.mean<-apply(cvscore,2,mean)
            if (is.null(dim(U.hat.reduced))){
              U.hat.reduced = matrix(U.hat.reduced, ncol=1)
            }
            U.hat.final<- U.hat.reduced[,min(which.max(CVscore.mean)), drop = FALSE] # OPTIMAL SOLUTION
            U.hat.final<-apply(U.hat.final,2, NORMALIZATION_UNIT)
            lambday_ALL[it,i.r]<-lambday.reduced[min(which.max(CVscore.mean))]
        
      } else {
            if (selection.criterion==2){ #Maximize test sample correlation
                  n = nrow(Y_data_st)
                  n.cv.sample<-trunc(n/n.cv)
                  whole.sample<-seq(1,n)
                  Y.data<-Y_data_st
                  ksi.data<-ksi
                  cvscore<-matrix(NA,ncol=length(lambday.reduced),nrow=n.cv)
                  for (i in 1:n.cv){
                    testing.sample<-whole.sample[((i-1)*n.cv.sample+1):(i*n.cv.sample)]
                    training.sample<-whole.sample[-(((i-1)*n.cv.sample+1):(i*n.cv.sample))]
                    Ycv = Y.data[training.sample, ]
                    ksicv = ksi[training.sample, ]
                    U.hat<-apply(lambday.reduced,2,UST,a=t(ksicv)%*%Ycv)
                    cvscore[i,]<-apply(U.hat,2,testsample.correlation,Xdata=Y.data[testing.sample, ],yscore=ksi[testing.sample, ] )  
                  } 
                  
                  CVscore.mean<-apply(cvscore,2,mean)
                  if (is.null(dim(U.hat.reduced))){
                    U.hat.reduced = matrix(U.hat.reduced, ncol=1)
                  }
                  U.hat.final<-U.hat.reduced[, min(which.max(CVscore.mean)), drop = FALSE] # OPTIMAL SOLUTION
                  U.hat.final<-apply(U.hat.final,2,NORMALIZATION_UNIT)
                  lambday_ALL[it,i.r]<-lambday.reduced[min(which.max(CVscore.mean))]
            } else {
              stop("selection.criterion needs to be equal to 1 (Difference test and training correlation) or 2 (Test correlation)")
            }
            
      }
            
      
      # Estimate canonical vector v
      V.hat<-apply(lambdaxseq,2,UST,a=t(omega)%*%X_data_st)
      V.hat.reduced<-V.hat[,which(apply(V.hat,2,sum)!=0), drop = FALSE]
      lambdax.reduced<-matrix(lambdaxseq[which(apply(V.hat,2,sum)!=0)],nrow=1)
      
      # Select tuning parameter
      if (selection.criterion==1){
              n = nrow(X_data_st)
              n.cv.sample<-trunc(n/n.cv)
              whole.sample<-seq(1,n)
              X.data<-X_data_st
              omega.data<-omega
              cvscore<-matrix(NA,ncol=length(lambdax.reduced),nrow=n.cv)
              for (i in 1:n.cv){
                testing.sample<-whole.sample[((i-1)*n.cv.sample+1):(i*n.cv.sample)]
                training.sample<-whole.sample[-(((i-1)*n.cv.sample+1):(i*n.cv.sample))]
                Xcv = X.data[training.sample, ]
                omegacv = omega[training.sample, ]
                V.hat<-apply(lambdax.reduced,2,UST,a=t(omegacv)%*%Xcv)
                
                cvscore[i,]<-apply(V.hat,2,delta.correlation,Xtest=X.data[testing.sample, ],Ytest=omega[testing.sample, ],Xtrain=X.data[training.sample, ],Ytrain=omega[training.sample, ] )  
                
              } 
              CVscore.mean<-apply(cvscore,2,mean)
              if (is.null(dim(V.hat.reduced))){
                V.hat.reduced = matrix(V.hat.reduced, ncol=1)
              }
              V.hat.final<-V.hat.reduced[,min(which.max(CVscore.mean)), drop = FALSE]# OPTIMAL SOLUTION
              V.hat.final<-apply(V.hat.final,2,NORMALIZATION_UNIT)
              lambdax_ALL[it,i.r]<-lambdax.reduced[min(which.max(CVscore.mean))]
        } else {
              if (selection.criterion==2){
                n = nrow(X_data_st)
                n.cv.sample<-trunc(n/n.cv)
                whole.sample<-seq(1,n)
                X.data<-X_data_st
                omega.data<-omega
                cvscore<-matrix(NA,ncol=length(lambdax.reduced),nrow=n.cv)
                for (i in 1:n.cv){
                  testing.sample<-whole.sample[((i-1)*n.cv.sample+1):(i*n.cv.sample)]
                  training.sample<-whole.sample[-(((i-1)*n.cv.sample+1):(i*n.cv.sample))]
                  Xcv = X.data[training.sample, ]
                  omegacv = omega[training.sample, ]
                  V.hat<-apply(lambdax.reduced,2,UST,a=t(omegacv)%*%Xcv)
                  
                  cvscore[i,]<-apply(V.hat,2,testsample.correlation,Xdata=X.data[testing.sample, ],yscore=omega[testing.sample, ] )  
                } 
                
                CVscore.mean<-apply(cvscore,2,mean)
                if (is.null(dim(V.hat.reduced))){
                  V.hat.reduced = matrix(V.hat.reduced, ncol=1)
                }
                V.hat.final<-matrix(V.hat.reduced[,min(which.max(CVscore.mean)), drop = FALSE],ncol=1) # OPTIMAL SOLUTION
                V.hat.final<-apply(V.hat.final,2,NORMALIZATION_UNIT)
                lambdax_ALL[it,i.r]<-lambdax.reduced[min(which.max(CVscore.mean))]
              } else {
                stop("selection.criterion needs to be equal to 1 (Difference test and training correlation) or 2 (Test correlation)")
              }
        }
              
        # Convergence measures
        diff.u<-max(abs(u.initial-U.hat.final))
        diff.v<-max(abs(v.initial-V.hat.final))
             
        # Updated starting values 
        ksi<-X_data_st%*%V.hat.final
        omega<-Y_data_st%*%U.hat.final    
        u.initial<-U.hat.final
        v.initial<-V.hat.final
      
        it<-it+1
      }
      
      # Number of ITERATIONS
      iterations[1,i.r]<-it
    
    # Final estimates of canonical vectors, variates and canonical correlation
    u_ALL[,i.r]<-U.hat.final
    v_ALL[,i.r]<-V.hat.final
    ksi_ALL[,i.r]<-ksi
    omega_ALL[,i.r]<-omega
    cancors[,i.r]<-abs(cor(ksi,omega))

    # Residual matrices
    gammahat<-matrix(MASS::ginv(t(ksi)%*%ksi)%*%t(ksi)%*%X_data,ncol=1)
    thetahat<-matrix(MASS::ginv(t(omega)%*%omega)%*%t(omega)%*%Y_data,ncol=1)

    X_data<-X_data-ksi%*%t(gammahat)
    Y_data<-Y_data-omega%*%t(thetahat)

  }  # END FOR-LOOP
    
    
  ##OUTPUT 
  out<-list(uhat=u_ALL,vhat=v_ALL,ksihat=ksi_ALL,omegahat=omega_ALL,cancors=cancors,lambdax_FINAL=lambdax_ALL,lambday_FINAL=lambday_ALL,it=iterations)
  
}

UST<-function(a,U){ # Univariate Soft Thresholding

  matrix((abs(a)-U/2 + abs(abs(a)-U/2))/2*sign(a),ncol=1)
}

delta.correlation<-function(U,Xtest,Ytest,Xtrain,Ytrain){ # Cross-validation criterium
  if (all(U==0)) {
    return(0)
  } else {
    traincor=abs(cor(Xtrain%*%U,Ytrain))
    testcor=abs(cor(Xtest%*%U,Ytest))
    return(abs(traincor-testcor))
  }
  
}
