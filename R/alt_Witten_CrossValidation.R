library(PMA)
#' Sparse CCA by Witten and Tibshirani (2009)
#'
#' @param X Matrix of predictors (n x p)
#' @param Y Matrix of responses (n x q)
#' @param n.cv Number of cross-validation folds (default is 5)
#' @param lambdax Vector of sparsity parameters for X (default is a sequence from 0 to 1 with step 0.1)
#' @param lambday Vector of sparsity parameters for Y (default is a sequence from 0 to 1 with step 0.1)
#' 
#'
#' @return the appropriate levels of regularisation
#' @export

Witten.CV<-function(X,Y,n.cv=5,
                    rank,
                    lambdax=matrix(seq(from=0,to=1,by=0.1),nrow=1),
                    lambday=matrix(seq(from=0,to=1,by=0.1),nrow=1),
                    standardize = TRUE){ 


  n = nrow(X)
  n.cv.sample<-trunc(n/n.cv)
  whole.sample<-seq(1,n)
  lambdax=matrix(lambdax,nrow=1)
  lambday=matrix(lambday,nrow=1)

  if (standardize){
    X = scale(X, center = TRUE, scale=TRUE)
    Y = scale(Y, center = TRUE, scale=TRUE)
  }
  
  cvscore<-array(NA,c(length(lambday),length(lambdax),n.cv)) #lambdax in columns, lambday in rows

  
  for (i in 1:n.cv){
    testing.sample<-whole.sample[((i-1)*n.cv.sample+1):(i*n.cv.sample)]
    training.sample<-whole.sample[-(((i-1)*n.cv.sample+1):(i*n.cv.sample))]
    Xcv = X[training.sample, ]
    Ycv = Y[training.sample, ]
    Xtest=X[testing.sample,]
    Ytest=Y[testing.sample,]

    cvscore[,,i]<-apply(lambdax,2, Witten.cv.lambdax, Xtrain=Xcv,Ytrain=Ycv,Xtest=Xtest,Ytest=Ytest,lambday=lambday)
  }
  cvscore.vec<-c(apply(cvscore,c(1,2),mean))

  
  LAMBDAX<-c(matrix(rep(lambdax,length(lambday)),byrow=T,nrow=length(lambday),ncol=length(lambdax)))
  LAMBDAY<-c(matrix(rep(lambday,length(lambdax)),byrow=F,nrow=length(lambday),ncol=length(lambdax)))
  
  lambdax.opt<-LAMBDAX[which.max(cvscore.vec)]
  lambday.opt<-LAMBDAY[which.max(cvscore.vec)]

  ### OUTPUT
  out<-list(lambdax.opt=lambdax.opt,lambday.opt=lambday.opt)
}


Witten.cv.lambdax<-function(U,Xtrain,Ytrain,Xtest,Ytest,lambday){ #AUXILIARY FUNCTION
  testcorrelations<-apply(lambday,2,Witten.cv.lambday,lambdaxfixed=U,Xtrain=Xtrain,Ytrain=Ytrain,Xtest=Xtest,Ytest=Ytest)
  return(testcorrelations)
}

Witten.cv.lambday<-function(V,Xtrain,Ytrain,Xtest,Ytest,lambdaxfixed){ #AUXILIARY FUNCTION
  #print(lambdaxfixed)
  Fit.Witten<-PMA::CCA(x=Xtrain,z=Ytrain,typex="standard",typez="standard",K=1,penaltyx=lambdaxfixed,penaltyz=V,trace=F)
  return(abs(cor(Xtest%*%Fit.Witten$u, Ytest%*%Fit.Witten$v)))
}  

