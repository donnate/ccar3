#' Function to perform Sparse CCA based on Waaijenborg et al. (2008)
#'  REFERENCE Parkhomenko et al. (2009), "Sparse Canonical Correlation Anlaysis with Application to Genomic Data Integration" in  Statistical Applications in Genetics and Molecular Biology, Volume 8, Issue 1, Article 1
#'
#' @param x.data Matrix of predictors (n x p)
#' @param y.data Matrix of responses (n x q)
#' @param n.cv Number of cross-validation folds (default is 5)
#' @param lambda.u.seq Vector of sparsity parameters for X (default is a sequence from 0 to 1 with step 0.1)
#' @param lambda.v.seq Vector of sparsity parameters for Y (default is a sequence from 0 to 1 with step 0.1)
#' @param Krank Number of canonical components to extract
#' @return A list with elements:
#' \describe{
#'   \item{U}{Canonical direction matrix for X (p x r)}
#'   \item{V}{Canonical direction matrix for Y (q x r)}
#'   \item{cor}{Canonical correlations}
#' }
#' @export

SCCA_Parkhomenko<-function(x.data,y.data,n.cv=5,lambda.v.seq=seq(0, 0.2, by=0.02), lambda.u.seq=seq(0, 0.2, by=0.02),Krank=1, standardize=TRUE){
  ### Function to perform Sparse CCA based on Parkhomenko et al. (2009)
  # REFERENCE Parkhomenko et al. (2009), "Sparse Canonical Correlation Anlaysis with Application to Genomic Data Integration" in  Statistical Applications in Genetics and Molecular Biology, Volume 8, Issue 1, Article 1
  
  ### INPUT
  # x.data:         1st datamatrix n times p
  # y.data:         2nd datamatrix n times q
  # n.cv:           n.cv-fold cross-validation
  # lambda.v.seq:   Possible values of sparseness parameters for data Y. Lower bouns should be 0, upper bound can be increased to 2.
  # lambda.u.seq:   Possible values of sparseness parameters for data X. Lower bouns should be 0, upper bound can be increased to 2.
  # Krank:          rank approximation singular value decompostion
  
  
  ### OUTPUT
  # a:              canonical vectors accosiated to matrix X
  # b:              Canonical vectors associated to matrix Y
  # cancor:         Canonical Correlations
  #lambda.uopt:     value of sparsity parameter lambdau
  #lambda.vopt      value of sparsity parameter lamdbav

  # START CODE
  
  # SCCA to obtain first pair of canonical vectors

  if (standardize){
    x.data <- scale(x.data, center = TRUE, scale=TRUE)
    y.data <-  scale(y.data, center = TRUE, scale=TRUE)
  }
  
  SCCA_FIT_K1 <- SCCA_Parkhomenko_K1(x.data=x.data,y.data=y.data,n.cv=n.cv,
                                     lambda.v.seq=lambda.v.seq, lambda.u.seq=lambda.u.seq)
  
  
  if (Krank==1){ #1 canonical vector pair
    ##OUTPUT
    out<-list(a=SCCA_FIT_K1$a,b=SCCA_FIT_K1$b,cancor=SCCA_FIT_K1$cancor,lambda.uopt=SCCA_FIT_K1$best.lambda.u,lambda.vopt=SCCA_FIT_K1$best.lambda.v)

  }   else{ #more than 1 canonical vector pair
    
    ### STORE RESULTS
    a_vectors<-matrix(data=0,ncol=Krank,nrow=ncol(x.data))
    b_vectors<-matrix(data=0,ncol=Krank,nrow=ncol(y.data))
    
    
    # First pair of canonical vectors
    a_vectors[,1]<- SCCA_FIT_K1$a
    b_vectors[,1]<- SCCA_FIT_K1$b 
    
    # Obtain higher-order pairs of canonical vectors
    for (irank in 2:Krank){ 
      if (irank==2){
        SCCA_FIT<-SCCA_Parkhomenko_Khigher(x.data=x.data,y.data=y.data,k=SCCA_FIT_K1$K_residual,best.lambda.u=max(SCCA_FIT_K1$best.lambda.u),best.lambda.v=max(SCCA_FIT_K1$best.lambda.v))
        a_vectors[,irank]<- SCCA_FIT$a 
        b_vectors[,irank]<- SCCA_FIT$b
        D<-diag(c(SCCA_FIT_K1$cancor,SCCA_FIT$cancor))
        
        # Work on residual space to extract higher order canoncical vectors
        K_hat<-a_vectors[,1:irank]%*%D%*%t(b_vectors[,1:irank])
        K_residual<-SCCA_FIT_K1$k-K_hat
      }   else {
        SCCA_FIT<-SCCA_Parkhomenko_Khigher(x.data=x.data,y.data=y.data,k=K_residual,best.lambda.u=max(SCCA_FIT_K1$best.lambda.u),best.lambda.v=max(SCCA_FIT_K1$best.lambda.v))
        a_vectors[,irank]<- SCCA_FIT$a 
        b_vectors[,irank]<- SCCA_FIT$b
        D<-diag(c(diag(D),SCCA_FIT$cancor))
        K_hat<-a_vectors[,1:irank]%*%D%*%t(b_vectors[,1:irank])
        K_residual<-SCCA_FIT_K1$k-K_hat
      } 
    }
    
    ##OUTPUT
    out <- list(uhat=a_vectors,vhat=b_vectors,cancor=c(diag(D)),lambda.uopt=SCCA_FIT_K1$best.lambda.u,lambda.vopt=SCCA_FIT_K1$best.lambda.v)
  } 
}


SCCA_Parkhomenko_K1<-function(x.data,y.data,n.cv=5,lambda.v.seq=seq(0, 0.2, by=0.02), 
                              lambda.u.seq=seq(0, 0.2, by=0.02)){ # Sparse CCA for first dimension
  # CODE taken from http://www.uhnres.utoronto.ca/labs/tritchler/
  
  ### INPUT
  # x.data:         1st datamatrix n times p
  # y.data:         2nd datamatrix n times q
  # n.cv:           number of cross-validation steps to select sparseness parameters
  # lambda.v.seq:   Possible values of sparseness parameters for data Y. Lower bouns should be 0, upper bound can be increased to 2.
  # lambda.u.seq:   Possible values of sparseness parameters for data X. Lower bouns should be 0, upper bound can be increased to 2.
  
  ### OUTPUT
  # v:              right singular vector associated to matrix Y
  # u:              left singular vector associated to matrix X
  # cancor:         first canonical correlation
  # a:              canonical vectors accosiated to matrix X
  # b:              Canonical vectors associated to matrix Y
  # K_hat:          rank 1 approximation of matrix K
  # K_residual:     K-K_hat
  # k:              matrix K (sample covariance matrix)
  # best.lambda.u:  optimal sparseness parameter for matrix X
  # best.lambda.v:  optimal sparseness parameter for matrix Y
  
  
  # ANALYSIS
  # Dimensions of input-dataset
  n.sample<-nrow(x.data) # Sample size
  p<-ncol(x.data) # Number of variables in X
  q<-ncol(y.data) # Number of variables in Y
  
  # Settings sparseness parameters
  n.lambdas.u <-  length(lambda.u.seq)
  n.lambdas.v <-  length(lambda.v.seq)
  lambda.v.matrix <- matrix(rep(lambda.v.seq, n.lambdas.u), nrow=n.lambdas.u, byrow=TRUE)
  lambda.u.matrix <- matrix(rep(lambda.u.seq, n.lambdas.v), nrow=n.lambdas.u, byrow=FALSE)
  ones.p <- rep(1, p)/p
  ones.q <- rep(1, q)/q
  
  
  # SCCA Analysis
  n.cv.sample <- trunc(n.sample/n.cv)
  whole.sample <- seq(1, n.sample)
  
  predict.corr.scca <- matrix(0, nrow=n.lambdas.u, ncol=n.lambdas.v)  # This matrix will contain average test sample correlation for each combination of sparseness parameters
  
  #_______Cross-validation to select optimal combination of sparseness parameters____________
  for (i.cv in 1:n.cv)
  {
    testing.sample <- whole.sample[((i.cv-1)*n.cv.sample+1):(i.cv*n.cv.sample)]
    training.sample <- whole.sample[!whole.sample%in%testing.sample]
    
    k <- sample.sigma12.function(x.data[training.sample, ], y.data[training.sample, ])
    
    # Get starting values for singular vectors
    # as column and row means from matrix K
    u.initial <- k %*% ones.q
    u.initial <- u.initial /sqrt(as.numeric(t(u.initial)%*%u.initial))
    v.initial <- t(k) %*% ones.p
    v.initial <- v.initial /sqrt(as.numeric(t(v.initial)%*%v.initial))
    
    # _______________Data for Predicted correlation (testing sample)_________________
    
    x.predict <- x.data[testing.sample, ]
    y.predict <- y.data[testing.sample, ]
    
    # Standardize data
    x.predict <- x.predict - mean(x.predict)
    y.predict <- y.predict - mean(y.predict)
    
    sigma11.predict <- var(x.predict)
    sigma22.predict <- var(y.predict)
    
    x.predict <- x.predict %*% diag( 1/sqrt(diag(sigma11.predict)) )
    y.predict <- y.predict %*% diag( 1/sqrt(diag(sigma22.predict)) )
    
    
    # ____________Loops for sparseness parameter combinations__________
    for(j.lambda.v in 1:n.lambdas.v)
    {
      
      flag.na <- 0
      
      for(j.lambda.u in 1:n.lambdas.u)
      {
        lambda.v <- lambda.v.seq[j.lambda.v]  # sparseness parameter for Y
        lambda.u <- lambda.u.seq[j.lambda.u]  # sparseness parameter for X
        
        if(flag.na==0)
        {
          uv <- scca.function(k, u.initial, v.initial, lambda.u, lambda.v)
          
          vj <- uv$v.new
          uj <- uv$u.new
          
          # Calculate predicted correlation for SCCA
          predict.corr.scca[j.lambda.u, j.lambda.v] <- predict.corr.scca[j.lambda.u, j.lambda.v] + abs(cor(x.predict%*%uj, y.predict%*%vj))
          if(is.na(predict.corr.scca[j.lambda.u, j.lambda.v])) flag.na <- 1
        }  # close if
        
        if(flag.na==1)
        {
          predict.corr.scca[j.lambda.u:n.lambdas.u, j.lambda.v] <- predict.corr.scca[j.lambda.u:n.lambdas.u, j.lambda.v] + NA
          break
        }
        
      }  # close loop on lambda.u
    }  # close loop on lambda.v
    
  }	# close cross-validation loop
  
  # ______________Identify optimal sparseness parameter combination___________  	
  
  predict.corr.scca[is.na(predict.corr.scca)] <- 0
  predict.corr.scca <- predict.corr.scca/n.cv
  
  best.predict.corr.scca <- max(abs(predict.corr.scca), na.rm=TRUE)
  best.lambda.v <- lambda.v.matrix[predict.corr.scca==best.predict.corr.scca]
  best.lambda.u <- lambda.u.matrix[predict.corr.scca==best.predict.corr.scca]  

  best.lambda.u <- best.lambda.u[1]
  best.lambda.v <- best.lambda.v[1]
  
  
  k <- sample.sigma12.function(x.data, y.data)
  
  # Get starting values for singular vectors
  # as column and row means from matrix K
  u.initial <- k %*% ones.q
  u.initial <- u.initial /sqrt(as.numeric(t(u.initial)%*%u.initial))
  v.initial <- t(k) %*% ones.p
  v.initial <- v.initial /sqrt(as.numeric(t(v.initial)%*%v.initial))
  
  uv <- scca.function(k, u.initial, v.initial, best.lambda.u, best.lambda.v)
  
  vj <- uv$v.new  # sparse singular vector (canonical vector for Y)
  uj <- uv$u.new	# sparse singular vector (canonical vector for X)
  
  corr.scca <- abs(cor(x.data%*%uj, y.data%*%vj))	# canonical correlation for X and Y data
  
  # Canoncial vectors
  a<-diag(1/apply(x.data,2,sd),p)%*%uj
  b<-diag(1/apply(y.data,2,sd),q)%*%vj
  
  # rank 1 approximation of K matrix
  K_hat<-a%*%corr.scca%*%t(b)
  K_residual<-k-K_hat
  
  out<-list(v=vj,u=uj,cancor=corr.scca,a=a,b=b,K_hat=K_hat,K_residual=K_residual,k=k,best.lambda.u=best.lambda.u,best.lambda.v=best.lambda.v)
}

SCCA_Parkhomenko_Khigher<-function(x.data,y.data,k,best.lambda.u,best.lambda.v){ # Sparse CCA for higher-order dimensions
  ### INPUT
  # x.data:         1st datamatrix n times p
  # y.data:         2nd datamatrix n times q
  # k:              K-K_hat
  # best.lambda.u:  optimal sparseness parameter for matrix X
  # best.lambda.v:  optimal sparseness parameter for matrix Y
  
  ### OUTPUT
  # v:              right singular vector associated to matrix Y
  # u:              left singular vector associated to matrix X
  # cancor:         first canonical correlation
  # a:              canonical vectors accosiated to matrix X
  # b:              Canonical vectors associated to matrix Y
  # k:              matrix K (sample covariance matrix)
  
  # ANALYSIS
  # Dimensions of input-dataset
  p<-ncol(x.data) # Number of variables in X
  q<-ncol(y.data) # Number of variables in Y
  ones.p <- rep(1, p)/p
  ones.q <- rep(1, q)/q
  
  # Get starting values for singular vectors
  # as column and row means from matrix K
  u.initial <- k %*% ones.q
  u.initial <- u.initial /sqrt(as.numeric(t(u.initial)%*%u.initial))
  v.initial <- t(k) %*% ones.p
  v.initial <- v.initial /sqrt(as.numeric(t(v.initial)%*%v.initial))
  
  uv <- scca.function(k, u.initial, v.initial, best.lambda.u, best.lambda.v)
  
  vj <- uv$v.new  # sparse singular vector (canonical vector for Y)
  uj <- uv$u.new  # sparse singular vector (canonical vector for X)
  
  corr.scca <- abs(cor(x.data%*%uj, y.data%*%vj))	# canonical correlation for X and Y data
  
  
  a<-diag(1/apply(x.data,2,sd),p)%*%uj
  b<-diag(1/apply(y.data,2,sd),q)%*%vj
  
  
  out<-list(v=vj,u=uj,cancor=corr.scca,a=a,b=b,k=k)
}


####################################################################################
####################### Code From Parkhomenko (2009) ###############################
####################################################################################
# Algorithm for Sparse Canonical Correlation Analysis
# Parkhomenko, Elena; Tritchler, David; and Beyene, Joseph (2008) "Sparse Canonical Correlation Analysis with Application to Genomic Data Integration," Statistical Applications in Genetics and Molecular Biology: Vol. 7 : Iss. 1, Article 1.  
# Available at: http://www.bepress.com/sagmb/vol7/iss1/1 


# This function computes sparse version of the first singular vectors of matrix K

# The number of selected variables, i.e. variables with non-zero entries 
# in computed singular vectors, is controlled by the sparseness parameters. 
# Increasing the sparseness parameter will decrease the number of selected 
# variables.


# Parameters:
# k is p by q covariance matrix for standardized data sets X and Y
# u.initial p by 1 vector of starting values for left singular vector
# v.initial q by 1 vector of starting values for right singular vector
# lambda.u is the sparseness parameter for left singular vector
# lambda.v is the sparseness parameter for right singular vector

scca.function <- function(k, u.initial, v.initial, lambda.u, lambda.v)
{
  i <- 0		# number of iterations used by SCCA
  eps <- 0.001	# convergence criterion 
  max.iter <- 50	# maximum nuber of iterations
  diff.u <- eps*10	
  diff.v <- eps*10
  
  while ((i < max.iter) & ((diff.u > eps) || (diff.v > eps)) )
  {
    i <- i+1
    
    # Update left singular vector
    
    vx <-  k %*% v.initial
    length.vx <- as.numeric(sqrt(t(vx)%*%vx))
    if(length.vx==0) length.vx <- 1
    vx <- vx / length.vx
    u.new <- abs(vx) - 0.5*lambda.u
    u.new <- (u.new + abs(u.new))/2
    u.new <- u.new*sign(vx)
    length.u.new <- as.numeric(sqrt(t(u.new)%*%u.new))
    if(length.u.new==0) length.u.new <- 1
    u.new <- u.new / length.u.new
    
    
    # Update right singular vector
    
    ux <- t(k) %*% u.new
    length.ux <- as.numeric(sqrt(t(ux)%*%ux))
    if(length.ux==0) length.ux <- 1
    ux <- ux / length.ux
    v.new <- abs(ux) - 0.5*lambda.v
    v.new <- (v.new + abs(v.new))/2
    v.new <- v.new * sign(ux)
    length.v.new <- as.numeric(sqrt(t(v.new)%*%v.new))
    if(length.v.new==0) length.v.new <- 1
    v.new <- v.new / length.v.new
    
    
    # Convergence measures
    
    diff.v <- max(abs(v.initial - v.new))
    diff.u <- max(abs(u.initial - u.new))
    
    v.initial <- v.new
    u.initial <- u.new
  }
  
  # Report the results:
  # u.new is computed left singular vector
  # v.new is computed right singular vector
  # i is the number of iterations used by SCCA
  
  list(u.new=u.new, v.new=v.new, i=i)
}


############  SOURCE FOR COVARIANCE MATRIX FUNCTION   #############################
# see source("sample_cov_function.R") above
###################################################################################
# Calculating sample covariance function
sample.sigma12.function <- function(x, y)
{
  x <- x - mean(x)
  y <- y - mean(y)
  
  # Sample variance-covariance matrices 
  sigma11 <- var(x)
  sigma22 <- var(y)
  
  x <- x %*% diag( 1/sqrt(diag(sigma11)) )
  y <- y %*% diag( 1/sqrt(diag(sigma22)) )
  
  return(cov(x,y))
  
  #   xstdize<-stdize(x)
  #   ystdize<-stdize(y)
  #   return(cov(xstdize,ystdize))
  
}