
scca_chao <- function(X, Y, rho = 1, lambda_max = .1, num_lambda = 10, r = 2, niter = 500,
                     nfold = 8, thresh = .01){
  p = dim(X)[2]
  q = dim(Y)[2]
  n = dim(X)[1]
  
  folds = createFolds(1:n, k = 3, list = T )
  
  result1 = Fantope(X[folds$Fold1, ], Y[folds$Fold1, ], r)
  
  result2_u = cv_group_lasso(X[folds$Fold2, ], Y[folds$Fold2,] %*% result1$v, rho, lambda_max, num_lambda, r, niter, nfold, thresh)
  result2_v = cv_group_lasso(Y[folds$Fold2, ], X[folds$Fold2,] %*% result1$u, rho, lambda_max, num_lambda, r, niter, nfold, thresh)
  
  u = result2_u %*% pracma::sqrtm( t(result2_u) %*% cov(X[folds$Fold3, ]) %*% result2_u)$Binv
  v = result2_v %*% pracma::sqrtm( t(result2_v) %*% cov(Y[folds$Fold3, ]) %*% result2_v)$Binv
  
  return(list(u = u, v = v))
  
  
  
}

Fantope <- function(X, Y, r = 2, rho=NULL){
  p = ncol(X)
  q = ncol(Y)
  n = nrow(X)
  if (is.null(rho)){
    rho = 0.5 * sqrt(log( p + q)/n)
  }
  

  Mask <- matrix(0, (p + q), (p + q))
  idx1 <- 1:p
  idx2 <- (p+ 1):(p + q)
  Mask[idx1, idx1] <- matrix(1,p, p)
  Mask[idx2, idx2] <- matrix(1, q, q)
  
  Data = cbind(X, Y)
  S = cov(Data)
  sigma0hat <- S * Mask
  
  ag <- sgca_init(A=S, B=sigma0hat, rho=rho,
                  K=r,  maxiter=1000, trace=FALSE)
  ainit <- init_process(ag$Pi, r) 
  test1 <- gca_to_cca(ainit, S, c(p, q))
  
  return(test1)
}

cv_Fantope <- function(X, Y, rho = 1, lambda_max = .1, num_lambda = 20, r = 2, niter = 500,
                           nfold = 8, thresh = .01){
  
  ## Create folds
  n = nrow(X)
  cv = createFolds(1:n, k = nfold, list = T )
  
  ## choose penalty lambda
  lambda_values <- rho * seq(from =0, to = 1, by = 1/num_lambda) * lambda_max
  cv_results <- data.frame(lambda = numeric(), mse = numeric(), std = numeric())
  
  ## Cross validation
  for(lambda in lambda_values) {
    mse_values = NULL
    
    ## Parallelly run cross validation
    results <- foreach(fold = cv, .export = c("group_lasso"), .packages = "SMUT") %dopar% {
      
      temp <- Fantope(X, Y, r = 2, rho = lambda)
      
      sum((X %*% temp$u - Y %*% temp$v)^2)
      
    }
    
    # Store average MSE for this lambda
    cv_results <- rbind(cv_results, data.frame(lambda = lambda, mse = mean(unlist(results)),
                                               std = sd(unlist(results)) ))
    # Break the cv process if lambda is so large that B = 0
    
    
  }
  
  min_ind = which(cv_results$mse == min(cv_results$mse) )
  lambda = cv_results$lambda[min_ind]
  print(lambda)
  B = group_lasso(X, Y,  lambda, rho, niter, thresh)
  return(B)
  
}



cv_group_lasso <- function(X, Y, rho = 1, lambda_max = .1, num_lambda = 20, r = 2, niter = 500,
                          nfold = 8, thresh = .01){

  ## Create folds
  n = nrow(X)
  cv = createFolds(1:n, k = nfold, list = T )
  
  ## choose penalty lambda
  lambda_values <- rho * seq(from =0, to = 1, by = 1/num_lambda) * lambda_max
  cv_results <- data.frame(lambda = numeric(), mse = numeric(), std = numeric())
  
  ## Cross validation
  for(lambda in lambda_values) {
    mse_values = NULL
    
    ## Parallelly run cross validation
    results <- foreach(fold = cv, .export = c("group_lasso"), .packages = "SMUT") %dopar% {
      
      soft_thresh <- function(A, lambda){
        result = sign(A) * pmax(abs(A) - lambda, 0)
        return(result)
      }
      soft_thresh2 <- function(A, lambda){
        result = A * pmax(1 - lambda/(sqrt(sum(A^2))), 0)
        return(result)
      }
      
      test_indices <- fold
      train_indices <- setdiff(1:n, test_indices)
      
      ## Fit lasso model
      B <- group_lasso(X[train_indices, ], Y[train_indices, ],  lambda, rho, niter, thresh)
      
      sum((X %*% B - Y)^2)
      
    }
    
    # Store average MSE for this lambda
    cv_results <- rbind(cv_results, data.frame(lambda = lambda, mse = mean(unlist(results)),
                                               std = sd(unlist(results)) ))
    # Break the cv process if lambda is so large that B = 0
    
    
  }
  
  min_ind = which(cv_results$mse == min(cv_results$mse) )
  lambda = cv_results$lambda[min_ind]
  print(lambda)
  B = group_lasso(X, Y,  lambda, rho, niter, thresh)
  return(B)
  
}

group_lasso <- function(X, Y_tilde, lambda = 0, rho = 1, niter = 500, thresh = 0.001){
  n = nrow(X)
  p = ncol(X)
  q = ncol(Y_tilde)
  U <- matrix(0, p, q)
  Z <- matrix(0, p, q)
  prod_xy = t(X) %*% Y_tilde/ n
  Sx = cov(X)
  invSx = solve(Sx + rho * diag(rep(1, p)))
  
  for (i in 1:niter){
    #print(paste0("iter is ", i))
    Uold = U
    Zold = Z
    B = invSx %*% (prod_xy  + rho * (Z - U))
    Bold = B
    Z = B + U
    norm_col = sapply(1:nrow(Z), function(i){sqrt(sum(Z[i,]^2))})
    index_0 = which(norm_col < lambda/rho)
    index_pos = which(norm_col  > lambda/rho)
    if(length(index_0)>0){
      Z[index_0,] = 0
    }
    if (length(index_pos) > 1){
      Z[index_pos,] = diag(sapply(index_pos, function(x){ 1- (lambda/rho)/norm_col[x]})) %*% Z[index_pos,] 
    }else{
      if (length(index_pos) == 1){
        Z[index_pos,] = ( 1- (lambda/rho)/norm_col[index_pos]) * Z[index_pos,]
      }
    }
    U = U + B - Z
    if (max(c(norm(Z - B)/ sqrt(p), norm(Zold - Z)/ sqrt(p))) <thresh){
      break
    }
  }
  
  return(Z)
  
}
