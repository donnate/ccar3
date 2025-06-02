library(foreach)
library(caret)
library(doParallel)
library(dplyr)
library(rARPACK)
library(SMUT)
library(crayon)
library(matrixStats)

soft_thresh = function(A, lambda){
  #A * pmax(1 - lambda / abs(A), 0)
  sign(A) * pmax(abs(A) - lambda, 0)
}

fnorm = function(A){
  sqrt(sum(A^2))
}

soft_thresh_group = function(A, lambda){
  if(fnorm(A) == 0) return(A)
  A * pmax(1 - lambda / fnorm(A), 0)
}

soft_thresh2 <- function(A, lambda){
  if(sum(A^2) == 0){
    return(A)
  }
  result = A * pmax(1 - lambda/(sqrt(sum(A^2))), 0)
  return(result)
}

matmul = function(A, B){
  eigenMapMatMult(A, B)
}

rmat = function(n, p){
  matrix(rnorm(n * p), n, p)
}


#' Sparse Canonical Correlation via Reduced-Rank Regression when both X and Y are high-dimensiona;
#'
#' Performs group-sparse reduced-rank regression for CCA using either ADMM or CVXR solvers.
#'
#' @param X Predictor matrix (n x p)
#' @param Y Response matrix (n x q)
#' @param groups List of index vectors defining groups of predictors
#' @param lambda Regularization parameter
#' @param r Target rank
#' @param do.scale Whether to scale variables
#' @param B0 Initial value for the coefficient matrix (optional)
#' @param eps Convergence threshold for ADMM
#' @param rho ADMM parameter
#' @param maxiter Maximum number of ADMM iterations
#' @param verbose Print diagnostics
#'
#' @return A list with elements:
#' \describe{
#'   \item{U}{Canonical direction matrix for X (p x r)}
#'   \item{V}{Canonical direction matrix for Y (q x r)}
#'   \item{cor}{Canonical covariances}
#'   \item{loss}{The prediction error 1/n * \| XU - YV\|^2}
#' }
#' @export
ecca = function(X, Y, lambda = 0, groups = NULL, r = 2,  
                rho = 1, B0 = NULL, eps = 1e-4, maxiter = 500, verbose = T){
  
  p = ncol(X)
  q = ncol(Y)
  n = nrow(X)
  ##### Need to make sure that they are centered
  X = scale(X, center = T, scale = F)
  Y = scale(Y, center = T, scale = F)
  
  Sxy = matmul(t(X), Y)/ n
  Sx = matmul(t(X), X) / n
  Sy = matmul(t(Y), Y) / n
  
  if(is.null(B0)) B = matrix(0, p, q)
  else B = B0

  
  EDx = eigen(Sx, symmetric = T)
  EDy = eigen(Sy, symmetric = T)
  
  Ux = EDx$vectors
  Lx = EDx$values
  Uy = EDy$vectors
  Ly = EDy$values
  
  Sx12 = matmul(matmul(Ux, diag(sqrt(pmax(Lx, 0)))),  t(Ux)) 
  Sy12 = matmul(matmul(Uy, diag(sqrt(pmax(Ly, 0)))),  t(Uy)) 
  
  b = outer(Lx, Ly) + rho
  B1 = matmul(matmul(t(Ux), Sxy),  Uy)
  U = list()
  V = list()
  
       
    
  # Step 1: ADMM
  H = matrix(0, p, q)
  Z = B
  iter = 0
  delta = Inf

  while(delta > eps && iter < maxiter){
        iter = iter + 1
        
        # Update B
        
        B0 = B
        Btilde = B1 + rho * matmul(matmul(t(Ux), Z - H),  Uy) 
        Btilde = Btilde / b
        B = matmul(matmul(Ux, Btilde), t(Uy))
        
        # Update Z
        
        Z = B + H
        if(is.null(groups)){
        Z = soft_thresh(Z, lambda / rho)
        }
        else{
        for (g in 1:length(groups)){
        Z[groups[[g]] ] =  soft_thresh2(Z[groups[[g]] ], sqrt(length(groups[[g]]) ) * lambda/rho)
        }
        }
        
        # Update H
        
        H = H + rho * (B - Z)
        
        if(fnorm(B0) > 0) delta = fnorm(B - B0)^2 / fnorm(B0)^2
        if(verbose && iter %% 10 == 0) cat("\niter:", iter, "fnorm:", fnorm(B0), "delta:", delta)
    }
    if(iter >= maxiter) cat(red("     ADMM did not converge!"))
    else cat(paste0(green("     ADMM converged in ", iter, " iterations")))

    # Step 2: map back

    B = Z
    C = matmul(matmul(Sx12, B), Sy12) 

    # SVD = svd(C)
    # U0 = SVD$u[, 1:r, drop = F]
    # V0 = SVD$v[, 1:r, drop = F]
    # L0 = SVD$d[1:r]

    SVD = svds(C, r)
    U0 = SVD$u
    V0 = SVD$v
    L0 = SVD$d

    if(max(L0) > 1e-8){
        U = matmul(matmul(matmul(B, Sy12),V0), diag(1 / L0, nrow = length(L0)))
        V = matmul(matmul(matmul(t(B), Sx12), U0), diag(1 / L0, nrow = length(L0)))
    } else{
        U = matrix(NA, p, r)
        V = matrix(NA, q, r)
    }
    
    return(list(U = U, V = V, loss = mean(apply((X %*% U - Y %*% V)^2,2,sum)), cor = diag(t(U) %*% matmul(Sxy, V))))
}





ecca_across_lambdas = function(X, Y, lambdas = 0, groups = NULL, r = 2,  
                rho = 1, B0 = NULL, eps = 1e-4, maxiter = 500, verbose = T){
  
  p = ncol(X)
  q = ncol(Y)
  n = nrow(X)
  
  Sxy = matmul(t(X), Y)/ n
  Sx = matmul(t(X), X) / n
  Sy = matmul(t(Y), Y) / n
  
  if(is.null(B0)) B = matrix(0, p, q)
  else B = B0

  
  EDx = eigen(Sx, symmetric = T)
  EDy = eigen(Sy, symmetric = T)
  
  Ux = EDx$vectors
  Lx = EDx$values
  Uy = EDy$vectors
  Ly = EDy$values
  
  Sx12 = matmul(matmul(Ux, diag(sqrt(pmax(Lx, 0)))),  t(Ux)) 
  Sy12 = matmul(matmul(Uy, diag(sqrt(pmax(Ly, 0)))),  t(Uy)) 
  
  b = outer(Lx, Ly) + rho
  B1 = matmul(matmul(t(Ux), Sxy),  Uy)
  U = list()
  V = list()
  
  for(i in 1:length(lambdas)) {
    lambda = lambdas[[i]]
   
    # Step 1: ADMM
    H = matrix(0, p, q)
    Z = B
    iter = 0
    delta = Inf
    
    while(delta > eps && iter < maxiter){
      iter = iter + 1
      
      # Update B
      
      B0 = B
      Btilde = B1 + rho * matmul(matmul(t(Ux), Z - H),  Uy) 
      Btilde = Btilde / b
      B = matmul(matmul(Ux, Btilde), t(Uy))
      
      # Update Z
      
      Z = B + H
      if(is.null(groups)){
        Z = soft_thresh(Z, lambda / rho)
      }
      else{
      for (g in 1:length(groups)){
        Z[groups[[g]] ] =  soft_thresh2(Z[groups[[g]] ], sqrt(length(groups[[g]]) ) * lambda/rho)
      }
      }
      
      # Update H
      
      H = H + rho * (B - Z)
      
      if(fnorm(B0) > 0) delta = fnorm(B - B0)^2 / fnorm(B0)^2
      if(verbose && iter %% 10 == 0) cat("\niter:", iter, "fnorm:", fnorm(B0), "delta:", delta)
    }
    if(iter >= maxiter) cat(red("     ADMM did not converge!"))
    else cat(paste0(green("     ADMM converged in ", iter, " iterations")))
    
    # Step 2: map back
    
    B = Z
    C = matmul(matmul(Sx12, B), Sy12) 
    
    # SVD = svd(C)
    # U0 = SVD$u[, 1:r, drop = F]
    # V0 = SVD$v[, 1:r, drop = F]
    # L0 = SVD$d[1:r]

    SVD = svds(C, r)
    U0 = SVD$u
    V0 = SVD$v
    L0 = SVD$d

    if(max(L0) > 1e-8){
      U[[i]] = matmul(matmul(matmul(B, Sy12),V0), diag(1 / L0, nrow = length(L0)))
      V[[i]] = matmul(matmul(matmul(t(B), Sx12), U0), diag(1 / L0, nrow = length(L0)))
    } else{
      U[[i]] = matrix(NA, p, r)
      V[[i]] = matrix(NA, q, r)
    }
  }
  if(length(lambdas) == 1) return(list(U = U[[1]], V = V[[1]]))
  else return(list(U = U, V = V))
}



ecca.eval = function(X, Y, lambdas = 0, groups = NULL, r = 2, 
                   rho = 1, B0 = NULL, nfold = 5, eps = 1e-4, maxiter = 500, verbose = T, parallel = T){
  p = ncol(X)
  q = ncol(Y)
  n = nrow(X)
  
  ## Create folds
  set.seed(0)
  folds = createFolds(1:n, k = nfold, list = T)
  
  ## Choose penalty lambda
  results = data.frame(lambda = numeric(), mse = numeric(), se = numeric())
  
  ## Cross validation
      if(parallel){
        registerDoParallel(makeCluster(detectCores() - 2, type = "PSOCK"))
        ## Parallel cross validation
        cv = foreach(fold = folds, .export = c("ecca", "matmul", "fnorm", "soft_thresh", "soft_thresh_group", "soft_thresh2"), .packages = c("SMUT", "rARPACK", "crayon")) %dopar% {
          
          ## Fit lasso model
          ECCA = ecca_across_lambdas(X[-fold,,drop = F], Y[-fold,,drop = F], lambdas, groups, r, rho, B0, eps, maxiter, verbose = verbose)
          
          ## Evaluate on test set
          scores = rep(0, length(lambdas))
          for(i in 1:length(lambdas)){
            if(length(lambdas) > 1){
              U = (ECCA$U)[[i]]
              V = (ECCA$V)[[i]]
            } else {
              U = ECCA$U
              V = ECCA$V
            }
            if(!is.null(U)) scores[i] = mean((X[fold,] %*% U - Y[fold,] %*% V)^2)
            else scores[i] = Inf
          }
          return(scores)
        }
        scores.cv = do.call(cbind, cv)
      } else {
        scores.cv = c()
        for(i in 1:length(folds)){
            
            if(verbose) print(paste0("\n\nfold:", i))
            fold = folds[[i]]
            
            ## Fit lasso model
            ECCA = ecca_across_lambdas(X[-fold,,drop = F], Y[-fold,,drop = F], lambdas, groups, r, rho, B0, eps, maxiter, verbose = verbose)
            
            ## Evaluate on test set
            scores = rep(0, length(lambdas))
            for(i in 1:length(lambdas)){
              if(length(lambdas) > 1){
                U = (ECCA$U)[[i]]
                V = (ECCA$V)[[i]]
              } else {
                U = ECCA$U
                V = ECCA$V
              }
              #if(is.na(U)) scores[i] = NA
              scores[i] = mean((X[fold,] %*% U - Y[fold,] %*% V)^2)
            }
            if(verbose) print(paste("\n\nMSEs:", scores))
            scores.cv = cbind(scores.cv, scores)
        }
    }
    scores = data.frame(lambda = lambdas, mse = rowMeans(scores.cv), se = rowSds(scores.cv)/sqrt(nfold))
    plt = scores %>%
      ggplot(aes(lambda, mse))+
      geom_point()+
      geom_line()+
      geom_errorbar(aes(ymin = mse - se, ymax = mse + se))+
      theme_bw()
    suppressWarnings(print(plt))
    lambda.min = scores %>% slice(which.min(mse)) %>% pull(lambda)
    upper = scores %>% slice(which.min(mse)) %>% mutate(upper = mse + se) %>% pull(upper)
    lambda.1se = scores %>% filter(mse <= upper) %>% pull(lambda) %>% max()
    return(list(scores = scores, lambda.min = lambda.min, lambda.1se = lambda.1se))
}

#' Sparse Canonical Correlation via Reduced-Rank Regression when both X and Y are high-dimensional, with Cross-Validation
#'
#' Performs group-sparse reduced-rank regression for CCA using either ADMM or CVXR solvers.
#'
#' @param X Predictor matrix (n x p)
#' @param Y Response matrix (n x q)
#' @param groups List of index vectors defining groups of predictors
#' @param lambdas Choice of regularization parameter 
#' @param r Target rank
#' @param nfold Number of cross-validation folds
#' @param select Which lambda to select: "lambda.min" or "lambda.1se"
#' @param do.scale Whether to scale variables
#' @param B0 Initial value for the coefficient matrix (optional)
#' @param eps Convergence threshold for ADMM
#' @param rho ADMM parameter
#' @param maxiter Maximum number of ADMM iterations
#' @param verbose Print diagnostics
#' @param parallel Whether to run cross-validation in parallel
#'
#' @return A list with elements:
#' \describe{
#'   \item{U}{Canonical direction matrix for X (p x r)}
#'   \item{V}{Canonical direction matrix for Y (q x r)}
#'   \item{cor}{Canonical covariances}
#'   \item{loss}{The prediction error 1/n * \| XU - YV\|^2}
#' }
#' @export
ecca.cv = function(X, Y, lambdas = 0, groups = NULL, r = 2, 
                 rho = 1, B0 = NULL, nfold = 5, select = "lambda.min", eps = 1e-4, maxiter = 500, verbose = F, parallel = F){
  p = ncol(X)
  q = ncol(Y)
  n = nrow(X)
  
  # Select lambda
  if(length(lambdas) > 1){
    eval = ecca.eval(X, Y, lambdas, groups, r, rho, B0, nfold, eps, maxiter, verbose, parallel)
    if(select == "lambda.1se") lambda.opt = eval$lambda.1se
    else lambda.opt = eval$lambda.min
  } else {
    lambda.opt = lambdas
  }
  cat("\n\nselected lambda:", lambda.opt)
  
  # Fit lasso
  ECCA = ecca_across_lambdas(X, Y, lambda.opt, groups, r, rho, B0, eps, maxiter, verbose = verbose)
  
  return(list(U = ECCA$U, V = ECCA$V, 
              cor = diag(t(ECCA$U) %*% matmul(t(X), Y) %*% ECCA$V), 
              loss = mean(apply((X %*% ECCA$U - Y %*% ECCA$V)^2,2,sum))))
}
