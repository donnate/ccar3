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
  norm_A <- sqrt(sum(A^2))
  if (norm_A == 0) return(A)
  A * pmax(0, 1 - lambda / norm_A)
}



soft_thresh2 <- function(A, lambda){
  norm_A <- sqrt(sum(A^2))
  if(norm_A == 0){
    return(A)
  }
  result = A * pmax(1 - lambda/(norm_A), 0)
  return(result)
}

matmul = function(A, B){
  SMUT::eigenMapMatMult(A, B)
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
#' @param standardize Whether to scale variables
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
ecca = function(X, Y, lambda = 0, groups = NULL, Sx = NULL,
                Sy = NULL, Sxy = NULL, r = 2,  standardize = F,
                rho = 1, B0 = NULL, eps = 1e-4, maxiter = 500, verbose = T){
  
  p = ncol(X)
  q = ncol(Y)
  n = nrow(X)
  ##### Need to make sure that they are centered
  if (standardize) {
    X = scale(X)
    Y = scale(Y)
  } 
  
  if (is.null(Sxy)) Sxy = matmul(t(X), Y)/ n
  if (is.null(Sx)) Sx = matmul(t(X), X) / n
  if (is.null(Sy)) {
    Sy = matmul(t(Y), Y) / n
  }
  
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
    Btilde = B1 + rho * (t(Ux) %*% ( Z - H) %*% Uy) 
    Btilde = Btilde / b
    B = ((Ux %*% Btilde) %*% t(Uy))
    
    # Update Z
    
    Z = B + H
    if(is.null(groups)){
      Z = soft_thresh(Z, lambda / rho)
    }
    else{
      for (g in seq_along(groups)) {
  
        # Get the indices for the current group
        current_indices <- groups[[g]]
        
        # 1. Subset Z only ONCE
        Z_subset <- Z[current_indices]
        
        # 2. Calculate the correctly scaled lambda for this group
        #    Using nrow() is correct for your 2-column coordinate matrix
        lambda_g <- sqrt(nrow(current_indices)) * lambda / rho
        
        # 3. Apply the single, robust thresholding function
        thresholded_values <- soft_thresh2(Z_subset, lambda_g)
        
        # 4. Assign the result back
        Z[current_indices] <- thresholded_values
      }


      # for (g in 1:length(groups)){
      #   Z[groups[[g]] ] =  soft_thresh2(Z[groups[[g]] ], sqrt(length(groups[[g]]) ) * lambda/rho)
      # }
    }
    
    # Update H
    
    H = H + rho * (B - Z)
    
    sB0 <- sum(B0^2)
    if(sB0 > 1e-20) { # Use a small tolerance instead of > 0 for numerical stability
      delta <- sum((B - B0)^2) / sB0
    } else {
      delta <- Inf
    }
    
    
    #if(fnorm(B0) > 0) delta = fnorm(B - B0)^2 / fnorm(B0)^2
    if(verbose && iter %% 10 == 0) cat("\niter:", iter, "delta:", delta)
  }
  if(iter >= maxiter) cat(crayon::red("     ADMM did not converge!"))
  else cat(paste0(crayon::green("     ADMM converged in ", iter, " iterations")))
  
  # Step 2: map back
  
  B = Z
  C = matmul(matmul(Sx12, B), Sy12) 
  
  # SVD = svd(C)
  # U0 = SVD$u[, 1:r, drop = F]
  # V0 = SVD$v[, 1:r, drop = F]
  # L0 = SVD$d[1:r]
  
  SVD = RSpectra::svds(C, r)
  U0 = SVD$u
  V0 = SVD$v
  L0 = SVD$d
  
  inv_L0 <- sapply(L0, function(d) ifelse(d > 1e-8, 1/d, 0))
  
  if(max(L0) > 1e-8){
    U = matmul(matmul(matmul(B, Sy12),V0), diag(inv_L0, nrow = length(L0)))
    V = matmul(matmul(matmul(t(B), Sx12), U0), diag(inv_L0, nrow = length(L0)))
    return(list(U = U, V = V, loss = mean(apply((matmul(X,U) - matmul(Y, V))^2,2,sum)), cor = diag(matmul(t(U), matmul(Sxy, V)))))
  } else{
    U = matrix(NA, p, r)
    V = matrix(NA, q, r)
    return(list(U = U, V = V, loss = mean(apply((matmul(X, U) - matmul(Y,V))^2,2,sum)), cor = rep(0, r)))
  }
  
  
}





ecca_across_lambdas = function(X, Y, lambdas = 0, groups = NULL, r = 2,  Sx = NULL,
                               Sy = NULL, Sxy = NULL, standardize = T, 
                               rho = 1, B0 = NULL, eps = 1e-4, maxiter = 500, verbose = T){
  
  p = ncol(X)
  q = ncol(Y)
  n = nrow(X)
  
  ##### Need to make sure that they are centered
  if (standardize) {
    X = scale(X)
    Y = scale(Y)
  } 
  
  
  
  if (is.null(Sxy)) Sxy = matmul(t(X), Y)/ n
  if (is.null(Sx)) Sx = matmul(t(X), X) / n
  if (is.null(Sy)) {
    Sy = matmul(t(Y), Y) / n
  }
  
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
  H = matrix(0, p, q)
  
  for(i in 1:length(lambdas)) {
    lambda = lambdas[[i]]
    
    # Step 1: ADMM
    
    Z = B
    iter = 0
    delta = Inf
    
    while(delta > eps && iter < maxiter){
      iter = iter + 1
      
      # Update B
      
      B0 = B
      Btilde = B1 + rho * (t(Ux) %*% ( Z - H) %*%  Uy) 
      Btilde = Btilde / b
      B = (Ux %*% Btilde) %*% t(Uy)
      
      # Update Z
      
      Z = B + H
      if(is.null(groups)){
        Z = soft_thresh(Z, lambda / rho)
      }
      else{
        for (g in seq_along(groups)) {
  
          # Get the indices for the current group
          current_indices <- groups[[g]]
          
          # 1. Subset Z only ONCE
          Z_subset <- Z[current_indices]
          
          # 2. Calculate the correctly scaled lambda for this group
          #    Using nrow() is correct for your 2-column coordinate matrix
          lambda_g <- sqrt(nrow(current_indices)) * lambda / rho
          
          # 3. Apply the single, robust thresholding function
          thresholded_values <- soft_thresh2(Z_subset, lambda_g)
          
          # 4. Assign the result back
          Z[current_indices] <- thresholded_values
        }

        # for (g in 1:length(groups)){
        #   Z[groups[[g]] ] =  soft_thresh2(Z[groups[[g]] ], sqrt(length(groups[[g]]) ) * lambda/rho)
        # }
      }
      
      # Update H
      
      H = H + rho * (B - Z)
      sB0 <- sum(B0^2) 
      if(sB0 > 1e-20) { # Use a small tolerance instead of > 0 for numerical stability
        delta <- sum((B - B0)^2) / sB0
      } else {
        delta <- Inf
      }
      if(verbose && iter %% 10 == 0) cat("\niter:", iter,  "delta:", delta)
    }
    if(iter >= maxiter) cat(crayon::red("     ADMM did not converge!"))
    else cat(paste0(crayon::green("     ADMM converged in ", iter, " iterations")))
    
    # Step 2: map back
    
    B = Z
    C = matmul(matmul(Sx12, B), Sy12) 
    
    # SVD = svd(C)
    # U0 = SVD$u[, 1:r, drop = F]
    # V0 = SVD$v[, 1:r, drop = F]
    # L0 = SVD$d[1:r]
    
    SVD = RSpectra::svds(C, r)
    U0 = SVD$u
    V0 = SVD$v
    L0 = SVD$d
    
    
    inv_L0 <- sapply(L0, function(d) ifelse(d > 1e-8, 1/d, 0))
    
    if(max(L0) > 1e-8){
      U[[i]] = matmul(matmul(matmul(B, Sy12),V0), diag(inv_L0, nrow = length(L0)))
      V[[i]] = matmul(matmul(matmul(t(B), Sx12), U0), diag(inv_L0, nrow = length(L0)))
    } else{
      U[[i]] = matrix(NA, p, r)
      V[[i]] = matrix(NA, q, r)
    }
  }
  if(length(lambdas) == 1) return(list(U = U[[1]], V = V[[1]]))
  else return(list(U = U, V = V))
}



ecca.eval = function(X, Y,  lambdas = 0, groups = NULL, r = 2, 
                     standardize = T, Sx = NULL, Sy = NULL, Sxy = NULL,
                     rho = 1, B0 = NULL, nfold = 5, eps = 1e-4,
                     maxiter = 500, verbose = T, parallel = T,
                     nb_cores = NULL, set_seed_cv=NULL){
  p = ncol(X)
  q = ncol(Y)
  n = nrow(X)
  if (n < 2 * nfold) {
    cat(crayon::yellow(paste0("\nWarning: Sample size (n=", n, ") is too small for ", nfold, "-fold CV. Skipping CV and fitting with the first lambda value.")))
    lambda.min = lambdas[1]
    lambda.1se =lambdas[1]
    scores = NULL
  } else {
    
    ##### Need to make sure that they are centered
    if (standardize) {
      X = scale(X)
      Y = scale(Y)
    } 
    if (is.null(Sxy)) Sxy = matmul(t(X), Y)/ n
    if (is.null(Sx)) Sx = matmul(t(X), X) /n
        if (is.null(Sy)) {
          Sy = matmul(t(Y), Y) / n
        }
    
    
    
    
    ## Create folds
    if (is.null(set_seed_cv) == FALSE){
       set.seed(set_seed_cv)
    }
   

    
    folds = caret::createFolds(1:n, k = nfold, list = T)
    n_success = nfold
    ## Choose penalty lambda
    results <- data.frame(lambda = numeric(), mse = numeric(), se = numeric())
    
    if (parallel) {
    # --- GRACEFUL PARALLEL SETUP ---
      cl <- setup_parallel_backend(nb_cores)
      
      if (!is.null(cl)) {
        # If the cluster was created successfully, register it and plan to stop it
        cat(crayon::green("Parallel backend successfully registered.\n"))
        doParallel::registerDoParallel(cl)
        on.exit(parallel::stopCluster(cl), add = TRUE)
      } else {
        # If setup_parallel_backend returned NULL, print a warning and proceed serially
        warning("All parallel setup attempts failed. Proceeding in serial mode.", immediate. = TRUE)
        parallel <- FALSE # Ensure %dopar% runs serially
    }
   }

   if (parallel) {
  
      ## Parallel cross validation
      cv = foreach(fold = folds, 
                   .export = c("ecca", "ecca_across_lambdas", "matmul", "fnorm", "soft_thresh", "soft_thresh_group", "soft_thresh2"), 
                   .packages = c("SMUT", "rARPACK", "crayon"),
                   .errorhandling = 'pass') %dopar% {
                     
                     n_full <- nrow(X)
                     X_train <- X[-fold, ]
                     Y_train <- Y[-fold, ]
                     X_val <- X[fold, ]
                     Y_val <- Y[fold, ]
                     
                     # We don't even need to create X_train and Y_train explicitly for the covariance calculation
                     n_train <- n_full - nrow(X_val)
                     
                     # --- DOWNDATE TRICK ---
                     Sx_train <- (n_full * Sx - crossprod(X_val)) / n_train
                     Sy_train <- (n_full * Sy - crossprod(Y_val)) / n_train
                     Sxy_train <- (n_full * Sxy - crossprod(X_val, Y_val)) / n_train
                     
                     ## Fit lasso model
                     ECCA = ecca_across_lambdas(X[-fold,,drop = F], Y[-fold,,drop = F], 
                                                Sx = Sx_train, Sy = Sy_train, Sxy = Sxy_train,
                                                standardize = F,
                                                lambdas=lambdas, groups=groups, r=r, rho=rho, 
                                                B0=B0, eps=eps, maxiter=maxiter, verbose = verbose)
                     
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
                       
                       X_val_mat <- X[fold, , drop = FALSE]
                       Y_val_mat <- Y[fold, , drop = FALSE]
                       
                       #if(is.na(U)) scores[i] = NA
                       if(!is.null(U) && !any(is.na(U))) {
                         scores[i] = mean((matmul(X_val_mat, U) - Y_val_mat %*%V)^2)
                       } else {
                         scores[i] = Inf
                       }
            
                     }
                     return(scores)
                   }
      
      
      # --- Post-processing to handle potential errors ---
      # Identify which folds resulted in an error
      is_error <- sapply(cv, function(x) inherits(x, "error"))
      
      if(any(is_error)){
        warning(paste(sum(is_error), "out of", nfold, "folds failed during cross-validation."))
        # Optional: print the actual error messages for debugging
        # print(cv_results_list[is_error])
      }
      
      # Filter out the errors and combine the successful results
      successful_results <- cv[!is_error]
      
      # If all folds failed, we can't proceed
      if(length(successful_results) == 0) {
        stop("All CV folds failed. Cannot select a lambda.")
      }
      
      scores.cv = do.call(cbind, successful_results)
      n_success <- length(successful_results)
    } else {
      scores.cv = c()
      for(i in 1:length(folds)){
        
        if(verbose) print(paste0("\n\nfold:", i))
        fold = folds[[i]]
        
        n_full <- nrow(X)
        X_train <- X[-fold, ]
        Y_train <- Y[-fold, ]
        X_val <- X[fold, ]
        Y_val <- Y[fold, ]
        
        # We don't even need to create X_train and Y_train explicitly for the covariance calculation
        n_train <- n_full - nrow(X_val)
        
        # --- DOWNDATE TRICK ---
        Sx_train <- (n_full * Sx - crossprod(X_val)) / n_train
        Sy_train <- (n_full * Sy - crossprod(Y_val)) / n_train
        Sxy_train <- (n_full * Sxy - crossprod(X_val, Y_val)) / n_train
        
        
        ## Fit lasso model
        ECCA = ecca_across_lambdas(X[-fold,,drop = F], Y[-fold,,drop = F], 
                                   Sx = Sx_train, Sy = Sy_train, Sxy = Sxy_train,
                                   standardize = F,
                                   lambdas = lambdas, groups = groups, r = r, 
                                   rho = rho, B0 = B0, eps= eps, maxiter = maxiter, verbose = verbose)
        
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
          #if(is.na(U)) scores[i] = NA
          X_val_mat <- X[fold, , drop = FALSE]
          Y_val_mat <- Y[fold, , drop = FALSE]
          
          if(!is.null(U) && !any(is.na(U))) {
            scores[i] = mean((matmul(X_val_mat, U) - Y_val_mat %*%V)^2)
          } else {
            scores[i] = Inf
          }
          
        }
        if(verbose) print(paste("\n\nMSEs:", scores))
        scores.cv = cbind(scores.cv, scores)
      }
    }
    scores = data.frame(lambda = lambdas, mse = rowMeans(scores.cv), 
                        se = matrixStats::rowSds(scores.cv)/sqrt(n_success))
    # plt = scores %>%
    #   ggplot(aes(lambda, mse))+
    #   geom_point()+
    #   geom_line()+
    #   geom_errorbar(aes(ymin = mse - se, ymax = mse + se))+
    #   theme_bw()
    # suppressWarnings(print(plt))
    lambda.min = scores %>% slice(which.min(mse)) %>% pull(lambda)
    upper = scores %>% slice(which.min(mse)) %>% mutate(upper = mse + se) %>% pull(upper)
    lambda.1se = scores %>% filter(mse <= upper) %>% pull(lambda) %>% max()
  }
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
#' @param standardize Whether to scale variables
#' @param B0 Initial value for the coefficient matrix (optional)
#' @param eps Convergence threshold for ADMM
#' @param rho ADMM parameter
#' @param maxiter Maximum number of ADMM iterations
#' @param verbose Print diagnostics
#' @param nb_cores Number of cores to use for parallel processing (default is NULL, which uses all available cores)
#' @param set_seed_cv Optional seed for reproducibility of cross-validation folds (default is NULL)
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
ecca.cv = function(X, Y, lambdas = 0, groups = NULL, r = 2, standardize = F,
                   rho = 1, B0 = NULL, nfold = 5, select = "lambda.min", eps = 1e-4, maxiter = 500, 
                   verbose = F, parallel = F,
                   nb_cores = NULL,
                   set_seed_cv=NULL){
  p = ncol(X)
  q = ncol(Y)
  n = nrow(X)
  
  if (standardize) {
    X = scale(X)
    Y = scale(Y)
  } 
  
  # Select lambda
  if(length(lambdas) > 1){
    eval = ecca.eval(X, Y, lambdas=lambdas, groups=groups, r=r, rho=rho,
                     standardize = F,
                     B0 = B0, nfold=nfold, eps=eps,  maxiter=maxiter, verbose=verbose, parallel= parallel,
                     nb_cores = nb_cores, set_seed_cv=set_seed_cv)
    if(select == "lambda.1se") lambda.opt = eval$lambda.1se
    else lambda.opt = eval$lambda.min
  } else {
    lambda.opt = lambdas
  }
  cat("\n\nselected lambda:", lambda.opt)
  
  # Fit lasso
  ECCA = ecca(X, Y, lambda=lambda.opt, groups = groups, r=r, rho=rho, B0=B0, eps=eps, 
              standardize = FALSE,
              maxiter = maxiter, verbose = verbose)
  
  return(list(U = ECCA$U, V = ECCA$V, 
              cor = diag(matmul(matmul(t(ECCA$U), matmul(t(X), Y)), ECCA$V)), 
              loss = mean(apply((matmul(X , ECCA$U) - matmul(Y , ECCA$V))^2,2,sum)),
              lambda.opt  = lambda.opt,
              cv.scores = eval$scores))
}
