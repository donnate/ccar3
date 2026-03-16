library(pracma)
library(caret)


generate_example_sparse_U <- function(n, p1, p2,
                                      r_pca = 3,
                                      nnzeros = 5,
                                      theta = diag(c(0.9,  0.8)),
                                      lambda_pca = 1,
                                      r = 2, overlapping_amount = 0,
                                      normalize_diagonal = TRUE,
                                      n_new = 50000) {
  ###
  # This generates a dataset (X,Y) from a multivariate normal where Sigma_XX  and 
  # Sigma_YY have a correlation structure, and Sigma_XY = U Lambda V^T is rank r on a set of nnzeros rows, 0 elsewhere.
  # The number of rows (resp. Columns) on which Sigma_XX (resp. Sigma_YY)  and Sigma_XY
  # overlap is controlled by overlapping_amount
  # n: number of samples
  # p1: nb of features for X 
  # p2: nb of features for Y
  # nnzeros: number of non zero rows of U and V
  # r: rank of Sigma_XY
  # theta : canonical correlation (must be a list of length r)
  # r_pca: rank of the PCA (used in the creation of Sigma_XX and SigmaYY)
  # lambda_pca: also used to create Sigma_XX as Sigma_{XX} = U_X Lambda_pca U_X^T, and setting diag(Sigma_XX) to 1
  # Returns:
  # S : empirical covariance matrix X^TY
  # Sigma: underlying (population) covariance Sigma_{XY}
  # u: ground truth for u
  # v: ground truth for v
  ###
  p_tot <- p1 + p2
  pp <- c(p1, p2)
  #print("We have:")
  #print(c(p_tot, nnzeros <= min(p1, p2), nnzeros))
  #print('--------------------------------------')
  #print('Generating data ...')
  Sigma <- diag(p1 + p2)
  s <- (1):(nnzeros)
  #print(paste0('Number of non zeros is: ', nnzeros))
  start <- ceiling((1 - overlapping_amount) * nnzeros)
  s_pca  <- (start + 1) : (start + nnzeros)
  s_pca2  <- 1:p2
  Lambda_pca <- rep(lambda_pca, r_pca)
  # generate vcovariance matrix for X and Y
  
  if (r_pca > 0){
    u1 = matrix(0, p1, r_pca)
    u1[s_pca, ] <- matrix(runif(n = nnzeros * r_pca, max = 3, min=0), nrow = nnzeros, ncol = r_pca)
    # Normalize u1
    u1[s_pca,] <- 0.5 * diag(1/ sqrt(apply(u1[s_pca, ]^2, 1, sum))) %*% u1[s_pca,]
    u1[s_pca, 1:r_pca] <- u1[s_pca,1:r_pca] #%*% (pracma::sqrtm(t(u1[s_pca,1:r_pca]) %*% u1[s_pca, 1:r_pca])$Binv)
    T1 <- u1 %*% diag(Lambda_pca) %*% t(u1)
    if (normalize_diagonal){
      diag(T1) <- 1
      Sigma[1:p1, 1:p1] <- T1
    }else{
      Sigma[1:p1, 1:p1] <- Sigma[1:p1, 1:p1] + T1 
    } 
    ### Same for Sigma_y
    u2 <- matrix(runif( p2 * r_pca, max = 1, min=0), nrow=p2)
    u2 <- 0.5 * diag(1/ sqrt(apply(u2^2, 1, sum))) %*% u2
    #u2[s_pca2, 1:r_pca] <- u2[s_pca2, 1:r_pca] #%*% (pracma::sqrtm(t(u2[s_pca2, 1:r_pca]) %*% u2[s_pca2, 1:r_pca])$Binv)
    T2 <- u2 %*% t(u2)
    if (normalize_diagonal){
      diag(T2) <- 1
      Sigma[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)] <- T2
    }else{
      Sigma[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)] <- T2 + Sigma[(p1 + 1):(p1 + p2), (p1 + 1) : (p1 + p2)]
    }
  }
  Sigmax <- Sigma[1:p1,1:p1];
  Sigmay <- Sigma[(p1+1):p_tot,(p1+1):p_tot];
  
  ### Generate cross covariance
  Tss <- Sigmax[s,s]
  prod <- matrix(0, pp[1], r)
  prod[s, 1:r] <- as.matrix(runif( nnzeros * r,max = 3, min=1), nrow=nnzeros)  * as.matrix(sample(c(-1,1), nnzeros*r, replace=TRUE), nrow=nnzeros)
  u <- prod %*% (pracma::sqrtm(t(prod[s, 1:r]) %*% Tss %*% prod[s, 1:r])$Binv)
  
  
  Tss_v <- Sigma[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)] 
  prod_v <- matrix(0, pp[2], r)
  prod_v[, 1:r] <- as.matrix(runif( p2 * r , max = 3, min=1), nrow=p2)  * as.matrix(sample(c(-1,1), p2 * r, replace=TRUE), nrow=p2)
  v <- prod_v %*% (pracma::sqrtm(t(prod_v[, 1:r]) %*% Tss_v %*% prod_v[, 1:r])$Binv)

  
  #Sigma[(p1 + 1) :(p1 + p2), 1:p1] <- Sigmax %*% u %*% theta %*% t(v) %*% Sigmay
  Sigma[1:p1, (p1 + 1):(p1 + p2)] <- Sigmax %*% u %*% theta %*% t(v) %*% Sigmay
  Sigma[(p1 + 1) :(p1 + p2), 1:p1]  <-  Sigmay %*% v%*% theta %*% t(u) %*% Sigmax
  
  #Generate Multivariate Normal Data According to Sigma
  sqrt_Sigma <- pracma::sqrtm(Sigma)$B
  Data <- matrix(rnorm(n = (n + n_new) * p_tot), ncol=p_tot) %*% sqrt_Sigma
  X <- Data[1:n, 1:p1]
  Y <- Data[1:n, (p1 + 1):(p1 + p2)]
  Xnew <- Data[(n+1):(n_new + n), 1:p1]
  Ynew <- Data[(n+1):(n_new + n), (p1 + 1):(p1 + p2)]
  #print('Data generated.')
  #print('--------------------------------------')
  Mask <- matrix(0, p_tot, p_tot)
  idx1 <- 1:pp[1]
  idx2 <- (pp[1] + 1):(pp[1] + pp[2])
  Mask[idx1, idx1] <- matrix(1, pp[1], pp[1])
  Mask[idx2, idx2] <- matrix(1, pp[2], pp[2])
  Sigma0 <- Sigma * Mask
  S <- cov(Data)
  sigma0hat <- S * Mask
  # Generate ground truth canonical vectors
  #Sigma_X_inv <- pracma::sqrtm(Sigma[1:p1, 1:p1])$Binv
  #Sigma_Y_inv <-  pracma::sqrtm(Sigma[(p1+1):(p_tot), (p1+1):(p_tot)])$Binv
  #GT = svd(Sigma_X_inv %*% Sigma[1:p1, (p1 + 1):p_tot] %*% Sigma_Y_inv, nu = r, nv = r)
  return(list(Sigma=Sigma, Sigma0=Sigma0,
              S = S, sigma0hat =  sigma0hat, Mask= Mask,
              X=X, Y = Y, Data=Data, 
              u=u,#pracma::sqrtm(Sigma[1:p1, 1:p1])$Binv %*% GT$u, 
              v=v, #pracma::sqrtm(Sigma[(p1+1):(p1+p2), (p1+1):(p1+p2)])$Binv %*% GT$v, 
              Xnew = Xnew, Ynew=Ynew,
              Sigmax=Sigmax, Sigmay=Sigmay
  ))
  
}



generate_example_group <- function(n, p1, p2, 
                                  r_pca = 3,
                                  nnzeros = 5,
                                  theta = diag(c(0.9,  0.8)),
                                  lambda_pca = 1,
                                  r = 2, overlapping_amount = 0,
                                  normalize_diagonal = TRUE,
                                  n_new = 5000) {
  ###
  # This generates a dataset (X,Y) from a multivariate normal where Sigma_XX  and 
  # Sigma_YY have a correlation structure, and Sigma_XY = U Lambda V^T is rank r on a set of nnzeros rows, 0 elsewhere.
  # The number of rows (resp. Columns) on which Sigma_XX (resp. Sigma_YY)  and Sigma_XY
  # overlap is controlled by overlapping_amount
  # n: number of samples
  # p1: nb of features for X 
  # p2: nb of features for Y
  # nnzeros: number of non zero rows of U and V
  # r: rank of Sigma_XY
  # theta : canonical correlation (must be a list of length r)
  # r_pca: rank of the PCA (used in the creation of Sigma_XX and SigmaYY)
  # lambda_pca: also used to create Sigma_XX as Sigma_{XX} = U_X Lambda_pca U_X^T, and setting diag(Sigma_XX) to 1
  # Returns:
  # S : empirical covariance matrix X^TY
  # Sigma: underlying (population) covariance Sigma_{XY}
  # u: ground truth for u 
  # v: ground truth for v
  ###
  p_tot <- p1 + p2
  pp <- c(p1, p2)
  #print("We have:")
  #print(c(p_tot, nnzeros <= min(p1, p2), nnzeros))
  #print('--------------------------------------')
  #print('Generating data ...')
  Sigma <- diag(p1 + p2)

  #print(paste0('Number of non zeros is: ', nnzeros))
  start = 1
  s_pca  <- 1:min(p2, 20)
  s_pca2  <- 1:min(p2, 20)
  #if (start + nnzeros > q){
  #  s_pca2 <- 1:nnzeros
  #}
  Lambda_pca <- rep(lambda_pca, r_pca)
  # generate covariance matrix for X and Y
  groups <- createFolds(1:p1, k = ceiling(p1/5), list = TRUE, returnTrain = FALSE) ### creates groups of equal size
  s <- c()
  for (i in 1:nnzeros){
    s <- c(s, groups[[i]])
  }
  s_q = 1:p2
  if (r_pca >0){
    u1 = matrix(0, p1, r_pca)
    u1[s_pca, ] <- matrix(runif(n = length(s_pca) * r_pca, max = 3, min=1), 
                          nrow = length(s_pca), ncol = r_pca) * matrix(sample(c(-1,1), length(s_pca) * r_pca, replace=TRUE), nrow=length(s_pca), ncol=r_pca)
    # Normalize u1
    #print(dim(u1))
    #print(c(s_pca, r_pca))
    u1[s_pca, 1:r_pca] <- u1[s_pca,1:r_pca] %*% (pracma::sqrtm(t(u1[s_pca,1:r_pca]) %*% u1[s_pca, 1:r_pca])$Binv)
    T1 <- u1 %*% diag(Lambda_pca) %*% t(u1)
    if (normalize_diagonal){
      diag(T1) <- 1
      Sigma[1:p1, 1:p1] <- T1
    }else{
      Sigma[1:p1, 1:p1] <- Sigma[1:p1, 1:p1] + T1 
    }
  
     ### Same for Sigma_y
    u2 <- matrix(0, pp[2], r_pca)
    u2[s_pca2, 1:r_pca] <- matrix(runif( nnzeros * r_pca,max = 3, min=1), nrow=nnzeros)  * matrix(sample(c(-1,1), nnzeros*r_pca, replace=TRUE), nrow=nnzeros, ncol=r_pca)
    u2[s_pca2, 1:r_pca] <- u2[s_pca2,1:r_pca] %*% (pracma::sqrtm(t(u2[s_pca2,1:r_pca]) %*% u2[s_pca2,1:r_pca])$Binv)
    T2 <- u2 %*% diag(Lambda_pca) %*% t(u2)
    if (normalize_diagonal){
      diag(T2) <- 1
      Sigma[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)] <- T2
    }else{
      Sigma[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)] <- T2 + Sigma[(p1 + 1):(p1 + p2), (p1 + 1) : (p1 + p2)]
    }
  }
  Sigmax = Sigma[1:p1,1:p1];
  Sigmay = Sigma[(p1+1):p_tot,(p1+1):p_tot];

  ### Generate cross covariance
  Tss <- Sigma[s,s]
  u <- matrix(0, pp[1], r)
  for (i in 1:nnzeros){
    u[groups[[i]], 1:r] <- as.matrix(runif(r * length(groups[[i]]),
                                           max = 3, min=1), 
                                     nrow=length(groups[[i]]))
  }
  
  u <- u %*% (pracma::sqrtm(t(u[s, 1:r]) %*% Tss %*% u[s, 1:r])$Binv)
  Tss_v <- Sigma[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)][s_q, s_q]
  v <- matrix(0, pp[2], r)
  v[s_q, 1:r] <- as.matrix(runif( length(s_q) * r,max = 3, min=1), nrow=nnzeros)  * as.matrix(sample(c(-1,1), length(s_q) *r, replace=TRUE), nrow=nnzeros)
  v <- v %*% (pracma::sqrtm(t(v[s_q, 1:r]) %*% Tss_v %*% v[s_q, 1:r])$Binv)
  Sigma[(p1 + 1) :(p1 + p2), 1:p1] <- Sigmay %*%  v  %*% theta %*% t(u) %*% Sigmax
  Sigma[1:p1, (p1 + 1):(p1 + p2)] <- Sigmax %*%  u  %*% theta %*% t(v) %*% Sigmay


  
  #Generate Multivariate Normal Data According to Sigma
  sqrt_Sigma <- pracma::sqrtm(Sigma)$B
  Data <- matrix(rnorm(n = (n + n_new) * p_tot), ncol=p_tot) %*% sqrt_Sigma
  X <- Data[1:n, 1:p1]
  Y <- Data[1:n, (p1 + 1):(p1 + p2)]
  Xnew <- Data[(n+1):(n_new + n), 1:p1]
  Ynew <- Data[(n+1):(n_new + n), (p1 + 1):(p1 + p2)]
  #print('Data generated.')
  #print('--------------------------------------')
  Mask <- matrix(0, p_tot, p_tot)
  idx1 <- 1:pp[1]
  idx2 <- (pp[1] + 1):(pp[1] + pp[2])
  Mask[idx1, idx1] <- matrix(1, pp[1], pp[1])
  Mask[idx2, idx2] <- matrix(1, pp[2], pp[2])
  Sigma0 <- Sigma * Mask
  S <- cov(Data)
  sigma0hat <- S * Mask
  # Generate ground truth canonical vectors
  Sigma_X_inv <- solve(Sigma[1:p1, 1:p1])
  Sigma_Y_inv <-  solve(Sigma[(p1 + 1):(p_tot), (p1 + 1):(p_tot)])
  #GT = svd(Sigma_X_inv %*% Sigma[1:p1, (p1 + 1):p_tot] %*% Sigma_Y_inv, nu = r, nv = r)
  return(list(Sigma=Sigma, Sigma0=Sigma0,
              S = S, sigma0hat =  sigma0hat, Mask= Mask,
              X=X, Y = Y, Data=Data, u=u, v=v, 
              groups = groups,
              Xnew=Xnew,
              Ynew=Ynew,
              Sigmax=Sigmax, Sigmay=Sigmay
  ))
}



generate_example_graph <- function(n, p1 = 10,
                                   type_graph="2d-grid",
                                            p2, order = 2,
                                            r_pca = 3,
                                            nnzeros = 5,
                                            do_plot = FALSE,
                                            theta = diag(c(0.9,  0.8)),
                                            lambda_pca = 1,
                                            nnzeros_pca = 10,
                                            r = 2, overlapping_amount = 0,
                                            normalize_diagonal = TRUE,
                                   gen.using.gamma = FALSE,
                                            n_new = 50000) {
  ###
  # This generates a dataset (X,Y) from a multivariate normal where Sigma_XX  and 
  # Sigma_YY have a correlation structure, and Sigma_XY = U Lambda V^T is rank r on a set of nnzeros rows, 0 elsewhere.
  # The number of rows (resp. Columns) on which Sigma_XX (resp. Sigma_YY)  and Sigma_XY
  # overlap is controlled by overlapping_amount
  # n: number of samples
  # p1: nb of features for X 
  # p2: nb of features for Y
  # nnzeros: number of non zero rows of U and V
  # r: rank of Sigma_XY
  # theta : canonical correlation (must be a list of length r)
  # r_pca: rank of the PCA (used in the creation of Sigma_XX and SigmaYY)
  # lambda_pca: also used to create Sigma_XX as Sigma_{XX} = U_X Lambda_pca U_X^T, and setting diag(Sigma_XX) to 1
  # Returns:
  # S : empirical covariance matrix X^TY
  # Sigma: underlying (population) covariance Sigma_{XY}
  # u: ground truth for u
  # v: ground truth for v
  ###

  if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("Package 'igraph' must be installed to use generate_example_graph().", call. = FALSE)
  }

  if (type_graph == "2d-grid"){
    g <- igraph::make_lattice(c(p1, p1))
  }else{
    g <- igraph::sample_pa(p1, power=1.2)
  }
  prop_g <- get_edge_incidence(g, weight = 1)
  p1<-length(igraph::V(g))
  p_tot <- p1 + p2
  
  pp <- c(p1, p2)
  #print("We have:")
  #print(c(p_tot, nnzeros <= min(p1, p2), nnzeros))
  #print('--------------------------------------')
  #print('Generating data ...')
  Sigma <- diag(p1 + p2)
  s <- sample(1:p1, nnzeros)
  
  #print(paste0('Number of non zeros is: ', nnzeros))
  start <- 1
  nnzeros_pca = 30
  s_pca  <- (start + 1) : (start + nnzeros_pca)
  s_pca2  <- 1:p2
  Lambda_pca <- rep(lambda_pca, r_pca)
  # generate vcovariance matrix for X and Y
  
  if (r_pca > 0){
    u1 = matrix(0, p1, r_pca)
    u1[s_pca, ] <- matrix(runif(n = nnzeros_pca * r_pca, max = 3, min=0), 
                          nrow = nnzeros_pca, ncol = r_pca)
    # Normalize u1
    u1[s_pca,] <- 0.5 * diag(1/ sqrt(apply(u1[s_pca, ]^2, 1, sum))) %*% u1[s_pca,]
    u1[s_pca, 1:r_pca] <- u1[s_pca,1:r_pca] #%*% (pracma::sqrtm(t(u1[s_pca,1:r_pca]) %*% u1[s_pca, 1:r_pca])$Binv)
    T1 <- u1 %*% diag(Lambda_pca) %*% t(u1)
    if (normalize_diagonal){
      diag(T1) <- 1
      Sigma[1:p1, 1:p1] <- T1
    }else{
      Sigma[1:p1, 1:p1] <- Sigma[1:p1, 1:p1] + T1 
    } 
    ### Same for Sigma_y
    u2 <- matrix(runif( p2 * r_pca, max = 1, min=0), nrow=p2)
    u2 <- 0.5 * diag(1/ sqrt(apply(u2^2, 1, sum))) %*% u2
    #u2[s_pca2, 1:r_pca] <- u2[s_pca2, 1:r_pca] #%*% (pracma::sqrtm(t(u2[s_pca2, 1:r_pca]) %*% u2[s_pca2, 1:r_pca])$Binv)
    T2 <- u2 %*% t(u2)
    if (normalize_diagonal){
      diag(T2) <- 1
      Sigma[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)] <- T2
    }else{
      Sigma[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)] <- T2 + Sigma[(p1 + 1):(p1 + p2), (p1 + 1) : (p1 + p2)]
    }
  }
  Sigmax <- Sigma[1:p1,1:p1];
  Sigmay <- Sigma[(p1+1):p_tot,(p1+1):p_tot];
  
  ### Generate cross covariance

  Tss <- Sigmax
  #svd_gamma = svd(prop_g)
  #Gamma_dagger = svd_gamma$v[,1:(length(svd_gamma$d) -1)] %*% diag(1/svd_gamma$d[1:(length(svd_gamma$d) -1)]) %*% t(svd_gamma$u[,1:(length(svd_gamma$d) -1)] )
  
  if (gen.using.gamma){
    prod <- matrix(0, nrow(prop_g), r)
    Gamma_dagger = pinv(prop_g)
    s <- sample(1:nrow(prop_g), nnzeros)
    prod[s, 1:r] <- as.matrix(runif( length(s) * r,max = 3, min=1), nrow=length(s))  * as.matrix(sample(c(-1,1), length(s) * r, replace=TRUE), nrow=length(s))
    prod =  Gamma_dagger %*% prod
  }else{
    prod <- matrix(0, p1, r)
    s <- sample(1:p1, nnzeros)
    for (ss in s){
      ### sample neighborhood
      ind <- as.numeric(igraph::neighborhood(g, order=order, nodes = ss)[[1]])
      for(rr in 1:r){
        prod[ind, rr] <- runif(1, max = 3, min=1)  *  sample(c(-1,1),1)
      }
      
    }
  }
  
  ####
  u <- prod %*% (pracma::sqrtm(t(prod[, 1:r]) %*% Tss %*% prod[, 1:r])$Binv)
  
  

  Tss_v <- Sigmay
  prod_v <- matrix(0, pp[2], r)
  prod_v[, 1:r] <- as.matrix(runif( p2 * r , max = 3, min=1), nrow=p2)  * as.matrix(sample(c(-1,1), p2 * r, replace=TRUE), nrow=p2)
  v <- prod_v %*% (pracma::sqrtm(t(prod_v[, 1:r]) %*% Tss_v %*% prod_v[, 1:r])$Binv)
  
  #Sigma[(p1 + 1) :(p1 + p2), 1:p1] <- Sigmax %*% u %*% theta %*% t(v) %*% Sigmay
  Sigma[1:p1, (p1 + 1):(p1 + p2)] <- Sigmax %*% u %*% theta %*% t(v) %*% Sigmay
  Sigma[(p1 + 1) :(p1 + p2), 1:p1]  <- t(Sigma[1:p1, (p1 + 1):(p1 + p2)] )
  
  #Generate Multivariate Normal Data According to Sigma
  sqrt_Sigma <- pracma::sqrtm(Sigma)$B
  Data <- matrix(rnorm(n = (n + n_new) * p_tot), ncol=p_tot) %*% sqrt_Sigma
  X <- Data[1:n, 1:p1]
  Y <- Data[1:n, (p1 + 1):(p1 + p2)]
  Xnew <- Data[(n+1):(n_new + n), 1:p1]
  Ynew <- Data[(n+1):(n_new + n), (p1 + 1):(p1 + p2)]
  #print('Data generated.')
  #print('--------------------------------------')
  Mask <- matrix(0, p_tot, p_tot)
  idx1 <- 1:pp[1]
  idx2 <- (pp[1] + 1):(pp[1] + pp[2])
  Mask[idx1, idx1] <- matrix(1, pp[1], pp[1])
  Mask[idx2, idx2] <- matrix(1, pp[2], pp[2])
  Sigma0 <- Sigma * Mask
  S <- cov(Data)
  sigma0hat <- S * Mask
  # Generate ground truth canonical vectors
  #Sigma_X_inv <- pracma::sqrtm(Sigma[1:p1, 1:p1])$Binv
  #Sigma_Y_inv <-  pracma::sqrtm(Sigma[(p1+1):(p_tot), (p1+1):(p_tot)])$Binv
  #GT = svd(Sigma_X_inv %*% Sigma[1:p1, (p1 + 1):p_tot] %*% Sigma_Y_inv, nu = r, nv = r)
  return(list(Sigma=Sigma, Sigma0=Sigma0,
              S = S, sigma0hat =  sigma0hat, Mask= Mask,
              X=X, Y = Y, Data=Data, 
              u=u,#pracma::sqrtm(Sigma[1:p1, 1:p1])$Binv %*% GT$u, 
              v=v, #pracma::sqrtm(Sigma[(p1+1):(p1+p2), (p1+1):(p1+p2)])$Binv %*% GT$v, 
              Xnew = Xnew, Ynew=Ynew,
              Sigmax=Sigmax, Sigmay=Sigmay,
              Gamma=prop_g
  ))
  
}



generate_example_none_trivial_pca <- function(n, p1, p2,
                                              r_pca = 3,
                                              nnzeros = 5,
                                              theta = diag(c(0.9,  0.8)),
                                              lambda_pca = 1,
                                              r = 2, overlapping_amount = 0,
                                              normalize_diagonal = TRUE) {
  ###
  # This generates a dataset (X,Y) from a multivariate normal where Sigma_XX  and 
  # Sigma_YY have a correlation structure, and Sigma_XY = U Lambda V^T is rank r on a set of nnzeros rows, 0 elsewhere.
  # The number of rows (resp. Columns) on which Sigma_XX (resp. Sigma_YY)  and Sigma_XY
  # overlap is controlled by overlapping_amount
  # n: number of samples
  # p1: nb of features for X 
  # p2: nb of features for Y
  # nnzeros: number of non zero rows of U and V
  # r: rank of Sigma_XY
  # theta : canonical correlation (must be a list of length r)
  # r_pca: rank of the PCA (used in the creation of Sigma_XX and SigmaYY)
  # lambda_pca: also used to create Sigma_XX as Sigma_{XX} = U_X Lambda_pca U_X^T, and setting diag(Sigma_XX) to 1
  # Returns:
  # S : empirical covariance matrix X^TY
  # Sigma: underlying (population) covariance Sigma_{XY}
  # u: ground truth for u 
  # v: ground truth for v
  ###
  p_tot <- p1 + p2
  pp <- c(p1, p2)
  print("We have:")
  print(c(p_tot, nnzeros <= min(p1, p2), nnzeros))
  print('--------------------------------------')
  print('Generating data ...')
  Sigma <- diag(p1 + p2)
  s <- (1):(nnzeros)
  print(paste0('Number of non zeros is: ', nnzeros))
  start <- ceiling((1 - overlapping_amount) * nnzeros)
  s_pca  <- (start + 1) : (start + nnzeros)
  s_pca2  <- s_pca
  Lambda_pca <- rep(lambda_pca, r_pca)
  # generate vcovariance matrix for X and Y
  u1 = matrix(0, p1, r_pca)
  u1[s_pca, ] <- matrix(runif(n = nnzeros * r_pca, max = 3, min=1), nrow = nnzeros, ncol = r_pca) * matrix(sample(c(-1,1), nnzeros * r_pca, replace=TRUE), nrow=nnzeros, ncol=r_pca)
  # Normalize u1
  u1[s_pca, 1:r_pca] <- u1[s_pca,1:r_pca] %*% (sqrtm(t(u1[s_pca,1:r_pca]) %*% u1[s_pca, 1:r_pca])$Binv)
  T1 <- u1 %*% diag(Lambda_pca) %*% t(u1)
  if (normalize_diagonal){
    diag(T1) <- 1
    Sigma[1:p1, 1:p1] <- T1
  }else{
    Sigma[1:p1, 1:p1] <- Sigma[1:p1, 1:p1] + T1 
  }

   ### Same for Sigma_y
  u2 <- matrix(0, pp[2], r_pca)
  u2[s_pca2, 1:r_pca] <- matrix(runif( nnzeros * r_pca,max = 3, min=1), nrow=nnzeros)  * matrix(sample(c(-1,1), nnzeros*r_pca, replace=TRUE), nrow=nnzeros, ncol=r_pca)
  u2[s_pca2, 1:r_pca] <- u2[s_pca2,1:r_pca] %*% (sqrtm(t(u2[s_pca2,1:r_pca]) %*% u2[s_pca2,1:r_pca])$Binv)
  T2 <- u2 %*% diag(Lambda_pca) %*% t(u2)
  if (normalize_diagonal){
    diag(T2) <- 1
    Sigma[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)] <- T2
  }else{
    Sigma[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)] <- T2 + Sigma[(p1 + 1):(p1 + p2), (p1 + 1) : (p1 + p2)]
  }
  Sigmax = Sigma[1:p1,1:p1];
  Sigmay = Sigma[(p1+1):p_tot,(p1+1):p_tot];

  ### Generate cross covariance
  Tss <- Sigma[s,s]
  u <- matrix(0, pp[1], r)
  u[s, 1:r] <- as.matrix(runif( nnzeros * r,max = 3, min=1), nrow=nnzeros)  * as.matrix(sample(c(-1,1), nnzeros*r, replace=TRUE), nrow=nnzeros)
  u <- u %*% (sqrtm(t(u[s, 1:r]) %*% Tss %*% u[s, 1:r])$Binv)
  Tss_v <- Sigma[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)][s, s]
  v <- matrix(0, pp[2], r)
  v[s, 1:r] <- as.matrix(runif( nnzeros * r,max = 3, min=1), nrow=nnzeros)  * as.matrix(sample(c(-1,1), nnzeros*r, replace=TRUE), nrow=nnzeros)
  v <- v %*% (sqrtm(t(v[s, 1:r]) %*% Tss_v %*% v[s, 1:r])$Binv)
  Sigma[(p1 + 1) :(p1 + p2), 1:p1] <- Sigmay %*%  v  %*% theta %*% t(u) %*% Sigmax
  Sigma[1:p1, (p1 + 1):(p1 + p2)] <- Sigmax %*%  u  %*% theta %*% t(v) %*% Sigmay


  
  #Generate Multivariate Normal Data According to Sigma
  sqrt_Sigma <- sqrtm(Sigma)$B
  Data <- matrix(rnorm(n = 2 * n * p_tot), ncol=p_tot) %*% sqrt_Sigma
  X <- Data[1:n, 1:p1]
  Y <- Data[1:n, (p1 + 1):(p1 + p2)]
  Xnew <- Data[(n+1):(2*n), 1:p1]
  Ynew <- Data[(n+1):(2*n), (p1 + 1):(p1 + p2)]
  print('Data generated.')
  print('--------------------------------------')
  Mask <- matrix(0, p_tot, p_tot)
  idx1 <- 1:pp[1]
  idx2 <- (pp[1] + 1):(pp[1] + pp[2])
  Mask[idx1, idx1] <- matrix(1, pp[1], pp[1])
  Mask[idx2, idx2] <- matrix(1, pp[2], pp[2])
  Sigma0 <- Sigma * Mask
  S <- cov(Data[1:n, ])
  sigma0hat <- S * Mask
  # Generate ground truth canonical vectors
  Sigma_X_inv <- solve(Sigma[1:p1, 1:p1])
  Sigma_Y_inv <-  solve(Sigma[(p1 + 1):(p_tot), (p1 + 1):(p_tot)])
  GT = svd(Sigma_X_inv %*% Sigma[1:p1, (p1 + 1):p_tot] %*% Sigma_Y_inv, nu = r, nv = r)
  return(list(Sigma=Sigma, Sigma0=Sigma0,
              S = S, sigma0hat =  sigma0hat, Mask= Mask,
              X=X, Y = Y, Data=Data[1:n, ], u=GT$u, v=GT$v,     Xnew = Xnew, Ynew=Ynew,
              Sigmax=Sigmax, Sigmay=Sigmay
  ))
}

