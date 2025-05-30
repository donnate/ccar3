
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