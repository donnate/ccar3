# LIBRARIES (Ensure these are loaded)
library(MASS) # For ginv()
library(caret) # For createFolds()

Waaijenborg <- function(X, Y,
                                   lambdaxseq = seq(from = 1, to = 0, by = -0.01),
                                   lambdayseq = seq(from = 1, to = 0, by = -0.01),
                                   rank,
                                   selection.criterion = 2, # Defaulting to criterion 2, which is more common
                                   n.cv = 5,
                                   max.iter = 20,
                                   conv = 1e-3,
                                   standardize = TRUE) {

  # ===================================================================
  # 1. INTERNAL HELPER FUNCTIONS
  # ===================================================================

  # Univariate Soft Thresholding (UST)
  UST <- function(a, U) {
    # This is equivalent to sign(a) * pmax(0, abs(a) - U/2)
    val <- (abs(a) - U / 2 + abs(abs(a) - U / 2)) / 2 * sign(a)
    return(matrix(val, ncol = 1))
  }

  # Normalization to unit length
  NORMALIZATION_UNIT <- function(vec) {
    norm_val <- sqrt(sum(vec^2))
    if (norm_val == 0) return(vec)
    return(vec / norm_val)
  }

  # Robust correlation calculation
  calculate_correlation <- function(v1, v2) {
    # Check for zero variance, which would lead to NA
    if (sd(v1) < 1e-6 || sd(v2) < 1e-6) return(0)
    
    corr_val <- cor(v1, v2, use = "pairwise.complete.obs")
    
    # If correlation is NA for any other reason, treat as 0
    return(ifelse(is.na(corr_val), 0, abs(corr_val)))
  }
  
  # Cross-validation criterion 1: Minimize difference between train/test correlation
  delta_correlation <- function(U, Xtest, Ytest, Xtrain, Ytrain) {
    if (all(U == 0)) return(0)
    train_cor <- calculate_correlation(Xtrain %*% U, Ytrain)
    test_cor <- calculate_correlation(Xtest %*% U, Ytest)
    return(abs(train_cor - test_cor))
  }

  # Cross-validation criterion 2: Maximize test sample correlation
  testsample_correlation <- function(U, Xdata, yscore) {
    if (all(U == 0)) return(0)
    return(calculate_correlation(Xdata %*% U, yscore))
  }
  
  # ---
  # The new function to find the optimal lambda for one canonical vector
  # ---
  find_optimal_lambda <- function(data_matrix, score_vector, lambda_seq, n_cv, criterion_func) {
    
    # First, get candidate vectors for all lambdas
    candidate_vectors <- apply(lambda_seq, 2, UST, a = t(score_vector) %*% data_matrix)
    
    # Filter out lambdas that result in all-zero vectors
    non_zero_indices <- which(colSums(abs(candidate_vectors)) > 1e-6)
    
    # === ROBUSTNESS CHECK 1: Handle no non-zero candidates ===
    if (length(non_zero_indices) == 0) {
      # If all lambdas produce a zero vector, there's no choice to make.
      # The "best" is the zero vector, achieved with the highest lambda.
      return(list(
        lambda_opt = max(lambda_seq),
        vector_opt = matrix(0, nrow = nrow(candidate_vectors), ncol = 1)
      ))
    }
    
    candidate_vectors_reduced <- candidate_vectors[, non_zero_indices, drop = FALSE]
    lambda_seq_reduced <- matrix(lambda_seq[non_zero_indices], nrow = 1)
    
    # --- Perform k-fold cross-validation ---
    folds <- caret::createFolds(1:nrow(data_matrix), k = n_cv, list = TRUE, returnTrain = TRUE)
    
    cv_scores <- matrix(NA, ncol = length(lambda_seq_reduced), nrow = n_cv)

    for (i in 1:n_cv) {
      train_indices <- folds[[i]]
      test_indices <- setdiff(1:nrow(data_matrix), train_indices)
      
      data_train <- data_matrix[train_indices, ]
      score_train <- score_vector[train_indices]
      
      data_test <- data_matrix[test_indices, ]
      score_test <- score_vector[test_indices]
      
      # Get candidate vectors based on THIS training fold
      fold_candidate_vectors <- apply(lambda_seq_reduced, 2, UST, a = t(score_train) %*% data_train)
      
      if (criterion_func == "delta_correlation") {
        cv_scores[i, ] <- apply(fold_candidate_vectors, 2, delta_correlation, Xtest=data_test, Ytest=score_test, Xtrain=data_train, Ytrain=score_train)
      } else { # testsample_correlation
        cv_scores[i, ] <- apply(fold_candidate_vectors, 2, testsample_correlation, Xdata=data_test, yscore=score_test)
      }
    }
    
    mean_cv_scores <- colMeans(cv_scores, na.rm = TRUE)
    
    # === ROBUSTNESS CHECK 2: Handle all-NA scores ===
    if (all(is.na(mean_cv_scores))) {
        # This can happen if all folds failed. Fallback to the most sparse solution.
        best_lambda_index <- length(lambda_seq_reduced)
    } else {
        # Select lambda that MAXIMIZES the criterion score.
        # min(which.max(...)) is a robust way to handle ties.
        best_lambda_index <- min(which(mean_cv_scores == max(mean_cv_scores, na.rm = TRUE)))
    }
    
    # Get the final optimal vector and lambda
    lambda_opt <- lambda_seq_reduced[best_lambda_index]
    # The final vector should be based on ALL the data, using the optimal lambda
    vector_opt <- UST(a = t(score_vector) %*% data_matrix, U = lambda_opt)
    
    return(list(lambda_opt = lambda_opt, vector_opt = vector_opt))
  }

  # ===================================================================
  # 2. INITIALIZATION
  # ===================================================================

  if (standardize) {
    X <- scale(X, center = TRUE, scale = TRUE)
    Y <- scale(Y, center = TRUE, scale = TRUE)
  }

  p <- ncol(Y)
  q <- ncol(X)
  
  # Store results
  u_ALL <- matrix(NA, ncol = rank, nrow = p)
  v_ALL <- matrix(NA, ncol = rank, nrow = q)
  ksi_ALL <- matrix(NA, ncol = rank, nrow = nrow(X))
  omega_ALL <- matrix(NA, ncol = rank, nrow = nrow(Y))
  cancors <- matrix(NA, ncol = rank, nrow = 1)
  lambdax_ALL <- numeric(rank)
  lambday_ALL <- numeric(rank)

  # Prepare data for the main loop
  X_data <- X
  Y_data <- Y
  
  criterion_func_name <- if (selection.criterion == 1) "delta_correlation" else "testsample_correlation"

  # ===================================================================
  # 3. MAIN LOOP: SEQUENTIALLY EXTRACT CANONICAL COMPONENTS
  # ===================================================================

  for (i.r in 1:rank) {
    
    cat(paste("\n--- Finding Component", i.r, "---\n"))
    
    # --- STEP 1: FIND OPTIMAL LAMBDAS FOR THIS COMPONENT ---
    # We do this once per component using a few iterations to stabilize the choice
    
    # Initial guess for the component
    u.temp <- matrix(1 / sqrt(p), nrow = p)
    v.temp <- matrix(1 / sqrt(q), nrow = q)
    
    for (k in 1:3) { # A few stabilization iterations
      ksi.temp <- X_data %*% v.temp
      lambda_y_res <- find_optimal_lambda(Y_data, ksi.temp, lambdayseq, n.cv, criterion_func_name)
      u.temp <- NORMALIZATION_UNIT(lambda_y_res$vector_opt)
      
      omega.temp <- Y_data %*% u.temp
      lambda_x_res <- find_optimal_lambda(X_data, omega.temp, lambdaxseq, n.cv, criterion_func_name)
      v.temp <- NORMALIZATION_UNIT(lambda_x_res$vector_opt)
    }
    
    lambda_y_opt <- lambda_y_res$lambda_opt
    lambda_x_opt <- lambda_x_res$lambda_opt
    
    lambday_ALL[i.r] <- lambda_y_opt
    lambdax_ALL[i.r] <- lambda_x_opt
    
    cat(paste("  Optimal lambday:", round(lambda_y_opt, 4), "| Optimal lambdax:", round(lambda_x_opt, 4), "\n"))

    # --- STEP 2: FIT THE COMPONENT USING THE CHOSEN LAMBDAS ---
    
    u.initial <- u.temp
    v.initial <- v.temp
    
    it <- 1
    diff.u <- conv * 10
    diff.v <- conv * 10
    
    # The main iterative loop, now much faster
    while (it < max.iter && (diff.u > conv || diff.v > conv)) {
      
      ksi <- X_data %*% v.initial
      U.hat.final <- UST(a = t(ksi) %*% Y_data, U = lambda_y_opt)
      U.hat.final <- NORMALIZATION_UNIT(U.hat.final)
      
      omega <- Y_data %*% u.initial
      V.hat.final <- UST(a = t(omega) %*% X_data, U = lambda_x_opt)
      V.hat.final <- NORMALIZATION_UNIT(V.hat.final)
      
      diff.u <- max(abs(u.initial - U.hat.final))
      diff.v <- max(abs(v.initial - V.hat.final))
      
      u.initial <- U.hat.final
      v.initial <- V.hat.final
      
      it <- it + 1
    }
    cat(paste("  Converged in", it, "iterations.\n"))
    
    # --- STEP 3: STORE RESULTS AND DEFLATE MATRICES ---
    
    ksi_final <- X_data %*% V.hat.final
    omega_final <- Y_data %*% U.hat.final
    
    u_ALL[, i.r] <- U.hat.final
    v_ALL[, i.r] <- V.hat.final
    ksi_ALL[, i.r] <- ksi_final
    omega_ALL[, i.r] <- omega_final
    cancors[, i.r] <- calculate_correlation(ksi_final, omega_final)
    
    # Deflate data matrices for the next component
    gammahat <- MASS::ginv(t(ksi_final) %*% ksi_final) %*% t(ksi_final) %*% X_data
    thetahat <- MASS::ginv(t(omega_final) %*% omega_final) %*% t(omega_final) %*% Y_data
    
    X_data <- X_data - ksi_final %*% gammahat
    Y_data <- Y_data - omega_final %*% thetahat
  }
  
  # ===================================================================
  # 4. FINAL OUTPUT
  # ===================================================================
  
  colnames(u_ALL) <- colnames(v_ALL) <- paste0("Component_", 1:rank)
  
  out <- list(
    U = u_ALL,
    V = v_ALL,
    ksi = ksi_ALL,
    omega = omega_ALL,
    cor = as.numeric(cancors),
    lambdax_optimal = lambdax_ALL,
    lambday_optimal = lambday_ALL
  )
  
  return(out)
}