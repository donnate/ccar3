
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


  #' Set up a parallel backend with graceful fallbacks.
#'
#' Attempts to create a parallel cluster, first trying the efficient FORK
#' method (on Unix-like systems), then falling back to PSOCK, and finally
#' returning NULL if all attempts fail.
#'
#' @param num_cores The number of cores to use. If NULL, it's determined automatically.
#' @return A cluster object `cl` on success, or `NULL` on failure.

setup_parallel_backend <- function(num_cores = NULL, verbose = FALSE) {

  # 1. Determine the number of cores to use
  if (is.null(num_cores)) {
    # For SLURM, read the allocated CPUs. Fallback for local use.
    n_cores_str <- Sys.getenv("SLURM_CPUS_PER_TASK")
    num_cores <- if (n_cores_str == "") (parallel::detectCores() - 1) else as.integer(n_cores_str)
  }
  # Ensure we have at least 1 core
  num_cores <- max(1, num_cores)
  if (verbose){
     cat(paste("\nAttempting to set up parallel backend with", num_cores, "cores.\n"))
  }
  
  
  cl <- NULL
  
  # Use a nested tryCatch to attempt different cluster types
  tryCatch({
    # The preferred method for Unix/macOS/Linux is FORK
    if (.Platform$OS.type == "unix") {
      if (verbose){
        cat("Unix-like system detected. Trying FORK backend (most efficient)...\n")
      }
      cl <- parallel::makeCluster(num_cores, type = "FORK")
    } else {
      # The only method for Windows is PSOCK
      if (verbose){
          cat("Windows system detected. Trying PSOCK backend...\n")
      }
      cl <- parallel::makeCluster(num_cores, type = "PSOCK")
    }
  }, error = function(e_fork) {
    # This block runs if the first attempt (e.g., FORK) failed
 
      cat(crayon::yellow("Initial backend setup failed with error: ", e_fork$message, "\n"))
      cat("Attempting fallback PSOCK backend...\n")
    
    
    # Second attempt: Try PSOCK, which works everywhere but may be blocked
    tryCatch({
      cl <<- parallel::makeCluster(num_cores, type = "PSOCK")
    }, error = function(e_psock) {
      # This block runs if PSOCK also fails
      cat(crayon::red("PSOCK setup also failed: ", e_psock$message, "\n"))
      # cl remains NULL
    })
  })
  
  return(cl)
}

