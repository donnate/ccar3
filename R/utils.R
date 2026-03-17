
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
#' @param verbose If TRUE, prints messages about the setup process.
#' @return A cluster object `cl` on success, or `NULL` on failure.

cleanup_parallel_backend <- function(cl = NULL, verbose = FALSE) {
  if (!is.null(cl)) {
    try(parallel::stopCluster(cl), silent = !verbose)
  }

  if (requireNamespace("foreach", quietly = TRUE)) {
    foreach::registerDoSEQ()
  }

  if (requireNamespace("doParallel", quietly = TRUE) &&
      exists("stopImplicitCluster", envir = asNamespace("doParallel"), inherits = FALSE)) {
    try(doParallel::stopImplicitCluster(), silent = !verbose)
  }

  invisible(NULL)
}

initialize_parallel_workers <- function(cl,
                                        pkg = "ccar3",
                                        pkg_path = NULL,
                                        verbose = FALSE) {
  if (is.null(cl)) {
    return(invisible(NULL))
  }

  parallel::clusterEvalQ(cl, {
    Sys.setenv(
      OMP_NUM_THREADS = 1,
      OPENBLAS_NUM_THREADS = 1,
      MKL_NUM_THREADS = 1,
      VECLIB_MAXIMUM_THREADS = 1
    )
    NULL
  })

  if (is.null(pkg_path)) {
    wd <- tryCatch(normalizePath(getwd(), mustWork = FALSE), error = function(e) "")
    if (nzchar(wd) &&
        file.exists(file.path(wd, "DESCRIPTION")) &&
        file.exists(file.path(wd, "NAMESPACE"))) {
      pkg_path <- wd
    }
  }

  if (!is.null(pkg_path) &&
      nzchar(pkg_path) &&
      requireNamespace("pkgload", quietly = TRUE)) {
    if (verbose) {
      cat("Initializing workers from source package at", pkg_path, "\n")
    }
    parallel::clusterCall(cl, function(path) {
      pkgload::load_all(
        path,
        quiet = TRUE,
        export_all = FALSE,
        helpers = FALSE,
        attach_testthat = FALSE
      )
      NULL
    }, pkg_path)
  } else {
    if (verbose) {
      cat("Initializing workers with installed package", pkg, "\n")
    }
    parallel::clusterCall(cl, function(pkg_name) {
      suppressPackageStartupMessages(library(pkg_name, character.only = TRUE))
      NULL
    }, pkg)
  }

  invisible(NULL)
}

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
  sysname <- Sys.info()[["sysname"]]
  prefer_fork <- .Platform$OS.type == "unix" && !identical(sysname, "Darwin")

  if (prefer_fork) {
    if (verbose) {
      cat("Unix-like system detected. Trying FORK backend first...\n")
    }
    cl <- tryCatch(
      parallel::makeCluster(num_cores, type = "FORK"),
      error = function(e_fork) {
        if (requireNamespace("crayon", quietly = TRUE)) {
          cat(crayon::yellow("Initial FORK backend setup failed: ", e_fork$message, "\n"))
        } else {
          cat("Initial FORK backend setup failed:", e_fork$message, "\n")
        }
        NULL
      }
    )
  } else if (verbose) {
    if (identical(sysname, "Darwin")) {
      cat("macOS detected. Using PSOCK backend to avoid orphaned FORK workers on interrupts...\n")
    } else {
      cat("Using PSOCK backend...\n")
    }
  }

  if (is.null(cl)) {
    if (verbose && prefer_fork) {
      cat("Attempting fallback PSOCK backend...\n")
    }
    cl <- tryCatch(
      parallel::makeCluster(num_cores, type = "PSOCK"),
      error = function(e_psock) {
        if (requireNamespace("crayon", quietly = TRUE)) {
          cat(crayon::red("PSOCK setup also failed: ", e_psock$message, "\n"))
        } else {
          cat("PSOCK setup also failed:", e_psock$message, "\n")
        }
        NULL
      }
    )
  }

  return(cl)
}
