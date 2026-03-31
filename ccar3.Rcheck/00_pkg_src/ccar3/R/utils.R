
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


# Set up a parallel backend with graceful fallbacks.
#
# Attempts to create a parallel cluster, first trying the efficient FORK
# method (on Unix-like systems), then falling back to PSOCK, and finally
# returning NULL if all attempts fail.

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

.safe_mean_na <- function(x) {
  out <- mean(x, na.rm = TRUE)
  if (is.nan(out)) NA_real_ else out
}

.safe_sd_na <- function(x) {
  out <- stats::sd(x, na.rm = TRUE)
  if (is.nan(out)) NA_real_ else out
}

.create_cv_folds <- function(n, k) {
  n <- as.integer(n)
  k <- as.integer(k)

  if (length(n) != 1L || is.na(n) || n < 2L) {
    stop("`n` must be a single integer >= 2.", call. = FALSE)
  }
  if (length(k) != 1L || is.na(k) || k < 2L) {
    stop("`k` must be a single integer >= 2.", call. = FALSE)
  }

  k <- min(k, n)
  fold_ids <- sample(rep(seq_len(k), length.out = n))
  unname(split(seq_len(n), fold_ids))
}

.find_local_global_names <- function(fun) {
  if (!requireNamespace("codetools", quietly = TRUE)) {
    return(character())
  }

  env <- environment(fun)
  if (is.null(env)) {
    return(character())
  }

  globals <- tryCatch(codetools::findGlobals(fun, merge = FALSE), error = function(e) NULL)
  if (is.null(globals)) {
    return(character())
  }

  names <- unique(c(globals$variables, globals$functions))
  names <- names[nzchar(names) & !grepl("^[[:punct:]]+$", names)]

  binding_env <- function(name, start) {
    cur <- start
    while (!identical(cur, emptyenv())) {
      if (exists(name, envir = cur, inherits = FALSE)) {
        return(cur)
      }
      cur <- parent.env(cur)
    }
    NULL
  }

  keep <- vapply(names, function(name) {
    where <- binding_env(name, env)
    if (is.null(where) || identical(where, baseenv())) {
      return(FALSE)
    }

    if ("ccar3" %in% loadedNamespaces()) {
      ccar3_ns <- asNamespace("ccar3")
      if (exists(name, envir = ccar3_ns, inherits = FALSE)) {
        return(FALSE)
      }
    }

    env_name <- environmentName(where)
    identical(where, globalenv()) || identical(env_name, "")
  }, logical(1))

  names[keep]
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
      KMP_SHM_DISABLE = 1,
      KMP_SHM_DIR = "/tmp",
      TMPDIR = "/tmp",
      OMP_NUM_THREADS = 1,
      OPENBLAS_NUM_THREADS = 1,
      MKL_NUM_THREADS = 1,
      VECLIB_MAXIMUM_THREADS = 1
    )
    NULL
  })

  if (is.null(pkg_path)) {
    ns_path <- tryCatch({
      if (requireNamespace("pkgload", quietly = TRUE) &&
          pkg %in% loadedNamespaces() &&
          pkgload::is_dev_package(pkg)) {
        getNamespaceInfo(asNamespace(pkg), "path")
      } else {
        ""
      }
    }, error = function(e) "")

    if (nzchar(ns_path) &&
        file.exists(file.path(ns_path, "DESCRIPTION")) &&
        file.exists(file.path(ns_path, "NAMESPACE"))) {
      pkg_path <- ns_path
    }
  }

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

run_lambda_jobs <- function(lambdas,
                            run_one,
                            parallelize = TRUE,
                            nb_cores = NULL,
                            packages = character(),
                            exports = NULL,
                            verbose = FALSE) {
  if (!parallelize) {
    return(lapply(lambdas, run_one))
  }

  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop("Package 'doParallel' must be installed to use the parallelization option.",
         call. = FALSE)
  }

  if (!requireNamespace("crayon", quietly = TRUE)) {
    stop("Package 'crayon' must be installed to use the parallelization option.",
         call. = FALSE)
  }

  cl <- setup_parallel_backend(nb_cores, verbose = verbose)
  if (is.null(cl)) {
    warning("Parallel setup failed; running CV serially.", immediate. = TRUE)
    return(lapply(lambdas, run_one))
  }

  initialize_parallel_workers(cl, verbose = verbose)
  local_env <- environment(run_one)
  local_names <- if (is.null(local_env)) character() else ls(envir = local_env, all.names = TRUE)
  local_names <- setdiff(local_names, "run_one")
  if (length(local_names) > 0) {
    parallel::clusterExport(cl, varlist = local_names, envir = local_env)
  }

  global_names <- unique(c(exports, .find_local_global_names(run_one)))
  global_names <- setdiff(global_names, c("run_one", local_names))
  global_names <- global_names[vapply(global_names, exists, logical(1), envir = globalenv(), inherits = FALSE)]
  if (length(global_names) > 0) {
    parallel::clusterExport(cl, varlist = global_names, envir = globalenv())
  }

  parallel::clusterExport(cl, varlist = "run_one", envir = parent.frame())
  doParallel::registerDoParallel(cl)
  on.exit(cleanup_parallel_backend(cl), add = TRUE)

  foreach::foreach(
    lambda = lambdas,
    .packages = unique(c(packages, "foreach"))
  ) %dopar% run_one(lambda)
}
