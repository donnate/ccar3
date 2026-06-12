cran_limits_cores <- function() {
  value <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = "false"))
  value %in% c("true", "1", "yes")
}

skip_if_cran_limits_cores <- function() {
  testthat::skip_if(
    cran_limits_cores(),
    "CRAN limits simultaneous parallel worker processes."
  )
}
