# Skip a test when RRatepol is not installed in a library. Parallel workers
#    load the installed package, and testthat::skip_if_not_installed() also
#    accepts the copy loaded by devtools::load_all().
skip_if_rratepol_not_installed <- function() {
  vec_path <-
    base::find.package(
      package = "RRatepol",
      lib.loc = base::.libPaths(),
      quiet = TRUE
    )

  testthat::skip_if(
    condition = base::length(vec_path) == 0,
    message = "RRatepol is not installed in a library"
  )
}
