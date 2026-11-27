#' @title Select the number of cores for parallel computation
#' @description
#' Translate the `use_parallel` argument of [estimate_roc()] into the
#' number of cores used for the cluster.
#' @param use_parallel
#' `logical` or `numeric`, already validated by [estimate_roc()].
#' \itemize{
#' \item `FALSE` - no parallel computation.
#' \item `TRUE` - number of cores detected by [parallel::detectCores()].
#' \item `1` - no parallel computation, same as `FALSE`.
#' \item A whole `numeric` above 1 - use that many cores.
#' }
#' @return
#' `NULL` when no cluster should be created, otherwise a single `integer`
#' with the number of cores.
#' @seealso [estimate_roc()]
#' @keywords internal
util_select_n_cores <- function(use_parallel) {
  if (
    isTRUE(use_parallel)
  ) {
    res <-
      base::as.integer(parallel::detectCores())

    return(res)
  }

  # isTRUE() is FALSE for any number, so the numeric case needs its own
  #    branch. One core runs without a cluster, so that its results equal
  #    those of use_parallel = FALSE for the same seed.
  if (
    is.numeric(use_parallel) && use_parallel > 1
  ) {
    res <-
      base::as.integer(use_parallel)

    return(res)
  }

  return(NULL)
}
