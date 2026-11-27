#' @title Set the random-number stream of a randomisation run
#' @description
#' Set the random-number generator of the current R session to a stream
#' created by [util_make_rng_streams()], so that the following random
#' draws of a randomisation run come from that stream.
#' @param rng_stream
#' `integer` vector of length 7, a value of `.Random.seed` for the
#' `"L'Ecuyer-CMRG"` generator. Any other value gives an error.
#' @return
#' `NULL`, invisibly. Called for its side effect on `.Random.seed` in the
#' global environment, which also switches the generator kind to
#' `"L'Ecuyer-CMRG"`. Callers that must preserve the user's generator
#' save and restore it themselves (see [estimate_roc()]).
#' @seealso [util_make_rng_streams()], [run_iteration()]
#' @keywords internal
util_set_rng_stream <- function(rng_stream) {
  util_check_class(rng_stream, "integer")

  assertthat::assert_that(
    base::length(rng_stream) == 7,
    msg = "'rng_stream' must be an integer vector of length 7"
  )

  base::assign(
    x = ".Random.seed",
    value = rng_stream,
    envir = base::globalenv()
  )

  return(base::invisible(NULL))
}
