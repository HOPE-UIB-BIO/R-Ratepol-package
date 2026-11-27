#' @title Create one random-number stream per randomisation run
#' @description
#' Create independent L'Ecuyer-CMRG random-number streams, one per
#' randomisation run of [estimate_roc()], so that each run draws the same
#' random numbers wherever it is computed (single core or any worker).
#' @param n_streams
#' `numeric`. Whole number of streams to create, at least 1.
#' @param seed
#' `numeric`. Single whole number that seeds the first stream.
#' @return
#' A `list` of length `n_streams`. Each element is an `integer` vector of
#' length 7, a value of `.Random.seed` for the `"L'Ecuyer-CMRG"`
#' generator.
#' @details
#' The first stream is the state of the `"L'Ecuyer-CMRG"` generator after
#' `set.seed(seed)`, and each following stream is
#' [parallel::nextRNGStream()] of the previous one. The caller's
#' random-number generator kind and state are left unchanged.
#' @seealso [util_set_rng_stream()], [estimate_roc()]
#' @keywords internal
util_make_rng_streams <- function(
  n_streams,
  seed
) {
  util_check_class(n_streams, c("numeric", "integer"))

  assertthat::assert_that(
    base::length(n_streams) == 1 && !base::is.na(n_streams) &&
      n_streams >= 1,
    msg = "'n_streams' must be a single number of at least 1"
  )

  util_check_if_integer(n_streams)

  util_check_class(seed, c("numeric", "integer"))

  assertthat::assert_that(
    base::length(seed) == 1 && !base::is.na(seed),
    msg = "'seed' must be a single number"
  )

  util_check_if_integer(seed)

  # set.seed() below replaces the caller's generator, so keep it and put
  #    it back on exit. Restoring .Random.seed also restores the kind.
  flag_seed_exists <-
    base::exists(
      x = ".Random.seed",
      envir = base::globalenv(),
      inherits = FALSE
    )

  vec_kind_user <-
    base::RNGkind()

  if (
    isTRUE(flag_seed_exists)
  ) {
    vec_seed_user <-
      base::get(".Random.seed", envir = base::globalenv())

    on.exit(
      base::assign(
        x = ".Random.seed",
        value = vec_seed_user,
        envir = base::globalenv()
      ),
      add = TRUE
    )
  } else {
    on.exit(
      {
        base::suppressWarnings(
          base::RNGkind(
            kind = vec_kind_user[1],
            normal.kind = vec_kind_user[2],
            sample.kind = vec_kind_user[3]
          )
        )

        base::rm(".Random.seed", envir = base::globalenv())
      },
      add = TRUE
    )
  }

  base::RNGkind(kind = "L'Ecuyer-CMRG")

  base::set.seed(seed)

  vec_stream_first <-
    base::get(".Random.seed", envir = base::globalenv())

  res <-
    purrr::accumulate(
      .x = base::seq_len(n_streams - 1),
      .f = ~ parallel::nextRNGStream(.x),
      .init = vec_stream_first
    )

  return(res)
}
