# ================================================================ #
# 1. Input validation                                              #
# ================================================================ #

testthat::test_that(
  "util_set_rng_stream() errors on an invalid stream",
  {
    purrr::walk(
      .x = list(NULL, "stream", 1:6, c(1, 2, 3, 4, 5, 6, 7)),
      .f = ~ testthat::expect_error(
        util_set_rng_stream(rng_stream = .x)
      )
    )
  }
)

# ================================================================ #
# 2. Correctness                                                   #
# ================================================================ #

testthat::test_that(
  "util_set_rng_stream() makes the following draws come from the stream",
  {
    base::set.seed(900723)

    vec_seed_user <-
      base::get(".Random.seed", envir = base::globalenv())

    list_streams <-
      util_make_rng_streams(
        n_streams = 2,
        seed = 900723
      )

    util_set_rng_stream(
      rng_stream = purrr::chuck(list_streams, 1)
    )

    vec_seed_set <-
      base::get(".Random.seed", envir = base::globalenv())

    vec_draws_first <-
      stats::runif(n = 3)

    util_set_rng_stream(
      rng_stream = purrr::chuck(list_streams, 1)
    )

    vec_draws_again <-
      stats::runif(n = 3)

    util_set_rng_stream(
      rng_stream = purrr::chuck(list_streams, 2)
    )

    vec_draws_other <-
      stats::runif(n = 3)

    # Restore the generator of the test session
    base::assign(
      x = ".Random.seed",
      value = vec_seed_user,
      envir = base::globalenv()
    )

    testthat::expect_identical(
      vec_seed_set,
      purrr::chuck(list_streams, 1)
    )

    testthat::expect_identical(vec_draws_first, vec_draws_again)

    testthat::expect_false(
      base::identical(vec_draws_first, vec_draws_other)
    )
  }
)
