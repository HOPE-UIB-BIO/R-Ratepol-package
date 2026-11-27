# ================================================================ #
# 1. Input validation                                              #
# ================================================================ #

testthat::test_that(
  "util_make_rng_streams() errors on invalid n_streams",
  {
    purrr::walk(
      .x = list("two", 0, -1, 1.5, c(1, 2)),
      .f = ~ testthat::expect_error(
        util_make_rng_streams(
          n_streams = .x,
          seed = 900723
        )
      )
    )
  }
)

testthat::test_that(
  "util_make_rng_streams() errors on invalid seed",
  {
    purrr::walk(
      .x = list("seed", NA_real_, 1.5, c(1, 2)),
      .f = ~ testthat::expect_error(
        util_make_rng_streams(
          n_streams = 2,
          seed = .x
        )
      )
    )
  }
)

# ================================================================ #
# 2. Output structure                                              #
# ================================================================ #

testthat::test_that(
  "util_make_rng_streams() returns one L'Ecuyer state per stream",
  {
    list_streams <-
      util_make_rng_streams(
        n_streams = 5,
        seed = 900723
      )

    testthat::expect_type(list_streams, "list")

    testthat::expect_length(list_streams, 5)

    purrr::walk(
      .x = list_streams,
      .f = ~ {
        testthat::expect_type(.x, "integer")

        testthat::expect_length(.x, 7)

        # The last two digits of the first element code the generator
        #    kind; 7 is "L'Ecuyer-CMRG"
        testthat::expect_identical(
          purrr::chuck(.x, 1) %% 100L,
          7L
        )
      }
    )
  }
)

# ================================================================ #
# 3. Correctness                                                   #
# ================================================================ #

testthat::test_that(
  "util_make_rng_streams() chains streams with nextRNGStream()",
  {
    list_streams <-
      util_make_rng_streams(
        n_streams = 4,
        seed = 900723
      )

    purrr::walk(
      .x = 1:3,
      .f = ~ testthat::expect_identical(
        purrr::chuck(list_streams, .x + 1),
        parallel::nextRNGStream(purrr::chuck(list_streams, .x))
      )
    )
  }
)

testthat::test_that(
  "util_make_rng_streams() is reproducible and depends on the seed",
  {
    list_streams_a <-
      util_make_rng_streams(
        n_streams = 3,
        seed = 900723
      )

    list_streams_b <-
      util_make_rng_streams(
        n_streams = 3,
        seed = 900723
      )

    list_streams_c <-
      util_make_rng_streams(
        n_streams = 3,
        seed = 123
      )

    testthat::expect_identical(list_streams_a, list_streams_b)

    testthat::expect_false(
      base::identical(list_streams_a, list_streams_c)
    )

    testthat::expect_length(
      base::unique(list_streams_a),
      3
    )
  }
)

testthat::test_that(
  "util_make_rng_streams() leaves the caller's generator unchanged",
  {
    base::set.seed(900723)

    vec_kind_before <-
      base::RNGkind()

    vec_seed_before <-
      base::get(".Random.seed", envir = base::globalenv())

    list_streams <-
      util_make_rng_streams(
        n_streams = 3,
        seed = 900723
      )

    testthat::expect_identical(base::RNGkind(), vec_kind_before)

    testthat::expect_identical(
      base::get(".Random.seed", envir = base::globalenv()),
      vec_seed_before
    )
  }
)
