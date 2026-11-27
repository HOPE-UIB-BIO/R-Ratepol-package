# ================================================================ #
# 1. Output: number of cores for each use_parallel mode            #
# ================================================================ #

testthat::test_that(
  "util_select_n_cores() returns NULL for FALSE",
  {
    testthat::expect_null(
      util_select_n_cores(use_parallel = FALSE)
    )
  }
)

testthat::test_that(
  "util_select_n_cores() returns NULL for a single core",
  {
    testthat::expect_null(
      util_select_n_cores(use_parallel = 1)
    )

    testthat::expect_null(
      util_select_n_cores(use_parallel = 1L)
    )
  }
)

testthat::test_that(
  "util_select_n_cores() uses available cores for TRUE",
  {
    n_available <-
      base::as.integer(parallelly::availableCores())

    testthat::skip_if(
      n_available < 2,
      message = "Fewer than 2 cores available"
    )

    res <-
      util_select_n_cores(use_parallel = TRUE)

    testthat::expect_type(res, "integer")

    testthat::expect_length(res, 1)

    testthat::expect_identical(res, n_available)
  }
)

testthat::test_that(
  "util_select_n_cores() respects the cores allocated to the session",
  {
    testthat::local_mocked_bindings(
      availableCores = function(...) c(slurm = 3L),
      .package = "parallelly"
    )

    testthat::expect_identical(
      util_select_n_cores(use_parallel = TRUE),
      3L
    )
  }
)

testthat::test_that(
  "util_select_n_cores() returns NULL for TRUE with one available core",
  {
    testthat::local_mocked_bindings(
      availableCores = function(...) c(system = 1L),
      .package = "parallelly"
    )

    testthat::expect_null(
      util_select_n_cores(use_parallel = TRUE)
    )
  }
)

testthat::test_that(
  "util_select_n_cores() returns a numeric value as integer",
  {
    purrr::walk(
      .x = c(2, 4),
      .f = ~ testthat::expect_identical(
        util_select_n_cores(use_parallel = .x),
        base::as.integer(.x)
      )
    )
  }
)

testthat::test_that(
  "util_select_n_cores() accepts an integer value",
  {
    testthat::expect_identical(
      util_select_n_cores(use_parallel = 3L),
      3L
    )
  }
)
