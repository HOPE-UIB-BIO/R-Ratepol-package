# ================================================================ #
# 1. Input validation                                              #
# ================================================================ #

testthat::test_that(
  "reduce_data_simple() errors on invalid data_source_reduce type",
  {
    purrr::walk(
      .x = list(NULL, "my_data", 123),
      .f = function(bad_input) {
        testthat::expect_error(
          reduce_data_simple(
            data_source_reduce = bad_input
          )
        )
      }
    )
  }
)

# ================================================================ #
# 2. Valid output                                                  #
# ================================================================ #

testthat::test_that(
  "reduce_data_simple() returns data.frame with community columns",
  {
    data_run <-
      make_run_data()

    data_subset <-
      subset_samples(
        data_source_subset = purrr::chuck(data_run, "data"),
        data_source_bins = purrr::chuck(data_run, "bins"),
        bin_selection = "first"
      )

    result <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )

    testthat::expect_s3_class(result, "data.frame")
    testthat::expect_true(base::nrow(result) > 0)
  }
)

testthat::test_that(
  "reduce_data_simple() errors on empty data.frame",
  {
    testthat::expect_error(
      reduce_data_simple(
        data_source_reduce = data.frame()
      ),
      "'data_source_reduce' must not be empty"
    )
  }
)

testthat::test_that(
  "reduce_data_simple() errors when check_taxa is not logical",
  {
    data_run <-
      make_run_data()

    data_subset <-
      subset_samples(
        data_source_subset = purrr::chuck(data_run, "data"),
        data_source_bins = purrr::chuck(data_run, "bins"),
        bin_selection = "first"
      )

    testthat::expect_error(
      reduce_data_simple(
        data_source_reduce = data_subset,
        check_taxa = 1
      )
    )
  }
)

testthat::test_that(
  "reduce_data_simple() errors when check_levels is not logical",
  {
    data_run <-
      make_run_data()

    data_subset <-
      subset_samples(
        data_source_subset = purrr::chuck(data_run, "data"),
        data_source_bins = purrr::chuck(data_run, "bins"),
        bin_selection = "first"
      )

    testthat::expect_error(
      reduce_data_simple(
        data_source_reduce = data_subset,
        check_levels = 1
      )
    )
  }
)
