# ================================================================ #
# 1. Input validation                                              #
# ================================================================ #

testthat::test_that(
  "transform_into_proportions() errors on invalid data type",
  {
    purrr::walk(
      .x = list(NULL, "my_data", 123, list()),
      .f = function(bad_input) {
        testthat::expect_error(
          transform_into_proportions(
            data_source_trans = bad_input
          ),
          "'data_source_trans' must be one of the following: 'data.frame'"
        )
      }
    )
  }
)

# ================================================================ #
# 2. Valid output: proportions                                     #
# ================================================================ #

testthat::test_that(
  "transform_into_proportions() returns proportions in [0, 1]",
  {
    data_run <-
      make_run_data()

    data_input <-
      reduce_data_simple(
        data_source_reduce =
          subset_samples(
            data_source_subset = purrr::chuck(data_run, "data"),
            data_source_bins = purrr::chuck(data_run, "bins"),
            bin_selection = "first"
          )
      )

    result <-
      transform_into_proportions(
        data_source_trans = data_input,
        sel_method = "proportions",
        silent = TRUE
      )

    testthat::expect_s3_class(result, "data.frame")

    taxa_values <-
      base::as.vector(
        base::as.matrix(
          subset_community(result)
        )
      )

    testthat::expect_true(
      base::all(taxa_values >= 0 & taxa_values <= 1)
    )
  }
)

testthat::test_that(
  "transform_into_proportions() returns percentages (values up to 100)",
  {
    data_run <-
      make_run_data()

    data_input <-
      reduce_data_simple(
        data_source_reduce = subset_samples(
          data_source_subset = purrr::chuck(data_run, "data"),
          data_source_bins = purrr::chuck(data_run, "bins"),
          bin_selection = "first"
        )
      )

    result <-
      transform_into_proportions(
        data_source_trans = data_input,
        sel_method = "percentages",
        silent = TRUE
      )

    testthat::expect_s3_class(result, "data.frame")

    taxa_max <-
      base::max(
        base::as.vector(
          base::as.matrix(
            subset_community(result)
          )
        )
      )

    testthat::expect_true(taxa_max <= 100)
  }
)

testthat::test_that(
  "transform_into_proportions() errors on empty data.frame",
  {
    testthat::expect_error(
      transform_into_proportions(
        data_source_trans = data.frame()
      ),
      "'data_source_trans' must not be empty"
    )
  }
)

testthat::test_that(
  "transform_into_proportions() errors when sel_method has length > 1",
  {
    data_run <-
      make_run_data()

    data_input <-
      reduce_data_simple(
        data_source_reduce =
          subset_samples(
            data_source_subset = purrr::chuck(data_run, "data"),
            data_source_bins = purrr::chuck(data_run, "bins"),
            bin_selection = "first"
          )
      )

    testthat::expect_error(
      transform_into_proportions(
        data_source_trans = data_input,
        sel_method = c("proportions", "percentages")
      ),
      "'sel_method' must be a single value"
    )
  }
)
