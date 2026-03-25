# ================================================================ #
# 1. Input validation                                              #
# ================================================================ #

testthat::test_that(
  "subset_community() errors on invalid data_source type",
  {
    purrr::walk(
      .x = list(NULL, "my_data", 123, list()),
      .f = function(bad_input) {
        testthat::expect_error(
          subset_community(data_source = bad_input)
        )
      }
    )
  }
)

# ================================================================ #
# 2. Valid output: drops meta columns                             #
# ================================================================ #

testthat::test_that(
  "subset_community() returns data.frame without meta columns",
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
      subset_community(data_source = data_subset)

    testthat::expect_s3_class(result, "data.frame")
    testthat::expect_false(
      "label" %in% base::colnames(result)
    )
    testthat::expect_false(
      "res_age" %in% base::colnames(result)
    )
  }
)
