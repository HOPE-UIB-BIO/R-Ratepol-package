# ================================================================ #
# 1. Input validation                                              #
# ================================================================ #

testthat::test_that(
  "subset_samples() errors with NULL data_source_bins",
  {
    testthat::expect_error(
      subset_samples(
        data_source_subset = NULL,
        data_source_bins = NULL,
        bin_selection = "first"
      )
    )
  }
)

# ================================================================ #
# 2. Valid output: returns data.frame with expected columns        #
# ================================================================ #

testthat::test_that(
  "subset_samples() returns data.frame with expected columns",
  {
    data_run <-
      make_run_data()

    result <-
      subset_samples(
        data_source_subset = purrr::chuck(data_run, "data"),
        data_source_bins = purrr::chuck(data_run, "bins"),
        bin_selection = "first"
      )

    testthat::expect_s3_class(result, "data.frame")
    testthat::expect_true(
      base::all(
        c("label", "res_age") %in% base::colnames(result)
      )
    )
    testthat::expect_true(base::nrow(result) > 0)
  }
)
