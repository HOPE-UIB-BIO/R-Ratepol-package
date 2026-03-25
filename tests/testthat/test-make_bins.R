# ================================================================ #
# 1. Input validation                                              #
# ================================================================ #

testthat::test_that(
  "make_bins() errors on invalid data_source_bins type",
  {
    purrr::walk(
      .x = list(NULL, "my_data", 123, data.frame()),
      .f = function(bad_input) {
        testthat::expect_error(
          make_bins(
            data_source_bins = bad_input,
            working_units = "levels"
          ),
          "'data_source_bins' must be one of the following: 'list'"
        )
      }
    )
  }
)

# ================================================================ #
# 2. Valid output: all working_units                               #
# ================================================================ #

testthat::test_that(
  "make_bins() returns data.frame for all working_units",
  {
    data_input <-
      make_extracted_data()

    purrr::walk(
      .x = c("levels", "bins", "MW"),
      .f = function(wu) {
        result <-
          make_bins(
            data_source_bins = data_input,
            working_units = wu,
            bin_size = 500,
            number_of_shifts = 1
          )

        testthat::expect_s3_class(result, "data.frame")
        testthat::expect_true(base::nrow(result) > 0)
      }
    )
  }
)

testthat::test_that(
  "make_bins() errors when working_units has length > 1",
  {
    data_input <-
      make_extracted_data()

    testthat::expect_error(
      make_bins(
        data_source_bins = data_input,
        working_units = c("levels", "bins")
      ),
      "'working_units' must be a single value"
    )
  }
)
