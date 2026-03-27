# ================================================================ #
# 1. Input validation                                              #
# ================================================================ #

testthat::test_that(
  "make_trend() errors on invalid data_source type",
  {
    purrr::walk(
      .x = list(NULL, "my_data", 123, list()),
      .f = function(bad_input) {
        testthat::expect_error(
          make_trend(
            data_source = bad_input
          ),
          "'data_source' must be one of the following: 'data.frame'"
        )
      }
    )
  }
)

# ================================================================ #
# 2. Valid output: both trend methods                              #
# ================================================================ #

testthat::test_that(
  "make_trend() returns numeric vector for both methods",
  {
    data_source <-
      make_roc_data()

    purrr::walk(
      .x = c("linear", "non_linear"),
      .f = function(method) {
        result <-
          make_trend(
            data_source = data_source,
            sel_method = method
          )

        testthat::expect_type(result, "double")
        testthat::expect_equal(
          base::length(result),
          base::nrow(data_source)
        )
      }
    )
  }
)

testthat::test_that(
  "make_trend() errors when sel_method has length > 1",
  {
    data_source <-
      make_roc_data()

    testthat::expect_error(
      make_trend(
        data_source = data_source,
        sel_method = c("linear", "non_linear")
      ),
      "'sel_method' must be a single value"
    )
  }
)

testthat::test_that(
  "make_trend() returns numeric (not array) for non_linear",
  {
    data_source <-
      make_roc_data()

    result <-
      make_trend(
        data_source = data_source,
        sel_method = "non_linear"
      )

    testthat::expect_type(result, "double")
    testthat::expect_false(
      base::inherits(result, "array")
    )
  }
)
