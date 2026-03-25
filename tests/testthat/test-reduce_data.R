# ================================================================ #
# 1. Input validation                                              #
# ================================================================ #

testthat::test_that(
  "reduce_data() errors on invalid data_source_reduce type",
  {
    purrr::walk(
      .x = list(NULL, "my_data", 123, data.frame()),
      .f = function(bad_input) {
        testthat::expect_error(
          reduce_data(
            data_source_reduce = bad_input
          ),
          "'data_source_reduce' must be one of the following: 'list'"
        )
      }
    )
  }
)

# ================================================================ #
# 2. Valid output                                                  #
# ================================================================ #

testthat::test_that(
  "reduce_data() returns list with expected structure",
  {
    result <-
      reduce_data(
        data_source_reduce = make_extracted_data()
      )

    testthat::expect_type(result, "list")
    testthat::expect_true(
      base::all(
        c("community", "age", "age_un") %in% base::names(result)
      )
    )
    testthat::expect_true(
      base::nrow(purrr::chuck(result, "community")) > 0
    )
  }
)


