# ================================================================ #
# 1. Input validation                                              #
# ================================================================ #

testthat::test_that(
  "prepare_data() errors on invalid data_source_prep type",
  {
    purrr::walk(
      .x = list(NULL, "my_data", 123, data.frame()),
      .f = function(bad_input) {
        testthat::expect_error(
          prepare_data(
            data_source_prep = bad_input,
            working_units = "levels",
            bin_size = 500,
            rand = 1
          ),
          "'data_source_prep' must be one of the following: 'list'"
        )
      }
    )
  }
)

# ================================================================ #
# 2. Valid output                                                  #
# ================================================================ #

testthat::test_that(
  "prepare_data() returns list for working_units = 'levels'",
  {
    result <-
      suppressWarnings(
        prepare_data(
          data_source_prep = make_extracted_data(),
          working_units = "levels",
          bin_size = 500,
          rand = 1
        )
      )

    testthat::expect_type(result, "list")
  }
)

testthat::test_that(
  "prepare_data() returns list for working_units = 'bins'",
  {
    result <-
      suppressWarnings(
        prepare_data(
          data_source_prep = make_extracted_data(),
          working_units = "bins",
          bin_size = 500,
          rand = 1
        )
      )

    testthat::expect_type(result, "list")
  }
)
