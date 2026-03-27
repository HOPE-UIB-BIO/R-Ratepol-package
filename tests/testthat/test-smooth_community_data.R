# ================================================================ #
# 1. Input validation                                              #
# ================================================================ #

testthat::test_that(
  "smooth_community_data() errors on invalid data_source_smooth type",
  {
    purrr::walk(
      .x = list(NULL, "my_data", 123, data.frame()),
      .f = function(bad_input) {
        testthat::expect_error(
          smooth_community_data(
            data_source_smooth = bad_input,
            smooth_method = "shep"
          ),
          "'data_source_smooth' must be one of the following: 'list'"
        )
      }
    )
  }
)

# ================================================================ #
# 2. Valid output                                                  #
# ================================================================ #

testthat::test_that(
  "smooth_community_data() returns list for smooth_method = 'shep'",
  {
    result <-
      suppressWarnings(
        smooth_community_data(
          data_source_smooth = make_extracted_data(),
          smooth_method = "shep"
        )
      )

    testthat::expect_type(result, "list")
    testthat::expect_true(
      base::nrow(purrr::chuck(result, "community")) > 0
    )
  }
)

testthat::test_that(
  "smooth_community_data() returns list for smooth_method = 'm.avg'",
  {
    result <-
      suppressWarnings(
        smooth_community_data(
          data_source_smooth = make_extracted_data(),
          smooth_method = "m.avg",
          smooth_n_points = 3
        )
      )

    testthat::expect_type(result, "list")
    testthat::expect_true(
      base::nrow(purrr::chuck(result, "community")) > 0
    )
  }
)
