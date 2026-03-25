# ================================================================ #
# 1. Required arguments: input validation                          #
# ================================================================ #

testthat::test_that(
  "extract_data() errors on missing data_community_extract",
  {
    testthat::expect_error(
      extract_data(
        data_community_extract = ,
        data_age_extract =
          RRatepol::example_data$sample_age[[1]],
        silent = TRUE
      ),
      'argument "data_community_extract" is missing, with no default'
    )
  }
)

testthat::test_that(
  "extract_data() errors on invalid data types",
  {
    purrr::walk(
      .x = list(NULL, "my_data", 123, NA),
      .f = function(bad_input) {
        testthat::expect_error(
          extract_data(
            data_community_extract = bad_input,
            data_age_extract =
              RRatepol::example_data$sample_age[[1]],
            silent = TRUE
          )
        )
      }
    )
  }
)

# ================================================================ #
# 2. Valid output: structure and column names                      #
# ================================================================ #

testthat::test_that(
  "extract_data() returns list with expected structure",
  {
    result <-
      suppressWarnings(
        extract_data(
          data_community_extract =
            RRatepol::example_data$pollen_data[[1]],
          data_age_extract =
            RRatepol::example_data$sample_age[[1]],
          silent = TRUE
        )
      )

    testthat::expect_type(result, "list")
    testthat::expect_true(
      base::all(
        c("community", "age", "age_un") %in% base::names(result)
      )
    )
    testthat::expect_s3_class(
      purrr::chuck(result, "community"), "data.frame"
    )
    testthat::expect_true(
      base::nrow(purrr::chuck(result, "community")) > 0
    )
  }
)
