# ================================================================ #
# 1. data_source: input validation                                 #
# ================================================================ #

testthat::test_that(
  "detect_peak_points() errors on missing data_source",
  {
    testthat::expect_error(
      detect_peak_points(
        data_source = ,
        sel_method = "trend_linear",
        sd_threshold = 2
      ),
      'argument "data_source" is missing, with no default'
    )
  }
)

testthat::test_that(
  "detect_peak_points() errors on invalid data_source type",
  {
    purrr::walk(
      .x = list(NULL, "my_data", 123, NA, list()),
      .f = function(bad_input) {
        testthat::expect_error(
          detect_peak_points(
            data_source = bad_input,
            sel_method = "trend_linear",
            sd_threshold = 2
          ),
          "'data_source' must be one of the following: 'data.frame'"
        )
      }
    )
  }
)

testthat::test_that(
  "detect_peak_points() errors when data_source lacks required columns",
  {
    testthat::expect_error(
      detect_peak_points(
        data_source = data.frame(),
        sel_method = "trend_linear",
        sd_threshold = 2
      ),
      "'data_source' must contains following columns: 'ROC'"
    )
  }
)

# ================================================================ #
# 2. sel_method: input validation                                  #
# ================================================================ #

testthat::test_that(
  "detect_peak_points() errors on invalid sel_method type",
  {
    data_source <-
      make_roc_data()

    purrr::walk(
      .x = list(NULL, 123, 0, NA),
      .f = function(bad_method) {
        testthat::expect_error(
          detect_peak_points(
            data_source = data_source,
            sel_method = bad_method,
            sd_threshold = 2
          ),
          "'sel_method' must be one of the following: 'character'"
        )
      }
    )
  }
)

testthat::test_that(
  "detect_peak_points() errors on unknown sel_method value",
  {
    data_source <-
      make_roc_data()

    testthat::expect_error(
      detect_peak_points(
        data_source = data_source,
        sel_method = "my_method",
        sd_threshold = 2
      ),
      "'sel_method' must contains one of the following values"
    )
  }
)

testthat::test_that(
  "detect_peak_points() errors on multiple sel_method values",
  {
    data_source <-
      make_roc_data()

    testthat::expect_error(
      detect_peak_points(
        data_source = data_source,
        sel_method = c("SNI", "threshold"),
        sd_threshold = 2
      ),
      "'arg' must be of length 1"
    )
  }
)

# ================================================================ #
# 3. sd_threshold: input validation                                #
# ================================================================ #

testthat::test_that(
  "detect_peak_points() errors on invalid sd_threshold",
  {
    data_source <-
      make_roc_data()

    purrr::walk(
      .x = list(NULL, "2", NA, list()),
      .f = function(bad_val) {
        testthat::expect_error(
          detect_peak_points(
            data_source = data_source,
            sel_method = "trend_linear",
            sd_threshold = bad_val
          )
        )
      }
    )
  }
)

# ================================================================ #
# 4. Valid inputs: all sel_methods return correct structure        #
# ================================================================ #

testthat::test_that(
  "detect_peak_points() returns data.frame with Peak column for all methods",
  {
    data_source <-
      make_roc_data()

    vec_methods <-
      c(
        "threshold",
        "trend_linear",
        "trend_non_linear",
        "GAM_deriv",
        "SNI"
      )

    purrr::walk(
      .x = vec_methods,
      .f = function(sel_method) {
        result <-
          detect_peak_points(
            data_source = data_source,
            sel_method = sel_method,
            sd_threshold = 2
          )

        testthat::expect_s3_class(result, "data.frame")
        testthat::expect_true("Peak" %in% base::names(result))
        testthat::expect_equal(
          base::nrow(data_source),
          base::nrow(result)
        )
      }
    )
  }
)
