# ================================================================ #
# 1. Required arguments: input validation                          #
# ================================================================ #

testthat::test_that(
  "run_iteration() errors on missing data_source_run",
  {
    testthat::expect_error(
      run_iteration(
        data_source_run = ,
        bin_selection = "first",
        standardise = FALSE,
        tranform_to_proportions = TRUE,
        dissimilarity_coefficient = "euc",
        time_standardisation = 500,
        silent = TRUE
      ),
      'argument "data_source_run" is missing, with no default'
    )
  }
)

testthat::test_that(
  "run_iteration() errors on invalid data_source_run type",
  {
    purrr::walk(
      .x = list(NULL, "my_data", 123),
      .f = function(bad_input) {
        testthat::expect_error(
          run_iteration(
            data_source_run = bad_input,
            bin_selection = "first",
            standardise = FALSE,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = "euc",
            time_standardisation = 500,
            silent = TRUE
          )
        )
      }
    )
  }
)

# ================================================================ #
# 2. Valid output: returns data.frame with expected columns        #
# ================================================================ #

testthat::test_that(
  "run_iteration() returns data.frame with label, res_age, roc",
  {
    data_run <-
      make_run_data()

    result <-
      run_iteration(
        data_source_run = data_run,
        bin_selection = "first",
        standardise = FALSE,
        tranform_to_proportions = TRUE,
        dissimilarity_coefficient = "euc",
        time_standardisation = 500,
        silent = TRUE
      )

    testthat::expect_s3_class(result, "data.frame")
    testthat::expect_true(
      base::all(
        c("label", "res_age", "roc") %in% base::colnames(result)
      )
    )
    testthat::expect_true(base::nrow(result) > 0)
  }
)

testthat::test_that(
  "run_iteration() works with all dissimilarity coefficients",
  {
    data_run <-
      make_run_data()

    purrr::walk(
      .x = c("euc", "euc.sd", "chord", "chisq", "gower", "bray"),
      .f = function(dc) {
        result <-
          run_iteration(
            data_source_run = data_run,
            bin_selection = "first",
            standardise = FALSE,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = dc,
            time_standardisation = 500,
            silent = TRUE
          )

        testthat::expect_s3_class(result, "data.frame")
        testthat::expect_true(base::nrow(result) > 0)
      }
    )
  }
)

testthat::test_that(
  "run_iteration() errors when time_standardisation is 0",
  {
    data_run <-
      make_run_data()

    testthat::expect_error(
      run_iteration(
        data_source_run = data_run,
        bin_selection = "first",
        standardise = FALSE,
        tranform_to_proportions = TRUE,
        dissimilarity_coefficient = "euc",
        time_standardisation = 0,
        silent = TRUE
      ),
      "'time_standardisation' must not be 0 or NA"
    )
  }
)

testthat::test_that(
  "run_iteration() errors when time_standardisation is NA",
  {
    data_run <-
      make_run_data()

    testthat::expect_error(
      run_iteration(
        data_source_run = data_run,
        bin_selection = "first",
        standardise = FALSE,
        tranform_to_proportions = TRUE,
        dissimilarity_coefficient = "euc",
        time_standardisation = NA_real_,
        silent = TRUE
      ),
      "'time_standardisation' must not be 0 or NA"
    )
  }
)

testthat::test_that(
  "run_iteration() errors when standardise is not logical",
  {
    data_run <-
      make_run_data()

    testthat::expect_error(
      run_iteration(
        data_source_run = data_run,
        bin_selection = "first",
        standardise = "yes",
        tranform_to_proportions = TRUE,
        dissimilarity_coefficient = "euc",
        time_standardisation = 500,
        silent = TRUE
      )
    )
  }
)
