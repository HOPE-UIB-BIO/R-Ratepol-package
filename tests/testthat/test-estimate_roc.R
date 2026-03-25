# ================================================================ #
# 1. Required arguments: input validation                          #
# ================================================================ #

testthat::test_that(
  "estimate_roc() errors on missing data_source_community",
  {
    testthat::expect_error(
      estimate_roc(
        data_source_community = ,
        data_source_age = RRatepol::example_data$sample_age[[1]],
        silent = TRUE
      ),
      "Object 'data_source_community' must be included as a 'data.frame'"
    )
  }
)

testthat::test_that(
  "estimate_roc() errors on invalid data_source_community type",
  {
    purrr::walk(
      .x = list(NULL, "data", 123, NA, list()),
      .f = function(bad_input) {
        testthat::expect_error(
          estimate_roc(
            data_source_community = bad_input,
            data_source_age =
              RRatepol::example_data$sample_age[[1]],
            silent = TRUE
          ),
          "'data_source_community' must be one of the following: 'data.frame'"
        )
      }
    )
  }
)

testthat::test_that(
  "estimate_roc() errors on invalid data_source_age type",
  {
    purrr::walk(
      .x = list(NULL, "data", 123, NA, list()),
      .f = function(bad_input) {
        testthat::expect_error(
          estimate_roc(
            data_source_community =
              RRatepol::example_data$pollen_data[[1]],
            data_source_age = bad_input,
            silent = TRUE
          ),
          "'data_source_age' must be one of the following: 'data.frame'"
        )
      }
    )
  }
)

# ================================================================ #
# 2. Character parameter validation                                #
# ================================================================ #

testthat::test_that(
  "estimate_roc() errors on invalid smooth_method",
  {
    purrr::walk(
      .x = list(NULL, 123, NA, "bad_method"),
      .f = function(bad_val) {
        testthat::expect_error(
          estimate_roc(
            data_source_community =
              RRatepol::example_data$pollen_data[[1]],
            data_source_age =
              RRatepol::example_data$sample_age[[1]],
            smooth_method = bad_val,
            silent = TRUE
          )
        )
      }
    )
  }
)

testthat::test_that(
  "estimate_roc() errors on invalid working_units",
  {
    purrr::walk(
      .x = list(NULL, 123, NA, "bad_units"),
      .f = function(bad_val) {
        testthat::expect_error(
          estimate_roc(
            data_source_community =
              RRatepol::example_data$pollen_data[[1]],
            data_source_age =
              RRatepol::example_data$sample_age[[1]],
            working_units = bad_val,
            silent = TRUE
          )
        )
      }
    )
  }
)

testthat::test_that(
  "estimate_roc() errors on invalid dissimilarity_coefficient",
  {
    purrr::walk(
      .x = list(NULL, 123, NA, "bad_dc"),
      .f = function(bad_val) {
        testthat::expect_error(
          estimate_roc(
            data_source_community =
              RRatepol::example_data$pollen_data[[1]],
            data_source_age =
              RRatepol::example_data$sample_age[[1]],
            dissimilarity_coefficient = bad_val,
            silent = TRUE
          )
        )
      }
    )
  }
)

testthat::test_that(
  "estimate_roc() errors on invalid bin_selection",
  {
    purrr::walk(
      .x = list(NULL, 123, NA, "bad_selection"),
      .f = function(bad_val) {
        testthat::expect_error(
          estimate_roc(
            data_source_community =
              RRatepol::example_data$pollen_data[[1]],
            data_source_age =
              RRatepol::example_data$sample_age[[1]],
            working_units = "bins",
            bin_size = 500,
            bin_selection = bad_val,
            silent = TRUE
          )
        )
      }
    )
  }
)

# ================================================================ #
# 3. Logical and numeric parameter validation                      #
# ================================================================ #

testthat::test_that(
  "estimate_roc() errors on invalid logical parameters",
  {
    vec_params <-
      c("standardise", "tranform_to_proportions")

    purrr::walk(
      .x = vec_params,
      .f = function(param) {
        base_args <-
          list(
            data_source_community =
              RRatepol::example_data$pollen_data[[1]],
            data_source_age =
              RRatepol::example_data$sample_age[[1]],
            silent = TRUE
          )
        base_args[[param]] <- "not_logical"

        testthat::expect_error(
          base::do.call(estimate_roc, base_args)
        )
      }
    )
  }
)

testthat::test_that(
  "estimate_roc() errors on invalid numeric parameters",
  {
    vec_params <-
      c("bin_size", "time_standardisation")

    purrr::walk(
      .x = vec_params,
      .f = function(param) {
        base_args <-
          list(
            data_source_community =
              RRatepol::example_data$pollen_data[[1]],
            data_source_age =
              RRatepol::example_data$sample_age[[1]],
            silent = TRUE
          )
        base_args[[param]] <- "not_a_number"

        testthat::expect_error(
          base::do.call(estimate_roc, base_args)
        )
      }
    )
  }
)

# ================================================================ #
# 4. Valid output: structure and column names                      #
# ================================================================ #

testthat::test_that(
  "estimate_roc() returns data.frame with expected columns",
  {
    result <-
      estimate_roc(
        data_source_community =
          RRatepol::example_data$pollen_data[[1]],
        data_source_age =
          RRatepol::example_data$sample_age[[1]],
        smooth_method = "none",
        working_units = "levels",
        bin_size = 500,
        bin_selection = "first",
        standardise = FALSE,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = NULL,
        use_parallel = FALSE,
        silent = TRUE
      )

    testthat::expect_s3_class(result, "data.frame")
    testthat::expect_true(
      base::all(
        c("Working_Unit", "Age", "ROC", "ROC_up", "ROC_dw") %in%
          base::colnames(result)
      )
    )
  }
)

testthat::test_that(
  "estimate_roc() returns data.frame for all working_units",
  {
    purrr::walk(
      .x = c("levels", "bins", "MW"),
      .f = function(wu) {
        result <-
          estimate_roc(
            data_source_community =
              RRatepol::example_data$pollen_data[[1]],
            data_source_age =
              RRatepol::example_data$sample_age[[1]],
            smooth_method = "none",
            working_units = wu,
            bin_size = 500,
            number_of_shifts = 1,
            bin_selection = "first",
            standardise = FALSE,
            dissimilarity_coefficient = "euc",
            tranform_to_proportions = TRUE,
            rand = NULL,
            use_parallel = FALSE,
            silent = TRUE
          )

        testthat::expect_s3_class(result, "data.frame")
        testthat::expect_true(base::nrow(result) > 0)
      }
    )
  }
)

testthat::test_that(
  "estimate_roc() returns data.frame for all dissimilarity coefficients",
  {
    purrr::walk(
      .x = c("euc", "euc.sd", "chord", "chisq", "gower", "bray"),
      .f = function(dc) {
        result <-
          estimate_roc(
            data_source_community =
              RRatepol::example_data$pollen_data[[1]],
            data_source_age =
              RRatepol::example_data$sample_age[[1]],
            smooth_method = "none",
            working_units = "levels",
            bin_size = 500,
            bin_selection = "first",
            standardise = FALSE,
            dissimilarity_coefficient = dc,
            tranform_to_proportions = TRUE,
            rand = NULL,
            use_parallel = FALSE,
            silent = TRUE
          )

        testthat::expect_s3_class(result, "data.frame")
        testthat::expect_true(base::nrow(result) > 0)
      }
    )
  }
)

testthat::test_that(
  "estimate_roc() returns data.frame for all smooth methods",
  {
    purrr::walk(
      .x = c("none", "m.avg", "grim", "age.w", "shep"),
      .f = function(sm) {
        result <-
          estimate_roc(
            data_source_community =
              RRatepol::example_data$pollen_data[[1]],
            data_source_age =
              RRatepol::example_data$sample_age[[1]],
            smooth_method = sm,
            smooth_n_points = 5,
            smooth_age_range = 500,
            smooth_n_max = 9,
            working_units = "levels",
            bin_size = 500,
            bin_selection = "first",
            standardise = FALSE,
            dissimilarity_coefficient = "euc",
            tranform_to_proportions = TRUE,
            rand = NULL,
            use_parallel = FALSE,
            silent = TRUE
          )

        testthat::expect_s3_class(result, "data.frame")
        testthat::expect_true(base::nrow(result) > 0)
      }
    )
  }
)

testthat::test_that(
  "estimate_roc() errors when standardise has length > 1",
  {
    testthat::expect_error(
      estimate_roc(
        data_source_community =
          RRatepol::example_data$pollen_data[[1]],
        data_source_age =
          RRatepol::example_data$sample_age[[1]],
        standardise = c(TRUE, FALSE),
        use_parallel = FALSE,
        silent = TRUE
      ),
      "'standardise' must be a single TRUE or FALSE"
    )
  }
)

testthat::test_that(
  "estimate_roc() errors when tranform_to_proportions has length > 1",
  {
    testthat::expect_error(
      estimate_roc(
        data_source_community =
          RRatepol::example_data$pollen_data[[1]],
        data_source_age =
          RRatepol::example_data$sample_age[[1]],
        tranform_to_proportions = c(TRUE, FALSE),
        use_parallel = FALSE,
        silent = TRUE
      ),
      "'tranform_to_proportions' must be a single TRUE or FALSE"
    )
  }
)

testthat::test_that(
  "estimate_roc() errors when verbose has length > 1",
  {
    testthat::expect_error(
      estimate_roc(
        data_source_community =
          RRatepol::example_data$pollen_data[[1]],
        data_source_age =
          RRatepol::example_data$sample_age[[1]],
        verbose = c(TRUE, FALSE),
        use_parallel = FALSE,
        silent = TRUE
      ),
      "'verbose' must be a single TRUE or FALSE"
    )
  }
)

testthat::test_that(
  "estimate_roc() errors when use_parallel is 0",
  {
    testthat::expect_error(
      estimate_roc(
        data_source_community =
          RRatepol::example_data$pollen_data[[1]],
        data_source_age =
          RRatepol::example_data$sample_age[[1]],
        use_parallel = 0,
        silent = TRUE
      ),
      "'use_parallel' must not be 0 or NA when numeric"
    )
  }
)

testthat::test_that(
  "estimate_roc() errors when interest_threshold has length > 1",
  {
    testthat::expect_error(
      estimate_roc(
        data_source_community =
          RRatepol::example_data$pollen_data[[1]],
        data_source_age =
          RRatepol::example_data$sample_age[[1]],
        interest_threshold = c(100, 200),
        use_parallel = FALSE,
        silent = TRUE
      ),
      "'interest_threshold' must be a single value"
    )
  }
)
