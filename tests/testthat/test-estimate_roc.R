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
  "estimate_roc() errors when use_parallel is negative",
  {
    testthat::expect_error(
      estimate_roc(
        data_source_community =
          purrr::chuck(RRatepol::example_data, "pollen_data", 1),
        data_source_age =
          purrr::chuck(RRatepol::example_data, "sample_age", 1),
        use_parallel = -2,
        silent = TRUE
      ),
      "'use_parallel' must be a positive number of cores when numeric"
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

# ================================================================ #
# 5. Parallel computation                                          #
# ================================================================ #

# Workers load the installed RRatepol, so these tests need the package
#    installed. CRAN limits checks to two cores.

testthat::test_that(
  "estimate_roc() runs with a numeric use_parallel",
  {
    testthat::skip_on_cran()
    skip_if_rratepol_not_installed()

    res_serial <-
      estimate_roc(
        data_source_community =
          purrr::chuck(RRatepol::example_data, "pollen_data", 1),
        data_source_age =
          purrr::chuck(RRatepol::example_data, "sample_age", 1),
        smooth_method = "none",
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 2,
        bin_selection = "random",
        rand = 4,
        use_parallel = FALSE,
        silent = TRUE
      )

    res_parallel <-
      estimate_roc(
        data_source_community =
          purrr::chuck(RRatepol::example_data, "pollen_data", 1),
        data_source_age =
          purrr::chuck(RRatepol::example_data, "sample_age", 1),
        smooth_method = "none",
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 2,
        bin_selection = "random",
        rand = 4,
        use_parallel = 2,
        silent = TRUE
      )

    testthat::expect_s3_class(res_parallel, "data.frame")

    testthat::expect_identical(
      base::colnames(res_parallel),
      base::colnames(res_serial)
    )

    testthat::expect_identical(
      base::nrow(res_parallel),
      base::nrow(res_serial)
    )
  }
)

# A numeric use_parallel that silently ran on one core gave results
#    identical to a serial run with the same seed. Worker RNG streams differ
#    from the main session, so a real cluster gives different results.
testthat::test_that(
  "estimate_roc() numeric use_parallel runs on a cluster",
  {
    testthat::skip_on_cran()
    skip_if_rratepol_not_installed()

    list_res <-
      purrr::map(
        .x = list(FALSE, 2),
        .f = ~ {
          base::set.seed(900723)

          estimate_roc(
            data_source_community =
              purrr::chuck(RRatepol::example_data, "pollen_data", 1),
            data_source_age =
              purrr::chuck(RRatepol::example_data, "sample_age", 1),
            smooth_method = "none",
            working_units = "MW",
            bin_size = 500,
            number_of_shifts = 2,
            bin_selection = "random",
            standardise = TRUE,
            rand = 4,
            use_parallel = .x,
            silent = TRUE
          )
        }
      )

    testthat::expect_false(
      base::identical(
        dplyr::pull(purrr::chuck(list_res, 1), "ROC"),
        dplyr::pull(purrr::chuck(list_res, 2), "ROC")
      )
    )
  }
)

testthat::test_that(
  "estimate_roc() use_parallel = 1 gives the single-core result",
  {
    list_res <-
      purrr::map(
        .x = list(FALSE, 1),
        .f = ~ {
          base::set.seed(900723)

          estimate_roc(
            data_source_community =
              purrr::chuck(RRatepol::example_data, "pollen_data", 1),
            data_source_age =
              purrr::chuck(RRatepol::example_data, "sample_age", 1),
            smooth_method = "none",
            working_units = "MW",
            bin_size = 500,
            number_of_shifts = 2,
            bin_selection = "random",
            standardise = TRUE,
            rand = 4,
            use_parallel = .x,
            silent = TRUE
          )
        }
      )

    testthat::expect_identical(
      purrr::chuck(list_res, 1),
      purrr::chuck(list_res, 2)
    )
  }
)

testthat::test_that(
  "estimate_roc() parallel results are reproducible with a seed",
  {
    testthat::skip_on_cran()
    skip_if_rratepol_not_installed()

    list_res <-
      purrr::map(
        .x = 1:2,
        .f = ~ {
          base::set.seed(900723)

          estimate_roc(
            data_source_community =
              purrr::chuck(RRatepol::example_data, "pollen_data", 1),
            data_source_age =
              purrr::chuck(RRatepol::example_data, "sample_age", 1),
            smooth_method = "none",
            working_units = "MW",
            bin_size = 500,
            number_of_shifts = 2,
            bin_selection = "random",
            rand = 4,
            use_parallel = 2,
            silent = TRUE
          )
        }
      )

    testthat::expect_identical(
      purrr::chuck(list_res, 1),
      purrr::chuck(list_res, 2)
    )
  }
)

testthat::test_that(
  "estimate_roc() parallel run leaves the main RNG kind unchanged",
  {
    testthat::skip_on_cran()
    skip_if_rratepol_not_installed()

    vec_rng_before <-
      base::RNGkind()

    res <-
      estimate_roc(
        data_source_community =
          purrr::chuck(RRatepol::example_data, "pollen_data", 1),
        data_source_age =
          purrr::chuck(RRatepol::example_data, "sample_age", 1),
        smooth_method = "none",
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 2,
        bin_selection = "random",
        rand = 4,
        use_parallel = 2,
        silent = TRUE
      )

    testthat::expect_identical(base::RNGkind(), vec_rng_before)
  }
)

# Opt-in timing test. Elapsed time depends on the load of the machine, and
#    starting the workers costs about a second, so a speed-up shows only
#    with enough runs. The test takes about 25 s and runs only when
#    the environment variable RRATEPOL_TEST_TIMING is "true":
#    Sys.setenv(RRATEPOL_TEST_TIMING = "true")
testthat::test_that(
  "estimate_roc() runs faster on 2 cores than on one",
  {
    testthat::skip_if_not(
      base::identical(base::Sys.getenv("RRATEPOL_TEST_TIMING"), "true"),
      message = "Set RRATEPOL_TEST_TIMING=true to run the timing test"
    )
    testthat::skip_on_cran()
    testthat::skip_on_ci()
    skip_if_rratepol_not_installed()

    n_cores <-
      parallel::detectCores()

    testthat::skip_if(
      base::is.na(n_cores) || n_cores < 2,
      message = "Fewer than 2 cores detected"
    )

    vec_elapsed <-
      purrr::map_dbl(
        .x = list(FALSE, 2),
        .f = ~ {
          base::set.seed(900723)

          time_run <-
            base::system.time(
              estimate_roc(
                data_source_community =
                  purrr::chuck(RRatepol::example_data, "pollen_data", 1),
                data_source_age =
                  purrr::chuck(RRatepol::example_data, "sample_age", 1),
                age_uncertainty =
                  purrr::chuck(RRatepol::example_data, "age_uncertainty", 1),
                smooth_method = "shep",
                working_units = "MW",
                bin_size = 500,
                number_of_shifts = 5,
                standardise = TRUE,
                rand = 100,
                use_parallel = .x,
                silent = TRUE
              )
            )

          purrr::chuck(time_run, "elapsed")
        }
      )

    testthat::expect_lt(
      purrr::chuck(vec_elapsed, 2),
      0.8 * purrr::chuck(vec_elapsed, 1)
    )
  }
)

# ================================================================ #
# 6. Session settings                                              #
# ================================================================ #

testthat::test_that(
  "estimate_roc() keeps the user's progress-bar setting",
  {
    list_pb_before <-
      pbapply::pboptions(type = "txt")

    res <-
      estimate_roc(
        data_source_community =
          purrr::chuck(RRatepol::example_data, "pollen_data", 1),
        data_source_age =
          purrr::chuck(RRatepol::example_data, "sample_age", 1),
        smooth_method = "none",
        working_units = "levels",
        rand = NULL,
        use_parallel = FALSE,
        silent = TRUE
      )

    type_after <-
      purrr::chuck(pbapply::pboptions(), "type")

    pbapply::pboptions(list_pb_before)

    testthat::expect_identical(type_after, "txt")
  }
)

testthat::test_that(
  "estimate_roc() restores the progress-bar setting after an error",
  {
    testthat::local_mocked_bindings(
      run_iteration = function(...) {
        base::stop("mocked run failure")
      }
    )

    list_pb_before <-
      pbapply::pboptions(type = "txt")

    testthat::expect_error(
      estimate_roc(
        data_source_community =
          purrr::chuck(RRatepol::example_data, "pollen_data", 1),
        data_source_age =
          purrr::chuck(RRatepol::example_data, "sample_age", 1),
        smooth_method = "none",
        working_units = "levels",
        rand = NULL,
        use_parallel = FALSE,
        silent = TRUE
      ),
      "mocked run failure"
    )

    type_after <-
      purrr::chuck(pbapply::pboptions(), "type")

    pbapply::pboptions(list_pb_before)

    testthat::expect_identical(type_after, "txt")
  }
)
