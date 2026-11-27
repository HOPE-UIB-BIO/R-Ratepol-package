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

# ================================================================ #
# 3. Random-number streams                                         #
# ================================================================ #

testthat::test_that(
  "run_iteration() gives the same result for the same rng_stream",
  {
    base::set.seed(900723)

    vec_seed_user <-
      base::get(".Random.seed", envir = base::globalenv())

    data_run <-
      make_run_data()

    data_run_stream <-
      c(
        data_run,
        list(
          rng_stream = purrr::chuck(
            util_make_rng_streams(n_streams = 1, seed = 900723),
            1
          )
        )
      )

    # Different generator states before each call: only the stream
    #    should determine the random draws of the run
    list_res <-
      purrr::map(
        .x = list(
          list(data = data_run, seed = 1),
          list(data = data_run, seed = 2),
          list(data = data_run_stream, seed = 1),
          list(data = data_run_stream, seed = 2)
        ),
        .f = ~ {
          base::set.seed(purrr::chuck(.x, "seed"))

          run_iteration(
            data_source_run = purrr::chuck(.x, "data"),
            bin_selection = "random",
            standardise = TRUE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = "euc",
            time_standardisation = 500,
            silent = TRUE
          )
        }
      )

    # Restore the generator of the test session
    base::assign(
      x = ".Random.seed",
      value = vec_seed_user,
      envir = base::globalenv()
    )

    # Without a stream, the result follows the generator state
    testthat::expect_false(
      base::identical(
        purrr::chuck(list_res, 1),
        purrr::chuck(list_res, 2)
      )
    )

    # With a stream, it does not
    testthat::expect_identical(
      purrr::chuck(list_res, 3),
      purrr::chuck(list_res, 4)
    )
  }
)
