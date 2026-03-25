# ================================================================ #
# 1. Input validation                                              #
# ================================================================ #

testthat::test_that(
  "standardise_community_data() errors on invalid data type",
  {
    purrr::walk(
      .x = list(NULL, "my_data", 123, list()),
      .f = function(bad_input) {
        testthat::expect_error(
          standardise_community_data(
            data_source_standard = bad_input
          ),
          "'data_source_standard' must be one of the following: 'data.frame'"
        )
      }
    )
  }
)

# ================================================================ #
# 2. Valid output                                                  #
# ================================================================ #

testthat::test_that(
  "standardise_community_data() returns data.frame with row sums == n",
  {
    data_run <-
      make_run_data()

    data_subset <-
      reduce_data_simple(
        data_source_reduce =
          subset_samples(
            data_source_subset = purrr::chuck(data_run, "data"),
            data_source_bins = purrr::chuck(data_run, "bins"),
            bin_selection = "first"
          )
      )

    n_ind <-
      base::min(
        c(
          base::rowSums(
            subset_community(data_subset),
            na.rm = TRUE
          ),
          150
        )
      )

    data_keep <-
      data_subset[
        base::rowSums(
          subset_community(data_subset),
          na.rm = TRUE
        ) >= n_ind,
      ]

    base::set.seed(900723)
    result <-
      standardise_community_data(
        data_source_standard = data_keep,
        n_individuals = n_ind
      )

    testthat::expect_s3_class(result, "data.frame")
    testthat::expect_true(
      base::all(
        base::rowSums(
          subset_community(result),
          na.rm = TRUE
        ) == n_ind
      )
    )
  }
)

testthat::test_that(
  "standardise_community_data() errors when n_individuals is 0",
  {
    data_run <-
      make_run_data()

    data_subset <-
      reduce_data_simple(
        data_source_reduce =
          subset_samples(
            data_source_subset = purrr::chuck(data_run, "data"),
            data_source_bins = purrr::chuck(data_run, "bins"),
            bin_selection = "first"
          )
      )

    testthat::expect_error(
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = 0
      ),
      "'n_individuals' must be a positive number"
    )
  }
)

testthat::test_that(
  "standardise_community_data() errors when n_individuals is negative",
  {
    data_run <-
      make_run_data()

    data_subset <-
      reduce_data_simple(
        data_source_reduce =
          subset_samples(
            data_source_subset = purrr::chuck(data_run, "data"),
            data_source_bins = purrr::chuck(data_run, "bins"),
            bin_selection = "first"
          )
      )

    testthat::expect_error(
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = -10
      ),
      "'n_individuals' must be a positive number"
    )
  }
)

testthat::test_that(
  "standardise_community_data() errors when n_individuals is NA",
  {
    data_run <-
      make_run_data()

    data_subset <-
      reduce_data_simple(
        data_source_reduce =
          subset_samples(
            data_source_subset = purrr::chuck(data_run, "data"),
            data_source_bins = purrr::chuck(data_run, "bins"),
            bin_selection = "first"
          )
      )

    testthat::expect_error(
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = NA_real_
      ),
      "'n_individuals' must be a positive number"
    )
  }
)
