# Shared test fixtures for RRatepol tests
# Called automatically by testthat before any test file

# Build minimal extracted data (output of extract_data)
make_extracted_data <- function() {
  suppressWarnings(
    extract_data(
      data_community_extract =
        RRatepol::example_data$pollen_data[[1]],
      data_age_extract =
        RRatepol::example_data$sample_age[[1]],
      silent = TRUE
    )
  )
}

# Build data_run[[1]] - list with $data and $bins
# for use in subset_samples, run_iteration, and related tests
make_run_data <- function() {
  data_prepared <-
    suppressWarnings(
      prepare_data(
        data_source_prep = make_extracted_data(),
        working_units = "bins",
        bin_size = 500,
        rand = 1
      )
    )

  purrr::chuck(
    RUtilpol::flatten_list_by_one(data_prepared),
    1
  )
}

# Build a minimal estimate_roc() result for detect_peak_points /
# detect_sni tests. Runs once per test session (results are
# reproducible because rand = 1 and a fixed seed is used).
make_roc_data <- function() {
  estimate_roc(
    data_source_community =
      RRatepol::example_data$pollen_data[[1]],
    data_source_age =
      RRatepol::example_data$sample_age[[1]],
    smooth_method = "shep",
    smooth_n_points = 5,
    working_units = "levels",
    bin_size = 500,
    number_of_shifts = 1,
    bin_selection = "first",
    standardise = FALSE,
    n_individuals = 150,
    dissimilarity_coefficient = "euc",
    tranform_to_proportions = TRUE,
    rand = 1,
    use_parallel = FALSE,
    silent = TRUE
  )
}

# Build pre-processed community data suitable as input to
# estimate_dissimilarity_coefficient().
make_dc_data <- function(sel_method = "proportions") {
  data_extracted <-
    suppressWarnings(
      extract_data(
        data_community_extract =
          RRatepol::example_data$pollen_data[[1]],
        data_age_extract =
          RRatepol::example_data$sample_age[[1]],
        age_uncertainty =
          RRatepol::example_data$age_uncertainty[[1]],
        silent = TRUE
      )
    )

  data_smooth <-
    suppressWarnings(
      smooth_community_data(
        data_source_smooth = data_extracted,
        smooth_method = "shep"
      )
    )

  data_reduced <-
    reduce_data(
      data_source_reduce = data_smooth,
      check_taxa = TRUE,
      check_levels = TRUE
    )

  data_prepared <-
    suppressWarnings(
      prepare_data(
        data_source_prep = data_reduced,
        working_units = "bins",
        bin_size = 500,
        rand = 1
      )
    )

  data_flat <-
    purrr::chuck(
      RUtilpol::flatten_list_by_one(data_prepared),
      1
    )

  data_subset <-
    subset_samples(
      data_source_subset = purrr::chuck(data_flat, "data"),
      data_source_bins = purrr::chuck(data_flat, "bins"),
      bin_selection = "first"
    ) |>
    (\(x) reduce_data_simple(data_source_reduce = x))()

  n_individuals <-
    min(
      c(
        rowSums(
          subset_community(data_subset),
          na.rm = TRUE
        ),
        150
      )
    )

  data_subset <-
    reduce_data_simple(
      data_source_reduce =
        data_subset[
          rowSums(
            subset_community(data_subset),
            na.rm = TRUE
          ) >= n_individuals,
        ]
    )

  base::set.seed(123)
  standardise_community_data(
    data_source_standard = data_subset,
    n_individuals = n_individuals
  ) |>
    (\(x) reduce_data_simple(data_source_reduce = x))() |>
    transform_into_proportions(
      sel_method = sel_method,
      silent = TRUE
    )
}

# Build a minimal estimate_roc() + detect_peak_points() result
# suitable as input to plot_roc(). Uses the lightweight make_roc_data()
# (rand = 1) so validation tests run quickly.
make_plot_roc_data <- function() {
  make_roc_data() |>
    detect_peak_points(
      sel_method = "trend_linear"
    )
}
