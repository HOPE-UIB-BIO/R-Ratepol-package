# Setup workflow within estimate_roc and run_iteration
# Prepare default data (with shep smoothing)

# Input validation (Error messages)
test_that(
  "transform_into_proportions throws error when data_source_trans parameter is missing for sel_method='proportions'",
  {
    expect_error(
      transform_into_proportions(
        sel_method = "proportions"
      ),
      'argument "data_source_trans" is missing, with no default'
    )
  }
)

test_that(
  "transform_into_proportions throws error when data_source_trans parameter is missing for sel_method='percentages'",
  {
    expect_error(
      transform_into_proportions(
        sel_method = "percentages"
      ),
      'argument "data_source_trans" is missing, with no default'
    )
  }
)

# NULL
test_that(
  "transform_into_proportions throws error when data_source_trans parameter is NULL for sel_method='proportions'",
  {
    expect_error(
      transform_into_proportions(
        data_source_trans = NULL,
        sel_method = "proportions"
      ),
      "'data_source_trans' must be one of the following: 'data.frame'"
    )
  }
)

test_that(
  "transform_into_proportions throws error when data_source_trans parameter is NULL for sel_method='percentages'",
  {
    expect_error(
      transform_into_proportions(
        data_source_trans = NULL,
        sel_method = "percentages"
      ),
      "'data_source_trans' must be one of the following: 'data.frame'"
    )
  }
)

# NA
test_that(
  "transform_into_proportions throws error when data_source_trans parameter is NA for sel_method='proportions'",
  {
    expect_error(
      transform_into_proportions(
        data_source_trans = NA,
        sel_method = "proportions"
      ),
      "'data_source_trans' must be one of the following: 'data.frame'"
    )
  }
)

test_that(
  "transform_into_proportions throws error when data_source_trans parameter is NA for sel_method='percentages'",
  {
    expect_error(
      transform_into_proportions(
        data_source_trans = NA,
        sel_method = "percentages"
      ),
      "'data_source_trans' must be one of the following: 'data.frame'"
    )
  }
)

# Character
test_that(
  "transform_into_proportions throws error when data_source_trans parameter is character for sel_method='proportions'",
  {
    expect_error(
      transform_into_proportions(
        data_source_trans = "my_data",
        sel_method = "proportions"
      ),
      "'data_source_trans' must be one of the following: 'data.frame'"
    )
  }
)

test_that(
  "transform_into_proportions throws error when data_source_trans parameter is character for sel_method='percentages'",
  {
    expect_error(
      transform_into_proportions(
        data_source_trans = "my_data",
        sel_method = "percentages"
      ),
      "'data_source_trans' must be one of the following: 'data.frame'"
    )
  }
)

# Numeric
test_that(
  "transform_into_proportions throws error when data_source_trans parameter is numeric for sel_method='proportions'",
  {
    expect_error(
      transform_into_proportions(
        data_source_trans = 123,
        sel_method = "proportions"
      ),
      "'data_source_trans' must be one of the following: 'data.frame'"
    )
  }
)

test_that(
  "transform_into_proportions throws error when data_source_trans parameter is numeric for sel_method='percentages'",
  {
    expect_error(
      transform_into_proportions(
        data_source_trans = 123,
        sel_method = "percentages"
      ),
      "'data_source_trans' must be one of the following: 'data.frame'"
    )
  }
)

# list
test_that(
  "transform_into_proportions throws error when data_source_trans parameter is list for sel_method='proportions'",
  {
    expect_error(
      transform_into_proportions(
        data_source_trans = list(),
        sel_method = "proportions"
      ),
      "'data_source_trans' must be one of the following: 'data.frame'"
    )
  }
)

test_that(
  "transform_into_proportions throws error when data_source_trans parameter is list for sel_method='percentages'",
  {
    expect_error(
      transform_into_proportions(
        data_source_trans = list(),
        sel_method = "percentages"
      ),
      "'data_source_trans' must be one of the following: 'data.frame'"
    )
  }
)

# Matrix
test_that(
  "transform_into_proportions throws error when data_source_trans parameter is matrix for sel_method='proportions'",
  {
    expect_error(
      transform_into_proportions(
        data_source_trans = matrix(),
        sel_method = "proportions"
      ),
      "'data_source_trans' must be one of the following: 'data.frame'"
    )
  }
)

test_that(
  "transform_into_proportions throws error when data_source_trans parameter is matrix for sel_method='percentages'",
  {
    expect_error(
      transform_into_proportions(
        data_source_trans = matrix(),
        sel_method = "percentages"
      ),
      "'data_source_trans' must be one of the following: 'data.frame'"
    )
  }
)

# empty dataframe (fails)
test_that(
  "transform_into_proportions handles empty data.frame for data_source_trans parameter with sel_method='proportions'",
  {
    expect_error(
      transform_into_proportions(
        data_source_trans = data.frame(),
        sel_method = "proportions"
      ),
      # none programmed into the function. returns empty result data.frame
    )
  }
)

test_that(
  "transform_into_proportions handles empty data.frame for data_source_trans parameter with sel_method='percentages'",
  {
    expect_error(
      transform_into_proportions(
        data_source_trans = data.frame(),
        sel_method = "percentages"
      ),
      # none programmed into the function. returns empty result data.frame
    )
  }
)

# 0
test_that(
  "transform_into_proportions throws error when data_source_trans parameter is zero for sel_method='proportions'",
  {
    expect_error(
      transform_into_proportions(
        data_source_trans = 0,
        sel_method = "proportions"
      ),
      "'data_source_trans' must be one of the following: 'data.frame'"
    )
  }
)

test_that(
  "transform_into_proportions throws error when data_source_trans parameter is zero for sel_method='percentages'",
  {
    expect_error(
      transform_into_proportions(
        data_source_trans = 0,
        sel_method = "percentages"
      ),
      "'data_source_trans' must be one of the following: 'data.frame'"
    )
  }
)


# Input validation for sel_method
# Invalid character
test_that(
  "transform_into_proportions throws error when sel_method parameter contains invalid character value",
  {
    suppressWarnings(
      data_to_run_bins <-
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
        ) %>%
        smooth_community_data(
          smooth_method = "shep"
        ) %>%
        reduce_data(
          check_taxa = TRUE,
          check_levels = TRUE
        ) %>%
        prepare_data(
          data_source_prep = .,
          working_units = "bins",
          bin_size = 500,
          rand = 1
        ) %>%
        RUtilpol::flatten_list_by_one() %>%
        .[[1]]
    )

    data_source_subset <-
      data_to_run_bins$data
    data_source_bins <-
      data_to_run_bins$bins

    data_subset <-
      subset_samples(
        data_source_subset = data_source_subset,
        data_source_bins = data_source_bins,
        bin_selection = "first"
      ) %>%
      reduce_data_simple()

    # standardisation
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      )
    data_sd <-
      reduce_data_simple(
        data_source_reduce = data_sd
      )
    # Test
    expect_error(
      transform_into_proportions(
        data_source_trans = data_sd,
        sel_method = "invalid_method"
      ),
      "'sel_method' must contains one of the following values: 'percentages', 'proportions'"
    )
  }
)

# Multiple methods
test_that(
  "transform_into_proportions throws error when sel_method parameter contains multiple values",
  {
    suppressWarnings(
      data_to_run_bins <-
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
        ) %>%
        smooth_community_data(
          smooth_method = "shep"
        ) %>%
        reduce_data(
          check_taxa = TRUE,
          check_levels = TRUE
        ) %>%
        prepare_data(
          data_source_prep = .,
          working_units = "bins",
          bin_size = 500,
          rand = 1
        ) %>%
        RUtilpol::flatten_list_by_one() %>%
        .[[1]]
    )

    data_source_subset <-
      data_to_run_bins$data
    data_source_bins <-
      data_to_run_bins$bins

    data_subset <-
      subset_samples(
        data_source_subset = data_source_subset,
        data_source_bins = data_source_bins,
        bin_selection = "first"
      ) %>%
      reduce_data_simple()

    # standardisation
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      )
    data_sd <-
      reduce_data_simple(
        data_source_reduce = data_sd
      )
    # Test
    expect_error(
      res <-
        transform_into_proportions(
          data_source_trans = data_sd,
          sel_method = c("proportions", "percentages"),
          verbose = TRUE
        ),
      # not programmed into the function yet - apparently
      "'arg' must be of length 1"
    )
  }
)

test_that(
  "transform_into_proportions throws error when sel_method parameter contains multiple values in different order",
  {
    suppressWarnings(
      data_to_run_bins <-
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
        ) %>%
        smooth_community_data(
          smooth_method = "shep"
        ) %>%
        reduce_data(
          check_taxa = TRUE,
          check_levels = TRUE
        ) %>%
        prepare_data(
          data_source_prep = .,
          working_units = "bins",
          bin_size = 500,
          rand = 1
        ) %>%
        RUtilpol::flatten_list_by_one() %>%
        .[[1]]
    )

    data_source_subset <-
      data_to_run_bins$data
    data_source_bins <-
      data_to_run_bins$bins

    data_subset <-
      subset_samples(
        data_source_subset = data_source_subset,
        data_source_bins = data_source_bins,
        bin_selection = "first"
      ) %>%
      reduce_data_simple()

    # standardisation
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      )
    data_sd <-
      reduce_data_simple(
        data_source_reduce = data_sd
      )
    # Test
    expect_error(
      res <-
        transform_into_proportions(
          data_source_trans = data_sd,
          sel_method = c("percentages", "proportions"),
          verbose = TRUE
        ),
      "'arg' must be of length 1"
    )
  }
)

# Numeric
test_that(
  "transform_into_proportions throws error when sel_method parameter is numeric instead of character",
  {
    suppressWarnings(
      data_to_run_bins <-
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
        ) %>%
        smooth_community_data(
          smooth_method = "shep"
        ) %>%
        reduce_data(
          check_taxa = TRUE,
          check_levels = TRUE
        ) %>%
        prepare_data(
          data_source_prep = .,
          working_units = "bins",
          bin_size = 500,
          rand = 1
        ) %>%
        RUtilpol::flatten_list_by_one() %>%
        .[[1]]
    )

    data_source_subset <-
      data_to_run_bins$data
    data_source_bins <-
      data_to_run_bins$bins

    data_subset <-
      subset_samples(
        data_source_subset = data_source_subset,
        data_source_bins = data_source_bins,
        bin_selection = "first"
      ) %>%
      reduce_data_simple()

    # standardisation
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      )
    data_sd <-
      reduce_data_simple(
        data_source_reduce = data_sd
      )
    # Test
    expect_error(
      res <-
        transform_into_proportions(
          data_source_trans = data_sd,
          sel_method = 123,
          verbose = TRUE
        ),
      "'sel_method' must be one of the following: 'character'"
    )
  }
)


# NULL
test_that(
  "transform_into_proportions throws error when sel_method parameter is NULL",
  {
    suppressWarnings(
      data_to_run_bins <-
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
        ) %>%
        smooth_community_data(
          smooth_method = "shep"
        ) %>%
        reduce_data(
          check_taxa = TRUE,
          check_levels = TRUE
        ) %>%
        prepare_data(
          data_source_prep = .,
          working_units = "bins",
          bin_size = 500,
          rand = 1
        ) %>%
        RUtilpol::flatten_list_by_one() %>%
        .[[1]]
    )

    data_source_subset <-
      data_to_run_bins$data
    data_source_bins <-
      data_to_run_bins$bins

    data_subset <-
      subset_samples(
        data_source_subset = data_source_subset,
        data_source_bins = data_source_bins,
        bin_selection = "first"
      ) %>%
      reduce_data_simple()

    # standardisation
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      )
    data_sd <-
      reduce_data_simple(
        data_source_reduce = data_sd
      )

    # transformation into proportions / percentages
    tranform_to_proportions <-
      TRUE
    expect_error(
      data_sd_prop <-
        transform_into_proportions(
          data_source_trans = data_sd,
          sel_method = NULL,
          verbose = FALSE
        ),
      "'sel_method' must be one of the following: 'character'"
    )
  }
)

# Empty
test_that(
  "transform_into_proportions throws warning when sel_method parameter is missing (uses default)",
  {
    suppressWarnings(
      data_to_run_bins <-
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
        ) %>%
        smooth_community_data(
          smooth_method = "shep"
        ) %>%
        reduce_data(
          check_taxa = TRUE,
          check_levels = TRUE
        ) %>%
        prepare_data(
          data_source_prep = .,
          working_units = "bins",
          bin_size = 500,
          rand = 1
        ) %>%
        RUtilpol::flatten_list_by_one() %>%
        .[[1]]
    )

    data_source_subset <-
      data_to_run_bins$data
    data_source_bins <-
      data_to_run_bins$bins

    data_subset <-
      subset_samples(
        data_source_subset = data_source_subset,
        data_source_bins = data_source_bins,
        bin_selection = "first"
      ) %>%
      reduce_data_simple()

    # standardisation
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      )
    data_sd <-
      reduce_data_simple(
        data_source_reduce = data_sd
      )

    # transformation into proportions / percentages
    tranform_to_proportions <-
      TRUE
    expect_warning(
      data_sd_prop <-
        transform_into_proportions(
          data_source_trans = data_sd,
          verbose = TRUE
        ),
      # none programmed into the function yet
      # e.g., 'Warning: No sel_method supplied. Using default "proportions"'
    )
  }
)

# Output validation:
# Valid data:
test_that(
  "transform_into_proportions returns values between 0-1 when sel_method='proportions' with valid data_source_trans",
  {
    suppressWarnings(
      data_to_run_bins <-
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
        ) %>%
        smooth_community_data(
          smooth_method = "shep"
        ) %>%
        reduce_data(
          check_taxa = TRUE,
          check_levels = TRUE
        ) %>%
        prepare_data(
          data_source_prep = .,
          working_units = "bins",
          bin_size = 500,
          rand = 1
        ) %>%
        RUtilpol::flatten_list_by_one() %>%
        .[[1]]
    )

    data_source_subset <-
      data_to_run_bins$data
    data_source_bins <-
      data_to_run_bins$bins

    data_subset <-
      subset_samples(
        data_source_subset = data_source_subset,
        data_source_bins = data_source_bins,
        bin_selection = "first"
      ) %>%
      reduce_data_simple()

    # standardisation
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      )
    data_sd <-
      reduce_data_simple(
        data_source_reduce = data_sd
      )

    # transformation into proportions / percentages
    tranform_to_proportions <-
      TRUE
    data_sd_prop <-
      transform_into_proportions(
        data_source_trans = data_sd,
        sel_method = "proportions", # or "percentages"
        verbose = FALSE
      )

    expect_true(
      all(
        data_sd_prop[, -c(1:3)] >= 0 &
          data_sd_prop[, -c(1:3)] <= 1
      )
    )
  }
)

test_that(
  "transform_into_proportions returns values between 0-100 when sel_method='percentages' with valid data_source_trans",
  {
    suppressWarnings(
      data_to_run_bins <-
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
        ) %>%
        smooth_community_data(
          smooth_method = "shep"
        ) %>%
        reduce_data(
          check_taxa = TRUE,
          check_levels = TRUE
        ) %>%
        prepare_data(
          data_source_prep = .,
          working_units = "bins",
          bin_size = 500,
          rand = 1
        ) %>%
        RUtilpol::flatten_list_by_one() %>%
        .[[1]]
    )

    data_source_subset <-
      data_to_run_bins$data
    data_source_bins <-
      data_to_run_bins$bins

    data_subset <-
      subset_samples(
        data_source_subset = data_source_subset,
        data_source_bins = data_source_bins,
        bin_selection = "first"
      ) %>%
      reduce_data_simple()

    # standardisation
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      )
    data_sd <-
      reduce_data_simple(
        data_source_reduce = data_sd
      )

    # transformation into proportions / percentages
    tranform_to_proportions <-
      TRUE
    data_sd_prop <-
      transform_into_proportions(
        data_source_trans = data_sd,
        sel_method = "percentages",
        verbose = FALSE
      )

    expect_true(
      all(
        data_sd_prop[, -c(1:3)] >= 0 &
          data_sd_prop[, -c(1:3)] <= 100
      )
    )
  }
)



# check mathematical accuracy of proportions and percentages
test_that(
  "transform_into_proportions calculates mathematically correct percentages when sel_method='percentages'",
  {
    suppressWarnings(
      data_to_run_bins <-
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
        ) %>%
        smooth_community_data(
          smooth_method = "shep"
        ) %>%
        reduce_data(
          check_taxa = TRUE,
          check_levels = TRUE
        ) %>%
        prepare_data(
          data_source_prep = .,
          working_units = "bins",
          bin_size = 500,
          rand = 1
        ) %>%
        RUtilpol::flatten_list_by_one() %>%
        .[[1]]
    )

    data_source_subset <-
      data_to_run_bins$data
    data_source_bins <-
      data_to_run_bins$bins

    data_subset <-
      subset_samples(
        data_source_subset = data_source_subset,
        data_source_bins = data_source_bins,
        bin_selection = "first"
      ) %>%
      reduce_data_simple()

    # standardisation
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      )
    data_sd <-
      reduce_data_simple(
        data_source_reduce = data_sd
      )

    # transformation into proportions / percentages
    tranform_to_proportions <-
      TRUE
    data_sd_perc <-
      transform_into_proportions(
        data_source_trans = data_sd,
        sel_method = "percentages",
        verbose = FALSE
      )

    # Control
    data_com_control <-
      subset_community(data_sd)
    data_rowsums <-
      rowSums(data_com_control, na.rm = TRUE)

    # percentage
    data_com_control_perc <-
      data_com_control / data_rowsums * 100
    res_perc <-
      data_sd
    res_perc[, names(
      data_com_control_perc
    )] <-
      data_com_control_perc

    expect_identical(
      res_perc,
      data_sd_perc
    )
  }
)

test_that(
  "transform_into_proportions calculates mathematically correct proportions when sel_method='proportions'",
  {
    suppressWarnings(
      data_to_run_bins <-
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
        ) %>%
        smooth_community_data(
          smooth_method = "shep"
        ) %>%
        reduce_data(
          check_taxa = TRUE,
          check_levels = TRUE
        ) %>%
        prepare_data(
          data_source_prep = .,
          working_units = "bins",
          bin_size = 500,
          rand = 1
        ) %>%
        RUtilpol::flatten_list_by_one() %>%
        .[[1]]
    )

    data_source_subset <-
      data_to_run_bins$data
    data_source_bins <-
      data_to_run_bins$bins

    data_subset <-
      subset_samples(
        data_source_subset = data_source_subset,
        data_source_bins = data_source_bins,
        bin_selection = "first"
      ) %>%
      reduce_data_simple()

    # standardisation
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      )
    data_sd <-
      reduce_data_simple(
        data_source_reduce = data_sd
      )

    # transformation into proportions / percentages
    tranform_to_proportions <-
      TRUE
    data_sd_prop <-
      transform_into_proportions(
        data_source_trans = data_sd,
        sel_method = "proportions",
        verbose = FALSE
      )
    # Control
    data_com_control <-
      subset_community(data_sd)
    data_rowsums <-
      rowSums(data_com_control, na.rm = TRUE)

    # proportions
    data_com_control_prop <-
      data_com_control / data_rowsums * 1
    res_prop <-
      data_sd
    res_prop[, names(
      data_com_control_prop
    )] <-
      data_com_control_prop

    expect_identical(
      res_prop,
      data_sd_prop
    )
  }
)

# check output types and names
test_that(
  "transform_into_proportions preserves required column names (label, res_age, age_diff) when sel_method='proportions'",
  {
    suppressWarnings(
      data_to_run_bins <-
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
        ) %>%
        smooth_community_data(
          smooth_method = "shep"
        ) %>%
        reduce_data(
          check_taxa = TRUE,
          check_levels = TRUE
        ) %>%
        prepare_data(
          data_source_prep = .,
          working_units = "bins",
          bin_size = 500,
          rand = 1
        ) %>%
        RUtilpol::flatten_list_by_one() %>%
        .[[1]]
    )

    data_source_subset <-
      data_to_run_bins$data
    data_source_bins <-
      data_to_run_bins$bins

    data_subset <-
      subset_samples(
        data_source_subset = data_source_subset,
        data_source_bins = data_source_bins,
        bin_selection = "first"
      ) %>%
      reduce_data_simple()

    # standardisation
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      )
    data_sd <-
      reduce_data_simple(
        data_source_reduce = data_sd
      )

    # transformation into proportions / percentages
    tranform_to_proportions <-
      TRUE
    data_sd_prop <-
      transform_into_proportions(
        data_source_trans = data_sd,
        sel_method = "proportions",
        verbose = FALSE
      )
    expect_true(
      all(
        c("label", "res_age", "age_diff") %in% colnames(data_sd_prop)
      )
    )
  }
)

test_that(
  "transform_into_proportions preserves required column names (label, res_age, age_diff) when sel_method='percentages'",
  {
    suppressWarnings(
      data_to_run_bins <-
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
        ) %>%
        smooth_community_data(
          smooth_method = "shep"
        ) %>%
        reduce_data(
          check_taxa = TRUE,
          check_levels = TRUE
        ) %>%
        prepare_data(
          data_source_prep = .,
          working_units = "bins",
          bin_size = 500,
          rand = 1
        ) %>%
        RUtilpol::flatten_list_by_one() %>%
        .[[1]]
    )

    data_source_subset <-
      data_to_run_bins$data
    data_source_bins <-
      data_to_run_bins$bins

    data_subset <-
      subset_samples(
        data_source_subset = data_source_subset,
        data_source_bins = data_source_bins,
        bin_selection = "first"
      ) %>%
      reduce_data_simple()

    # standardisation
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      )
    data_sd <-
      reduce_data_simple(
        data_source_reduce = data_sd
      )

    # transformation into proportions / percentages
    tranform_to_proportions <-
      TRUE
    data_sd_prop <-
      transform_into_proportions(
        data_source_trans = data_sd,
        sel_method = "percentages",
        verbose = FALSE
      )
    expect_true(
      all(
        c("label", "res_age", "age_diff") %in% colnames(data_sd_prop)
      )
    )
  }
)

# s3 type
test_that(
  "transform_into_proportions returns data.frame S3 class when sel_method='proportions'",
  {
    suppressWarnings(
      data_to_run_bins <-
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
        ) %>%
        smooth_community_data(
          smooth_method = "shep"
        ) %>%
        reduce_data(
          check_taxa = TRUE,
          check_levels = TRUE
        ) %>%
        prepare_data(
          data_source_prep = .,
          working_units = "bins",
          bin_size = 500,
          rand = 1
        ) %>%
        RUtilpol::flatten_list_by_one() %>%
        .[[1]]
    )

    data_source_subset <-
      data_to_run_bins$data
    data_source_bins <-
      data_to_run_bins$bins

    data_subset <-
      subset_samples(
        data_source_subset = data_source_subset,
        data_source_bins = data_source_bins,
        bin_selection = "first"
      ) %>%
      reduce_data_simple()

    # standardisation
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      )
    data_sd <-
      reduce_data_simple(
        data_source_reduce = data_sd
      )

    # transformation into proportions / percentages
    tranform_to_proportions <-
      TRUE
    data_sd_prop <-
      transform_into_proportions(
        data_source_trans = data_sd,
        sel_method = "proportions",
        verbose = FALSE
      )
    expect_s3_class(
      data_sd_prop,
      "data.frame"
    )
  }
)

test_that(
  "transform_into_proportions returns data.frame S3 class when sel_method='percentages'",
  {
    suppressWarnings(
      data_to_run_bins <-
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
        ) %>%
        smooth_community_data(
          smooth_method = "shep"
        ) %>%
        reduce_data(
          check_taxa = TRUE,
          check_levels = TRUE
        ) %>%
        prepare_data(
          data_source_prep = .,
          working_units = "bins",
          bin_size = 500,
          rand = 1
        ) %>%
        RUtilpol::flatten_list_by_one() %>%
        .[[1]]
    )

    data_source_subset <-
      data_to_run_bins$data
    data_source_bins <-
      data_to_run_bins$bins

    data_subset <-
      subset_samples(
        data_source_subset = data_source_subset,
        data_source_bins = data_source_bins,
        bin_selection = "first"
      ) %>%
      reduce_data_simple()

    # standardisation
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      )
    data_sd <-
      reduce_data_simple(
        data_source_reduce = data_sd
      )

    # transformation into proportions / percentages
    tranform_to_proportions <-
      TRUE
    data_sd_prop <-
      transform_into_proportions(
        data_source_trans = data_sd,
        sel_method = "percentages",
        verbose = FALSE
      )
    expect_s3_class(
      data_sd_prop,
      "data.frame"
    )
  }
)

# Numeric values inside result
test_that(
  "transform_into_proportions returns numeric values in community columns when sel_method='proportions'",
  {
    suppressWarnings(
      data_to_run_bins <-
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
        ) %>%
        smooth_community_data(
          smooth_method = "shep"
        ) %>%
        reduce_data(
          check_taxa = TRUE,
          check_levels = TRUE
        ) %>%
        prepare_data(
          data_source_prep = .,
          working_units = "bins",
          bin_size = 500,
          rand = 1
        ) %>%
        RUtilpol::flatten_list_by_one() %>%
        .[[1]]
    )

    data_source_subset <-
      data_to_run_bins$data
    data_source_bins <-
      data_to_run_bins$bins

    data_subset <-
      subset_samples(
        data_source_subset = data_source_subset,
        data_source_bins = data_source_bins,
        bin_selection = "first"
      ) %>%
      reduce_data_simple()

    # standardisation
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      )
    data_sd <-
      reduce_data_simple(
        data_source_reduce = data_sd
      )

    # transformation into proportions / percentages
    tranform_to_proportions <-
      TRUE
    data_sd_prop <-
      transform_into_proportions(
        data_source_trans = data_sd,
        sel_method = "proportions",
        verbose = FALSE
      )

    expect_true(
      all(
        sapply(
          data_sd_prop[, -c(1:3)],
          is.numeric
        )
      )
    )
  }
)

test_that(
  "transform_into_proportions returns numeric values in community columns when sel_method='percentages'",
  {
    suppressWarnings(
      data_to_run_bins <-
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
        ) %>%
        smooth_community_data(
          smooth_method = "shep"
        ) %>%
        reduce_data(
          check_taxa = TRUE,
          check_levels = TRUE
        ) %>%
        prepare_data(
          data_source_prep = .,
          working_units = "bins",
          bin_size = 500,
          rand = 1
        ) %>%
        RUtilpol::flatten_list_by_one() %>%
        .[[1]]
    )

    data_source_subset <-
      data_to_run_bins$data
    data_source_bins <-
      data_to_run_bins$bins

    data_subset <-
      subset_samples(
        data_source_subset = data_source_subset,
        data_source_bins = data_source_bins,
        bin_selection = "first"
      ) %>%
      reduce_data_simple()

    # standardisation
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      )
    data_sd <-
      reduce_data_simple(
        data_source_reduce = data_sd
      )

    # transformation into proportions / percentages
    tranform_to_proportions <-
      TRUE
    data_sd_prop <-
      transform_into_proportions(
        data_source_trans = data_sd,
        sel_method = "percentages",
        verbose = FALSE
      )

    expect_true(
      all(
        sapply(
          data_sd_prop[, -c(1:3)],
          is.numeric
        )
      )
    )
  }
)
