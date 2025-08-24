# Setup workflow within estimate_roc and run_iteration
# Prepare default data (without smoothing)

# Input validation (Error messages)
test_that("transform_into_proportions throws error with missing data input for 'proportions'", {
  expect_error(
    transform_into_proportions(
      sel_method = "proportions"
    ),
    'argument "data_source_trans" is missing, with no default'
  )
})

test_that("transform_into_proportions throws error with missing data input for 'percentages'", {
  expect_error(
    transform_into_proportions(
      sel_method = "percentages"
    ),
    'argument "data_source_trans" is missing, with no default'
  )
})

# NULL
test_that("transform_into_proportions throws error with NULL data input for 'proportions'", {
  expect_error(
    transform_into_proportions(
      data_source_trans = NULL,
      sel_method = "proportions"
    ),
    "'data_source_trans' must be one of the following: 'data.frame'"
  )
})

test_that("transform_into_proportions throws error with NULL data input for 'percentages'", {
  expect_error(
    transform_into_proportions(
      data_source_trans = NULL,
      sel_method = "percentages"
    ),
    "'data_source_trans' must be one of the following: 'data.frame'"
  )
})

# NA
test_that("transform_into_proportions throws error with NA data input for 'proportions'", {
  expect_error(
    transform_into_proportions(
      data_source_trans = NA,
      sel_method = "proportions"
    ),
    "'data_source_trans' must be one of the following: 'data.frame'"
  )
})

test_that("transform_into_proportions throws error with NA data input for 'percentages'", {
  expect_error(
    transform_into_proportions(
      data_source_trans = NA,
      sel_method = "percentages"
    ),
    "'data_source_trans' must be one of the following: 'data.frame'"
  )
})

# Character
test_that("transform_into_proportions throws error with character data input for 'proportions'", {
  expect_error(
    transform_into_proportions(
      data_source_trans = "my_data",
      sel_method = "proportions"
    ),
    "'data_source_trans' must be one of the following: 'data.frame'"
  )
})

test_that("transform_into_proportions throws error with character data input for 'percentages'", {
  expect_error(
    transform_into_proportions(
      data_source_trans = "my_data",
      sel_method = "percentages"
    ),
    "'data_source_trans' must be one of the following: 'data.frame'"
  )
})

# Numeric
test_that("transform_into_proportions throws error with numeric data input for 'proportions'", {
  expect_error(
    transform_into_proportions(
      data_source_trans = 123,
      sel_method = "proportions"
    ),
    "'data_source_trans' must be one of the following: 'data.frame'"
  )
})

test_that("transform_into_proportions throws error with numeric data input for 'percentages'", {
  expect_error(
    transform_into_proportions(
      data_source_trans = 123,
      sel_method = "percentages"
    ),
    "'data_source_trans' must be one of the following: 'data.frame'"
  )
})

# list
test_that("transform_into_proportions throws error with list data input for 'proportions'", {
  expect_error(
    transform_into_proportions(
      data_source_trans = list(),
      sel_method = "proportions"
    ),
    "'data_source_trans' must be one of the following: 'data.frame'"
  )
})

test_that("transform_into_proportions throws error with list data input for 'percentages'", {
  expect_error(
    transform_into_proportions(
      data_source_trans = list(),
      sel_method = "percentages"
    ),
    "'data_source_trans' must be one of the following: 'data.frame'"
  )
})

# Matrix
test_that("transform_into_proportions throws error with matrix data input for 'proportions'", {
  expect_error(
    transform_into_proportions(
      data_source_trans = matrix(),
      sel_method = "proportions"
    ),
    "'data_source_trans' must be one of the following: 'data.frame'"
  )
})

test_that("transform_into_proportions throws error with matrix data input for 'percentages'", {
  expect_error(
    transform_into_proportions(
      data_source_trans = matrix(),
      sel_method = "percentages"
    ),
    "'data_source_trans' must be one of the following: 'data.frame'"
  )
})

# empty dataframe
test_that("transform_into_proportions throws error with empty dataframe for 'proportions'", {
  expect_error(
    transform_into_proportions(
      data_source_trans = data.frame(),
      sel_method = "proportions"
    ),
    # none programmed into the function. returns empty result data.frame
  )
})

test_that("transform_into_proportions throws error with empty dataframe for 'percentages'", {
  expect_error(
    transform_into_proportions(
      data_source_trans = data.frame(),
      sel_method = "percentages"
    ),
    # none programmed into the function. returns empty result data.frame
  )
})

# 0
test_that("transform_into_proportions throws error with 0 data input for 'proportions'", {
  expect_error(
    transform_into_proportions(
      data_source_trans = 0,
      sel_method = "proportions"
    ),
    "'data_source_trans' must be one of the following: 'data.frame'"
  )
})

test_that("transform_into_proportions throws error with 0 data input for 'percentages'", {
  expect_error(
    transform_into_proportions(
      data_source_trans = 0,
      sel_method = "percentages"
    ),
    "'data_source_trans' must be one of the following: 'data.frame'"
  )
})


# Input validation for sel_method
# Invalid character
test_that("transform_into_proportions throws error with invalid character supplied to sel_method", {
  data_to_run_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
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
        150
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
      n_individuals = 150
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
})

# Multiple methods
test_that("transform_into_proportions throws error with multiple supplied to sel_method", {
  data_to_run_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
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
        150
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
      n_individuals = 150
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
})

test_that("transform_into_proportions throws error with invalid character supplied to sel_method", {
  data_to_run_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
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
        150
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
      n_individuals = 150
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
})

# Numeric
test_that("transform_into_proportions throws error with numeric supplied to sel_method", {
  data_to_run_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
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
        150
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
      n_individuals = 150
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
})


# NULL
test_that("transform_into_proportions throws error with NULL supplied to sel_method", {
  data_to_run_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
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
        150
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
      n_individuals = 150
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
})

# Empty
test_that("transform_into_proportions throws warning if no input supplied to sel_method", {
  data_to_run_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
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
        150
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
      n_individuals = 150
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
})

# Output validation:
# Valid data:
test_that("transform_into_proportions returns correct output for 'proportions' with valid input", {
  data_to_run_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
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
        150
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
      n_individuals = 150
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
})

test_that("transform_into_proportions returns correct output for 'percentages' with valid input", {
  data_to_run_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
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
        150
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
      n_individuals = 150
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
})


# To Do:
# check mathematical accuracy of proportions and percentages
# minimal case scenarios with 1 row or 1 column and both
# check output types and names
# check: default parameters within run iteration output
