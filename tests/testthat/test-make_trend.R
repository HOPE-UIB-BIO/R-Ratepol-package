# 2 parameters:
# data_source
# sel_method = c("linear", "non_linear")

# 1. linear
# empty
test_that("make_trend errors when data_source is missing and sel_method is 'linear'", {
  expect_error(
    make_trend(
      data_source = ,
      sel_method = "linear"
    ),
    'argument "data_source" is missing, with no default'
  )
})

# NULL
test_that("make_trend errors when data_source is NULL and sel_method is 'linear'", {
  expect_error(
    make_trend(
      data_source = NULL,
      sel_method = "linear"
    ),
    "'data_source' must be one of the following: 'data.frame'"
  )
})

# character
test_that("make_trend errors when data_source is character and sel_method is 'linear'", {
  expect_error(
    make_trend(
      data_source = "invalid",
      sel_method = "linear"
    ),
    "'data_source' must be one of the following: 'data.frame'"
  )
})
# numeric
test_that("make_trend errors when data_source is numeric and sel_method is 'linear'", {
  expect_error(
    make_trend(
      data_source = 123,
      sel_method = "linear"
    ),
    "'data_source' must be one of the following: 'data.frame'"
  )
})

# empty list
test_that("make_trend errors when data_source is an empty list and sel_method is 'linear'", {
  expect_error(
    make_trend(
      data_source = list(),
      sel_method = "linear"
    ),
    "'data_source' must be one of the following: 'data.frame'"
  )
})
# empty data.frame
test_that("make_trend errors when data_source is an empty data.frame and sel_method is 'linear'", {
  expect_error(
    make_trend(
      data_source = data.frame(),
      sel_method = "linear"
    ),
    "'data_source' must contains following columns: 'ROC', 'Age'"
  )
})

# zero
test_that("make_trend errors when data_source is zero and sel_method is 'linear'", {
  expect_error(
    make_trend(
      data_source = 0,
      sel_method = "linear"
    ),
    "'data_source' must be one of the following: 'data.frame'"
  )
})

# NA
test_that("make_trend errors when data_source is NA and sel_method is 'linear'", {
  expect_error(
    make_trend(
      data_source = NA,
      sel_method = "linear"
    ),
    "'data_source' must be one of the following: 'data.frame'"
  )
})

# 2. Non-linear
# empty
test_that("make_trend errors when data_source is missing and sel_method is 'non_linear'", {
  expect_error(
    make_trend(
      data_source = ,
      sel_method = "non_linear"
    ),
    'argument "data_source" is missing, with no default'
  )
})

# NULL
test_that("make_trend errors when data_source is NULL and sel_method is 'non_linear'", {
  expect_error(
    make_trend(
      data_source = NULL,
      sel_method = "non_linear"
    ),
    "'data_source' must be one of the following: 'data.frame'"
  )
})

# character
test_that("make_trend errors when data_source is character and sel_method is 'non_linear'", {
  expect_error(
    make_trend(
      data_source = "invalid",
      sel_method = "non_linear"
    ),
    "'data_source' must be one of the following: 'data.frame'"
  )
})
# numeric
test_that("make_trend errors when data_source is numeric and sel_method is 'non_linear'", {
  expect_error(
    make_trend(
      data_source = 123,
      sel_method = "non_linear"
    ),
    "'data_source' must be one of the following: 'data.frame'"
  )
})

# empty list
test_that("make_trend errors when data_source is an empty list and sel_method is 'non_linear'", {
  expect_error(
    make_trend(
      data_source = list(),
      sel_method = "non_linear"
    ),
    "'data_source' must be one of the following: 'data.frame'"
  )
})
# empty data.frame
test_that("make_trend errors when data_source is an empty data.frame and sel_method is 'non_linear'", {
  expect_error(
    make_trend(
      data_source = data.frame(),
      sel_method = "non_linear"
    ),
    "'data_source' must contains following columns: 'ROC', 'Age'"
  )
})

# zero
test_that("make_trend errors when data_source is zero and sel_method is 'non_linear'", {
  expect_error(
    make_trend(
      data_source = 0,
      sel_method = "non_linear"
    ),
    "'data_source' must be one of the following: 'data.frame'"
  )
})

# NA
test_that("make_trend errors when data_source is NA and sel_method is 'non_linear'", {
  expect_error(
    make_trend(
      data_source = NA,
      sel_method = "non_linear"
    ),
    "'data_source' must be one of the following: 'data.frame'"
  )
})

# Valid data - invalid sel_method

# empty
test_that("make_trend uses default 'linear' when sel_method is missing and data_source is valid", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      # age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "shep",
      smooth_n_points = 5,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      silent = TRUE
    )

  expect_no_error(
    make_trend(
      data_source = data_source,
      sel_method = # silently uses 'linear' as default
      )
  )
})

# NULL
test_that("make_trend errors when sel_method is NULL and data_source is valid", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      # age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "shep",
      smooth_n_points = 5,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      silent = TRUE
    )

  expect_error(
    make_trend(
      data_source = data_source,
      sel_method = NULL
    ),
    "'sel_method' must be one of the following: 'character'"
  )
})

# Character
test_that("make_trend errors when sel_method is invalid character and data_source is valid", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      # age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "shep",
      smooth_n_points = 5,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      silent = TRUE
    )

  expect_error(
    make_trend(
      data_source = data_source,
      sel_method = "invalid"
    ),
    "'sel_method' must contains one of the following values: 'linear', 'non_linear'"
  )
})

# Numeric
test_that("make_trend errors when sel_method is numeric and data_source is valid", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      # age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "shep",
      smooth_n_points = 5,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      silent = TRUE
    )

  expect_error(
    make_trend(
      data_source = data_source,
      sel_method = 123
    ),
    "'sel_method' must be one of the following: 'character'"
  )
})

# Zero
test_that("make_trend errors when sel_method is zero and data_source is valid", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      # age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "shep",
      smooth_n_points = 5,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      silent = TRUE
    )

  expect_error(
    make_trend(
      data_source = data_source,
      sel_method = 0
    ),
    "'sel_method' must be one of the following: 'character'"
  )
})

# NA

test_that("make_trend errors when sel_method is NA and data_source is valid", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      # age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "shep",
      smooth_n_points = 5,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      silent = TRUE
    )

  expect_error(
    make_trend(
      data_source = data_source,
      sel_method = NA
    ),
    "'sel_method' must be one of the following: 'character'"
  )
})

# multiple

test_that("make_trend errors when sel_method is length > 1 and data_source is valid", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      # age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "shep",
      smooth_n_points = 5,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      silent = TRUE
    )

  expect_error(
    make_trend(
      data_source = data_source,
      sel_method = c("linear", "non_linear")
    ),
    # not programmed yet
    "'sel_method' must be of length 1"
  )
})

# empty list
test_that("make_trend errors when sel_method is an empty list and data_source is valid", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      # age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "shep",
      smooth_n_points = 5,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      silent = TRUE
    )

  expect_error(
    make_trend(
      data_source = data_source,
      sel_method = list()
    ),
    "'sel_method' must be one of the following: 'character'"
  )
})

# empty dataframe

test_that("make_trend errors when sel_method is an empty data.frame and data_source is valid", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      # age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "shep",
      smooth_n_points = 5,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      silent = TRUE
    )

  expect_error(
    make_trend(
      data_source = data_source,
      sel_method = data.frame()
    ),
    "'sel_method' must be one of the following: 'character'"
  )
})

# Functionality
# test that both methods produce different results
test_that("make_trend produces different results for 'linear' and 'non_linear' sel_method", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      # age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "shep",
      smooth_n_points = 5,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      silent = TRUE
    )

  res_linear <-
    make_trend(
      data_source = data_source,
      sel_method = "linear"
    )

  res_non_linear <-
    make_trend(
      data_source = data_source,
      sel_method = "non_linear"
    )

  expect_false(
    all(
      res_linear == res_non_linear
    )
  )
})


# Output validation
test_that("make_trend returns a numeric vector when data_source is valid and sel_method is 'linear'", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      # age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "shep",
      smooth_n_points = 5,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      silent = TRUE
    )

  res <-
    make_trend(
      data_source = data_source,
      sel_method = "linear"
    )

  expect_true(
    is.numeric(res) & is.vector(res)
  )
})

test_that("make_trend returns a numeric vector when data_source is valid and sel_method is 'non_linear'", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      # age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "shep",
      smooth_n_points = 5,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      silent = TRUE
    )

  res <-
    make_trend(
      data_source = data_source,
      sel_method = "non_linear"
    )

  expect_true(
    is.numeric(res) & is.vector(res) # is array instead.
  )
})


## Additional tests for edge cases / functionality:
# Test for data frame with NA values in ROC column
test_that("make_trend with 'linear' handles NA values in ROC column appropriately (predict NA))", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      # age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "shep",
      smooth_n_points = 5,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      silent = TRUE
    )

  data_source$ROC[2] <- NA # Introduce NA value in ROC column

  expect_false(
    any(
      is.na(
        make_trend(
          data_source,
          sel_method = "linear"
        )
      )
    )
  )
})

test_that("make_trend with 'non_linear' handles NA values in ROC column appropriately (predict NA))", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      # age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "shep",
      smooth_n_points = 5,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      silent = TRUE
    )

  data_source$ROC[2] <- NA # Introduce NA value in ROC column

  expect_false(
    any(
      is.na(
        make_trend(
          data_source,
          sel_method = "non_linear"
        )
      )
    )
  )
})

# Test for data frame with NA values in Age column
test_that("make_trend with 'linear' handles NA values in Age column appropriately", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      # age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "shep",
      smooth_n_points = 5,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      silent = TRUE
    )

  data_source$Age[2] <- NA # Introduce NA value in Age column

  expect_false(
    any(
      is.na(
        make_trend(
          data_source,
          sel_method = "linear"
        )
      )
    )
  )
})

test_that("make_trend with 'non_linear' handles NA values in Age column appropriately", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      # age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "shep",
      smooth_n_points = 5,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      silent = TRUE
    )

  data_source$Age[2] <- NA # Introduce NA value in Age column

  expect_false(
    any(
      is.na(
        make_trend(
          data_source,
          sel_method = "non_linear"
        )
      )
    )
  )
})

# Test for data frame with only one row (edge case)
test_that("make_trend with 'linear' throws warning with single row data frame", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      # age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "shep",
      smooth_n_points = 5,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      silent = TRUE
    )
  data_source <-
    data_source[1, ] # Keep only the first row

  expect_warning(
    res <-
      make_trend(
        data_source,
        sel_method = "linear"
      ),
    "NaNs produced"
  )
})

test_that("make_trend with 'non_linear' throws error with single row data frame", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      # age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "shep",
      smooth_n_points = 5,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      silent = TRUE
    )
  data_source <-
    data_source[1, ] # Keep only the first row
  # Non-linear might need more data points
  expect_error(
    make_trend(
      data_source,
      sel_method = "non_linear"
    ),
    regexp = "Not enough \\(non-NA\\) data to do anything meaningful"
  )
})


# Test for data frame with negative ROC values (edge case)
test_that("make_trend with 'linear' rejects negative ROC values", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      # age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "shep",
      smooth_n_points = 5,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      silent = TRUE
    )

  data_source$ROC <-
    data_source$ROC * -1 # Introduce negative values in ROC column

  expect_warning(
    expect_error(
      make_trend(
        data_source,
        sel_method = "linear"
      ),
      "missing value where TRUE/FALSE needed"
    ),
    "NaNs produced"
  )
})

test_that("make_trend with 'non_linear' rejects negative ROC values", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      # age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "shep",
      smooth_n_points = 5,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      silent = TRUE
    )

  data_source$ROC <-
    data_source$ROC * -1 # Introduce negative values in ROC column

  expect_warning(
    expect_error(
      make_trend(
        data_source,
        sel_method = "non_linear"
      ),
      "missing value where TRUE/FALSE needed"
    ),
    "NaNs produced"
  )
})

# Test that output length matches input data frame rows
test_that("make_trend output length matches input rows with 'linear'", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      # age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "shep",
      smooth_n_points = 5,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      silent = TRUE
    )

  res <-
    make_trend(
      data_source,
      sel_method = "linear"
    )

  expect_equal(
    length(res),
    nrow(data_source)
  )
})

test_that("make_trend output length matches input rows with 'non_linear'", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      # age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "shep",
      smooth_n_points = 5,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      silent = TRUE
    )

  res <-
    make_trend(
      data_source,
      sel_method = "non_linear"
    )

  expect_equal(
    length(res),
    nrow(data_source)
  )
})
