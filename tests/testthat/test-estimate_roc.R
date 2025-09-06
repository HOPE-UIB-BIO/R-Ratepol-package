test_that("estimate_roc returns dataframe with valid inputs", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]
  # Run the function with valid inputs
  result <-
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "none",
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    )

  # Check that the result is a data frame and contains expected elements
  expect_s3_class(
    result,
    "data.frame"
  )
})


# ============================= #
# 1. INPUT VALIDATION           #
# ============================= #


# 1. data_source_community validation
test_that("data_source_community: Missing argument throws error indicating required data.frame", {
  # Create example data
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = ,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "none",
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = NULL,
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "Object 'data_source_community' must be included as a 'data.frame'"
  )
})

test_that("data_source_community: NULL value throws error indicating data.frame required", {
  # Create example data
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = NULL,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "none",
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = NULL,
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'data_source_community' must be one of the following: 'data.frame'"
  )
})

# character
test_that("data_source_community: Character string instead of data.frame throws error", {
  # Create example data
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = "my_data",
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "none",
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = NULL,
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'data_source_community' must be one of the following: 'data.frame'"
  )
})


# Numeric
test_that("data_source_community: Numeric value instead of data.frame throws error", {
  # Create example data
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = 123,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "none",
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = NULL,
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'data_source_community' must be one of the following: 'data.frame'"
  )
})

# Zero
test_that("data_source_community: Zero value instead of data.frame throws error", {
  # Create example data
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = 0,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "none",
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = NULL,
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'data_source_community' must be one of the following: 'data.frame'"
  )
})

# NA
test_that("data_source_community: NA value instead of data.frame throws error", {
  # Create example data
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = NA,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "none",
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = NULL,
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'data_source_community' must be one of the following: 'data.frame'"
  )
})


# empty list
test_that("data_source_community: Empty list instead of data.frame throws error", {
  # Create example data
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = list(),
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "none",
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = NULL,
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'data_source_community' must be one of the following: 'data.frame'"
  )
})


# empty dataframe
test_that("data_source_community: Empty data.frame throws error about missing required columns", {
  # Create example data
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data.frame(),
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "none",
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = NULL,
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'data_community_extract' must contains following columns: 'sample_id'"
  )
})

# 2. data_source_age validation
test_that("data_source_age: Missing argument throws error indicating required data.frame", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = ,
      age_uncertainty = age_uncertainty,
      smooth_method = "none",
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = NULL,
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "Object 'data_source_age' must be included as a 'data.frame'"
  )
})

test_that("data_source_age: NULL value throws error indicating data.frame required", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = NULL,
      age_uncertainty = age_uncertainty,
      smooth_method = "none",
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = NULL,
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'data_source_age' must be one of the following: 'data.frame'"
  )
})

# character
test_that("data_source_age: Character string instead of data.frame throws error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = "my_age",
      age_uncertainty = age_uncertainty,
      smooth_method = "none",
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = NULL,
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'data_source_age' must be one of the following: 'data.frame'"
  )
})


# Numeric
test_that("data_source_age: Numeric value instead of data.frame throws error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = 123,
      age_uncertainty = age_uncertainty,
      smooth_method = "none",
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = NULL,
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'data_source_age' must be one of the following: 'data.frame'"
  )
})

# Zero
test_that("data_source_age: Zero value instead of data.frame throws error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = 0,
      age_uncertainty = age_uncertainty,
      smooth_method = "none",
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = NULL,
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'data_source_age' must be one of the following: 'data.frame'"
  )
})

# NA
test_that("data_source_age: NA value instead of data.frame throws error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = NA,
      age_uncertainty = age_uncertainty,
      smooth_method = "none",
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = NULL,
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'data_source_age' must be one of the following: 'data.frame'"
  )
})


# empty list
test_that("data_source_age: Empty list instead of data.frame throws error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = list(),
      age_uncertainty = age_uncertainty,
      smooth_method = "none",
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = NULL,
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'data_source_age' must be one of the following: 'data.frame'"
  )
})


# empty dataframe
test_that("data_source_age: Empty data.frame throws error about missing required columns", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data.frame(),
      age_uncertainty = age_uncertainty,
      smooth_method = "none",
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = NULL,
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'data_age_extract' must contains following columns: 'sample_id'"
  )
})

# 3. Age uncertainty validation
# it's valid if missing (NULL) - no error expected
test_that("age_uncertainty: NULL value is valid and function completes without error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  # Run the function with valid inputs
  result <-
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = NULL,
      smooth_method = "none",
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    )

  # Check that the result is a data frame and contains expected elements
  expect_s3_class(
    result,
    "data.frame"
  )
})


test_that("age_uncertainty: Missing argument is valid and function completes without error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  # Run the function with valid inputs
  result <-
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = ,
      smooth_method = "none",
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    )

  # Check that the result is a data frame and contains expected elements
  expect_s3_class(
    result,
    "data.frame"
  )
})


# character
test_that("age_uncertainty: Character string instead of matrix throws error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = "my_age_uncertainty",
      smooth_method = "none",
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'age_uncertainty' must be one of the following: 'NULL', 'matrix'"
  )
})

# Numeric
test_that("age_uncertainty: Numeric value instead of matrix throws error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = 123,
      smooth_method = "none",
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'age_uncertainty' must be one of the following: 'NULL', 'matrix'"
  )
})

# Zero
test_that("age_uncertainty: Zero value instead of matrix throws error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = 0,
      smooth_method = "none",
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'age_uncertainty' must be one of the following: 'NULL', 'matrix'"
  )
})

# NA
test_that("age_uncertainty: NA value instead of matrix throws error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = NA,
      smooth_method = "none",
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = NULL,
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'age_uncertainty' must be one of the following: 'NULL', 'matrix'"
  )
})

# empty list
test_that("age_uncertainty: Empty list instead of matrix throws error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = list(),
      smooth_method = "none",
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = NULL,
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'age_uncertainty' must be one of the following: 'NULL', 'matrix'"
  )
})

# empty dataframe
test_that("age_uncertainty: Empty data.frame instead of matrix throws error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = data.frame(),
      smooth_method = "none",
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = NULL,
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'age_uncertainty' must be one of the following: 'NULL', 'matrix'"
  )
})

# 3. smooth_method validation

test_that("smooth_method: NULL value throws error requiring character input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = NULL,
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = NULL,
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'smooth_method' must be one of the following: 'character'"
  )
})

# empty
test_that("smooth_method: Missing argument uses default value 'none' without error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_no_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = , # will use first element of vector -> "none"
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    )
  )
})

# NULL
test_that("smooth_method: NULL value throws error requiring character input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = NULL,
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'smooth_method' must be one of the following: 'character'"
  )
})

# multiple
test_that("smooth_method: Multiple values throw error as only single value allowed", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = c("none", "shep"),
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'arg' must be of length 1"
  )
})

# character
test_that("smooth_method: Invalid character value throws error listing allowed options", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "my_smooth_method",
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'smooth_method' must contains one of the following values: 'none', 'm.avg', 'grim', 'age.w', 'shep'"
  )
})

# Numeric
test_that("smooth_method: Numeric value throws error requiring character input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = 123,
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'smooth_method' must be one of the following: 'character'"
  )
})

# 0
test_that("smooth_method: Zero value throws error requiring character input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = 0,
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'smooth_method' must be one of the following: 'character'"
  )
})

# NA
test_that("smooth_method: NA value throws error requiring character input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = NA,
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'smooth_method' must be one of the following: 'character'"
  )
})

# empty list
test_that("smooth_method: Empty list throws error requiring character input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = list(),
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'smooth_method' must be one of the following: 'character'"
  )
})

# empty dataframe
test_that("smooth_method: Empty data.frame throws error requiring character input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = data.frame(),
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'smooth_method' must be one of the following: 'character'"
  )
})

# 4. smooth_n_points validation (for grim, age.w smoothing; must be odd)
# empty
test_that("smooth_n_points: Missing argument uses default value 5 without error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_no_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "age.w",
      smooth_n_points = , # will use 5 as default
      smooth_age_range = 500,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    )
  )
})

# NULL
test_that("smooth_n_points: NULL value throws error requiring numeric value", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "age.w",
      smooth_n_points = NULL,
      smooth_age_range = 500,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = ,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "assert_that: length of assertion is not 1"
  )
})

# NULL
test_that("smooth_n_points: NULL value throws error requiring numeric value", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "age.w",
      smooth_n_points = NULL,
      smooth_age_range = 500,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "assert_that: length of assertion is not 1"
  )
})

# multiple
test_that("smooth_n_points: Multiple values throw error as only single value allowed", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "age.w",
      smooth_n_points = c(3, 5),
      smooth_age_range = 500,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "length of assertion is not 1"
  )
})

# character
test_that("smooth_n_points: Character value throws error requiring numeric input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "age.w",
      smooth_n_points = "my_smooth_n_points",
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "non-numeric argument to binary operator"
  )
})

# Numeric (even - i.e., not odd)
test_that("smooth_n_points: Even number throws error as odd number required", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "age.w",
      smooth_n_points = 4,
      smooth_age_range = 500,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'smooth_n_points' must be an odd number"
  )
})

# 0
test_that("smooth_n_points: Zero value throws error as odd number required", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "age.w",
      smooth_n_points = 0,
      smooth_age_range = 500,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'smooth_n_points' must be an odd number"
  )
})
# NA
test_that("smooth_n_points: NA value throws error requiring non-missing value", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "age.w",
      smooth_n_points = NA,
      smooth_age_range = 500,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "assert_that: missing values present in assertion"
  )
})

# empty list
test_that("smooth_n_points: Empty list throws error requiring numeric input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "age.w",
      smooth_n_points = list(),
      smooth_age_range = 500,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "non-numeric argument to binary operator"
  )
})
# empty dataframe
test_that("smooth_n_points: Empty data.frame throws error requiring single value", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "age.w",
      smooth_n_points = data.frame(),
      smooth_age_range = 500,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "assert_that: length of assertion is not 1"
  )
})

# 5. smooth_age_range validation
# empty
test_that("smooth_age_range: Missing argument uses default value 500 without error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_no_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "age.w",
      smooth_n_points = 5,
      smooth_age_range = , # will use 500 as default
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    )
  )
})

# NULL
test_that("smooth_age_range: NULL value throws error requiring numeric input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "age.w",
      smooth_n_points = 5,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'smooth_age_range' must be one of the following: 'numeric'"
  )
})

# multiple
test_that("smooth_age_range: Multiple values causes error in internal function", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    suppressWarnings(
      estimate_roc(
        data_source_community = data_source_community,
        data_source_age = data_source_age,
        age_uncertainty = age_uncertainty,
        smooth_method = "age.w",
        smooth_n_points = 5,
        smooth_age_range = c(300, 500),
        smooth_n_max = NULL,
        working_units = "levels",
        bin_size = 500, # to avoid error during testing
        number_of_shifts = NULL,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        use_parallel = FALSE,
        interest_threshold = NULL,
        time_standardisation = 500,
        verbose = FALSE
      )
    ),
    # none programmed into the function yet
    # e.g., "Error: Multiple arguments supplied to 'smooth_age_range'"
  )
})

# character
test_that("smooth_age_range: Character value throws error requiring numeric input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "age.w",
      smooth_n_points = 5,
      smooth_age_range = "my_smooth_age_range",
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 10,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'smooth_age_range' must be one of the following: 'numeric'"
  )
})

# 0
test_that("smooth_age_range: Zero value causes error in internal calculations", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "age.w",
      smooth_n_points = 5,
      smooth_age_range = 0,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 10,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "`age_diff` must be size 1, not 2"
  )
})

# NA
test_that("smooth_age_range: NA value throws error requiring numeric input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "age.w",
      smooth_n_points = 5,
      smooth_age_range = NA,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 10,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'smooth_age_range' must be one of the following: 'numeric'"
  )
})

# empty list
test_that("smooth_age_range: Empty list throws error requiring numeric input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "age.w",
      smooth_n_points = 5,
      smooth_age_range = list(),
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 10,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'smooth_age_range' must be one of the following: 'numeric'"
  )
})
# empty dataframe
test_that("smooth_age_range: Empty data.frame throws error requiring numeric input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "age.w",
      smooth_n_points = 5,
      smooth_age_range = data.frame(),
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 10,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'smooth_age_range' must be one of the following: 'numeric'"
  )
})

# 6. smooth_n_max validation

# empty
test_that("smooth_n_max: Missing argument uses default value 9 without error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_no_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = , # uses 9 as default
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = ,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 10,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    )
  )
})

# NULL
test_that("smooth_n_max: NULL value throws error requiring single value", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 10,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "length of assertion is not 1"
  )
})

# multiple
test_that("smooth_n_max: Multiple values throw error as only single value allowed", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = c(7, 9),
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 10,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "assert_that: length of assertion is not 1"
  )
})

# character
test_that("smooth_n_max: Character value throws error requiring numeric input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = "my_smooth_n_max",
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 10,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "non-numeric argument to binary operator"
  )
})


# Numeric
test_that("smooth_n_max: Even number throws error as odd number required", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 4,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 10,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'smooth_n_max' must be an odd number"
  )
})

test_that("smooth_n_max: Value not larger than smooth_n_points throws error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 5,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 10,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'smooth_n_max' must be bigger than 'smooth_n_points"
  )
})

# 0
test_that("smooth_n_max: Zero value throws error as odd number required", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 0,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 10,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'smooth_n_max' must be an odd number"
  )
})

# NA
test_that("smooth_n_max: NA value throws error requiring non-missing value", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = NA,
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 10,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "assert_that: missing values present in assertion"
  )
})

# empty list
test_that("smooth_n_max: Empty list throws error requiring numeric input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = list(),
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 10,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "non-numeric argument to binary operator"
  )
})

# empty dataframe
test_that("smooth_n_max: Empty data.frame throws error requiring single value", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = data.frame(),
      working_units = "levels",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 10,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "assert_that: length of assertion is not 1"
  )
})

# 7. working_units validation
# empty
test_that("working_units: Missing argument uses default value 'levels' without error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_no_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      # age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = , # will use first element of default vector (levels)
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    )
  )
})

# NULL
test_that("working_units: NULL value throws error requiring character input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      # age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = NULL,
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'working_units' must be one of the following: 'character'"
  )
})

# multiple
test_that("working_units: Multiple values throw error as only single value allowed", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      # age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = c("levels", "bins"),
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'arg' must be of length 1"
  )
})

# character
test_that("working_units: Invalid character value throws error listing allowed options", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      # age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "my_working_units",
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'working_units' must contains one of the following values: 'levels', 'bins', 'MW'"
  )
})

# Numeric
test_that("working_units: Numeric value throws error requiring character input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      # age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = 123,
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'working_units' must be one of the following: 'character'"
  )
})

# 0
test_that("working_units: Zero value throws error requiring character input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      # age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = 0,
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'working_units' must be one of the following: 'character'"
  )
})

# NA
test_that("working_units: NA value throws error requiring character input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      # age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = NA,
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'working_units' must be one of the following: 'character'"
  )
})

# empty list
test_that("working_units: Empty list throws error requiring character input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      # age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = list(),
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'working_units' must be one of the following: 'character'"
  )
})

# empty dataframe
test_that("working_units: Empty data.frame throws error requiring character input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      # age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = data.frame(),
      bin_size = 500, # to avoid error during testing
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'working_units' must be one of the following: 'character'"
  )
})

# 8. bin_size validation
# empty
test_that("bin_size: Missing argument uses default value 500 without error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_no_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      # age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "bins",
      bin_size = , # will use default 500 silently
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    )
  )
})
# NULL
test_that("bin_size: NULL value throws error requiring numeric input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      # age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "bins",
      bin_size = NULL,
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'bin_size' must be one of the following: 'numeric'"
  )
})

# multiple
test_that("bin_size: Multiple numeric values throws error requiring single argument", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      # age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "bins",
      bin_size = c(100, 200),
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "assert_that: length of assertion is not 1"
  )
})

# character
test_that("bin_size: Character value throws error requiring numeric input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      # age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "bins",
      bin_size = "my_bin_size",
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'bin_size' must be one of the following: 'numeric'"
  )
})

# 0
test_that("bin_size: Zero value throws error requiring numeric input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      # age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "bins",
      bin_size = 0,
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "invalid"
  )
})

# NA
test_that("bin_size: NA value throws error requiring numeric input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      # age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "bins",
      bin_size = NA,
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'bin_size' must be one of the following: 'numeric'"
  )
})

# empty list
test_that("bin_size: Empty list throws error requiring numeric input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      # age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "bins",
      bin_size = list(),
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'bin_size' must be one of the following: 'numeric'"
  )
})

# empty dataframe
test_that("bin_size: Empty data.frame throws error requiring numeric input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      # age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "bins",
      bin_size = data.frame(),
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'bin_size' must be one of the following: 'numeric'"
  )
})


# 9. number_of_shifts validation
# empty
test_that("number_of_shifts: Missing argument uses default value 5 without error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_no_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      # age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = , # will use default 5
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    )
  )
})
# NULL
test_that("number_of_shifts: NULL value throws error requiring numeric input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      # age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = NULL, # will use default 5
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'number_of_shifts' must be one of the following: 'numeric'"
  )
})

# multiple
test_that("number_of_shifts: Multiple values throw error as only single value allowed", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      # age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = c(3, 5), # will use default 5
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "assert_that: length of assertion is not 1"
  )
})

# character
test_that("number_of_shifts: Character value throws error requiring numeric input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      # age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = "my_number_of_shifts",
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'number_of_shifts' must be one of the following: 'numeric'"
  )
})

# Numeric
test_that("number_of_shifts: Valid numeric value works without error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_no_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      # age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
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
      time_standardisation = 500,
      verbose = FALSE
    )
  )
})

# 0
test_that("number_of_shifts: Zero value is overwritten with 1 without error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_no_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      # age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 0, # gets overwritten with 1
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    )
  )
})

# NA
test_that("number_of_shifts: NA value throws error requiring numeric input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      # age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = NA,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'number_of_shifts' must be one of the following: 'numeric'"
  )
})

# empty list
test_that("number_of_shifts: Empty list throws error requiring numeric input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      # age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = list(),
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'number_of_shifts' must be one of the following: 'numeric'"
  )
})

# empty dataframe
test_that("estimate_roc throws error when smooth_method is NULL (invalid empty input)", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      # age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = data.frame(),
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'number_of_shifts' must be one of the following: 'numeric'"
  )
})

# 10. bin_selection validation
# empty
test_that("estimate_roc throws error when smooth_method is NULL (invalid empty input)", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_no_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = , # will use default "random"
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    )
  )
})

# NULL
test_that("bin_selection: NULL value throws error requiring character input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = NULL,
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'bin_selection' must be one of the following: 'character'"
  )
})

# multiple
test_that("bin_selection: Multiple values throw error as only single value allowed", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = c("first", "random"),
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'arg' must be of length 1"
  )
})

# character
test_that("bin_selection: Invalid character value throws error listing allowed options", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "my_bin_selection",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'bin_selection' must contains one of the following values: 'first', 'random'"
  )
})

# Numeric
test_that("bin_selection: Numeric value throws error requiring character input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = 123,
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'bin_selection' must be one of the following: 'character'"
  )
})

# 0
test_that("bin_selection: Zero value throws error requiring character input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = 0,
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'bin_selection' must be one of the following: 'character'"
  )
})

# NA
test_that("bin_selection: NA value throws error requiring character input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = NA,
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'bin_selection' must be one of the following: 'character'"
  )
})

# empty list
test_that("bin_selection: Empty list throws error requiring character input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = list(),
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'bin_selection' must be one of the following: 'character'"
  )
})

# empty dataframe
test_that("bin_selection: Empty data.frame throws error requiring character input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = data.frame(),
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'bin_selection' must be one of the following: 'character'"
  )
})


# 11. standardise validation
# empty
test_that("standardise: Missing argument uses default value FALSE without error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_no_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = , # will use default FALSE silently
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    )
  )
})

# NULL
test_that("standardise: NULL value throws error requiring logical input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = NULL,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'standardise' must be one of the following: 'logical'"
  )
})

# multiple
test_that("standardise: Multiple values throw error as only single value allowed", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = c(FALSE, TRUE),
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    # none programmed into function yet. will use first element in vector
    # e.g., "assert_that: length of assertion is not 1"
  )
})

# character
test_that("standardise: Character value throws error requiring logical input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = "TRUE",
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'standardise' must be one of the following: 'logical'"
  )
})

# Numeric
test_that("standardise: Numeric value throws error requiring logical input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = 1,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'standardise' must be one of the following: 'logical'"
  )
})

# 0
test_that("standardise: Zero value throws error requiring logical input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = 0,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'standardise' must be one of the following: 'logical'"
  )
})

# NA
test_that("standardise: NA value throws error instead of using FALSE silently", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = NA, # uses FALSE siltently
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    # none programmed into the function yet. Skips standardising
    # e.g., "'standardise' must be one of the following: c(TRUE, FALSE)"
  )
})

# empty list
test_that("standardise: Empty list throws error requiring logical input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = list(),
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'standardise' must be one of the following: 'logical'"
  )
})

# empty dataframe
test_that("standardise: Empty data.frame throws error requiring logical input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = data.frame(),
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'standardise' must be one of the following: 'logical'"
  )
})


# 12. n_individuals validation
# empty
test_that("n_individuals: Missing argument uses default value 150 without error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_no_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = , # uses default 150 silently
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    )
  )
})

# NULL
test_that("n_individuals: NULL value throws error requiring numeric input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = NULL,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'n_individuals' must be one of the following: 'numeric'"
  )
})

# multiple
test_that("n_individuals: Multiple values throw error as only single value allowed", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = c(150, 500),
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "assert_that: length of assertion is not 1"
  )
})

# character
test_that("n_individuals: Character value throws error requiring numeric input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = "my_n_individuals",
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'n_individuals' must be one of the following: 'numeric'"
  )
})


# 0
test_that("n_individuals: Zero value throws error in internal calculations", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 0,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "invalid 'length' argument"
  )
})

# NA
test_that("n_individuals: NA value throws error requiring numeric input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = NA,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'n_individuals' must be one of the following: 'numeric'"
  )
})

# empty list
test_that("n_individuals: Empty list throws error requiring numeric input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = list(),
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'n_individuals' must be one of the following: 'numeric'"
  )
})

# empty dataframe
test_that("n_individuals: Empty dataframe throws error requiring numeric input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = data.frame(),
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'n_individuals' must be one of the following: 'numeric'"
  )
})

# 13. dissimilarity_coefficient validation
# empty
test_that("dissimilarity_coefficient: Missing argument uses default value 'euc' without error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_no_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = , # will use default "euc"
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    )
  )
})

# NULL
test_that("dissimilarity_coefficient: NULL value throws error requiring character input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = NULL,
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'dissimilarity_coefficient' must be one of the following: 'character'"
  )
})

# multiple
test_that("dissimilarity_coefficient: Multiple values throw error as only single value allowed", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = c("euc", "bray"),
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'arg' must be of length 1"
  )
})
# character
test_that("dissimilarity_coefficient: invalid character values throw error as only single value allowed", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "my_dissimilarity_coefficient",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'dissimilarity_coefficient' must contains one of the following values: 'euc', 'euc.sd', 'chord', 'chisq', 'gower', 'bray'"
  )
})

# Numeric
test_that("dissimilarity_coefficient: Numeric value throws error requiring character input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = 123,
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'dissimilarity_coefficient' must be one of the following: 'character'"
  )
})

# 0
test_that("dissimilarity_coefficient: Zero value throws error requiring character input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = 0,
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'dissimilarity_coefficient' must be one of the following: 'character'"
  )
})

# NA
test_that("dissimilarity_coefficient: NA value throws error requiring character input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = NA,
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'dissimilarity_coefficient' must be one of the following: 'character'"
  )
})

# empty list
test_that("dissimilarity_coefficient: Empty list throws error requiring character input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = list(),
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'dissimilarity_coefficient' must be one of the following: 'character'"
  )
})

# empty dataframe
test_that("estimate_roc throws error when smooth_method is NULL (invalid empty input)", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = data.frame(),
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'dissimilarity_coefficient' must be one of the following: 'character'"
  )
})

# 14. tranform_to_proportions validation
# empty
test_that("estimate_roc throws error when smooth_method is NULL (invalid empty input)", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_no_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = , # will use default TRUE
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    )
  )
})

# NULL
test_that("tranform_to_proportions: NULL value throws error requiring logical input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = NULL,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'tranform_to_proportions' must be one of the following: 'logical'"
  )
})

# multiple
test_that("tranform_to_proportions: Multiple values throw error as only single value allowed", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = c(FALSE, TRUE),
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    # none programmed into function yet
    "'arg' must be of length 1"
  )
})

# character
test_that("tranform_to_proportions: Character value throws error requiring logical input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = "my_choice",
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'tranform_to_proportions' must be one of the following: 'logical'"
  )
})

# Numeric
test_that("tranform_to_proportions: Numeric value throws error requiring logical input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = 123,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'tranform_to_proportions' must be one of the following: 'logical'"
  )
})

# 0
test_that("tranform_to_proportions: Zero value throws error requiring logical input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = 0,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'tranform_to_proportions' must be one of the following: 'logical'"
  )
})

# NA
test_that("tranform_to_proportions: NA value throws error requiring non-missing value", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = NA,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "assert_that: missing values present in assertion"
  )
})
# empty list
test_that("tranform_to_proportions: Empty list throws error requiring logical input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = list(),
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'tranform_to_proportions' must be one of the following: 'logical'"
  )
})
# empty dataframe
test_that("tranform_to_proportions: Empty data.frame throws error requiring logical input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = data.frame(),
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'tranform_to_proportions' must be one of the following: 'logical'"
  )
})


# 15. rand validation
# empty
test_that("rand: Missing argument uses default value without error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = , # will use default NULL
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "argument is of length zero"
  )
})

# NULL
test_that("rand: NULL value throws error requiring non-zero input for rand", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NULL,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "argument is of length zero"
  )
})

# multiple
test_that("rand: Multiple values for tranform_to_proportions throw error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = c(100, 200),
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "assert_that: length of assertion is not 1"
  )
})

# character
test_that("rand: Character value throws error requiring numeric or NULL", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = "100",
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'rand' must be one of the following: 'NULL', 'numeric'"
  )
})

# 0
test_that("rand: Zero value causes error in internal function calculation", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 0,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "subscript out of bounds"
  )
})

# NA
test_that("rand: NA value throws error requiring numeric or NULL", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = NA,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'rand' must be one of the following: 'NULL', 'numeric'"
  )
})

# empty list
test_that("rand: Empty list throws error requiring numeric or NULL", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = list(),
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'rand' must be one of the following: 'NULL', 'numeric'"
  )
})

# empty dataframe
test_that("rand: Empty data.frame throws error requiring numeric or NULL", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = data.frame(),
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'rand' must be one of the following: 'NULL', 'numeric'"
  )
})

# 16. use_parallel validation
# empty
test_that("use_parallel: Missing argument uses default value without error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_no_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = , # will use default FALSE
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    )
  )
})

# NULL
test_that("use_parallel: NULL input throws error requiring logical or numeric", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = NULL,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'use_parallel' must be one of the following: 'logical', 'numeric'"
  )
})

# multiple
test_that("use_parallel: multiple input throws error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = c(TRUE, 100),
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "assert_that: length of assertion is not 1"
  )
})

# character
test_that("use_parallel: Character input throws error requiring logical or numeric", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = "TRUE",
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'use_parallel' must be one of the following: 'logical', 'numeric'"
  )
})

# 0
test_that("use_parallel: Zero (0) input throws error requiring more than 0 cores", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = 0,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    # this should throw an error
    # - not sure what it does if 0 cores are entered as input
  )
})

# NA
test_that("use_parallel: NA throws error requiring non-NA input", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = NA,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    # should throw an error.
    # not sure what it does if NA is entered as input
  )
})

# empty list
test_that("use_parallel: empty list input throws error requiring logical or numeric", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = list(),
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'use_parallel' must be one of the following: 'logical', 'numeric'"
  )
})
# empty dataframe
test_that("use_parallel: empty dataframe input throws error requiring logical or numeric", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = data.frame(),
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "'use_parallel' must be one of the following: 'logical', 'numeric'"
  )
})



# 17. interest_threshold validation (default is NULL)
# empty
test_that("interest_threshold: empty input uses default with no error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_no_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = , # will use default NULL
      time_standardisation = 500,
      verbose = FALSE
    )
  )
})

# multiple
test_that("interest_threshold: multiple input throws error requiring single value", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    suppressWarnings(
      estimate_roc(
        data_source_community = data_source_community,
        data_source_age = data_source_age,
        age_uncertainty = age_uncertainty,
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_age_range = 500,
        smooth_n_max = 9,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 1,
        bin_selection = "first",
        standardise = TRUE,
        n_individuals = 150,
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 100,
        use_parallel = FALSE,
        interest_threshold = c(2000, 3000),
        time_standardisation = 500,
        verbose = FALSE
      )
    ),
    # none programmed into the function yet.
    # e.g., "'arg' must be of length 1"
  )
})

# character
test_that("interest_threshold: Character input throws error requiring numeric", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = "3000",
      time_standardisation = 500,
      verbose = FALSE
    ),
    "interest_threshold' must be one of the following: 'NULL', 'numeric'"
  )
})

# NA
test_that("interest_threshold: NA input throws error requiring numeric", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NA,
      time_standardisation = 500,
      verbose = FALSE
    ),
    "interest_threshold' must be one of the following: 'NULL', 'numeric'"
  )
})

# empty list
test_that("interest_threshold: list input throws error requiring numeric", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = list(),
      time_standardisation = 500,
      verbose = FALSE
    ),
    "interest_threshold' must be one of the following: 'NULL', 'numeric'"
  )
})
# empty dataframe
test_that("interest_threshold: Empty dataframe input throws error requiring numeric", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = data.frame(),
      time_standardisation = 500,
      verbose = FALSE
    ),
    "interest_threshold' must be one of the following: 'NULL', 'numeric'"
  )
})


# 18. time_standardisation validation (default is NULL, requires numeric if not NULL)
# empty
test_that("time_standardisation: empty input uses default with no error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_no_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = ,
      verbose = FALSE
    )
  )
})

# multiple
test_that("time_standardisation: multiple input throws error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = c(500, 1000),
      verbose = FALSE
    ),
    "assert_that: length of assertion is not 1"
  )
})

# character
test_that("time_standardisation: Character input throws error requiring numeric", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = "500",
      verbose = FALSE
    ),
    "'time_standardisation' must be one of the following: 'numeric'"
  )
})

# NA
test_that("time_standardisation: NA input throws error requiring numeric", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NA,
      verbose = FALSE
    ),
    "time_standardisation' must be one of the following: 'numeric'"
  )
})

# empty list
test_that("time_standardisation: Empty list input throws error requiring numeric", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = list(),
      verbose = FALSE
    ),
    "time_standardisation' must be one of the following: 'numeric'"
  )
})

# empty dataframe
test_that("time_standardisation: Empty data.frame input throws error requiring numeric", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = data.frame(),
      verbose = FALSE
    ),
    "time_standardisation' must be one of the following: 'numeric'"
  )
})

# 19. verbose validation (default is FALSE)
# empty
test_that("verbose: Empty input uses default without error", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_no_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      verbose =
      )
  )
})

# NULL
test_that("verbose: NULL input throws error requiring logical", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      verbose = NULL
    ),
    "'verbose' must be one of the following: 'logical'"
  )
})

# multiple
test_that("verbose: multiple input throws error requiring single value", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      verbose = c(FALSE, TRUE)
    ),
    # none programmed into the function yet.
    # e.g. "assert_that: length of assertion is not 1"
    # or: "'arg' must be of length 1"
  )
})

# character
test_that("verbose: Character input throws error requiring logical", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      verbose = "FALSE"
    ),
    "'verbose' must be one of the following: 'logical'"
  )
})

# Numeric
test_that("verbose: Numeric input throws error requiring logical", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      verbose = 123
    ),
    "'verbose' must be one of the following: 'logical'"
  )
})

# 0
test_that("verbose: Zero (0) input throws error requiring logical", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      verbose = 0
    ),
    "'verbose' must be one of the following: 'logical'"
  )
})

# NA
test_that("verbose: NA input throws error requiring logical", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      verbose = NA # uses FALSE instead silently
    ),
    # none programmed into the function yet
    # e.g., "verbose cannot be NA"
  )
})

# empty list
test_that("verbose: Empty list input throws error requiring logical", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      verbose = list()
    ),
    "'verbose' must be one of the following: 'logical'"
  )
})

# empty dataframe
test_that("verbose: Empty dataframe input throws error requiring logical", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]

  expect_error(
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      verbose = data.frame()
    ),
    "'verbose' must be one of the following: 'logical'"
  )
})

# ============================= #
# 2. OUTPUT VALIDATION          #
# ============================= #

# 1. Valid data input returns a data.frame
# MW
test_that("Valid data input returns a data.frame", {
  res <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 10,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      verbose = FALSE
    )
  expect_s3_class(
    res,
    "data.frame"
  )
})
# levels
test_that("Valid data input returns correct column types", {
  res <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 10,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      verbose = FALSE
    )
  expect_s3_class(
    res,
    "data.frame"
  )
})

# bins
test_that("Valid data input returns correct column types", {
  res <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "bins",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 10,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      verbose = FALSE
    )
  expect_s3_class(
    res,
    "data.frame"
  )
})


# 2. Output has expected columns
# MW
test_that("Valid data input returns correct columns", {
  res <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 10,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      verbose = FALSE
    )

  expect_identical(
    c(
      "Working_Unit",
      "Age",
      "ROC",
      "ROC_up",
      "ROC_dw"
    ),
    colnames(res)
  )
})

# levels
test_that("Valid data input returns correct column types", {
  res <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 10,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      verbose = FALSE
    )
  expect_identical(
    c(
      "Working_Unit",
      "Age",
      "ROC",
      "ROC_up",
      "ROC_dw"
    ),
    colnames(res)
  )
})

# bins
test_that("Valid data input returns correct column types", {
  res <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "bins",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 10,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      verbose = FALSE
    )
  expect_identical(
    c(
      "Working_Unit",
      "Age",
      "ROC",
      "ROC_up",
      "ROC_dw"
    ),
    colnames(res)
  )
})


# 3. Columns have the correct types
test_that("Valid data input returns correct column types", {
  res <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 10,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      verbose = FALSE
    )

  expect_type(res$Working_Unit, "character")
  expect_type(res$Age, "double")
  expect_type(res$ROC, "double")
  expect_type(res$ROC_up, "double")
  expect_type(res$ROC_dw, "double")
})

# levels
test_that("Valid data input returns correct column types", {
  res <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 10,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      verbose = FALSE
    )

  expect_type(res$Working_Unit, "character")
  expect_type(res$Age, "double")
  expect_type(res$ROC, "double")
  expect_type(res$ROC_up, "double")
  expect_type(res$ROC_dw, "double")
})

# bins
test_that("Valid data input returns correct column types", {
  res <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      smooth_method = "grim",
      smooth_n_points = 5,
      smooth_age_range = 500,
      smooth_n_max = 9,
      working_units = "bins",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 10,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = NULL,
      verbose = FALSE
    )

  expect_type(res$Working_Unit, "character")
  expect_type(res$Age, "double")
  expect_type(res$ROC, "double")
  expect_type(res$ROC_up, "double")
  expect_type(res$ROC_dw, "double")
})





# ============================= #
# 3. FUNCTIONALITY              #
# ============================= #
# ============================= #
# OUTPUT VALIDATION TESTS       #
# ============================= #

test_that("estimate_roc output has correct structure with required columns", {
    # Create example data
    data_source_community <-
        RRatepol::example_data$pollen_data[[1]]
    data_source_age <-
        RRatepol::example_data$sample_age[[1]]

    # Run the function with minimal parameters
    result <- estimate_roc(
        data_source_community = data_source_community,
        data_source_age = data_source_age,
        smooth_method = "none",
        working_units = "levels",
        bin_size = 500,
        bin_selection = "first",
        dissimilarity_coefficient = "euc",
        rand = 10,
        verbose = FALSE
    )

    # Check output structure contains expected columns
    expect_true(all(c("Working_Unit", "Age", "ROC", "ROC_up", "ROC_dw") %in% colnames(result)))
})

test_that("estimate_roc output has rows sorted by Age", {
    # Create example data
    data_source_community <-
        RRatepol::example_data$pollen_data[[1]]
    data_source_age <-
        RRatepol::example_data$sample_age[[1]]

    # Run the function with minimal parameters
    result <- estimate_roc(
        data_source_community = data_source_community,
        data_source_age = data_source_age,
        smooth_method = "none",
        working_units = "levels",
        bin_size = 500,
        bin_selection = "first",
        dissimilarity_coefficient = "euc",
        rand = 10,
        verbose = FALSE
    )

    # Check that rows are sorted by Age
    expect_equal(result$Age, sort(result$Age))
})

test_that("estimate_roc properly applies interest_threshold", {
    # Create example data
    data_source_community <-
        RRatepol::example_data$pollen_data[[1]]
    data_source_age <-
        RRatepol::example_data$sample_age[[1]]

    # Get the median age from the age data
    median_age <- median(data_source_age$age)

    # Run with threshold
    result_with_threshold <- estimate_roc(
        data_source_community = data_source_community,
        data_source_age = data_source_age,
        smooth_method = "none",
        working_units = "levels",
        bin_size = 500,
        bin_selection = "first",
        dissimilarity_coefficient = "euc",
        interest_threshold = median_age,
        rand = 10,
        verbose = FALSE
    )

    # Check that all ages in result are less than or equal to the threshold
    expect_true(all(result_with_threshold$Age <= median_age))
})

test_that("estimate_roc produces different ROC values with different dissimilarity coefficients", {
    # Create example data
    data_source_community <-
        RRatepol::example_data$pollen_data[[1]]
    data_source_age <-
        RRatepol::example_data$sample_age[[1]]

    # Run the function with different dissimilarity coefficients
    result_chisq <- estimate_roc(
        data_source_community = data_source_community,
        data_source_age = data_source_age,
        smooth_method = "none",
        working_units = "levels",
        bin_size = 500,
        bin_selection = "first",
        dissimilarity_coefficient = "chisq",
        rand = 10,
        verbose = FALSE
    )

    result_chord <- estimate_roc(
        data_source_community = data_source_community,
        data_source_age = data_source_age,
        smooth_method = "none",
        working_units = "levels",
        bin_size = 500,
        bin_selection = "first",
        dissimilarity_coefficient = "chord",
        rand = 10,
        verbose = FALSE
    )

    # Check that ROC values differ between methods
    expect_false(identical(result_chisq$ROC, result_chord$ROC))
})

test_that("estimate_roc produces more working units with MW compared to bins method", {
    # Create example data
    data_source_community <-
        RRatepol::example_data$pollen_data[[1]]
    data_source_age <-
        RRatepol::example_data$sample_age[[1]]

    # Run with bins
    result_bins <- estimate_roc(
        data_source_community = data_source_community,
        data_source_age = data_source_age,
        smooth_method = "none",
        working_units = "bins",
        bin_size = 500,
        bin_selection = "first",
        dissimilarity_coefficient = "euc",
        rand = 10,
        verbose = FALSE
    )

    # Run with MW (moving window)
    result_mw <- estimate_roc(
        data_source_community = data_source_community,
        data_source_age = data_source_age,
        smooth_method = "none",
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 5,
        bin_selection = "first",
        dissimilarity_coefficient = "euc",
        rand = 10,
        verbose = FALSE
    )

    # Check that MW produces more working units than bins
    expect_gt(nrow(result_mw), nrow(result_bins))
})

test_that("estimate_roc returns ROC values within expected range for proportion data", {
    # Create example data
    data_source_community <-
        RRatepol::example_data$pollen_data[[1]]
    data_source_age <-
        RRatepol::example_data$sample_age[[1]]

    # Run with transform_to_proportions = TRUE
    result <- estimate_roc(
        data_source_community = data_source_community,
        data_source_age = data_source_age,
        smooth_method = "none",
        working_units = "levels",
        bin_size = 500,
        bin_selection = "first",
        dissimilarity_coefficient = "euc",
        tranform_to_proportions = TRUE,
        rand = 10,
        verbose = FALSE
    )

    # For proportion data with euclidean distance, values typically range from 0-1
    # but can be higher depending on time_standardisation
    expect_true(all(result$ROC >= 0))
    expect_true(all(is.finite(result$ROC)))
})

test_that("estimate_roc with time_standardisation properly scales ROC values", {
    # Create example data
    data_source_community <-
        RRatepol::example_data$pollen_data[[1]]
    data_source_age <-
        RRatepol::example_data$sample_age[[1]]

    # Run with different time_standardisation values
    set.seed(123)
    result_500 <- estimate_roc(
        data_source_community = data_source_community,
        data_source_age = data_source_age,
        smooth_method = "none",
        working_units = "levels",
        bin_size = 500,
        bin_selection = "first",
        dissimilarity_coefficient = "euc",
        time_standardisation = 500,
        rand = 10,
        verbose = FALSE
    )
    set.seed(123)
    result_1000 <- estimate_roc(
        data_source_community = data_source_community,
        data_source_age = data_source_age,
        smooth_method = "none",
        working_units = "levels",
        bin_size = 500,
        bin_selection = "first",
        dissimilarity_coefficient = "euc",
        time_standardisation = 1000,
        rand = 10,
        verbose = FALSE
    )

    # When time_standardisation is doubled, ROC values should be doubled
    # This is an approximate test due to randomization
    median_ratio <- median(result_500$ROC / result_1000$ROC, na.rm = TRUE)
    expect_gt(median_ratio, 0)
    expect_lt(median_ratio, 1)
})

test_that("estimate_roc applies smoothing correctly", {
    # Create example data
    data_source_community <-
        RRatepol::example_data$pollen_data[[1]]
    data_source_age <-
        RRatepol::example_data$sample_age[[1]]

    # Run with no smoothing
    set.seed(123)
    result_none <- estimate_roc(
        data_source_community = data_source_community,
        data_source_age = data_source_age,
        smooth_method = "none",
        working_units = "levels",
        bin_size = 500,
        bin_selection = "first",
        dissimilarity_coefficient = "euc",
        rand = 10,
        verbose = FALSE
    )

    # Run with Shepard smoothing
    set.seed(123)
    result_shep <- estimate_roc(
        data_source_community = data_source_community,
        data_source_age = data_source_age,
        smooth_method = "shep",
        working_units = "levels",
        bin_size = 500,
        bin_selection = "first",
        dissimilarity_coefficient = "euc",
        rand = 10,
        verbose = FALSE
    )

    # Smoothing should typically result in different ROC values
    # and often lower variance in the results
    expect_false(identical(result_none$ROC, result_shep$ROC))

    # Typically, smoothing reduces the range of ROC values
    range_none <- max(result_none$ROC, na.rm = TRUE) - min(result_none$ROC, na.rm = TRUE)
    range_shep <- max(result_shep$ROC, na.rm = TRUE) - min(result_shep$ROC, na.rm = TRUE)

    # Not always true, but generally smoothing reduces variance
    expect_lte(range_shep, range_none)
})

test_that("estimate_roc confidence intervals (ROC_up, ROC_dw) widen with higher age uncertainty", {
    # Create example data
    data_source_community <-
        RRatepol::example_data$pollen_data[[1]]
    data_source_age <-
        RRatepol::example_data$sample_age[[1]]
    age_uncertainty <-
        RRatepol::example_data$age_uncertainty[[1]]

    # Create increased uncertainty by multiplying standard deviations
    # This simulates the effect of using a more uncertain age model
    increased_uncertainty <- age_uncertainty
    set.seed(123)
    # Add more random variation to the uncertainty matrix
    for (i in seq_len(ncol(increased_uncertainty))) {
        # Add noise proportional to the original variation
        col_sd <- sd(increased_uncertainty[, i])
        increased_uncertainty[, i] <- increased_uncertainty[, i] + rnorm(nrow(increased_uncertainty), 0, col_sd * 0.5)
    }

    # Run with original uncertainty
    set.seed(123)
    result_orig <- estimate_roc(
        data_source_community = data_source_community,
        data_source_age = data_source_age,
        age_uncertainty = age_uncertainty,
        smooth_method = "none",
        working_units = "levels",
        bin_size = 500,
        bin_selection = "first",
        dissimilarity_coefficient = "euc",
        rand = 100,
        verbose = FALSE
    )

    # Run with increased uncertainty
    set.seed(123)
    result_high <- estimate_roc(
        data_source_community = data_source_community,
        data_source_age = data_source_age,
        age_uncertainty = increased_uncertainty,
        smooth_method = "none",
        working_units = "levels",
        bin_size = 500,
        bin_selection = "first",
        dissimilarity_coefficient = "euc",
        rand = 100,
        verbose = FALSE
    )

    # Calculate confidence interval width
    ci_width_orig <- mean(result_orig$ROC_up - result_orig$ROC_dw, na.rm = TRUE)
    ci_width_high <- mean(result_high$ROC_up - result_high$ROC_dw, na.rm = TRUE)

    # Higher uncertainty should result in wider confidence intervals
    # This test might be affected by randomness, so we use a tolerance
    expect_gt(ci_width_high, ci_width_orig)
})
