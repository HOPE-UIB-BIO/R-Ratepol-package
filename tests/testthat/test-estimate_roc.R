# Parameters in the function:

# # data_source_community,
# # data_source_age,
# # age_uncertainty = NULL,
# # smooth_method = c("none", "m.avg", "grim", "age.w", "shep"),
# # smooth_n_points = 5,
# # smooth_age_range = 500,
# # smooth_n_max = 9,
# # working_units = c("levels", "bins", "MW"),
# # bin_size = 500,
# # number_of_shifts = 5,
# # bin_selection = c("random", "first"),
# # standardise = FALSE,
# # n_individuals = 150,
# # dissimilarity_coefficient = c("euc", "euc.sd", "chord", "chisq", "gower", "bray"),
# # tranform_to_proportions = TRUE,
# # rand = NULL,
# # use_parallel = FALSE,
# # interest_threshold = NULL,
# # time_standardisation = NULL,
# # verbose = FALSE


test_that("estimate_roc works with valid inputs", {
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
test_that("estimate_roc throws error when data_source_community argument is missing (no default value provided)", {
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

test_that("estimate_roc throws error when data_source_community is NULL (invalid empty input)", {
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
test_that("estimate_roc throws error when data_source_community is a character string ('my_data') instead of required data structure", {
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
test_that("estimate_roc throws error when data_source_community is a numeric value (123) instead of required data structure", {
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
test_that("estimate_roc throws error when data_source_community is zero (numeric 0) instead of required data structure", {
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
test_that("estimate_roc throws error when data_source_community is NA (missing value) instead of required data structure", {
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
test_that("estimate_roc throws error when data_source_community is an empty list (list() with no elements)", {
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
test_that("estimate_roc throws error when data_source_community is an empty data frame (data.frame() with no rows/columns)", {
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
test_that("estimate_roc throws error when data_source_age argument is missing (no default value provided)", {
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

test_that("estimate_roc throws error when data_source_age is NULL (invalid empty input)", {
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
test_that("estimate_roc throws error when data_source_age is a character string ('my_data') instead of required data structure", {
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
test_that("estimate_roc throws error when data_source_age is a numeric value (123) instead of required data structure", {
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
test_that("estimate_roc throws error when data_source_age is zero (numeric 0) instead of required data structure", {
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
test_that("estimate_roc throws error when data_source_age is NA (missing value) instead of required data structure", {
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
test_that("estimate_roc throws error when data_source_age is an empty list (list() with no elements)", {
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
test_that("estimate_roc throws error when data_source_age is an empty data frame (data.frame() with no rows/columns)", {
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
test_that("estimate_roc works when age_uncertainty is NULL (no uncertainty randomizations)", {
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


test_that("estimate_roc works when age_uncertainty is NULL (no uncertainty randomizations)", {
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
test_that("estimate_roc works when age_uncertainty is NULL (no uncertainty randomizations)", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  # Run the function with valid inputs
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
test_that("estimate_roc works when age_uncertainty is NULL (no uncertainty randomizations)", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  # Run the function with valid inputs
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
test_that("estimate_roc works when age_uncertainty is NULL (no uncertainty randomizations)", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  # Run the function with valid inputs
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
test_that("estimate_roc works when age_uncertainty is NULL (no uncertainty randomizations)", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  # Run the function with valid inputs
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

# empty list
test_that("estimate_roc works when age_uncertainty is NULL (no uncertainty randomizations)", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  # Run the function with valid inputs
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

# empty dataframe
test_that("estimate_roc works when age_uncertainty is NULL (no uncertainty randomizations)", {
  # Create example data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]

  # Run the function with valid inputs
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

# 3. smooth_method validation

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
    "length of assertion is not 1"
  )
})

# multiple
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
test_that("estimate_roc throws error when smooth_method is NULL (invalid empty input)", {
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
      smooth_n_max = , # uses 9 as default
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
  )
})

# NULL
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
    "assert_that: length of assertion is not 1"
  )
})

# multiple
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
# NULL
# multiple
# character
# Numeric
# 0
# NA
# empty list
# empty dataframe

# 8. bin_size validation
# empty
# NULL
# multiple
# character
# Numeric
# 0
# NA
# empty list
# empty dataframe

# 9. number_of_shifts validation
# empty
# NULL
# multiple
# character
# Numeric
# 0
# NA
# empty list
# empty dataframe

# 10. bin_selection validation
# empty
# NULL
# multiple
# character
# Numeric
# 0
# NA
# empty list
# empty dataframe

# 11. standardise validation
# empty
# NULL
# multiple
# character
# Numeric
# 0
# NA
# empty list
# empty dataframe

# 12. n_individuals validation
# empty
# NULL
# multiple
# character
# Numeric
# 0
# NA
# empty list
# empty dataframe

# 13. dissimilarity_coefficient validation
# empty
# NULL
# multiple
# character
# Numeric
# 0
# NA
# empty list
# empty dataframe

# 14. tranform_to_proportions validation
# empty
# NULL
# multiple
# character
# Numeric
# 0
# NA
# empty list
# empty dataframe

# 15. rand validation
# empty
# NULL
# multiple
# character
# Numeric
# 0
# NA
# empty list
# empty dataframe

# 16. use_parallel validation
# empty
# NULL
# multiple
# character
# Numeric
# 0
# NA
# empty list
# empty dataframe

# 17. interest_threshold validation
# empty
# NULL
# multiple
# character
# Numeric
# 0
# NA
# empty list
# empty dataframe

# 18. time_standardisation validation
# empty
# NULL
# multiple
# character
# Numeric
# 0
# NA
# empty list
# empty dataframe

# 19. verbose validation
# empty
# NULL
# multiple
# character
# Numeric
# 0
# NA
# empty list
# empty dataframe



# ============================= #
# 2. OUTPUT VALIDATION          #
# ============================= #



# ============================= #
# 3. FUNCTIONALITY              #
# ============================= #

