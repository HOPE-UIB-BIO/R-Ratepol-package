test_that("detect_peak_points() errors if 'data_source' argument is missing: expects error about missing required argument", {
  expect_error(
    detect_peak_points(
      data_source = ,
      sel_method = "trend_linear",
      sd_threshold = 2
    ),
    'argument "data_source" is missing, with no default'
  )
})

test_that("detect_peak_points() errors if 'data_source' is NULL: expects type error requiring data.frame", {
  expect_error(
    detect_peak_points(
      data_source = NULL,
      sel_method = "trend_linear",
      sd_threshold = 2
    ),
    "'data_source' must be one of the following: 'data.frame'"
  )
})

test_that("detect_peak_points() errors if 'data_source' is character: expects type error requiring data.frame", {
  expect_error(
    detect_peak_points(
      data_source = "my_data",
      sel_method = "trend_linear",
      sd_threshold = 2
    ),
    "'data_source' must be one of the following: 'data.frame'"
  )
})

test_that("detect_peak_points() errors if 'data_source' is numeric: expects type error requiring data.frame", {
  expect_error(
    detect_peak_points(
      data_source = 123,
      sel_method = "trend_linear",
      sd_threshold = 2
    ),
    "'data_source' must be one of the following: 'data.frame'"
  )
})

test_that("detect_peak_points() errors if 'data_source' is NA: expects type error requiring data.frame", {
  expect_error(
    detect_peak_points(
      data_source = NA,
      sel_method = "trend_linear",
      sd_threshold = 2
    ),
    "'data_source' must be one of the following: 'data.frame'"
  )
})

test_that("detect_peak_points() errors if 'data_source' is an empty list: expects type error requiring data.frame", {
  expect_error(
    detect_peak_points(
      data_source = list(),
      sel_method = "trend_linear",
      sd_threshold = 2
    ),
    "'data_source' must be one of the following: 'data.frame'"
  )
})

test_that("detect_peak_points() errors if 'data_source' is an empty data.frame: expects missing required columns error", {
  expect_error(
    detect_peak_points(
      data_source = data.frame(),
      sel_method = "trend_linear",
      sd_threshold = 2
    ),
    "'data_source' must contains following columns: 'ROC'"
  )
})

test_that("detect_peak_points() errors if 'data_source' has empty ROC and Age columns: expects model fitting error", {
  expect_error(
    suppressWarnings(
      detect_peak_points(
        data_source = data.frame(ROC = numeric(), Age = numeric()),
        sel_method = "trend_linear",
        sd_threshold = 2
      )
    ),
    "object 'fit' not found"
  )
})

test_that("detect_peak_points() errors if 'data_source' has NA in ROC and Age columns: expects model fitting error", {
  expect_error(
    suppressWarnings(
      detect_peak_points(
        data_source = data.frame(ROC = NA, Age = NA),
        sel_method = "trend_linear",
        sd_threshold = 2
      )
    ),
    "object 'fit' not found"
  )
})

test_that("detect_peak_points() uses default sel_method 'trend_linear' when not supplied: expects no error and default method used", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_no_error(
    detect_peak_points(
      data_source = data_source,
      sel_method = , # uses 'trend_linear' as default
      sd_threshold = 2
    )
  )
})

test_that("detect_peak_points() errors if sel_method is NULL: expects character type error for method parameter", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_error(
    detect_peak_points(
      data_source = data_source,
      sel_method = NULL,
      sd_threshold = 2
    ),
    "'sel_method' must be one of the following: 'character'"
  )
})

test_that("detect_peak_points() errors if sel_method is invalid character: expects allowed values error for method parameter", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_error(
    detect_peak_points(
      data_source = data_source,
      sel_method = "my_method",
      sd_threshold = 2
    ),
    "'sel_method' must contains one of the following values: 'trend_linear', 'trend_non_linear', 'threshold', 'GAM_deriv', 'SNI'"
  )
})

test_that("detect_peak_points() errors if sel_method is a vector of multiple values: expects length 1 error for method parameter", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_error(
    detect_peak_points(
      data_source = data_source,
      sel_method = c("SNI", "threshold"),
      sd_threshold = 2
    ),
    "'arg' must be of length 1"
  )
})

test_that("detect_peak_points() errors if sel_method is numeric: expects character type error for method parameter", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_error(
    detect_peak_points(
      data_source = data_source,
      sel_method = 123,
      sd_threshold = 2
    ),
    "'sel_method' must be one of the following: 'character'"
  )
})

test_that("detect_peak_points() errors if sel_method is zero: expects character type error for method parameter", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_error(
    detect_peak_points(
      data_source = data_source,
      sel_method = 0,
      sd_threshold = 2
    ),
    "'sel_method' must be one of the following: 'character'"
  )
})

test_that("detect_peak_points() errors if sel_method is NA: expects character type error for method parameter", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_error(
    detect_peak_points(
      data_source = data_source,
      sel_method = NA,
      sd_threshold = 2
    ),
    "'sel_method' must be one of the following: 'character'"
  )
})

test_that("detect_peak_points() errors if sel_method is an empty list: expects character type error for method parameter", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_error(
    detect_peak_points(
      data_source = data_source,
      sel_method = list(),
      sd_threshold = 2
    ),
    "'sel_method' must be one of the following: 'character'"
  )
})

test_that("detect_peak_points() errors if sel_method is an empty data.frame: expects character type error for method parameter", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_error(
    detect_peak_points(
      data_source = data_source,
      sel_method = data.frame(),
      sd_threshold = 2
    ),
    "'sel_method' must be one of the following: 'character'"
  )
})

test_that("detect_peak_points() uses default sd_threshold=2 when not supplied: expects no error and default threshold used", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_no_error(
    detect_peak_points(
      data_source = data_source,
      sel_method = "trend_linear",
      sd_threshold = # uses 2 as default silently
      )
  )
})

test_that("detect_peak_points() errors if sd_threshold is NULL: expects numeric type error for threshold parameter", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_error(
    detect_peak_points(
      data_source = data_source,
      sel_method = "trend_linear",
      sd_threshold = NULL
    ),
    "'sd_threshold' must be one of the following: 'numeric'"
  )
})

test_that("detect_peak_points() errors if sd_threshold is character: expects numeric type error for threshold parameter", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_error(
    detect_peak_points(
      data_source = data_source,
      sel_method = "trend_linear",
      sd_threshold = "my_threshold"
    ),
    "'sd_threshold' must be one of the following: 'numeric'"
  )
})

test_that("detect_peak_points() errors if sd_threshold is zero: expects value greater than zero error for threshold parameter", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_error(
    detect_peak_points(
      data_source = data_source,
      sel_method = "trend_linear",
      sd_threshold = 0
    ),
    "'sd_threshold' must be bigger than 0"
  )
})

test_that("detect_peak_points() errors if sd_threshold is NA: expects numeric type error for threshold parameter", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_error(
    detect_peak_points(
      data_source = data_source,
      sel_method = "trend_linear",
      sd_threshold = NA
    ),
    "'sd_threshold' must be one of the following: 'numeric'"
  )
})

test_that("detect_peak_points() errors if sd_threshold is an empty list: expects numeric type error for threshold parameter", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_error(
    detect_peak_points(
      data_source = data_source,
      sel_method = "trend_linear",
      sd_threshold = list()
    ),
    "'sd_threshold' must be one of the following: 'numeric'"
  )
})

test_that("detect_peak_points() errors if sd_threshold is an empty data.frame: expects numeric type error for threshold parameter", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_error(
    detect_peak_points(
      data_source = data_source,
      sel_method = "trend_linear",
      sd_threshold = data.frame()
    ),
    "'sd_threshold' must be one of the following: 'numeric'"
  )
})

test_that("detect_peak_points() errors if sd_threshold is a vector of multiple values: expects length 1 assertion error for threshold parameter", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_error(
    detect_peak_points(
      data_source = data_source,
      sel_method = "trend_linear",
      sd_threshold = c(2, 6)
    ),
    "assert_that: length of assertion is not 1"
  )
})

# Output validation
test_that("detect_peak_points() returns a data.frame output with valid input and parameters", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  res <-
    detect_peak_points(
      data_source = data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )

  expect_s3_class(
    res,
    "data.frame"
  )
})

test_that("detect_peak_points() returns output with expected column names for valid input and parameters", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  res <-
    detect_peak_points(
      data_source = data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )

  expect_named(
    res,
    c("Working_Unit", "Age", "ROC", "ROC_up", "ROC_dw", "Peak")
  )
})

# test for different sel_method values
# and check if output is as expected
test_that("detect_peak_points() returns output with expected columns for sel_method='trend_non_linear'", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  res <-
    detect_peak_points(
      data_source = data_source,
      sel_method = "trend_non_linear",
      sd_threshold = 2
    )

  expect_named(
    res,
    c("Working_Unit", "Age", "ROC", "ROC_up", "ROC_dw", "Peak")
  )
})

test_that("detect_peak_points() returns output with expected columns for sel_method='threshold'", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  res <-
    detect_peak_points(
      data_source = data_source,
      sel_method = "threshold",
      sd_threshold = 2
    )

  expect_named(
    res,
    c("Working_Unit", "Age", "ROC", "ROC_up", "ROC_dw", "Peak")
  )
})

test_that("detect_peak_points() returns output with expected columns for sel_method='GAM_deriv'", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  res <-
    detect_peak_points(
      data_source = data_source,
      sel_method = "GAM_deriv",
      sd_threshold = 2
    )

  expect_named(
    res,
    c("Working_Unit", "Age", "ROC", "ROC_up", "ROC_dw", "Peak")
  )
})

test_that("detect_peak_points() returns output with expected columns for sel_method='SNI'", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  res <-
    detect_peak_points(
      data_source = data_source,
      sel_method = "SNI",
      sd_threshold = 2
    )

  expect_named(
    res,
    c("Working_Unit", "Age", "ROC", "ROC_up", "ROC_dw", "Peak")
  )
})

# Testing the correct functionality

test_that("detect_peak_points() with method 'threshold' correctly identifies peaks", {
  # Prepare test data with known properties
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  # Create a version with artificially high ROC_dw values to ensure peaks
  high_roc_data <-
    data_source
  median_roc <-
    median(high_roc_data$ROC)
  # Set the first 3 ROC_dw values to be very high (ensuring they'll be detected as peaks)
  high_roc_data$ROC_dw[1:3] <-
    median_roc * 10

  result <-
    detect_peak_points(
      high_roc_data,
      sel_method = "threshold"
    )

  # Check that at least the first 3 values are detected as peaks
  expect_true(
    all(
      result$Peak[1:3]
    )
  )
})

test_that("detect_peak_points() with method 'trend_linear' correctly identifies peaks", {
  # Prepare test data with known properties
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  # Create a version with artificially high ROC values to ensure peaks
  high_roc_data <-
    data_source
  # Find the mean and SD of ROC values
  mean_roc <-
    mean(high_roc_data$ROC)
  sd_roc <-
    sd(high_roc_data$ROC)
  # Set the first 3 ROC values to be very high (ensuring they'll be detected as peaks)
  high_roc_data$ROC[1:3] <-
    mean_roc + (sd_roc * 5)

  result <-
    detect_peak_points(
      high_roc_data,
      sel_method = "trend_linear",
      sd_threshold = 2
    )

  # Check that at least the first 3 values are detected as peaks
  expect_true(
    all(
      result$Peak[1:3]
    )
  )
})

test_that("detect_peak_points() with method 'trend_non_linear' correctly identifies peaks", {
  # Prepare test data with known properties
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  # Create a version with artificially high ROC values to ensure peaks
  high_roc_data <-
    data_source
  # Find the mean and SD of ROC values
  mean_roc <-
    mean(high_roc_data$ROC)
  sd_roc <-
    sd(high_roc_data$ROC)
  # Set the first 3 ROC values to be very high (ensuring they'll be detected as peaks)
  high_roc_data$ROC[1:3] <-
    mean_roc + (sd_roc * 5)

  result <-
    detect_peak_points(
      high_roc_data,
      sel_method = "trend_non_linear",
      sd_threshold = 2
    )

  # Check that the result has a Peak column and it's logical
  expect_type(
    result$Peak,
    "logical"
  )
  # Check that at least the first 3 values are detected as peaks
  expect_true(
    all(
      result$Peak[1:3]
    )
  )
})

test_that("detect_peak_points() with method 'GAM_deriv' correctly processes the data", {
  # Prepare test data with known properties
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  result <-
    detect_peak_points(
      data_source,
      sel_method = "GAM_deriv"
    )

  # Check that the result has a Peak column and it's logical
  expect_type(
    result$Peak,
    "logical"
  )
})

test_that("detect_peak_points() with method 'SNI' correctly processes the data", {
  # Prepare test data with known properties
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  result <-
    detect_peak_points(
      data_source,
      sel_method = "SNI"
    )

  # Check that the result has a Peak column and it's logical
  expect_type(
    result$Peak,
    "logical"
  )
})

test_that("detect_peak_points() correctly responds to different sd_threshold values", {
  # Prepare test data with known properties
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  # Test with different sd_threshold values
  result_sd1 <-
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 1
    )
  result_sd3 <-
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 3
    )

  # Higher threshold should detect fewer peaks
  expect_gt(
    sum(result_sd1$Peak, na.rm = TRUE),
    sum(result_sd3$Peak, na.rm = TRUE)
  )
})

test_that("detect_peak_points() methods give different results", {
  # Prepare test data with known properties
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  # Test with different methods
  methods <-
    c(
      "threshold",
      "trend_linear",
      "trend_non_linear",
      "SNI"
    )
  results <-
    list()

  for (method in methods) {
    results[[method]] <-
      detect_peak_points(
        data_source,
        sel_method = method
      )
  }

  # Check that at least some methods give different results
  peak_sums <-
    sapply(
      results,
      function(df) sum(df$Peak, na.rm = TRUE)
    )

  expect_true(
    length(
      unique(peak_sums)
    ) > 1
  )
})

test_that("detect_peak_points() preserves original data structure with Peak column added", {
  # Prepare test data with known properties
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  original_cols <-
    names(data_source)

  # Check result
  result <-
    detect_peak_points(
      data_source,
      sel_method = "trend_linear"
    )

  # The Peak column should be added
  expect_true(
    all(
      c("Peak", original_cols) %in% colnames(result)
    )
  )
  # Row count should be preserved
  expect_equal(
    nrow(data_source),
    nrow(result)
  )
})

# Testing different input data scenarios created using estimate_roc()
# Check if no error:
## 1. bins
test_that("detect_peak_points() works without error for non-smoothed data with working_units='bins', using sel_method='trend_linear' and sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "none",
      smooth_n_points = 5,
      working_units = "bins",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_no_error(
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )
  )
})

test_that("detect_peak_points() works without error for m.avg-smoothed data with working_units='bins', using sel_method='trend_linear' and sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "m.avg",
      smooth_n_points = 5,
      working_units = "bins",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_no_error(
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )
  )
})

test_that("detect_peak_points() works without error for grim-smoothed data with working_units='bins', using sel_method='trend_linear' and sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "grim",
      smooth_n_points = 5,
      working_units = "bins",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_no_error(
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )
  )
})

test_that("detect_peak_points() works without error for age.w-smoothed data with working_units='bins', using sel_method='trend_linear' and sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "age.w",
      smooth_n_points = 5,
      working_units = "bins",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_no_error(
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )
  )
})

test_that("detect_peak_points() works without error for shep-smoothed data with working_units='bins', using sel_method='trend_linear' and sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "shep",
      smooth_n_points = 5,
      working_units = "bins",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_no_error(
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )
  )
})

## 2. levels
test_that("detect_peak_points() works without error for non-smoothed data with working_units='levels', using sel_method='trend_linear' and sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "none",
      smooth_n_points = 5,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_no_error(
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )
  )
})

test_that("detect_peak_points() works without error for m.avg-smoothed data with working_units='levels', using sel_method='trend_linear' and sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "m.avg",
      smooth_n_points = 5,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_no_error(
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )
  )
})

test_that("detect_peak_points() works without error for grim-smoothed data with working_units='levels', using sel_method='trend_linear' and sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "grim",
      smooth_n_points = 5,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_no_error(
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )
  )
})

test_that("detect_peak_points() works without error for age.w-smoothed data with working_units='levels', using sel_method='trend_linear' and sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "age.w",
      smooth_n_points = 5,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_no_error(
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )
  )
})

test_that("detect_peak_points() works without error for shep-smoothed data with working_units='levels', using sel_method='trend_linear' and sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_no_error(
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )
  )
})

## 3. MW
test_that("detect_peak_points() works without error for non-smoothed data with working_units='MW', using sel_method='trend_linear' and sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "none",
      smooth_n_points = 5,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_no_error(
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )
  )
})

test_that("detect_peak_points() works without error for m.avg-smoothed data with working_units='MW', using sel_method='trend_linear' and sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "m.avg",
      smooth_n_points = 5,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_no_error(
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )
  )
})

test_that("detect_peak_points() works without error for grim-smoothed data with working_units='MW', using sel_method='trend_linear' and sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "grim",
      smooth_n_points = 5,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_no_error(
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )
  )
})

test_that("detect_peak_points() works without error for age.w-smoothed data with working_units='MW', using sel_method='trend_linear' and sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "age.w",
      smooth_n_points = 5,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_no_error(
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )
  )
})

test_that("detect_peak_points() works without error for shep-smoothed data with working_units='MW', using sel_method='trend_linear' and sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "shep",
      smooth_n_points = 5,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  expect_no_error(
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )
  )
})

# Check if output structure is as expected (data.frame with correct columns):
## 1. bins
test_that("detect_peak_points() returns data.frame with Peak column and preserved row count for non-smoothed input, working_units='bins', sel_method='trend_linear', sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "none",
      smooth_n_points = 5,
      working_units = "bins",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  res <-
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )

  expect_s3_class(
    res, "data.frame"
  )
  expect_true(
    "Peak" %in% names(res)
  )
  expect_equal(
    nrow(data_source),
    nrow(res)
  )
})

test_that("detect_peak_points() returns data.frame with Peak column and preserved row count for m.avg-smoothed input, working_units='bins', sel_method='trend_linear', sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "m.avg",
      smooth_n_points = 5,
      working_units = "bins",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  res <-
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )

  expect_s3_class(
    res, "data.frame"
  )
  expect_true(
    "Peak" %in% names(res)
  )
  expect_equal(
    nrow(data_source),
    nrow(res)
  )
})

test_that("detect_peak_points() returns data.frame with Peak column and preserved row count for grim-smoothed input, working_units='bins', sel_method='trend_linear', sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "grim",
      smooth_n_points = 5,
      working_units = "bins",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  res <-
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )

  expect_s3_class(
    res, "data.frame"
  )
  expect_true(
    "Peak" %in% names(res)
  )
  expect_equal(
    nrow(data_source),
    nrow(res)
  )
})

test_that("detect_peak_points() returns data.frame with Peak column and preserved row count for age.w-smoothed input, working_units='bins', sel_method='trend_linear', sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "age.w",
      smooth_n_points = 5,
      working_units = "bins",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  res <-
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )

  expect_s3_class(
    res, "data.frame"
  )
  expect_true(
    "Peak" %in% names(res)
  )
  expect_equal(
    nrow(data_source),
    nrow(res)
  )
})

test_that("detect_peak_points() returns data.frame with Peak column and preserved row count for shep-smoothed input, working_units='bins', sel_method='trend_linear', sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "shep",
      smooth_n_points = 5,
      working_units = "bins",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  res <-
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )

  expect_s3_class(
    res, "data.frame"
  )
  expect_true(
    "Peak" %in% names(res)
  )
  expect_equal(
    nrow(data_source),
    nrow(res)
  )
})

## 2. levels
test_that("detect_peak_points() returns data.frame with Peak column and preserved row count for non-smoothed input, working_units='levels', sel_method='trend_linear', sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "none",
      smooth_n_points = 5,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  res <-
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )

  expect_s3_class(
    res, "data.frame"
  )
  expect_true(
    "Peak" %in% names(res)
  )
  expect_equal(
    nrow(data_source),
    nrow(res)
  )
})

test_that("detect_peak_points() returns data.frame with Peak column and preserved row count for m.avg-smoothed input, working_units='levels', sel_method='trend_linear', sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "m.avg",
      smooth_n_points = 5,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  res <-
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )

  expect_s3_class(
    res, "data.frame"
  )
  expect_true(
    "Peak" %in% names(res)
  )
  expect_equal(
    nrow(data_source),
    nrow(res)
  )
})

test_that("detect_peak_points() returns data.frame with Peak column and preserved row count for grim-smoothed input, working_units='levels', sel_method='trend_linear', sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "grim",
      smooth_n_points = 5,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  res <-
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )

  expect_s3_class(
    res, "data.frame"
  )
  expect_true(
    "Peak" %in% names(res)
  )
  expect_equal(
    nrow(data_source),
    nrow(res)
  )
})

test_that("detect_peak_points() returns data.frame with Peak column and preserved row count for age.w-smoothed input, working_units='levels', sel_method='trend_linear', sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "age.w",
      smooth_n_points = 5,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  res <-
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )

  expect_s3_class(
    res, "data.frame"
  )
  expect_true(
    "Peak" %in% names(res)
  )
  expect_equal(
    nrow(data_source),
    nrow(res)
  )
})

test_that("detect_peak_points() returns data.frame with Peak column and preserved row count for shep-smoothed input, working_units='levels', sel_method='trend_linear', sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
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
      use_parallel = FALSE,
      verbose = FALSE
    )

  res <-
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )

  expect_s3_class(
    res, "data.frame"
  )
  expect_true(
    "Peak" %in% names(res)
  )
  expect_equal(
    nrow(data_source),
    nrow(res)
  )
})

## 3. MW
test_that("detect_peak_points() returns data.frame with Peak column and preserved row count for non-smoothed input, working_units='MW', sel_method='trend_linear', sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "none",
      smooth_n_points = 5,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  res <-
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )

  expect_s3_class(
    res, "data.frame"
  )
  expect_true(
    "Peak" %in% names(res)
  )
  expect_equal(
    nrow(data_source),
    nrow(res)
  )
})

test_that("detect_peak_points() returns data.frame with Peak column and preserved row count for m.avg-smoothed input, working_units='MW', sel_method='trend_linear', sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "m.avg",
      smooth_n_points = 5,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  res <-
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )

  expect_s3_class(
    res, "data.frame"
  )
  expect_true(
    "Peak" %in% names(res)
  )
  expect_equal(
    nrow(data_source),
    nrow(res)
  )
})

test_that("detect_peak_points() returns data.frame with Peak column and preserved row count for grim-smoothed input, working_units='MW', sel_method='trend_linear', sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "grim",
      smooth_n_points = 5,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  res <-
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )

  expect_s3_class(
    res, "data.frame"
  )
  expect_true(
    "Peak" %in% names(res)
  )
  expect_equal(
    nrow(data_source),
    nrow(res)
  )
})

test_that("detect_peak_points() returns data.frame with Peak column and preserved row count for age.w-smoothed input, working_units='MW', sel_method='trend_linear', sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "age.w",
      smooth_n_points = 5,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  res <-
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )

  expect_s3_class(
    res, "data.frame"
  )
  expect_true(
    "Peak" %in% names(res)
  )
  expect_equal(
    nrow(data_source),
    nrow(res)
  )
})

test_that("detect_peak_points() returns data.frame with Peak column and preserved row count for shep-smoothed input, working_units='MW', sel_method='trend_linear', sd_threshold=2", {
  data_source <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "shep",
      smooth_n_points = 5,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      bin_selection = "first",
      standardise = FALSE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE
    )

  res <-
    detect_peak_points(
      data_source,
      sel_method = "trend_linear",
      sd_threshold = 2
    )

  expect_s3_class(
    res, "data.frame"
  )
  expect_true(
    "Peak" %in% names(res)
  )
  expect_equal(
    nrow(data_source),
    nrow(res)
  )
})

test_that("detect_peak_points() works with standardised data", {
  # Test with standardisation enabled
  data_standardised <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "shep",
      working_units = "levels",
      standardise = TRUE,
      n_individuals = 100,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE,
      rand = 100
    )

  # Test with standardisation disabled
  data_unstandardised <-
    estimate_roc(
      data_source_community = RRatepol::example_data$pollen_data[[1]],
      data_source_age = RRatepol::example_data$sample_age[[1]],
      smooth_method = "shep",
      working_units = "levels",
      standardise = FALSE,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      use_parallel = FALSE,
      verbose = FALSE,
      rand = 100
    )

  # Run detect_peak_points on both datasets
  result_standardised <-
    detect_peak_points(
      data_standardised,
      sel_method = "trend_linear",
      sd_threshold = 2
    )

  result_unstandardised <-
    detect_peak_points(
      data_unstandardised,
      sel_method = "trend_linear",
      sd_threshold = 2
    )

  # Check that both results have the Peak column
  expect_true(
    "Peak" %in% names(result_standardised)
  )
  expect_true(
    "Peak" %in% names(result_unstandardised)
  )
})
