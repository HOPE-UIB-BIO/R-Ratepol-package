test_that("detect_peak_points() throws error when 'data_source' parameter is missing: expects error about missing argument", {
  expect_error(
    detect_peak_points(
      data_source = ,
      sel_method = "trend_linear",
      sd_threshold = 2
    ),
    'argument "data_source" is missing, with no default'
  )
})

test_that("detect_peak_points() throws error when 'data_source' is NULL: expects type error requiring data.frame", {
  expect_error(
    detect_peak_points(
      data_source = NULL,
      sel_method = "trend_linear",
      sd_threshold = 2
    ),
    "'data_source' must be one of the following: 'data.frame'"
  )
})

test_that("detect_peak_points() throws error when 'data_source' is character: expects type error requiring data.frame", {
  expect_error(
    detect_peak_points(
      data_source = "my_data",
      sel_method = "trend_linear",
      sd_threshold = 2
    ),
    "'data_source' must be one of the following: 'data.frame'"
  )
})

test_that("detect_peak_points() throws error when 'data_source' is numeric: expects type error requiring data.frame", {
  expect_error(
    detect_peak_points(
      data_source = 123,
      sel_method = "trend_linear",
      sd_threshold = 2
    ),
    "'data_source' must be one of the following: 'data.frame'"
  )
})

test_that("detect_peak_points() throws error when 'data_source' is NA: expects type error requiring data.frame", {
  expect_error(
    detect_peak_points(
      data_source = NA,
      sel_method = "trend_linear",
      sd_threshold = 2
    ),
    "'data_source' must be one of the following: 'data.frame'"
  )
})

test_that("detect_peak_points() throws error when 'data_source' is an empty list: expects type error requiring data.frame", {
  expect_error(
    detect_peak_points(
      data_source = list(),
      sel_method = "trend_linear",
      sd_threshold = 2
    ),
    "'data_source' must be one of the following: 'data.frame'"
  )
})

test_that("detect_peak_points() throws error when 'data_source' is an empty data.frame: expects missing required columns error", {
  expect_error(
    detect_peak_points(
      data_source = data.frame(),
      sel_method = "trend_linear",
      sd_threshold = 2
    ),
    "'data_source' must contains following columns: 'ROC'"
  )
})

test_that("detect_peak_points() throws error when 'data_source' has empty ROC and Age columns: expects model fitting error", {
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

test_that("detect_peak_points() throws error when 'data_source' has NA in ROC and Age columns: expects model fitting error", {
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

test_that("detect_peak_points() uses default sel_method 'trend_linear' when not supplied: expects no error", {
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

test_that("detect_peak_points() throws error when sel_method is NULL: expects character type error", {
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

test_that("detect_peak_points() throws error when sel_method is invalid character: expects allowed values error", {
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

test_that("detect_peak_points() throws error when sel_method is a vector of multiple values: expects length 1 error", {
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

test_that("detect_peak_points() throws error when sel_method is numeric: expects character type error", {
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

test_that("detect_peak_points() throws error when sel_method is zero: expects character type error", {
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

test_that("detect_peak_points() throws error when sel_method is NA: expects character type error", {
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

test_that("detect_peak_points() throws error when sel_method is an empty list: expects character type error", {
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

test_that("detect_peak_points() throws error when sel_method is an empty data.frame: expects character type error", {
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

test_that("detect_peak_points() uses default sd_threshold=2 when not supplied: expects no error", {
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

test_that("detect_peak_points() throws error when sd_threshold is NULL: expects numeric type error", {
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

test_that("detect_peak_points() throws error when sd_threshold is character: expects numeric type error", {
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

test_that("detect_peak_points() throws error when sd_threshold is zero: expects value greater than zero error", {
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

test_that("detect_peak_points() throws error when sd_threshold is NA: expects numeric type error", {
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

test_that("detect_peak_points() throws error when sd_threshold is an empty list: expects numeric type error", {
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

test_that("detect_peak_points() throws error when sd_threshold is an empty data.frame: expects numeric type error", {
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

test_that("detect_peak_points() throws error when sd_threshold is a vector of multiple values: expects length 1 assertion error", {
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
test_that("detect_peak_points() returns valid output class with valid input", {
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

test_that("detect_peak_points() returns valid output class with valid input", {
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

# test for differernt sel_method values
# and check if output is as expected
test_that("detect_peak_points() returns valid output class with valid input", {
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

test_that("detect_peak_points() returns valid output class with valid input", {
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

test_that("detect_peak_points() returns valid output class with valid input", {
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

test_that("detect_peak_points() returns valid output class with valid input", {
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
