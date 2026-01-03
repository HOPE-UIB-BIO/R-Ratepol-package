# ------------------------------------------------------------------------
# Input validation tests for CharData parameter
# ------------------------------------------------------------------------

# empty
test_that("detect_sni rejects empty CharData input", {
  expect_error(
    detect_sni(
      CharData = ,
      BandWidth = 5
    ),
    'argument "CharData" is missing, with no default'
  )
})

# NULL
test_that("detect_sni rejects NULL CharData input", {
  expect_error(
    suppressWarnings(
      detect_sni(
        CharData = NULL,
        BandWidth = 5
      )
    ),
    "missing value where TRUE/FALSE needed"
  )
})

# character
test_that("detect_sni rejects character CharData input", {
  expect_error(
    suppressWarnings(
      detect_sni(
        CharData = "my_data",
        BandWidth = 5
      )
    ),
    "incorrect number of dimensions"
  )
})

# numeric
test_that("detect_sni rejects numeric CharData input", {
  expect_error(
    suppressWarnings(
      detect_sni(
        CharData = 123,
        BandWidth = 5
      )
    ),
    "incorrect number of dimensions"
  )
})

# zero
test_that("detect_sni rejects zero (0) CharData input", {
  expect_error(
    suppressWarnings(
      detect_sni(
        CharData = 0,
        BandWidth = 5
      )
    ),
    "incorrect number of dimensions"
  )
})

# NA
test_that("detect_sni rejects NA  CharData input", {
  expect_error(
    suppressWarnings(
      detect_sni(
        CharData = NA,
        BandWidth = 5
      )
    ),
    "incorrect number of dimensions"
  )
})

# list()
test_that("detect_sni rejects empty CharData input", {
  expect_error(
    detect_sni(
      CharData = list(),
      BandWidth = 5
    ),
    "incorrect number of dimensions"
  )
})

# data.frame()
test_that("detect_sni rejects empty CharData input", {
  expect_error(
    detect_sni(
      CharData = data.frame(),
      BandWidth = 5
    ),
    "undefined columns selected"
  )
})

# matrix()
test_that("detect_sni rejects matrix CharData input", {
  expect_error(
    detect_sni(
      CharData = matrix(),
      BandWidth = 5
    ),
    "subscript out of bounds"
  )
})


# ------------------------------------------------------------------------
# Input validation tests for BandWidth parameter
# ------------------------------------------------------------------------

# Test missing BandWidth parameter
test_that("detect_sni rejects missing BandWidth input", {
  # Setup test data
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

  pred_gam <-
    make_trend(
      data_source = data_source,
      sel_method = "non_linear"
    )

  # Create input for detect_sni
  char_data <- data.frame(
    data_source$Age,
    data_source$ROC,
    pred_gam
  )

  # Test expectation
  expect_error(
    detect_sni(
      CharData = char_data,
      BandWidth =
      ),
    'argument "BandWidth" is missing, with no default'
  )
})

# Test NULL BandWidth parameter
test_that("detect_sni rejects NULL BandWidth input", {
  # Setup test data
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

  pred_gam <-
    make_trend(
      data_source = data_source,
      sel_method = "non_linear"
    )

  # Create input for detect_sni
  char_data <- data.frame(
    data_source$Age,
    data_source$ROC,
    pred_gam
  )

  # Test expectation
  expect_error(
    detect_sni(
      CharData = char_data,
      BandWidth = NULL
    ),
    "argument is of length zero"
  )
})

# Test character BandWidth parameter
test_that("detect_sni rejects character BandWidth input", {
  # Setup test data
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

  pred_gam <-
    make_trend(
      data_source = data_source,
      sel_method = "non_linear"
    )

  # Create input for detect_sni
  char_data <- data.frame(
    data_source$Age,
    data_source$ROC,
    pred_gam
  )

  # Test expectation
  expect_error(
    detect_sni(
      CharData = char_data,
      BandWidth = "5"
    ),
    "non-numeric argument to binary operator"
  )
})

# Test negative BandWidth parameter
test_that("detect_sni throws error with negative BandWidth input", {
  # Setup test data
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

  pred_gam <-
    make_trend(
      data_source = data_source,
      sel_method = "non_linear"
    )

  # Create input for detect_sni
  char_data <- data.frame(
    data_source$Age,
    data_source$ROC,
    pred_gam
  )

  expect_error(
    detect_sni(
      CharData = char_data,
      BandWidth = -5
    ),
    "'f' must be finite and > 0"
  )
})

# Test zero BandWidth parameter
test_that("detect_sni handles zero BandWidth input", {
  # Setup test data
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

  pred_gam <-
    make_trend(
      data_source = data_source,
      sel_method = "non_linear"
    )

  # Create input for detect_sni
  char_data <- data.frame(
    data_source$Age,
    data_source$ROC,
    pred_gam
  )

  # Test expectation
  expect_error(
    detect_sni(
      CharData = char_data,
      BandWidth = 0
    ),
    "'f' must be finite and > 0"
  )
})

# ------------------------------------------------------------------------
# CharData format validation tests
# ------------------------------------------------------------------------

# Test for incorrect number of columns in CharData
test_that("detect_sni rejects CharData with incorrect number of columns", {
  # Setup test data
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

  # Create input for detect_sni with only 2 columns instead of required 3
  char_data <- data.frame(
    data_source$Age,
    data_source$ROC
  )

  # Test expectation
  expect_error(
    detect_sni(
      CharData = char_data,
      BandWidth = 5
    ),
    "undefined columns selected"
  )
})

# Test CharData with NAs in first column (ages)
test_that("detect_sni handles CharData with NAs in ages column", {
  # Setup test data
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

  pred_gam <-
    make_trend(
      data_source = data_source,
      sel_method = "non_linear"
    )

  # Create input for detect_sni with NA in age column
  char_data <- data.frame(
    data_source$Age,
    data_source$ROC,
    pred_gam
  )

  # Insert NA into first column (ages)
  char_data[5, 1] <- NA

  # Test expectation - this should cause an error or warning
  expect_error(
    suppressWarnings(
      detect_sni(
        CharData = char_data,
        BandWidth = 5 * mean(diff(data_source$Age), na.rm = TRUE)
      )
    ),
    "missing value where TRUE/FALSE needed"
  )
})

# Test CharData with NAs in second column (ROC values)
test_that("detect_sni handles CharData with NAs in ROC column", {
  # Setup test data
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

  pred_gam <-
    make_trend(
      data_source = data_source,
      sel_method = "non_linear"
    )

  # Create input for detect_sni with NA in ROC column
  char_data <- data.frame(
    data_source$Age,
    data_source$ROC,
    pred_gam
  )

  # Insert NA into second column (ROC)
  char_data[5, 2] <- NA

  # Test expectation - this should cause an error or warning
  expect_error(
    detect_sni(
      CharData = char_data,
      BandWidth = 5 * mean(diff(data_source$Age))
    ),
    # none programmed into function yet.
    # returns NA values in output instead of error
    # e.g., "Error: NA detected in ROC data"
  )
})

# Test CharData with NAs in third column (threshold values)
test_that("detect_sni handles CharData with NAs in threshold column", {
  # Setup test data
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

  pred_gam <-
    make_trend(
      data_source = data_source,
      sel_method = "non_linear"
    )

  # Create input for detect_sni with NA in threshold column
  char_data <- data.frame(
    data_source$Age,
    data_source$ROC,
    pred_gam
  )

  # Insert NA into third column (threshold)
  char_data[5, 3] <- NA

  # Test expectation - this should cause an error or warning
  expect_error(
    detect_sni(
      CharData = char_data,
      BandWidth = 5 * mean(diff(data_source$Age))
    ),
    # none programmed into function yet.
    # returns unexpected results (empty, NAs, numbers)
    # e.g., "Error: NA detected in threshold data"
  )
})

# ------------------------------------------------------------------------
# Output validation tests
# ------------------------------------------------------------------------

# Test output structure - list with expected elements
test_that("detect_sni returns a list with expected elements", {
  # Setup test data
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

  pred_gam <-
    make_trend(
      data_source = data_source,
      sel_method = "non_linear"
    )

  # Create input for detect_sni
  char_data <-
    data.frame(
      data_source$Age,
      data_source$ROC,
      pred_gam
    )

  # Get result
  result <-
    detect_sni(
      CharData = char_data,
      BandWidth = 5 * mean(diff(data_source$Age))
    )

  # Test expectation
  expect_named(
    result,
    c("SNI", "winInd", "popN", "popS", "meanN", "stdN", "CF")
  )
})

# Test that SNI is numeric vector of proper length
test_that("detect_sni returns SNI as numeric vector with proper length", {
  # Setup test data
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

  pred_gam <-
    make_trend(
      data_source = data_source,
      sel_method = "non_linear"
    )

  # Create input for detect_sni
  char_data <-
    data.frame(
      data_source$Age,
      data_source$ROC,
      pred_gam
    )

  # Get result
  result <-
    detect_sni(
      CharData = char_data,
      BandWidth = 5 * mean(diff(data_source$Age))
    )

  # Test expectation
  expect_type(
    result$SNI,
    "double"
  )
  expect_equal(
    length(result$SNI),
    nrow(char_data)
  )
})

# Test that winInd is a matrix with expected dimensions
test_that("detect_sni returns winInd as matrix with proper dimensions", {
  # Setup test data
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

  pred_gam <-
    make_trend(
      data_source = data_source,
      sel_method = "non_linear"
    )

  # Create input for detect_sni
  char_data <- data.frame(
    data_source$Age,
    data_source$ROC,
    pred_gam
  )

  # Get result
  result <-
    detect_sni(
      CharData = char_data,
      BandWidth = 5 * mean(diff(data_source$Age))
    )

  # Test expectation
  expect_true(
    is.matrix(result$winInd)
  )
  expect_equal(
    dim(result$winInd),
    c(nrow(char_data), 2)
  )
})

# Test that popN and popS are lists with expected length
test_that("detect_sni returns popN and popS as lists with proper length", {
  # Setup test data
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

  pred_gam <-
    make_trend(
      data_source = data_source,
      sel_method = "non_linear"
    )

  # Create input for detect_sni
  char_data <-
    data.frame(
      data_source$Age,
      data_source$ROC,
      pred_gam
    )

  # Get result
  result <-
    detect_sni(
      CharData = char_data,
      BandWidth = 5 * mean(diff(data_source$Age))
    )

  # Test expectation
  expect_type(
    result$popN, "list"
  )
  expect_type(
    result$popS, "list"
  )
  expect_equal(
    length(result$popN), nrow(char_data)
  )
  expect_equal(
    length(result$popS), nrow(char_data)
  )
})

# Test that meanN, stdN, and CF are numeric vectors with expected length
test_that("detect_sni returns meanN, stdN, CF as numeric vectors with proper length", {
  # Setup test data
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

  pred_gam <-
    make_trend(
      data_source = data_source,
      sel_method = "non_linear"
    )

  # Create input for detect_sni
  char_data <-
    data.frame(
      data_source$Age,
      data_source$ROC,
      pred_gam
    )

  # Get result
  result <-
    detect_sni(
      CharData = char_data,
      BandWidth = 5 * mean(diff(data_source$Age))
    )

  # Test expectation
  expect_type(
    result$meanN,
    "double"
  )
  expect_type(
    result$stdN,
    "double"
  )
  expect_type(
    result$CF,
    "double"
  )
  expect_equal(
    length(result$meanN),
    nrow(char_data)
  )
  expect_equal(
    length(result$stdN),
    nrow(char_data)
  )
  expect_equal(
    length(result$CF),
    nrow(char_data)
  )
})

# ------------------------------------------------------------------------
# Functionality tests with valid data
# ------------------------------------------------------------------------

# Test SNI values are within expected range
test_that("detect_sni produces SNI values within expected range", {
  # Setup test data
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

  pred_gam <-
    make_trend(
      data_source = data_source,
      sel_method = "non_linear"
    )

  # Create input for detect_sni
  char_data <-
    data.frame(
      data_source$Age,
      data_source$ROC,
      pred_gam
    )

  # Get result
  result <-
    detect_sni(
      CharData = char_data,
      BandWidth = 5 * mean(diff(data_source$Age))
    )

  # Test expectation - SNI should be finite numbers
  expect_true(
    all(
      is.finite(
        result$SNI
      )
    )
  )
})

# Test that changing bandwidth affects the result
test_that("detect_sni produces different results with different bandwidth", {
  # Setup test data
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

  pred_gam <-
    make_trend(
      data_source = data_source,
      sel_method = "non_linear"
    )

  # Create input for detect_sni
  char_data <-
    data.frame(
      data_source$Age,
      data_source$ROC,
      pred_gam
    )

  # Get results with two different bandwidths
  result1 <-
    detect_sni(
      CharData = char_data,
      BandWidth = 5 * mean(diff(data_source$Age))
    )

  result2 <- detect_sni(
    CharData = char_data,
    BandWidth = 10 * mean(diff(data_source$Age))
  )

  # Test expectation - results should be different
  expect_false(
    identical(
      result1$SNI,
      result2$SNI
    )
  )
})

# Valid input produces output is a list
test_that("detect_sni produces a list output", {
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
      verbose = FALSE
    )

  pred_gam <-
    make_trend(
      data_source = data_source,
      sel_method = "non_linear"
    )

  SNI_calc <-
    detect_sni(
      CharData =
        data.frame(
          data_source$Age,
          data_source$ROC,
          pred_gam
        ),
      BandWidth =
        5 * mean(diff(data_source$Age))
    )

  expect_type(
    SNI_calc,
    "list"
  )
})


test_that("detect_sni returns correct names in list elements", {
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
      verbose = FALSE
    )

  pred_gam <-
    make_trend(
      data_source = data_source,
      sel_method = "non_linear"
    )

  SNI_calc <-
    detect_sni(
      CharData =
        data.frame(
          data_source$Age,
          data_source$ROC,
          pred_gam
        ),
      BandWidth =
        5 * mean(diff(data_source$Age))
    )

  expect_named(
    SNI_calc,
    c("SNI", "winInd", "popN", "popS", "meanN", "stdN", "CF")
  )
})

test_that("detect_sni produces output with no NA values", {
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
      verbose = FALSE
    )

  pred_gam <-
    make_trend(
      data_source = data_source,
      sel_method = "non_linear"
    )

  SNI_calc <-
    detect_sni(
      CharData =
        data.frame(
          data_source$Age,
          data_source$ROC,
          pred_gam
        ),
      BandWidth =
        5 * mean(diff(data_source$Age))
    )

  expect_false(
    any(
      is.na(
        SNI_calc
      )
    )
  )
})
