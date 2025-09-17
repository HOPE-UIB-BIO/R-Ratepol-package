# ----------------------------------------------------------- #
# 1. INPUT PARAMETER TYPE VALIDATION -----
# ----------------------------------------------------------- #

## 1.1 Test required parameter validation - wrong type -----
test_that(
  "smooth_community_data() throws error if data_source_smooth is not list",
  {
    expect_error(
      smooth_community_data(
        data_source_smooth = "not_a_list",
        smooth_method = "m.avg",
        smooth_n_points = 5,
        verbose = FALSE
      ),
      "'data_source_smooth' must be a list"
    )
  }
)

## 1.2 Test required parameter validation - NULL input -----
test_that(
  "smooth_community_data() throws error if data_source_smooth is NULL",
  {
    expect_error(
      smooth_community_data(
        data_source_smooth = NULL,
        smooth_method = "m.avg",
        smooth_n_points = 5,
        verbose = FALSE
      ),
      "'data_source_smooth' must be a list"
    )
  }
)

## 1.3 Test optional parameter validation - wrong type -----
test_that(
  "smooth_community_data() throws error if smooth_method is not character",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    expect_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = 123,
        smooth_n_points = 5,
        verbose = FALSE
      ),
      "'smooth_method' must be character"
    )
  }
)

## 1.4 Test smooth_n_points type validation -----
test_that(
  "smooth_community_data() throws error if smooth_n_points is not numeric",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    expect_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "m.avg",
        smooth_n_points = "not_numeric",
        verbose = FALSE
      ),
      "'smooth_n_points' must be numeric"
    )
  }
)

## 1.5 Test smooth_n_max type validation -----
test_that(
  "smooth_community_data() throws error if smooth_n_max is not numeric",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    expect_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_n_max = "not_numeric",
        smooth_age_range = 500,
        verbose = FALSE
      ),
      "'smooth_n_max' must be numeric"
    )
  }
)

## 1.6 Test smooth_age_range type validation for grim -----
test_that(
  "smooth_community_data() throws error if smooth_age_range is not numeric for grim",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    expect_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_n_max = 9,
        smooth_age_range = "not_numeric",
        verbose = FALSE
      ),
      "'smooth_age_range' must be numeric"
    )
  }
)

## 1.7 Test smooth_age_range type validation for age.w -----
test_that(
  "smooth_community_data() throws error if smooth_age_range is not numeric for age.w",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    expect_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "age.w",
        smooth_n_points = 5,
        smooth_age_range = "not_numeric",
        verbose = FALSE
      ),
      "'smooth_age_range' must be numeric"
    )
  }
)

## 1.8 Test round_results type validation -----
test_that(
  "smooth_community_data() throws error if round_results is not logical",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    expect_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "m.avg",
        smooth_n_points = 5,
        round_results = "TRUE",
        verbose = FALSE
      ),
      "'round_results' must be logical and either TRUE or FALSE"
    )
  }
)

## 1.9 Test verbose type validation -----
test_that(
  "smooth_community_data() throws error if verbose is not logical",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    expect_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "m.avg",
        smooth_n_points = 5,
        verbose = "not_logical"
      ),
      "'verbose' must be logical and either TRUE or FALSE"
    )
  }
)

# ----------------------------------------------------------- #
# 2. INPUT PARAMETER VALUE VALIDATION -----
# ----------------------------------------------------------- #

## 2.1 Test parameter bounds - invalid smooth_method -----
test_that(
  "smooth_community_data() throws error if smooth_method is invalid",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    expect_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "invalid_method",
        smooth_n_points = 5,
        verbose = FALSE
      ),
      'should be one of "m.avg", "grim", "age.w", "shep"'
    )
  }
)

# ----------------------------------------------------------- #
# 3. INPUT DATA STRUCTURE VALIDATION -----
# ----------------------------------------------------------- #

## 3.1 Test required components presence -----
test_that(
  "smooth_community_data() throws error if community component is missing",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    invalid_data <-
      data_source_smooth
    invalid_data$community <-
      NULL

    expect_error(
      smooth_community_data(
        data_source_smooth = invalid_data,
        smooth_method = "m.avg",
        smooth_n_points = 5,
        verbose = FALSE
      ),
      "assert_that: length of assertion is not 1"
    )
  }
)

## 3.2 Test component data types -----
test_that(
  "smooth_community_data() throws warning if community component is a character matrix",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    invalid_data <-
      data_source_smooth
    invalid_data$community <-
      as.character(
        invalid_data$community
      )

    # Capture all warnings to verify the expected warning occurs
    expect_error(
      smooth_community_data(
        data_source_smooth = invalid_data,
        smooth_method = "m.avg",
        smooth_n_points = 5,
        verbose = FALSE
      ),
      "assert_that: length of assertion is not 1"
    )
  }
)

## 3.3 Test data dimensions - empty data -----
test_that(
  "smooth_community_data() handles empty community data",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    empty_data <-
      data_source_smooth
    empty_data$community <-
      data.frame()
    empty_data$age <-
      data.frame()

    expect_error(
      smooth_community_data(
        data_source_smooth = empty_data,
        smooth_method = "m.avg",
        smooth_n_points = 5,
        verbose = FALSE
      ),
      "'community' must have at least one row"
    )
  }
)

## 3.4 Test data dimensions - single row -----
test_that(
  "smooth_community_data() handles single row data",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    single_row_data <-
      data_source_smooth
    single_row_data$community <-
      single_row_data$community[1, , drop = FALSE]
    single_row_data$age <-
      single_row_data$age[1, , drop = FALSE]

    expect_no_error(
      smooth_community_data(
        data_source_smooth = single_row_data,
        smooth_method = "m.avg",
        smooth_n_points = 1,
        verbose = FALSE
      )
    )
  }
)

## 3.5 Test relationship between data inputs - dimension mismatch -----
test_that(
  "smooth_community_data() throws error when community and age have incompatible dimensions",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    mismatched_data <-
      data_source_smooth
    mismatched_data$age <-
      mismatched_data$age[1:5, , drop = FALSE]

    expect_error(
      smooth_community_data(
        data_source_smooth = mismatched_data,
        smooth_method = "m.avg",
        smooth_n_points = 5,
        verbose = FALSE
      ),
      "'community' and 'age' must have identical rownames"
    )
  }
)

# ----------------------------------------------------------- #
# 4. INPUT DATA CONTENT VALIDATION -----
# ----------------------------------------------------------- #

## 4.4 Test extreme values handling - very large
test_that(
  "smooth_community_data() handles very large values in community data",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    extreme_data <-
      data_source_smooth
    extreme_data$community[1, 1] <-
      1e6

    expect_no_error(
      smooth_community_data(
        data_source_smooth = extreme_data,
        smooth_method = "m.avg",
        smooth_n_points = 5,
        verbose = FALSE
      )
    )
  }
)

## 4.5 Test extreme values handling - very small -----
test_that(
  "smooth_community_data() handles very small values in community data",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    extreme_data <-
      data_source_smooth
    extreme_data$community[1, 1] <-
      1e-6

    expect_no_error(
      smooth_community_data(
        data_source_smooth = extreme_data,
        smooth_method = "m.avg",
        smooth_n_points = 5,
        verbose = FALSE
      )
    )
  }
)

# ----------------------------------------------------------- #
# 6. BOOLEAN/LOGICAL PARAMETER TESTING
# ----------------------------------------------------------- #

## 6.1 Test boolean parameter - TRUE behavior
test_that(
  "smooth_community_data() behaves correctly when round_results = TRUE",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    result <-
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "m.avg",
        smooth_n_points = 5,
        round_results = TRUE,
        verbose = FALSE
      )

    expect_true(
      all(
        result$community == round(
          result$community
        )
      )
    )
  }
)

## 6.2 Test boolean parameter - FALSE behavior -----
test_that(
  "smooth_community_data() behaves correctly when round_results = FALSE",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    result <-
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "m.avg",
        smooth_n_points = 5,
        round_results = FALSE,
        verbose = FALSE
      )
    expect_false(
      all(
        result$community == round(
          result$community
        )
      )
    )
  }
)

## 6.3 Test verbose parameter - produces expected message -----
test_that(
  "smooth_community_data() produces expected message when verbose = TRUE",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    expect_message(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "m.avg",
        smooth_n_points = 5,
        verbose = TRUE
      ),
      "Data will be smoothed by 'moving average'"
    )
  }
)

## 6.4 Test verbose parameter - silent when FALSE -----
test_that(
  "smooth_community_data() produces no messages when verbose = FALSE",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    expect_silent(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "m.avg",
        smooth_n_points = 5,
        verbose = FALSE
      )
    )
  }
)

# ----------------------------------------------------------- #
# 7. OUTPUT STRUCTURE VALIDATION ----
# ----------------------------------------------------------- #

## 7.1 Test output type -----
test_that(
  "smooth_community_data() returns correct output type",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    result <-
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "m.avg",
        smooth_n_points = 5,
        verbose = FALSE
      )

    expect_type(
      result,
      "list"
    )
  }
)

## 7.2 Test output class -----
test_that(
  "smooth_community_data() returns community component as data.frame",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    result <-
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "m.avg",
        smooth_n_points = 5,
        verbose = FALSE
      )

    expect_s3_class(
      result$community,
      "data.frame"
    )
  }
)

## 7.3 Test output names/structure -----
test_that(
  "smooth_community_data() returns correctly named output",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    result <-
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "m.avg",
        smooth_n_points = 5,
        verbose = FALSE
      )

    expect_named(
      result,
      c(
        "community",
        "age",
        "age_un"
      )
    )
  }
)

## 7.4 Test output dimensions -----
test_that(
  "smooth_community_data() returns output with correct dimensions",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    result <-
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "m.avg",
        smooth_n_points = 5,
        verbose = FALSE
      )

    expect_equal(
      nrow(
        result$community
      ),
      nrow(data_source_smooth$community)
    )
  }
)

test_that(
  "smooth_community_data() returns output with correct column count",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    result <-
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "m.avg",
        smooth_n_points = 5,
        verbose = FALSE
      )

    expect_equal(
      ncol(
        result$community
      ),
      ncol(data_source_smooth$community)
    )
  }
)

## 7.5 Test output component types -----
test_that(
  "smooth_community_data() returns age component with correct type",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    result <-
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "m.avg",
        smooth_n_points = 5,
        verbose = FALSE
      )

    expect_s3_class(
      result$age,
      "data.frame"
    )
  }
)

# ----------------------------------------------------------- #
# 8. OUTPUT CONTENT VALIDATION -----
# ----------------------------------------------------------- #

## 8.1 Test output values - value ranges -----
test_that(
  "smooth_community_data() produces non-negative values",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    result <-
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "m.avg",
        smooth_n_points = 5,
        verbose = FALSE
      )

    expect_true(
      all(
        result$community >= 0
      )
    )
  }
)

## 8.2 Test output values - no missing values -----
test_that(
  "smooth_community_data() produces no missing values in result",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    result <-
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "m.avg",
        smooth_n_points = 5,
        verbose = FALSE
      )

    expect_true(
      all(
        !is.na(
          result$community
        )
      )
    )
  }
)

## 8.3 Test output consistency - reproducible results -----
test_that(
  "smooth_community_data() produces consistent results across calls",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    result1 <-
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "m.avg",
        smooth_n_points = 5,
        verbose = FALSE
      )

    result2 <-
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "m.avg",
        smooth_n_points = 5,
        verbose = FALSE
      )

    expect_true(
      identical(
        result1,
        result2
      )
    )
  }
)

# ----------------------------------------------------------- #
# 9. AGE UNCERTAINTY INTEGRATION -----
# ----------------------------------------------------------- #

## 9.1 Test with age uncertainty preservation -----
test_that(
  "smooth_community_data() preserves age uncertainty structure",
  {
    extracted_data <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
      )

    result <-
      smooth_community_data(
        data_source_smooth = extracted_data,
        smooth_method = "m.avg",
        smooth_n_points = 5,
        verbose = FALSE
      )

    # Verify age uncertainty is preserved
    expect_false(
      is.null(
        result$age_un
      )
    )
  }
)
# ----------------------------------------------------------- #
# 10. MOVING AVERAGE METHOD TESTS -----
# ----------------------------------------------------------- #

## 10.1 Test m.avg method basic functionality -----
test_that(
  "smooth_community_data() works with m.avg method",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    expect_no_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "m.avg",
        smooth_n_points = 5,
        verbose = FALSE
      )
    )
  }
)

## 10.2 Test m.avg parameter bounds - even smooth_n_points -----
test_that(
  "smooth_community_data() throws error if smooth_n_points is even for m.avg",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    expect_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "m.avg",
        smooth_n_points = 4,
        verbose = FALSE
      ),
      "'smooth_n_points' must be odd"
    )
  }
)

## 10.3 Test m.avg parameter bounds - at boundaries (valid) -----
test_that(
  "smooth_community_data() accepts smooth_n_points = 1 for m.avg",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    expect_no_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "m.avg",
        smooth_n_points = 1,
        verbose = FALSE
      )
    )
  }
)

## 10.4 Test m.avg with uniform data -----
test_that(
  "smooth_community_data() produces identical values for uniform data with m.avg",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    # Create uniform data
    uniform_data <-
      data_source_smooth
    uniform_data$community[, 1] <-
      10

    result <-
      smooth_community_data(
        data_source_smooth = uniform_data,
        smooth_method = "m.avg",
        smooth_n_points = 5,
        verbose = FALSE
      )

    expect_true(
      all(
        result$community[, 1] == 10
      )
    )
  }
)

## 10.5 Test m.avg ignores optional parameters -----
test_that(
  "smooth_community_data() works with m.avg ignoring optional parameters",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    expect_no_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "m.avg",
        smooth_n_points = 5,
        smooth_n_max = 9,
        smooth_age_range = 500,
        verbose = FALSE
      )
    )
  }
)

## 10.6 Test warning when smooth_n_points exceeds available data -----
test_that(
  "smooth_community_data() throws error if smooth_n_points larger than available samples",
  {
    # Create minimal dataset with only 3 samples
    minimal_data <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]][1:3, ],
        RRatepol::example_data$sample_age[[1]][1:3, ]
      )

    expect_error(
      smooth_community_data(
        data_source_smooth = minimal_data,
        smooth_method = "m.avg",
        smooth_n_points = 7, # More than available samples
        verbose = FALSE
      ),
      "'smooth_n_points' must be <= number of samples in 'community'"
    )
  }
)

# ----------------------------------------------------------- #
# 11. GRIM METHOD TESTS -----
# ----------------------------------------------------------- #

## 11.1 Test grim method basic functionality -----
test_that(
  "smooth_community_data() works with grim method (no error)",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    expect_no_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_n_max = 9,
        smooth_age_range = 500,
        verbose = FALSE
      )
    )
  }
)

## 11.2 Test grim parameter bounds - even smooth_n_points -----
test_that(
  "smooth_community_data() throws error if smooth_n_points is even for grim",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    expect_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "grim",
        smooth_n_points = 4,
        smooth_n_max = 8,
        smooth_age_range = 500,
        verbose = FALSE
      ),
      "'smooth_n_points' must be odd"
    )
  }
)

## 11.4 Test grim smooth_n_max must be bigger than smooth_n_points -----
test_that(
  "smooth_community_data() throws error if grim smooth_n_max <= smooth_n_points",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    expect_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_n_max = 5,
        smooth_age_range = 500,
        verbose = FALSE
      ),
      "'smooth_n_points' must be < 'smooth_n_max' for 'grim' smoothing"
    )
  }
)

## 11.5 Test grim with very small age range -----
test_that(
  "smooth_community_data() handles very small smooth_age_range for grim",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    expect_no_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "grim",
        smooth_n_points = 3,
        smooth_n_max = 5,
        smooth_age_range = 1, # Very small range
        verbose = FALSE
      )
    )
  }
)

## 11.6 Test with very large smooth_n_max -----
test_that(
  "smooth_community_data() throws error with large smooth_n_max values",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    expect_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_n_max = 65, # Large value
        smooth_age_range = 5000,
        verbose = FALSE
      ),
      "'smooth_n_max' must be <= number of samples in 'community'"
    )
  }
)

## 11.7 Test grim minimum data requirements -----
test_that(
  "smooth_community_data() works with grim minimum data requirements",
  {
    minimal_5 <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]][1:5, ],
        RRatepol::example_data$sample_age[[1]][1:5, ]
      )

    expect_no_error(
      smooth_community_data(
        data_source_smooth = minimal_5,
        smooth_method = "grim",
        smooth_n_points = 3,
        smooth_n_max = 5,
        smooth_age_range = 100,
        verbose = FALSE
      )
    )
  }
)

## 11.8 Test exact error messages - grim specific validation -----
test_that(
  "smooth_community_data() produces exact expected error message for grim smooth_n_max",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    expect_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "grim",
        smooth_n_points = 5,
        smooth_n_max = 5,
        smooth_age_range = 500,
        verbose = FALSE
      ),
      "'smooth_n_points' must be < 'smooth_n_max' for 'grim' smoothing"
    )
  }
)

# ----------------------------------------------------------- #
# 12. AGE WEIGHTED METHOD TESTS -----
# ----------------------------------------------------------- #

## 12.1 Test age.w method basic functionality -----
test_that(
  "smooth_community_data() works with age.w method (no error)",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    expect_no_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "age.w",
        smooth_n_points = 5,
        smooth_age_range = 500,
        verbose = FALSE
      )
    )
  }
)

## 12.2 Test age.w parameter bounds - even smooth_n_points -----
test_that(
  "smooth_community_data() throws error if smooth_n_points is even for age.w",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    expect_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "age.w",
        smooth_n_points = 4,
        smooth_age_range = 500,
        verbose = FALSE
      ),
      "'smooth_n_points' must be odd"
    )
  }
)

## 12.3 Test age.w with non-sequential age data -----
test_that(
  "smooth_community_data() handles non-sequential age data for age.w method",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]][1:10, ],
        RRatepol::example_data$sample_age[[1]][1:10, ]
      )

    # Scramble age order
    scrambled_data <-
      data_source_smooth
    scrambled_data$age$age <-
      sample(
        scrambled_data$age$age
      )

    expect_error(
      smooth_community_data(
        data_source_smooth = scrambled_data,
        smooth_method = "age.w",
        smooth_n_points = 5,
        smooth_age_range = 500,
        verbose = FALSE
      ),
      "'age' must be sorted in increasing order"
    )
  }
)


# ----------------------------------------------------------- #
# 13. SHEPARD METHOD TESTS -----
# ----------------------------------------------------------- #

## 13.1 Test shep method basic functionality -----
test_that(
  "smooth_community_data() works with shep method (no error)",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    expect_no_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "shep",
        smooth_n_points = 5,
        verbose = FALSE
      )
    )
  }
)

## 13.2 Test shep with even numbers -----
test_that(
  "smooth_community_data() accepts even smooth_n_points for shep",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    expect_no_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "shep",
        smooth_n_points = 4,
        verbose = FALSE
      )
    )
  }
)

## 13.3 Test shep minimum valid smooth_n_points -----
test_that(
  "smooth_community_data() handles shep method with minimum smooth_n_points requirements",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    # Test that shep fails appropriately with smooth_n_points = 2
    expect_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "shep",
        smooth_n_points = 2,
        verbose = FALSE
      ),
      "'smooth_n_points' must be > 2 for 'shep' smoothing"
    )

    # Test that shep works with smooth_n_points >= 3
    expect_no_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "shep",
        smooth_n_points = 3,
        verbose = FALSE
      )
    )
  }
)

## 13.4 Test shep works without optional parameters -----
test_that(
  "smooth_community_data() works with shep without additional parameters",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    expect_no_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "shep",
        smooth_n_points = 5
      )
    )
  }
)


## 13.6 Test shep minimum data requirements
test_that(
  "smooth_community_data() works with shep minimum data requirements",
  {
    minimal_3 <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]][1:3, ],
        RRatepol::example_data$sample_age[[1]][1:3, ]
      )

    expect_no_error(
      smooth_community_data(
        data_source_smooth = minimal_3,
        smooth_method = "shep",
        smooth_n_points = 3,
        verbose = FALSE
      )
    )
  }
)

# ----------------------------------------------------------- #
# 14. EDGE CASES AND BOUNDARY CONDITIONS -----
# ----------------------------------------------------------- #

## 14.1 Test minimum valid input -----
test_that(
  "smooth_community_data() handles minimum valid input",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]][1:3, ],
        RRatepol::example_data$sample_age[[1]][1:3, ]
      )

    expect_no_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "m.avg",
        smooth_n_points = 1,
        verbose = FALSE
      )
    )
  }
)

## 14.2 Test exact error messages - parameter validation -----
test_that(
  "smooth_community_data() produces exact expected error message for even smooth_n_points",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]]
      )

    expect_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = "grim",
        smooth_n_points = 4,
        smooth_n_max = 8,
        smooth_age_range = 500,
        verbose = FALSE
      ),
      "'smooth_n_points' must be odd",
      fixed = TRUE
    )
  }
)

# multiple smooth_methods
test_that(
  "smooth_community_data fails if more than one smooth_method is supplied",
  {
    data_source_smooth <-
      extract_data(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$sample_age[[1]],
        RRatepol::example_data$age_uncertainty[[1]]
      )

    expect_error(
      smooth_community_data(
        data_source_smooth = data_source_smooth,
        smooth_method = c("shep", "m.avg"),
        smooth_n_points = 5
      ),
      "'arg' must be of length 1"
    )
  }
)
