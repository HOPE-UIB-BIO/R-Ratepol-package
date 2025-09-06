# Comprehensive Test Suite for extract_data()
# v3 - Complete edge case coverage

# ----------------------------------------------------------- #
# 1. VERBOSE PARAMETER TESTING
# ----------------------------------------------------------- #

## 1.1a Test verbose parameter validation - character input
test_that(
  "extract_data() throws error if verbose is character",
  {
    expect_error(
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        verbose = "TRUE"
      ),
      "'verbose' must be one of the following: 'logical'"
    )
  }
)

## 1.1b Test verbose parameter validation - numeric input
test_that(
  "extract_data() throws error if verbose is numeric",
  {
    expect_error(
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        verbose = 1
      ),
      "'verbose' must be one of the following: 'logical'"
    )
  }
)

## 1.1c Test verbose parameter validation - NULL input
test_that(
  "extract_data() throws error if verbose is NULL",
  {
    expect_error(
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        verbose = NULL
      ),
      "'verbose' must be one of the following: 'logical'"
    )
  }
)

## 1.2a Test verbose = TRUE produces start message
test_that(
  "extract_data() produces start message when verbose = TRUE",
  {
    expect_message(
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        verbose = TRUE
      ),
      "Data extraction started"
    )
  }
)

## 1.2b Test verbose = TRUE produces completion message
test_that(
  "extract_data() produces completion message when verbose = TRUE",
  {
    expect_message(
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        verbose = TRUE
      ),
      "Data extraction completed"
    )
  }
)

## 1.3 Test verbose = FALSE produces no messages (when no issues)
test_that(
  "extract_data() produces no messages when verbose = FALSE and no data issues",
  {
    expect_silent(
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        verbose = FALSE
      )
    )
  }
)

# ----------------------------------------------------------- #
# 2. SAMPLE_ID DATA TYPE VALIDATION
# ----------------------------------------------------------- #

## 2.1 Test numeric sample_id throws error in community data
test_that(
  "extract_data() throws error if sample_id is numeric in community data",
  {
    community_numeric_id <-
      RRatepol::example_data$pollen_data[[1]]
    community_numeric_id$sample_id <-
      as.numeric(factor(community_numeric_id$sample_id))

    expect_error(
      extract_data(
        data_community_extract = community_numeric_id,
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        verbose = FALSE
      ),
      "Variable 'sample_id' in 'data_community' must.*be a 'character'"
    )
  }
)

## 2.2 Test factor sample_id throws error in community data
test_that(
  "extract_data() throws error if sample_id is factor in community data",
  {
    community_factor_id <-
      RRatepol::example_data$pollen_data[[1]]
    community_factor_id$sample_id <-
      as.factor(community_factor_id$sample_id)

    expect_error(
      extract_data(
        data_community_extract = community_factor_id,
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        verbose = FALSE
      ),
      "Variable 'sample_id' in 'data_community' must.*be a 'character'"
    )
  }
)

## 2.3 Test numeric sample_id in age data (should cause mismatch error)
test_that(
  "extract_data() throws error if sample_id is numeric in age data",
  {
    age_numeric_id <-
      RRatepol::example_data$sample_age[[1]]
    age_numeric_id$sample_id <-
      as.numeric(factor(age_numeric_id$sample_id))

    expect_error(
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = age_numeric_id,
        verbose = FALSE
      ),
      "Variable 'sample_id' must have same values in.*'data_age' and 'data_community'"
    )
  }
)

# ----------------------------------------------------------- #
# 3. AGE COLUMN VALIDATION
# ----------------------------------------------------------- #

## 3.1 Test missing age column in age data
test_that(
  "extract_data() throws error if age column is missing in age data",
  {
    age_no_age_col <-
      RRatepol::example_data$sample_age[[1]] %>%
      dplyr::select(-age)

    expect_error(
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = age_no_age_col,
        verbose = FALSE
      ),
      "Variable 'age' in 'data_source_age' must be a 'numeric'"
    )
  }
)

## 3.2a Test non-numeric age column - character
test_that(
  "extract_data() throws error if age column is character",
  {
    age_char_age <-
      RRatepol::example_data$sample_age[[1]]
    age_char_age$age <-
      as.character(age_char_age$age)

    expect_error(
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = age_char_age,
        verbose = FALSE
      ),
      "Variable 'age' in 'data_source_age' must be a 'numeric'"
    )
  }
)

## 3.2b Test non-numeric age column - factor
test_that(
  "extract_data() throws error if age column is factor",
  {
    age_factor_age <-
      RRatepol::example_data$sample_age[[1]]
    age_factor_age$age <-
      as.factor(age_factor_age$age)

    expect_error(
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = age_factor_age,
        verbose = FALSE
      ),
      "Variable 'age' in 'data_source_age' must be a 'numeric'"
    )
  }
)

# ----------------------------------------------------------- #
# 4. AGE UNCERTAINTY MATRIX VALIDATION
# ----------------------------------------------------------- #

## 4.1 Test age_uncertainty wrong dimensions (fewer columns than samples)
test_that(
  "extract_data() throws error if age_uncertainty has wrong number of columns",
  {
    age_un_wrong_size <-
      RRatepol::example_data$age_uncertainty[[1]][, 1:5] # reduce columns

    expect_error(
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = age_un_wrong_size,
        verbose = FALSE
      ),
      "Object 'data_source_age' and 'age_uncertainty' must have.*the same number of levels"
    )
  }
)

## 4.2 Test age_uncertainty wrong dimensions (more columns than samples)
test_that(
  "extract_data() throws error if age_uncertainty has too many columns",
  {
    age_un_extra <-
      cbind(
        RRatepol::example_data$age_uncertainty[[1]],
        extra_col = rnorm(nrow(RRatepol::example_data$age_uncertainty[[1]]))
      )

    expect_error(
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = age_un_extra,
        verbose = FALSE
      ),
      "Object 'data_source_age' and 'age_uncertainty' must have.*the same number of levels"
    )
  }
)

## 4.3 Test age_uncertainty as data.frame instead of matrix
test_that(
  "extract_data() throws error if age_uncertainty is data.frame instead of matrix",
  {
    age_un_df <-
      as.data.frame(RRatepol::example_data$age_uncertainty[[1]])

    expect_error(
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = age_un_df,
        verbose = FALSE
      ),
      "'age_uncertainty' must be one of the following: 'NULL', 'matrix'"
    )
  }
)

## 4.4 Test age_uncertainty as character vector
test_that(
  "extract_data() throws error if age_uncertainty is character vector",
  {
    expect_error(
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = "uncertainty_data",
        verbose = FALSE
      ),
      "'age_uncertainty' must be one of the following: 'NULL', 'matrix'"
    )
  }
)

# ----------------------------------------------------------- #
# 5. DATA DIMENSION MISMATCHES
# ----------------------------------------------------------- #

## 5.1 Test community data with more rows than age data
test_that(
  "extract_data() throws error if community has additional samples (rows)",
  {
    extra_sample_in_community <-
      rbind(
        RRatepol::example_data$pollen_data[[1]],
        RRatepol::example_data$pollen_data[[1]][1, ] # duplicate first row
      )

    extra_sample_in_community$sample_id[nrow(
      extra_sample_in_community
    )] <-
      "EXTRA_SAMPLE"

    expect_error(
      extract_data(
        data_community_extract = extra_sample_in_community,
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        verbose = FALSE
      ),
      "Variable 'sample_id' must have same values in.*'data_age' and 'data_community'"
    )
  }
)

## 5.2 Test age data with more rows than community data
test_that(
  "extract_data() throws error if age has more rows than community data",
  {
    extra_age <-
      rbind(
        RRatepol::example_data$sample_age[[1]],
        data.frame(sample_id = "EXTRA_SAMPLE", depth = NA, age = 999)
      )

    expect_error(
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = extra_age,
        verbose = FALSE
      ),
      "Variable 'sample_id' must have same values in.*'data_age' and 'data_community'"
    )
  }
)

# ----------------------------------------------------------- #
# 6. EMPTY AND MINIMAL DATA TESTING
# ----------------------------------------------------------- #

test_that(
  "extract_data() throws no error if data is empty",
  {
    zero_community <-
      RRatepol::example_data$pollen_data[[1]]
    zero_community[,-1] <- 0

    expect_no_error(
      extract_data(
        data_community_extract = zero_community,
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        verbose = FALSE
      )
    )
  }
)

## 6.1 Test empty community data frame throws no error
test_that(
  "extract_data() throws no error if data is empty",
  {
    empty_community <-
      data.frame(sample_id = character(0))
    empty_age <-
      data.frame(sample_id = character(0), depth = numeric(0), age = numeric(0))

    expect_no_error(
      extract_data(
        data_community_extract = empty_community,
        data_age_extract = empty_age,
        verbose = FALSE
      )
    )
  }
)

## 6.2a Test single row data returns correct type
test_that(
  "extract_data() handles single row data - returns list",
  {
    single_community <-
      RRatepol::example_data$pollen_data[[1]][1, ]
    single_age <-
      RRatepol::example_data$sample_age[[1]][1, ]

    result <-
      extract_data(
        data_community_extract = single_community,
        data_age_extract = single_age,
        verbose = FALSE
      )

    expect_type(result, "list")
  }
)

## 6.2b Test single row data returns correct names
test_that(
  "extract_data() handles single row data - returns named list",
  {
    single_community <-
      RRatepol::example_data$pollen_data[[1]][1, ]
    single_age <-
      RRatepol::example_data$sample_age[[1]][1, ]

    result <-
      extract_data(
        data_community_extract = single_community,
        data_age_extract = single_age,
        verbose = FALSE
      )

    expect_named(result, c("community", "age", "age_un"))
  }
)

## 6.2c Test single row data community has correct dimensions
test_that(
  "extract_data() handles single row data - community has 1 row",
  {
    single_community <-
      RRatepol::example_data$pollen_data[[1]][1, ]
    single_age <-
      RRatepol::example_data$sample_age[[1]][1, ]

    result <-
      extract_data(
        data_community_extract = single_community,
        data_age_extract = single_age,
        verbose = FALSE
      )

    expect_equal(nrow(result$community), 1)
  }
)

## 6.2d Test single row data age has correct dimensions
test_that(
  "extract_data() handles single row data - age has 1 row",
  {
    single_community <-
      RRatepol::example_data$pollen_data[[1]][1, ]
    single_age <-
      RRatepol::example_data$sample_age[[1]][1, ]

    result <-
      extract_data(
        data_community_extract = single_community,
        data_age_extract = single_age,
        verbose = FALSE
      )

    expect_equal(nrow(result$age), 1)
  }
)

## 6.3a Test community data with only sample_id column returns list
test_that(
  "extract_data() handles minimal community data - returns list",
  {
    minimal_community <-
      RRatepol::example_data$pollen_data[[1]][, 1, drop = FALSE]

    result <-
      extract_data(
        data_community_extract = minimal_community,
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        verbose = FALSE
      )

    expect_type(result, "list")
  }
)

## 6.3b Test community data with only sample_id column has zero columns after processing
test_that(
  "extract_data() handles minimal community data - community has 0 columns",
  {
    minimal_community <-
      RRatepol::example_data$pollen_data[[1]][, 1, drop = FALSE]

    result <-
      extract_data(
        data_community_extract = minimal_community,
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        verbose = FALSE
      )

    expect_equal(
      ncol(
        result$community
      ), 0
    ) # after removing sample_id column
  }
)

# ----------------------------------------------------------- #
# 7. DUPLICATE SAMPLE_ID TESTING
# ----------------------------------------------------------- #

## 7.1 Test duplicate sample_id in community data
test_that(
  "extract_data() handles duplicate sample_id in community data",
  {
    dup_community <-
      RRatepol::example_data$pollen_data[[1]]
    dup_community$sample_id[2] <-
      dup_community$sample_id[1] # create duplicate

    # This should cause issues when converting to rownames
    expect_error(
      extract_data(
        data_community_extract = dup_community,
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        verbose = FALSE
      )
    )
  }
)

## 7.2 Test duplicate sample_id in age data
test_that(
  "extract_data() handles duplicate sample_id in age data",
  {
    dup_age <-
      RRatepol::example_data$sample_age[[1]]
    dup_age$sample_id[2] <-
      dup_age$sample_id[1] # create duplicate

    # This should cause issues when converting to rownames
    expect_error(
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = dup_age,
        verbose = FALSE
      )
    )
  }
)

# ----------------------------------------------------------- #
# 8. AGE SORTING BEHAVIOR TESTING
# ----------------------------------------------------------- #

## 8.1 Test unsorted age data gets sorted
test_that(
  "extract_data() sorts unsorted age data",
  {
    unsorted_age <-
      RRatepol::example_data$sample_age[[1]]

    unsorted_age <-
      unsorted_age[order(
        runif(
          nrow(
            unsorted_age
          )
        )
      ), ] # randomize order

    rownames(unsorted_age) <-
      NULL

    result <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = unsorted_age,
        verbose = FALSE
      )

    # Check that result ages are sorted
    expect_true(is.unsorted(result$age$age) == FALSE)
  }
)

## 8.2 Test already sorted age data remains sorted
test_that(
  "extract_data() maintains sorted age data",
  {
    sorted_age <-
      RRatepol::example_data$sample_age[[1]] %>%
      dplyr::arrange(age)

    result <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = sorted_age,
        verbose = FALSE
      )

    expect_true(is.unsorted(result$age$age) == FALSE)
  }
)

# ----------------------------------------------------------- #
# 9. COMMUNITY DATA CONTENT VALIDATION
# ----------------------------------------------------------- #

## 9.1 Test community data with character columns returns error
test_that(
  "extract_data() handles community data with character columns - returns error",
  {
    char_community <-
      RRatepol::example_data$pollen_data[[1]] %>%
      sapply(., as.character) %>%
      as.data.frame()

    expect_error(
      extract_data(
        data_community_extract = char_community,
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        verbose = FALSE
      ),
      "'x' must be numeric"
    )
  }
)

## 9.2 Test community data with factor columns returns list
test_that(
  "extract_data() handles community data with factor columns - returns error",
  {
    factor_community <-
      RRatepol::example_data$pollen_data[[1]]
    factor_community$habitat_type <-
      as.factor(sample(c("forest", "grassland"), nrow(factor_community), replace = TRUE))

    expect_error(
      extract_data(
        data_community_extract = factor_community,
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        verbose = FALSE
      ),
      "'x' must be numeric"
    )
  }
)


# ----------------------------------------------------------- #
# 10. NA HANDLING EDGE CASES
# ----------------------------------------------------------- #

## 10.1a Test all community data is NA produces warning
test_that(
  "extract_data() handles all NA community data - produces message",
  {
    all_na_community <-
      RRatepol::example_data$pollen_data[[1]]

    all_na_community[, -1] <-
      NA # make all species columns NA, keep sample_id

    expect_message(
      extract_data(
        data_community_extract = all_na_community,
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        verbose = FALSE
      ),
      "Missing data have been detected in community data"
    )
  }
)

## 10.1b Test all community data is NA replaces with zeros and returns empty result
test_that(
  "extract_data() handles all NA community data - replaces with zeros & returns empty result",
  {
    all_na_community <-
      RRatepol::example_data$pollen_data[[1]]

    all_na_community[, -1] <-
      NA # make all species columns NA, keep sample_id

    result <-
      extract_data(
        data_community_extract = all_na_community,
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        verbose = FALSE
      )


    # Is result empty?
    expect_equal(nrow(result$community), 0)
    expect_equal(ncol(result$community), 0)
    expect_equal(nrow(result$age), 0)
    expect_equal(ncol(result$age), 1)
  }
)

## 10.3a Test mixed NA in age data produces warning
test_that(
  "extract_data() handles partial NA in age data - produces message",
  {
    partial_na_age <-
      RRatepol::example_data$sample_age[[1]]
    partial_na_age$age[1:3] <-
      NA

    expect_message(
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = partial_na_age,
        verbose = FALSE
      ),
      "Missing 'age' values"
    )
  }
)

## 10.3b Test for some NAs in age data - filters out NA rows
test_that(
  "extract_data() handles some NA in age data - filters out NA rows",
  {
    partial_na_age <-
      RRatepol::example_data$sample_age[[1]]

    partial_na_age$age[1:3] <-
      NA

    result <-
      suppressWarnings(
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = partial_na_age,
          verbose = FALSE
        )
      )

    n_na <-
      sum(is.na(partial_na_age$age))
    expect_equal(nrow(result$age), nrow(partial_na_age) - n_na)
  }
)

## 10.3c Test mixed NA in age data has no remaining NAs
test_that(
  "extract_data() handles partial NA in age data - no remaining NAs",
  {
    partial_na_age <-
      RRatepol::example_data$sample_age[[1]]

    partial_na_age$age[1:3] <-
      NA

    result <-
      suppressWarnings(
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = partial_na_age,
          verbose = FALSE
        )
      )

    expect_true(all(!is.na(result$age$age)))
  }
)

# ----------------------------------------------------------- #
# 11. TIBBLE VS DATA.FRAME HANDLING
# ----------------------------------------------------------- #

## 11.1a Test function works with tibbles - returns list
test_that(
  "extract_data() handles tibble input - returns list",
  {
    tibble_community <-
      tibble::as_tibble(RRatepol::example_data$pollen_data[[1]])
    tibble_age <-
      tibble::as_tibble(RRatepol::example_data$sample_age[[1]])

    result <-
      extract_data(
        data_community_extract = tibble_community,
        data_age_extract = tibble_age,
        verbose = FALSE
      )

    expect_type(result, "list")
  }
)

## 11.1b Test function works with tibbles - community is data.frame
test_that(
  "extract_data() handles tibble input - community is data.frame",
  {
    tibble_community <-
      tibble::as_tibble(RRatepol::example_data$pollen_data[[1]])
    tibble_age <-
      tibble::as_tibble(RRatepol::example_data$sample_age[[1]])

    result <-
      extract_data(
        data_community_extract = tibble_community,
        data_age_extract = tibble_age,
        verbose = FALSE
      )

    expect_s3_class(result$community, "data.frame")
  }
)

## 11.1c Test function works with tibbles - age is data.frame
test_that(
  "extract_data() handles tibble input - age is data.frame",
  {
    tibble_community <-
      tibble::as_tibble(RRatepol::example_data$pollen_data[[1]])
    tibble_age <-
      tibble::as_tibble(RRatepol::example_data$sample_age[[1]])

    result <-
      extract_data(
        data_community_extract = tibble_community,
        data_age_extract = tibble_age,
        verbose = FALSE
      )

    expect_s3_class(result$age, "data.frame")
  }
)

## 11.1d Test function works with tibbles - community not tibble
test_that(
  "extract_data() handles tibble input - community not tibble",
  {
    tibble_community <-
      tibble::as_tibble(RRatepol::example_data$pollen_data[[1]])
    tibble_age <-
      tibble::as_tibble(RRatepol::example_data$sample_age[[1]])

    result <-
      extract_data(
        data_community_extract = tibble_community,
        data_age_extract = tibble_age,
        verbose = FALSE
      )

    expect_false(inherits(result$community, "tbl"))
  }
)

## 11.1e Test function works with tibbles - age not tibble
test_that(
  "extract_data() handles tibble input - age not tibble",
  {
    tibble_community <-
      tibble::as_tibble(RRatepol::example_data$pollen_data[[1]])
    tibble_age <-
      tibble::as_tibble(RRatepol::example_data$sample_age[[1]])

    result <-
      extract_data(
        data_community_extract = tibble_community,
        data_age_extract = tibble_age,
        verbose = FALSE
      )

    expect_false(inherits(result$age, "tbl"))
  }
)

# ----------------------------------------------------------- #
# 12. OUTPUT STRUCTURE VALIDATION FOR EDGE CASES
# ----------------------------------------------------------- #

## 12.1a Test output structure with minimal valid data - returns list
test_that(
  "extract_data() returns correct structure with minimal data - returns list",
  {
    minimal_community <-
      RRatepol::example_data$pollen_data[[1]][3, 1:2]

    minimal_age <-
      RRatepol::example_data$sample_age[[1]][3, 1:3]

    result <-
      extract_data(
        data_community_extract = minimal_community,
        data_age_extract = minimal_age,
        verbose = FALSE
      )

    expect_type(result, "list")
  }
)

## 12.1b Test output structure with minimal valid data - named correctly
test_that(
  "extract_data() returns correct structure with minimal data - named correctly",
  {
    minimal_community <-
      RRatepol::example_data$pollen_data[[1]][3, 1:2]

    minimal_age <-
      RRatepol::example_data$sample_age[[1]][3, 1:3]

    result <-
      extract_data(
        data_community_extract = minimal_community,
        data_age_extract = minimal_age,
        verbose = FALSE
      )

    expect_named(result, c("community", "age", "age_un"))
  }
)

## 12.1c Test output structure with minimal valid data - community is data.frame
test_that(
  "extract_data() returns correct structure with minimal data - community is data.frame",
  {
    minimal_community <-
      RRatepol::example_data$pollen_data[[1]][3, 1:2]

    minimal_age <-
      RRatepol::example_data$sample_age[[1]][3, 1:3]

    result <-
      extract_data(
        data_community_extract = minimal_community,
        data_age_extract = minimal_age,
        verbose = FALSE
      )

    expect_s3_class(result$community, "data.frame")
  }
)

## 12.1d Test output structure with minimal valid data - age is data.frame
test_that(
  "extract_data() returns correct structure with minimal data - age is data.frame",
  {
    minimal_community <-
      RRatepol::example_data$pollen_data[[1]][3, 1:2]

    minimal_age <-
      RRatepol::example_data$sample_age[[1]][3, 1:3]

    result <-
      extract_data(
        data_community_extract = minimal_community,
        data_age_extract = minimal_age,
        verbose = FALSE
      )

    expect_s3_class(result$age, "data.frame")
  }
)

## 12.1e Test output structure with minimal valid data - age_un is NULL
test_that(
  "extract_data() returns correct structure with minimal data - age_un is NULL",
  {
    minimal_community <-
      RRatepol::example_data$pollen_data[[1]][3, 1:2]

    minimal_age <-
      RRatepol::example_data$sample_age[[1]][3, 1:3]

    result <-
      extract_data(
        data_community_extract = minimal_community,
        data_age_extract = minimal_age,
        verbose = FALSE
      )

    expect_null(result$age_un)
  }
)

## 12.1f Test output structure with minimal valid data - community rownames correct
test_that(
  "extract_data() returns correct structure with minimal data - community rownames",
  {
    minimal_community <-
      RRatepol::example_data$pollen_data[[1]][3, 1:2]

    minimal_age <-
      RRatepol::example_data$sample_age[[1]][3, 1:3]

    result <-
      extract_data(
        data_community_extract = minimal_community,
        data_age_extract = minimal_age,
        verbose = FALSE
      )

    expect_equal(rownames(result$community), minimal_community$sample_id)
  }
)

## 12.1g Test output structure with minimal valid data - age rownames correct
test_that(
  "extract_data() returns correct structure with minimal data - age rownames",
  {
    minimal_community <-
      RRatepol::example_data$pollen_data[[1]][3, 1:2]

    minimal_age <-
      RRatepol::example_data$sample_age[[1]][3, 1:3]

    result <-
      extract_data(
        data_community_extract = minimal_community,
        data_age_extract = minimal_age,
        verbose = FALSE
      )

    expect_equal(rownames(result$age), minimal_age$sample_id)
  }
)

## 12.1h Test output structure with minimal valid data - community values correct
test_that(
  "extract_data() returns correct structure with minimal data - community values",
  {
    minimal_community <-
      RRatepol::example_data$pollen_data[[1]][3, 1:2]

    minimal_age <-
      RRatepol::example_data$sample_age[[1]][3, 1:3]

    result <-
      extract_data(
        data_community_extract = minimal_community,
        data_age_extract = minimal_age,
        verbose = FALSE
      )

    expect_equal(result$community$`Chenopodiaceae/Amaranthaceae`, minimal_community[2]$`Chenopodiaceae/Amaranthaceae`)
  }
)

## 12.1i Test output structure with minimal valid data - age values correct
test_that(
  "extract_data() returns correct structure with minimal data - age values",
  {
    minimal_community <-
      RRatepol::example_data$pollen_data[[1]][3, 1:2]

    minimal_age <-
      RRatepol::example_data$sample_age[[1]][3, 1:3]

    result <-
      extract_data(
        data_community_extract = minimal_community,
        data_age_extract = minimal_age,
        verbose = FALSE
      )

    expect_equal(result$age$age, minimal_age$age)
  }
)


# ----------------------------------------------------------- #
# 14. SPECIFIC ERROR MESSAGE VALIDATION
# ----------------------------------------------------------- #

## 14.1a Test exact error messages for verbose parameter
test_that(
  "extract_data() produces exact expected error message for verbose parameter",
  {
    expect_error(
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        verbose = "yes"
      ),
      "'verbose' must be one of the following: 'logical'",
      fixed = TRUE
    )
  }
)

## 14.1b Test exact error messages for data_community_extract parameter
test_that(
  "extract_data() produces exact expected error message for data_community_extract parameter",
  {
    expect_error(
      extract_data(
        data_community_extract = "not_a_dataframe",
        data_age_extract = RRatepol::example_data$sample_age[[1]]
      ),
      "'data_community_extract' must be one of the following: 'data.frame'",
      fixed = TRUE
    )
  }
)

## 14.1c Test exact error messages for data_age_extract parameter
test_that(
  "extract_data() produces exact expected error message for data_age_extract parameter",
  {
    expect_error(
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = "not_a_dataframe"
      ),
      "'data_age_extract' must be one of the following: 'data.frame'",
      fixed = TRUE
    )
  }
)

## 14.2a Test column requirement error message for missing sample_id in community
test_that(
  "extract_data() produces exact error message for missing sample_id in community",
  {
    no_sample_id_community <-
      RRatepol::example_data$pollen_data[[1]] %>%
      dplyr::select(-sample_id)

    expect_error(
      extract_data(
        data_community_extract = no_sample_id_community,
        data_age_extract = RRatepol::example_data$sample_age[[1]]
      ),
      "'data_community_extract' must contains following columns: 'sample_id'",
      fixed = TRUE
    )
  }
)

## 14.2b Test column requirement error message for missing sample_id in age
test_that(
  "extract_data() produces exact error message for missing sample_id in age",
  {
    no_sample_id_age <-
      RRatepol::example_data$sample_age[[1]] %>%
      dplyr::select(-sample_id)

    expect_error(
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = no_sample_id_age
      ),
      "'data_age_extract' must contains following columns: 'sample_id'",
      fixed = TRUE
    )
  }
)

# ----------------------------------------------------------- #
# 15. EDGE CASES WITH AGE UNCERTAINTY
# ----------------------------------------------------------- #

## 15.1a Test age uncertainty with NA values returns list
test_that(
  "extract_data() handles age uncertainty matrix with NA values - returns list",
  {
    age_un_with_na <-
      RRatepol::example_data$age_uncertainty[[1]]
    age_un_with_na[1:5, 1:3] <-
      NA

    result <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = age_un_with_na,
        verbose = FALSE
      )

    expect_type(result, "list")
  }
)

## 15.1b Test age uncertainty with NA values returns data.frame
test_that(
  "extract_data() handles age uncertainty matrix with NA values - returns data.frame",
  {
    age_un_with_na <-
      RRatepol::example_data$age_uncertainty[[1]]
    age_un_with_na[1:5, 1:3] <-
      NA

    result <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = age_un_with_na,
        verbose = FALSE
      )

    expect_s3_class(result$age_un, "data.frame")
  }
)

## 15.1c Test age uncertainty with NA values preserves NAs
test_that(
  "extract_data() handles age uncertainty matrix with NA values - preserves NAs",
  {
    age_un_with_na <-
      RRatepol::example_data$age_uncertainty[[1]]
    age_un_with_na[1:5, 1:3] <-
      NA

    result <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = age_un_with_na,
        verbose = FALSE
      )

    expect_true(any(is.na(result$age_un)))
  }
)

## 15.2a Test age uncertainty with single row returns list
test_that(
  "extract_data() handles age uncertainty matrix with single row - returns list",
  {
    single_row_uncertainty <-
      matrix(
        rnorm(
          ncol(
            RRatepol::example_data$age_uncertainty[[1]]
          )
        ),
        nrow = 1
      )

    result <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = single_row_uncertainty,
        verbose = FALSE
      )

    expect_type(result, "list")
  }
)

## 15.2b Test age uncertainty with single row returns data.frame
test_that(
  "extract_data() handles age uncertainty matrix with single row - returns data.frame",
  {
    single_row_uncertainty <-
      matrix(
        rnorm(
          ncol(
            RRatepol::example_data$age_uncertainty[[1]]
          )
        ),
        nrow = 1
      )

    result <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = single_row_uncertainty,
        verbose = FALSE
      )

    expect_s3_class(result$age_un, "data.frame")
  }
)

## 15.2c Test age uncertainty with single row has correct dimensions
test_that(
  "extract_data() handles age uncertainty matrix with single row - correct dimensions",
  {
    single_row_uncertainty <-
      matrix(
        rnorm(
          ncol(
            RRatepol::example_data$age_uncertainty[[1]]
          )
        ),
        nrow = 1
      )

    result <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = single_row_uncertainty,
        verbose = FALSE
      )

    expect_equal(nrow(result$age_un), 1)
  }
)

## 15.3 Test age uncertainty column names assignment
test_that(
  "extract_data() correctly assigns column names to age uncertainty",
  {
    result <-
      extract_data(
        data_community_extract = RRatepol::example_data$pollen_data[[1]],
        data_age_extract = RRatepol::example_data$sample_age[[1]],
        age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
        verbose = FALSE
      )

    expected_names <-
      RRatepol::example_data$sample_age[[1]]$sample_id
    expect_identical(names(result$age_un), expected_names)
  }
)
