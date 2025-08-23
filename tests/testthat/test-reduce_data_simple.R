# uses only 1 parameter by default function call within run_iteration:
# result from subset_samples() "data_subset"
# other parameters:
## omit_vars = c("label", "res_age", "age_diff")
## check_taxa = TRUE
## check_levels = TRUE


test_that("reduce_data_simple throws error without input data", {
  expect_error(
    reduce_data_simple(
      data_source_reduce = NULL
    ),
    "no applicable method for 'select'"
  )
})

test_that("reduce_data_simple throws error without input data", {
  expect_error(
    reduce_data_simple(
      data_source_reduce = NULL,
      ommit_vars = c("label", "res_age", "age_diff"),
      check_taxa = TRUE,
      check_levels = TRUE
    ),
    "no applicable method for 'select'"
  )
})

test_that("reduce_data_simple throws error with NULL input data", {
  expect_error(
    reduce_data_simple(
      data_source_reduce = NULL,
      ommit_vars = c("label", "res_age", "age_diff"),
      check_taxa = TRUE,
      check_levels = TRUE
    ),
    "no applicable method for 'select'"
  )
})

test_that("reduce_data_simple throws error with character input data", {
  expect_error(
    reduce_data_simple(
      data_source_reduce = "my_data",
      ommit_vars = c("label", "res_age", "age_diff"),
      check_taxa = TRUE,
      check_levels = TRUE
    ),
    "no applicable method for 'select'"
  )
})

test_that("reduce_data_simple throws error with numeric input data", {
  expect_error(
    reduce_data_simple(
      data_source_reduce = 123,
      ommit_vars = c("label", "res_age", "age_diff"),
      check_taxa = TRUE,
      check_levels = TRUE
    ),
    "no applicable method for 'select'"
  )
})

test_that("reduce_data_simple throws error with list input data", {
  expect_error(
    reduce_data_simple(
      data_source_reduce = list(),
      ommit_vars = c("label", "res_age", "age_diff"),
      check_taxa = TRUE,
      check_levels = TRUE
    ),
    "doesn't handle lists"
  )
})

test_that("reduce_data_simple throws error with data.frame input data", {
  expect_error(
    reduce_data_simple(
      data_source_reduce = data.frame(),
      ommit_vars = c("label", "res_age", "age_diff"),
      check_taxa = TRUE,
      check_levels = TRUE
    ),
    # none programmed into the function yet
  )
})


# Check_taxa input validation

test_that("reduce_data_simple throws error with invalid check_taxa argument", {
  data_to_run_levels <-
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
      working_units = "levels",
      rand = 1
    ) %>%
    RUtilpol::flatten_list_by_one() %>%
    .[[1]]

  data_source_reduce <-
    subset_samples(
      data_source_subset = data_to_run_levels$data,
      data_source_bins = data_to_run_levels$bins,
      bin_selection = NULL
    )

  expect_error(
    reduce_data_simple(
      data_source_reduce = data_source_reduce,
      ommit_vars = c("label", "res_age", "age_diff"),
      check_taxa = "invalid",
      check_levels = TRUE
    ),
    # none programmed into the function yet
  )
})


# Check_levels input validation

test_that("reduce_data_simple throws error with invalid check_levels argument", {
  data_to_run_levels <-
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
      working_units = "levels",
      rand = 1
    ) %>%
    RUtilpol::flatten_list_by_one() %>%
    .[[1]]

  data_source_reduce <-
    subset_samples(
      data_source_subset = data_to_run_levels$data,
      data_source_bins = data_to_run_levels$bins,
      bin_selection = NULL
    )

  expect_error(
    reduce_data_simple(
      data_source_reduce = data_source_reduce,
      ommit_vars = c("label", "res_age", "age_diff"),
      check_taxa = TRUE,
      check_levels = "invalid"
    ),
    # none programmed into the function yet
  )
})


# Output validation:
# With valid data
test_that("reduce_data_simple works with valid data", {
  data_to_run_levels <-
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
      working_units = "levels",
      rand = 1
    ) %>%
    RUtilpol::flatten_list_by_one() %>%
    .[[1]]

  data_source_reduce <-
    subset_samples(
      data_source_subset = data_to_run_levels$data,
      data_source_bins = data_to_run_levels$bins,
      bin_selection = NULL
    )

  res <-
    reduce_data_simple(
      data_source_reduce = data_source_reduce,
      ommit_vars = c("label", "res_age", "age_diff"),
      check_taxa = TRUE,
      check_levels = TRUE
    )

  # ensure important columns are added back to data
  expect_true(
    all(
      c("label", "res_age", "age_diff") %in% names(res)
    )
  )
})


# Test that NA-rows/samples are dropped if check_levels = TRUE
test_that("reduce_data_simple drops samples without observations if check_levels = TRUE", {
  ## --- first bin --- ##
  set.seed(123)
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

  data_to_reduce <-
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = "first"
    )

  res <-
    reduce_data_simple(
      data_to_reduce,
      check_levels = TRUE,
      check_taxa = TRUE
    )

  expect_lt(
    nrow(res),
    nrow(data_to_reduce)
  )

  expect_true(
    all(
      rowSums(is.na(res)) == 0
    )
  )
})

# Test that NA-rows/samples are dropped if check_levels = FALSE
test_that("reduce_data_simple does not drop samples without observations if check_levels = FALSE", {
  ## --- first bin --- ##
  set.seed(123)
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

  data_to_reduce <-
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = "first"
    )

  res <-
    reduce_data_simple(
      data_to_reduce,
      check_levels = FALSE,
      check_taxa = TRUE
    )

  expect_equal(
    nrow(res),
    nrow(data_to_reduce)
  )

  # expect that there are still NAs
  expect_false(
    all(
      rowSums(is.na(res)) == 0
    )
  )
})

# Check_taxa = TRUE
test_that("reduce_data_simple drops taxa without observations if check_taxa = TRUE", {
  ## --- first bin --- ##
  set.seed(123)
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

  data_to_reduce <-
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = "first"
    )

  zero_taxa <-
    colSums(data_to_reduce[, -c(1:3)], na.rm = T)
  dropped_taxa <-
    names(zero_taxa[zero_taxa == 0])


  res <-
    reduce_data_simple(
      data_to_reduce,
      check_levels = FALSE,
      check_taxa = TRUE
    )

  # no more zero-taxa columns
  expect_false(
    any(
      dropped_taxa %in% names(res)
    )
  )

  # still NAs in rows
  expect_true(
    any(
      rowSums(is.na(res)) > 0
    )
  )
})

# Test that 0 taxa are dropped if check_taxa = FALSE
test_that("reduce_data_simple does not drop taxa without observations if check_taxa = FALSE", {
  ## --- first bin --- ##
  set.seed(123)
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

  data_to_reduce <-
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = "first"
    )
  zero_taxa <-
    colSums(data_to_reduce[, -c(1:3)], na.rm = T)
  dropped_taxa <-
    names(zero_taxa[zero_taxa == 0])

  res <-
    reduce_data_simple(
      data_to_reduce,
      check_levels = FALSE,
      check_taxa = FALSE
    )

  # expect that they are still there
  expect_true(
    all(
      dropped_taxa %in% names(res)
    )
  )
})


# --------------------------------------------------- #
#     Additional tests for reduce_data_simple         #
# --------------------------------------------------- #

# Additional Input Validation Tests

test_that("reduce_data_simple throws error with invalid ommit_vars argument", {
  data_to_run_levels <-
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
      working_units = "levels",
      rand = 1
    ) %>%
    RUtilpol::flatten_list_by_one() %>%
    .[[1]]

  data_source_reduce <-
    subset_samples(
      data_source_subset = data_to_run_levels$data,
      data_source_bins = data_to_run_levels$bins,
      bin_selection = NULL
    )

  expect_error(
    reduce_data_simple(
      data_source_reduce = data_source_reduce,
      ommit_vars = 123, # should be character vector
      check_taxa = TRUE,
      check_levels = TRUE
    )
  )
})

test_that("reduce_data_simple handles wrong ommit_vars", {
  data_to_run_levels <-
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
      working_units = "levels",
      rand = 1
    ) %>%
    RUtilpol::flatten_list_by_one() %>%
    .[[1]]

  data_source_reduce <-
    subset_samples(
      data_source_subset = data_to_run_levels$data,
      data_source_bins = data_to_run_levels$bins,
      bin_selection = NULL
    )

  expect_error(
    reduce_data_simple(
      data_source_reduce = data_source_reduce,
      ommit_vars = c("x", "y"),
      check_taxa = TRUE,
      check_levels = TRUE
    ),
    "'x' must be numeric"
  )
})


# Edge Cases and Corner Cases

test_that("reduce_data_simple throws warning if all taxa are zero", {
  data_to_run_levels <-
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
      working_units = "levels",
      rand = 1
    ) %>%
    RUtilpol::flatten_list_by_one() %>%
    .[[1]]

  data_source_reduce <-
    subset_samples(
      data_source_subset = data_to_run_levels$data,
      data_source_bins = data_to_run_levels$bins,
      bin_selection = NULL
    )

  # Set all community data to zero
  taxa_cols <-
    setdiff(names(data_source_reduce), c("label", "res_age", "age_diff"))
  data_source_reduce[, taxa_cols] <-
    0
  expect_warning(
    res <-
      reduce_data_simple(
        data_source_reduce = data_source_reduce,
        check_taxa = TRUE,
        check_levels = TRUE
      ),
    # none programmed into the function yet
    "Warning: Community data is all-zero. Return empty result."
  )
})

test_that("reduce_data_simple handles single row data", {
  data_to_run_levels <-
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
      working_units = "levels",
      rand = 1
    ) %>%
    RUtilpol::flatten_list_by_one() %>%
    .[[1]]

  data_source_reduce <-
    subset_samples(
      data_source_subset = data_to_run_levels$data,
      data_source_bins = data_to_run_levels$bins,
      bin_selection = NULL
    )[1, ] # Take only first row

  res <-
    reduce_data_simple(
      data_source_reduce = data_source_reduce,
      check_taxa = TRUE,
      check_levels = TRUE
    )

  expect_true(is.data.frame(res))
  expect_gte(nrow(res), 0) # Should handle single row gracefully
})

test_that("reduce_data_simple throws warning if only NA values in community data", {
  data_to_run_levels <-
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
      working_units = "levels",
      rand = 1
    ) %>%
    RUtilpol::flatten_list_by_one() %>%
    .[[1]]

  data_source_reduce <-
    subset_samples(
      data_source_subset = data_to_run_levels$data,
      data_source_bins = data_to_run_levels$bins,
      bin_selection = NULL
    )

  # Set all community data to NA
  taxa_cols <-
    setdiff(names(data_source_reduce), c("label", "res_age", "age_diff"))
  data_source_reduce[, taxa_cols] <-
    NA
  expect_warning(
    res <-
      reduce_data_simple(
        data_source_reduce = data_source_reduce,
        check_taxa = TRUE,
        check_levels = TRUE
      ),
    # none programmed into the function yet
    # e.g., "Warning: community data is all-NA. Returning empty result."
  )
})

# Parameter Combination Tests

test_that("reduce_data_simple works with both check_taxa and check_levels FALSE", {
  data_to_run_levels <-
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
      working_units = "levels",
      rand = 1
    ) %>%
    RUtilpol::flatten_list_by_one() %>%
    .[[1]]

  data_source_reduce <-
    subset_samples(
      data_source_subset = data_to_run_levels$data,
      data_source_bins = data_to_run_levels$bins,
      bin_selection = NULL
    )

  res <-
    reduce_data_simple(
      data_source_reduce = data_source_reduce,
      check_taxa = FALSE,
      check_levels = FALSE
    )

  # Should return identical to input
  expect_identical(
    res,
    data_source_reduce
  )
})

# Output Structure Validation
test_that("reduce_data_simple preserves data types", {
  data_to_run_levels <-
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
      working_units = "levels",
      rand = 1
    ) %>%
    RUtilpol::flatten_list_by_one() %>%
    .[[1]]

  data_source_reduce <-
    subset_samples(
      data_source_subset = data_to_run_levels$data,
      data_source_bins = data_to_run_levels$bins,
      bin_selection = NULL
    )

  res <-
    reduce_data_simple(
      data_source_reduce = data_source_reduce,
      check_taxa = TRUE,
      check_levels = TRUE
    )

  expect_true(is.data.frame(res))

  # Check that metadata columns maintain their types
  if ("label" %in% names(res)) {
    expect_true(is.character(res$label) || is.factor(res$label))
  }
  if ("res_age" %in% names(res)) {
    expect_true(is.numeric(res$res_age))
  }
  if ("age_diff" %in% names(res)) {
    expect_true(is.numeric(res$age_diff))
  }
})

# Integration with run_iteration context
test_that("reduce_data_simple works within run_iteration workflow", {
  # Test that the function works as expected within the context of run_iteration
  data_to_run_levels <-
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
      working_units = "levels",
      rand = 1
    ) %>%
    RUtilpol::flatten_list_by_one() %>%
    .[[1]]

  # Simulate the subset_samples call from run_iteration
  data_subset <-
    subset_samples(
      data_source_subset = data_to_run_levels$data,
      data_source_bins = data_to_run_levels$bins,
      bin_selection = "first"
    )

  # Call reduce_data_simple with default parameters as in run_iteration
  res <-
    reduce_data_simple(
      data_source_reduce = data_subset
    )

  expect_true(
    is.data.frame(res)
  )
  expect_true(
    all(
      c("label", "res_age", "age_diff") %in% names(res)
    )
  )
})
