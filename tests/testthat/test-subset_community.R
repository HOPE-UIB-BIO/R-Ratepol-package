# Input validation
test_that("subset_community fails with NULL arguments", {
  expect_error(
    subset_community(
      data_source = NULL,
      ommit_vars = NULL
    ),
    'to an object of class "NULL"'
  )
})

test_that("subset_community fails with NULL data_source input", {
  expect_error(
    subset_community(
      data_source = NULL,
      ommit_vars = c("label", "res_age", "age_diff", "age")
    ),
    'to an object of class "NULL"'
  )
})

test_that("subset_community fails with character data_source input", {
  expect_error(
    subset_community(
      data_source = "my_data",
      ommit_vars = c("label", "res_age", "age_diff", "age")
    ),
    'to an object of class "character"'
  )
})


test_that("subset_community fails with numeric data_source input", {
  expect_error(
    subset_community(
      data_source = 123,
      ommit_vars = c("label", "res_age", "age_diff", "age")
    ),
    "to an object of class" # double/numeric
  )
})

test_that("subset_community throws error with empty data_source input", {
  expect_warning(
    subset_community(
      data_source = data.frame(),
      ommit_vars = c("label", "res_age", "age_diff", "age")
    ),
    # returns empty result
    # no warning implemented yet
    # e.g.., "Warning: Empty data source supplied, returns empty result"
  )
})

test_that("subset_community fails with numeric data_source input", {
  expect_error(
    subset_community(
      data_source = list(),
      ommit_vars = c("label", "res_age", "age_diff", "age")
    ),
    "doesn't handle lists."
  )
})

test_that("subset_community fails with numeric data_source input", {
  expect_error(
    subset_community(
      data_source = matrix(),
      ommit_vars = c("label", "res_age", "age_diff", "age")
    ),
    "no applicable method for 'select' applied to an object of class."
  )
})

# ommit_vars modifications
test_that("subset_community returns result without dropping vars with NULL ommit_vars", {
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
    subset_community(
      data_source_reduce,
      ommit_vars = NULL
    )

  expect_named(
    res,
    names(data_source_reduce)
  )
})


test_that("subset_community correctly drops label if it's in ommit_vars", {
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
    subset_community(
      data_source_reduce,
      ommit_vars = c("label")
    )

  expect_false(
    "label" %in% names(res)
  )
})

test_that("subset_community throws error if invalid ommit_vars is entered", {
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
    res <-
      subset_community(
        data_source_reduce,
        ommit_vars = c("invalid_label")
      ),
    # none programmed into the function yet
  )
})

test_that("subset_community throws error if invalid ommit_vars is entered", {
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
    res <-
      subset_community(
        data_source_reduce,
        ommit_vars = 123
      ),
    # none programmed into the function yet
  )
})


# Output validation
test_that("subset_community works with valid data", {
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
    subset_community(
      data_source_reduce,
      ommit_vars = c("label", "res_age", "age_diff")
    )

  expect_named(
    res,
    setdiff(
      names(data_source_reduce),
      c("label", "res_age", "age_diff")
    )
  )
})


test_that("subset_community works with valid data", {
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
    subset_community(
      data_source_reduce
    )

  expect_named(
    res,
    setdiff(
      names(data_source_reduce),
      c("label", "res_age", "age_diff", "age")
    )
  )
})
