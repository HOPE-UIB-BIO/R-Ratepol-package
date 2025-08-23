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



