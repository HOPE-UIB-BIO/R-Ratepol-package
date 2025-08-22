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
