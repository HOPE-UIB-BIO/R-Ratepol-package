# with input data for all working_unit methods.
# working_units = c("levels", "bins", "MW")

# different workflow for "levels" than for "bins" and "MW"
# levels ignores bin_selection.

# I will test without smoothing to reduce the amount of code and focus on
# issues that are more likely.... (?)


# ---------------------------------------------------- #
# Input validation : Default parameters #
# ---------------------------------------------------- #
test_that("subset_samples throws error with no data input", {
  expect_error(
    subset_samples(),
    'argument "data_source_bins" is missing, with no default'
  )
})

test_that("subset_samples throws error with NULL data input", {
  expect_error(
    subset_samples(
      data_source_subset = NULL,
      data_source_bins = NULL,
      bin_selection = NULL
    ),
    "argument of length 0"
  )
})

# ---- levels ---- #
# subset-data NULL
test_that("subset_samples and levels throws error with NULL data_subset input", {
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

  data_source_bins <-
    data_to_run_levels$bins

  expect_error(
    subset_samples(
      data_source_subset = NULL,
      data_source_bins = data_source_bins,
      bin_selection = NULL
    ),
    "is not TRUE"
  )
})

test_that("subset_samples and levels throws error with numeric data_subset input", {
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

  data_source_bins <-
    data_to_run_levels$bins

  expect_error(
    subset_samples(
      data_source_subset = 123,
      data_source_bins = data_source_bins,
      bin_selection = NULL
    ),
    "is not TRUE"
  )
})

test_that("subset_samples and levels throws error with character data_subset input", {
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

  data_source_bins <-
    data_to_run_levels$bins

  expect_error(
    subset_samples(
      data_source_subset = "my_data",
      data_source_bins = data_source_bins,
      bin_selection = NULL
    ),
    "is not TRUE"
  )
})

test_that("subset_samples and levels throws error with list data_subset input", {
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

  data_source_subset <-
    data_to_run_levels$data

  data_source_bins <-
    data_to_run_levels$bins

  expect_error(
    subset_samples(
      data_source_subset = list(data_source_subset),
      data_source_bins = data_source_bins,
      bin_selection = NULL
    ),
    "is not TRUE"
  )
})

# bins-data = NULL
test_that("subset_samples and levels throws error with NULL data_bins input", {
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

  data_source_subset <-
    data_to_run_levels$data

  expect_error(
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = NULL,
      bin_selection = NULL
    ),
    "non-numeric matrix extent"
  )
})

test_that("subset_samples and levels throws error with numeric data_bins input", {
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

  data_source_subset <-
    data_to_run_levels$data

  expect_error(
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = 123,
      bin_selection = NULL
    ),
    "operator is invalid for atomic vectors"
  )
})

test_that("subset_samples and levels throws error with character data_bins input", {
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

  data_source_subset <-
    data_to_run_levels$data

  expect_error(
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = "my_bins",
      bin_selection = NULL
    ),
    "operator is invalid for atomic vectors"
  )
})

test_that("subset_samples and levels throws error with list data_bins input", {
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

  data_source_subset <-
    data_to_run_levels$data

  data_source_bins <-
    data_to_run_levels$bins

  expect_error(
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = list(data_source_bins),
      bin_selection = NULL
    ),
    "non-numeric matrix extent"
  )
})

# ---- bins ---- #
# subset-data NULL
test_that("subset_samples and bins throws error with NULL data_subset input", {
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

  data_source_bins <-
    data_to_run_bins$bins

  expect_error(
    subset_samples(
      data_source_subset = NULL,
      data_source_bins = data_source_bins,
      bin_selection = NULL
    ),
    "argument of length 0"
  )
})

test_that("subset_samples and bins throws error with numeric data_subset input", {
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

  data_source_bins <-
    data_to_run_bins$bins

  expect_error(
    subset_samples(
      data_source_subset = 123,
      data_source_bins = data_source_bins,
      bin_selection = NULL
    ),
    "argument of length 0"
  )
})

test_that("subset_samples and bins throws error with character data_subset input", {
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

  data_source_bins <-
    data_to_run_bins$bins

  expect_error(
    subset_samples(
      data_source_subset = "my_data",
      data_source_bins = data_source_bins,
      bin_selection = NULL
    ),
    "argument of length 0"
  )
})

test_that("subset_samples and bins throws error with list data_subset input", {
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

  data_source_bins <-
    data_to_run_bins$bins

  data_source_subset <-
    data_to_run_bins$data

  expect_error(
    subset_samples(
      data_source_subset = list(data_source_subset),
      data_source_bins = data_source_bins,
      bin_selection = NULL
    ),
    "argument of length 0"
  )
})

# bins-data = NULL
test_that("subset_samples and bins throws error with NULL data_bins input", {
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

  data_source_bins <-
    data_to_run_bins$bins

  data_source_subset <-
    data_to_run_bins$data

  expect_error(
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = NULL,
      bin_selection = NULL
    ),
    "non-numeric matrix extent"
  )
})

test_that("subset_samples and bins throws error with numeric data_bins input", {
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

  data_source_bins <-
    data_to_run_bins$bins

  data_source_subset <-
    data_to_run_bins$data

  expect_error(
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = 123,
      bin_selection = NULL
    ),
    "operator is invalid for atomic vectors"
  )
})

test_that("subset_samples and bins throws error with character data_bins input", {
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

  data_source_bins <-
    data_to_run_bins$bins

  data_source_subset <-
    data_to_run_bins$data

  expect_error(
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = "my_bins",
      bin_selection = NULL
    ),
    "operator is invalid for atomic vectors"
  )
})

test_that("subset_samples and bins throws error with list data_bins input", {
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

  data_source_bins <-
    data_to_run_bins$bins

  data_source_subset <-
    data_to_run_bins$data

  expect_error(
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = list(data_source_bins),
      bin_selection = NULL
    ),
    "non-numeric matrix extent"
  )
})

# ---- MW ---- #
# subset-data NULL
test_that("subset_samples and MW throws error with NULL data_subset input", {
  data_to_run_MW <-
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
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 5,
      rand = 1
    ) %>%
    RUtilpol::flatten_list_by_one() %>%
    .[[1]]

  data_source_bins <-
    data_to_run_MW$bins

  expect_error(
    subset_samples(
      data_source_subset = NULL,
      data_source_bins = data_source_bins,
      bin_selection = NULL
    ),
    "argument of length 0"
  )
})

test_that("subset_samples and MW throws error with numeric data_subset input", {
  data_to_run_MW <-
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
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 5,
      rand = 1
    ) %>%
    RUtilpol::flatten_list_by_one() %>%
    .[[1]]

  data_source_bins <-
    data_to_run_MW$bins

  expect_error(
    subset_samples(
      data_source_subset = 123,
      data_source_bins = data_source_bins,
      bin_selection = NULL
    ),
    "argument of length 0"
  )
})

test_that("subset_samples and MW throws error with character data_subset input", {
  data_to_run_MW <-
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
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 5,
      rand = 1
    ) %>%
    RUtilpol::flatten_list_by_one() %>%
    .[[1]]

  data_source_bins <-
    data_to_run_MW$bins

  expect_error(
    subset_samples(
      data_source_subset = "my_data",
      data_source_bins = data_source_bins,
      bin_selection = NULL
    ),
    "argument of length 0"
  )
})

test_that("subset_samples and MW throws error with list data_subset input", {
  data_to_run_MW <-
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
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 5,
      rand = 1
    ) %>%
    RUtilpol::flatten_list_by_one() %>%
    .[[1]]

  data_source_bins <-
    data_to_run_MW$bins

  data_source_subset <-
    data_to_run_MW$data

  expect_error(
    subset_samples(
      data_source_subset = list(data_source_subset),
      data_source_bins = data_source_bins,
      bin_selection = NULL
    ),
    "argument of length 0"
  )
})

# bins-data = NULL
test_that("subset_samples and MW throws error with NULL data_bins input", {
  data_to_run_MW <-
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
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 5,
      rand = 1
    ) %>%
    RUtilpol::flatten_list_by_one() %>%
    .[[1]]

  data_source_bins <-
    data_to_run_MW$bins

  data_source_subset <-
    data_to_run_MW$data

  expect_error(
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = NULL,
      bin_selection = NULL
    ),
    "non-numeric matrix extent"
  )
})

test_that("subset_samples and MW throws error with numeric data_bins input", {
  data_to_run_MW <-
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
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 5,
      rand = 1
    ) %>%
    RUtilpol::flatten_list_by_one() %>%
    .[[1]]

  data_source_bins <-
    data_to_run_MW$bins

  data_source_subset <-
    data_to_run_MW$data

  expect_error(
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = 123,
      bin_selection = NULL
    ),
    "operator is invalid for atomic vectors"
  )
})

test_that("subset_samples and MW throws error with character data_bins input", {
  data_to_run_MW <-
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
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 5,
      rand = 1
    ) %>%
    RUtilpol::flatten_list_by_one() %>%
    .[[1]]

  data_source_bins <-
    data_to_run_MW$bins

  data_source_subset <-
    data_to_run_MW$data

  expect_error(
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = "my_bins",
      bin_selection = NULL
    ),
    "operator is invalid for atomic vectors"
  )
})

test_that("subset_samples and MW throws error with list data_bins input", {
  data_to_run_MW <-
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
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 5,
      rand = 1
    ) %>%
    RUtilpol::flatten_list_by_one() %>%
    .[[1]]

  data_source_bins <-
    data_to_run_MW$bins

  data_source_subset <-
    data_to_run_MW$data

  expect_error(
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = list(data_source_bins),
      bin_selection = NULL
    ),
    "non-numeric matrix extent"
  )
})


# ---------------------------------------- #
# Bin_selection input validation
# ---------------------------------------- #

# ---- levels ---- #

test_that(
  "subset_samples and levels returns correct result without input bin_selection",
  {
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

    data_source_subset <-
      data_to_run_levels$data
    data_source_bins <-
      data_to_run_levels$bins

    res_fn <-
      subset_samples(
        data_source_subset = data_source_subset,
        data_source_bins = data_source_bins,
        bin_selection = NULL # bin_selection is ignored for levels
      )
    res_raw <-
      data_source_bins %>%
      dplyr::select("label", "age_diff", "res_age", "start") %>%
      dplyr::inner_join(
        data_source_subset %>%
          tibble::rownames_to_column("start"),
        by = "start"
      ) %>%
      dplyr::mutate(
        age_diff = c(
          diff(
            .data$age
          ), Inf
        ),
        age_diff = ifelse(.data$age_diff == 0, 0.1, .data$age_diff),
        res_age = .data$age
      ) %>%
      dplyr::select(-c("start", "age"))

    expect_identical(
      res_fn,
      res_raw
    )
  }
)

test_that(
  "subset_samples and levels ignores bin_selection parameter",
  {
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

    data_source_subset <-
      data_to_run_levels$data
    data_source_bins <-
      data_to_run_levels$bins

    # Should produce identical results regardless of bin_selection
    res_first <-
      subset_samples(
        data_source_subset = data_source_subset,
        data_source_bins = data_source_bins,
        bin_selection = "first"
      )

    res_random <-
      subset_samples(
        data_source_subset = data_source_subset,
        data_source_bins = data_source_bins,
        bin_selection = "random"
      )

    expect_identical(
      res_first,
      res_random
    )
  }
)

test_that("subset_samples and levels ignores invalid bin_selection", {
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

  data_source_subset <-
    data_to_run_levels$data
  data_source_bins <-
    data_to_run_levels$bins

  res_fn <-
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = "invalid" # bin_selection is ignored for levels
    )
  res_raw <-
    data_source_bins %>%
    dplyr::select("label", "age_diff", "res_age", "start") %>%
    dplyr::inner_join(
      data_source_subset %>%
        tibble::rownames_to_column("start"),
      by = "start"
    ) %>%
    dplyr::mutate(
      age_diff = c(
        diff(
          .data$age
        ), Inf
      ),
      age_diff = ifelse(.data$age_diff == 0, 0.1, .data$age_diff),
      res_age = .data$age
    ) %>%
    dplyr::select(-c("start", "age"))

  expect_identical(
    res_fn,
    res_raw
  )
})

test_that("subset_samples and levels ignores invalid bin_selection", {
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

  data_source_subset <-
    data_to_run_levels$data
  data_source_bins <-
    data_to_run_levels$bins

  res_fn <-
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = 123 # bin_selection is ignored for levels
    )
  res_raw <-
    data_source_bins %>%
    dplyr::select("label", "age_diff", "res_age", "start") %>%
    dplyr::inner_join(
      data_source_subset %>%
        tibble::rownames_to_column("start"),
      by = "start"
    ) %>%
    dplyr::mutate(
      age_diff = c(
        diff(
          .data$age
        ), Inf
      ),
      age_diff = ifelse(.data$age_diff == 0, 0.1, .data$age_diff),
      res_age = .data$age
    ) %>%
    dplyr::select(-c("start", "age"))

  expect_identical(
    res_fn,
    res_raw
  )
})

test_that(
  "subset_samples and levels produces res_age == original age",
  {
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
    data_source_subset <-
      data_to_run_levels$data
    data_source_bins <-
      data_to_run_levels$bins

    res_fn <-
      subset_samples(
        data_source_subset = data_source_subset,
        data_source_bins = data_source_bins
      )

    # res_age should equal the original age from data_source_subset
    expect_equal(
      res_fn$res_age,
      data_source_subset$age
    )
  }
)

test_that(
  "subset_samples and levels preserves all community columns (- age)",
  {
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
    data_source_subset <-
      data_to_run_levels$data
    data_source_bins <-
      data_to_run_levels$bins

    res_fn <-
      subset_samples(
        data_source_subset = data_source_subset,
        data_source_bins = data_source_bins
      )

    # Should have all original community columns plus metadata columns
    community_cols <-
      names(
        data_source_subset
      )[!names(
        data_source_subset
      ) %in% "age"]
    expected_cols <-
      c("label", "age_diff", "res_age", community_cols)

    expect_true(
      all(
        expected_cols %in% names(res_fn)
      )
    )
    expect_equal(
      nrow(res_fn),
      nrow(data_source_bins)
    )
  }
)

# ---- bins ---- #
test_that("subset_samples, bins, NULL bin_selection", {
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

  expect_error(
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = NULL
    ),
    "argument is of length zero"
  )
})


test_that("subset_samples, bins, character bin_selection", {
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

  expect_error(
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = "invalid"
    ),
    # none programmed into function yet (returns all-NA data)
  )
})


test_that("subset_samples, bins, numeric bin_selection", {
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

  expect_error(
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = 123
    ),
    # none programmed into function yet (returns all-NA data)
  )
})


test_that("subset_samples, bins, multiple bin_selection", {
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

  expect_error(
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = c("first", "random")
    ),
    "the condition has length > 1"
  )
})


test_that("subset_samples, bins and no input bin_selection", {
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

  expect_warning(
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins
    ),
    # none programmed into function yet.
    # e.g., "Warning: no bin_selection method supplied. Using default bin_selection = "first"."
  )
})

## --- first bin --- ##
test_that("subset_samples, bins and first bin", {
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

  res_fn <-
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = "first"
    )

  expect_false(
    any(
      is.na(res_fn$label)
    )
  )
  expect_false(
    any(
      is.na(res_fn$res_age)
    )
  )

  # expect no NAs rows (samples for taxa) (?)
  expect_true(
    all(
      rowSums(is.na(res_fn)) == 0
    )
  )
})

## --- random bin --- ##
test_that("subset_samples, bins and random bin", {
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

  res_fn <-
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = "random"
    )

  expect_false(
    any(
      is.na(res_fn$label)
    )
  )
  expect_false(
    any(
      is.na(res_fn$res_age)
    )
  )

  # expect no NAs rows (samples for taxa) (?)
  expect_true(
    all(
      rowSums(is.na(res_fn)) == 0
    )
  )
})

# ---- MW ---- #
test_that("subset_samples, MW, NULL bin_selection", {
  data_to_run_MW <-
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
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 5,
      rand = 1
    ) %>%
    RUtilpol::flatten_list_by_one() %>%
    .[[1]]
  
  data_source_bins <-
    data_to_run_MW$bins
  
  data_source_subset <-
    data_to_run_MW$data
  
  expect_error(
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = NULL
    ),
    "argument is of length zero"
  )
})


test_that("subset_samples, MW, character bin_selection", {
  data_to_run_MW <-
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
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 5,
      rand = 1
    ) %>%
    RUtilpol::flatten_list_by_one() %>%
    .[[1]]
  
  data_source_bins <-
    data_to_run_MW$bins
  
  data_source_subset <-
    data_to_run_MW$data
  
  expect_error(
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = "invalid"
    ),
    # none programmed into function yet (returns all-NA data)
  )
})


test_that("subset_samples, MW, numeric bin_selection", {
  data_to_run_MW <-
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
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 5,
      rand = 1
    ) %>%
    RUtilpol::flatten_list_by_one() %>%
    .[[1]]
  
  data_source_bins <-
    data_to_run_MW$bins
  
  data_source_subset <-
    data_to_run_MW$data
  
  expect_error(
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = 123
    ),
    # none programmed into function yet (returns all-NA data)
  )
})


test_that("subset_samples, MW, multiple bin_selection", {
  data_to_run_MW <-
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
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 5,
      rand = 1
    ) %>%
    RUtilpol::flatten_list_by_one() %>%
    .[[1]]
  
  data_source_bins <-
    data_to_run_MW$bins
  
  data_source_subset <-
    data_to_run_MW$data
  
  expect_error(
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = c("first", "random")
    ),
    "the condition has length > 1"
  )
})


test_that("subset_samples, MW and no input bin_selection", {
  data_to_run_MW <-
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
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 5,
      rand = 1
    ) %>%
    RUtilpol::flatten_list_by_one() %>%
    .[[1]]
  
  data_source_bins <-
    data_to_run_MW$bins
  
  data_source_subset <-
    data_to_run_MW$data
  
  expect_warning(
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins
    ),
    # none programmed into function yet.
    # e.g., "Warning: no bin_selection method supplied. Using default bin_selection = "first"."
  )
})

## --- first bin --- ##
test_that("subset_samples, MW and first bin", {
  data_to_run_MW <-
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
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 5,
      rand = 1
    ) %>%
    RUtilpol::flatten_list_by_one() %>%
    .[[1]]
  
  data_source_bins <-
    data_to_run_MW$bins
  
  data_source_subset <-
    data_to_run_MW$data
  
  res_fn <-
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = "first"
    )
  
  expect_false(
    any(
      is.na(res_fn$label)
    )
  )
  expect_false(
    any(
      is.na(res_fn$res_age)
    )
  )
  
  # expect no NAs rows (samples for taxa) (?)
  expect_true(
    all(
      rowSums(is.na(res_fn)) == 0
    )
  )
})

## --- random bin --- ##
test_that("subset_samples, MW and random bin", {
  data_to_run_MW <-
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
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 5,
      rand = 1
    ) %>%
    RUtilpol::flatten_list_by_one() %>%
    .[[1]]
  
  data_source_bins <-
    data_to_run_MW$bins
  
  data_source_subset <-
    data_to_run_MW$data
  
  res_fn <-
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = "random"
    )
  
  expect_false(
    any(
      is.na(res_fn$label)
    )
  )
  expect_false(
    any(
      is.na(res_fn$res_age)
    )
  )
  
  # expect no NAs rows (samples for taxa) (?)
  expect_true(
    all(
      rowSums(is.na(res_fn)) == 0
    )
  )
})
