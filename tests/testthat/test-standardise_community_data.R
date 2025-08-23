# False input validation
test_that(
  "standardise_community_data throws error if wrong input data is supploed", {
  expect_error(
    standardise_community_data(
      n_individuals = 150
    ),
    '"data_source_standard" is missing, with no default'
  )
})

test_that(
  "standardise_community_data throws error if wrong input data is supploed", {
  expect_error(
    standardise_community_data(
      data_source_standard = NULL,
      n_individuals = 150
    ),
    "no applicable method for 'select' applied"
  )
})

test_that(
  "standardise_community_data throws error if wrong input data is supploed", {
  expect_error(
    standardise_community_data(
      data_source_standard = 123,
      n_individuals = 150
    ),
    "no applicable method for 'select' applied"
  )
})

test_that(
  "standardise_community_data throws error if wrong input data is supploed", {
  expect_error(
    standardise_community_data(
      data_source_standard = "my_data",
      n_individuals = 150
    ),
    "no applicable method for 'select' applied"
  )
})

test_that(
  "standardise_community_data throws error if wrong input data is supploed", {
  expect_error(
    standardise_community_data(
      data_source_standard = NA,
      n_individuals = 150
    ),
    "no applicable method for 'select' applied"
  )
})

test_that(
  "standardise_community_data throws error if wrong input data is supploed", {
  expect_error(
    standardise_community_data(
      data_source_standard = list(),
      n_individuals = 150
    ),
    "doesn't handle lists."
  )
})

test_that(
  "standardise_community_data throws error if wrong input data is supploed", {
  expect_error(
    standardise_community_data(
      data_source_standard = data.frame(),
      n_individuals = 150
    ),
    "undefined columns selected"
  )
})

test_that(
  "standardise_community_data throws error if wrong input data is supploed", {
  expect_error(
    standardise_community_data(
      data_source_standard = matrix(),
      n_individuals = 150
    ),
    "no applicable method for 'select' applied"
  )
})

# workflow within run_iteration:
## Wrong n_individuals input
### NULL
test_that(
  "standardise_community_data throws warning with n_individuals = NULL", {
  n_individuals <-
    NULL
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

  data_subset <-
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = "first"
    ) %>%
    reduce_data_simple()

  com_data_sums <-
    rowSums(
      subset_community(
        data_source = data_subset
      ),
      na.rm = TRUE
    )

  # adjust the value to a minimal of presented values
  n_individuals <-
    min(
      c(
        com_data_sums,
        n_individuals
      )
    )

  # check if all samples has n_individuals of individuals
  data_subset <-
    data_subset[com_data_sums >= n_individuals, ]

  data_subset <-
    reduce_data_simple(
      data_source_reduce = data_subset
    )

  # standardisation
  set.seed(123)
  expect_warning(
    data_sd <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ),
    # none programmed into function yet
    # e.g.,
    # "Warning: n_individuals = NULL will use min number of observations from the data for standardisation."
  )
})

### Character
test_that(
  "standardise_community_data returns identical results with character or numeric n_individuals", {
  n_individuals_char <-
    "10"

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

  data_subset <-
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = "first"
    ) %>%
    reduce_data_simple()

  com_data_sums <-
    rowSums(
      subset_community(
        data_source = data_subset
      ),
      na.rm = TRUE
    )

  # adjust the value to a minimal of presented values
  n_individuals_char <-
    min(
      c(
        com_data_sums,
        n_individuals_char
      )
    )

  # check if all samples has n_individuals of individuals:
  data_subset_char <-
    data_subset[com_data_sums >= n_individuals_char, ]

  data_subset_char <-
    reduce_data_simple(
      data_source_reduce = data_subset_char
    )

  # standardisation
  set.seed(123)
  expect_no_error(
    data_sd_char <-
      standardise_community_data(
        data_source_standard = data_subset_char,
        n_individuals = n_individuals_char
      )
    #
  )

  expect_true(
    all(
      rowSums(
    data_sd_char[, -c(
    1:3)]) == 10
    )
  )

  # Control:
  n_individuals_num <-
    10
  
  # adjust the value to a minimal of presented values
  n_individuals_num <-
    min(
      c(
        com_data_sums,
        n_individuals_num
      )
    )

  # check if all samples has n_individuals of individuals
  data_subset_num <-
    data_subset[com_data_sums >= n_individuals_num, ]

  data_subset_num <-
    reduce_data_simple(
      data_source_reduce = data_subset_num
    )

  expect_identical(
    data_subset_char,
    data_subset_num
  )

  # standardisation
  set.seed(123)
  data_sd_num <-
    standardise_community_data(
      data_source_standard = data_subset_num,
      n_individuals = n_individuals_num
    )

  expect_identical(
    data_sd_char,
    data_sd_num
  )
})


### high n_individuals
test_that(
  "standardise_community_data returns n_individuals observations in samples with high n_individuals (within run_iteration workflow)", {
  n_individuals <-
    100000
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

  data_subset <-
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = "first"
    ) %>%
    reduce_data_simple()

  com_data_sums <-
    rowSums(
      subset_community(
        data_source = data_subset
      ),
      na.rm = TRUE
    )

  # adjust the value to a minimal of presented values
  n_individuals <-
    min(
      c(
        com_data_sums,
        n_individuals
      )
    )

  # check if all samples has n_individuals of individuals
  data_subset <-
    data_subset[com_data_sums >= n_individuals, ]

  data_subset <-
    reduce_data_simple(
      data_source_reduce = data_subset
    )

  # standardisation
  set.seed(123)
  data_sd <-
    standardise_community_data(
      data_source_standard = data_subset,
      n_individuals = n_individuals
    )

  expect_false(
    identical(
      data_subset,
      data_sd
    )
  )

  expect_true(
    all(rowSums(data_sd[, -c(1:3)]) == n_individuals)
  )
})

### high
test_that(
  "standardise_community_data returns standardised data with high n_individuals (within run_iteration workflow)", {
  n_individuals <-
    100000
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

  data_subset <-
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = "first"
    ) %>%
    reduce_data_simple()

  com_data_sums <-
    rowSums(
      subset_community(
        data_source = data_subset
      ),
      na.rm = TRUE
    )

  # adjust the value to a minimal of presented values
  n_individuals <-
    min(
      c(
        com_data_sums,
        n_individuals
      )
    )

  # check if all samples has n_individuals of individuals
  data_subset <-
    data_subset[com_data_sums >= n_individuals, ]

  data_subset <-
    reduce_data_simple(
      data_source_reduce = data_subset
    )

  # standardisation
  set.seed(123)
  data_sd <-
    standardise_community_data(
      data_source_standard = data_subset,
      n_individuals = n_individuals
    )

  expect_false(
    identical(
      data_subset,
      data_sd
    )
  )
})

### negative 
test_that(
  "standardise_community_data fails with negative n_individuals", {
    n_individuals <-
      -150
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
    
    data_subset <-
      subset_samples(
        data_source_subset = data_source_subset,
        data_source_bins = data_source_bins,
        bin_selection = "first"
      ) %>%
      reduce_data_simple()
    
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    
    # adjust the value to a minimal of presented values
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    
    # check if all samples has n_individuals of individuals
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    
    # standardisation
    set.seed(123)
    expect_error(
      data_sd <-
        standardise_community_data(
          data_source_standard = data_subset,
          n_individuals = n_individuals
        ),
      "invalid 'size' argument"
    )
  })

### n_individuals = 1
test_that(
  "standardise_community_data works with n_individuals = 1", {
    n_individuals <-
      1
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
    
    data_subset <-
      subset_samples(
        data_source_subset = data_source_subset,
        data_source_bins = data_source_bins,
        bin_selection = "first"
      ) %>%
      reduce_data_simple()
    
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    
    # adjust the value to a minimal of presented values
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    
    # check if all samples has n_individuals of individuals
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    
    # standardisation
    set.seed(123)
    data_sd <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      )
    
    expect_true(
      all(
        rowSums(
          data_sd[, -c(
            1:3)]) == n_individuals
      )
    )
  })

# n_individuals = 0
test_that(
  "standardise_community_data fails with zero n_individuals", {
    n_individuals <-
      0
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
    
    data_subset <-
      subset_samples(
        data_source_subset = data_source_subset,
        data_source_bins = data_source_bins,
        bin_selection = "first"
      ) %>%
      reduce_data_simple()
    
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    
    # adjust the value to a minimal of presented values
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    
    # check if all samples has n_individuals of individuals
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    
    # standardisation
    set.seed(123)
    expect_error(
      data_sd <-
        standardise_community_data(
          data_source_standard = data_subset,
          n_individuals = n_individuals
        ),
      # none programmed into the function yet.
      # e.g.,
      # "Error: n_individuals must be > 0"
    )
  })

# without run_iteration workflow:
### high
test_that(
  "standardise_community_data fails if n_individuals >> than n observations in sample (outside of run_iteration workflow)", {
  n_individuals <-
    100000
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

  data_subset <-
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = "first"
    ) %>%
    reduce_data_simple()

  # standardisation
  set.seed(123)

  expect_error(
    standardise_community_data(
      data_source_standard = data_subset,
      n_individuals = n_individuals
    ),
    "cannot take a sample larger than the population when 'replace = FALSE'"
  )
})

# Output validation
# Valid data:
test_that(
  "standardise_community_data functions correctly (rowSums = n_individuals) with valid/default data", {
  n_individuals <-
    150
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

  data_subset <-
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = "first"
    ) %>%
    reduce_data_simple()

  com_data_sums <-
    rowSums(
      subset_community(
        data_source = data_subset
      ),
      na.rm = TRUE
    )

  # adjust the value to a minimal of presented values
  n_individuals <-
    min(
      c(
        com_data_sums,
        n_individuals
      )
    )

  # check if all samples has n_individuals of individuals
  data_subset <-
    data_subset[com_data_sums >= n_individuals, ]

  data_subset <-
    reduce_data_simple(
      data_source_reduce = data_subset
    )

  # standardisation
  set.seed(123)
  data_sd <-
    standardise_community_data(
      data_source_standard = data_subset,
      n_individuals = n_individuals
    )

  expect_false(
    identical(
      data_subset,
      data_sd
    )
  )
  expect_true(
    all(
      rowSums(
    data_sd[, -c(
    1:3)]) == n_individuals
    )
  )
})

test_that(
  "standardise_community_data returns data.frame with valid data", {
  n_individuals <-
    150
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

  data_subset <-
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = "first"
    ) %>%
    reduce_data_simple()

  com_data_sums <-
    rowSums(
      subset_community(
        data_source = data_subset
      ),
      na.rm = TRUE
    )

  # adjust the value to a minimal of presented values
  n_individuals <-
    min(
      c(
        com_data_sums,
        n_individuals
      )
    )

  # check if all samples has n_individuals of individuals
  data_subset <-
    data_subset[com_data_sums >= n_individuals, ]

  data_subset <-
    reduce_data_simple(
      data_source_reduce = data_subset
    )

  # standardisation
  set.seed(123)
  data_sd <-
    standardise_community_data(
      data_source_standard = data_subset,
      n_individuals = n_individuals
    )

  expect_s3_class(
    data_sd, "data.frame"
  )
})

test_that(
  "standardise_community_data returns colnames correctly with valid data", {
  n_individuals <-
    150
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

  data_subset <-
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = "first"
    ) %>%
    reduce_data_simple()

  com_data_sums <-
    rowSums(
      subset_community(
        data_source = data_subset
      ),
      na.rm = TRUE
    )

  # adjust the value to a minimal of presented values
  n_individuals <-
    min(
      c(
        com_data_sums,
        n_individuals
      )
    )

  # check if all samples has n_individuals of individuals
  data_subset <-
    data_subset[com_data_sums >= n_individuals, ]

  data_subset <-
    reduce_data_simple(
      data_source_reduce = data_subset
    )

  # standardisation
  set.seed(123)
  data_sd <-
    standardise_community_data(
      data_source_standard = data_subset,
      n_individuals = n_individuals
    )

  expect_true(
    all(
      c("label", "res_age", "age_diff") %in% names(data_sd)
    )
  )
})

# Reproducibility tests
test_that(
  "standardise_community_data returns consistent results with set seed", {
  n_individuals <-
    150
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

  data_subset <-
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = "first"
    ) %>%
    reduce_data_simple()

  com_data_sums <-
    rowSums(
      subset_community(
        data_source = data_subset
      ),
      na.rm = TRUE
    )

  # adjust the value to a minimal of presented values
  n_individuals <-
    min(
      c(
        com_data_sums,
        n_individuals
      )
    )

  # check if all samples has n_individuals of individuals
  data_subset <-
    data_subset[com_data_sums >= n_individuals, ]

  data_subset <-
    reduce_data_simple(
      data_source_reduce = data_subset
    )

  set.seed(123)
  data_sd_1 <-
    standardise_community_data(
      data_source_standard = data_subset,
      n_individuals = n_individuals
    )

  set.seed(123)
  data_sd_2 <-
    standardise_community_data(
      data_source_standard = data_subset,
      n_individuals = n_individuals
    )

  expect_identical(
    data_sd_1,
    data_sd_2
  )
})

# Different seed produces different results
test_that(
  "standardise_community_data produces different results with different seeds", {
  n_individuals <-
    150
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

  data_subset <-
    subset_samples(
      data_source_subset = data_source_subset,
      data_source_bins = data_source_bins,
      bin_selection = "first"
    ) %>%
    reduce_data_simple()

  com_data_sums <-
    rowSums(
      subset_community(
        data_source = data_subset
      ),
      na.rm = TRUE
    )

  # adjust the value to a minimal of presented values
  n_individuals <-
    min(
      c(
        com_data_sums,
        n_individuals
      )
    )

  # check if all samples has n_individuals of individuals
  data_subset <-
    data_subset[com_data_sums >= n_individuals, ]

  data_subset <-
    reduce_data_simple(
      data_source_reduce = data_subset
    )
  set.seed(123)
  result1 <-
    standardise_community_data(
    data_subset,
    n_individuals
  )

  set.seed(456)
  result2 <-
    standardise_community_data(
    data_subset,
    n_individuals
  )

  expect_false(identical(result1, result2))
})
