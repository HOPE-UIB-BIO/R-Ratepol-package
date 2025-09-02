# Unit tests for run_iteration.
# Parameters:
# data_source_run =
# bin_selection = "first",
# standardise = FALSE,
# n_individuals = 150,
# tranform_to_proportions = TRUE,
# dissimilarity_coefficient = "euc",
# time_standardisation = 500,
# verbose = FALSE

test_that("run_iteration throws error when data_source_run argument is missing", {
    expect_error(
        run_iteration(
            data_source_run = ,
            bin_selection = "first",
            standardise = FALSE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = "euc",
            time_standardisation = 500,
            verbose = FALSE
        ),
        'argument "data_source_run" is missing, with no default'
    )
})

test_that("run_iteration throws error when data_source_run is NULL", {
    expect_error(
        run_iteration(
            data_source_run = NULL,
            bin_selection = "first",
            standardise = FALSE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = "euc",
            time_standardisation = 500,
            verbose = FALSE
        ),
        "argument of length 0"
    )
})

# character
test_that("run_iteration throws error when data_source_run is a character string ('my_data')", {
    expect_error(
        run_iteration(
            data_source_run = "my_data",
            bin_selection = "first",
            standardise = FALSE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = "euc",
            time_standardisation = 500,
            verbose = FALSE
        ),
        "invalid for atomic vectors"
    )
})

# Numeric
test_that("run_iteration throws error when data_source_run is a numeric value (123)", {
    expect_error(
        run_iteration(
            data_source_run = 123,
            bin_selection = "first",
            standardise = FALSE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = "euc",
            time_standardisation = 500,
            verbose = FALSE
        ),
        "invalid for atomic vectors"
    )
})

# Zero
test_that("run_iteration throws error when data_source_run is zero", {
    expect_error(
        run_iteration(
            data_source_run = 0,
            bin_selection = "first",
            standardise = FALSE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = "euc",
            time_standardisation = 500,
            verbose = FALSE
        ),
        "invalid for atomic vectors"
    )
})

# NA
test_that("run_iteration throws error when data_source_run is NA", {
    expect_error(
        run_iteration(
            data_source_run = NA,
            bin_selection = "first",
            standardise = FALSE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = "euc",
            time_standardisation = 500,
            verbose = FALSE
        ),
        "invalid for atomic vectors"
    )
})

# empty list
test_that("run_iteration throws error when data_source_run is an empty list", {
    expect_error(
        run_iteration(
            data_source_run = list(),
            bin_selection = "first",
            standardise = FALSE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = "euc",
            time_standardisation = 500,
            verbose = FALSE
        ),
        "argument of length 0"
    )
})

# empty dataframe
test_that("run_iteration throws error when data_source_run is an empty data frame", {
    expect_error(
        run_iteration(
            data_source_run = data.frame(),
            bin_selection = "first",
            standardise = FALSE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = "euc",
            time_standardisation = 500,
            verbose = FALSE
        ),
        "argument of length 0"
    )
})


# Valid data
test_that("run_iteration throws error when bin_selection is an invalid string ('my_choice')", {
    suppressWarnings(
        data_to_run_bins <-
            extract_data(
                data_community_extract = RRatepol::example_data$pollen_data[[1]],
                data_age_extract = RRatepol::example_data$sample_age[[1]],
                age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
            ) %>%
            smooth_community_data(
                smooth_method = "shep"
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
    )

    expect_error(
        run_iteration(
            data_source_run = data_to_run_bins,
            bin_selection = "my_choice",
            standardise = FALSE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = "euc",
            time_standardisation = 500,
            verbose = FALSE
        ),
        "invalid 'length' argument"
    )
})

test_that("run_iteration throws error when bin_selection is a numeric value (1)", {
    suppressWarnings(
        data_to_run_bins <-
            extract_data(
                data_community_extract = RRatepol::example_data$pollen_data[[1]],
                data_age_extract = RRatepol::example_data$sample_age[[1]],
                age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
            ) %>%
            smooth_community_data(
                smooth_method = "shep"
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
    )

    expect_error(
        run_iteration(
            data_source_run = data_to_run_bins,
            bin_selection = 1,
            standardise = FALSE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = "euc",
            time_standardisation = 500,
            verbose = FALSE
        ),
        "invalid 'length' argument"
    )
})

# multiple bin_selection
test_that("run_iteration throws error when bin_selection contains multiple values ('first', 'random')", {
    suppressWarnings(
        data_to_run_bins <-
            extract_data(
                data_community_extract = RRatepol::example_data$pollen_data[[1]],
                data_age_extract = RRatepol::example_data$sample_age[[1]],
                age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
            ) %>%
            smooth_community_data(
                smooth_method = "shep"
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
    )

    expect_error(
        run_iteration(
            data_source_run = data_to_run_bins,
            bin_selection = c("first", "random"),
            standardise = FALSE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = "euc",
            time_standardisation = 500,
            verbose = FALSE
        ),
        "the condition has length > 1"
    )
})

test_that("run_iteration throws error when bin_selection contains multiple values in different order ('random', 'first')", {
    suppressWarnings(
        data_to_run_bins <-
            extract_data(
                data_community_extract = RRatepol::example_data$pollen_data[[1]],
                data_age_extract = RRatepol::example_data$sample_age[[1]],
                age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
            ) %>%
            smooth_community_data(
                smooth_method = "shep"
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
    )

    expect_error(
        run_iteration(
            data_source_run = data_to_run_bins,
            bin_selection = c("random", "first"),
            standardise = FALSE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = "euc",
            time_standardisation = 500,
            verbose = FALSE
        ),
        "the condition has length > 1"
    )
})

test_that("run_iteration issues warning when bin_selection is missing and uses default value", {
    suppressWarnings(
        data_to_run_bins <-
            extract_data(
                data_community_extract = RRatepol::example_data$pollen_data[[1]],
                data_age_extract = RRatepol::example_data$sample_age[[1]],
                age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
            ) %>%
            smooth_community_data(
                smooth_method = "shep"
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
    )

    expect_warning(
        run_iteration(
            data_source_run = data_to_run_bins,
            bin_selection = ,
            standardise = FALSE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = "euc",
            time_standardisation = 500,
            verbose = FALSE
        ),
        # none programmed into the function yet.
        # e.g., "Warning: No bin_selection supplied.Using default 'first' instead"
    )
})

test_that("run_iteration throws error when bin_selection is NULL", {
    suppressWarnings(
        data_to_run_bins <-
            extract_data(
                data_community_extract = RRatepol::example_data$pollen_data[[1]],
                data_age_extract = RRatepol::example_data$sample_age[[1]],
                age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
            ) %>%
            smooth_community_data(
                smooth_method = "shep"
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
    )

    expect_error(
        run_iteration(
            data_source_run = data_to_run_bins,
            bin_selection = NULL,
            standardise = FALSE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = "euc",
            time_standardisation = 500,
            verbose = FALSE
        ),
        "argument is of length zero"
    )
})



# standardize
test_that("run_iteration handles NULL in standardise parameter without error", {
    suppressWarnings(
        data_to_run_bins <-
            extract_data(
                data_community_extract = RRatepol::example_data$pollen_data[[1]],
                data_age_extract = RRatepol::example_data$sample_age[[1]],
                age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
            ) %>%
            smooth_community_data(
                smooth_method = "shep"
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
    )

    expect_error(
        run_iteration(
            data_source_run = data_to_run_bins,
            bin_selection = "first",
            standardise = NULL,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = "euc",
            time_standardisation = 500,
            verbose = FALSE
        ),
        # No error programmed into the fuction for standardise != TRUE
        # gets ignored - i.e., equivalent to standardise = FALSE
    )
})

# n_individuals:
test_that("run_iteration throws error when n_individuals is zero with standardise=TRUE", {
    suppressWarnings(
        data_to_run_bins <-
            extract_data(
                data_community_extract = RRatepol::example_data$pollen_data[[1]],
                data_age_extract = RRatepol::example_data$sample_age[[1]],
                age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
            ) %>%
            smooth_community_data(
                smooth_method = "shep"
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
    )

    expect_error(
        run_iteration(
            data_source_run = data_to_run_bins,
            bin_selection = "first",
            standardise = TRUE,
            n_individuals = 0,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = "euc",
            time_standardisation = 500,
            verbose = FALSE
        ),
        "invalid 'length' argument"
    )
})

test_that("run_iteration issues warning when n_individuals is too large and uses minimum available count", {
    suppressWarnings(
        data_to_run_bins <-
            extract_data(
                data_community_extract = RRatepol::example_data$pollen_data[[1]],
                data_age_extract = RRatepol::example_data$sample_age[[1]],
                age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
            ) %>%
            smooth_community_data(
                smooth_method = "shep"
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
    )

    expect_warning(
        run_iteration(
            data_source_run = data_to_run_bins,
            bin_selection = "first",
            standardise = TRUE,
            n_individuals = 10000000,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = "euc",
            time_standardisation = 500,
            verbose = TRUE
        ),
        # None programmed into the function yet.
        # If n_individuals is larger than the smallest number of available observations accross samples,
        # it choses min N observations over n_individuals.
        # eg., "Warning: n_individuals was chosen too large.
        # Defaulting to the min number of observations instead."
    )
})

# tranform_to_proportions is only checked if TRUE. Anything else will be regarded as FALSE.
test_that("run_iteration treats numeric zero for tranform_to_proportions as FALSE without error", {
    suppressWarnings(
        data_to_run_bins <-
            extract_data(
                data_community_extract = RRatepol::example_data$pollen_data[[1]],
                data_age_extract = RRatepol::example_data$sample_age[[1]],
                age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
            ) %>%
            smooth_community_data(
                smooth_method = "shep"
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
    )

    expect_error(
        run_iteration(
            data_source_run = data_to_run_bins,
            bin_selection = "first",
            standardise = TRUE,
            n_individuals = 150,
            tranform_to_proportions = 0,
            dissimilarity_coefficient = "euc",
            time_standardisation = 500,
            verbose = FALSE
        )
        # None programmed into the function yet
        # e.g., "invalid argument supplied to 'tranform_to_proportions'"
    )
})

# Dissimilarity_coefficient
# Empty
test_that("run_iteration issues warning when dissimilarity_coefficient is missing and uses default", {
    suppressWarnings(
        data_to_run_bins <-
            extract_data(
                data_community_extract = RRatepol::example_data$pollen_data[[1]],
                data_age_extract = RRatepol::example_data$sample_age[[1]],
                age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
            ) %>%
            smooth_community_data(
                smooth_method = "shep"
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
    )

    expect_warning(
        run_iteration(
            data_source_run = data_to_run_bins,
            bin_selection = "first",
            standardise = TRUE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = ,
            time_standardisation = 500,
            verbose = FALSE
        ),
        # none programmed into the function yet.
        # e.g., "No dissimilarity_coefficient supplied. Defaulting to to 'euc'."
    )
})
# NULL
test_that("run_iteration throws error when dissimilarity_coefficient is NULL", {
    suppressWarnings(
        data_to_run_bins <-
            extract_data(
                data_community_extract = RRatepol::example_data$pollen_data[[1]],
                data_age_extract = RRatepol::example_data$sample_age[[1]],
                age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
            ) %>%
            smooth_community_data(
                smooth_method = "shep"
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
    )

    expect_error(
        run_iteration(
            data_source_run = data_to_run_bins,
            bin_selection = "first",
            standardise = TRUE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = NULL,
            time_standardisation = 500,
            verbose = FALSE
        ),
        "argument is of length zero"
    )
})

# Character
test_that("run_iteration throws error when dissimilarity_coefficient is an invalid string ('my_choice')", {
    suppressWarnings(
        data_to_run_bins <-
            extract_data(
                data_community_extract = RRatepol::example_data$pollen_data[[1]],
                data_age_extract = RRatepol::example_data$sample_age[[1]],
                age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
            ) %>%
            smooth_community_data(
                smooth_method = "shep"
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
    )

    expect_error(
        run_iteration(
            data_source_run = data_to_run_bins,
            bin_selection = "first",
            standardise = TRUE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = "my_choice",
            time_standardisation = 500,
            verbose = FALSE
        ),
        "object 'corrmat' not found"
    )
})

# Numeric
test_that("run_iteration throws error when dissimilarity_coefficient is a numeric value (123)", {
    suppressWarnings(
        data_to_run_bins <-
            extract_data(
                data_community_extract = RRatepol::example_data$pollen_data[[1]],
                data_age_extract = RRatepol::example_data$sample_age[[1]],
                age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
            ) %>%
            smooth_community_data(
                smooth_method = "shep"
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
    )

    expect_error(
        run_iteration(
            data_source_run = data_to_run_bins,
            bin_selection = "first",
            standardise = TRUE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = 123,
            time_standardisation = 500,
            verbose = FALSE
        ),
        "object 'corrmat' not found"
    )
})

# Multiple
test_that("run_iteration throws error when dissimilarity_coefficient contains multiple values ('euc', 'euc.sd')", {
    suppressWarnings(
        data_to_run_bins <-
            extract_data(
                data_community_extract = RRatepol::example_data$pollen_data[[1]],
                data_age_extract = RRatepol::example_data$sample_age[[1]],
                age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
            ) %>%
            smooth_community_data(
                smooth_method = "shep"
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
    )

    expect_error(
        run_iteration(
            data_source_run = data_to_run_bins,
            bin_selection = "first",
            standardise = TRUE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = c("euc", "euc.sd"),
            time_standardisation = 500,
            verbose = FALSE
        ),
        "the condition has length > 1"
    )
})

# Zero
test_that("run_iteration throws error when dissimilarity_coefficient is zero", {
    suppressWarnings(
        data_to_run_bins <-
            extract_data(
                data_community_extract = RRatepol::example_data$pollen_data[[1]],
                data_age_extract = RRatepol::example_data$sample_age[[1]],
                age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
            ) %>%
            smooth_community_data(
                smooth_method = "shep"
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
    )

    expect_error(
        run_iteration(
            data_source_run = data_to_run_bins,
            bin_selection = "first",
            standardise = TRUE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = 0,
            time_standardisation = 500,
            verbose = FALSE
        ),
        "object 'corrmat' not found"
    )
})

# NA
test_that("run_iteration throws error when dissimilarity_coefficient is NA", {
    suppressWarnings(
        data_to_run_bins <-
            extract_data(
                data_community_extract = RRatepol::example_data$pollen_data[[1]],
                data_age_extract = RRatepol::example_data$sample_age[[1]],
                age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
            ) %>%
            smooth_community_data(
                smooth_method = "shep"
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
    )

    expect_error(
        run_iteration(
            data_source_run = data_to_run_bins,
            bin_selection = "first",
            standardise = TRUE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = NA,
            time_standardisation = 500,
            verbose = FALSE
        ),
        "missing value where TRUE/FALSE needed"
    )
})

# time standardiisation:
# empty
test_that("run_iteration issues warning when time_standardisation is missing and uses default value", {
    suppressWarnings(
        data_to_run_bins <-
            extract_data(
                data_community_extract = RRatepol::example_data$pollen_data[[1]],
                data_age_extract = RRatepol::example_data$sample_age[[1]],
                age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
            ) %>%
            smooth_community_data(
                smooth_method = "shep"
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
    )

    expect_warning(
        run_iteration(
            data_source_run = data_to_run_bins,
            bin_selection = "first",
            standardise = TRUE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = "euc",
            time_standardisation = ,
            verbose = FALSE
        ),
        # none programmed into the function yet.
        # e.g., "No time_standardization supplied.
        # Defaulting to time_standardisation = 500 years"
    )
})

# NULL
test_that("run_iteration throws error when time_standardisation is NULL", {
    suppressWarnings(
        data_to_run_bins <-
            extract_data(
                data_community_extract = RRatepol::example_data$pollen_data[[1]],
                data_age_extract = RRatepol::example_data$sample_age[[1]],
                age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
            ) %>%
            smooth_community_data(
                smooth_method = "shep"
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
    )

    expect_error(
        run_iteration(
            data_source_run = data_to_run_bins,
            bin_selection = "first",
            standardise = TRUE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = "euc",
            time_standardisation = NULL,
            verbose = FALSE
        ),
        "`age_diff_st` must be size" # 16 or 1, not 0.
    )
})

# character
test_that("run_iteration throws error when time_standardisation is a character string ('123')", {
    suppressWarnings(
        data_to_run_bins <-
            extract_data(
                data_community_extract = RRatepol::example_data$pollen_data[[1]],
                data_age_extract = RRatepol::example_data$sample_age[[1]],
                age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
            ) %>%
            smooth_community_data(
                smooth_method = "shep"
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
    )

    expect_error(
        run_iteration(
            data_source_run = data_to_run_bins,
            bin_selection = "first",
            standardise = TRUE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = "euc",
            time_standardisation = "123",
            verbose = FALSE
        ),
        "non-numeric argument to binary operator"
        # none programmed into the function yet.
        # e.g., "No time_standardization supplied. Defaulting to "bin_size" standardisation."
    )
})

# 0
test_that("run_iteration throws error when time_standardisation is zero", {
    suppressWarnings(
        data_to_run_bins <-
            extract_data(
                data_community_extract = RRatepol::example_data$pollen_data[[1]],
                data_age_extract = RRatepol::example_data$sample_age[[1]],
                age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
            ) %>%
            smooth_community_data(
                smooth_method = "shep"
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
    )

    expect_error(
        run_iteration(
            data_source_run = data_to_run_bins,
            bin_selection = "first",
            standardise = TRUE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = "euc",
            time_standardisation = 0,
            verbose = FALSE
        ),
        # none programmed into the function yet.
        # e.g., "Error: time_standardisation = 0 results in zero roc."
    )
})

test_that("run_iteration throws error when time_standardisation is NA", {
    suppressWarnings(
        data_to_run_bins <-
            extract_data(
                data_community_extract = RRatepol::example_data$pollen_data[[1]],
                data_age_extract = RRatepol::example_data$sample_age[[1]],
                age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
            ) %>%
            smooth_community_data(
                smooth_method = "shep"
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
    )

    expect_error(
        run_iteration(
            data_source_run = data_to_run_bins,
            bin_selection = "first",
            standardise = TRUE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = "euc",
            time_standardisation = NA,
            verbose = FALSE
        ),
        # none programmed into the function yet.
        # e.g.,
        # "Error: time_standardisation = NA results in NA roc."
    )
})

## Edge cases
# all zero data

test_that("run_iteration handles community data with all zeros when using euc dissimilarity coefficient", {
    community <- RRatepol::example_data$pollen_data[[1]]
    community[,-1 ] <- 0
    age <- RRatepol::example_data$sample_age[[1]]
    age_un <- RRatepol::example_data$age_uncertainty[[1]]
    suppressWarnings(
        data_to_run_bins <-
            extract_data(
                data_community_extract = community,
                data_age_extract = age,
                age_uncertainty = age_un
            ) %>%
            smooth_community_data(
                smooth_method = "shep"
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
    )

    expect_error(
        run_iteration(
            data_source_run = data_to_run_bins,
            bin_selection = "first",
            standardise = TRUE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = "euc",
            time_standardisation = TRUE,
            verbose = FALSE
        ),
        # none programmed into the function yet.
        # e.g.,
        # "Error: time_standardisation = NA results in NA roc."
    )
})

test_that("run_iteration handles community data with all zeros when using euc.sd dissimilarity coefficient", {
    community <- RRatepol::example_data$pollen_data[[1]]
    community[,-1 ] <- 0
    age <- RRatepol::example_data$sample_age[[1]]
    age_un <- RRatepol::example_data$age_uncertainty[[1]]
    suppressWarnings(
        data_to_run_bins <-
            extract_data(
                data_community_extract = community,
                data_age_extract = age,
                age_uncertainty = age_un
            ) %>%
            smooth_community_data(
                smooth_method = "shep"
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
    )

    expect_error(
        run_iteration(
            data_source_run = data_to_run_bins,
            bin_selection = "first",
            standardise = TRUE,
            n_individuals = 150,
            tranform_to_proportions = TRUE,
            dissimilarity_coefficient = "euc.sd",
            time_standardisation = NA,
            verbose = FALSE
        ),
        # none programmed into the function yet.
        # e.g.,
        # "Error: time_standardisation = NA results in NA roc."
    )
})
