# ================================ #
# INPUT VALIDATION TESTS:
# ================================ #


# --------------------------------- #
#   1. data_source_run validation   #
# --------------------------------- #

test_that(
    "run_iteration throws error when data_source_run argument is missing (no default value provided)",
    {
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
    }
)

test_that(
    "run_iteration throws error when data_source_run is NULL (invalid empty input)",
    {
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
    }
)

# character
test_that(
    "run_iteration throws error when data_source_run is a character string ('my_data') instead of required data structure",
    {
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
    }
)

# Numeric
test_that(
    "run_iteration throws error when data_source_run is a numeric value (123) instead of required data structure",
    {
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
    }
)

# Zero
test_that(
    "run_iteration throws error when data_source_run is zero (numeric 0) instead of required data structure",
    {
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
    }
)

# NA
test_that(
    "run_iteration throws error when data_source_run is NA (missing value) instead of required data structure",
    {
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
    }
)

# empty list
test_that(
    "run_iteration throws error when data_source_run is an empty list (list() with no elements)",
    {
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
    }
)

# empty dataframe
test_that(
    "run_iteration throws error when data_source_run is an empty data frame (data.frame() with no rows/columns)",
    {
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
    }
)


# --------------------------------- #
#     2. bin_selection validation   #
# --------------------------------- #


# bin_selection = string
test_that(
    "run_iteration throws error when bin_selection is an invalid string ('my_choice') not among allowed options ('first', 'random', 'last')",
    {
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
    }
)

# bin_selection numeric
test_that(
    "run_iteration throws error when bin_selection is a numeric value (1) instead of required character string",
    {
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
    }
)

# multiple bin_selection
test_that(
    "run_iteration throws error when bin_selection contains multiple values ('first', 'random') instead of a single option",
    {
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
    }
)

# multiple bin selection in reversed order
test_that(
    "run_iteration throws error when bin_selection contains multiple values in different order ('random', 'first') instead of a single option",
    {
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
    }
)

# bin_selection missing
test_that(
    "run_iteration issues warning when bin_selection is missing and uses default value ('first')",
    {
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
    }
)

# bin_selection = NULL
test_that(
    "run_iteration throws error when bin_selection is NULL instead of a valid selection option",
    {
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
    }
)


# --------------------------------- #
#    3. 'standardise' validation    #
# --------------------------------- #

# invalid standardise parameter:
# (only TRUE is checked, anything else is FALSE)

# standardise = NULL
test_that(
    "run_iteration throws error with NULL in standardise parameter",
    {
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
    }
)

# output validation for (invalid) standardise parameter (1)
test_that(
    "run_iteration treats non-boolean value for standardise as FALSE",
    {
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

        # This should work like standardise=FALSE since isTRUE(1) is FALSE
        result <-
            run_iteration(
                data_source_run = data_to_run_bins,
                bin_selection = "first",
                standardise = 1,
                n_individuals = 150,
                tranform_to_proportions = TRUE,
                dissimilarity_coefficient = "euc",
                time_standardisation = 500,
                verbose = FALSE
            )

        expect_s3_class(result, "data.frame")
    }
)


# --------------------------------- #
#  4. 'n_individuals' validation    #
# --------------------------------- #

# valid standardise but invalid n_individuals:
test_that(
    "run_iteration throws error when n_individuals is zero with standardise=TRUE (cannot standardize with zero individuals)",
    {
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
    }
)

# Test for combinations of parameters
test_that(
    "run_iteration throws warning for invalid combinations of standardise=TRUE with n_individuals=NULL",
    {
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
                n_individuals = NULL,
                tranform_to_proportions = TRUE,
                dissimilarity_coefficient = "euc",
                time_standardisation = 500,
                verbose = FALSE
            ),
            # this scenario silently uses the minimum count available
            # instead of user supplied n_individuals
        )
    }
)


# --------------------------------- #
#    5.tranform_to_proportions      #
# --------------------------------- #


# tranform_to_proportions is only checked if TRUE. Anything else will be regarded as FALSE.
test_that(
    "run_iteration throws error for non-boolian tranform_to_proportions = 0",
    {
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
            # None programmed into the function yet since it only recognizes TRUE
            # gets ignored - i.e., equivalent to tranform_to_proportions = FALSE
            # e.g., "invalid argument supplied to 'tranform_to_proportions'"
        )
    }
)


# --------------------------------- #
#  6. 'dissimilarity_coefficient'   #
# --------------------------------- #


# Dissimilarity_coefficient
# Empty
test_that(
    "run_iteration issues warning when dissimilarity_coefficient is missing and uses default ('euc')",
    {
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
    }
)
# NULL
test_that(
    "run_iteration throws error when dissimilarity_coefficient is NULL instead of valid coefficient option",
    {
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
    }
)

# Character
test_that(
    "run_iteration throws error when dissimilarity_coefficient is an invalid string ('my_choice') not among allowed options",
    {
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
    }
)

# Numeric
test_that(
    "run_iteration throws error when dissimilarity_coefficient is a numeric value (123) instead of valid character string",
    {
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
    }
)

# Multiple
test_that(
    "run_iteration throws error when dissimilarity_coefficient contains multiple values ('euc', 'euc.sd') instead of single option",
    {
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
    }
)

# Zero
test_that(
    "run_iteration throws error when dissimilarity_coefficient is zero (numeric 0) instead of valid character string",
    {
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
    }
)

# NA
test_that(
    "run_iteration throws error when dissimilarity_coefficient is NA (missing value) instead of valid coefficient option",
    {
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
    }
)


# --------------------------------- #
#   7. 'time_standardisation'       #
# --------------------------------- #

# time standardiisation:
# empty
test_that(
    "run_iteration issues warning when time_standardisation is missing and uses default value",
    {
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
    }
)

# NULL
test_that(
    "run_iteration throws error when time_standardisation is NULL",
    {
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
    }
)

# character
test_that(
    "run_iteration throws error when time_standardisation is a character string ('123')",
    {
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
    }
)

# 0
test_that(
    "run_iteration throws error when time_standardisation is zero",
    {
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
    }
)

test_that(
    "run_iteration throws error when time_standardisation is NA",
    {
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
    }
)


# --------------------------------- #
#   8. data edge cases              #
# --------------------------------- #
# all zero data

test_that(
    "run_iteration handles community data with all zeros when using euc dissimilarity coefficient",
    {
        community <-
            RRatepol::example_data$pollen_data[[1]]
        community[, -1] <-
            0
        age <-
            RRatepol::example_data$sample_age[[1]]
        age_un <-
            RRatepol::example_data$age_uncertainty[[1]]

        expect_error(
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
                    .[[1]] %>%
                    run_iteration(
                        data_source_run = .,
                        bin_selection = "first",
                        standardise = TRUE,
                        n_individuals = 150,
                        tranform_to_proportions = TRUE,
                        dissimilarity_coefficient = "euc",
                        time_standardisation = TRUE,
                        verbose = FALSE
                    )
            ),
            "subscript out of bounds"
        )
    }
)

test_that(
    "run_iteration handles community data with all zeros when using euc.sd dissimilarity coefficient",
    {
        community <-
            RRatepol::example_data$pollen_data[[1]]
        community[, -1] <-
            0
        age <-
            RRatepol::example_data$sample_age[[1]]
        age_un <-
            RRatepol::example_data$age_uncertainty[[1]]
        expect_error(
            suppressWarnings(
                data_to_run_bins <-
                    extract_data(
                        data_community_extract = community,
                        data_age_extract = age,
                        age_uncertainty = age_un
                    ) %>%
                    # it fails during smoothing, here!
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
                    .[[1]] %>%
                    run_iteration(
                        data_source_run = .,
                        bin_selection = "first",
                        standardise = TRUE,
                        n_individuals = 150,
                        tranform_to_proportions = TRUE,
                        dissimilarity_coefficient = "euc.sd",
                        time_standardisation = TRUE,
                        verbose = FALSE
                    )
            ),
            "subscript out of bounds"
        )
    }
)



# Similarly for tranform_to_proportions
test_that(
    "run_iteration treats non-boolean value for tranform_to_proportions as FALSE",
    {
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

        # This should work like tranform_to_proportions=FALSE since isTRUE(1) is FALSE
        result <-
            run_iteration(
                data_source_run = data_to_run_bins,
                bin_selection = "first",
                standardise = FALSE,
                n_individuals = 150,
                tranform_to_proportions = 1,
                dissimilarity_coefficient = "euc",
                time_standardisation = 500,
                verbose = FALSE
            )

        expect_s3_class(result, "data.frame")
    }
)

# Test for the combination of time_standardisation with TRUE
test_that(
    "run_iteration handles boolean TRUE for time_standardisation as 1",
    {
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

        result <-
            run_iteration(
                data_source_run = data_to_run_bins,
                bin_selection = "first",
                standardise = FALSE,
                n_individuals = 150,
                tranform_to_proportions = TRUE,
                dissimilarity_coefficient = "euc",
                time_standardisation = TRUE, # TRUE is coerced to 1
                verbose = FALSE
            )

        expect_s3_class(result, "data.frame")
    }
)

# Test for warning when n_individuals is large
## Not sure what's going on here. sometimes the test passes, sometimes it fails.
test_that(
    "run_iteration throws error if standardisation failed if verbose = TRUE",
    {
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

        # this works:
        set.seed(123)
        expect_condition(
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
            "Data standardisation was unsuccesfull, try 'standardise' = FALSE"
        )
    }
)
# to do: what is the expectation?
# with verbose = FALSE
test_that(
  "run_iteration throws no message if standardisation failed if verbose = FALSE",
  {
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


        set.seed(123)
        expect_no_condition(
            res <-
                run_iteration(
                    data_source_run = data_to_run_bins,
                    bin_selection = "first",
                    standardise = TRUE,
                    n_individuals = 10000000,
                    tranform_to_proportions = TRUE,
                    dissimilarity_coefficient = "euc",
                    time_standardisation = 500,
                    verbose = FALSE
                )
        )

        
    }
)

# tranform_to_proportions
# tranform_to_proportions is only checked if TRUE. Anything else will be regarded as FALSE.
test_that(
    "run_iteration treats numeric zero for tranform_to_proportions as FALSE without error",
    {
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
    }
)

# Dissimilarity_coefficient
# Empty
test_that(
    "run_iteration issues warning when dissimilarity_coefficient is missing and uses default",
    {
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
    }
)

# NULL
test_that(
    "run_iteration throws error when dissimilarity_coefficient is NULL",
    {
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
    }
)

# Character
test_that(
    "run_iteration throws error when dissimilarity_coefficient is an invalid string ('my_choice')",
    {
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
    }
)

# Numeric
test_that(
    "run_iteration throws error when dissimilarity_coefficient is a numeric value (123)",
    {
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
    }
)

# Multiple
test_that(
    "run_iteration throws error when dissimilarity_coefficient contains multiple values ('euc', 'euc.sd')",
    {
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
    }
)

# Zero
test_that(
    "run_iteration throws error when dissimilarity_coefficient is zero",
    {
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
    }
)

# NA
test_that(
    "run_iteration throws error when dissimilarity_coefficient is NA",
    {
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
    }
)

# 9. 'verbose' validation


test_that(
    "run_iteration issues a warning if NULL is supplied for verbose",
    {
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
                time_standardisation = 500,
                verbose = NULL
            )
            # a message warning that verbose = NULL will be treated as verbose = FALSE
        )
    }
)

test_that(
    "run_iteration issues a warning if NA is supplied for verbose",
    {
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
                time_standardisation = 500,
                verbose = NA
            )
        )
        # a message warning that verbose = NA will be treated as verbose = FALSE
    }
)


# ================================ #
# OUTPUT VALIDATION TESTS:
# ================================ #


test_that(
    "run_iteration output has expected structure and column names",
    {
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

        set.seed(123)
        result <-
            run_iteration(
                data_source_run = data_to_run_bins,
                bin_selection = "first",
                standardise = FALSE,
                n_individuals = 150,
                tranform_to_proportions = TRUE,
                dissimilarity_coefficient = "euc",
                time_standardisation = 500,
                verbose = FALSE
            )

        # Check output structure
        expect_s3_class(result, "data.frame")
        expect_true(all(c("label", "res_age", "roc") %in% names(result)))
        expect_true(
            nrow(
                result
            ) > 0
        ) # Should have at least one row
    }
)

test_that(
    "run_iteration output has non-NA values in expected columns",
    {
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

        set.seed(123)
        result <-
            run_iteration(
                data_source_run = data_to_run_bins,
                bin_selection = "first",
                standardise = FALSE,
                n_individuals = 150,
                tranform_to_proportions = TRUE,
                dissimilarity_coefficient = "euc",
                time_standardisation = 500,
                verbose = FALSE
            )

        # Check for NA values
        expect_false(any(is.na(result$label)))
        expect_false(any(is.na(result$res_age)))
        expect_false(any(is.na(result$roc)))
    }
)

test_that(
    "run_iteration RoC values are numeric and positive",
    {
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

        set.seed(123)
        result <-
            run_iteration(
                data_source_run = data_to_run_bins,
                bin_selection = "first",
                standardise = FALSE,
                n_individuals = 150,
                tranform_to_proportions = TRUE,
                dissimilarity_coefficient = "euc",
                time_standardisation = 500,
                verbose = FALSE
            )

        # Check RoC values
        expect_type(result$roc, "double")
        expect_true(
            all(
                result$roc >= 0
            )
        ) # RoC should be non-negative
        expect_true(
            all(
                is.finite(
                    result$roc
                )
            )
        ) # No Inf or -Inf values
    }
)

test_that(
    "run_iteration with different time_standardisation gives proportional RoC values",
    {
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

        set.seed(123)
        result1 <-
            run_iteration(
                data_source_run = data_to_run_bins,
                bin_selection = "first",
                standardise = FALSE,
                n_individuals = 150,
                tranform_to_proportions = TRUE,
                dissimilarity_coefficient = "euc",
                time_standardisation = 500,
                verbose = FALSE
            )

        set.seed(123)
        result2 <-
            run_iteration(
                data_source_run = data_to_run_bins,
                bin_selection = "first",
                standardise = FALSE,
                n_individuals = 150,
                tranform_to_proportions = TRUE,
                dissimilarity_coefficient = "euc",
                time_standardisation = 1000, # Double the time standardization
                verbose = FALSE
            )

        # With double time_standardisation, RoC should be doubled
        expect_equal(result1$roc * 2, result2$roc, tolerance = 1e-6)
    }
)


test_that(
    "run_iteration with standardise=TRUE produces valid output",
    {
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

        set.seed(123)
        result <-
            run_iteration(
                data_source_run = data_to_run_bins,
                bin_selection = "first",
                standardise = TRUE,
                n_individuals = 150,
                tranform_to_proportions = TRUE,
                dissimilarity_coefficient = "euc",
                time_standardisation = 500,
                verbose = FALSE
            )

        # Check output with standardization
        expect_s3_class(result, "data.frame")
        expect_true(all(c("label", "res_age", "roc") %in% names(result)))
        expect_false(any(is.na(result$roc)))
        expect_true(all(result$roc >= 0))
    }
)

test_that(
    "run_iteration with different dissimilarity_coefficient produces different RoC values",
    {
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

        set.seed(123)
        result_euc <-
            run_iteration(
                data_source_run = data_to_run_bins,
                bin_selection = "first",
                standardise = FALSE,
                n_individuals = 150,
                tranform_to_proportions = TRUE,
                dissimilarity_coefficient = "euc",
                time_standardisation = 500,
                verbose = FALSE
            )

        set.seed(123)
        result_euc_sd <-
            run_iteration(
                data_source_run = data_to_run_bins,
                bin_selection = "first",
                standardise = FALSE,
                n_individuals = 150,
                tranform_to_proportions = TRUE,
                dissimilarity_coefficient = "euc.sd",
                time_standardisation = 500,
                verbose = FALSE
            )

        # Different dissimilarity coefficients should produce different results
        expect_false(identical(result_euc$roc, result_euc_sd$roc))
    }
)

test_that(
    "run_iteration output is deterministic with fixed seed for random bin_selection",
    {
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

        set.seed(123)
        result1 <-
            run_iteration(
                data_source_run = data_to_run_bins,
                bin_selection = "random",
                standardise = FALSE,
                n_individuals = 150,
                tranform_to_proportions = TRUE,
                dissimilarity_coefficient = "euc",
                time_standardisation = 500,
                verbose = FALSE
            )

        set.seed(123)
        result2 <-
            run_iteration(
                data_source_run = data_to_run_bins,
                bin_selection = "random",
                standardise = FALSE,
                n_individuals = 150,
                tranform_to_proportions = TRUE,
                dissimilarity_coefficient = "euc",
                time_standardisation = 500,
                verbose = FALSE
            )

        # Same seed should produce identical results with random bin selection
        expect_identical(result1, result2)
    }
)

test_that(
    "run_iteration output has valid label column",
    {
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

        set.seed(123)
        result <-
            run_iteration(
                data_source_run = data_to_run_bins,
                bin_selection = "first",
                standardise = FALSE,
                n_individuals = 150,
                tranform_to_proportions = TRUE,
                dissimilarity_coefficient = "euc",
                time_standardisation = 500,
                verbose = FALSE
            )

        # Check that labels are valid
        expect_type(result$label, "character")
        expect_true(all(nchar(result$label) > 0))
        expect_equal(
            length(
                unique(
                    result$label
                )
            ), nrow(
                result
            )
        ) # All labels should be unique
    }
)

test_that(
    "run_iteration handles very small time_standardisation correctly",
    {
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

        # Use a very small time_standardisation to test edge behavior
        result <-
            run_iteration(
                data_source_run = data_to_run_bins,
                bin_selection = "first",
                standardise = FALSE,
                n_individuals = 150,
                tranform_to_proportions = TRUE,
                dissimilarity_coefficient = "euc",
                time_standardisation = 0.001, # Very small value to test division behavior
                verbose = FALSE
            )

        # Check that results are still valid (not Inf)
        expect_true(all(is.finite(result$roc)))
    }
)

test_that(
    "run_iteration output age values match input data range",
    {
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

        set.seed(123)
        result <-
            run_iteration(
                data_source_run = data_to_run_bins,
                bin_selection = "first",
                standardise = FALSE,
                n_individuals = 150,
                tranform_to_proportions = TRUE,
                dissimilarity_coefficient = "euc",
                time_standardisation = 500,
                verbose = FALSE
            )

        # Check age range consistency
        min_input_age <-
            min(data_to_run_bins$data$age)
        max_input_age <-
            max(data_to_run_bins$data$age)

        # Result ages should be within or close to input data range
        expect_true(min(result$res_age) >= min_input_age || abs(min(result$res_age) - min_input_age) < 500)
        expect_true(max(result$res_age) <= max_input_age || abs(max(result$res_age) - max_input_age) < 500)
    }
)

test_that(
    "run_iteration output has reasonable RoC values for typical inputs",
    {
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

        set.seed(123)
        result <-
            run_iteration(
                data_source_run = data_to_run_bins,
                bin_selection = "first",
                standardise = FALSE,
                n_individuals = 150,
                tranform_to_proportions = TRUE,
                dissimilarity_coefficient = "euc",
                time_standardisation = 500,
                verbose = FALSE
            )

        # RoC values should be in a reasonable range for ecological data
        # Most common range is between 0 and 10 for standardized data
        mean_roc <-
            mean(result$roc)
        expect_true(mean_roc > 0 && mean_roc < 100)
    }
)

# Test for successful run with all valid parameters
test_that(
    "run_iteration executes successfully with all valid parameters",
    {
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

        result <-
            run_iteration(
                data_source_run = data_to_run_bins,
                bin_selection = "first",
                standardise = FALSE,
                n_individuals = 150,
                tranform_to_proportions = TRUE,
                dissimilarity_coefficient = "euc",
                time_standardisation = 500,
                verbose = FALSE
            )

        expect_s3_class(result, "data.frame")
        expect_true("roc" %in% names(result))
    }
)

