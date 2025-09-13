# Input validation (Errors)
## 1. Empty argument tests
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc argument is missing for euc coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = ,
        dissimilarity_coefficient = "euc"
      ),
      'argument "data_source_dc" is missing, with no default'
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc argument is missing for euc.sd coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = ,
        dissimilarity_coefficient = "euc.sd"
      ),
      'argument "data_source_dc" is missing, with no default'
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc argument is missing for chord coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = ,
        dissimilarity_coefficient = "chord"
      ),
      'argument "data_source_dc" is missing, with no default'
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc argument is missing for chisq coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = ,
        dissimilarity_coefficient = "chisq"
      ),
      'argument "data_source_dc" is missing, with no default'
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc argument is missing for gower coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = ,
        dissimilarity_coefficient = "gower"
      ),
      'argument "data_source_dc" is missing, with no default'
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc argument is missing for bray coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = ,
        dissimilarity_coefficient = "bray"
      ),
      'argument "data_source_dc" is missing, with no default'
    )
  }
)

## 2. NULL argument tests
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is NULL for euc coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = NULL,
        dissimilarity_coefficient = "euc"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is NULL for euc.sd coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = NULL,
        dissimilarity_coefficient = "euc.sd"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is NULL for chord coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = NULL,
        dissimilarity_coefficient = "chord"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is NULL for chisq coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = NULL,
        dissimilarity_coefficient = "chisq"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is NULL for gower coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = NULL,
        dissimilarity_coefficient = "gower"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is NULL for bray coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = NULL,
        dissimilarity_coefficient = "bray"
      ),
      "invalid 'length' argument"
    )
  }
)
# Character string input tests
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is character string for euc coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = "my_data",
        dissimilarity_coefficient = "euc"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is character string for euc.sd coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = "my_data",
        dissimilarity_coefficient = "euc.sd"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is character string for chord coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = "my_data",
        dissimilarity_coefficient = "chord"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is character string for chisq coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = "my_data",
        dissimilarity_coefficient = "chisq"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is character string for gower coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = "my_data",
        dissimilarity_coefficient = "gower"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is character string for bray coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = "my_data",
        dissimilarity_coefficient = "bray"
      ),
      "invalid 'length' argument"
    )
  }
)
# Numeric input tests
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is numeric value for euc coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = 123,
        dissimilarity_coefficient = "euc"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is numeric value for euc.sd coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = 123,
        dissimilarity_coefficient = "euc.sd"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is numeric value for chord coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = 123,
        dissimilarity_coefficient = "chord"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is numeric value for chisq coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = 123,
        dissimilarity_coefficient = "chisq"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is numeric value for gower coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = 123,
        dissimilarity_coefficient = "gower"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is numeric value for bray coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = 123,
        dissimilarity_coefficient = "bray"
      ),
      "invalid 'length' argument"
    )
  }
)

# Empty list input tests
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is empty list for euc coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = list(),
        dissimilarity_coefficient = "euc"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is empty list for euc.sd coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = list(),
        dissimilarity_coefficient = "euc.sd"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is empty list for chord coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = list(),
        dissimilarity_coefficient = "chord"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is empty list for chisq coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = list(),
        dissimilarity_coefficient = "chisq"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is empty list for gower coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = list(),
        dissimilarity_coefficient = "gower"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is empty list for bray coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = list(),
        dissimilarity_coefficient = "bray"
      ),
      "invalid 'length' argument"
    )
  }
)
# Empty dataframe input tests
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is empty dataframe for euc coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = data.frame(),
        dissimilarity_coefficient = "euc"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is empty dataframe for euc.sd coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = data.frame(),
        dissimilarity_coefficient = "euc.sd"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is empty dataframe for chord coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = data.frame(),
        dissimilarity_coefficient = "chord"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is empty dataframe for chisq coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = data.frame(),
        dissimilarity_coefficient = "chisq"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is empty dataframe for gower coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = data.frame(),
        dissimilarity_coefficient = "gower"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is empty dataframe for bray coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = data.frame(),
        dissimilarity_coefficient = "bray"
      ),
      "invalid 'length' argument"
    )
  }
)
# Zero value input tests
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is zero for euc coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = 0,
        dissimilarity_coefficient = "euc"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is zero for euc.sd coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = 0,
        dissimilarity_coefficient = "euc.sd"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is zero for chord coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = 0,
        dissimilarity_coefficient = "chord"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is zero for chisq coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = 0,
        dissimilarity_coefficient = "chisq"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is zero for gower coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = 0,
        dissimilarity_coefficient = "gower"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is zero for bray coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = 0,
        dissimilarity_coefficient = "bray"
      ),
      "invalid 'length' argument"
    )
  }
)
# NA value input tests
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is NA for euc coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = NA,
        dissimilarity_coefficient = "euc"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is NA for euc.sd coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = NA,
        dissimilarity_coefficient = "euc.sd"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is NA for chord coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = NA,
        dissimilarity_coefficient = "chord"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is NA for chisq coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = NA,
        dissimilarity_coefficient = "chisq"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is NA for gower coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = NA,
        dissimilarity_coefficient = "gower"
      ),
      "invalid 'length' argument"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when data_source_dc is NA for bray coefficient",
  {
    expect_error(
      estimate_dissimilarity_coefficient(
        data_source_dc = NA,
        dissimilarity_coefficient = "bray"
      ),
      "invalid 'length' argument"
    )
  }
)
## Input validation: dissimilarity_coefficient parameter tests
# # Empty dissimilarity_coefficient tests
# test_that(
#   "estimate_dissimilarity_coefficient() warns when dissimilarity_coefficient is empty with proportions data",
#   {
#     suppressWarnings(
#       data_to_run_bins <-
#         extract_data(
#           data_community_extract = RRatepol::example_data$pollen_data[[1]],
#           data_age_extract = RRatepol::example_data$sample_age[[1]],
#           age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
#         ) %>%
#         smooth_community_data(
#           smooth_method = "shep"
#         ) %>%
#         reduce_data(
#           check_taxa = TRUE,
#           check_levels = TRUE
#         ) %>%
#         prepare_data(
#           data_source_prep = .,
#           working_units = "bins",
#           bin_size = 500,
#           rand = 1
#         ) %>%
#         RUtilpol::flatten_list_by_one() %>%
#         .[[1]]
#     )
#     data_source_subset <-
#       data_to_run_bins$data
#     data_source_bins <-
#       data_to_run_bins$bins
#     data_subset <-
#       subset_samples(
#         data_source_subset = data_source_subset,
#         data_source_bins = data_source_bins,
#         bin_selection = "first"
#       ) %>%
#       reduce_data_simple()
#     # standardisation
#     standardise <-
#       TRUE
#     n_individuals <-
#       150
#     com_data_sums <-
#       rowSums(
#         subset_community(
#           data_source = data_subset
#         ),
#         na.rm = TRUE
#       )
#     n_individuals <-
#       min(
#         c(
#           com_data_sums,
#           n_individuals
#         )
#       )
#     data_subset <-
#       data_subset[com_data_sums >= n_individuals, ]
#     data_subset <-
#       reduce_data_simple(
#         data_source_reduce = data_subset
#       )
#     set.seed(123)
#     data_sd_prop <-
#       standardise_community_data(
#         data_source_standard = data_subset,
#         n_individuals = n_individuals
#       ) %>%
#       reduce_data_simple(
#         data_source_reduce = .
#       ) %>%
#       transform_into_proportions(
#         data_source_trans = .,
#         sel_method = "proportions",
#         verbose = FALSE
#       )
#     expect_warning(
#       dc_res <-
#         estimate_dissimilarity_coefficient(
#           data_source_dc = data_sd_prop,
#           dissimilarity_coefficient =
#           ), # none programmed into function yet
#       # e.g., "Warning: no dissimilarity coefficient supplied. Using default 'chord' instead."
#     )
#   }
# )
# test_that(
#   "estimate_dissimilarity_coefficient() warns when dissimilarity_coefficient is empty with percentages data",
#   {
#     suppressWarnings(
#       data_to_run_bins <-
#         extract_data(
#           data_community_extract = RRatepol::example_data$pollen_data[[1]],
#           data_age_extract = RRatepol::example_data$sample_age[[1]],
#           age_uncertainty = RRatepol::example_data$age_uncertainty[[1]]
#         ) %>%
#         smooth_community_data(
#           smooth_method = "shep"
#         ) %>%
#         reduce_data(
#           check_taxa = TRUE,
#           check_levels = TRUE
#         ) %>%
#         prepare_data(
#           data_source_prep = .,
#           working_units = "bins",
#           bin_size = 500,
#           rand = 1
#         ) %>%
#         RUtilpol::flatten_list_by_one() %>%
#         .[[1]]
#     )
#     data_source_subset <-
#       data_to_run_bins$data
#     data_source_bins <-
#       data_to_run_bins$bins
#     data_subset <-
#       subset_samples(
#         data_source_subset = data_source_subset,
#         data_source_bins = data_source_bins,
#         bin_selection = "first"
#       ) %>%
#       reduce_data_simple()
#     # standardisation
#     standardise <-
#       TRUE
#     n_individuals <-
#       150
#     com_data_sums <-
#       rowSums(
#         subset_community(
#           data_source = data_subset
#         ),
#         na.rm = TRUE
#       )
#     n_individuals <-
#       min(
#         c(
#           com_data_sums,
#           n_individuals
#         )
#       )
#     data_subset <-
#       data_subset[com_data_sums >= n_individuals, ]
#     data_subset <-
#       reduce_data_simple(
#         data_source_reduce = data_subset
#       )
#     set.seed(123)
#     data_sd_prop <-
#       standardise_community_data(
#         data_source_standard = data_subset,
#         n_individuals = n_individuals
#       ) %>%
#       reduce_data_simple(
#         data_source_reduce = .
#       ) %>%
#       transform_into_proportions(
#         data_source_trans = .,
#         sel_method = "percentages",
#         verbose = FALSE
#       )
#     expect_warning(
#       dc_res <-
#         estimate_dissimilarity_coefficient(
#           data_source_dc = data_sd_prop,
#           dissimilarity_coefficient =
#           ), # none programmed into function yet
#       # e.g., "Warning: no dissimilarity coefficient supplied. Using default 'chord' instead."
#     )
#   }
# )
# NULL dissimilarity_coefficient tests
test_that(
  "estimate_dissimilarity_coefficient() throws error when dissimilarity_coefficient is NULL with proportions data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "proportions",
        verbose = FALSE
      )
    expect_error(
      dc_res <-
        estimate_dissimilarity_coefficient(
          data_source_dc = data_sd_prop,
          dissimilarity_coefficient = NULL
        ), # none programmed into function yet
      "argument is of length zero"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when dissimilarity_coefficient is NULL with percentages data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "percentages",
        verbose = FALSE
      )
    expect_error(
      dc_res <-
        estimate_dissimilarity_coefficient(
          data_source_dc = data_sd_prop,
          dissimilarity_coefficient = NULL
        ), #
      "argument is of length zero"
    )
  }
)
# Invalid character dissimilarity_coefficient tests
test_that(
  "estimate_dissimilarity_coefficient() throws error when dissimilarity_coefficient is invalid character with proportions data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "proportions",
        verbose = FALSE
      )
    expect_error(
      dc_res <-
        estimate_dissimilarity_coefficient(
          data_source_dc = data_sd_prop,
          dissimilarity_coefficient = "my_coefficient"
        ), # none programmed into function yet
      "object 'corrmat' not found"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when dissimilarity_coefficient is invalid character with percentages data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "percentages",
        verbose = FALSE
      )
    expect_error(
      dc_res <-
        estimate_dissimilarity_coefficient(
          data_source_dc = data_sd_prop,
          dissimilarity_coefficient = "my_coefficient"
        ), 
      "object 'corrmat' not found"
    )
  }
)
# Numeric dissimilarity_coefficient tests
test_that(
  "estimate_dissimilarity_coefficient() throws error when dissimilarity_coefficient is numeric with proportions data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "proportions",
        verbose = FALSE
      )
    expect_error(
      dc_res <-
        estimate_dissimilarity_coefficient(
          data_source_dc = data_sd_prop,
          dissimilarity_coefficient = 123
        ), # none programmed into function yet
      "object 'corrmat' not found"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when dissimilarity_coefficient is numeric with percentages data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "percentages",
        verbose = FALSE
      )
    expect_error(
      dc_res <-
        estimate_dissimilarity_coefficient(
          data_source_dc = data_sd_prop,
          dissimilarity_coefficient = 123
        ), 
      "object 'corrmat' not found"
    )
  }
)
# Multiple values dissimilarity_coefficient tests
test_that(
  "estimate_dissimilarity_coefficient() throws error when dissimilarity_coefficient has multiple values with proportions data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "proportions",
        verbose = FALSE
      )
    expect_error(
      dc_res <-
        estimate_dissimilarity_coefficient(
          data_source_dc = data_sd_prop,
          dissimilarity_coefficient = c("euc", "chord")
        ), # none programmed into function yet
      "the condition has length > 1"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when dissimilarity_coefficient has multiple values with percentages data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "percentages",
        verbose = FALSE
      )
    expect_error(
      dc_res <-
        estimate_dissimilarity_coefficient(
          data_source_dc = data_sd_prop,
          dissimilarity_coefficient = c("euc", "chord")
        ), #
      "the condition has length > 1"
    )
  }
)
# Zero dissimilarity_coefficient tests
test_that(
  "estimate_dissimilarity_coefficient() throws error when dissimilarity_coefficient is zero with proportions data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "proportions",
        verbose = FALSE
      )
    expect_error(
      dc_res <-
        estimate_dissimilarity_coefficient(
          data_source_dc = data_sd_prop,
          dissimilarity_coefficient = 0
        ), # none programmed into function yet
      "object 'corrmat' not found"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when dissimilarity_coefficient is zero with percentages data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "percentages",
        verbose = FALSE
      )
    expect_error(
      dc_res <-
        estimate_dissimilarity_coefficient(
          data_source_dc = data_sd_prop,
          dissimilarity_coefficient = 0
        ), #
      "object 'corrmat' not found"
    )
  }
)
# NA dissimilarity_coefficient tests
test_that(
  "estimate_dissimilarity_coefficient() throws error when dissimilarity_coefficient is NA with proportions data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "proportions",
        verbose = FALSE
      )
    expect_error(
      dc_res <-
        estimate_dissimilarity_coefficient(
          data_source_dc = data_sd_prop,
          dissimilarity_coefficient = NA
        ), # none programmed into function yet
      "missing value where TRUE/FALSE needed"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() throws error when dissimilarity_coefficient is NA with percentages data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "percentages",
        verbose = FALSE
      )
    expect_error(
      dc_res <-
        estimate_dissimilarity_coefficient(
          data_source_dc = data_sd_prop,
          dissimilarity_coefficient = NA
        ), #
      "missing value where TRUE/FALSE needed"
    )
  }
)

# Valid parameter input tests
# 1. Euclidean distance tests
test_that(
  "estimate_dissimilarity_coefficient() correctly calculates euclidean distance with proportions data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "proportions",
        verbose = FALSE
      )
    dc_res <-
      estimate_dissimilarity_coefficient(
        data_source_dc = data_sd_prop,
        dissimilarity_coefficient = "euc"
      )
    # Control euc:
    n_res <-
      nrow(
        data_sd_prop
      ) - 1
    dat_res <-
      vector(
        mode = "numeric",
        length = n_res
      )
    data_com <-
      subset_community(data_sd_prop)
    corrmat <-
      as.matrix(
        vegan::vegdist(
          data_com,
          method = "euclidean"
        )
      )
    control_res <-
      corrmat[row(corrmat) == col(corrmat) + 1]

    expect_identical(
      dc_res,
      control_res
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() correctly calculates euclidean distance with percentages data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "percentages",
        verbose = FALSE
      )
    dc_res <-
      estimate_dissimilarity_coefficient(
        data_source_dc = data_sd_prop,
        dissimilarity_coefficient = "euc"
      )
    # Control euc:
    n_res <-
      nrow(
        data_sd_prop
      ) - 1
    dat_res <-
      vector(
        mode = "numeric",
        length = n_res
      )
    data_com <-
      subset_community(data_sd_prop)
    corrmat <-
      as.matrix(
        vegan::vegdist(
          data_com,
          method = "euclidean"
        )
      )
    control_res <-
      corrmat[row(corrmat) == col(corrmat) + 1]

    expect_identical(
      dc_res,
      control_res
    )
  }
)
# 2. Euclidean distance tests
test_that(
  "estimate_dissimilarity_coefficient() correctly calculates euclidean distance with proportions data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "proportions",
        verbose = FALSE
      )
    dc_res <-
      estimate_dissimilarity_coefficient(
        data_source_dc = data_sd_prop,
        dissimilarity_coefficient = "euc.sd"
      )
    # Control euc.sd:
    n_res <-
      nrow(
        data_sd_prop
      ) - 1
    control_res <-
      vector(
        mode = "numeric",
        length = n_res
      )
    data_com <-
      subset_community(data_sd_prop)
    # calculate the SD for each species
    df_sp_supp <-
      apply(data_com, 2, stats::sd)
    # calculation of the dissimilarity_coefficient
    # for each sample (except the last)
    for (i in 1:n_res) {
      # select only 2 samples (observed + 1 after)
      df_work <-
        data_com[c(
          i, i + 1
        ), , drop = FALSE]
      # get rid of "empty species" in data & in sp.std
      df_sp_supp_work <-
        df_sp_supp[colSums(
          df_work,
          na.rm = TRUE
        ) > 0]
      df_work <-
        as.data.frame(df_work[, colSums(df_work, , na.rm = TRUE) > 0])
      # vector for result for each species
      vector_work <-
        vector(
          mode = "numeric",
          length = ncol(df_work)
        )
      # for each species
      for (j in 1:ncol(
        df_work
      )) {
        # check if the standard deviation is not equal zero
        if (
          df_sp_supp_work[j] != 0
        ) {
          a <-
            .subset2(
              df_work, j
            )[1]
          b <-
            .subset2(
              df_work, j
            )[2]
          # calculate the difference
          vector_work[j] <-
            ((a - b) / df_sp_supp_work[j])**2
        }
      }
      # save the square root of sum of all differece
      control_res[i] <-
        sqrt(sum(vector_work))
    }
    expect_identical(
      dc_res,
      control_res
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() correctly calculates euclidean standardized distance with percentages data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "percentages",
        verbose = FALSE
      )
    dc_res <-
      estimate_dissimilarity_coefficient(
        data_source_dc = data_sd_prop,
        dissimilarity_coefficient = "euc.sd"
      )
    # Control euc.sd:
    n_res <-
      nrow(
        data_sd_prop
      ) - 1
    control_res <-
      vector(
        mode = "numeric",
        length = n_res
      )
    data_com <-
      subset_community(data_sd_prop)
    # calculate the SD for each species
    df_sp_supp <-
      apply(data_com, 2, stats::sd)
    # calculation of the dissimilarity_coefficient
    # for each sample (except the last)
    for (i in 1:n_res) {
      # select only 2 samples (observed + 1 after)
      df_work <-
        data_com[c(
          i, i + 1
        ), , drop = FALSE]
      # get rid of "empty species" in data & in sp.std
      df_sp_supp_work <-
        df_sp_supp[colSums(
          df_work,
          na.rm = TRUE
        ) > 0]
      df_work <-
        as.data.frame(df_work[, colSums(df_work, , na.rm = TRUE) > 0])
      # vector for result for each species
      vector_work <-
        vector(
          mode = "numeric",
          length = ncol(df_work)
        )
      # for each species
      for (j in 1:ncol(
        df_work
      )) {
        # check if the standard deviation is not equal zero
        if (
          df_sp_supp_work[j] != 0
        ) {
          a <-
            .subset2(
              df_work, j
            )[1]
          b <-
            .subset2(
              df_work, j
            )[2]
          # calculate the difference
          vector_work[j] <-
            ((a - b) / df_sp_supp_work[j])**2
        }
      }
      # save the square root of sum of all differece
      control_res[i] <-
        sqrt(sum(vector_work))
    }
    expect_identical(
      dc_res,
      control_res
    )
  }
)
# 3. Chord distance tests
test_that(
  "estimate_dissimilarity_coefficient() correctly calculates chord distance with proportions data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "proportions",
        verbose = FALSE
      )
    dc_res <-
      estimate_dissimilarity_coefficient(
        data_source_dc = data_sd_prop,
        dissimilarity_coefficient = "chord"
      )
    # Control chord:
    n_res <-
      nrow(
        data_sd_prop
      ) - 1
    dat_res <-
      vector(
        mode = "numeric",
        length = n_res
      )
    data_com <-
      subset_community(data_sd_prop)
    corrmat <-
      as.matrix(
        vegan::vegdist(
          data_com,
          method = "chord"
        )
      )
    control_res <-
      corrmat[row(corrmat) == col(corrmat) + 1]

    expect_identical(
      dc_res,
      control_res
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() correctly calculates chord distance with percentages data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "percentages",
        verbose = FALSE
      )
    dc_res <-
      estimate_dissimilarity_coefficient(
        data_source_dc = data_sd_prop,
        dissimilarity_coefficient = "chord"
      )
    # Control chord:
    n_res <-
      nrow(
        data_sd_prop
      ) - 1
    dat_res <-
      vector(
        mode = "numeric",
        length = n_res
      )
    data_com <-
      subset_community(data_sd_prop)
    corrmat <-
      as.matrix(
        vegan::vegdist(
          data_com,
          method = "chord"
        )
      )
    control_res <-
      corrmat[row(corrmat) == col(corrmat) + 1]

    expect_identical(
      dc_res,
      control_res
    )
  }
)
# 4. Chi-square distance tests
test_that(
  "estimate_dissimilarity_coefficient() correctly calculates chi-square distance with proportions data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "proportions",
        verbose = FALSE
      )
    dc_res <-
      estimate_dissimilarity_coefficient(
        data_source_dc = data_sd_prop,
        dissimilarity_coefficient = "chisq"
      )
    # Control chisq:
    n_res <-
      nrow(
        data_sd_prop
      ) - 1
    dat_res <-
      vector(
        mode = "numeric",
        length = n_res
      )
    data_com <-
      subset_community(data_sd_prop)
    corrmat <-
      as.matrix(
        vegan::vegdist(
          data_com,
          method = "chisq"
        )
      )
    control_res <-
      corrmat[row(corrmat) == col(corrmat) + 1]

    expect_identical(
      dc_res,
      control_res
    )
  }
)

test_that(
  "estimate_dissimilarity_coefficient() correctly calculates chi-square distance with percentages data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "percentages",
        verbose = FALSE
      )
    dc_res <-
      estimate_dissimilarity_coefficient(
        data_source_dc = data_sd_prop,
        dissimilarity_coefficient = "chisq"
      )
    # Control chisq:
    n_res <-
      nrow(
        data_sd_prop
      ) - 1
    dat_res <-
      vector(
        mode = "numeric",
        length = n_res
      )
    data_com <-
      subset_community(data_sd_prop)
    corrmat <-
      as.matrix(
        vegan::vegdist(
          data_com,
          method = "chisq"
        )
      )
    control_res <-
      corrmat[row(corrmat) == col(corrmat) + 1]

    expect_identical(
      dc_res,
      control_res
    )
  }
)
# 5. Gower distance tests
test_that(
  "estimate_dissimilarity_coefficient() correctly calculates gower distance with proportions data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "proportions",
        verbose = FALSE
      )
    dc_res <-
      estimate_dissimilarity_coefficient(
        data_source_dc = data_sd_prop,
        dissimilarity_coefficient = "gower"
      )
    # Control gower:
    n_res <-
      nrow(
        data_sd_prop
      ) - 1
    dat_res <-
      vector(
        mode = "numeric",
        length = n_res
      )
    data_com <-
      subset_community(data_sd_prop)
    corrmat <-
      as.matrix(
        vegan::vegdist(
          data_com,
          method = "gower"
        )
      )
    control_res <-
      corrmat[row(corrmat) == col(corrmat) + 1]
    expect_identical(
      dc_res,
      control_res
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() correctly calculates gower distance with percentages data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "percentages",
        verbose = FALSE
      )
    dc_res <-
      estimate_dissimilarity_coefficient(
        data_source_dc = data_sd_prop,
        dissimilarity_coefficient = "gower"
      )
    # Control gower:
    n_res <-
      nrow(
        data_sd_prop
      ) - 1
    dat_res <-
      vector(
        mode = "numeric",
        length = n_res
      )
    data_com <-
      subset_community(data_sd_prop)
    corrmat <-
      as.matrix(
        vegan::vegdist(
          data_com,
          method = "gower"
        )
      )
    control_res <-
      corrmat[row(corrmat) == col(corrmat) + 1]
    expect_identical(
      dc_res,
      control_res
    )
  }
)
# 6. Bray-Curtis distance tests
test_that(
  "estimate_dissimilarity_coefficient() correctly calculates bray-curtis distance with proportions data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "proportions",
        verbose = FALSE
      )
    dc_res <-
      estimate_dissimilarity_coefficient(
        data_source_dc = data_sd_prop,
        dissimilarity_coefficient = "bray"
      )
    # Control bray:
    n_res <-
      nrow(
        data_sd_prop
      ) - 1
    dat_res <-
      vector(
        mode = "numeric",
        length = n_res
      )
    data_com <-
      subset_community(
        data_sd_prop
      )
    corrmat <-
      as.matrix(
        vegan::vegdist(
          data_com,
          method = "bray"
        )
      )
    control_res <-
      corrmat[row(corrmat) == col(corrmat) + 1]
    expect_identical(
      dc_res,
      control_res
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() correctly calculates bray-curtis distance with percentages data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "percentages",
        verbose = FALSE
      )
    dc_res <-
      estimate_dissimilarity_coefficient(
        data_source_dc = data_sd_prop,
        dissimilarity_coefficient = "bray"
      )
    # Control bray:
    n_res <-
      nrow(
        data_sd_prop
      ) - 1
    dat_res <-
      vector(
        mode = "numeric",
        length = n_res
      )
    data_com <-
      subset_community(
        data_sd_prop
      )
    corrmat <-
      as.matrix(
        vegan::vegdist(
          data_com,
          method = "bray"
        )
      )
    control_res <-
      corrmat[row(corrmat) == col(corrmat) + 1]
    expect_identical(
      dc_res,
      control_res
    )
  }
)

# Output test: vector:
# 1. Euclidean distance tests
test_that(
  "estimate_dissimilarity_coefficient() correctly calculates euclidean distance with proportions data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "proportions",
        verbose = FALSE
      )
    dc_res <-
      estimate_dissimilarity_coefficient(
        data_source_dc = data_sd_prop,
        dissimilarity_coefficient = "euc"
      )
    expect_type(
      dc_res,
      "double"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() correctly calculates euclidean distance with percentages data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "percentages",
        verbose = FALSE
      )
    dc_res <-
      estimate_dissimilarity_coefficient(
        data_source_dc = data_sd_prop,
        dissimilarity_coefficient = "euc"
      )
    expect_type(
      dc_res,
      "double"
    )
  }
)
# 2. Euclidean distance tests
test_that(
  "estimate_dissimilarity_coefficient() correctly calculates euclidean distance with proportions data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "proportions",
        verbose = FALSE
      )
    dc_res <-
      estimate_dissimilarity_coefficient(
        data_source_dc = data_sd_prop,
        dissimilarity_coefficient = "euc.sd"
      )
    expect_type(
      dc_res,
      "double"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() correctly calculates euclidean standardized distance with percentages data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "percentages",
        verbose = FALSE
      )
    dc_res <-
      estimate_dissimilarity_coefficient(
        data_source_dc = data_sd_prop,
        dissimilarity_coefficient = "euc.sd"
      )
    expect_type(
      dc_res,
      "double"
    )
  }
)
# 3. Chord distance tests
test_that(
  "estimate_dissimilarity_coefficient() correctly calculates chord distance with proportions data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "proportions",
        verbose = FALSE
      )
    dc_res <-
      estimate_dissimilarity_coefficient(
        data_source_dc = data_sd_prop,
        dissimilarity_coefficient = "chord"
      )
    expect_type(
      dc_res,
      "double"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() correctly calculates chord distance with percentages data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "percentages",
        verbose = FALSE
      )
    dc_res <-
      estimate_dissimilarity_coefficient(
        data_source_dc = data_sd_prop,
        dissimilarity_coefficient = "chord"
      )
    expect_type(
      dc_res,
      "double"
    )
  }
)
# 4. Chi-square distance tests
test_that(
  "estimate_dissimilarity_coefficient() correctly calculates chi-square distance with proportions data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "proportions",
        verbose = FALSE
      )
    dc_res <-
      estimate_dissimilarity_coefficient(
        data_source_dc = data_sd_prop,
        dissimilarity_coefficient = "chisq"
      )
    expect_type(
      dc_res,
      "double"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() correctly calculates chi-square distance with percentages data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "percentages",
        verbose = FALSE
      )
    dc_res <-
      estimate_dissimilarity_coefficient(
        data_source_dc = data_sd_prop,
        dissimilarity_coefficient = "chisq"
      )
    expect_type(
      dc_res,
      "double"
    )
  }
)
# 5. Gower distance tests
test_that(
  "estimate_dissimilarity_coefficient() correctly calculates gower distance with proportions data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "proportions",
        verbose = FALSE
      )
    dc_res <-
      estimate_dissimilarity_coefficient(
        data_source_dc = data_sd_prop,
        dissimilarity_coefficient = "gower"
      )
    expect_type(
      dc_res,
      "double"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() correctly calculates gower distance with percentages data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "percentages",
        verbose = FALSE
      )
    dc_res <-
      estimate_dissimilarity_coefficient(
        data_source_dc = data_sd_prop,
        dissimilarity_coefficient = "gower"
      )
    expect_type(
      dc_res,
      "double"
    )
  }
)
# 6. Bray-Curtis distance tests
test_that(
  "estimate_dissimilarity_coefficient() correctly calculates bray-curtis distance with proportions data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "proportions",
        verbose = FALSE
      )
    dc_res <-
      estimate_dissimilarity_coefficient(
        data_source_dc = data_sd_prop,
        dissimilarity_coefficient = "bray"
      )
    expect_type(
      dc_res,
      "double"
    )
  }
)
test_that(
  "estimate_dissimilarity_coefficient() correctly calculates bray-curtis distance with percentages data",
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
    standardise <-
      TRUE
    n_individuals <-
      150
    com_data_sums <-
      rowSums(
        subset_community(
          data_source = data_subset
        ),
        na.rm = TRUE
      )
    n_individuals <-
      min(
        c(
          com_data_sums,
          n_individuals
        )
      )
    data_subset <-
      data_subset[com_data_sums >= n_individuals, ]
    data_subset <-
      reduce_data_simple(
        data_source_reduce = data_subset
      )
    set.seed(123)
    data_sd_prop <-
      standardise_community_data(
        data_source_standard = data_subset,
        n_individuals = n_individuals
      ) %>%
      reduce_data_simple(
        data_source_reduce = .
      ) %>%
      transform_into_proportions(
        data_source_trans = .,
        sel_method = "percentages",
        verbose = FALSE
      )
    dc_res <-
      estimate_dissimilarity_coefficient(
        data_source_dc = data_sd_prop,
        dissimilarity_coefficient = "bray"
      )

    expect_type(
      dc_res,
      "double"
    )
  }
)

