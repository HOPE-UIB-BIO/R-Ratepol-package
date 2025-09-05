# Parameters in the function:

# # data_source_community,
# # data_source_age,
# # age_uncertainty = NULL,
# # smooth_method = c("none", "m.avg", "grim", "age.w", "shep"),
# # smooth_n_points = 5,
# # smooth_age_range = 500,
# # smooth_n_max = 9,
# # working_units = c("levels", "bins", "MW"),
# # bin_size = 500,
# # number_of_shifts = 5,
# # bin_selection = c("random", "first"),
# # standardise = FALSE,
# # n_individuals = 150,
# # dissimilarity_coefficient = c("euc", "euc.sd", "chord", "chisq", "gower", "bray"),
# # tranform_to_proportions = TRUE,
# # rand = NULL,
# # use_parallel = FALSE,
# # interest_threshold = NULL,
# # time_standardisation = NULL,
# # verbose = FALSE


test_that("estimate_roc works with valid inputs", {
  # Create mock data
  data_source_community <-
    RRatepol::example_data$pollen_data[[1]]
  data_source_age <-
    RRatepol::example_data$sample_age[[1]]
  age_uncertainty <-
    RRatepol::example_data$age_uncertainty[[1]]
  # Run the function with valid inputs
  result <-
    estimate_roc(
      data_source_community = data_source_community,
      data_source_age = data_source_age,
      age_uncertainty = age_uncertainty,
      smooth_method = "none",
      smooth_n_points = NULL,
      smooth_age_range = NULL,
      smooth_n_max = NULL,
      working_units = "levels",
      bin_size = NULL,
      number_of_shifts = NULL,
      bin_selection = "first",
      standardise = TRUE,
      n_individuals = 150,
      dissimilarity_coefficient = "euc",
      tranform_to_proportions = TRUE,
      rand = 100,
      use_parallel = FALSE,
      interest_threshold = NULL,
      time_standardisation = 500,
      verbose = FALSE
    )

  # Check that the result is a list and contains expected elements
  expect_s3_class(
    result,
    "data.frame"
  )
})
