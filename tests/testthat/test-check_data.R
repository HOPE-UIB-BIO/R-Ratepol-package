# Extract example data for testing
example_community <-
  RRatepol::example_data$pollen_data[[1]]

example_age <-
  RRatepol::example_data$sample_age[[1]]

example_uncertainty <-
  RRatepol::example_data$age_uncertainty[[1]]

result <- extract_data(example_community,
                       example_age,
                       example_uncertainty)


test_that("returns NA in min, max, mean, median age if there is NA in the data", {
  # introduce NAs:
  example_community_NA <-
    example_community # copy of community
  
  example_community_NA[5:10, 3] <-
    NA
  
  example_age_NA <-
    example_age
  
  example_age_NA[5:10, 3] <-
    NA
  
  example_uncertainty_NA <-
    example_uncertainty
  
  example_uncertainty_NA[5:10, 5] <-
    NA
  
  result_NA <- extract_data(example_community_NA,
                            example_age_NA,
                            example_uncertainty_NA)
  
  expect_message(
    check_data(result_NA),
    regexp = "Age data has values of min NA, max NA, mean NA, and median NA"
  )
})
