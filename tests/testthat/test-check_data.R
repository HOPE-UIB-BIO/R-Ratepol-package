
test_that("does not return NA in min, max, mean, median age if there is NA in the data", {
  # introduce NAs:
  example_community_NA <-
    RRatepol::example_data$pollen_data[[1]]

  example_community_NA[5:10, 3] <-
    NA

  example_age_NA <-
    RRatepol::example_data$sample_age[[1]]

  example_age_NA[5:10, 3] <-
    NA

  example_uncertainty_NA <-
    RRatepol::example_data$age_uncertainty[[1]]

  example_uncertainty_NA[5:10, 5] <-
    NA

  result_NA <-
    extract_data(
      example_community_NA,
      example_age_NA,
      example_uncertainty_NA
    )

  msg <- 
    capture.output(
      check_data(result_NA),
      type = "message"
      )
  
  expect_false(
    any(
      grepl(
        "Age data has values of min NA, max NA, mean NA, and median NA", 
        msg
        )
      )
    )
})

