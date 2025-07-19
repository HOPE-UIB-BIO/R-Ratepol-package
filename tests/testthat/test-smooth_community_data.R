data_source_smooth <-
  extract_data(
    RRatepol::example_data$pollen_data[[1]],
    RRatepol::example_data$sample_age[[1]]
  )

result <-
  smooth_community_data(
    data_source_smooth,
    smooth_method = c("m.avg", "grim", "age.w", "shep"),
    smooth_n_points = 5,
    smooth_n_max = 9,
    smooth_age_range = 500,
    round_results = FALSE,
    verbose = FALSE
  )

str(result)


# -------------------------------- #
# 1. Test output structure:
# -------------------------------- #

test_that("smoothing returns expected output structure", {
  # m.avg
  result <-
    smooth_community_data(
      data_source_smooth,
      smooth_method = c("m.avg"),
      smooth_n_points = 5,
      smooth_n_max = 9,
      smooth_age_range = 500,
      round_results = FALSE,
      verbose = FALSE
    )
  
  expect_type(
    result,
    "list"
  )
  
  expect_named(
    result,
    c("community", "age", "age_un")
  )
  
  expect_s3_class(
    result$community,
    "data.frame"
  )
  
  expect_s3_class(
    result$age,
    "data.frame"
  )
  
  if (!is.null(result$age_un)) expect_s3_class(result$age_un, "data.frame")
  
  # grim
  result <-
    smooth_community_data(
      data_source_smooth,
      smooth_method = c("grim"),
      smooth_n_points = 5,
      smooth_n_max = 9,
      smooth_age_range = 500,
      round_results = FALSE,
      verbose = FALSE
    )
  
  expect_type(
    result,
    "list"
  )
  
  expect_named(
    result,
    c("community", "age", "age_un")
  )
  
  expect_s3_class(
    result$community,
    "data.frame"
  )
  
  expect_s3_class(
    result$age,
    "data.frame"
  )
  
  if (!is.null(result$age_un)) expect_s3_class(result$age_un, "data.frame")
  
  # age.w
  result <-
    smooth_community_data(
      data_source_smooth,
      smooth_method = c("age.w"),
      smooth_n_points = 5,
      smooth_n_max = 9,
      smooth_age_range = 500,
      round_results = FALSE,
      verbose = FALSE
    )
  
  expect_type(
    result,
    "list"
  )
  
  expect_named(
    result,
    c("community", "age", "age_un")
  )
  
  expect_s3_class(
    result$community,
    "data.frame"
  )
  
  expect_s3_class(
    result$age,
    "data.frame"
  )
  
  if (!is.null(result$age_un)) expect_s3_class(result$age_un, "data.frame")
  
  # shep
  result <-
    smooth_community_data(
      data_source_smooth,
      smooth_method = c("shep"),
      smooth_n_points = 5,
      smooth_n_max = 9,
      smooth_age_range = 500,
      round_results = FALSE,
      verbose = FALSE
    )
  
  expect_type(
    result,
    "list"
  )
  
  expect_named(
    result,
    c("community", "age", "age_un")
  )
  
  expect_s3_class(
    result$community,
    "data.frame"
  )
  
  expect_s3_class(
    result$age,
    "data.frame"
  )
  
  if (!is.null(result$age_un)) expect_s3_class(result$age_un, "data.frame")
})