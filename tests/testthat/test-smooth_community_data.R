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


# -------------------------------- #
# 2. Test NA behaviour
# -------------------------------- #

test_that("NAs are handled correctly during smoothing", {
  data_NA <-
    example_data$pollen_data[[1]]
  
  data_NA[c(1:10), 2] <- NA
  
  
  age_NA <-
    RRatepol::example_data$sample_age[[1]]
  
  age_NA[c(1:4), 3] <- NA
  
  data_source_smooth_NA <-
    extract_data(
      data_NA,
      age_NA
    )
  
  
  # m.avg
  result_1 <-
    smooth_community_data(
      data_source_smooth_NA,
      smooth_method = c("m.avg"),
      smooth_n_points = 5,
      smooth_n_max = 9,
      smooth_age_range = 500,
      round_results = FALSE,
      verbose = FALSE
    )
  
  expect_false(
    identical(
      data_source_smooth_NA,
      result_1
    )
  )
  
  
  expect_type(
    result_1,
    "list"
  )
  
  expect_named(
    result_1,
    c("community", "age", "age_un")
  )
  
  expect_s3_class(
    result_1$community,
    "data.frame"
  )
  
  expect_s3_class(
    result_1$age,
    "data.frame"
  )
  
  if (!is.null(result_1$age_un)) expect_s3_class(result_1$age_un, "data.frame")
  
  # grim
  result_2 <-
    smooth_community_data(
      data_source_smooth_NA,
      smooth_method = c("grim"),
      smooth_n_points = 5,
      smooth_n_max = 9,
      smooth_age_range = 500,
      round_results = FALSE,
      verbose = FALSE
    )
  
  expect_false(
    identical(
      data_source_smooth_NA,
      result_2
    )
  )
  
  expect_type(
    result_2,
    "list"
  )
  
  expect_named(
    result_2,
    c("community", "age", "age_un")
  )
  
  expect_s3_class(
    result_2$community,
    "data.frame"
  )
  
  expect_s3_class(
    result_2$age,
    "data.frame"
  )
  
  if (!is.null(result_2$age_un)) expect_s3_class(result_2$age_un, "data.frame")
  
  
  # age.w
  result_3 <-
    smooth_community_data(
      data_source_smooth_NA,
      smooth_method = c("age.w"),
      smooth_n_points = 5,
      smooth_n_max = 9,
      smooth_age_range = 500,
      round_results = FALSE,
      verbose = FALSE
    )
  
  expect_false(
    identical(
      data_source_smooth_NA,
      result_3
    )
  )
  expect_type(
    result_3,
    "list"
  )
  
  expect_named(
    result_3,
    c("community", "age", "age_un")
  )
  
  expect_s3_class(
    result_3$community,
    "data.frame"
  )
  
  expect_s3_class(
    result_3$age,
    "data.frame"
  )
  
  if (!is.null(result_3$age_un)) expect_s3_class(result_3$age_un, "data.frame")
  
  
  # shep
  result_4 <-
    smooth_community_data(
      data_source_smooth_NA,
      smooth_method = c("shep"),
      smooth_n_points = 5,
      smooth_n_max = 9,
      smooth_age_range = 500,
      round_results = FALSE,
      verbose = FALSE
    )
  
  expect_false(
    identical(
      data_source_smooth_NA,
      result_4
    )
  )
  
  expect_type(
    result_4,
    "list"
  )
  
  expect_named(
    result_4,
    c("community", "age", "age_un")
  )
  
  expect_s3_class(
    result_4$community,
    "data.frame"
  )
  
  expect_s3_class(
    result_4$age,
    "data.frame"
  )
  
  if (!is.null(result_4$age_un)) expect_s3_class(result_4$age_un, "data.frame")
})

