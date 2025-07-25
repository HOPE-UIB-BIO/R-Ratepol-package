data_source_smooth <-
  extract_data(
    RRatepol::example_data$pollen_data[[1]],
    RRatepol::example_data$sample_age[[1]]
  )

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

  # Community
  expect_false(
    is.null(
      result$community
    )
  )

  expect_s3_class(
    result$community,
    "data.frame"
  )

  expect_true(
    all(
      !is.na(
        result$community
      )
    )
  )

  # Age
  expect_false(
    is.null(
      result$age
    )
  )

  expect_s3_class(
    result$age,
    "data.frame"
  )

  expect_true(
    all(
      !is.na(
        result$age
      )
    )
  )

  # Age uncertainty
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

  # Community
  expect_false(
    is.null(
      result$community
    )
  )

  expect_s3_class(
    result$community,
    "data.frame"
  )

  expect_true(
    all(
      !is.na(
        result$community
      )
    )
  )

  # Age
  expect_false(
    is.null(
      result$age
    )
  )

  expect_s3_class(
    result$age,
    "data.frame"
  )

  expect_true(
    all(
      !is.na(
        result$age
      )
    )
  )

  # Age uncertainty
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

  # Community
  expect_false(
    is.null(
      result$community
    )
  )

  expect_s3_class(
    result$community,
    "data.frame"
  )

  expect_true(
    all(
      !is.na(
        result$community
      )
    )
  )

  # Age
  expect_false(
    is.null(
      result$age
    )
  )

  expect_s3_class(
    result$age,
    "data.frame"
  )

  expect_true(
    all(
      !is.na(
        result$age
      )
    )
  )

  # Age uncertainty
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

  # Community
  expect_false(
    is.null(
      result$community
    )
  )

  expect_s3_class(
    result$community,
    "data.frame"
  )

  expect_true(
    all(
      !is.na(
        result$community
      )
    )
  )

  # Age
  expect_false(
    is.null(
      result$age
    )
  )

  expect_s3_class(
    result$age,
    "data.frame"
  )

  expect_true(
    all(
      !is.na(
        result$age
      )
    )
  )

  # Age uncertainty
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
    is.null(
      result_1$community
    )
  )

  expect_false(
    is.null(
      result_1$age
    )
  )

  expect_true(
    all(
      !is.na(
        result_1$community
      )
    )
  )

  expect_true(
    all(
      !is.na(
        result_1$age
      )
    )
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
    is.null(
      result_2$community
    )
  )

  expect_false(
    is.null(
      result_2$age
    )
  )

  expect_true(
    all(
      !is.na(
        result_2$community
      )
    )
  )

  expect_true(
    all(
      !is.na(
        result_2$age
      )
    )
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
    is.null(
      result_3$community
    )
  )

  expect_false(
    is.null(
      result_3$age
    )
  )

  expect_true(
    all(
      !is.na(
        result_3$community
      )
    )
  )

  expect_true(
    all(
      !is.na(
        result_3$age
      )
    )
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
    is.null(
      result_4$community
    )
  )

  expect_false(
    is.null(
      result_4$age
    )
  )

  expect_true(
    all(
      !is.na(
        result_4$community
      )
    )
  )

  expect_true(
    all(
      !is.na(
        result_4$age
      )
    )
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


# -------------------------------- #
# 3. Test Error messages
# -------------------------------- #

# ------------------------------------------ #
# 3.0 shep works without parameter values for
#     smooth_n_max, smooth_age_range
# ------------------------------------------ #

test_that("shep works without additional parameters", {
  res_shep <-
    smooth_community_data(
      data_source_smooth,
      smooth_method = c("shep"),
      smooth_n_points = 5,
      smooth_n_max = NULL, # no n max
      smooth_age_range = NULL, # no age range
      round_results = FALSE,
      verbose = FALSE
    )

  expect_false(
    is.null(
      res_shep$community
    )
  )

  expect_false(
    is.null(
      res_shep$age
    )
  )

  expect_true(
    all(
      !is.na(
        res_shep$community
      )
    )
  )

  expect_true(
    all(
      !is.na(
        res_shep$age
      )
    )
  )

  expect_false(
    identical(
      data_source_smooth,
      res_shep
    )
  )

  expect_type(
    res_shep,
    "list"
  )
})


# ------------------------------------------ #
# 3.1. shep requires smooth_n_points to be
#      >= 2
# ------------------------------------------ #
test_that("shep throws error message when smooth_n_points <= 2", {
  # this should throw an error that is not programmed as
  # error message into the function (yet)
  # It is a placeholder test that fails until a less-cryptic error message
  # is added to the smooth_community_data() function

  smooth_community_data(
    data_source_smooth,
    smooth_method = c("shep"),
    smooth_n_points = 2, # even number
    smooth_n_max = 9,
    smooth_age_range = 500,
    round_results = FALSE,
    verbose = FALSE
  )
})

# ------------------------------------------ #
# 3.2. all but shep require
#      smooth_n_points to be odd
# ------------------------------------------ #
test_that("Error messages are thrown when even smooth_n_points is supplied", {
  # no error for shep
  res_shep <-
    smooth_community_data(
      data_source_smooth,
      smooth_method = c("shep"),
      smooth_n_points = 4, # even number
      smooth_n_max = 9,
      smooth_age_range = 500,
      round_results = FALSE,
      verbose = FALSE
    )

  expect_false(
    is.null(
      res_shep$community
    )
  )

  expect_false(
    is.null(
      res_shep$age
    )
  )

  expect_true(
    all(
      !is.na(
        res_shep$community
      )
    )
  )

  expect_true(
    all(
      !is.na(
        res_shep$age
      )
    )
  )

  expect_false(
    identical(
      data_source_smooth,
      res_shep
    )
  )

  expect_type(
    res_shep,
    "list"
  )

  # throws error for m.avg:
  expect_error(
    smooth_community_data(
      data_source_smooth,
      smooth_method = c("m.avg"),
      smooth_n_points = 4, # even
      smooth_n_max = 9,
      smooth_age_range = 500,
      round_results = FALSE,
      verbose = FALSE
    ),
    "'smooth_n_points' must be an odd number"
  )

  # throws error for grim:
  expect_error(
    smooth_community_data(
      data_source_smooth,
      smooth_method = c("grim"),
      smooth_n_points = 4, # even
      smooth_n_max = 9,
      smooth_age_range = 500,
      round_results = FALSE,
      verbose = FALSE
    ),
    "'smooth_n_points' must be an odd number"
  )

  # throws error for age.w
  expect_error(
    smooth_community_data(
      data_source_smooth,
      smooth_method = c("age.w"),
      smooth_n_points = 4, # even
      smooth_n_max = 9,
      smooth_age_range = 500,
      round_results = FALSE,
      verbose = FALSE
    ),
    "'smooth_n_points' must be an odd number"
  )
})


# ------------------------------------------ #
# 3.3. all but m.avg & shep require
#      smooth_age_range to be numeric
# ------------------------------------------ #
# test that smooth_age_range = "A" throws error
test_that("Error is thrown if incorrect smooth_age_range is supplied", {
  # no error for m.avg & smooth_age_range = "A"
  res_mavg <-
    smooth_community_data(
      data_source_smooth,
      smooth_method = c("m.avg"),
      smooth_n_points = 5,
      smooth_n_max = 9,
      smooth_age_range = "A",
      round_results = FALSE,
      verbose = FALSE
    )

  expect_false(
    is.null(
      res_mavg$community
    )
  )

  expect_false(
    is.null(
      res_mavg$age
    )
  )

  expect_true(
    all(
      !is.na(
        res_mavg$community
      )
    )
  )

  expect_true(
    all(
      !is.na(
        res_mavg$age
      )
    )
  )

  expect_false(
    identical(
      data_source_smooth,
      res_mavg
    )
  )

  expect_type(
    res_mavg,
    "list"
  )

  # no error for shep & smooth_age_range = "A"
  res_shep <- smooth_community_data(
    data_source_smooth,
    smooth_method = c("shep"),
    smooth_n_points = 5,
    smooth_n_max = 9,
    smooth_age_range = "A",
    round_results = FALSE,
    verbose = FALSE
  )

  expect_false(
    is.null(
      res_shep$community
    )
  )

  expect_false(
    is.null(
      res_shep$age
    )
  )

  expect_true(
    all(
      !is.na(
        res_shep$community
      )
    )
  )

  expect_true(
    all(
      !is.na(
        res_shep$age
      )
    )
  )

  expect_false(
    identical(
      data_source_smooth_NA,
      res_shep
    )
  )

  expect_type(
    res_shep,
    "list"
  )

  ## age.w throws error when smooth_age_range = "A"
  expect_error(
    smooth_community_data(
      data_source_smooth,
      smooth_method = c("age.w"),
      smooth_n_points = 5,
      smooth_n_max = 9,
      smooth_age_range = "A",
      round_results = FALSE,
      verbose = FALSE
    ),
    "'smooth_age_range' must be one of the following: 'numeric'"
  )

  ## grim throws error when smooth_age_range = "A"
  expect_error(
    smooth_community_data(
      data_source_smooth,
      smooth_method = c("grim"),
      smooth_n_points = 5,
      smooth_n_max = 9,
      smooth_age_range = "A",
      round_results = FALSE,
      verbose = FALSE
    ),
    "'smooth_age_range' must be one of the following: 'numeric'"
  )
})


# ------------------------------------------ #
# 3.4 grim requires
#     smooth_n_max to be odd
# ------------------------------------------ #
test_that("grim throws error if smoth_n_max is even", {
  expect_error(
    smooth_community_data(
      data_source_smooth,
      smooth_method = c("grim"),
      smooth_n_points = 5,
      smooth_n_max = 8,
      smooth_age_range = 500,
      round_results = FALSE,
      verbose = FALSE
    ),
    "'smooth_n_max' must be an odd number"
  )
})

# ------------------------------------------ #
# 3.5 grim requires
#     smooth_n_max to be > than smooth_n_points
# ------------------------------------------ #
test_that("grim throws error if smooth_n_points is > than smooth_n_max", {
  expect_error(
    smooth_community_data(
      data_source_smooth,
      smooth_method = c("grim"),
      smooth_n_points = 5,
      smooth_n_max = 3,
      smooth_age_range = 500,
      round_results = FALSE,
      verbose = FALSE
    ),
    "'smooth_n_max' must be bigger than 'smooth_n_points"
  )
})
