# ---------------------------------------------------------- #
#               Data_source_prep Input Tests                 #
# ---------------------------------------------------------- #

test_that(
  "prepare_data rejects NULL data_source", {
  expect_error(
    prepare_data(
      data_source_prep = NULL,
      working_units = c("levels", "bins", "MW"),
      bin_size = 500,
      number_of_shifts = 5,
      rand = NULL
    ), "'data_source_prep' must be one of the following: 'list'"
  )
})

test_that(
  "prepare_data rejects incomplete list(data.frame) in data_source", {
  expect_error(
    prepare_data(
      data_source_prep = RRatepol::example_data$pollen_data[[1]],
      working_units = c("levels", "bins", "MW"),
      bin_size = 500,
      number_of_shifts = 5,
      rand = NULL
    ), "'data_source_prep' must be one of the following: 'list'"
  )
})

test_that(
  "prepare_data fails if no community in data_source", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  example_data$community <-
    NULL
  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = c("levels", "bins", "MW"),
      bin_size = 500,
      number_of_shifts = 5,
      rand = NULL
    ), "is not TRUE"
  )
})


test_that(
  "prepare_data fails if no age in data_source", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  example_data$age <-
    NULL
  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = c("levels", "bins", "MW"),
      bin_size = 500,
      number_of_shifts = 5,
      rand = NULL
    ), "`names` must be a character vector"
  )
})

test_that(
  "prepare_data works if no age_uncertainty in data_source", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  example_data$age_un <-
    NULL
  expect_no_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = c("levels", "bins", "MW"),
      bin_size = 500,
      number_of_shifts = 5,
      rand = NULL
    )
  )
})




## Community data
# Zeros
test_that(
  "prepare_data with working_units='levels' works if there are 0s in community data", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  example_data$community[1:3] <-
    0

  expect_no_error(
    prepare_data(
      example_data,
      working_units = "levels"
    )
  )
})

test_that(
  "prepare_data with working_units='bins' works if there are 0s in community data", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  example_data$community[1:3] <-
    0

  expect_no_error(
    prepare_data(
      example_data,
      working_units = "bins",
      bin_size = 500
    )
  )
})


test_that(
  "prepare_data with working_units='MW' works if there are 0s in community data", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  example_data$community[1:3] <-
    0

  expect_no_error(
    prepare_data(
      example_data,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 5
    )
  )
})

# NAs in Age
test_that(
  "prepare_data with working_units='levels' works if there are NAs in community data", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  example_data$community[1:3] <-
    NA

  expect_no_error(
    prepare_data(
      example_data,
      working_units = "levels"
    )
  )
})

test_that(
  "prepare_data with working_units='bins' works if there are NAs in community data", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  example_data$community[1:3] <-
    NA

  expect_no_error(
    prepare_data(
      example_data,
      working_units = "bins",
      bin_size = 500
    )
  )
})


test_that(
  "prepare_data with working_units='MW' works if there are NAs in community data", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  example_data$community[1:3] <-
    NA

  expect_no_error(
    prepare_data(
      example_data,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 5
    )
  )
})


## Age data
# Zeros
test_that(
  "prepare_data with working_units='levels' works if there are 0s in age data", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  example_data$age$age[1:3] <-
    0

  expect_no_error(
    prepare_data(
      example_data,
      working_units = "levels"
    )
  )
})

test_that(
  "prepare_data with working_units='bins' works if there are 0s in age data", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  example_data$age$age[1:3] <-
    0

  expect_no_error(
    prepare_data(
      example_data,
      working_units = "bins",
      bin_size = 500
    )
  )
})


test_that(
  "prepare_data with working_units='MW' works if there are 0s in age data", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  example_data$age$age[1:3] <-
    0

  expect_no_error(
    prepare_data(
      example_data,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 5
    )
  )
})

# NAs in Age
test_that(
  "prepare_data with working_units='levels' works if there are NAs in age data", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  example_data$age$age[1:3] <-
    NA

  expect_no_error(
    prepare_data(
      example_data,
      working_units = "levels"
    )
  )
})

test_that(
  "prepare_data with working_units='bins' fails if there are NAs in age data", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  example_data$age$age[1:3] <-
    NA

  expect_error(
    prepare_data(
      example_data,
      working_units = "bins",
      bin_size = 500
    ),
    "'from' must be a finite number"
  )
})


test_that(
  "prepare_data with working_units='MW' fails if there are NAs in age data", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  example_data$age$age[1:3] <-
    NA

  expect_error(
    prepare_data(
      example_data,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 5
    ),
    "'from' must be a finite number"
  )
})


## Age_uncertainty data
# Zeros
test_that(
  "prepare_data with working_units='levels' and rand=1 works if there are 0s in the age uncertainty", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  example_data$age_un[1:3] <-
    0

  expect_no_error(
    res <-
      prepare_data(
      example_data,
      working_units = "levels",
      rand = 1
    )
  )

  expect_false(
    identical(
      res[[1]][[1]]$data$age,
      example_data$age$age
    )
  )
})

test_that(
  "prepare_data with working_units='bins' and rand=1 works if there are 0s in the age uncertainty", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  example_data$age_un[1:3] <-
    0

  expect_no_error(
    res <-
      prepare_data(
      example_data,
      working_units = "bins",
      bin_size = 500,
      rand = 1
    )
  )

  expect_false(
    identical(
      res[[1]][[1]]$data$age,
      example_data$age$age
    )
  )
})


test_that(
  "prepare_data with working_units='MW' and rand=1 works if there are 0s in the age uncertainty", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  example_data$age_un[1:3] <-
    0

  expect_no_error(
    res <-
      prepare_data(
      example_data,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 5,
      rand = 1
    )
  )

  expect_false(
    identical(
      res[[1]][[1]]$data$age,
      example_data$age$age
    )
  )
})


# NAs in Age_uncertainty
test_that(
  "prepare_data with working_units='levels' works if there are NAs in age_un data", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  example_data$age_un[1:3] <-
    NA

  expect_no_error(
    prepare_data(
      example_data,
      working_units = "levels",
      rand = 1
    )
  )
})

test_that(
  "prepare_data with working_units='bins' works if there are NAs in age_un data", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  example_data$age_un[1:3] <-
    NA

  expect_no_error(
    prepare_data(
      example_data,
      working_units = "bins",
      bin_size = 500,
      rand = 1
    )
  )
})


test_that(
  "prepare_data with working_units='MW' returns NA in data if there is NA in age_uncertainty", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  example_data$age_un[1:3] <-
    NA

  expect_no_error(
    res <-
      prepare_data(
        example_data,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 5,
        rand = 1
      )
  )

  # returns NA in output$age
  expect_true(any(is.na(res[[1]][[1]]$data$age)))
})

# ---------------------------------------------------------- #
#                  Working Units Input Tests                 #
# ---------------------------------------------------------- #

test_that(
  "prepare_data validates working_units is one of the three valid options", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "invalid",
      bin_size = 500,
      number_of_shifts = 5,
      rand = NULL
    ),
    "'working_units' must contains one of the following values: 'levels', 'bins', 'MW'"
  )
})

test_that(
  "prepare_data validates working_units parameter is not numeric", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = 123,
      bin_size = 500,
      number_of_shifts = 5,
      rand = NULL
    ),
    "'working_units' must be one of the following: 'character'"
  )
})

test_that(
  "prepare_data validates working_units parameter is not numeric", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = NA,
      bin_size = 500,
      number_of_shifts = 5,
      rand = NULL
    ),
    "'working_units' must be one of the following: 'character'"
  )
})


test_that(
  "prepare_data validates working_units parameter is not a vector of parameters", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = c("MW", "bins", "levels"),
      rand = NULL
    ),
    "'arg' must be of length 1"
  )
})

test_that(
  "prepare_data validates working_units parameter is not a vector of parameters", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = c("MW", "bins")
    ),
    "'arg' must be of length 1"
  )
})


# test fails because function silently uses the default parameters
test_that(
  "prepare_data rejects no working_unit as input", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    prepare_data(
      data_source_bins
    )
    # no error programmed into function yet
  )
})

test_that(
  "prepare_data validates working_units parameter is not 0", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = 0
    ),
    "'working_units' must be one of the following: 'character'"
  )
})


test_that(
  "prepare_data validates working_units parameter is not NULL", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = NULL
    ),
    "'working_units' must be one of the following: 'character'"
  )
})



### Bin Size ###

# ---------------------------------------------------------- #
#                   Bin_size Input Tests                     #
# ---------------------------------------------------------- #

# 1.1 "levels"

test_that(
  "prepare_data validates bin_size is not NULL with working_unit='levels'", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "levels",
      bin_size = NULL,
      number_of_shifts = NULL,
      rand = NULL
    ),
    "'bin_size' must be one of the following: 'numeric'"
  )
})

test_that(
  "prepare_data validates bin_size is not Inf with working_unit='levels'", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "levels",
      bin_size = Inf,
      number_of_shifts = NULL,
      rand = NULL
    ),
    "'to' must be a finite number"
  )
})


test_that(
  "prepare_data validates bin_size is not negative with working_unit='levels'", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "levels",
      bin_size = -500,
      number_of_shifts = NULL,
      rand = NULL
    ),
    "wrong sign in 'by' argument"
  )
})

test_that(
  "prepare_data validates bin_size is not NA with working_unit='levels'", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "levels",
      bin_size = NA,
      number_of_shifts = NULL,
      rand = NULL
    ),
    "'bin_size' must be one of the following: 'numeric'"
  )
})

test_that(
  "prepare_data validates bin_size is integer with working_unit='levels'", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "levels",
      bin_size = 500.5,
      number_of_shifts = NULL,
      rand = NULL
    ),
    "' bin_size ' must be a an integer"
  )
})

test_that(
  "prepare_data validates bin_size is not 0 with working_unit='levels'", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "levels",
      bin_size = 0,
      number_of_shifts = NULL,
      rand = NULL
    ),
    "invalid"
  )
})


test_that(
  "prepare_data works with minimum bin_size (1) with working_unit='levels'", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  expect_no_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "levels",
      bin_size = 1,
      number_of_shifts = NULL,
      rand = NULL
    )
  )
})

# 1.2 "bins"

test_that(
  "prepare_data validates bin_size is numeric", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "bins",
      bin_size = "500",
      number_of_shifts = NULL,
      rand = NULL
    ),
    "'bin_size' must be one of the following: 'numeric'"
  )
})

test_that(
  "prepare_data validates bin_size is not Inf", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "bins",
      bin_size = Inf,
      number_of_shifts = NULL,
      rand = NULL
    ),
    "'to' must be a finite number"
  )
})


test_that(
  "prepare_data validates bin_size is not negative", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "bins",
      bin_size = -500,
      number_of_shifts = NULL,
      rand = NULL
    ),
    "wrong sign in 'by' argument"
  )
})

test_that(
  "prepare_data validates bin_size is not NA", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "bins",
      bin_size = NA,
      number_of_shifts = NULL,
      rand = NULL
    ),
    "'bin_size' must be one of the following: 'numeric'"
  )
})

test_that(
  "prepare_data validates bin_size is integer", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "bins",
      bin_size = 500.5,
      number_of_shifts = NULL,
      rand = NULL
    ),
    "' bin_size ' must be a an integer"
  )
})

test_that(
  "prepare_data validates bin_size is not 0", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "bins",
      bin_size = 0,
      number_of_shifts = NULL,
      rand = NULL
    ),
    "invalid"
  )
})


test_that(
  "prepare_data works with minimum bin_size (1)", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  expect_no_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "bins",
      bin_size = 1,
      number_of_shifts = NULL,
      rand = NULL
    )
  )
})


# 1.3 "MW"

test_that(
  "prepare_data validates bin_size is numeric (working_units=WM)", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "MW",
      bin_size = "500",
      number_of_shifts = 5,
      rand = NULL
    ),
    "'bin_size' must be one of the following: 'numeric'"
  )
})

test_that(
  "prepare_data validates bin_size is not Inf (working_units=WM)", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "MW",
      bin_size = Inf,
      number_of_shifts = 5,
      rand = NULL
    ),
    "'to' must be a finite number"
  )
})


test_that(
  "prepare_data validates bin_size is not negative (working_units=WM)", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "MW",
      bin_size = -500,
      number_of_shifts = 5,
      rand = NULL
    ),
    "wrong sign in 'by' argument"
  )
})

test_that(
  "prepare_data validates bin_size is not NA (working_units=WM)", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "MW",
      bin_size = NA,
      number_of_shifts = 5,
      rand = NULL
    ),
    "'bin_size' must be one of the following: 'numeric'"
  )
})


test_that(
  "prepare_data validates bin_size is integer (working_units=WM)", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "MW",
      bin_size = 500.5,
      number_of_shifts = 5,
      rand = NULL
    ),
    "' bin_size ' must be a an integer"
  )
})


test_that(
  "prepare_data validates bin_size is not 0 (working_units=WM)", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "MW",
      bin_size = 0,
      number_of_shifts = 5,
      rand = NULL
    ),
    "invalid"
  )
})


test_that(
  "prepare_data works with minimum bin_size (1) and MW", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  expect_no_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "MW",
      bin_size = 1,
      number_of_shifts = 5,
      rand = NULL
    )
  )
})



# ---------------------------------------------------------- #
#                Number of shifts Input Tests                #
# ---------------------------------------------------------- #

test_that(
  "prepare_data validates number_of_shifts is present (working_units=WM)", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = NULL,
      rand = NULL
    ),
    "'number_of_shifts' must be one of the following: 'numeric'"
  )
})

test_that(
  "prepare_data validates number_of_shifts is not Inf (working_units=WM)", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = Inf,
      rand = NULL
    ),
    "result would be too long a vector"
  )
})


test_that(
  "prepare_data validates number_of_shifts is not negative (working_units=WM)", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = -5,
      rand = NULL
    ),
    "invalid 'times' argument"
  )
})

test_that(
  "prepare_data validates number_of_shifts is not NA (working_units=WM)", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = NA,
      rand = NULL
    ),
    "'number_of_shifts' must be one of the following: 'numeric'"
  )
})


test_that(
  "prepare_data validates number_of_shifts is integer (working_units=WM)", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 5.5,
      rand = NULL
    ),
    "' number_of_shifts ' must be a an integer"
  )
})


test_that(
  "prepare_data works with number_of_shifts = 0 (working_units=WM)", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  expect_no_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 0,
      rand = NULL
    ),
  )
})


# ---------------------------------------------------------- #
#                 Rand Input Tests                           #
# ---------------------------------------------------------- #

test_that(
  "prepare_data validates rand parameter when provided", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 5,
      rand = "invalid"
    ),
    "'rand' must be one of the following: 'NULL', 'numeric'"
  )
})

test_that(
  "prepare_data validates rand is integer when numeric", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    prepare_data(
      data_source_prep = example_data,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 5,
      rand = 5.5
    ),
    "' rand ' must be a an integer"
  )
})

test_that(
  "prepare_data produces N = rand random samples", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_no_error(
    dat_5 <-
      prepare_data(
        data_source_prep = example_data,
        working_units = "levels",
        bin_size = 500,
        number_of_shifts = 5,
        rand = 5
      )
  )

  dat_3 <-
    prepare_data(
      data_source_prep = example_data,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 5,
      rand = 3
    )

  expect_equal(length(dat_3), 3)
  expect_equal(length(dat_5), 5)
})

test_that(
  "prepare_data produces _randomized_ samples with every trial", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  dat_5 <-
    prepare_data(
      data_source_prep = example_data,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 5,
      rand = 5
    )

  dat_3 <-
    prepare_data(
      data_source_prep = example_data,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 5,
      rand = 3
    )
  expect_false(identical(dat_3, dat_5[1:3]))
})


test_that(
  "prepare_data samples reproducibly if rand > 1 and seed set manually", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  set.seed(123)
  res1 <-
    prepare_data(
      data_source_prep = example_data,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 5,
      rand = 5
    )

  set.seed(123)
  res2 <-
    prepare_data(
      data_source_prep = example_data,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 5,
      rand = 5
    )

  expect_identical(res1, res2)
})

#fails currently
test_that(
  "prepare_data samples reproducibly if rand > 1 and no seed set manually", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  res1 <-
    prepare_data(
      data_source_prep = example_data,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 5,
      rand = 1
    )

  res2 <-
    prepare_data(
      data_source_prep = example_data,
      working_units = "levels",
      bin_size = 500,
      number_of_shifts = 5,
      rand = 1
    )

  expect_identical(res1, res2)
})

# fails currentöy
test_that(
  "prepare_data does not change age if no age_uncertainty in data", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = NULL,
      verbose = FALSE
    )

  res <-
    prepare_data(
      data_source_prep = example_data,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1,
      rand = 1
    )

  expect_identical(
    res[[1]][[1]]$bins$res_age,
    example_data$age$age
  )
})



# WIP:
# ---------------------------------------------------------- #
#                  Output structure Tests                    #
# ---------------------------------------------------------- #

# with valid input:
test_that(
  "prepare_data returns the correct output structure with default parameters (levels)", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  res <-
    1
  
    prepare_data(
      data_source_prep = example_data,
      working_units = "levels",
      bin_size = NULL,
      number_of_shifts = NULL,
      rand = NULL
    )

  # General output structure tests:
  expect_type(
    res, "list"
  )

  # Named list elements
  expect_identical(
    names(
    res[[1]][[1]]),
    c("data", "bins")
  )

  # Named data object
  expect_identical(
    names(
    res[[1]][[1]]$data),
    c("age", names(example_data$community))
  )

  # Named bins object
  expect_true(
    all(
      c("name", "shift", "age_diff", "start", "end", "res_age", "label")
      %in%
        colnames(res[[1]][[1]]$bins)
    )
  )

  # $data is not all 0 or empty
  expect_false(all(res[[1]][[1]]$data == 0))
  expect_false(length(res[[1]][[1]]$data) == 0)
  # $bins is not all 0 or empty
  expect_false(all(res[[1]][[1]]$bins == 0))
  expect_false(nrow(res[[1]][[1]]$bins) == 0)
})



test_that(
  "prepare_data returns the correct output structure with default parameters (MW)", {
  example_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  res <-
    prepare_data(
      data_source_prep = example_data,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 5,
      rand = NULL
    )

  # General output structure tests:
  expect_type(
    res, "list"
  )

  # Named list elements
  expect_identical(
    names(
    res[[1]][[1]]),
    c("data", "bins")
  )

  # Named data object
  expect_identical(
    names(
    res[[1]][[1]]$data),
    c("age", names(example_data$community))
  )

  # Named bins object
  expect_true(
    all(
      c("name", "shift", "age_diff", "start", "end", "res_age", "label")
      %in%
        colnames(res[[1]][[1]]$bins)
    )
  )

  # $data is not all 0 or empty
  expect_false(all(res[[1]][[1]]$data == 0))
  expect_false(length(res[[1]][[1]]$data) == 0)
  # $bins is not all 0 or empty
  expect_false(all(res[[1]][[1]]$bins == 0))
  expect_false(nrow(res[[1]][[1]]$bins) == 0)
})




















# General output structure tests:
expect_type(
  res, "list"
)

# Named list elements
expect_identical(
  names(
    res[[1]][[1]]),
  c("data", "bins")
)

# Named data object
expect_identical(
  names(
    res[[1]][[1]]$data),
  c("age", names(example_data$community))
)

# Named bins object
expect_true(
  all(
    c("name", "shift", "age_diff", "start", "end", "res_age", "label")
    %in%
      colnames(res[[1]][[1]]$bins)
  )
)

# $data is not all 0 or empty
expect_false(all(res[[1]][[1]]$data == 0))
expect_false(length(res[[1]][[1]]$data) == 0)
# $bins is not all 0 or empty
expect_false(all(res[[1]][[1]]$bins == 0))
expect_false(nrow(res[[1]][[1]]$bins) == 0)
