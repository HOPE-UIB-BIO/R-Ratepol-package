# ---------------------------------------------------------- #
#  1. Input validation (Errors)
# ---------------------------------------------------------- #

# ---------------------------------------------------------- #
#               data_source_bins validation                  #
# ---------------------------------------------------------- #

test_that("make_bins rejects missing age data in data_source_bins", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  data_source_bins$age <-
    NULL

  expect_error(
    make_bins(
      data_source_bins,
      working_units = c("levels", "bins", "MW"),
      bin_size = 500,
      number_of_shifts = 5
    ),
    "`names` must be a character vector"
  )
})

test_that("make_bins validates data input class == list", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  class(data_source_bins) <-
    "data.frame"

  expect_error(
    make_bins(
      data_source_bins,
      working_units = c("levels", "bins", "MW"),
      bin_size = 500,
      number_of_shifts = 5
    ),
    "'data_source_bins' must be one of the following: 'list'"
  )
})

test_that("make_bins rejects data input class == character", {
  expect_error(
    make_bins(
      "data",
      working_units = c("levels", "bins", "MW"),
      bin_size = 500,
      number_of_shifts = 5
    ),
    "'data_source_bins' must be one of the following: 'list'"
  )
})

test_that("make_bins rejects data input class == numeric", {
  expect_error(
    make_bins(
      123,
      working_units = c("levels", "bins", "MW"),
      bin_size = 500,
      number_of_shifts = 5
    ),
    "'data_source_bins' must be one of the following: 'list'"
  )
})

# ---------------------------------------------------------- #
#               working_unit validation                      #
# ---------------------------------------------------------- #

test_that("make_bins validates that there is a working_unit as input", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    make_bins(
      data_source_bins,
      working_units = NULL
    ), "'working_units' must be one of the following: 'character'"
  )
})

## (fails - it silently uses "bins" and default bin_size = 500
test_that("make_bins rejects no working_unit as input", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    make_bins(
      data_source_bins,
    ), "'working_units' must be one of the following: 'character'"
  )
})

test_that("make_bins validates that working_unit is not NULL", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    make_bins(
      data_source_bins,
      working_units = NULL
    ), "'working_units' must be one of the following: 'character'"
  )
})

test_that("make_bins validates correct working_unit value as input", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    make_bins(
      data_source_bins,
      working_units = 123
    ), "'working_units' must be one of the following: 'character'"
  )
})

test_that("make_bins validates correct working_unit value as input", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    make_bins(
      data_source_bins,
      working_units = "invalid input"
    ), "'working_units' must contains one of the following values: 'levels', 'bins', 'MW'"
  )
})

test_that("make_bins validates correct working_unit value as input", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    make_bins(
      data_source_bins,
      working_units = c("levels", "bins")
    ), "'arg' must be of length 1"
  )
})

test_that("make_bins validates correct working_unit value as input", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    make_bins(
      data_source_bins,
      working_units = c("levels", "bins"),
      bin_size = 500,
      number_of_shifts = 5
    ), "'arg' must be of length 1"
  )
})

## (fails right now - it will take always the first method in the vector without error or warning)
test_that("make_bins rejects multiple working_units as input", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  expect_error(
    bins <-
      make_bins(
        data_source_bins,
        working_units = c("levels", "bins", "MW"),
        bin_size = 500,
        number_of_shifts = 5
      )
    # no error/warning programmed into the function yet
  )
})

test_that("make_bins rejects NULL values for all parameters (but data_source_bins)", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      verbose = FALSE
    )

  expect_error(
    make_bins(
      data_source_bins = data_source_bins,
      working_units = NULL,
      bin_size = NULL,
      number_of_shifts = NULL
    ),
    "'working_units' must be one of the following: 'character'"
  )
})





# ---------------------------------------------------------- #
#               working_units = "levels"                     #
# ---------------------------------------------------------- #

test_that("make_bins with working_units='levels' accepts data without age uncertainty", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      verbose = FALSE
    )
  expect_no_error(
    make_bins(
      data_source_bins,
      working_units = "levels",
      bin_size = NULL,
      number_of_shifts = NULL
    )
  )
})

test_that("make_bins with working_units='levels' works if there are 0s in age data", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  data_source_bins$age$age[1:3] <- 0

  expect_no_error(
    bins <-
      make_bins(
        data_source_bins,
        working_units = "levels"
      )
  )
})

## (test fails - function works and produces invalid output)
test_that("make_bins with working_units='levels' fails if there are no rownames in age data", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  rownames(data_source_bins$age) <- NULL

  expect_error(
    bins <-
      make_bins(
        data_source_bins,
        working_units = "levels",
        bin_size = NULL,
        number_of_shifts = NULL
      ),
    # no error/warning programmed into the function yet
  )

  expect_false(bins$label[1] == "...1-...2")
})

test_that("make_bins with working_units='levels' works if rownames in age are not numeric", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  length_levels <-
    length(rownames(data_source_bins$age))
  rownames(data_source_bins$age) <-
    paste0("ABC", 1:length_levels)

  expect_no_error(
    bins <-
      make_bins(
        data_source_bins,
        working_units = "levels",
        bin_size = NULL,
        number_of_shifts = NULL
      )
  )
})

test_that("make_bins with working_units='levels' ignores NAs in age data", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  data_source_bins$age$age[1:3] <- NA

  expect_no_error(
    bins <-
      make_bins(
        data_source_bins,
        working_units = "levels",
        bin_size = NULL,
        number_of_shifts = NULL
      )
  )
})


# ---------------------------------------------------------- #
#               working_units = "bins"                     #
# ---------------------------------------------------------- #

test_that("make_bins with working_units='bins' accepts data without age uncertainty", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      verbose = FALSE
    )
  expect_no_error(
    make_bins(
      data_source_bins,
      working_units = "bins",
      bin_size = 500,
      number_of_shifts = NULL
    )
  )
})

test_that("make_bins with working_units='bins' works if there are 0s in age data", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  data_source_bins$age$age[1:3] <- 0

  expect_no_error(
    bins <-
      make_bins(
        data_source_bins,
        working_units = "bins",
        bin_size = 500,
      )
  )
})

test_that("make_bins fails with working_units='bins' and no bin_size", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    bins <-
      make_bins(
        data_source_bins,
        working_units = "bins",
        bin_size = NULL,
        number_of_shifts = NULL
      ),
    "'bin_size' must be one of the following: 'numeric'"
  )
})

test_that("make_bins rejects non-numeric bin_size if working_units='bins'", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    bins <-
      make_bins(
        data_source_bins,
        working_units = "bins",
        bin_size = NA,
        number_of_shifts = NULL
      ),
    "'bin_size' must be one of the following: 'numeric'"
  )
})



test_that("make_bins works with minimum valid bin_size if working_units='bins'", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_no_error(
    bins <-
      make_bins(
        data_source_bins,
        working_units = "bins",
        bin_size = 1,
        number_of_shifts = NULL
      )
  )
})


test_that("make_bins works with maximum valid bin_size if working_units='bins'", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_no_error(
    bins <-
      make_bins(
        data_source_bins,
        working_units = "bins",
        bin_size = 2000000,
        number_of_shifts = NULL
      )
  )
})


test_that("make_bins rejects Inf bin_size if working_units='bins'", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    bins <-
      make_bins(
        data_source_bins,
        working_units = "bins",
        bin_size = Inf,
        number_of_shifts = NULL
      ),
    "'to' must be a finite number"
  )
})

test_that("make_bins rejects negative bin_size if working_units='bins'", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    bins <-
      make_bins(
        data_source_bins,
        working_units = "bins",
        bin_size = -500,
        number_of_shifts = NULL
      ),
    "wrong sign in 'by' argument"
  )
})

test_that("make_bins rejects non-integer bin_size if working_units='bins'", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    bins <-
      make_bins(
        data_source_bins,
        working_units = "bins",
        bin_size = 500.5,
        number_of_shifts = NULL
      ),
    "' bin_size ' must be a an integer"
  )
})


test_that("make_bins rejects vector of multiple bin_sizes if working_units='bins'", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    bins <-
      make_bins(
        data_source_bins,
        working_units = "bins",
        bin_size = c(100:500),
        number_of_shifts = NULL
      ),
    "'bin_size' must be one of the following: 'numeric'"
  )
})

test_that("make_bins throws error if 'bins' and there are NAs in age data", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  data_source_bins$age$age[1:3] <- NA

  expect_error(
    bins <-
      make_bins(
        data_source_bins,
        working_units = "bins",
        bin_size = 500,
        number_of_shifts = NULL
      ),
    "'from' must be a finite number"
  )
})

test_that("make_bins throws error if working_units='bins' and bin_size < 1", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    make_bins(
      data_source_bins,
      working_units = "bins",
      bin_size = 0,
    ),
    "invalid" # invalid '(to - from)/by'
  )
})


# ---------------------------------------------------------- #
#               working_units = "MW"                     #
# ---------------------------------------------------------- #

test_that("make_bins with working_units='MW' accepts data without age uncertainty", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      verbose = FALSE
    )

  expect_no_error(
    make_bins(
      data_source_bins,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1
    )
  )
})

test_that("make_bins with working_units='MW' throws error if there are NAs in age data", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  data_source_bins$age$age[1:3] <- NA

  expect_error(
    bins <-
      make_bins(
        data_source_bins,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 5
      ),
    "'from' must be a finite number"
  )
})

test_that("make_bins with working_units='MW' works if there are 0s in age data", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  data_source_bins$age$age[1:3] <- 0

  expect_no_error(
    bins <-
      make_bins(
        data_source_bins,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = 5
      )
  )
})


test_that("make_bins throws error if working_units='MW' and bin_size is missing", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )


  expect_error(
    bins <-
      make_bins(
        data_source_bins,
        working_units = "MW",
        bin_size = NULL,
        number_of_shifts = NULL
      ),
    "'bin_size' must be one of the following: 'numeric'"
  )
})


test_that("make_bins throws error if working_units='MW' and number_of_shifts is missing", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    bins <-
      make_bins(
        data_source_bins,
        working_units = "MW",
        bin_size = 500,
        number_of_shifts = NULL
      ),
    "'number_of_shifts' must be one of the following: 'numeric'"
  )
})

test_that("make_bins throws error if working_units='MW' and bin_size is NA", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )


  expect_error(
    bins <-
      make_bins(
        data_source_bins,
        working_units = "MW",
        bin_size = NA,
        number_of_shifts = 1
      ),
    "'bin_size' must be one of the following: 'numeric'"
  )
})

test_that("make_bins uses default parameters if working_units='MW' and no input bin_size and number_of_shifts", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  bins_res <-
    make_bins(
      data_source_bins,
      working_units = "MW"
    )

  bins_default <-
    make_bins(
      data_source_bins,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 5
    )

  expect_identical(bins_default, bins_res)
})

test_that("make_bins throws error if working_units='MW' and non-integer bin_size", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    make_bins(
      data_source_bins,
      working_units = "MW",
      bin_size = 500.5,
      number_of_shifts = 5
    ),
    "' bin_size ' must be a an integer"
  )
})

test_that("make_bins throws error if working_units='MW' and bin_size < 1", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    make_bins(
      data_source_bins,
      working_units = "MW",
      bin_size = 0,
      number_of_shifts = 5
    ),
    "invalid" # invalid '(to - from)/by'
  )
})

test_that("make_bins throws error if working_units='MW' and number_of_shifts < 1", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    make_bins(
      data_source_bins,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 0
    ),
    "arguments imply differing number of rows: 0, 36" # arguments imply differing number of rows: 0, 36
  )
})

test_that("make_bins throws no error if working_units='MW' and number_of_shifts = 1", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_no_error(
    make_bins(
      data_source_bins,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 1
    )
  )
})

test_that("make_bins throws error if working_units='MW' and number_of_shifts = Inf", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    make_bins(
      data_source_bins,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = Inf
    ),
    "result would be too long a vector"
  )
})

test_that("make_bins throws error if working_units='MW' and number_of_shifts = negative", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    make_bins(
      data_source_bins,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = -5
    ),
    "invalid 'times' argument"
  )
})

test_that("make_bins throws error if working_units='MW' and number_of_shifts = NA", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    make_bins(
      data_source_bins,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = NA
    ),
    "'number_of_shifts' must be one of the following: 'numeric'"
  )
})


test_that("make_bins throws error if working_units='MW' and number_of_shifts = non-numeric", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_error(
    make_bins(
      data_source_bins,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = "5"
    ),
    "'number_of_shifts' must be one of the following: 'numeric'"
  )
})

# ---------------------------------------------------------- #
# 2. Output validation
# ---------------------------------------------------------- #

test_that("make_bins returns the correct output", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  bins <-
    make_bins(data_source_bins = data_source_bins)

  # General output structure tests:
  expect_s3_class(
    bins, "data.frame"
  )

  expect_identical(
    bins$name,
    rownames(data_source_bins$age)
  )

  expect_true(
    all(
      c("name", "shift", "age_diff", "start", "end", "res_age", "label")
      %in%
        colnames(bins)
    )
  )

  expect_true(
    class(bins$name) == "character"
  )
  expect_true(
    class(bins$shift) == "numeric"
  )
  expect_true(
    class(bins$age_diff) == "numeric"
  )
  expect_true(
    class(bins$start) == "character"
  )
  expect_true(
    class(bins$end) == "character"
  )
  expect_true(
    class(bins$res_age) == "numeric"
  )
  expect_true(
    class(bins$label) == "character"
  )
})




test_that("make_bins returns the correct output without age_uncertainty", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      verbose = FALSE
    )

  bins <-
    make_bins(data_source_bins = data_source_bins)

  # General output structure tests:
  expect_s3_class(
    bins, "data.frame"
  )

  expect_identical(
    bins$name,
    rownames(data_source_bins$age)
  )

  expect_true(
    all(
      c("name", "shift", "age_diff", "start", "end", "res_age", "label")
      %in%
        colnames(bins)
    )
  )

  expect_true(
    class(bins$name) == "character"
  )
  expect_true(
    class(bins$shift) == "numeric"
  )
  expect_true(
    class(bins$age_diff) == "numeric"
  )
  expect_true(
    class(bins$start) == "character"
  )
  expect_true(
    class(bins$end) == "character"
  )
  expect_true(
    class(bins$res_age) == "numeric"
  )
  expect_true(
    class(bins$label) == "character"
  )
})


test_that("make_bins returns the correct output with default parameters for 'levels'", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  bins <-
    make_bins(
      data_source_bins = data_source_bins,
      working_units = "levels"
    )

  # General output structure tests:
  expect_s3_class(
    bins, "data.frame"
  )

  expect_identical(
    bins$name,
    rownames(data_source_bins$age)
  )

  expect_true(
    all(
      c("name", "shift", "age_diff", "start", "end", "res_age", "label")
      %in%
        colnames(bins)
    )
  )

  expect_true(
    class(bins$name) == "character"
  )
  expect_true(
    class(bins$shift) == "numeric"
  )
  expect_true(
    class(bins$age_diff) == "numeric"
  )
  expect_true(
    class(bins$start) == "character"
  )
  expect_true(
    class(bins$end) == "character"
  )
  expect_true(
    class(bins$res_age) == "numeric"
  )
  expect_true(
    class(bins$label) == "character"
  )
})

test_that("make_bins with working_units='levels' ignores NAs in age data and returns valid output", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  data_source_bins$age$age[1:3] <- NA

  bins <-
    make_bins(
      data_source_bins,
      working_units = "levels",
      bin_size = NULL,
      number_of_shifts = NULL
    )

  # General output structure tests:
  expect_s3_class(
    bins, "data.frame"
  )

  expect_identical(
    bins$name,
    rownames(data_source_bins$age)
  )

  expect_true(
    all(
      c("name", "shift", "age_diff", "start", "end", "res_age", "label")
      %in%
        colnames(bins)
    )
  )

  expect_true(
    class(bins$name) == "character"
  )
  expect_true(
    class(bins$shift) == "numeric"
  )
  expect_true(
    class(bins$age_diff) == "numeric"
  )
  expect_true(
    class(bins$start) == "character"
  )
  expect_true(
    class(bins$end) == "character"
  )
  expect_true(
    class(bins$res_age) == "numeric"
  )
  expect_true(
    class(bins$label) == "character"
  )
})



test_that("make_bins returns the correct output with default parameters for 'bins'", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  bins <-
    make_bins(
      data_source_bins = data_source_bins,
      working_units = "bins",
      bin_size = 500
    )


  # General output structure tests:
  expect_s3_class(
    bins, "data.frame"
  )

  expect_true(
    all(
      c("name", "shift", "age_diff", "start", "end", "res_age", "label")
      %in%
        colnames(bins)
    )
  )

  expect_true(
    class(bins$name) == "numeric"
  )
  expect_true(
    class(bins$shift) == "numeric"
  )
  expect_true(
    class(bins$age_diff) == "numeric"
  )
  expect_true(
    class(bins$start) == "numeric"
  )
  expect_true(
    class(bins$end) == "numeric"
  )
  expect_true(
    class(bins$res_age) == "numeric"
  )
  expect_true(
    class(bins$label) == "character"
  )
})

test_that("make_bins returns the correct output with default parameters for 'MW'", {
  data_source_bins <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  bins <-
    make_bins(
      data_source_bins = data_source_bins,
      working_units = "MW",
      bin_size = 500,
      number_of_shifts = 5
    )


  # General output structure tests:
  expect_s3_class(
    bins, "data.frame"
  )

  expect_true(
    all(
      c("name", "shift", "age_diff", "start", "end", "res_age", "label")
      %in%
        colnames(bins)
    )
  )

  expect_true(
    class(bins$name) == "numeric"
  )
  expect_true(
    class(bins$shift) == "integer"
  )
  expect_true(
    class(bins$age_diff) == "numeric"
  )
  expect_true(
    class(bins$start) == "numeric"
  )
  expect_true(
    class(bins$end) == "numeric"
  )
  expect_true(
    class(bins$res_age) == "numeric"
  )
  expect_true(
    class(bins$label) == "character"
  )
})
