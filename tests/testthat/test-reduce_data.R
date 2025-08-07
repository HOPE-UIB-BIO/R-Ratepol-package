# ==================================================================== #
#                        TESTS FOR reduce_data()                      #
# ==================================================================== #
#
# This file contains comprehensive tests for the reduce_data() function
# which filters taxa and/or levels based on zero-sum criteria.
#
# Test structure:
# 1. Input Validation Tests (Errors)
# 2. No Filtering Tests (check_taxa = FALSE, check_levels = FALSE)
# 3. Taxa-Only Filtering Tests (check_taxa=TRUE, check_levels=FALSE)
# 4. Levels-Only Filtering Tests (check_taxa=FALSE, check_levels=TRUE)
# 5. Both Filtering Tests (check_taxa=TRUE, check_levels=TRUE)
# 6. All-zero / All-NA / Empty data cases
# 7. Output validation
#
# ==================================================================== #

# --------------------------------------------------- #
# 1. INPUT VALIDATION TESTS - ERRORS
# --------------------------------------------------- #

# 1.1 data_source_reduce validation
test_that(
  "reduce_data rejects NULL data_source_reduce", {
  expect_error(
    reduce_data(
    data_source_reduce = NULL),
    "data_source_reduce.*list"
  )
})

test_that(
  "reduce_data rejects string data_source_reduce", {
  expect_error(
    reduce_data(
    data_source_reduce = "invalid"),
    "data_source_reduce.*list"
  )
})

test_that(
  "reduce_data rejects numeric data_source_reduce", {
  expect_error(
    reduce_data(
    data_source_reduce = 123),
    "data_source_reduce.*list"
  )
})

test_that(
  "reduce_data rejects data.frame data_source_reduce", {
  expect_error(
    reduce_data(
    data_source_reduce = data.frame(
    x = 1)),
    "data_source_reduce.*list"
  )
})

test_that(
  "reduce_data rejects empty list data_source_reduce", {
  expect_error(
    reduce_data(
    data_source_reduce = list(
    )),
    "'x' must be an array of at least two dimensions"
  )
})

# 1.2 check_taxa validation
test_that(
  "reduce_data rejects invalid check_taxa argument", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  expect_error(
    reduce_data(
    data_source_reduce = raw_data, check_taxa = "TRUE"),
    "check_taxa.*logical"
  )
})

test_that(
  "reduce_data rejects invalid check_taxa argument", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  expect_error(
    reduce_data(
    data_source_reduce = raw_data, check_taxa = 123),
    "check_taxa.*logical"
  )
})

test_that(
  "reduce_data rejects invalid check_taxa argument", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  expect_error(
    reduce_data(
    data_source_reduce = raw_data, check_taxa = NULL),
    "check_taxa.*logical"
  )
})

# 1.3 check_levels validation
test_that(
  "reduce_data rejects invalid check_levels argument", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  expect_error(
    reduce_data(
    data_source_reduce = raw_data, check_levels = "TRUE"),
    "check_levels.*logical"
  )
})

test_that(
  "reduce_data rejects invalid check_levels argument", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  expect_error(
    reduce_data(
    data_source_reduce = raw_data, check_levels = 123),
    "check_levels.*logical"
  )
})

test_that(
  "reduce_data rejects invalid check_levels argument", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  expect_error(
    reduce_data(
      data_source_reduce = raw_data,
      check_levels = NULL
    ),
    "check_levels.*logical"
  )
})

# 1.4 data_source_reduce structure validation
test_that(
  "reduce_data rejects missing community component", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  incomplete_data <-
    list(age = raw_data$age, age_un = raw_data$age_un)
  expect_error(
    reduce_data(
    data_source_reduce = incomplete_data),
    "'x' must be an array of at least two dimensions"
  )
})

test_that(
  "reduce_data rejects missing age component", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  incomplete_data <-
    list(community = raw_data$community, age_un = raw_data$age_un)
  expect_error(
    reduce_data(
    data_source_reduce = incomplete_data),
    "Empty name found at location 1"
  )
})

test_that(
  "reduce_data rejects community as matrix", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  raw_data$community <-
    as.matrix(raw_data$community)
  expect_error(
    reduce_data(
    data_source_reduce = raw_data),
    "no applicable method for 'select'"
  )
})

test_that(
  "reduce_data rejects age as matrix", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  raw_data$age <-
    as.matrix(raw_data$age)

  expect_error(
    reduce_data(
    data_source_reduce = raw_data),
    "is not TRUE"
  )
})

test_that(
  "reduce_data rejects community as list", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  raw_data$community <-
    list(community = raw_data$community)

  expect_error(
    reduce_data(
    data_source_reduce = raw_data),
    "no applicable method for 'select'"
  )
})

test_that(
  "reduce_data rejects age as list", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  raw_data$age <-
    list(age = raw_data$age)

  expect_error(
    reduce_data(
    data_source_reduce = raw_data),
    "no applicable method for 'select'"
  )
})

test_that(
  "reduce_data validates age_un is not a list", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  raw_data$age_un <-
    list(age_un = raw_data$age_un)

  expect_error(
    reduce_data(
    data_source_reduce = raw_data),
    "incorrect number of dimensions"
  )
})


test_that(
  "reduce_data rejects community with non-numeric columns", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  raw_data$community$text_col <-
    "NA"
  expect_error(
    reduce_data(
    data_source_reduce = raw_data),
    "'x' must be numeric"
  )
})

test_that(
  "reduce_data validates community must have colnames", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  colnames(raw_data$community) <-
    NULL

  expect_error(
    reduce_data(
    data_source_reduce = raw_data),
    "Can't select within an unnamed vector."
  )
})

test_that(
  "reduce_data throws warning if community without rownames", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  rownames(raw_data$community) <-
    NULL
  expect_condition(
    result <-
      reduce_data(data_source_reduce = raw_data)
  )
})

test_that(
  "reduce_data throws warning if age without rownames", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  rownames(raw_data$age) <-
    NULL
  expect_warning(
    reduce_data(data_source_reduce = raw_data)
  )
})

# --------------------------------------------------- #
# 2. NO FILTERING TESTS (check_taxa=FALSE, check_levels=FALSE)
# --------------------------------------------------- #

test_that(
  "reduce_data with no filtering returns identical data", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  # one zero column in community
  raw_data$community[, 1] <-
    0
  # one zero row in community
  raw_data$community[1, ] <-
    0

  result <-
    reduce_data(
      data_source_reduce = raw_data,
      check_taxa = FALSE,
      check_levels = FALSE
    )

  expect_identical(raw_data, result)
})

# --------------------------------------------------- #
# 3. TAXA-ONLY FILTERING TESTS (check_taxa=TRUE, check_levels=FALSE)
# --------------------------------------------------- #

test_that(
  "taxa filtering drops taxa with zero column sums", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  # one zero column/taxa in community
  raw_data$community[, 1] <-
    0

  result <-
    reduce_data(
      data_source_reduce = raw_data,
      check_taxa = TRUE,
      check_levels = FALSE
    )

  expect_true(all(colSums(result$community, na.rm = TRUE) > 0))
  expect_equal(ncol(result$community), ncol(raw_data$community) - 1)
})


test_that(
  "check_levels = FALSE preserves rownames / sample ids", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  # one zero column in community
  raw_data$community[, 1] <-
    0
  # one zero row in community
  raw_data$community[1, ] <-
    0

  result <-
    reduce_data(
      data_source_reduce = raw_data,
      check_taxa = TRUE,
      check_levels = FALSE
    )
  expect_equal(nrow(raw_data$community), nrow(result$community))
  expect_identical(rownames(raw_data$community), rownames(result$community))
  expect_identical(rownames(raw_data$age), rownames(result$age))
  expect_identical(colnames(raw_data$age_un), colnames(result$age_un))
})

test_that(
  "taxa filtering removes all-NA taxa", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  raw_data$community$na_taxon <-
    NA

  result <-
    reduce_data(
    data_source_reduce = raw_data,
    check_taxa = TRUE,
    check_levels = FALSE
  )

  expect_false("na_taxon" %in% colnames(result$community))
})

# --------------------------------------------------- #
# 4. LEVELS-ONLY FILTERING TESTS (check_taxa=FALSE, check_levels=TRUE)
# --------------------------------------------------- #

test_that(
  "levels filtering removes zero-sum levels from community and age", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  # one zero column in community
  raw_data$community[, 1] <-
    0
  # one zero row in community
  raw_data$community[1, ] <-
    0

  first_level <-
    rownames(raw_data$community)[1]
  raw_data$community[1, ] <-
    0

  result <-
    reduce_data(
    data_source_reduce = raw_data,
    check_taxa = FALSE,
    check_levels = TRUE
  )

  expect_true(all(rowSums(result$community, na.rm = TRUE) > 0))
  expect_false(first_level %in% rownames(result$community))
  expect_false(first_level %in% rownames(result$age))
  expect_false(first_level %in% colnames(result$age_un))
})

test_that(
  "levels filtering preserves community colnames/taxa", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  # one zero column in community
  raw_data$community[, 1] <-
    0
  # one zero row in community
  raw_data$community[1, ] <-
    0

  result <-
    reduce_data(
    data_source_reduce = raw_data,
    check_taxa = FALSE,
    check_levels = TRUE
  )

  expect_identical(colnames(raw_data$community), colnames(result$community))
})

test_that(
  "levels filtering maintains matching rownames between community and age and age_un", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  # one zero column in community
  raw_data$community[, 1] <-
    0
  # one zero row in community
  raw_data$community[1, ] <-
    0

  result <-
    reduce_data(
    data_source_reduce = raw_data,
    check_taxa = FALSE,
    check_levels = TRUE
  )

  expect_identical(rownames(result$community), rownames(result$age))
  expect_identical(rownames(result$community), colnames(result$age_un))
})


test_that(
  "levels filtering works with NULL age_un", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      verbose = FALSE
    )

  raw_data$community[1, ] <-
    0

  expect_no_error(
    reduce_data(
      data_source_reduce = raw_data,
      check_taxa = FALSE,
      check_levels = TRUE
    )
  )
})

test_that(
  "levels filtering removes all-NA sample", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  first_level <-
    rownames(raw_data$community[1, ])
  raw_data$community[1, ] <-
    NA

  result <-
    reduce_data(
    data_source_reduce = raw_data,
    check_taxa = FALSE,
    check_levels = TRUE
  )

  expect_false(first_level %in% rownames(result$community))
})

# --------------------------------------------------- #
# 5. BOTH FILTERING TESTS (check_taxa=TRUE, check_levels=TRUE)
# --------------------------------------------------- #

test_that(
  "both filtering removes zero-sum taxa", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  # one zero column in community
  raw_data$community[, 1] <-
    0
  # one zero row in community
  raw_data$community[1, ] <-
    0

  result <-
    reduce_data(
    data_source_reduce = raw_data,
    check_taxa = TRUE,
    check_levels = TRUE
  )

  expect_false("zero_taxon" %in% colnames(result$community))
  expect_true(all(colSums(result$community, na.rm = TRUE) > 0))
})

test_that(
  "both filtering removes zero-sum levels", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  first_level <-
    rownames(raw_data$community)[1]

  # one zero column in community
  raw_data$community[, 1] <-
    0
  # one zero row in community
  raw_data$community[1, ] <-
    0

  result <-
    reduce_data(
    data_source_reduce = raw_data,
    check_taxa = TRUE,
    check_levels = TRUE
  )

  expect_false(first_level %in% rownames(result$community))
  expect_true(all(rowSums(result$community, na.rm = TRUE) > 0))
  expect_false(first_level %in% rownames(result$community))
  expect_false(first_level %in% rownames(result$age))
  expect_false(first_level %in% colnames(result$age_un))
})

test_that(
  "both filtering removes samples with zero taxa", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  # one zero column in community
  raw_data$community[, 1] <-
    0
  # one zero row in community
  raw_data$community[1, ] <-
    0

  result <-
    reduce_data(
    data_source_reduce = raw_data,
    check_taxa = TRUE,
    check_levels = TRUE
  )

  expect_true(all(rowSums(result$community, na.rm = TRUE) > 0))
})

test_that(
  "both filtering maintains matching rownames between community and age / age_un", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  # one zero column in community
  raw_data$community[, 1] <-
    0
  # one zero row in community
  raw_data$community[1, ] <-
    0

  result <-
    reduce_data(
    data_source_reduce = raw_data,
    check_taxa = TRUE,
    check_levels = TRUE
  )

  expect_identical(rownames(result$community), rownames(result$age))
  expect_identical(rownames(result$community), colnames(result$age_un))
})

test_that(
  "both filtering handles NULL age_un", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  # one zero column in community
  raw_data$community$zero_taxon <-
    0

  # one zero row in community
  raw_data$community[1, ] <-
    0

  # no uncertainty data
  raw_data$age_un <-
    NULL

  expect_no_error(
    result <-
      reduce_data(
      data_source_reduce = raw_data,
      check_taxa = TRUE,
      check_levels = TRUE
    )
  )

  expect_type(result, "list")
  expect_named(
    result,
    c("community", "age")
  )
})


# WIP below


test_that(
  "reduce_data works with a single sample in community", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  raw_data$community <-
    raw_data$community[1, ]

  result <-
    reduce_data(raw_data)

  expect_equal(nrow(result$community), 1)
  expect_equal(nrow(result$age), 1)
  expect_equal(ncol(result$age_un), 1)
})


test_that(
  "reduce_data works with a single taxon in community", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]][1:2],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  expect_no_error(
    result <-
      reduce_data(raw_data)
  )
})



# --------------------------------------------------- #
# 6. All Zero Handling
# --------------------------------------------------- #

test_that(
  "taxa filtering with all-zero community data returns empty community result", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  raw_data$community[, ] <-
    0

  result <-
    reduce_data(
      data_source_reduce = raw_data,
      check_taxa = TRUE,
      check_levels = FALSE
    )

  expect_equal(ncol(result$community), 0)
  expect_false(nrow(result$age) == 0)
  expect_false(ncol(result$age_un) == 0)
})

test_that(
  "levels filtering with all zero levels returns empty result", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  raw_data$community[, ] <-
    0

  result <-
    reduce_data(
    data_source_reduce = raw_data,
    check_taxa = FALSE,
    check_levels = TRUE
  )

  expect_equal(nrow(result$community), 0)
  expect_equal(nrow(result$age), 0)
  expect_equal(ncol(result$age_un), 0)
})

test_that(
  "both filtering with all zero data returns empty result", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  # all zero in community
  raw_data$community[, ] <-
    0


  result <-
    reduce_data(
    data_source_reduce = raw_data,
    check_taxa = TRUE,
    check_levels = TRUE
  )

  expect_equal(ncol(result$community), 0)
  expect_equal(nrow(result$community), 0)
  expect_equal(nrow(result$age), 0)
  expect_equal(ncol(result$age_un), 0)
})


test_that(
  "reduce_data validates community has columns", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  # 0 columns in community
  raw_data$community <-
    raw_data$community[, 0]

  result <-
    reduce_data(data_source_reduce = raw_data)

  expect_equal(ncol(result$community), 0)
  expect_equal(nrow(result$community), 0)
  expect_equal(nrow(result$age), 0)
  expect_equal(ncol(result$age_un), 0)
})


test_that(
  "both filtering returns empty result if all data are NA", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )


  raw_data$community[, ] <-
    NA

  result <-
    reduce_data(
    data_source_reduce = raw_data,
    check_taxa = TRUE,
    check_levels = TRUE
  )

  expect_equal(nrow(result$community), 0)
  expect_equal(nrow(result$age), 0)
  expect_equal(ncol(result$age_un), 0)
})


# --------------------------------------------------- #
# 7. OUTPUT VALIDATION
# --------------------------------------------------- #


test_that(
  "reduce_data produces valid output structure", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  # one zero column in community
  raw_data$community[, 1] <-
    0
  # one zero row in community
  raw_data$community[1, ] <-
    0

  result <-
    reduce_data(
    data_source_reduce = raw_data,
    check_taxa = TRUE,
    check_levels = TRUE
  )

  expect_type(
    result, "list"
  )

  expect_named(
    result,
    c("community", "age", "age_un")
  )
})

test_that(
  "reduce_data produces valid community structure", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  # one zero column in community
  raw_data$community[, 1] <-
    0
  # one zero row in community
  raw_data$community[1, ] <-
    0

  result <-
    reduce_data(
    data_source_reduce = raw_data,
    check_taxa = TRUE,
    check_levels = TRUE
  )

  expect_s3_class(
    result$community,
    "data.frame"
  )

  if (ncol(
    result$community) > 0) {
    expect_true(all(vapply(result$community, is.numeric, logical(1))))
  }

  expect_true(
    !is.null(
      rownames(
        result$community
      )
    )
  )

  expect_true(
    !is.null(
      colnames(
        result$community
      )
    )
  )
})


test_that(
  "reduce_data produces valid age structure", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  # one zero column in community
  raw_data$community[, 1] <-
    0
  # one zero row in community
  raw_data$community[1, ] <-
    0

  result <-
    reduce_data(
    data_source_reduce = raw_data,
    check_taxa = TRUE,
    check_levels = TRUE
  )

  expect_s3_class(
    result$age,
    "data.frame"
  )

  expect_true(
    "age" %in% colnames(result$age)
  )

  expect_true(
    !is.null(
      rownames(result$age)
    )
  )

  expect_equal(
    rownames(
    result$age),
    rownames(result$community)
  )
})

test_that(
  "reduce_data produces valid age_un structure", {
  raw_data <-
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  # one zero column in community
  raw_data$community[, 1] <-
    0
  # one zero row in community
  raw_data$community[1, ] <-
    0

  result <-
    reduce_data(
    data_source_reduce = raw_data,
    check_taxa = TRUE,
    check_levels = TRUE
  )

  if (!is.null(
    result$age_un)) {
    expect_true(
      is.matrix(result$age_un) || is.data.frame(result$age_un)
    )
  }

  expect_equal(
    colnames(result$age_un), rownames(result$age)
  )
})
