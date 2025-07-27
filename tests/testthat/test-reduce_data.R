# ==================================================================== #
#                        TESTS FOR reduce_data()                      #
# ==================================================================== #
#
# This file contains comprehensive tests for the reduce_data() function
# which filters taxa and/or levels based on zero-sum criteria.
#
# Test structure:
# 1. Basic Functionality Tests
# 2. Edge Case Tests  
# 3. Input Validation Tests
#
# ==================================================================== #

# --------------------------------------------------- #
# 1. BASIC FUNCTIONALITY TESTS
# --------------------------------------------------- #

# --------------------------------------------------- #
# 1.1 Taxa Reduction Only
# --------------------------------------------------- #

test_that("reduce_data preserves rownames when filtering taxa only", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  out_taxa <-
    reduce_data(
      data_source_reduce = data_smooth,
      check_taxa = TRUE,
      check_levels = FALSE
    )

  expect_true(
    all(rownames(data_smooth$community) == rownames(out_taxa$community))
  )
})

test_that("reduce_data removes zero-sum taxa when filtering taxa only", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  out_taxa <-
    reduce_data(
      data_source_reduce = data_smooth,
      check_taxa = TRUE,
      check_levels = FALSE
    )

  expect_true(
    all(colSums(out_taxa$community, na.rm = TRUE) > 0)
  )
})

test_that("reduce_data reduces or maintains number of taxa when filtering taxa only", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  full_taxa <- ncol(data_smooth$community)

  out_taxa <-
    reduce_data(
      data_source_reduce = data_smooth,
      check_taxa = TRUE,
      check_levels = FALSE
    )

  expect_lte(
    ncol(out_taxa$community), full_taxa
  )
})

test_that("reduce_data returns valid structure when filtering taxa only", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  out_taxa <-
    reduce_data(
      data_source_reduce = data_smooth,
      check_taxa = TRUE,
      check_levels = FALSE
    )

  expect_type(
    out_taxa, "list"
  )
})

# --------------------------------------------------- #
# 1.2 Level Reduction Only
# --------------------------------------------------- #

test_that("reduce_data preserves colnames when filtering levels only", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  out_levels <-
    reduce_data(
      data_source_reduce = data_smooth,
      check_taxa = FALSE,
      check_levels = TRUE
    )

  expect_true(
    all(colnames(data_smooth$community) == colnames(out_levels$community))
  )
})

test_that("reduce_data removes zero-sum levels when filtering levels only", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  out_levels <-
    reduce_data(
      data_source_reduce = data_smooth,
      check_taxa = FALSE,
      check_levels = TRUE
    )

  expect_true(
    all(rowSums(out_levels$community, na.rm = TRUE) > 0)
  )
})

test_that("reduce_data reduces or maintains number of levels when filtering levels only", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  full_levels <- nrow(data_smooth$community)

  out_levels <-
    reduce_data(
      data_source_reduce = data_smooth,
      check_taxa = FALSE,
      check_levels = TRUE
    )

  expect_lte(
    nrow(out_levels$community), full_levels
  )
})

test_that("reduce_data maintains matching rownames between age and community when filtering levels only", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  out_levels <-
    reduce_data(
      data_source_reduce = data_smooth,
      check_taxa = FALSE,
      check_levels = TRUE
    )

  expect_equal(
    rownames(out_levels$community), rownames(out_levels$age)
  )
})

test_that("reduce_data maintains matching colnames between age_un and community when filtering levels only", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  out_levels <-
    reduce_data(
      data_source_reduce = data_smooth,
      check_taxa = FALSE,
      check_levels = TRUE
    )

  expect_equal(
    colnames(out_levels$age_un), rownames(out_levels$community)
  )
})

# --------------------------------------------------- #
# 1.3 Both Taxa and Levels Filtering
# --------------------------------------------------- #

test_that("reduce_data removes zero-sum levels when filtering both taxa and levels", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  out_both <-
    reduce_data(
      data_source_reduce = data_smooth,
      check_taxa = TRUE,
      check_levels = TRUE
    )

  expect_true(
    all(rowSums(out_both$community, na.rm = TRUE) > 0)
  )
})

test_that("reduce_data removes zero-sum taxa when filtering both taxa and levels", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  out_both <-
    reduce_data(
      data_source_reduce = data_smooth,
      check_taxa = TRUE,
      check_levels = TRUE
    )

  expect_true(
    all(colSums(out_both$community, na.rm = TRUE) > 0)
  )
})

test_that("reduce_data maintains matching rownames between age and community when filtering both", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  out_both <-
    reduce_data(
      data_source_reduce = data_smooth,
      check_taxa = TRUE,
      check_levels = TRUE
    )

  expect_equal(
    rownames(out_both$community), rownames(out_both$age)
  )
})

test_that("reduce_data maintains matching colnames between age_un and community when filtering both", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  out_both <-
    reduce_data(
      data_source_reduce = data_smooth,
      check_taxa = TRUE,
      check_levels = TRUE
    )

  expect_equal(
    colnames(out_both$age_un), rownames(out_both$community)
  )
})

# --------------------------------------------------- #
# 1.4 No Filtering
# --------------------------------------------------- #

test_that("reduce_data with both FALSE doesn't change the data", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  no_filter <-
    reduce_data(
      data_source_reduce = data_smooth,
      check_taxa = FALSE,
      check_levels = FALSE
    )

  expect_equal(
    data_smooth, no_filter
  )
})

# --------------------------------------------------- #
# 2. EDGE CASE TESTS
# --------------------------------------------------- #

# --------------------------------------------------- #
# 2.1 NULL age_un Handling
# --------------------------------------------------- #

test_that("reduce_data handles NULL age_un correctly when filtering both", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  data_null_age_un <- data_smooth

  data_null_age_un$age_un <- NULL

  expect_no_error(
    reduce_data(
      data_source_reduce = data_null_age_un,
      check_taxa = TRUE,
      check_levels = TRUE
    )
  )
})

# --------------------------------------------------- #
# 2.2 Zero Taxa Handling
# --------------------------------------------------- #

test_that("reduce_data preserves rownames when removing zero taxa", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  test_data_zeros <- data_smooth

  test_data_zeros$community$fake_taxon <- 0

  out_zeros <-
    reduce_data(
      data_source_reduce = test_data_zeros,
      check_taxa = TRUE,
      check_levels = FALSE
    )

  expect_equal(
    rownames(test_data_zeros$community), rownames(out_zeros$community)
  )
})

test_that("reduce_data reduces number of columns when removing zero taxa", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  test_data_zeros <- data_smooth

  test_data_zeros$community$fake_taxon <- 0

  out_zeros <-
    reduce_data(
      data_source_reduce = test_data_zeros,
      check_taxa = TRUE,
      check_levels = FALSE
    )

  expect_lt(
    ncol(out_zeros$community), ncol(test_data_zeros$community)
  )
})

# --------------------------------------------------- #
# 2.3 Zero Levels Handling
# --------------------------------------------------- #

test_that("reduce_data removes fake level from community when filtering levels", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  test_data_empty_levels <- data_smooth

  test_data_empty_levels$community["fake_level", ] <- 0

  test_data_empty_levels$age["fake_level", ] <- 9999

  test_data_empty_levels$age_un <-
    cbind(
      test_data_empty_levels$age_un,
      fake_level = rep(9999, nrow(test_data_empty_levels$age_un))
    )

  out_empty_levels <-
    reduce_data(
      data_source_reduce = test_data_empty_levels,
      check_taxa = FALSE,
      check_levels = TRUE
    )

  expect_false(
    "fake_level" %in% rownames(out_empty_levels$community)
  )
})

test_that("reduce_data removes fake level from age when filtering levels", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  test_data_empty_levels <- data_smooth

  test_data_empty_levels$community["fake_level", ] <- 0

  test_data_empty_levels$age["fake_level", ] <- 9999

  test_data_empty_levels$age_un <-
    cbind(
      test_data_empty_levels$age_un,
      fake_level = rep(9999, nrow(test_data_empty_levels$age_un))
    )

  out_empty_levels <-
    reduce_data(
      data_source_reduce = test_data_empty_levels,
      check_taxa = FALSE,
      check_levels = TRUE
    )

  expect_false(
    "fake_level" %in% rownames(out_empty_levels$age)
  )
})

test_that("reduce_data removes fake level from age_un when filtering levels", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  test_data_empty_levels <- data_smooth

  test_data_empty_levels$community["fake_level", ] <- 0

  test_data_empty_levels$age["fake_level", ] <- 9999

  test_data_empty_levels$age_un <-
    cbind(
      test_data_empty_levels$age_un,
      fake_level = rep(9999, nrow(test_data_empty_levels$age_un))
    )

  out_empty_levels <-
    reduce_data(
      data_source_reduce = test_data_empty_levels,
      check_taxa = FALSE,
      check_levels = TRUE
    )

  expect_false(
    "fake_level" %in% colnames(out_empty_levels$age_un)
  )
})

# --------------------------------------------------- #
# 2.4 NA Taxa Handling
# --------------------------------------------------- #

test_that("reduce_data removes NA taxon from community", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  test_data_NA_taxon <- data_smooth

  test_data_NA_taxon$community$na_taxon <- NA

  out_NA_taxon <-
    reduce_data(
      data_source_reduce = test_data_NA_taxon,
      check_taxa = TRUE,
      check_levels = FALSE
    )

  expect_false(
    "na_taxon" %in% colnames(out_NA_taxon$community)
  )
})

# --------------------------------------------------- #
# 2.5 All Taxa Zero - Taxa Filtering Only
# --------------------------------------------------- #

test_that("reduce_data returns zero columns when all taxa are zero with taxa filtering only", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  test_all_taxa_zero <- data_smooth

  test_all_taxa_zero$community[, ] <- 0

  out_all_zero_taxa <-
    reduce_data(
      data_source_reduce = test_all_taxa_zero,
      check_taxa = TRUE,
      check_levels = FALSE
    )

  expect_true(
    ncol(out_all_zero_taxa$community) == 0
  )
})

test_that("reduce_data preserves age rows when all taxa are zero with taxa filtering only", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  test_all_taxa_zero <- data_smooth

  test_all_taxa_zero$community[, ] <- 0

  out_all_zero_taxa <-
    reduce_data(
      data_source_reduce = test_all_taxa_zero,
      check_taxa = TRUE,
      check_levels = FALSE
    )

  expect_true(
    nrow(out_all_zero_taxa$age) == nrow(data_smooth$age)
  )
})

test_that("reduce_data preserves age_un columns when all taxa are zero with taxa filtering only", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  test_all_taxa_zero <- data_smooth

  test_all_taxa_zero$community[, ] <- 0

  out_all_zero_taxa <-
    reduce_data(
      data_source_reduce = test_all_taxa_zero,
      check_taxa = TRUE,
      check_levels = FALSE
    )

  expect_false(
    is.null(out_all_zero_taxa$age_un)
  )
})

# --------------------------------------------------- #
# 2.6 All Taxa Zero - Both Taxa and Levels Filtering
# --------------------------------------------------- #

test_that("reduce_data returns zero columns when all taxa are zero with both filtering", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  test_all_taxa_zero <- data_smooth

  test_all_taxa_zero$community[, ] <- 0

  out_all_zero_taxa <-
    reduce_data(
      data_source_reduce = test_all_taxa_zero,
      check_taxa = TRUE,
      check_levels = TRUE
    )

  expect_true(
    ncol(out_all_zero_taxa$community) == 0
  )
})

test_that("reduce_data returns zero rows when all taxa are zero with both filtering", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  test_all_taxa_zero <- data_smooth

  test_all_taxa_zero$community[, ] <- 0

  out_all_zero_taxa <-
    reduce_data(
      data_source_reduce = test_all_taxa_zero,
      check_taxa = TRUE,
      check_levels = TRUE
    )

  expect_true(
    nrow(out_all_zero_taxa$community) == 0
  )
})

test_that("reduce_data returns zero age_un columns when all taxa are zero with both filtering", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  test_all_taxa_zero <- data_smooth

  test_all_taxa_zero$community[, ] <- 0

  out_all_zero_taxa <-
    reduce_data(
      data_source_reduce = test_all_taxa_zero,
      check_taxa = TRUE,
      check_levels = TRUE
    )

  expect_true(
    is.null(out_all_zero_taxa$age_un) || ncol(out_all_zero_taxa$age_un) == 0
  )
})

# --------------------------------------------------- #
# 2.7 All Levels Zero Handling
# --------------------------------------------------- #

test_that("reduce_data returns zero community rows when all levels are zero", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  test_all_levels_zero <- data_smooth

  test_all_levels_zero$community[, ] <- 0

  out_all_zero_levels <-
    reduce_data(
      data_source_reduce = test_all_levels_zero,
      check_taxa = FALSE,
      check_levels = TRUE
    )

  expect_equal(
    nrow(out_all_zero_levels$community), 0
  )
})

test_that("reduce_data returns zero age rows when all levels are zero", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  test_all_levels_zero <- data_smooth

  test_all_levels_zero$community[, ] <- 0

  out_all_zero_levels <-
    reduce_data(
      data_source_reduce = test_all_levels_zero,
      check_taxa = FALSE,
      check_levels = TRUE
    )

  expect_equal(
    nrow(out_all_zero_levels$age), 0
  )
})

test_that("reduce_data returns zero age_un columns when all levels are zero", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  test_all_levels_zero <- data_smooth

  test_all_levels_zero$community[, ] <- 0

  out_all_zero_levels <-
    reduce_data(
      data_source_reduce = test_all_levels_zero,
      check_taxa = FALSE,
      check_levels = TRUE
    )

  expect_equal(
    ncol(out_all_zero_levels$age_un), 0
  )
})

# --------------------------------------------------- #
# 3. INPUT VALIDATION TESTS
# --------------------------------------------------- #

# --------------------------------------------------- #
# 3.1 Data Source Validation
# --------------------------------------------------- #

test_that("reduce_data validates NULL data", {
  expect_error(
    reduce_data(data_source_reduce = NULL)
  )
})

test_that("reduce_data validates string data input", {
  expect_error(
    reduce_data(data_source_reduce = "invalid")
  )
})

test_that("reduce_data validates numeric data input", {
  expect_error(
    reduce_data(data_source_reduce = 123)
  )
})

test_that("reduce_data validates data.frame data input", {
  expect_error(
    reduce_data(data_source_reduce = data.frame(x = 1))
  )
})

test_that("reduce_data validates empty list", {
  expect_error(
    reduce_data(data_source_reduce = list())
  )
})

# --------------------------------------------------- #
# 3.2 Required Components Validation
# --------------------------------------------------- #

test_that("reduce_data validates missing community component", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  incomplete_data <- list(age = data_smooth$age, age_un = data_smooth$age_un)
  expect_error(
    reduce_data(data_source_reduce = incomplete_data)
  )
})

test_that("reduce_data validates missing age component", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  incomplete_data <- list(community = data_smooth$community, age_un = data_smooth$age_un)
  expect_error(
    reduce_data(data_source_reduce = incomplete_data)
  )
})

# --------------------------------------------------- #
# 3.3 Community Component Validation
# --------------------------------------------------- #

test_that("reduce_data validates community is not a matrix", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  data_smooth$community <- as.matrix(data_smooth$community)

  expect_no_error(
    reduce_data(data_source_reduce = data_smooth)
  )
})

test_that("reduce_data validates community is not a list", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  data_smooth$community <- list(x = 1)

  expect_error(
    reduce_data(data_source_reduce = data_smooth)
  )
})

test_that("reduce_data validates community has rows", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  data_smooth$community <- data_smooth$community[0, ]

  expect_error(
    reduce_data(data_source_reduce = data_smooth)
  )
})

test_that("reduce_data validates community has columns", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  data_smooth$community <- data_smooth$community[, 0]

  expect_no_error(
    reduce_data(data_source_reduce = data_smooth)
  )
})

test_that("reduce_data validates community contains only numeric columns", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  data_smooth$community$text_col <- "text"

  expect_error(
    reduce_data(data_source_reduce = data_smooth)
  )
})

test_that("reduce_data validates community must have rownames", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  rownames(data_smooth$community) <- NULL

  expect_error(
    reduce_data(data_source_reduce = data_smooth)
  )
})

test_that("reduce_data validates community must have colnames", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  colnames(data_smooth$community) <- NULL

  expect_error(
    reduce_data(data_source_reduce = data_smooth)
  )
})

test_that("reduce_data validates community cannot contain Inf values", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  data_smooth$community[1, 1] <- Inf

  expect_error(
    reduce_data(data_source_reduce = data_smooth)
  )
})

test_that("reduce_data validates community cannot contain -Inf values", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  data_smooth$community[1, 1] <- -Inf

  expect_error(
    reduce_data(data_source_reduce = data_smooth)
  )
})

# --------------------------------------------------- #
# 3.4 Age Component Validation
# --------------------------------------------------- #

test_that("reduce_data validates age is a data.frame", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  data_smooth$age <- as.matrix(data_smooth$age)

  expect_no_error(
    reduce_data(data_source_reduce = data_smooth)
  )
})

test_that("reduce_data validates age contains required age column", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  data_smooth$age$age <- NULL

  expect_error(
    reduce_data(data_source_reduce = data_smooth)
  )
})

test_that("reduce_data validates age column is numeric", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  data_smooth$age$age <- as.character(data_smooth$age$age)

  expect_error(
    reduce_data(data_source_reduce = data_smooth)
  )
})

test_that("reduce_data validates age must have rownames", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  rownames(data_smooth$age) <- NULL

  expect_error(
    reduce_data(data_source_reduce = data_smooth)
  )
})

test_that("reduce_data validates age column cannot contain NA values", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  data_smooth$age$age[1] <- NA

  expect_error(
    reduce_data(data_source_reduce = data_smooth)
  )
})

test_that("reduce_data validates age column cannot contain Inf values", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  data_smooth$age$age[1] <- Inf

  expect_error(
    reduce_data(data_source_reduce = data_smooth)
  )
})

test_that("reduce_data allows negative values in age", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  data_smooth$age$age[1] <- -1000

  expect_no_error(
    reduce_data(data_source_reduce = data_smooth)
  )
})

# --------------------------------------------------- #
# 3.5 Age Uncertainty Component Validation
# --------------------------------------------------- #

test_that("reduce_data validates age_un is not a list", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  data_smooth$age_un <- list(x = 1)

  expect_error(
    reduce_data(data_source_reduce = data_smooth)
  )
})

test_that("reduce_data validates age_un is not a string", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  data_smooth$age_un <- "invalid"

  expect_error(
    reduce_data(data_source_reduce = data_smooth)
  )
})

# --------------------------------------------------- #
# 3.6 Cross-Component Validation
# --------------------------------------------------- #

test_that("reduce_data validates community and age have same row counts", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  data_smooth$age <- data_smooth$age[-1, ]

  expect_error(
    reduce_data(data_source_reduce = data_smooth)
  )
})

test_that("reduce_data validates community and age have matching row names", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  rownames(data_smooth$age)[1] <- "different_name"

  expect_error(
    reduce_data(data_source_reduce = data_smooth)
  )
})

test_that("reduce_data validates age_un dimensions match community", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  data_smooth$age_un <- data_smooth$age_un[, -1]

  expect_error(
    reduce_data(data_source_reduce = data_smooth)
  )
})

test_that("reduce_data validates age_un column names match community row names", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  colnames(data_smooth$age_un)[1] <- "different_name"

  expect_error(
    reduce_data(data_source_reduce = data_smooth)
  )
})

# --------------------------------------------------- #
# 3.7 Parameter Validation
# --------------------------------------------------- #

test_that("reduce_data validates check_taxa is not a string", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  expect_error(
    reduce_data(data_source_reduce = data_smooth, check_taxa = "invalid")
  )
})

test_that("reduce_data validates check_taxa is not numeric", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  expect_error(
    reduce_data(data_source_reduce = data_smooth, check_taxa = 123)
  )
})

test_that("reduce_data validates check_taxa is not NULL", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  expect_error(
    reduce_data(data_source_reduce = data_smooth, check_taxa = NULL)
  )
})

test_that("reduce_data validates check_taxa has length 1", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  expect_error(
    reduce_data(data_source_reduce = data_smooth, check_taxa = c(TRUE, FALSE))
  )
})

test_that("reduce_data validates check_levels is not a string", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  expect_error(
    reduce_data(data_source_reduce = data_smooth, check_levels = "invalid")
  )
})

test_that("reduce_data validates check_levels is not numeric", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  expect_error(
    reduce_data(data_source_reduce = data_smooth, check_levels = 123)
  )
})

test_that("reduce_data validates check_levels is not NULL", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  expect_error(
    reduce_data(data_source_reduce = data_smooth, check_levels = NULL)
  )
})

test_that("reduce_data validates check_levels has length 1", {
  data_smooth <-
    smooth_community_data(
      data_source_smooth =
        extract_data(
          data_community_extract = RRatepol::example_data$pollen_data[[1]],
          data_age_extract = RRatepol::example_data$sample_age[[1]],
          age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
          verbose = FALSE
        ),
      smooth_method = "m.avg",
      verbose = FALSE
    )

  expect_error(
    reduce_data(data_source_reduce = data_smooth, check_levels = c(TRUE, FALSE))
  )
})
