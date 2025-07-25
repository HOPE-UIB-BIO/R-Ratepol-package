# --------------------------------------------------- #
# Load and prepare example smoothed data
# --------------------------------------------------- #

data_smooth <-
  extract_data(
    RRatepol::example_data$pollen_data[[1]],
    RRatepol::example_data$sample_age[[1]],
    RRatepol::example_data$age_uncertainty[[1]]
  ) %>%
  smooth_community_data(
    smooth_method = c("m.avg"),
    smooth_n_points = 5,
    smooth_n_max = 9,
    smooth_age_range = 500,
    round_results = FALSE,
    verbose = FALSE
  )


# Store full taxa and levels
full_taxa <-
  colnames(data_smooth$community)

full_levels <-
  rownames(data_smooth$community)

# --------------------------------------------------- #
# Helper function to validate output structure
# --------------------------------------------------- #

check_reduce_data_structure <- function(out) {
  # ---- Top-level structure ----
  expect_type(
    out, "list"
  )
  
  expect_named(
    out, 
    c("community", "age", "age_un")
  )
  
  # ---- community ----
  expect_s3_class(
    out$community, 
    "data.frame"
  )
  
  if (ncol(out$community) > 0) {
    expect_true(all(vapply(out$community, is.numeric, logical(1))))
  }
  
  expect_true(
    !is.null(
      rownames(
        out$community)
    )
  )
  
  expect_true(
    !is.null(
      colnames(
        out$community)
    )
  )
  
  # ---- age ----
  expect_s3_class(
    out$age, 
    "data.frame"
  )
  
  expect_true(
    "age" %in% colnames(out$age)
  )
  
  expect_true(
    !is.null(
      rownames(out$age)
    )
  )
  
  expect_equal(
    rownames(out$age), 
    rownames(out$community)
  )
  
  # ---- age_un ----
  if (!is.null(out$age_un)) {
    expect_true(
      is.matrix(out$age_un) || is.data.frame(out$age_un)
    )
    
    expect_equal(
      colnames(out$age_un), rownames(out$community)
    )
  }
}

# --------------------------------------------------- #
# Modified data for edge cases
# --------------------------------------------------- #

# 1 taxon with zeros
test_data_zeros <-
  data_smooth

test_data_zeros$community$fake_taxon <-
  0

# all taxa zero
test_all_taxa_zero <-
  data_smooth

test_all_taxa_zero$community[, ] <-
  0

# 1 taxon with NA in community
test_data_NA_taxon <-
  data_smooth

test_data_NA_taxon$community$na_taxon <-
  NA

# 1 level with zeros
test_data_empty_levels <-
  data_smooth

test_data_empty_levels$community["fake_level", ] <-
  0

test_data_empty_levels$age["fake_level", ] <-
  9999

test_data_empty_levels$age_un <-
  cbind(
    test_data_empty_levels$age_un,
    "fake_level" = rep(0, nrow(test_data_empty_levels$age_un))
  )

# all levels zero
test_all_levels_zero <-
  data_smooth

test_all_levels_zero$community[, ] <-
  0

# without age_un matrix
data_null_age_un <-
  data_smooth

data_null_age_un$age_un <-
  NULL







# --------------------------------------------------- #
# Test 1: test taxa reduction only
# --------------------------------------------------- #

test_that("reduce_data filters taxa correctly using smoothed example data", {
  out_taxa <-
    reduce_data(
      data_smooth,
      check_taxa = TRUE,
      check_levels = FALSE
    )

  expect_true(
    all(
      rownames(out_taxa$community) == rownames(data_smooth$community)
    )
  )

  expect_true(
    all(
      colSums(out_taxa$community, na.rm = TRUE) > 0
    )
  )

  expect_lte(
    ncol(out_taxa$community), length(full_taxa)
  )
  
  check_reduce_data_structure(out_taxa)
  
})


# --------------------------------------------------- #
# Test 2: test level reduction only
# --------------------------------------------------- #
test_that("reduce_data filters levels correctly using smoothed example data", {
  out_levels <-
    reduce_data(
      data_smooth,
      check_taxa = FALSE,
      check_levels = TRUE
    )

  expect_true(
    all(
      colnames(out_levels$community) == colnames(data_smooth$community)
    )
  )

  expect_true(
    all(rowSums(out_levels$community, na.rm = TRUE) > 0)
  )

  expect_lte(
    nrow(out_levels$community), length(full_levels)
  )

  expect_equal(
    rownames(out_levels$age), rownames(out_levels$community)
  )

  expect_equal(
    colnames(out_levels$age_un), rownames(out_levels$community)
  )
  
  check_reduce_data_structure(out_levels)
  
})

# --------------------------------------------------- #
# Test 3: test both at the same time
# --------------------------------------------------- #
test_that("reduce_data filters taxa and levels correctly using smoothed example data", {
  out_both <-
    reduce_data(
      data_smooth,
      check_taxa = TRUE,
      check_levels = TRUE
    )

  expect_true(
    all(
      rowSums(out_both$community, na.rm = TRUE) > 0
    )
  )

  expect_true(
    all(
      colSums(out_both$community, na.rm = TRUE) > 0
    )
  )

  expect_equal(
    rownames(out_both$age), rownames(out_both$community)
  )

  expect_equal(
    colnames(out_both$age_un), rownames(out_both$community)
  )

  # Test: NULL age_un is handled
  expect_no_error(
    reduce_data(
      data_null_age_un
    )
  )
  
  check_reduce_data_structure(out_both)
})




# --------------------------------------------------- #
# Test 4: FALSE input doesn't reduce the data
# --------------------------------------------------- #
test_that("reduce_data with both = FALSE doesn't change the data", {
  # Control: no filtering needed
  no_filter <-
    reduce_data(
      data_smooth,
      check_taxa = FALSE,
      check_levels = FALSE
    )

  expect_equal(
    no_filter,
    data_smooth
  )
  
  check_reduce_data_structure(no_filter)
})



# --------------------------------------------------- #
# Test 5: all-zero taxa
# --------------------------------------------------- #
test_that("reduce_data handles taxa with all-zeros in samples correctly", {
  out_zeros <-
    reduce_data(
      test_data_zeros,
      check_taxa = TRUE,
      check_levels = FALSE
    )

  expect_equal(
    rownames(out_zeros$community),
    rownames(test_data_zeros$community)
  )

  expect_lt(
    ncol(out_zeros$community),
    ncol(test_data_zeros$community)
  )
  
  check_reduce_data_structure(out_zeros)
  
})


# --------------------------------------------------- #
# Test 6: all-zero levels
# --------------------------------------------------- #
test_that("reduce_data handles levels with all-zeros in samples correctly", {
  out_empty_levels <-
    reduce_data(
      test_data_empty_levels,
      check_taxa = FALSE,
      check_levels = TRUE
    )

  expect_false(
    "fake_level" %in% rownames(out_empty_levels$community)
  )

  expect_false(
    "fake_level" %in% rownames(out_empty_levels$age)
  )

  expect_false(
    "fake_level" %in% colnames(out_empty_levels$age_un)
  )
  
  check_reduce_data_structure(out_empty_levels)
  
  
})
# --------------------------------------------------- #
# Test 7: NA taxa
# --------------------------------------------------- #

test_that("reduce_data handles NA in taxon (community) correctly", {
  out_NA_taxon <-
    reduce_data(
      test_data_NA_taxon,
      check_taxa = TRUE,
      check_levels = FALSE
    )

  expect_false(
    "na_taxon" %in% colnames(out_NA_taxon$community)
  )
  
  check_reduce_data_structure(out_NA_taxon)
  
})

# --------------------------------------------------- #
# Test 8: if all Taxa are zero
# --------------------------------------------------- #
test_that("reduce_data handles all-zero community data correctly", {
  # Test: all taxa zero — everything from community data should be removed
  # only taxa (no level reduction)
  out_all_zero_taxa <-
    reduce_data(
      test_all_taxa_zero,
      check_taxa = TRUE,
      check_levels = FALSE
    )

  expect_true(
    ncol(
      out_all_zero_taxa$community
    ) == 0
  )

  expect_true(
    nrow(
      out_all_zero_taxa$age
    ) != 0
  )


  expect_false(
    ncol(
      out_all_zero_taxa$age_un
    ) == 0
  )


  # with levels
  out_all_zero_taxa <-
    reduce_data(
      test_all_taxa_zero,
      check_taxa = TRUE,
      check_levels = TRUE
    )

  expect_true(
    ncol(
      out_all_zero_taxa$community
    ) == 0
  )

  expect_true(
    nrow(
      out_all_zero_taxa$age
    ) == 0
  )

  expect_true(
    ncol(
      out_all_zero_taxa$age_un
    ) == 0
  )
  
  check_reduce_data_structure(out_all_zero_taxa)
  
})

# --------------------------------------------------- #
# Test 9: if all levels are zero
# --------------------------------------------------- #
test_that("reduce_data handles all-zero levels correctly", {
  out_all_zero_levels <-
    reduce_data(
      test_all_levels_zero,
      check_taxa = FALSE,
      check_levels = TRUE
    )

  expect_equal(
    nrow(
      out_all_zero_levels$community
    ), 0
  )

  expect_equal(
    nrow(
      out_all_zero_levels$age
    ), 0
  )

  expect_equal(
    ncol(
      out_all_zero_levels$age_un
    ), 0
  )
  
  check_reduce_data_structure(out_all_zero_levels)
  
})

# Test 10a: Correct output structure
test_that("reduce_data returns correct structure and format", {

  out <- 
    reduce_data(
      data_smooth, 
      check_taxa = TRUE, 
      check_levels = TRUE)
  
  check_reduce_data_structure(out)
  
})

# Test 10b: Correct output structure
test_that("reduce_data returns correct structure and format", {
  
  out <- 
    reduce_data(
      data_smooth, 
      check_taxa = FALSE, 
      check_levels = TRUE)
  
  check_reduce_data_structure(out)
  
})

# Test 10c: Correct output structure
test_that("reduce_data returns correct structure and format", {
  
  out <- 
    reduce_data(
      data_smooth, 
      check_taxa = TRUE, 
      check_levels = FALSE)
  
  check_reduce_data_structure(out)
  
})

# Test 10d: Correct output structure
test_that("reduce_data returns correct structure and format", {
  
  out <- 
    reduce_data(
      data_smooth, 
      check_taxa = FALSE, 
      check_levels = FALSE)
  
check_reduce_data_structure(out)
})





