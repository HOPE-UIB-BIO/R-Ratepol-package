# ==================================================================== #
#                        TESTS FOR reduce_data()                      #
# ==================================================================== #
#
# This file contains comprehensive tests for the reduce_data() function
# which filters taxa and/or levels based on zero-sum criteria.
#
# Test structure:
# 1. Input Validation Tests
# 2. No Filtering Tests (both FALSE)
# 3. Taxa-Only Filtering Tests (check_taxa=TRUE, check_levels=FALSE)
# 4. Levels-Only Filtering Tests (check_taxa=FALSE, check_levels=TRUE)
# 5. Both Filtering Tests (check_taxa=TRUE, check_levels=TRUE)
#
# ==================================================================== #

# --------------------------------------------------- #
# 1. INPUT VALIDATION TESTS - ERRORS
# --------------------------------------------------- #

# 1.1 data_source_reduce validation
test_that("reduce_data rejects NULL data_source_reduce", {
  expect_error(
    reduce_data(data_source_reduce = NULL),
    "data_source_reduce.*list"
  )
})

test_that("reduce_data rejects string data_source_reduce", {
  expect_error(
    reduce_data(data_source_reduce = "invalid"),
    "data_source_reduce.*list"
  )
})

test_that("reduce_data rejects numeric data_source_reduce", {
  expect_error(
    reduce_data(data_source_reduce = 123),
    "data_source_reduce.*list"
  )
})

test_that("reduce_data rejects data.frame data_source_reduce", {
  expect_error(
    reduce_data(data_source_reduce = data.frame(x = 1)),
    "data_source_reduce.*list"
  )
})

test_that("reduce_data rejects empty list data_source_reduce", {
  expect_error(
    reduce_data(data_source_reduce = list())
  )
})

# 1.2 check_taxa validation
test_that("reduce_data rejects invalid check_taxa argument", {
  raw_data <- 
  extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  expect_error(
    reduce_data(data_source_reduce = raw_data, check_taxa = "TRUE"),
    "check_taxa.*logical"
  )
})

test_that("reduce_data rejects invalid check_taxa argument", {
  raw_data <- 
  extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  expect_error(
    reduce_data(data_source_reduce = raw_data, check_taxa = 123),
    "check_taxa.*logical"
  )
})

test_that("reduce_data rejects invalid check_taxa argument", {
  raw_data <- 
  extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  expect_error(
    reduce_data(data_source_reduce = raw_data, check_taxa = NULL),
    "check_taxa.*logical"
  )
})

# 1.3 check_levels validation
test_that("reduce_data rejects invalid check_levels argument", {
  raw_data <- 
  extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  expect_error(
    reduce_data(data_source_reduce = raw_data, check_levels = "TRUE"),
    "check_levels.*logical"
  )
})

test_that("reduce_data rejects invalid check_levels argument", {
  raw_data <- 
  extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  expect_error(
    reduce_data(data_source_reduce = raw_data, check_levels = 123),
    "check_levels.*logical"
  )
})

test_that("reduce_data rejects invalid check_levels argument", {
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
test_that("reduce_data rejects missing community component", {
  raw_data <- 
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  incomplete_data <- list(age = raw_data$age, age_un = raw_data$age_un)
  expect_error(
    reduce_data(data_source_reduce = incomplete_data)
  )
})

test_that("reduce_data rejects missing age component", {
  raw_data <- 
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  incomplete_data <- list(community = raw_data$community, age_un = raw_data$age_un)
  expect_error(
    reduce_data(data_source_reduce = incomplete_data)
  )
})

test_that("reduce_data rejects matrix community", {
  raw_data <- 
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  raw_data$community <- as.matrix(raw_data$community)
  expect_error(
    reduce_data(data_source_reduce = raw_data)
  )
})

test_that("reduce_data rejects community without rownames", {
  raw_data <- 
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  rownames(raw_data$community) <- NULL
  expect_error(
    reduce_data(data_source_reduce = raw_data)
  )
})

test_that("reduce_data rejects community with non-numeric columns", {
  raw_data <- 
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  raw_data$community$text_col <- "NA"
  expect_error(
    reduce_data(data_source_reduce = raw_data),
    "'x' must be numeric"
  )
})

test_that("reduce_data rejects age without rownames", {
  raw_data <- 
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  rownames(raw_data$age) <- NULL
  expect_error(
    reduce_data(data_source_reduce = raw_data)
  )
})


# --------------------------------------------------- #
# 2. NO FILTERING TESTS (check_taxa=FALSE, check_levels=FALSE)
# --------------------------------------------------- #

test_that("reduce_data with no filtering returns identical data", {
  raw_data <- 
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  
  # one zero column in community
  raw_data$community[,1] <- 0  
  # one zero row in community
  raw_data$community[1,] <- 0 

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
test_that("taxa filtering drops taxa with zero column sums", {
  raw_data <- 
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  
  # one zero column in community
  raw_data$community[,1] <- 0  

  result <- 
    reduce_data(
      data_source_reduce = raw_data,
      check_taxa = TRUE,
      check_levels = FALSE
    )

  expect_true(all(colSums(result$community, na.rm = TRUE) > 0))
  expect_equal(ncol(result$community), ncol(raw_data$community) - 1)
})


test_that("check_levels = FALSE preserves rownames / sample ids", {
  raw_data <- 
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  
  # one zero column in community
  raw_data$community[,1] <- 0  
  # one zero row in community
  raw_data$community[1,] <- 0 

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

test_that("taxa filtering removes all-NA taxa", {
  raw_data <- 
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  raw_data$community$na_taxon <- NA

  result <- reduce_data(
    data_source_reduce = raw_data,
    check_taxa = TRUE,
    check_levels = FALSE
  )
  
  expect_false("na_taxon" %in% colnames(result$community))
})

test_that("taxa filtering with all-zero community data returns empty community result", {
  raw_data <- 
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  raw_data$community[,] <- 0

  result <- 
  reduce_data(
    data_source_reduce = raw_data,
    check_taxa = TRUE,
    check_levels = FALSE
  )
  
  expect_equal(ncol(result$community), 0)
})

# --------------------------------------------------- #
# 4. LEVELS-ONLY FILTERING TESTS (check_taxa=FALSE, check_levels=TRUE)
# --------------------------------------------------- #

test_that("levels filtering removes zero-sum levels from community", {
  raw_data <- 
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  
  # one zero column in community
  raw_data$community[,1] <- 0  
  # one zero row in community
  raw_data$community[1,] <- 0

  first_level <- rownames(raw_data$community)[1]
  raw_data$community[1,] <- 0

  result <- reduce_data(
    data_source_reduce = raw_data,
    check_taxa = FALSE,
    check_levels = TRUE
  )
  
  expect_true(all(rowSums(result$community, na.rm = TRUE) > 0)) 
  expect_false(first_level %in% rownames(result$community))
  expect_false(first_level %in% rownames(result$age))
  expect_false(first_level %in% colnames(result$age_un))
})

test_that("levels filtering preserves community colnames", {
  raw_data <- 
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  
  # one zero column in community
  raw_data$community[,1] <- 0  
  # one zero row in community
  raw_data$community[1,] <- 0 
  
  result <- reduce_data(
    data_source_reduce = raw_data,
    check_taxa = FALSE,
    check_levels = TRUE
  )

  expect_identical(colnames(raw_data$community), colnames(result$community))
})

test_that("levels filtering maintains matching rownames between community and age", {
  raw_data <- 
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  
  # one zero column in community
  raw_data$community[,1] <- 0  
  # one zero row in community
  raw_data$community[1,] <- 0 
  
  result <- reduce_data(
    data_source_reduce = data_smooth,
    check_taxa = FALSE,
    check_levels = TRUE
  )
  
  expect_identical(rownames(result$community), rownames(result$age))
  expect_identical(rownames(result$community), colnames(result$age_un))
})



test_that("levels filtering with all zero levels returns empty result", {
  raw_data <- 
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  raw_data$community[,] <- 0

  result <- reduce_data(
    data_source_reduce = raw_data,
    check_taxa = FALSE,
    check_levels = TRUE
  )
  
  expect_equal(nrow(result$community), 0)
  expect_equal(nrow(result$age), 0)
  expect_equal(ncol(result$age_un), 0)
})


test_that("levels filtering handles NULL age_un", {
  raw_data <- 
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      verbose = FALSE
    )

  raw_data$community[1,] <- 0
  
  expect_no_error(
    reduce_data(
      data_source_reduce = raw_data,
      check_taxa = FALSE,
      check_levels = TRUE
    )
  )
})

# --------------------------------------------------- #
# 5. BOTH FILTERING TESTS (check_taxa=TRUE, check_levels=TRUE)
# --------------------------------------------------- #

test_that("both filtering removes zero-sum taxa", {
  raw_data <- 
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  
  # one zero column in community
  raw_data$community[,1] <- 0  
  # one zero row in community
  raw_data$community[1,] <- 0 
  
  result <- reduce_data(
    data_source_reduce = raw_data,
    check_taxa = TRUE,
    check_levels = TRUE
  )
  
  expect_false("zero_taxon" %in% colnames(result$community))
})

test_that("both filtering removes zero-sum levels", {
  raw_data <- 
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )
  
  first_level <- rownames(raw_data$community)[1]
  raw_data$community[1,] <- 0

  result <- reduce_data(
    data_source_reduce = raw_data,
    check_taxa = TRUE,
    check_levels = TRUE
  )
  
  expect_false(first_level %in% rownames(result$community))
})

test_that("both filtering removes taxa with zero observations", {
  raw_data <- 
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  raw_data$community$zero_taxon <- 0

  result <- reduce_data(
    data_source_reduce = raw_data,
    check_taxa = TRUE,
    check_levels = TRUE
  )
  
  expect_true(all(colSums(result$community, na.rm = TRUE) > 0))
})

test_that("both filtering removes samples with zero taxa", {
  raw_data <- 
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  raw_data$community[1,] <- 0

  result <- reduce_data(
    data_source_reduce = raw_data,
    check_taxa = TRUE,
    check_levels = TRUE
  )
  
  expect_true(all(rowSums(result$community, na.rm = TRUE) > 0))
})

test_that("both filtering maintains matching rownames between community and age", {
  raw_data <- 
    extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    )

  raw_data$community$zero_taxon <- 0
  raw_data$community[1,] <- 0

  result <- reduce_data(
    data_source_reduce = raw_data,
    check_taxa = TRUE,
    check_levels = TRUE
  )
  
  expect_identical(rownames(result$community), rownames(result$age))
  expect_identical(colnames(result$age_un), rownames(result$community))
})

















# WIP below


test_that("both filtering maintains matching colnames between age_un and community", {
  data_smooth <- smooth_community_data(
    data_source_smooth = extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    ),
    smooth_method = "m.avg",
    verbose = FALSE
  )
  
  data_smooth$community$zero_taxon <- 0
  data_smooth$community[1,] <- 0
  
  result <- reduce_data(
    data_source_reduce = data_smooth,
    check_taxa = TRUE,
    check_levels = TRUE
  )
  
  expect_equal(colnames(result$age_un), rownames(result$community))
})

test_that("both filtering with all zero data returns empty community", {
  data_smooth <- smooth_community_data(
    data_source_smooth = extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    ),
    smooth_method = "m.avg",
    verbose = FALSE
  )
  
  data_smooth$community[,] <- 0
  
  result <- reduce_data(
    data_source_reduce = data_smooth,
    check_taxa = TRUE,
    check_levels = TRUE
  )
  
  expect_equal(ncol(result$community), 0)
  expect_equal(nrow(result$community), 0)
})

test_that("both filtering with all zero data returns empty age", {
  data_smooth <- smooth_community_data(
    data_source_smooth = extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    ),
    smooth_method = "m.avg",
    verbose = FALSE
  )
  
  data_smooth$community[,] <- 0
  
  result <- reduce_data(
    data_source_reduce = data_smooth,
    check_taxa = TRUE,
    check_levels = TRUE
  )
  
  expect_equal(nrow(result$age), 0)
})

test_that("both filtering with all zero data handles age_un correctly", {
  data_smooth <- smooth_community_data(
    data_source_smooth = extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      age_uncertainty = RRatepol::example_data$age_uncertainty[[1]],
      verbose = FALSE
    ),
    smooth_method = "m.avg",
    verbose = FALSE
  )
  
  data_smooth$community[,] <- 0
  
  result <- reduce_data(
    data_source_reduce = data_smooth,
    check_taxa = TRUE,
    check_levels = TRUE
  )
  
  expect_true(is.null(result$age_un) || ncol(result$age_un) == 0)
})

test_that("both filtering handles NULL age_un", {
  data_smooth <- smooth_community_data(
    data_source_smooth = extract_data(
      data_community_extract = RRatepol::example_data$pollen_data[[1]],
      data_age_extract = RRatepol::example_data$sample_age[[1]],
      verbose = FALSE
    ),
    smooth_method = "m.avg",
    verbose = FALSE
  )
  
  data_smooth$age_un <- NULL
  data_smooth$community$zero_taxon <- 0
  data_smooth$community[1,] <- 0
  
  expect_no_error(
    reduce_data(
      data_source_reduce = data_smooth,
      check_taxa = TRUE,
      check_levels = TRUE
    )
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
  
  expect_error(
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

test_that("reduce_data validates community must have sampleID/levels/rownames", {
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

test_that("reduce_data validates age must have sampleID/levels/rownames", {
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

test_that("reduce_data validates community and age have matching sampleID/levels/rownames", {
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

test_that("reduce_data validates age_un column names match community sampleID/levels/rownames", {
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
