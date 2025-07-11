
# extract example data for testing
example_community <- RRatepol::example_data$pollen_data[[1]]
example_age <- RRatepol::example_data$sample_age[[1]]
example_uncertainty <- RRatepol::example_data$age_uncertainty[[1]]


# 1. Test output structure:
test_that("extract_data returns expected structure", {
  result <- extract_data(
    data_community_extract = example_community,
    data_age_extract = example_age,
    age_uncertainty = NULL,
    verbose = FALSE
  )
  
  expect_type(result, "list")
  expect_named(result, c("community", "age", "age_un"))
  expect_s3_class(result$community, "data.frame")
  expect_s3_class(result$age, "data.frame")
  expect_null(result$age_un)
})

# 2. Test for correct output when community data has NAs:
test_that("missing values in community are replaced with 0", {
  example_community_NA <- example_community # copy of community
  example_community_NA[5,3] <- NA # introduce NAs
  
  result <- extract_data(example_community_NA, example_age, verbose = FALSE)
  expect_true(all(!is.na(result$community)))
  expect_true(all(result$community["B", "sp2"] == 0))
})

# 3. Test if "sample.id" column works too 
# if not provided as "sample_id":
test_that("sample.id column is automatically renamed", {
  comm <- example_community
  names(comm)[1] <- "sample.id"
  age <- example_age
  names(age)[1] <- "sample.id"
  
  result <- extract_data(comm, age, verbose = FALSE)

  expect_equal(rownames(result$community), example_community$sample_id)
})


# 4. Test error messages if sample_ids are not the same in age and community
test_that("throws error if sample_ids do not match", {
  broken_age <- example_age
  broken_age$sample_id[3] <- "D"
  expect_error(
    extract_data(example_community, broken_age),
    "Variable 'sample_id' must have same values in
    'data_age' and 'data_community'"
  )
})

# 5. Test if age_uncertainty has sample_ids as rownames
test_that("age_uncertainty is handled correctly", {
  result <- extract_data(
    data_community_extract = example_community,
    data_age_extract = example_age,
    age_uncertainty = example_uncertainty,
    verbose = FALSE
  )
  
  expect_s3_class(result$age_un, "data.frame")
  expect_named(result$age_un, example_age$sample_id)
})
