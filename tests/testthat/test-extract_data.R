
# Extract example data for testing
example_community <- 
  RRatepol::example_data$pollen_data[[1]]

example_age <- 
  RRatepol::example_data$sample_age[[1]]

example_uncertainty <- 
  RRatepol::example_data$age_uncertainty[[1]]



# 1. Test output structure:
test_that("returns expected output structure", {
  # run extract_data()
  result <-
    extract_data(
      example_community,
      example_age,
      example_uncertainty,
      verbose = FALSE
    )

  expect_type(result, "list")
  expect_named(result, c("community", "age", "age_un"))
  
  expect_s3_class(result$community, "data.frame")
  expect_s3_class(result$age, "data.frame")
  expect_s3_class(result$age_un, "data.frame")
  
  expect_equal(row.names(result$community), example_community$sample_id)
  expect_equal(row.names(result$age), example_age$sample_id)
  expect_named(result$age_un, example_age$sample_id)
  
})


# 2. Test NA behaviour and related output structure:
test_that("returns expected NA behaviour", {
  # introduce NAs:
  example_community_NA <-
    example_community # copy of community

  example_community_NA[5:10, 3] <-
    NA

  example_age_NA <-
    example_age

  example_age_NA[5:10, 3] <-
    NA

  example_uncertainty_NA <-
    example_uncertainty

  example_uncertainty_NA[5:10, 5] <-
    NA

  # run extract_data()
  result <-
    extract_data(
      example_community_NA,
      example_age_NA,
      example_uncertainty_NA,
      verbose = TRUE
    )

  expect_true(all(!is.na(result$community)))
  expect_true(all(result$community[5:10, 2] == 0))

  expect_type(result, "list")
  expect_named(result, c("community", "age", "age_un"))
  
  expect_s3_class(result$community, "data.frame")
  expect_s3_class(result$age, "data.frame")
  expect_s3_class(result$age_un, "data.frame")
  
  expect_equal(row.names(result$community), example_community_NA$sample_id)
  expect_equal(row.names(result$age), example_age_NA$sample_id)
  expect_named(result$age_un, example_age_NA$sample_id)
  
})



# 3. Test if "sample.id" column works too 
# if not provided as "sample_id":
test_that("sample.id column is automatically renamed", {
  # rename sample_id to sample.id
  comm <- example_community
  names(comm)[1] <- "sample.id"
  age <- example_age
  names(age)[1] <- "sample.id"
  
  result <- extract_data(comm, age, verbose = FALSE)

  expect_equal(rownames(result$community), example_community$sample_id)
  
})



# 4. Test error messages if sample_ids are not the same in age and community
test_that("throws error if sample_ids do not match", {
  # test wrong sample IDs
  broken_age <- example_age
  broken_age$sample_id[3] <- "D"
  expect_error(
    extract_data(example_community, broken_age),
    "Variable 'sample_id' must have same values in
    'data_age' and 'data_community'"
  )
  
})



# 5. Test that age_uncertainty has sample_ids as row names
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

