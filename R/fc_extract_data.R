fc_extract_data <-  function (data_community_extract,
                              data_age_extract,
                              age_uncertainty = FALSE,
                              Debug = FALSE)
{
  if (Debug == T) {
    cat("", fill = T)
    cat(paste("Data extraction started", Sys.time()), fill = TRUE)
  }
  
  # Initial tests ----
  
  # community data
  assertthat::assert_that(
    'data.frame' %in% class(data_community_extract),
    msg = "Object `data_source_community` must be a `data.frame`")
  
  n_samples_com <- nrow(data_community_extract)
  n_species_com <- ncol(data_community_extract)
  
  assertthat::assert_that(
    'sample.id' %in% names(data_community_extract),
    msg = "Variable `sample.id` must be present in `data_source_community`")
  
  assertthat::assert_that(
    'character' %in% class(data_community_extract$sample.id),
    msg = "Variable `sample.id` in `data_source_community` must be a `character`")
  
  
  # age data
  assertthat::assert_that(
    'data.frame' %in% class(data_age_extract),
    msg = "Object `data_source_age` must be a `data.frame`")
  
  n_samples_age <- nrow(data_age_extract)
  
  assertthat::assert_that(
    'sample.id' %in% names(data_age_extract),
    msg = "Variable `sample.id` must be present in `data_source_age`")
  
  assertthat::assert_that(
    'character' %in% class(data_age_extract$sample.id),
    msg = "Variable `sample.id` in `data_source_age` must be a `character`")
  
  assertthat::assert_that(
    'age' %in% names(data_age_extract),
    msg = "Variable `age` must be present in `data_source_age`")
  
  assertthat::assert_that(
    'numeric' %in% class(data_age_extract$age),
    msg = "Variable `age` in `data_source_age` must be a `numeric`")
  
  assertthat::assert_that(
    is.unsorted(data_age_extract$age) == FALSE,
    msg = "Variable `age` in `data_source_age` must be ordered by age")
  
  # order of the age
  if(data_age_extract$age[1] > data_age_extract$age[n_samples_age]){
    
    data_age_extract$age <- data_age_extract[order(data_age_extract$age),]
    
    
    if (Debug == T){
      cat("Variable `age` in `data_source_age` was stored in decreesing format,
          changed accordingly for analyses", fill = TRUE)
    }
    
  } else {
    if (Debug == T){
      cat("Variable `age` in `data_source_age` is stored in increesing format", fill = TRUE)
    }
    
  }
  
  # both
  assertthat::assert_that(
    n_samples_com == n_samples_age,
    msg = "Object `data_source_community` and `data_source_age` 
    must have the same number of levels")
  
  assertthat::assert_that(
    all(data_community_extract$sample.id == data_age_extract$sample.id),
    msg = "Variable `sample.id` must have same values in 
    `data_source_age` and `data_source_community`")
  
  # extract both important tables if stored as tibbles
  age <- as.data.frame(data_age_extract)
  dat_community <- as.data.frame(data_community_extract)
  
  # exclude any other varibales than age and sample id
  age <-  age[ ,names(age) %in% c("sample.id", "age")]
  
  # if age_uncertainty is used
  if (!all(age_uncertainty == F)){
    assertthat::assert_that(
      'matrix' %in% class(age_uncertainty),
      msg = "Object `age_uncertainty` must be a `matrix`")
    
    n_samples_un <- ncol(age_uncertainty)
    
    assertthat::assert_that(
      n_samples_un == n_samples_un,
      msg = "Object `data_source_age` and `age_uncertainty` must have 
      the same number of levels. `age_uncertainty` must have samples stored as 
      columns")
    
    # save as dataframe
    age.un <- data.frame(age_uncertainty)
    
  } else {
    age.un <- data.frame(t(matrix(rep(age$age, 10), ncol = 10)))
  }
  
  # sample.id ----

  # add row.names to commity, age, and uncertainty data
  row.names(dat_community) <-  dat_community$sample.id
  row.names(age) <-  age$sample.id
  names(age.un) <- age$sample.id
  
  # remove the sample.id
  dat_community <- dat_community[ , !names(dat_community) %in% "sample.id"]
  
  # create a new variable that would be used all latter analyses
  # Newage is a value of interpolated time (time which start with 0)
  age$newage <-age$age
  
  # dim.val are values of the size of the dataset
  dim.val <- vector(mode = "integer", length = 3)
  names(dim.val) <- c("N Species","N samples community","N samples Age")
  
  # create list  class of 4 variables Ppllen, age, age.un, dim.val
  dat.merge <- RRatepolList(Community = dat_community, Age=age, Age.un = age.un, Dim.val = dim.val)
  
  # perform check = round number of species and samples and exclude "empty" ones
  dat.merge <- fc_check_data(dat.merge, proportion = F, Debug = Debug)
  
  if (Debug==T)
  {
    cat("",fill = T)
    cat(paste("Data extraction completed",Sys.time()),fill = T)
    cat("",fill = T)
  }
  
  return(dat.merge)
}
