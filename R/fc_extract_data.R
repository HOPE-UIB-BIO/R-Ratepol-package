fc_extract_data <-  function (data_community_extract,
                              data_age_extract,
                              age_uncertainty = FALSE,
                              Debug = FALSE)
{
  if (Debug == T) {
    cat("", fill = T)
    cat(paste("Data extraction started", Sys.time()), fill = TRUE)
  }
  
  # 1. Initial tests -----

  n_samples_com <- nrow(data_community_extract)
  n_samples_age <- nrow(data_age_extract)
  
  # 1.1 Data types -----
  
  assertthat::assert_that(
    'data.frame' %in% class(data_community_extract),
    msg = "Object 'data_source_community' must be a 'data.frame'")
  
  assertthat::assert_that(
    'data.frame' %in% class(data_age_extract),
    msg = "Object 'data_source_age' must be a 'data.frame'")
  
  # 1.2. Sample id -----
  
  assertthat::assert_that(
    'sample.id' %in% names(data_community_extract) | 
      'sample_id' %in% names(data_community_extract)  ,
    msg = "Variable 'sample.id' or 'sample_id' must be present in the 
    'data_source_community'")
  
  assertthat::assert_that(
    'sample.id' %in% names(data_age_extract) |
      'sample_id' %in% names(data_age_extract) ,
    msg = "Variable 'sample.id' or 'sample_id' must be present in the 
    'data_source_age'")

  # community
  if('sample_id' %in% names(data_community_extract)){
    data_community_extract <- 
      data_community_extract %>% 
      dplyr::rename(sample.id  = sample_id)
  }
  
  assertthat::assert_that(
    'character' %in% class(data_community_extract$sample.id),
    msg = "Variable 'sample.id' or 'sample_id' in 'data_source_community' must
    be a 'character'")
  
  #age
  if('sample_id' %in% names(data_age_extract)){
    data_age_extract <- 
      data_age_extract %>% 
      dplyr::rename(sample.id  = sample_id)
  }
  
  assertthat::assert_that(
    'character' %in% class(data_age_extract$sample.id),
    msg = "Variable 'sample.id' in or 'sample_id' 'data_source_age' must be
    a 'character'")
  
  assertthat::assert_that(
    all(data_community_extract$sample.id == data_age_extract$sample.id),
    msg = "Variable 'sample.id'/'sample_id' must have same values in 
    'data_source_age' and 'data_source_community'")
  
  # 1.3. Age test -----
  
  assertthat::assert_that(
    'age' %in% names(data_age_extract),
    msg = "Variable 'age' must be present in 'data_source_age'")
  
  assertthat::assert_that(
    'numeric' %in% class(data_age_extract$age),
    msg = "Variable 'age' in 'data_source_age' must be a 'numeric'")
  
  assertthat::assert_that(
    is.unsorted(data_age_extract$age) == FALSE,
    msg = "Variable 'age' in 'data_source_age' must be ordered by age")
  
  # order of the age
  if(data_age_extract$age[1] > data_age_extract$age[n_samples_age]){
    
    data_age_extract$age <- data_age_extract[order(data_age_extract$age),]
    
    
    if (Debug == T){
      cat("Variable 'age' in 'data_source_age' was stored in decreesing format,
          changed accordingly for analyses", fill = TRUE)
    }
    
  } else {
    if (Debug == T){
      cat("Variable 'age' in 'data_source_age' is stored in increesing format", fill = TRUE)
    }
    
  }
  
  # 1.4. Size test ----- 
  
  assertthat::assert_that(
    n_samples_com == n_samples_age,
    msg = "Object 'data_source_community' and 'data_source_age' 
    must have the same number of levels")

 
  # 2. Processes data ----- 
  
  # 2.1 Extract data  ----- 
  
  # extract both important tables if stored as tibbles
  age <- as.data.frame(data_age_extract)
  dat_community <- as.data.frame(data_community_extract)
  
  # exclude any other varibales than age and sample id
  age <-  age[ ,names(age) %in% c("sample.id", "age")]
  
  # 2.2 Age uncertainty  ----- 
  
  # if age_uncertainty is used
  if (!all(age_uncertainty == FALSE)){
    assertthat::assert_that(
      'matrix' %in% class(age_uncertainty),
      msg = "Object 'age_uncertainty' must be a 'matrix'")
    
    n_samples_un <- ncol(age_uncertainty)
    
    assertthat::assert_that(
      n_samples_age == n_samples_un,
      msg = "Object 'data_source_age' and 'age_uncertainty' must have 
      the same number of levels. 'age_uncertainty' must have samples stored as 
      columns")
    
    # save as dataframe
    age.un <- data.frame(age_uncertainty)
    
  } else {
    age.un <- data.frame(t(matrix(rep(age$age, 10), ncol = 10)))
  }
  
  # 2.3 Row.names  ----- 
  
  # add row.names to commity, age, and uncertainty data
  row.names(dat_community) <-  dat_community$sample.id
  row.names(age) <-  age$sample.id
  names(age.un) <- age$sample.id
  
  # remove the sample.id
  dat_community <- dat_community[ , !names(dat_community) %in% "sample.id"]
  
  # create a new variable that would be used all latter analyses
  # Newage is a value of interpolated time (time which start with 0)
  age$newage <-age$age
  
  # 2.4 Summary  -----   
  
  # dim_val are values of the size of the dataset
  dim_val <- vector(mode = "integer", length = 3)
  names(dim_val) <- c("N Species","N samples community","N samples Age")
  
  # create list  class of 4 variables Ppllen, age, age.un, dim_val
  dat_merge <- 
    RRatepolList(
      Community = dat_community,
      Age = age,
      Age.un = age.un,
      Dim.val = dim_val)
  
  # perform check = round number of species and samples and exclude "empty" ones
  dat_merge <- fc_check_data(dat_merge, proportion = F, Debug = Debug)
  
  if (Debug==T)
  {
    cat("",fill = T)
    cat(paste("Data extraction completed",Sys.time()),fill = T)
    cat("",fill = T)
  }
  
  return(dat_merge)
}
