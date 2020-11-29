fc_extract_data <-  function (data_community_extract,
                              data_age_extract,
                              age_uncertainty = F,
                              Debug = F)
{
  if (Debug==T) {
    cat("",fill = T)
    cat(paste("Data extraction started",Sys.time()),fill = T)
  }
  
  # extract both important tables a) age data, b) community data
  age <- data_age_extract
  dat_community <- data_community_extract
  
  if (all(age_uncertainty == F)){
    age.un <- data.frame(t(matrix(rep(age$age,1000), ncol = 1000)))
  } else {
    if( is.matrix(age_uncertainty)){
      age.un <- data.frame(age_uncertainty)
    } else {
      if(is.data.frame(age_uncertainty)){
        age.un <- age_uncertainty
      } else {
        stop("Age uncertainty has to be either FALSE or matrix")
      }
    }
  }
  
  names(age.un) <- age$sample.id
  
  if(names(dat_community) %in% "sample.id"){
    row.names(dat_community) <-  dat_community$sample.id
  } else {
    row.names(dat_community) <-  age$sample.id
  }
  
  # remove the sample ID
  if (is.numeric(unlist(dat_community[,1]))==F){
    dat_community <- dat_community[,-1]
  }
  
  # create a new variable that would be used all latter analyses
  # Newage is a value of interpolated time (time which start with 0)
  age$newage <-age$age
  
  # dim.val are values of the size of the dataset
  dim.val <- vector(mode = "integer", length = 3)
  names(dim.val) <- c("N Species","N samples community","N samples Age")
  
  # create list  class of 4 variables Ppllen, age, age.un, dim.val
  dat.merge <- list(Community=dat_community, Age=age, Age.un = age.un, Dim.val = dim.val)
  class(dat.merge) <- "RRatepolList"
  
  # perform check = round number of species and samples and exclude "empty" ones
  dat.merge <- fc_check_data(dat.merge, proportion = F, Debug = Debug)
  
  suppressWarnings(row.names(dat.merge$Age) <- dat.merge$Age$sample.id)
  suppressWarnings(row.names(dat.merge$Community) <- dat.merge$Age$sample.id)
  
  if (dat.merge$Dim.val[3]!=dat.merge$Dim.val[2]) # check number of rows
    stop("Community and Age data have different number of samples")
  
  if (dat.merge$Dim.val[3]!=ncol(dat.merge$Age.un))
    stop("Age uncertainty data does not have appropriete number of samples")
  
  # check if are sample iD saved as characters and tranform if necesery
  if(is.character(dat.merge$Age$sample.id)==F) {dat.merge$Age$sample.id <- as.character(dat.merge$Age$sample.id)}
  if(is.character(row.names(dat.merge$Community))==F) {row.names(dat.merge$Community) <- as.character(row.names(dat.merge$Community))}
  
  if (any(row.names(dat.merge$Community)!=dat.merge$Age$sample.id)) # check if samples have same name
    stop("Samples code for community data and age data have different names")
  
  if(is.unsorted(dat.merge$Age$age)==T) # check if is age of samples in order
    stop("Age data is not sorted")
  
  if(dat.merge$Age$age[1]>dat.merge$Age$age[dat.merge$Dim.val[3]]){
    if (Debug==T)
    {
      cat("Age data was in decreasing format, changed accordingly",fill = T)
      dat.merge$Age <- dat.merge$Age[order(dat.merge$Age$age),]
    }
    
  }
  if(dat.merge$Age$age[1]<dat.merge$Age$age[dat.merge$Dim.val[3]]){
    if (Debug==T)
    {
      cat("Age data is in increasing format",fill = T)
    }
    
  }
  
  if (Debug==T)
  {
    cat("",fill = T)
    cat(paste("Data extraction completed",Sys.time()),fill = T)
    cat("",fill = T)
  }
  
  return(dat.merge)
}
