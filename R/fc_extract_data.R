fc_extract_data <-  function (data_pollen_extract,
                              data_age_extract,
                              age_uncertainty = F,
                              Debug = F)
{
  # # data.source.pollen = pollen data with species as collumns and levels as rows, level ID as row names
  # data.source.age = list of 2:
  #                       $ages = datframe
  #                                   $sample.id = unique ID of each level
  #                                   #age = age of level
  #                       $age_position = matrix with number of collumns as number of levels. Each column is one level
  #                                         each row is one age sequence from bchron
  #
  # result of function is list length 3
  # [1] Pollen data
  # [2] Age data with samples ordered by age
  # [3] Number of pollen species and number of samples for pollen & age data
  # [4] Dataframe with all age uncertainties

  if (Debug==T) {
    cat("",fill = T)
    cat(paste("Data extraction started",Sys.time()),fill = T)
  }

  # extract both important tables a) age data, b) pollen data
  age <- data_age_extract$ages
  p.counts <- data_pollen_extract

  if (age_uncertainty == T){
    age.un <- data.frame(data_age_extract$age_position)
  } else {
    age.un <- data.frame(t(matrix(rep(age$age,1000), ncol = 1000)))
  }

  names(age.un) <- age$sample.id

  # remove the sample ID
  if (is.numeric(unlist(p.counts[,1]))==F){
    p.counts <- p.counts[,-1]
  }

  # create a new variable that would be used all latter analyses
  # Newage is a value of interpolated time (time which start with 0)
  age$newage <-age$age

  # dim.val are values of the size of the dataset
  dim.val <- vector(mode = "integer", length = 3)
  names(dim.val) <- c("N Species","N samples pollen","N samples Age")

  # create list  class of 4 variables Ppllen, age, age.un, dim.val
  dat.merge <- list(Pollen=p.counts, Age=age, Age.un = age.un, Dim.val = dim.val)
  class(dat.merge) <- "RRatepolList"

  # perform check = round number of species and samples and exclude "empty" ones
  dat.merge <- fc_check_data(dat.merge, proportion = F, Debug = Debug)

  row.names(dat.merge$Age) <- dat.merge$Age$sample.id
  row.names(dat.merge$Pollen) <- dat.merge$Age$sample.id

  if (dat.merge$Dim.val[3]!=dat.merge$Dim.val[2]) # check number of rows
    stop("Pollen and Age data have different number of samples")

  if (dat.merge$Dim.val[3]!=ncol(dat.merge$Age.un))
    stop("Age uncertainty data does not have appropriete number of samples")

  # check if are sample iD saved as characters and tranform if necesery
  if(is.character(dat.merge$Age$sample.id)==F) {dat.merge$Age$sample.id <- as.character(dat.merge$Age$sample.id)}
  if(is.character(row.names(dat.merge$Pollen))==F) {row.names(dat.merge$Pollen) <- as.character(row.names(dat.merge$Pollen))}

  if (any(row.names(dat.merge$Pollen)!=dat.merge$Age$sample.id)) # check if samples have same name
    stop("Samples code for pollen data and age data have different names")

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
