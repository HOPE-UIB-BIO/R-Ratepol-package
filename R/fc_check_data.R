fc_check_data <- function (data_source_check, proportion = F, Species = TRUE, Samples = TRUE, Debug = F)
{
  # check if there is a sample that do not have a individuals data and delete it
  # & 
  # check if there are any specie without individuals record and delete them
  
  data_source_check <-  
    RRatepol:::fc_kill_all(
      data_source_check, 
      Species = Species, 
      Samples = Samples) 
  
  if(Debug == TRUE){
    
    cat("", fill = TRUE)
    cat(
      paste(
        "Community data have", data_source_check@Dim.val[1], "species with records and",
        data_source_check@Dim.val[2], "samples. Age data have", data_source_check@Dim.val[3], "samples"),
      fill = TRUE)
    
    cat("", fill = TRUE)
    cat(
      paste(
        "Age data has values of min", min(data_source_check@Age$age),
        ", max",max(data_source_check@Age$age),
        ", mean", round(mean(data_source_check@Age$age), 2),
        ", and median", round(median(data_source_check@Age$age), 2)),
      fill = TRUE)
    
    cat("", fill = TRUE)
  }
  
  # check if all values is new age are in positive values and interpolate if necesery
  if(any(data_source_check@Age$newage<0, na.rm = TRUE))
  {
    data_source_check@Age$newage <- 
      data_source_check@Age$newage + min(data_source_check@Age$newage)*(-1)
  }
  
  if (proportion == TRUE)
  {
    if (Debug == TRUE){
      cat("Community data values are being converted to proportions", fill = TRUE)}
    
    # convert the values community data to proportion of sum of each sample
    p.counts.row.sums <- apply(data_source_check@Community, 1, sum)
    
    data_source_check@Community <- 
      as.data.frame(lapply(data_source_check@Community, function(x) x/p.counts.row.sums))
    
    data_source_check <- 
      Ratepol:::fc_kill_all(
        data_source_check,
        Species = Species,
        Samples = Samples)
  }
  
  return(data_source_check)
}