fc_check_data <- function (data_source_check, proportion = F, Species = T, Samples = T, Debug = F)
{
  # check if there is a sample that do not have a individuals data and delete it
  # & 
  # check if there are any specie without individuals record and delete them
  
  kill.all <- function(data_source_check.kill, Species = Species, Samples = Samples){
    if(any(rowSums(data_source_check.kill$Community, na.rm = T) == 0) & Samples == T){ # if there are some samples without individuals
      
      data_source_check.kill$Age <-  
        data_source_check.kill$Age[rowSums(data_source_check.kill$Community, na.rm = T) > 0, ]
      
      data_source_check.kill$Age.un <- 
        data_source_check.kill$Age.un[ ,rowSums(data_source_check.kill$Community, na.rm = T) > 0]
      
      data_source_check.kill$Community <- 
        data_source_check.kill$Community[rowSums(data_source_check.kill$Community, na.rm = T) > 0, ]
    }
    
    if(any(colSums(data_source_check.kill$Community, na.rm=T)==0) & Species == T){ # if there are some species without individuals
      
      data_source_check.kill$Community<- data_source_check.kill$Community[ ,colSums(data_source_check.kill$Community, na.rm = T) > 0]
    }
    
    # count the species and samples
    data_source_check.kill$Dim.val[1] <-  ncol(data_source_check.kill$Community)
    data_source_check.kill$Dim.val[2] <-  nrow(data_source_check.kill$Community)
    data_source_check.kill$Dim.val[3] <-  nrow(data_source_check.kill$Age)
    
    return(data_source_check.kill)
  }
  
  data_source_check <-  
    kill.all(
      data_source_check, 
      Species = Species, 
      Samples = Samples) 
  
  if(Debug==T){
    
    cat("", fill=T)
    cat(paste("Community data have", data_source_check$Dim.val[1], "species with records and",
              data_source_check$Dim.val[2], "samples. Age data have", data_source_check$Dim.val[3], "samples"),
        fill=T)
    
    cat("", fill=T)
    cat(paste("Age data has values of min", min(data_source_check$Age$age),
              ", max",max(data_source_check$Age$age),
              ",mean", round(mean(data_source_check$Age$age),2),
              ",and median", round(median(data_source_check$Age$age),2)),
        fill=T)
    
    cat("", fill=T)
  }
  
  # check if all values is new age are in positive values and interpolate if necesery
  if(any(data_source_check$Age$newage<0, na.rm = T))
  {
    data_source_check$Age$newage <- data_source_check$Age$newage + min(data_source_check$Age$newage)*(-1)
  }
  
  if (proportion == T)
  {
    if (Debug==T){ cat("Community data values are being converted to proportions", fill=T)}
    
    # convert the values community data to proportion of sum of each sample
    p.counts.row.sums <- apply(data_source_check$Community, 1, sum)
    data_source_check$Community <- as.data.frame(lapply(data_source_check$Community, function(x) x/p.counts.row.sums))
    data_source_check <- kill.all(data_source_check, Species = Species, Samples = Samples)
  }
  
  return(data_source_check)
}