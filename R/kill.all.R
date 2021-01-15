kill.all <- function(data_source_check.kill, Species = Species, Samples = Samples){
  
  #' @Description: check if there is a sample that do not have a individuals data and delete it
  #     and check if there are any specie without individuals record and delete them

  
  if(any(rowSums(data_source_check.kill@Community, na.rm = TRUE) == 0) & Samples == TRUE){ # if there are some samples without individuals
    
    data_source_check.kill@Age <-  
      data_source_check.kill@Age[rowSums(data_source_check.kill@Community, na.rm = TRUE) > 0, ]
    
    data_source_check.kill@Age.un <- 
      data_source_check.kill@Age.un[ ,rowSums(data_source_check.kill@Community, na.rm = TRUE) > 0]
    
    data_source_check.kill@Community <- 
      data_source_check.kill@Community[rowSums(data_source_check.kill@Community, na.rm = TRUE) > 0, ]
  }
  
  if(any(colSums(data_source_check.kill@Community, na.rm = TRUE) == 0) & Species == TRUE){ # if there are some species without individuals
    
    data_source_check.kill@Community<- data_source_check.kill@Community[ ,colSums(data_source_check.kill@Community, na.rm = TRUE) > 0]
  }
  
  # count the species and samples
  data_source_check.kill@Dim.val[1] <-  ncol(data_source_check.kill@Community)
  data_source_check.kill@Dim.val[2] <-  nrow(data_source_check.kill@Community)
  data_source_check.kill@Dim.val[3] <-  nrow(data_source_check.kill@Age)
  
  return(data_source_check.kill)
}
