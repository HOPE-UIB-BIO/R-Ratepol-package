fc_standardise_pollen_data <- function (data.source, fc_standar.S.value, Debug=F)
{
  # data.source = imput data
  # fc_standar.S.value = values to standardise number of pollen grains
  
  if(Debug==T){cat(paste("Data standardization started",Sys.time()), fill=T)}
  
  # pollen randomization
  Samples <- 1:nrow(data.source$Pollen) 
  Samples<- Samples[!is.na(data.source$Age$sample.id)]
  
  for(i in Samples)  # for each row(sample)
  {
    select.row <- data.source$Pollen[i,] # selected row
    
    n1 <- 1:ncol(data.source$Pollen) #number for each species name in pollen data
    ab1 <- as.vector(select.row) #frequencies of species or pollen in each sample
    
    vec1 <- NULL    #a vector for the species or pollen pool
    
    # create a vector with species numbers replicatet X times, where X is number of pollen grains
    for(j in 1:length(n1))  
      { #1: number of species or pollen types
      v1 <- rep(n1[j], ab1[j]) #repeat species names 
      vec1 <- c(vec1, v1) #a vector that repeat the species or pollen types for the occurrences
      }
    
    # sample species X time (150 times)
    rsample <- sample(vec1, size = fc_standar.S.value, replace = FALSE)
    
    # replace all values in pollen data by 0
    data.source$Pollen[i,]<- rep(0,length(select.row))
    
    # replace pollen by new randomised values
    data.source$Pollen[i,as.numeric(names(table(rsample)))] <- as.numeric(table(rsample))
    
  }
  
  if(Debug==T){cat(paste("Data standardization finished",Sys.time()),fill=T)}
  
  return (data.source) 
}