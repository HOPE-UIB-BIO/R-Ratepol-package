fc_standardise_community_data <- function (data_source, N_individuals = 150, Debug = F)
{
  
  if(Debug == T){
    cat(paste("Data standardization started", Sys.time()), fill = T)}
  
  Samples <-  1:nrow(data_source$Community) 
  Samples <-  Samples[!is.na(data_source$Age$sample.id)]
  
  for(i in Samples){  # for each row(sample)
    
    select_row <-  data_source$Community[i, ] # selected row
    
    n1 <-  1:ncol(data_source$Community) #number for each species name in community data
    ab1 <-  as.vector(select_row) #frequencies of species in each sample
    
    vec1 <-  NULL    #a vector for the species pool
    
    # create a vector with species numbers replicated X times, where X is number of individuals
    for(j in 1:length(n1))  {
      
      #1: number of species 
      v1 <-  rep(n1[j], ab1[j]) #repeat species names 
      vec1 <-  c(vec1, v1) #a vector that repeat the species for the occurrences
      }
    
    # sample species X time (150 times)
    rsample <-  sample(vec1, size = N_individuals, replace = FALSE)
    
    # replace all values in community data by 0
    data_source$Community[i, ]<-  rep(0,length(select_row))
    
    # replace individuals by new randomised values
    data_source$Community[i,as.numeric(names(table(rsample)))] <- as.numeric(table(rsample))
    
  }
  
  if(Debug == T){cat(paste("Data standardization finished", Sys.time()), fill = T)}
  
  return (data_source) 
}