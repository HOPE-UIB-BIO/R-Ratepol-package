fc_calculate_DC <- function (data_source_DC, DC = "chord", Debug = F){

  dat_res <-  vector(
    mode="numeric",
    length = data_source_DC$Dim.val[2]-1)

  #----------------------------------------------------------#
  # Euclidan distance ----- 
  #----------------------------------------------------------#

  if (DC == "euc"){
    
    if (Debug == T){
      cat("Euclidan distance will be used as DC", fill=T)}
    
    for (i in 1:(data_source_DC$Dim.val[2] - 1)){ # for each sample (except the last)
    
      df_work <-  data_source_DC$Community[c(i, i + 1), ] # select only 2 samples (observed + 1 after)

      df_work <-  as.data.frame(df_work[ ,colSums(df_work) > 0]) # get rid of "empty species"

      vector_work <-  vector( # vector for result for each species
        mode="numeric",
        length = ncol(df_work)) 

      for( j in 1:ncol(df_work)){ # for each species
      
        a <-  .subset2(df_work, j)[1]
        b <-  .subset2(df_work, j)[2]
        vector_work[j] <-  (a-b)**2 # calculate the diference
      }

      dat_res[i] <-  sqrt(sum(vector_work)) # save the square root of sum of all dufereces

    }

  }

  #----------------------------------------------------------#
  # Standardised euclidan distace -----
  #----------------------------------------------------------#

  if (DC == "euc.sd"){
    
    if(Debug == T){
      cat("Standardised Euclidan distance will be used as DC", fill=T)}

    # calculation of standard deviation for each species

    # vector for standar deviation for each species
    df_sp_supp <-  vector(
      mode="numeric",
      length = data_source_DC$Dim.val[1])
    
    # calculate the SD for each species
    df_sp_supp <-  apply(data_source_DC$Community, 2 , sd)

    # calculation of the DC
    for (i in 1:(data_source_DC$Dim.val[2] - 1)){ # for each sample (except the last)
      
      df_work <-  data_source_DC$Community[c(i, i + 1), ] # select only 2 samples (observed + 1 after)

      # get rid of "empty species" in data & in sp.std
      df_sp_supp_work <-  df_sp_supp[colSums(df_work) > 0]
      df_work <-  as.data.frame(df_work[ ,colSums(df_work) > 0])

      vector_work <-  vector( # vector for result for each species
        mode="numeric",
        length = ncol(df_work)) 

      for( j in 1:ncol(df_work)){ # for each species
        
        if (df_sp_supp_work[j] != 0) { # check if the standard deviation is not equal zero
          a <-  .subset2(df_work, j)[1]
          b <-  .subset2(df_work, j)[2]
          vector_work[j] <-  ((a - b) / df_sp_supp_work[j])**2 # calculate the diference
        }

      }

      dat_res[i] <-  sqrt(sum(vector_work)) # save the square root of sum of all duferece
    }
  }


  #----------------------------------------------------------#
  # Chord's distance -----
  #----------------------------------------------------------#

  if (DC == "chord"){
    
    if(Debug == T){
      cat("Chord distance will be used as DC", fill = T)}

    for (i in 1:(data_source_DC$Dim.val[2] - 1)){ # for each sample (except the last)
    
      df_work <-  data_source_DC$Community[c(i, i + 1), ] # select only 2 samples (observed + 1 after)

      df_work <-  as.data.frame(df_work[ ,colSums(df_work) > 0]) # get rid of "empty species"

      vector_work <-  vector( # vector for result for each species
        mode="numeric",
        length = ncol(df_work)) 

      for( j in 1:ncol(df_work)){ # for each species
      
        a <-  .subset2(df_work, j)[1]
        b <-  .subset2(df_work, j)[2]
        vector_work[j] <-  ( sqrt(a) - sqrt(b) )**2 # calculate the diference
      }

      dat_res[i] <-  sqrt(sum(vector_work)) # save the square root of sum of all dufereces

    }

  }


  #----------------------------------------------------------#
  # Chi-squared coeficient ----- 
  #----------------------------------------------------------#

  if (DC == "chisq"){
    
    if(Debug == T){
      cat("Chi-squared coeficient will be used as DC", fill = T)}

    for (i in 1:(data_source_DC$Dim.val[2] - 1)){ # for each sample (except the last)
   
      df_work <-  data_source_DC$Community[c(i, i + 1), ] # select only 2 samples (observed + 1 after)

      df_work <-  as.data.frame(df_work[ ,colSums(df_work) > 0]) # get rid of "empty species"

      vector_work <-  vector( # vector for result for each species
        mode="numeric",
        length = ncol(df_work)) 

      for( j in 1:ncol(df_work)){ # for each species 
      
        a <-  .subset2(df_work, j)[1]
        b <-  .subset2(df_work, j)[2]
        vector_work[j] <-  ( (a - b)**2 ) / (a + b) # calculate the diference
      }

      dat_res[i] <-  sqrt(sum(vector_work)) # save the square root of sum of all dufereces

    }

  }

  #----------------------------------------------------------#
  # Gower's disatnce -----
  #----------------------------------------------------------#
  
  if (DC == "gower"){
    
    if (Debug == T){
      cat("Gower's distance will be used as DC", fill = T)}
    
    # use premade function from StatMatch package to calculate correlation between all samples
    corrmat <-  StatMatch::gower.dist(data_source_DC$Community)
    
    # only include calculation between subsequent samples 
    dat_res <-  corrmat[row(corrmat) == col(corrmat) + 1]
  }
  
  #----------------------------------------------------------#
  # Return result -----
  #----------------------------------------------------------#

  return(dat_res)


}
