fc_calculate_DC <- function (data_source_DC, DC = "chord", Debug = FALSE){
  
  
  assertthat::assert_that(
    any(DC == c("euc", "euc.sd", "chord", "chisq", "gower", "bray")),
    msg = "'DC' must be one of the following:
    'euc', 'euc.sd', 'chord', 'chisq', 'gower', 'bray'")
  
  
  # pre-allocate some space
  dat_res <-  
    vector(
      mode = "numeric",
      length = data_source_DC@Dim.val[2] - 1)
  
  #----------------------------------------------------------#
  # Euclidan distance ----- 
  #----------------------------------------------------------#
  
  if(DC == "euc"){
    
    if (Debug == TRUE){
      cat("Euclidan distance will be used as DC", fill = TRUE)
      }
    
    # for each sample (except the last)
    for(i in 1:(data_source_DC@Dim.val[2] - 1)){ 
      
      # select only 2 samples (observed + 1 after)
      df_work <-  data_source_DC@Community[c(i, i + 1), ] 
      
      # get rid of "empty species"
      df_work <-  as.data.frame(df_work[ ,colSums(df_work) > 0]) 
      
      # vector for result for each species
      vector_work <-  
        vector( 
          mode="numeric",
          length = ncol(df_work)) 
      
      # for each species
      for(j in 1:ncol(df_work)){ 
        a <-  .subset2(df_work, j)[1]
        b <-  .subset2(df_work, j)[2]
        
        # calculate the diference
        vector_work[j] <-  (a-b)**2 
      }
      
      # save the square root of sum of all differeces
      dat_res[i] <-  sqrt(sum(vector_work)) 
    }
    
  }
  
  #----------------------------------------------------------#
  # Standardised euclidan distace -----
  #----------------------------------------------------------#
  
  if(DC == "euc.sd"){
    
    if(Debug == TRUE){
      cat("Standardised Euclidan distance will be used as DC", fill = TRUE)
      }
    
    # calculation of standard deviation for each species
    
    # vector for standard deviation for each species
    df_sp_supp <-  
      vector(
      mode = "numeric",
      length = data_source_DC@Dim.val[1])
    
    # calculate the SD for each species
    df_sp_supp <-  apply(data_source_DC@Community, 2 , stats::sd)
    
    # calculation of the DC
    # for each sample (except the last)
    for(i in 1:(data_source_DC@Dim.val[2] - 1)){ 
      
      # select only 2 samples (observed + 1 after)
      df_work <-  data_source_DC@Community[c(i, i + 1), , drop = FALSE] 
      
      # get rid of "empty species" in data & in sp.std
      df_sp_supp_work <-  df_sp_supp[colSums(df_work) > 0]
      df_work <-  as.data.frame(df_work[ ,colSums(df_work) > 0])
      
      # vector for result for each species
      vector_work <-  vector( 
        mode = "numeric",
        length = ncol(df_work)) 
      
      # for each species
      for(j in 1:ncol(df_work)){ 
        
        # check if the standard deviation is not equal zero
        if (df_sp_supp_work[j] != 0) { 
          a <-  .subset2(df_work, j)[1]
          b <-  .subset2(df_work, j)[2]
         
          # calculate the difference
          vector_work[j] <-  ((a - b) / df_sp_supp_work[j])**2 
        }
      }
      
      # save the square root of sum of all differece
      dat_res[i] <-  sqrt(sum(vector_work)) 
    }
  }
  
  
  #----------------------------------------------------------#
  # Chord's distance -----
  #----------------------------------------------------------#
  
  if(DC == "chord"){
    
    if(Debug == TRUE){
      cat("Chord distance will be used as DC", fill = TRUE)
    }
    
    # for each sample (except the last)
    for(i in 1:(data_source_DC@Dim.val[2] - 1)){ 
      
      # select only 2 samples (observed + 1 after)
      df_work <-  data_source_DC@Community[c(i, i + 1), , drop = FALSE]
      
      # get rid of "empty species"
      df_work <-  as.data.frame(df_work[ ,colSums(df_work) > 0]) 
      
      # vector for result for each species
      vector_work <-  
        vector( 
          mode = "numeric",
          length = ncol(df_work)) 
      
      # for each species
      for(j in 1:ncol(df_work)){ 
        a <-  .subset2(df_work, j)[1]
        b <-  .subset2(df_work, j)[2]
        
        # calculate the diference
        vector_work[j] <-  ( sqrt(a) - sqrt(b) )**2 
      }
      
      # save the square root of sum of all differences
      dat_res[i] <-  sqrt(sum(vector_work)) 
    }
  }
  
  
  #----------------------------------------------------------#
  # Chi-squared coeficient ----- 
  #----------------------------------------------------------#
  
  if(DC == "chisq"){
    
    if(Debug == TRUE){
      cat("Chi-squared coeficient will be used as DC", fill = TRUE)}
    
    # for each sample (except the last)
    for(i in 1:(data_source_DC@Dim.val[2] - 1)){ 
      
      # select only 2 samples (observed + 1 after)
      df_work <-  data_source_DC@Community[c(i, i + 1), , drop = FALSE] 
      
      # get rid of "empty species"
      df_work <-  as.data.frame(df_work[ ,colSums(df_work) > 0]) 
      
      # vector for result for each species
      vector_work <-  vector( 
        mode = "numeric",
        length = ncol(df_work)) 
      
      # for each species 
      for(j in 1:ncol(df_work)){ 
        a <-  .subset2(df_work, j)[1]
        b <-  .subset2(df_work, j)[2]
        
        # calculate the diference
        vector_work[j] <-  ( (a - b)**2 ) / (a + b) 
      }
      
      # save the square root of sum of all differeces
      dat_res[i] <-  sqrt(sum(vector_work)) 
    }
  }
  
  
  #----------------------------------------------------------#
  # Gower's distance -----
  #----------------------------------------------------------#
  
  if(DC == "gower"){
    
    if(Debug == T){
      cat("Gower's distance will be used as DC", fill = T)}
    
    # use pre-made function from vegan package to 
    #   calculate correlation between all samples
    corrmat <-  as.matrix(vegan::vegdist(data_source_DC@Community, method = "gower"))
    
    # only include calculation between subsequent samples 
    dat_res <-  corrmat[row(corrmat) == col(corrmat) + 1]
  }
  
  
  #----------------------------------------------------------#
  # Bray–Curtis distance -----
  #----------------------------------------------------------#
  
  if(DC == "bray"){
    
    if(Debug == T){
      cat("Bray–Curtis distance will be used as DC", fill = T)}
    
    # use pre-made function from vegan package to 
    #   calculate correlation between all samples
    corrmat <-  as.matrix(vegan::vegdist(data_source_DC@Community, method = "bray"))
    
    # only include calculation between subsequent samples 
    dat_res <-  corrmat[row(corrmat) == col(corrmat) + 1]
  }
  
  
  #----------------------------------------------------------#
  # Return result -----
  #----------------------------------------------------------#
  
  return(dat_res)
}
