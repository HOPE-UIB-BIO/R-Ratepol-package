fc_subset_samples <- function(data_subset, bins, WU, bin_selection){
  
  # create empty template with length of number of bins for both age and pollen
  
  data_result_age <-  
    data.frame(
      matrix(
        ncol = ncol(data_subset@Age),
        nrow = nrow(bins)))
  
  names(data_result_age) <-  names(data_subset@Age)
  
  data_result_community <-  
    data.frame(
      matrix(
        ncol = ncol(data_subset@Community),
        nrow = nrow(bins)))
  
  names(data_result_community) <-  names(data_subset@Community)
  
  row.names(data_result_age) <-  bins$name
  row.names(data_result_community) <-  bins$name
  
  # calculate bin size
  bin_size <-  bins$name[2] - bins$name[1]
  
  for(i in 1:nrow(bins)){ # for each bin

    selected_bin <- bins$name[i] #select the bin
    
    # subset age data so it selected all samples which has higher values than the BIN itself but
    # still small then selected bin + calculated BIN size
    subset_w <-  data_subset@Age[data_subset@Age$newage < bins$name[i] + bin_size &
                                  data_subset@Age$newage > bins$name[i], ]
    
    if (nrow(subset_w) > 0){ # If selected subset has at least one sample
  
      # for binning
      if (WU != "levels" & bin_selection == "random"){
        # select one random sample from the bin
        random_row <- sample(1:nrow(subset_w), 1)
        suppressWarnings(data_result_age[i, ] <-  subset_w[random_row, ] )
        data_result_age$age[i] <- selected_bin + (bin_size/2)
        data_result_age$newage[i] <- data_result_age$age[i]
        
        data_result_community[i, ]<-  data_subset@Community[row.names(data_subset@Community) %in% data_result_age$sample.id[i], ]  
      }
      
      # for moving window
      if (WU != "levels" & bin_selection == "first"){
        # select the sample which is the closest to the beggining of the bin
        subset_w$diff <-  abs(subset_w$newage - selected_bin)
        suppressWarnings(data_result_age[i, ] <-  subset_w[subset_w$diff == min(subset_w$diff), ] )
        data_result_age$age[i] <- selected_bin + (bin_size/2)
        data_result_age$newage[i] <- data_result_age$age[i]
        
        data_result_community[i, ]<-  data_subset@Community[row.names(data_subset@Community) %in% data_result_age$sample.id[i], ]  
      }
      
    }
  }
  
  # reduce the empty collumns 
  keep_row <-
    ifelse(is.na(data_result_age$newage),FALSE,TRUE)

  data_result_age <-
    data_result_age[keep_row, ] 
  
  data_result_community <- 
    data_result_community[keep_row, ] 
  
  list.res <-
    RRatepolList(
      Community = data_result_community,
      Age = data_result_age,
      Age.un = data_subset@Age.un,
      Dim.val = data_subset@Dim.val )
  
  list.res <- fc_check_data(list.res, proportion = F, verbose = F)

  
  return(list.res)
}
