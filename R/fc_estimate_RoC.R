fc_estimate_RoC <- function(data_source_community,
                            data_source_age,
                            age_uncertainty = FALSE,
                            smooth_method = "none",
                            smooth_N_points = 5,
                            smooth_age_range = 500,
                            smooth_N_max = 9,
                            Working_Units = "levels",
                            bin_size = 500,
                            Number_of_shifts = 1,
                            rand = 1,
                            treads = FALSE,
                            standardise = FALSE,
                            N_individuals = 150,
                            tranform_to_proportions = TRUE,
                            DC = "euc",
                            interest_threshold = FALSE,
                            Debug = FALSE){
  
  # Start of the code
  start.time <- Sys.time()
  if (Debug == TRUE){
    cat(paste("RATEPOL started", start.time), fill=TRUE)
  }
  
  #----------------------------------------------------------# 
  # 1. Data extraction -----
  #----------------------------------------------------------#
  
  # extract data into working format
  # already include data check
  data_extract <- 
    RRatepol:::fc_extract_data(
      data_source_community,
      data_source_age,
      age_uncertainty = age_uncertainty,
      Debug = Debug)
  
  #----------------------------------------------------------#
  # 2. Data smoothing ----- 
  #----------------------------------------------------------#
  data_smooth <- data_extract
  
  # smoothdata by selected smoothing type
  data_smooth <- 
    RRatepol:::fc_smooth_community_data(
      data_smooth,
      smooth_method = smooth_method,
      smooth_N_points = smooth_N_points,
      smooth_N_max = smooth_N_max,
      smooth_age_range = smooth_age_range,
      round_results = standardise,
      Debug = Debug)
  
  #check data and reduce data dimentions
  data_work <- 
    RRatepol:::fc_check_data(
      data_smooth, 
      proportion = FALSE,
      Debug = Debug)
  
  #----------------------------------------------------------#
  # 3. Working Unit selection -----
  #----------------------------------------------------------#
  
  if(Working_Units == "levels"){
    Number_of_shifts <-  1
  } else if(Working_Units == "bins"){
    Number_of_shifts <-  1
    
    assertthat::assert_that(
      is.numeric(bin_size),
      msg = " `bin_size` must be a `numeric`")
    
    assertthat::assert_that(
      round(bin_size) == bin_size,
      msg = " `bin_size` must be a whole number")
    
    bin_sizes <- 
      RRatepol:::fc_create_bins(
        data_work,
        shift_value = bin_size,
        Number_of_shifts = 1)
  } else if(Working_Units == "MW"){
    
    assertthat::assert_that(
      is.numeric(bin_size),
      msg = " `bin_size` must be a `numeric`")
    
    assertthat::assert_that(
      round(bin_size) == bin_size,
      msg = " `bin_size` must be a whole number")
    
    assertthat::assert_that(
      is.numeric(Number_of_shifts),
      msg = " `Number_of_shifts` must be a `numeric`")
    
    assertthat::assert_that(
      round(Number_of_shifts) == Number_of_shifts,
      msg = " `Number_of_shifts` must be a whole number")
    
    shift_value <-  bin_size/Number_of_shifts
    
    bin_sizes <- 
      RRatepol:::fc_create_bins(
        data_work,
        shift_value = shift_value,
        Number_of_shifts = Number_of_shifts)
  }
  
  
  #----------------------------------------------------------#
  # 4. Randomisation ----- 
  #----------------------------------------------------------#
  
  # select the prefetred number of cores for of cores for parallel computation
  if(class(treads) == "numeric"){
    Ncores <-  treads # set value
  } else {
    if (treads == TRUE) {
      Ncores <-  parallel::detectCores() # detect number
    } else {
      Ncores <-  1
    }
  }
  
  # create cluster
  cl <-  parallel::makeCluster(Ncores)
  doSNOW::registerDoSNOW(cl)
  
  # add all functions to the cluster
  envir <-  environment(fc_estimate_RoC)
  parallel::clusterExport(cl, varlist = c(ls(envir)))
  
  # create progress bar based os the number of replication
  pb <-  utils::txtProgressBar(max = rand, style = 3)
  progress <-  function(n) setTxtProgressBar(pb, n)
  opts <-  list(progress = progress)
  
  # create template for result tibble
  shift_tibble_template <-  tibble::tibble()
  
  result_tibble <-  
    foreach::`%dopar%`(foreach::foreach(
      l = 1:rand,
      .combine = rbind,
      .options.snow=opts), {
        
        # TIME SAMPLING
        # sample random time sequence from time uncern.
        data_work@Age$newage <- 
          as.numeric(data_work@Age.un[sample(c(1:max(1, nrow(data_work@Age.un))), 1), ])
        
        # create result tible
        shift_tibble <-  shift_tibble_template
        
        # repeat for number of shifts
        for(k in 1:Number_of_shifts){
          
          #----------------------------------------------------------#
          # 4.1 Data subsetting ----- 
          #----------------------------------------------------------#
          data_subset <-  data_work
          
          # select one sample for each bin based on the age of the samples. Sample is chones if it is the closes one to the upper end of the bin
          if (Working_Units != "levels"){
            
            # select bin for this shift
            selected_bins <-  bin_sizes[bin_sizes$shift == k, ]
            
            #subset data
            data_subset <-  RRatepol:::fc_subset_samples(data_subset, selected_bins, Working_Units)
            
            data_subset <-  RRatepol:::fc_check_data(data_subset, proportion = FALSE)
            
          }
          
          #----------------------------------------------------------#
          # 4.2 Data Standardisation ----- 
          #----------------------------------------------------------#
          data_sd <-  data_subset
          
          # standardisation of community data to X(N_individuals) number of individuals
          if(standardise == TRUE){
            
            assertthat::assert_that(
              is.numeric(N_individuals),
              msg = " `N_individuals` must be a `numeric`")
            
            assertthat::assert_that(
              round(N_individuals) == N_individuals,
              msg = " `N_individuals` must be a whole number")
            
            # adjust the value by the minimal Community or to a minimal of presected values
            N_individuals <-  min(c(rowSums(data_sd@Community), N_individuals) )
            
            # check if all samples has N_individuals of individuals
            data_sd@Age <-  data_sd@Age[rowSums(data_sd@Community, na.rm = TRUE) >= N_individuals, ]
            data_sd@Age.un <-  data_sd@Age.un[ ,rowSums(data_sd@Community, na.rm = TRUE) >= N_individuals]
            data_sd@Community <-  data_sd@Community[rowSums(data_sd@Community, na.rm = TRUE) >= N_individuals, ]
            data_sd <-  RRatepol:::fc_check_data(data_sd, proportion = FALSE, Samples = TRUE, Debug = Debug)
            
            # standardisation
            data_sd <-  RRatepol:::fc_standardise_community_data(data_sd, N_individuals, Debug = Debug)
            
            assertthat::assert_that(
              any(rowSums(data_sd@Community, na.rm = TRUE) == N_individuals),
              msg = "Data standardisation was unsuccesfull, try `standardise` = FALSE")
          }
          
          # data check with proportioning
          data_sd_check <-  
            RRatepol:::fc_check_data(data_sd,
                                     proportion = tranform_to_proportions,
                                     Samples = FALSE, Debug = Debug)
          
          #----------------------------------------------------------#
          # 4.3 DC Calculation ----- 
          #----------------------------------------------------------#
          
          # calculate DC between each subsequent samples/bins
          DC_res <-  fc_calculate_DC(data_sd_check, DC = DC, Debug = Debug)
          
          #----------------------------------------------------------#
          # 4.4 Age Standardisation ----- 
          #----------------------------------------------------------#
          
          # create empty vector with size = number of samples-1
          sample_size_work <-  data_sd_check$Dim.val[2]-1
          
          # create empty vectors for age difference calcualtion
          age_diff <-  vector(mode = "numeric", length = sample_size_work)
          age_diff_names <-  vector(mode = "character", length = sample_size_work)
          age_mean <-  age_diff
          
          for (i in 1:sample_size_work){ # for each RoC
            
            # calcualte the age difference between subsequesnt samples
            age_diff[i] <-  data_sd_check$Age$newage[i + 1] - data_sd_check$Age$newage[i]
            
            # Set age difference as 1, if age difference between samples is smaller than 1
            if(age_diff[i] < 1){ age_diff[i] <-  1}
            
            #calculate the average position of RoC
            age_mean[i] <-  mean(c(data_sd_check$Age$age[i+1],
                                   data_sd_check$Age$age[i]))
            
            # create vector with bin names
            age_diff_names[i] <-
              paste(
                row.names(data_sd_check$Age)[i],
                "-",
                row.names(data_sd_check$Age)[i + 1])
          }
          
          if (Debug == TRUE){
            cat("", fill = TRUE)
            cat(paste("The time standardisation unit (TSU) is",round(mean(age_diff),2)), fill=TRUE)
          }
          
          #  calculate DC standardise by time
          DC_res_s <-  vector(mode = "numeric", length = sample_size_work)
          DC_res_s <-  (DC_res * mean(age_diff)) / age_diff
          
          #----------------------------------------------------------#
          # 4.5 Result of a single window shift -----
          #----------------------------------------------------------#
          
          # add the results from this shift into the result tibble
          shift_tibble <- 
            rbind(
              shift_tibble,
              data.frame(
                bin = age_diff_names,
                DC = DC_res,
                age_position = age_mean,
                age_diff = age_diff,
                RoC = DC_res_s,
                shift = rep(k, sample_size_work))
            )
        }
        
        #----------------------------------------------------------#
        # 4.6 Result of a single randomisation run -----
        #----------------------------------------------------------#
        
        # save result from single randomisation into data.frame with number of randomisation as ID.
        data_result_temp <- 
          as.data.frame(
            list(
              ID = l,
              shift_tibble))
        
        return(data_result_temp)
      })# end of the randomization
  
  # close progress bar and cluster
  close(pb)
  parallel::stopCluster(cl)
  
  
  #----------------------------------------------------------#
  # 5. Results Summary -----
  #----------------------------------------------------------#
  
  # create new dataframe with summary of randomisation results
  
  # extract results and match them by bin
  results_full <- 
    dplyr::right_join(
      fc_extract_result(
        result_tibble,
        "RoC",
        rand),
      fc_extract_result(
        result_tibble,
        "age_position",
        rand),
      by = c("sample_id","shift"))
  
  
  # reduce results by the focus age time
  if (interest_threshold != FALSE){
    results_full <- 
      dplyr::filter(
        results_full,
        age_position <= interest_threshold)
  }
  
  # sort samples by age and add smoothing to avoid "waiving" from different shifts
  results_full <-  results_full[order(results_full$age_position), ]
  results_full$RoC_sm <-  stats::lowess(results_full$age_position, results_full$RoC, f = .1, iter = 100)$y
  results_full$RoC_95q_sm <-  stats::lowess(results_full$age_position, results_full$RoC_95q, f = .1, iter = 100)$y
  results_full$RoC_05q_sm <-  stats::lowess(results_full$age_position, results_full$RoC_05q, f = .1, iter = 100)$y
  results_full$RoC_sm <-  ifelse(results_full$RoC_sm <= 0, 0.0001, results_full$RoC_sm)
  results_full$RoC_05q_sm <-  ifelse(results_full$RoC_05q_sm < 0, 0.0001, results_full$RoC_05q_sm)
  
  
  # final tibble
  results_full_fin <- 
    dplyr::select(
      results_full, 
      sample_id,
      age_position,
      RoC_sm,
      RoC_95q_sm,
      RoC_05q_sm)
  
  names(results_full_fin) <-  c("Working_Unit","Age", "ROC", "ROC_up", "ROC_dw")
  
  end.time <-  Sys.time()
  time.length <-  end.time - start.time
  
  if(Debug == TRUE){
    cat("", fill=TRUE)
    cat(paste(
      "RATEPOL finished", end.time, "taking", time.length, units(time.length)), fill=TRUE)
  }
  
  
  return(results_full_fin)
  
}
# end of code
