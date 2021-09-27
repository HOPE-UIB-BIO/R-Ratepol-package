fc_estimate_RoC <- function(data_source_community,
                            data_source_age,
                            age_uncertainty,
                            smooth_method = c("none", "m.avg", "grim", "age.w", "shep"),
                            smooth_N_points = 5,
                            smooth_age_range = 500,
                            smooth_N_max = 9,
                            Working_Units = c("levels", "bins", "MW"),
                            bin_size = 500,
                            Number_of_shifts = 5,
                            bin_selection = c( "first", "random"),
                            rand = 1,
                            treads = TRUE,
                            standardise = FALSE,
                            N_individuals = 150,
                            tranform_to_proportions = TRUE,
                            DC = c("euc", "euc.sd", "chord", "chisq", "gower"),
                            interest_threshold = FALSE,
                            only_subsequent = TRUE,
                            time_standardisation = 100,
                            Debug = FALSE){
  
  # Start of the code
  
  #----------------------------------------------------------# 
  # 0. Arguments check -----
  #----------------------------------------------------------#
  
  assertthat::assert_that(
    !missing(data_source_community),
    msg = "Object 'data_source_community' must be included as a 'data.frame'")
  
  assertthat::assert_that(
    !missing(data_source_age),
    msg = "Object 'data_source_age' must be included as a 'data.frame'")
  
  if(missing(age_uncertainty)){
    age_uncertainty <- FALSE
  }
  
  assertthat::assert_that(
    is.numeric(rand),
    msg = "'rand' must be a 'numeric'")
  
  assertthat::assert_that(
    round(rand) == rand,
    msg = "'rand' must be a whole number")
  
  Working_Units <- match.arg(Working_Units)
  
  assertthat::assert_that(
    Working_Units %in% c("levels", "bins", "MW"),
    msg = "'Working_Units' must be a 'levels' 'bins' or 'MW'")
  
  assertthat::assert_that(
    is.numeric(time_standardisation) | time_standardisation == "auto",
    msg = "'time_standardisation' must be a 'numeric' or 'auto'")
  
  if( is.numeric(time_standardisation)){
    assertthat::assert_that(
      round(time_standardisation) == time_standardisation,
      msg = "'time_standardisation' must be a whole number")  
  }
  
  if(Working_Units != "levels"){
    
    assertthat::assert_that(
      is.numeric(bin_size),
      msg = "'bin_size' must be a 'numeric'")
    
    assertthat::assert_that(
      round(bin_size) == bin_size,
      msg = "'bin_size' must be a whole number")
    
    bin_selection <- match.arg(bin_selection)
    
    assertthat::assert_that(
      bin_selection == "first" | bin_selection == "random",
      msg = "'bin_selection' must be a 'first' or 'random'")
    
    assertthat::assert_that(
      is.logical(only_subsequent),
      msg = "'only_subsequent' must be a 'TRUE' or 'FALSE'")
    
    if (Working_Units == "MW"){
      assertthat::assert_that(
        is.numeric(Number_of_shifts),
        msg = "'Number_of_shifts' must be a 'numeric'")
      
      assertthat::assert_that(
        round(Number_of_shifts) == Number_of_shifts,
        msg = "'Number_of_shifts' must be a whole number")
    }
  }  
  
  assertthat::assert_that(
    is.logical(standardise),
    msg = "'standardise' must be a 'TRUE' or 'FALSE'")
  
  if(standardise == TRUE){
    assertthat::assert_that(
      is.numeric(N_individuals),
      msg = "'N_individuals' must be a 'numeric'")
    
    assertthat::assert_that(
      round(N_individuals) == N_individuals,
      msg = "'N_individuals' must be a whole number")
  }
  
  assertthat::assert_that(
    is.logical(tranform_to_proportions),
    msg = "'tranform_to_proportions' must be a 'TRUE' or 'FALSE'")
  
  assertthat::assert_that(
    is.numeric(interest_threshold) | interest_threshold == FALSE,
    msg = "'interest_threshold' must be a 'numeric' or 'FALSE'")
  
  smooth_method <-  match.arg(smooth_method)
  
  assertthat::assert_that(
    any(smooth_method == c("none", "m.avg", "grim", "age.w", "shep")),
    msg = "'smooth_method' must be one of the following:
    'none', 'm.avg', 'grim', 'age.w', 'shep'")
  
  if(!smooth_method %in% c("none", "shep")){
    
    assertthat::assert_that(
      smooth_N_points%%2 != 0,
      msg = "'smooth_N_points' must be an odd number")
    
    if(smooth_method != "m.avg"){
      assertthat::assert_that(
        is.numeric(smooth_age_range),
        msg = "'smooth_age_range' must be 'numeric")
      
      if(smooth_method == "grim"){
        assertthat::assert_that(
          smooth_N_max%%2 != 0,
          msg = "'smooth_N_max' must be an odd number")
        
        assertthat::assert_that(
          smooth_N_points < smooth_N_max,
          msg = "'smooth_N_max' must be bigger than 'smooth_N_points")
      }
    }
  }
  
  DC <- match.arg(DC)
  
  assertthat::assert_that(
    any(DC == c("euc", "euc.sd", "chord", "chisq", "gower")),
    msg = "'DC' must be one of the following:
    'euc', 'euc.sd', 'chord', 'chisq', 'gower'")
  
  assertthat::assert_that(
    is.numeric(treads) | is.logical(treads),
    msg = "'treads' must be a 'numeric' or 'TRUE'/'FALSE'")
  
  if(is.numeric(treads)){
    assertthat::assert_that(
      round(treads) == treads,
      msg = "'treads' must be a whole number")
  }
  
  
  #--------------------------------------------------#
  # 0.1. Report to user -----
  #--------------------------------------------------#
  cat("\n")
  
  start_time <- Sys.time()
  cat(paste("R-RATEPOL started", start_time),
      "\n", fill = TRUE)
  
  if (!all(age_uncertainty == FALSE)){
    cat(
      "'age_uncertainty' will be used for in the RoC estimation",
      "\n", fill = TRUE)
    
    if (rand < 100){
      cat(
        paste(
          "'age_uncertainty' was selected to be used with low number",
          "of replication. Recommend to increase 'rand'"),
        "\n", fill = TRUE)
    }
  }
  
  if(smooth_method == "m.avg"){
    cat(
      paste(
        "Data will be smoothed by 'moving average' over", smooth_N_points,
        "points"),
      "\n",fill = TRUE)  
  }
  
  if(smooth_method == "grim"){
    cat(
      paste(
        "Data will be smoothed by 'Grimm method' with min samples", smooth_N_points,
        "max samples", smooth_N_max, "and max age range of", smooth_age_range),
      "\n", fill = TRUE)  
  }
  
  if(smooth_method == "age.w"){
    cat(
      paste(
        "Data will be smoothed by 'age-weighed average' over", smooth_N_points,
        "points with a threshold of", smooth_age_range),
      "\n", fill = TRUE)  
  }
  
  if(smooth_method == "shep"){
    cat(
      paste(
        "Data will be smoothed by 'Shepard's 5-term filter'"),
      "\n", fill = TRUE)  
  }
  
  if(Working_Units == "levels"){
    cat(
      "RoC will be estimated between individual subsequent levels",
      "\n", fill = TRUE)  
  }
  
  if(Working_Units != "levels"){
    
    if(Working_Units == "bins"){
      cat(
        paste(
          "RoC will be estimated using selective binning with", bin_size,
          "yr time bin"),
        "\n", fill = TRUE)  
    }
    
    if(Working_Units == "MW"){
      cat(
        paste(
          "RoC will be estimated using 'binning with the mowing window' of",
          bin_size, "yr time bin over", Number_of_shifts, "number of window shifts"),
        "\n", fill = TRUE)  
    }
    
    if(bin_selection == "random"){
      cat(
        "Sample will randomly selected for each bin",
        "\n", fill = TRUE)  
      
      if(rand < 100){
        cat(
          paste(
            "'bin_selection' was selected as 'random' with low number",
            "of replication. Recommend to increase 'rand'"),
          "\n", fill = TRUE)  
      }
      
    } else {
      cat(
        "First sample of each time bin will selected",
        "\n", fill = TRUE)  
    }
    
    if(only_subsequent == FALSE & Working_Units == "MW"){
      cat(
        paste(
          "'only_subsequent == FALSE' and 'Working_Units == MW'.",
          "This is not a recommended setting.",
          "Please use 'only_subsequent == TRUE' for 'Working_Units == MW'",
          "see '?fc_estimate_ROC' for more information"), 
        "\n", fill = TRUE)
    } else if(only_subsequent == FALSE){
      cat(
        paste(
          "'only_subsequent' was selected as 'FALSE'.",
          "This is not a recommended setting. Results will be affected",
          "see '?fc_estimate_ROC' for more information"),
        "\n", fill = TRUE)
    }
    
  }
  
  if(is.numeric(time_standardisation)){
    cat(
      paste(
        "'time_standardisation' =", time_standardisation,":",
        "RoC values will be reported as disimilarity per", time_standardisation,
        "years."),
      "\n", fill = TRUE)
    
    if(Working_Units != "levels" & time_standardisation != bin_size){
      cat(
        paste(
          "RoC values will be reported in different units than size of bin.",
          "Recommend to keep 'time_standardisation'",
          "and 'bin_size' as same values"),
        "\n", fill = TRUE)
    }
  }
  
  if(time_standardisation == "auto"){
    cat(
      paste(
        "'time_standardisation' = 'auto' is not recomended setting.", 
        "RoC values will be reported as standardised by the average distance", 
        "between Working Units (levels/ bins)"),
      "\n", fill = TRUE)
  }
  
  if(standardise == TRUE ){
    cat(
      paste("Data will be standardise in each Working unit to", N_individuals,
            "or the lowest number detected in dataset"),
      "\n", fill=TRUE)
    
    if(rand < 100){
      cat(
        paste(
          "'standardise' was selected as 'TRUE' with low number of replication.",
          "Recommend to increase 'rand'"),
        "\n", fill = TRUE) 
    }
  }
  
  
  #----------------------------------------------------------# 
  # 1. Data extraction -----
  #----------------------------------------------------------#
  
  # extract data into working format
  # already include data check
  data_extract <- 
    fc_extract_data(
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
    fc_smooth_community_data(
      data_smooth,
      smooth_method = smooth_method,
      smooth_N_points = smooth_N_points,
      smooth_N_max = smooth_N_max,
      smooth_age_range = smooth_age_range,
      round_results = standardise,
      Debug = Debug)
  
  #check data and reduce data dimentions
  data_work <- 
    fc_check_data(
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
    
    bin_sizes <- 
      fc_create_bins(
        data_work,
        shift_value = bin_size,
        Number_of_shifts = 1)
  } else if(Working_Units == "MW"){
    
    shift_value <-  bin_size/Number_of_shifts
    
    bin_sizes <- 
      fc_create_bins(
        data_work,
        shift_value = shift_value,
        Number_of_shifts = Number_of_shifts)
  }
  
  
  #----------------------------------------------------------#
  # 4. Randomisation ----- 
  #----------------------------------------------------------#
  
  # create template for result tibble
  shift_tibble_template <-  tibble::tibble()
  
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
  doParallel::registerDoParallel(cl)
  parallel::clusterEvalQ(cl, {
    library("RRatepol")
    library("tidyverse")
  })
  
  if(rand > 1){
    cat(
      paste(
        "Starting the randomisation. Number of randomisations set to", rand),
      "\n", fill = TRUE)  
  }
  
  # progress bar
  cat("Progress bar is not working in the current version, please wait",
      "\n", fill=TRUE)
  
  result_tibble <-  
    foreach::`%dopar%`(foreach::foreach(
      l = 1:rand,
      .combine = rbind), {
        
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
          
          # select one sample for each bin based on the age of the samples. 
          # Sample is choses if it is the closes one to the upper end of the bin
          if (Working_Units != "levels"){
            
            # select bin for this shift
            selected_bins <-  bin_sizes[bin_sizes$shift == k, ]
            
            #subset data
            data_subset <-  
              fc_subset_samples(
                data_subset,
                bins = selected_bins,
                WU = Working_Units,
                bin_selection = bin_selection)
            
            data_subset <-  
              fc_check_data(
                data_subset,
                proportion = FALSE)
          }
          
          #----------------------------------------------------------#
          # 4.2 Data Standardisation ----- 
          #----------------------------------------------------------#
          data_sd <-  data_subset
          
          # standardisation of community data to X(N_individuals) number of individuals
          if(standardise == TRUE){
            
            # adjust the value by the minimal Community or to a minimal of presected values
            N_individuals <-  min(c(rowSums(data_sd@Community), N_individuals) )
            
            # check if all samples has N_individuals of individuals
            data_sd@Age <-  
              data_sd@Age[rowSums(data_sd@Community, na.rm = TRUE) >= N_individuals, ]
            
            data_sd@Age.un <-  
              data_sd@Age.un[ ,rowSums(data_sd@Community, na.rm = TRUE) >= N_individuals]
            
            data_sd@Community <-  
              data_sd@Community[rowSums(data_sd@Community, na.rm = TRUE) >= N_individuals, ]
            
            data_sd <-  
              fc_check_data(
                data_sd,
                proportion = FALSE,
                Samples = TRUE,
                Debug = Debug)
            
            # standardisation
            data_sd <-  
              fc_standardise_community_data(
                data_sd
                , N_individuals,
                Debug = Debug)
            
            assertthat::assert_that(
              any(rowSums(data_sd@Community, na.rm = TRUE) == N_individuals),
              msg = "Data standardisation was unsuccesfull, try 'standardise' = FALSE")
          }
          
          # data check with proportioning
          data_sd_check <-  
            fc_check_data(
              data_sd,
              proportion = tranform_to_proportions,
              Samples = FALSE, Debug = Debug)
          
          
          #----------------------------------------------------------#
          # 4.3 DC Calculation ----- 
          #----------------------------------------------------------#
          
          # calculate DC between each subsequent samples/bins
          DC_res <-  
            fc_calculate_DC(
              data_sd_check,
              DC = DC,
              Debug = Debug)
          
          
          #----------------------------------------------------------#
          # 4.4 Age Standardisation ----- 
          #----------------------------------------------------------#
          
          # create empty tible with size = number of samples-1
          shift_tibble_res <- 
            tibble::tibble(
              DC = DC_res) 
          
          # create empty vectors for age difference calcualtion
          shift_tibble_res$age_diff <-  
            vector(
              mode = "numeric",
              length = nrow(shift_tibble_res))
          
          shift_tibble_res$bin <-  
            vector(
              mode = "character",
              length = nrow(shift_tibble_res))
          
          shift_tibble_res$age_distance <- 
            vector(
              mode = "numeric",
              length = nrow(shift_tibble_res))
          
          shift_tibble_res$age_position <-  
            vector(
              mode = "numeric",
              length = nrow(shift_tibble_res))
          
          for (i in 1:nrow(shift_tibble_res)){ # for each RoC
            
            # calcualte the age difference between subsequesnt samples
            shift_tibble_res$age_diff[i] <-  
              data_sd_check@Age$newage[i + 1] - data_sd_check@Age$newage[i]
            
            # Set age difference as 1, if age difference between samples is 
            #   smaller than 1
            if(shift_tibble_res$age_diff[i] < 1){shift_tibble_res$age_diff[i] <-  1}
            
            #calculate the average position of RoC
            shift_tibble_res$age_position[i] <-  
              mean(c(data_sd_check@Age$age[i+1],
                     data_sd_check@Age$age[i]))
            
            # create vector with WU names
            shift_tibble_res$bin[i] <-
              paste(
                row.names(data_sd_check@Age)[i],
                "-",
                row.names(data_sd_check@Age)[i + 1])
            
            if(Working_Units != "levels"){
              shift_tibble_res$age_distance[i] <- 
                as.numeric(row.names(data_sd_check@Age)[i + 1]) -
                as.numeric(row.names(data_sd_check@Age)[i])
            } else {
              shift_tibble_res$age_distance[i] <- NA
            }
          }
          
          # remove the non-subsequent levels.
          if(Working_Units != "levels" & only_subsequent == TRUE){
            shift_tibble_res <-
              shift_tibble_res %>% 
              dplyr::filter(age_distance <= bin_size)
          }
          
          if(time_standardisation == "auto"){
            if(Working_Units != "levels"){
              time_standardisation_unit <- bin_size
            } else {
              time_standardisation_unit <- mean(shift_tibble_res$age_diff)
            }
          } else {
            time_standardisation_unit <- time_standardisation
          }
          
          if (Debug == TRUE){
            cat("", fill = TRUE)
            cat(paste("The time standardisation unit (TSU) is",
                      round(time_standardisation_unit,2)), fill=TRUE)
          }
          
          #  calculate DC standardise by time
          shift_tibble_res <-
            shift_tibble_res %>% 
            dplyr::mutate(
              age_diff_st = age_diff / time_standardisation_unit,
              RoC = DC / age_diff_st,
              shift = k)
          
          
          #----------------------------------------------------------#
          # 4.5 Result of a single window shift -----
          #----------------------------------------------------------#
          
          # add the results from this shift into the result tibble
          shift_tibble <- 
            rbind(
              shift_tibble,
              shift_tibble_res)
        }
        
        
        #----------------------------------------------------------#
        # 4.6 Result of a single randomisation run -----
        #----------------------------------------------------------#
        
        if(nrow(shift_tibble)<1 & Working_Units != "levels" & only_subsequent == TRUE){
          stop("Estimation not succesfull, try increase the bin size")
        }
        
        shift_tibble <- 
          shift_tibble %>% 
          dplyr::mutate(
            ID = l)
        
        return(shift_tibble)
      })# end of the randomization
  
  # close progress bar and cluster
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
      by = c("sample_id","shift","age_distance"))
  
  # reduce results by the focus age time
  if(interest_threshold != FALSE){
    results_full <- 
      dplyr::filter(
        results_full,
        age_position <= interest_threshold)
  }
  
  # final tibble (sort samples by age and select variables)
  results_full_fin <- 
    results_full %>% 
    dplyr::arrange(age_position) %>% 
    dplyr::select( 
      sample_id,
      age_position,
      RoC,
      RoC_95q,
      RoC_05q)
  
  names(results_full_fin) <-  c("Working_Unit","Age", "ROC", "ROC_up", "ROC_dw")
  
  end_time <-  Sys.time()
  time_duration <-  end_time - start_time
  cat(paste(
    "R-RATEPOL finished", end_time, "taking", time_duration, units(time_duration)),
    fill = TRUE)
  
  return(results_full_fin)
  
}
# end of code
