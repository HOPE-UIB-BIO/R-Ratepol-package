fc_estimate_RoC <- function(data_source_community,
                            data_source_age,
                            age_uncertainty = F,
                            smooth_method = "none",
                            smooth_N_points = 5,
                            smooth_age_range = 500,
                            smooth_N_max = 9,
                            Working_Units = "levels",
                            bin_size = 500,
                            Number_of_shifts = 1,
                            rand = 1,
                            treads = F,
                            standardise = F,
                            N_individuals = 150,
                            tranform_to_proportions = T,
                            DC = "euc",
                            interest_threshold = F,
                            Peak_detection = "GAM",
                            Debug = F){
  
  # Start of the code
  start.time <- Sys.time()
  if (Debug == T){
    cat(paste("RATEPOL started", start.time), fill=T)
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
      Debug=Debug)

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
      Debug=Debug)

  #check data and reduce data dimentions
  data_work <- 
    fc_check_data(
      data_smooth, 
      proportion = F,
      Debug=Debug)

  #----------------------------------------------------------#
  # 3. Working Unit selection -----
  #----------------------------------------------------------#

  if(Working_Units == "levels"){
    Number_of_shifts <-  1
  }

  if(Working_Units == "selected_bins"){
    Number_of_shifts <-  1
    
    bin_sizes <- 
      fc_create_bins(
        data_work,
        shift_value = bin_size,
        Number_of_shifts = 1)
  }

  if(Working_Units == "MW"){
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

  # select the prefetred number of cores for of cores for parallel computation
  if(class(treads) == "numeric"){
    Ncores <-  treads # set value
  } else {
    if (treads == T) {
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
      l=1:rand,
      .combine = rbind,
      .options.snow=opts), {

    # TIME SAMPLING
    # sample random time sequence from time uncern.
    data_work$Age$newage <- 
      as.numeric(data_work$Age.un[sample(c(1:max(1,nrow(data_work$Age.un))),1),])

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
        data_subset <-  fc_subset_samples(data_subset, selected_bins)
      }

      #----------------------------------------------------------#
      # 4.2 Data Standardisation ----- 
      #----------------------------------------------------------#
      data_sd <-  data_subset

      # standardisation of community data to X(N_individuals) number of individuals
      if(standardise == T){

        # adjust the value by the minimal Community or to a minimal of presected values
        N_individuals <-  min(c( rowSums(data_subset$Community), N_individuals) )

        # check if all samples has N_individuals of individuals
        data_sd$Age <-  data_sd$Age[rowSums(data_sd$Community, na.rm = T) >= N_individuals, ]
        data_sd$Age.un <-  data_sd$Age.un[ ,rowSums(data_sd$Community, na.rm = T) >= N_individuals]
        data_sd$Community <-  data_sd$Community[rowSums(data_sd$Community, na.rm = T) >= N_individuals, ]
        data_sd <-  fc_check_data(data_sd, proportion = F, Samples = T, Debug=Debug)

        # standardisation
        data_sd <-  fc_standardise_community_data(data_sd, N_individuals, Debug = Debug)

        if(any(rowSums(data_sd$Community, na.rm = T)!=N_individuals))
          stop("standardisation was unsuccesfull")
      }

      # data check with proportioning
      data_sd_check <-  fc_check_data(data_sd, proportion = tranform_to_proportions, Samples = F, Debug = Debug)

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

      if (Debug == T){
        cat("", fill = T)
        cat(paste("The time standardisation unit (TSU) is",round(mean(age_diff),2)), fill=T)
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
  # 5. Results Summary
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
  if (interest_threshold != F){
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


  #----------------------------------------------------------#
  # 6. Peak Deetction
  #----------------------------------------------------------#

  # Median peak threshold
  if(Peak_detection == "Threshold"){
    
    # threshold for RoC peaks is set as median of all RoC in dataset
    r_threshold <- median(results_full$RoC_sm)
    
    # mark peaks which have 95% quantile above the threshold as Peak
    results_full$Peak <-  results_full$RoC_05q_sm > r_threshold
  }

  # Median peak threshold
  if(Peak_detection == "linear"){
    
    # mark points that are abowe the linear model (exactly 1.5 SD higher than prediction)
    results_full$pred_linear <-
      predict.lm(
        lm(RoC_sm~age_position,
          data = results_full),
        type="response")
    
    results_full$pred_linear_diff <-  results_full$RoC_sm - results_full$pred_linear
    results_full$Peak <-  (results_full$pred_linear_diff) > 1.5 * sd(results_full$pred_linear_diff)
  }
  
  # GAM
  if(Peak_detection == "GAM"){
    
    # mark points that are abowe the GAM model (exactly 1.5 SD higher than GAM prediction)
    results_full$pred_gam <-  
      mgcv::predict.gam(
        mgcv::gam(
          RoC_sm~s(age_position,k=3),
          data = results_full,
          family = "Gamma()",
          method = "REML"),
        type="response")
    
    results_full$pred_gam_diff <-  results_full$RoC_sm - results_full$pred_gam
    results_full$Peak <-  (results_full$pred_gam_diff) > 1.5 * sd(results_full$pred_gam_diff)

  }
  
  # GAM
  if(Peak_detection == "GAM_deriv"){
    
    # fit gam well smother gam model and use first derivative of the function to detect signifiant increases in the function 
    gam_model <-  
        mgcv::gam(
          RoC_sm~s(age_position),
          data = results_full,
          family = "Gamma()",
          method = "REML")
    
    new_data  <-
      tibble(age_position = results_full$age_position)
    
    gam_deriv <-  gratia::fderiv(gam_model, newdata = new_data , n = 1000)
    
    gam_deriv_conf <-
      with(new_data, cbind(
        confint(
          gam_deriv,
          nsim = 1000,
          type = "simultaneous",
          transform	= T
        ),
        age_position = age_position
      )) %>%
      mutate(Peak = upper < 0)
    
    results_full$Peak <-  gam_deriv_conf$Peak
    
  }

  # SNI
  if (Peak_detection == "SNI"){
    
    # set moving window of 5 times higher than average distance between samples
    mean_age_window <- 5 * mean( diff(results_full$age_position) )
    
    # create GAM
    pred_gam <-  
      mgcv::predict.gam(
        mgcv::gam(
          RoC_sm~s(age_position, k=3),
          data = results_full,
          family = "Gamma()",
          method = "REML"),
        type="response")
    
    # calculate SNI (singal to noise ratio)
    SNI_calc <- 
      fc_CharSNI(
        data.frame(
          results_full$age_position,
          results_full$RoC_sm,
          pred_gam),
        mean_age_window)
    
    # mark points with SNI higher than 3
    results_full$Peak <-  SNI_calc$SNI > 3 & results_full$RoC_sm > pred_gam
  }

  # outro
  results_full_fin <- 
    dplyr::select(
      results_full, 
      sample_id,
      RoC_sm,
      RoC_95q_sm,
      RoC_05q_sm,
      age_position,
      Peak)
  
  names(results_full_fin) <-  c("Working_Unit", "ROC", "ROC_up", "ROC_dw", "Age", "Peak")

  end.time <-  Sys.time()
  time.length <-  end.time - start.time
  
  if(Debug == T){
    cat("", fill=T)
    cat(paste(
      "RATEPOL finished", end.time, "taking", time.length, units(time.length)), fill=T)
  }


  return(results_full_fin)

}
# end of code
