fc_detect_peak_points <-  function(data_source, method = "trend_linear"){
  
  
  if(any(method == c("threshold", "trend_linear", "trend_non_linear", "GAM_deriv", "SNI")) == F)
    stop("Method has to be one of the following 'threshold', 'trend_linear', 'trend_non_linear', 'GAM_deriv', 'SNI'")
  
  
  #----------------------------------------------------------#
  # 1. Median peak threshold ----- 
  #----------------------------------------------------------#
  if(method == "threshold"){
    
    # threshold for RoC peaks is set as median of all RoC in dataset
    r_threshold <- median(data_source$ROC)
    
    # mark peaks which have 95% quantile above the threshold as Peak
    data_source$Peak <-  data_source$ROC_dw > r_threshold
  }
  
  #----------------------------------------------------------#
  # 2. Linear trend  ----- 
  #----------------------------------------------------------#
  if(method == "trend_linear"){
    
    # mark points that are abowe the linear model (exactly 1.5 SD higher than prediction)
    data_source$pred_linear <-
      predict.glm(
        glm(ROC ~ Age,
           data = data_source,
           family = stats::Gamma()),
        type = "response")
    
    data_source$pred_linear_diff <-  data_source$ROC - data_source$pred_linear
    data_source$Peak <-  (data_source$pred_linear_diff) > 1.5 * sd(data_source$pred_linear_diff)
  }
  
  #----------------------------------------------------------#
  # 3. Non-linear trend  ----- 
  #----------------------------------------------------------#
  if(method == "trend_non_linear"){
    
    # mark points that are abowe the GAM model (exactly 1.5 SD higher than GAM prediction)
    data_source$pred_gam <-  
      mgcv::predict.gam(
        mgcv::gam(
          ROC ~ s(Age, k = 3),
          data = data_source,
          family = stats::Gamma(),
          method = "REML"),
        type="response")
    
    data_source$pred_gam_diff <-  data_source$ROC - data_source$pred_gam
    data_source$Peak <-  (data_source$pred_gam_diff) > 1.5 * sd(data_source$pred_gam_diff)
    
  }
  
  #----------------------------------------------------------#
  # 4. Firts derivative of GAM model  ----- 
  #----------------------------------------------------------#
  if(method == "GAM_deriv"){
    
    # fit gam well smother gam model and use first derivative of the function to detect signifiant increases in the function 
    gam_model <-  
      mgcv::gam(
        ROC ~ s(Age),
        data = data_source,
        family = stats::Gamma(),
        method = "REML")
    
    new_data  <-
      tibble::tibble(Age = data_source$Age)
    
    gam_deriv <-  gratia::fderiv(gam_model, newdata = new_data , n = 1000)
    
    gam_deriv_conf <-
      with(new_data, cbind(
        confint(
          gam_deriv,
          nsim = 1000,
          type = "simultaneous",
          transform	= T
        ),
        Age = Age
      )) 
    data_source$Peak <-  gam_deriv_conf$upper < 0
    
  }
  
  #----------------------------------------------------------#
  # 5. Signal-to-Noise-ratio Index  ----- 
  #----------------------------------------------------------#
  if (method == "SNI"){
    
    # set moving window of 5 times higher than average distance between samples
    mean_age_window <- 5 * mean( diff(data_source$Age) )
    
    # create GAM
    pred_gam <-  
      mgcv::predict.gam(
        mgcv::gam(
          ROC ~ s(Age, k = 3),
          data = data_source,
          family = stats::Gamma(),
          method = "REML"),
        type = "response")
    
    # calculate SNI (singal to noise ratio)
    SNI_calc <- 
      fc_CharSNI(
        data.frame(
          data_source$Age,
          data_source$ROC,
          pred_gam),
        mean_age_window)
    
    # mark points with SNI higher than 3
    data_source$Peak <-  SNI_calc$SNI > 3 & data_source$ROC > pred_gam
  }
  
  #----------------------------------------------------------#
  # 6. save result  ----- 
  #----------------------------------------------------------#
  
  data_result <- 
    dplyr::select(
      data_source, 
      Working_Unit,
      Age,
      ROC,
      ROC_up,
      ROC_dw,
      Peak)
  
  return(data_result)
  
}