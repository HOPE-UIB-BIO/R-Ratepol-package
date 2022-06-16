#' @title Detect significant peak points
#'
#' @param data_source Data.frame. Output of `fc_estimate_RoC` function 
#' @param method 
#' Character. A method to use for peak-poit detection:
#' \itemize{
#' \item `"threshold"` - Each point in the RoC sequence is compared to a median 
#' of all RoC scores from the whole sequence (i.e. threshold value). The RoC 
#' value for a point is considered significant if the 95th quantile of the RoC 
#' scores from all calculations is higher than the threshold value.
#' \item `"trend_linear"` - A linear model is fitted between the RoC values and 
#' their ages. Differences between the model and each point are calculated (residuals).
#' The standard deviation (SD) is calculated from all the residuals. A peak is considered 
#' significant if it is 2 SD higher than the model (`sd_threshold` = 2).
#' \item `"trend_non_linear"` - A conservative generalised additive model (GAM) 
#' is fitted through the RoC scores and their ages (GAM = `RoC ~ s(age, k = 3)` 
#' using the `mgcv` package (Wood, 2011). The distance between each point and 
#' the fitted value is calculated (residuals). The standard deviation (SD) is 
#' calculated from all the residuals. A peak is considered significant if it 
#' is 2 SD higher than the model (`sd_threshold` = 2). 
#' }
#' @param sd_threshold 
#' Numeric. Threshold that SD of residuals are compared to, to determine 
#' peak-point (default = 2)
#'
#' @description Detect points of sudden increase of Rate-of-change values
#' @details 
#' A rapid change in composition or relative abundances of variables within the 
#' sequence can provide a means of comparing RoC between sequences and interpreting 
#' the potential drivers of assemblage change. To detect such significant peak-points 
#' of RoC scores in each sequence, each point is tested to see if it represents 
#' a significant increase in RoC values. There are various ways to detect 
#' peak-points in a time series and R-Ratepol is able to detect peak-points 
#' using five methods:
#' \itemize{
#' \item Threshold (`method` = `"threshold"`) - Each point in the RoC sequence is 
#' compared to a median of all RoC scores from the whole sequence (i.e. threshold value).
#'  The ROC value for a point is considered significant if the 95th quantile of 
#'  the RoC scores from all calculations is higher than the threshold value. 
#' \item Linear trend (`method` = `"trend_linear"`) - A linear model is fitted 
#' between the RoC values and their ages. Differences between the model and each 
#' point are calculated (residuals). The standard deviation (SD) is calculated 
#' from all the residuals. A peak is considered significant if it is 2 SD higher 
#' than the model (`sd_threshold` = 2).
#' \item Non-linear trend (`method` = `"trend_non_linear"`) - A conservative 
#' generalised additive model (GAM) is fitted through the RoC scores and their 
#' ages (GAM = `RoC ~ s(age, k = 3)` using the `mgcv` package (Wood, 2011). 
#' The distance between each point and the fitted value is calculated (residuals). 
#' The standard deviation (SD) is calculated from all the residuals. A peak is 
#' considered significant if it is 2 SD higher than the model (`sd_threshold` = 2). 
#' \item F-deriv GAM  (`method` = `"GAM_deriv"`) - A smooth GAM model is fitted 
#' to the RoC scores and their ages (GAM = `RoC ~ s(age)`. The first derivative 
#' as well as continuous confidence intervals are calculated from the model 
#' using the `gratia` package (Simpson, 2019). A peak is considered significant 
#' if the confidence intervals of the first derivative differ from 0 
#' (for more information see Simpson, 2018).
#' \item Signal-to-noise method (`method` = `"SNI"`) - We adapted SNI from 
#' Kelly et al. (2011), which was developed to detect changes in charcoal 
#' stratigraphical records. SNI is calculated for the whole RoC sequence and 
#' a peak-point is considered significant if it has an SNI value higher than 3. 
#' }
#' @references 
#' Kelly, R.F., Higuera, P.E., Barrett, C.M., Feng Sheng, H., 2011. A signal-to-noise 
#' index to quantify the potential for peak detection in sediment-charcoal records. 
#' Quat. Res. 75, 11–17. https://doi.org/10.1016/j.yqres.2010.07.011
#' 
#' Simpson, G.L., 2019. gratia: graceful'ggplot'–based graphics and other functions 
#' for GAMs fitted using “mgcv.” R Packag. version 0.2–1.
#' 
#' Simpson, G.L., 2018. Modelling palaeoecological time series using generalised 
#' additive models. Front. Ecol. Evol. 6, 1–21. https://doi.org/10.3389/fevo.2018.00149
#' 
#' Wood, S.N., 2011. Fast stable restricted maximum likelihood and marginal 
#' likelihood estimation of semiparametric generalized linear models. 
#' J. R. Stat. Soc. Ser. B Stat. Methodol. 73, 3–36. 
#' https://doi.org/10.1111/j.1467-9868.2010.00749.x
#' @export
#'
#' @examples
#' \dontrun{
#' example_data <-  RRatepol::example_data
#' 
#' sequence_01 <- 
#' fc_estimate_RoC(
#' data_source_community = example_data$pollen_data[[1]],
#' data_source_age = example_data$sample_age[[1]],
#' age_uncertainty = FALSE,
#' smooth_method = "shep",
#' Working_Units = "MW",
#' rand = 1e3, 
#' treads = TRUE,
#' DC = "chisq")
#' 
#' sequence_01_peak <-
#' fc_detect_peak_points(
#' sequence_01,
#' method = "trend_non_linear",
#' sd_threshold = 2)
#' 
#' fc_plot_RoC_sequence(
#' sequence_01_peak,
#' age_threshold = 8e3,
#' Roc_threshold = 2,
#' Peaks = TRUE,
#' trend = "trend_non_linear")
#' }
fc_detect_peak_points <-  
  function(
    data_source,
    method = "trend_linear",
    sd_threshold = 2){
    
    
    if(any(method == c("threshold", "trend_linear", "trend_non_linear", "GAM_deriv", "SNI")) == F)
      stop("Method has to be one of the following 'threshold', 'trend_linear', 'trend_non_linear', 'GAM_deriv', 'SNI'")
    
    
    #----------------------------------------------------------#
    # 1. Median peak threshold ----- 
    #----------------------------------------------------------#
    if(method == "threshold"){
      
      # threshold for RoC peaks is set as median of all RoC in dataset
      r_threshold <- stats::median(data_source$ROC)
      
      # mark peaks which have 95% quantile above the threshold as Peak
      data_source$Peak <-  data_source$ROC_dw > r_threshold
    }
    
    #----------------------------------------------------------#
    # 2. Linear trend  ----- 
    #----------------------------------------------------------#
    if(method == "trend_linear"){
      
      # mark points that are abowe the linear model (exactly sd_threshold SD higher than prediction)
      data_source$pred_linear <-
        stats::predict.glm(
          stats::glm(ROC ~ Age,
                     data = data_source,
                     family = mgcv::Tweedie(p=2)),
          type = "response")
      
      data_source$pred_linear_diff <-  data_source$ROC - data_source$pred_linear
      data_source$Peak <-  (data_source$pred_linear_diff) > sd_threshold * stats::sd(data_source$pred_linear_diff)
    }
    
    #----------------------------------------------------------#
    # 3. Non-linear trend  ----- 
    #----------------------------------------------------------#
    if(method == "trend_non_linear"){
      
      # mark points that are abowe the GAM model (exactly sd_threshold SD higher than GAM prediction)
      data_source$pred_gam <-  
        mgcv::predict.gam(
          mgcv::gam(
            ROC ~ s(Age, k = 3),
            data = data_source,
            family = mgcv::Tweedie(p=2),
            method = "REML"),
          type="response")
      
      data_source$pred_gam_diff <-  data_source$ROC - data_source$pred_gam
      data_source$Peak <-  (data_source$pred_gam_diff) > sd_threshold * stats::sd(data_source$pred_gam_diff)
      
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
          family = mgcv::Tweedie(p=2),
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
      data_source$Peak <-  gam_deriv_conf$lower > 0
      
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
            family = mgcv::Tweedie(p=2),
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