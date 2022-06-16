#' @title Function to plot the Rate-of-Change sequence
#'
#' @param data_source 
#' Data.frame. Output of `fc_estimate_RoC` function 
#' @param age_threshold 
#' Numeric. Cut-off value used as maximum age.
#' @param Roc_threshold 
#' Numeric Cut-off value used as maximum RoC value.
#' @param Peaks 
#' Logical. If peak-points are presented in the dataset and `Peaks` == `TRUE`, 
#' then peak points will be displayed  
#' @param trend 
#' If peak-points are presented in the dataset and `Peaks` == `TRUE`,
#'  then one of the three method can be used to visualise the process of peak detection:
#'  \itemize{
#'  \item `"threshold"` - Each point in the RoC sequence is compared to a median 
#'  of all RoC scores from the whole sequence (i.e. threshold value). 
#'  The ROC value for a point is considered significant if the 95th quantile of 
#'  the RoC scores from all calculations is higher than the threshold value.
#'  \item `"trend_linear"` -  A linear model is fitted between the RoC values 
#'  and their ages. Differences between the model and each point are calculated 
#'  (residuals). The standard deviation (SD) is calculated from all the residuals. 
#'  A peak is considered significant if it is 2 SD higher than the model.
#'  \item `"trend_non_linear"` - A conservative generalised additive model (GAM) 
#'  is fitted through the RoC scores and their ages (GAM = `RoC ~ s(age,k=3)` using 
#'  the `mgcv` package (Wood, 2011). The distance between each point and the 
#'  fitted value is calculated (residuals). The standard deviation (SD) is 
#'  calculated from all the residuals. A peak is considered significant if it 
#'  is 2 SD higher than the model. 
#'  }
#' @description Plot Rate-of-Change sequence with a error estimate and trend 
#' and/or peak-points in present.
#' @export
#'
#' @examples
#' \dontrun{
#' 
#' example_data <-  RRatepol::example_data
#' 
#' # example 1
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
#' fc_plot_RoC_sequence(
#' sequence_01,
#' age_threshold = 8e3,
#' Roc_threshold = 1)
#' 
#' # example 2
#' sequence_02 <- 
#' fc_estimate_RoC(
#' data_source_community = example_data$pollen_data[[2]],
#' data_source_age = example_data$sample_age[[2]],
#' age_uncertainty = FALSE,
#' smooth_method = "shep",
#' Working_Units = "MW",
#' rand = 1e3,
#' treads = TRUE,
#' DC = "chisq")
#' 
#' sequence_02_peak <-
#' fc_detect_peak_points(sequence_01, method = "trend_non_linear")
#' 
#' fc_plot_RoC_sequence(
#' sequence_02_peak,
#' age_threshold = 8e3,
#' Roc_threshold = 2,
#' Peaks = TRUE,
#' trend = "trend_non_linear")
#' }
fc_plot_RoC_sequence <- 
  function(
    data_source,
    age_threshold = FALSE,
    Roc_threshold = FALSE,
    Peaks = FALSE,
    trend = FALSE){
    
    # age_threshold
    assertthat::assert_that(
      age_threshold == FALSE | is.numeric(age_threshold),
      msg = "'age_threshold' must be 'FALSE' or 'numeric'")
    
    if(age_threshold == FALSE){
      age_threshold <- max(data_source$Age)
    } 
    
    data_source_filter <-  
      dplyr::filter(data_source, Age <= age_threshold)
    
    
    #Roc_threshold
    assertthat::assert_that(
      Roc_threshold == FALSE | is.numeric(Roc_threshold),
      msg = "'Roc_threshold' must be 'FALSE' or 'numeric'")
    
    if(Roc_threshold == FALSE){
      Roc_threshold <- max(data_source$ROC_up)
    }
    
    
    p.fin <- 
      ggplot2::ggplot(data_source_filter, mapping =  ggplot2::aes( y=ROC,x= Age)) +
      ggplot2::theme_classic() +
      ggplot2::scale_x_continuous(trans = "reverse") +
      ggplot2::geom_vline(xintercept = seq(0, age_threshold, 2e3), colour = "gray90", size = 0.1)+
      ggplot2::coord_flip(xlim = c(age_threshold, 0), ylim = c(0, Roc_threshold)) +
      ggplot2::geom_ribbon(mapping = ggplot2::aes(ymin=ROC_up, ymax=ROC_dw), fill="gray90") +
      ggplot2::geom_line(alpha=1, size=1, color="gray30") +
      ggplot2::geom_hline(yintercept = 0, color="gray30", lty = 3) +
      ggplot2::labs(x = "Age (cal yr BP)",
                    y = "Rate of change score")
    
    if (Peaks == TRUE){
      
      assertthat::assert_that(
        "Peak" %in% names(data_source),
        msg = "'Peak' are not detected in the data and cannot be plotted.
      Select 'Peaks = FALSE'")
      
      p.fin <- 
        p.fin + 
        ggplot2::geom_point(data = . %>% dplyr::filter(Peak==T),
                            color="green",
                            alpha=1,
                            size=1)
    }
    
    if (trend != F){
      
      assertthat::assert_that(
        trend %in% c("threshold","trend_linear","trend_non_linear"),
        msg = "'trend' method for peak detection have to specified: 'threshold', 'trend_linear','trend_non_linear'"
      )
      
      if (trend == "threshold"){
        p.fin <- 
          p.fin+
          ggplot2::geom_hline(yintercept = stats::median(data_source_filter$ROC),
                              color = "blue", size = 1)
      }
      
      
      if (trend == "trend_linear"){
        p.fin <- 
          p.fin +
          ggplot2::geom_line(
            data=data.frame(ROC = stats::predict.glm(stats::glm(ROC~Age,
                                                                data = data_source_filter,
                                                                family = stats::Gamma()),
                                                     type="response"),
                            Age = data_source_filter$Age),
            color = "blue", size = 1)
      }
      
      if (trend == "trend_non_linear"){
        p.fin <- 
          p.fin +
          ggplot2::geom_line(
            data=data.frame(ROC = mgcv::predict.gam(mgcv::gam(ROC~s(Age, k=3),
                                                              data = data_source_filter,
                                                              family = stats::Gamma(),
                                                              method = "REML"),
                                                    type="response"),
                            Age = data_source_filter$Age),
            color = "blue", size = 1)
      }
    }
    
    return(p.fin)
  }
