fc_plot_RoC_sequence <- function (data_source, age_treshold = 15000, Roc_threshold=2, Peaks = T, method = "GAM"){
  
  data_source_filter <-  dplyr::filter(data_source, Age <= age_treshold)
  
  p.fin <- 
    ggplot2::ggplot(data_source_filter, mapping =  ggplot2::aes( y=ROC,x= Age)) +
    ggplot2::theme_classic() +
    ggplot2::scale_x_continuous(trans = "reverse") +
    ggplot2::geom_vline(xintercept = seq(0, age_treshold, 2e3), colour = "gray90", size = 0.1)+
    ggplot2::coord_flip(xlim = c(age_treshold, 0), ylim = c(0, Roc_threshold)) +
    ggplot2::geom_ribbon(mapping = ggplot2::aes(ymin=ROC_up, ymax=ROC_dw), fill="gray90") +
    ggplot2::geom_line(alpha=1, size=1, color="gray30") +
    ggplot2::geom_hline(yintercept = 0, color="gray30", lty = 3) +
    ggplot2::labs(x = "Age (cal yr BP)",
                  y = "Rate of change score")

  if (Peaks == T){
    p.fin <- 
      p.fin + 
      ggplot2::geom_point(data = . %>% dplyr::filter(Peak==T),
                          color="green",
                          alpha=1,
                          size=1)
    
    if (! method %in% c("GAM","Threshold")) 
      stop("Used method for peak detection have to specified 'GAM' / 'Threshold' ")
    
    if (method == "GAM"){
      p.fin <- 
        p.fin +
        ggplot2::geom_line(
          data=data.frame(ROC = mgcv::predict.gam(mgcv::gam(ROC~s(Age, k=3),
                                                            data = data_source_filter,
                                                            family = "Gamma()",
                                                            method = "REML"),
                                                  type="response"),
                                         Age = data_source_filter$Age),
                         color = "blue", size = 1)
    }
    
    if (method == "Threshold"){
      p.fin <- 
        p.fin+
        ggplot2::geom_hline(yintercept = median(data_source_filter$ROC), color = "blue", size = 1)
    }
    
  }

  return(p.fin)
}
