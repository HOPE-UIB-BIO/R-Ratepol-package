fc_plot_RoC_sequence <- function (data_source, age_treshold = 15000, Roc_threshold=2, Peaks = T, method = "GAM")
{
  data_source_filter <- dplyr::filter(data_source, AGE <= age_treshold)
  
  p.fin <- ggplot2::ggplot(data_source_filter, mapping =  ggplot2::aes( y=ROC,x= AGE)) +
    ggplot2::theme_classic() +
    ggplot2::scale_x_continuous(trans = "reverse") +
    ggplot2::geom_vline(xintercept = seq(0,age_treshold, 2e3), colour = "gray90", size = 0.1)+
    ggplot2::coord_flip(xlim=c(age_treshold,0), ylim = c(0,Roc_threshold)) +
    ggplot2::geom_ribbon(mapping = ggplot2::aes(ymin=ROC.up, ymax=ROC.dw), fill="gray90") +
    ggplot2::geom_line(alpha=1, size=1, color="gray30") +
    ggplot2::geom_hline(yintercept = 0, color="red") +
    ggplot2::labs(x = "Age (cal yr BP)",
                  y = "Rate of change score")

  if (Peaks == T){
    p.fin <- p.fin+ 
      ggplot2::geom_point(data = . %>% dplyr::filter(PEAK==T),color="green", alpha=1, size=1)
    
    if (! method %in% c("GAM","Threshold")) 
      stop("Used method for peak detection have to specified 'GAM' / 'Threshold' ")
    
    if (method == "GAM"){
      p.fin <- p.fin+
      ggplot2::geom_line(data=data.frame(ROC = mgcv::predict.gam(mgcv::gam(ROC~s(AGE, k=3), data = data_source_filter,
                                                                           family = "Gamma()",
                                                                           method = "REML"), type="response"),
                                         AGE = data_source_filter$AGE),
                         color="blue", size=1)
    }
    
    if (method == "Threshold"){
      p.fin <- p.fin
      ggplot2::geom_hline(yintercept = median(data_source_filter$ROC), color="blue", size=1)
    }
    
  }

  return(p.fin)
}
