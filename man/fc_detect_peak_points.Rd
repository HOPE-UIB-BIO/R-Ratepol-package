\name{fc_detect_peak_points}
\alias{fc_detect_peak_points}
\title{
Detect significant peak points
}
\description{
Detect points of sudden increae of Rate-of-change values
}
\usage{
fc_detect_peak_points(data_source, method = "trend_linear", sd_threshold = 1.5)
}

\arguments{
  \item{data_source}{
data from fc_estimate_RoC function.
}
  \item{method}{
method of peak detection:
  
"threshold"        = Threshold value is set for whole dataset (after subsetting for interest_threshold) as median of all RoC values. Peak is  consider significant if 95 percent quantile (gain from randomisation) is higher than threshold

"trend_linear"     = linear model is fitted with RoC and Age. Differences between the model and each point is calculated. SD is calculated from all the differences.Peak is considered significat if it is sd_threshold (default = 1.5) SD higher than the model
  
"trend_non_linear" = very crude GAM model (k = 3) is fitted with RoC and Age. Differences between GAM and each point is calculated. SD is calculated from all the differences. Peak is considered significat if it is sd_threshold (default = 1.5) SD higher than GAM
            
"GAM_deriv"        = smooth GAM model is fitted with RoC and Age. First derivative as welll as continuous confidence intervals are calculated from the model. Peak is considered significat if confidence intervals of first derivative is higher than 0.

"SNI"              = Signal-to-Noise Index, following adapted the SNI from Kelly et al. (2011) written to detect changes in charcoal records. The SNI is calculated for the whole RoC sequence and point is consider significant if has SNI value higher than 3 (following suggesting from Kelly et al. (2011)).
}
\item{sd_threshold}{
Threshold that SD of residuals are compared to, to determine peak point (default = 1.5).
}

}
\details{
A rapid change in composition or relative abundances of taxa within the sequence can provide a means of comparing RoC between sequences and interpreting the potential drivers of assemblage change. To detect such significant peak-points of RoC scores in each sequence, each point is tested to see if it represents a significant increase in RoC values. There are various ways to detect peak-points in a time series and R-Ratepol is able detect peak-points using five methods: 

    - Threshold: Each point in the RoC sequence is compared to a median of all RoC scores from the whole sequence (i.e. threshold value). The point is considered significant if the 95th quantile of the RoC scores from all calculations is higher than the threshold value. 
    
    - Linear trend: linear model is fitted between the RoC values and their ages. Differences between the model and each point are calculated (residuals). The standard deviation (SD) is calculated from all the residuals. A peak is considered significant if it is 1.5 SD higher than the model.
    
    - Non-linear trend: A conservative generalised additive model (GAM) is fitted through the RoC scores and their ages (GAM= RoC ~ s(age,k=3) using the mgcv package (Wood, 2011). The distance between each point and the fitted value is calculated (residuals). The standard deviation (SD) is calculated from all the residuals. A peak is considered significant if it is 1.5 SD higher than the model. 
    
    - F-deriv GAM: smooth GAM model is fitted to the RoC scores and their ages (GAM= RoC ~ s(age). The first derivative as well as continuous confidence intervals are calculated from the model using the gratia package (Simpson, 2019). A peak is considered significant if the confidence intervals of the first derivative differ from 0 (for more information see Simpson, 2018).
    
    - Signal-to-noise (SNI) method: We adapted SNI from Kelly et al. (2011), which was developed to detect changes in charcoal stratigraphical records. SNI is calculated for the whole RoC sequence and a peak-point is considered significant if it has an SNI value higher than 3. 
}

\references{
Kelly, R.F., Higuera, P.E., Barrett, C.M., Feng Sheng, H., 2011. A signal-to-noise index to quantify the potential for peak detection in sediment-charcoal records. Quat. Res. 75, 11–17. https://doi.org/10.1016/j.yqres.2010.07.011

Simpson, G.L., 2019. gratia: graceful’ggplot’‐based graphics and other functions for GAMs fitted using “mgcv.” R Packag. version 0.2‐1.

Simpson, G.L., 2018. Modelling palaeoecological time series using generalised additive models. Front. Ecol. Evol. 6, 1–21. https://doi.org/10.3389/fevo.2018.00149

Wood, S.N., 2011. Fast stable restricted maximum likelihood and marginal likelihood estimation of semiparametric generalized linear models. J. R. Stat. Soc. Ser. B Stat. Methodol. 73, 3–36. https://doi.org/10.1111/j.1467-9868.2010.00749.x



}
\author{
Ondrej Mottl
}

\examples{

# Toy example with low number of randomisations
# We higly recommend to increase the randomisations to 10.000 

example_data <-  RRatepol::example_data

sequence_01 <- 
  fc_estimate_RoC(
    data_source_community = example_data$pollen_data[[1]],
    data_source_age = example_data$sample_age[[1]],
    age_uncertainty = FALSE,
    smooth_method = "shep",
    Working_Units = "MW",
    rand = 100, 
    treads = FALSE,
    DC = "chisq")

sequence_01_peak <-
  fc_detect_peak_points(sequence_01, method = "trend_non_linear", sd_threshold = 1.5)
  
fc_plot_RoC_sequence(
  sequence_01_peak,
  age_threshold = 8e3,
  Roc_threshold = 2,
  Peaks = TRUE,
  trend = "trend_non_linear")
}

