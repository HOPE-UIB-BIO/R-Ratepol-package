\name{fc_plot_RoC_sequence}
\alias{fc_plot_RoC_sequence}

\title{
Function to plot the Rate-of-Change sequence.
}
\description{
Plot Rate-of-Change sequence with a error estimate and trend and/or peak-points in present.
}
\usage{
fc_plot_RoC_sequence(
  data_source,
  age_threshold = FALSE,
  Roc_threshold = FALSE,
  Peaks = FALSE,
  trend = FALSE)
}

\arguments{
  \item{data_source}{
Output of fc_estimate_RoC (tibble)
}
  \item{age_threshold}{
Cut-off value used as maximum age
}
  \item{Roc_threshold}{
Cut-off value used as maximum RoC value
}
  \item{Peaks}{
If peak-points are presented in the dataset and \code{Peaks} == \code{TRUE}, then peak points will be displayed  
}
  \item{trend}{
If peak-points are presented in the dataset and \code{Peaks} == \code{TRUE}, then one of the three method can be used to visualise the process of peak detection:

\code{"threshold"}         = Each point in the RoC sequence is compared to a median of all RoC scores from the whole sequence (i.e. threshold value). The ROC value for a point is considered significant if the 95th quantile of the RoC scores from all calculations is higher than the threshold value.

\code{"trend_linear"}      = A linear model is fitted between the RoC values and their ages. Differences between the model and each point are calculated (residuals). The standard deviation (SD) is calculated from all the residuals. A peak is considered significant if it is 1.5 SD higher than the model.
  
\code{"trend_non_linear"}  = A conservative generalised additive model (GAM) is fitted through the RoC scores and their ages (GAM = RoC ~ s(age,k=3) using the mgcv package (Wood, 2011). The distance between each point and the fitted value is calculated (residuals). The standard deviation (SD) is calculated from all the residuals. A peak is considered significant if it is 1.5 SD higher than the model. 
}
}

\author{
Ondrej Mottl
}

\examples{

\dontrun{
example_data <-  RRatepol::example_data

# example 1
sequence_01 <- 
  fc_estimate_RoC(
    data_source_community = example_data$pollen_data[[1]],
    data_source_age = example_data$sample_age[[1]],
    age_uncertainty = FALSE,
    smooth_method = "shep",
    Working_Units = "MW",
    rand = 1e3,
    treads = TRUE,
    DC = "chisq")

fc_plot_RoC_sequence(
  sequence_01,
  age_threshold = 8e3,
  Roc_threshold = 1)

# example 2
sequence_02 <- 
  fc_estimate_RoC(
    data_source_community = example_data$pollen_data[[2]],
    data_source_age = example_data$sample_age[[2]],
    age_uncertainty = FALSE,
    smooth_method = "shep",
    Working_Units = "MW",
    rand = 1e3,
    treads = TRUE,
    DC = "chisq")

sequence_02_peak <-
  fc_detect_peak_points(sequence_01, method = "trend_non_linear")
  
fc_plot_RoC_sequence(
  sequence_02_peak,
  age_threshold = 8e3,
  Roc_threshold = 2,
  Peaks = TRUE,
  trend = "trend_non_linear")
}

}

