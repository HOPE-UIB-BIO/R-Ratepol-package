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

"threshold"         = Threshold value is set for whole dataset (after subsetting for interest_threshold) as median of all RoC values. Peak is  consider significant if 95 percent quantile (gain from randomisation) is higher than threshold

"trend_linear"      = linear model is fitted with RoC and Age. Differences between the model and each point is calculated. SD is calculated from all the differences.Peak is considered significant if it is 1.5 SD higher than the model

"trend_non_linear"  = very crude GAM model (k = 3) is fitted with RoC and Age. Differences between GAM and each point is calculated. SD is calculated from all the differences. Peak is considered significant if it is 1.5 SD higher than GAM
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
    rand = 10e3,
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
    rand = 10e3,
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

