% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_roc.R
\name{plot_roc}
\alias{plot_roc}
\title{Function to plot the Rate-of-Change sequence}
\usage{
plot_roc(
  data_source,
  age_threshold = NULL,
  roc_threshold = NULL,
  peaks = FALSE,
  trend = NULL
)
}
\arguments{
\item{data_source}{Data.frame. Output of \code{estimate_roc} function}

\item{age_threshold}{Numeric. Cut-off value used as maximum age.}

\item{roc_threshold}{Numeric Cut-off value used as maximum RoC value.}

\item{peaks}{Logical. If peak-points are presented in the dataset and \code{peaks} == \code{TRUE},
then peak points will be displayed}

\item{trend}{If peak-points are presented in the dataset and \code{peaks} == \code{TRUE},
then one of the three method can be used to visualise the process of peak detection:
\itemize{
\item \code{"threshold"} - Each point in the RoC sequence is compared to a median
of all RoC scores from the whole sequence (i.e. threshold value).
The ROC value for a point is considered significant if the 95th quantile of
the RoC scores from all calculations is higher than the threshold value.
\item \code{"trend_linear"} -  A linear model is fitted between the RoC values
and their ages. Differences between the model and each point are calculated
(residuals). The standard deviation (SD) is calculated from all the residuals.
A peak is considered significant if it is 2 SD higher than the model.
\item \code{"trend_non_linear"} - A conservative generalised additive model (GAM)
is fitted through the RoC scores and their ages (GAM = \code{RoC ~ s(age, k = 3)} using
the \code{mgcv} package (Wood, 2011). The distance between each point and the
fitted value is calculated (residuals). The standard deviation (SD) is
calculated from all the residuals. A peak is considered significant if it
is 2 SD higher than the model.
}}
}
\description{
Plot Rate-of-Change sequence with a error estimate and trend
and/or peak-points in present.
}
\examples{
\dontrun{

example_data <- RRatepol::example_data

# example 1
sequence_01 <-
  estimate_roc(
    data_source_community = example_data$pollen_data[[1]],
    data_source_age = example_data$sample_age[[1]],
    age_uncertainty = FALSE,
    smooth_method = "shep",
    working_units = "MW",
    rand = 1e3,
    treads = TRUE,
    dissimilarity_coefficient = "chisq"
  )

plot_roc(
  sequence_01,
  age_threshold = 8e3,
  roc_threshold = 1
)

# example 2
sequence_02 <-
  estimate_roc(
    data_source_community = example_data$pollen_data[[2]],
    data_source_age = example_data$sample_age[[2]],
    age_uncertainty = FALSE,
    smooth_method = "shep",
    working_units = "MW",
    rand = 1e3,
    treads = TRUE,
    dissimilarity_coefficient = "chisq"
  )

sequence_02_peak <-
  detect_peak_points(sequence_01, sel_method = "trend_non_linear")

plot_roc(
  sequence_02_peak,
  age_threshold = 8e3,
  roc_threshold = 2,
  peaks = TRUE,
  trend = "trend_non_linear"
)
}
}
