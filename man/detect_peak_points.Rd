% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detect_peak_points.R
\name{detect_peak_points}
\alias{detect_peak_points}
\title{Detect significant peak points}
\usage{
detect_peak_points(
  data_source,
  sel_method = c("trend_linear", "trend_non_linear", "threshold", "GAM_deriv", "SNI"),
  sd_threshold = 2
)
}
\arguments{
\item{data_source}{Data.frame. Output of \code{estimate_roc} function}

\item{sel_method}{Character. A method to use for peak-poit detection:
\itemize{
\item \code{"threshold"} - Each point in the RoC sequence is compared to a median
of all RoC scores from the whole sequence (i.e. threshold value). The RoC
value for a point is considered significant if the 95th quantile of the RoC
scores from all calculations is higher than the threshold value.
\item \code{"trend_linear"} - A linear model is fitted between the RoC values and
their ages. Differences between the model and each point are calculated (residuals).
The standard deviation (SD) is calculated from all the residuals. A peak is considered
significant if it is 2 SD higher than the model (\code{sd_threshold} = 2).
\item \code{"trend_non_linear"} - A conservative generalised additive model (GAM)
is fitted through the RoC scores and their ages (GAM = \code{RoC ~ s(age, k = 3)})
using the \code{mgcv} package (Wood, 2011). The distance between each point and
the fitted value is calculated (residuals). The standard deviation (SD) is
calculated from all the residuals. A peak is considered significant if it
is 2 SD higher than the model (\code{sd_threshold} = 2).
}}

\item{sd_threshold}{Numeric. Threshold that SD of residuals are compared to, to determine
peak-point (default = 2)}
}
\description{
Detect points of sudden increase of Rate-of-change values
}
\details{
A rapid change in composition or relative abundances of variables within the
sequence can provide a means of comparing RoC between sequences and interpreting
the potential drivers of assemblage change. To detect such significant peak-points
of RoC scores in each sequence, each point is tested to see if it represents
a significant increase in RoC values. There are various ways to detect
peak-points in a time series and R-Ratepol is able to detect peak-points
using five methods:
\itemize{
\item Threshold (\code{sel_method} = \code{"threshold"}) - Each point in the RoC sequence is
compared to a median of all RoC scores from the whole sequence (i.e. threshold value).
The ROC value for a point is considered significant if the 95th quantile of
the RoC scores from all calculations is higher than the threshold value.
\item Linear trend (\code{sel_method} = \code{"trend_linear"}) - A linear model is fitted
between the RoC values and their ages. Differences between the model and each
point are calculated (residuals). The standard deviation (SD) is calculated
from all the residuals. A peak is considered significant if it is 2 SD higher
than the model (\code{sd_threshold} = 2).
\item Non-linear trend (\code{sel_method} = \code{"trend_non_linear"}) - A conservative
generalised additive model (GAM) is fitted through the RoC scores and their
ages (GAM = \code{RoC ~ s(age, k = 3)} using the \code{mgcv} package (Wood, 2011).
The distance between each point and the fitted value is calculated (residuals).
The standard deviation (SD) is calculated from all the residuals. A peak is
considered significant if it is 2 SD higher than the model (\code{sd_threshold} = 2).
\item F-deriv GAM  (\code{sel_method} = \code{"GAM_deriv"}) - A smooth GAM model is fitted
to the RoC scores and their ages (GAM = \code{RoC ~ s(age)}). The first derivative
as well as continuous confidence intervals are calculated from the model
using the \code{gratia} package (Simpson, 2019). A peak is considered significant
if the confidence intervals of the first derivative differ from 0
(for more information see Simpson, 2018).
\item Signal-to-noise method (\code{sel_method} = \code{"SNI"}) - We adapted SNI from
Kelly et al. (2011), which was developed to detect changes in charcoal
stratigraphical records. SNI is calculated for the whole RoC sequence and
a peak-point is considered significant if it has an SNI value higher than 3.
}
}
\examples{
\dontrun{
example_data <- RRatepol::example_data

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

sequence_01_peak <-
  detect_peak_points(
    sequence_01,
    sel_method = "trend_non_linear",
    sd_threshold = 2
  )

plot_roc(
  sequence_01_peak,
  age_threshold = 8e3,
  roc_threshold = 2,
  peaks = TRUE,
  trend = "trend_non_linear"
)
}
}
\references{
Kelly, R.F., Higuera, P.E., Barrett, C.M., Feng Sheng, H., 2011. A signal-to-noise
index to quantify the potential for peak detection in sediment-charcoal records.
Quat. Res. 75, 11-17.

Simpson, G.L., 2019. gratia: graceful'ggplot'-based graphics and other functions
for GAMs fitted using 'mgcv' R Packag. version 0.2-1.

Simpson, G.L., 2018. Modelling palaeoecological time series using generalised
additive models. Front. Ecol. Evol. 6, 1-21.

Wood, S.N., 2011. Fast stable restricted maximum likelihood and marginal
likelihood estimation of semiparametric generalized linear models.
J. R. Stat. Soc. Ser. B Stat. Methodol. 73, 3-36.
}
