% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_iteration.R
\name{run_iteration}
\alias{run_iteration}
\title{Run a single interation of RoC estimation}
\usage{
run_iteration(
  data_source_run,
  bin_selection = "first",
  standardise = FALSE,
  n_individuals = 150,
  tranform_to_proportions = TRUE,
  dissimilarity_coefficient = "euc",
  time_standardisation = 500,
  verbose = FALSE
)
}
\arguments{
\item{data_source_run}{List with \code{data} and \code{bins} prepared by \code{prepare_data}}

\item{bin_selection}{Character. Setting determining the the process of selection of
samples from bins.
\itemize{
\item \code{"first"} - sample closest to the beginning of the bin is selected
as a representation.
\item \code{"random"} - a random sample is selected as a representation.
}}

\item{standardise}{Logical. If \code{standardise} == \code{TRUE}, then standardise
each Working Unit to certain number of individuals (using random resampling
without repetition)}

\item{n_individuals}{Numeric. Number of grain to perform standardisation to.
The \code{N_individual} is automatically adjusted to the smallest number
of pollen grains in sequence.}

\item{tranform_to_proportions}{Logical. Should the community data be transformed to a
proportion during calculations?}

\item{dissimilarity_coefficient}{Character. Dissimilarity coefficient. Type of calculation of differences
between Working Units. See \code{vegan::vegdist} for more details.
\itemize{
\item \code{"euc"} - Euclidean distance
\item \code{"euc.sd"} - Standardised Euclidean distance
\item \code{"chord"} - Chord distance
\item \code{"chisq"} - Chi-squared coefficient
\item \code{"gower"} - Gower's distance
\item \code{"bray"} - Bray-Curtis distance
}}

\item{time_standardisation}{Numeric. Units scaling for result RoC values. For example,
if \code{time_standardisation} = 100, the RoC will be reported as
dissimilarity per 100 yr.}

\item{verbose}{Logical. If \code{TRUE}, function will output messages about internal processes}
}
\description{
A single run is computed following the simple steps:
\itemize{
\item Subsetting levels in each bin: Here the working units (WU) are defined
\item Standardisation of assemblage data in each WU
\item Calculation of calculated as the dissimilarity coefficient (dissimilarity_coefficient)
\item Calculation of RoC between WUs: RoC is calculated as (dissimilarity_coefficient)
standardised by age differences between WUs.
}
}
\seealso{
\code{\link[=estimate_roc]{estimate_roc()}}
}
\keyword{internal}
