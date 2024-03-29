% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_roc.R
\name{estimate_roc}
\alias{estimate_roc}
\title{RRatepol: Estimate rate of change}
\usage{
estimate_roc(
  data_source_community,
  data_source_age,
  age_uncertainty = NULL,
  smooth_method = c("none", "m.avg", "grim", "age.w", "shep"),
  smooth_n_points = 5,
  smooth_age_range = 500,
  smooth_n_max = 9,
  working_units = c("levels", "bins", "MW"),
  bin_size = 500,
  number_of_shifts = 5,
  bin_selection = c("random", "first"),
  standardise = FALSE,
  n_individuals = 150,
  dissimilarity_coefficient = c("euc", "euc.sd", "chord", "chisq", "gower", "bray"),
  tranform_to_proportions = TRUE,
  rand = NULL,
  use_parallel = FALSE,
  interest_threshold = NULL,
  time_standardisation = NULL,
  verbose = FALSE
)
}
\arguments{
\item{data_source_community}{Data.frame. Community data with species as columns and
levels (samples) as rows. First column should be \code{sample_id} (character).}

\item{data_source_age}{Data.frame with two columns:
\itemize{
\item \code{sample_id} - unique ID of each level (character)
\item \code{age} - age of level (numeric)
}}

\item{age_uncertainty}{Usage of age uncertainty form Age-depth models. Either:
\itemize{
\item matrix with number of columns as number of samples. Each column is one sample,
each row is one age sequence from age-depth model. Age sequence is randomly
sampled from age-depth model uncertainties at the beginning of each run.
\item \code{NULL} - Age uncertainties are not available and, therefore, will not be used.
}}

\item{smooth_method}{Character. type of smoothing applied for the each of the pollen type
\itemize{
\item \code{"none"} - Pollen data is not smoothed
\item \code{"m.avg"} - Moving average
\item \code{"grim"} - Grimm's smoothing
\item \verb{"age.w""} - Age-weighted average
\item \code{"shep"} - Shepard's 5-term filter
}}

\item{smooth_n_points}{Numeric. Number of points for used for moving average,
Grimm and Age-Weighted smoothing (odd number)}

\item{smooth_age_range}{Numeric. Maximal age range for both Grimm and Age-weight smoothing}

\item{smooth_n_max}{Numeric. Maximal number of samples to look in Grimm smoothing}

\item{working_units}{Character. Selection of units that the dissimilarity_coefficient will be calculated between.
\itemize{
\item \code{"levels"} - individual levels are going to be used
\item \code{"bins"} - samples in predefined bins will be pooled together and one sample
will be selected from each time bin as a representation.
\item \code{"MW"} - Bins of selected size are created, starting from the beginning of the core.
This is repeated many times, with each time bin (window) shifting by Z years forward.
This is repeated X times, where X = bin size / Z.
}}

\item{bin_size}{Numeric. Size of the time bin (in years)}

\item{number_of_shifts}{Numeric. Value determining the number of shifts of window used
in Moving window method}

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

\item{tranform_to_proportions}{Logical. Should the community data be transformed to a
proportion during calculations?}

\item{rand}{Numeric. Number of runs used in randomisation.}

\item{use_parallel}{Preference of usage of parallel computation of randomisation
\itemize{
\item \verb{[value]} - selected number of cores
\item \code{TRUE} - automatically selected number of cores
\item \code{FALSE} - does not use parallel computation (only single core)
}}

\item{interest_threshold}{Numeric. Optional. Age, after which all results of RoC are excluded.}

\item{time_standardisation}{Numeric. Units scaling for result RoC values. For example,
if \code{time_standardisation} = 100, the RoC will be reported as
dissimilarity per 100 yr.}

\item{verbose}{Logical. If \code{TRUE}, function will output messages about internal processes}
}
\description{
A function to estimate Rate of change in community data in time series
}
\details{
R-Ratepol is written as an R package and includes a range of
possible settings including a novel method to evaluate RoC in a single
stratigraphical sequence using assemblage data and age uncertainties for
each level. There are multiple built-in dissimilarity coefficients (dissimilarity_coefficient) for
different types of assemblage data, and various levels of data smoothing
that can be applied depending on the type and variance of the data.
In addition, R-Ratepol can use randomisation, accompanied by use of age
uncertainties of each level and taxon standardisation to detect RoC patterns
in datasets with high data noise or variability (i.e. numerous rapid changes
in composition or sedimentation rates).

The computation of RoC in R-Ratepol is performed using the following steps:
\enumerate{
\item Assemblage and age-model data are extracted from the original source and
should be compiled together, i.e. depth, age, variable (taxon) 1, variable (taxon) 2, etc.
\item (optional) Smoothing of assemblage data: Each variable within the
assemblage data is smoothed using one of five in-built smoothing methods:
\itemize{
\item none (\code{smooth_method} = \code{"none"})
\item Shepard's 5-term filter (\code{smooth_method} = \code{"shep"}; Davis, 1986; Wilkinson, 2005)
\item moving average (\code{smooth_method} = \verb{"m.avg"\}})
\item age-weighted average (\code{smooth_method} = \code{"age.w"})
\item Grimm's smoothing (\code{smooth_method} = \code{"grim"}; Grimm and Jacobson, 1992)
}
\item Creation of time bins: A template for all time bins in all window movements is created.
\item A single run (an individual loop) is computed:
\itemize{
\item (optional) Selection of one time series from age uncertainties (see section on randomisation)
\item Subsetting levels in each bin: Here the working units (WU) are defined
\item (optional) Standardisation of assemblage data in each WU
\item The summary of a single run is produced based on all moving windows
\item Calculation of RoC between WUs: RoC is calculated as the dissimilarity
coefficient (dissimilarity_coefficient) standardised by age differences between WUs. Five in-built
dissimilarity coefficients are available:
\itemize{
\item Euclidean distance (\code{dissimilarity_coefficient} = \code{"euc"})
\item standardised Euclidean distance (\code{dissimilarity_coefficient} = \code{"euc.sd"})
\item Chord distance (\code{dissimilarity_coefficient} = \code{"chord"})
\item Chi-squared coefficient (\code{dissimilarity_coefficient} = \code{"chisq"})
\item Gower's distance (\code{dissimilarity_coefficient} = \code{"gower"})
\item Bray-Curtis distance (\code{dissimilarity_coefficient} = \code{"bray"})
}
}
\item Step 4 is repeated multiple times (e.g. 10,000 times).
\item Validation and summary of results from all runs of RoC calculation are produced.
\item (Optional) Data beyond a certain age can be excluded.
}
\subsection{Selection of working units (WU; Step 3)}{

RoC is calculated between consecutive Working Units (WU). Traditionally,
these WUs represent individual stratigraphical levels. However, changes in
sedimentation rates and sampling strategies can result in an uneven temporal
distribution of levels within a time sequence, which in turn makes
the comparison of RoC between sequences problematic. There are various methods
that attempt to minimise such problems. The first is interpolation of levels
to evenly spaced time intervals, and the use of the interpolated data as WUs.
This can lead to a loss of information when the density of levels is high.
Second is binning of levels: assemblage data are pooled into age brackets
of various size (i.e. time bins) and these serve as WUs. Here, the issue
is a lower resolution of WUs and their uneven size in terms of total
assemblage count (bins with more levels have higher assemblage counts).
Third is selective binning: like classical binning, bins of selected size
are created, but instead of pooling assemblage data together, only one
level per time bin is selected as representative of each bin. This results
in an even number of WUs in bins with a similar count size in the assemblage.
However, the issue of low resolution remains.
Therefore, we propose a new method of binning with a moving window,
which is a compromise between using individual levels and selective binning.
This method follows a simple sequence: time bins are created,
levels are selected as in selective binning, and RoC between bins is calculated.
However, the brackets of the time bin (window) are then moved forward by a
selected amount of time (Z), levels are selected again (subset into bins),
and RoC calculated for the new set of WUs. This is repeated X times
(where X is the bin size divided by Z) while retaining all the results.

R-Ratepol currently provides several options for selecting WU, namely as i
ndividual levels (\code{working_units} = \code{"levels"}), selective binning of levels
(\code{working_units} = \code{"bins"}), and our new method of binning with a moving
window (\code{working_units} = \code{"MW"})
}

\subsection{Randomisation}{

Due to the inherent statistical errors in uncertainties in the age estimates
from age-depth and the assemblage datasets (e.g. pollen counts in each level;
Birks and Gordon, 1985), R-Ratepol can be run several times and the results
summarised (Steps 5-6). Therefore, two optional settings are available by
using age uncertainties and assemblage data standardisation.
}

\subsection{Age uncertainties}{

For each run, a single age sequence from the age uncertainties is randomly
selected. The calculation between two consecutive WUs (i.e. one working-unit
combination) results in a RoC score and a time position (which is calculated
as the mean age position of the two WUs). However, due to random sampling
of the age sequence, each WU combination will result in multiple RoC values.
The final RoC value for a single WU combination is calculated as the median
of the scores from all randomisations. In addition, the 95th quantile from all
randomisations is calculated as an error estimate.
}

\subsection{Data standardisation (Step 4b)}{

Taxa in the assemblage dataset can be standardised to a certain count
(e.g. number of pollen grains in each WU) by rarefaction. Random sampling
without replacement is used to draw a selected number of individuals from
each WU (e.g. 150 pollen grains).
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
    use_parallel = TRUE,
    dissimilarity_coefficient = "chisq"
  )

plot_roc(
  sequence_01,
  age_threshold = 8e3,
  roc_threshold = 1
)
}
}
\references{
Birks, H.J.B., Gordon, A.D., 1985. Numerical Methods in Quaternary Pollen
Analysis. Academic Press, London.

Davis, J.C., 1986. Statistics and Data Analysis in Geology, 2nd edn. ed.
J. Wiley & Sons, New York.

Grimm, E.C., Jacobson, G.L., 1992. Fossil-pollen evidence for abrupt
climate changes during the past 18000 years in eastern North America.
Clim. Dyn. 6, 179-184.

Wilkinson, L., 2005. The Grammar of Graphics. Springer-Verlag, New York,
USA 37.
}
