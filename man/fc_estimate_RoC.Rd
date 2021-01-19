\name{fc_estimate_RoC}
\alias{fc_estimate_RoC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
RRatepol: Estimate Rate of change
%%  ~~function to do ... ~~
}
\description{
A function to estimate Rate of change in community data in time series 
}
\usage{
fc_estimate_RoC(data_source_community,
                data_source_age,
                age_uncertainty = FALSE,
                smooth_method = "shep",
                smooth_N_points = 5,
                smooth_age_range = 500,
                smooth_N_max = 9,
                Working_Units = "MW",
                bin_size = 500,
                Number_of_shifts = 5,
                rand = 1,
                treads = TRUE,
                standardise = FALSE,
                N_individuals = 150,
                tranform_to_proportions = TRUE,
                DC = "chord",
                interest_threshold = FALSE,
                Debug = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_source_community}{
Community data with species as collumns and levels (samples) as rows. First column should be sample.id. (character)
}
  \item{data_source_age}{
Dataframe with two columns:
sample.id = unique ID of each level (character)
age = age of level (numeric)
}
 \item{age_uncertainty}{
Usage of age uncertainty form Age-depth models. Either:

[1] matrix with number of columns as number of samples. Each column is one sample, each row is one age sequence from age-depth model. Age sequence is randomly sampled from age-depth model uncertainties at the begining of each run.

[2] FALSE. Age uncertainties are not available and, therefore, will not be used.
}

 \item{smooth_method}{
type of smoothing applied for the each of the pollen type

"none"    = None: Pollen data is not smoothed

"m.avg"   = Moving average:

"grim"    = 	Grimm’s smoothing

"age.w"   = 	Age-weighted average

"shep"    =   Shepard's 5-term filter
}
 \item{smooth_N_points}{
Number of points for used for moving average, Grimm and Age-Weighted smoothing (odd number)
}
 \item{smooth_N_max}{
Maximal number of samples to look in Grimm smoothing
}
 \item{smooth_age_range}{
Maximal age range for both Grimm and Age-weight smoothing
}
 \item{Working_Units}{
Selection of units that the DC will be calculated between.

"levels" = individual levels are ogoing to be used

"bins"   = samples in predifined bins will be pooled together and one random sample will be selected from each time bin as a representation.

"MW"     = Bins of selected size are created, starting from zero. This is repeated many times, with each time bin (window) shifting by Z years forward. This is repeated X times, where X = bin size / Z. Sample, whitch is the closest to the begining of the bin is selected asa representation of the bin.  

}
 \item{bin_size}{
bin_size = size of the time bin (in years)
}
 \item{Number_of_shifts}{
Value determining the number of shifts of window used in Moving window method
}
 \item{rand}{
Number of runs used in randomisation.
}
 \item{treads}{
Preference of usage of parallel computation of randomisation
[value] = selected number of cores
T = automaticly selected number of cores
F = does not use parallel computation
}

 \item{standardise}{
T = standardise each Working Unit to cetrain number of individuals (using random resampling without repetition)
}
 \item{N_individuals}{
Number of grain to perform standardisation to. The N_individuals is automatically adjusted to the smallest  number of pollen grains in sequence. 
}
 \item{tranform_to_proportions}{
Should the community data be transformed to a proportion during calculations?
}
 \item{DC}{
Disimilarity coeficient. Type of calculation of differences between Working Units

"euc"     = 	Euclidean distance

"euc.sd"  = 	Standardised Euclidean distance

"chord"   = 	Chord distance

"chisq    = 	Chi-squared coefficient

"gower"   =   Gower's distance
}
 \item{interest_threshold}{
[value] = age, after which all results of RoC are excluded before detection of peak points 
}
 \item{Debug}{
T = show messages from internal processes
}

}
\details{
The general process of computation RoC in R-Ratepol follow simple
sequence::

1.  Both community and age data are extracted and compiled together.

2.  (optional) Community data is smoothed. Each taxon is smoothed using
    one of five in-build smoothing methods: none, Shepard’s 5-term
    filter (Davis, 1986; Wilkinson, 2005)., moving average, age-weighted
    average, Grimm’s smoothing (Grimm and Jacobson, 1992).

3.  Working Units (WU) for computation are selected.

4.  Single run (an individual loop) is computed:

    -   (optional) A single age sequence is randomly select from age
        uncertainties for all levels.

    -   (optional) Community data is standardised, i.e. community data
        is subsampled in each WU to a selected total count of
        individuals (e.g. 150 pollen grains).

    -   RoC between adjacent WUs is calculated as the dissimilarity
        coefficient (DC) standardised by age differences between WUs.
        Four in-build dissimilarity coefficients: Euclidean distance,
        standardised Euclidean distance, Chord distance, Chi-squared
        coefficient (Prentice, 1980), Gower’s distance (Gower, 1971).

5.  Single repetition (a loop) is repeated multiple times (e.g. 10,000
    times).

6.  Results from all runs are summarised.

7.  Significant peak-points are detected and validated. Five in-build
    methods of testing for the significance of individual RoC points:
    Threshold, Linear trend, Non-linear trend, first derivative of
    generalised additive model (f-deriv GAM; Simpson, 2018), and
    Signal-to-Noise Index (SNI; Kelly et al., 2011).

Selection of working units

RoC is calculated between consecutive Working Units (WU). Traditionally,
these WUs represent individual stratigraphical levels. However, changes
in sedimentation rates and sampling strategies can result in an uneven
distribution of levels within a time sequence, which in turn makes the
comparison of RoC between sequences problematic. There are various
methods that attempt to minimise such problems. The first is
interpolation of levels to evenly spaced time intervals, and the use of
the interpolated data as WUs. This can lead to a loss of information
when the density of levels is high. Second is binning of levels: pollen
data are pooled into age brackets of various size (bins) and these serve
as WUs. Here, the issue is a lower resolution of WUs and their uneven
size in terms of total pollen count (bins with more levels have higher
pollen counts). Third is selective binning: like classical binning, bins
of selected size are created, but instead of pooling pollen data
together, only one level is selected as representative of each bin. This
results in an even number of WUs in bins with a similar number of pollen
grains. However, the issue of low resolution remains. We propose a new
method of binning with a moving window, which is a compromise between
using individual levels and selective binning. This method follows a
simple sequence: bins are created, levels are selected as in selective
binning, and RoC between bins is calculated. However, the brackets of
the time bin (window) are then moved forward by a selected amount of
time (Z), levels are selected again, and RoC calculated for the new set
of WUs. This is repeated X times (where X is the bin size divided by Z)
while retaining all the results. 

Randomisation

Due to the inherent statistical errors in a community datasets
(e.g. pollen count in each level; Birks and Gordon, 1985) and
uncertainties in the age estimates from age-depth modelling, R-Ratepol
can be run several times and the results summarised. Therefore, two
optional settings can be used: usage of age uncertainties and community
data standardisation.

Community data standardisation

Taxa in community dataset can standardised to a certain abundance
(e.g. number of pollen grains in each WU) by rarefaction. Random
sampling without replacement is used to a draw selected number of
individuals in each WU (e.g. 150 pollen).

Age uncertainties

For each run, a single age sequence from age uncertainties is randomly
selected. The calculation between two consecutive WUs (i.e. one
working-unit combination) results in a RoC score and a time position
(which is calculated as the mean age position of the two WUs). However,
due to random sampling of the age sequence, each WU combination will
result in multiple RoC values and age positions. R Ratepol assigns the
age position of each WU combination as the median time position from all
calculations. The final RoC value for a single WU combination is
calculated as the median of the scores from all randomisations. In
addition, 95th quantile from all randomisations is calculated as error
estimate.

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Birks, H.J.B., Gordon, A.D., 1985. Numerical Methods in Quaternary Pollen Analysis. Academic Press, London.

Davis, J.C., 1986. Statistics and Data Analysis in Geology, 2nd edn. ed. J. Wiley & Sons, New York.

Gower, J.C., 1971. A general coefficient of similarity and some of its properties. Biometrics 27, 857–871.

Grimm, E.C., Jacobson, G.L., 1992. Fossil-pollen evidence for abrupt climate changes during the past 18000 years in eastern North America. Clim. Dyn. 6, 179–184.

Kelly, R.F., Higuera, P.E., Barrett, C.M., Feng Sheng, H., 2011. A signal-to-noise index to quantify the potential for peak detection in sediment-charcoal records. Quat. Res. 75, 11–17. https://doi.org/10.1016/j.yqres.2010.07.011

Prentice, I.C., 1980. Multidimensional scaling as a research tool in Quaternary palynology: A review of theory and methods. Rev. Palaeobot. Palynol. 31, 71–104. https://doi.org/10.1016/0034-6667(80)90023-8

Simpson, G.L., 2018. Modelling palaeoecological time series using generalised additive models. Front. Ecol. Evol. 6, 1–21. https://doi.org/10.3389/fevo.2018.00149
Wilkinson, L., 2005. The Grammar of Graphics. Springer-Verlag, New York, USA 37. https://doi.org/10.2307/2669493
}
\author{
Ondřej Mottl
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
sequence_01 <- 
  fc_estimate_RoC(
    data_source_community = example_data$pollen_data[[1]],
    data_source_age = example_data$sample_age[[1]])

fc_plot_RoC_sequence(
  sequence_01,
  age_treshold = 8e3,
  Roc_threshold = 1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
