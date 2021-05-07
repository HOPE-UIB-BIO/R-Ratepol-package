\name{fc_estimate_RoC}
\alias{fc_estimate_RoC}
\title{
RRatepol: Estimate Rate of change
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
                DC = "chisq",
                interest_threshold = FALSE,
                Debug = FALSE)
}

\arguments{
  \item{data_source_community}{
Community data with species as columns and levels (samples) as rows. First column should be sample.id. (character)
}
  \item{data_source_age}{
Dataframe with two columns:
sample.id = unique ID of each level (character)
age = age of level (numeric)
}
 \item{age_uncertainty}{
Usage of age uncertainty form Age-depth models. Either:

[1] matrix with number of columns as number of samples. Each column is one sample, each row is one age sequence from age-depth model. Age sequence is randomly sampled from age-depth model uncertainties at the beginning of each run.

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

"levels" = individual levels are going to be used

"bins"   = samples in predefined bins will be pooled together and one random sample will be selected from each time bin as a representation.

"MW"     = Bins of selected size are created, starting from zero. This is repeated many times, with each time bin (window) shifting by Z years forward. This is repeated X times, where X = bin size / Z. Sample, whitch is the closest to the beginning of the bin is selected asa representation of the bin.  

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
T = automatically selected number of cores
F = does not use parallel computation
}

 \item{standardise}{
T = standardise each Working Unit to certrain number of individuals (using random resampling without repetition)
}
 \item{N_individuals}{
Number of grain to perform standardisation to. The N_individuals is automatically adjusted to the smallest  number of pollen grains in sequence. 
}
 \item{tranform_to_proportions}{
Should the community data be transformed to a proportion during calculations?
}
 \item{DC}{
Dissimilarity coefficient. Type of calculation of differences between Working Units

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
R-Ratepol is written as an R package and includes a range of possible settings including a novel method to evaluate RoC in a single stratigraphical sequence using assemblage data and age uncertainties for each level. There are multiple built-in dissimilarity coefficients (DC) for different types of data, and various levels of data smoothing that can be applied depending on the type and variance of the data. In addition, R-Ratepol can use randomisation, accompanied by use of age uncertainties of each level and taxon standardisation to detect RoC patterns in datasets with high data noise or variability (i.e. numerous rapid changes in composition or sedimentation rates). 


The computation of RoC in R-Ratepol is performed using the following steps::

1.  1.	Assemblage and age-model data are extracted from the original source and is compiled together, i.e. depth, age, taxon 1, taxon 2, etc.

2. (optional) Smoothing of assemblage data: Each variable within the assemblage data is smoothed using one of five in-built smoothing methods: none, Shepard's 5-term filter (Davis, 1986; Wilkinson, 2005), moving average, age-weighted average, Grimm’s smoothing (Grimm and Jacobson, 1992). 

3.  Creation of time bins: A template for all time bins in all window movements is created.

4.  A single run (an individual loop) is computed:

    -   (optional) Selection of one time series from age uncertainties
    
    -   Subsetting levels in each bin: Here the working units are created 

    -   Standardisation of assemblage data in each WU.

    -   Calculation of RoC between adjacent WUs: RoC is calculated as the dissimilarity coefficient (DC) standardised by age differences between WUs. Five in-built dissimilarity coefficients are available: Euclidean distance, standardised Euclidean distance, Chord distance, Chi-squared coefficient (Prentice, 1980), Gower’s distance (Gower, 1971). The choice of DC depends on the type of data.
    
    -   Summary of a single run is produced based on all moving windows

5.  Step 4 is repeated multiple times (e.g. 10,000 times).

6.  Validation and summary of results from all runs of RoC calculation are produced.

7.  (Optional) Data beyond a certain age can be excluded.


Selection of working units (WU; Step 3)

RoC is calculated between consecutive Working Units (WU). Traditionally, these WUs represent individual stratigraphical levels. However, changes in sedimentation rates and sampling strategies can result in an uneven temporal distribution of levels within a time sequence, which in turn makes the comparison of RoC between sequences problematic. There are various methods that attempt to minimise such problems. The first is interpolation of levels to evenly spaced time intervals, and the use of the interpolated data as WUs. This can lead to a loss of information when the density of levels is high. Second is binning of levels: assemblage data are pooled into age brackets of various size (bins) and these serve as WUs. Here, the issue is a lower resolution of WUs and their uneven size in terms of total assemblage count (bins with more levels have higher assemblage counts). Third is selective binning: like classical binning, bins of selected size are created, but instead of pooling assemblage data together, only one level is selected as representative of each bin. This results in an even number of WUs in bins with a similar count size in the assemblage. However, the issue of low resolution remains. 
Therefore, we propose a new method of binning with a moving window, which is a compromise between using individual levels and selective binning. This method follows a simple sequence: bins are created, levels are selected as in selective binning, and RoC between bins is calculated. However, the brackets of the time bin (window) are then moved forward by a selected amount of time (Z), levels are selected again, and RoC calculated for the new set of WUs. This is repeated X times (where X is the bin size divided by Z) while retaining all the results. R-Ratepol currently provides several options for selecting WU, namely as individual levels, binning of levels, selective binning, and our new method of binning with a moving window.


Randomisation

Due to the inherent statistical errors in uncertainties in the age estimates from age-depth and the assemblage datasets (e.g. pollen counts in each level; Birks and Gordon, 1985), R-Ratepol can be run several times and the results summarised (Steps 5-6). Therefore, two optional settings are available by using age uncertainties and assemblage data standardisation. 

Age uncertainties 

For each run, a single age sequence from the age uncertainties is randomly selected. The calculation between two consecutive WUs (i.e. one working-unit combination) results in a RoC score and a time position (which is calculated as the mean age position of the two WUs). However, due to random sampling of the age sequence, each WU combination will result in multiple RoC values and age positions. R Ratepol assigns the age position of each WU combination as the median time position from all calculations (Step 6). 

The final RoC value for a single WU combination is calculated as the median of the scores from all randomisations. In addition, the 95th quantile from all randomisations is calculated as an error estimate.


Data standardisation (Step 4b)

Taxa in the assemblage dataset can be standardised to a certain count (e.g. number of pollen grains in each WU) by rarefaction. Random sampling without replacement is used to draw a selected number of individuals from each WU (e.g. 150 pollen grains).

}

\references{
Birks, H.J.B., Gordon, A.D., 1985. Numerical Methods in Quaternary Pollen Analysis. Academic Press, London.

Davis, J.C., 1986. Statistics and Data Analysis in Geology, 2nd edn. ed. J. Wiley & Sons, New York.

Gower, J.C., 1971. A general coefficient of similarity and some of its properties. Biometrics 27, 857–871.

Grimm, E.C., Jacobson, G.L., 1992. Fossil-pollen evidence for abrupt climate changes during the past 18000 years in eastern North America. Clim. Dyn. 6, 179–184.

Prentice, I.C., 1980. Multidimensional scaling as a research tool in Quaternary palynology: A review of theory and methods. Rev. Palaeobot. Palynol. 31, 71–104. https://doi.org/10.1016/0034-6667(80)90023-8

Wilkinson, L., 2005. The Grammar of Graphics. Springer-Verlag, New York, USA 37. https://doi.org/10.2307/2669493

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

fc_plot_RoC_sequence(
  sequence_01,
  age_threshold = 8e3,
  Roc_threshold = 1)

}

