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
                age_uncertainty = F,
                smooth_method = "none",
                smooth_N_points = 5,
                smooth_age_range = 500,
                smooth_N_max = 9,
                Working_Units = "levels",
                bin_size = 500,
                Number_of_shifts = 1,
                rand = 1,
                treads = F,
                standardise = F,
                N_individuals = 150,
                tranform_to_proportions = T,
                DC = "euc",
                interest_threshold = F,
                Peak_detection = "GAM",
                Debug = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_source_community}{
Community data with species as collumns and levels (samples) as rows. First column should be ID (row names)
}
  \item{data_source_age}{
Dataframe with two columns:
sample.id = unique ID of each level
age = age of level
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
Selection of units that the DC will be calculated between
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
 \item{Peak_detection}{
 method of peak detection:
  
"Threshold" = Treshold value is set for whole dataset (after subseting for interest_threshold) as median of all RoC values.                   
                Peak is  consider significant if 95 percent quantile (gain from randomisation) is higher than treshold

"linear" = linear model is fitted with RoC and Age. Differences between the model and each point is calculated. SD is calculated from all the                differences.
           Peak is considered significat if it is 1.5 SD higher than the model
  
"GAM" = very crude GAM model (k = 3) is fitted with RoC and Age. Differences between GAM and each point is calculated. SD is calculated from               all the differences
            Peak is considered significat if it is 1.5 SD higher than GAM
            
"GAM_deriv" = smooth GAM model is fitted with RoC and Age. First derivative as welll as continuous confidence intervals are calculted from the 
                model. 
                Peak is considered significat if confidence intervals of first derivative is higher than 0.

"SNI" = Signal-to-Noise Index, following adapted the SNI from Kelly et al. (2011)
                 written to detect changes in charcoal records. We SNI I calculated for the whole RoC sequence
                 and point is consider significant if has SNI value higher than 3 (following suggesting
                 from Kelly et al. (2011)).
                 
}

 \item{Debug}{
T = show messages from internal processes
}

}
\details{
Randomisation
Due to the inherent statistical errors in a pollen count in each level (Birks & Gordon, 1985) and uncertainties in the age estimates from age-depth models, R-Ratepol is run several times and the results summarised. Each run, one age prediction   is randomly selected from the pool of age sequences, and pollen taxa are standardised to a certain number of pollen grains in each working unit (WU) by rarefaction (using random sampling without replacement to a default of 150 pollen grains). The calculation between two consecutive WUs (i.e. one working-unit combination) results in a RoC score and a time position (calculated as mean age position of the WUs). However, due to random sampling of the age sequence, each WU combination will result in multiple RoC values and age positions. R Ratepol assigns the age position of each WU combination as the median time position from all calculations. The final RoC values are calculated as the median and 95th quantile of the scores from all calculations. In addition, because of excluding empty bins in the selective binning with a moving window approach, it is possible that some working-unit combinations will not be present in some calculations. Therefore, R-Ratepol only includes working-unit combinations that are present in at least 10 percent of all calculations.  

Description of individual steps:

1) Pollen data smoothing 
Each pollen taxon is smoothed by one of five smoothing functions in R-Ratepol. 
	None: Pollen data are not smoothed
	Moving average (Wilkinson, 2005): Each focus value is calculated as the average over N number of levels (preferably ½ N before and ½ after the focus level; values are adjusted at the beginning and end of the core. N must be odd number). Note that each calculation is done from scratch and results are saved separately in order to avoid cumulative rounding errors. Default N is set to 5.
	Grimm’s smoothing (Grimm & Jacobson, 1992): Similar to a moving average but N is not fixed. For each level, N is selected as an odd number between N_a and N_b  (N ⊂{N_a,N_b}) while maintaining the maximum age difference from the selected levels as 〖Age〗_MAX. Default values are set as N_a= 5, N_b = 9, 〖Age〗_MAX = 500. 
	Age-weighted average (Wilkinson, 2005): Similar to the moving average but the average is weighted by the age difference from the focus level and multiplied by 1/CONSTANT (CONSTANT is selected by user). To avoid up-weighting levels, if CONSTANT/AGEDIFF exceeds 1, it is saved as 1. This means that levels closer than the CONSTANT to the target age are given full weighting, but those farther away are downweighed by an amount increasing with age difference. CONSTANT has a default setting of 500 years.
	Shepard's 5-term filter (Davis, 1986; Wilkinson, 2005): Smoothing over 5 points following the equation below 
V_NEW=(17*V + 12*(V_((+1) )  + V_((-1) ) )-3*(V_((+2) )  + V_((-2) )))/35 ,
where V is the focal level value. All result values smaller than zero are saved as zero.

2) Creation of template for all time bins in all window movement
Bins of selected size are created, starting from zero. This is repeated many times, with each time bin (window) shifting by Z years forward. This is repeated X times, where X = bin size / Z. 

3a) Random time sequence sampling from age-depth model
One age sequence is randomly selected from the 1000 sequences generated from the Bchron age-depth models. Each of the levels is given one age value and the correct temporal order is checked in the sequence. 

3b) Creation of working units (WUs) and data standardisation
RoC can be calculated between individual levels or bins. In terms of calculation between levels, each level becomes one working unit. If bins are selected, the level closest to the start of each bin is selected as representative of that bin. All other levels are excluded. Bins with valid levels become the new WUs, and empty bins are excluded. Each WU is standardised to the same number of pollen grains (default value = 150), using random sampling of pollen grains without replacement. The dataset is reduced by excluding all WUs without any pollen and all pollen taxa without any counts.

3c) Calculate RoC between adjacent WUs
First, the dissimilarity coefficient between adjacent WUs is calculated using one of four dissimilarity coefficients (Prentice, 1980):
	Euclidean distance: 
√(∑_(i=1)^m▒〖(A_i-B_i)〗^2 ), where A_i and B_i are pollen values for WUs A and B of taxon i, and there are m taxa.

	Standardised Euclidean distance:
√(∑_(i=1)^m▒〖((A_i-B_i)/〖SD〗_i )〗^2 ) , where A_i and B_i are pollen values for WUs A and B of taxon i,〖SD〗_i is the standard deviation for taxon i calculated from the whole sequence, and there are m taxa.

	Chord distance:
√(∑_(i=1)^m▒〖(√(A_i )-√(B_i ))〗^2 ), where A_i and B_i are pollen values for WUs A and B of taxon i, and there are m taxa.

	Chi-squared coefficient:
√(∑_(i=1)^m▒〖(A_i-B_i)〗^2/((A_i+B_i))), where A_i and B_i are pollen values for WUs A and B of taxon i , and there are m taxa.

 Gower's disatnce:
Usabe for data with NA's and also community data in difefrent units

Next, the age difference between subsequent WUs is calculated. Finally, the RoC score between WUs is calculated as:
〖RoC〗_(A-B)=(〖DC〗_(A-B)*〖Age.diff〗_mean)/〖Age.diff〗_(A-B)  , where 〖DC〗_(A-B) is the dissimilarity between WUs A and B, 〖Age.diff〗_mean is the average age difference between WUs in the whole sequence, and 〖Age.diff〗_(A-B) is the age difference between WUs A and B. 

3d) Summarisation of results from all moving windows
Bins (windows) are moved forward by a selected amount of time Z, the steps 3a–3c are repeated for a new set of WUs (age sampling, levels selection, RoC calculation). This is repeated X times, where X=  (BIN size)/Z. Results from all window positions are merged and saved as a single value.

4) Step 3 is repeated X times, with default number of 10.000 repetitions

5) Validation and summary of results from all randomisations
Results from all randomisations are merged while keeping the identity of each working-unit combination. Due to random selection of the age sequence, each of the working-unit combinations will result in multiple RoC values and time positions. The time position of each working-unit combination is assigned as the median of the time positions from all randomisations. Final RoC values are calculated as the median score of all randomisations. In addition, due to excluding empty bins in selective binning with a moving window, there is a chance that some WU combinations will be present only in some randomisations. Therefore, only working-unit combinations that are present in at least 10 percent of all randomisations are saved.  

5) Data beyond the selected age are excluded (8000 cal yr BP in this study)

6) Detection of significant peak points
R-Ratepol detects peak points using one of 5 methods: 
	Threshold method: Each point in the RoC sequence is compared to a median of all RoC scores from the whole (i.e. threshold value). The point is considered significant if the 95th quantile of the RoC scores from all calculations is higher than the threshold value. 
	GAM method: A GAM is fitted through the RoC scores and time positions (GAM= RoC ~ s(time,k=3). The distance between each point and the fitted value is calculated. A peak point is considered significant if the distance between the point and the fitted GAM is greater than 1.5 * the standard deviation of all differences. 
Signal-to-Noise Index (SNI) method: We adapted SNI from Kelly et al. (2011) developed to detect changes in charcoal records. The SNI is calculated for the whole RoC sequence and a peak point is considered significant if it has an SNI value higher than 3 with the GAM set as a reference value and with a time window of 500 yrs.   

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
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
    example_data$filtered.counts[[1]],
    example_data$list_ages[[1]]$age,
    rand = 10000)

fc_plot_RoC_sequence(sequence_01)

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
