R-Ratepol package
=================

R-Ratepol is an R package for estimating rate of change (RoC) from
community data in time series.

Reference: Ondřej Mottl, John-Arvid Grytnes, Alistair W.R. Seddon,
Manuel J. Steinbauer, Kuber P. Bhatta, Vivian A. Felde, Suzette G.A.
Flantua, H. John B. Birks. Rate-of-change analysis in palaeoecology
revisited: a new approach bioRxiv 2020.12.16.422943; doi:
<a href="https://doi.org/10.1101/2020.12.16.422943" class="uri">https://doi.org/10.1101/2020.12.16.422943</a>

R-Ratepol is written as an R package (R Core Team 2018) and include
plethora of possible setting including a novel statistical approach to
evaluate RoC in a single stratigraphical sequence using community data
and age uncertainties for each level. There are multiple build-in
dissimilarity coefficients (DC), for different types of dataset, and
various levels of data smoothing (from none to Grimm’s) can be applied
depending on the variance of the data. In addition, R-Ratepol is able to
use randomisation, accompanied with taxa standardisation and/or usage of
uncertainty for age of each level, to detect RoC patterns in dataset
with more noise in data.

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

    -   (optional) Community data is standardised, i.e. community data
        is subsampled in each WU to a selected total count of
        individuals (e.g. 150 pollen grains).

    -   RoC between adjacent WUs is calculated as the dissimilarity
        coefficient (DC) standardised by age differences between WUs.
        Four in-build dissimilarity coefficients: Euclidean distance,
        standardised Euclidean distance, Chord distance, Chi-squared
        coefficient (Prentice, 1980), Gower’s distance (Gower, 1971).

5.  Single repetition (a loop) is repeated multiple times (e.g. 10,000
    times).

6.  Results from all runs are summarised.

7.  Significant peak-points are detected and validated. Five in-build
    methods of testing for the significance of individual RoC points:
    Threshold, Linear trend, Non-linear trend, first derivative of
    generalised additive model (f-deriv GAM; Simpson, 2018), and
    Signal-to-Noise Index (SNI; Kelly et al., 2011).

Selection of working units
--------------------------

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
while retaining all the results. The method is summarised in figure.

![](README_files/figure-source/Fig_Scheme.png)

Randomisation
-------------

Due to the inherent statistical errors in a community datasets
(e.g. pollen count in each level; Birks and Gordon, 1985) and
uncertainties in the age estimates from age-depth modelling, R-Ratepol
can be run several times and the results summarised. Therefore, two
optional settings can be used: usage of age uncertainties and community
data standardisation.

### Community data standardisation

Taxa in community dataset can standardised to a certain abundance
(e.g. number of pollen grains in each WU) by rarefaction. Random
sampling without replacement is used to a draw selected number of
individuals in each WU (e.g. 150 pollen).

### Age uncertainties

For each run, a single age sequence from age uncertainties is randomly
selected. The calculation between two consecutive WUs (i.e. one
working-unit combination) results in a RoC score and a time position
(which is calculated as the mean age position of the two WUs). However,
due to random sampling of the age sequence, each WU combination will
result in multiple RoC values and age positions. R Ratepol assigns the
age position of each WU combination as the median time position from all
calculations. The final RoC value for a single WU combination is
calculated as the median of the scores from all randomisations. In
addition, 95th quantile from all randomisations is calculated as error
estimate.

Detection of peak-points in RoC sequence
----------------------------------------

A rapid change in taxonomic composition or relative abundances of taxa
within the sequence can provides a means of comparing RoC between
sequences and interpreting the potential drivers of assemblage change.
To detect such significant peak-points of RoC scores in each sequence,
each point is tested to see if it represents a significant increase in
RoC values. There are various ways of detection peak points in time
series and R-Ratepol is able detect peak-points using five methods:

1.  Threshold: Each point in the RoC sequence is compared to a median of
    all RoC scores from the whole (i.e. threshold value). The point is
    considered significant if the 95th quantile of the RoC scores from
    all calculations is higher than the threshold value.

2.  Linear trend: linear model is fitted between RoC values and their
    ages. Differences between the model and each point is calculated
    (residuals). Standard deviation (SD) is calculated from all the
    residuals. Peak is considered significant if it is 1.5 SD higher
    than the model.

3.  Non-linear trend: A conservative generalised additive model (GAM) is
    fitted through the RoC scores and their ages (GAM= RoC ~ s(age,k=3)
    using mgcv package (Wood, 2011). The distance between each point and
    the fitted value is calculated (residuals). Standard deviation (SD)
    is calculated from all the residuals. Peak is considered significant
    if it is 1.5 SD higher than the model.

4.  F-deriv GAM: smooth GAM model is fitted the RoC scores and their
    ages (GAM= RoC ~ s(age). First derivative as well as continuous
    confidence intervals are calculated from the model using gratia
    package (Simpson, 2019). Peak is considered significant if
    confidence intervals of first derivative is higher than 0 (for more
    information see Simpson, 2018). SNI method: We adapted SNI from
    Kelly et al. (2011), which was developed to detect changes in
    charcoal stratigraphical records. SNI is calculated for the whole
    RoC sequence and a peak-point is considered significant if it has an
    SNI value higher than 3.

Examples
--------

### Instaling package

    devtools::install_github("HOPE-UIB-BIO/R-Ratepol-package")

### Example data

Pollen data from four European sequences the *Neotoma database* (Goring
et al., 2015). Taxa were standardised to the taxonomically highest
pollen morphotype (Level = MHVar2) using the pollen harmonisation table
in Giesecke et al. (2019).

Age-depth models were developed using the pre-selected radiometric
control points provided in Giesecke et al. (2014) and calibrated the
radiocarbon dates using the IntCal13 Northern Hemisphere calibration
curve (Reimer et al., 2013). For each sequence, an age-depth model was
constructed using the *Bchron R package* (Haslett & Parnell, 2008) to
generate 1000 possible age predictions (i.e. age uncertainties) for all
levels. We calculated the median of all the uncertainties for each level
to give the most probable age (default age) in calibrated years before
present (cal yr BP, where 0 = 1950 CE).

In each sequence, we excluded all levels that contained less than 150
pollen grain counts of the terrestrial taxa, and all levels beyond a
3000-years extrapolation of the oldest chronological control point. In
addition, we excluded all levels with an age older than 8500 cal yr BP
to focus on the period of most substantial human impact.

    library(RRatepol)
    library(tidyverse)

    example_data <-  RRatepol::example_data

    glimpse(example_data)

    ## Rows: 4
    ## Columns: 7
    ## $ dataset.id        <chr> "4012", "40951", "45314", "17334"
    ## $ collection.handle <chr> "DALLICAN", "STEERMOS", "KILOALA", "GL"
    ## $ lat               <dbl> 60.38736, 47.80567, 67.96611, 53.00735
    ## $ long              <dbl> -1.096480, 8.200150, 20.460278, -6.348035
    ## $ pollen_data       <list> [<tbl_df[63 x 51]>, <tbl_df[273 x 104]>, <tbl_df...
    ## $ sample_age        <named list> [<data.frame[63 x 3]>, <data.frame[273 x 3...
    ## $ age_uncertainty   <named list> [<matrix[1000 x 63]>, <matrix[1000 x 273]>...

    example_data %>%
      ggplot(aes(x=long, y=lat))+
      borders(fill = "gray90", colour = NA)+
      geom_point(shape = 0, size = 2)+
      geom_point(shape = 20, size = 2)+
      coord_quickmap(xlim = c(-10,25), ylim= c(47,70))+
      labs(x= "Longitude",y="Latitude")+
      theme_classic()

![](README_files/figure-markdown_strict/plot%20data-1.png)

### Example 1

Estimate RoC values for *Dallican Water* site using *Age-weighed
smoothing* of the data and *Chord dissimilarity* coefficient. Pollen
data will not standardised to a certain pollen count and age
uncertainties from *Bchron* will not be used.

    sequence_01 <- 
      fc_estimate_RoC(
        data_source_community = example_data$pollen_data[[1]],
        data_source_age = example_data$sample_age[[1]],
        age_uncertainty = F,
        smooth_method = "age.w",
        DC = "chord",
        Working_Units = "levels",
        standardise = F,
        rand = 1)

    fc_plot_RoC_sequence(sequence_01, age_treshold = 8e3, Roc_threshold = 0.5, Peaks = F, trend = F)

![](README_files/figure-markdown_strict/plot%201-1.png)

### Example 2

Now try to standardise pollen data in each sample to a maximum of 150
pollen grains and use age uncertainties from *age-depth model*. Process
will be repeated 1000 times on multiple cores using *parallel package*.
This will produce error *wrapper* showing 95th percent quantile.

    sequence_02 <-
      fc_estimate_RoC(
        data_source_community = example_data$pollen_data[[1]],
        data_source_age = example_data$sample_age[[1]],
        age_uncertainty = example_data$age_uncertainty[[1]],
        smooth_method = "age.w",
        DC = "chord",
        Working_Units = "levels",
        standardise = T,
        N_individuals = 150,
        rand = 1000,
        treads = T)

    fc_plot_RoC_sequence(sequence_02, age_treshold= 8e3, Roc_threshold = 2.5, Peaks = F, trend = F)

![](README_files/figure-markdown_strict/plot%202-1.png)

### Example 3

Use *Binning with the mowing window* approach with `bin_size` = 500 and
`Number_of_shifts` = 5.

    sequence_03 <-
        fc_estimate_RoC(
        data_source_community = example_data$pollen_data[[1]],
        data_source_age = example_data$sample_age[[1]],
        age_uncertainty = example_data$age_uncertainty[[1]],
        smooth_method = "age.w",
        DC = "chord",
        Working_Units = "MW",
        bin_size = 500,
        Number_of_shifts  = 5,
        standardise = T,
        N_individuals = 150,
        rand = 1000,
        treads = T)

    fc_plot_RoC_sequence(sequence_03, age_treshold= 8e3, Roc_threshold = 1.5, Peaks = F, trend = F)

![](README_files/figure-markdown_strict/plot%203-1.png)

### Example 4

Detect the *peak points* using *trend\_non\_linear* method.

    sequence_03_peaks <-
      fc_detect_peak_points(sequence_03, method = "trend_non_linear")

    fc_plot_RoC_sequence(sequence_03_peaks, age_treshold= 8e3, Roc_threshold = 1, Peaks = T, trend = "trend_non_linear")

![](README_files/figure-markdown_strict/plot%204-1.png)
