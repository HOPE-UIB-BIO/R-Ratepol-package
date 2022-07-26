# RRatepol package <img src="man/figures/RRatepol_logo.png" align="right" width="200" />

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/RRatepol)](https://CRAN.R-project.org/package=RRatepol)
[![R-CMD-check](https://github.com/HOPE-UIB-BIO/R-Ratepol-package/workflows/R-CMD-check/badge.svg)](https://github.com/HOPE-UIB-BIO/R-Ratepol-package/actions)
<!-- badges: end -->

## Current version: 1.0.0

What is new in the package? See
[NEWS](https://hope-uib-bio.github.io/R-Ratepol-package/news/index.html)

## Description

RRatepol is an R package for estimating rate of change (RoC) from
community data in time series.

RRatepol is written as an R package and includes a range of possible
settings including a novel method to evaluate RoC in a single
stratigraphical sequence using assemblage data and age uncertainties for
each level. There are multiple built-in dissimilarity coefficients (DC)
for different types of assemblage data, and various levels of data
smoothing that can be applied depending on the type and variance of the
data. In addition, RRatepol can use randomisation, accompanied by use of
age uncertainties of each level and taxon standardisation to detect RoC
patterns in datasets with high data noise or variability (i.e. numerous
rapid changes in composition or sedimentation rates).

## Instaling package

    devtools::install_github("HOPE-UIB-BIO/R-Ratepol-package")

## Cite as

    citation(package = "RRatepol")

Ondřej Mottl, John-Arvid Grytnes, Alistair W.R. Seddon, Manuel J.
Steinbauer, Kuber P. Bhatta, Vivian A. Felde, Suzette G.A. Flantua, H.
John B. Birks. Rate-of-change analysis in palaeoecology revisited: a new
approach Review of Palaeobotany and Palynology 293, doi:
[![](https://img.shields.io/badge/doi-10.1016/j.revpalbo.2021.104483-yellow.svg)](https://doi.org/10.1016/j.revpalbo.2021.104483)

## Package website

More detailed information can be found on [RRatepol package
website](https://hope-uib-bio.github.io/R-Ratepol-package/)

This include description of the individual steps for RoC estimation
[Package
Description](https://hope-uib-bio.github.io/R-Ratepol-package/articles/package-description.html)

## Examples

### Workflow

Example of workflow showing full strength of RRatepol package, with as
step by step guidance starting from downloading dataset from Neotoma,
building age-depth models, to estimating rate-of-change using age
uncertainty. [Example of full
workflow](https://hope-uib-bio.github.io/R-Ratepol-package/articles/workflow-example.html)

### Seminar Day of PalaeoSIG

For additional examples of RRatepol setting, see [Materials for Seminar
Day of PalaeoSIG: Data Manipulation & Statistical Analysis in
Palaeoecology](https://github.com/OndrejMottl/palaeoR_R-Ratepol_workshop)

### Build-in example

Pollen data from four European sequences the *Neotoma database* (Goring
et al., 2015) were obtained. Taxa were standardised to the taxonomically
highest pollen morphotype (Level = MHVar2) using the pollen
harmonisation table in Giesecke et al. (2019).

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

    example_data <-  
      RRatepol::example_data

    dplyr::glimpse(example_data)
    #> Rows: 4
    #> Columns: 7
    #> $ dataset_id        <chr> "4012", "40951", "45314", "17334"
    #> $ collection_handle <chr> "DALLICAN", "STEERMOS", "KILOALA", "GL"
    #> $ lat               <dbl> 60.38736, 47.80567, 67.96611, 53.00735
    #> $ long              <dbl> -1.096480, 8.200150, 20.460278, -6.348035
    #> $ pollen_data       <list> [<tbl_df[63 x 51]>], [<tbl_df[273 x 104]>], [<tbl_df~
    #> $ sample_age        <named list> [<data.frame[63 x 3]>], [<data.frame[273 x 3]>], [<d~
    #> $ age_uncertainty   <named list> <<matrix[1000 x 63]>>, <<matrix[1000 x 273]>>, <<mat~

    example_data %>%
      ggplot2::ggplot(
        ggplot2::aes(
          x = long,
          y = lat)) +
      ggplot2::borders(
        fill = "gray90",
        colour = NA) +
      ggplot2::geom_point(
        shape = 0,
        size = 2) +
      ggplot2::geom_point(
        shape = 20,
        size = 2) +
      ggplot2::coord_quickmap(
        xlim = c(-10, 25),
        ylim = c(47, 70)) +
      ggplot2::labs(
        x = "Longitude",
        y = "Latitude") +
      ggplot2::theme_classic()

![](man/figures/README-plot_data-1.png)

#### Example 1

Estimate RoC values for *Dallican Water* site using *Age-weighed
smoothing* of the data and *Chord dissimilarity* coefficient. Pollen
data will not standardised to a certain pollen count and age
uncertainties from *Bchron* will not be used.

    sequence_01 <- 
      RRatepol::fc_estimate_RoC(
        data_source_community = example_data$pollen_data[[1]],
        data_source_age = example_data$sample_age[[1]],
        smooth_method = "shep",
        DC = "chisq",
        Working_Units = "levels"
        )
    #> #----------------------------------------------------------#
    #> i RRatepol started 2022-07-26 15:44:23
    #> #----------------------------------------------------------#
    #> i RoC will be estimated between individual subsequent levels
    #> i 'time_standardisation' = 500 : RoC values will be reported as disimilarity per 500 years.
    #> #----------------------------------------------------------#
    #> i RRatepol finished 2022-07-26 15:44:25 taking 2.04 secs
    #> #----------------------------------------------------------#

    RRatepol::fc_plot_RoC_sequence(
      data_source = sequence_01)

![](man/figures/README-plot_1-1.png)

#### Example 2

Now try to standardise pollen data in each sample to a maximum of 150
pollen grains and use age uncertainties from *age-depth model*. Process
will be repeated 1000 times on multiple cores using parallel
computation. This will produce error *wrapper* showing 95th percent
quantile.

    sequence_02 <-
      RRatepol::fc_estimate_RoC(
        data_source_community = example_data$pollen_data[[1]],
        data_source_age = example_data$sample_age[[1]],
        age_uncertainty = example_data$age_uncertainty[[1]],
        smooth_method = "shep",
        DC = "chisq",
        Working_Units = "levels",
        standardise = TRUE,
        N_individuals = 150,
        rand = 1000,
        use_parallel = TRUE)
    #> #----------------------------------------------------------#
    #> i RRatepol started 2022-07-26 15:44:26
    #> #----------------------------------------------------------#
    #> i RoC will be estimated between individual subsequent levels
    #> i 'time_standardisation' = 500 : RoC values will be reported as disimilarity per 500 years.
    #> i Data will be standardise in each Working unit to 150 or the lowest number detected in dataset
    #> #----------------------------------------------------------#
    #> i RRatepol finished 2022-07-26 15:45:56 taking 1.51 mins
    #> #----------------------------------------------------------#

    RRatepol::fc_plot_RoC_sequence(
      data_source = sequence_02)

![](man/figures/README-plot_2-1.png)

#### Example 3

Use *Binning with the mowing window* approach with `bin_size` = 500 and
`Number_of_shifts` = 5.

    sequence_03 <-
        RRatepol::fc_estimate_RoC(
        data_source_community = example_data$pollen_data[[1]],
        data_source_age = example_data$sample_age[[1]],
        age_uncertainty = example_data$age_uncertainty[[1]],
        smooth_method = "shep",
        DC = "chisq",
        Working_Units = "MW",
        bin_size = 500,
        Number_of_shifts  = 5,
        standardise = TRUE,
        N_individuals = 150,
        rand = 1000,
        use_parallel = TRUE)
    #> #----------------------------------------------------------#
    #> i RRatepol started 2022-07-26 15:45:56
    #> #----------------------------------------------------------#
    #> i RoC will be estimated using 'binning with the mowing window' of 500 yr time bin over 5 number of window shifts
    #> i Sample will randomly selected for each bin
    #> i 'time_standardisation' = 500 : RoC values will be reported as disimilarity per 500 years.
    #> i Data will be standardise in each Working unit to 150 or the lowest number detected in dataset
    #> #----------------------------------------------------------#
    #> i RRatepol finished 2022-07-26 15:49:38 taking 3.7 mins
    #> #----------------------------------------------------------#

    RRatepol::fc_plot_RoC_sequence(
      data_source = sequence_03)

![](man/figures/README-plot_3-1.png)

#### Example 4

Detect the *peak points* using *trend\_non\_linear* method.

    sequence_03_with_peaks <-
      RRatepol::fc_detect_peak_points(
        data_source = sequence_03,
        sel_method = "trend_non_linear")

    RRatepol::fc_plot_RoC_sequence(
      data_source = sequence_03_with_peaks,
      Peaks = TRUE,
      trend = "trend_non_linear")

![](man/figures/README-plot_4-1.png)
