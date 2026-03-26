

# RRatepol package <img src="man/figures/RRatepol_logo.png" align="right" width="200" alt="RRatepol package logo" />

<!-- badges: start -->

[![CRAN status](https://www.r-pkg.org/badges/version/RRatepol.png)](https://CRAN.R-project.org/package=RRatepol) [![R-CMD-check](https://github.com/HOPE-UIB-BIO/R-Ratepol-package/workflows/R-CMD-check/badge.svg)](https://github.com/HOPE-UIB-BIO/R-Ratepol-package/actions) [![Codecov-test-coverage](https://codecov.io/gh//HOPE-UIB-BIO/R-Ratepol-package/branch/main/graph/badge.svg)](https://app.codecov.io/gh//HOPE-UIB-BIO/R-Ratepol-package?branch=main)

<!-- badges: end -->

What is new in the package? See [NEWS](https://hope-uib-bio.github.io/R-Ratepol-package/news/index.html)

### Package logo

The original sketch for logo was done by [Vanesa Surtkova](https://www.instagram.com/vavatattoo/). Check her out!

## Description

RRatepol is an R package for estimating rate of change (RoC) from community data in time series.

RRatepol is written as an R package and includes a range of possible settings including a novel method to evaluate RoC in a single stratigraphical sequence using assemblage data and age uncertainties for each level. There are multiple built-in dissimilarity coefficients (dissimilarity_coefficient) for different types of assemblage data, and various levels of data smoothing that can be applied depending on the type and variance of the data. In addition, RRatepol can use randomisation, accompanied by use of age uncertainties of each level and taxon standardisation to detect RoC patterns in datasets with high data noise or variability (i.e. numerous rapid changes in composition or sedimentation rates).

## Installing package

``` r
devtools::install_github("HOPE-UIB-BIO/R-Ratepol-package")
```

## Cite as

``` r
citation(package = "RRatepol")
```

Ondrej Mottl, John-Arvid Grytnes, Alistair W.R. Seddon, Manuel J. Steinbauer, Kuber P. Bhatta, Vivian A. Felde, Suzette G.A. Flantua, H. John B. Birks. Rate-of-change analysis in palaeoecology revisited: a new approach Review of Palaeobotany and Palynology 293, doi: [![DOI badge](https://img.shields.io/badge/doi-10.1016/j.revpalbo.2021.104483-yellow.svg)](https://doi.org/10.1016/j.revpalbo.2021.104483)

## Package website

More detailed information can be found on [RRatepol package website](https://hope-uib-bio.github.io/R-Ratepol-package/)

This include description of the individual steps for RoC estimation [Package Description](https://hope-uib-bio.github.io/R-Ratepol-package/articles/package-description.html)

## FOSSILPOL

If you are interested in estimating rate of change for several records, please see [**FOSSILPOL**](https://hope-uib-bio.github.io/FOSSILPOL-website/), an R-based modular workflow to process multiple fossil pollen records to create a comprehensive, standardized dataset compilation, ready for multi-record and multi-proxy analyses at various spatial and temporal scales.

## Examples

### Workflow

Example of workflow showing full strength of RRatepol package, with as step by step guidance starting from downloading dataset from Neotoma, building age-depth models, to estimating rate-of-change using age uncertainty. [Example of full workflow](https://hope-uib-bio.github.io/R-Ratepol-package/articles/workflow-example.html)

### APD R-Ratepol workshop

For additional examples of RRatepol setting, see [Materials for R-Ratepol workshop with an African focus (APD data users)](https://ondrejmottl.github.io/APD_R-Ratepol_workshop/)

### OCCR R-Ratepol workshop

For examples using other data types than fosssil pollen, see [Oeschger Centre for Climate Change Research Workshop](https://ondrejmottl.github.io/OCCR_R-Ratepol_workshop/) for workflows using geochemistry and XRF data.

### Build-in example

Pollen data from four European sequences the *Neotoma database* (Goring et al., 2015) were obtained. Taxa were standardised to the taxonomically highest pollen morphotype (Level = MHVar2) using the pollen harmonisation table in Giesecke et al. (2019).

Age-depth models were developed using the pre-selected radiometric control points provided in Giesecke et al. (2014) and calibrated the radiocarbon dates using the IntCal13 Northern Hemisphere calibration curve (Reimer et al., 2013). For each sequence, an age-depth model was constructed using the *Bchron R package* (Haslett & Parnell, 2008) to generate 1000 possible age predictions (i.e. age uncertainties) for all levels. We calculated the median of all the uncertainties for each level to give the most probable age (default age) in calibrated years before present (cal yr BP, where 0 = 1950 CE).

In each sequence, we excluded all levels that contained less than 150 pollen grain counts of the terrestrial taxa, and all levels beyond a 3000-years extrapolation of the oldest chronological control point. In addition, we excluded all levels with an age older than 8500 cal yr BP to focus on the period of most substantial human impact.

``` r
library(RRatepol)
library(tidyverse)
```

``` r
example_data <-
  RRatepol::example_data

dplyr::glimpse(example_data)
#> Rows: 4
#> Columns: 7
#> $ dataset_id        <chr> "4012", "40951", "45314", "17334"
#> $ collection_handle <chr> "DALLICAN", "STEERMOS", "KILOALA", "GL"
#> $ lat               <dbl> 60.38736, 47.80567, 67.96611, 53.00735
#> $ long              <dbl> -1.096480, 8.200150, 20.460278, -6.348035
#> $ pollen_data       <list> [<tbl_df[63 x 51]>], [<tbl_df[273 x 104]>], [<tbl_df…
#> $ sample_age        <named list> [<data.frame[63 x 3]>], [<data.frame[273 x 3]>], [<d…
#> $ age_uncertainty   <named list> <<matrix[1000 x 63]>>, <<matrix[1000 x 273]>>, <<mat…
```

``` r
example_data %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = long,
      y = lat
    )
  ) +
  ggplot2::borders(
    fill = "gray90",
    colour = NA
  ) +
  ggplot2::geom_point(
    shape = 0,
    size = 2
  ) +
  ggplot2::geom_point(
    shape = 20,
    size = 2
  ) +
  ggplot2::coord_quickmap(
    xlim = c(-10, 25),
    ylim = c(47, 70)
  ) +
  ggplot2::labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  ggplot2::theme_classic()
#> Warning: `borders()` was deprecated in ggplot2 4.0.0.
#> 
