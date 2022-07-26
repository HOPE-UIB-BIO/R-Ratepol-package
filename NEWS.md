# RRatepol 1.0.0

## General
* package has been re-factored (most internal processes have been updated)
* several arguments of functions have been changed and marked with `lifecycle`
* add several more checks while data sourcing

## Dissimilarity calculation
* use `vegan` package to calculate all dissimilarities except of standardised euclidean.

## parallel computation
* use the `pbapply` package 
* data for iteration is prepared separately
* added progress bar

## example data
* names for example columns as been adjusted to avoid using `.`


# RRatepol 0.6.1

## Dissimilarity calculation
* use `vegan` package to calculate Gower and Bray-Curtis (new) dissimilarity

## Univariate RoC estimation
* fix issues when using only single-taxon dataset 

## Code cosmetics changes
* TRUE/FALSE written as full
* all arguments are written for each function
* spelling fixes


# RRatepol 0.6.0

## Overhaul change of the estimation of RoC
* RoC score is now returned in units selected by user (time_standardisation)
* Method of sample selection from bin is now optional (bin_selection)
* RoC is estimated only for subsequent bins (only_subsequent == TRUE)
* Bins are now created form the beginning of the core (instead of 0)

## Console outputs
* General overhaul of console outputs presented to user while running RRatepol
* more information and warning messages are presented with each run

## doParallel
have been curated by you during the last years packages `doSNOW` and `snow` are replaced by `doParallel`
* progress bar for randomisation is currently not present 

## Other
* Added various checks for correct argument selection
* all hyphens are replaced with dashes (`–`) in the whole package
* examples are wrapped into `dontrun`


# RRatepol 0.5.6

## README
* Figures are now saved in the new folder (`man` folder) with the new names.
* cosmetics changes in the code (spaces, new lines, etc)

## Vignette
* stop using pre-saved data
* decrease the `i_multiplier` to 0.5 to speed up the vignette building
* change the `treads` arguments in `fc_estimate_RoC` to FALSE so the vignette can be build in machines without multiple cores.
* rename r chunks to remove empty spaces in names

## Other
* remove `RRatepol` namespace inside package for internal functions
* fix couple of typos in function descriptions
* add cran-comments


# RRatepol 0.5.5
* `sample.id` or `sample_id` can be used in as the sample identification
* RRatepol::fc_plot_RoC_sequence: fix typo in argument `age_threshold` 
* citation(package = "RRatepol") now show correct citation
* update the Figure in README file
