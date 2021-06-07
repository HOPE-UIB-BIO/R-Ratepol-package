# R-Ratepol package 0.6.0.
### Overhaul change of the estimation of RoC
 - RoC score is now returned in units selected by user (time_standardisation)
 - Method of sample selection from bin is now selectable (bin_selection)
 - RoC is estimated only for subsequent bins (only_subsequent == TRUE)
 - Bins are now created form the beginning of the core (instead of 0)

### doPararell
 - packages 'doSNOW' and 'snow' are replaced by 'doParallel'
 - progress bar for randomisation is currently not present 

### Other
 - Added various checks for correct argumet selection
 - all hyphens are replaced with dashes ('–') in the whole package
 - examples are wrapped into 'dontrun'

# R-Ratepol package 0.5.6.
### README
 - Figures are now saved in the new folder ('man' folder) with the new names.
 - cosmetics changes in the code (spaces, new lines, etc)

### Vignette
 - stop using pre-saved data
 - decrease the 'i_multiplier' to 0.5 to speed up the vignette building
 - change the 'treads' arguments in 'fc_estimate_RoC' to FALSE so the vignette can be build in machines without multiple cores.
 - rename r chunks to remove empty spaces in names

### Other
 - remove 'RRatepol' namespace inside package for internal functions
 - fix couple of typos in function descriptions
 - add cran-comments

# R-Ratepol package 0.5.5.
 - 'sample.id' or 'sample_id' can be used in as the sample identification
 - RRatepol::fc_plot_RoC_sequence: fix typo in argument 'age_threshold' 
 - citation(package = "RRatepol") now show correct citation
 - update the Figure in README file