# R-Ratepol package 0.5.6.
### README
 - Figures are now saved in the new folder ('man' folder) with the new names.
 - cosmetics changes in the code (spaces, new lines, etc)

### Vignette
 - stop using pre-saved data
 - decrease the 'i_multiplier' to 0.5 to speed up the vignette building
 - change the 'treads' arguments in 'fc_estimate_RoC' to FALSE so the vignette can be build in machines without multiple cores.

### Other
 - remove 'RRatepol' namespace inside package for internal functions

# R-Ratepol package 0.5.5.
 - 'sample.id' or 'sample_id' can be used in as the sample identification
 - RRatepol::fc_plot_RoC_sequence: fix typo in argument 'age_threshold' 
 - citation(package = "RRatepol") now show correct citation
 - update the Figure in README file