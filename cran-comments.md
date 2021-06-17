## Test environments
* local R istallation, R 4.0.1
* Fedora Linux, R-devel, clang, gfortran (R-hub)
* Ubuntu Linux 20.04.1 LTS, R-release, GCC (R-Hub)
* win-builder (devel)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 2 NOTEs:

* Possibly mis-spelled words in DESCRIPTION:
  Ratepol (3:10, 6:16)
  RoC (6:71)

Those are not mis-spelled words. 'Ratepol' is part of the name of the package 'R–Ratepol'. RoC stands for Rate of Change.

* checking R code for possible problems ... NOTE
Undefined global functions or variables:
  . Age ID Peak ROC ROC_dw ROC_up RoC_05q_sm RoC_95q_sm RoC_sm
  Working_Unit age_position bin bin_order l sample_id shift

All of those variables are new variables defined within the functions.
