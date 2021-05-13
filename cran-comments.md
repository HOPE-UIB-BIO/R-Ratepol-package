## Test environments
* local R istallation, R 4.0.1
* Fedora Linux, R-devel, clang, gfortran (R-hub)
* Ubuntu Linux 20.04.1 LTS, R-release, GCC (R-Hub)
* win-builder (devel)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 3 NOTEs:

* Possibly mis-spelled words in DESCRIPTION:
  Ratepol (3:10, 6:16)
  RoC (6:71)

Those are not mis-spelled words. 'Ratepol' is part of the name of the package 'R–Ratepol'. RoC stands for Rate of Change.

* checking R code for possible problems ... NOTE
fc_detect_peak_points: no visible binding for global variable
  ‘Working_Unit’
fc_detect_peak_points: no visible binding for global variable ‘Age’
fc_detect_peak_points: no visible binding for global variable ‘ROC’
fc_detect_peak_points: no visible binding for global variable ‘ROC_up’
fc_detect_peak_points: no visible binding for global variable ‘ROC_dw’
fc_detect_peak_points: no visible binding for global variable ‘Peak’
fc_estimate_RoC: no visible binding for global variable ‘l’
fc_estimate_RoC: no visible binding for global variable ‘age_position’
fc_estimate_RoC: no visible binding for global variable ‘sample_id’
fc_estimate_RoC: no visible binding for global variable ‘RoC_sm’
fc_estimate_RoC: no visible binding for global variable ‘RoC_95q_sm’
fc_estimate_RoC: no visible binding for global variable ‘RoC_05q_sm’
fc_extract_data: no visible binding for global variable ‘sample_id’
fc_extract_result: no visible binding for global variable ‘shift’
fc_extract_result: no visible binding for global variable ‘bin’
fc_extract_result: no visible binding for global variable ‘ID’
fc_extract_result: no visible binding for global variable ‘.’
fc_extract_result: no visible binding for global variable ‘bin_order’
fc_extract_result: no visible binding for global variable ‘sample_id’
fc_plot_RoC_sequence: no visible binding for global variable ‘Age’
fc_plot_RoC_sequence: no visible binding for global variable ‘ROC’
fc_plot_RoC_sequence: no visible binding for global variable ‘ROC_up’
fc_plot_RoC_sequence: no visible binding for global variable ‘ROC_dw’
fc_plot_RoC_sequence: no visible binding for global variable ‘.’
fc_plot_RoC_sequence: no visible binding for global variable ‘Peak’
Undefined global functions or variables:
  . Age ID Peak ROC ROC_dw ROC_up RoC_05q_sm RoC_95q_sm RoC_sm
  Working_Unit age_position bin bin_order l sample_id shift



All of those variables are new variables defined within the functions.

* checking examples ... NOTE
Examples with CPU (user + system) or elapsed time > 5s
                      user system elapsed
fc_plot_RoC_sequence 0.837  0.012    6.75

Our package works with lost of randomization and has ability to use parallel calculation, which is not used in the example.
