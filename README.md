R-Ratepol package
-----------------

R-Ratepol is an R package for estimating rate of change form community
data in time series

    library(RRatepol)

    example_data = readRDS("example_data/ex_data.rda")

    sequence_01 = fc_estimate_RoC(example_data$filtered.counts[[1]],example_data$list_ages[[1]], rand = 100)

### Plots

You can also embed plots, for example:

    fc_plot_RoC_sequence(sequence_01)

![](README_files/figure-markdown_strict/plot%201-1.png)
