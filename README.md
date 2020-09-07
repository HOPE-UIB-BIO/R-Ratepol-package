---
title: "RRatepol"
author: "Ondrej Mottl"
date: "07/09/2020"
output:
  html_document:
    keep_md: true
---



## R-Ratepol package

R-Ratepol is an R package for estimating rate of change form community data in time series



```r
library(RRatepol)

example_data = readRDS("example_data/ex_data.rda")

sequence_01 = fc_estimate_RoC(example_data$filtered.counts[[1]],example_data$list_ages[[1]], rand = 100)
```

### Plots

You can also embed plots, for example:


```r
fc_plot_RoC_sequence(sequence_01)
```

![](README_files/figure-html/plot 1-1.png)<!-- -->


