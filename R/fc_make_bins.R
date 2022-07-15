fc_make_bins <-
  function(data_source_bin,
           bin_size = 500) {
    util_check_class("data_source_bin", "list")

    util_check_class("bin_size", "numeric")

    util_check_if_integer("bin_size")

    bin_oldest <-
      ceiling(
        max(data_source_bin$age$age)
      )

    bin_youngest <-
      floor(
        min(data_source_bin$age$age)
      )

    bin_breaks <- seq(
      from = bin_youngest,
      to = bin_oldest,
      by = bin_size
    )

    data.frame(
      name = bin_breaks
    ) %>%
      return()
  }
