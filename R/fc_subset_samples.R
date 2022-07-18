fc_subset_samples <-
  function(data_subset,
           bins,
           bin_selection = "first") {
    if (
      is.character(bins$start)
    ) {
      bins %>%
        dplyr::select(label, res_age, start) %>%
        dplyr::inner_join(
          data_subset %>%
            tibble::rownames_to_column("start"),
          by = "start"
        ) %>%
        dplyr::select(-start) %>%
        return()
    }

    res_com <-
      as.data.frame(
        matrix(
          nrow = nrow(bins),
          ncol = ncol(data_subset) - 1,
          dimnames = list(
            NULL,
            names(data_subset)[2:ncol(data_subset)]
          )
        )
      )

    # for each bin
    for (i in 1:nrow(bins)) {

      # subset age data so it selected all samples which has higher values
      # than the BIN itself but
      # still small then selected bin + calculated BIN size
      subset_w <-
        data_subset[data_subset$age >= bins$start[i] &
          data_subset$age < bins$end[i], ]

      # If selected subset has at least one sample
      if (nrow(subset_w) > 0) {
        if (
          bin_selection == "random"
        ) {
          # select one random sample from the bin
          random_row <-
            sample(1:nrow(subset_w), 1)

          res_com[i, ] <-
            subset_w[random_row, -1]
        }

        if (
          bin_selection == "first"
        ) {
          # select the sample which is the closest to the beggining of the bin
          res_com[i, ] <-
            subset_w[1, -1]
        }
      }
    }

    res <-
      dplyr::bind_cols(
        bins %>%
          dplyr::select(label, res_age),
        res_com
      )

    return(res)
  }
