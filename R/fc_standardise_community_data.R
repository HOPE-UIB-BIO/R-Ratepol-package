fc_standardise_community_data <-
  function(data_source_standard,
           N_individuals = 150) {
    data_community <-
      util_subset_community(data_source_standard)  %>% 
      round()

    n_taxa <-
      ncol(data_community)

    # for each row(sample)
    for (i in 1:nrow(data_community)) {

      # selected row
      select_row <-
        data_community[i, ]

      # a vector for the species pool
      vec1 <- vector(, length = 0)

      # create a vector with species numbers replicated X times,
      #  where X is number of individuals
      for (j in 1:n_taxa) {

        # repeat species names
        v1 <-
          rep(
            names(select_row)[j],
            select_row[j]
          )
        vec1 <-
          c(vec1, v1) # a vector that repeat the species for the occurrences
      }

      # sample species X time
      rsample <-
        sample(vec1,
          size = N_individuals,
          replace = FALSE
        )

      sel_names <-
        table(rsample)

      # replace all values in community data by 0
      data_community[i, ] <-
        rep(0, n_taxa)

      # replace individuals by new randomised values
      data_community[i, names(sel_names)] <-
        as.numeric(sel_names)
    }

    data_source_standard[, names(data_community)] <-
      data_community

    return(data_source_standard)
  }
