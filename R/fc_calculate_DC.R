#' @title Calculate the dissimilarity coeficient
#'
#' @inheritParams fc_estimate_RoC
#' @param data_source_DC Data.frame with taxons as columns
#' @details
#' Five in-built dissimilarity coefficients are available:
#' \itemize{
#' \item Euclidean distance (`DC` = `"euc"`)
#' \item standardised Euclidean distance (`DC` = `"euc.sd"`)
#' \item Chord distance (`DC` = `"chord"`)
#' \item Chi-squared coefficient (`DC` = `"chisq"`)
#' \item Gower's distance (`DC` = `"gower"`)
#' \item Bray-Curtis distance (`DC` = `"bray"`)
#' }
#' The choice of DC depends on the type of assemblage data. In addition, RoC
#' between WUs be calculated using every consecutive WU (`only_subsequent` = `FALSE`),
#' or alternatively, calculation of RoC can be restricted to only directly
#' adjacent WUs (`only_subsequent` = `TRUE`). Using the former increases
#' the number of samples for which RoC can be calculated within a sequence,
#' which varies in terms of sample resolution, but may still introduce
#' biases related to the RoC estimation as a result of the varying
#' inter-sample distances.
#' @seealso [vegan::vegdist()]
fc_calculate_DC <-
  function(data_source_DC,
           DC = "chord",
           verbose = FALSE) {
    n_res <-
      nrow(data_source_DC) - 1

    # pre-allocate some space
    dat_res <-
      vector(
        mode = "numeric",
        length = n_res
      )

    data_com <-
      util_subset_community(data_source_DC)

    #----------------------------------------------------------#
    # Standardised euclidan distace -----
    #----------------------------------------------------------#

    # for euc.sd use custom made calculation

    if (
      DC == "euc.sd"
    ) {
      if (
        verbose == TRUE
      ) {
        RUtilpol::output_comment(
          "Standardised Euclidan distance will be used as DC"
        )
      }

      # calculation of standard deviation for each species

      # calculate the SD for each species
      df_sp_supp <-
        apply(data_com, 2, stats::sd)

      # calculation of the DC
      # for each sample (except the last)
      for (i in 1:n_res) {

        # select only 2 samples (observed + 1 after)
        df_work <-
          data_com[c(i, i + 1), , drop = FALSE]

        # get rid of "empty species" in data & in sp.std
        df_sp_supp_work <-
          df_sp_supp[colSums(df_work, na.rm = TRUE) > 0]
        df_work <-
          as.data.frame(df_work[, colSums(df_work, , na.rm = TRUE) > 0])

        # vector for result for each species
        vector_work <- vector(
          mode = "numeric",
          length = ncol(df_work)
        )

        # for each species
        for (j in 1:ncol(df_work)) {

          # check if the standard deviation is not equal zero
          if (df_sp_supp_work[j] != 0) {
            a <-
              .subset2(df_work, j)[1]
            b <-
              .subset2(df_work, j)[2]

            # calculate the difference
            vector_work[j] <-
              ((a - b) / df_sp_supp_work[j])**2
          }
        }

        # save the square root of sum of all differece
        dat_res[i] <- sqrt(sum(vector_work))
      }

      return(dat_res)
    }

    # for all other DC u se pre-made function from vegan package to
    #   calculate correlation between all samples and then
    #   only include calculation between subsequent samples

    #----------------------------------------------------------#
    # Euclidan distance -----
    #----------------------------------------------------------#

    if (
      DC == "euc"
    ) {
      if (
        verbose == TRUE
      ) {
        RUtilpol::output_comment(
          "Euclidan distance will be used as DC"
        )
      }

      corrmat <-
        as.matrix(
          vegan::vegdist(
            data_com,
            method = "euclidean"
          )
        )
    }
    #----------------------------------------------------------#
    # Chord's distance -----
    #----------------------------------------------------------#

    if (
      DC == "chord"
    ) {
      if (
        verbose == TRUE
      ) {
        RUtilpol::output_comment(
          "Chord distance will be used as DC"
        )
      }

      # use pre-made function from vegan package to
      #   calculate correlation between all samples
      corrmat <-
        as.matrix(
          vegan::vegdist(
            data_com,
            method = "chord"
          )
        )
    }


    #----------------------------------------------------------#
    # Chi-squared coeficient -----
    #----------------------------------------------------------#

    if (
      DC == "chisq"
    ) {
      if (
        verbose == TRUE
      ) {
        RUtilpol::output_comment(
          "Chi-squared coeficient will be used as DC"
        )
      }

      # use pre-made function from vegan package to
      #   calculate correlation between all samples
      corrmat <-
        as.matrix(
          vegan::vegdist(
            data_com,
            method = "chord"
          )
        )
    }


    #----------------------------------------------------------#
    # Gower's distance -----
    #----------------------------------------------------------#

    if (
      DC == "gower"
    ) {
      if (
        verbose == TRUE
      ) {
        RUtilpol::output_comment(
          "Gower's distance will be used as DC"
        )
      }

      # use pre-made function from vegan package to
      #   calculate correlation between all samples
      corrmat <- as.matrix(
        vegan::vegdist(
          data_com,
          method = "gower"
        )
      )
    }


    #----------------------------------------------------------#
    # Bray–Curtis distance -----
    #----------------------------------------------------------#

    if (
      DC == "bray"
    ) {
      if (
        verbose == TRUE
      ) {
        RUtilpol::output_comment(
          "Bray-Curtis distance will be used as DC"
        )
      }

      # use pre-made function from vegan package to
      #   calculate correlation between all samples
      corrmat <-
        as.matrix(
          vegan::vegdist(
            data_com,
            method = "bray"
          )
        )
    }

    # only include calculation between subsequent samples
    dat_res <-
      corrmat[row(corrmat) == col(corrmat) + 1]

    #----------------------------------------------------------#
    # Return result -----
    #----------------------------------------------------------#

    return(dat_res)
  }
