#' @title Calculate the dissimilarity coeficient
#'
#' @inheritParams estimate_roc
#' @param data_source_dc Data.frame with taxons as columns
#' @details
#' Five in-built dissimilarity coefficients are available:
#' \itemize{
#' \item Euclidean distance (`dissimilarity_coefficient` = `"euc"`)
#' \item standardised Euclidean distance (`dissimilarity_coefficient` = `"euc.sd"`)
#' \item Chord distance (`dissimilarity_coefficient` = `"chord"`)
#' \item Chi-squared coefficient (`dissimilarity_coefficient` = `"chisq"`)
#' \item Gower's distance (`dissimilarity_coefficient` = `"gower"`)
#' \item Bray-Curtis distance (`dissimilarity_coefficient` = `"bray"`)
#' }
#' The choice of dissimilarity_coefficient depends on the type of assemblage data. In addition, RoC
#' between WUs be calculated using every consecutive WU (`only_subsequent` = `FALSE`),
#' or alternatively, calculation of RoC can be restricted to only directly
#' adjacent WUs (`only_subsequent` = `TRUE`). Using the former increases
#' the number of samples for which RoC can be calculated within a sequence,
#' which varies in terms of sample resolution, but may still introduce
#' biases related to the RoC estimation as a result of the varying
#' inter-sample distances.
#' @seealso [vegan::vegdist()]
estimate_dissimilarity_coefficient <-
  function(data_source_dc,
           dissimilarity_coefficient = "chord",
           verbose = FALSE) {
    n_res <-
      nrow(data_source_dc) - 1

    # pre-allocate some space
    dat_res <-
      vector(
        mode = "numeric",
        length = n_res
      )

    data_com <-
      subset_community(data_source_dc)

    #----------------------------------------------------------#
    # Standardised euclidan distace -----
    #----------------------------------------------------------#

    # for euc.sd use custom made calculation

    if (
      dissimilarity_coefficient == "euc.sd"
    ) {
      if (
        isTRUE(verbose)
      ) {
        RUtilpol::output_comment(
          "Standardised Euclidan distance will be used as dissimilarity_coefficient"
        )
      }

      # calculation of standard deviation for each species

      # calculate the SD for each species
      df_sp_supp <-
        apply(data_com, 2, stats::sd)

      # calculation of the dissimilarity_coefficient
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
          if (
            df_sp_supp_work[j] != 0
          ) {
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

    # for all other dissimilarity_coefficient u se pre-made function from vegan package to
    #   calculate correlation between all samples and then
    #   only include calculation between subsequent samples

    #----------------------------------------------------------#
    # Euclidan distance -----
    #----------------------------------------------------------#

    if (
      dissimilarity_coefficient == "euc"
    ) {
      if (
        isTRUE(verbose)
      ) {
        RUtilpol::output_comment(
          "Euclidan distance will be used as dissimilarity_coefficient"
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
      dissimilarity_coefficient == "chord"
    ) {
      if (
        isTRUE(verbose)
      ) {
        RUtilpol::output_comment(
          "Chord distance will be used as dissimilarity_coefficient"
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
      dissimilarity_coefficient == "chisq"
    ) {
      if (
        isTRUE(verbose)
      ) {
        RUtilpol::output_comment(
          "Chi-squared coeficient will be used as dissimilarity_coefficient"
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
      dissimilarity_coefficient == "gower"
    ) {
      if (
        isTRUE(verbose)
      ) {
        RUtilpol::output_comment(
          "Gower's distance will be used as dissimilarity_coefficient"
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
      dissimilarity_coefficient == "bray"
    ) {
      if (
        isTRUE(verbose)
      ) {
        RUtilpol::output_comment(
          "Bray-Curtis distance will be used as dissimilarity_coefficient"
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
