#' @title Smooth the community data
#'
#' @param data_source_smooth
#' List with `community`, and `age`
#' @param round_results
#' Logical. Should smoothed values be rounded to integers?
#' @inheritParams estimate_roc
#' @description
#' A function to apply one of the 4 smoothers.
#' @details
#' Smoothing of assemblage data: Each variable within the
#' assemblage data is smoothed using one of five in-built smoothing methods:
#' \itemize{
#' \item Shepard's 5-term filter (`smooth_method` = `"shep"`;
#'  Davis, 1986; Wilkinson, 2005)
#' \item moving average (`smooth_method` = `"m.avg"}`)
#' \item age-weighted average (`smooth_method` = `"age.w"`)
#' \item Grimm's smoothing (`smooth_method` = `"grim"`;
#'  Grimm and Jacobson, 1992)
#' }
#' @seealso [estimate_roc()]
#' @references
#' Davis, J.C., 1986. Statistics and Data Analysis in Geology, 2nd edn. ed.
#' J. Wiley & Sons, New York.
#'
#' Grimm, E.C., Jacobson, G.L., 1992. Fossil-pollen evidence for abrupt
#' climate changes during the past 18000 years in eastern North America.
#' Clim. Dyn. 6, 179-184.
#'
#' Wilkinson, L., 2005. The Grammar of Graphics. Springer-Verlag, New York,
#' USA 37.
smooth_community_data <-
  function(data_source_smooth,
           smooth_method = c("m.avg", "grim", "age.w", "shep"),
           smooth_n_points = 5,
           smooth_n_max = 9,
           smooth_age_range = 500,
           round_results = FALSE,
           verbose = FALSE) {
    # ----------------------------------------------
    # SETUP -----
    # ----------------------------------------------
    # Mandatory assertions for all methods:
    # Assert data_source_smooth is a list
    assertthat::assert_that(
      is.list(data_source_smooth),
      msg = "'data_source_smooth' must be a list"
    )
    # Assert that community is not empty
    assertthat::assert_that(
      nrow(data_source_smooth$community) > 0,
      msg = "'community' must have at least one row"
    )
    assertthat::assert_that(
      ncol(data_source_smooth$community) > 0,
      msg = "'community' must have at least one column"
    )
    # Assert community and age have the same rownames
    assertthat::assert_that(
      identical(
        rownames(data_source_smooth$community),
        rownames(data_source_smooth$age)
      ),
      msg = "'community' and 'age' must have identical rownames"
    )
    # Assert all values in community are numeric
    assertthat::assert_that(
      all(
        sapply(data_source_smooth$community, is.numeric)
      ),
      msg = "All values in 'community' must be numeric"
    )
    # Assert that there are no NAs
    assertthat::assert_that(
      !any(is.na(data_source_smooth$community)),
      msg = "'community' must not contain any NAs"
    )
    assertthat::assert_that(
      !any(is.na(data_source_smooth$age)),
      msg = "'age' must not contain any NAs"
    )

    # Assert age is sorted
    assertthat::assert_that(
      is.unsorted(data_source_smooth$age$age) == FALSE,
      msg = "'age' must be sorted in increasing order"
    )
    # Assert smooth_method is character and valid
    assertthat::assert_that(
      is.character(smooth_method),
      msg = "'smooth_method' must be character"
    )
    smooth_method <-
      match.arg(
        smooth_method,
        choices = c("m.avg", "grim", "age.w", "shep")
      )
    # Assert smooth_n_points is only 1 argument, numeric and smaller than N samples
    assertthat::assert_that(
      length(smooth_n_points) == 1,
      msg = "'smooth_n_points' must be length 1"
    )
    assertthat::assert_that(
      is.numeric(smooth_n_points),
      msg = "'smooth_n_points' must be numeric"
    )
    assertthat::assert_that(
      smooth_n_points <= nrow(data_source_smooth$community),
      msg = "'smooth_n_points' must be <= number of samples in 'community'"
    )
    # Assert that logical parameters are logical and either TRUE or FALSE
    assertthat::assert_that(
      is.logical(round_results) && length(round_results) == 1 && (round_results == TRUE || round_results == FALSE),
      msg = "'round_results' must be logical and either TRUE or FALSE"
    )
    assertthat::assert_that(
      is.logical(verbose) && length(verbose) == 1 && (verbose == TRUE || verbose == FALSE),
      msg = "'verbose' must be logical and either TRUE or FALSE"
    )


    # Method-specific assertions:
    if (smooth_method == "shep") {
      # must be > 2
      assertthat::assert_that(
        smooth_n_points > 2,
        msg = "'smooth_n_points' must be > 2 for 'shep' smoothing"
      )
    }

    if (smooth_method %in% c("m.avg", "age.w", "grim")) {
      # must be odd
      assertthat::assert_that(
        smooth_n_points %% 2 != 0,
        msg = "'smooth_n_points' must be odd"
      )

      if (smooth_method %in% c("age.w", "grim")) {
        # smooth_age_range must be length 1
        assertthat::assert_that(
          length(smooth_age_range) == 1,
          msg = "'smooth_age_range' must be length 1"
        )
        # smooth_age_range must be numeric
        assertthat::assert_that(
          is.numeric(smooth_age_range),
          msg = "'smooth_age_range' must be numeric"
        )

        if (smooth_method == "grim") {
          # smooth_n_max must be length 1
          assertthat::assert_that(
            is.numeric(smooth_n_max),
            msg = "'smooth_n_max' must be numeric"
          )
          assertthat::assert_that(
            smooth_n_max %% 2 != 0,
            msg = "'smooth_n_max' must be odd"
          )
          assertthat::assert_that(
            length(smooth_n_max) == 1,
            msg = "'smooth_n_max' must be length 1"
          )
          assertthat::assert_that(
            smooth_n_max <= nrow(data_source_smooth$community),
            msg = "'smooth_n_max' must be <= number of samples in 'community'"
          )
          # smooth_n_points must be < smooth_n_max
          assertthat::assert_that(
            smooth_n_points < smooth_n_max,
            msg = "'smooth_n_points' must be < 'smooth_n_max' for 'grim' smoothing"
          )
        }
      }
    }

    # ----------------------------------------------
    # Additional information -----
    # ----------------------------------------------

    if (
      isTRUE(verbose)
    ) {
      switch(smooth_method,
        "m.avg" = {
          RUtilpol::output_comment(
            paste(
              "Data will be smoothed by 'moving average' over", smooth_n_points,
              "points"
            )
          )
        },
        "grim" = {
          RUtilpol::output_comment(
            paste(
              "Data will be smoothed by 'Grimm method' with min samples",
              smooth_n_points,
              "max samples", smooth_n_max, "and max age range of",
              smooth_age_range
            )
          )
        },
        "age.w" = {
          RUtilpol::output_comment(
            paste(
              "Data will be smoothed by 'age-weighed average' over",
              smooth_n_points,
              "points with a threshold of", smooth_age_range
            )
          )
        },
        "shep" = {
          RUtilpol::output_comment(
            paste(
              "Data will be smoothed by 'Shepard's 5-term filter'"
            )
          )
        }
      )
    }

    # ----------------------------------------------
    # CALCULATION -----
    # ----------------------------------------------

    # split data into 2 datasets
    dat_community <- as.data.frame(data_source_smooth$community)
    dat_age <- as.data.frame(data_source_smooth$age)

    # pre-allocate some space
    focus_par <- matrix(data = NA, nrow = nrow(dat_age), ncol = 2)

    # for every species
    for (j in 1:ncol(dat_community)) {
      # select the species
      col_work <- .subset2(dat_community, j)

      # create empty vector of same lengt for values to be saved
      col_res <- rep(0, length(col_work))

      for (i in 1:nrow(dat_community)) { # for each sample

        # ----------------------------------------------
        # MOVING AVERAGE SMOOTHING -----
        # ----------------------------------------------
        if (
          smooth_method == "m.avg"
        ) {
          # Samples near beginning (moving window truncated)
          if (
            i < round(0.5 * (smooth_n_points)) + 1
          ) {
            focus_par[i, ] <- c(1, (i + round(0.5 * (smooth_n_points))))
          } else {
            # Samples near end
            if (
              i > nrow(dat_age) - round(0.5 * (smooth_n_points))
            ) {
              focus_par[i, ] <-
                c(
                  (i - round(0.5 * (smooth_n_points))),
                  nrow(dat_age)
                )
            } else {
              focus_par[i, ] <-
                c(
                  (i - round(0.5 * (smooth_n_points))),
                  (i + round(0.5 * (smooth_n_points)))
                )
            }
          }
          col_res[i] <-
            mean(col_work[focus_par[i, 1]:focus_par[i, 2]])
        }

        # ----------------------------------------------
        # GRIMMM SMOOTHING -----
        # ----------------------------------------------
        if (
          smooth_method == "grim"
        ) {
          # Samples near beginning (moving window truncated)
          if (
            i < round(0.5 * (smooth_n_max)) + 1
          ) {
            focus_par[i, 1] <- 1

            focus_par[i, 2] <-
              (i + round(0.5 * (smooth_n_points)))

            focus_par[i, ] <-
              util_search_parameter(
                focus_par[i, 1],
                focus_par[i, 2],
                smooth_age_range,
                smooth_n_max,
                smooth_n_points,
                dat_age,
                dat_community
              )
          } else {
            # Samples near end
            if (
              i > nrow(dat_age) - round(0.5 * (smooth_n_points))
            ) {
              focus_par[i, 1] <-
                (i - round(0.5 * (smooth_n_points)))

              focus_par[i, 2] <-
                nrow(dat_age)

              focus_par[i, ] <-
                util_search_parameter(
                  focus_par[i, 1],
                  focus_par[i, 2],
                  smooth_age_range,
                  smooth_n_max,
                  smooth_n_points,
                  dat_age,
                  dat_community
                )
            } else {
              focus_par[i, 1] <-
                (i - round(0.5 * (smooth_n_points)))

              focus_par[i, 2] <-
                (i + round(0.5 * (smooth_n_points)))

              focus_par[i, ] <-
                util_search_parameter(
                  focus_par[i, 1],
                  focus_par[i, 2],
                  smooth_age_range,
                  smooth_n_max,
                  smooth_n_points,
                  dat_age,
                  dat_community
                )
            }
          }
          col_res[i] <-
            mean(col_work[focus_par[i, 1]:focus_par[i, 2]])
        }

        # ----------------------------------------------
        # AGE-WEIGHTED SMOOTHING -----
        # ----------------------------------------------
        if (
          smooth_method == "age.w"
        ) {
          # Samples near beginning (moving window truncated)
          if (
            i < round(0.5 * (smooth_n_points)) + 1
          ) {
            focus_par[i, ] <-
              c(
                1,
                (i + round(0.5 * (smooth_n_points)))
              )
          } else {
            # Samples near end
            if (
              i > nrow(dat_age) - round(0.5 * (smooth_n_points))
            ) {
              focus_par[i, ] <-
                c(
                  (i - round(0.5 * (smooth_n_points))),
                  nrow(dat_age)
                )
            } else {
              focus_par[i, ] <-
                c(
                  (i - round(0.5 * (smooth_n_points))),
                  (i + round(0.5 * (smooth_n_points)))
                )
            }
          }

          # create small df with values around observed sample
          #   (in range of offset)
          df_work <-
            data.frame(
              values = col_work[focus_par[i, 1]:focus_par[i, 2]],
              age = dat_age$age[focus_par[i, 1]:focus_par[i, 2]],
              weight = 1
            )

          # Weith of points is calculated as smooth_age_range / distance
          #   bewtween oldest and youngest points.
          # If cannot be smaller than 1. Values very far away from the point
          F_age_dist <- abs(df_work$age - dat_age$age[i])

          const <- smooth_age_range / F_age_dist

          const[const > 1] <- 1

          df_work$weight <- const

          col_res[i] <-
            stats::weighted.mean(df_work$values, df_work$weight)
        }

        # ----------------------------------------------
        # Shepard's 5-term filter -----
        # ----------------------------------------------
        if (
          smooth_method == "shep"
        ) {
          if (
            i < round(0.5 * (smooth_n_points)) + 1
          ) {
            col_res[i] <- col_work[i]
          } else {
            if (
              i > nrow(dat_age) - round(0.5 * (smooth_n_points))
            ) {
              col_res[i] <- col_work[i]
            } else {
              w.value <-
                (17 * .subset(col_work, i) +
                  12 * (.subset(col_work, i + 1) + .subset(col_work, i - 1)) -
                  3 * (.subset(col_work, i + 2) + .subset(col_work, i - 2))) / 35
              if (
                w.value < 0
              ) {
                w.value <- 0
              }
              col_res[i] <- w.value
            }
          }
        }
      }
      dat_community[, j] <- col_res
    }

    if (
      isTRUE(round_results)
    ) {
      dat_community <- round(dat_community)
    }

    final_list <-
      list(
        community = dat_community,
        age = dat_age,
        age_un = data_source_smooth$age_un
      )

    return(final_list)
  }
