#' @title Smooth the community data
#'
#' @param data_source_smooth
#' List with `community`, and `age`
#' @param round_results
#' Logical. Should smoothed values be rounded to integers?
#' @inheritParams fc_estimate_RoC
#' @description 
#' A function to apply one of the 4 smoothers.
#' @details 
#' Smoothing of assemblage data: Each variable within the
#' assemblage data is smoothed using one of five in-built smoothing methods:
#' \itemize{
#' \item Shepard's 5-term filter (`smooth_method` = `"shep"`;
#'  Davis, 1986; Wilkinson, 2005)
#' \item moving average ´(`smooth_method` = `"m.avg"}`)
#' \item age-weighted average (`smooth_method` = `"age.w"`)
#' \item Grimm's smoothing (`smooth_method` = `"grim"`;
#'  Grimm and Jacobson, 1992)
#' }
#' @seealso [fc_estimate_RoC()]
#' @references 
#' Davis, J.C., 1986. Statistics and Data Analysis in Geology, 2nd edn. ed.
#' J. Wiley & Sons, New York.
#'
#' Grimm, E.C., Jacobson, G.L., 1992. Fossil-pollen evidence for abrupt
#' climate changes during the past 18000 years in eastern North America.
#' Clim. Dyn. 6, 179–184.
#'
#' Wilkinson, L., 2005. The Grammar of Graphics. Springer-Verlag, New York,
#' USA 37.
fc_smooth_community_data <-
  function(data_source_smooth,
           smooth_method = c("m.avg", "grim", "age.w", "shep"),
           smooth_N_points = 5,
           smooth_N_max = 9,
           smooth_age_range = 500,
           round_results = FALSE,
           verbose = FALSE) {

    # ----------------------------------------------
    # SETUP -----
    # ----------------------------------------------

    RUtilpol::check_class("data_source_smooth", "list")

    RUtilpol::check_class("smooth_method", "character")

    RUtilpol::check_vector_values(
      "smooth_method",
      c("m.avg", "grim", "age.w", "shep")
    )

    smooth_method <- match.arg(smooth_method)

    if (
      !smooth_method != "shep"
    ) {
      assertthat::assert_that(
        smooth_N_points %% 2 != 0,
        msg = "'smooth_N_points' must be an odd number"
      )

      if (
        smooth_method != "m.avg"
      ) {
        RUtilpol::check_class("smooth_age_range", "numeric")

        if (
          smooth_method == "grim"
        ) {
          assertthat::assert_that(
            smooth_N_max %% 2 != 0,
            msg = "'smooth_N_max' must be an odd number"
          )

          assertthat::assert_that(
            smooth_N_points < smooth_N_max,
            msg = "'smooth_N_max' must be bigger than 'smooth_N_points"
          )
        }
      }
    }

    # crete helper function for GRIMM smoothing
    # test if this increase does not invalidate rules:
    #   1) seach parameter cannot go outside of the sample size
    #     (up or down)
    #   2) seach parameter cannot be biger than selected maximum sample
    #     sizes
    #   3) the age difference between samples selected by the seach
    #     paramated cannot be higher than defined max age range if all
    #     of those ARE TRUE then increase the real search parameter

    util_search_parameter <-
      function(A, B, smooth_age_range) {
        for (k in 1:(smooth_N_max - smooth_N_points)) {
          # create new search parameter that is lower by 1
          A_test <- A - 1
          if (
            A_test > 0 && B - A_test < smooth_N_max
          ) { # i+N.active.test < nrow(dat_community) &
            if (
              abs(dat_age$age[A_test] - dat_age$age[B]) < smooth_age_range
            ) {
              A <- A_test
            }
          }

          # create new search parameter that higher by 1
          B_test <- B + 1
          if (
            B_test < nrow(dat_community) && B - A_test < smooth_N_max
          ) {
            if (
              abs(dat_age$age[A] - dat_age$age[B_test]) < smooth_age_range
            ) {
              B <- B_test
            }
          }
        }
        return(c(A, B))
      }

    # ----------------------------------------------
    # Additional information -----
    # ----------------------------------------------

    if (verbose == TRUE) {
      switch(smooth_method,
        "m.avg" = {
          RUtilpol::output_comment(
            paste(
              "Data will be smoothed by 'moving average' over", smooth_N_points,
              "points"
            )
          )
        },
        "grim" = {
          RUtilpol::output_comment(
            paste(
              "Data will be smoothed by 'Grimm method' with min samples",
              smooth_N_points,
              "max samples", smooth_N_max, "and max age range of",
              smooth_age_range
            )
          )
        },
        "age.w" = {
          RUtilpol::output_comment(
            paste(
              "Data will be smoothed by 'age-weighed average' over",
              smooth_N_points,
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
            i < round(0.5 * (smooth_N_points)) + 1
          ) {
            focus_par[i, ] <- c(1, (i + round(0.5 * (smooth_N_points))))
          } else {
            # Samples near end
            if (
              i > nrow(dat_age) - round(0.5 * (smooth_N_points))
            ) {
              focus_par[i, ] <-
                c(
                  (i - round(0.5 * (smooth_N_points))),
                  nrow(dat_age)
                )
            } else {
              focus_par[i, ] <-
                c(
                  (i - round(0.5 * (smooth_N_points))),
                  (i + round(0.5 * (smooth_N_points)))
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
            i < round(0.5 * (smooth_N_max)) + 1
          ) {
            focus_par[i, 1] <- 1

            focus_par[i, 2] <-
              (i + round(0.5 * (smooth_N_points)))

            focus_par[i, ] <-
              util_search_parameter(
                focus_par[i, 1],
                focus_par[i, 2],
                smooth_age_range
              )
          } else {
            # Samples near end
            if (
              i > nrow(dat_age) - round(0.5 * (smooth_N_points))
            ) {
              focus_par[i, 1] <-
                (i - round(0.5 * (smooth_N_points)))

              focus_par[i, 2] <-
                nrow(dat_age)

              focus_par[i, ] <-
                util_search_parameter(
                  focus_par[i, 1],
                  focus_par[i, 2],
                  smooth_age_range
                )
            } else {
              focus_par[i, 1] <-
                (i - round(0.5 * (smooth_N_points)))

              focus_par[i, 2] <-
                (i + round(0.5 * (smooth_N_points)))

              focus_par[i, ] <-
                util_search_parameter(
                  focus_par[i, 1],
                  focus_par[i, 2],
                  smooth_age_range
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
            i < round(0.5 * (smooth_N_points)) + 1
          ) {
            focus_par[i, ] <-
              c(
                1,
                (i + round(0.5 * (smooth_N_points)))
              )
          } else {
            # Samples near end
            if (
              i > nrow(dat_age) - round(0.5 * (smooth_N_points))
            ) {
              focus_par[i, ] <-
                c(
                  (i - round(0.5 * (smooth_N_points))),
                  nrow(dat_age)
                )
            } else {
              focus_par[i, ] <-
                c(
                  (i - round(0.5 * (smooth_N_points))),
                  (i + round(0.5 * (smooth_N_points)))
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
            i < round(0.5 * (smooth_N_points)) + 1
          ) {
            col_res[i] <- col_work[i]
          } else {
            if (
              i > nrow(dat_age) - round(0.5 * (smooth_N_points))
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
      round_results == TRUE
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
