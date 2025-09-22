#' Window growth helper for GRIMM smoothing
#'
#' Expands the indices [A, B] while respecting (a) dataset bounds,
#' (b) maximum window size, and (c) maximum age range.
#'
#' @keywords internal
#' @noRd
util_search_parameter <- function(A, B, smooth_age_range,
                                  smooth_n_max, smooth_n_points,
                                  dat_age, dat_community) {
            # create new search parameter that is lower by 1
            A_test <- A - 1
            if (
                A_test > 0 && B - A_test < smooth_n_max
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
                B_test < nrow(dat_community) && B - A_test < smooth_n_max
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
