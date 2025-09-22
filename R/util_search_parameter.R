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
        for (k in 1:(smooth_n_max - smooth_n_points)) {
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
