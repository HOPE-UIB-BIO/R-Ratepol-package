#' @title Predict data for a fitted trend
#' @description
#' Fit a linear or GAM model to RoC scores and return the fitted values.
#' @param data_source
#' `data.frame` with columns `ROC` and `Age`.
#' @param sel_method
#' `character`. Which trend to fit:
#' \itemize{
#' \item `"linear"` - a linear model fitted between RoC and age.
#' \item `"non_linear"` - a conservative GAM
#' (`RoC ~ s(age, k = 3)`) using the `mgcv` package.
#' }
#' @return
#' `numeric` vector of fitted values, one per row of `data_source`.
#' @seealso [detect_peak_points()]
#' @keywords internal
make_trend <-
    function(data_source,
             sel_method = c("linear", "non_linear")) {
        util_check_class(data_source, "data.frame")

        util_check_col_names(data_source, c("ROC", "Age"))

        util_check_class(sel_method, "character")

        util_check_vector_values(sel_method, c("linear", "non_linear"))

        assertthat::assert_that(
            base::length(sel_method) == 1,
            msg = "'sel_method' must be a single value"
        )

        sel_method <- match.arg(sel_method)

        if (
            sel_method == "non_linear"
        ) {
            res <-
                mgcv::predict.gam(
                    mgcv::gam(
                        ROC ~ s(Age, k = 3),
                        data = data_source,
                        family = mgcv::Tweedie(p = 2),
                        method = "REML"
                    ),
                    type = "response"
                )
        } else {
            res <-
                stats::predict.glm(
                    stats::glm(ROC ~ Age,
                        data = data_source,
                        family = mgcv::Tweedie(p = 2)
                    ),
                    type = "response"
                )
        }

        return(base::as.numeric(res))
    }
