#' @title Predict data for certain trend
#'
#' @param data_source
#' Data.frame with `ROC` and `Age`
#' @param sel_method
#' Which trend should be used:
#' \itemize{
#' \item `"linear"` - A linear model is fitted between the RoC values and
#' their ages.
#' \item `"non_linear"` - A conservative generalised additive model (GAM)
#' is fitted through the RoC scores and their ages (GAM = `RoC ~ s(age, k = 3)`
#' using the `mgcv` package (Wood, 2011).
#' }
#' @seealso [fc_detect_peak_points()]
util_make_trend <-
    function(data_source,
             sel_method = c("linear", "non_linear")) {
        RUtilpol::check_class("data_source", "data.frame")

        RUtilpol::check_col_names("data_source", c("ROC", "Age"))

        RUtilpol::check_class("sel_method", "character")

        RUtilpol::check_vector_values("sel_method", c("linear", "non_linear"))

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

        return(res)
    }
