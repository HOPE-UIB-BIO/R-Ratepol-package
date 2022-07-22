
util_make_trend <-
    function(data_source,
             sel_method = c("linear", "non_linear")) {
        util_check_class("data_source", "data.frame")

        util_check_col_names("data_source", c("ROC", "Age"))

        util_check_class("sel_method", "character")

        util_check_vector_values("sel_method", c("linear", "non_linear"))

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
