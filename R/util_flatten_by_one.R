util_flatten_by_one <-
    function(x) {
        len <- length(x)
        names_high <- names(x)
        y <- vector("list", 0)

        for (i in 1:len) {
            sel_item <- x[[i]]
            names_low <- names(sel_item)
            names(sel_item) <- paste(names_high[i], names_low, sep = "-")

            y <- c(y, sel_item)
        }
        return(y)
    }