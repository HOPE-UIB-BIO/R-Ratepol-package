.onAttach <- function(lib, pkg) {
    packageStartupMessage(
        paste(
            "R-Ratepol version",
            utils::packageDescription("RRatepol",
                fields = "Version"
            ),
            "\n"
        ),
        appendLF = TRUE
    )
}