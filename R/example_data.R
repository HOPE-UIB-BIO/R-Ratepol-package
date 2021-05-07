#' Pollen data from four European sequences the Neotoma database (Goring et al., 2015)
#'
#' Taxa were standardised to the taxonomically highest pollen morphotype (Level = MHVar2) using the pollen harmonisation table in Giesecke et al. (2019).
#' Age-depth models were developed using the pre-selected radiometric control points provided in Giesecke et al. (2014)
#'  and calibrated the radiocarbon dates using the IntCal13 Northern Hemisphere calibration curve (Reimer et al., 2013). 
#'  For each sequence, an age-depth model was constructed using the Bchron R package (Haslett & Parnell, 2008) to 
#'  generate 1000 possible age predictions (i.e. age uncertainties) for all levels. We calculated the median of all the uncertainties for each level
#'   to give the most probable age (default age) in calibrated years before present (cal yr BP, where 0 = 1950 CE).
#' In each sequence, we excluded all levels that contained less than 150 pollen grain counts of the terrestrial taxa, 
#' and all levels beyond a 3000-years extrapolation of the oldest chronological control point. In addition, 
#' we excluded all levels with an age older than 8500 cal yr BP to focus on the period of most substantial human impact.
#'
#' @format A tibble with columns:
#' \describe{
#'  \item{dataset.id}{Neotoma unique number for each sequence}
#'  \item{collection.handle}{Neotoma sequence name abberation} 
#'  \item{lat}{Latitude}
#'  \item{long}{Longitude}
#'  \item{pollen_data}{Dataframe contain pollen counts of taxa in each level. }
#'  \item{sample_age}{Dataframe including estimated ages for each level}
#'  \item{age_uncertainty}{Matrix with age uncertainties from age-depth model}
#' }
#' @examples
#' \dontrun{
#'  example_data
#' }
"example_data"