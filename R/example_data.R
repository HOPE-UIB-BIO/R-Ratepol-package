#' @title Pollen data from four European sequences the Neotoma database (Goring et al., 2015)
#'
#' @docType data
#' @description
#' Taxa were standardised to the taxonomically highest pollen morphotype (Level = MHVar2) using the pollen harmonisation table in Giesecke et al. (2019).
#' Age-depth models were developed using the pre-selected radiometric control points provided in Giesecke et al. (2014)
#' and calibrated the radiocarbon dates using the IntCal13 Northern Hemisphere calibration curve (Reimer et al., 2013).
#' For each sequence, an age-depth model was constructed using the Bchron R package (Haslett & Parnell, 2008) to
#' generate 1000 possible age predictions (i.e. age uncertainties) for all levels. We calculated the median of all the uncertainties for each level
#'  to give the most probable age (default age) in calibrated years before present (cal yr BP, where 0 = 1950 CE).
#' In each sequence, we excluded all levels that contained less than 150 pollen grain counts of the terrestrial taxa,
#' and all levels beyond a 3000-years extrapolation of the oldest chronological control point. In addition,
#' we excluded all levels with an age older than 8500 cal yr BP to focus on the period of most substantial human impact.
#' @format A data.frame with columns:
#' \itemize{
#'  \item `dataset.id` - Neotoma unique number for each sequence
#'  \item `collection.handle` - Neotoma sequence name aberration
#'  \item `lat` - Latitude (degrees)
#'  \item `long` - Longitude (degrees)
#'  \item `pollen_data` - Data.frame contain pollen counts of taxa in each level.
#'  \item `sample_age` - Data.frame including estimated ages for each level
#'  \item `age_uncertainty` - Matrix with age uncertainties from age-depth model
#' }
#' @keywords datasets
#' @source https://www.neotomadb.org/
#' @details
#' We obtained pollen data from the Neotoma database (Williams, Grimm, et al., 2018)
#' using the Neotoma R package (Goring et al., 2015). We chose four European sequences (A - D).
#' I each sequence, taxa were standardised to the taxonomically highest pollen morphotype (Level = MHVar2)
#' using the pollen harmonisation table in Giesecke et al. (2019).
#' To develop age-depth models, we used the pre-selected radiometric control points
#' provided in Giesecke et al. (2014) and calibrated the radiocarbon dates using the
#' IntCal13 Northern Hemisphere calibration curve (Reimer et al., 2013). For each sequence,
#' we constructed an age-depth model using the Bchron R package (Haslett and Parnell, 2008) to
#' generate 1000 possible age estimates for all sample depths at the original
#' sampling resolution of the original pollen sequences. We used these 1000
#' draws to build posterior estimates of age uncertainty. We calculated the
#' median age estimate for each sample depth to obtain the default age used
#' in following analyses.
#' In each sequence, we excluded all levels that contained less than 150 pollen
#' counts of terrestrial taxa, and all levels beyond a 3000-years extrapolation
#' of the oldest chronological control point. In addition, we excluded all
#' levels with an age older than 8500 cal yr BP to ensure focus on the period
#' with substantial human impact.
#' @references
#' Giesecke, T., Davis, B., Brewer, S., Finsinger, W., Wolters, S., Blaauw, M.,
#' de Beaulieu, J.L., Binney, H., Fyfe, R.M., Gaillard, M.J., Gil-Romera, G.,
#' van der Knaap, W.O., Kuneš, P., Kühl, N., van Leeuwen, J.F.N., Leydet, M.,
#' Lotter, A.F., Ortu, E., Semmler, M., Bradshaw, R.H.W., 2013. Towards mapping
#' the late Quaternary vegetation change of Europe. Veg. Hist. Archaeobot. 23,
#' 75–86.
#'
#' Giesecke, T., Wolters, S., van Leeuwen, J.F.N., van der Knaap, P.W.O.,
#' Leydet, M., Brewer, S., 2019. Postglacial change of the floristic diversity
#' gradient in Europe. Nat. Commun. 10.
#'
#' Goring, S., Dawson, A., Simpson, G.L., Ram, K., Graham, R.W., Grimm, E.C.,
#' Williams, J.W., 2015. Neotoma: A programmatic interface to the Neotoma
#' paleoecological database. Open Quat. 1, 1–17.
#'
#' Haslett, J., Parnell, A., 2008. A simple monotone process with application
#' to radiocarbon-dated depth chronologies. J. R. Stat. Soc. Ser. C Appl. Stat.
#' 57, 399–418.
#'
#' Reimer, P.J., Bard, E., Bayliss, A., Beck, J.W., Blackwell, P.G., Ramsey,
#' C.B., Buck, C.E., Cheng, H., Edwards, R.L., Friedrich, M., Grootes, P.M.,
#' Guilderson, T.P., Haflidason, H., Hajdas, I., Hatté, C., Heaton, T.J.,
#' Hoffmann, D.L., Hogg, A.G., Hughen, K.A.,
#'
#' Kaiser, K.F., Kromer, B., Manning, S.W., Niu, M., Reimer, R.W., Richards,
#' D.A., Scott, E.M., Southon, J.R., Staff, R.A., Turney, C.S.M.,
#' van der Plicht, J., 2013. IntCal13 and Marine13 Radiocarbon Age Calibration
#' Curves 0–50,000 years cal BP. Radiocarbon 55, 1869–1887.
#' @examples
#' \dontrun{
#' data(example_data)
#' }
"example_data"
