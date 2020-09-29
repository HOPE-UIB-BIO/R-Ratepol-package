\name{example_data}
\alias{example_data}
\docType{data}
\title{
Pollen data from four European sequences
}
\description{
Tibble cointaining four European sequences obtained the Neotoma database (Goring et al., 2015) with age-depth models from Bchron (Haslett & Parnell, 2008) and hamronised pollen taxa.
}
\usage{data("example_data")}
\format{
  A data frame with 4 observations on the following 6 variables.
  \describe{
    \item{\code{dataset.id}}{a character vector containing Neotoma unique number for each sequence}
    \item{\code{collection.handle}}{a character vector containing Neotoma sequence name abberation}
    \item{\code{lat}}{a numeric vector of Latitude}
    \item{\code{long}}{a numeric vector of Longitude}
    \item{\code{list_ages}}{a list vector including estimated ages for each level \code{ages} and all age uncertainties from age-depth model \code{age_position} }
    \item{\code{filtered.counts}}{a list vector cointaining tibbles with pollen counts of taxa in each level}
  }
}
\details{
Taxa were standardised to the taxonomically highest pollen morphotype (Level = MHVar2) using the pollen harmonisation table in Giesecke et al. (2019).
Age-depth models were developed using the pre-selected radiometric control points provided in Giesecke et al. (2014)
and calibrated the radiocarbon dates using the IntCal13 Northern Hemisphere calibration curve (Reimer et al., 2013). 
For each sequence, an age-depth model was constructed using the Bchron R package (Haslett & Parnell, 2008) to 
generate 1000 possible age predictions (i.e. age uncertainties) for all levels. We calculated the median of all the uncertainties for each level
to give the most probable age (default age) in calibrated years before present (cal yr BP, where 0 = 1950 CE).
In each sequence, we excluded all levels that contained less than 150 pollen grain counts of the terrestrial taxa, 
and all levels beyond a 3000-years extrapolation of the oldest chronological control point. In addition, 
we excluded all levels with an age older than 8500 cal yr BP to focus on the period of most substantial human impact.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(example_data)
}
\keyword{datasets}
