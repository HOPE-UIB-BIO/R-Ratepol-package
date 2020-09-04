\name{fc_estimate_RoC}
\alias{fc_estimate_RoC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
RRatepol
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fc_estimate_RoC <- function (data_source_pollen,
                             data_source_age,
                             smooth_method = "grim",
                             smooth_N_points = 5,
                             smooth_age_range = 500,
                             smooth_N_max = 9,
                             Working_Units = "MW",
                             bin_size = 500,
                             Number_of_shifts = 5,
                             rand = 1000,
                             treads = T,
                             age_uncertainty = T,
                             standardise = T,
                             N_pollen_grains = 150,
                             DC = "Chord",
                             interest_threshold = F,
                             Peak_detection = "GAM",
                             Debug = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.source.pollen}{
pollen data with species as collumns and levels as rows, level ID as row names
}
  \item{data.source.age}{
list of 2:

[1] ages  dataframe:
sample.id = unique ID of each level
age = age of level

[2] age_position = matrix with number of collumns as number of levels. Each column is one level each row is one age sequence from bchron
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
