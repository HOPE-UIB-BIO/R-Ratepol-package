\name{fc_extract_result}
\alias{fc_extract_result}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract results
}
\description{
Extract results for a single variable from all randomisation
}
\usage{
fc_extract_result(x, sel_var, rand)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Tibble with all the results from all randomisation and window shifts 
}
  \item{sel_var}{
Selected varibale that should been evaluated
}
  \item{rand}{
The number of total randomisations. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
tibble including the median and 95th quantile of the selected variable in each time bin
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
Because of excluding empty bins in the selective binning with a moving window approach, it is possible that some working-unit combinations will not be present in some calculations. Therefore, R-Ratepol only includes working-unit combinations that are present in at least 10 percent of all calculations.  
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

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
