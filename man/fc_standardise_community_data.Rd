\name{fc_standardise_community_data}
\alias{fc_standardise_community_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Standardies community data
}
\description{
Each Working Unit (sample) is standardised to the same number of individuals (default value = 150), using random sampling of individuals without replacement. The dataset is reduced by excluding all WUs without any individuals and all taxa without any counts.
}
\usage{
fc_standardise_community_data(data_source, N_individuals = 150, Debug = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_source}{
Data in the RRatepolList class list
}
  \item{N_individuals}{
Number of individuals to perform standardisation to. The N_individuals is automatically adjusted to the smallest  number of individuals in dataset. 
}
  \item{Debug}{
Show messages from internal processes
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
Ondřej Mottl
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}

