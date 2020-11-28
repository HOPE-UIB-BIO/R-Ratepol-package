\name{fc_subset_samples}
\alias{fc_subset_samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subsample data in selected time bins
}
\description{
Select a single sample from each time bin either a) randomly (WU="bins"), or b) that is the closest to the start of the time bin (WU="MW").
}
\usage{
fc_subset_samples(data_subset, bins)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_subset}{
Data in the RRatepolList class list
}
  \item{bins}{
Vector with selected time bins
}
 \item{WU}{
Working unit selection method
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
