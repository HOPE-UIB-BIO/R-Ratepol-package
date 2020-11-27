\name{fc_create_bins}
\alias{fc_create_bins}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Crete a vector with time bins 
}
\description{
Function for pre-alocation of space for time bins of various size with a posibility of mowing window shift
}
\usage{
fc_create_bins(data_source_bin, shift_value, Number_of_shifts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_source_bin}{
Data in the RRatepolList class list
}
  \item{shift_value}{
Size of a time bin divided by the number of shifts (bin_size/Number_of_shifts)
}
  \item{Number_of_shifts}{
Number of mowing window shifts
}
}
\details{
Bins of selected size are created, starting from zero. This is repeated many times, with each time bin (window) shifting by Z years forward. This is repeated X times, where X = bin size / Z. 
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

}

