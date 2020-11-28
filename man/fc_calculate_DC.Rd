\name{fc_calculate_DC}
\alias{fc_calculate_DC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate disimilarity between Working Units
}
\description{
The dissimilarity between adjacent Working Units is calculated using one of four dissimilarity coefficients
}
\usage{
fc_calculate_DC(data_source_DC, DC = "chord", Debug = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_source_DC}{
Data in the RRatepolList class list
}
  \item{DC}{
Disimilarity coeficient. Type of calculation of differences between Working Units

"euc"     = 	Euclidean distance

"euc.sd"  = 	Standardised Euclidean distance

"chord"   = 	Chord distance

"chisq    = 	Chi-squared coefficient

"gower"   =   Gower's distance
}
  \item{Debug}{
Show messages from internal processes
}
}
\details{
	Euclidean distance: 
√(∑_(i=1)^m▒〖(A_i-B_i)〗^2 ), where A_i and B_i are pollen values for WUs A and B of taxon i, and there are m taxa.

	Standardised Euclidean distance:
√(∑_(i=1)^m▒〖((A_i-B_i)/〖SD〗_i )〗^2 ) , where A_i and B_i are pollen values for WUs A and B of taxon i,〖SD〗_i is the standard deviation for taxon i calculated from the whole sequence, and there are m taxa.

	Chord distance:
√(∑_(i=1)^m▒〖(√(A_i )-√(B_i ))〗^2 ), where A_i and B_i are pollen values for WUs A and B of taxon i, and there are m taxa.

	Chi-squared coefficient:
√(∑_(i=1)^m▒〖(A_i-B_i)〗^2/((A_i+B_i))), where A_i and B_i are pollen values for WUs A and B of taxon i , and there are m taxa.

  Gower's disatnce:
Usabe for data with NA's and also community data in difefrent units

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

