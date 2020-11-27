\name{fc_smooth_community_data}
\alias{fc_smooth_community_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Smoothing of data
}
\description{
Function to smooth each taxa in the dataset.
}
\usage{
fc_smooth_community_data(data_source, smooth_method = "age.w", smooth_N_points = 5, smooth_N_max = 9, smooth_age_range = 300, round_results = T, Debug = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_source}{
Data in the RRatepolList class list
}
  \item{smooth_method}{
type of smoothing applied for the each of the pollen type

"none"    = None: Pollen data is not smoothed

"m.avg"   = Moving average:

"grim"    = 	Grimm’s smoothing

"age.w"   = 	Age-weighted average

"shep"    =   Shepard's 5-term filter
}
  \item{smooth_N_points}{
Number of points for used for moving average, Grimm and Age-Weighted smoothing (odd number)
}
  \item{smooth_N_max}{
Maximal number of samples to look in Grimm smoothing
}
  \item{smooth_age_range}{
Maximal age range for both Grimm and Age-weight smoothing
}
  \item{round_results}{
Round result values ?
}
  \item{Debug}{
Show messages from internal processes
}
}
\details{
Types of smoothing:

	None: Pollen data are not smoothed

	Moving average (Wilkinson, 2005): Each focus value is calculated as the average over N number of levels (smooth_N_points; preferably ½ N before and ½ after the focus level; values are adjusted at the beginning and end of the core. N must be odd number). Note that each calculation is done from scratch and results are saved separately in order to avoid cumulative rounding errors. Default N is set to 5.

	Grimm’s smoothing (Grimm & Jacobson, 1992): Similar to a moving average but N is not fixed. For each level, N is selected as an odd number between smooth_N_points and smooth_N_max while maintaining the maximum age difference from the selected levels as smooth_age_range. Default value is smooth_age_range = 500. 

	Age-weighted average (Wilkinson, 2005): Similar to the moving average but the average is weighted by the age difference from the focus level and multiplied by 1/smooth_age_range. To avoid up-weighting levels, if smooth_age_range/AGEDIFF exceeds 1, it is saved as 1. This means that levels closer than the smooth_age_range to the target age are given full weighting, but those farther away are downweighed by an amount increasing with age difference. smooth_age_range has a default setting of 500 years.

	Shepard's 5-term filter (Davis, 1986; Wilkinson, 2005): Smoothing over 5 points following the equation below 
V_NEW=(17*V + 12*(V_((+1) )  + V_((-1) ) )-3*(V_((+2) )  + V_((-2) )))/35 ,
where V is the focal level value. All result values smaller than zero are saved as zero.
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
