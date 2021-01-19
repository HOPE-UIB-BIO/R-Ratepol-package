\name{fc_plot_RoC_sequence}
\alias{fc_plot_RoC_sequence}

\title{
Function to plot the Rate-of-Change sequence.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fc_plot_RoC_sequence(
  data_source,
  age_treshold = 15000,
  Roc_threshold = 2,
  Peaks = F,
  trend = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_source}{
Output of fc_estimate_RoC (tibble)
}
  \item{age_treshold}{
Cut-off value used as maximum age
}
  \item{Roc_threshold}{
Cut-off value used as maximum RoC value
}
  \item{Peaks}{
If peak-points are presented in the datasset and \code{Peaks} == \code{TRUE}, then peak points will be displayed  
}
  \item{trend}{
If peak-points are presented in the datasset and \code{Peaks} == \code{TRUE}, then one of the three method can be used to visulise the process of peak detection:

"threshold"         = Treshold value is set for whole dataset (after subseting for interest_threshold) as median of all RoC values. Peak is  consider significant if 95 percent quantile (gain from randomisation) is higher than treshold

"trend_linear"      = linear model is fitted with RoC and Age. Differences between the model and each point is calculated. SD is calculated from all the differences.Peak is considered significat if it is 1.5 SD higher than the model

"trend_non_linear"  = very crude GAM model (k = 3) is fitted with RoC and Age. Differences between GAM and each point is calculated. SD is calculated from all the differences. Peak is considered significat if it is 1.5 SD higher than GAM
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
Ondrej Mottl
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# example 1
sequence_01 <- 
  fc_estimate_RoC(
    data_source_community = example_data$pollen_data[[1]],
    data_source_age = example_data$sample_age[[1]])

fc_plot_RoC_sequence(
  sequence_01,
  age_treshold = 8e3,
  Roc_threshold = 1)

# example 2
sequence_02 <- 
  fc_estimate_RoC(
    data_source_community = example_data$pollen_data[[1]],
    data_source_age = example_data$sample_age[[1]])

sequence_02_peak <-
  fc_detect_peak_points(sequence_01, method = "trend_non_linear")
  
fc_plot_RoC_sequence(
  sequence_02_peak,
  age_treshold = 8e3,
  Roc_threshold = 2,
  Peaks = T,
  trend = "trend_non_linear")


}

