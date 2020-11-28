\name{fc_detect_peak_points}
\alias{fc_detect_peak_points}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Detect significant peak points
}
\description{
Detect points of suddent increae of Rate-of-change values
}
\usage{
fc_detect_peak_points(data_source, method = "trend_linear")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_source}{
data from fc_estimae_RoC function.
}
  \item{method}{
method of peak detection:
  
"threshold"        = Treshold value is set for whole dataset (after subseting for interest_threshold) as median of all RoC values.              
                      Peak is  consider significant if 95 percent quantile (gain from randomisation) is higher than treshold

"trend_linear"     = linear model is fitted with RoC and Age. Differences between the model and each point is calculated. SD is calculated                          from all the differences.Peak is considered significat if it is 1.5 SD higher than the model
  
"trend_non_linear" = very crude GAM model (k = 3) is fitted with RoC and Age. Differences between GAM and each point is calculated. SD is                             calculated from all the differences. Peak is considered significat if it is 1.5 SD higher than GAM
            
"GAM_deriv"        = smooth GAM model is fitted with RoC and Age. First derivative as welll as continuous confidence intervals are calculted                          from the model. Peak is considered significat if confidence intervals of first derivative is higher than 0.

"SNI"              = Signal-to-Noise Index, following adapted the SNI from Kelly et al. (2011) written to detect changes in charcoal records. 
                      We SNI I calculated for the whole RoC sequence
                      and point is consider significant if has SNI value higher than 3 (following suggesting from Kelly et al. (2011)).
}
}
\details{
R-Ratepol detects peak points using one of 5 methods: 
	
	Threshold method: Each point in the RoC sequence is compared to a median of all RoC scores from the whole (i.e. threshold value). The point is considered significant if the 95th quantile of the RoC scores from all calculations is higher than the threshold value. 
	
	Linear trend  method: A generalised linear model is fitted through the RoC scores and time positions (GLM = RoC ~ time), family = mgcv::Tweedie(p=2)). The distance between each point and the fitted value is calculated. A peak point is considered significant if the distance between the point and the fitted GLM is greater than 1.5 * the standard deviation of all differences. 
	
	Non-linear trend  method: A GAM is fitted through the RoC scores and time positions (GAM = RoC ~ s(time, k=3), family = mgcv::Tweedie(p=2)). The distance between each point and the fitted value is calculated. A peak point is considered significant if the distance between the point and the fitted GAM is greater than 1.5 * the standard deviation of all differences. 
	
	First derivative of GAM method: A GAM is fitted through the RoC scores and time positions (GAM = RoC ~ s(time), family = mgcv::Tweedie(p=2)).First derivative as welll as continuous confidence intervals are calculted from the model. Peak is considered significat if confidence intervals of first derivative is higher than 0. See Simpson 2018 for more information. 
	
  Signal-to-Noise Index (SNI) method: We adapted SNI from Kelly et al. (2011) developed to detect changes in charcoal records. The SNI is calculated for the whole RoC sequence and a peak point is considered significant if it has an SNI value higher than 3 with the GAM set as a reference value and with a time window of 500 yrs.   
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Simpson GL. 2018 Modelling palaeoecological time series using generalised additive models. Front. Ecol. Evol. 6, 1–29. (doi:10.3389/fevo.2018.00149)

Kelly RF, Higuera PE, Barrett CM, Feng Sheng H. 2011 A signal-to-noise index to quantify the potential for peak detection in sediment-charcoal records. Quat. Res. 75, 11–17. (doi:10.1016/j.yqres.2010.07.011)

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

