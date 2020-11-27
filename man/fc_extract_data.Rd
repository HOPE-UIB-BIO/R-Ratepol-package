\name{fc_extract_data}
\alias{fc_extract_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extracting data to correct format
}
\description{
Function to extract data and save it to the prefered format.
}
\usage{
fc_extract_data(data_community_extract, data_age_extract, age_uncertainty = F, Debug = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_community_extract}{
Community data with species as collumns and samples as rows. First column should be ID (row names) 
}
  \item{data_age_extract}{
Dataframe with two columns:
sample.id = unique ID of each level
age = age of level
}

  \item{age_uncertainty}{
Usage of age uncertainty form Age-depth models. Either:

[1] matrix with number of columns as number of samples. Each column is one sample, each row is one age sequence from age-depth model. Age sequence is randomly sampled from age-depth model uncertainties at the begining of each run.

[2] FALSE. Age uncertainties are not available and, therefore, will not be used.

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
