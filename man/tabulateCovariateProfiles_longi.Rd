\name{tabulateCovariateProfiles_longi}
\alias{tabulateCovariateProfiles_longi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make table of covariate profiles with information for plotting
}
\description{
Used by \code{\link{plotProfilesByCluster_longi}}
}
\usage{
tabulateCovariateProfiles_longi(riskProfObj, whichCovariates = NULL, rhoMinimum = NULL, useProfileStar = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{riskProfObj}{
Object of type \code{riskProfObj}, output of
\code{calcAvgRiskAndProfile} or \code{calcAvgRiskAndProfile_longi}.
}
  \item{whichCovariates}{
A vector of indices or a vector of strings
corresponding to the covariates that are to be displayed.
}
  \item{rhoMinimum}{
Threshold value for rho: all covariates with mean rho values higher than rhoMinimum will be selected.
}
  \item{useProfileStar}{
The definition of the star
profile is given in Liverani, S., Hastie, D. I. and Richardson,
S. (2013) PReMiuM: An R package for Bayesian profile regression.
Ignored silently if variable selection was not used.
}
}
\details{
From the tabulateCovariateProfiles function in the R premiumPlots package (https://github.com/simisc/premiumPlots).
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
Rob Johson and Anais Rouanet
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
