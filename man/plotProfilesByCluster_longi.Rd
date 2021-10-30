\name{plotProfilesByCluster_longi}
\alias{plotProfilesByCluster_longi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot longifutinal and covariate profiles by cluster
}
\description{
Plots the longitudinal outcome (yModel="MVN" or "Longitudinal") and (discrete) covariate profiles, facetted vertically by cluster and (optionally) horizontally by some grouping (split) of the covariates. Covariate categories that are more common within a cluster than they are overall are shaded dark; those that are not are shaded light.
}
\usage{
plotProfilesByCluster_longi(riskProfObj, whichCovariates = NULL, rhoMinimum = NULL, useProfileStar = TRUE, covariate_info = list(title = "Covariate\ncategory", levels = NULL, labels = NULL, split = NULL))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{riskProfObj}{
Object of type \code{riskProfObj}, output of calcAvgRiskAndProfile or calcAvgRiskAndProfile_longi.
}
  \item{whichCovariates}{
A vector of indices or a vector of strings
corresponding to the covariates that are to be displayed.
}
  \item{rhoMinimum}{
Threshold value for rho: all covariates with mean rho values higher than rhoMinimum will be selected.
}
  \item{useProfileStar}{
The definition of the star profile is given in Liverani, S., Hastie, D. I. and Richardson, S. (2013) PReMiuM: An R package for Bayesian profile regression.
}
\item{covariate_info}{Optional list of details about the covariates,
with (some of) the following named elements:
\describe{
    \item{\code{title}:}{String describing the covariate, default
    to "Covariate category".}
    \item{\code{levels}:}{Vector of integer values taken by the
    (discrete) covariates, specifying the order of the facets.}
    \item{\code{labels}:}{Vector of strings giving labels for each
    level of the covariate, in the same order as \code{levels}.}
    \item{\code{split}:}{A character vector specifying the covariate
    grouping for horizontal facetting. The first element appears
    in a subset of covariates to be plotted in a separate horizontal
    facet. The second and third elements, if present, are facet labels
    for covariates that do and do not meet this criterion, respectively,
    e.g. \code{c("_citr", "Citrullinated", "Not citrullinated")}.}
    }
  }
}
\details{
From the plotProfilesByCluster function in the R premiumPlots package  (https://github.com/simisc/premiumPlots).
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
Rob Johnson and Anais Rouanet
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
