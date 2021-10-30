\name{plotRiskProfile_longi}
\alias{plotRiskProfile_longi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the Risk Profiles for longitudinal outcome
}
\description{
Plots the risk profiles for longitudinal outcome (yModel = "MVN" or "Longitudinal")
}
\usage{
plotRiskProfile_longi(riskProfObj, outFile, showRelativeRisk = F, orderBy = NULL, whichClusters = NULL, whichCovariates = NULL, useProfileStar = F, riskLim = NULL, bycol = FALSE, profile_X = NULL, timevar = NULL, double_plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{riskProfObj}{An object of type riskProfObj.}
\item{outFile}{Path and file name to save the plot.}
\item{showRelativeRisk}{Whether to show the relative risk (with respect to the risk of the first cluster). This option is not available for Normal outcome. For Survival outcomes it computed proportional hazards, but only if the option proportionalHazards=T was used in the function calcAvgRiskAndProfile().}
\item{orderBy}{Order by which the clusters are to be displayed. It can take values "Empirical", "ClusterSize" and "Risk" (the latter only if the outcome is provided). It can also take the name of a covariate to order the clusters, in which case the clusters are ordered.}
\item{whichClusters}{Either a vector of indeces that corresponds to the clusters that are to be displayed. The length of this vector must be greater than 1. The default is that all clusters are shown.}
\item{whichCovariates}{Either a vector of indeces or a vector of strings that corresponds to the covariates that are to be displayed. The length of this vector must be greater than 1. The default is that all covariates are shown.}
\item{useProfileStar}{To be set equal to TRUE only if a variable selection procedure has been run. The definition of the star profile is given in Liverani, S., Hastie, D. I. and Richardson, S. (2013) PReMiuM: An R package for Bayesian profile regression.}
\item{riskLim}{Limits of the y-axis for the plot of the boxplots for the response variable. The default is NULL. If the riskLim are provided, they should be a vector of length 2. }
  \item{bycol}{
Boolean. If TRUE, individual outcome trajectories are coloured by assigned cluster. FALSE by default.
}
  \item{profile_X}{
Profile of fixed effects if yModel="LME" (in development). By default, NULL.
}
  \item{timevar}{
A vector of time points to compute predictions if yModel="LME" (in development). By default, NULL.
}
  \item{double_plot}{
Bolean to make extra plots. See value. FALSE by default.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Plots the cluster-specific outcome profiles and covariates profiles. If double_plot=TRUE, makes extra plots:

\item{-trajectories-data.png}{
Individual trajectories  by clusters and cluster-specific estimated trajectories with 90\% credible intervals}

\item{-trajectories.png}{
Cluster-specific estimated trajectories }

\item{-all_trajectories.png}{
Cluster-specific estimated trajectories  and covariate profiles in terms of selected covariates }

\item{-all_trajectories-data.png}{
Individual trajectories by clusters  and covariate profiles in terms of selected covariates }

\item{-all.png}{
Individual trajectories by clusters  and covariate profiles }
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
