\name{calcAvgRiskAndProfile_longi}
\alias{calcAvgRiskAndProfile_longi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation of the average risks and profiles for longitudinal outcome data.
}
\description{
Calculation of the average risks and profiles for longitudinal outcome data.
}
\usage{
calcAvgRiskAndProfile_longi(clusObj, includeFixedEffects = F, proportionalHazards = F, nSweeps)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clusObj}{
Object of type \code{clusObj}.
}
  \item{includeFixedEffects}{
By default this is set to FALSE. If it is set to FALSE then the risk profile is computed with the parameters beta of the fixed effects assumed equal to zero. If it is set to TRUE, then risk profile at each sweep is computed adjusting for the sample of the beta parameter at that sweep.
}
  \item{proportionalHazards}{
Whether the risk matrix should include lambda only for the yModel="Survival" case so that the proportional hazards can be computed in the plotting function. The default is the average survival time.
}
  \item{nSweeps}{
Number of MCMC iterations.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with the following components. This is an object of type riskProfileObj.
\item{riskProfClusObj}{The object of type clusObj as given in the input of this function.}
\item{risk}{A matrix that has a column for each cluster and a row for each sweep. Each element of the matrix represents the estimated risk at each sweep for each cluster.}
\item{profile}{An array whose first dimension is the number of sweeps, the second is the number of clusters, the third is the number of discrete covariates and the fourth is the number of categories of each of the covariates. Each element of the array represents the covariate profile at each sweep for each cluster. The fourth dimension does not exists if the covariate type is Normal. If the covariate type is mixed, then instead of this element, the two elements below are defined, 'profilePhi' and 'profileMu'.}
\item{profileStar}{This is NULL if there has not been any variable selection. Otherwise it contains the selected profile.}
\item{empiricals}{A vector of length of the optimal number of clusters, where each value is the empirical mean of the outcome for each cluster.}
\item{profileStdDev}{An array whose first dimension is the number of sweeps, the second is the number of clusters, the third and the fourth are the number of continuous covariates. Each square matrix identified by the first and second dimension of the array represents the standard deviation at each sweep for each cluster. This element is only available if the covariate type is continuous or mixed.}
\item{profilePhi}{This array is the equivalent of the 'profile' above for discrete covariates in case of mixed covariates.}
\item{profileStarPhi}{This array is defined as profile and profilePhi, but the values are computed only if a variable selection procedure has been run. The definition of the star profile is given in Liverani, S., Hastie, D. I. and Richardson, S. (2013) PReMiuM: An R package for Bayesian profile regression.}
\item{profileMu}{This array is the equivalent of the 'profile' above for Normal covariates in case of mixed covariates.}
\item{profileStarMu}{This array is defined as profile and profileMu, but the values are computed only if a variable selection procedure has been run. The definition of the star profile is given in Liverani, S., Hastie, D. I. and Richardson, S. (2013) PReMiuM: An R package for Bayesian profile regression.}
\item{nuArray}{For yModel=Survival when weibullFixedShape=FALSE this array contains the sampled values of the shape parameter nu. The first dimension is the number of sweeps, the second is the number of clusters.}
\item{LArray}{For yModel=Longitudinal. This array contains the sampled values of the L parameters. The first dimension is the number of sweeps, the second is the number of clusters, the third has the L parameters in the order signal variance, lengthscale and noise variance.}
\item{MVNmuArray}{For yModel=MVN. This array contains the sampled values of the MVNmu parameters.}
\item{MVNSigmaArray}{For yModel=MVN. This array contains the sampled values of the lower triangular part of the matrix MVNSigma, ordered by row.}
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
\keyword{postprocessing}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
