\name{setHyperparams}
\alias{setHyperparams}
\title{Definition of characteristics of sample datasets for profile regression}
\description{
Hyperparameters for the priors can be specified here and passed as an argument to profRegr.

The user can specify some or all hyperparameters.
Those hyperparameters not specified will take their default values. Where the
file is not provided, all hyperparameters will take their default values.
}
\usage{
setHyperparams(shapeAlpha = NULL, rateAlpha = NULL, aPhi = NULL, mu0
                    = NULL, Tau0 = NULL, R0 = NULL, kappa0 = NULL, nu0 =
                    NULL, muTheta = NULL, sigmaTheta = NULL, dofTheta =
                    NULL, muBeta = NULL, sigmaBeta = NULL, dofBeta = NULL,
                    shapeTauEpsilon = NULL, rateTauEpsilon = NULL, aRho =
                    NULL, bRho = NULL, atomRho = NULL, shapeSigmaSqY =
                    NULL, scaleSigmaSqY = NULL, rSlice = NULL,
                    truncationEps = NULL, shapeTauCAR = NULL, rateTauCAR =
                    NULL, shapeNu = NULL, scaleNu = NULL, initAlloc =
                    NULL, initL = NULL, muLSignal = 0, sigmaLSignal = 0,
                    muLLengthscale = 0, sigmaLLengthscale = 0, muLNoise =
                    0, sigmaLNoise = 0, MVNmu0 = NULL, MVNTau0 = NULL,
                    MVNR0 = NULL, MVNkappa0 = NULL, MVNnu0 = NULL)
}
\arguments{

\item{shapeAlpha}{The shape parameter for Gamma prior on alpha (default=2)}
\item{rateAlpha}{The inverse-scale (rate) parameter for the Gamma prior on alpha (default=1)}
\item{aPhi}{The vector of parameters for the Dirichlet prior on phi_j. Element j corresponds to covariate j which then has a prior Dirichlet(aPhi[j],aPhi[j],....,aPhi[j]). Only used in discrete case, default=(1 1 1 ... 1).}
\item{mu0}{The mean vector for mu_c in the Normal covariate case (only used in Normal covariate case, default=empirical covariate means)}
\item{Tau0}{The precision matrix for mu_c in the Normal covariate case (only used in Normal covariate case, default=inverse of diagonal matrix with elements equal to squareof empirical range for each covariate)}
\item{R0}{The matrix parameter for the Wishart distribution for Tau_c (only used in Normal covariate case, default=1/nCovariates * inverse of empirical covariance matrix)}
\item{kappa0}{The degrees of freedom parameter for the Wishart distribution for Tau_c (only used in Normal covariate case, default=nCovariates).}
\item{nu0}{Hyperparameter for the conjugate Normal inverse Wishart prior for Normal covariates. The Normal distribution of mu_c has covariance Sigma_c/nu0. The default value is 0.01. The other hyperparameters for this parametrisation are re-used from the independant priors. This hyperparameter is only useful when the option useNormInvWishPrior=TRUE in the function profRegr().}
\item{muTheta}{The location parameter for the t-Distribution for theta_c (only used if response included in model, default=0)}
\item{sigmaTheta}{The scale parameter for the t-Distribution for theta_c (only used if response included in model, default=2.5)}
\item{dofTheta}{The degrees of freedom parameter for the t-Distribution for theta_c (only used if response included in model, default=7)}
\item{muBeta}{The location parameter for the t-Distribution for beta (only used when fixed effects present, default=0)}
\item{sigmaBeta}{The scale parameter for the t-Distribution for beta (only used when fixed effects present, default=2.5)}
\item{dofBeta}{The dof parameter for the t-Distribution for beta (only used when fixed effects present, default=7)}
\item{shapeTauEpsilon}{Shape parameter for gamma distribution for prior for precision tau of extra variation errors epsilon (only used if extra variation is used i.e. extraYVar argument is included, default=5.0)}
\item{rateTauEpsilon}{Inverse-scale (rate) parameter for gamma distribution for prior for precision tau of extra variation errors epsilon (only used if extra variation is used i.e. extraYVar argument is used, default=0.5)}
\item{aRho}{Parameter for beta distribution for prior on rho in variable selection (default=0.5)}
\item{bRho}{Parameter for beta distribution for prior on rho in variable selection (default=0.5)}
\item{atomRho}{Parameter for the probability for the atom at zero, i.e. the 0.5 probability in w_j distributed Bernoulli(0.5) in the formulation of the sparsity inducing prior (default=0.5). This parameter must be in the interval (0,1], where atomRho=1 corresponds to the case where the prior for rho is a Beta(aRho,bRho).}
\item{shapeSigmaSqY}{Shape parameter of inverse-gamma prior for sigma_Y^2 (only used in the Normal response model, default =2.5)}
\item{scaleSigmaSqY}{Scale parameter of inverse-gamma prior for sigma_Y^2 (only used in the Normal response model, default =2.5)}
\item{rSlice}{Slice parameter for independent slice sampler such that xi_c = (1-rSlice)*rSlice^c for c=0,1,2,... (only used for slice independent sampler i.e. sampler=SliceIndependent, default 0.75).}
\item{truncationEps}{Parameter for determining the truncation level of the finite Dirichlet process (only used for truncated sampler i.e. sampler=Truncated}
\item{shapeTauCAR}{Shape parameter for gamma distribution for precision TauCAR of spatial CAR term (only used if a spatial term is included i.e. includeCAR argument is TRUE, default=0.001)}
\item{rateTauCAR}{Inverse-scale (rate) parameter for gamma distribution for precision TauCAR of spatial CAR term (only used if a spatial term is included i.e. includeCAR argument is TRUE, default=0.001)}
\item{shapeNu}{Shape parameter of Gamma prior for the shape parameter of the Weibull for survival response (only used in the Survival response model, default = 2.5)}
\item{scaleNu}{Scale parameter of Gamma prior for the shape parameter of the Weibull for survival response  (only used in the Survival response model, default = 1)}
\item{initAlloc}{Vector of the initial allocation of the individuals to clusters. This is NULL by default, which implies a random start. Useful for starting the MCMC from a specific partition. Note that if this overwrites the option nClusInit in the function profRegr: nClusInit is set equal to the maximum value in initAlloc. }
\item{initL}{Default to NULL.}
\item{muLSignal}{Signal-variance mean for Gaussian processes to model longitudinal outcomes.}
\item{sigmaLSignal}{Signal-variance variance for Gaussian processes to model longitudinal outcomes.}
\item{muLLengthscale}{Lengthscale mean for Gaussian processes to model longitudinal outcomes.}
\item{sigmaLLengthscale}{Lengthscale variance for Gaussian processes to model longitudinal outcomes.}
\item{muLNoise}{Noise-variance mean for Gaussian processes to model longitudinal outcomes.}
\item{sigmaLNoise}{Noise-variance variance for Gaussian processes to model longitudinal outcomes.}
\item{MVNmu0}{The mean vector for MVNmu_c in the MVN outcome case (default=empirical covariate means)}
\item{MVNTau0}{The precision matrix for MVNmu_c in the MVN outcome case (default=inverse of diagonal matrix with elements equal to square of empirical range for each outcome)}
\item{MVNR0}{The matrix parameter for the Wishart distribution for MVNTau_c (default=nSubjects/nOutcomes * inverse of empirical covariance matrix)}
\item{MVNnu0}{The degrees of freedom parameter for the Wishart distribution for MVNTau_c (default=nOutcomes).}
\item{MVNkappa0}{Hyperparameter for the conjugate Normal inverse Wishart prior for MVN outcome case. The Normal distribution of MVNmu_c has covariance MVNSigma_c/MVNkappa0. The default value is 0.01. The other hyperparameters for this parametrisation are re-used from the independent priors.}
}

\value{
The output of this function is a list with the components defined as above.
}
\section{Authors}{
David Hastie, Department of Epidemiology and Biostatistics, Imperial College London, UK

Silvia Liverani, Department of Epidemiology and Biostatistics, Imperial College London and MRC Biostatistics Unit, Cambridge, UK

Maintainer: Silvia Liverani <liveranis@gmail.com>
}
\references{

Silvia Liverani, David I. Hastie, Lamiae Azizi, Michail Papathomas, Sylvia Richardson (2015). PReMiuM: An R Package for Profile Regression Mixture Models Using Dirichlet Processes. Journal of Statistical Software, 64(7), 1-30. URL http://www.jstatsoft.org/v64/i07/.

}
\examples{

hyp <- setHyperparams(shapeAlpha=3,rateAlpha=2,mu0=c(30,13),R0=3.2*diag(2))

inputs <- generateSampleDataFile(clusSummaryPoissonNormal())
runInfoObj<-profRegr(yModel=inputs$yModel,
    xModel=inputs$xModel, nSweeps=2, nClusInit=15,
    nBurn=2, data=inputs$inputData, output="output",
    covNames = inputs$covNames, outcomeT = inputs$outcomeT,
    fixedEffectsNames = inputs$fixedEffectNames,
    hyper=hyp)

}
\keyword{hyperparameters}
