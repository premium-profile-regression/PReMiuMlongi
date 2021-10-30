\name{GP_post}
\alias{GP_post}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes the posterior mean and variance covariance matrix of a GP
}
\description{
A function to compute the posterior mean and variance covariance matrix of a GP
}
\usage{
GP_post(t0, ts, y, L, kernelType)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t0}{
A vector of times corresponding to outcome vector y (e.g. vector of times for subjects in a given cluster)
}
  \item{ts}{
A vector of times on which the mean and covariance matrix are computed.
}
  \item{y}{
A vector of outcome values (e.g. vector of outcomes for subjects in a given cluster)
}
  \item{L}{
A vector of hyperparameters for the covariance function (e.g. variance signal, lengthscale and measurement error variance for kernelType="SQexponential")
}
  \item{kernelType}{
Kernel function (yModel=Longitudinal) characterizing the covariance of the Gaussian Process. By default, we use a square exponential function ("SQexponential").

}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{mu }{Posterior mean vector on times ts}
  \item{GPSigma }{Posterior variance covariance matrix for time vector ts}
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
