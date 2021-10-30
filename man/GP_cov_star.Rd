\name{GP_cov_star}
\alias{GP_cov_star}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes the covariance matrix of a GP (no measurement error) for two different time vectors
}
\description{
A function that computes the covariance matrix of a GP (no measurement error) for two different time vectors
}
\usage{
GP_cov_star(X1, X2, L, kernelType)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X1}{
A vector of times
}
  \item{X2}{
A vector of times
}
  \item{L}{
A vector of hyperparameters for the covariance function (e.g. variance signal, lengthscale for kernelType="SQexponential")
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
\item{Mat }{The covariance matrix of X1 and X2. Not necessarily squared matrix.}
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
