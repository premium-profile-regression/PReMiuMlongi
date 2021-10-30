\name{mypheatDissMat}
\alias{mypheatDissMat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the heatmap of the dissimilarity matrix and compare to known clustering
}
\description{
Plot the heatmap of the dissimilarity matrix and compare to known clustering (for example true clustering in simulation studies.)
}
\usage{
mypheatDissMat(dissimObj, true_clusters = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dissimObj}{
An object of class dissimObj.
}
  \item{true_clusters}{
A vector containing the known clusters.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Plot of the heatmap of the dissimilary matrix. This functions uses the function 'heatmap' of package 'stats'. Note that this function has not been optimised for large datasets.
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
