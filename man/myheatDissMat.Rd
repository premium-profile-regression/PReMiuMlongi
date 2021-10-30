\name{myheatDissMat}
\alias{myheatDissMat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the heatmap of the dissimilarity matrix with row reordering
}
\description{
Function to plot the heatmap of the dissimilarity matrix with row reordering
}
\usage{
myheatDissMat(dissimObj, main = NULL, xlab = NULL, ylab = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{dissimObj}{An object of class dissimObj.}
\item{main}{The usual plot option, to be passed to the heatmap function.}
\item{ylab}{The usual plot option, to be passed to the heatmap function.}
\item{xlab}{The usual plot option, to be passed to the heatmap function.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Plot of the heatmap of the dissimilary matrix. This functions uses the function 'heatmap' of package 'stats', with the parameter Rowv=NULL (Rowv=NA in heatDissMat to suppress any row dendrogram and reordering). Note that this function has not been optimised for large datasets.
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
