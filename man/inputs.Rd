\name{inputs}
\alias{inputs}
\docType{data}
\title{
Inputs data for toy example
}
\description{
Inputs data for toy example
}
\usage{data("inputs")}
\format{
  The format is:
List of 6
 $ inputData  :'data.frame':	100 obs. of  9 variables:
  ..$ ID       : num [1:100] 1 2 3 4 5 6 7 8 9 10 ...
  ..$ outcome1 : num [1:100] -0.0115 0.8149 0.382 0.8379 3.462 ...
  ..$ outcome2 : num [1:100] -0.5 1.374 -0.104 0.848 2.379 ...
  ..$ outcome3 : num [1:100] 0.181 1.582 -0.597 0.787 -0.68 ...
  ..$ outcome4 : num [1:100] 1.152 0.826 -0.221 2.006 -0.589 ...
  ..$ Variable1: num [1:100] 0 2 2 2 2 0 0 0 2 2 ...
  ..$ Variable2: num [1:100] 1 0 0 0 0 1 1 1 0 0 ...
  ..$ Variable3: num [1:100] 0 0 0 0 0 0 0 0 0 0 ...
  ..$ Variable4: num [1:100] 1 1 1 1 1 1 1 1 1 0 ...
 $ covNames   : chr [1:4] "Variable1" "Variable2" "Variable3" "Variable4"
 $ outcome    : chr [1:4] "outcome1" "outcome2" "outcome3" "outcome4"
 $ xModel     : chr "Discrete"
 $ yModel     : chr "MVN"
 $ nCovariates: num 4
}
\examples{
data(inputs)
## maybe str(inputs) ; plot(inputs) ...
}
\keyword{datasets}
