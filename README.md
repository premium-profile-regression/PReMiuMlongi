
# PReMiuMlongi

<!-- badges: start -->
[![R-CMD-check](https://github.com/premium-profile-regression/PReMiuMlongi/workflows/R-CMD-check/badge.svg)](https://github.com/premium-profile-regression/PReMiuMlongi/actions)
<!-- badges: end -->

`PReMiuMlongi` is an extension of the R `PReMiuM` package (Liverani et al., 2015), implementing Bayesian profile regression for longitudinal Gaussian response and covariates.

## Installation

You can install `PReMiuMlongi` from [GitHub](https://CRAN.R-project.org) with the following commands:
``` r
if(!requireNamespace("remotes", quietly=TRUE)){
  install.packages(remotes)
}
if(!requireNamespace("premiumPlots", quietly=TRUE)){
  remotes::install_github("simisc/premiumPlots")
}
remotes::install_github("premium-profile-regression/PReMiuMlongi")
```

**NB:** the `premiumPlots` package which `PReMiuMlongi` depends on is only available from GitHub.

