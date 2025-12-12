
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mantar - Missingness Alleviation for NeTwork Analysis in R <img src="man/figures/sticker.png" align="right" height="138" alt="mantar sticker" /></a>

<!-- badges: start -->

[![CRAN_Status_Badge](https://r-pkg.org/badges/version/mantar)](https://cran.r-project.org/package=mantar)
[![Download_Badge](https://cranlogs.r-pkg.org/badges/grand-total/mantar)](https://www.datasciencemeta.com/rpackages)
[![R-CMD-check](https://github.com/kai-nehler/mantar/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kai-nehler/mantar/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`mantar` provides several methods for estimating psychological networks
with and without missing data. For network estimation, two main
approaches are implemented:

- Neighborhood selection with information-criteria-based model selection
- Regularized network estimation (e.g., glasso) combined with
  information criteria

For missing data handling, the preferred two strategies are:

- stacked multiple imputation
- two-step EM algorithm

The workflow is designed to support typical use cases in psychological
research, including varying sample sizes, missingness levels, and
variable types.

## Installation

The current stable version (0.2.0) is [available on
CRAN](https://cran.r-project.org/package=stuart) and can be installed
using the usual approach:

``` r
install.packages("mantar")
```

You can install the development version of `mantar` from
[GitHub](https://github.com/kai-nehler/mantar). To do so, you need the
`remotes` package.

``` r
# install.packages("remotes")
remotes::install_github("kai-nehler/mantar@develop")
#> Downloading GitHub repo kai-nehler/mantar@develop
#> 
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#> * checking for file ‘/tmp/RtmpfoYoh5/remotes6378386c42c7/kai-nehler-mantar-f23df72/DESCRIPTION’ ... OK
#> * preparing ‘mantar’:
#> * checking DESCRIPTION meta-information ... OK
#> * installing the package to process help pages
#> Loading required namespace: mantar
#> * saving partial Rd database
#> * checking for LF line-endings in source and make files and shell scripts
#> * checking for empty or unneeded directories
#> * building ‘mantar_0.2.0.tar.gz’
#> Installing package into '/tmp/Rtmp3nbteu/temp_libpath3fe5733c09c'
#> (as 'lib' is unspecified)
```

The suffix `@develop` ensures you install the latest development version
with new features and updates.

## Basic Example

`mantar` offers a few dummy data sets for demonstration and testing
purposes. Here, we illustrate how to estimate a neighborhood selection
network on a continuous data set with missing values using stacked
multiple imputation.

``` r
library(mantar)
data(mantar_dummy_mis_cont)
```

Network analysis based on neighborhood selection can be performed using
the `neighborhood_net()` function. By default, it uses the Bayesian
Information Criterion (BIC) for model selection.

``` r
result <- neighborhood_net(mantar_dummy_full_cont,
                           missing_handling = "stacked-mi")
#> No missing values in data. Sample size for each variable is equal to the number of rows in the data.
summary(result)
#> The density of the estimated network is 0.250
#> 
#> Network was estimated using neighborhood selection with the information criterion: BIC
#> and the 'and' rule for the inclusion of edges based on a full data set.
#> 
#> The sample sizes used for the nodewise regressions were as follows:
#> EmoReactivity     TendWorry    StressSens     SelfAware     Moodiness 
#>           400           400           400           400           400 
#>      Cautious ThoughtFuture RespCriticism 
#>           400           400           400
result$pcor
#>               EmoReactivity TendWorry StressSens SelfAware Moodiness  Cautious
#> EmoReactivity     0.0000000 0.2617524   0.130019 0.0000000 0.0000000 0.0000000
#> TendWorry         0.2617524 0.0000000   0.000000 0.2431947 0.0000000 0.0000000
#> StressSens        0.1300190 0.0000000   0.000000 0.0000000 0.0000000 0.0000000
#> SelfAware         0.0000000 0.2431947   0.000000 0.0000000 0.0000000 0.0000000
#> Moodiness         0.0000000 0.0000000   0.000000 0.0000000 0.0000000 0.4377322
#> Cautious          0.0000000 0.0000000   0.000000 0.0000000 0.4377322 0.0000000
#> ThoughtFuture     0.0000000 0.2595917   0.000000 0.0000000 0.0000000 0.0000000
#> RespCriticism     0.0000000 0.0000000   0.000000 0.0000000 0.2762595 0.2523658
#>               ThoughtFuture RespCriticism
#> EmoReactivity     0.0000000     0.0000000
#> TendWorry         0.2595917     0.0000000
#> StressSens        0.0000000     0.0000000
#> SelfAware         0.0000000     0.0000000
#> Moodiness         0.0000000     0.2762595
#> Cautious          0.0000000     0.2523658
#> ThoughtFuture     0.0000000     0.0000000
#> RespCriticism     0.0000000     0.0000000
```

The summary shows the density of the estimated network and details about
the estimation procedure. The partial correlation matrix can be accessed
via `result$pcor`.

## Further Documentation

The package vignette provides:

- network estimation workflows
- missing-data recommendations
- examples for both estimation approaches
