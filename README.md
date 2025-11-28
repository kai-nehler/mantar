
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mantar - Missingness Alleviation for NeTwork Analysis in R <img src="man/figures/sticker.png" align="right" height="138" alt="mantar sticker" /></a>

<!-- badges: start -->

[![CRAN_Status_Badge](https://r-pkg.org/badges/version/mantar)](https://cran.r-project.org/package=mantar)
[![Download_Badge](https://cranlogs.r-pkg.org/badges/grand-total/mantar)](https://www.datasciencemeta.com/rpackages)
[![R-CMD-check](https://github.com/kai-nehler/mantar/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kai-nehler/mantar/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`mantar` provides users with several methods for handling missing data
in the context of network analysis. The package focuses on estimating
psychological networks using **neighborhood selection** with
**information criteria** for model selection in node-wise regressions.
Furthemore, as of version 0.2.0, the package also supports network
estimation based on **regularization techniques** in combination with
**information criteria** for model selection. Both functionalities are
available for datasets with and without missing values.

## Installation

The current stable version (0.2.0) is [available on
CRAN](https://cran.r-project.org/package=stuart) and can be installed
using the usual approach:

``` r
install.packages("mantar")
```

You can install the development version of `mantar` from
[GitHub](https://github.com/kai-nehler/mantar). Before installing, make
sure your system is able to compile packages from source. You can verify
this using:

``` r
pkgbuild::check_build_tools()
#> Your system is ready to build packages!
```

If build tools are missing or outdated, you may need to install/update
appropriate development tools (e.g. Rtools for Windows) for your
operating system.

Once your system is ready, install the development branch via:

``` r
# install.packages("pak")
pak::pak("kai-nehler/mantar@develop")
#> ! Using bundled GitHub PAT. Please add your own PAT using `gitcreds::gitcreds_set()`.
#> 
#> ✔ Updated metadata database: 3.25 MB in 8 files.
#> 
#> ℹ Updating metadata database
#> ✔ Updating metadata database ... done
#> 
#> 
#> → Will install 4 packages.
#> → Will update 1 package.
#> → Will download 4 CRAN packages (2.62 MB).
#> → Will download 1 package with unknown size.
#> + glassoFast         1.0.1 [bld][cmp][dl] (5.43 kB)
#> + mantar     0.2.0 → 0.2.0 [bld][cmp][dl] (GitHub: bd97012)
#> + mathjaxr           1.8-0 [bld][cmp][dl] (1.06 MB)
#> + rbibutils          2.4   [bld][cmp][dl] (1.17 MB)
#> + Rdpack             2.6.4 [bld][dl] (379.14 kB)
#> 
#> ℹ Getting 4 pkgs (2.62 MB) and 1 pkg with unknown size
#> ✔ Got glassoFast 1.0.1 (source) (5.43 kB)
#> ✔ Got rbibutils 2.4 (source) (1.17 MB)
#> ✔ Got Rdpack 2.6.4 (source) (379.14 kB)
#> ✔ Got mathjaxr 1.8-0 (source) (1.06 MB)
#> ✔ Got mantar 0.2.0 (source) (240.00 kB)
#> ℹ Building glassoFast 1.0.1
#> ℹ Building mathjaxr 1.8-0
#> ℹ Building rbibutils 2.4
#> ✔ Built glassoFast 1.0.1 (1.1s)
#> ✔ Installed glassoFast 1.0.1  (1s)
#> ✔ Built mathjaxr 1.8-0 (2.2s)
#> ✔ Installed mathjaxr 1.8-0  (34ms)
#> ✔ Built rbibutils 2.4 (10.4s)
#> ✔ Installed rbibutils 2.4  (38ms)
#> ℹ Building Rdpack 2.6.4
#> ✔ Built Rdpack 2.6.4 (2.3s)
#> ✔ Installed Rdpack 2.6.4  (19ms)
#> ℹ Packaging mantar 0.2.0
#> ✔ Packaged mantar 0.2.0 (2.1s)
#> ℹ Building mantar 0.2.0
#> ✔ Built mantar 0.2.0 (1.1s)
#> ✔ Installed mantar 0.2.0 (github::kai-nehler/mantar@bd97012) (1s)
#> ✔ 1 pkg + 6 deps: kept 1, upd 1, added 4, dld 5 (NA B) [25.3s]
```

The extension `@develop` ensures that you get the latest development
version of the package, which may include new features and bug fixes not
yet available in the stable release on CRAN. Exluding this extension
will install the same version as the one on CRAN.

After installation the easiest way to get an overview of functions and
capabilities is to use `help(package = "mantar")` to open the package
help-file. You could also read the rest of this README for an
introduction and some examples.

## Features

As described above, the package currently implements two of the many
approaches for estimating network structures:

- **Neighborhood selection**: This approach is based on node-wise
  regressions with information criteria for model selection. Each node
  is used once as a dependent variable, while all other nodes serve as
  potential predictors. For each node, the best set of predictors is
  selected using an information criterion. The resulting regression
  weights are then combined to obtain partial correlations between
  nodes.
- **Regularization**: This approach relies on penalized maximum
  likelihood estimation of the (inverse) covariance matrix. Through
  appropriate penalty parameters, some entries of the estimated inverse
  covariance matrix are shrunk to zero, which in turn yields zeros in
  the partial correlation matrix and therefore a sparse network
  structure. Multiple values of the penalty parameter are evaluated, and
  information criteria are used to select the optimal penalization.

For data sets with missing values, two modern missing approaches are
implemented:

- **Two-step Expectation-Maximization (EM):** A fast method that
  estimates the correlation matrix via an EM algorithm using the
  `lavaan` package. It performs well when the sample size is very large
  relative to the amount of missingness and the complexity of the
  network.
- **Stacked Multiple Imputation (MI):** A more robust approach across a
  wider range of sample sizes. Multiple imputation is performed using
  predictive mean matching (PMM) with the `mice` package. The imputed
  data sets are stacked into a single data set, and a correlation matrix
  is estimated from this combined data.

Both methods produce a correlation matrix that is then used for network
estimation. It is also possible to compute the correlation matrix using
pairwise or listwise deletion. However, these methods are generally not
recommended, except in specific cases (e.g., when data are missing
completely at random and the proportion of missingness is very small).
By default, correlations are computed using Pearson correlations.
However, with complete data, listwise deletion, or the stacked MI
approach, users may choose to treat variables as ordered categorical, in
which case polychoric and polyserial correlations are computed where
appropriate. This option is particularly advisable when variables have a
low number of categories or exhibit noticeable non-normality. At the
same time, estimating polychoric and polyserial correlations requires a
sufficiently large number of observations relative to the number of
variables to ensure stable and reliable estimates.

In addition to network estimation, the package also supports **stepwise
regression search** based on information criteria for a **single
dependent variable**. This regression search is available for both
complete and incomplete data and relies on the same two-step EM or
stacked MI procedures to handle missing values as the network analysis.
While both methods to handle missingness are expected to perform well in
this context, no specific simulation study has been conducted to compare
their effectiveness for single regression modeling, and thus their
relative strengths remain an open question.

## Examples

The package includes dummy datasets that resemble a typical
psychological dataset, where the number of observations is considerably
larger than the number of variables. Although the variables have
descriptive names, these are included solely to make the examples more
engaging - the data themselves are **fully synthetic**.

Three data sets without missing values are included: -
`mantar_dummy_full_cont`: Fully observed data (no missing values) -
`mantar_dummy_full_cat`: Fully observed data with ordered categorical
variables - `mantar_dummy_full_mix`: Fully observed data with a mix of
continuous and ordered categorical variables

Additionally, three data sets with missing values are provided: -
`mantar_dummy_mis_cont`: Data with approximately 30% missing values in
each continuous variable - `mantar_dummy_mis_cat`: Data with
approximately 25% missing values in each ordered categorical variable -
`mantar_dummy_mis_mix`: Data with approximately 25% missing values in
each variable, with a mix of continuous and ordered categorical
variables

These data sets are intended for examples and testing only.

``` r
library(mantar)

# Load some example data sets for the ReadMe
data(mantar_dummy_full_cont)
data(mantar_dummy_full_cat)
data(mantar_dummy_mis_cont)

# Preview the first few rows of these data sets
head(mantar_dummy_full_cont)
#>   EmoReactivity  TendWorry StressSens  SelfAware  Moodiness    Cautious
#> 1   -0.08824641 -0.2659269 -1.2036137 -2.3499259  0.6693700  0.04102854
#> 2   -0.44657803 -0.4588384 -0.2431794 -0.1656722 -0.3361568  0.88919849
#> 3   -1.06934325 -1.5050242 -0.8986388 -1.0857552  0.2249633  0.77060142
#> 4    0.58282029 -0.5036316 -1.6020000  1.0820676 -0.1858346 -0.03462852
#> 5    0.58791759  0.5972580 -0.5882332  1.7461103  0.7160714  1.58280444
#> 6    0.10224725  0.1494428 -1.0877812 -1.7886107  1.3522197 -0.25494638
#>   ThoughtFuture RespCriticism
#> 1     0.6484939   -0.77992262
#> 2     0.2949630   -0.91747608
#> 3    -1.3519007    0.56000763
#> 4    -0.4702988    0.34653985
#> 5     0.9503597    0.82981174
#> 6    -0.8938618   -0.01593388
head(mantar_dummy_full_cat)
#>   EmoReactivity TendWorry StressSens SelfAware Moodiness Cautious ThoughtFuture
#> 1             3         3          2         1         4        3             4
#> 2             3         3          3         3         3        4             4
#> 3             2         2          3         2         4        4             2
#> 4             4         3          2         5         3        3             3
#> 5             4         4          3         5         4        5             4
#> 6             4         4          2         2         5        3             3
#>   RespCriticism
#> 1             3
#> 2             3
#> 3             4
#> 4             4
#> 5             4
#> 6             3
head(mantar_dummy_mis_cont)
#>   EmoReactivity  TendWorry StressSens   SelfAware  Moodiness    Cautious
#> 1    -1.7551632 -0.4376210 -0.5774722  0.10562820  0.6614044          NA
#> 2    -1.7551688 -0.7039623  0.9070330  0.03418623  0.6140406  0.83879818
#> 3     2.0493638         NA         NA          NA -0.8872971  0.04830719
#> 4     0.1056282         NA         NA -1.24779117 -0.7298623 -0.62263184
#> 5    -0.6338512  0.4361078 -0.5564631 -0.01032403         NA -0.09690612
#> 6     0.1054382  0.6935808  2.6557231          NA         NA -0.04358574
#>   ThoughtFuture RespCriticism
#> 1     0.7710993    0.37233355
#> 2    -1.5588119   -0.55079199
#> 3            NA   -0.90103222
#> 4    -0.7100126    0.80773402
#> 5     1.0583312    0.20820252
#> 6            NA   -0.03915726
```

First, we will present an example of network estimation using the
`neighborhood_net()` function for a data set without missing values.
Following that, we will illustrate how to estimate a network when the
data contains missing values. Then, we will briefly demonstrate how to
perform a network analysis using regularization techniques with the
`regularization_net()` function.

### Network Estimation via Neighborhood Selection

The `neighborhood_net()` function estimates a network structure based on
neighborhood selection using information criteria for model selection in
node-wise regressions. The function can either be provided with raw data
(data frame or matrix) or a correlation matrix along with sample sizes
for each variable. The examples will use raw data, as this is the more
complex case. The following arguments are particularly relevant for
controlling the network estimation process (with fully observed data):

#### Information Criteria

The `k` argument controls the penalty used in model selection for
node-wise regressions. It reflects the penalty per parameter (i.e.,
number of predictors + 1):

- `k = "log(n)"` (default): corresponds to the **Bayesian Information
  Criterion (BIC)**
- `k = "2"`: corresponds to the **Akaike Information Criterion (AIC)**

#### Estimation of Partial Correlation

The `pcor_merge_rule` argument determines how partial correlations are
estimated based on the regression results between two nodes:

- `"and"` (default): a partial correlation is estimated **only if both**
  regression weights (from node A to B and from B to A) are non-zero.
- `"or"`: a partial correlation is estimated if **at least one** of the
  two regression weights is non-zero.

Although both options are available, current simulation evidence
suggests that the `"and"` rule yields more accurate partial correlation
estimates than the `"or"` rule. Therefore, changing this default is
**not recommended** unless you have a specific reason.

#### Type of Correlation

The `ordered` argument specifies how variables are treated when
estimating correlations from raw data.

- **Global specification:**
  - `ordered = TRUE`: all variables are treated as ordered categorical  
  - `ordered = FALSE`: all variables are treated as continuous  
- **Variable-specific specification:**
  - A logical vector of length equal to the number of variables can be
    supplied to indicate which variables are treated as ordered
    categorical (e.g., `ordered = c(TRUE, FALSE, FALSE, TRUE)`).  
- If a global specification is used, the function automatically creates
  a vector of this value with the same length as the number of
  variables.

Based on these specifications, the function applies the appropriate
correlation type for each pair of variables:  
- both `FALSE`: **Pearson correlation**  
- one `TRUE` and one `FALSE`: **polyserial correlation**  
- both `TRUE`: **polychoric correlation**

#### Estimation with Continuous Complete Data

After discussing the key arguments, we can now illustrate how to
estimate a network structure using the `neighborhood_net()` function
with a data set without missing values.

``` r
# Estimate network from full data set using BIC, the and rule as well as treating the data as continuous
result_full_cont <- neighborhood_net(data = mantar_dummy_full_cont, 
                                             k = "log(n)", 
                           pcor_merge_rule = "and",
                           ordered = FALSE)
#> No missing values in data. Sample size for each variable is equal to the number of rows in the data.

# View estimated partial correlations
result_full_cont
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

# Create and view a summary of the network estimation
sum_result_full_cont <- summary(result_full_cont)
sum_result_full_cont
#> The density of the estimated network is 0.250
#> 
#> Network was estimated using neighborhood selection with a penalty term of log(n)
#> and the 'and' rule for the inclusion of edges based on a full data set.
#> 
#> The sample sizes used for the nodewise regressions were as follows:
#> EmoReactivity     TendWorry    StressSens     SelfAware     Moodiness 
#>           400           400           400           400           400 
#>      Cautious ThoughtFuture RespCriticism 
#>           400           400           400
```

#### Estimation with Ordered Categorical Complete Data

We can also estimate a network structure when some variables are ordered
categorical. In the following example, we treat all variables as ordered
categorical.

``` r
# Estimate network from full data set using BIC, the and rule as well as treating the
# data as ordered categorical
result_full_cat <- neighborhood_net(data = mantar_dummy_full_cat,
                                    k = "log(n)", 
                                    pcor_merge_rule = "and",
                                    ordered = TRUE)
#> No missing values in data. Sample size for each variable is equal to the number of rows in the data.

# View estimated partial correlations
result_full_cat
#>               EmoReactivity TendWorry StressSens SelfAware Moodiness  Cautious
#> EmoReactivity     0.0000000 0.2742356   0.136029 0.0000000 0.0000000 0.0000000
#> TendWorry         0.2742356 0.0000000   0.000000 0.2679285 0.0000000 0.0000000
#> StressSens        0.1360290 0.0000000   0.000000 0.0000000 0.0000000 0.0000000
#> SelfAware         0.0000000 0.2679285   0.000000 0.0000000 0.0000000 0.0000000
#> Moodiness         0.0000000 0.0000000   0.000000 0.0000000 0.0000000 0.4398609
#> Cautious          0.0000000 0.0000000   0.000000 0.0000000 0.4398609 0.0000000
#> ThoughtFuture     0.0000000 0.2224662   0.000000 0.0000000 0.0000000 0.0000000
#> RespCriticism     0.0000000 0.0000000   0.000000 0.0000000 0.2752566 0.2687388
#>               ThoughtFuture RespCriticism
#> EmoReactivity     0.0000000     0.0000000
#> TendWorry         0.2224662     0.0000000
#> StressSens        0.0000000     0.0000000
#> SelfAware         0.0000000     0.0000000
#> Moodiness         0.0000000     0.2752566
#> Cautious          0.0000000     0.2687388
#> ThoughtFuture     0.0000000     0.0000000
#> RespCriticism     0.0000000     0.0000000

# Create and view a summary of the network estimation
sum_result_full_cat <- summary(result_full_cat)
sum_result_full_cat
#> The density of the estimated network is 0.250
#> 
#> Network was estimated using neighborhood selection with a penalty term of log(n)
#> and the 'and' rule for the inclusion of edges based on a full data set.
#> 
#> The sample sizes used for the nodewise regressions were as follows:
#> EmoReactivity     TendWorry    StressSens     SelfAware     Moodiness 
#>           400           400           400           400           400 
#>      Cautious ThoughtFuture RespCriticism 
#>           400           400           400
```

The standard output does not differ between the continuous and ordered
categorical cases.

In the case of missing data, the `neighborhood_net()` function offers
several additional arguments that control how sample size and
missingness are handled.

#### Calculation of Sample Size

The `n_calc` argument specifies how the sample size is calculated for
each node-wise regression. This affects the penalty term used in model
selection.

The available options are:

- `"individual"` *(default)*: Uses the number of non-missing
  observations for each individual variable. This is the recommended
  approach.
- `"average"`: Uses the average number of non-missing observations
  across all variables.
- `"max"`: Uses the maximum number of non-missing observations across
  all variables.
- `"total"`: Uses the total number of observations in the data set
  (i.e., the number of rows).

#### Handling Missing Data

The `missing_handling` argument specifies how the correlation matrix is
estimated when the input data contains missing values. Two approaches
are supported:

- `"two-step-em"`: Applies a standard **Expectation-Maximization (EM)**
  algorithm to estimate the covariance matrix. This method is the
  default as it is computationally efficient. However, it only performs
  well when the sample size is large relative to the amount of edges in
  the network and the proportion of missingness.
- `"stacked-mi"`: Applies **multiple imputation** to create several
  completed datasets, which are then stacked into a single dataset. A
  correlation matrix is computed from this stacked data.

As described previously, deletion techniques (listwise and pairwise) are
also available, but their use is not recommended. When `"two-step-em"`
is selected, the correlation matrix is always based on Pearson
correlations, regardless of the `ordered` argument. In contrast, when
`"stacked-mi"` is used, the `ordered` argument determines how variables
are treated (continuous vs. ordered categorical) during the correlation
estimation.

If `"stacked-mi"` is used, the `nimp` argument controls the number of
imputations (default: `20`), while `imp_method` specifies the imputation
method (default: `"pmm"` for predictive mean matching).

#### Estimation with Missing Data

We can now illustrate how to estimate a network structure using the
`neighborhood_net()` function with a data set that contains missing
values. All variables are continuous in this example.

``` r
# Estimate network for data set with missing values
result_mis_cont <- neighborhood_net(data = mantar_dummy_mis_cont,
                                    n_calc = "individual", 
                                    missing_handling = "two-step-em",
                                    pcor_merge_rule = "and")

# View estimated partial correlations
result_mis_cont
#>               EmoReactivity TendWorry StressSens SelfAware Moodiness  Cautious
#> EmoReactivity     0.0000000 0.0000000  0.1947561 0.0000000 0.0000000 0.0000000
#> TendWorry         0.0000000 0.0000000  0.0000000 0.2769383 0.0000000 0.1487042
#> StressSens        0.1947561 0.0000000  0.0000000 0.0000000 0.0000000 0.0000000
#> SelfAware         0.0000000 0.2769383  0.0000000 0.0000000 0.0000000 0.0000000
#> Moodiness         0.0000000 0.0000000  0.0000000 0.0000000 0.0000000 0.4328186
#> Cautious          0.0000000 0.1487042  0.0000000 0.0000000 0.4328186 0.0000000
#> ThoughtFuture     0.2210913 0.2603475  0.0000000 0.0000000 0.0000000 0.0000000
#> RespCriticism     0.0000000 0.0000000  0.0000000 0.2991248 0.2782289 0.1809261
#>               ThoughtFuture RespCriticism
#> EmoReactivity     0.2210913     0.0000000
#> TendWorry         0.2603475     0.0000000
#> StressSens        0.0000000     0.0000000
#> SelfAware         0.0000000     0.2991248
#> Moodiness         0.0000000     0.2782289
#> Cautious          0.0000000     0.1809261
#> ThoughtFuture     0.0000000     0.0000000
#> RespCriticism     0.0000000     0.0000000

# Create and view a summary of the network estimation
sum_result_mis_cont <- summary(result_mis_cont)
sum_result_mis_cont
#> The density of the estimated network is 0.321
#> 
#> Network was estimated using neighborhood selection on data with missing values.
#> Missing data were handled using 'two-step-em'.
#> The penalty term was log(n) and the 'and' rule was used for edge inclusion.
#> 
#> The sample sizes used for the nodewise regressions were as follows:
#> EmoReactivity     TendWorry    StressSens     SelfAware     Moodiness 
#>           273           276           281           284           287 
#>      Cautious ThoughtFuture RespCriticism 
#>           281           281           282
```

The output does not differ much from the complete data cased, only
offering additional information about how missingness was handled.
Additionally, due to the `n_calc` argument being set to `individual`,
the sample sizes for each variable may differ.

### Network Estimation via Regularization

The `regularization_net()` function estimates a network structure based
on regularization techniques using information criteria for model
selection. Similar to `neighborhood_net()`, this function can either be
provided with raw data (data frame or matrix) or a correlation matrix
along with sample sizes for each variable. The examples will use raw
data, as this is the more complex case. The following arguments are
particularly relevant for controlling the network estimation process
(with fully observed data):

#### Type of Regularization Penalty and corresponding Parameters

The `penalty` argument controls the type of regularization used in the
network estimation. The recommended options are using the graphical
lasso (`"glasso"`) as a convex penalty or `"atan"` as a non-convex
penalty.

For `glasso`, the `lambda_min_ratio` and `n_lambdas` arguments control
the range and number of penalty parameters evaluated during model
selection. The default values are generally appropriate. For all
nonconvex penalties (e.g., `"atan"`), there is the option to specify an
additional parameter via the `gamma` argument. However, just using one
default value (default value is different for different penalty types)
for `gamma`, by setting the argument `vary = "lambda"` is sufficient.

The last argument controlling the regularization process is `pen_diag`
which specifies whether the diagonal elements of the covariance matrix
should be penalized (`TRUE`) or not (`FALSE`, default).

#### Information Criteria

The `k` argument controls the penalty used in model selection for the
regularization process. It reflects the penalty per parameter (in this
case for every included edge / nonzero partial correlation). The default
value is `k = "log(n)"`, which corresponds to the **Bayesian Information
Criterion (BIC)**. Another option is `k = "2"`, which corresponds to the
**Akaike Information Criterion (AIC)**. The function also offers using
the extended version of the information criteria via `extended = TRUE`
(especially helpful to specify extended BIC). In this case, an
additional parameter `extended_gamma` can be specified (default: `0.5`).
The default of `extended` depends on the penalty type as the `"glasso"`
penalty generally does generally use the extended version, while
non-convex penalties do not.

#### Type of Correlation

The `ordered` argument specifies how variables are treated when
estimating correlations from raw data.

- **Global specification:**
  - `ordered = TRUE`: all variables are treated as ordered categorical  
  - `ordered = FALSE`: all variables are treated as continuous  
- **Variable-specific specification:**
  - A logical vector of length equal to the number of variables can be
    supplied to indicate which variables are treated as ordered
    categorical (e.g., `ordered = c(TRUE, FALSE, FALSE, TRUE)`).  
- If a global specification is used, the function automatically creates
  a vector of this value with the same length as the number of
  variables.

Based on these specifications, the function applies the appropriate
correlation type for each pair of variables:  
- both `FALSE`: **Pearson correlation**  
- one `TRUE` and one `FALSE`: **polyserial correlation**  
- both `TRUE`: **polychoric correlation**

#### Estimation with Continuous Complete Data

After discussing the key arguments, we can now illustrate how to
estimate a network structure using the `regularization_net()` function
with a data set without missing values.

``` r
# Estimate network from full data set using BIC and the glasso penalty
result_full_cont <- regularization_net(data = mantar_dummy_full_cont,
                                      penalty = "glasso",
                                      vary = "lambda",
                                      n_lambda = 100,
                                      lambda_min_ratio = 0.1,
                                      k = "log(n)",
                                      extended = TRUE, 
                                      extended_gamma = 0.5,
                                      pcor_merge_rule = "and",
                                      ordered = FALSE)

# View estimated partial correlations
result_full_cont
#>               EmoReactivity TendWorry  StressSens  SelfAware  Moodiness
#> EmoReactivity    0.00000000 0.1947651  0.05726695 0.00000000 0.00000000
#> TendWorry        0.19476514 0.0000000  0.00000000 0.17447288 0.00000000
#> StressSens       0.05726695 0.0000000  0.00000000 0.00000000 0.00000000
#> SelfAware        0.00000000 0.1744729  0.00000000 0.00000000 0.00000000
#> Moodiness        0.00000000 0.0000000  0.00000000 0.00000000 0.00000000
#> Cautious         0.06142038 0.0000000  0.00000000 0.00000000 0.38632297
#> ThoughtFuture    0.07236467 0.2030423  0.00000000 0.00000000 0.01130336
#> RespCriticism    0.00000000 0.0000000 -0.01643745 0.04990513 0.24301446
#>                 Cautious ThoughtFuture RespCriticism
#> EmoReactivity 0.06142038    0.07236467    0.00000000
#> TendWorry     0.00000000    0.20304232    0.00000000
#> StressSens    0.00000000    0.00000000   -0.01643745
#> SelfAware     0.00000000    0.00000000    0.04990513
#> Moodiness     0.38632297    0.01130336    0.24301446
#> Cautious      0.00000000    0.05748129    0.22087095
#> ThoughtFuture 0.05748129    0.00000000    0.00000000
#> RespCriticism 0.22087095    0.00000000    0.00000000

# Create and view a summary of the network estimation
sum_result_full_cont <- summary(result_full_cont)
sum_result_full_cont
#> The density of the estimated network is 0.464
#> 
#> Network was estimated using regularization with the glasso penalty. The observed-data log-likelihood was used in the model selection.
#> The effective sample size used for the information criteria
#>                 computation was 400.
```

With missing data, the `regularization_net()` function offers several
additional arguments that control how sample size, information criteria
computation and missingness are handled.

#### Calculation of Sample Size for Information Criteria

The `n_calc` argument specifies how the sample size is calculated for
the information criteria computation. Only one value is needed here, as
the regularization approach does not rely on node-wise regressions. The
default input is `"average"`, which uses the average number of
non-missing observations for all estimated correlations - this includes
the correlations of variables with themselves. Ignoring these
correlations (i.e., using the average number of non-missing observations
across different variables only) is also possible with setting
`count_diagonal` to `FALSE`.

Within the information criteria computation, the likelihood for the
candidate models has to be computed. The `likelihood` argument controls
how this is done:

- `"mat_based"` (default): The likelihood is computed based on the
  sample correlation matrix.
- `"obs_based"`: The likelihood is computed based on the observed data.
  This option is only available when the raw input data contains no
  ordered categorical variables. In these cases, the observed data
  log-likelihood is recommended as it is a better representation of the
  sample data than using the sample correlation matrix.

These options to compute the likelihood are also available with full
data. However, they return the exact same results in this case.

#### Handling Missing Data

The `missing_handling` argument specifies how the correlation matrix is
estimated when the input data contains missing values. Two approaches
are supported:

- `"two-step-em"`: Applies a standard **Expectation-Maximization (EM)**
  algorithm to estimate the covariance matrix. This method is the
  default as it is computationally efficient. However, it only performs
  well when the sample size is large relative to the amount of edges in
  the network and the proportion of missingness.
- `"stacked-mi"`: Applies **multiple imputation** to create several
  completed datasets, which are then stacked into a single data set. A
  correlation matrix is computed from this stacked data.

As described previously, deletion techniques (listwise and pairwise) are
also available, but their use is not recommended. When `"two-step-em"`
is selected, the correlation matrix is always based on Pearson
correlations, regardless of the `ordered` argument. In contrast, when
`"stacked-mi"` is used, the `ordered` argument determines how variables
are treated (continuous vs. ordered categorical) during the correlation
estimation.

If `"stacked-mi"` is used, the `nimp` argument controls the number of
imputations (default: `20`), while `imp_method` specifies the imputation
method (default: `"pmm"` for predictive mean matching).

#### Estimation with Missing Data

We can now illustrate how to estimate a network structure using the
`regularization_net()` function with a data set that contains missing
values. All variables are continuous in this example.

``` r
# Estimate network for data set with missing values
result_mis_cont <- regularization_net(data = mantar_dummy_mis_cont,
                                      likelihood = "obs_based",
                                     penalty = "glasso",
                                     vary = "lambda",
                                     n_lambda = 100,
                                     lambda_min_ratio = 0.1,
                                     k = "log(n)",
                                     extended = TRUE,
                                     extended_gamma = 0.5,
                                     n_calc = "average",
                                     missing_handling = "two-step-em",
                                     pcor_merge_rule = "and",
                                     ordered = FALSE)

# View estimated partial correlations
result_mis_cont
#>               EmoReactivity  TendWorry StressSens  SelfAware Moodiness
#> EmoReactivity    0.00000000 0.00000000  0.0866789 0.03103689 0.0000000
#> TendWorry        0.00000000 0.00000000  0.0000000 0.19973839 0.0000000
#> StressSens       0.08667890 0.00000000  0.0000000 0.00000000 0.0000000
#> SelfAware        0.03103689 0.19973839  0.0000000 0.00000000 0.0000000
#> Moodiness        0.00000000 0.00000000  0.0000000 0.00000000 0.0000000
#> Cautious         0.00000000 0.07321694  0.0000000 0.00000000 0.3645097
#> ThoughtFuture    0.12858049 0.18479199  0.0000000 0.03402553 0.0000000
#> RespCriticism    0.00000000 0.00000000  0.0000000 0.22080386 0.2346446
#>                 Cautious ThoughtFuture RespCriticism
#> EmoReactivity 0.00000000    0.12858049     0.0000000
#> TendWorry     0.07321694    0.18479199     0.0000000
#> StressSens    0.00000000    0.00000000     0.0000000
#> SelfAware     0.00000000    0.03402553     0.2208039
#> Moodiness     0.36450971    0.00000000     0.2346446
#> Cautious      0.00000000    0.00000000     0.1504089
#> ThoughtFuture 0.00000000    0.00000000     0.0000000
#> RespCriticism 0.15040891    0.00000000     0.0000000

# Create and view a summary of the network estimation
sum_result_mis_cont <- summary(result_mis_cont)
sum_result_mis_cont
#> The density of the estimated network is 0.393
#> 
#> Network was estimated using regularization with the glasso penalty. The observed-data log-likelihood was used in the model selection.
#> The effective sample size used for the information criteria
#>                 computation was 207.140625. 
#> Missing data were handled using 'two-step-em'.
```
