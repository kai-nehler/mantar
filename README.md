
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mantar - Missingness Alleviation for NeTwork Analysis in R

<!-- badges: start -->
<!-- badges: end -->

`mantar` provides users with several methods for handling missing data
in the context of network analysis. Currently, these methods are
specifically implemented for network estimation via neighborhood
selection using the Bayesian Information Criterion (BIC).

## Installation

You can install the development version of mantar from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("kai-nehler/mantar")
```

After installation the easiest way to get an overview of functions and
capabilities is to use `?mantar` to open the package help-file. You
could also read the rest of this README for an introduction and some
examples.

## Features

As already described, the package currently focuses on network
estimation using **neighborhood selection** with **information
criteria** for model selection in node-wise regressions. This
functionality is available for both complete and incomplete data.

For datasets with missing values, two modern missing approaches are
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

Both methods produce a correlation matrix that is then used to estimate
the network via node-wise regressions. It is also possible to compute
the correlation matrix using pairwise or listwise deletion. However,
these methods are generally not recommended, except in specific cases,
such as when data are missing completely at random and the proportion of
missingness is very small.

In addition to network estimation, the package also supports **stepwise
regression search** based on information criteria for a **single
dependent variable**. This regression search is available for both
complete and incomplete data and relies on the same two-step EM or
stacked MI procedures to handle missing values as the network analysis.
While both methods to handle missingness are expected to perform well in
this context, no specific simulation study has been conducted to compare
their effectiveness for single regression modeling, and thus their
relative strengths remain an open question.

## Example

The package includes two dummy datasets that resemble a typical
psychological dataset, where the number of observations is considerably
larger than the number of variables. Although the variables have
descriptive names, these are included solely to make the examples more
engaging - the data themselves are **fully synthetic**.

- `mantar_dummy_full`: Fully observed data (no missing values)
- `mantar_dummy_mis`: Data with missing values

These data sets are intended for examples and testing only.

``` r
library(mantar)

# Load example data
data(mantar_dummy_full)
data(mantar_dummy_mis)

# Preview the first few rows
head(mantar_dummy_full)
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
head(mantar_dummy_mis)
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

The main function for estimating a network is `neighborhood_net()`. In
the case of fully observed data, the function takes the dataset as input
and estimates a network structure using **neighborhood selection**
guided by **information criteria**. With default arguments, only the
dataset needs to be provided.

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

The `cor_method` argument specifies the type of correlation to be
estimated. Currently, only two types are supported in this context, as
these were evaluated in the given setting.

- `"pearson"`: computes the **Pearson correlation**.
- `"polychoric"`: computes the **polychoric correlation**.

The default for this argument is `"adapted"`, which automatically
selects an appropriate correlation type based on the data. If all
variables are continuous, `"pearson"` is used. If at least one variable
is ordinal, the procedure considers the number of variables, the number
of categories, and the sample size to decide whether to apply
`"pearson"` or `"polychoric"`. If you want to force a specific type of
correlation, you can set `cor_method` to either `"pearson"` or
`"polychoric"`.

### Example of Network Estimation without Missing Data

``` r
# Estimate network from full data set using BIC, the and rule as well as an adapted correlation type
result <- neighborhood_net(data = mantar_dummy_full, 
                           k = "log(n)", 
                           pcor_merge_rule = "and",
                           cor_method = "adapted")
#> Using the 'pearson' correlation method because the number of distinct values suggests continuous variables.
#> No missing values in data. Sample size for each variable is equal to the number of rows in the data.
# View estimated partial correlations
result
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
sum_result <- summary(result)
sum_result
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
- `"total"`: Uses the total number of observations in the dataset (i.e.,
  the number of rows).

#### Handling Missing Data

The `missing_handling` argument specifies how the correlation matrix is
estimated when the input data contains missing values. Two approaches
are supported:

- `"two-step-em"`: Applies a standard **Expectation-Maximization (EM)**
  algorithm to estimate the covariance matrix.
- `"stacked-mi"`: Applies **multiple imputation** to create several
  completed datasets, which are then stacked into a single dataset. A
  correlation matrix is computed from this stacked data.

As described previously, deletion techniques (listwise and pairwise) are
also available, but their use is not recommended. When `"two-step-em"`
is selected, the correlation matrix is always based on Pearson
correlations.

If `"stacked-mi"` is used, the `nimp` argument controls the number of
imputations (default: `20`), while `imp_method` specifies the imputation
method (default: `"pmm"` for predictive mean matching).

### Example of Network Estimation with Missing Data

``` r
# Estimate network for data set with missing values
result_mis <- neighborhood_net(data = mantar_dummy_mis, 
                                n_calc = "individual", 
                                missing_handling = "two-step-em",
                                pcor_merge_rule = "and",
                                cor_method = "adapted")
#> Using the 'pearson' correlation method because the number of distinct values suggests continuous variables.
# View estimated partial correlations
result_mis
#>               EmoReactivity TendWorry StressSens SelfAware Moodiness  Cautious
#> EmoReactivity     0.0000000 0.1295824   0.230612 0.0000000 0.0000000 0.0000000
#> TendWorry         0.1295824 0.0000000   0.000000 0.2515697 0.0000000 0.0000000
#> StressSens        0.2306120 0.0000000   0.000000 0.0000000 0.0000000 0.0000000
#> SelfAware         0.0000000 0.2515697   0.000000 0.0000000 0.0000000 0.0000000
#> Moodiness         0.0000000 0.0000000   0.000000 0.0000000 0.0000000 0.4768098
#> Cautious          0.0000000 0.0000000   0.000000 0.0000000 0.4768098 0.0000000
#> ThoughtFuture     0.1446363 0.2991518   0.000000 0.0000000 0.0000000 0.0000000
#> RespCriticism     0.0000000 0.0000000   0.000000 0.3008107 0.1930326 0.2210164
#>               ThoughtFuture RespCriticism
#> EmoReactivity     0.1446363     0.0000000
#> TendWorry         0.2991518     0.0000000
#> StressSens        0.0000000     0.0000000
#> SelfAware         0.0000000     0.3008107
#> Moodiness         0.0000000     0.1930326
#> Cautious          0.0000000     0.2210164
#> ThoughtFuture     0.0000000     0.0000000
#> RespCriticism     0.0000000     0.0000000

# Create and view a summary of the network estimation
sum_result_mis <- summary(result_mis)
sum_result_mis
#> The density of the estimated network is 0.321
#> 
#> Network was estimated using neighborhood selection on data with missing values.
#> Missing data were handled using 'two-step-em'.
#> The penalty term was log(n) and the 'and' rule was used for edge inclusion.
#> 
#> The sample sizes used for the nodewise regressions were as follows:
#> EmoReactivity     TendWorry    StressSens     SelfAware     Moodiness 
#>           427           426           425           428           424 
#>      Cautious ThoughtFuture RespCriticism 
#>           423           422           420
```
