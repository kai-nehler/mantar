# mantar 0.3.1

## Refinements
* Do not include the README in the installation.
* Numerical tolerances in unit tests adjusted to avoid platform-specific test failures.
* Regenerating documentation was with the current roxygen2 version.

# mantar 0.3.0

## Bug fixes
* Fixed issue where off-diagonal elements were counted twice when computing the average sample size.

## New defaults
* When both `mat` and `data` are provided, `mat` is used for network estimation, while `data` is only used for additional calculations (e.g., for compatibility with `bootnet`).
* `cor_calc()` now automatically scales raw data when no variables are treated as ordered. This was already the case when called indirectly, but now also applies when calling it directly. Only affects the returned means, not the estimated correlations.

## New arguments and logic
* Added `means` argument to `regularization_net()`. Not intended for typical use, but needed for integration with `bootnet`.
* `ns` no longer accepts vectors in `regularization_net()`, instead accepts a matrix.
* Added `network_vars` and `auxiliary_vars` arguments to `cor_calc()`, `neighborhood_net()`, and `regularization_net()`. These arguments accept both variable names and numeric indices referring to the provided data.

## New features
* Functions now return imputed data as a `mids` object when using stacked multiple imputation.
* With the new `network_vars` and `auxiliary_vars` arguments, users can now specify which variables are used for network estimation and which additional variables are included for correlation estimation in the presence of missing data.

## Refinements
* Using only a single value for one of the regularization parameters in nonconvex regularization now results in a message rather than a warning.
* Updated all dummy data sets, including the complete versions, and introduced varying amounts of missingness across variables in the incomplete data sets.

# mantar 0.2.0

## Change of argument names and input
* The argument `k`, which controlled the penalty term in information-criterion calculations, has been removed for security reasons. Instead, the penalty type is now specified via the argument `ic_type` (see the corresponding help pages).

## New defaults 
* Changed the default handling of the `ridge` penalty in the multiple-imputation `pmm` workflow when looking for donors through regressions (`mice`). Instead of forcing `ridge = 0`, the function now uses the default value defined by `mice`, ensuring consistent and method-appropriate regularization.

## New features
* Added function `reg_network()` for network estimation using regularization, supporting both convex and non-convex penalties as well as multiple options for computing the likelihood in the information criterion when missing values are present.
* Added function `ordered_suggest()`, a heuristic procedure for identifying variables that may be treated as ordered categorical based on their distribution and available information.
* Added dummy data sets `mantar_dummy_full_cat` and `mantar_dummy_mis_cat`, containing only ordered categorical variables (with and without missing values).
* Added dummy data sets `mantar_dummy_full_mix` and `mantar_dummy_mis_mix`, containing mixtures of ordered categorical and continuous variables (with and without missing values).

## Improvements
* Added support for treating variables as ordered categorical in the estimation of correlations.
* Renamed dummy data sets to `mantar_dummy_full_cont` and `mantar_dummy_mis_cont` to better reflect that they contain only continuous variables.
* Improved documentation for several functions and updated the README to reflect new functionality.

# mantar 0.1.0

* Initial release
