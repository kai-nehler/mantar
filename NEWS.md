# mantar 0.2.0

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
