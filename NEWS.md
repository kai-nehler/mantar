# mantar 0.2.0

## Changes to defaults

* The default correlation method has changed: In version 0.1.0 all correlations were computed as **Pearson correlations**.  
  Starting with version 0.2.0, the default is to select the correlation type adaptively.  
  As a result, results obtained with version 0.1.0 may differ from those obtained with version 0.2.0.

## New features

* Allow different methods for correlation calculation based on the argument `ordered` across different functions, defining whether variables should be treated as ordered categorical or continuous.
* Add dummy data set `mantar_dummy_full_cat` with only ordered categorical variables.

# mantar 0.1.0

* Initial release
