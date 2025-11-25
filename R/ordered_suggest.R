#' @title
#' Heuristic procedure for identifying ordered categorical variables
#'
#' @description
#' Suggests which variables in a data set may be treated as ordered categorical
#' based on their number of unique categories and the amount of available
#' information for estimating the network structure. This function provides a
#' preliminary, non-binding recommendation and should be interpreted as a
#' beta-level heuristic.
#'
#' @param data Raw data matrix or data frame containing the variables to be
#' included in the network. May include missing values.
#' @param max_categories Maximum number of categories a variable may have to be
#' treated as ordered (default: 7).
#'
#' @details
#' While polychoric correlations are generally more appropriate for ordered categorical
#' data \insertCite{@foldness.2022}{mantar}, they may encounter estimation problems if
#' the number of available observations is small relative to the number of estimated
#' parameters \insertCite{@see e.g., @johal.2023}{mantar}. Our preliminary simulations
#' suggest that in such cases Pearson correlations may introduce less bias, an effect
#' that becomes even more pronounced when data are missing.
#'
#' This helper function provides a recommendation on which variables to treat as
#' ordered. In general, variables with more than \code{max_categories} categories are
#' recommended to be treated as continuous, whereas for variables with fewer categories
#' the procedure evaluates whether the amount of available information is too limited
#' to justify polychoric estimation, in which case Pearson correlations are recommended
#' instead. This procedure is only a helper, is still under early development, and may
#' be refined in future versions.
#'
#' @references
#' \insertAllCited{}
#'
#' @returns
#' A logical vector of length \code{ncol(data)} indicating, for each variable, whether
#' it is recommended to be treated as ordered (\code{TRUE}) or continuous
#' (\code{FALSE}). Additionally, a message is printed to the console summarizing the
#' recommendation in terms of which correlation methods to use.
#'
#' @export
#'
#' @examples
#' # Suggest ordered variables in a data set with mixed variable types
#' # (400 observations for 8 variables)
#' ordered_suggest(data = mantar_dummy_full_mix, max_categories = 7)
#'
ordered_suggest <- function(data, max_categories = 7) {

  # Number of variables
  p <- ncol(data)

  # Check ratio of available information to estimated coefficients
  # if this is too small, suggest pearson for all variables
  if ((sum(!is.na(data)) / 80) < (p * (p - 1) / 2)) {
    ordered <- rep(FALSE, p)
    names(ordered) <- colnames(data)
    text <- "Using the 'pearson' correlation method is advised due to the small ratio of available information to estimated correlations."
  } else {

    # if the ratio of available information to estimated coefficients is sufficient,
    # suggest ordered for variables with <= max_categories distinct values
    ordered <- vapply(
      data,
      function(x) {
        length(unique(x[!is.na(x)])) <= max_categories
      },
      logical(1)
    )
    names(ordered) <- colnames(data)

    # Add messages according to the suggested treatment as ordered
    if (any(ordered)) {
      if (all(ordered)) {
        text <- "Using only 'polychoric' correlations is advised because all variables are ordered and there is sufficient information."
      } else {
        text <- paste(
          "Using a mix of 'pearson', 'polychoric', and 'polyserial' correlations is advised because some variables are ordered but not all, and there is sufficient information. Ordered variables are",
          paste(names(ordered)[ordered], collapse = ", ")
        )
      }
    } else {
      text <- "Using the 'pearson' correlation method is advised because the number of distinct values suggests continuous variables."
    }
  }

  # print message and return suggestion vector
  message(text)
  return(ordered)

}
