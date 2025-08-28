#' Calculate correlation matrix with different methods and missing data handling
#'
#' @param data A data frame or matrix containing the variables for which the correlation matrix is to be calculated.
#' @param ordered Specifies whether variables should be treated as ordered categorical when determining correlations.
#' Options are `TRUE`, `FALSE`, or `"adapted"`. The argument can be provided as a single value
#' (applied to all variables) or as a vector of length equal to the number of variables (using only `TRUE` and `FALSE`),
#' allowing mixed specifications. With the default `"adapted"`, the treatment of each variable is determined according to
#' guidelines from preliminary simulations (considering the number of cases, number of variables,
#' and number of categories).
#' @param missing_handling Method for estimating the correlation matrix in the presence of missing data.
#' `"tow-step-em"` uses a classic EM algorithm to estimate the correlation matrix from the data.
#' `"stacked-mi"` uses multiple imputation to estimate the correlation matrix from the data.
#' `"pairwise"` uses pairwise deletion to estimate the correlation matrix from the data.
#' `"listwise"` uses listwise deletion to estimate the correlation matrix from the data.
#' @param nimp Number of multiple imputations to perform when using multiple imputation for missing data (default: 20).
#' @param imp_method Method for multiple imputation when using `"stacked-mi"` for missing data handling. Default is `"pmm"` (predictive mean matching).
#' @param max_categories Maximum number of categories for a variable to be considered ordered if `cor_method = "adapted"` (default: 7).
#'
#' @details
#' While polychoric correlations are generally more appropriate for ordered categorical data, they rely on the assumption of multivariate normality
#' and may encounter estimation problems if the number of available observations is small relative to the number of estimated parameters. Preliminary
#' simulations suggest that in such cases Pearson correlations may introduce less bias, an effect that becomes even more pronounced when data are missing.
#'
#' We therefore allow users either to specify the correlation type directly or to use an adaptive option that selects the method based on the simulation
#' results. In general, variables with more than seven categories can be treated as continuous, whereas for variables with fewer categories the procedure
#' evaluates whether the amount of available information is too limited to justify polychoric estimation, in which
#' case Pearson correlations are used instead. The method is applied pairwise: polychoric correlations are computed if both variables are treated as
#' ordinal, polyserial correlations if one variable is ordinal and the other continuous, and Pearson correlations if both are continuous. The adaptive
#' procedure is still under development and may be refined in future versions.
#'
#' @returns A list with the following elements:
#' \describe{
#' \item{mat}{Estimated correlation matrix.}
#' \item{missing_handling}{Method used for handling missing data (if applicable).}
#' \item{nimp}{Number of imputations used (if applicable).}
#' \item{cor_method}{Matrix indicating the correlation method used for each pair of variables.}
#' }
#' @export
#'
#' @examples
#' # Estimate correlation matrix from full data set
#' result <- cor_calc(data = mantar_dummy_full, ordered = FALSE)
#'
#' # View estimated correlation matrix and methods used
#' result$mat
#' result$cor_method
#'
#' # Estimate correlation matrix for data set with missings
#' result_mis <- cor_calc(data = mantar_dummy_mis,
#'                       ordered = "adapted",
#'                       missing_handling = "two-step-em")
#'
#' # View estimated correlation matrix and methods used
#' result_mis$mat
#' result_mis$cor_method
cor_calc <- function(data, ordered = "adapted",
                     missing_handling = "two-step-em",
                     nimp = 20, imp_method = "pmm", max_categories = 7) {

  missing_handling <- match.arg(tolower(missing_handling), choices = c("two-step-em", "stacked-mi", "pairwise", "listwise"))

  checker(data = data, mat = NULL)

  if (length(ordered) == 1L) {
    if (is.logical(ordered)) {
      if (is.na(ordered)) {
        stop("`ordered` must not contain NA.")
      }
      ordered <- rep(ordered, ncol(data))
    } else if (is.character(ordered)) {
      ordered <- tolower(ordered)
      if (!identical(ordered, "adapted")) {
        stop('If `ordered` is a single character, it must be "adapted".')
      }
    } else {
      stop('If `ordered` has length 1, it must be logical (TRUE/FALSE) or "adapted".')
    }
  } else { # length > 1
    if (!is.logical(ordered)) {
      stop("If `ordered` has length > 1, it must be a logical vector (TRUE/FALSE).")
    }
    if (length(ordered) != ncol(data)) {
      stop("If `ordered` has length > 1, its length must equal the number of variables (p = ", ncol(data), ").")
    }
    if (any(is.na(ordered))) {
      stop("`ordered` must not contain NA.")
    }
  }

  if (unique(ordered) == "adapted") {
    if ((sum(!is.na(data))/80) < (ncol(data) * (ncol(data) - 1) / 2)) {
      ordered <- rep(FALSE, ncol(data))
      text <- "Using 'pearson' correlation method due to small ratio of available information to estimated correlations."
    } else {
      ordered <- vapply(data, function(x) {
        length(unique(x[!is.na(x)])) <= max_categories
      }, logical(1))
      if(any(ordered)) {
        if (sum(ordered) == ncol(data)) {
          text <- "Using only 'polychoric' correlations because all variables are ordered."
        } else {
          text <- paste("Using mix of 'pearson', 'polychoric' and 'polyserial' correlations because some variables are ordered, but not all. Ordered Variable are ", names(data)[ordered], collapse = " ")
        }
      } else {
        text <- "Using the 'pearson' correlation method because the number of distinct values suggests continuous variables."
      }
    }
    message(text)
  }

  # Check if package is installed if polychoric method is selected
  if (any(ordered) && !requireNamespace("lavaan", quietly = TRUE)) {
    stop(
      "Package \"lavaan\" must be installed to compute polychoric or polyserial correlations.",
      call. = FALSE
    )
  }

  if (anyNA(data)){
    if (missing_handling == "two-step-em"){

      if (!requireNamespace("lavaan", quietly = TRUE)) {
        stop(
          "Package \"lavaan\" must be installed to use this function.",
          call. = FALSE
        )
      }

      lavobject <- suppressWarnings(try(lavaan::lavCor(data,
                                                       missing = "ml", se = "none", meanstructure = TRUE,
                                                       estimator = "ML", output = "fit")))
      if (inherits(lavobject, "try-error")) stop("lavaan::lavCor failed. Check your data.")
      mat <- stats::cov2cor(lavaan::inspect(lavobject, "cov.ov"))
      if (any(ordered)){
        warning("Using a specific two-step EM algorithm based on polychoric or polyserial correlations is not implemented at the moment.
                  Will proceed with the two-step EM approach based on pearson correlations.")
      }
      nimp <- NULL

    } else if (missing_handling == "stacked-mi"){

      if (!requireNamespace("mice", quietly = TRUE)) {
        stop(
          "Package \"mice\" must be installed to use this function.",
          call. = FALSE
        )
      }
      original_names <- colnames(data)
      colnames(data) <- make.names(original_names)

      imputed_data <- suppressMessages(mice::mice(data, m = nimp, maxit = 10, method = imp_method,
                                                  ridge = 0, donors = 5, ls.method = "qr"))
      stacked_data <- mice::complete(imputed_data, c(1:nimp))
      colnames(stacked_data) <- original_names  # restore original column names
      if (any(ordered)){
        mat <- suppressWarnings(try(lavaan::lavCor(data, ordered = names(data)[ordered],
                                                   se = "none", output = "cor")))
        if (inherits(mat, "try-error")) stop("lavaan::lavCor failed. Check your data.")
      }else {
        mat <- stats::cor(stacked_data)
      }

    } else if (missing_handling == "pairwise"){
      if (any(ordered)){
        warning("Using pairwise deletion based on polychoric or polyserial correlations is not implemented at the moment.
                Will proceed with pairwise deletion based on pearson correlations.")
      }
      mat <- stats::cov2cor(stats::cov(data, use = "pairwise.complete.obs"))
      nimp <- NULL
    } else if (missing_handling == "listwise"){
      if (any(ordered)){
        mat <- suppressWarnings(try(lavaan::lavCor(data, ordered = names(data)[ordered],
          missing = "listwise", se = "none", output = "cor")))
        if (inherits(mat, "try-error")) stop("lavaan::lavCor failed. Check your data.")
      } else {
        mat <- stats::cov2cor(stats::cov(data, use = "complete.obs"))
      }
      nimp <- NULL
    }
  } else {
    if (any(ordered)){
      mat <- suppressWarnings(try(lavaan::lavCor(data, ordered = names(data)[ordered],
                                                 se = "none", output = "cor")))
      if (inherits(mat, "try-error")) stop("lavaan::lavCor failed. Check your data.")
    } else {
      mat <- stats::cov2cor(stats::cov(data))
    }
    missing_handling <- NULL
    nimp <- NULL
  }

  cor_method <- matrix(NA_character_, ncol(data), ncol(data), dimnames = list(names(data), names(data)))

  cor_method[lower.tri(cor_method)] <- ifelse(
    outer(ordered, ordered, "&")[lower.tri(cor_method)], "polychoric",
    ifelse(
      outer(ordered, ordered, "!=")[lower.tri(cor_method)], "polyserial",
      "pearson"
    )
  )


  return(list(mat = mat, missing_handling = missing_handling, nimp = nimp, cor_method = cor_method))
}
