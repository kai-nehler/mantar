#' @title
#' Correlation Matrix Estimation with Support for Multiple Correlation Types
#'
#' @description
#' Computes a correlation matrix from raw data while accounting for missing
#' values through several missing-data handling strategies. Supports different
#' correlation types based on whether variables are treated as ordered.
#'
#' @param data Data frame or matrix containing the variables for which the
#' correlation matrix is to be computed. May include missing values.
#' @param ordered Logical vector indicating whether each variable in `data`
#' should be treated as ordered categorical when computing the correlation
#' matrix. If a single logical value is supplied, it is recycled to all
#' variables.
#' @param missing_handling Character string specifying how the correlation
#' matrix is estimated from `data` in the presence of missing values. Possible
#' values are:
#' \describe{
#'   \item{`"two-step-em"`}{Uses a classical EM algorithm to estimate the
#'   correlation matrix from `data`.}
#'   \item{`"stacked-mi"`}{Uses stacked multiple imputation to estimate the
#'   correlation matrix from `data`.}
#'   \item{`"pairwise"`}{Uses pairwise deletion to compute correlations from
#'   `data`.}
#'   \item{`"listwise"`}{Uses listwise deletion to compute correlations from
#'   `data`.}
#' }
#' @param nimp Number of imputations (default: 20) to be used when
#' `missing_handling = "stacked-mi"`.
#' @param imp_method Character string specifying the imputation method to be
#' used when `missing_handling = "stacked-mi"` (default: `"pmm"` - predictive
#' mean matching).
#' @param maxit Maximum number of iterations for the imputation algorithm when
#' `missing_handling = "stacked-mi"` (default: 10).
#' @param ... Further arguments passed to internal functions.
#'
#' @details
#' Correlations are computed pairwise:
#' * Polychoric correlations for two ordered variables,
#' * Polyserial correlations for one ordered and one continuous variable,
#' * Pearson correlations for two continuous variables.
#'
#' Treating variables as ordered requires the missing handling method to be either
#' `"stacked-mi"` or `"listwise"`
#'
#' Means are computed whenever Pearson correlations are used. If any variable
#' is treated as ordered, `means` is returned as NULL.
#'
#' @returns A list containing:
#' \describe{
#'   \item{mat}{Estimated correlation matrix.}
#'   \item{means}{Vector of estimated means. If any variable is treated as ordered,
#'                means is returned as NULL.}
#'   \item{cor_method}{A matrix indicating the correlation method used for each
#'                     variable pair.}
#'   \item{args}{List of settings used in the correlation estimation.}
#' }
#' @export
#'
#' @examples
#' # Estimate correlation matrix from full data set
#' result <- cor_calc(data = mantar_dummy_full_cont,
#'                    ordered = FALSE)
#'
#' # View estimated correlation matrix and methods used
#' result$mat
#' result$cor_method
#'
#' # Estimate correlation matrix for data set with missings
#' result_mis <- cor_calc(data = mantar_dummy_mis_cont,
#'                       ordered = FALSE,
#'                       missing_handling = "two-step-em")
#'
#' # View estimated correlation matrix and methods used
#' result_mis$mat
#' result_mis$cor_method
cor_calc <- function(data, ordered = FALSE,
                     missing_handling = "two-step-em",
                     nimp = 20, imp_method = "pmm",
                     maxit = 10, ...) {

  # Capture additional arguments
  dots <- list(...)

  # Match and validate missing handling method
  missing_handling <- match.arg(tolower(missing_handling),
                                choices = c("two-step-em", "stacked-mi", "pairwise", "listwise"))
  # Check data input - setting mat = NULL important to perform the correct check
  checker(data = data, mat = NULL)

  # Check ordered input to be logical, without NAs and have correct length
  if (any(is.na(ordered))) {
    stop("`ordered` must not contain NA.")
  }
  if (any(!is.logical(ordered))){
    stop("All entries of `ordered` must be logical (TRUE/FALSE).")
  }
  if (length(ordered) == 1L) {
      # Recycle to all variables
      ordered <- rep(ordered, ncol(data))
    } else {
      if (length(ordered) != ncol(data)) {
        # Incorrect length - abort the estimation
        stop("If `ordered` has length > 1, its length must equal the number of variables (p = ", ncol(data), ").")
      }
    }

  # Check if package is installed if polychoric method is selected
  any_ord <- any(ordered)
  if (any_ord && !requireNamespace("lavaan", quietly = TRUE)) {
    stop(
      "Package \"lavaan\" must be installed to compute polychoric or polyserial correlations.",
      call. = FALSE
    )
  }

  # Prepare means vector
  means <- NULL

  # Handle missing data
  if (anyNA(data)){
    if (missing_handling == "two-step-em"){

      # EM algorithm cannot treat variables as ordered
      if (any_ord) {
        stop(
          "Cannot use 'two-step-em' when any variables are ordered Either set ordered = FALSE or choose a
           different missing data method."
        )
      }

      # Check if lavaan is installed for the EM algorithm
      if (!requireNamespace("lavaan", quietly = TRUE)) {
        stop(
          "Package \"lavaan\" must be installed to use the EM algorithm to handle missingness.",
          call. = FALSE
        )
      }

      # Estimate correlation matrix with EM algorithm
      lavobject <- suppressWarnings(try(lavaan::lavCor(data,
                                                       missing = "ml", se = "none", meanstructure = TRUE,
                                                       estimator = "ML", output = "fit")))
      # If lavcor fails - stop the estimation
      if (inherits(lavobject, "try-error")) stop("lavaan::lavCor failed. Check your data.")
      # Extract estimated correlation matrix and means
      mat <- stats::cov2cor(lavaan::inspect(lavobject, "cov.ov"))
      means <- lavaan::inspect(lavobject, "mean.ov")

      # Set nimp to NULL as no multiple imputation was used
      nimp <- NULL

    } else if (missing_handling == "stacked-mi"){

      # Check if mice is installed for multiple imputation
      if (!requireNamespace("mice", quietly = TRUE)) {
        stop(
          "Package \"mice\" must be installed to use this function.",
          call. = FALSE
        )
      }

      # Keep column names safe during imputation
      original_names <- colnames(data)
      colnames(data) <- make.names(original_names)

      # Perform multiple imputation and stack the data
      imputed_data <- suppressMessages(mice::mice(data, m = nimp, maxit = 10,
                                                  method = imp_method, ...))
      stacked_data <- mice::complete(imputed_data, c(1:nimp))
      colnames(stacked_data) <- original_names  # restore original column names
      if (any_ord) {
        mat <- suppressWarnings(
          try(lavaan::lavCor(stacked_data, ordered = names(data)[ordered],
                             se = "none", output = "cor"))
        )
        if (inherits(mat, "try-error"))
          stop("lavaan::lavCor failed. Check your data.")
      } else {
        mat <- stats::cor(stacked_data)
        means <- colMeans(stacked_data)
      }

    } else if (missing_handling == "pairwise"){
      if (any_ord) {
        stop(
          "Cannot use 'pairwise' when any variables are treated as ordered. Either set ordered = FALSE or choose a
           different missing data method."
        )
      }
      mat <- stats::cov2cor(stats::cov(data, use = "pairwise.complete.obs"))
      means <- colMeans(data, na.rm = TRUE)
      nimp <- imp_method <- maxit <- NULL

    } else if (missing_handling == "listwise") {

      data_complete <- stats::na.omit(data)

      if (any_ord) {
        mat <- suppressWarnings(
          try(lavaan::lavCor(data_complete, ordered = names(data)[ordered],
                             missing = "listwise", se = "none", output = "cor"))
        )
        if (inherits(mat, "try-error"))
          stop("lavaan::lavCor failed. Check your data.")
      } else {
        mat <- stats::cov2cor(stats::cov(data_complete, use = "complete.obs"))
        means <- colMeans(data_complete)
      }
      nimp <- imp_method <- maxit <- NULL
    }
  } else {
    if (any(ordered)){
      mat <- suppressWarnings(try(lavaan::lavCor(data, ordered = names(data)[ordered],
                                                 se = "none", output = "cor")))
      if (inherits(mat, "try-error")) stop("lavaan::lavCor failed. Check your data.")
    } else {
      mat <- stats::cov2cor(stats::cov(data))
      means <- colMeans(data)
    }
    missing_handling <- NULL
    nimp <- imp_method <- maxit <- NULL
  }

  cor_method <- matrix("", ncol(data), ncol(data), dimnames = list(names(data), names(data)))
  diag(cor_method) <- "-"

  cor_method[lower.tri(cor_method)] <- ifelse(
    outer(ordered, ordered, "&")[lower.tri(cor_method)], "polychoric",
    ifelse(
      outer(ordered, ordered, "!=")[lower.tri(cor_method)], "polyserial",
      "pearson"
    )
  )


  return(list(
    mat = mat,
    means = means,
    cor_method = cor_method,
    args = list(
      missing_handling = missing_handling,
      nimp = nimp)
  ))
}


transform_inverse <- function(inverse_cov){
  rows <- nrow(inverse_cov)
  cols <- ncol(inverse_cov)
  cov <- solve(inverse_cov)
  inverse_cor <- matrix(NA, rows, cols)

  # it is important to do it this way, cause using solve, cov2cor and solve again introduces small errors
  for (i in 1:rows){
    for (j in 1:cols){
      inverse_cor[i,j] <- inverse_cov[i,j] * sqrt(diag(cov)[i]) * sqrt(diag(cov)[j])
    }
  }
  return(inverse_cor)
}

#' @title Convert Inverse Covariance Matrix to Network Matrix
#'
#' @param theta Inverse covariance matrix.
#'
#' @returns Network matrix with partial correlations.
#' @noRd
inv_to_net <- function(theta){
  net <- -stats::cov2cor(theta)
  diag(net) <- 0
  net <- Matrix::forceSymmetric(net) |> as.matrix()
  return(net)
}


