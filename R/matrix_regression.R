#' @title Multiple Regression based on a Covariance/Correlation Matrix
#'
#' @param mat Covariance or correlation matrix.
#' @param dep_ind Index of the dependent variable in `mat`.
#' @param pred_ind Vector of indices for the predictor variables in `mat`.
#'
#' @returns A list containing:
#' \describe{
#'  \item{beta}{Vector of regression coefficients for the predictors in
#'  `pred_ind`. If no predictors are specified, this is `NULL`.}
#'  \item{resid_var}{Residual variance of the regression model.}
#'  }
#' @noRd
matrix_regression <- function(mat, dep_ind, pred_ind) {

  # Check that exactly one dependent variable is specified
  if (length(dep_ind) != 1) {
    stop("dep_ind must have exactly one index as the function supports only a single dependent variable.")
  }

  # If there are no predictors, return variance of the dependent variable only
  if (length(pred_ind) == 0) {
    resid_var <- mat[dep_ind, dep_ind]
    return(list(beta = NULL, resid_var = resid_var))
  }

  # Compute regression coefficients and residual variance
  XX <- mat[pred_ind, pred_ind, drop = FALSE]
  Xy <- mat[pred_ind, dep_ind, drop = FALSE]
  beta <- solve(XX) %*% Xy
  resid_var <- mat[dep_ind, dep_ind] - drop(crossprod(beta, Xy))

  return(list(beta = beta,
              resid_var = resid_var))
}
