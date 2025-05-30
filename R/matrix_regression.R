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

  return(list(beta = beta, resid_var = resid_var))
}
