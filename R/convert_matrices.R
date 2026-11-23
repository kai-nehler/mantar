#' Convert a correlation matrix to a covariance matrix using standard deviations
#' @noRd
cor_to_cov <- function(cor_mat, sds) {
  cor_mat * outer(sds, sds)
}

#' Convert temporal VAR coefficients to partial directed correlations
#' @noRd
calculate_pdc <- function(beta, sigma, kappa) {
  t(beta / sqrt(diag(sigma) %o% diag(kappa) + beta^2))
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
