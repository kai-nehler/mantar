#' @title Calculate sample size for the regression models
#'
#' @param data Data frame or matrix with data
#' @param n_calc Character string indicating how sample size should be calculated:
#'   \itemize{
#'     \item `"individual"` – number of non-missing values per variable.
#'     \item `"average"` – average number of non-missing observations across variables.
#'     \item `"max"` – maximum number of non-missing observations across variables.
#'     \item `"total"` – total number of rows in the dataset.
#'   }
#'
#' @returns A numeric vector with sample sizes for each variable
#' @noRd
reg_calculate_sample_size <- function(data, n_calc = c("individual", "average", "max", "total")){

  # keep user input consistent
  n_calc <- tolower(n_calc[1])
  n_calc <- match.arg(n_calc, c("individual", "average", "max", "total"))

  # if no missing values, return total n for all variables
  if (!anyNA(data)) {
    ns <- rep(nrow(data), ncol(data))
    message("No missing values in data. Sample size for each variable is equal to the number of rows in the data.")
  } else {
    if (n_calc == "individual"){
      # each variable gets its own sample size
      ns <- colSums(!is.na(data))
      } else if (n_calc == "average"){
        # calculate mean number of non-missing observations across variables
        ns <- rep(mean(colSums(!is.na(data))), ncol(data))
        } else if (n_calc == "max"){
          # number of observation of "best" case variable
          ns <- rep(max(colSums(!is.na(data))), ncol(data))
          } else if (n_calc == "total"){
            # total number of rows in data - disregarding missings completely
            ns <- rep(nrow(data), ncol(data))
            }
  }
  return(ns)
}




#' @title Calculate information criterion for regression model based on residual variance
#'
#' @param resid_var Residual variance of the regression model
#' @param n Integer; number of observations used to fit the model
#' @param n_preds Integer; number of predictors in the model
#' @param ic_type Character; type of information criterion to compute.
#' Options: "AIC", "BIC", "AICc"
#'
#' @return A single numeric value with the information criterion
#' @noRd
reg_ic_calc <- function(resid_var, n, n_preds, ic_type = "bic"){

  ic_type <- match.arg(tolower(ic_type), choices =c("aic", "bic", "aicc"))

  # Determine penalty k based on the type
  # add 1 to n_preds for the intercept in aicc calculation
  k <- switch(ic_type,
              aic  = 2,
              bic  = log(n),
              aicc = 2 + (2 * (n_preds + 1) * (n_preds + 2)) / (n - n_preds - 2)
  )

  # calculate log-likelihood
  LL <- -(n/2) * (log(resid_var) + log(2*pi) + 1)

  # calculate information criterion; add 1 to the for the intercept
  IC <- -2 * LL + k * (n_preds + 1)

  return(IC)
}



#' @title Calculate effective sample size for the graphical models
#'
#' @param data Data frame or matrix with data
#' @param n_calc Character string indicating how sample size should be calculated:
#'       \itemize{
#'       \item `"average"` – average number of non-missing observations across variable pairs.
#'       \item `"max"` – maximum number of non-missing observations across variable pairs.
#'       \item `"total"` – total number of rows in the dataset.
#'       }
#' @param count_diagonal Logical; whether to include diagonal elements in the average sample size calculation
#'
#' @returns A numeric value with the effective sample size
#' @noRd
mat_calculate_sample_size <- function(data, n_calc = c("average", "max", "total"),
                                      count_diagonal = TRUE) {

  # keep user input consistent
  n_calc <- tolower(n_calc[1])
  n_calc <- match.arg(n_calc)

  # create non-missingness indicator matrix
  R <- !is.na(data)
  # calculate number of non-missing observations for each variable pair
  nonmisMatrix <- t(R) %*% R

  if (n_calc == "average") {

    if (count_diagonal) {
      # include diagonal elements in the average calculation
      n <- mean(nonmisMatrix)
    } else {
      # exclude diagonal elements in the average calculation
      n <- mean(nonmisMatrix[upper.tri(nonmisMatrix)])
    }

  } else if (n_calc == "max") {
    # maximum number of non-missing observations across variable pairs (without diagonal)
    n <- max(nonmisMatrix[upper.tri(nonmisMatrix)])
  } else if (n_calc == "total") {
    # total number of rows in data - disregarding missings completely
    n <- nrow(data)
  }

  return(n)
}

#' @title Computes information criterion for graphical models
#'
#' @param data Optional data matrix or data frame. Required if
#' `likelihood = "obs_based"`.
#' @param sample_cor Optional sample correlation matrix. Required if
#' `data` is not provided.
#' @param theta Precision matrix of the graphical model
#' @param mu Optional vector of variable means. Required if
#' `likelihood = "obs_based"`.
#' @param n Integer; effective sample size for the model.
#' @param ic_type Type of information criterion to compute. Options are: `"bic"`,
#' `"ebic"` and `"aic"`.
#' @param extended_gamma Numeric; tuning parameter for the extended information
#' criterion
#' @param likelihood Character string; which likelihood to use for the
#' IC calculation. Possible values are:
#'      \itemize{
#'      \item `"obs_based"` – uses the observed-data based likelihood.
#'      \item `"mat_based"` – uses the matrix-based likelihood, which only
#'      requires a sample correlation matrix and an effective sample size.
#'      }
#' @param ... Further arguments passed to internal functions.
#'
#' @returns A single numeric value with the information criterion
#' @noRd
mat_ic_calc <- function(data = NULL, sample_cor = NULL, theta, mu = NULL, n,
                        ic_type = "ebic", extended_gamma = 0.5, likelihood = "obs_based", ...){

  # ensure ic_type is consistent
  ic_type <- match.arg(tolower(ic_type), choices = c("ebic", "bic", "aic"))

  # number of variables
  p <- ncol(theta)

  if (likelihood == "obs_based"){

    # define matrix with missing checks
    R <- !is.na(data)

    # find all unique rows in R
    unique_patterns <- unique(R)
    # find number of unique missing data patterns
    n_patterns <- nrow(unique_patterns)

    # match all persons to a pattern using the missing check in R
    matched_indices <- apply(unique_patterns, 1, function(pat) {
      which(colSums(pat == t(R)) == ncol(R))
    }, simplify = FALSE
    )

    ## make sure that the inverse is an inverse correlation matrix
    theta <- transform_inverse(theta)
    # Compute covariance matrix
    sigma <- solve(theta)


    if (likelihood == "obs_based"){
      # theta is an inverse correlation matrix (so the entries of sigma on the diagonal are 1)
      # need to transform the data in a way, where the SAMPLE covariance is 1
      cova <- sigma * (nrow(data)-1)/nrow(data)
      vars <- t(rbind(replicate(nrow(data), sqrt(diag(cova)), simplify = TRUE)))
      data <- data / vars
    }


    # starting likelihood calculation
    # count number of missing patterns
    ll_patterns <- numeric(n_patterns)

    # calculate loglikelihood for each pattern
    for (i in 1:n_patterns){
      ll_patterns[i] <- pattern_ll(data = data, rows = matched_indices[[i]],
                                   cols = unique_patterns[i,], sigma = sigma, theta = theta, mu = mu)
    }

    # total loglikelihood is sum over all patterns
    loglikelihood <- sum(ll_patterns)

  }  else if (likelihood == "mat_based"){
    # matrix-based loglikelihood calculation
    const1 <- -0.5 * n * p * log(2 * pi)
    loglikelihood <- const1 + n/2 * (log(det(theta)) - sum(diag(theta %*% sample_cor)))
  }

  # find number of non-zero edges
  edges <- sum(theta[lower.tri(theta, diag = FALSE)] != 0)

  if (ic_type == "ebic"){
    IC <- -2 * loglikelihood + edges * log(n) + 4 * edges * extended_gamma * log(p)
  } else if(ic_type == "bic") {
    IC <- -2 * loglikelihood + edges * log(n)
  } else if(ic_type == "aic") {
    IC <- -2 * loglikelihood + edges * 2
  }

  return(IC)
}


#' @title Log-likelihood contribution for a single specific missing data pattern
#'
#' @param data Data frame or matrix with observed data
#' @param rows Indices of rows belonging to this missing data pattern
#' @param cols Logical vector indicating observed (TRUE) and missing (FALSE)
#' variables for the specific pattern
#' @param sigma Covariance matrix of the graphical model
#' @param theta Inverse covariance matrix of the graphical model
#' @param mu Vector of variable means
#'
#' @returns A single numeric value giving the log-likelihood contribution of
#' all observations belongig to the specific missing data pattern
#' @noRd
pattern_ll <- function (data, rows, cols, sigma, theta, mu)
{

  # reduce X to observations belonging to this pattern
  observations <- data[rows, cols, drop = F]

  # check thatt observations is not empty
  if (ncol(observations) == 0){
    warning("There are observations with no observed variables. Their log-likelihood contribution is \n set to zero. Consider rerunning the analysis without these observations as they still count \n towards the number of observations.")
    return(0)
  }

  # define submatrix of sigma with only observed variables
  sigma_oo <- sigma[cols, cols, drop = FALSE]

  # use lavaan to get submatrix of the inverse without needing to invert
  na_idx <- which(!cols)
  theta_oo <- mat_inverse_update(theta = theta,
                                 rm_idx = na_idx,
                                 logdet = FALSE)

  # pre term involving number of observed variables
  pre_log <- ncol(observations) * log(2 * pi)

  # log determinant of the correlation matrix
  logdet <- log(det(sigma_oo))

  # mahalanobis distance seperately for every row
  # first, center the observations by subtracting the mean
  observations_cen <- as.matrix(sweep(observations, 2, mu[cols], FUN = "-"))
   # now calculate the mahalanobis distance
  mahalanobis_dist <- rowSums((observations_cen %*% theta_oo) * observations_cen)

  # loglikelihood contributions for all persons of this pattern
  loglik <-  sum(-(pre_log + logdet + mahalanobis_dist)/2)
  loglik
}



#' @title Update inverse of a matrix after removing rows/columns
#'
#' @param theta Symmetric, positive definite inverse covariance matrix
#' @param rm_idx Indices of rows/columns to be removed
#' @param logdet Logical; whether to update the log-determinant attribute
#' @param mat_logdet Log-determinant of the original covariance matrix (if `logdet = TRUE`)
#'
#' @details
#' This function imitates the behavior of the `lavaan` version 0.6.20
#' \insertCite{lavaan.2025}{mantar}internal function
#' \code{lav_matrix_symmetric_inverse_update} as this functionality was
#' not exported in that version.
#'
#' @returns Updated inverse covariance matrix with removed rows/columns.
#' If `logdet = TRUE`, the updated log-determinant is added as an attribute.
#' @noRd
mat_inverse_update <- function(theta, rm_idx, logdet, mat_logdet){

  # number of variables
  p <- ncol(theta)

  # ensure rm_idx contains unique indices
  if (any(duplicated(rm_idx))) {
    stop("`rm_idx` must contain unique indices.")
  }
  # ensure rm_idx is inside matrix bounds
  if (any(rm_idx < 1) || any(rm_idx > p)) {
    stop("All entries of `rm_idx` must be between 1 and ", p, ".")
  }
  # count number of removed indices
  n_del <- length(rm_idx)

  # if logdet is requested, require the log-determinant input
  if (logdet && is.null(mat_logdet)) {
    stop("`logdet = TRUE`, but `mat_logdet` is NULL. ",
         "Provide the log-determinant of the original covariance matrix.")
  }

  if (n_del == 0) {
    # if nothing is removed, return original matrix
    out <- theta
    if (logdet) {
      attr(out, "logdet") <- mat_logdet
      }
    } else if (n_del == 1){
      # if one index is removed, avoid costly matrix inversion
      h <- theta[rm_idx, rm_idx]
      a <- theta[-rm_idx, rm_idx, drop = FALSE]/sqrt(h)
      out <- theta[-rm_idx, -rm_idx, drop = FALSE] - tcrossprod(a)
      if (logdet) {
        attr(out, "logdet") <- mat_logdet + log(h)
        }
      } else if (n_del < p) {
        # if multiple indices are removed, perform full update
        H <- theta[rm_idx, rm_idx, drop = FALSE]
        A <- theta[rm_idx, -rm_idx, drop = FALSE]
        out <- theta[-rm_idx, -rm_idx, drop = FALSE] - crossprod(A, solve.default(H,A))
        if (logdet) {
          attr(out, "logdet") <- mat_logdet + log(det(H))
          }
        } else if (n_del == p) {
          # if all indices are removed, return empty matrix
          out <- matrix(NA, nrow = 0, ncol = 0)
          if (logdet) {
            attr(out, "logdet") <- 0
          }
          }

  return(out)

}
