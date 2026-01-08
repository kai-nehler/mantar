#' @title Network Estimation via Neighborhood Selection using Information Criteria
#'
#' @description
#' Estimates a network structure through node-wise regression models, where each
#' regression is selected via an information-criterion–based stepwise procedure.
#' The selected regression coefficients are subsequently combined into partial
#' correlations to form the final network.
#'
#' @param data Optional raw data matrix or data frame containing the variables
#' to be included in the network. May include missing values. If `data` is not
#' provided (`NULL`), a covariance or correlation matrix must be supplied in `mat`.
#' @param ns Optional numeric sample size specification. Can be a single value
#' (same sample size is used for all regressions) or a vector (e.g., variable-wise sample
#' sizes). When `data` is provided and `ns` is `NULL`, sample sizes are derived
#' automatically from `data`. When `mat` is supplied instead of raw data,
#' `ns` must be provided and should reflect the sample size underlying `mat`.
#' @param mat Optional covariance or correlation matrix for the variables to be
#' included in the network. Used only when `data` is `NULL`. If both `data` and
#' `mat` are supplied, `mat` is ignored. When `mat` is used, `ns` must also be
#' provided.
#' @param n_calc Character string specifying how per-variable sample sizes for
#' node-wise regression models are computed when `ns` is not supplied. If `ns`
#' is provided, its values are used directly and `n_calc` is ignored. Possible
#' values are:
#' \describe{
#'   \item{`"individual"`}{For each variable, uses the number of non-missing
#'   observations for that variable.}
#'   \item{`"average"`}{Computes the average number of non-missing observations
#'   across all variables and uses this average as the sample size for every
#'   variable.}
#'   \item{`"max"`}{Computes the maximum number of non-missing observations
#'   across all variables and uses this maximum as the sample size for every
#'   variable.}
#'   \item{`"total"`}{Uses the total number of rows in `data` as the sample size
#'   for every variable.}
#' }
#' @param ic_type Type of information criterion to compute for model selection in
#' the node-wise regression models. Options are `bic` (default), `aic`, `aicc`.
#' @param ordered Logical vector indicating whether each variable in `data`
#' should be treated as ordered categorical. Only used when `data` is provided.
#' If a single logical value is supplied, it is recycled to all variables.
#' @param pcor_merge_rule Character string specifying how regression weights
#' from the node-wise models are merged into partial correlations. Possible
#' values are:
#' \describe{
#'   \item{`"and"`}{Estimates a partial correlation only if the regression
#'   weights in both directions (e.g., from node 1 to 2 and from node 2 to 1)
#'   are non-zero in the final models.}
#'   \item{`"or"`}{Uses the available regression weight from one direction as
#'   the partial correlation if the corresponding regression in the other
#'   direction is not included in the final model.}
#' }
#' @param missing_handling Character string specifying how correlations are
#' estimated from the `data` input in the presence of missing values. Possible
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
#' @param ... Further arguments passed to internal functions.
#'
#' @details
#' This function estimates a network structure using neighborhood selection guided by information criteria.
#' Simulations by \insertCite{williams.2019;textual}{mantar} indicated that using the `"and"` rule for merging regression weights tends to yield more accurate partial correlation estimates than the `"or"` rule.
#'
#' The argument `ic_type` specifies which information criterion is computed.
#' All criteria are computed based on the log-likelihood of the maximum
#' likelihood estimated regression model, where the residual variance
#' determines the likelihood. The following options are available:
#'
#' \describe{
#'
#'   \item{\code{"aic"}:}{
#'     Akaike Information Criterion \insertCite{akaike.1974}{mantar}; defined as
#'     \mjseqn{\mathrm{AIC} = -2 \ell + 2k},
#'     where \eqn{\ell} is the log-likelihood of the model and \eqn{k} is the
#'     number of estimated parameters (including the intercept).
#'   }
#'
#'   \item{\code{"bic"}:}{
#'     Bayesian Information Criterion \insertCite{schwarz.1978}{mantar}; defined as
#'    \mjseqn{\mathrm{BIC} = -2 \ell + k \log(n)}, where \eqn{\ell} is
#'     the log-likelihood of the model, \eqn{k} is the
#'     number of estimated parameters (including the intercept)
#'     and \eqn{n} is the sample size.
#'   }
#'
#'   \item{\code{"aicc"}:}{
#'     Corrected Akaike Information Criterion \insertCite{hurvich.1989}{mantar};
#'     particularly useful in small samples where AIC tends to be biased.
#'     Defined as
#'      \mjseqn{\mathrm{AIC_c} = \mathrm{AIC} + \frac{2k(k+1)}{n - k - 1}},
#'     where \eqn{k} is the number of estimated parameters (including
#'     the intercept) and \eqn{n} is the sample size.
#'   }
#'
#' }
#'
#' \strong{Missing Handling}
#'
#' To handle missing data, the function offers two approaches: a two-step expectation-maximization (EM) algorithm and stacked multiple imputation.
#' According to simulations by \insertCite{nehler.2024;textual}{mantar}, stacked multiple imputation performs reliably across a range of sample sizes.
#' In contrast, the two-step EM algorithm provides accurate results primarily when the sample size is large relative to the amount of missingness and network complexity - but may still be preferred in such cases due to its much faster runtime.
#'
#' Currently, the function only supports variables that are directly included in the network analysis; auxiliary variables for missing handling are not yet supported.
#' During imputation, all variables are imputed by default using predictive mean matching \insertCite{@see e.g., @vanbuuren.2018}{mantar}, with all other variables in the data set serving as predictors.
#'
#' @import Rdpack
#'
#' @references
#' \insertAllCited{}
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{pcor}{Partial correlation matrix estimated from the node-wise regressions.}
#'   \item{betas}{Matrix of regression coefficients from the final regression models.}
#'   \item{ns}{Sample sizes used for each variable in the node-wise regressions.}
#'   \item{args}{List of settings used in the network estimation.}
#' }
#' @export
#'
#' @examples
#' # Estimate network from full data set
#' # Using Akaike information criterion
#' result <- neighborhood_net(data = mantar_dummy_full_cont,
#' ic_type = "aic")
#'
#' # View estimated partial correlations
#' result$pcor
#'
#' # Estimate network for data set with missings
#' # Using Bayesian Information Criterion, individual sample sizes, and two-step EM
#' result_mis <- neighborhood_net(data = mantar_dummy_mis_cont,
#' n_calc = "individual",
#' missing_handling = "two-step-em",
#' ic_type = "bic")
#'
#' # View estimated partial correlations
#' result_mis$pcor
neighborhood_net <- function(data = NULL, ns = NULL, mat = NULL, n_calc = "individual", ic_type = "bic",
                             ordered = FALSE, pcor_merge_rule = "and",
                             missing_handling = "two-step-em",
                              nimp = 20, imp_method = "pmm", ...){

  # Match arguments
  ic_type <- match.arg(tolower(ic_type), choices = c("bic", "aic", "aicc"))
  n_calc <- match.arg(tolower(n_calc), choices =c("average", "individual", "max", "total"))
  missing_handling <- match.arg(tolower(missing_handling), choices = c("two-step-em", "stacked-mi", "pairwise", "listwise"))
  pcor_merge_rule <- match.arg(tolower(pcor_merge_rule), choices = c("and", "or"))

  # Check: Which input is provided?
  checker(data = data, mat = mat)

  # If data is provided, compute correlation matrix
  if (!is.null(data)) {

    cor_out <- cor_calc(data = data, missing_handling = missing_handling,
                        ordered = ordered, nimp = nimp, imp_method = imp_method, ...)
    list2env(cor_out, envir = environment())
    list2env(cor_out$args, envir = environment())

    # if numbers of observations are not provided, calculate them
    if (is.null(ns)){
      if (!(n_calc %in% c("individual", "average", "max", "total"))){
        stop("Invalid n_calc value. Choose from 'individual', 'average', 'max', or 'total'.")
      }
      ns <- reg_calculate_sample_size(data = data, n_calc = n_calc)
    } else {
      # validate ns input if provided; if it is of length 1 recycle to all variables
      checker(ns = ns, data = data)
      if (length(ns) == 1){
        ns <- rep(ns, ncol(data))
      }
    }
  } else if (!is.null(mat)) {
    if (is.null(ns)) {
      stop("If 'mat' is provided, 'ns' must also be specified.")
    }
    # Check if combined input is valid
    checker(ns = ns, mat = mat)
    # If ns is of length 1, recycle to all variables
    if (length(ns) == 1){
      ns <- rep(ns, ncol(mat))
    }
    # Ensure that mat is a correlation matrix
    mat <- stats::cov2cor(mat)
    nimp <- missing_handling <- cor_method <- imp_method <- NULL
  }

  # Call helper function to perform neighborhood selection
  mod <- neighborhood_sel(mat = mat, ns = ns, ic_type = ic_type, pcor_merge_rule = pcor_merge_rule)

  # Prepare result
  result <- list(
    pcor = mod$partials,
    betas = mod$beta_mat,
    ns = ns,
    args = list(ic_type = ic_type, cor_method = cor_method, pcor_merge_rule = pcor_merge_rule,
                missing_handling = missing_handling, nimp = nimp, imp_method = imp_method)
  )

  # Set class and return
  class(result) <- c("mantar_neighborhood", "mantar_network")
  return(result)

}



#' @title Helper Function for Neighborhood Selection
#'
#' @description Internal helpfer function to perform neighborhood selection for
#' each variable in the correlation matrix and combine the results into partial
#' correlations.
#'
#' @param mat Correlation matrix
#' @param ns Sample sizes per variable
#' @param ic_type Type of information criterion for model selection
#' @param pcor_merge_rule Rule for merging regression weights into partial correlations
#'
#' @returns
#' A list with elements:
#' \describe{
#'   \item{partials}{Symmetric matrix of estimated partial correlations.}
#'   \item{beta_mat}{Matrix of regression coefficients from the nodewise
#'     regressions.}
#' }
#'
#' @noRd
neighborhood_sel <- function(mat, ns, ic_type, pcor_merge_rule){

  # be sure that mat is of class matrix and inputs match is valid
  class(mat) <- "matrix"
  checker(ns = ns, mat = mat)

  # number of variables
  p <- ncol(mat)
  # initialize matrix for regression coefficients
  beta_mat <- matrix(NA, p, p)
  colnames(beta_mat) <- rownames(beta_mat) <- colnames(mat)

  # process each variable as dependent variable
  for (dep in 1:p){
    # possible predictors and sample size for dependent variable
    possible_preds <- (1:p)[-dep]
    n <- ns[dep]
    # perform neighorhood selection for dependent node
    mod <- pred_search(mat = mat, dep_ind = dep, possible_pred_ind = possible_preds,
                  n = n, ic_type = ic_type)
    # store regression coefficients in beta matrix
    beta_mat[dep, mod$actual_preds] <- mod$actual_betas
  }

  # compute partial correlations from the resulting beta matrix
  partials <- compute_partials(betas = beta_mat, rule = pcor_merge_rule)

  # return results
  return(list(partials = partials$wadj,
              beta_mat = beta_mat))

}


#' @title Compute partial correlations from regression coefficients
#'
#' @description
#' Converts a matrix of regression coefficients (from neighborhood selection)
#' into a matrix of partial correlations using a specified symmetry rule
#'
#' @param betas Square matrix of regression coefficients
#' @param rule Rule for merging regression weights into partial correlations. Currently
#' supported:
#' \describe{
#'   \item{`"or"`}{Use the available coefficient if at least one directional
#'   regression provides an estimate. Missing coefficients on one side of the
#'   diagonal are imputed from the opposite side.}
#'   \item{`"and"`} (default) Only estimate partial correlations if both directional
#'   regressions provide non-zero coefficients.
#' }
#'
#' @returns A list with the following elements:
#' \describe{
#'   \item{wadj}{A matrix of partial correlations (weighted adjacency matrix).}
#'   \item{adj}{A binary adjacency matrix indicating edge presence (1) or
#'              absence (0).}
#' }
#' @noRd
#'
compute_partials <- function(betas, rule){

  if (rule == "or"){
    # replace missing entry with estimated entry on the opposing side of the diagonal
    betas[upper.tri(betas)][is.na(betas[upper.tri(betas)])]  <- betas[lower.tri(betas)][is.na(betas[upper.tri(betas)])]
    betas[lower.tri(betas)][is.na(betas[lower.tri(betas)])]  <- betas[upper.tri(betas)][is.na(betas[lower.tri(betas)])]
  }

  # values for the partial correlations
  Dummy <- betas * t(betas)

  # look for entries with differing signs of beta; bounded between 0 and 1
  Dummy[Dummy < 0] <- 0
  Dummy[Dummy > 1] <- 1

  # add the corresponding sign to the partial correlation
  wadj <- sign(betas) * sqrt(Dummy)
  # replace NA with 0
  wadj[is.na(wadj)] <- 0
  adj <- ifelse(wadj == 0, 0, 1)

  # set row and column names
  colnames(wadj) <- rownames(wadj) <- colnames(adj) <- rownames(adj) <- colnames(betas)

  # return results
  return(list(wadj = wadj,
              adj = adj))
}
