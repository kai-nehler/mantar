#' Estimate Network using Neighborhood Selection based on Information Criteria
#'
#' @param data Raw data containing only the variables to be included in the network. May include missing values.
#' @param ns Numeric vector specifying the sample size for each variable in the data.
#' If not provided, it will be computed based on the data.
#' Must be provided if a correlation matrix (`mat`) is supplied instead of raw data.
#' @param mat Optional covariance or correlation matrix for the variables to be included in the network.
#' Used only if \code{data} is \code{NULL}.
#' @param n_calc Method for calculating the sample size for node-wise regression models. Can be one of:
#' `"individual"` (sample size for each variable is the number of non-missing observations for that variable),
#' `"average"` (sample size is the average number of non-missing observations across all variables),
#' `"max"` (sample size is the maximum number of non-missing observations across all variables),
#' `"total"` (sample size is the total number of observations across in the data set / number of rows).
#' @param k Penalty per parameter (number of predictor + 1) to be used in node-wise regressions; the default `k = "log(n)"` (number of observations for the dependent variable) is the classical BIC. Alternatively, classical AIC would be `k = "2"`.
#' @param ordered Specifies whether variables should be treated as ordered categorical when determining correlations.
#' Options are `TRUE`, `FALSE`, or `"adapted"`. The argument can be provided as a single value
#' (applied to all variables) or as a vector of length equal to the number of variables (using only `TRUE` and `FALSE`),
#' allowing mixed specifications. With the default `"adapted"`, the treatment of each variable is determined according to
#' guidelines from preliminary simulations (considering the number of cases, number of variables,
#' and number of categories). See the *Details* section of [cor_calc()] for further elaboration.
#' @param pcor_merge_rule Rule for merging regression weights into partial correlations.
#' `"and"` estimates a partial correlation only if regression weights in both directions (e.g., from node 1 to 2 and from 2 to 1) are non-zero in the final models.
#' `"or"` uses the available regression weight from one direction as partial correlation if the other is not included in the final model.
#' @param missing_handling Method for estimating the correlation matrix in the presence of missing data.
#' `"tow-step-em"` uses a classic EM algorithm to estimate the correlation matrix from the data.
#' `"stacked-mi"` uses multiple imputation to estimate the correlation matrix from the data.
#' `"pairwise"` uses pairwise deletion to estimate the correlation matrix from the data.
#' `"listwise"` uses listwise deletion to estimate the correlation matrix from the data.
#' @param nimp Number of multiple imputations to perform when using multiple imputation for missing data (default: 20).
#' @param imp_method Method for multiple imputation when using `"stacked-mi"` for missing data handling. Default is `"pmm"` (predictive mean matching).
#'
#' @details
#' This function estimates a network structure using neighborhood selection guided by information criteria.
#' Simulations by Williams et al. (2019) indicated that using the `"and"` rule for merging regression weights tends to yield more accurate partial correlation estimates than the `"or"` rule.
#' Both the Akaike Information Criterion (AIC) and the Bayesian Information Criterion (BIC) are supported and have been shown to produce valid network structures.
#'
#' To handle missing data, the function offers two approaches: a two-step expectation-maximization (EM) algorithm and stacked multiple imputation.
#' According to simulations by Nehler and Schultze (2024), stacked multiple imputation performs reliably across a range of sample sizes.
#' In contrast, the two-step EM algorithm provides accurate results primarily when the sample size is large relative to the amount of missingness and network complexity—but may still be preferred in such cases due to its much faster runtime.
#'
#' Currently, the function only supports variables that are directly included in the network analysis; auxiliary variables for missing handling are not yet supported.
#' During imputation, all variables are imputed by default using predictive mean matching (see, e.g., van Buuren, 2018), with all other variables in the dataset serving as predictors.
#'
#' @references
#' Nehler, K. J., & Schultze, M. (2024). *Handling missing values when using neighborhood selection for network analysis*. https://doi.org/10.31234/osf.io/qpj35
#'
#' van Buuren, S. (2018). *Flexible Imputation of Missing Data* (2nd ed.). CRC Press.
#'
#' Williams, D. R., Rhemtulla, M., Wysocki, A. C., & Rast, P. (2019). On nonregularized estimation of psychological networks. *Multivariate Behavioral Research, 54*(5), 719–750. https://doi.org/10.1080/00273171.2019.1575716

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
#' result <- neighborhood_net(data = mantar_dummy_full,
#' k = "2")
#'
#' # View estimated partial correlations
#' result$pcor
#'
#' # Estimate network for data set with missings
#' # Using Bayesian Information Criterion, individual sample sizes, and two-step EM
#' result_mis <- neighborhood_net(data = mantar_dummy_mis,
#' n_calc = "individual",
#' missing_handling = "two-step-em")
#'
#' # View estimated partial correlations
#' result_mis$pcor
neighborhood_net <- function(data = NULL, ns = NULL, mat = NULL, n_calc = "individual", k = "log(n)",
                             ordered = "adapted", pcor_merge_rule = "and",
                             missing_handling = "two-step-em",
                              nimp = 20, imp_method = "pmm"){

  n_calc <- match.arg(tolower(n_calc), choices =c("average", "individual", "max", "total"))
  missing_handling <- match.arg(tolower(missing_handling), choices = c("two-step-em", "stacked-mi", "pairwise", "listwise"))
  pcor_merge_rule <- match.arg(tolower(pcor_merge_rule), choices = c("and", "or"))

  # Check: Which input is provided?
  checker(data = data, mat = mat)

  if (!is.null(data)) {

    cor_out <- cor_calc(data = data, missing_handling = missing_handling,
                        ordered = ordered, nimp = nimp, imp_method = imp_method)
    list2env(cor_out, envir = environment())

    if (is.null(ns)){
      if (!(n_calc %in% c("individual", "average", "max", "total"))){
        stop("Invalid n_calc value. Choose from 'individual', 'average', 'max', or 'total'.")
      }
      ns <- calculate_sample_size(data = data, n_calc = n_calc)
    } else checker(ns = ns, data = data)
  } else if (!is.null(mat)) {
    if (is.null(ns)) {
      stop("If 'mat' is provided, 'ns' must also be specified.")
    }
    checker(ns = ns, mat = mat)
    mat <- stats::cov2cor(mat)
    nimp <- missing_handling <- cor_method <- NULL
  }

  mod <- neighborhood_sel(mat = mat, ns = ns, k = k, pcor_merge_rule = pcor_merge_rule)

  result <- list(
    pcor = mod$partials,
    betas = mod$beta_mat,
    ns = ns,
    args = list(k = k, cor_method = cor_method, pcor_merge_rule = pcor_merge_rule,
                missing_handling = missing_handling, nimp = nimp, imp_method = imp_method)
  )

  class(result) <- c("mantar_network")
  return(result)

}



neighborhood_sel <- function(mat, ns, k, pcor_merge_rule){

  class(mat) <- "matrix"              # be sure that mat is of class matrix
  checker(ns = ns, mat = mat)


  p <- ncol(mat)                      # number of variables
  beta_mat <- matrix(NA, p, p)        # initialize matrix for regression coefficients
  colnames(beta_mat) <- rownames(beta_mat) <- colnames(mat)

  # calculate sample size for each variable


  for (dep in 1:p){
    possible_preds <- (1:p)[-dep]  # possible predictors for dependent variable
    n <- ns[dep]                # sample size for dependent variable
    # perform neighorhood selection for dependent node
    mod <- pred_search(mat = mat, dep_ind = dep, possible_pred_ind = possible_preds,
                  n = n, k = k)
    # store regression coefficients in beta matrix
    beta_mat[dep, mod$actual_preds] <- mod$actual_betas
  }

  # compute partial correlations from the resulting beta matrix
  partials <- compute_partials(betas = beta_mat, rule = pcor_merge_rule)

  return(list(partials = partials$wadj, beta_mat = beta_mat))

}


compute_partials <- function(betas, rule){

  if (rule == "or"){
    # replace missing entry with estimated entry on the opposing side of the diagonal
    betas[upper.tri(betas)][is.na(betas[upper.tri(betas)])]  <- betas[lower.tri(betas)][is.na(betas[upper.tri(betas)])]
    betas[lower.tri(betas)][is.na(betas[lower.tri(betas)])]  <- betas[upper.tri(betas)][is.na(betas[lower.tri(betas)])]
  }
  Dummy <- betas * t(betas) # values for the partial correlations

  Dummy[Dummy < 0] <- 0  # look for entries with differing signs of beta
  Dummy[Dummy > 1] <- 1  # bounded between 0 and 1

  wadj <- sign(betas) * sqrt(Dummy) # add the corresponding sign to the partial correlation
  wadj[is.na(wadj)] <- 0 # replace NA with 0
  adj <- ifelse(wadj == 0, 0, 1)

  colnames(wadj) <- rownames(wadj) <- colnames(adj) <- rownames(adj) <- colnames(betas)

  return(list(wadj = wadj, adj = adj))
}
