#' @title
#' Regularized Network Estimation
#'
#' @description
#' Estimate cross-sectional network structures using regularization. The function first
#' computes the correlations (if needed),
#' constructs a grid of tuning parameters tailored to the chosen penalty,
#' and then selects the final network by minimizing a user‑specified
#' information criterion.
#'
#' @param data Optional raw data matrix or data frame containing the variables
#' to be included in the network. May include missing values. If `data` is not
#' provided (`NULL`), a covariance or correlation matrix must be supplied in `mat`.
#' @param ns Optional numeric sample size specification. Can be a single value
#' (one sample size for all variables) or a vector (e.g., variable-wise sample
#' sizes). When `data` is provided and `ns` is `NULL`, sample sizes are derived
#' automatically from `data`. When `mat` is supplied instead of raw data,
#' `ns` must be provided and should reflect the sample size underlying `mat`.
#' @param mat Optional covariance or correlation matrix for the variables to be
#' included in the network. Used only when `data` is `NULL`. If both `data` and
#' `mat` are supplied, `mat` is ignored. When `mat` is used, `ns` must also be
#' provided.
#' @param likelihood Character string specifying how the log-likelihood
#' is computed. Possible values are:
#' \describe{
#'   \item{`"obs_based"`}{Uses the observed-data log-likelihood.}
#'   \item{`"mat_based"`}{Uses log-likelihood based on the sample correlation matrix.}
#' }
#' @param n_calc Character string specifying how the effective sample size is
#' determined. When `data` are provided, it controls how the observation counts
#' across variables are aggregated. When `ns` is a vector, it controls how the
#' entries of `ns` are combined. If both `data` and `ns` are supplied, the
#' values in `ns` take precedence. This argument is ignored when `ns` is a
#' single numeric value. Possible values are:
#' \describe{
#'   \item{`"average"`}{Uses the average sample size across variables or across
#'   the entries of `ns`.}
#'   \item{`"max"`}{Uses the maximum sample size across variables or across
#'   the entries of `ns`.}
#'   \item{`"total"`}{Uses the total number of observations. Applicable only when
#'   `ns` is not provided by the user.}
#' }
#' @param count_diagonal Logical; should observations contributing to the
#' diagonal elements be included when computing the sample size? Only relevant
#' when `data` is provided and `n_calc = "average"`.
#' @param k Penalty term per parameter (number of non-zero partial correlations)
#' used in the network information criterion. The default `k = "log(n)"`,
#' where `n` is the sample size, corresponds to the classical BIC. Setting
#' `k = "2"` yields the classical AIC.
#' @param extended Logical; should the penalty in the information criterion be
#' extended (e.g., using the EBIC instead of the BIC)? If `NULL`, the default is
#' `TRUE` when `penalty = "glasso"` and `FALSE` otherwise.
#' @param extended_gamma Numeric gamma parameter used in the extended information
#' criterion calculation. Only relevant when `extended = TRUE`
#' @param penalty Character string indicating the type of penalty used for
#' regularization. Available options are described in the Details section.
#' @param vary Character string specifying which penalty parameter(s) are varied
#' during regularization to determine the optimal network. Possible values are
#' `"lambda"`, `"gamma"`, or `"both"`.
#' @param n_lambda Number of lambda values to be evaluated. If not specified,
#' the default is 100 when `penalty = "glasso"` and 50 otherwise.
#' @param lambda_min_ratio Ratio of the smallest to the largest lambda value.
#' @param n_gamma Number of gamma values to be evaluated.
#' @param pen_diag Logical; should the diagonal elements be penalized in the
#' regularization process?
#' @param lambda Optional user-specified vector of lambda values.
#' @param gamma Optional user-specified vector of gamma values.
#' @param ordered Logical vector indicating which variables in `data` are
#' treated as ordered (ordinal). Only used when `data` is provided. If a single
#' logical value is supplied, it is recycled to the length of `data`.
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
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{pcor}{Estimated partial correlation matrix corresponding to the
#'   selected (optimal) network.}
#'   \item{n}{Effective sample size used, either supplied
#'   directly via `n` or derived based on `n_calc`.}
#'   \item{cor_method}{Correlation estimation method used for each
#'   variable pair.}
#'   \item{full_results}{Full set of results returned by the model selection
#'   procedure, including all evaluated networks and their fit statistics.}
#'   \item{args}{A list of settings used in the estimation procedure.}
#' }
#'
#' @import mathjaxr
#'
#' @details
#' \loadmathjax
#' \strong{Penalties}
#'
#' This function supports a range of convex and nonconvex penalties for
#' regularized network estimation.
#'
#' For convex penalties, the graphical lasso can be used via
#' `penalty = "glasso"` \insertCite{friedman.2008}{mantar}.
#'
#' Another option is the adaptive lasso, specified with
#' `penalty = "adapt"`.
#' By default, it employs \mjseqn{\gamma = 0.5} \insertCite{zou.2008}{mantar}.
#' Smaller values of \mjseqn{\gamma} (i.e., \mjseqn{\gamma \to 0}) correspond
#' to stronger penalization, whereas \mjseqn{\gamma = 1} yields standard
#' \mjseqn{\ell_1} regularization.
#'
#' The available nonconvex penalties follow the work of
#' \insertCite{williams.2020;textual}{mantar}, who identified the
#' atan penalty as particularly promising. It serves as the default in this
#' implementation because it has desirable theoretical properties, including
#' consistency in recovering the true model as \eqn{n \to \infty}.
#' Additional nonconvex penalties are included for completeness. These were
#' originally implemented in the now–deprecated `R` package \pkg{GGMncv} \insertCite{williams.2021}{mantar},
#' and the implementation in \pkg{mantar} is based on the corresponding methods
#' from that package.
#'
#' Several algorithms exist for nonconvex regularized network estimation.
#' In \pkg{mantar}, we use the one-step estimator of \insertCite{zou.2008;textual}{mantar}
#' because of its computational efficiency and its good performance in settings
#' where \eqn{n > p}, which is typically the case in psychological research.
#'
#'
#' \itemize{
#'
#' \item \strong{Atan}: `penalty = "atan"` \insertCite{wang.2016}{mantar}.
#'  This is currently the default.
#'
#'
#'   \item \strong{Exponential}: `penalty = "exp"`
#'   \insertCite{wang.2018}{mantar}.
#'
#'   \item \strong{Log}: `penalty = "log"`
#'   \insertCite{mazumder.2011}{mantar}.
#'
#'   \item \strong{MCP}: `penalty = "mcp"`
#'   \insertCite{zhang.2010}{mantar}.
#'
#'   \item \strong{SCAD}: `penalty = "scad"`
#'   \insertCite{fan.2001}{mantar}.
#'
#'   \item \strong{Seamless } \mjseqn{\ell_0}: `penalty = "selo"`
#'   \insertCite{dicker.2013}{mantar}.
#'
#'   \item \strong{SICA}: `penalty = "sica"`
#'   \insertCite{lv.2009}{mantar}.
#'
#' }
#'
#' \strong{Conditional Defaults}
#'
#' By default, some tuning parameters depend on the chosen penalty.
#' Specifically, when `penalty = "glasso"`, the number of lambda
#' values `n_lambda` defaults to `100` and `extended`
#' defaults to `TRUE`. For all other penalties, the defaults are
#' `n_lambda = 50` and `extended = FALSE`. These defaults can
#' be overridden by specifying `n_lambda` and/or `extended`
#' explicitly.
#'
#' \strong{Missing Handling}
#'
#' To handle missing data, the function offers two approaches: a two-step expectation-maximization
#' (EM) algorithm and stacked multiple imputation. According to simulations by \insertCite{nehler.2025;textual}{mantar},
#'  stacked multiple imputation performs reliably across a range of sample sizes. In contrast,
#' the two-step EM algorithm provides accurate results primarily when the sample size is large relative
#' to the amount of missingness and network complexity - but may still be preferred in such cases due
#' to its much faster runtime.
#' Currently, the function only supports variables that are directly included in the network analysis;
#' auxiliary variables for missing handling are not yet supported. During imputation, all variables
#' are imputed by default using predictive mean matching \insertCite{@see e.g., @vanbuuren.2018}{mantar},
#' with all other variables in the data set serving as predictors.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' # Estimate regularized network from full data set
#' # Using observed-data loglikelihood and atan penalty
#' result <- regularization_net(mantar_dummy_full_cont,
#'                             likelihood = "obs_based",
#'                             penalty = "atan")
#'
#' # View estimated partial correlation network
#' result$pcor
#'
#' # Estimate regularized network from data set with missings
#' # Using correlation-matrix-based loglikelihood, glasso penalty,
#' # and stacked multiple imputation to handle missings
#' result <- regularization_net(mantar_dummy_mis_mix,
#'                            likelihood = "mat_based",
#'                            penalty = "glasso",
#'                            missing_handling = "stacked-mi",
#'                            ordered = c(FALSE,FALSE,TRUE,TRUE,
#'                                        FALSE,FALSE,TRUE,TRUE))
#'
#' # View used correlation method and effective sample size
#' result$cor_method
#' result$n
#' # View estimated partial correlation network
#' result$pcor
#'
#'
#' @export
#'
regularization_net <- function(data = NULL, ns = NULL, mat = NULL,
                               likelihood = "obs_based", n_calc = "average", count_diagonal = TRUE,
                               k = "log(n)", extended = NULL,
                               extended_gamma = 0.5,
                               penalty = "glasso",
                               vary = "lambda", n_lambda = NULL,
                               lambda_min_ratio = 0.01,
                               n_gamma = 1, pen_diag = FALSE,
                               lambda = NULL, gamma = NULL,
                               ordered = FALSE,
                               missing_handling = "two-step-em",
                               nimp = 20, imp_method = "pmm", ...) {


  # Check: Valid penalty type, likelihood
  if (!penalty %in% c("atan", "mcp", "scad", "exp", "selo",
                      "log", "glasso", "sica", "lq", "adapt")){
    stop("Invalid penalty type. Choose from 'atan', 'mcp', 'scad', 'exp', 'selo',
         'log', 'glasso', 'sica', 'lq', or 'adapt'.")
  }
  likelihood <- match.arg(tolower(likelihood),
                          choices = c("obs_based", "mat_based"))

  # Set conditional defaults for the extended and n_lambda arguments
  if (is.null(extended)) {
    extended <- identical(penalty, "glasso")
  }
  if (is.null(n_lambda)) {
    n_lambda <- if (identical(penalty, "glasso")) 100 else 50
  }

  # Check: Which input is provided?
  checker(data = data, mat = mat)

  # If observed-data loglikelihood is used, scale data to have mean 0 and sd 1; data must be provided
  # and all variables must be treated as continuous
  if (likelihood == "obs_based"){
    if (any(ordered)){
    stop("Calculation of the observed data loglikelihood is only implemented for data
         treated as continuous. ")
    } else if (is.null(data)){
      stop("Calculation of the observed data loglikelihood is only implemented when data is provided.")
    }
    data <- as.data.frame(scale(data))
  }

  # set up means object
  means <- NULL

  # Calculate correlation matrix if data is provided
  if (!is.null(data)) {

    cor_out <- cor_calc(data = data, missing_handling = missing_handling,
                        ordered = ordered, nimp = nimp, imp_method = imp_method, ...)
    list2env(cor_out, envir = environment())
    list2env(cor_out$args, envir = environment())

    # Determine effective sample size if not provided; check that ns is valid if provided
    if (is.null(ns)){
      if (!(n_calc %in% c("average", "max", "total"))){
        stop("Invalid n_calc value. Choose from 'average', 'max', or 'total'.")
      }
      n <- mat_calculate_sample_size(data = data, n_calc = n_calc, count_diagonal = count_diagonal)
    } else checker(ns = ns, data = data)
  } else if (!is.null(mat)) {
    if (is.null(ns)) {
      stop("If 'mat' is provided, 'ns' must also be specified.")
    }
    # calculate n based on ns and n_calc
    if (n_calc == "average"){
      n <- mean(ns)
    } else if (n_calc == "max"){
      n <- max(ns)
    } else {
      stop("When 'mat' is provided, 'n_calc' must be one of 'average' or 'max'.")
    }
    # Ensure mat is a correlation matrix
    mat <- stats::cov2cor(mat)
    nimp <- missing_handling <- cor_method <- NULL
  }

  # Call helper function to perform regularization and model selection
  mod <- regularization_sel(
    mat = mat,
    data = data,
    means = means,
    n = n,
    k = k,
    likelihood = likelihood,
    extended = extended,
    extended_gamma = extended_gamma,
    penalty = penalty,
    vary = vary,
    n_lambda = n_lambda,
    lambda_min_ratio = lambda_min_ratio,
    n_gamma = n_gamma,
    pen_diag = pen_diag,
    lambda = lambda,
    gamma = gamma
  )

  # Prepare output
  result <- list(
    pcor = mod$opt_net,
    n = n,
    cor_method = cor_method,
    full_results = mod$full_results,
    args = list(likelihood = likelihood, n_calc = n_calc, count_diagonal = count_diagonal,
                k = k, extended = extended, extended_gamma = extended_gamma,
                penalty = penalty, vary = vary, n_lambda = n_lambda,
                lambda_min_ratio = lambda_min_ratio, n_gamma = n_gamma,
                pen_diag = pen_diag, lambda = lambda, gamma = gamma,
                ordered = ordered,
                missing_handling = missing_handling,
                nimp = nimp, imp_method = imp_method)
  )

  # Set class and return
  class(result) <- c("mantar_regularization", "mantar_network")
  return(result)

}


#' @title Regularized network analysis helper function
#'
#' @param mat Covariance or correlation matrix.
#' @param data Optional raw data matrix or data frame containing the variables
#' @param means Optional vector of variable means.
#' @param n Effective sample size.
#' @param k Penalty term per parameter used in the information criterion.
#' @param likelihood Character string specifying how the log-likelihood is computed.
#' @param extended Logical; should the penalty in the information criterion be extended?
#' @param extended_gamma Numeric gamma parameter used in the extended information criterion calculation.
#' @param penalty Character string indicating the type of penalty used for regularization.
#' @param vary Character string specifying which penalty parameter(s) are varied during regularization.
#' @param n_lambda Number of lambda values to be evaluated.
#' @param lambda_min_ratio Ratio of the smallest to the largest lambda value.
#' @param n_gamma Number of gamma values to be evaluated.
#' @param pen_diag Logical; should the diagonal elements be penalized in the regularization process?
#' @param lambda Optional user-specified vector of lambda values.
#' @param gamma Optional user-specified vector of gamma values.
#'
#' @returns
#' \describe{
#'   \item{opt_net}{Optimal network based on the information criterion.}
#'   \item{full_results}{A list of fitted models (one per tuning-parameter
#'                       combination), each containing the estimated
#'                       precision matrix `wi`, the covariance matrix
#'                       `w`, the information criterion `ic`, and
#'                       the corresponding penalty matrix `rho_mat`.}
#' }
#' @noRd
regularization_sel <- function(mat, data = NULL, means = NULL, n, k,
                               likelihood, extended, extended_gamma,
                               penalty = "glasso", vary = "lambda",
                               n_lambda = 50,
                               lambda_min_ratio = 0.01,
                               n_gamma = 1, pen_diag = FALSE,
                               lambda = NULL, gamma = NULL) {


  # be sure that mat is of class matrix
  class(mat) <- "matrix"
  # store variable names
  varnames <- colnames(mat)
  # check combination of inputes to be valid
  checker(n = n, mat = mat)
  # check that data is provided if obs_based likelihood is used
  if (is.null(data) & likelihood == "obs_based"){
    stop("Calculation of the observed data loglikelihood is only implemented when data is provided.")
  }

  # be sure to only handle correlation matrices
  mat <- stats::cov2cor(mat)

  # define penalty matrices based on specified penalty and tuning parameters
  pen_results <- def_pen_mats(
    mat = mat,
    penalty = penalty,
    vary = vary,
    n_lambda = n_lambda,
    lambda_min_ratio = lambda_min_ratio,
    n_gamma = n_gamma,
    n = n,
    pen_diag = pen_diag,
    lambda = lambda,
    gamma = gamma
  )

  # process separately for all penalty matrices
  results <- lapply(pen_results$pen_mats, function(rho_mat) {
    # estimate network for given penalty matrix
    cand_net <- glassoFast::glassoFast(S = mat, rho = rho_mat)
    dimnames(cand_net$wi) <- list(varnames, varnames)
    dimnames(cand_net$w)  <- list(varnames, varnames)

    # calculate information criterion for the candidate network
    cand_net$ic <- mat_ic_calc(data = data,
                               sample_cor = mat,
                               theta = cand_net$wi,
                               mu = means,
                               n = n,
                               k = k,
                               extended = extended,
                               extended_gamma = extended_gamma,
                               likelihood = likelihood)
    # store penalty matrix
    cand_net$rho_mat <- rho_mat
    return(cand_net)
  })

  # select optimal network based on information criterion
  opt_wi <- results[[ which.min(
    vapply(results, function(x) x$ic, numeric(1L))
  ) ]]$wi
  opt_net <- inv_to_net(theta = opt_wi)

  # return optimal network and full results
  return(list(
    opt_net = opt_net,
    full_results = results
  ))

}


#' @title Create penalty matrix for regularization based on specified penalty
#'
#' @description Construct a grid of penalty parameters (`lambda`, `gamma`)
#' and the corresponding penalty matrices for a variety of penalties used
#' in graphical model estimation.
#'
#' @param mat Correlation or covariance matrix.
#' @param penalty penalty type used for regularization. Supported optons are
#' "atan", "mcp", "scad", "exp", "selo", "log", "glasso", "sica", "lq", and "adapt".
#' @param vary Character string specifying which penalty parameter(s) are varied
#' during regularization to determine the optimal network. Possible values are
#' `"lambda"`, `"gamma"`, or `"both"`.
#' @param n_lambda Number of lambda values to be evaluated when varying lambda.
#' @param lambda_min_ratio Ratio of the smallest to the largest lambda value.
#' @param n_gamma Number of gamma values to be evaluated when varying gamma.
#' @param n Sample size used for calculating default lambda when not varying lambda.
#' @param pen_diag Logical; should the diagonal elements be penalized in the
#' regularization process?
#' @param lambda Optional user-specified vector of lambda values.
#' @param gamma Optional user-specified vector of gamma values.
#'
#' @returns
#' A list with two elements:
#' \describe{
#'   \item{grid}{A data frame with all combinations of `lambda` and `gamma`
#'               used to construct penalty matrices.}
#'   \item{pen_mats}{A list of penalty matrices, one for each row in `grid`.}
#' }
#' @noRd
def_pen_mats <- function(mat,
                         penalty = "glasso",
                         vary =  "lambda",
                         n_lambda = 50,
                         lambda_min_ratio = 0.01,
                         n_gamma = 1,
                         n = NULL,
                         pen_diag = FALSE,
                         lambda = NULL,
                         gamma  = NULL) {

  # Be sure to only handle correlation matrices - set up precision and identity matrices
  mat <- stats::cov2cor(mat)
  theta <- solve(mat)
  p <- ncol(theta)
  I_p <- diag(p)



  # create values for lambda when varying lambda and it is not user-specified
  if (!is.null(lambda)) {
    lambda_vec <- lambda
    message("Using user-specified lambda values.")
  } else if (vary %in% c("lambda", "both")) {

    if (vary == "lambda" & n_gamma > 1) {
      warning("Varying 'lambda' only, but n_gamma > 1. Using only 1 gamma value.")
    }

    # largest lambda value is the largest off-diagonal absolute value in the correlation matrix
    lambda_max <- max(max(mat - I_p), -min(mat - I_p))
    # minimal lambda value at the specified ration
    lambda_min <- lambda_min_ratio * lambda_max
    # space lambda values on log scale
    lambda_vec <- exp(seq(log(lambda_max), log(lambda_min), length.out = n_lambda)) |> sort()

  } else {
    # if lambda is not varied and not user-specified, set to default value
    lambda_vec <- sqrt(log(p) / n)
  }

  # create values for gamma when varying gamma and it is not user-specified
  if (!is.null(gamma)) {
    gamma_vec <- gamma
    message("Using user-specified gamma values.")
  } else if (vary %in% c("gamma", "both")) {

    if (vary == "gamma" & n_lambda > 1) {
      warning("Varying 'gamma' only, but n_lambda > 1. Using only 1 lambda value.")
    }
    # Defaults for multiple gamma values based on penalty type
    if (penalty == "scad") {
      gamma_vec <- seq(2.001, 5, n_gamma)
    } else if (penalty == "mcp") {
      gamma_vec <- seq(1.001, 4, n_gamma)
    } else if (penalty == "adapt") {
      gamma_vec <- seq(0.1, 1, n_gamma)
    } else {
      gamma_vec <- seq(0.001, 0.1, length.out =  n_gamma)
    }
    if (penalty == "glasso") {
      warning("The glasso penalty does not use a gamma parameter. The gamma grid will be ignored.")
    }

  } else {
    # if gamma is not varied and not user-specified, set to default value  based on penalty type
    if (penalty == "scad") {
      gamma_vec <- 3.7
    } else if (penalty == "mcp") {
      gamma_vec <- 2
    } else if (penalty == "adapt") {
      gamma_vec <- 0.5
    } else {
      gamma_vec <- 0.01
    }
  }

 # create grid of (lambda, gamma) combinations
  grid <- expand.grid(
    lambda = lambda_vec,
    gamma  = gamma_vec,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  # create penalty matrices for each (lambda, gamma) combination
  pen_mats <- lapply(
    X = seq_len(nrow(grid)),
    FUN = function(i, theta) {
      lam <- grid$lambda[i]
      gam <- grid$gamma[i]
      p   <- ncol(theta)

      if (penalty == "glasso") {

        # for glasso, penalty matrix is simply lambda times a matrix of ones
        pen_mat <- matrix(lam, p, p)

      # other penalties differ in their definition
      } else if (penalty == "atan") {

        theta_safe <- abs(theta)
        pen_mat <- lam * ((gam * (gam + 2/pi)) / (gam^2 + theta_safe^2))

      } else if (penalty == "scad") {

        theta_safe <- abs(theta)
        pen_mat <- lam * (theta_safe <= lam) +
          (pmax(gam * lam - theta_safe, 0) / (gam - 1)) * (theta_safe > lam)

      } else if (penalty == "mcp") {

        theta_safe <- abs(theta)
        pen_mat <- (theta_safe <= gam * lam) * (lam - theta_safe / gam)

      } else if (penalty == "exp") {

        theta_safe <- abs(theta)
        pen_mat <- (lam / gam) * exp(-(theta_safe / gam))

      } else if (penalty == "selo") {

        if (!requireNamespace("numDeriv", quietly = TRUE)) {
          stop("Package 'numDeriv' must be installed to use selo penalty.", call. = FALSE)
        }

        theta_safe <- abs(theta)
        theta_safe <- ifelse(theta == 0, 1e-5, theta_safe)

        selo_pen <- function(x, lam, gam = 0.01) {
          x <- abs(x)
          (lam / log(2)) * log((x / (x + gam)) + 1)
        }

        pen_mat <- matrix(
          numDeriv::grad(selo_pen, x = theta_safe, lam = lam, gam = gam),
          nrow = p, ncol = p
        )

      } else if (penalty == "sica") {

        if (!requireNamespace("numDeriv", quietly = TRUE)) {
          stop("Package 'numDeriv' must be installed to use sica penalty.", call. = FALSE)
        }

        theta_safe <- abs(theta + 1e-4)

        sica_pen <- function(x, lam, gam){
          x <- abs(x)
          lam * (((gam + 1) * x) /(x + gam))
        }

        pen_mat <- matrix(
          numDeriv::grad(sica_pen, x = theta_safe, lam = lam, gam = gam),
          nrow = p, ncol = p
        )

      } else if (penalty == "log") {

        if (!requireNamespace("numDeriv", quietly = TRUE)) {
          stop("Package 'numDeriv' must be installed to use log penalty.", call. = FALSE)
        }

        theta_safe <- abs(theta + 1e-4)

        log_pen <- function(x, lam, gam = 0.01){
          x <- abs(x)
          gam_inv <- 1/gam
          (lam / log(gam_inv + 1)) * log(gam_inv * x + 1)
        }

        pen_mat <- matrix(
          numDeriv::grad(log_pen, x = as.numeric(theta_safe), lam = lam, gam = gam),
          nrow = p, ncol = p
        )

      } else if (penalty == "lq") {

        theta_safe <- abs(theta)
        epsilon <- 1e-4
        pen_mat <- lam * gam * (theta_safe + epsilon)^(gam - 1)

      } else if (penalty == "adapt") {

        theta_safe <- abs(theta + 1e-4)
        pen_mat <- lam * (theta_safe)^(-(1 - gam))
      }

      pen_mat
    },
    theta = theta
  )


  # set penalty for diagonal to 0 if pen_diag = FALSE
  if (!pen_diag) {
    pen_mats <- lapply(pen_mats, function(mat) {
      diag(mat) <- 0
      return(mat)
    })
  }

  list(
    grid     = grid,
    pen_mats = pen_mats
  )
}





