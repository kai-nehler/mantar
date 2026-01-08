#' @title Stepwise Multiple Regression Model Search based on Information Criteria
#'
#' @description
#' Performs stepwise model selection for multiple regression using information
#' criteria to identify the optimal regression model.
#'
#' @param data Raw data matrix or data frame containing the variables to be
#' included in the regression models. May include missing values. If `data` is
#' `NULL`, a covariance or correlation matrix must be supplied in `mat`.
#' @param n Numeric value specifying the sample size used in calculating the
#' information criteria. If not provided, it is derived from `data`. When `mat`
#' is supplied instead of raw data, `n` must be provided.
#' @param mat Optional covariance or correlation matrix for the variables to be
#' included in the regression. Used only when `data` is `NULL`.
#' @param dep_ind Index of the column in `data` to be used as the dependent
#' variable in the regression model.
#' @param n_calc Character string specifying how the sample size is calculated
#' when `n` is not provided. Possible values are:
#' \describe{
#'   \item{`"individual"`}{Uses the number of non-missing observations for the
#'   variable used as the dependent variable.}
#'   \item{`"average"`}{Uses the average number of non-missing observations
#'   across all variables.}
#'   \item{`"max"`}{Uses the maximum number of non-missing observations across
#'   all variables.}
#'   \item{`"total"`}{Uses the total number of rows in `data`.}
#' }
#' @param ic_type Type of information criterion to compute for model selection.
#' Options are `bic` (default), `aic`, `aicc`.
#' @param ordered Logical vector indicating whether each variable in `data`
#' should be treated as ordered categorical when computing the correlation
#' matrix. If a single logical value is supplied, it is recycled to all
#' variables. Only used when `data` is provided.
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
#' @param ... Further arguments passed to internal functions.
#'
#' @return A list with the following elements:
#' \describe{
#'  \item{regression}{Named vector of regression coefficients for the dependent variable.}
#'  \item{R2}{R-squared value of the regression model.}
#'  \item{n}{Sample size used in the regression model.}
#'  \item{args}{List of settings used in the regression model.}
#'  }
#'
#' @details
#' \loadmathjax
#'
#' This function performs stepwise model selection for multiple regression
#' using information criteria. It was originally developed as a component of
#' the neighborhood selection framework for network estimation
#' \insertCite{nehler.2024}{mantar}, where each node-wise regression model is
#' selected individually. However, the procedure can also be used as a
#' standalone tool for exploratory regression model search, particularly in
#' settings with missing data. Unlike standard stepwise regression functions,
#' this implementation explicitly supports missing-data handling strategies,
#' making it suitable for situations in which classical methods fail or produce
#' biased results.
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
#' @references
#' \insertAllCited{}
#'
#' @export
#'
#' @examples
#' # For full data using AIC
#' # First variable of the data set as dependent variable
#' result <- regression_opt(
#'   data = mantar_dummy_full_cont,
#'   dep_ind = 1,
#'   ic_type = "aic"
#' )
#'
#' # View regression coefficients and R-squared
#' result$regression
#' result$R2
#'
#' # For data with missingess using BIC
#' # Second variable of the data set as dependent variable
#' # Using individual sample size of the dependent variable and stacked Multiple Imputation
#'
#' result_mis <- regression_opt(
#'  data = mantar_dummy_mis_cont,
#'  dep_ind = 2,
#'  n_calc = "individual",
#'  missing_handling = "two-step-em",
#'  ic_type = "bic"
#'  )
#'
#'  # View regression coefficients and R-squared
#'  result_mis$regression
#'  result_mis$R2
regression_opt <- function(data = NULL, n = NULL, mat = NULL, dep_ind,
                           n_calc = "individual", ic_type = "bic", ordered = FALSE,
                           missing_handling = "stacked-mi", nimp = 20, imp_method = "pmm", ...) {

  # Argument checks
  ic_type <- match.arg(tolower(ic_type), choices = c("aic", "bic", "aicc"))
  n_calc <- match.arg(tolower(n_calc), choices =c("average", "individual", "max", "total"))
  missing_handling <- match.arg(tolower(missing_handling), choices = c("two-step-em", "stacked-mi", "pairwise", "listwise"))

  # Check: Which input is provided?
  checker(data = data, mat = mat)

  # Determine column names
  col_names <- if (!is.null(data)) colnames(data) else colnames(mat)

  # Resolve dep_ind to a column index
  # if it is given as a number, check if it is valid
  if (is.character(dep_ind)) {
    if (is.null(col_names)) {
      stop("dep_ind was given as a name, but the selected input (data or mat) has no column names.")
    }
    if (!(dep_ind %in% col_names)) {
      stop("dep_ind must be a valid column name.")
    }
    dep_ind <- which(col_names == dep_ind)
  } else if (is.numeric(dep_ind)) {
    max_cols <- if (!is.null(col_names)) length(col_names) else if (!is.null(data)) ncol(data) else ncol(mat)
    if (dep_ind < 1 || dep_ind > max_cols) {
      stop("dep_ind must be a valid column index.")
    }
  } else {
    stop("dep_ind must be either a numeric index or a character string.")
  }

  # Calculate correlation matrix and sample size if data is provided
  if (!is.null(data)) {

    cor_out <- cor_calc(data = data, ordered = ordered,
                        missing_handling = missing_handling, nimp = nimp, imp_method = imp_method, ...)
    list2env(cor_out, envir = environment())

    if (is.null(n)){
      # Calculate sample sizes
      ns <- reg_calculate_sample_size(data = data, n_calc = n_calc)
      # use the sample size in the dep_ind place
      n <- ns[dep_ind]
    } else {checker(n)}

  # if mat is provided, transform to correlation matrix and check for n
  } else if (!is.null(mat)) {
    if (is.null(n)) {
      stop("If 'mat' is provided, 'n' must also be specified.")
    }
    if (length(n) != 1) {
      stop("Length of 'n' must be 1.")
    }
    mat <- stats::cov2cor(mat)
    nimp <- missing_handling <- cor_method <- imp_method <- NULL
  }

  # Search for the optimal regression model
  mod <- pred_search(mat = mat, dep_ind = dep_ind,
                         possible_pred_ind = setdiff(1:ncol(mat), dep_ind),
                     n = n, ic_type = ic_type)

  # Retrive the estimated beta weights and assign names
  regression <- mod$actual_betas
  names(regression) <- colnames(mat)[mod$actual_preds]

  # Prepare output
  result <- list(
    regression = regression,
    R2 = 1 - mod$actual_resid_var,
    n = n,
    args = list(ic_type = ic_type, cor_method = cor_method, missing_handling = missing_handling,
                nimp = nimp, imp_method = imp_method)
  )

  # Assign class and return
  class(result) <- c("mantar_regression")
  return(result)
}



#' @title Helper for stepwise predictor search based on information criteria
#'
#' @description Internal helper that performs a greedy stepwise search for
#' predictors of a given dependent variable, using an information criterion
#' (e.g., AIC or BIC) as selection rule. Predictors are added or removed one
#' at a time based on whether they improve the information criterion.
#'
#' @param mat Correlation matrix
#' @param dep_ind Index of the dependent variable in of `mat`
#' @param possible_pred_ind Vector of indices for the possible predictor variables in `mat`
#' @param n Sample size for information criteria calculation
#' @param ic_type Type of information criterion to compute.
#' Options: "AIC", "BIC", "AICc"
#'
#' @returns
#' A list with the following elements:
#' \describe{
#'   \item{actual_preds}{Integer vector with the indices of the selected predictors.}
#'   \item{actual_betas}{Numeric vector with the corresponding regression coefficients.}
#'   \item{actual_resid_var}{Residual variance of the final model.}
#'   \item{best_IC}{Information criterion value of the final model.}
#' }
#' @noRd
pred_search <- function(mat, dep_ind, possible_pred_ind, n, ic_type = "bic"){

  # Compute Null Model and the corresponding Information Criterion
  null_mod <- matrix_regression(mat = mat, dep_ind = dep_ind, pred_ind = NULL)
  null_IC <- reg_ic_calc(resid_var = null_mod$resid_var, n = n, n_preds = 0, ic_type = ic_type)

  # initialize vector for actual predictors
  actual_preds <- c()
  # initialize vector for actual betas
  actual_betas <- c()
  # initialize object for actual residual variance using null model
  actual_resid_var <- null_mod$resid_var
  # initialize best IC with null model
  best_IC <- null_IC

  while(TRUE){
    mods <- lapply(possible_pred_ind, function(pred){
      # if the predictor is not already in the model
      if (!(pred %in% actual_preds)){
        # add the predictor to the model
        mod_preds <- sort(c(pred, actual_preds))
        # calculate regression coefficients
        mod <- matrix_regression(mat = mat, dep_ind = dep_ind, pred_ind = mod_preds)
        mod_betas <- mod$beta
        mod_resid_var <- mod$resid_var
      }
      # if the predictor is already in the model
      else{
        # remove the predictor from the model
        mod_preds <- sort(actual_preds[actual_preds != pred])
        # if there are no predictors left - empty model
        if (length(mod_preds) == 0){
          mod_preds <- NULL
          mod_betas <- NULL
        } else{
          # if there are predictors left - calculate regression coefficients
          mod <- matrix_regression(mat = mat, dep_ind = dep_ind, pred_ind = mod_preds)
          mod_betas <- mod$beta
          mod_resid_var <- mod$resid_var
        }
      }

      # calculate information criteria for the model
      # if there are no predictors, the information criteria is the one from the null model
      if(is.null(mod_betas)) {
        mod_resid_var <- null_mod$resid_var
        IC <- null_IC
      } else {
        # if there are predictors, calculate the information criteria for the model
        IC <- reg_ic_calc(resid_var = mod_resid_var, n = n, n_preds = length(mod_preds), ic_type = ic_type)
      }

      return(list(IC = IC, mod_preds = mod_preds, mod_betas = mod_betas, mod_resid_var = mod_resid_var))
    })

    # extract information criteria from models
    ICs <-  lapply(mods, "[[", "IC") |> unlist()

    # if there is a model with a lower information criteria than the best model so far
    if (any(ICs< best_IC)){
      # extract predictors and regression coefficients from all models
      mod_betas <-  lapply(mods, "[[", "mod_betas")
      mod_preds <-  lapply(mods, "[[", "mod_preds")
      mod_resid_vars <-  lapply(mods, "[[", "mod_resid_var")

      # find the model with the lowest information criteria
      mod_num <- which.min(ICs)
      best_IC <- ICs[mod_num]

      # extract predictors and regression coefficients from the best model
      actual_preds <- mod_preds[mod_num] |> unlist()
      actual_betas <- mod_betas[mod_num] |> unlist()
      actual_resid_var <- mod_resid_vars[mod_num] |> unlist()
    } else{
      # if there is no model with a lower information criteria than the best model so far
      break
    }
  }

  # return model information for the best model
  return(list(actual_preds = actual_preds,
              actual_betas = actual_betas,
              actual_resid_var = actual_resid_var,
              best_IC = best_IC))

}

