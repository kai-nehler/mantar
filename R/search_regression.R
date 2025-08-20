#' Stepwise Multiple Regression Search based on Information Criteria
#'
#' @param data Raw data containing only the variables to be tested within the multiple regression as dependent or independent variable. May include missing values.
#' @param n Numeric value specifying the sample size used in calculating information criteria for model search.
#' If not provided, it will be computed based on the data.
#' If a correlation matrix (`mat`) is supplied instead of raw data, `n` must be provided.
#' @param mat Optional covariance or correlation matrix for the variables to be used within the multiple regression.
#' #' Used only if \code{data} is \code{NULL}.
#' @param dep_ind Index of the column within a data set to be used as dependent variable within in the regression model.
#' @param n_calc Method for calculating the sample size for node-wise regression models. Can be one of:
#' `"individual"` (sample size for each variable is the number of non-missing observations for that variable),
#' `"average"` (sample size is the average number of non-missing observations across all variables),
#' `"max"` (sample size is the maximum number of non-missing observations across all variables),
#' `"total"` (sample size is the total number of observations across in the data set / number of rows).
#' @param k Penalty per parameter (number of predictors + 1) to be used in node-wise regressions; the default log(n) (number of observations observation) is the classical BIC. Alternatively, classical AIC would be `k = 2`.
#' @param cor_method Two correlation methods are currently implemented, namely polychoric and pearson. The default setting for is `"adapted"`. In this mode, an appropriate method is automatically selected based on the data (number of cases, number of variables, number of categories). This automatic selection can be overridden by explicitly specifying a method (`"polychoric"` or `"pearson"`).
#' @param missing_handling Method for estimating the correlation matrix in the presence of missing data.
#' `"tow-step-em"` uses a classic EM algorithm to estimate the covariance matrix from the data.
#' `"stacked-mi"` uses multiple imputation to estimate the covariance matrix from the data.
#' `"pairwise"` uses pairwise deletion to estimate the covariance matrix from the data.
#' `"listwise"` uses listwise deletion to estimate the covariance matrix from the data.
#' @param nimp Number of multiple imputations to perform when using multiple imputation for missing data (default: 20).
#' @param imp_method Method for multiple imputation when using `"stacked-mi"` for missing data handling. Default is `"pmm"` (predictive mean matching).
#'
#' @return A list with the following elements:
#' \describe{
#'  \item{regression}{Named vector of regression coefficients for the dependent variable.}
#'  \item{R2}{R-squared value of the regression model.}
#'  \item{n}{Sample size used in the regression model.}
#'  \item{args}{List of settings used in the regression model.}
#'  }
#'
#' @export
#'
#' @examples
#' # For full data using AIC
#' # First variable of the data set as dependent variable
#' result <- regression_opt(
#'   data = mantar_dummy_full,
#'   dep_ind = 1,
#'   k = "2"
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
#'  data = mantar_dummy_mis,
#'  dep_ind = 2,
#'  n_calc = "individual",
#'  missing_handling = "two-step-em",
#'  )
#'
#'  # View regression coefficients and R-squared
#'  result_mis$regression
#'  result_mis$R2
regression_opt <- function(data = NULL, n = NULL, mat = NULL, dep_ind,
                           n_calc = "individual", k = "log(n)", cor_method = "adapted",
                           missing_handling = "stacked-mi", nimp = 20, imp_method = "pmm") {

  n_calc <- match.arg(tolower(n_calc), choices =c("average", "individual", "max", "total"))
  missing_handling <- match.arg(tolower(missing_handling), choices = c("two-step-em", "stacked-mi", "pairwise", "listwise"))
  cor_method <- match.arg(tolower(cor_method), choices = c("adapted", "pearson", "polychoric"))

  # Check: Which input is provided?
  checker(data = data, mat = mat)

  # Determine column names relevant for dep_ind
  col_names <- if (!is.null(data)) colnames(data) else colnames(mat)

  # Resolve dep_ind to a column index
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

  if (!is.null(data)) {

    cor_out <- cor_calc(data = data, cor_method = cor_method,
                        missing_handling = missing_handling, nimp = nimp, imp_method = imp_method)
    list2env(cor_out, envir = environment())

    if (is.null(n)){
      ns <- calculate_sample_size(data = data, n_calc = n_calc)
      n <- ns[dep_ind]  # use the sample size for the dependent variable
    } else if (length(n) != 1){
      stop("Length of 'n' must be 1.")
    }

  } else if (!is.null(mat)) {
    if (is.null(n)) {
      stop("If 'mat' is provided, 'n' must also be specified.")
    }
    if (length(n) != 1) {
      stop("Length of 'n' must be 1.")
    }
    mat <- stats::cov2cor(mat)
    nimp <- missing_handling <- NULL
  }


  mod <- pred_search(mat = mat, dep_ind = dep_ind,
                         possible_pred_ind = setdiff(1:ncol(mat), dep_ind),
                     n = n, k = k)

  regression <- mod$actual_betas
  names(regression) <- colnames(mat)[mod$actual_preds]

  result <- list(
    regression = regression,
    R2 = 1 - mod$actual_resid_var,
    n = n,
    args = list(k = k, cor_method = cor_method, missing_handling = missing_handling,
                nimp = nimp, imp_method = imp_method)
  )

  class(result) <- c("mantar_regression")
  return(result)
}



pred_search <- function(mat, dep_ind, possible_pred_ind, n, k = log(n)){

  null_mod <- matrix_regression(mat = mat, dep_ind = dep_ind, pred_ind = NULL)
  null_IC <- reg_ic_calc(resid_var = null_mod$resid_var, n = n, n_preds = 0, k = k)

  actual_preds <- c()     # initialize vector for actual predictors
  actual_betas <- c()     # initialize vector for actual betas
  actual_resid_var <- null_mod$resid_var # initialize object for actual residual variance

  best_IC <- null_IC  # initialize best IC with null model

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
        IC <- reg_ic_calc(resid_var = mod_resid_var, n = n, n_preds = length(mod_preds), k = k)
      }

      return(list(IC = IC, mod_preds = mod_preds, mod_betas = mod_betas, mod_resid_var = mod_resid_var))
    })

    ICs <-  lapply(mods, "[[", "IC") |> unlist()  # extract information criteria from models

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

  return(list(actual_preds = actual_preds, actual_betas = actual_betas, actual_resid_var = actual_resid_var, best_IC = best_IC))

}

