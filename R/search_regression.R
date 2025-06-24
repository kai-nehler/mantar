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
#' @param missing_handling Method for estimating the correlation matrix in the presence of missing data.
#' `"tow-step-em"` uses a classic EM algorithm to estimate the covariance matrix from the data.
#' `"stacked-mi"` uses multiple imputation to estimate the covariance matrix from the data.
#' `"pairwise"` uses pairwise deletion to estimate the covariance matrix from the data.
#' `"listwise"` uses listwise deletion to estimate the covariance matrix from the data.
#' @param k Penalty per parameter (number of predictors + 1) to be used in node-wise regressions; the default log(n) (number of observations observation) is the classical BIC. Alternatively, classical AIC would be `k = 2`.
#' @param nimp Number of multiple imputations to perform when using multiple imputation for missing data (default: 20).
#'
#' @return A list with the following elements:
#' \describe{
#'  \item{regression}{Named vector of regression coefficients for the dependent variable.}
#'  \item{R2}{R-squared value of the regression model.}
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
#'  missing_handling = "stacked-mi",
#'  )
#'
#'  # View regression coefficients and R-squared
#'  result_mis$regression
#'  result_mis$R2
regression_opt <- function(data = NULL, n = NULL, mat = NULL, dep_ind,
                           n_calc = "individual",
                           missing_handling = "two-step-em",
                           k = "log(n)", nimp = 20) {

  n_calc <- match.arg(tolower(n_calc), choices =c("average", "individual", "max", "total"))
  missing_handling <- match.arg(tolower(missing_handling), choices = c("two-step-em", "stacked-mi"))

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
    if (anyNA(data)){
      if (missing_handling == "two-step-em"){

        if (!requireNamespace("lavaan", quietly = TRUE)) {
          stop(
            "Package \"lavaan\" must be installed to use this function.",
            call. = FALSE
          )
        }

        lavobject <- suppressWarnings(try(lavaan::lavCor(data,
                                                         missing = "ml", se = "none", meanstructure = TRUE,
                                                         estimator = "ML", output = "fit")))
        if (inherits(lavobject, "try-error")) stop("lavaan::lavCor failed. Check your data.")
        mat <- try(stats::cov2cor(lavaan::inspect(lavobject, "cov.ov")))

      } else if (missing_handling == "stacked-mi"){

        if (!requireNamespace("mice", quietly = TRUE)) {
          stop(
            "Package \"mice\" must be installed to use this function.",
            call. = FALSE
          )
        }
        original_names <- colnames(data)
        colnames(data) <- make.names(original_names)

        imputed_data <- suppressMessages(mice::mice(data, m = nimp, maxit = 10, method = 'pmm',
                                                    ridge = 0, donors = 5, ls.method = "qr"))
        stacked_data <- mice::complete(imputed_data, c(1:nimp))
        colnames(stacked_data) <- original_names  # restore original column names
        mat <- stats::cor(stacked_data)
      } else if (missing_handling == "pairwise"){
        mat <- stats::cor(data, use = "pairwise.complete.obs")
      } else if (missing_handling == "listwise"){
        mat <- stats::cor(data, use = "complete.obs")
      }
    } else {
      mat <- stats::cor(data)
    }

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
  }

  mod <- pred_search(mat = mat, dep_ind = dep_ind,
                         possible_pred_ind = setdiff(1:ncol(mat), dep_ind),
                     n = n, k = k)

  regression <- mod$actual_betas
  names(regression) <- mod$actual_preds
  R2 <- 1 - mod$actual_resid_var

  return(list(regression = regression,  R2 = R2))
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

