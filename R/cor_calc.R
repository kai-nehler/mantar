#' @title
#' Correlation Matrix Estimation with Support for Multiple Correlation Types
#'
#' @description
#' Computes a correlation matrix from raw data while accounting for missing
#' values through several missing-data handling strategies. Supports different
#' correlation types based on whether variables are treated as ordered.
#'
#' @param data Data frame or matrix containing the variables for which the
#' correlation matrix is to be computed. May include missing values. Variables
#' selected in `network_vars` must be numeric or ordered factors; variables
#' selected in `auxiliary_vars` may also be unordered factors or character
#' vectors (see Details).
#' @param ordered Logical vector indicating whether each variable in
#' `network_vars` should be treated as ordered categorical when computing the
#' correlation matrix. If a single logical value is supplied, it is recycled to
#' all network variables. `ordered` only affects the correlation type, not the
#' imputation; how variables are imputed depends on their class and on
#' `imp_method` (see Details).
#' @param network_vars Optional character or numeric vector specifying the variables
#' for which the correlation matrix is returned, typically the variables used in
#' subsequent network estimation. All selected variables must be numeric or
#' ordered factors. Default is `NULL`, which means that all variables in the
#' provided data are included.
#' @param auxiliary_vars Optional character or numeric vector specifying variables
#' in the data set that are additionally used for missing-data handling, but are
#' not included in the returned correlation matrix. Can only be used together
#' with `network_vars`, when at least one network variable contains missing
#' values, and with `missing_handling = "two-step-em"` or `"stacked-mi"`.
#' Default is `NULL`, which means that no auxiliary variables are used.
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
#' @param imp_method Optional character vector specifying the imputation method
#' for each variable when `missing_handling = "stacked-mi"`. Default is `NULL`,
#' which uses the defaults of [mice::mice()]: predictive mean matching (`"pmm"`)
#' for numeric variables (including ordered factors, which are converted to
#' integer scores), logistic regression (`"logreg"`) for factors with two
#' levels and multinomial logistic regression (`"polyreg"`) for factors with
#' more than two levels. If supplied, `imp_method` must contain one entry per
#' variable in `c(network_vars, auxiliary_vars)`, in this order. Entries for
#' fully observed variables are ignored. To impute an ordinal variable (numeric
#' or ordered factor) as ordered categorical instead of with `"pmm"`, set its
#' entry to `"polr"` (proportional odds model). See
#' [mice::mice()] for further methods. Note that \pkg{mice} does not check every
#' method against the variable type; for example, `"pmm"` applied to an
#' unordered factor treats its categories as ordered.
#' @param maxit Maximum number of iterations for the imputation algorithm when
#' `missing_handling = "stacked-mi"` (default: 10).
#' @param ... Further arguments passed to internal functions. When
#' `missing_handling = "stacked-mi"`, these are passed on to [mice::mice()]
#' (e.g., `seed` or `printFlag`).
#'
#' @details
#' Correlations between the network variables are computed pairwise, depending
#' on `ordered`:
#' * Polychoric correlations for two ordered variables,
#' * Polyserial correlations for one ordered and one continuous variable,
#' * Pearson correlations for two continuous variables.
#'
#' Ordered factors (network or auxiliary variables) are converted to integer
#' scores 1, ..., k following the order of their levels, so the levels must be
#' in the intended order. Unordered categorical variables cannot be used as
#' network variables.
#'
#' Binary network variables can be treated as ordered, which yields tetrachoric
#' and biserial correlations. As partial correlation networks assume (latent)
#' continuous variables, a warning is issued for binary network variables; for
#' purely binary data, an Ising model is usually more appropriate.
#'
#' Setting `ordered = TRUE` for any network variable requires the missing
#' handling method to be either `"stacked-mi"` or `"listwise"`, unless the
#' network variables are fully observed. Whether a network variable is treated
#' as ordered depends only on `ordered`, not on its class: an ordered factor
#' with `ordered = FALSE` enters the correlation matrix as continuous scores.
#'
#' Means are computed whenever all network variables are treated as continuous.
#' If any network variable is treated as ordered, `means` is returned as `NULL`.
#'
#' Auxiliary variables are only used to improve the handling of missing values
#' in the network variables and do not enter the returned correlation matrix.
#' They are not standardized. With `missing_handling = "two-step-em"`, they are
#' included in the joint EM estimation and must therefore be numeric or ordered
#' factors; both are treated as continuous. With `missing_handling = "stacked-mi"`,
#' they are used as predictors in the imputation models and may additionally be
#' unordered categorical variables stored as factors or character vectors (converted to factors).
#'
#' @returns A list containing:
#' \describe{
#'   \item{mat}{Estimated correlation matrix of the network variables.}
#'   \item{means}{Vector of estimated means of the network variables. If any
#'                network variable is treated as ordered, `means` is returned
#'                as `NULL`.}
#'   \item{cor_method}{A matrix indicating the correlation method used for each
#'                     pair of network variables.}
#'   \item{imputed_data}{Imputed datasets used for estimation. For
#'                       `missing_handling = "stacked-mi"`, the imputed data
#'                       are returned as a `mids` object from the \pkg{mice}
#'                       package (including the auxiliary variables); otherwise
#'                       `NULL`. The imputation methods actually used are
#'                       stored in `imputed_data$method`.}
#'   \item{args}{List of settings used in the correlation estimation.}
#' }
#' @export
#'
#' @examples
#' # Estimate correlation matrix from full data set
#' result <- cor_calc(data = mantar_dummy_full_cont,
#'                    ordered = FALSE)
#'
#' # View estimated correlation matrix and methods used
#' result$mat
#' result$cor_method
#'
#' # Estimate correlation matrix for data set with missings
#' result_mis <- cor_calc(data = mantar_dummy_mis_cont,
#'                       ordered = FALSE,
#'                       missing_handling = "two-step-em")
#'
#' # View estimated correlation matrix and methods used
#' result_mis$mat
#' result_mis$cor_method
cor_calc <- function(data, ordered = FALSE,
                     network_vars = NULL, auxiliary_vars = NULL,
                     missing_handling = "two-step-em",
                     nimp = 20, imp_method = NULL,
                     maxit = 10, ...) {

  # Capture additional arguments
  dots <- list(...)

  # Set variable names if not provided
  if (!is.null(data) && is.null(colnames(data))) {
    colnames(data) <- paste0("V", seq_len(ncol(data)))
  }
  # Transform the indicators for the used variables to character and perform some checks
  checker(network_vars = network_vars, auxiliary_vars = auxiliary_vars, called_from = "before_network_vars_check")
  network_vars <- chosen_vars_handling(vars = network_vars, data = data, type = "network")
  auxiliary_vars <- chosen_vars_handling(vars = auxiliary_vars, data = data, type = "auxiliary")
  checker(network_vars = network_vars, auxiliary_vars = auxiliary_vars, data = data, called_from = "after_network_vars_check")

  # reduce data set to variables used for network estimation and auxiliary variables
  # (network variables first, followed by auxiliary variables)
  data <- data[, c(network_vars, auxiliary_vars), drop = FALSE]

  # Match and validate missing handling method
  missing_handling <- match.arg(tolower(missing_handling),
                                choices = c("two-step-em", "stacked-mi", "pairwise", "listwise"))
  # Check data input - setting mat = NULL important to perform the correct check
  checker(data = data, mat = NULL)
  # Check usage of auxiliary variables and their type
  checker(
    data = data,
    network_vars = network_vars,
    auxiliary_vars = auxiliary_vars,
    missing_handling = missing_handling
  )

  # Check ordered input: logical, without NAs and one entry per network variable
  if (any(is.na(ordered))) {
    stop("`ordered` must not contain NA.")
  }
  if (any(!is.logical(ordered))){
    stop("All entries of `ordered` must be logical (TRUE/FALSE).")
  }
  p_net <- length(network_vars)
  p_aux <- length(auxiliary_vars)
  if (length(ordered) == 1L) {
    # Recycle to all network variables
    ordered <- rep(ordered, p_net)
  } else if (p_aux > 0L && length(ordered) == p_net + p_aux) {
    # Transition rule: earlier versions expected one entry per network and auxiliary variable
    warning(
      "`ordered` now only refers to the variables in `network_vars`. ",
      "The entries for `auxiliary_vars` are ignored, as auxiliary variables ",
      "do not enter the correlation matrix. Use `imp_method` to choose how ",
      "they are imputed.",
      call. = FALSE
    )
    ordered <- ordered[seq_len(p_net)]
  } else if (length(ordered) != p_net) {
    stop(
      "If `ordered` has length > 1, its length must equal the number of ",
      "network variables (p = ", p_net, ").",
      call. = FALSE
    )
  }
  names(ordered) <- network_vars
  any_ord <- any(ordered)

  # Ordinal network correlations require lavaan
  if (any_ord && !requireNamespace("lavaan", quietly = TRUE)) {
    stop(
      "Package 'lavaan' must be installed to compute polychoric ",
      "or polyserial correlations.",
      call. = FALSE
    )
  }


  data <- as.data.frame(data)
  # Convert ordered factors to integer scores 1, ..., k (following the order of the levels)
  ord_vars <- vapply(data, is.ordered, logical(1))
  data[ord_vars] <- lapply(data[ord_vars], as.integer)
  # Store character variables (only possible for auxiliary variables) as factors,
  # otherwise mice treats them as constant and silently drops them as predictors
  chr_vars <- vapply(data, is.character, logical(1))
  data[chr_vars] <- lapply(data[chr_vars], factor)

  # Warn about binary network variables
  binary_vars <- network_vars[vapply(data[network_vars], function(x) {
    length(unique(x[!is.na(x)])) == 2L
  }, logical(1))]
  if (length(binary_vars) > 0L) {
    msg <- paste0(
      "Binary network variable(s): ", paste(binary_vars, collapse = ", "), ". ",
      "The network estimation assumes (latent) continuous variables. For binary ",
      "data, consider an Ising model instead."
    )
    if (any(!ordered[binary_vars])) {
      msg <- paste0(
        msg, " Binary variables with ordered = FALSE enter the correlation ",
        "matrix with Pearson (phi) correlations; consider ordered = TRUE."
      )
    }
    warning(msg, call. = FALSE)
  }

  # Standardize the network variables when all of them are continuous;
  # auxiliary variables are kept on their original scale
  if (!any_ord) {
    data[network_vars] <- lapply(data[network_vars], function(x) {
      as.numeric(scale(x))
    })
  }

  # Prepare means vector and imputed data output
  means <- NULL
  imputed_data <- NULL

  # Handle missing data
  if (anyNA(data)){
    if (missing_handling == "two-step-em"){

      # EM algorithm cannot treat variables as ordered
      if (any_ord) {
        stop(
          "'two-step-em' treats all variables as continuous and cannot handle ",
          "ordered variables. Either set ordered = FALSE or choose a different ",
          "missing data method.",
          call. = FALSE
        )
      }

      # Check if lavaan is installed for the EM algorithm
      if (!requireNamespace("lavaan", quietly = TRUE)) {
        stop(
          "Package \"lavaan\" must be installed to use the EM algorithm to handle missingness.",
          call. = FALSE
        )
      }

      # Estimate correlation matrix with EM algorithm
      lavobject <- suppressWarnings(try(lavaan::lavCor(data,
                                                       missing = "ml", se = "none", meanstructure = TRUE,
                                                       estimator = "ML", output = "fit")))
      # If lavcor fails - stop the estimation
      if (inherits(lavobject, "try-error")) stop("lavaan::lavCor failed. Check your data.")
      # Extract estimated correlation matrix and means
      mat <- stats::cov2cor(lavaan::inspect(lavobject, "cov.ov"))
      means <- lavaan::inspect(lavobject, "mean.ov")

      # Set nimp to NULL as no multiple imputation was used
      nimp <- NULL

    } else if (missing_handling == "stacked-mi"){

      # Check if mice is installed for multiple imputation
      if (!requireNamespace("mice", quietly = TRUE)) {
        stop(
          "Package \"mice\" must be installed to use this function.",
          call. = FALSE
        )
      }

      # check multiple imputation to have a valid entry
      checker(nimp = nimp)

      # Keep column names safe during imputation
      original_names <- colnames(data)
      colnames(data) <- make.names(original_names)

      # Imputation methods: mice defaults based on the variable class (pmm for
      # numeric incl. converted ordered factors, logreg for binary and polyreg
      # for unordered factors) or the user input
      if (is.null(imp_method)) {
        imp_method <- mice::make.method(data)
      } else if (length(imp_method) != ncol(data)) {
        stop(
          "'imp_method' must be NULL or contain one entry per variable in ",
          "'network_vars' and 'auxiliary_vars' (p = ", ncol(data), ").",
          call. = FALSE
        )
      } else {
        imp_method <- unname(imp_method)
      }
      # Numeric variables imputed with "polr" are temporarily stored as ordered
      # factors; their observed values are kept to convert them back afterwards
      polr_vars <- which(imp_method == "polr" & vapply(data, is.numeric, logical(1)))
      polr_values <- lapply(data[polr_vars], function(x) sort(unique(x[!is.na(x)])))
      data[polr_vars] <- Map(function(x, v) factor(x, levels = v, ordered = TRUE),
                             data[polr_vars], polr_values)

      # Perform multiple imputation and stack the data
      mice_args <- list(data = data, m = nimp, maxit = maxit,
                        method = imp_method, printFlag = FALSE)
      if (length(dots) > 0L) mice_args[names(dots)] <- dots
      imputed_data <- suppressMessages(do.call(mice::mice, mice_args))
      stacked_data <- mice::complete(imputed_data, c(1:nimp))
      #  convert variables imputed with "polr" back to their original numeric values
      stacked_data[polr_vars] <- Map(function(x, v) v[as.integer(x)],
                                     stacked_data[polr_vars], polr_values)
      colnames(stacked_data) <- original_names  # restore original column names
      colnames(data) <- original_names
      stacked_data <- stacked_data[, network_vars, drop = FALSE] # reduce only to the network variables

      if (any_ord) {
        mat <- suppressWarnings(
          try(lavaan::lavCor(stacked_data, ordered = network_vars[ordered],
                             se = "none", output = "cor"))
        )
        if (inherits(mat, "try-error"))
          stop("lavaan::lavCor failed. Check your data.")
      } else {
        mat <- stats::cor(stacked_data)
        means <- colMeans(stacked_data)
      }

    } else if (missing_handling == "pairwise"){

      # auxiliary variables are not allowed here, so data only contains network variables
      if (any_ord) {
        stop(
          "'pairwise' computes Pearson correlations and cannot handle ordered ",
          "variables. Either set ordered = FALSE or choose a different missing ",
          "data method.",
          call. = FALSE
        )
      }
      mat <- stats::cov2cor(stats::cov(data, use = "pairwise.complete.obs"))
      means <- colMeans(data, na.rm = TRUE)
      nimp <- imp_method <- maxit <- NULL

    } else if (missing_handling == "listwise") {

      # auxiliary variables are not allowed here, so data only contains network variables
      data_complete <- stats::na.omit(data)

      if (any_ord) {
        mat <- suppressWarnings(
          try(lavaan::lavCor(data_complete, ordered = network_vars[ordered],
                             missing = "listwise", se = "none", output = "cor"))
        )
        if (inherits(mat, "try-error"))
          stop("lavaan::lavCor failed. Check your data.")
      } else {
        mat <- stats::cov2cor(stats::cov(data_complete, use = "complete.obs"))
        means <- colMeans(data_complete)
      }
      nimp <- imp_method <- maxit <- NULL
    }
  } else {

    # no missing values: auxiliary variables are not allowed, so data only contains network variables
    if (any_ord){
      mat <- suppressWarnings(try(lavaan::lavCor(data, ordered = network_vars[ordered],
                                                 se = "none", output = "cor")))
      if (inherits(mat, "try-error")) stop("lavaan::lavCor failed. Check your data.")
    } else {
      mat <- stats::cov2cor(stats::cov(data))
      means <- colMeans(data)
    }
    missing_handling <- NULL
    nimp <- imp_method <- maxit <- NULL
  }

  cor_method <- matrix("", p_net, p_net, dimnames = list(network_vars, network_vars))
  diag(cor_method) <- "-"

  cor_method[lower.tri(cor_method)] <- ifelse(
    outer(ordered, ordered, "&")[lower.tri(cor_method)], "polychoric",
    ifelse(
      outer(ordered, ordered, "!=")[lower.tri(cor_method)], "polyserial",
      "pearson"
    )
  )

  # reduce matrix and means to the network variables
  mat <- mat[network_vars, network_vars, drop = FALSE]
  if (!is.null(means)) means <- means[network_vars]


  return(list(
    mat = mat,
    means = means,
    cor_method = cor_method,
    imputed_data = imputed_data,
    args = list(
      missing_handling = missing_handling,
      nimp = nimp)
  ))
}


transform_inverse <- function(inverse_cov){
  rows <- nrow(inverse_cov)
  cols <- ncol(inverse_cov)
  cov <- solve(inverse_cov)
  inverse_cor <- matrix(NA, rows, cols)

  # it is important to do it this way, cause using solve, cov2cor and solve again introduces small errors
  for (i in 1:rows){
    for (j in 1:cols){
      inverse_cor[i,j] <- inverse_cov[i,j] * sqrt(diag(cov)[i]) * sqrt(diag(cov)[j])
    }
  }
  return(inverse_cor)
}



#' @title Handle Chosen Variables for Network and Auxiliary Variables
#'
#' @param vars Vector of the variables to be used.
#' @param data Data frame containing the variables.
#' @param mat Correlation matrix containing the variables.
#' @param type Character string indicating the type of variables being handled.
#'
#' @returns Numeric indices of the selected variables.
#' @noRd
chosen_vars_handling <- function(vars, data = NULL, mat = NULL,
                                 type = c("network", "auxiliary")) {

  type <- match.arg(type)
  arg_name <- paste0(type, "_vars")

  if (!is.null(mat) && type == "auxiliary" && !is.null(vars)) {
    stop(
      "'auxiliary_vars' can only be used when raw data are provided via 'data', not when a correlation matrix is supplied via 'mat'."
    )
  }

  x <- if (!is.null(data)) data else mat
  p <- ncol(x)
  var_names <- colnames(x)


  if (is.null(vars)) {
    if (type == "network") return(var_names)
    if (type == "auxiliary") return(NULL)
  }

  if (is.numeric(vars)) {
    if (any(vars < 1 | vars > p)) {
      stop(
        "'", arg_name, "' contains indices outside the valid range 1 to ",
        p, "."
      )
    }
    return(var_names[vars])
  }

  if (is.character(vars)) {
    unknown <- setdiff(vars, var_names)

    if (length(unknown) > 0) {
      stop(
        "Unknown variable name(s) in '", arg_name, "': ",
        paste(unknown, collapse = ", ")
      )
    }

    return(vars)
  }

  stop(
    "'", arg_name, "' must be numeric indices, character variable names, or NULL."
  )
}
