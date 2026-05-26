#' @title Checker function to validate inputs
#'
#' @param ... Any arguments to be checked.
#'
#' @returns Either stops execution with an error message if checks fail, or
#'          returns NULL invisibly if all checks pass.
#' @noRd
checker <- function(called_from = NULL, ...) {
  # Capture arguments
  args <- list(...)
  # Extract specific arguments (with default NULL if not provided)
  data <- args$data
  mat <- args$mat
  ns <- args$ns
  n <- args$n
  penalty <- args$penalty
  vary <- args$vary
  network_vars <- args$network_vars
  auxiliary_vars <- args$auxiliary_vars

  # This section deals with checking data and mat inputs
  if (all(c("data", "mat") %in% names(args))) {
    if (is.null(data) & is.null(mat)) {
      stop("Either 'data' or 'mat' must be provided.")
    }
    if (!is.null(data)) {
      if (!is.data.frame(data) & !is.matrix(data)) {
        stop("'data' must be a data frame or matrix.")
      }
      if (!all(sapply(data, is.numeric))) {
        stop("All variables in 'data' must be numeric.
           If you have an ordered categorical variable stored as a factor,
           please convert it to numeric first. This requirement helps prevent
           accidentally including unordered categorical variables in the calculation.")
      }
    }
    if (!is.null(mat)) {
      if (!is.matrix(mat)) {
        stop("'mat' must be a matrix.")
      }
      if (!isSymmetric(mat)) {
        stop("'mat' must be a symmetric matrix.")
      }
      if (!is.numeric(mat)) {
        stop("All entries in 'mat' must be numeric.")
      }
      if (any(is.na(mat))) {
        stop("'mat' must not contain missing values.")
      }
      if (!is.null(rownames(mat)) & !is.null(colnames(mat))) {
        if (!identical(rownames(mat), colnames(mat))) {
          stop("Row and column names of 'mat' must be identical.")
        }
      }
    }
    if (!is.null(data) & !is.null(mat)) {
      if (ncol(data) != ncol(mat)) {
        stop("The number of columns used for network analysis in 'data' must match the number of columns in 'mat'.")
      }
      if (!is.null(colnames(data)) & !is.null(colnames(mat))) {
        if (!identical(colnames(data), colnames(mat))) {
          stop("Column names of 'data' and 'mat' must match and be in the same order.")
        }
      }
    }
  }

  if (!is.null(called_from)){
    if(called_from == "neighborhood"){
    if ("ns" %in% names(args)){
      if ("mat" %in% names(args)){
        if (length(ns) != ncol(mat) & length(ns) != 1){
          stop("'ns' must be either a single value or a vector with one entry per column used for network estimation in 'mat' (after optional selection via 'network_vars').")
        }
      }
      if ("data" %in% names(args)){
        if (length(ns) != ncol(data) & length(ns) != 1){
          stop("'ns' must be either a single value or a vector with one entry per column used for network estimation in 'data' (after optional selection via 'network_vars').")
        }
      }
    }
  }

  if (called_from == "regularization"){
  # This section deals with checking ns for multiple regression model
  if ("ns" %in% names(args)) {
    if ("mat" %in% names(args)) {
      if (!(length(ns) == 1 || (is.matrix(ns) && nrow(ns) == ncol(mat) && ncol(ns) == ncol(mat)))) {
        stop("'ns' must be either a single value or a matrix with dimensions matching the matrix used for network estimation in 'mat' (after optional selection via 'network_vars').")
      }
    }
    if ("data" %in% names(args)) {
      if (!(length(ns) == 1 || (is.matrix(ns) && nrow(ns) == ncol(data) && ncol(ns) == ncol(data)))) {
        stop("'ns' must be either a single value or a matrix with dimensions matching the matrix used for network estimation in 'data' (after optional selection via 'network_vars').")
      }
    }
  }
  }}

  # This section deals with checking n for a single regression model
  if ("n" %in% names(args)) {
    if (length(n) != 1){
      stop("Length of 'n' must be 1.")
    } else if (!is.numeric(n) || n <= 0) {
      stop("'n' must be a positive numeric value.")
    }
  }

  # this section deals with selected network and auxiliary variables

  if (!is.null(called_from)){
  if (called_from == "before_network_vars_check") {
    if (all(c("network_vars", "auxiliary_vars") %in% names(args))) {
      if (is.null(network_vars) & !is.null(auxiliary_vars)) {
      stop( "'auxiliary_vars' can only be used when 'network_vars' is specified. If 'network_vars' is NULL, all variables are used for network estimation, so no separate auxiliary variables can be defined.")
    }
    }
  }
  if (called_from == "after_network_vars_check") {
    # Check duplicate selected variables
    if (!is.null(network_vars)) {
      if (anyDuplicated(network_vars)) {
        stop("'network_vars' must not contain duplicate variables.")
      }
    }

    if (!is.null(auxiliary_vars)) {
      if (anyDuplicated(auxiliary_vars)) {
        stop("'auxiliary_vars' must not contain duplicate variables.")
      }
    }

    # Check overlap between network and auxiliary variables
    if (!is.null(network_vars) && !is.null(auxiliary_vars)) {
      if (length(intersect(network_vars, auxiliary_vars)) > 0) {
        stop(
          "'network_vars' and 'auxiliary_vars' must not contain overlapping variables."
        )
      }
    }
  }
  }



  # Check lambda to be varying if the penalty is glasso
  if (all(c("penalty", "vary") %in% names(args))) {
    if (penalty == "glasso" & vary != "lambda") {
      stop("For 'glasso' penalty, 'vary' must be set to 'lambda' as this is the only penalty parameter. If you want to provide your own lambda values you can do this in the corresponding argument 'lambda' but you still have to set 'vary' to 'lambda'.")
    }
  }

}

