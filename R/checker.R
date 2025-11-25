#' @title Checker function to validate inputs
#'
#' @param ... Any arguments to be checked.
#'
#' @returns Either stops execution with an error message if checks fail, or
#'          returns NULL invisibly if all checks pass.
#' @noRd
checker <- function(...) {
  # Capture arguments
  args <- list(...)
  # Extract specific arguments (with default NULL if not provided)
  data <- args$data
  mat <- args$mat
  ns <- args$ns
  n <- args$n

  # This section deals with checking data and mat inputs
  if (all(c("data", "mat") %in% names(args))) {
    if (is.null(data) && is.null(mat)) {
      stop("Either 'data' or 'mat' must be provided.")
      } else if (!is.null(data) && !is.null(mat)) {
    message("Both 'data' and 'mat' provided. 'data' will be used.")
        if (!is.data.frame(data) & !is.matrix(data)) {
          stop("'data' must be a data frame or matrix.")
        }
        if (!all(sapply(data, is.numeric))){
          stop("All variables in 'data' must be numeric.")
        }
      } else if (!is.null(data)) {
        if (!all(sapply(data, is.numeric))){
          stop("All variables in 'data' must be numeric.
               If you have an ordered categorical variable stored as a factor, please convert it to numeric first.
               This requirement helps prevent accidentally including unordered categorical variables in the calculation.")
        }
      } else if (!is.null(mat)) {
        if (!is.matrix(mat)) {
          stop("'mat' must be a matrix.")
        } else if(!isSymmetric(mat)) {
          stop("'mat' must be a symmetric matrix.")
        }
        if (!is.numeric(mat)) {
          stop("All entries in 'mat' must be numeric")
        }
        if (any(is.na(mat))) {
          stop("'mat' must not contain missing values.")
        }
      }
  }

  # This section deals with checking ns for multiple regression model
  if ("ns" %in% names(args)) {
    if ("mat" %in% names(args)){
      if (length(ns) != ncol(mat) & length(ns) != 1) {
        stop("'ns' must be either a single value or a vector with one entry per column in 'mat'.")
      }
    }
    if ("data" %in% names(args)){
      if (length(ns) != ncol(data) & length(ns) != 1) {
        stop("'ns' must be either a single value or a vector with one entry per column in 'data'.")
      }
    }
  }

  # This section deals with checking n for a single regression model
  if ("n" %in% names(args)) {
    if (length(n) != 1){
      stop("Length of 'n' must be 1.")
    } else if (!is.numeric(n) || n <= 0) {
      stop("'n' must be a positive numeric value.")
    }
  }
}

