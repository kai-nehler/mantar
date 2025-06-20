checker <- function(...) {
  args <- list(...)
  data <- args$data
  mat <- args$mat
  ns <- args$ns

  # Check: Which input is provided?
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
          stop("All variables in 'data' must be numeric.")
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

  # Check length of ns
  if ("ns" %in% names(args)) {
    if ("mat" %in% names(args)){
      if (length(ns) != ncol(mat)) {
        stop("Length of 'ns' must match the number of columns in 'mat'.")
      }
    }
    if ("data" %in% names(args)){
      if (length(ns) != ncol(data)) {
        stop("Length of 'ns' must match the number of columns in 'data'.")
      }
    }
  }
}
