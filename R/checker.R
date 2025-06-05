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
      }
  }

  # Check length of ns
  if (all(c("data", "mat") %in% names(args))) {
    if (length(ns) != ncol(mat)) {
      stop("Length of 'ns' must match the number of columns in 'mat'.")
    }
  }
}
