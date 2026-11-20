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

  # This section deals with checking data and mat inputs
  if (all(c("data", "mat") %in% names(args))) {
    if (is.null(args$data) & is.null(args$mat)) {
      stop("Either 'data' or 'mat' must be provided.")
    }

    if (!is.null(args$data)) {
      if (!is.data.frame(args$data) & !is.matrix(args$data)) {
        stop("'data' must be a data frame or matrix.")
      }
      if (!all(sapply(args$data, is.numeric))) {
        stop("All variables in 'data' must be numeric.
           If you have an ordered categorical variable stored as a factor,
           please convert it to numeric first. This requirement helps prevent
           accidentally including unordered categorical variables in the calculation.")
      }
    }

    if (!is.null(args$mat)) {
      if (!is.matrix(args$mat)) {
        stop("'mat' must be a matrix.")
      }
      if (!isSymmetric(args$mat)) {
        stop("'mat' must be a symmetric matrix.")
      }
      if (!is.numeric(args$mat)) {
        stop("All entries in 'mat' must be numeric.")
      }
      if (any(is.na(args$mat))) {
        stop("'mat' must not contain missing values.")
      }
      if (!is.null(rownames(args$mat)) & !is.null(colnames(args$mat))) {
        if (!identical(rownames(args$mat), colnames(args$mat))) {
          stop("Row and column names of 'mat' must be identical.")
        }
      }
    }

    if (!is.null(args$data) & !is.null(args$mat)) {
      if (ncol(args$data) != ncol(args$mat)) {
        stop("The number of columns used for network analysis in 'data' must match the number of columns in 'mat'.")
      }
      if (!is.null(colnames(args$data)) & !is.null(colnames(args$mat))) {
        if (!identical(colnames(args$data), colnames(args$mat))) {
          stop("Column names of 'data' and 'mat' must match and be in the same order.")
        }
      }
    }
  }

  if (!is.null(called_from)) {
    if (called_from == "neighborhood") {
      if ("ns" %in% names(args)) {
        if ("mat" %in% names(args)) {
          if (length(args$ns) != ncol(args$mat) & length(args$ns) != 1) {
            stop("'ns' must be either a single value or a vector with one entry per column used for network estimation in 'mat' (after optional selection via 'network_vars').")
          }
        }
        if ("data" %in% names(args)) {
          if (length(args$ns) != ncol(args$data) & length(args$ns) != 1) {
            stop("'ns' must be either a single value or a vector with one entry per column used for network estimation in 'data' (after optional selection via 'network_vars').")
          }
        }
      }
    }

    if (called_from == "regularization") {
      # This section deals with checking ns for multiple regression model
      if ("ns" %in% names(args)) {
        if ("mat" %in% names(args)) {
          if (!(length(args$ns) == 1 ||
                (is.matrix(args$ns) &&
                 nrow(args$ns) == ncol(args$mat) &&
                 ncol(args$ns) == ncol(args$mat)))) {
            stop("'ns' must be either a single value or a matrix with dimensions matching the matrix used for network estimation in 'mat' (after optional selection via 'network_vars').")
          }
        }
        if ("data" %in% names(args)) {
          if (!(length(args$ns) == 1 ||
                (is.matrix(args$ns) &&
                 nrow(args$ns) == ncol(args$data) &&
                 ncol(args$ns) == ncol(args$data)))) {
            stop("'ns' must be either a single value or a matrix with dimensions matching the matrix used for network estimation in 'data' (after optional selection via 'network_vars').")
          }
        }
      }
    }
  }

  # This section deals with checking n for a single regression model
  if ("n" %in% names(args)) {
    if (length(args$n) != 1) {
      stop("Length of 'n' must be 1.")
    } else if (!is.numeric(args$n) || args$n <= 0) {
      stop("'n' must be a positive numeric value.")
    }
  }

  # Check number of imputations
  if ("nimp" %in% names(args)) {
    if (!is.numeric(args$nimp) || length(args$nimp) != 1L ||
        !is.finite(args$nimp) || args$nimp < 1 || args$nimp %% 1 != 0) {
      stop("'nimp' must be a single positive integer.")
    }
  }

  if ("conf_level" %in% names(args)) {
    if (!is.numeric(args$conf_level) ||
        length(args$conf_level) != 1L ||
        is.na(args$conf_level) ||
        args$conf_level <= 0 ||
        args$conf_level >= 1) {
      stop("`conf.level` must be a single number between 0 and 1.")
    }
  }

  # This section deals with selected network and auxiliary variables
  if (!is.null(called_from)) {
    if (called_from == "before_network_vars_check") {
      if (all(c("network_vars", "auxiliary_vars") %in% names(args))) {
        if (is.null(args$network_vars) & !is.null(args$auxiliary_vars)) {
          stop("'auxiliary_vars' can only be used when 'network_vars' is specified. If 'network_vars' is NULL, all variables are used for network estimation, so no separate auxiliary variables can be defined.")
        }
      }
    }

    if (called_from == "after_network_vars_check") {
      # Check duplicate selected variables
      if (!is.null(args$network_vars)) {
        if (anyDuplicated(args$network_vars)) {
          stop("'network_vars' must not contain duplicate variables.")
        }
      }

      if (!is.null(args$auxiliary_vars)) {
        if (anyDuplicated(args$auxiliary_vars)) {
          stop("'auxiliary_vars' must not contain duplicate variables.")
        }
      }

      # Check overlap between network and auxiliary variables
      if (!is.null(args$network_vars) && !is.null(args$auxiliary_vars)) {
        if (length(intersect(args$network_vars, args$auxiliary_vars)) > 0) {
          stop("'network_vars' and 'auxiliary_vars' must not contain overlapping variables.")
        }
      }
    }
  }

  # Check lambda to be varying if the penalty is glasso
  if (all(c("penalty", "vary") %in% names(args))) {
    if (args$penalty == "glasso" & args$vary != "lambda") {
      stop("For 'glasso' penalty, 'vary' must be set to 'lambda' as this is the only penalty parameter. If you want to provide your own lambda values you can do this in the corresponding argument 'lambda' but you still have to set 'vary' to 'lambda'.")
    }
  }
}
