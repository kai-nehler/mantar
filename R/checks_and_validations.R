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

  # Check scalar counts
  count_args <- c("nimp", "nperson", "nnode", "max_shrink_attempts")
  for (name in intersect(count_args, names(args))) {
    x <- args[[name]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < 1 || x != floor(x)) {
      stop("'", name, "' must be a single positive integer.")
    }
  }

  if ("conf_level" %in% names(args)) {
    validate_proportion(args$conf_level, "conf_level", inclusive = FALSE)
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

  # mlVAR Simulation dimensions
  if ("ntime" %in% names(args)) {
    if (!is.numeric(args$ntime) ||
        !(length(args$ntime) %in% c(1L, args$nperson)) ||
        any(!is.finite(args$ntime)) || any(args$ntime < 1) ||
        any(args$ntime != floor(args$ntime))) {
      stop("ntime must be a positive integer or one positive integer per person.")
    }
  }


  # Check that each supplied flag is a single TRUE or FALSE
  for (name in intersect(c("only_nets", "shrink_temporal_fixed"), names(args))) {
    x <- args[[name]]
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
      stop(name, " must be TRUE or FALSE.")
    }
  }

  # mlVAR_sim check Temporal population network
  if ("beta_fixed" %in% names(args)) {
    validate_population_matrix(args$beta_fixed, "beta_fixed", args$nnode)

    if (is.null(args$beta_fixed)) {
      validate_proportion(args$temporal_sparsity, "temporal_sparsity")
      validate_proportion(args$temporal_prop_positive, "temporal_prop_positive")
    }
  }

  # Temporal random effects are generated even when beta_fixed is supplied.
  if ("temporal_beta_re_sd_range" %in% names(args)) {
    validate_sd_range(args$temporal_beta_re_sd_range, "temporal_beta_re_sd_range")
  }
  if ("temporal_re_correlation_sparsity" %in% names(args)) {
    validate_proportion(
      args$temporal_re_correlation_sparsity, "temporal_re_correlation_sparsity"
    )
  }

  # Contemporaneous network
  if ("theta_fixed" %in% names(args)) {
    validate_population_matrix(
      args$theta_fixed, "theta_fixed", args$nnode, positive_definite = TRUE
    )

    if (args$contemporaneous == "random_ggm" && !is.null(args$theta_fixed)) {
      stop("theta_fixed cannot be specified with contemporaneous = 'random_ggm': ",
           "the population innovation covariance is computed from the ",
           "individual covariance matrices.")
    }

    if (is.null(args$theta_fixed)) {
      validate_proportion(args$contemporaneous_sparsity, "contemporaneous_sparsity")

      if (!is.numeric(args$theta_var) ||
          length(args$theta_var) != args$nnode ||
          any(!is.finite(args$theta_var)) || any(args$theta_var <= 0)) {
        stop("theta_var must contain one positive variance per node.")
      }
    }

    if (args$contemporaneous == "orthogonal_re") {
      validate_sd_range(
        args$contemporaneous_pcor_re_sd_range, "contemporaneous_pcor_re_sd_range"
      )
    }
  }

  # Between-person network
  if ("between_cov" %in% names(args)) {
    validate_population_matrix(
      args$between_cov, "between_cov", args$nnode, positive_definite = TRUE
    )

    if (is.null(args$between_cov)) {
      validate_proportion(args$between_sparsity, "between_sparsity")
      validate_sd_range(args$random_intercept_sd_range, "random_intercept_sd_range")
    }
  }

  if ("fixed_mean_sd" %in% names(args)) {
    if (!is.numeric(args$fixed_mean_sd) || length(args$fixed_mean_sd) != 1L ||
        !is.finite(args$fixed_mean_sd) || args$fixed_mean_sd < 0) {
      stop("fixed_mean_sd must be one finite, non-negative number.")
    }
  }

  # mlVAR sim Shrinkage factors
  shrink_args <- c(
    "shrink_temporal_fixed_factor",
    "shrink_temporal_re_factor",
    "shrink_contemporaneous_re_factor"
  )
  for (name in intersect(shrink_args, names(args))) {
    validate_proportion(args[[name]], name, inclusive = FALSE)
  }

  invisible(NULL)
}




#' Check that a single numeric value lies in [0, 1] or (0, 1)
#' @noRd
validate_proportion <- function(x, name, inclusive = TRUE) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
      x < 0 || x > 1) {
    stop(name, " must be a single number between 0 and 1.")
  }

  if (!inclusive && (x == 0 || x == 1)) {
    stop(name, " must be strictly between 0 and 1.")
  }
}

#' Check a range of two finite, non-negative SD bounds
#' @noRd
validate_sd_range <- function(x, name) {
  if (length(x) != 2L || !is.numeric(x) || any(!is.finite(x)) ||
      any(x < 0) || x[1] > x[2]) {
    stop(name, " must be c(lower, upper), with 0 <= lower <= upper.")
  }
}

#' Check the dimensions and entries of an optional population matrix
#' @noRd
validate_population_matrix <- function(x, name, size, positive_definite = FALSE) {
  if (is.null(x)) return(invisible(NULL))

  if (!is.matrix(x) || !is.numeric(x) ||
      !identical(dim(x), c(as.integer(size), as.integer(size))) ||
      any(!is.finite(x))) {
    stop(name, " must be a finite numeric ", size, " x ", size, " matrix.")
  }

  if (positive_definite) {
    if (!isTRUE(all.equal(unname(x), unname(t(x)), tolerance = 1e-10))) {
      stop(name, " must be symmetric.")
    }

    eig <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
    if (min(eig) <= max(eig) * 1e-12) {
      stop(name, " must be positive definite.")
    }
  }

  invisible(NULL)
}
