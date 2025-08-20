cor_calc <- function(data, cor_method = "adapted",
                     missing_handling = "two-step-em",
                     nimp = 20, imp_method = "pmm", max_categories = 7) {

  cor_method <- match.arg(tolower(cor_method), choices = c("adapted", "pearson", "polychoric"))
  missing_handling <- match.arg(tolower(missing_handling), choices = c("two-step-em", "stacked-mi", "pairwise", "listwise"))

  if (cor_method == "adapted") {
    if ((sum(!is.na(data))/80) < (ncol(data) * (ncol(data) - 1) / 2)) {
      ordered <- rep(FALSE, ncol(data))
      text <- "Using 'pearson' correlation method due to small ratio of available information to estimated correlations."
    } else {
      ordered <- vapply(data, function(x) {
        length(unique(x[!is.na(x)])) <= max_categories
      }, logical(1))
      if(any(ordered)) {
        if (sum(ordered) == ncol(data)) {
          text <- "Using only 'polychoric' correlations because all variables are ordered."
        } else {
          text <- paste("Using mix of 'pearson', 'polychoric' and 'polyserial' correlations because some variables are ordered, but not all. Ordered Variable are ", names(data)[ordered], collapse = " ")
        }
      } else {
        text <- "Using the 'pearson' correlation method because the number of distinct values suggests continuous variables."
      }
    }
    message(text)
  } else if (cor_method == "polychoric"){
    ordered <- rep(TRUE, ncol(data))
  } else if (cor_method == "pearson") {
    ordered <- rep(FALSE, ncol(data))
  }

  # Check if package is installed if polychoric method is selected
  if (any(ordered) && !requireNamespace("lavaan", quietly = TRUE)) {
    stop(
      "Package \"lavaan\" must be installed to compute polychoric or polyserial correlations.",
      call. = FALSE
    )
  }

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
      mat <- stats::cov2cor(lavaan::inspect(lavobject, "cov.ov"))
      if (cor_method == "polychoric" | any(ordered)){
        warning("Using a specific two-step EM algorithm based on polychoric or polyserial correlations is not implemented at the moment.
                  Will proceed with the two-step EM approach based on pearson correlations.")
      }
      nimp <- NULL

    } else if (missing_handling == "stacked-mi"){

      if (!requireNamespace("mice", quietly = TRUE)) {
        stop(
          "Package \"mice\" must be installed to use this function.",
          call. = FALSE
        )
      }
      original_names <- colnames(data)
      colnames(data) <- make.names(original_names)

      imputed_data <- suppressMessages(mice::mice(data, m = nimp, maxit = 10, method = imp_method,
                                                  ridge = 0, donors = 5, ls.method = "qr"))
      stacked_data <- mice::complete(imputed_data, c(1:nimp))
      colnames(stacked_data) <- original_names  # restore original column names
      if (any(ordered)){
        mat <- suppressWarnings(try(lavaan::lavCor(data, ordered = names(data)[ordered],
                                                   se = "none", output = "cor")))
        if (inherits(mat, "try-error")) stop("lavaan::lavCor failed. Check your data.")
      }else {
        mat <- stats::cor(stacked_data)
      }

    } else if (missing_handling == "pairwise"){
      if (any(ordered)){
        warning("Using pairwise deletion based on polychoric or polyserial correlations is not implemented at the moment.
                Will proceed with pairwise deletion based on pearson correlations.")
      }
      mat <- stats::cov2cor(stats::cov(data, use = "pairwise.complete.obs"))
      nimp <- NULL
    } else if (missing_handling == "listwise"){
      if (any(ordered)){
        mat <- suppressWarnings(try(lavaan::lavCor(data, ordered = names(data)[ordered],
          missing = "listwise", se = "none", output = "cor")))
        if (inherits(mat, "try-error")) stop("lavaan::lavCor failed. Check your data.")
      } else {
        mat <- stats::cov2cor(stats::cov(data, use = "complete.obs"))
      }
      nimp <- NULL
    }
  } else {
    if (any(ordered)){
      mat <- suppressWarnings(try(lavaan::lavCor(data, ordered = names(data)[ordered],
                                                 se = "none", output = "cor")))
      if (inherits(mat, "try-error")) stop("lavaan::lavCor failed. Check your data.")
    } else {
      mat <- stats::cov2cor(stats::cov(data))
    }
    missing_handling <- NULL
    nimp <- NULL
  }

  cor_method <- matrix(NA_character_, ncol(data), ncol(data), dimnames = list(names(data), names(data)))

  cor_method[lower.tri(cor_method)] <- ifelse(
    outer(ordered, ordered, "&")[lower.tri(cor_method)], "polychoric",
    ifelse(
      outer(ordered, ordered, "!=")[lower.tri(cor_method)], "polyserial",
      "pearson"
    )
  )


  return(list(mat = mat, missing_handling = missing_handling, nimp = nimp, cor_method = cor_method))
}
