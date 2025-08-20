cor_calc <- function(data, cor_method = "adapted",
                     missing_handling = "two-step-em",
                     nimp = 20, imp_method = "pmm") {

  cor_method <- match.arg(tolower(cor_method), choices = c("adapted", "pearson", "polychoric"))
  missing_handling <- match.arg(tolower(missing_handling), choices = c("two-step-em", "stacked-mi", "pairwise", "listwise"))

  if (cor_method == "adapted") {
    cor_method <- "pearson"  # Default to Pearson if adapted is selected
  }

  # Check if package is installed if polychoric method is selected
  if (cor_method == "polychoric" && !requireNamespace("psych", quietly = TRUE)) {
    stop(
      "Package \"psych\" must be installed to compute .",
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
      mat <- try(stats::cov2cor(lavaan::inspect(lavobject, "cov.ov")))
      if (cor_method == "polychoric"){
        warning("Using a specific two-step EM algorithm based on polychoric correlations is not implemented at the moment.
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
      if (cor_method == "polychoric"){
        mat <- psych::polychoric(stacked_data)$rho
      }else {
        mat <- stats::cor(stacked_data)
      }

    } else if (missing_handling == "pairwise"){
      mat <- stats::cov2cor(stats::cov(data, use = "pairwise.complete.obs"))
      if (cor_method == "polychoric"){
        warning("Using pairwise deletion based on polychoric correlations is not implemented at the moment.
                Will proceed with pairwise deletion based on pearson correlations.")
      }
      nimp <- NULL
    } else if (missing_handling == "listwise"){
      if (cor_method == "polychoric"){
        mat <- psych::polychoric(data)$rho
      } else {
        mat <- stats::cov2cor(stats::cov(data, use = "complete.obs"))
      }
      nimp <- NULL
    }
  } else {
    if (cor_method == "polychoric"){
      mat <- psych::polychoric(data)$rho
    } else {
      mat <- stats::cov2cor(stats::cov(data))
    }
    missing_handling <- NULL
    nimp <- NULL
  }


  return(list(mat = mat, missing_handling = missing_handling, nimp = nimp, cor_method = cor_method))
}
