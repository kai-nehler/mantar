#' @export
print.mantar_network <- function(x, ...) {
  print(x$pcor)

  invisible(x)
}


#' @export
summary.mantar_neighborhood <- function(object, ...) {
  stopifnot(inherits(object, "mantar_neighborhood"))

  # Compute network density
  pcor <- object$pcor
  n_nodes <- ncol(pcor)
  n_edges <- sum(pcor[upper.tri(pcor)] != 0)
  density <- n_edges / choose(n_nodes, 2)

  result <- list(
    density = density,
    args = object$args,
    ns = object$ns,
    varnames = colnames(pcor)
  )

  class(result) <- "summary.mantar_neighborhood"
  return(result)
}

#' @export
print.summary.mantar_neighborhood <- function(x, ...) {
  cat(sprintf("The density of the estimated network is %.3f\n\n", x$density))

  # Extract arguments
  ic_type <- to.upper(x$args$ic_type)
  pcor_rule <- x$args$pcor_merge_rule
  missing_handling <- x$args$missing_handling
  nimp <- x$args$nimp

  if (is.null(missing_handling)) {
    cat(sprintf("Network was estimated using neighborhood selection with the information criterion: %s\n", ic_type))
    cat(sprintf("and the '%s' rule for the inclusion of edges based on a full data set.\n\n", pcor_rule))
  } else {
    cat("Network was estimated using neighborhood selection on data with missing values.\n")
    cat(sprintf("Missing data were handled using '%s'.\n", missing_handling))
    if (!is.null(nimp)) {
      cat(sprintf("Stacked multiple imputation was performed with %d imputations.\n", nimp))
    }
    cat(sprintf("The information criterion was %s and the '%s' rule was used for edge inclusion.\n\n", ic_type, pcor_rule))
  }

  cat("The sample sizes used for the nodewise regressions were as follows:\n")
  ns_named <- stats::setNames(x$ns, x$varnames)
  print(ns_named)


  invisible(x)
}

#' @export
summary.mantar_regularization <- function(object, ...) {
  stopifnot(inherits(object, "mantar_regularization"))

  # Compute network density
  pcor <- object$pcor
  n_nodes <- ncol(pcor)
  n_edges <- sum(pcor[upper.tri(pcor)] != 0)
  density <- n_edges / choose(n_nodes, 2)

  result <- list(
    density = density,
    args = object$args,
    n = object$n,
    varnames = colnames(pcor)
  )

  class(result) <- "summary.mantar_regularization"
  return(result)
}

#' @export
print.summary.mantar_regularization <- function(x, ...) {
  cat(sprintf("The density of the estimated network is %.3f\n\n", x$density))

  # Extract arguments
  penalty <- x$args$penalty
  missing_handling <- x$args$missing_handling
  nimp <- x$args$nimp
  likelihood <- x$args$likelihood

  cat(sprintf("Network was estimated using regularization with the %s penalty.", penalty))
  if (likelihood == "obs_based") {
    cat(" The observed-data log-likelihood was used in the model selection.\n")
    cat(sprintf("The effective sample size used for the information criteria
                computation was %s. \n", x$n))
    } else if (likelihood == "mat_based") {
      cat("Log-Likelihood calculation was based on the sample correlation matrix. \n")
      cat(sprintf("The effective sample size used for the information criteria
                  and log-liklihood computation was %s. \n", x$n))
    }

    if (!is.null(missing_handling)) {
      cat(sprintf("Missing data were handled using '%s'.\n", missing_handling))
      if (!is.null(nimp)) {
        cat(sprintf("Stacked multiple imputation was performed with %d imputations.\n", nimp))
    }}

  invisible(x)
}


#' @export
plot.mantar_network <- function(x, layout = "spring", ...) {
  if (!requireNamespace("qgraph", quietly = TRUE)) {
    stop(
      "The 'qgraph' package must be installed to plot a mantar_network object.\n",
      "You can install it with: install.packages('qgraph')"
    )
  }

  # Call qgraph::qgraph using the partial correlation matrix and user arguments
  qgraph::qgraph(x$pcor, layout = layout, ...)
}

#' @export
print.mantar_regression <- function(x, ...) {
  print(x$regression)

  invisible(x)
}

#' @export
summary.mantar_regression <- function(object, ...) {
  stopifnot(inherits(object, "mantar_regression"))

  result <- list(
    regression = object$regression,
    R2 = object$R2,
    n = object$n,
    args = object$args
  )

  class(result) <- "summary.mantar_regression"
  return(result)
}

#' @export
print.summary.mantar_regression <- function(x, ...) {
  cat("Regression results:\n")
  print(x$regression)

  cat(sprintf("\nR-squared: %.3f\n", x$R2))
  cat(sprintf("The sample size used in the predictor selection: %d\n", x$n))

  invisible(x)
}
