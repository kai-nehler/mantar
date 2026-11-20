#' @title
#' Calculate network density
#'
#' @description
#' Computes the proportion of possible edges whose absolute weight exceeds
#' a specified threshold. Network sparsity is one minus network density.
#'
#' @param x A square numeric matrix of edge weights.
#' @param directed Logical. If `FALSE`, each node pair is counted once using
#'   the upper triangle of `x`. If `TRUE`, both directions are counted.
#' @param include_diag Logical. If `TRUE`, diagonal entries are included as
#'   possible self-loops. This can be used to include autoregressive effects
#'   in a temporal network.
#' @param tol A finite, non-negative threshold. Edge weights with an absolute
#'   value less than or equal to `tol` count as zero. Defaults to `0`.
#'
#' @returns A single numeric value between 0 and 1. Returns `NA_real_` if
#'   there are no possible edges or if a selected edge weight is missing.
#'
#' @examples
#' # Undirected network: two of three possible edges are present
#' pcor <- matrix(c(
#'   1.0, 0.4, 0.0,
#'   0.4, 1.0, 0.2,
#'   0.0, 0.2, 1.0
#' ), nrow = 3, byrow = TRUE)
#' network_density(pcor)
#'
#' # Temporal network: three of four possible effects are present
#' beta <- matrix(c(
#'   0.3, -0.2,
#'   0.0,  0.1
#' ), nrow = 2, byrow = TRUE)
#' network_density(beta, directed = TRUE, include_diag = TRUE)
#'
#' @export
network_density <- function(x, directed = FALSE, include_diag = FALSE, tol = 0) {
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) != ncol(x)) {
    stop("x must be a square numeric matrix.", call. = FALSE)
  }
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) {
    stop("directed must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(include_diag) || length(include_diag) != 1L ||
      is.na(include_diag)) {
    stop("include_diag must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0) {
    stop("tol must be one finite, non-negative number.", call. = FALSE)
  }

  if (directed) {
    weights <- if (include_diag) x else x[row(x) != col(x)]
  } else {
    weights <- x[upper.tri(x, diag = include_diag)]
  }

  if (length(weights) == 0L || anyNA(weights)) {
    return(NA_real_)
  }

  mean(abs(weights) > tol)
}
