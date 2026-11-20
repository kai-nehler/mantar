mlvar_pool <- function(fits, conf_level = 0.95, dfcom = Inf) {
  checker(conf_level = conf_level)

  if (fits$fits[[1]]$input$estimator != "lmer") {
    stop("mlvar_pool() supports only estimator = 'lmer'.", call. = FALSE)
  }
  if (!identical(fits$fits[[1]]$input$AR, FALSE)) {
    stop("mlvar_pool() supports only AR = FALSE.")
  }

  mu <- lapply(fits$fits, function(x) x$results$mu)

  # Gemeinsamen Intercept für jeden Knoten nach Rubin poolen
  mu_table <- pool_rubin(
    estimates = lapply(mu, function(x) x$mean),
    ses = lapply(mu, function(x) x$SE),
    dfcom = dfcom
  )

  mu_summary <- summarize_pooled(mu_table, conf_level = conf_level)

  mu_pooled <- update_mlvar_array(
    x = mu[[1]],
    pooled = mu_summary,
    index = seq_along(mu[[1]]$mean)
  )

  # Varianz der zufälligen Intercepts mitteln und als SD ausgeben
  mu_pooled$SD <- sqrt(pool_mean(lapply(mu, function(x) x$SD^2)))

  # Für jede Person ihren Intercept-Vektor über die Imputationen mitteln
  mu_pooled$subject <- lapply(seq_along(mu[[1]]$subject), function(i) {
    pool_mean(lapply(mu, function(x) x$subject[[i]]))
  })

  mu_pooled$imp_infos_table <- mu_table

  gamma_omega <- lapply(fits$fits, function(x) x$results$Gamma_Omega_mu)
  gamma_omega_pooled <- pool_mlvar_gamma(gamma_omega, conf_level, dfcom)


  # Between-Person-Matrizen aus den gepoolten Parametern berechnen
  omega_derived <- derive_omega_mu(
    gamma = gamma_omega_pooled$mean,
    mu_sd = mu_pooled$SD
  )

  # mlVAR-Struktur übernehmen und die vier Punktschätzungen einsetzen
  omega_pooled <- fits$fits[[1]]$results$Omega_mu

  for (component in names(omega_derived)) {
    omega_pooled[[component]]$mean[] <- omega_derived[[component]]$mean
  }

  beta <- lapply(
    fits$fits,
    function(x) x$results$Beta
  )

  # Alle Beta-Einträge poolen, einschließlich der autoregressiven Diagonale
  beta_table <- pool_rubin(
    estimates = lapply(beta, function(x) x$mean),
    ses = lapply(beta, function(x) x$SE),
    dfcom = dfcom
  )

  beta_summary <- summarize_pooled(
    beta_table,
    conf_level = conf_level
  )

  beta_pooled <- update_mlvar_array(
    x = beta[[1]],
    pooled = beta_summary,
    index = rep(TRUE, length(beta[[1]]$mean))
  )

  beta_pooled$SD <- sqrt(pool_mean(lapply(beta, function(x) x$SD^2)))

  if (fits$fits[[1]]$input$temporal != "fixed") {
    beta_pooled$subject <- lapply(seq_along(beta[[1]]$subject), function(i) {
      pool_mean(lapply(beta, function(x) x$subject[[i]]))
    })
  }

  beta_pooled$imp_infos_table <- beta_table

  results <- list(
    mu = mu_pooled,
    Gamma_Omega_mu = gamma_omega_pooled,
    Omega_mu = omega_pooled,
    Beta = beta_pooled
  )

  theta <- lapply(fits$fits, function(x) x$results$Theta)
  theta_pooled <- theta[[1]]
  is_fixed_theta <- FALSE

  if (is.null(fits$fits[[1]]$results$Gamma_Theta)) {
    # "fixed" oder "unique"
    theta_cov <- theta_pooled$cov
    is_fixed_theta <- all(vapply(
      theta_cov$subject, identical, logical(1), y = theta_cov$mean
    ))

    if (is_fixed_theta) {
      group_cor <- pool_mean(lapply(theta, function(x) x$cor$mean))
    } else {
      # "unique": zuerst die Kovarianz jeder Person poolen
      subject_derived <- lapply(seq_along(theta_cov$subject), function(i) {
        cov_i <- pool_mean(lapply(theta, function(x) x$cov$subject[[i]]))
        convert_mlvar_matrix(cov_i, type = "cov")
      })

      subject_cov <- lapply(subject_derived, function(x) x$cov$mean)
      group_cor <- stats::cov2cor(pool_mean(subject_cov))
    }

    group_var <- pool_mean(lapply(theta, function(x) diag(x$cov$mean)))
    group_sd <- sqrt(group_var)
    group_cov <- group_cor * outer(group_sd, group_sd)
    group_derived <- convert_mlvar_matrix(group_cov, type = "cov")

  } else {
    # "correlated" oder "orthogonal": gerichtete Regressionsgewichte poolen
    gamma_theta <- lapply(fits$fits, function(x) x$results$Gamma_Theta)
    gamma_theta_pooled <- pool_mlvar_gamma(gamma_theta, conf_level, dfcom)

    # Random-Slope-Varianzen aus allen Modellblöcken jeder Imputation lesen
    slope_vars <- lapply(fits$fits, function(fit) {
      vars <- fit$results$Gamma_Theta$SD
      vars[] <- 0

      for (i in seq_len(nrow(vars))) {
        model <- fit$output$contemporaneous[[i]]
        vc <- as.data.frame(lme4::VarCorr(model))
        slopes <- vc[!is.na(vc$var1) & is.na(vc$var2), ]
        vars[i, as.character(slopes$var1)] <- slopes$vcov
      }

      vars
    })

    gamma_theta_pooled$SD[] <- sqrt(pool_mean(slope_vars))

    # Jede Personenmatrix enthält bereits Fixed Effects + Random Effects
    for (i in seq_along(gamma_theta_pooled$subject)) {
      gamma_theta_pooled$subject[[i]][] <- pool_mean(
        lapply(gamma_theta, function(x) x$subject[[i]]))
    }

    # Eine gepoolte Residualvarianz je Knoten aus der zweiten Modellstufe
    resid_var <- pool_mean(lapply(fits$fits, function(x) {
      models <- x$output$contemporaneous
      vapply(models, stats::sigma, numeric(1))^2
    }))

    # Gruppe und Personen verwenden dieselben knotenweisen Residualvarianzen
    group_derived <- derive_mlvar_network(gamma_theta_pooled$mean, resid_var)
    subject_derived <- lapply(
      gamma_theta_pooled$subject, derive_mlvar_network, resid_var = resid_var
    )

    idvar <- fits$fits[[1]]$input$idvar
    nodes <- rownames(gamma_theta_pooled$mean)
    model <- fits$fits[[1]]$output$contemporaneous[[1]]
    ids <- rownames(lme4::ranef(model, condVar = FALSE)[[idvar]])

    cutoff <- 10 * sum(diag(group_derived$cov$mean))

    for (i in seq_along(subject_derived)) {
      if (sum(diag(subject_derived[[i]]$cov$mean)) <= cutoff) next

      posthoc_cov <- lapply(fits$fits, function(fit) {
        rows <- as.character(fit$data[[idvar]]) == ids[i]
        res <- fit$step1_residuals[rows, nodes, drop = FALSE]
        stats::cov(res, use = "pairwise.complete.obs")
      })

      subject_derived[[i]] <- convert_mlvar_matrix(pool_mean(posthoc_cov), "cov")
    }

    results$Gamma_Theta <- gamma_theta_pooled
  }

  # Alle Varianten nutzen dieselbe mlVAR-Ausgabestruktur
  for (component in names(group_derived)) {
    entry <- theta_pooled[[component]]
    entry$mean[] <- group_derived[[component]]$mean

    if (is_fixed_theta) {
      entry$subject[] <- list(entry$mean)
      entry$SD[] <- 0
    } else {
      for (i in seq_along(entry$subject)) {
        entry$subject[[i]][] <- subject_derived[[i]][[component]]$mean
      }
      entry$SD[] <- apply(simplify2array(entry$subject), c(1, 2), stats::sd)
    }

    theta_pooled[[component]] <- entry
  }

  results$Theta <- theta_pooled

  list(results = results)
}


update_mlvar_array <- function(x, pooled, index) {
  x$mean[index] <- pooled$estimate
  x$SE[index] <- pooled$SE
  x$lower[index] <- pooled$lower
  x$upper[index] <- pooled$upper
  x$P[index] <- pooled$p
  x
}

#' Pool the directed off-diagonal coefficients of an mlVAR gamma matrix
#' @noRd
pool_mlvar_gamma <- function(gamma, conf_level, dfcom) {
  template <- gamma[[1]]
  mat <- template$mean
  off_diag <- row(mat) != col(mat)

  terms <- outer(rownames(mat), colnames(mat), paste, sep = ":")[off_diag]

  gamma_table <- pool_rubin(
    estimates = lapply(gamma, function(x) {
      setNames(x$mean[off_diag], terms)
    }),
    ses = lapply(gamma, function(x) {
      setNames(x$SE[off_diag], terms)
    }),
    dfcom = dfcom
  )

  pooled <- update_mlvar_array(
    x = template,
    pooled = summarize_pooled(gamma_table, conf_level = conf_level),
    index = off_diag
  )

  pooled$imp_infos_table <- gamma_table
  pooled
}

#' Derive mlVAR network matrices from nodewise regression coefficients
#' @noRd
derive_mlvar_network <- function(gamma, resid_var) {
  prec <- sweep(diag(nrow(gamma)) - gamma, 1L, resid_var, "/")
  prec <- correct_mlvar_matrix(prec)

  convert_mlvar_matrix(prec, type = "prec")
}

#' Derive the between-person network from pooled model parameters
#' @noRd
derive_omega_mu <- function(gamma, mu_sd) {
  if (any(mu_sd == 0)) {
    warning(
      "Zero pooled random-intercept SD for: ",
      paste(names(mu_sd)[mu_sd == 0], collapse = ", "),
      ". The between-person network could not be derived.",
      call. = FALSE
    )

    # Dimensionen und Knotennamen für alle vier Ausgaben erhalten
    gamma[] <- NA_real_
    missing <- list(mean = gamma)
    return(list(cov = missing, cor = missing, prec = missing, pcor = missing))
  }

  derive_mlvar_network(gamma, resid_var = mu_sd^2)
}

#' Convert a covariance or precision matrix to mlVAR matrix representations
#' @noRd
convert_mlvar_matrix <- function(mat, type = c("prec", "cov")) {
  type <- match.arg(type)
  node_names <- dimnames(mat)

  if (min(eigen(mat, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    warning("Pooled ", type, " matrix corrected with nearPD().", call. = FALSE)
    mat <- as.matrix(Matrix::nearPD(mat, keepDiag = TRUE)$mat)
    dimnames(mat) <- node_names
  }

  if (type == "prec") {
    prec <- mat
    cov <- solve(prec)
  } else {
    cov <- mat
    prec <- solve(cov)
  }

  pcor <- inv_to_net(prec)
  diag(pcor) <- 1
  dimnames(pcor) <- node_names

  list(
    cov  = list(mean = cov),
    cor  = list(mean = stats::cov2cor(cov)),
    prec = list(mean = prec),
    pcor = list(mean = pcor)
  )
}


#' Apply the eigenvalue correction used in mlVAR
#' @noRd
correct_mlvar_matrix <- function(mat) {
  mat <- (mat + t(mat)) / 2
  min_ev <- min(eigen(mat, symmetric = TRUE, only.values = TRUE)$values)

  if (min_ev < 0) {
    diag(mat) <- diag(mat) - min_ev
    mat <- mat + 0.001
  }

  mat
}

#' @title Add Inference Statistics to Pooled Estimates
#'
#' @description
#' Adds standard errors, confidence intervals, test statistics, and
#' p-values to the data frame returned by `pool_rubin()`.
#'
#' @param pooled A data frame returned by `pool_rubin()`, containing
#'   `estimate`, `t`, and `df` columns.
#' @param conf_level Confidence level for the interval. Defaults to 0.95.
#'
#' @return The input data frame with additional `SE`, `lower`, `upper`,
#'   `statistic`, and `p` columns.
#'
#' @details
#' The pooled standard error is calculated as `sqrt(t)`. Confidence
#' intervals and p-values use the pooled degrees of freedom in `df`.
#' When `df` is infinite, the normal distribution is used.
#'
#' @noRd
summarize_pooled <- function(pooled, conf_level = 0.95) {
  pooled$SE <- sqrt(pooled$t)
  pooled$statistic <- pooled$estimate / pooled$SE

  critical_value <- stats::qt(
    1 - (1 - conf_level) / 2,
    df = pooled$df
  )

  infinite_df <- is.infinite(pooled$df)
  critical_value[infinite_df] <- stats::qnorm(
    1 - (1 - conf_level) / 2
  )

  pooled$lower <- pooled$estimate - critical_value * pooled$SE
  pooled$upper <- pooled$estimate + critical_value * pooled$SE

  pooled$p <- 2 * stats::pt(
    -abs(pooled$statistic),
    df = pooled$df
  )
  pooled$p[infinite_df] <- 2 * stats::pnorm(
    -abs(pooled$statistic[infinite_df])
  )

  pooled
}

#' @title Pool Structured Estimates Across Multiple Imputations
#'
#' @description
#' Pools corresponding estimates and standard errors across multiple
#' imputed datasets using Rubin's rules. Each scalar parameter receives
#' one row in the returned data frame.
#'
#' @param estimates A list of numeric vectors, matrices, or arrays.
#'   Each element contains estimates from one imputed dataset.
#' @param ses A list of numeric vectors, matrices, or arrays containing
#'   the corresponding standard errors.
#' @param dfcom Degrees of freedom for the complete-data analysis. A
#'   single value is applied to all parameters. Alternatively, provide a
#'   numeric object with the same structure as each element of
#'   `estimates`. Defaults to `Inf`, which assumes large complete-data
#'   degrees of freedom.
#'
#' @return A data frame with one row per scalar parameter and columns
#'   for `term`, `nimp`, `estimate`, `ubar`, `b`, `t`, `dfcom`, `df`,
#'   `riv`, `lambda`, and `fmi`. For matrices and arrays, additional
#'   columns (`dim1`, `dim2`, ...) identify each parameter's position.
#'
#' @details
#' All estimate and standard-error objects must have the same dimensions
#' and names. Parameters are pooled independently. Non-finite estimates
#' or standard errors are not silently removed.
#'
#' @noRd
pool_rubin <- function(estimates, ses, dfcom = Inf) {

  # Estimates and standard errors must be supplied as lists.
  if (!is.list(estimates) || !is.list(ses)) {
    stop("`estimates` and `ses` must be lists.", call. = FALSE)
  }

  # Each estimate list must have a corresponding standard-error list.
  if (length(estimates) != length(ses)) {
    stop("`estimates` and `ses` must have the same length.", call. = FALSE)
  }

  # At least two imputations are needed to calculate between-imputation variance.
  if (length(estimates) < 2L) {
    stop("Rubin pooling requires at least two imputations.", call. = FALSE)
  }

  # If both lists are named, require the imputations to be in matching order.
  if (
    !is.null(names(estimates)) &&
    !is.null(names(ses)) &&
    !identical(names(estimates), names(ses))
  ) {
    stop(
      "When named, `estimates` and `ses` must have matching names and order.",
      call. = FALSE
    )
  }

  template <- estimates[[1]]

  if (!is.numeric(template) || length(template) == 0L) {
    stop("Each input must be a non-empty numeric object.", call. = FALSE)
  }

  # Verify that all estimates and SEs have the same structure as the template.
  for (i in seq_along(estimates)) {
    estimate_i <- estimates[[i]]
    se_i <- ses[[i]]

    if (!is.numeric(estimate_i) || !is.numeric(se_i)) {
      stop("All estimates and standard errors must be numeric.", call. = FALSE)
    }

    estimate_matches <- identical(dim(estimate_i), dim(template)) &&
      identical(dimnames(estimate_i), dimnames(template)) &&
      (
        !is.null(dim(template)) ||
          identical(names(estimate_i), names(template))
      )

    se_matches <- identical(dim(se_i), dim(template)) &&
      identical(dimnames(se_i), dimnames(template)) &&
      (
        !is.null(dim(template)) ||
          identical(names(se_i), names(template))
      )

    if (!estimate_matches || !se_matches) {
      stop(
        "All estimates and standard errors must have the same dimensions and names.",
        call. = FALSE
      )
    }
  }

  n_parameters <- length(template)

  # dfcom may be one value for all parameters or match the template structure.
  if (!is.numeric(dfcom)) {
    stop("`dfcom` must be numeric.", call. = FALSE)
  }

  if (length(dfcom) == 1L) {
    dfcom_values <- rep(dfcom, n_parameters)
  } else {
    dfcom_matches <- identical(dim(dfcom), dim(template)) &&
      identical(dimnames(dfcom), dimnames(template)) &&
      (
        !is.null(dim(template)) ||
          identical(names(dfcom), names(template))
      )

    if (!dfcom_matches) {
      stop(
        "`dfcom` must be a single value or have the same structure as the estimates.",
        call. = FALSE
      )
    }

    dfcom_values <- as.vector(dfcom)
  }

  # Flatten each imputation to a column:
  # rows correspond to parameters, columns correspond to imputations.
  Q <- do.call(cbind, lapply(estimates, as.vector))
  SE <- do.call(cbind, lapply(ses, as.vector))

  # Pool each parameter, storing parameters in rows and statistics in columns.
  first_result <- pool_rubin_scalar(
    q = Q[1, ],
    se = SE[1, ],
    dfcom = dfcom_values[1]
  )

  pooled <- matrix(
    NA_real_,
    nrow = n_parameters,
    ncol = length(first_result),
    dimnames = list(NULL, names(first_result))
  )

  pooled[1, ] <- first_result

  if (n_parameters > 1L) {
    for (j in 2:n_parameters) {
      pooled[j, ] <- pool_rubin_scalar(
        q = Q[j, ],
        se = SE[j, ],
        dfcom = dfcom_values[j]
      )
    }
  }

  # Give each row a readable parameter label.
  if (is.null(dim(template))) {
    term <- names(template)

    if (is.null(term)) {
      term <- paste0("parameter_", seq_len(n_parameters))
    }
  } else {
    # arrayInd() provides the coordinates of each flattened parameter.
    coordinates <- arrayInd(seq_len(n_parameters), .dim = dim(template))
    term <- rep("", n_parameters)

    for (d in seq_len(ncol(coordinates))) {
      dimension_names <- dimnames(template)[[d]]

      if (is.null(dimension_names)) {
        labels <- as.character(coordinates[, d])
      } else {
        labels <- dimension_names[coordinates[, d]]
      }

      if (d == 1L) {
        term <- labels
      } else {
        term <- paste(term, labels, sep = ":")
      }
    }
  }

  # Assemble the MICE-like table: one row per parameter.
  out <- data.frame(
    term = term,
    pooled,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # Keep the separate coordinates too, so matrix/array parameters remain identifiable.
  if (!is.null(dim(template))) {
    for (d in seq_len(ncol(coordinates))) {
      dimension_names <- dimnames(template)[[d]]

      if (is.null(dimension_names)) {
        labels <- as.character(coordinates[, d])
      } else {
        labels <- dimension_names[coordinates[, d]]
      }

      out[[paste0("dim", d)]] <- labels
    }
  }

  out
}

#' @title Pool a Scalar Estimate Across Multiple Imputations
#'
#' @description
#' Pools the estimates and standard errors of one scalar parameter across
#' multiple imputed datasets using Rubin's rules.
#'
#' @param q Numeric vector containing the estimate of the same parameter
#'   from each imputed dataset.
#' @param se Numeric vector containing the corresponding standard error
#'   from each imputed dataset.
#' @param dfcom Degrees of freedom for the complete-data analysis. A
#'   finite value applies the Barnard-Rubin small-sample adjustment.
#'   The default, \code{Inf}, assumes large complete-data degrees of
#'   freedom.
#'
#' @return A named numeric vector with:
#' \describe{
#'   \item{nimp}{Number of imputations.}
#'   \item{estimate}{Pooled estimate.}
#'   \item{ubar}{Mean within-imputation variance.}
#'   \item{b}{Between-imputation variance.}
#'   \item{t}{Total pooled variance.}
#'   \item{dfcom}{Complete-data degrees of freedom used.}
#'   \item{df}{Degrees of freedom for inference.}
#'   \item{riv}{Relative increase in variance due to missing data.}
#'   \item{lambda}{Proportion of total variance attributable to missing data.}
#'   \item{fmi}{Fraction of missing information.}
#' }
#'
#' @details
#' The total variance is calculated as
#' \deqn{t = \bar{U} + (1 + 1/nimp)B,}
#' where \eqn{\bar{U}} is the mean within-imputation variance and \eqn{B}
#' is the between-imputation variance. The degrees of freedom use Rubin's
#' large-sample formula when \code{dfcom = Inf}, and the Barnard-Rubin
#' adjustment when \code{dfcom} is finite.
#'
#' All imputations must estimate the same parameter, and each estimate
#' must have a corresponding valid standard error. Non-finite values are
#' treated as errors.
#'
#' @noRd
pool_rubin_scalar <- function(q, se, dfcom = Inf) {



  # check that estimates and standard errors are numeric vectors
  if (!is.numeric(q) || !is.null(dim(q)) || !is.numeric(se) || !is.null(dim(se))) {
    stop("`q` and `se` must be numeric vectors.", call. = FALSE)
    }

  # each estimate must have a matching standard error
  if (length(q) != length(se)) {
    stop("`q` and `se` must have the same length.", call. = FALSE)
  }

  nimp <- length(q) # number of imputations
  # at least two imputations are needed to estimate between-imputation variance
  if (nimp < 2L) {
    stop("Rubin pooling requires at least two imputations.", call. = FALSE)
  }

  # do not silently drop failed or invalid estimates
  if (any(!is.finite(q)) || any(!is.finite(se)) || any(se < 0)) {
    stop("`q` and `se` must be finite, and `se` must be non-negative.",
      call. = FALSE)
  }

  # dfcom can be Inf, but otherwise must be a positive scalar.
  if (!is.numeric(dfcom) || length(dfcom) != 1L || is.na(dfcom) || dfcom <= 0) {
    stop("`dfcom` must be a positive number or `Inf`.", call. = FALSE)
  }

  # average of the estimates across imputations as estimate
  estimate <- mean(q)

  # average squared standard error as within-imputation variance:
  ubar <- mean(se^2)

  # variability of the estimates across imputations as between-imputation variance
  b <- stats::var(q)

  # Total variance combines within- and between-imputation uncertainty.
  t <- ubar + (1 + 1 / nimp) * b

  if (ubar == 0 || t == 0) {
    stop("Pooling is undefined when the within- or total variance is zero.",
      call. = FALSE
    )
  }

  # relative increase in variance due to missing data
  riv <- (1 + 1 / nimp) * b / ubar

  # proportion of the total variance attributable to missing data
  lambda <- (1 + 1 / nimp) * b / t

  # Rubin's degrees of freedom under the large-sample assumption
  df_old <- if (riv == 0) {
    Inf
  } else {
    (nimp - 1) * (1 + 1 / riv)^2
  }

  # Observed-data degrees of freedom for the Barnard-Rubin adjustment.
  # When dfcom is infinite, this correction does not reduce df_old.
  df_obs <- if (is.infinite(dfcom)) {
    Inf
  } else {
    ((dfcom + 1) / (dfcom + 3)) * dfcom * (1 - lambda)
  }

  # Combine the two degrees-of-freedom components.
  df <- 1 / (1 / df_old + 1 / df_obs)

  # Fraction of missing information, including the df adjustment.
  fmi <- if (is.infinite(riv)) {
    1
  } else {
    (riv + 2 / (df + 3)) / (riv + 1)
  }

  # return pooled quantities
  c(
    nimp = nimp,
    estimate = estimate,
    ubar = ubar,
    b = b,
    t = t,
    dfcom = dfcom,
    df = df,
    riv = riv,
    lambda = lambda,
    fmi = fmi
  )
}

#' @title Average Estimates Across Imputations
#'
#' @description
#' Computes the element-wise arithmetic mean of a list of numeric objects
#' with identical structure.
#'
#' @param x A list of numeric vectors, matrices, or arrays with identical
#'   dimensions.
#'
#' @return A numeric object with the same structure as the elements of `x`,
#'   containing the element-wise mean across imputations.
#'
#' @noRd
pool_mean <- function(x) {
  Reduce(`+`, x) / length(x)
}

