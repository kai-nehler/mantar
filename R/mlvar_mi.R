mlvar_mi <- function(
    data,
    idvar,
    network_vars,
    nimp,
    maxit,
    imp_format = c("long", "wide"),
    imp_method,
    seed = NULL,
    estimator = "lmer",
    AR = FALSE,
    temporal = c(
      "default", "correlated", "orthogonal", "fixed", "unique"
    ),
    contemporaneous = c(
      "default", "correlated", "orthogonal", "fixed", "unique"
    ),
    keep_imputations = FALSE,
    keep_fits = FALSE,
    mice_args = list(),
    mlvar_args = list()) {

  if (estimator != "lmer"){
    stop("mlvar_mi() supports only estimator = 'lmer'.")
  }
  if (!identical(AR, FALSE)) {
    stop("mlvar_mi() supports only AR = FALSE.")
  }

  imp_format <- match.arg(tolower(imp_format), choices = c("long", "wide"))


  imputations <- mlvar_impute(
    data = data,
    nimp = nimp,
    maxit = maxit,
    imp_format = imp_format,
    imp_method = imp_method,
    seed = seed,
    mice_args = mice_args
  )

  fits <- mlvar_fits(
    imputations = imputations,
    idvar = idvar,
    network_vars = network_vars,
    mlvar_args = mlvar_args,
    estimator = estimator,
    AR = AR,
    temporal = temporal,
    contemporaneous = contemporaneous
  )

  pooled <- mlvar_pool(
    fits = fits
  )

  out <- list(
    call = match.call(),
    results = pooled$results
  )

  if (keep_imputations) {
    out$imputations <- imputations
  }

  if (keep_fits) {
    out$fits <- fits
  }

  class(out) <- "mlvar_mi"
  out
}


mlvar_impute <- function(
    data,
    nimp,
    imp_format = c("long", "wide"),
    imp_method = NULL,
    seed = NULL,
    maxit = 5L,
    mice_args = list()) {

  imp_format <- match.arg(tolower(imp_format), choices = c("long", "wide"))

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  checker(nimp = nimp)

  args <- merge_args(
    fixed = list(
      data = data,
      m = nimp,
      method = imp_method,
      maxit = maxit,
      seed = seed,
      printFlag = FALSE
    ),
    extra = mice_args,
    arg_name = "mice_args"
  )

  mids <- do.call(mice::mice, args)

  result <- list(
    datasets = mice::complete(mids, action = "all"),
    mids = mids,
    nimp = nimp,
    imp_format = imp_format,
    imp_method = mids$method,
    seed = seed,
    call = match.call()
  )

  result
}


mlvar_fits <- function(
    imputations,
    network_vars,
    idvar,
    lags = 1,
    AR,
    estimator = c("default", "lmer", "lm", "Mplus"),
    temporal = c(
      "default", "correlated", "orthogonal", "fixed", "unique"
    ),
    contemporaneous = c(
      "default", "correlated", "orthogonal", "fixed", "unique"
    ),
    mlvar_args = list()) {
  estimator <- match.arg(estimator)
  temporal <- match.arg(temporal)
  contemporaneous <- match.arg(contemporaneous)

  args <- merge_args(
    fixed = list(
      data = NULL,
      vars = network_vars,
      idvar = idvar,
      lags = lags,
      estimator = estimator,
      AR = AR,
      temporal = temporal,
      contemporaneous = contemporaneous
    ),
    extra = mlvar_args,
    arg_name = "mlvar_args"
  )

  fits <- lapply(imputations$datasets, function(dat) {
    args$data <- dat
    do.call(mlVAR::mlVAR, args)
  })

  names(fits) <- paste0("imp", seq_along(fits))

  result <- list(
    fits = fits,
    nimp = length(fits),
    call = match.call()
  )

  result
}




merge_args <- function(fixed, extra, arg_name) {
  if (length(extra) > 0L &&
      (is.null(names(extra)) || anyNA(names(extra)) || any(names(extra) == ""))) {
    stop("All entries in `", arg_name, "` must be named.", call. = FALSE)
  }

  duplicates <- intersect(names(extra), names(fixed))

  if (length(duplicates) > 0L) {
    warning(
      "Ignoring ", paste(sprintf("`%s`", duplicates), collapse = ", "),
      " in `", arg_name, "` because the wrapper already supplies ",
      "these arguments.",
      call. = FALSE
    )
    extra <- extra[!names(extra) %in% duplicates]
  }

  c(fixed, extra)
}
