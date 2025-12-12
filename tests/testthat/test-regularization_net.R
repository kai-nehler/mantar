#### Tests for regularization_net() ####

test_that("regularization_net() works with full data and glasso defaults", {

  expect_warning(
  res <- regularization_net(
    data = mantar_dummy_full_cont,
    penalty = "glasso"
  ),
  regexp = "Varying 'lambda' only, n_gamma is set to 1."
  )

  # structure and classes
  expect_type(res, "list")
  expect_named(res, c("pcor", "n", "cor_method", "full_results", "args"))
  expect_s3_class(res, "mantar_regularization")
  expect_s3_class(res, "mantar_network")

  p <- ncol(mantar_dummy_full_cont)
  expect_true(is.matrix(res$pcor))
  expect_equal(dim(res$pcor), c(p, p))

  # Defaults for glasso
  expect_equal(res$args$penalty, "glasso")
  expect_equal(res$args$likelihood, "obs_based")
  expect_equal(res$args$n_lambda, 100)
  expect_equal(res$args$ic_type, "ebic")
})

test_that("regularization_net() sets conditional defaults for extended and n_lambda", {

  expect_warning(
  res_glasso <- regularization_net(
    data   = mantar_dummy_full_cont,
    penalty = "glasso"
  ),
  regexp = "Varying 'lambda' only, n_gamma is set to 1."
  )

  expect_warning(
  res_atan <- regularization_net(
    data   = mantar_dummy_full_cont,
    penalty = "atan"
  ),
  regexp = "Varying 'lambda' only, n_gamma is set to 1."
  )

  # glasso
  expect_equal(res_glasso$args$ic_type, "ebic")
  expect_equal(res_glasso$args$n_lambda, 100)

  # atan
  expect_equal(res_atan$args$ic_type, "bic")
  expect_equal(res_atan$args$n_lambda, 50)
})

test_that("regularization_net() works with mat + ns input", {

  mat <- stats::cov(mantar_dummy_full_cont)
  ns  <- rep(nrow(mantar_dummy_full_cont) - 5, ncol(mat))

  expect_warning(
  res <- regularization_net(
    mat = mat,
    ns  = ns,
    n_calc = "average",
    penalty = "glasso",
    likelihood = "mat_based"
  ),
  regexp = "Varying 'lambda' only, n_gamma is set to 1."
  )

  expect_true(is.matrix(res$pcor))
  expect_equal(res$n, mean(ns))
  expect_null(res$cor_method)
})


test_that("regularization_net() errors for invalid penalty type", {
  expect_error(
    regularization_net(
      data = mantar_dummy_full_cont,
      penalty = "foo"
    ),
    "Invalid penalty type",
    fixed = FALSE
  )
})

test_that("regularization_net() errors for obs_based likelihood without data", {
  expect_error(
    regularization_net(
      mat = stats::cov(mantar_dummy_full_cont),
      likelihood = "obs_based",
      penalty = "glasso",
      ns = nrow(mantar_dummy_full_cont)
    ),
    "observed data loglikelihood is only implemented when data is provided",
    fixed = FALSE
  )
})

test_that("regularization_net() errors for obs_based likelihood with ordered variables", {
  expect_error(
    regularization_net(
      data = mantar_dummy_full_cont,
      ordered = TRUE,
      likelihood = "obs_based",
      penalty = "glasso"
    ),
    "only implemented for data\\s+treated as continuous",
    fixed = FALSE
  )
})

test_that("regularization_net() requires ns when mat is provided", {
  expect_error(
    regularization_net(
      mat = stats::cov(mantar_dummy_full_cont),
      penalty = "glasso",
      likelihood = "mat_based"
    ),
    "If 'mat' is provided, 'ns' must also be specified",
    fixed = FALSE
  )
})


#### Tests for regularization_sel() ####

test_that("regularization_sel() runs and returns expected structure for glasso", {

  mat <- stats::cov2cor(stats::cov(mantar_dummy_full_cont))
  n   <- nrow(mantar_dummy_full_cont)

  res <- regularization_sel(
    mat        = mat,
    data       = NULL,
    means      = NULL,
    n          = n,
    likelihood = "mat_based",
    ic_type   = "bic",
    extended_gamma = 0.5,
    penalty    = "glasso",
    vary       = "lambda",
    n_lambda   = 5,
    lambda_min_ratio = 0.1,
    n_gamma    = 1,
    pen_diag   = FALSE
  )

  expect_type(res, "list")
  expect_named(res, c("opt_net", "full_results"))
  expect_true(is.matrix(res$opt_net))

  # full_results has to be of length n_lambda * n_gamma
  expect_equal(length(res$full_results), 5)

  # every result consitists of wi, w, ic, rho_mat
  first <- res$full_results[[1]]
  expect_true(is.matrix(first$wi))
  expect_true(is.matrix(first$w))
  expect_type(first$ic, "double")
  expect_true(is.matrix(first$rho_mat))
})


#### Tests for def_pen_mats() ####

test_that("def_pen_mats() constructs correct grid and penalty matrices for glasso", {
  mat <- stats::cov2cor(stats::cov(mantar_dummy_full_cont))
  p   <- ncol(mat)

  pm <- def_pen_mats(
    mat = mat,
    penalty = "glasso",
    vary = "lambda",
    n_lambda = 10,
    n_gamma  = 1,
    n = nrow(mantar_dummy_full_cont),
    pen_diag = FALSE
  )

  expect_named(pm, c("grid", "pen_mats"))

  # Grid: 10 combinations
  expect_s3_class(pm$grid, "data.frame")
  expect_equal(nrow(pm$grid), 10)

  # Penalty-matrices - one row per grid combination
  expect_equal(length(pm$pen_mats), 10)
  expect_true(all(vapply(pm$pen_mats, is.matrix, logical(1L))))
  expect_true(all(vapply(pm$pen_mats, function(m) all(dim(m) == p), logical(1L))))

  # Diagonal = 0 with pen_diag = FALSE
  expect_true(all(vapply(pm$pen_mats, function(m) all(diag(m) == 0), logical(1L))))
})


test_that("def_pen_mats() uses user-specified lambda and gamma", {
  mat <- diag(3)

  expect_message(
    pm <- def_pen_mats(
      mat = mat,
      penalty = "glasso",
      vary = "lambda",
      lambda = c(0.1, 0.2),
      gamma = c(0.5, 1)
    ),
    "Using user-specified lambda values",
    fixed = FALSE
  )

  expect_equal(nrow(pm$grid), 4)
  expect_equal(pm$grid$lambda, rep(c(0.1, 0.2), times = 2))
  expect_equal(pm$grid$gamma, rep(c(0.5, 1), each = 2))
})


test_that("def_pen_mats() works for atan penalty", {
  mat <- stats::cov2cor(stats::cov(mantar_dummy_full_cont))
  p   <- ncol(mat)

  pm <- def_pen_mats(
    mat = mat,
    penalty = "atan",
    vary = "both",
    n_lambda = 3,
    n_gamma  = 2,
    n = nrow(mantar_dummy_full_cont),
    pen_diag = FALSE
  )

  expect_equal(nrow(pm$grid), 3 * 2)
  expect_equal(length(pm$pen_mats), 3 * 2)
  expect_true(all(vapply(pm$pen_mats, function(m) all(dim(m) == p), logical(1L))))
})
