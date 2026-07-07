#### Tests for regularization_net() ####

test_that("regularization_net() works with full data and glasso defaults", {


  res <- regularization_net(
    data = mantar_dummy_full_cont,
    penalty = "glasso"
  )

  # structure and classes
  expect_type(res, "list")
  expect_named(res, c("pcor", "n", "cor_method", "imputed_data", "full_results", "args"))
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

  res_glasso <- regularization_net(
    data   = mantar_dummy_full_cont,
    penalty = "glasso"
  )

  expect_message(
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
  ns  <- matrix(data = nrow(mantar_dummy_full_cont), nrow = ncol(mat), ncol = ncol(mat))

  res <- regularization_net(
    mat = mat,
    ns  = ns,
    n_calc = "average",
    penalty = "glasso",
    likelihood = "mat_based"
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

test_that("regularization_net() errors with glasso and vary not as lambda", {
  expect_error(
    regularization_net(
      data = mantar_dummy_full_cont,
      likelihood = "obs_based",
      penalty = "glasso",
      vary = "gamma"
    ),
    "For 'glasso' penalty, 'vary' must be set to 'lambda' as this is the only penalty parameter. If you want to provide your own lambda values you can do this in the corresponding argument 'lambda' but you still have to set 'vary' to 'lambda'.",
    fixed = FALSE
  )

  expect_error(
    regularization_net(
      data = mantar_dummy_full_cont,
      likelihood = "obs_based",
      penalty = "glasso",
      vary = "both"
    ),
    "For 'glasso' penalty, 'vary' must be set to 'lambda' as this is the only penalty parameter. If you want to provide your own lambda values you can do this in the corresponding argument 'lambda' but you still have to set 'vary' to 'lambda'.",
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

test_that("regularization_net() fails false definition of ns", {
  expect_error(
    regularization_net(
      mat = stats::cov(mantar_dummy_full_cont),
      penalty = "glasso",
      likelihood = "mat_based",
      ns = colSums(!is.na(mantar_dummy_full_cont))
    ),
    "'ns' must be either a single value or a matrix with dimensions matching the matrix used for network estimation in 'mat' (after optional selection via 'network_vars').",
    fixed = TRUE
  )
})


test_that("regularization_net() works with ns as matrix", {

  res <- regularization_net(
      mat = stats::cov(mantar_dummy_full_cont),
      penalty = "glasso",
      likelihood = "mat_based",
      ns = matrix(100, nrow = ncol(mantar_dummy_full_cont), ncol = ncol(mantar_dummy_full_cont))
    )

  expect_equal(res$n, 100)

})


test_that("regularization_net() requires means when mat and data is provided with obs_based likelihood", {
  expect_error(
    regularization_net(
      mat = stats::cor(mantar_dummy_full_cont),
      data = mantar_dummy_full_cont,
      penalty = "glasso",
      likelihood = "obs_based",
      ns = nrow(mantar_dummy_full_cont)
    ),
    "When likelihood = 'obs_based' and a user-supplied matrix is used for regularization, 'means' must be provided, as no estimation from raw data is performed.",
    fixed = FALSE
  )
})

test_that("regularization_net() does not require means when mat and data is provided with mat_based likelihood", {
  expect_message(
    regularization_net(
      mat = stats::cor(mantar_dummy_full_cont),
      data = mantar_dummy_full_cont,
      penalty = "glasso",
      likelihood = "mat_based"
    ),
    "Both 'data' and 'mat' are provided. 'mat' will be used for regularization and likelihood computation. 'data' is ignored, and 'ns' will be computed from 'data'.",
    fixed = FALSE
  )

  expect_message(
    regularization_net(
      mat = stats::cor(mantar_dummy_full_cont),
      data = mantar_dummy_full_cont,
      penalty = "glasso",
      likelihood = "mat_based",
      ns = nrow(mantar_dummy_full_cont)
    ),
    "Both 'data' and 'mat' are provided. 'mat' will be used for regularization.",
    fixed = FALSE
  )

})

test_that("regularization_net() works when means, mat and data is provided with obs_based likelihood", {

  expect_message(
    regularization_net(
      mat = stats::cor(mantar_dummy_full_cont),
      data = mantar_dummy_full_cont,
      penalty = "glasso",
      likelihood = "obs_based",
      means = colMeans(mantar_dummy_full_cont),
      ns = nrow(mantar_dummy_full_cont)
    ),
    "Both 'data' and 'mat' are provided. 'mat' will be used for regularization, while 'data' will be used in the calculation of the observed-data loglikelihood.",
    fixed = FALSE
  )
})



test_that("regularization_net() works identically for obs_based likelihood while using cor_calc or providing values", {

  mantar_dummy_mis_cont_half <- mantar_dummy_mis_cont[1:(nrow(mantar_dummy_mis_cont)/2), ]

  expect_message(
    res_provided <- regularization_net(
      mat = stats::cor(mantar_dummy_mis_cont_half, use = "complete"),
      data = mantar_dummy_mis_cont_half,
      penalty = "glasso",
      likelihood = "obs_based",
      means = colMeans(na.omit(scale(mantar_dummy_mis_cont_half))),
      ns = mean({m <- t(!is.na(mantar_dummy_mis_cont_half)) %*% !is.na(mantar_dummy_mis_cont_half); m[upper.tri(m, diag = TRUE)]})
    ),
    "Both 'data' and 'mat' are provided. 'mat' will be used for regularization, while 'data' will be used in the calculation of the observed-data loglikelihood.",
    fixed = FALSE
  )

  res_aut <- regularization_net(
    data = mantar_dummy_mis_cont_half,
    penalty = "glasso",
    missing_handling = "listwise")

  expect_equal(res_provided$pcor, res_aut$pcor)

})

test_that("regularization_net() gives identical results for obs_based and mat_based likelihood when using cor_calc", {

  data_obs <- regularization_net(data = mantar_dummy_full_cont,
                                 likelihood = "obs_based",
                                 n_calc = "average",
                                 penalty = "glasso",
                                 n_lambda = 60)

  mat_data_obs <- regularization_net(data = mantar_dummy_full_cont,
                                     mat = cor(mantar_dummy_full_cont),
                                     means = colMeans(mantar_dummy_full_cont),
                                     likelihood = "obs_based",
                                     n_calc = "average",
                                     penalty = "glasso",
                                     n_lambda = 60)

  mat_mat <- regularization_net(mat = cor(mantar_dummy_full_cont),
                                ns = nrow(mantar_dummy_full_cont),
                                likelihood = "mat_based",
                                n_calc = "average",
                                penalty = "glasso",
                                n_lambda = 60)

  expect_equal(data_obs$pcor, mat_data_obs$pcor)
  expect_equal(data_obs$pcor, mat_mat$pcor)

  expect_equal(
    data_obs$pcor,
    mat_data_obs$pcor,
    tolerance = 1e-4
  )
  expect_equal(
    data_obs$pcor,
    mat_mat$pcor,
    tolerance = 1e-4
  )

})

test_that("regularization_net() gives identical results for obs_based and mat_based likelihood when using cor_calc with missing data handling", {

  data_obs <- regularization_net(data = mantar_dummy_mis_cont,
                                 likelihood = "obs_based",
                                 n_calc = "average",
                                 penalty = "glasso",
                                 n_lambda = 50,
                                 missing_handling = "two-step-em")

  mis_handling <- cor_calc(data = mantar_dummy_mis_cont)

  mat_data_obs <- regularization_net(data = mantar_dummy_mis_cont,
                                     mat = mis_handling$mat,
                                     means = mis_handling$means,
                                     likelihood = "obs_based",
                                     n_calc = "average",
                                     penalty = "glasso",
                                     n_lambda = 50,
                                     missing_handling = "two-step-em")

  expect_equal(data_obs$pcor, mat_data_obs$pcor)
})


test_that("regularization_net() mat based differs from observation based but not from matbased with data (data ignored)", {

  mis_handling <- cor_calc(data = mantar_dummy_mis_cont[1:180,])

  mat_data_obs <- regularization_net(data = mantar_dummy_mis_cont[1:180,],
                                     mat = mis_handling$mat,
                                     means = mis_handling$means,
                                     likelihood = "obs_based",
                                     n_calc = "average",
                                     count_diagonal = FALSE,
                                     penalty = "glasso",
                                     n_lambda = 100,
                                     missing_handling = "two-step-em")

  mat_mat <- regularization_net(mat = mis_handling$mat,
                                ns = mat_calculate_sample_size(mantar_dummy_mis_cont[1:180,], n_calc = "average", count_diagonal = FALSE),
                                likelihood = "mat_based",
                                n_calc = "average",
                                penalty = "glasso",
                                n_lambda = 100,
                                missing_handling = "two-step-em")

  mat_mat_dataignore <- regularization_net(data = mantar_dummy_mis_cont[1:180,],
                                           mat = mis_handling$mat,
                                           ns = mat_calculate_sample_size(mantar_dummy_mis_cont[1:180,], n_calc = "average", count_diagonal = FALSE),
                                           likelihood = "mat_based",
                                           n_calc = "average",
                                           penalty = "glasso",
                                           n_lambda = 100)

  expect_failure(expect_equal(mat_data_obs$pcor, mat_mat$pcor))
  expect_equal(mat_mat$pcor, mat_mat_dataignore$pcor)
})



test_that("regularization_net() fails when auxiliary_vars are provided without network_vars", {
  expect_error(
    regularization_net(
      data = mantar_dummy_full_cont,
      likelihood = "mat_based",
      penalty = "glasso",
      auxiliary_vars = 1
    ),
    "'auxiliary_vars' can only be used when 'network_vars' is specified.",
    fixed = TRUE
  )
})


test_that("regularization_net() fails when network_vars contain duplicates", {
  expect_error(
    regularization_net(
      data = mantar_dummy_full_cont,
      likelihood = "mat_based",
      penalty = "glasso",
      network_vars = c(1, 1)
    ),
    "'network_vars' must not contain duplicate variables.",
    fixed = TRUE
  )
})

test_that("regularization_net() fails when auxiliary_vars contain duplicates", {
  expect_error(
    regularization_net(
      data = mantar_dummy_mis_cont,
      likelihood = "mat_based",
      penalty = "glasso",
      network_vars = 1:3,
      auxiliary_vars = c(4, 4)
    ),
    "'auxiliary_vars' must not contain duplicate variables.",
    fixed = TRUE
  )
})

test_that("regularization_net() fails when network_vars and auxiliary_vars overlap", {
  expect_error(
    regularization_net(
      data = mantar_dummy_mis_cont,
      likelihood = "mat_based",
      penalty = "glasso",
      network_vars = 1:4,
      auxiliary_vars = c(4, 5)
    ),
    "'network_vars' and 'auxiliary_vars' must not contain overlapping variables.",
    fixed = TRUE
  )
})

test_that("regularization_net() works with network_vars only", {
  res <- regularization_net(
    data = mantar_dummy_full_cont,
    likelihood = "mat_based",
    penalty = "glasso",
    network_vars = 1:4
  )

  expect_s3_class(res, "mantar_regularization")
  expect_equal(dim(res$pcor), c(4L, 4L))
})

test_that("regularization_net() works with network_vars and auxiliary_vars", {
  res <- regularization_net(
    data = mantar_dummy_mis_cont,
    likelihood = "mat_based",
    penalty = "glasso",
    network_vars = 1:3,
    auxiliary_vars = 4:5,
    missing_handling = "two-step-em"
  )

  expect_s3_class(res, "mantar_regularization")
  expect_equal(dim(res$pcor), c(3L, 3L))
})

test_that("regularization_net() gives same result for preselected data and network_vars", {
  res_network_vars <- regularization_net(
    data = mantar_dummy_full_cont,
    likelihood = "mat_based",
    penalty = "glasso",
    network_vars = c(1, 2, 3)
  )

  res_reduced_data <- regularization_net(
    data = mantar_dummy_full_cont[, c(1, 2, 3)],
    likelihood = "mat_based",
    penalty = "glasso"
  )

  expect_equal(res_network_vars$pcor, res_reduced_data$pcor)
})


test_that("regularization_net() can differ when auxiliary_vars are included", {
  res_reduced_data <- regularization_net(
    data = mantar_dummy_mis_cont[, c(1, 2, 3, 4)],
    likelihood = "mat_based",
    penalty = "glasso",
    missing_handling = "two-step-em"
  )

  res_auxiliary <- regularization_net(
    data = mantar_dummy_mis_cont,
    likelihood = "mat_based",
    penalty = "glasso",
    network_vars = c(1, 2, 3, 4),
    auxiliary_vars = c(5, 6),
    missing_handling = "two-step-em"
  )

  expect_false(isTRUE(all.equal(res_reduced_data$pcor, res_auxiliary$pcor)))
})


test_that("regularization_net() respects the order of network_vars", {
  res_12 <- regularization_net(
    data = mantar_dummy_full_cont,
    likelihood = "mat_based",
    penalty = "glasso",
    network_vars = c(1, 2)
  )

  res_21 <- regularization_net(
    data = mantar_dummy_full_cont,
    likelihood = "mat_based",
    penalty = "glasso",
    network_vars = c(2, 1)
  )

  expect_equal(
    res_12$pcor,
    res_21$pcor[c(2, 1), c(2, 1)]
  )
})


test_that("regularization_net() returns imputed data with stacked-mi and auxiliary_vars", {
  res_aux <- regularization_net(
    data = mantar_dummy_mis_cont,
    likelihood = "mat_based",
    penalty = "glasso",
    network_vars = 1:6,
    auxiliary_vars = 7:8,
    missing_handling = "stacked-mi",
    nimp = 2,
    imp_method = "pmm"
  )

  res_no_aux <- regularization_net(
    data = mantar_dummy_mis_cont[, 1:6],
    likelihood = "mat_based",
    penalty = "glasso",
    missing_handling = "stacked-mi",
    nimp = 2,
    imp_method = "pmm"
  )

  expect_s3_class(res_aux, "mantar_regularization")
  expect_s3_class(res_aux$imputed_data, "mids")

  expect_equal(res_aux$imputed_data$m, 2)
  expect_equal(ncol(res_aux$imputed_data$data), 8L)
  expect_equal(dim(res_aux$pcor), c(6L, 6L))

  expect_false(isTRUE(all.equal(res_aux$pcor, res_no_aux$pcor)))
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
      penalty = "atan",
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


test_that("def_pen_mats() sends correct warning with glasso and gamma", {
  mat <- diag(3)

  expect_warning(
    pm <- def_pen_mats(
      mat = mat,
      penalty = "glasso",
      vary = "lambda",
      lambda = c(0.1, 0.2),
      gamma = c(0.5, 1)
    ),
    "Gamma values are not used for the glasso penalty and will be ignored.",
    fixed = FALSE
  )

  expect_equal(nrow(pm$grid), 2)
  expect_equal(pm$grid$lambda, rep(c(0.1, 0.2)))
  expect_equal(pm$grid$gamma, rep(NA, each = 2))
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
