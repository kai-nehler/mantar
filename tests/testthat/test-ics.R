test_that("reg_sample_size_full_works", {
  expect_equal(reg_calculate_sample_size(data = mantar_dummy_full_cont), rep(400, 8))
})

test_that("reg_sample_size_mis_works", {
  expect_equal(reg_calculate_sample_size(data = mantar_dummy_mis_cont, n_calc = "Total"), rep(400, 8))
})

test_that("reg_error_for_different_ncalc", {
  expect_error(
    reg_calculate_sample_size(data = mantar_dummy_mis_cont, n_calc = "abc"),
    class = "error"
  )
})

test_that("ic_calculated as number",{
  expect_equal(
    class(reg_ic_calc(resid_var = 1, n = 400, n_preds = 2, ic_type = "bic")),
    "numeric"
  )
})


test_that("mat_ic_calc: obs_based and mat_based agree for complete data", {

  # complete data
  data <- mantar_dummy_full_mix
  data <- scale(data)  # standardize for identity model
  n <- nrow(data)
  p <- ncol(data)

  # Model-implied mean (here: sample mean)
  mu <- rep(0, p)
  # Correlation Matrix
  S <- cor(data)
  # precision matrix
  theta <- solve(S)

  # IC from observation-based likelihood
  IC_obs <- mat_ic_calc(
    data      = data,
    sample_cor = S,
    theta     = theta,
    mu        = mu,
    n         = n,
    ic_type   = "bic",
    likelihood = "obs_based"
  )

  # IC from matrix-based likelihood
  IC_mat <- mat_ic_calc(
    data       = data,
    sample_cor = S,
    theta      = theta,
    mu         = mu,
    n          = n,
    k          = "log(n)",
    extended   = FALSE,
    likelihood = "mat_based"
  )

  # For the identity model and consistent scaling, both formulations
  # should give (numerically) the same IC
  expect_equal(IC_obs, IC_mat, tolerance = 1e-10)
})




test_that("pattern_ll matches closed-form log-likelihood for complete data", {

  data <- matrix(
    c(0.5, 1.0,
      1.5, 2.0,
      0.0, -1.0),
    ncol = 2,
    byrow = TRUE
  )
  data <- scale(data)  # standardize for identity model

  # Model-implied mean and covariance
  mu <- c(0, 0)
  sigma <- matrix(c(1, 0.3,
                    0.3, 1), ncol = 2, byrow = TRUE)
  theta <- solve(sigma)

  cova <- sigma*(nrow(data)-1)/nrow(data)
  var <-  t(rbind(replicate(nrow(data), sqrt(diag(cova)), simplify = TRUE)))
  data_st <- data / var

  rows <- 1:nrow(data_st)
  cols <- rep(TRUE, 2)

  # value from pattern_ll()
  ll_pattern <- pattern_ll(
    data  = data_st,
    rows  = rows,
    cols  = cols,
    sigma = sigma,
    theta = theta,
    mu    = mu
  )

  # Closed-form log-likelihood (manually)
  observations   <- data[rows, cols, drop = FALSE]
  sigma_oo       <- sigma[cols, cols, drop = FALSE]
  theta_oo       <- solve(sigma_oo)
  pre_log        <- ncol(observations) * log(2 * pi)
  logdet         <- log(det(sigma_oo))
  observations_c <- sweep(observations, 2, mu[cols], FUN = "-")
  cova <- stats::cov(observations_c)*(nrow(observations_c)-1)/nrow(observations_c)
  var <-  t(rbind(replicate(nrow(observations_c), sqrt(diag(cova)), simplify = TRUE)))
  observations_c <- observations_c / var
  mahalanobis    <- rowSums((observations_c %*% theta_oo) * observations_c)
  ll_manual      <- sum(-(pre_log + logdet + mahalanobis) / 2)

  expect_equal(ll_pattern, ll_manual, tolerance = 1e-10)
})


test_that("pattern_ll handles missing variable pattern correctly", {

  # 3-dimensional example, but only variables 1 and 3 observed
  data <- matrix(
    c( 1,  2,  3,
       -1,  0,  1,
       0, -2,  4),
    ncol = 3,
    byrow = TRUE
  )

  mu <- c(0, 0, 0)

  # Positive definite covariance & its inverse
  sigma <- matrix(c(1.0, 0.2, 0.1,
                    0.2, 1.0, 0.3,
                    0.1, 0.3, 1.0),
                  ncol = 3, byrow = TRUE)
  theta <- solve(sigma)

  cova <- sigma*(nrow(data)-1)/nrow(data)
  var <-  t(rbind(replicate(nrow(data), sqrt(diag(cova)), simplify = TRUE)))
  data_st <- data / var

  rows <- 1:nrow(data_st)
  cols <- c(TRUE, FALSE, TRUE)  # variable 2 is unobserved for this pattern

  # value from pattern_ll()
  ll_pattern <- pattern_ll(
    data  = data_st,
    rows  = rows,
    cols  = cols,
    sigma = sigma,
    theta = theta,
    mu    = mu
  )

  # Manual calculation: only use observed variables 1 and 3
  cova <- sigma*(nrow(data)-1)/nrow(data)
  var <-  t(rbind(replicate(nrow(data), sqrt(diag(cova)), simplify = TRUE)))
  observations <- data / var
  observations   <- observations[rows, cols, drop = FALSE]
  sigma_oo       <- sigma[cols, cols, drop = FALSE]
  theta_oo       <- solve(sigma_oo)  # direct inverse of submatrix
  pre_log        <- ncol(observations) * log(2 * pi)
  logdet         <- log(det(sigma_oo))
  observations_c <- sweep(observations, 2, mu[cols], FUN = "-")
  mahalanobis    <- rowSums((observations_c %*% theta_oo) * observations_c)
  ll_manual      <- sum(-(pre_log + logdet + mahalanobis) / 2)

  expect_equal(ll_pattern, ll_manual, tolerance = 1e-10)
})


test_that("pattern_ll returns zero contribution when no variables observed and a warning", {

  data <- matrix(
    c(1, 2, 3,
      4, 5, 6),
    ncol = 3,
    byrow = TRUE
  )

  mu <- c(0, 0, 0)
  sigma <- diag(3)
  theta <- diag(3)

  rows <- 1:nrow(data)
  cols <- c(FALSE, FALSE, FALSE)  # nothing observed

  expect_warning(
    ll_pattern <- pattern_ll(
      data  = data,
      rows  = rows,
      cols  = cols,
      sigma = sigma,
      theta = theta,
      mu    = mu
    ),
    "There are observations with no observed variables",
    fixed = FALSE
  )

  expect_equal(ll_pattern, 0)
})

test_that("mat_inverse_update matches lavaan for typical deletion patterns", {


  S <- crossprod(matrix(c(1, 0.5, 0.3, 0.2,
                             0.5, 1, 0.4, 0.3,
                             0.3, 0.4, 1, 0.5,
                             0.2, 0.3, 0.5, 1), nrow = 4))
  theta <- solve(S)
  S_logdet <- as.numeric(determinant(S, logarithm = TRUE)$modulus)

  # verschiedene rm_idx-Kombinationen: nichts, eine, mehrere (aber nicht alle)
  rm_list <- list(
    integer(0L),  # nothing removed
    2L,           # single index
    c(1L, 4L)     # multiple indices
  )

  for (rm_idx in rm_list) {
    lav_res <- lavaan:::lav_matrix_symmetric_inverse_update(
      S.inv   = theta,
      rm.idx  = rm_idx,
      logdet  = TRUE,
      S.logdet = S_logdet
    )

    our_res <- mat_inverse_update(
      theta      = theta,
      rm_idx     = rm_idx,
      logdet     = TRUE,
      mat_logdet = S_logdet
    )

    # Matrixinhalte vergleichen
    expect_true(
      isTRUE(all.equal(as.matrix(our_res), as.matrix(lav_res), tolerance = 1e-8)),
      info = paste("Mismatch for rm_idx =", paste(rm_idx, collapse = ","))
    )

    # logdet-Attribut vergleichen
    expect_equal(attr(our_res, "logdet"), attr(lav_res, "logdet"))
  }
})



test_that("mat_inverse_update handles full deletion case correctly", {
  S <- crossprod(matrix(c(1, 0.5, 0.3,
                             0.5, 1, 0.4,
                             0.3, 0.4, 1), nrow = 3))
  theta <- solve(S)
  p <- ncol(theta)

  res <- mat_inverse_update(
    theta      = theta,
    rm_idx     = seq_len(p),
    logdet     = FALSE
  )

  # should be no matrix left
  expect_equal(dim(res), c(0L, 0L))
})

test_that("mat_inverse_update checks inputs and throws informative errors", {
  S <- crossprod(matrix(c(1, 0.5, 0.3,
                             0.5, 1, 0.4,
                             0.3, 0.4, 1), nrow = 3))
  theta <- solve(S)

  # Test invalid inputs - double index
  expect_error(
    mat_inverse_update(theta = theta, rm_idx = c(1, 1), logdet = FALSE, mat_logdet = NULL),
    "`rm_idx` must contain unique indices."
  )

  # Test invalid inputs - out of bounds index
  expect_error(
    mat_inverse_update(theta = theta, rm_idx = 4, logdet = FALSE, mat_logdet = NULL),
    "All entries of `rm_idx` must be between 1 and 3."
  )

  # Test invalid inputs - no logdet provided but an update is requested
  expect_error(
    mat_inverse_update(theta = theta, rm_idx = 1, logdet = TRUE, mat_logdet = NULL),
    "`logdet = TRUE`, but `mat_logdet` is NULL"
  )
})
