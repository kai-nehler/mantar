#### Tests for regression_opt() ----

test_that("regression_opt() works with full data and AIC", {
  result <- regression_opt(
    data = mantar_dummy_full_cont,
    dep_ind = 1,
    ic_type = "aic"
  )

  expect_type(result, "list")
  expect_named(result, c("regression", "R2", "n", "args"))
  expect_type(result$regression, "double")
  expect_type(result$R2, "double")
  expect_true(length(result$regression) >= 0)
})

test_that("regression_opt() works with data with two-step EM missingness and BIC", {
  result_mis <- regression_opt(
   data = mantar_dummy_mis_cont,
   dep_ind = 2,
   n_calc = "individual",
   missing_handling = "two-step-em",
   )

  expect_type(result_mis, "list")
  expect_named(result_mis, c("regression", "R2", "n", "args"))
  expect_type(result_mis$regression, "double")
  expect_type(result_mis$R2, "double")
})


test_that("regression_opt() works with data with stacked MI missingness and BIC", {
  result_mis <- regression_opt(
    data = mantar_dummy_mis_cont,
    dep_ind = 2,
    n_calc = "individual",
    missing_handling = "stacked-mi",
  )

  expect_type(result_mis, "list")
  expect_named(result_mis, c("regression", "R2", "n", "args"))
  expect_type(result_mis$regression, "double")
  expect_type(result_mis$R2, "double")
})


test_that("regression_opt() works when mat and n are supplied", {
  mat <- stats::cov(mantar_dummy_full_cont)
  n   <- nrow(mantar_dummy_full_cont)

  res <- regression_opt(
    mat = mat,
    n   = n,
    dep_ind = 1
  )

  expect_s3_class(res, "mantar_regression")
  expect_named(res, c("regression", "R2", "n", "args"))
  expect_equal(res$n, n)

  expect_type(res$regression, "double")
  expect_true(length(res$regression) >= 0)

  # check args
  expect_equal(res$args$ic_type, "bic")
  expect_null(res$args$cor_method)
  expect_null(res$args$missing_handling)
  expect_null(res$args$nimp)
  expect_null(res$args$imp_method)
})




test_that("regression_opt() accepts dep_ind as column name", {
  dep_name <- colnames(mantar_dummy_full_cont)[1]

  res_index <- regression_opt(
    data   = mantar_dummy_full_cont,
    dep_ind = 1,
    ic_type = "aic"
  )

  res_name <- regression_opt(
    data   = mantar_dummy_full_cont,
    dep_ind = dep_name,
    ic_type = "aic"
  )

  # gleiche Regressionsergebnisse erwartet
  expect_equal(res_index$regression, res_name$regression)
  expect_equal(res_index$R2, res_name$R2)
  expect_equal(res_index$n, res_name$n)
})

test_that("regression_opt() errors for invalid dep_ind specifications", {
  expect_error(
    regression_opt(
      data   = mantar_dummy_full_cont,
      dep_ind = "not_a_column"
    ),
    "dep_ind must be a valid column name.",
    fixed = TRUE
  )

  mat <- stats::cov(mantar_dummy_full_cont)
  n   <- nrow(mantar_dummy_full_cont)

  expect_error(
    regression_opt(
      mat = mat,
      n   = n,
      dep_ind = 0
    ),
    "dep_ind must be a valid column index.",
    fixed = TRUE
  )
})

test_that("errors in regression optimization work for mat input", {
  # Test with incorrect mat type
  expect_error(regression_opt(
    mat =  matrix(c("a", "b", "b", "a"), nrow = 2),  # non-numeric matrix
    dep_ind = 1,
    n = 100,
    ic_type = "bic"
  ), "All entries in 'mat' must be numeric")

  expect_error(regression_opt(
    mat =  matrix(c(1, 2, NA, 3), nrow = 2),  # non-numeric matrix
    dep_ind = 1,
    n = 100,
    ic_type = "bic"
  ), "'mat' must be a symmetric matrix.")

  expect_error(regression_opt(
    mat =  matrix(c(1, NA, NA, 1), nrow = 2),  # non-numeric matrix
    dep_ind = 1,
    n = 100,
    ic_type = "bic"
  ), "'mat' must not contain missing values.")
})

test_that("errors in regression optimization work for data input", {
  expect_error(regression_opt(
    data =  cbind(mantar_dummy_full_cont, rep("a", nrow(mantar_dummy_full_cont))),  # non-numeric matrix
    dep_ind = 1,
    n = 100,
    ic_type = "bic"
  ), "All variables in 'data' must be numeric.")

  expect_error(regression_opt(
    data =  cbind(mantar_dummy_full_cont, rep("a", nrow(mantar_dummy_full_cont))),  # non-numeric matrix
    dep_ind = 1,
    mat = matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow =3),
    n = 100,
    ic_type = "bic"
  ), "All variables in 'data' must be numeric.")
})


#### Tests for pred_search() ----
test_that("pred search works for correct specification", {
  mat <- cor(mantar_dummy_full_cont)
  n <- nrow(mantar_dummy_full_cont)  # sample size for IC calculation

  # Set dependent variable to be the first column
  dep_ind <- 1
  possible_pred_ind <- setdiff(seq_len(ncol(mat)), dep_ind)

  # Run pred_search
  res <- pred_search(
    mat = mat,
    dep_ind = dep_ind,
    possible_pred_ind = possible_pred_ind,
    n = n,
    ic_type = "bic"  # BIC penalty
  )

  # Basic structure checks
  expect_type(res, "list")
  expect_named(res, c("actual_preds", "actual_betas", "actual_resid_var", "best_IC"), ignore.order = TRUE)

  # Predictors and betas must match in length
  expect_equal(length(res$actual_preds), length(res$actual_betas))

  # IC must be numeric
  expect_type(res$best_IC, "double")
})


test_that("pred_search() keeps null model when predictors do not improve IC", {
  mat <- diag(3)
  n   <- 100

  res <- pred_search(
    mat = mat,
    dep_ind = 1,
    possible_pred_ind = 2:3,
    n = n,
    ic_type = "bic"
  )

  # No predictors in the model
  expect_equal(length(res$actual_preds), 0)
  expect_equal(length(res$actual_betas), 0)

  # Residual variance should be 1 (since dep var is standardized)
  expect_equal(res$actual_resid_var, 1)

  # IC shold be of type double
  expect_type(res$best_IC, "double")
})
