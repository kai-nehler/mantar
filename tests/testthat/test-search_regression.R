test_that("regression_opt() works with full data and AIC", {
  result <- regression_opt(
    data = mantar_dummy_full,
    dep_ind = 1,
    k = "2"
  )

  expect_type(result, "list")
  expect_named(result, c("regression", "R2"))
  expect_type(result$regression, "double")
  expect_type(result$R2, "double")
  expect_true(length(result$regression) >= 0)
})

test_that("regression_opt() works with data with two-step EM missingness and BIC", {
  result_mis <- regression_opt(
   data = mantar_dummy_mis,
   dep_ind = 2,
   n_calc = "individual",
   missing_handling = "two-step-em",
   )

  expect_type(result_mis, "list")
  expect_named(result_mis, c("regression", "R2"))
  expect_type(result_mis$regression, "double")
  expect_type(result_mis$R2, "double")
})


test_that("regression_opt() works with data with stacked MI missingness and BIC", {
  result_mis <- regression_opt(
    data = mantar_dummy_mis,
    dep_ind = 2,
    n_calc = "individual",
    missing_handling = "stacked-mi",
  )

  expect_type(result_mis, "list")
  expect_named(result_mis, c("regression", "R2"))
  expect_type(result_mis$regression, "double")
  expect_type(result_mis$R2, "double")
})

test_that("pred search works for correct specification", {
  mat <- cor(mantar_dummy_full)
  n <- nrow(mantar_dummy_full)  # sample size for IC calculation

  # Set dependent variable to be the first column
  dep_ind <- 1
  possible_pred_ind <- setdiff(seq_len(ncol(mat)), dep_ind)

  # Run pred_search
  res <- pred_search(
    mat = mat,
    dep_ind = dep_ind,
    possible_pred_ind = possible_pred_ind,
    n = n,
    k = "log(n)"  # BIC penalty
  )

  # Basic structure checks
  expect_type(res, "list")
  expect_named(res, c("actual_preds", "actual_betas", "actual_resid_var", "best_IC"), ignore.order = TRUE)

  # Predictors and betas must match in length
  expect_equal(length(res$actual_preds), length(res$actual_betas))

  # IC must be numeric
  expect_type(res$best_IC, "double")
})


test_that("errors in regression optimization work for mat input", {
  # Test with incorrect mat type
  expect_error(regression_opt(
    mat =  matrix(c("a", "b", "b", "a"), nrow = 2),  # non-numeric matrix
    dep_ind = 1,
    n = 100,
    k = "log(n)"
  ), "All entries in 'mat' must be numeric")

  expect_error(regression_opt(
    mat =  matrix(c(1, 2, NA, 3), nrow = 2),  # non-numeric matrix
    dep_ind = 1,
    n = 100,
    k = "log(n)"
  ), "'mat' must be a symmetric matrix.")

  expect_error(regression_opt(
    mat =  matrix(c(1, NA, NA, 1), nrow = 2),  # non-numeric matrix
    dep_ind = 1,
    n = 100,
    k = "log(n)"
  ), "'mat' must not contain missing values.")
})

test_that("errors in regression optimization work for data input", {
  expect_error(regression_opt(
    data =  cbind(mantar_dummy_full, rep("a", nrow(mantar_dummy_full))),  # non-numeric matrix
    dep_ind = 1,
    n = 100,
    k = "log(n)"
  ), "All variables in 'data' must be numeric.")

  expect_error(regression_opt(
    data =  cbind(mantar_dummy_full, rep("a", nrow(mantar_dummy_full))),  # non-numeric matrix
    dep_ind = 1,
    mat = matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow =3),
    n = 100,
    k = "log(n)"
  ), "All variables in 'data' must be numeric.")
})
