example_matrix <- matrix(c(1, 0.3, 0.5,
                            0.3, 1, 0.2,
                            0.5, 0.2, 1), nrow = 3, ncol = 3)

test_that("matrix regression returns correct dimensions for beta weights", {
  expect_equal(length(matrix_regression(mat = example_matrix, dep_ind = 1, pred_ind = c(2,3))$beta), 2)
  expect_equal(length(matrix_regression(mat = example_matrix, dep_ind = 1, pred_ind = NULL)$beta), 0)
})

test_that("matrix regression returns correct dimensions for residual variance", {
  expect_equal(length(matrix_regression(mat = example_matrix, dep_ind = 1, pred_ind = c(2,3))$resid_var), 1)
  expect_equal(length(matrix_regression(mat = example_matrix, dep_ind = 1, pred_ind = NULL)$resid_var), 1)
})

test_that("matrix_regression() returns a list with expected elements and types", {
  res <- matrix_regression(mat = example_matrix, dep_ind = 1, pred_ind = c(2, 3))

  expect_type(res, "list")
  expect_named(res, c("beta", "resid_var"))

  expect_true(is.matrix(res$beta))
  expect_true(is.numeric(res$resid_var))
  expect_length(res$resid_var, 1)
})


test_that("matrix_regression() with no predictors returns variance of dependent variable", {
  res <- matrix_regression(mat = example_matrix, dep_ind = 2, pred_ind = integer(0))

  expect_null(res$beta)
  expect_equal(res$resid_var, example_matrix[2, 2])
})

test_that("matrix_regression() errors if dep_ind has length != 1", {
  expect_error(
    matrix_regression(mat = example_matrix, dep_ind = c(1, 2), pred_ind = 3),
    "dep_ind must have exactly one index",
    fixed = TRUE
  )

  expect_error(
    matrix_regression(mat = example_matrix, dep_ind = integer(0), pred_ind = 1),
    "dep_ind must have exactly one index",
    fixed = TRUE
  )
})

