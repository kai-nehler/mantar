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
