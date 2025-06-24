test_that("neighborhood_net() works with complete data and default arguments", {

  result <- neighborhood_net(data = mantar_dummy_full)

  # Check return type
  expect_type(result, "list")
  expect_named(result, c("pcor", "betas", "ns", "args"))

  # Check matrix dimensions
  expect_true(is.matrix(result$pcor))
  expect_true(is.matrix(result$betas))
  expect_equal(ncol(result$pcor), ncol(mantar_dummy_full))
  expect_equal(nrow(result$pcor), ncol(mantar_dummy_full))
})




test_that("errors in neighborhood network selection work for mat input", {
  # Test with incorrect mat type
  expect_error(neighborhood_net(
    mat =  matrix(c("a", "b", "b", "a"), nrow = 2),  # non-numeric matrix
    k = "log(n)"
  ), "All entries in 'mat' must be numeric")

  expect_error(neighborhood_net(
    mat =  matrix(c(1, 2, NA, 3), nrow = 2),  # non-numeric matrix
    k = "log(n)"
  ), "'mat' must be a symmetric matrix.")

  expect_error(neighborhood_net(
    mat =  matrix(c(1, NA, NA, 1), nrow = 2),  # non-numeric matrix
    k = "log(n)"
  ), "'mat' must not contain missing values.")
})

test_that("errors in neighborhood network selection work for data input", {
  expect_error(neighborhood_net(
    data =  cbind(mantar_dummy_full, rep("a", nrow(mantar_dummy_full))),  # non-numeric matrix
    k = "log(n)"
  ), "All variables in 'data' must be numeric.")

  expect_error(neighborhood_net(
    data =  cbind(mantar_dummy_full, rep("a", nrow(mantar_dummy_full))),  # non-numeric matrix
    mat = matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow =3),
    k = "log(n)"
  ), "All variables in 'data' must be numeric.")
})
