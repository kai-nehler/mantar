#### Tests for cor_calc() ####

test_that("cor_calc() works with treating variables as continuous", {
  res <- cor_calc(data = mantar_dummy_full_cat, ordered = FALSE)
  expect_type(res, "list")
  expect_true(all(res$cor_method[lower.tri(res$cor_method)] == "pearson"))
  expect_equal(res$means, colMeans(mantar_dummy_full_cat))
  expect_null(res$args$missing_handling)
  expect_null(res$args$nimp)
})

test_that("cor_calc() works with treating variables as ordered categorical", {
  result <- cor_calc(data = mantar_dummy_full_cat, ordered = rep(T, 8))
  expect_type(result, "list")
  expect_true(is.numeric(result$mat))
  expect_true("matrix" %in% class(result$mat))
  expect_true(all(result$cor_method[lower.tri(result$cor_method)] == "polychoric"))
})

test_that("cor_calc() recycles ordered if length 1", {
  res <- cor_calc(data = mantar_dummy_full_cat, ordered = TRUE)
  expect_true(all(res$cor_method[lower.tri(res$cor_method)] == "polychoric"))
})


test_that("errors in cor_calc work", {
  expect_error(
    cor_calc(data = mantar_dummy_full_cat, ordered = c(FALSE, TRUE)),
    "If `ordered` has length > 1, its length must equal the number of variables (p = 8).",
    fixed = TRUE
  )

  expect_error(
    cor_calc(data = mantar_dummy_full_cat, ordered = NA),
    "`ordered` must not contain NA.",
    fixed = TRUE
  )

  expect_error(
    cor_calc(data = mantar_dummy_full_cat, ordered = "yes"),
    'All entries of `ordered` must be logical (TRUE/FALSE).',
    fixed = TRUE
  )
})

test_that("cor_calc() works with mixed variable types", {
  result <- cor_calc(data = mantar_dummy_full_mix, ordered = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE))
  expect_type(result, "list")
  expect_true(is.numeric(result$mat))
  expect_true("matrix" %in% class(result$mat))
  expect_true(all(c("pearson", "polychoric", "polyserial") %in% result$cor_method[lower.tri(result$cor_method)]))
})

test_that("cor_calc() errors for two-step-em with ordered variables", {
  dat <- mantar_dummy_full_cat
  dat[1, 1] <- NA

  expect_error(
    cor_calc(
      data = dat,
      ordered = TRUE,
      missing_handling = "two-step-em"
    ),
    "Cannot use 'two-step-em' when any variables are ordered",
    fixed = FALSE
  )
})

#### Tests for transform_inverse() ####
test_that("transform_inverse() returns identity for identity inverse covariance", {
  inv_cov <- diag(3)
  res <- transform_inverse(inv_cov)

  expect_true(is.matrix(res))
  expect_equal(res, diag(3))
})


test_that("transform_inverse() preserves dimensions and symmetry", {
  inv_cov <- matrix(c(
    2, 0.5, 0.1,
    0.5, 1, 0.3,
    0.1, 0.3, 0.7
  ), nrow = 3, byrow = TRUE)

  res <- transform_inverse(inv_cov)

  expect_equal(dim(res), dim(inv_cov))
  expect_true(all(abs(res - t(res)) < 1e-10))
})

#### Tests for inv_to_net() ####
test_that("inv_to_net() returns symmetric matrix with zero diagonal", {
  theta <- matrix(c(
    1,  0.3, -0.2,
    0.3, 1,   0.1,
    -0.2, 0.1, 1
  ), nrow = 3, byrow = TRUE)

  net <- inv_to_net(theta)

  expect_true(is.matrix(net) || methods::is(net, "Matrix"))
  expect_equal(dim(net), dim(theta))
  expect_true(all(diag(net) == 0))

  expect_true(all(abs(net - t(net)) < 1e-10))
})
