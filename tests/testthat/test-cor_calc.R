test_that("cor_calc() works with treating variables as continuous", {
  result <- cor_calc(data = mantar_dummy_full_cat, ordered = FALSE)
  expect_type(result, "list")
  expect_true(all(result$cor_method[lower.tri(result$cor_method)] == "pearson"))
})

test_that("cor_calc() works with treating variables as ordered categorical", {
  result <- cor_calc(data = mantar_dummy_full_cat, ordered = rep(T, 8))
  expect_type(result, "list")
  expect_true(is.numeric(result$mat))
  expect_true("matrix" %in% class(result$mat))
  expect_true(all(result$cor_method[lower.tri(result$cor_method)] == "polychoric"))
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
    'If `ordered` is a single character, it must be "adapted".',
    fixed = TRUE
  )
})

test_that("cor_calc() works with mixed variable types", {
  result <- cor_calc(data = mantar_dummy_full, ordered = c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE))
  expect_type(result, "list")
  expect_true(is.numeric(result$mat))
  expect_true("matrix" %in% class(result$mat))
  expect_true(all(c("pearson", "polychoric", "polyserial") %in% result$cor_method[lower.tri(result$cor_method)]))
})

test_that("cor_calc() works with adapted variable types", {
  expect_message(result <- cor_calc(data = mantar_dummy_full_cat, ordered = "adapted"))
  expect_type(result, "list")
  expect_true(is.numeric(result$mat))
  expect_true("matrix" %in% class(result$mat))
  expect_true(all(c("polychoric", "polyserial") %in% result$cor_method[lower.tri(result$cor_method)]))
})
