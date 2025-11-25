test_that("ordered_suggest() treats continuous data as non-ordered and gives appropriate message", {
  expect_message(
    res <- ordered_suggest(mantar_dummy_full_cont, max_categories = 7),
    "Using the 'pearson' correlation method is advised because the number of distinct values suggests continuous variables.",
    fixed = TRUE
  )

  # alle Variablen als kontinuierlich (nicht ordered)
  expect_false(any(res))

  # Länge passt
  expect_equal(length(res), ncol(mantar_dummy_full_cont))

  # Namen passen zu den Spaltennamen
  expect_identical(names(res), colnames(mantar_dummy_full_cont))
})


test_that("ordered_suggest() treats ordered categorical data as ordered and gives appropriate message", {
  expect_message(
    res <- ordered_suggest(mantar_dummy_full_cat, max_categories = 7),
    "Using only 'polychoric' correlations is advised because all variables are ordered and there is sufficient information.",
    fixed = TRUE
  )

  # alle Variablen als ordered
  expect_true(all(res))

  # Länge passt
  expect_equal(length(res), ncol(mantar_dummy_full_cat))

  # Namen passen zu den Spaltennamen
  expect_identical(names(res), colnames(mantar_dummy_full_cat))
})


test_that("ordered_suggest() treats mixed data correctly and gives appropriate message", {
  expect_message(
    res <- ordered_suggest(mantar_dummy_full_mix, max_categories = 7),
    "Using a mix of 'pearson', 'polychoric', and 'polyserial' correlations is advised because some variables are ordered but not all, and there is sufficient information. Ordered variables are StressSens, SelfAware, ThoughtFuture, RespCriticism",
    fixed = TRUE
  )

  # alle Variablen als ordered
  expect_true(class(res) == "logical")

  # Länge passt
  expect_equal(length(res), ncol(mantar_dummy_full_mix))

  # Namen passen zu den Spaltennamen
  expect_identical(names(res), colnames(mantar_dummy_full_mix))
})
