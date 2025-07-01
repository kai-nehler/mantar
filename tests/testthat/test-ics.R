test_that("sample_size_full_works", {
  expect_equal(calculate_sample_size(data = mantar_dummy_full), rep(400, 8))
})

test_that("sample_size_mis_works", {
  expect_equal(calculate_sample_size(data = mantar_dummy_mis, n_calc = "Total"), rep(600, 8))
})

test_that("error_for_different_ncalc", {
  expect_error(
    calculate_sample_size(data = mantar_dummy_mis, n_calc = "abc"),
    class = "error"
  )
})

test_that("ic_calculated as number",{
  expect_equal(
    class(reg_ic_calc(resid_var = 1, n = 400, n_preds = 2, k = "log(n)")),
    "numeric"
  )
})
