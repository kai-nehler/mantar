test_that("print.mantar_network() prints pcor and returns invisibly", {
  pcor <- matrix(c(1, 0.3,
                   0.3, 1), nrow = 2)
  obj <- list(pcor = pcor)
  class(obj) <- "mantar_network"

  expect_output(
    out <- print(obj),
    "0.3"  # one entry in the matrix
  )
  expect_identical(out, obj)
})


test_that("summary.mantar_neighborhood() computes density and stores meta info", {
  pcor <- matrix(c(
    0, 0.5, 0,
    0.5, 0, 0.2,
    0, 0.2, 0
  ), nrow = 3, byrow = TRUE)
  colnames(pcor) <- c("A", "B", "C")

  obj <- list(
    pcor = pcor,
    args = list(
      ic_type = "bic",
      pcor_merge_rule = "and",
      missing_handling = "stacked-mi",
      nimp = 20
    ),
    ns = c(100, 110, 120)
  )
  class(obj) <- "mantar_neighborhood"

  s <- summary(obj)

  expect_s3_class(s, "summary.mantar_neighborhood")
  expect_named(s, c("density", "args", "ns", "varnames"))

  # density = number of edges / number of possible edges
  expect_equal(s$density, 2 / choose(3, 2))

  expect_equal(s$args, obj$args)
  expect_equal(s$ns, obj$ns)
  expect_equal(s$varnames, colnames(pcor))
})

test_that("print.summary.mantar_neighborhood() reports density and settings", {
  pcor <- matrix(c(
    0, 0.5,
    0.5, 0
  ), nrow = 2)
  colnames(pcor) <- c("A", "B")

  obj <- list(
    pcor = pcor,
    args = list(
      ic_type = "bic",
      pcor_merge_rule = "or",
      missing_handling = NULL,
      nimp = NULL
    ),
    ns = c(100, 120)
  )
  class(obj) <- "mantar_neighborhood"

  s <- summary(obj)

  expect_output(
    print(s),
    "density of the estimated network is",
  )
  expect_output(
    print(s),
    "the 'or' rule for the inclusion of edges",
    ignore.case = TRUE
  )
  expect_output(
    print(s),
    "The sample sizes used for the nodewise regressions were as follows:",
  )
})


test_that("summary.mantar_regularization() computes density and stores meta info", {
  pcor <- matrix(c(
    0, 0.4,
    0.4, 0
  ), nrow = 2)
  colnames(pcor) <- c("X1", "X2")

  obj <- list(
    pcor = pcor,
    args = list(
      penalty = "glasso",
      missing_handling = "stacked-mi",
      nimp = 10,
      likelihood = "obs_based"
    ),
    n = 150
  )
  class(obj) <- "mantar_regularization"

  s <- summary(obj)

  expect_s3_class(s, "summary.mantar_regularization")
  expect_equal(s$density, 1 / choose(2, 2)) # nur eine Kante bei 2 Knoten → 1/1
  expect_equal(s$varnames, colnames(pcor))
  expect_equal(s$n, 150)
  expect_equal(s$args$penalty, "glasso")
})

test_that("print.summary.mantar_regularization() reports likelihood, n and missing handling", {
  pcor <- diag(2)
  colnames(pcor) <- c("X1", "X2")

  obj <- list(
    pcor = pcor,
    args = list(
      penalty = "glasso",
      missing_handling = "two-step-em",
      nimp = NULL,
      likelihood = "mat_based"
    ),
    n = 200
  )
  class(obj) <- "mantar_regularization"

  s <- summary(obj)

  expect_output(
    print(s),
    "density of the estimated network is",
  )
  expect_output(
    print(s),
    "regularization with the glasso penalty",
    ignore.case = TRUE
  )
  expect_output(
    print(s),
    "Log-Likelihood calculation was based on the sample correlation matrix",
    ignore.case = TRUE
  )
  expect_output(
    print(s),
    "Missing data were handled using 'two-step-em'",
    ignore.case = TRUE
  )
})


test_that("plot.mantar_network() calls qgraph when available", {
  testthat::skip_if_not_installed("qgraph")

  pcor <- matrix(c(
    0, 0.2,
    0.2, 0
  ), nrow = 2)
  colnames(pcor) <- c("A", "B")

  obj <- list(pcor = pcor)
  class(obj) <- "mantar_network"

  expect_silent(
    plot(obj)
  )
})


test_that("plot.mantar_network() calls qgraph when available", {
  testthat::skip_if_not_installed("qgraph")

  pcor <- matrix(c(
    0, 0.2,
    0.2, 0
  ), nrow = 2)
  colnames(pcor) <- c("A", "B")

  obj <- list(pcor = pcor)
  class(obj) <- "mantar_network"

  expect_silent(
    plot(obj)  # plot
  )
})
