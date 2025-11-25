#### Tests for neighborhood_net() ####

test_that("neighborhood_net() works with complete data and default arguments", {

  result <- neighborhood_net(data = mantar_dummy_full_cont)

  # Check return type
  expect_type(result, "list")
  expect_named(result, c("pcor", "betas", "ns", "args"))

  # Check matrix dimensions
  expect_true(is.matrix(result$pcor))
  expect_true(is.matrix(result$betas))
  expect_equal(ncol(result$pcor), ncol(mantar_dummy_full_cont))
  expect_equal(nrow(result$pcor), ncol(mantar_dummy_full_cont))
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
    data =  cbind(mantar_dummy_full_cont, rep("a", nrow(mantar_dummy_full_cont))),  # non-numeric matrix
    k = "log(n)"
  ), "All variables in 'data' must be numeric.")

  expect_error(neighborhood_net(
    data =  cbind(mantar_dummy_full_cont, rep("a", nrow(mantar_dummy_full_cont))),  # non-numeric matrix
    mat = matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow =3),
    k = "log(n)"
  ), "All variables in 'data' must be numeric.")
})


test_that("neighborhood_net() works with mat + ns input", {
  mat <- stats::cov(mantar_dummy_full_cont)
  ns  <- nrow(mantar_dummy_full_cont)

  result <- neighborhood_net(
    mat = mat,
    ns  = ns
  )

  # strcuture
  expect_type(result, "list")
  expect_named(result, c("pcor", "betas", "ns", "args"))

  # classes
  expect_s3_class(result, "mantar_neighborhood")
  expect_s3_class(result, "mantar_network")

  # dimensions
  p <- ncol(mantar_dummy_full_cont)
  expect_true(is.matrix(result$pcor))
  expect_true(is.matrix(result$betas))
  expect_equal(dim(result$pcor), c(p, p))
  expect_equal(dim(result$betas), c(p, p))

  # ns-Recycling
  expect_equal(result$ns, rep(ns, p))

  # args: missing infos should be NULL
  expect_equal(result$args$k, "log(n)")
  expect_null(result$args$cor_method)
  expect_null(result$args$missing_handling)
  expect_null(result$args$nimp)
  expect_null(result$args$imp_method)
})


test_that("neighborhood_net() recycles scalar ns for data input", {
  ns_scalar <- 50

  result <- neighborhood_net(
    data = mantar_dummy_full_cont,
    ns   = ns_scalar
  )

  p <- ncol(mantar_dummy_full_cont)
  expect_equal(result$ns, rep(ns_scalar, p))
})



#### Tests for neihborhood_sel() ####
test_that("neighborhood_sel() returns zero network for identity correlation matrix", {
  mat <- diag(3)
  colnames(mat) <- paste0("V", 1:3)
  ns  <- rep(100, 3)
  k   <- log(100)

  res <- neighborhood_sel(mat = mat, ns = ns, k = k, pcor_merge_rule = "and")

  expect_type(res, "list")
  expect_named(res, c("partials", "beta_mat"))

  # dimensions and names
  expect_equal(dim(res$partials), c(3, 3))
  expect_equal(dim(res$beta_mat), c(3, 3))
  expect_equal(colnames(res$partials), colnames(mat))
  expect_equal(rownames(res$partials), colnames(mat))

  # symmetry & zeros
  expect_true(all(abs(res$partials - t(res$partials)) < 1e-10))
  expect_true(all(res$partials == 0))
})


#### Tests for compute_partials() ####

test_that("compute_partials() with 'and' rule uses only edges with non-zero betas in both directions and same sign", {
  betas <- matrix(0, nrow = 3, ncol = 3)
  colnames(betas) <- rownames(betas) <- c("A", "B", "C")

  # A <-> B: both positive
  betas["A", "B"] <- 0.5
  betas["B", "A"] <- 0.4

  # A <-> C: only one direction
  betas["A", "C"] <- 0.3
  betas["C", "A"] <- 0

  # B <-> C: different signs
  betas["B", "C"] <- 0.2
  betas["C", "B"] <- -0.1

  res <- compute_partials(betas, rule = "and")

  wadj <- res$wadj
  adj  <- res$adj

  # A-B: sqrt(0.5*0.4) > 0
  expect_gt(wadj["A", "B"], 0)
  expect_equal(wadj["A", "B"], wadj["B", "A"])
  expect_equal(adj["A", "B"], 1)
  expect_equal(adj["B", "A"], 1)

  # A-C: only one direction → 0
  expect_equal(wadj["A", "C"], 0)
  expect_equal(wadj["C", "A"], 0)
  expect_equal(adj["A", "C"], 0)
  expect_equal(adj["C", "A"], 0)

  # B-C: different signs
  expect_equal(wadj["B", "C"], 0)
  expect_equal(wadj["C", "B"], 0)
  expect_equal(adj["B", "C"], 0)
  expect_equal(adj["C", "B"], 0)

  # Diagonal = 0
  expect_true(all(diag(wadj) == 0))
})
