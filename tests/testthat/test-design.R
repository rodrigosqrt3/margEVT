test_that("intercept-only matrices are correct", {
  df <- data.frame(x = 1:5, y = rnorm(5))
  dm <- build_design_matrices(df)
  expect_equal(ncol(dm$X_mu),    1L)
  expect_equal(ncol(dm$X_sigma), 1L)
  expect_equal(ncol(dm$X_xi),    1L)
  expect_equal(length(dm$idx_pen_mu),    0L)
  expect_equal(length(dm$idx_pen_sigma), 0L)
  expect_equal(length(dm$idx_pen_xi),    0L)
})

test_that("covariates appear in correct columns", {
  df <- data.frame(a = 1:5, b = rnorm(5), c = rnorm(5))
  dm <- build_design_matrices(df, loc_vars = c("a", "b"), shape_vars = "c")
  expect_equal(colnames(dm$X_mu),    c("(Intercept)", "a", "b"))
  expect_equal(colnames(dm$X_xi),    c("(Intercept)", "c"))
  expect_equal(dm$idx_pen_mu,        c(2L, 3L))
  expect_equal(dm$idx_pen_xi,        2L)
})

test_that("free_vars are excluded from penalized indices", {
  df <- data.frame(cos1 = rnorm(5), x = rnorm(5))
  dm <- build_design_matrices(df,
                              loc_vars  = c("cos1", "x"),
                              free_vars = "cos1")
  # cos1 is column 2, x is column 3 — only x should be penalized
  expect_equal(dm$idx_pen_mu, 3L)
})

test_that("missing column throws informative error", {
  df <- data.frame(a = 1:5)
  expect_error(
    build_design_matrices(df, loc_vars = c("a", "z")),
    regexp = "columns not found"
  )
})
