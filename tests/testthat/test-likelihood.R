# Helper: minimal stationary setup (no covariates, no penalty)
make_stationary_setup <- function(n = 200, seed = 1L) {
  set.seed(seed)
  y  <- c(rexp(n), runif(20, 5, 10))   # some exceedances above u=4
  df <- data.frame(y = y)
  dm <- build_design_matrices(df)
  list(y = y, dm = dm,
       par = c(mu.Intercept = 6, sigma.Intercept = log(1.5), xi.Intercept = 0.1))
}

test_that("pp_nllh returns a finite scalar", {
  s <- make_stationary_setup()
  val <- pp_nllh(s$par, s$dm, s$y, threshold = 4,
                 lambda = 0, obs_per_year = 365.25)
  expect_true(is.finite(val))
  expect_length(val, 1L)
})

test_that("pp_nllh returns 1e9 when support is violated", {
  s   <- make_stationary_setup()
  bad <- s$par
  bad["xi.Intercept"] <- 10   # forces z_u <= 0 for many obs
  val <- pp_nllh(bad, s$dm, s$y, threshold = 4,
                 lambda = 0, obs_per_year = 365.25)
  expect_equal(val, 1e9)
})

test_that("pp_grad has correct length", {
  s <- make_stationary_setup()
  g <- pp_grad(s$par, s$dm, s$y, threshold = 4,
               lambda = 0, obs_per_year = 365.25)
  expect_length(g, length(s$par))
})

test_that("pp_grad is numerically close to finite differences", {
  s   <- make_stationary_setup()
  eps <- 1e-6
  g_an <- pp_grad(s$par, s$dm, s$y, threshold = 4,
                  lambda = 0, obs_per_year = 365.25)
  g_fd <- numericDeriv_simple <- sapply(seq_along(s$par), function(i) {
    p1 <- p2 <- s$par
    p1[i] <- p1[i] + eps
    p2[i] <- p2[i] - eps
    (pp_nllh(p1, s$dm, s$y, 4, lambda = 0, obs_per_year = 365.25) -
        pp_nllh(p2, s$dm, s$y, 4, lambda = 0, obs_per_year = 365.25)) / (2 * eps)
  })
  expect_equal(g_an, g_fd, tolerance = 1e-4)
})

test_that("penalty increases nllh when lambda > 0", {
  s    <- make_stationary_setup()
  dm_c <- build_design_matrices(
    data.frame(y = s$y, x = rnorm(length(s$y))),
    loc_vars = "x"
  )
  par_c <- c(mu.Intercept = 6, mu.x = 0.5,
             sigma.Intercept = log(1.5), xi.Intercept = 0.1)
  v0 <- pp_nllh(par_c, dm_c, s$y, 4, lambda = 0)
  v1 <- pp_nllh(par_c, dm_c, s$y, 4, lambda = 1)
  expect_gt(v1, v0)
})
