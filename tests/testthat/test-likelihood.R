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

test_that("fast likelihood path agrees with the reference implementation", {
  set.seed(3L)
  n <- 600L
  df <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  df$y <- 20 + 2 * df$x1 + exp(1 + 0.2 * df$x2) *
    (((50 * -log(stats::runif(n)))^0.1 - 1) / -0.1)
  dm <- build_design_matrices(df, loc_vars = "x1", scale_vars = "x2")
  u  <- unname(stats::quantile(df$y, 0.9))
  fit <- fit_nhpp(df, threshold = u, loc_vars = "x1", scale_vars = "x2",
                  penalty = "none", obs_per_year = 50, verbose = FALSE)
  set.seed(9L)
  n_inside <- 0L
  for (i in seq_len(300L)) {
    par <- fit$par + stats::rnorm(length(fit$par), sd = c(1, 0.5, 0.3, 0.2, 0.3))
    if (i %% 50L == 0L) par[length(par)] <- 1e-8        # Gumbel limit
    a <- margEVT:::.nhpp_eval_ref(par, dm, df$y, u, 50, TRUE, -0.5)
    b <- margEVT:::.nhpp_eval(par, dm, df$y, u, 50, TRUE, -0.5)
    v <- margEVT:::.nhpp_eval(par, dm, df$y, u, 50, FALSE, -0.5)
    expect_equal(b$value, a$value, tolerance = 1e-10)
    expect_equal(v$value, a$value, tolerance = 1e-10)
    expect_equal(b$grad, a$grad, tolerance = 1e-9)
    if (a$value < 1e9) n_inside <- n_inside + 1L
  }
  expect_gt(n_inside, 50L)
})

# -----------------------------------------------------------------------------
# Coordinate penalties, objective and KKT conditions
# -----------------------------------------------------------------------------

test_that(".resolve_block handles scalar, named and unnamed vectors", {
  expect_equal(margEVT:::.resolve_block(2, "sigma"), 2)
  expect_equal(margEVT:::.resolve_block(c(xi = 3, mu = 1, sigma = 2), "sigma"), 2)
  expect_equal(margEVT:::.resolve_block(c(1, 2, 3), "xi"), 3)
})

test_that(".penalty_vectors assigns block lambdas and weights by coordinate", {
  df <- make_pen_df(200L)
  dm <- build_design_matrices(df, loc_vars = "x1", scale_vars = "x2",
                              shape_vars = "x3")
  pv <- margEVT:::.penalty_vectors(dm, c(mu = 1, sigma = 2, xi = 3), 1, TRUE)
  expect_equal(pv$lam, c(0, 1, 0, 2, 0, 3))
  pv2 <- margEVT:::.penalty_vectors(dm, 1, 0.5, FALSE,
                                    penalty_weights = c(9, 2, 9, 3, 9, 4))
  expect_equal(pv2$lam, c(0, 2, 0, 3, 0, 0))
  expect_equal(pv2$alp, c(0, 0.5, 0, 0.5, 0, 0))
})

test_that("weighted pp_nllh and pp_grad agree with finite differences", {
  df <- make_pen_df(400L)
  dm <- build_design_matrices(df, loc_vars = c("x1", "x2"))
  par <- c(15, 1, 0.2, log(3), 0.05)
  w <- c(0, 2, 0.5, 0, 0)
  g <- pp_grad(par, dm, df$y, 16, lambda = 1.5, alpha = 0.7,
               penalty_weights = w)
  fd <- vapply(seq_along(par), function(k) {
    h <- 1e-6; a <- b <- par; a[k] <- a[k] + h; b[k] <- b[k] - h
    (pp_nllh(a, dm, df$y, 16, 1.5, 0.7, penalty_weights = w) -
       pp_nllh(b, dm, df$y, 16, 1.5, 0.7, penalty_weights = w)) / (2 * h)
  }, numeric(1L))
  expect_equal(g, fd, tolerance = 1e-4)
})

test_that(".objective_factory caches and flags the support barrier", {
  df <- make_pen_df(300L)
  dm <- build_design_matrices(df)
  pv <- margEVT:::.penalty_vectors(dm, 0, 1, TRUE)
  obj <- margEVT:::.objective_factory(dm, df$y, 16, pv)
  par <- c(15, log(3), 0.05)
  expect_identical(obj$fn(par), obj$fn(par))
  expect_length(obj$gr(par), 3L)
  expect_equal(obj$fn(c(15, log(3), -5)), 1e9)
})

test_that(".kkt_violation uses the subgradient at zero and rejects non-finite input", {
  pv <- list(lam = c(0, 2), alp = c(0, 1))
  expect_equal(margEVT:::.kkt_violation(c(1, 0), c(0.1, 1.5), pv, 1e-2), 0.1)
  expect_equal(margEVT:::.kkt_violation(c(1, 0), c(0, 3), pv, 1e-2), 1)
  expect_identical(margEVT:::.kkt_violation(c(1, 0), c(NA, 0), pv, 1e-2), Inf)
})

# -----------------------------------------------------------------------------
# [.]_+ convention of the integrated intensity and the shape bound
# -----------------------------------------------------------------------------

make_endpoint_setup <- function(n = 500L, seed = 11L) {
  set.seed(seed)
  x <- stats::rnorm(n)
  y <- rep(0, n)
  y[x > 1] <- 11.2                     # exceedances where the endpoint is high
  df <- data.frame(y = y, x = x)
  list(df = df, dm = build_design_matrices(df, loc_vars = "x"),
       par = c(10, 3, log(1), -0.4), u = 11)
}

test_that("periods with upper endpoint below u contribute zero to Lambda", {
  s <- make_endpoint_setup()
  mu <- 10 + 3 * s$df$x; sig <- 1; xi <- -0.4
  z_u <- 1 + xi * (s$u - mu) / sig
  expect_true(any(z_u <= 0))           # the case the old likelihood rejected
  tau <- ifelse(z_u > 0, pmax(z_u, 0)^(-1 / xi), 0)
  exc <- s$df$y > s$u
  z_e <- 1 + xi * (s$df$y[exc] - mu[exc]) / sig
  manual <- sum(tau) / 365.25 - sum(-log(sig) - (1 + 1 / xi) * log(z_e))
  val <- pp_nllh(s$par, s$dm, s$df$y, s$u, lambda = 0)
  expect_equal(val, manual, tolerance = 1e-10)
})

test_that("gradient matches finite differences when some endpoints lie below u", {
  s <- make_endpoint_setup()
  g <- pp_grad(s$par, s$dm, s$df$y, s$u, lambda = 0)
  fd <- vapply(seq_along(s$par), function(k) {
    h <- 1e-6; a <- b <- s$par; a[k] <- a[k] + h; b[k] <- b[k] - h
    (pp_nllh(a, s$dm, s$df$y, s$u, lambda = 0) -
       pp_nllh(b, s$dm, s$df$y, s$u, lambda = 0)) / (2 * h)
  }, numeric(1L))
  expect_equal(g, fd, tolerance = 1e-4)
})

test_that("support violations and the shape bound still return 1e9", {
  s <- make_endpoint_setup()
  pos <- s$par; pos[4L] <- 0.4; pos[1L] <- 30   # u below the lower endpoint
  expect_equal(pp_nllh(pos, s$dm, s$df$y, s$u, lambda = 0), 1e9)
  beyond <- s$par; beyond[1L] <- 5; beyond[2L] <- 0   # exceedance beyond endpoint
  expect_equal(pp_nllh(beyond, s$dm, s$df$y, s$u, lambda = 0), 1e9)
  expect_equal(pp_nllh(s$par, s$dm, s$df$y, s$u, lambda = 0,
                       shape_lower = -0.3), 1e9)
  expect_true(is.finite(pp_nllh(s$par, s$dm, s$df$y, s$u, lambda = 0,
                                shape_lower = -0.5)))
})

test_that(".shape_bound_map handles intercept, one covariate, and other designs", {
  df <- data.frame(y = 1:10, t = (1:10) / 10, z = stats::rnorm(10))
  dm0 <- build_design_matrices(df)
  m0 <- margEVT:::.shape_bound_map(dm0, -0.5)
  expect_equal(unname(m0$A), matrix(1))
  expect_equal(m0$idx, 3L)
  dm1 <- build_design_matrices(df, shape_vars = "t")
  m1 <- margEVT:::.shape_bound_map(dm1, -0.5)
  expect_equal(unname(m1$A), rbind(c(1, 0.1), c(1, 1)))
  expect_null(margEVT:::.shape_bound_map(dm1, -Inf))
  dm2 <- build_design_matrices(df, shape_vars = c("t", "z"))
  expect_null(margEVT:::.shape_bound_map(dm2, -0.5))
  # a constant shape covariate gives a singular map
  dmc <- build_design_matrices(data.frame(y = 1:10, t = rep(0.3, 10)),
                               shape_vars = "t")
  expect_null(margEVT:::.shape_bound_map(dmc, -0.5))
})

test_that(".kkt_violation_bounded uses multiplier signs at the bound", {
  df <- data.frame(y = 1:5)
  dm <- build_design_matrices(df)
  bmap <- margEVT:::.shape_bound_map(dm, -0.5)
  pv <- list(lam = c(0, 0, 0), alp = c(0, 0, 0))
  par <- c(1, 0, -0.5)
  # at the bound, a positive derivative is admissible, a negative one is not
  expect_equal(margEVT:::.kkt_violation_bounded(par, c(0, 0, 2), pv, 1e-2, bmap), 0)
  expect_equal(margEVT:::.kkt_violation_bounded(par, c(0, 0, -2), pv, 1e-2, bmap), 2)
  # interior shape: ordinary stationarity
  expect_equal(margEVT:::.kkt_violation_bounded(c(1, 0, -0.2), c(0, 0, 2), pv,
                                                1e-2, bmap), 2)
  # non-finite gradient
  expect_identical(margEVT:::.kkt_violation_bounded(par, c(0, NA, 0), pv,
                                                    1e-2, bmap), Inf)
})
