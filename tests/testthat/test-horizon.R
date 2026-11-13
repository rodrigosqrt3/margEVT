# Tests for risk over a planning horizon: horizon_risk() and
# design_life_level().

test_that("a one-year horizon reproduces the marginal return level", {
  s <- make_bt_df()
  set.seed(2L)
  mc <- lapply(1:200, function(i)
    s$df[s$df$year == sample(unique(s$df$year), 1L), "x", drop = FALSE])
  rl <- marginalize(s$fit, s$df, TRs = c(10, 50), approaches = "B",
                    mc_sample = mc)
  dll <- design_life_level(s$fit, lapply(mc, list), prob = c(0.1, 0.02))
  expect_equal(dll$DLL, rl$RL, tolerance = 1e-3)
  expect_equal(dll$se_mc, rl$se_mc, tolerance = 1e-3)
  expect_true(all(dll$horizon == 1L))
})

test_that("constant covariates give the closed forms of independent years", {
  s <- make_bt_df()
  block <- data.frame(x = rep(0.5, s$n_obs))
  L <- 15L
  traj <- rep(list(block), L)
  z <- c(12, 14)
  G <- annual_prob(s$fit, z, list(block))
  hr <- horizon_risk(s$fit, z, traj)
  expect_equal(hr$prob_exceed, 1 - G^L)
  expect_equal(hr$expected_years, L * (1 - G))
  expect_equal(hr$wait, (1 - G^L) / (1 - G))
  expect_true(all(is.na(hr$se_mc)))
  cum <- attr(hr, "cumulative")
  expect_identical(dim(cum), c(L, 2L))
  expect_equal(cum[L, ], hr$prob_exceed, ignore_attr = TRUE)
  expect_equal(unname(cum[1, ]), 1 - G)
  expect_identical(attr(hr, "horizon"), L)
  dll <- design_life_level(s$fit, traj, prob = 0.2)
  expect_equal(horizon_risk(s$fit, dll$DLL, traj)$prob_exceed, 0.2,
               tolerance = 1e-4)
})

test_that("rising covariates raise the risk; ensembles get a standard error", {
  s <- make_bt_df()
  flat <- rep(list(data.frame(x = rep(0, s$n_obs))), 10L)
  rise <- lapply(seq(0, 2, length.out = 10L),
                 function(v) data.frame(x = rep(v, s$n_obs)))
  expect_gt(horizon_risk(s$fit, 13, rise)$prob_exceed,
            horizon_risk(s$fit, 13, flat)$prob_exceed)
  expect_gt(design_life_level(s$fit, rise, prob = 0.1)$DLL,
            design_life_level(s$fit, flat, prob = 0.1)$DLL)
  set.seed(4L)
  ens <- replicate(50L, lapply(1:10, function(i)
    data.frame(x = stats::rnorm(s$n_obs))), simplify = FALSE)
  hr <- horizon_risk(s$fit, 13, ens)
  expect_true(is.finite(hr$se_mc) && hr$se_mc > 0)
  expect_true(is.finite(design_life_level(s$fit, ens, prob = 0.1)$se_mc))
})

test_that("fixed values and the default upper bound are used", {
  s <- make_bt_df()
  traj <- rep(list(data.frame(w = rep(0, s$n_obs))), 5L)
  a <- horizon_risk(s$fit, 13, traj, fixed = list(x = 1))
  b <- horizon_risk(s$fit, 13, rep(list(data.frame(x = rep(1, s$n_obs))), 5L))
  expect_equal(a$prob_exceed, b$prob_exceed)
  legacy <- s$fit
  legacy$y <- NULL
  expect_error(design_life_level(legacy, traj, fixed = list(x = 1)), "z_hi")
  expect_true(all(is.finite(design_life_level(legacy, traj, fixed = list(x = 1),
                                              z_hi = 60)$DLL)))
})

test_that("horizon functions validate their input", {
  s <- make_bt_df()
  one <- data.frame(x = rep(0, s$n_obs))
  expect_error(horizon_risk(list(), 13, list(one)), "nhpp_fit")
  expect_error(horizon_risk(s$fit, NA, list(one)), "finite")
  expect_error(horizon_risk(s$fit, 13, one), "trajectories")
  expect_error(horizon_risk(s$fit, 13, list(list(one), list(one, one))),
               "same number")
  expect_error(horizon_risk(s$fit, 13, list(data.frame(w = 1))),
               "missing required")
  expect_error(horizon_risk(s$fit, 13, list(one), n_obs = 0), "n_obs")
  expect_error(design_life_level(s$fit, list(one), prob = 1), "prob")
  expect_error(design_life_level(s$fit, list(one), z_hi = 1), "z_hi")
})
