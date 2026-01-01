# =============================================================================
# Test Suite: Module 03 - Conditional Moments of K_J | alpha
# =============================================================================
#
# Author: JoonHo Lee
# Date: December 2025
# =============================================================================

# -----------------------------------------------------------------------------
# Test: Polygamma vs Summation Formula Equivalence
# -----------------------------------------------------------------------------

test_that("Conditional mean matches summation formula", {
  test_cases <- list(
    c(J = 50, alpha = 0.5),
    c(J = 50, alpha = 1.0),
    c(J = 50, alpha = 2.0),
    c(J = 50, alpha = 5.0),
    c(J = 100, alpha = 2.0),
    c(J = 300, alpha = 2.0)
  )

  for (tc in test_cases) {
    J <- tc["J"]
    alpha <- tc["alpha"]
    expect_equal(
      mean_K_given_alpha(J, alpha),
      mean_K_given_alpha_sum(J, alpha),
      tolerance = 1e-10,
      info = sprintf("J=%d, alpha=%.1f", J, alpha)
    )
  }
})

test_that("Conditional variance matches summation formula", {
  test_cases <- list(
    c(J = 50, alpha = 0.5),
    c(J = 50, alpha = 1.0),
    c(J = 50, alpha = 2.0),
    c(J = 50, alpha = 5.0),
    c(J = 100, alpha = 2.0),
    c(J = 300, alpha = 2.0)
  )

  for (tc in test_cases) {
    J <- tc["J"]
    alpha <- tc["alpha"]
    expect_equal(
      var_K_given_alpha(J, alpha),
      var_K_given_alpha_sum(J, alpha),
      tolerance = 1e-10,
      info = sprintf("J=%d, alpha=%.1f", J, alpha)
    )
  }
})

# -----------------------------------------------------------------------------
# Test: Underdispersion Inequality
# -----------------------------------------------------------------------------

test_that("Underdispersion inequality holds: 0 < Var(K) < E[K]", {
  J <- 50
  for (alpha in c(0.1, 0.5, 1, 2, 5, 10, 20)) {
    mu <- mean_K_given_alpha(J, alpha)
    v <- var_K_given_alpha(J, alpha)
    expect_true(
      v > 0 && v < mu,
      info = sprintf("alpha=%.1f: mu=%.4f, v=%.4f", alpha, mu, v)
    )
  }
})

test_that("Dispersion index is always < 1", {
  J <- 50
  for (alpha in c(0.1, 0.5, 1, 2, 5, 10, 20)) {
    D <- dispersion_K_given_alpha(J, alpha)
    expect_true(D > 0 && D < 1, info = sprintf("alpha=%.1f: D=%.4f", alpha, D))
  }
})

# -----------------------------------------------------------------------------
# Test: Limiting Behavior
# -----------------------------------------------------------------------------

test_that("Limiting behavior as alpha -> 0+", {
  J <- 50

  # E[K] -> 1
  expect_equal(mean_K_given_alpha(J, 1e-10), 1.0, tolerance = 1e-5)

  # Var(K) -> 0
  expect_equal(var_K_given_alpha(J, 1e-10), 0.0, tolerance = 1e-5)
})

test_that("Limiting behavior as alpha -> infinity", {
  J <- 50

  # E[K] -> J
  expect_equal(mean_K_given_alpha(J, 1e6), J, tolerance = 0.01)

  # Var(K) -> 0
  expect_lt(var_K_given_alpha(J, 1e6), 0.01)
})

# -----------------------------------------------------------------------------
# Test: Derivative Verification
# -----------------------------------------------------------------------------

test_that("Derivative of mean matches finite difference", {
  test_cases <- list(
    c(J = 50, alpha = 0.5),
    c(J = 50, alpha = 2.0),
    c(J = 50, alpha = 10.0),
    c(J = 100, alpha = 2.0)
  )

  eps <- 1e-6

  for (tc in test_cases) {
    J <- tc["J"]
    alpha <- tc["alpha"]

    deriv_analytic <- dmean_dalpha(J, alpha)
    deriv_fd <- (mean_K_given_alpha(J, alpha + eps) -
                   mean_K_given_alpha(J, alpha - eps)) / (2 * eps)

    expect_equal(
      deriv_analytic, deriv_fd,
      tolerance = 1e-5,
      info = sprintf("J=%d, alpha=%.1f", J, alpha)
    )
  }
})

test_that("Derivative of variance matches finite difference", {
  J <- 50
  alpha <- 2.0
  eps <- 1e-6

  deriv_analytic <- dvar_dalpha(J, alpha)
  deriv_fd <- (var_K_given_alpha(J, alpha + eps) -
                 var_K_given_alpha(J, alpha - eps)) / (2 * eps)

  expect_equal(deriv_analytic, deriv_fd, tolerance = 1e-5)
})

test_that("Derivative is always positive (monotonicity)", {
  J <- 50
  alpha_seq <- seq(0.1, 10, by = 0.1)
  derivs <- dmean_dalpha(J, alpha_seq)

  expect_true(all(derivs > 0))
})

# -----------------------------------------------------------------------------
# Test: Vectorization
# -----------------------------------------------------------------------------

test_that("mean_K_given_alpha is properly vectorized", {
  J <- 50
  alpha_vec <- c(0.5, 1, 2, 5)

  means <- mean_K_given_alpha(J, alpha_vec)

  expect_length(means, 4)
  expect_equal(means[1], mean_K_given_alpha(J, 0.5))
  expect_equal(means[2], mean_K_given_alpha(J, 1))
  expect_equal(means[3], mean_K_given_alpha(J, 2))
  expect_equal(means[4], mean_K_given_alpha(J, 5))
})

test_that("var_K_given_alpha is properly vectorized", {
  J <- 50
  alpha_vec <- c(0.5, 1, 2, 5)

  vars <- var_K_given_alpha(J, alpha_vec)

  expect_length(vars, 4)
  expect_equal(vars[3], var_K_given_alpha(J, 2))
})

# -----------------------------------------------------------------------------
# Test: Monotonicity
# -----------------------------------------------------------------------------

test_that("E[K] is strictly increasing in alpha", {
  J <- 50
  alpha_seq <- seq(0.1, 10, by = 0.1)
  means <- mean_K_given_alpha(J, alpha_seq)

  expect_true(all(diff(means) > 0))
})

# -----------------------------------------------------------------------------
# Test: Golden Data Values
# -----------------------------------------------------------------------------

test_that("Golden test data matches", {
  # Values from golden_moments_cond.csv
  expect_equal(mean_K_given_alpha(50, 0.5), 2.93777, tolerance = 1e-4)
  expect_equal(var_K_given_alpha(50, 0.5), 1.70907, tolerance = 1e-4)

  expect_equal(mean_K_given_alpha(50, 1.0), 4.49921, tolerance = 1e-4)
  expect_equal(var_K_given_alpha(50, 1.0), 2.87407, tolerance = 1e-4)

  expect_equal(mean_K_given_alpha(50, 2.0), 7.03763, tolerance = 1e-4)
  expect_equal(var_K_given_alpha(50, 2.0), 4.53556, tolerance = 1e-4)

  expect_equal(mean_K_given_alpha(50, 5.0), 12.46049, tolerance = 1e-4)
  expect_equal(var_K_given_alpha(50, 5.0), 7.38611, tolerance = 1e-4)

  expect_equal(mean_K_given_alpha(100, 2.0), 8.39456, tolerance = 1e-4)
  expect_equal(var_K_given_alpha(100, 2.0), 5.85423, tolerance = 1e-4)

  expect_equal(mean_K_given_alpha(300, 2.0), 10.57197, tolerance = 1e-4)
  expect_equal(var_K_given_alpha(300, 2.0), 8.00550, tolerance = 1e-4)
})

# -----------------------------------------------------------------------------
# Test: Convenience Functions
# -----------------------------------------------------------------------------

test_that("moments_K_given_alpha returns correct structure", {
  m <- moments_K_given_alpha(50, 2.0)

  expect_named(m, c("mean", "var"))
  expect_equal(unname(m["mean"]), mean_K_given_alpha(50, 2.0))
  expect_equal(unname(m["var"]), var_K_given_alpha(50, 2.0))
})

test_that("moments_K_given_alpha handles vector alpha correctly", {
  # Vector alpha should return a matrix with 2 columns (mean, var)
  alpha_vec <- c(0.5, 1, 2, 5)
  m <- moments_K_given_alpha(50, alpha_vec)

  expect_true(is.matrix(m))
  expect_equal(nrow(m), length(alpha_vec))
  expect_equal(ncol(m), 2)
  expect_equal(colnames(m), c("mean", "var"))

  # Verify values match individual calls (use unname to ignore attribute differences)
  for (i in seq_along(alpha_vec)) {
    expect_equal(unname(m[i, "mean"]), unname(mean_K_given_alpha(50, alpha_vec[i])))
    expect_equal(unname(m[i, "var"]), unname(var_K_given_alpha(50, alpha_vec[i])))
  }
})

test_that("cv_K_given_alpha computes correctly", {
  J <- 50
  alpha <- 2.0

  cv <- cv_K_given_alpha(J, alpha)
  expected_cv <- sqrt(var_K_given_alpha(J, alpha)) / mean_K_given_alpha(J, alpha)

  expect_equal(cv, expected_cv)
})

test_that("summary_K_given_alpha returns complete summary", {
  s <- summary_K_given_alpha(50, 2.0)

  expect_true(is.list(s))
  expect_named(s, c("J", "alpha", "mean", "var", "sd", "cv", "dispersion", "dmean_dalpha"))
  expect_equal(s$J, 50)
  expect_equal(s$alpha, 2.0)
  expect_equal(s$sd, sqrt(s$var))
  expect_equal(s$cv, s$sd / s$mean)
  expect_equal(s$dispersion, s$var / s$mean)
})

# -----------------------------------------------------------------------------
# Test: Input Validation
# -----------------------------------------------------------------------------

test_that("mean_K_given_alpha validates inputs", {
  expect_error(mean_K_given_alpha(0, 2), "positive integer")
  expect_error(mean_K_given_alpha(-5, 2), "positive integer")
  expect_error(mean_K_given_alpha(50, 0), "positive")
  expect_error(mean_K_given_alpha(50, -1), "positive")
  expect_error(mean_K_given_alpha(1000, 2), "exceeds maximum")
})

test_that("var_K_given_alpha validates inputs", {
  expect_error(var_K_given_alpha(0, 2), "positive integer")
  expect_error(var_K_given_alpha(50, 0), "positive")
})

# -----------------------------------------------------------------------------
# Test: Edge Cases
# -----------------------------------------------------------------------------

test_that("Functions work at J boundary values", {
  # Minimum J
  expect_equal(mean_K_given_alpha(1, 2), 1)  # K can only be 1 when J=1
  expect_equal(var_K_given_alpha(1, 2), 0)   # No variance when K must be 1

  # Maximum J (uses max from constants)
  expect_silent(mean_K_given_alpha(500, 2))
})

test_that("Variance of K at J=1 is zero", {
  # When J=1, there's only one observation -> exactly 1 cluster
  for (alpha in c(0.5, 1, 2, 5)) {
    expect_equal(var_K_given_alpha(1, alpha), 0, tolerance = 1e-15)
  }
})

# -----------------------------------------------------------------------------
# Test: Cross-validation with RN-01 Verification Code
# -----------------------------------------------------------------------------

test_that("Verification functions return TRUE for valid inputs", {
  expect_true(validate_moments_conditional(50, 2.0, verbose = FALSE))
  expect_true(verify_underdispersion(50, verbose = FALSE))
  expect_true(verify_derivative(50, 2.0, verbose = FALSE))
})
