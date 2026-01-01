# =============================================================================
# Unit Tests: Module 08 - First Stick-Breaking Weight (w₁) Distribution
# =============================================================================
#
# This file contains comprehensive tests for the w₁ distribution functions.
# Golden values are computed directly from the closed-form formulas and
# verified by the inverse identity CDF(Q(u)) = u.
#
# Author: JoonHo Lee
# Date: December 2025
# Part of: DPprior R Package
# =============================================================================

library(testthat)

# =============================================================================
# Test 1: CDF and Quantile are Inverses
# =============================================================================

test_that("CDF and quantile are inverses", {
  for (a in c(0.5, 1, 2, 5)) {
    for (b in c(0.5, 1, 2)) {
      for (u in c(0.1, 0.5, 0.9)) {
        q <- quantile_w1(u, a, b)

        # For some (a, b, u) combinations, the true quantile is so close to 1
        # that it rounds to 1.0 in double precision. Skip strict inversion
        # checks in that numerically-degenerate regime.
        if (q >= 1 - 1e-15) next

        p <- cdf_w1(q, a, b)
        expect_equal(p, u, tolerance = 1e-10,
                     info = sprintf("a=%.1f, b=%.1f, u=%.1f", a, b, u))
      }
    }
  }
})


# =============================================================================
# Test 2: CDF is Monotonically Increasing
# =============================================================================

test_that("CDF is monotonically increasing", {
  a <- 2; b <- 1
  x_seq <- seq(0.01, 0.99, by = 0.01)
  cdf_vals <- cdf_w1(x_seq, a, b)
  expect_true(all(diff(cdf_vals) > 0))
})


# =============================================================================
# Test 3: Density Integrates to 1
# =============================================================================

test_that("Density integrates to 1", {
  a <- 2; b <- 1
  # Use upper limit very close to 1 for better coverage of the tail
  integral <- integrate(function(x) density_w1(x, a, b), 1e-6, 1 - 1e-16)$value
  expect_equal(integral, 1, tolerance = 1e-3)
})


# =============================================================================
# Test 4: Survival Function Equals 1 - CDF
# =============================================================================

test_that("prob_w1_exceeds matches 1 - CDF", {
  a <- 1.6; b <- 1.22
  t <- 0.5
  expect_equal(prob_w1_exceeds(t, a, b), 1 - cdf_w1(t, a, b), tolerance = 1e-12)
})


# =============================================================================
# Test 5: Conditional Expectation E[w₁|α] = 1/(1+α)
# =============================================================================

test_that("Conditional expectation E[w1|alpha] = 1/(1+alpha)", {
  alpha <- 2
  cond_mean <- integrate(function(x) x * alpha * (1 - x)^(alpha - 1), 0, 1)$value
  expect_equal(cond_mean, 1 / (1 + alpha), tolerance = 1e-12)
})


# =============================================================================
# Test 6: Golden Values (High-Precision, Verified by Inverse Identity)
# =============================================================================

test_that("w1 quantiles match golden values (closed form)", {

  # These values are computed from the closed-form quantile function and
  # verified by the inverse identity F(Q(u)) = u.
  # Note: Using tolerance = 1e-9 to account for floating point variations
  # across different platforms and R versions.
  golden <- data.frame(
    a = c(1.6, 2.0, 0.5),
    b = c(1.22, 1.0, 1.0),
    q50 = c(0.4839219199376489, 0.3391401985931721, 0.9502129316321360),
    prob_gt_50 = c(0.4868311039310716, 0.3488273884247458, 0.7685155230091048)
  )

  for (i in seq_len(nrow(golden))) {
    actual_q50 <- quantile_w1(0.5, golden$a[i], golden$b[i])
    actual_prob <- prob_w1_exceeds(0.5, golden$a[i], golden$b[i])

    expect_equal(actual_q50, golden$q50[i], tolerance = 1e-9,
                 info = sprintf("q50 for a=%.1f, b=%.2f", golden$a[i], golden$b[i]))
    expect_equal(actual_prob, golden$prob_gt_50[i], tolerance = 1e-9,
                 info = sprintf("P(w1>0.5) for a=%.1f, b=%.2f",
                                golden$a[i], golden$b[i]))
  }
})


# =============================================================================
# Test 7: Boundary Handling
# =============================================================================

test_that("CDF handles boundaries correctly", {
  a <- 2; b <- 1

  # At or below 0
  expect_equal(cdf_w1(0, a, b), 0)
  expect_equal(cdf_w1(-1, a, b), 0)
  expect_equal(cdf_w1(-Inf, a, b), 0)

  # At or above 1
  expect_equal(cdf_w1(1, a, b), 1)

  expect_equal(cdf_w1(2, a, b), 1)
  expect_equal(cdf_w1(Inf, a, b), 1)

  # NA handling
  expect_true(is.na(cdf_w1(NA, a, b)))
  expect_true(is.na(cdf_w1(NaN, a, b)))
})


test_that("prob_w1_exceeds handles boundaries correctly", {
  a <- 2; b <- 1

  # At or below 0
  expect_equal(prob_w1_exceeds(0, a, b), 1)
  expect_equal(prob_w1_exceeds(-1, a, b), 1)

  # At or above 1
  expect_equal(prob_w1_exceeds(1, a, b), 0)
  expect_equal(prob_w1_exceeds(2, a, b), 0)
})


test_that("quantile_w1 handles boundaries correctly", {
  a <- 2; b <- 1

  expect_equal(quantile_w1(0, a, b), 0)
  expect_equal(quantile_w1(1, a, b), 1)
})


test_that("density_w1 handles boundaries correctly", {
  a <- 2; b <- 1

  # Outside (0, 1) -> 0 (or -Inf on log scale)
  expect_equal(density_w1(0, a, b), 0)
  expect_equal(density_w1(1, a, b), 0)
  expect_equal(density_w1(-1, a, b), 0)
  expect_equal(density_w1(2, a, b), 0)

  expect_equal(density_w1(0, a, b, log = TRUE), -Inf)
  expect_equal(density_w1(1, a, b, log = TRUE), -Inf)
})


# =============================================================================
# Test 8: Input Validation
# =============================================================================

test_that("input validation works correctly", {
  # Negative parameters
  expect_error(cdf_w1(0.5, a = -1, b = 1))
  expect_error(cdf_w1(0.5, a = 1, b = -1))
  expect_error(cdf_w1(0.5, a = 0, b = 1))
  expect_error(cdf_w1(0.5, a = 1, b = 0))

  # Invalid probability
  expect_error(quantile_w1(-0.1, a = 1, b = 1))
  expect_error(quantile_w1(1.1, a = 1, b = 1))
})


# =============================================================================
# Test 9: Vectorization
# =============================================================================

test_that("functions are properly vectorized", {
  a <- 2; b <- 1

  # CDF
  x_vec <- c(0.1, 0.3, 0.5, 0.7)
  cdf_vals <- cdf_w1(x_vec, a, b)
  expect_length(cdf_vals, 4)
  expect_true(all(cdf_vals >= 0 & cdf_vals <= 1))

  # Quantile
  u_vec <- c(0.1, 0.5, 0.9)
  q_vals <- quantile_w1(u_vec, a, b)
  expect_length(q_vals, 3)
  expect_true(all(q_vals >= 0 & q_vals <= 1))

  # Density
  dens_vals <- density_w1(x_vec, a, b)
  expect_length(dens_vals, 4)
  expect_true(all(dens_vals > 0))
})


# =============================================================================
# Test 10: Mean and Variance via Monte Carlo
# =============================================================================

test_that("mean_w1 matches Monte Carlo estimate", {
  skip_on_cran()  # Skip on CRAN due to Monte Carlo variance

  set.seed(42)
  a <- 2; b <- 1
  n_samples <- 50000

  samples <- rw1(n_samples, a, b)
  mc_mean <- mean(samples)
  quad_mean <- mean_w1(a, b)

  # Allow 5% relative error for MC
  expect_equal(quad_mean, mc_mean, tolerance = 0.05)
})


test_that("var_w1 matches Monte Carlo estimate", {
  skip_on_cran()

  set.seed(42)
  a <- 2; b <- 1
  n_samples <- 50000

  samples <- rw1(n_samples, a, b)
  mc_var <- var(samples)
  quad_var <- var_w1(a, b)

  # Allow 10% relative error for MC variance estimate
  expect_equal(quad_var, mc_var, tolerance = 0.1)
})


# =============================================================================
# Test 11: summary_w1 Output Structure
# =============================================================================

test_that("summary_w1 returns correct structure", {
  s <- summary_w1(a = 2, b = 1)

  expect_s3_class(s, "w1_summary")
  expect_named(s, c("mean", "var", "sd", "median", "quantiles",
                    "prob_gt_50", "prob_gt_90", "params"))

  expect_true(s$mean > 0 && s$mean < 1)
  expect_true(s$var > 0)
  expect_equal(s$sd, sqrt(s$var))
  expect_true(s$median > 0 && s$median < 1)
  expect_length(s$quantiles, 5)  # Default 5 quantiles
  expect_true(s$prob_gt_50 >= 0 && s$prob_gt_50 <= 1)
  expect_true(s$prob_gt_90 >= 0 && s$prob_gt_90 <= 1)
  expect_equal(s$params, list(a = 2, b = 1))
})


# =============================================================================
# Test 12: rw1 Random Generation
# =============================================================================

test_that("rw1 generates valid samples", {
  set.seed(42)
  samples <- rw1(1000, a = 2, b = 1)

  expect_length(samples, 1000)
  expect_true(all(samples > 0))
  expect_true(all(samples < 1))
})


test_that("rw1 validates input", {
  expect_error(rw1(-1, a = 2, b = 1))
  expect_error(rw1(0, a = 2, b = 1))
  expect_error(rw1(1.5, a = 2, b = 1))
  expect_error(rw1(100, a = -1, b = 1))
})


# =============================================================================
# Test 13: Log-Density Consistency
# =============================================================================

test_that("log-density is consistent with density", {
  a <- 2; b <- 1
  x <- c(0.1, 0.3, 0.5, 0.7)

  log_dens <- density_w1(x, a, b, log = TRUE)
  dens <- density_w1(x, a, b, log = FALSE)

  expect_equal(log_dens, log(dens), tolerance = 1e-10)
})


# =============================================================================
# Test 14: Numerical Stability - expm1 Test
# =============================================================================

test_that("CDF is stable for very small x (expm1 test)", {
  a <- 2; b <- 1

  # For very small x, CDF should be close to 0
  # Using expm1 is critical here
  x_small <- c(1e-10, 1e-8, 1e-6, 1e-4)
  cdf_vals <- cdf_w1(x_small, a, b)

  # All should be positive and increasing
  expect_true(all(cdf_vals > 0))
  expect_true(all(diff(cdf_vals) > 0))

  # Check against naive calculation for moderate x
  # where both methods should agree
  x_moderate <- 0.1
  cdf_stable <- cdf_w1(x_moderate, a, b)

  # Naive calculation (less stable)
  denom <- b - log1p(-x_moderate)
  cdf_naive <- 1 - (b / denom)^a

  expect_equal(cdf_stable, cdf_naive, tolerance = 1e-10)
})


# =============================================================================
# Test 15: Consistency Check - CDF Derivation
# =============================================================================

test_that("CDF matches numerical integration of conditional", {
  # Verify: F(x) = ∫ [1 - (1-x)^α] × p(α) dα
  test_cases <- list(
    c(a = 2.0, b = 1.0, x = 0.3),
    c(a = 1.6, b = 1.22, x = 0.5),
    c(a = 0.5, b = 1.0, x = 0.3)
  )

  for (case in test_cases) {
    a <- case["a"]
    b <- case["b"]
    x <- case["x"]

    # Numerical integration
    cdf_numerical <- integrate(
      function(alpha) (1 - (1 - x)^alpha) * dgamma(alpha, shape = a, rate = b),
      lower = 0, upper = 100
    )$value

    # Closed form
    cdf_closed <- cdf_w1(x, a, b)

    expect_equal(cdf_numerical, cdf_closed, tolerance = 1e-6,
                 info = sprintf("Derivation check: a=%.1f, b=%.2f, x=%.1f",
                                a, b, x))
  }
})
