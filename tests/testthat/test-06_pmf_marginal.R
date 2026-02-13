# =============================================================================
# Test Suite: Module 06 - Marginal PMF of K_J
# =============================================================================
#
# Comprehensive tests for pmf_K_marginal() and related functions.
#
# Run with: devtools::test() or testthat::test_file("tests/testthat/test-06_pmf_marginal.R")
# =============================================================================

# -----------------------------------------------------------------------------
# Setup: Pre-compute Stirling numbers for all tests
# -----------------------------------------------------------------------------

logS_50 <- compute_log_stirling(50)
logS_100 <- compute_log_stirling(100)

# -----------------------------------------------------------------------------
# 1. Basic PMF Properties
# -----------------------------------------------------------------------------

test_that("Marginal PMF sums to 1", {
  pmf <- pmf_K_marginal(50, 1.5, 0.5, logS_50)
  expect_equal(sum(pmf), 1, tolerance = 1e-10)
})

test_that("P(K=0) = 0 for marginal PMF", {
  pmf <- pmf_K_marginal(50, 1.5, 0.5, logS_50)
  expect_equal(pmf[1], 0)
})

test_that("All probabilities are non-negative", {
  pmf <- pmf_K_marginal(50, 1.5, 0.5, logS_50)
  expect_true(all(pmf >= 0))
})

test_that("PMF has correct length", {
  pmf <- pmf_K_marginal(50, 1.5, 0.5, logS_50)
  expect_length(pmf, 51)  # k = 0, 1, ..., 50
})

test_that("Marginal PMF sums to 1 across parameter grid", {
  for (J in c(10, 50)) {
    for (a in c(1, 2)) {
      for (b in c(0.5, 1)) {
        pmf <- pmf_K_marginal(J, a, b, logS_100)
        expect_equal(sum(pmf), 1, tolerance = 1e-10,
                     info = sprintf("J=%d, a=%.2f, b=%.2f", J, a, b))
      }
    }
  }
})

# -----------------------------------------------------------------------------
# 2. Log-Space Implementation
# -----------------------------------------------------------------------------

test_that("log_pmf_K_marginal returns correct structure", {
  logp <- log_pmf_K_marginal(50, 1.5, 0.5, logS_50)
  expect_length(logp, 51)
  expect_equal(logp[1], -Inf)  # log(0) = -Inf
  expect_true(all(is.finite(logp[-1])))  # k=1..J should be finite
})

test_that("log_pmf and pmf are consistent", {
  logp <- log_pmf_K_marginal(50, 1.5, 0.5, logS_50)
  pmf <- pmf_K_marginal(50, 1.5, 0.5, logS_50)

  # Compare non-zero elements
  expect_equal(exp(logp[-1]), pmf[-1], tolerance = 1e-12)
})

# -----------------------------------------------------------------------------
# 3. Moments Consistency (PMF vs exact_K_moments)
# -----------------------------------------------------------------------------

test_that("Mean from PMF matches exact_K_moments", {
  pmf <- pmf_K_marginal(50, 1.5, 0.5, logS_50)
  mean_pmf <- sum((0:50) * pmf)
  exact <- exact_K_moments(50, 1.5, 0.5)
  expect_equal(mean_pmf, exact$mean, tolerance = 1e-6)
})

test_that("Variance from PMF matches exact_K_moments", {
  pmf <- pmf_K_marginal(50, 1.5, 0.5, logS_50)
  k <- 0:50
  mean_pmf <- sum(k * pmf)
  var_pmf <- sum(k^2 * pmf) - mean_pmf^2
  exact <- exact_K_moments(50, 1.5, 0.5)
  expect_equal(var_pmf, exact$var, tolerance = 1e-5)
})

test_that("Moments match across multiple parameter sets", {
  test_cases <- list(
    list(J = 50, a = 1.5, b = 0.5),
    list(J = 50, a = 2.0, b = 1.0),
    list(J = 30, a = 1.0, b = 0.5)
  )

  for (tc in test_cases) {
    pmf <- pmf_K_marginal(tc$J, tc$a, tc$b, logS_100)
    k <- 0:tc$J
    mean_pmf <- sum(k * pmf)
    exact <- exact_K_moments(tc$J, tc$a, tc$b)

    expect_equal(mean_pmf, exact$mean, tolerance = 1e-6,
                 info = sprintf("J=%d, a=%.2f, b=%.2f", tc$J, tc$a, tc$b))
  }
})

# -----------------------------------------------------------------------------
# 4. CDF Properties
# -----------------------------------------------------------------------------

test_that("CDF is non-decreasing", {
  cdf <- cdf_K_marginal(50, 1.5, 0.5, logS_50)
  expect_true(all(diff(cdf) >= -1e-15))
})

test_that("CDF starts at 0 and ends at 1", {
  cdf <- cdf_K_marginal(50, 1.5, 0.5, logS_50)
  expect_equal(cdf[1], 0)  # F(0) = P(K=0) = 0
  expect_equal(cdf[51], 1, tolerance = 1e-10)  # F(50) = 1
})

test_that("CDF has correct length", {
  cdf <- cdf_K_marginal(50, 1.5, 0.5, logS_50)
  expect_length(cdf, 51)
})

# -----------------------------------------------------------------------------
# 5. Quantile Function
# -----------------------------------------------------------------------------

test_that("Quantile function returns integers", {
  qs <- quantile_K_marginal(c(0.1, 0.5, 0.9), 50, 1.5, 0.5, logS_50)
  expect_type(qs, "integer")
})

test_that("Quantile function is vectorized", {
  probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  qs <- quantile_K_marginal(probs, 50, 1.5, 0.5, logS_50)
  expect_length(qs, 5)
})

test_that("Quantile-CDF consistency", {
  probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  qs <- quantile_K_marginal(probs, 50, 1.5, 0.5, logS_50)
  cdf <- cdf_K_marginal(50, 1.5, 0.5, logS_50)

  for (i in seq_along(probs)) {
    # CDF at quantile should be >= probability level
    expect_true(cdf[qs[i] + 1] >= probs[i],
                info = sprintf("p=%.2f, q=%d", probs[i], qs[i]))
  }
})

test_that("Quantile function handles edge cases", {
  q0 <- quantile_K_marginal(0, 50, 1.5, 0.5, logS_50)
  q1 <- quantile_K_marginal(1, 50, 1.5, 0.5, logS_50)

  expect_equal(q0, 0L)
  expect_equal(q1, 50L)
})

test_that("Quantiles are non-decreasing", {
  probs <- seq(0.1, 0.9, by = 0.1)
  qs <- quantile_K_marginal(probs, 50, 1.5, 0.5, logS_50)
  expect_true(all(diff(qs) >= 0))
})

# -----------------------------------------------------------------------------
# 6. Mode Function
# -----------------------------------------------------------------------------

test_that("Mode is an integer >= 1", {
  mode_val <- mode_K_marginal(50, 1.5, 0.5, logS_50)
  expect_type(mode_val, "integer")
  expect_true(mode_val >= 1)
})

test_that("Mode maximizes PMF", {
  pmf <- pmf_K_marginal(50, 1.5, 0.5, logS_50)
  mode_val <- mode_K_marginal(50, 1.5, 0.5, logS_50)

  # PMF at mode should be the maximum
  expect_equal(pmf[mode_val + 1], max(pmf))
})

# -----------------------------------------------------------------------------
# 7. Summary Function
# -----------------------------------------------------------------------------

test_that("Summary function returns correct structure", {
  summary_stats <- summary_K_marginal(50, 1.5, 0.5, logS_50)

  expect_true(is.list(summary_stats))
  expect_true("mean" %in% names(summary_stats))
  expect_true("var" %in% names(summary_stats))
  expect_true("sd" %in% names(summary_stats))
  expect_true("mode" %in% names(summary_stats))
  expect_true("median" %in% names(summary_stats))
  expect_true("quantiles" %in% names(summary_stats))
  expect_true("pmf" %in% names(summary_stats))
  expect_true("cdf" %in% names(summary_stats))
})

test_that("Summary statistics are consistent with exact moments", {
  summary_stats <- summary_K_marginal(50, 1.5, 0.5, logS_50)
  exact <- exact_K_moments(50, 1.5, 0.5)

  expect_equal(summary_stats$mean, exact$mean, tolerance = 1e-6)
  expect_equal(summary_stats$var, exact$var, tolerance = 1e-5)
})

test_that("Summary respects custom probs argument", {
  custom_probs <- c(0.025, 0.5, 0.975)
  summary_stats <- summary_K_marginal(50, 1.5, 0.5, logS_50, probs = custom_probs)

  expect_length(summary_stats$quantiles, 3)
  expect_equal(names(summary_stats$quantiles), c("q2", "q50", "q98"))
})

test_that("Summary mode and median are valid", {
  summary_stats <- summary_K_marginal(50, 1.5, 0.5, logS_50)

  expect_true(summary_stats$mode >= 1)
  expect_true(summary_stats$mode <= 50)
  expect_true(summary_stats$median >= 1)
  expect_true(summary_stats$median <= 50)
})

# -----------------------------------------------------------------------------
# 8. Golden Test Values (Python-verified)
# -----------------------------------------------------------------------------

test_that("Golden test: J=50, a=1.5, b=0.5", {
  summary_stats <- summary_K_marginal(50, 1.5, 0.5, logS_50)

  expect_equal(summary_stats$mean, 8.355487, tolerance = 1e-4)
  expect_equal(summary_stats$var, 22.768950, tolerance = 1e-3)
  expect_equal(summary_stats$mode, 6L)
  expect_equal(summary_stats$median, 8L)
})

test_that("Golden test: J=50, a=2.0, b=1.0", {
  summary_stats <- summary_K_marginal(50, 2.0, 1.0, logS_50)

  expect_equal(summary_stats$mean, 6.639693, tolerance = 1e-4)
  expect_equal(summary_stats$var, 12.954502, tolerance = 1e-3)
  expect_equal(summary_stats$mode, 5L)
  expect_equal(summary_stats$median, 6L)
})

test_that("Golden test: J=100, a=1.5, b=0.5", {
  summary_stats <- summary_K_marginal(100, 1.5, 0.5, logS_100)

  expect_equal(summary_stats$mean, 10.311916, tolerance = 1e-4)
  expect_equal(summary_stats$var, 39.260433, tolerance = 1e-3)
  expect_equal(summary_stats$mode, 7L)
  expect_equal(summary_stats$median, 9L)
})

test_that("Golden test: J=30, a=1.0, b=0.5", {
  summary_stats <- summary_K_marginal(30, 1.0, 0.5, logS_100)

  expect_equal(summary_stats$mean, 5.378524, tolerance = 1e-4)
  expect_equal(summary_stats$var, 12.352028, tolerance = 1e-3)
  expect_equal(summary_stats$mode, 2L)
  expect_equal(summary_stats$median, 5L)
})

test_that("Golden test: J=50, a=3.0, b=1.5", {
  summary_stats <- summary_K_marginal(50, 3.0, 1.5, logS_50)

  expect_equal(summary_stats$mean, 6.762616, tolerance = 1e-4)
  expect_equal(summary_stats$var, 10.442325, tolerance = 1e-3)
  expect_equal(summary_stats$mode, 6L)
  expect_equal(summary_stats$median, 6L)
})

# -----------------------------------------------------------------------------
# 9. Quadrature Convergence
# -----------------------------------------------------------------------------

test_that("PMF converges with increasing quadrature nodes", {
  M_vals <- c(40, 80, 120)
  pmfs <- lapply(M_vals, function(M) pmf_K_marginal(50, 1.5, 0.5, logS_50, M = M))

  # L1 difference should decrease
  diff1 <- sum(abs(pmfs[[2]] - pmfs[[1]]))
  diff2 <- sum(abs(pmfs[[3]] - pmfs[[2]]))

  expect_true(diff2 < diff1)
})

test_that("Mean stabilizes with increasing M", {
  M_vals <- c(40, 80, 120)
  means <- sapply(M_vals, function(M) {
    pmf <- pmf_K_marginal(50, 1.5, 0.5, logS_50, M = M)
    sum((0:50) * pmf)
  })

  # Changes should decrease
  change1 <- abs(means[2] - means[1])
  change2 <- abs(means[3] - means[2])

  expect_true(change2 < change1)
})

# -----------------------------------------------------------------------------
# 10. Input Validation
# -----------------------------------------------------------------------------

test_that("Invalid J is rejected", {
  expect_error(pmf_K_marginal(0, 1.5, 0.5, logS_50))
  expect_error(pmf_K_marginal(-1, 1.5, 0.5, logS_50))
})

test_that("Invalid a is rejected", {
  expect_error(pmf_K_marginal(50, 0, 0.5, logS_50))
  expect_error(pmf_K_marginal(50, -1, 0.5, logS_50))
})

test_that("Invalid b is rejected", {
  expect_error(pmf_K_marginal(50, 1.5, 0, logS_50))
  expect_error(pmf_K_marginal(50, 1.5, -1, logS_50))
})

test_that("Invalid M is rejected", {
  expect_error(pmf_K_marginal(50, 1.5, 0.5, logS_50, M = 0))
  expect_error(pmf_K_marginal(50, 1.5, 0.5, logS_50, M = -1))
  expect_error(pmf_K_marginal(50, 1.5, 0.5, logS_50, M = 1.5))
})

test_that("Invalid probability in quantile_K_marginal is rejected", {
  expect_error(quantile_K_marginal(-0.1, 50, 1.5, 0.5, logS_50))
  expect_error(quantile_K_marginal(1.1, 50, 1.5, 0.5, logS_50))
})

# -----------------------------------------------------------------------------
# 11. Variance Inflation Property
# -----------------------------------------------------------------------------

test_that("Marginal variance exceeds conditional variance at E[alpha]", {
  J <- 50
  a <- 1.5
  b <- 0.5

  # Marginal variance
  summary_marg <- summary_K_marginal(J, a, b, logS_50)
  var_marginal <- summary_marg$var

  # Conditional variance at E[alpha] = a/b
  E_alpha <- a / b
  cond_summary <- summary_K_given_alpha(J, E_alpha)
  var_conditional <- cond_summary$var


  # Marginal variance should be larger due to variance inflation
  # Var(K|a,b) = E[Var(K|alpha)] + Var[E(K|alpha)]
  expect_true(var_marginal > var_conditional)
})

# -----------------------------------------------------------------------------
# 12. Convenience Functions
# -----------------------------------------------------------------------------

test_that("mean_K_from_marginal_pmf works correctly", {
  mean1 <- mean_K_from_marginal_pmf(50, 1.5, 0.5, logS_50)
  exact <- exact_K_moments(50, 1.5, 0.5)
  expect_equal(mean1, exact$mean, tolerance = 1e-6)
})

test_that("var_K_from_marginal_pmf works correctly", {
  var1 <- var_K_from_marginal_pmf(50, 1.5, 0.5, logS_50)
  exact <- exact_K_moments(50, 1.5, 0.5)
  expect_equal(var1, exact$var, tolerance = 1e-5)
})

# -----------------------------------------------------------------------------
# 13. Verification Functions
# -----------------------------------------------------------------------------

test_that("verify_pmf_marginal_properties returns TRUE for valid input", {
  result <- verify_pmf_marginal_properties(50, 1.5, 0.5, logS_50, verbose = FALSE)
  expect_true(result)
})

test_that("verify_pmf_marginal_moments returns TRUE for valid input", {
  result <- verify_pmf_marginal_moments(50, 1.5, 0.5, logS_50, verbose = FALSE)
  expect_true(result)
})

test_that("verify_pmf_marginal_all returns TRUE", {
  result <- verify_pmf_marginal_all(verbose = FALSE)
  expect_true(result)
})
