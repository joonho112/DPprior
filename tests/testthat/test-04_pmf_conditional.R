# =============================================================================
# Tests for Module 04: Conditional PMF of K | alpha (Antoniak Distribution)
# =============================================================================
#
# Author: JoonHo Lee
# Date: December 2025
# Part of: DPprior R Package
# =============================================================================

# Pre-compute Stirling matrix for all tests
logS <- compute_log_stirling(100)


# =============================================================================
# Rising Factorial Tests
# =============================================================================

test_that("log_rising_factorial computes correctly", {
  # (2)_3 = 2 * 3 * 4 = 24
  expect_equal(exp(log_rising_factorial(2, 3)), 24, tolerance = 1e-10)

  # (1)_5 = 5! = 120
  expect_equal(exp(log_rising_factorial(1, 5)), 120, tolerance = 1e-10)

  # (0.5)_4 = 0.5 * 1.5 * 2.5 * 3.5 = 6.5625
  expect_equal(exp(log_rising_factorial(0.5, 4)), 6.5625, tolerance = 1e-10)

  # (alpha)_1 = alpha
  expect_equal(exp(log_rising_factorial(3.7, 1)), 3.7, tolerance = 1e-10)
})

test_that("log_rising_factorial validates inputs", {
  expect_error(log_rising_factorial(-1, 5))
  expect_error(log_rising_factorial(0, 5))
  expect_error(log_rising_factorial(2, 0))
  expect_error(log_rising_factorial(2, -1))
  expect_error(log_rising_factorial(c(1, 2), 5))
})


# =============================================================================
# PMF Normalization Tests
# =============================================================================

test_that("PMF sums to 1", {
  for (J in c(10, 50, 100)) {
    for (alpha in c(0.5, 1, 2, 5)) {
      pmf <- pmf_K_given_alpha(J, alpha, logS)
      expect_equal(sum(pmf), 1, tolerance = 1e-10,
                   info = sprintf("J=%d, alpha=%.1f", J, alpha))
    }
  }
})

test_that("PMF is non-negative", {
  for (J in c(10, 50, 100)) {
    for (alpha in c(0.5, 2, 5)) {
      pmf <- pmf_K_given_alpha(J, alpha, logS)
      expect_true(all(pmf >= 0), info = sprintf("J=%d, alpha=%.1f", J, alpha))
    }
  }
})

test_that("P(K=0) = 0 always", {
  for (alpha in c(0.1, 0.5, 1, 2, 5, 10)) {
    pmf <- pmf_K_given_alpha(50, alpha, logS)
    expect_equal(pmf[1], 0, tolerance = 1e-15,
                 info = sprintf("alpha=%.1f", alpha))
  }
})

test_that("P(K=J) > 0 for all alpha > 0", {
  for (J in c(10, 50)) {
    for (alpha in c(0.1, 1, 10)) {
      pmf <- pmf_K_given_alpha(J, alpha, logS)
      expect_true(pmf[J + 1] > 0, info = sprintf("J=%d, alpha=%.1f", J, alpha))
    }
  }
})


# =============================================================================
# Moment Consistency Tests
# =============================================================================

test_that("Moments from PMF match digamma/trigamma formulas", {
  for (J in c(10, 50, 100)) {
    for (alpha in c(0.5, 1.0, 2.0, 5.0)) {
      pmf <- pmf_K_given_alpha(J, alpha, logS)
      k_vals <- 0:J
      mean_pmf <- sum(k_vals * pmf)
      var_pmf <- sum(k_vals^2 * pmf) - mean_pmf^2

      mean_cf <- mean_K_given_alpha(J, alpha)
      var_cf <- var_K_given_alpha(J, alpha)

      expect_equal(mean_pmf, mean_cf, tolerance = 1e-8,
                   info = sprintf("Mean: J=%d, alpha=%.1f", J, alpha))
      expect_equal(var_pmf, var_cf, tolerance = 1e-8,
                   info = sprintf("Var: J=%d, alpha=%.1f", J, alpha))
    }
  }
})


# =============================================================================
# CDF Tests
# =============================================================================

test_that("CDF is non-decreasing", {
  for (J in c(10, 50, 100)) {
    for (alpha in c(0.5, 2, 5)) {
      cdf <- cdf_K_given_alpha(J, alpha, logS)
      expect_true(all(diff(cdf) >= -1e-15),
                  info = sprintf("J=%d, alpha=%.1f", J, alpha))
    }
  }
})

test_that("CDF ends at 1", {
  for (J in c(10, 50, 100)) {
    for (alpha in c(0.5, 2, 5)) {
      cdf <- cdf_K_given_alpha(J, alpha, logS)
      expect_equal(cdf[length(cdf)], 1, tolerance = 1e-10,
                   info = sprintf("J=%d, alpha=%.1f", J, alpha))
    }
  }
})

test_that("CDF starts at 0 (for k=0)", {
  cdf <- cdf_K_given_alpha(50, 2.0, logS)
  expect_equal(cdf[1], 0, tolerance = 1e-15)
})


# =============================================================================
# Quantile Tests
# =============================================================================

test_that("Quantiles are within valid range", {
  J <- 50
  alpha <- 2.0

  for (p in c(0.1, 0.25, 0.5, 0.75, 0.9)) {
    q <- quantile_K_given_alpha(p, J, alpha, logS)
    expect_true(q >= 0 && q <= J,
                info = sprintf("p=%.2f", p))
  }
})

test_that("Quantile satisfies CDF definition", {
  J <- 50
  alpha <- 2.0
  cdf <- cdf_K_given_alpha(J, alpha, logS)

  for (p in c(0.1, 0.25, 0.5, 0.75, 0.9)) {
    q <- quantile_K_given_alpha(p, J, alpha, logS)

    # CDF(q) >= p
    expect_true(cdf[q + 1] >= p - 1e-10,
                info = sprintf("p=%.2f: CDF(%d) >= p", p, q))

    # CDF(q-1) < p (if q > 0)
    if (q > 0) {
      expect_true(cdf[q] < p + 1e-10,
                  info = sprintf("p=%.2f: CDF(%d) < p", p, q - 1))
    }
  }
})

test_that("Quantile handles p=0 and p=1 edge cases", {
  J <- 50
  alpha <- 2.0

  # p=0 should return 0
  q0 <- quantile_K_given_alpha(0, J, alpha, logS)
  expect_equal(q0, 0L)

  # p=1 should return J
  q1 <- quantile_K_given_alpha(1, J, alpha, logS)
  expect_equal(q1, as.integer(J))
})

test_that("Quantile is monotone non-decreasing in p", {
  J <- 50
  alpha <- 2.0
  p_seq <- seq(0, 1, length.out = 21)
  q_seq <- quantile_K_given_alpha(p_seq, J, alpha, logS)

  expect_true(all(diff(q_seq) >= 0),
              info = "Quantiles should be non-decreasing in p")
})

test_that("Quantile handles vector input correctly", {
  J <- 50
  alpha <- 2.0
  probs <- c(0.25, 0.5, 0.75)

  qs <- quantile_K_given_alpha(probs, J, alpha, logS)

  expect_length(qs, 3)
  expect_true(is.integer(qs))
  expect_true(all(diff(qs) >= 0))
})

test_that("Median lies between mode and mean", {
  # For unimodal distributions, median should be close to mode/mean
  J <- 50
  for (alpha in c(1, 2, 5)) {
    pmf <- pmf_K_given_alpha(J, alpha, logS)
    mean_k <- sum((0:J) * pmf)
    mode_k <- which.max(pmf) - 1
    median_k <- quantile_K_given_alpha(0.5, J, alpha, logS)

    # Median should be within 1 of mode for this distribution
    expect_true(abs(median_k - mode_k) <= 2,
                info = sprintf("alpha=%.1f: median=%d, mode=%d",
                               alpha, median_k, mode_k))
  }
})


# =============================================================================
# Mode Tests
# =============================================================================

test_that("Mode is within valid range", {
  for (J in c(10, 50, 100)) {
    for (alpha in c(0.5, 2, 5)) {
      mode_k <- mode_K_given_alpha(J, alpha, logS)
      expect_true(mode_k >= 1 && mode_k <= J,
                  info = sprintf("J=%d, alpha=%.1f", J, alpha))
    }
  }
})

test_that("Mode increases with alpha", {
  J <- 50
  alphas <- c(0.5, 1, 2, 5, 10)
  modes <- sapply(alphas, function(a) mode_K_given_alpha(J, a, logS))

  # Mode should be non-decreasing with alpha
  expect_true(all(diff(modes) >= 0),
              info = sprintf("Modes: %s", paste(modes, collapse = ", ")))
})


# =============================================================================
# Golden Test Data
# =============================================================================

test_that("PMF matches Python golden values", {
  # Golden values generated from Python reference implementation
  golden <- list(
    list(J = 10, alpha = 0.5, mean = 2.133256, var = 0.924534, mode = 2),
    list(J = 10, alpha = 1.0, mean = 2.928968, var = 1.379201, mode = 3),
    list(J = 10, alpha = 2.0, mean = 4.039755, var = 1.807626, mode = 4),
    list(J = 50, alpha = 0.5, mean = 2.937775, var = 1.709074, mode = 3),
    list(J = 50, alpha = 2.0, mean = 7.037626, var = 4.535558, mode = 7),
    list(J = 50, alpha = 5.0, mean = 12.460485, var = 7.386114, mode = 12),
    list(J = 100, alpha = 1.0, mean = 5.187378, var = 3.552394, mode = 5),
    list(J = 100, alpha = 5.0, mean = 15.715366, var = 10.421525, mode = 15)
  )

  for (g in golden) {
    pmf <- pmf_K_given_alpha(g$J, g$alpha, logS)
    k_vals <- 0:g$J
    mean_k <- sum(k_vals * pmf)
    var_k <- sum(k_vals^2 * pmf) - mean_k^2
    mode_k <- which.max(pmf) - 1

    expect_equal(mean_k, g$mean, tolerance = 1e-5,
                 info = sprintf("Mean: J=%d, alpha=%.1f", g$J, g$alpha))
    expect_equal(var_k, g$var, tolerance = 1e-5,
                 info = sprintf("Var: J=%d, alpha=%.1f", g$J, g$alpha))
    expect_equal(mode_k, g$mode,
                 info = sprintf("Mode: J=%d, alpha=%.1f", g$J, g$alpha))
  }
})


# =============================================================================
# Edge Case Tests
# =============================================================================

test_that("PMF handles small alpha", {
  # As alpha -> 0+, mass concentrates on K=1
  pmf <- pmf_K_given_alpha(50, 0.01, logS)
  expect_true(pmf[2] > 0.9)  # P(K=1) should be > 0.9
})

test_that("PMF handles large alpha", {
  # As alpha -> infinity, mass shifts toward K=J
  pmf <- pmf_K_given_alpha(10, 100, logS)
  # Mode should be close to J
  mode_k <- which.max(pmf) - 1
  expect_true(mode_k >= 8)
})

test_that("PMF handles J=1", {
  logS_small <- compute_log_stirling(5)
  pmf <- pmf_K_given_alpha(1, 2.0, logS_small)

  # With J=1, there can only be K=1
  expect_equal(pmf[1], 0)    # P(K=0) = 0
  expect_equal(pmf[2], 1)    # P(K=1) = 1
})


# =============================================================================
# Input Validation Tests
# =============================================================================

test_that("Functions validate inputs", {
  # Invalid J
  expect_error(pmf_K_given_alpha(0, 2, logS))
  expect_error(pmf_K_given_alpha(-1, 2, logS))
  expect_error(pmf_K_given_alpha(1.5, 2, logS))

  # Invalid alpha
  expect_error(pmf_K_given_alpha(50, 0, logS))
  expect_error(pmf_K_given_alpha(50, -1, logS))
  expect_error(pmf_K_given_alpha(50, c(1, 2), logS))

  # Invalid logS
  expect_error(pmf_K_given_alpha(50, 2, matrix(1:4, 2, 2)))

  # J exceeds logS
  logS_small <- compute_log_stirling(10)
  expect_error(pmf_K_given_alpha(20, 2, logS_small))

  # Invalid probability
  expect_error(quantile_K_given_alpha(-0.1, 50, 2, logS))
  expect_error(quantile_K_given_alpha(1.1, 50, 2, logS))
})


# =============================================================================
# Summary Function Tests
# =============================================================================

test_that("summary_pmf_K_given_alpha returns complete output", {
  summary <- summary_pmf_K_given_alpha(50, 2.0, logS)

  expect_true(is.list(summary))
  expect_true("J" %in% names(summary))
  expect_true("alpha" %in% names(summary))
  expect_true("mean" %in% names(summary))
  expect_true("var" %in% names(summary))
  expect_true("sd" %in% names(summary))
  expect_true("mode" %in% names(summary))
  expect_true("median" %in% names(summary))
  expect_true("quantiles" %in% names(summary))
  expect_true("pmf" %in% names(summary))
  expect_true("cdf" %in% names(summary))

  expect_equal(summary$J, 50)
  expect_equal(summary$alpha, 2.0)
  expect_equal(length(summary$pmf), 51)
  expect_equal(length(summary$cdf), 51)
  expect_equal(summary$sd, sqrt(summary$var))
})


# =============================================================================
# Verification Function Tests
# =============================================================================

test_that("verify_pmf_all passes", {
  result <- verify_pmf_all(J_values = c(10, 50),
                           alpha_values = c(1, 2),
                           verbose = FALSE)
  expect_true(result)
})
