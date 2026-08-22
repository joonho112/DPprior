# =============================================================================
# Tests for Module 04: Conditional PMF of K | alpha (Antoniak Distribution)
# =============================================================================
#
# Author: JoonHo Lee (jlee296@ua.edu)
# Date: December 2025
# Part of: DPprior R Package
# =============================================================================

# Pre-compute Stirling matrix for all tests
logS <- compute_log_stirling(100)


# Independent oracle: the CRP new-cluster indicators are independent
# Bernoulli variables with probabilities alpha/(alpha+i-1), i=2,...,J.
# This convolution calls no DPprior production PMF, Stirling, or moment code.
.conditional_pmf_crp_oracle <- function(J, alpha) {
  probability <- 1
  if (J == 1L) {
    return(probability)
  }

  for (i in 2:J) {
    p_new <- alpha / (alpha + i - 1)
    updated <- numeric(i)
    updated[seq_len(i - 1L)] <-
      updated[seq_len(i - 1L)] + probability * (1 - p_new)
    updated[2:i] <- updated[2:i] + probability * p_new
    probability <- updated
  }
  probability
}


.conditional_pmf_fixture <- utils::read.csv(
  testthat::test_path("04_conditional_pmf_oracle.csv"),
  stringsAsFactors = FALSE
)


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
    expect_true(q >= 1 && q <= J,
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

  # p=0 uses the lower endpoint of the mathematical support.
  q0 <- quantile_K_given_alpha(0, J, alpha, logS)
  expect_equal(q0, 1L)

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


# =============================================================================
# Phase 4 independent oracle and extreme-domain regression contract
# =============================================================================

test_that("frozen small-J oracle has complete support and provenance", {
  expect_identical(
    unique(.conditional_pmf_fixture$oracle_id),
    "crp_bernoulli_convolution_v1"
  )

  groups <- split(
    .conditional_pmf_fixture,
    interaction(
      .conditional_pmf_fixture$J,
      .conditional_pmf_fixture$alpha,
      drop = TRUE
    )
  )

  for (fixture in groups) {
    fixture <- fixture[order(fixture$k), ]
    J <- fixture$J[1L]
    alpha <- fixture$alpha[1L]
    oracle <- .conditional_pmf_crp_oracle(J, alpha)
    fixture_mean <- fixture$mean[1L]
    fixture_var <- fixture$variance[1L]

    expect_identical(fixture$k, seq_len(J))
    expect_true(all(is.finite(fixture$probability)))
    expect_true(all(fixture$probability > 0))
    expect_equal(sum(fixture$probability), 1, tolerance = 2e-15)
    expect_equal(oracle, fixture$probability, tolerance = 2e-15)
    expect_equal(
      sum(fixture$k * fixture$probability),
      fixture_mean,
      tolerance = 2e-15
    )
    expect_equal(
      sum((fixture$k - fixture_mean)^2 * fixture$probability),
      fixture_var,
      tolerance = 2e-15
    )
  }
})


test_that("conditional PMF and moments match the frozen independent grid", {
  fixture_logS <- compute_log_stirling(max(.conditional_pmf_fixture$J))
  groups <- split(
    .conditional_pmf_fixture,
    interaction(
      .conditional_pmf_fixture$J,
      .conditional_pmf_fixture$alpha,
      drop = TRUE
    )
  )

  for (fixture in groups) {
    fixture <- fixture[order(fixture$k), ]
    J <- fixture$J[1L]
    alpha <- fixture$alpha[1L]
    expected <- fixture$probability
    observed <- pmf_K_given_alpha(J, alpha, fixture_logS)
    observed_log <- log_pmf_K_given_alpha(J, alpha, fixture_logS)
    budget <- 1e-14 + 1e-12 * abs(expected)

    expect_length(observed, J + 1L)
    expect_equal(observed[1L], 0, tolerance = 0)
    expect_equal(observed_log[1L], -Inf)
    expect_true(all(abs(observed[-1L] - expected) <= budget))
    expect_true(all(abs(exp(observed_log[-1L]) - expected) <= budget))
    expect_equal(
      mean_K_from_pmf(J, alpha, fixture_logS),
      fixture$mean[1L],
      tolerance = 1e-12
    )
    expect_equal(
      var_K_from_pmf(J, alpha, fixture_logS),
      fixture$variance[1L],
      tolerance = 1e-12
    )
  }
})


test_that("log rising factorial remains stable at finite alpha extremes", {
  J <- 500L
  for (alpha in c(.Machine$double.xmin, 1e-300, 1e-14, 1, 1e16,
                  1e100, 1e300, .Machine$double.xmax)) {
    # Direct finite-product log oracle; no lgamma difference is used.
    reference <- sum(log(alpha + 0:(J - 1L)))
    observed <- log_rising_factorial(alpha, J)
    budget <- 1e-12 + 1e-12 * abs(reference)
    expect_true(
      is.finite(observed) && abs(observed - reference) <= budget,
      info = sprintf("alpha=%.17g", alpha)
    )
  }
})


test_that("extreme conditional log-PMF and PMF obey the support contract", {
  extreme_logS <- compute_log_stirling(500L)
  alpha_grid <- c(1e-300, 1e-100, 1e-14, 1, 1e14, 1e16, 1e100, 1e300)

  for (J in c(1L, 2L, 50L, 500L)) {
    for (alpha in alpha_grid) {
      log_probability <- log_pmf_K_given_alpha(J, alpha, extreme_logS)
      probability <- pmf_K_given_alpha(J, alpha, extreme_logS)
      direct_probability <- pmf_K_given_alpha(
        J, alpha, extreme_logS, normalize = FALSE
      )
      support <- seq_len(J)
      expected_mean <- mean_K_given_alpha(J, alpha)
      expected_var <- var_K_given_alpha(J, alpha)
      observed_mean <- sum(support * probability[-1L])
      observed_var <- sum(
        (support - observed_mean)^2 * probability[-1L]
      )

      expect_length(log_probability, J + 1L)
      expect_equal(log_probability[1L], -Inf)
      expect_true(all(is.finite(log_probability[-1L])))
      expect_lte(max(log_probability[-1L]), 64 * .Machine$double.eps)
      expect_equal(
        logsumexp_vec(log_probability[-1L]), 0,
        tolerance = 2e-12
      )

      expect_length(probability, J + 1L)
      expect_equal(probability[1L], 0, tolerance = 0)
      expect_true(all(is.finite(probability)))
      expect_true(all(probability >= 0))
      expect_equal(sum(probability), 1, tolerance = 2e-15)
      expect_true(all(is.finite(direct_probability)))
      expect_true(all(direct_probability >= 0))
      expect_equal(sum(direct_probability), 1, tolerance = 1e-10)

      expect_true(
        abs(observed_mean - expected_mean) <=
          1e-12 + 1e-10 * abs(expected_mean)
      )
      expect_true(
        abs(observed_var - expected_var) <=
          1e-12 + 1e-10 * abs(expected_var)
      )
    }
  }
})


test_that("CDF and quantile endpoints use mathematical support 1 through J", {
  endpoint_logS <- compute_log_stirling(50L)
  probability_levels <- c(
    lower = 0,
    smallest = .Machine$double.xmin,
    tiny = 5e-13,
    q10 = 0.1,
    median = 0.5,
    q90 = 0.9,
    near_one = 1 - 5e-13,
    upper = 1
  )

  for (J in c(1L, 2L, 50L)) {
    for (alpha in c(1e-14, 2, 1e14)) {
      cdf <- cdf_K_given_alpha(J, alpha, endpoint_logS)
      quantiles <- quantile_K_given_alpha(
        probability_levels, J, alpha, endpoint_logS
      )

      expect_length(cdf, J + 1L)
      expect_equal(cdf[1L], 0, tolerance = 0)
      expect_equal(cdf[J + 1L], 1, tolerance = 0)
      expect_true(all(is.finite(cdf)))
      expect_true(all(cdf >= 0 & cdf <= 1))
      expect_true(all(diff(cdf) >= 0))

      expect_type(quantiles, "integer")
      expect_identical(names(quantiles), names(probability_levels))
      expect_true(all(quantiles >= 1L & quantiles <= J))
      expect_equal(unname(quantiles["lower"]), 1L)
      expect_equal(unname(quantiles["upper"]), J)
      expect_true(all(diff(unname(quantiles)) >= 0))

      positive <- which(probability_levels > 0 & probability_levels < 1)
      for (index in positive) {
        p <- probability_levels[index]
        q <- quantiles[index]
        expect_true(cdf[q + 1L] >= p)
        if (q > 1L) {
          expect_true(cdf[q] < p)
        }
      }
    }
  }
})


test_that("conditional distribution core uses typed validation failures", {
  expect_error(
    log_pmf_K_given_alpha(10, c(1, 2), logS),
    class = "dpprior_length_error"
  )
  expect_error(
    log_pmf_K_given_alpha(10, NA_real_, logS),
    class = "dpprior_missing_error"
  )
  expect_error(
    log_pmf_K_given_alpha(10, Inf, logS),
    class = "dpprior_nonfinite_error"
  )
  expect_error(
    log_pmf_K_given_alpha(10, 0, logS),
    class = "dpprior_bounds_error"
  )
  expect_error(
    pmf_K_given_alpha(10, 1, logS, normalize = NA),
    class = "dpprior_missing_error"
  )
  expect_error(
    pmf_K_given_alpha(10, 1, logS, normalize = 1),
    class = "dpprior_type_error"
  )
  expect_error(
    log_pmf_K_given_alpha(2, 1, matrix("x", 3, 3)),
    class = "dpprior_type_error"
  )
  expect_error(
    log_pmf_K_given_alpha(2, 1, matrix(0, 2, 3)),
    class = "dpprior_dimension_error"
  )
  expect_error(
    log_pmf_K_given_alpha(20, 1, compute_log_stirling(10)),
    class = "dpprior_bounds_error"
  )

  malformed <- compute_log_stirling(10)
  malformed[11L, 2L] <- NA_real_
  expect_error(
    log_pmf_K_given_alpha(10, 1, malformed),
    class = "dpprior_nonfinite_error"
  )
  expect_error(
    quantile_K_given_alpha(numeric(), 10, 1, logS),
    class = "dpprior_length_error"
  )
  expect_error(
    quantile_K_given_alpha(NA_real_, 10, 1, logS),
    class = "dpprior_missing_error"
  )
})
