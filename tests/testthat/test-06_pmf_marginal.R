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

  expect_equal(as.vector(q0), 1L)
  expect_equal(as.vector(q1), 50L)
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
  expect_identical(
    attr(mode_val, "marginal_metadata", exact = TRUE)$status,
    "approximate"
  )
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

test_that("PMF selected/refined discrepancy meets an explicit L1 budget", {
  checked <- pmf_K_marginal(
    50, 1.5, 0.5, logS_50,
    M = 80L, M_verify = 160L, abs_tol = 1e-10, rel_tol = 1e-7
  )
  metadata <- attr(checked, "marginal_metadata", exact = TRUE)

  expect_identical(metadata$status, "converged")
  expect_lte(
    metadata$verification$l1_difference,
    metadata$verification$tolerance
  )
  expect_true(metadata$verification$l1_passed)
  expect_true(metadata$verification$moment_consistency$passed)
  expect_true(
    metadata$verification$moment_consistency$selected$mean$passed
  )
  expect_true(
    metadata$verification$moment_consistency$selected$variance$passed
  )
  expect_true(
    metadata$verification$moment_consistency$verification$mean$passed
  )
  expect_true(
    metadata$verification$moment_consistency$verification$variance$passed
  )
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
# 11. Law of Total Variance
# -----------------------------------------------------------------------------

test_that("PMF variance agrees with the total-variance decomposition", {
  summary_marg <- summary_K_marginal(50, 1.5, 0.5, logS_50)
  mixed <- exact_K_moments(50, 1.5, 0.5)

  expect_equal(summary_marg$var, mixed$var, tolerance = 1e-10)
  expect_equal(
    summary_marg$var,
    mixed$decomposition$within_alpha + mixed$decomposition$between_alpha,
    tolerance = 1e-10
  )
  expect_gte(mixed$decomposition$between_alpha, 0)
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


# -----------------------------------------------------------------------------
# 14. Phase 4 marginal distribution contract
# -----------------------------------------------------------------------------

test_that("marginal PMF matches adaptive integration cell by cell", {
  J <- 8L
  a <- 0.5
  b <- 0.5
  logS <- compute_log_stirling(J)
  observed <- pmf_K_marginal(J, a, b, logS, M = 320L)

  reference <- c(0, vapply(seq_len(J), function(k) {
    stats::integrate(
      function(alpha) {
        conditional <- vapply(alpha, function(value) {
          pmf_K_given_alpha(J, value, logS)[k + 1L]
        }, numeric(1))
        conditional * stats::dgamma(alpha, shape = a, rate = b)
      },
      lower = 0,
      upper = Inf,
      abs.tol = 1e-11,
      rel.tol = 1e-10,
      subdivisions = 1000L,
      stop.on.error = TRUE
    )$value
  }, numeric(1)))

  expect_equal(sum(reference), 1, tolerance = 1e-10)
  expect_equal(observed, reference, tolerance = 2e-9, ignore_attr = TRUE)
})


test_that("PMF-derived moments agree with mixed moment formulas on a grid", {
  grid <- data.frame(
    J = c(2L, 10L, 50L, 100L),
    a = c(0.5, 0.5, 1.5, 2),
    b = c(2, 0.5, 0.5, 1)
  )
  logS <- compute_log_stirling(max(grid$J))

  for (i in seq_len(nrow(grid))) {
    case <- grid[i, ]
    pmf <- pmf_K_marginal(case$J, case$a, case$b, logS, M = 160L)
    k <- 0:case$J
    mean_pmf <- sum(k * pmf)
    var_pmf <- sum((k - mean_pmf)^2 * pmf)
    moments <- exact_K_moments(case$J, case$a, case$b, M = 160L)

    expect_equal(mean_pmf, moments$mean, tolerance = 2e-10)
    expect_equal(var_pmf, moments$var, tolerance = 2e-9)
  }
})


test_that("marginal PMF exposes support normalization and truncation metadata", {
  pmf <- pmf_K_marginal(50, 2, 1, logS_50, M = 80L)
  metadata <- attr(pmf, "marginal_metadata", exact = TRUE)

  expect_type(metadata, "list")
  expect_identical(metadata$status, "approximate")
  expect_identical(metadata$reason, "fixed_order_unverified")
  expect_false(metadata$verification$performed)
  expect_true(is.na(metadata$verification$passed))
  expect_identical(metadata$support, c(lower = 1L, upper = 50L))
  expect_identical(metadata$returned_support, c(lower = 0L, upper = 50L))
  expect_false(metadata$truncation$truncated)
  expect_identical(metadata$truncation$requested_mass, 1)
  expect_equal(metadata$truncation$achieved_mass, 1, tolerance = 1e-14)
  expect_identical(metadata$truncation$omitted_mass, 0)
  expect_equal(metadata$normalization$probability_sum, 1, tolerance = 1e-14)
})


test_that("higher-order PMF verification distinguishes convergence", {
  converged <- pmf_K_marginal(
    50, 2, 1, logS_50, M = 80L, M_verify = 160L,
    abs_tol = 1e-10, rel_tol = 1e-8
  )
  converged_metadata <- attr(converged, "marginal_metadata", exact = TRUE)
  expect_identical(converged_metadata$status, "converged")
  expect_true(converged_metadata$verification$passed)
  expect_identical(converged_metadata$verification$M_verification, 160L)

  approximate <- pmf_K_marginal(
    50, 0.5, 0.2, logS_50, M = 80L, M_verify = 160L,
    abs_tol = 1e-12, rel_tol = 1e-10
  )
  approximate_metadata <- attr(approximate, "marginal_metadata", exact = TRUE)
  selected <- pmf_K_marginal(
    50, 0.5, 0.2, logS_50, M = 80L
  )
  expect_identical(approximate_metadata$status, "approximate")
  expect_identical(approximate_metadata$reason, "higher_order_disagreement")
  expect_false(approximate_metadata$verification$passed)
  expect_gt(approximate_metadata$verification$l1_difference,
            approximate_metadata$verification$tolerance)
  expect_identical(as.vector(approximate), as.vector(selected))

  expect_error(
    pmf_K_marginal(
      50, 0.5, 0.2, logS_50, M = 80L, M_verify = 160L,
      abs_tol = 1e-12, rel_tol = 1e-10, strict = TRUE
    ),
    class = "dpprior_marginal_convergence_error"
  )
})


test_that("marginal PMF verification controls retain typed bounds failures", {
  expect_error(
    pmf_K_marginal(50, 2, 1, logS_50, M = 80L, M_verify = 80L),
    class = "dpprior_marginal_verification_error"
  )
  expect_error(
    pmf_K_marginal(50, 2, 1, logS_50, M = 80L, M_verify = 80L),
    class = "dpprior_bounds_error"
  )
  expect_error(
    pmf_K_marginal(50, 2, 1, logS_50, M = 80L, M_verify = 513L),
    class = "dpprior_bounds_error"
  )
  expect_error(
    pmf_K_marginal(50, 2, 1, logS_50, abs_tol = -1),
    class = "dpprior_control_error"
  )
})


test_that("invalid conditional mixtures fail instead of substituting a PMF", {
  testthat::local_mocked_bindings(
    log_pmf_K_given_alpha = function(J, alpha, logS) {
      rep(-Inf, as.integer(J) + 1L)
    },
    .package = "DPprior"
  )

  expect_error(
    pmf_K_marginal(10, 2, 1, compute_log_stirling(10)),
    class = "dpprior_marginal_pmf_error"
  )
})


test_that("PMF convergence rejects a mutually agreeing malicious mixture", {
  testthat::local_mocked_bindings(
    log_pmf_K_given_alpha = function(J, alpha, logS) {
      c(-Inf, rep(-log(as.integer(J)), as.integer(J)))
    },
    .package = "DPprior"
  )

  checked <- pmf_K_marginal(
    20, 2, 1, compute_log_stirling(20),
    M = 80L, M_verify = 160L
  )
  metadata <- attr(checked, "marginal_metadata", exact = TRUE)

  # The fabricated PMFs agree exactly across orders, but disagree with the
  # independently computed marginal moments and therefore cannot converge.
  expect_true(metadata$verification$l1_passed)
  expect_false(metadata$verification$moment_consistency$passed)
  expect_false(metadata$verification$passed)
  expect_identical(metadata$status, "approximate")
  expect_identical(metadata$reason, "pmf_moment_disagreement")
  expect_false(
    metadata$verification$moment_consistency$selected$mean$passed
  )
  expect_false(
    metadata$verification$moment_consistency$verification$variance$passed
  )

  expect_error(
    pmf_K_marginal(
      20, 2, 1, compute_log_stirling(20),
      M = 80L, M_verify = 160L, strict = TRUE
    ),
    class = "dpprior_marginal_convergence_error"
  )
})


test_that("verification PMF remains private state for downstream auditors", {
  selected <- pmf_K_marginal(
    20, 2, 1, compute_log_stirling(20),
    M = 80L, M_verify = 160L
  )
  verification <- attr(
    selected, ".marginal_verification_pmf", exact = TRUE
  )
  independent <- pmf_K_marginal(
    20, 2, 1, compute_log_stirling(20), M = 160L
  )

  expect_equal(verification, independent, tolerance = 0, ignore_attr = TRUE)
  expect_false(identical(as.vector(selected), as.vector(verification)))
})


test_that("marginal quantiles stay on support and satisfy discrete inversion", {
  J <- 50L
  probs <- c(0, .Machine$double.eps, 0.01, 0.25, 0.5, 0.9, 1)
  cdf <- cdf_K_marginal(J, 1.5, 0.5, logS_50)
  quantiles <- quantile_K_marginal(probs, J, 1.5, 0.5, logS_50)

  expect_identical(quantiles[[1L]], 1L)
  expect_identical(quantiles[[length(quantiles)]], J)
  expect_true(all(quantiles >= 1L & quantiles <= J))
  expect_true(all(diff(quantiles) >= 0L))

  for (i in seq_along(probs)) {
    p <- probs[[i]]
    q <- quantiles[[i]]
    if (p > 0) {
      expect_gte(cdf[[q + 1L]], p)
      if (q > 1L && p < 1) {
        expect_lt(cdf[[q]], p)
      }
    }
  }
})


test_that("marginal quantiles preserve probability names", {
  probabilities <- c(lower = 0, median = 0.5, upper = 1)
  quantiles <- quantile_K_marginal(
    probabilities, 50, 1.5, 0.5, logS_50
  )
  expect_identical(names(quantiles), names(probabilities))
})


test_that("marginal quantiles expose disagreement and retain selected order", {
  probability <- 0.41052028739476654
  selected <- quantile_K_marginal(
    probability, 50, 0.5, 0.2, logS_50, M = 80L
  )
  higher_order <- quantile_K_marginal(
    probability, 50, 0.5, 0.2, logS_50, M = 160L
  )
  checked <- quantile_K_marginal(
    probability, 50, 0.5, 0.2, logS_50,
    M = 80L, M_verify = 160L, abs_tol = 1e-12, rel_tol = 1e-10
  )
  metadata <- attr(checked, "marginal_metadata", exact = TRUE)

  expect_identical(as.vector(checked), as.vector(selected))
  expect_false(identical(as.vector(checked), as.vector(higher_order)))
  expect_identical(metadata$status, "approximate")
  expect_identical(metadata$reason, "higher_order_disagreement")
})


test_that("midpoint auditor catches an unstable discrete quantile", {
  J <- 20L
  logS <- compute_log_stirling(J)
  selected_cdf <- cdf_K_marginal(J, 0.1, 0.5, logS, M = 80L)
  verification_cdf <- cdf_K_marginal(J, 0.1, 0.5, logS, M = 160L)
  probability <- (selected_cdf[[2L]] + verification_cdf[[2L]]) / 2

  selected <- quantile_K_marginal(
    probability, J, 0.1, 0.5, logS, M = 80L
  )
  verification <- quantile_K_marginal(
    probability, J, 0.1, 0.5, logS, M = 160L
  )
  checked_pmf <- pmf_K_marginal(
    J, 0.1, 0.5, logS, M = 80L, M_verify = 160L
  )
  checked <- quantile_K_marginal(
    probability, J, 0.1, 0.5, logS,
    M = 80L, M_verify = 160L
  )
  metadata <- attr(checked, "marginal_metadata", exact = TRUE)

  expect_equal(probability, 0.81527353170979544, tolerance = 1e-15)
  expect_identical(as.vector(selected), 2L)
  expect_identical(as.vector(verification), 1L)
  expect_identical(
    attr(checked_pmf, "marginal_metadata", exact = TRUE)$status,
    "converged"
  )
  expect_identical(as.vector(checked), as.vector(selected))
  expect_identical(metadata$status, "approximate")
  expect_identical(metadata$reason, "discrete_quantile_disagreement")
  expect_false(metadata$discrete_verification$quantile$passed)
  expect_identical(
    metadata$discrete_verification$quantile$verification, 1L
  )

  expect_error(
    quantile_K_marginal(
      probability, J, 0.1, 0.5, logS,
      M = 80L, M_verify = 160L, strict = TRUE
    ),
    class = "dpprior_marginal_convergence_error"
  )
})


test_that("summary audits quantiles median and mode without substitution", {
  stable <- summary_K_marginal(
    50, 2, 1, logS_50, M = 80L, M_verify = 160L
  )
  expect_identical(stable$metadata$status, "converged")
  expect_true(stable$metadata$discrete_verification$quantile$passed)
  expect_true(stable$metadata$discrete_verification$median$passed)
  expect_true(stable$metadata$discrete_verification$mode$passed)

  # With a deliberately loose distribution tolerance, the continuous PMF
  # contract passes but the median still changes from 2 to 1.
  median_case <- summary_K_marginal(
    10, 0.1, 0.02, compute_log_stirling(10),
    M = 20L, M_verify = 60L, abs_tol = 0.5
  )
  expect_identical(median_case$median, 2L)
  expect_identical(median_case$metadata$status, "approximate")
  expect_false(median_case$metadata$discrete_verification$median$passed)
  expect_identical(
    median_case$metadata$discrete_verification$median$verification, 1L
  )

  # A separate grid cell changes only the mode among the primary summaries.
  mode_case <- mode_K_marginal(
    100, 1, 0.1, logS_100,
    M = 80L, M_verify = 160L, abs_tol = 0.01
  )
  mode_metadata <- attr(mode_case, "marginal_metadata", exact = TRUE)
  expect_identical(as.vector(mode_case), 11L)
  expect_identical(mode_metadata$status, "approximate")
  expect_false(mode_metadata$discrete_verification$mode$passed)
  expect_identical(mode_metadata$discrete_verification$mode$verification, 13L)

  expect_error(
    summary_K_marginal(
      10, 0.1, 0.02, compute_log_stirling(10),
      M = 20L, M_verify = 60L, abs_tol = 0.5, strict = TRUE
    ),
    class = "dpprior_marginal_convergence_error"
  )
  expect_error(
    mode_K_marginal(
      100, 1, 0.1, logS_100,
      M = 80L, M_verify = 160L, abs_tol = 0.01, strict = TRUE
    ),
    class = "dpprior_marginal_convergence_error"
  )
})


test_that("CDF has exact endpoints and carries marginal metadata", {
  cdf <- cdf_K_marginal(
    50, 2, 1, logS_50, M = 80L, M_verify = 160L
  )
  metadata <- attr(cdf, "marginal_metadata", exact = TRUE)
  expect_identical(cdf[[1L]], 0)
  expect_identical(cdf[[length(cdf)]], 1)
  expect_true(all(diff(cdf) >= 0))
  expect_identical(metadata$status, "converged")
  expect_identical(metadata$support, c(lower = 1L, upper = 50L))

  # Regression: the unscaled cumulative sum can end at 1 + 2^-52 here.
  # Replacing only its final value by one used to create a tiny decrease.
  rounding_cdf <- cdf_K_marginal(
    100, 1, 1, compute_log_stirling(100),
    M = 80L, M_verify = 160L
  )
  expect_identical(rounding_cdf[[1L]], 0)
  expect_identical(rounding_cdf[[length(rounding_cdf)]], 1)
  expect_true(all(diff(rounding_cdf) >= 0))
})


test_that("marginal summary reports achieved quantile mass and no truncation", {
  probs <- c(0, 0.05, 0.5, 0.95, 1)
  result <- summary_K_marginal(
    50, 2, 1, logS_50, M = 80L, probs = probs,
    M_verify = 160L
  )

  expect_true("metadata" %in% names(result))
  expect_identical(result$metadata$status, "converged")
  expect_false(result$metadata$truncation$truncated)
  expect_identical(result$metadata$truncation$requested_mass, 1)
  expect_equal(result$metadata$truncation$achieved_mass, 1, tolerance = 1e-14)
  expect_identical(result$metadata$truncation$omitted_mass, 0)
  expect_equal(result$metadata$quantiles$probabilities, probs)
  expect_equal(result$metadata$quantiles$values, unname(result$quantiles))
  expect_true(all(result$metadata$quantiles$achieved_cdf >= probs))
  expect_identical(result$metadata$quantiles$previous_cdf[[1L]], 0)
  interior <- probs > 0 & probs < 1
  expect_true(all(
    result$metadata$quantiles$previous_cdf[interior] < probs[interior]
  ))
  expect_lte(result$metadata$quantiles$previous_cdf[[length(probs)]], 1)
  expect_identical(
    attr(result$pmf, "marginal_metadata", exact = TRUE), result$metadata
  )
  expect_identical(
    attr(result$cdf, "marginal_metadata", exact = TRUE), result$metadata
  )
  expect_true(all(result$mode >= 1L & result$mode <= 50L))
  expect_true(all(result$median >= 1L & result$median <= 50L))
})


test_that("J=1 marginal distribution is deterministic on support one", {
  logS <- compute_log_stirling(1)
  pmf <- pmf_K_marginal(1, 0.2, 4, logS)
  cdf <- cdf_K_marginal(1, 0.2, 4, logS)
  quantiles <- quantile_K_marginal(c(0, 0.5, 1), 1, 0.2, 4, logS)
  summary <- summary_K_marginal(1, 0.2, 4, logS, probs = c(0, 0.5, 1))

  expect_identical(as.vector(pmf), c(0, 1))
  expect_identical(as.vector(cdf), c(0, 1))
  expect_identical(as.vector(quantiles), rep(1L, 3L))
  expect_identical(summary$mean, 1)
  expect_identical(summary$var, 0)
  expect_identical(summary$mode, 1L)
  expect_identical(summary$median, 1L)
})
