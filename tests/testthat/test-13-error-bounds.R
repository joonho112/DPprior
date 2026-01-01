# =============================================================================
# Unit Tests for Module 13: Error Bounds (Revised)
# =============================================================================
#
# Tests for error quantification functions using the correct Chen-Stein
# prefactor (1 - exp(-lambda))/lambda.
#
# Author: JoonHo Lee
# Date: December 2025
# =============================================================================

# Note: context() is deprecated in testthat 3rd edition

# =============================================================================
# Golden Values (from Python verification with correct prefactor)
# =============================================================================

golden_50_2_1 <- list(
  J = 50L, a = 2, b = 1,
  exact_mean = 6.6396928911,
  exact_var = 12.9545022869,
  a1_mean = 8.8240460109,
  a1_var = 38.4318940009,
  pois_raw_at_E_alpha = 1.5020688045,
  pois_bound_at_E_alpha = 0.2481907538,
  lin_bound_at_E_alpha = 0.3328077886,
  total_tv_at_E_alpha = 0.5809985424,
  marginal_tv_bound = 0.5213110800,
  error_mean_rel = 32.898406,
  error_var_rel = 196.668241
)

golden_100_1_1 <- list(
  J = 100L, a = 1, b = 1,
  exact_mean = 4.8373970578,
  exact_var = 13.2154435967,
  a1_mean = 5.6051701860,
  a1_var = 25.8127626279,
  pois_raw_at_E_alpha = 0.6349839002,
  pois_bound_at_E_alpha = 0.1493395318,
  lin_bound_at_E_alpha = 0.0988742108,
  total_tv_at_E_alpha = 0.2482137426,
  marginal_tv_bound = 0.2465068320,
  error_mean_rel = 15.871617,
  error_var_rel = 95.322710
)

golden_50_1_1 <- list(
  J = 50L, a = 1, b = 1,
  exact_mean = 4.1586595424,
  exact_var = 8.8287378070,
  a1_mean = 4.9120230054,
  a1_var = 19.2159470004,
  pois_raw_at_E_alpha = 0.6251327336,
  pois_bound_at_E_alpha = 0.1732508656,
  lin_bound_at_E_alpha = 0.1062796064,
  total_tv_at_E_alpha = 0.2795304720,
  marginal_tv_bound = 0.2728911762,
  error_mean_rel = 18.115536,
  error_var_rel = 117.652256
)

golden_25_2_2 <- list(
  J = 25L, a = 2, b = 2,
  exact_mean = 3.6366862654,
  exact_var = 4.0830325590,
  a1_mean = 4.2188758249,
  a1_var = 8.3994566128,
  pois_raw_at_E_alpha = 0.6057234036,
  pois_bound_at_E_alpha = 0.2022304413,
  lin_bound_at_E_alpha = 0.1147621810,
  total_tv_at_E_alpha = 0.3169926223,
  marginal_tv_bound = 0.3112972680,
  error_mean_rel = 16.008793,
  error_var_rel = 105.716131
)


# =============================================================================
# Test: Poissonization Bound Properties
# =============================================================================

test_that("Poissonization bound with Chen-Stein prefactor is in [0, 1]", {
  for (J in c(25, 50, 100, 200)) {
    for (alpha in c(0.5, 1, 2, 5, 10)) {
      bound <- compute_poissonization_bound(J, alpha)
      expect_true(is.finite(bound))
      expect_true(bound >= 0)
      expect_true(bound <= 1)
    }
  }
})

test_that("Raw poissonization bound (sum_p_sq) >= Chen-Stein bound", {
  for (J in c(25, 50, 100)) {
    for (alpha in c(0.5, 1, 2, 5)) {
      raw <- compute_poissonization_bound(J, alpha, raw = TRUE)
      full <- compute_poissonization_bound(J, alpha, raw = FALSE)
      # Chen-Stein prefactor is always <= 1, so full <= raw (before capping)
      expect_true(raw >= full)
    }
  }
})

test_that("Chen-Stein prefactor is in (0, 1]", {
  # Test the prefactor behavior implicitly
  J <- 50
  alpha <- 2

  lambda <- mean_K_given_alpha(J, alpha) - 1
  raw <- compute_poissonization_bound(J, alpha, raw = TRUE)
  full <- compute_poissonization_bound(J, alpha, raw = FALSE)

  # Implied prefactor
  prefactor <- full / raw
  expect_gt(prefactor, 0)
  expect_lte(prefactor, 1)
})


# =============================================================================
# Test: Linearization Bound Edge Cases
# =============================================================================

test_that("Linearization bound handles edge cases correctly", {
  J <- 50

  # Normal case: bound in (0, 1)
  bound_normal <- compute_linearization_bound(J, alpha = 2)
  expect_gt(bound_normal, 0)
  expect_lt(bound_normal, 1)

  # When cJ = 0 and lambda > 0: bound = 1 (infinite KL)
  bound_cj0 <- compute_linearization_bound(J, alpha = 2, cJ = 0)
  expect_equal(bound_cj0, 1)

  # Small alpha: bound should still be valid
  bound_small <- compute_linearization_bound(J, alpha = 0.01)
  expect_true(is.finite(bound_small))
  expect_gte(bound_small, 0)
  expect_lte(bound_small, 1)
})

test_that("Linearization bound decreases with J (for fixed alpha)", {
  alpha <- 2
  bounds <- sapply(c(25, 50, 100, 200, 500), function(J) {
    compute_linearization_bound(J, alpha)
  })

  # Should generally decrease
  expect_lt(bounds[5], bounds[1])
  expect_lt(bounds[4], bounds[2])
})

test_that("Linearization bound validates cJ input", {
  expect_error(compute_linearization_bound(50, 1, cJ = -1))
  expect_error(compute_linearization_bound(50, 1, cJ = NA))
  expect_error(compute_linearization_bound(50, 1, cJ = c(1, 2)))
})


# =============================================================================
# Test: Total TV Bound
# =============================================================================

test_that("Total TV bound equals sum of components (capped at 1)", {
  J <- 50
  alpha <- 2
  cJ <- log(J)

  B_pois <- compute_poissonization_bound(J, alpha)
  B_lin <- compute_linearization_bound(J, alpha, cJ)
  B_total <- compute_total_tv_bound(J, alpha, cJ)

  expect_equal(B_total, min(1, B_pois + B_lin))
})

test_that("Total TV bound is always in [0, 1]", {
  for (J in c(10, 50, 200)) {
    for (alpha in c(0.5, 2, 10)) {
      bound <- compute_total_tv_bound(J, alpha)
      expect_gte(bound, 0)
      expect_lte(bound, 1)
    }
  }
})


# =============================================================================
# Test: A1 Moment Error
# =============================================================================

test_that("A1 moment error function returns complete results", {
  errors <- a1_moment_error(50, 2, 1)

  expect_true(all(c("exact_mean", "exact_var", "a1_mean", "a1_var",
                    "error_mean_abs", "error_var_abs",
                    "error_mean_rel", "error_var_rel") %in% names(errors)))

  # All values should be positive
  expect_true(all(sapply(errors, function(x) x > 0)))
})

test_that("A1 moment error improves with J", {
  a <- 2
  b <- 1

  errors_50 <- a1_moment_error(50, a, b)
  errors_200 <- a1_moment_error(200, a, b)

  expect_lt(errors_200$error_mean_rel, errors_50$error_mean_rel)
})

test_that("A1 moment error accepts cJ parameter", {
  J <- 50
  a <- 2
  b <- 1

  errors_log <- a1_moment_error(J, a, b, cJ = log(J))
  errors_H <- a1_moment_error(J, a, b, cJ = digamma(J) + 0.5772)

  # Results should differ slightly
  expect_false(identical(errors_log$a1_mean, errors_H$a1_mean))
})


# =============================================================================
# Test: Expected TV Bound
# =============================================================================

test_that("Expected TV bound is in [0, 1]", {
  for (J in c(25, 50, 100)) {
    for (a in c(1, 2)) {
      bound <- expected_tv_bound(J, a, b = 1)
      expect_true(is.finite(bound))
      expect_gte(bound, 0)
      expect_lte(bound, 1)
    }
  }
})

test_that("Expected TV bound decreases with J", {
  bounds <- sapply(c(25, 50, 100, 200), function(J) {
    expected_tv_bound(J, a = 2, b = 1)
  })

  # Should generally decrease
  expect_lt(bounds[4], bounds[1])
})


# =============================================================================
# Test: DPprior_error_bounds Main Function
# =============================================================================

test_that("DPprior_error_bounds returns correct structure", {
  bounds <- DPprior_error_bounds(50, 2, 1)

  expect_s3_class(bounds, "DPprior_error_bounds")
  expect_true(all(c("J", "a", "b", "cJ", "moment_errors", "tv_bounds",
                    "recommendation", "threshold_J") %in% names(bounds)))

  expect_true(bounds$recommendation %in% c("A1_sufficient", "A2_recommended"))
})

test_that("DPprior_error_bounds print method works", {
  bounds <- DPprior_error_bounds(50, 1.6, 1.2)

  expect_output(print(bounds), "DPprior A1 Approximation Error Analysis")
  expect_output(print(bounds), "Moment Errors")
  expect_output(print(bounds), "Recommendation")
})

test_that("DPprior_error_bounds summary method works", {
  bounds <- DPprior_error_bounds(50, 1.6, 1.2)

  expect_output(summary(bounds), "Conditional TV Bounds")
})


# =============================================================================
# Test: Golden Value Verification
# =============================================================================

test_that("Golden values: J=50, a=2, b=1 match Python reference", {
  g <- golden_50_2_1
  errors <- a1_moment_error(g$J, g$a, g$b)

  expect_equal(errors$exact_mean, g$exact_mean, tolerance = 1e-6)
  expect_equal(errors$a1_mean, g$a1_mean, tolerance = 1e-6)

  alpha <- g$a / g$b
  pois_raw <- compute_poissonization_bound(g$J, alpha, raw = TRUE)
  pois_bound <- compute_poissonization_bound(g$J, alpha, raw = FALSE)

  expect_equal(pois_raw, g$pois_raw_at_E_alpha, tolerance = 1e-6)
  expect_equal(pois_bound, g$pois_bound_at_E_alpha, tolerance = 1e-6)
})

test_that("Golden values: J=100, a=1, b=1 match Python reference", {
  g <- golden_100_1_1
  errors <- a1_moment_error(g$J, g$a, g$b)

  expect_equal(errors$exact_mean, g$exact_mean, tolerance = 1e-6)

  total_tv <- compute_total_tv_bound(g$J, g$a / g$b, log(g$J))
  expect_equal(total_tv, g$total_tv_at_E_alpha, tolerance = 1e-4)
})

test_that("Golden values: Marginal TV bounds match Python reference", {
  g <- golden_50_1_1
  marginal <- expected_tv_bound(g$J, g$a, g$b)

  expect_equal(marginal, g$marginal_tv_bound, tolerance = 0.01)
})


# =============================================================================
# Test: Vectorization
# =============================================================================

test_that("Poissonization bound is vectorized over alpha", {
  J <- 50
  alpha_vec <- c(0.5, 1, 2, 5)

  bounds <- compute_poissonization_bound(J, alpha_vec)

  expect_length(bounds, 4)
  expect_true(all(bounds >= 0 & bounds <= 1))
})

test_that("Linearization bound is vectorized over alpha", {
  J <- 50
  alpha_vec <- c(0.5, 1, 2, 5)

  bounds <- compute_linearization_bound(J, alpha_vec)

  expect_length(bounds, 4)
  expect_true(all(bounds >= 0 & bounds <= 1))
})


# =============================================================================
# Test: compute_error_landscape
# =============================================================================

test_that("compute_error_landscape returns correct structure", {
  landscape <- compute_error_landscape(
    J_seq = c(25, 50),
    alpha_seq = c(1, 2)
  )

  expect_equal(nrow(landscape), 4)  # 2 J x 2 alpha
  expect_true(all(c("J", "alpha", "lambda_exact", "lambda_approx",
                    "pois_raw", "pois_bound", "lin_bound", "total_tv")
                  %in% names(landscape)))
})


# =============================================================================
# Test: Input Validation
# =============================================================================

test_that("Functions reject invalid inputs", {
  # Invalid J
  expect_error(compute_poissonization_bound(0, 1))
  expect_error(compute_poissonization_bound(-5, 1))

  # Invalid alpha
  expect_error(compute_poissonization_bound(50, -1))
  expect_error(compute_poissonization_bound(50, 0))

  # Invalid a, b
  expect_error(a1_moment_error(50, -1, 1))
  expect_error(a1_moment_error(50, 1, 0))
})
