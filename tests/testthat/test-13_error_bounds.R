# =============================================================================
# Unit Tests for Module 13: Error Bounds (Revised)
# =============================================================================
#
# Tests for error quantification functions using the correct Chen-Stein
# prefactor (1 - exp(-lambda))/lambda.
#
# Author: JoonHo Lee (jlee296@ua.edu)
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

.test_poisson_binomial_pmf <- function(probabilities) {
  pmf <- 1
  for (probability in probabilities) {
    pmf <- c(pmf * (1 - probability), 0) +
      c(0, pmf * probability)
  }
  pmf
}

.test_tv_finite_vs_poisson <- function(finite_pmf, lambda) {
  support <- 0:(length(finite_pmf) - 1L)
  poisson_mass <- stats::dpois(support, lambda)
  poisson_tail <- stats::ppois(
    max(support), lambda, lower.tail = FALSE
  )
  0.5 * (sum(abs(finite_pmf - poisson_mass)) + poisson_tail)
}

.test_tv_poisson <- function(lambda, lambda_prime) {
  upper <- max(
    stats::qpois(1 - 1e-14, lambda),
    stats::qpois(1 - 1e-14, lambda_prime),
    20
  )
  support <- 0:upper
  mass_difference <- abs(
    stats::dpois(support, lambda) -
      stats::dpois(support, lambda_prime)
  )
  tail_difference <- abs(
    stats::ppois(upper, lambda, lower.tail = FALSE) -
      stats::ppois(upper, lambda_prime, lower.tail = FALSE)
  )
  0.5 * (sum(mass_difference) + tail_difference)
}


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

test_that("poissonization theorem crosswalk and empirical TV are noncontradictory", {
  grid <- expand.grid(
    J = c(2L, 10L, 50L),
    alpha = c(0.01, 0.1, 1, 10, 100)
  )

  for (i in seq_len(nrow(grid))) {
    J <- grid$J[[i]]
    alpha <- grid$alpha[[i]]
    probabilities <- alpha / (alpha + seq_len(J - 1L))
    lambda <- sum(probabilities)
    raw <- sum(probabilities^2)
    exact_pmf <- .test_poisson_binomial_pmf(probabilities)
    actual_tv <- .test_tv_finite_vs_poisson(exact_pmf, lambda)
    implemented <- compute_poissonization_bound(J, alpha)
    manuscript_d10 <- min(1, 1 / lambda) * raw

    expect_lte(actual_tv, implemented + 2e-13)
    expect_lte(implemented, manuscript_d10 + 2e-13)
    expect_equal(
      DPprior:::compute_sum_p_squared(J, alpha),
      raw,
      tolerance = 1e-14
    )
  }
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

test_that("Poisson KL is stable near equality and linearization bounds exact TV", {
  lambda <- c(0.01, 1, 100)
  lambda_prime <- lambda * (1 + 1e-10)
  kl <- DPprior:::poisson_kl_divergence(lambda, lambda_prime)
  second_order <- lambda * (1e-10)^2 / 2

  expect_true(all(is.finite(kl)))
  expect_true(all(kl >= 0))
  expect_equal(kl, second_order, tolerance = 1e-7)

  for (J in c(2L, 10L, 50L)) {
    for (alpha in c(0.01, 0.1, 1, 10)) {
      exact_lambda <- sum(alpha / (alpha + seq_len(J - 1L)))
      proxy_lambda <- alpha * log(J)
      actual_tv <- .test_tv_poisson(exact_lambda, proxy_lambda)
      expect_lte(
        actual_tv,
        compute_linearization_bound(J, alpha) + 2e-12
      )
    }
  }
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

test_that("conditional total bound contains empirical exact discrepancy", {
  for (J in c(2L, 10L, 50L)) {
    for (alpha in c(0.01, 0.1, 1, 10, 100)) {
      probabilities <- alpha / (alpha + seq_len(J - 1L))
      exact_pmf <- .test_poisson_binomial_pmf(probabilities)
      actual_tv <- .test_tv_finite_vs_poisson(
        exact_pmf, alpha * log(J)
      )
      poissonization <- compute_poissonization_bound(J, alpha)
      linearization <- compute_linearization_bound(J, alpha)
      total <- compute_total_tv_bound(J, alpha)

      expect_lte(actual_tv, total + 2e-12)
      expect_gte(total + 2e-14, poissonization)
      expect_gte(total + 2e-14, linearization)
    }
  }
})

test_that("all TV bounds remain finite at accepted extreme alpha values", {
  alpha <- c(1e-300, 1e-200, 1e-16, 1, 1e100, 1e308)
  for (J in c(1L, 2L, 50L, 500L)) {
    poissonization <- compute_poissonization_bound(J, alpha)
    linearization <- compute_linearization_bound(J, alpha)
    total <- compute_total_tv_bound(J, alpha)
    bounds <- c(poissonization, linearization, total)

    expect_true(all(is.finite(bounds)))
    expect_true(all(bounds >= 0 & bounds <= 1))
    expect_true(all(total + 2e-14 >= poissonization))
    expect_true(all(total + 2e-14 >= linearization))
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

test_that("A1 proxy moments avoid finite underflow and overflow intermediates", {
  J <- 50L
  cJ <- log(J)

  # b^2 overflows, so the legacy product/quotient returned zero although the
  # variance is positive and representable.
  a <- 100
  b <- 1e160
  expect_identical(a * cJ * (b + cJ) / (b * b), 0)
  reference <- exp(log(a) + log(cJ) + log(b + cJ) - 2 * log(b))
  errors <- a1_moment_error(J, a, b)
  expect_equal(errors$a1_var, reference, tolerance = 1e-13)
  expect_gt(errors$a1_var, 0)

  # The proxy calculation itself also remains finite when b^2 underflows.
  tiny_shape <- DPprior:::.a1_proxy_moments(1e-300, 1e-170, cJ)
  tiny_shape_reference <- exp(
    log(1e-300) + log(cJ) + log(cJ + 1e-170) - 2 * log(1e-170)
  )
  expect_true(is.finite(tiny_shape$variance))
  expect_equal(tiny_shape$variance, tiny_shape_reference, tolerance = 1e-13)
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

test_that("expected TV bound verifies adaptively and exposes GL disagreement", {
  bound <- expected_tv_bound(
    J = 50, a = 2, b = 1, M = 80, M_verify = 160
  )
  metadata <- attr(bound, "tv_bound_metadata", exact = TRUE)

  expect_true(is.finite(bound))
  expect_gte(as.numeric(bound), 0)
  expect_lte(as.numeric(bound), 1)
  expect_identical(metadata$status, "converged")
  expect_true(metadata$verified)
  expect_identical(metadata$reason, "adaptive_integration_agreement")
  expect_equal(metadata$M_selected, 80L)
  expect_equal(metadata$M_verification, 160L)
  expect_equal(as.numeric(bound), metadata$selected)
  expect_gte(as.numeric(bound), metadata$selected_estimate)
  expect_gt(metadata$selected_numerical_error_bound, 0)
  expect_identical(metadata$quadrature_audit$status, "approximate")
  expect_false(metadata$quadrature_audit$verified)
  expect_identical(
    metadata$quadrature_audit$reason,
    "higher_order_disagreement"
  )

  strict_bound <- expected_tv_bound(
    J = 50, a = 2, b = 1, M = 80, M_verify = 160,
    strict = TRUE
  )
  expect_equal(as.numeric(strict_bound), as.numeric(bound), tolerance = 1e-12)

  expect_error(
    expected_tv_bound(
      J = 50, a = 2, b = 1, M = 80, M_verify = 160,
      abs_tol = 0, rel_tol = 0, strict = TRUE
    ),
    class = "dpprior_tv_bound_convergence_error"
  )
})

test_that("marginal TV bound does not contradict exact PMF discrepancy", {
  cases <- list(
    c(J = 10, a = 0.5, b = 1),
    c(J = 25, a = 2, b = 1),
    c(J = 50, a = 1, b = 1)
  )
  for (case in cases) {
    J <- as.integer(case[["J"]])
    a <- case[["a"]]
    b <- case[["b"]]
    exact <- pmf_K_marginal(
      J, a, b, compute_log_stirling(J), M = 160
    )[-1L]
    probability <- b / (b + log(J))
    proxy <- stats::dnbinom(0:(J - 1L), size = a, prob = probability)
    proxy_tail <- stats::pnbinom(
      J - 1L, size = a, prob = probability, lower.tail = FALSE
    )
    actual_tv <- 0.5 * (sum(abs(exact - proxy)) + proxy_tail)
    bound <- expected_tv_bound(J, a, b, M = 80, M_verify = 160)

    expect_lte(actual_tv, as.numeric(bound) + 2e-10)
  }
})


# =============================================================================
# Test: DPprior_error_bounds Main Function
# =============================================================================

test_that("A1 threshold scan distinguishes verified and provisional results", {
  provisional <- DPprior:::find_a1_threshold_J(
    2, 1, J_min = 10, J_max = 30, step = 10
  )
  provisional_metadata <- attr(
    provisional, "threshold_metadata", exact = TRUE
  )
  expect_identical(provisional_metadata$status, "approximate")
  expect_false(provisional_metadata$verified)
  expect_true(provisional_metadata$candidate_is_provisional)
  expect_true(all(
    provisional_metadata$scan$exact_moment_status == "approximate"
  ))

  verified <- DPprior:::find_a1_threshold_J(
    2, 1, J_min = 10, J_max = 30, step = 10,
    M = 80, M_verify = 160
  )
  verified_metadata <- attr(
    verified, "threshold_metadata", exact = TRUE
  )
  expect_identical(verified_metadata$status, "converged")
  expect_true(verified_metadata$verified)
  expect_false(verified_metadata$candidate_is_provisional)
  expect_true(all(
    verified_metadata$scan$exact_moment_status == "converged"
  ))
  expect_identical(as.integer(provisional), as.integer(verified))
})

test_that("DPprior_error_bounds returns correct structure", {
  bounds <- DPprior_error_bounds(50, 2, 1)

  expect_s3_class(bounds, "DPprior_error_bounds")
  expect_true(all(c("J", "a", "b", "cJ", "moment_errors", "tv_bounds",
                    "recommendation", "threshold_J", "status", "verified",
                    "usable", "verification", "provenance") %in% names(bounds)))

  expect_true(bounds$recommendation %in% c("A1_sufficient", "A2_recommended"))
  expect_identical(
    bounds$verified,
    identical(bounds$status, "converged")
  )
  expect_true(bounds$status %in% c("converged", "approximate"))
  expect_identical(
    bounds$verification$component_status[["adequacy_threshold"]],
    bounds$verification$adequacy_threshold$status
  )
  expect_true(bounds$verification$adequacy_threshold$verified)
  expect_true(
    bounds$provenance$theorem_crosswalk$
      mixing_is_not_an_additional_approximation
  )
})

test_that("DPprior_error_bounds print method works", {
  bounds <- DPprior_error_bounds(50, 1.6, 1.2)

  expect_output(print(bounds), "DPprior A1 Approximation Error Analysis")
  expect_output(print(bounds), "Status:")
  expect_output(print(bounds), "Verified:")
  expect_output(print(bounds), "Method used:")
  expect_output(print(bounds), "Moment Errors")
  expect_output(print(bounds), "Recommendation")
  expect_output(print(bounds), "verified scan-grid|No scanned J")
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

  expect_equal(as.numeric(marginal), g$marginal_tv_bound, tolerance = 0.01)
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
