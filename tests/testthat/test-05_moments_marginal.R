# =============================================================================
# Test File: Module 05 - Marginal Moments
# =============================================================================
#
# testthat tests for exact_K_moments() and related functions.
# Tolerances match Python reference precision.
#
# Author: JoonHo Lee (jlee296@ua.edu)
# Date: December 2025
# Part of: DPprior R Package
# =============================================================================

# Note: context() removed - deprecated in testthat 3rd edition

# =============================================================================
# Golden Test Data (Python-verified via Gauss-Laguerre quadrature M=80)
# =============================================================================
# High-precision values verified by:
# 1. Python reference implementation (verify_moments_marginal.py)
# 2. Monte Carlo cross-validation (< 1% relative error)
#
# Note: RN-01 Appendix B contains different values (10.23, 18.76 for first row).
# Those are preserved in golden_moments_marg_rn01_table.csv for traceability,
# but the exact digamma/trigamma formulas yield the values below.

golden_data <- data.frame(
  J = c(50, 50, 100, 10, 10, 100, 200, 300),
  a = c(1.5, 2.0, 2.0, 0.5, 5.0, 1.0, 2.0, 1.5),
  b = c(0.5, 1.0, 1.0, 0.5, 2.0, 0.5, 1.0, 0.5),
  mean_K = c(8.35548676, 6.63969289, 7.97855126, 2.47193732,
             4.31096274, 7.64046385, 9.34050887, 13.52153920),
  var_K = c(22.76895012, 12.95450229, 20.40684333, 2.90510882,
            2.51257669, 31.61663148, 29.87384434, 77.48541362)
)


# =============================================================================
# Test: Golden Test Data Matching (Tight Tolerance)
# =============================================================================

test_that("Marginal moments match Python-verified golden data (tol = 1e-6)", {
  for (i in seq_len(nrow(golden_data))) {
    J <- golden_data$J[i]
    a <- golden_data$a[i]
    b <- golden_data$b[i]
    expected_mean <- golden_data$mean_K[i]
    expected_var <- golden_data$var_K[i]

    result <- exact_K_moments(J, a, b, M = 80)

    expect_equal(
      result$mean, expected_mean,
      tolerance = 1e-6,
      info = sprintf("Mean mismatch at J=%d, a=%.1f, b=%.1f", J, a, b)
    )
    expect_equal(
      result$var, expected_var,
      tolerance = 1e-6,
      info = sprintf("Variance mismatch at J=%d, a=%.1f, b=%.1f", J, a, b)
    )
  }
})


# =============================================================================
# Test: Independent Reference Match (Canonical Case)
# =============================================================================

test_that("independent reference values match exactly (J=50, a=1.5, b=0.5)", {
  # Reference values computed independently with high-order quadrature (M=80).
  moments <- exact_K_moments(50, 1.5, 0.5, M = 80)

  # Match to 8 decimal places (independent-reference precision).
  expect_equal(moments$mean, 8.3554867566, tolerance = 1e-8)
  expect_equal(moments$var, 22.7689501247, tolerance = 1e-8)
})


# =============================================================================
# Test: Mean Bounds Property
# =============================================================================

test_that("Mean is bounded by 1 and J", {
  test_cases <- expand.grid(
    J = c(10, 50, 100),
    a = c(0.5, 1.0, 2.0),
    b = c(0.5, 1.0, 2.0)
  )

  for (i in seq_len(nrow(test_cases))) {
    J <- test_cases$J[i]
    a <- test_cases$a[i]
    b <- test_cases$b[i]

    moments <- exact_K_moments(J, a, b)

    expect_true(
      moments$mean >= 1 && moments$mean <= J,
      info = sprintf("Mean bounds violated at J=%d, a=%.1f, b=%.1f: E[K]=%.4f",
                     J, a, b, moments$mean)
    )
  }
})


# =============================================================================
# Test: Law of Total Variance
# =============================================================================

test_that("Marginal variance obeys the law of total variance", {
  test_cases <- list(
    list(J = 50, a = 2.0, b = 1.0),
    list(J = 100, a = 1.5, b = 0.5),
    list(J = 50, a = 0.5, b = 0.5),
    list(J = 2, a = 0.1, b = 0.177828)
  )

  for (tc in test_cases) {
    marginal <- exact_K_moments(tc$J, tc$a, tc$b, M = 80)
    within <- marginal$decomposition$within_alpha
    between <- marginal$decomposition$between_alpha

    expect_equal(
      marginal$var, within + between, tolerance = 1e-13,
      info = sprintf("total-variance identity at J=%d", tc$J)
    )
    expect_gte(within, 0)
    expect_gte(between, 0)
    expect_gte(marginal$var, within)
  }
})


# =============================================================================
# Test: Non-negative Variance
# =============================================================================

test_that("Variance is always non-negative", {
  test_cases <- expand.grid(
    J = c(5, 10, 50, 100),
    a = c(0.1, 0.5, 1.0, 5.0),
    b = c(0.1, 0.5, 1.0, 5.0)
  )

  for (i in seq_len(nrow(test_cases))) {
    moments <- exact_K_moments(test_cases$J[i], test_cases$a[i], test_cases$b[i])
    expect_true(moments$var >= 0)
  }
})


# =============================================================================
# Test: Convenience Wrapper
# =============================================================================

test_that("K_moments returns correct structure", {
  result <- K_moments(50, 2.0, 1.0)

  expect_true(is.numeric(result))
  expect_equal(length(result), 2)
  expect_true("mean" %in% names(result))
  expect_true("var" %in% names(result))
})

test_that("K_moments matches exact_K_moments", {
  J <- 50; a <- 1.5; b <- 0.5

  full <- exact_K_moments(J, a, b)
  quick <- K_moments(J, a, b)

  # Use unname() to compare values without name attributes
  expect_equal(unname(quick["mean"]), full$mean, tolerance = 1e-10)
  expect_equal(unname(quick["var"]), full$var, tolerance = 1e-10)
})

test_that("K_moments propagates approximation and verification status", {
  unverified <- K_moments(50, 2, 1)
  unverified_metadata <- attr(
    unverified, "marginal_metadata", exact = TRUE
  )
  expect_identical(unverified_metadata$status, "approximate")
  expect_identical(
    unverified_metadata$quadrature$reason, "fixed_order_unverified"
  )

  verified <- K_moments(50, 2, 1, M = 80L, M_verify = 160L)
  verified_metadata <- attr(verified, "marginal_metadata", exact = TRUE)
  expect_identical(verified_metadata$status, "converged")
  expect_true(verified_metadata$quadrature$verification_passed)
})


# =============================================================================
# Test: Quadrature Convergence
# =============================================================================

test_that("Moment selected/refined discrepancy meets an explicit budget", {
  checked <- exact_K_moments(
    50, 1.5, 0.5, M = 80L, M_verify = 160L,
    abs_tol = 1e-10, rel_tol = 1e-8
  )

  expect_identical(checked$status, "converged")
  expect_lte(
    checked$quadrature$mean_difference,
    checked$quadrature$mean_tolerance
  )
  expect_lte(
    checked$quadrature$variance_difference,
    checked$quadrature$variance_tolerance
  )
})


# =============================================================================
# Test: Jacobian Computation
# =============================================================================

test_that("Jacobian has correct structure", {
  result <- marginal_moments_with_jacobian(50, 2.0, 1.0)

  expect_true("jacobian" %in% names(result))
  expect_true(is.matrix(result$jacobian))
  expect_equal(dim(result$jacobian), c(2, 2))
  expect_equal(rownames(result$jacobian), c("mean", "var"))
  expect_equal(colnames(result$jacobian), c("a", "b"))
})

test_that("Jacobian values are finite", {
  result <- marginal_moments_with_jacobian(50, 2.0, 1.0)
  expect_true(all(is.finite(result$jacobian)))
})

test_that("legacy marginal Jacobian helper delegates to canonical contract", {
  cases <- data.frame(
    J = c(1L, 50L, 100L),
    a = c(0.2, 2, 0.1),
    b = c(4, 1, 0.2)
  )

  for (i in seq_len(nrow(cases))) {
    case <- cases[i, ]
    adapter <- marginal_moments_with_jacobian(
      case$J, case$a, case$b, M = 80L
    )
    canonical <- moments_with_jacobian(
      case$J, case$a, case$b, M = 80L
    )

    expect_identical(adapter$mean, canonical$mean)
    expect_identical(adapter$var, canonical$var)
    expect_identical(
      unname(adapter$jacobian), unname(canonical$jacobian)
    )
    expect_identical(
      adapter$derivative_diagnostics, canonical$derivative_diagnostics
    )
    expect_identical(adapter$conditioning, canonical$conditioning)
  }
})

test_that("Jacobian matches numerical derivatives", {
  J <- 50; a <- 2.0; b <- 1.0
  h <- 1e-5

  result <- marginal_moments_with_jacobian(J, a, b)

  # Numerical derivatives via central difference
  mom_a_plus <- exact_K_moments(J, a + h, b)
  mom_a_minus <- exact_K_moments(J, a - h, b)
  mom_b_plus <- exact_K_moments(J, a, b + h)
  mom_b_minus <- exact_K_moments(J, a, b - h)

  num_dmean_da <- (mom_a_plus$mean - mom_a_minus$mean) / (2 * h)
  num_dmean_db <- (mom_b_plus$mean - mom_b_minus$mean) / (2 * h)
  num_dvar_da <- (mom_a_plus$var - mom_a_minus$var) / (2 * h)
  num_dvar_db <- (mom_b_plus$var - mom_b_minus$var) / (2 * h)

  # Compare with relaxed tolerance for numerical differentiation error
  expect_equal(result$jacobian["mean", "a"], num_dmean_da, tolerance = 1e-3)
  expect_equal(result$jacobian["mean", "b"], num_dmean_db, tolerance = 1e-3)
  expect_equal(result$jacobian["var", "a"], num_dvar_da, tolerance = 1e-3)
  expect_equal(result$jacobian["var", "b"], num_dvar_db, tolerance = 1e-3)
})


# =============================================================================
# Test: Input Validation
# =============================================================================

test_that("Invalid inputs throw errors", {
  expect_error(exact_K_moments(-1, 1.5, 0.5))
  expect_error(exact_K_moments(0, 1.5, 0.5))
  expect_error(exact_K_moments(50, -1, 0.5))
  expect_error(exact_K_moments(50, 0, 0.5))
  expect_error(exact_K_moments(50, 1.5, -1))
  expect_error(exact_K_moments(50, 1.5, 0))
})


# =============================================================================
# Test: Edge Cases
# =============================================================================

test_that("Works for small J", {
  result <- exact_K_moments(2, 1.0, 1.0)
  expect_true(result$mean >= 1 && result$mean <= 2)
  expect_true(result$var >= 0)
})

test_that("Works for large J", {
  result <- exact_K_moments(500, 2.0, 1.0)
  expect_true(result$mean >= 1 && result$mean <= 500)
  expect_true(result$var >= 0)
})

test_that("Works for small alpha parameters", {
  result <- exact_K_moments(50, 0.1, 0.1)
  expect_true(is.finite(result$mean))
  expect_true(is.finite(result$var))
})

test_that("Works for large alpha parameters", {
  result <- exact_K_moments(50, 10, 5)
  expect_true(is.finite(result$mean))
  expect_true(is.finite(result$var))
})


# =============================================================================
# Test: NegBin Comparison Structure
# =============================================================================

test_that("compare_to_negbin returns correct structure", {
  result <- compare_to_negbin(50, 1.5, 0.5)

  # The comparison helper returns nested exact and approximation summaries
  expect_true("exact" %in% names(result))
  expect_true("negbin" %in% names(result))
  expect_true("abs_error" %in% names(result))
  expect_true("rel_error" %in% names(result))

  # Check nested structure
  expect_true("mean" %in% names(result$exact))
  expect_true("var" %in% names(result$exact))
  expect_true("mean" %in% names(result$negbin))
  expect_true("var" %in% names(result$negbin))
})

test_that("NegBin approximation overestimates for small J", {
  result <- compare_to_negbin(50, 1.5, 0.5)

  # A1 NegBin typically overestimates both mean and variance
  expect_true(result$rel_error$mean > 0)
  expect_true(result$rel_error$var > 0)
})

test_that("NegBin error decreases with J", {
  result_50 <- compare_to_negbin(50, 1.5, 0.5)
  result_300 <- compare_to_negbin(300, 1.5, 0.5)

  # Relative error should decrease with larger J
  expect_true(abs(result_300$rel_error$mean) < abs(result_50$rel_error$mean))
})


# =============================================================================
# Phase 4 independent marginal-moment contract
# =============================================================================

.adaptive_gamma_expectation <- function(f, a, b) {
  stats::integrate(
    function(alpha) f(alpha) * stats::dgamma(alpha, shape = a, rate = b),
    lower = 0,
    upper = Inf,
    abs.tol = 1e-11,
    rel.tol = 1e-10,
    subdivisions = 1000L,
    stop.on.error = TRUE
  )$value
}


.adaptive_marginal_moments <- function(J, a, b) {
  mean_ref <- .adaptive_gamma_expectation(
    function(alpha) mean_K_given_alpha(J, alpha), a, b
  )
  within_ref <- .adaptive_gamma_expectation(
    function(alpha) var_K_given_alpha(J, alpha), a, b
  )
  between_ref <- .adaptive_gamma_expectation(
    function(alpha) (mean_K_given_alpha(J, alpha) - mean_ref)^2,
    a, b
  )
  c(
    mean = mean_ref,
    var = within_ref + between_ref,
    within = within_ref,
    between = between_ref
  )
}


test_that("Gamma-mixed moments match adaptive integration on a stratified grid", {
  grid <- data.frame(
    J = c(2L, 10L, 50L, 100L, 50L),
    a = c(0.5, 0.5, 1.5, 2.0, 0.5),
    b = c(2.0, 0.5, 0.5, 1.0, 0.2)
  )

  for (i in seq_len(nrow(grid))) {
    case <- grid[i, ]
    reference <- .adaptive_marginal_moments(case$J, case$a, case$b)
    observed <- exact_K_moments(case$J, case$a, case$b, M = 320L)
    mean_budget <- 1e-10 + 2e-9 * abs(reference[["mean"]])
    var_budget <- 1e-10 + 2e-9 * abs(reference[["var"]])

    expect_true(
      abs(observed$mean - reference[["mean"]]) <= mean_budget,
      info = sprintf("adaptive mean J=%d a=%g b=%g", case$J, case$a, case$b)
    )
    expect_true(
      abs(observed$var - reference[["var"]]) <= var_budget,
      info = sprintf("adaptive variance J=%d a=%g b=%g", case$J, case$a, case$b)
    )
  }
})


test_that("marginal moments expose law-of-total-variance and verification metadata", {
  result <- exact_K_moments(
    50, 2, 1, M = 80L, M_verify = 160L,
    abs_tol = 1e-10, rel_tol = 1e-8
  )

  expect_true(all(c("status", "decomposition", "quadrature") %in% names(result)))
  expect_identical(result$status, "converged")
  expect_equal(
    result$decomposition$within_alpha + result$decomposition$between_alpha,
    result$var,
    tolerance = 1e-13
  )
  expect_gte(result$decomposition$within_alpha, 0)
  expect_gte(result$decomposition$between_alpha, 0)
  expect_identical(result$quadrature$M_selected, 80L)
  expect_identical(result$quadrature$M_verification, 160L)
  expect_identical(result$quadrature$M_verification_required, 160L)
  expect_true(result$quadrature$verification_available)
  expect_true(result$quadrature$verification_passed)
})


test_that("unverified and discrepant moment rules are explicitly approximate", {
  unverified <- exact_K_moments(50, 2, 1, M = 80L)
  expect_identical(unverified$status, "approximate")
  expect_identical(unverified$quadrature$reason, "fixed_order_unverified")
  expect_false(unverified$quadrature$verification_performed)
  expect_true(is.na(unverified$quadrature$verification_passed))
  expect_identical(unverified$quadrature$M_verification_required, 160L)
  expect_true(unverified$quadrature$verification_available)

  selected <- exact_K_moments(50, 0.5, 0.2, M = 80L)
  discrepant <- exact_K_moments(
    50, 0.5, 0.2, M = 80L, M_verify = 160L,
    abs_tol = 1e-12, rel_tol = 1e-10
  )
  expect_identical(discrepant$status, "approximate")
  expect_identical(discrepant$quadrature$reason, "higher_order_disagreement")
  expect_false(discrepant$quadrature$verification_passed)
  expect_identical(discrepant$mean, selected$mean)
  expect_identical(discrepant$var, selected$var)
  expect_identical(discrepant$decomposition, selected$decomposition)

  expect_error(
    exact_K_moments(
      50, 0.5, 0.2, M = 80L, M_verify = 160L,
      abs_tol = 1e-12, rel_tol = 1e-10, strict = TRUE
    ),
    class = "dpprior_marginal_convergence_error"
  )
  expect_error(
    exact_K_moments(50, 2, 1, strict = TRUE),
    class = "dpprior_marginal_convergence_error"
  )
})


test_that("J=1 marginal moments retain the deterministic boundary", {
  result <- exact_K_moments(1, 0.2, 4, M = 40L, M_verify = 80L)
  expect_identical(result$mean, 1)
  expect_identical(result$var, 0)
  expect_identical(result$sd, 0)
  expect_identical(result$cv, 0)
  expect_identical(result$decomposition$within_alpha, 0)
  expect_identical(result$decomposition$between_alpha, 0)
  expect_identical(result$status, "converged")
})


test_that("marginal verification controls fail with typed conditions", {
  expect_error(
    exact_K_moments(50, 2, 1, M = 80L, M_verify = 80L),
    class = "dpprior_marginal_verification_error"
  )
  expect_error(
    exact_K_moments(50, 2, 1, M = 80L, M_verify = 80L),
    class = "dpprior_bounds_error"
  )
  expect_error(
    exact_K_moments(50, 2, 1, M = 80L, M_verify = 513L),
    class = "dpprior_marginal_verification_error"
  )
  expect_error(
    exact_K_moments(50, 2, 1, M = 80L, M_verify = 513L),
    class = "dpprior_bounds_error"
  )
  expect_error(
    exact_K_moments(50, 2, 1, abs_tol = -1),
    class = "dpprior_control_error"
  )
})


test_that("marginal verification order honors required and unavailable ranges", {
  lower_boundary <- .marginal_verification_controls(
    1L, M_verify = 41L
  )
  expect_identical(lower_boundary$M_verification_required, 41L)
  expect_true(lower_boundary$verification_available)

  ceiling_boundary <- .marginal_verification_controls(
    256L, M_verify = 512L
  )
  expect_identical(ceiling_boundary$M_verification_required, 512L)
  expect_true(ceiling_boundary$verification_available)

  expect_error(
    .marginal_verification_controls(80L, M_verify = 159L),
    "at least 160", class = "dpprior_marginal_verification_error"
  )
  expect_error(
    .marginal_verification_controls(80L, M_verify = 159L),
    class = "dpprior_bounds_error"
  )

  unavailable <- exact_K_moments(2L, 1, 1, M = 257L)
  expect_identical(unavailable$status, "approximate")
  expect_identical(
    unavailable$quadrature$M_verification_required, 514L
  )
  expect_false(unavailable$quadrature$verification_available)
  expect_false(unavailable$quadrature$verification_performed)
  expect_true(is.na(unavailable$quadrature$M_verification))

  expect_error(
    exact_K_moments(2L, 1, 1, M = 257L, M_verify = 512L),
    "verification is unavailable",
    class = "dpprior_marginal_verification_error"
  )
  expect_error(
    exact_K_moments(2L, 1, 1, M = 257L, strict = TRUE),
    "verification is unavailable",
    class = "dpprior_marginal_convergence_error"
  )
})


test_that("invalid variance decomposition fails with typed numerical error", {
  testthat::local_mocked_bindings(
    var_K_given_alpha = function(J, alpha) rep(-1e-8, length(alpha)),
    .package = "DPprior"
  )

  expect_error(
    exact_K_moments(50, 2, 1),
    class = "dpprior_marginal_moment_error"
  )
  expect_error(
    exact_K_moments(50, 2, 1),
    class = "dpprior_numerical_error"
  )
})
