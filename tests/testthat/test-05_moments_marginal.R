# =============================================================================
# Test File: Module 05 - Marginal Moments
# =============================================================================
#
# testthat tests for exact_K_moments() and related functions.
# Tolerances match Python reference precision.
#
# Author: JoonHo Lee
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
# Test: Python Reference Match (Canonical Case)
# =============================================================================

test_that("Python reference values match exactly (J=50, a=1.5, b=0.5)", {
  # Reference computed by dev/verify_moments_marginal.py with M=80
  moments <- exact_K_moments(50, 1.5, 0.5, M = 80)

  # Match to 8 decimal places (Python self-test precision)
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
# Test: Variance Inflation Property
# =============================================================================

test_that("Marginal variance exceeds conditional variance at E[alpha]", {
  test_cases <- list(
    list(J = 50, a = 2.0, b = 1.0),
    list(J = 100, a = 1.5, b = 0.5),
    list(J = 50, a = 0.5, b = 0.5)
  )

  for (tc in test_cases) {
    marginal <- exact_K_moments(tc$J, tc$a, tc$b, M = 80)
    alpha_mean <- tc$a / tc$b
    cond_var <- var_K_given_alpha(tc$J, alpha_mean)

    expect_true(
      marginal$var > cond_var,
      info = sprintf("Variance inflation violated at J=%d: marginal=%.4f, cond=%.4f",
                     tc$J, marginal$var, cond_var)
    )
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


# =============================================================================
# Test: Quadrature Convergence
# =============================================================================

test_that("Moments converge with increasing quadrature nodes", {
  J <- 50; a <- 1.5; b <- 0.5

  M_values <- c(20, 40, 60, 80, 100)
  means <- numeric(length(M_values))
  vars <- numeric(length(M_values))

  for (i in seq_along(M_values)) {
    result <- exact_K_moments(J, a, b, M = M_values[i])
    means[i] <- result$mean
    vars[i] <- result$var
  }

  # Check convergence: differences should decrease
  mean_diffs <- abs(diff(means))
  var_diffs <- abs(diff(vars))

  expect_true(all(mean_diffs < 1e-3))
  expect_true(all(var_diffs < 1e-3))

  # Final values should be very close
  expect_true(abs(means[5] - means[4]) < 1e-8)
  expect_true(abs(vars[5] - vars[4]) < 1e-8)
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
# Test: NegBin Comparison Structure (Claude's nested list structure)
# =============================================================================

test_that("compare_to_negbin returns correct structure", {
  result <- compare_to_negbin(50, 1.5, 0.5)

  # Claude's implementation uses nested lists
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
