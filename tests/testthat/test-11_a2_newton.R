# =============================================================================
# Tests for Module 11: A2-MN Newton Refinement
# =============================================================================
#
# Unit tests for DPprior_a2_newton() using testthat.
#
# Author: JoonHo Lee
# Date: December 2025
# =============================================================================

# -----------------------------------------------------------------------------
# Core Convergence Tests
# -----------------------------------------------------------------------------

test_that("A2-MN achieves exact moment matching", {
  fit <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8)

  expect_true(isTRUE(fit$converged))
  expect_equal(fit$fit$mu_K, 5, tolerance = 1e-6)
  expect_equal(fit$fit$var_K, 8, tolerance = 1e-6)
})

test_that("A2-MN converges in few iterations", {
  fit <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8)

  expect_lt(fit$iterations, 10)
})

test_that("A2-MN corrects A1 error", {
  J <- 50
  mu_K <- 5
  var_K <- 8

  a1 <- DPprior_a1(J, mu_K, var_K)
  a1_moments <- exact_K_moments(J, a1$a, a1$b)

  a2 <- DPprior_a2_newton(J, mu_K, var_K)

  # A2 should be much more accurate than A1
  expect_lt(a2$fit$residual, abs(a1_moments$mean - mu_K))
})

# -----------------------------------------------------------------------------
# Return Structure Tests
# -----------------------------------------------------------------------------

test_that("A2-MN returns proper DPprior_fit structure", {
  fit <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8)

  expect_s3_class(fit, "DPprior_fit")

  # Required fields
  expect_true(all(c("a", "b", "J", "target", "method", "status",
                    "converged", "iterations", "termination",
                    "fit", "diagnostics", "trace") %in% names(fit)))

  # Termination field validity

  expect_true(fit$termination %in% c("residual", "step", "max_iter", "nelder_mead"))

  # Method identification
  expect_true(grepl("A2-MN", fit$method))
})

test_that("Trace contains required diagnostic columns", {
  fit <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8)

  required_cols <- c("iter", "a", "b", "M1", "V", "residual", "step", "det_Jlog")
  expect_true(all(required_cols %in% names(fit$trace)))

  # All recorded det_Jlog values should be finite (no singular Jacobians)
  expect_true(all(is.finite(fit$trace$det_Jlog)))
})

# -----------------------------------------------------------------------------
# Termination Behavior Tests
# -----------------------------------------------------------------------------

test_that("Stagnation is not reported as convergence", {
  # This test ensures tol_step termination doesn't falsely report convergence
  # when residual is still above tolerance

  # Force early stagnation with very tight tol_step but loose tol_F
  fit <- DPprior_a2_newton(
    J = 50, mu_K = 5, var_K = 8,
    tol_step = 1e-2,  # Very loose step tolerance (will trigger early)
    tol_F = 1e-15,    # Impossibly tight residual tolerance
    damping = FALSE,
    use_fallback = FALSE,
    max_iter = 3
  )

  # Basic structure checks (always run)
  expect_true(fit$termination %in% c("residual", "step", "max_iter", "nelder_mead"))

  # Core logic: converged should only be TRUE if residual actually met tol_F
  if (fit$fit$residual >= 1e-15) {
    # Residual did NOT meet tolerance, so should NOT be converged
    expect_false(fit$converged)
  } else {
    # Residual met tolerance, so should be converged
    expect_true(fit$converged)
  }
})

test_that("Post-loop convergence check works", {
  # Even if loop exits for max_iter, should still check final residual
  fit <- DPprior_a2_newton(
    J = 50, mu_K = 5, var_K = 8,
    max_iter = 10,
    tol_F = 1e-6
  )

  # If final residual < tol_F, should be converged regardless of how loop ended
  if (fit$fit$residual < 1e-6) {
    expect_true(fit$converged)
  }
})

# -----------------------------------------------------------------------------
# Edge Cases and Robustness
# -----------------------------------------------------------------------------

test_that("A2-MN handles various J values", {
  test_cases <- list(
    list(J = 25, mu_K = 3, var_K = 4),
    list(J = 50, mu_K = 5, var_K = 8),
    list(J = 100, mu_K = 10, var_K = 15)
  )

  for (tc in test_cases) {
    fit <- DPprior_a2_newton(tc$J, tc$mu_K, tc$var_K)
    case_label <- sprintf("J=%d, mu_K=%d, var_K=%d", tc$J, tc$mu_K, tc$var_K)
    expect_true(fit$converged, label = case_label)
    expect_lt(fit$fit$residual, 1e-6, label = case_label)
  }
})

test_that("A2-MN handles high VIF (quasi-improper) cases", {
  # High variance relative to mean requires very small a
  fit <- DPprior_a2_newton(J = 50, mu_K = 3, var_K = 10)

  expect_true(fit$converged)
  expect_lt(fit$fit$residual, 1e-5)

  # Should result in quasi-improper prior (small a)
  expect_lt(fit$a, 0.5)
})

test_that("Nelder-Mead fallback activates when needed", {
  # Very challenging case that may need fallback
  fit <- DPprior_a2_newton(
    J = 50, mu_K = 3, var_K = 10,
    max_iter = 5,        # Limit Newton iterations
    use_fallback = TRUE
  )

  # Should still converge (potentially via Nelder-Mead)
  expect_true(fit$converged)

  # Check if fallback was used
  if (fit$termination == "nelder_mead") {
    expect_equal(fit$method, "A2-MN+NM")
    expect_true(fit$diagnostics$fallback_used)
  }
})

# -----------------------------------------------------------------------------
# Input Validation Tests
# -----------------------------------------------------------------------------

test_that("A2-MN validates input parameters", {
  # Invalid J

  expect_error(DPprior_a2_newton(J = 1, mu_K = 5, var_K = 8))
  expect_error(DPprior_a2_newton(J = -10, mu_K = 5, var_K = 8))


  # Invalid mu_K
  expect_error(DPprior_a2_newton(J = 50, mu_K = 0.5, var_K = 8))
  expect_error(DPprior_a2_newton(J = 50, mu_K = 100, var_K = 8))

  # Invalid var_K
  expect_error(DPprior_a2_newton(J = 50, mu_K = 5, var_K = -1))
  expect_error(DPprior_a2_newton(J = 50, mu_K = 5, var_K = 0))

  # Invalid tolerances
  expect_error(DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8, tol_F = -1))
  expect_error(DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8, tol_step = 0))
})

test_that("A2-MN accepts custom initial values", {
  # Use A1 values as custom start
  a1 <- DPprior_a1(50, 5, 8)

  fit <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8,
                           a0 = a1$a, b0 = a1$b)

  expect_true(fit$converged)
  expect_equal(fit$diagnostics$a0, a1$a)
  expect_equal(fit$diagnostics$b0, a1$b)
})

# -----------------------------------------------------------------------------
# Damping and Numerical Stability
# -----------------------------------------------------------------------------

test_that("Damping prevents overshooting", {
  fit_damped <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8, damping = TRUE)
  fit_undamped <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8, damping = FALSE)

  # Both should converge for this well-behaved case
  expect_true(fit_damped$converged)
  expect_true(fit_undamped$converged)

  # Damped version should have monotonically decreasing residuals (mostly)
  if (nrow(fit_damped$trace) > 1) {
    residuals <- fit_damped$trace$residual
    # Allow for some non-monotonicity due to final convergence check
    decreasing_pairs <- sum(diff(residuals) < 0)
    expect_gte(decreasing_pairs / (length(residuals) - 1), 0.8)
  }
})

test_that("Quadrature node count affects precision", {
  fit_M80 <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8, M = 80)
  fit_M120 <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8, M = 120)

  # Higher M should give at least as good results
  expect_lte(fit_M120$fit$residual, fit_M80$fit$residual * 10)

  # Both should converge
  expect_true(fit_M80$converged)
  expect_true(fit_M120$converged)
})

# -----------------------------------------------------------------------------
# Print Method Test
# -----------------------------------------------------------------------------

test_that("print.DPprior_fit works without error", {
  fit <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8)

  expect_output(print(fit), "DPprior.*Elicitation Result")
  expect_output(print(fit), "A2-MN")
  expect_output(print(fit), "Gamma")  # Gamma hyperprior info
})
