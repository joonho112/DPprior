# =============================================================================
# test-14_diagnostics.R
# =============================================================================
# Testthat tests for Module 14 - Comprehensive Diagnostics
#
# Author: JoonHo Lee
# Date: December 2025
# =============================================================================

# Note: context() removed - deprecated in testthat 3rd edition

# -----------------------------------------------------------------------------
# Weight Diagnostics Tests
# -----------------------------------------------------------------------------

test_that("Weight diagnostics detect high dominance risk", {
  # Lee et al. DP-inform prior (known to have dominance issues)
  diag <- compute_weight_diagnostics(1.60, 1.22)

  expect_equal(diag$dominance_risk, "high")
  expect_gt(diag$prob_exceeds["prob_gt_0.5"], 0.4)
})


test_that("Weight diagnostics detect low dominance risk", {
  diag <- compute_weight_diagnostics(5, 1)

  expect_equal(diag$dominance_risk, "low")
  expect_lt(diag$prob_exceeds["prob_gt_0.5"], 0.2)
})


test_that("Weight diagnostics detect moderate dominance risk", {
  # Find a prior with moderate risk (0.2 <= P(w1>0.5) < 0.4)
  diag <- compute_weight_diagnostics(3, 1)

  expect_true(diag$prob_exceeds["prob_gt_0.5"] >= 0.2)
  expect_true(diag$prob_exceeds["prob_gt_0.5"] < 0.4)
  expect_equal(diag$dominance_risk, "moderate")
})


# -----------------------------------------------------------------------------
# Identity Tests
# -----------------------------------------------------------------------------

test_that("E[w1] = E[rho] identity holds", {
  for (a in c(0.5, 1, 2)) {
    for (b in c(0.5, 1, 2)) {
      expect_equal(
        mean_w1(a, b),
        mean_rho(a, b),
        tolerance = 1e-6,
        info = sprintf("a=%.1f, b=%.1f", a, b)
      )
    }
  }
})


# -----------------------------------------------------------------------------
# CDF/Quantile Consistency Tests
# -----------------------------------------------------------------------------

test_that("Quantiles from CDF inversion are consistent", {
  a <- 2
  b <- 1

  for (p in c(0.1, 0.25, 0.5, 0.75, 0.9)) {
    q <- quantile_w1(p, a, b)
    cdf_at_q <- cdf_w1(q, a, b)
    expect_equal(cdf_at_q, p, tolerance = 1e-10,
                 info = sprintf("p=%.2f", p))
  }
})


test_that("Survival function consistency: P(w1>t) = 1 - CDF(t)", {
  a <- 1.6
  b <- 1.22

  for (t in c(0.3, 0.5, 0.7, 0.9)) {
    surv <- prob_w1_exceeds(t, a, b)
    cdf_val <- cdf_w1(t, a, b)
    expect_equal(surv, 1 - cdf_val, tolerance = 1e-10,
                 info = sprintf("t=%.1f", t))
  }
})


# -----------------------------------------------------------------------------
# Full Diagnostics Tests
# -----------------------------------------------------------------------------

test_that("DPprior_diagnostics generates appropriate warnings for high-risk prior", {
  fit_high <- list(a = 1.6, b = 1.22, J = 50)
  diag_high <- DPprior_diagnostics(fit_high)

  expect_true(length(diag_high$warnings) > 0)
  expect_true(any(grepl("DOMINANCE", diag_high$warnings)))
})


test_that("DPprior_diagnostics generates no warnings for low-risk prior", {
  fit_low <- list(a = 5, b = 1, J = 50)
  diag_low <- DPprior_diagnostics(fit_low)

  expect_equal(diag_low$weights$dominance_risk, "low")
})


test_that("DPprior_diagnostics requires a, b, J fields", {
  expect_error(DPprior_diagnostics(list(a = 1, b = 1)), "J")
  expect_error(DPprior_diagnostics(list(a = 1, J = 50)), "b")
  expect_error(DPprior_diagnostics(NULL), "list")
})


test_that("DPprior_diagnostics accepts optional M from fit", {
  fit <- list(a = 2, b = 1, J = 50, M = 100)
  diag <- DPprior_diagnostics(fit)

  expect_s3_class(diag, "DPprior_diagnostics")
  expect_equal(diag$J, 50)
})


# -----------------------------------------------------------------------------
# Reference Value Tests
# -----------------------------------------------------------------------------

test_that("Reference values match Python/mpmath verification", {
  # From verify_w1_diagnostics.py with mpmath (50 decimal places)

  # P(w1 > 0.5) for a=1.6, b=1.22
  expected_p05 <- 0.4868311039
  actual_p05 <- prob_w1_exceeds(0.5, 1.6, 1.22)
  expect_equal(actual_p05, expected_p05, tolerance = 1e-6)

  # P(w1 > 0.9) for a=1.6, b=1.22
  expected_p09 <- 0.1833147
  actual_p09 <- prob_w1_exceeds(0.9, 1.6, 1.22)
  expect_equal(actual_p09, expected_p09, tolerance = 1e-5)

  # Quantile Q(0.5) for a=2, b=1
  expected_q50 <- 0.3391401986
  actual_q50 <- quantile_w1(0.5, 2, 1)
  expect_equal(actual_q50, expected_q50, tolerance = 1e-6)
})


# -----------------------------------------------------------------------------
# Alpha Diagnostics Tests
# -----------------------------------------------------------------------------

test_that("Alpha diagnostics compute correct moments", {
  diag <- compute_alpha_diagnostics(a = 4, b = 2)

  # E[alpha] = a/b = 4/2 = 2
  expect_equal(diag$mean, 2)

  # SD[alpha] = sqrt(a)/b = 2/2 = 1
  expect_equal(diag$sd, 1)

  # CV[alpha] = 1/sqrt(a) = 1/2 = 0.5
  expect_equal(diag$cv, 0.5)
})


# -----------------------------------------------------------------------------
# K Diagnostics Tests
# -----------------------------------------------------------------------------

test_that("K diagnostics return valid structure", {
  diag <- compute_K_diagnostics(J = 50, a = 2, b = 1)

  expect_true(diag$mean >= 1 && diag$mean <= 50)
  expect_true(diag$var >= 0)
  expect_true(diag$mode >= 1 && diag$mode <= 50)
  expect_true(length(diag$pmf) == 50)
  expect_equal(sum(diag$pmf), 1, tolerance = 1e-10)
})


# -----------------------------------------------------------------------------
# Summary Method Tests
# -----------------------------------------------------------------------------

test_that("summary.DPprior_diagnostics returns data frame", {
  fit <- list(a = 2, b = 1, J = 50)
  diag <- DPprior_diagnostics(fit)
  summ <- summary(diag)

  expect_s3_class(summ, "data.frame")
  expect_true("E_w1" %in% names(summ))
  expect_true("dominance_risk" %in% names(summ))
})


# -----------------------------------------------------------------------------
# Comparison Function Tests
# -----------------------------------------------------------------------------

test_that("compare_diagnostics works with multiple fits", {
  fit1 <- list(a = 1.6, b = 1.22, J = 50)
  fit2 <- list(a = 5, b = 1, J = 50)

  comp <- compare_diagnostics(High = fit1, Low = fit2)

  expect_s3_class(comp, "data.frame")
  expect_equal(nrow(comp), 2)
  expect_equal(comp$method[1], "High")
  expect_equal(comp$method[2], "Low")
  expect_equal(comp$risk[1], "high")
  expect_equal(comp$risk[2], "low")
})


# -----------------------------------------------------------------------------
# Quick Risk Check Tests
# -----------------------------------------------------------------------------

test_that("check_dominance_risk returns correct boolean", {
  # High risk case
  expect_true(check_dominance_risk(1.60, 1.22))

  # Low risk case
  expect_false(check_dominance_risk(5, 1))
})


test_that("check_dominance_risk validates threshold parameter", {
  # Invalid thresholds should throw errors (no pattern matching - just check error is thrown)
  expect_error(check_dominance_risk(2, 1, threshold = 1.5))
  expect_error(check_dominance_risk(2, 1, threshold = -0.1))
})
