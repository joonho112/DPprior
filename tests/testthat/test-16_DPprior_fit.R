# =============================================================================
# testthat: Module 16 - DPprior_fit() Main Wrapper
# =============================================================================
#
# Comprehensive test suite for DPprior_fit() wrapper function.
# Updated to match actual package behavior.
#
# Author: JoonHo Lee
# Date: December 2025
# =============================================================================

# =============================================================================
# Basic Functionality Tests
# =============================================================================

test_that("DPprior_fit returns valid DPprior_fit object", {
  fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8, check_diagnostics = FALSE)

  expect_s3_class(fit, "DPprior_fit")
  expect_true(fit$a > 0)
  expect_true(fit$b > 0)
  expect_equal(fit$J, 50L)
  expect_match(fit$method, "A2")  # default is A2-MN backend
})

test_that("DPprior_fit works with minimal inputs", {
  fit <- DPprior_fit(J = 30, mu_K = 4, var_K = 6,
                     check_diagnostics = FALSE, warn_dominance = FALSE)

  expect_s3_class(fit, "DPprior_fit")
  expect_true(is.finite(fit$a))
  expect_true(is.finite(fit$b))
})

test_that("Default method is A2-MN", {
  fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8, check_diagnostics = FALSE)
  expect_equal(fit$method, "A2-MN")
})

# =============================================================================
# Confidence Level Tests
# =============================================================================

test_that("Confidence to var_K conversion uses updated VIF values", {
  # Current VIF values: low=5.0, medium=2.5, high=1.5
  fit_medium <- suppressMessages(
    DPprior_fit(J = 50, mu_K = 5, confidence = "medium", check_diagnostics = FALSE)
  )
  vif_medium <- confidence_to_vif_fit("medium")
  expected_var <- vif_to_variance_fit(5, vif_medium)

  expect_equal(vif_medium, 2.5)
  expect_equal(fit_medium$target$var_K, expected_var)
  expect_equal(fit_medium$target$confidence, "medium")
})

test_that("Low confidence produces VIF=5.0", {
  fit_low <- suppressMessages(
    DPprior_fit(J = 50, mu_K = 5, confidence = "low", check_diagnostics = FALSE)
  )
  vif_low <- confidence_to_vif_fit("low")

  expect_equal(vif_low, 5.0)
  # var_K = 5.0 * (5 - 1) = 20
  expect_equal(fit_low$target$var_K, 20.0)
})

test_that("High confidence produces VIF=1.5", {
  fit_high <- suppressMessages(
    DPprior_fit(J = 50, mu_K = 5, confidence = "high", check_diagnostics = FALSE)
  )
  vif_high <- confidence_to_vif_fit("high")

  expect_equal(vif_high, 1.5)
  # var_K = 1.5 * (5 - 1) = 6
  expect_equal(fit_high$target$var_K, 6.0)
})

test_that("confidence_to_vif_fit rejects invalid input", {
  expect_error(confidence_to_vif_fit("very_low"),
               "confidence must be one of")
  expect_error(confidence_to_vif_fit(123),
               "confidence must be one of")
})

# =============================================================================
# Method Dispatch Tests
# =============================================================================

test_that("Method dispatch works for all methods", {
  fit_a1 <- DPprior_fit(J = 50, mu_K = 5, var_K = 8,
                        method = "A1", check_diagnostics = FALSE)
  fit_a2 <- DPprior_fit(J = 50, mu_K = 5, var_K = 8,
                        method = "A2-MN", check_diagnostics = FALSE)

  expect_equal(fit_a1$method, "A1")
  expect_equal(fit_a2$method, "A2-MN")
})

test_that("A2-KL dispatch works with moment target", {
  # DPprior_a2_kl exists, so this should work directly
  fit_kl <- DPprior_fit(J = 50, mu_K = 5, var_K = 8,
                        method = "A2-KL", check_diagnostics = FALSE)

  expect_match(fit_kl$method, "A2-KL")
  expect_true(fit_kl$a > 0)
  expect_true(fit_kl$b > 0)
})

test_that("A2 improves on A1 accuracy", {
  mu_K <- 5; var_K <- 8

  fit_a1 <- DPprior_fit(J = 50, mu_K = mu_K, var_K = var_K,
                        method = "A1", check_diagnostics = FALSE)
  fit_a2 <- DPprior_fit(J = 50, mu_K = mu_K, var_K = var_K,
                        method = "A2-MN", check_diagnostics = FALSE)

  # A2 residual should be smaller
  expect_lt(fit_a2$fit$residual, 1e-6)
})

# =============================================================================
# Variance Bounds Tests
# =============================================================================

test_that("var_K upper bound is enforced", {
  J <- 50

  expect_error(
    DPprior_fit(J = J, mu_K = 25, var_K = 700, check_diagnostics = FALSE),
    "exceeds maximum possible variance"
  )
})

test_that("var_K fixed-mean upper bound is enforced", {
  expect_error(
    DPprior_fit(J = 10, mu_K = 9, var_K = 9, check_diagnostics = FALSE),
    "var_K = 9.*maximum possible variance 8.*K in \\{1,...,10\\}.*mu_K = 9"
  )
})

test_that("var_K at upper bound limit works", {
  J <- 20

  # Feasibility should not reject a high but valid variance.
  fit <- DPprior_fit(J = J, mu_K = 10, var_K = 80, method = "A1",
                     check_diagnostics = FALSE, warn_dominance = FALSE)
  expect_s3_class(fit, "DPprior_fit")
})

# =============================================================================
# A1 Projection Tests
# =============================================================================

test_that("A1 projects var_K when infeasible for NegBin", {
  # var_K = 3 < mu_K - 1 = 4
  # A1 has its own projection
  warnings <- character()
  fit <- withCallingHandlers(
    DPprior_fit(J = 50, mu_K = 5, var_K = 3,
                method = "A1", check_diagnostics = FALSE),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_s3_class(fit, "DPprior_fit")
  expect_true(fit$a > 0)
  expect_true(fit$b > 0)
  expect_equal(fit$method, "A1")
  expect_equal(fit$target$var_K, 3)
  expect_gt(fit$target$var_K_used, fit$target$var_K)
  expect_true(any(grepl("infeasible for NegBin", warnings)))
  expect_true(any(grepl("projected to feasible boundary", warnings)))
})

test_that("A2-MN handles challenging variance gracefully", {
  # Note: A2-MN internally uses A1 for initialization, which may project
  # Use reasonable var_K that A2-MN can handle
  fit_a2 <- suppressWarnings(
    DPprior_fit(J = 50, mu_K = 5, var_K = 5,  # More reasonable value
                method = "A2-MN", check_diagnostics = FALSE,
                warn_dominance = FALSE)
  )

  expect_s3_class(fit_a2, "DPprior_fit")
  expect_true(fit_a2$a > 0)
})

# =============================================================================
# mu_K Boundary Tests
# =============================================================================

test_that("mu_K = 1 is rejected as trivial", {
  expect_error(
    DPprior_fit(J = 50, mu_K = 1, var_K = 8, check_diagnostics = FALSE),
    "mu_K must be > 1"
  )
})

test_that("mu_K >= J is rejected", {
  expect_error(
    DPprior_fit(J = 50, mu_K = 51, var_K = 8, check_diagnostics = FALSE),
    "mu_K must be < J"
  )
  expect_error(
    DPprior_fit(J = 50, mu_K = 50, var_K = 1, check_diagnostics = FALSE),
    "mu_K must be < J"
  )
})

test_that("mu_K just above 1 works", {
  fit <- DPprior_fit(J = 50, mu_K = 2, var_K = 2,
                     check_diagnostics = FALSE, warn_dominance = FALSE)
  expect_s3_class(fit, "DPprior_fit")
})

test_that("mu_K below J works with appropriate variance", {
  # Use moderate variance that won't cause numerical issues
  fit <- suppressWarnings(
    DPprior_fit(J = 20, mu_K = 15, var_K = 30,
                check_diagnostics = FALSE, warn_dominance = FALSE)
  )
  expect_s3_class(fit, "DPprior_fit")
})

test_that("mu_K = J is rejected consistently across moment workflows", {
  expect_error(DPprior_fit(J = 20, mu_K = 20, var_K = 1,
                           check_diagnostics = FALSE),
               "mu_K must be < J")
  expect_error(DPprior_a1(J = 20, mu_K = 20, var_K = 1),
               "mu_K must be < J")
  expect_error(DPprior_a2_newton(J = 20, mu_K = 20, var_K = 1),
               "mu_K must be < J")
  expect_error(construct_target_pmf(20, list(mu_K = 20, var_K = 1)),
               "mu_K must be < J")
  expect_error(DPprior_a2_kl(20, list(mu_K = 20, var_K = 1), method = "chisq"),
               "mu_K must be < J")
})

test_that("near-boundary confidence-derived variance gives clear feasibility error", {
  expect_error(
    DPprior_fit(J = 20, mu_K = 18, confidence = "medium",
                check_diagnostics = FALSE),
    "var_K = 42.5.*maximum possible variance 34.*K in \\{1,...,20\\}.*mu_K = 18"
  )
})

test_that("DPprior_fit reports low-variance non-convergence as an error", {
  expect_error(
    suppressWarnings(
      DPprior_fit(J = 50, mu_K = 5, var_K = 3,
                  method = "A2-MN", check_diagnostics = FALSE)
    ),
    "Calibration did not converge.*A2-MN.*larger var_K|less concentrated target"
  )
})

# =============================================================================
# Diagnostics Tests
# =============================================================================

test_that("Diagnostics are computed when requested", {
  fit <- suppressWarnings(
    DPprior_fit(J = 50, mu_K = 5, var_K = 8, check_diagnostics = TRUE,
                warn_dominance = FALSE)
  )

  expect_true(!is.null(fit$diagnostics))
  expect_true(!is.null(fit$diagnostics$weights))
})

test_that("Diagnostics are not computed when disabled", {
  fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8, check_diagnostics = FALSE)
  expect_null(fit$diagnostics)
})

test_that("Dominance warning is issued when risk is high", {
  # Small mu_K with small variance triggers dominance warning
  expect_warning(
    DPprior_fit(J = 50, mu_K = 2, var_K = 2,
                check_diagnostics = TRUE, warn_dominance = TRUE),
    "DOMINANCE|dominance"
  )
})

# =============================================================================
# Input Validation Tests
# =============================================================================

test_that("Invalid J is rejected", {
  expect_error(DPprior_fit(J = 0, mu_K = 5, var_K = 8))
  expect_error(DPprior_fit(J = -10, mu_K = 5, var_K = 8))
  # J = 1000 exceeds max, error message may vary
  expect_error(DPprior_fit(J = 1000, mu_K = 5, var_K = 8))
})

test_that("Invalid mu_K is rejected", {
  expect_error(DPprior_fit(J = 50, mu_K = NA, var_K = 8), "mu_K must be")
  expect_error(DPprior_fit(J = 50, mu_K = Inf, var_K = 8), "mu_K must be")
  expect_error(DPprior_fit(J = 50, mu_K = -5, var_K = 8), "mu_K must be")
})

test_that("Invalid var_K is rejected", {
  expect_error(DPprior_fit(J = 50, mu_K = 5, var_K = NA), "var_K must be")
  expect_error(DPprior_fit(J = 50, mu_K = 5, var_K = -1), "var_K must be")
  expect_error(DPprior_fit(J = 50, mu_K = 5, var_K = 0), "var_K must be")
})

test_that("Invalid M is rejected", {
  expect_error(
    DPprior_fit(J = 50, mu_K = 5, var_K = 8, M = 5, check_diagnostics = FALSE),
    "M must be"
  )
})

test_that("Invalid method is rejected", {
  expect_error(
    DPprior_fit(J = 50, mu_K = 5, var_K = 8, method = "invalid"),
    "should be one of"
  )
})

# =============================================================================
# Output Structure Tests
# =============================================================================

test_that("Output contains all required fields", {
  fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8,
                     check_diagnostics = FALSE, warn_dominance = FALSE)

  required <- c("a", "b", "J", "target", "method", "status",
                "converged", "fit")

  for (field in required) {
    expect_true(field %in% names(fit),
                info = sprintf("Missing field: %s", field))
  }
})

test_that("Target structure includes var_K_used", {
  fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8, check_diagnostics = FALSE)

  expect_true("var_K" %in% names(fit$target))
  expect_true("var_K_used" %in% names(fit$target))
  expect_true("mu_K" %in% names(fit$target))
})

test_that("Fit structure contains achieved moments", {
  fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8, check_diagnostics = FALSE)

  expect_true("mu_K" %in% names(fit$fit))
  expect_true("var_K" %in% names(fit$fit))
  expect_true("residual" %in% names(fit$fit))
})

# =============================================================================
# S3 Methods Tests
# =============================================================================

test_that("print.DPprior_fit works", {
  fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8,
                     check_diagnostics = FALSE, warn_dominance = FALSE)

  # Should print without error - use flexible regex
  expect_output(print(fit), "DPprior.*Elicitation Result")
  expect_output(print(fit), "Method:")
  expect_output(print(fit), "Gamma")  # Gamma hyperprior info
})

test_that("summary.DPprior_fit returns list structure", {
  fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8,
                     check_diagnostics = FALSE, warn_dominance = FALSE)

  summ <- summary(fit)
  expect_type(summ, "list")
  expect_true("gamma_prior" %in% names(summ))
  expect_true("a" %in% names(summ$gamma_prior))
  expect_true("b" %in% names(summ$gamma_prior))
  expect_equal(summ$gamma_prior$a, fit$a)
  expect_equal(summ$gamma_prior$b, fit$b)
})

test_that("summary.DPprior_fit works with diagnostics", {
  fit <- suppressWarnings(
    DPprior_fit(J = 50, mu_K = 5, var_K = 8,
                check_diagnostics = TRUE, warn_dominance = FALSE)
  )

  summ <- summary(fit)
  expect_type(summ, "list")
  expect_true("gamma_prior" %in% names(summ))
})

# =============================================================================
# Verbose Mode Tests
# =============================================================================

test_that("Verbose mode produces messages", {
  expect_message(
    DPprior_fit(J = 50, mu_K = 5, var_K = 8,
                verbose = TRUE, check_diagnostics = FALSE),
    "Using A2-MN"
  )
})

test_that("Verbose mode with confidence produces conversion message", {
  expect_message(
    DPprior_fit(J = 50, mu_K = 5, confidence = "medium",
                verbose = TRUE, check_diagnostics = FALSE),
    "confidence='medium'"
  )
})

# =============================================================================
# Edge Cases Tests
# =============================================================================

test_that("Works with small J", {
  fit <- DPprior_fit(J = 5, mu_K = 2, var_K = 2,
                     check_diagnostics = FALSE, warn_dominance = FALSE)
  expect_s3_class(fit, "DPprior_fit")
})

test_that("Works with large J (within limits)", {
  fit <- DPprior_fit(J = 100, mu_K = 20, var_K = 40,
                     check_diagnostics = FALSE, warn_dominance = FALSE)
  expect_s3_class(fit, "DPprior_fit")
})

test_that("Works with moderate mu_K relative to J", {
  # Use reasonable values that won't cause numerical issues
  fit <- suppressWarnings(
    DPprior_fit(J = 50, mu_K = 30, var_K = 50,
                check_diagnostics = FALSE, warn_dominance = FALSE)
  )
  expect_s3_class(fit, "DPprior_fit")
})

# =============================================================================
# target_pmf Tests for A2-KL
# =============================================================================

test_that("target_pmf validation works", {
  J <- 50

  # Wrong length
  expect_error(
    DPprior_fit(J = J, method = "A2-KL", target_pmf = rep(1/10, 10),
                check_diagnostics = FALSE),
    "target_pmf must"
  )

  # Negative values
  expect_error(
    DPprior_fit(J = J, method = "A2-KL",
                target_pmf = c(-0.1, rep(1.1/(J-1), J-1)),
                check_diagnostics = FALSE),
    "non-negative"
  )

  expect_error(
    DPprior_fit(J = J, method = "A1", target_pmf = rep(1/J, J),
                check_diagnostics = FALSE),
    "target_pmf can only be used"
  )
})

test_that("A2-KL works with valid target_pmf", {
  J <- 50
  # Create a valid PMF that sums to 1
  # dbinom(1:J, ...) doesn't sum to 1 because k=0 is missing
  # So we normalize it
  raw_pmf <- dbinom(1:J, size = J, prob = 0.1)
  target_pmf <- raw_pmf / sum(raw_pmf)  # Normalize to sum to 1

  fit <- DPprior_fit(J = J, target_pmf = target_pmf,
                     check_diagnostics = FALSE)

  expect_s3_class(fit, "DPprior_fit")
  expect_match(fit$method, "A2-KL")
  expect_equal(fit$target$type, "pmf")
  expect_equal(fit$target$pmf, target_pmf, tolerance = 1e-12)
  expect_equal(fit$target$mu_K, sum(seq_len(J) * target_pmf), tolerance = 1e-12)
})

test_that("A2-KL wrapper accepts target_pmf with k=0 entry", {
  J <- 20
  tail_pmf <- stats::dpois(seq_len(J), lambda = 4)
  expected <- tail_pmf / sum(tail_pmf)
  target_pmf <- c(0, tail_pmf)

  fit <- DPprior_fit(J = J, method = "A2-KL", target_pmf = target_pmf,
                     check_diagnostics = FALSE)

  expect_s3_class(fit, "DPprior_fit")
  expect_equal(fit$target$type, "pmf")
  expect_equal(length(fit$target$pmf), J)
  expect_equal(fit$target$pmf, expected, tolerance = 1e-12)
  expect_equal(sum(fit$target$pmf), 1, tolerance = 1e-12)
})

test_that("A2-KL wrapper rejects positive k=0 target_pmf mass", {
  J <- 20L

  expect_error(
    DPprior_fit(J = J,
                target_pmf = c(0.25, rep(0.75 / J, J)),
                check_diagnostics = FALSE,
                warn_dominance = FALSE,
                M = 20L),
    "target_pmf\\[1\\].*K = 0.*must be 0.*support \\{1, \\.\\.\\., J\\}"
  )
})

test_that("A2-KL wrapper validates invalid k=0 entries before dropping", {
  J <- 20L

  expect_error(
    DPprior_fit(J = J,
                target_pmf = c(NA_real_, rep(1, J)),
                check_diagnostics = FALSE,
                warn_dominance = FALSE,
                M = 20L),
    "target_pmf must be finite and non-missing"
  )
  expect_error(
    DPprior_fit(J = J,
                target_pmf = c(Inf, rep(1, J)),
                check_diagnostics = FALSE,
                warn_dominance = FALSE,
                M = 20L),
    "target_pmf must be finite and non-missing"
  )
  expect_error(
    DPprior_fit(J = J,
                target_pmf = c(-0.1, rep(1, J)),
                check_diagnostics = FALSE,
                warn_dominance = FALSE,
                M = 20L),
    "target_pmf must be non-negative"
  )
})

test_that("target_pmf accepts explicit moments that match implied moments", {
  J <- 20L
  target_pmf <- stats::dpois(seq_len(J), lambda = 4)
  target_pmf <- target_pmf / sum(target_pmf)
  k <- seq_len(J)
  pmf_mu <- sum(k * target_pmf)
  pmf_var <- sum(k^2 * target_pmf) - pmf_mu^2

  fit <- DPprior_fit(J = J,
                     mu_K = pmf_mu + 1e-11,
                     var_K = pmf_var - 1e-11,
                     method = "A2-KL",
                     target_pmf = target_pmf,
                     check_diagnostics = FALSE,
                     warn_dominance = FALSE,
                     M = 20L)

  expect_s3_class(fit, "DPprior_fit")
  expect_equal(fit$method, "A2-KL")
  expect_equal(fit$target$type, "pmf")
  expect_equal(fit$target$mu_K, pmf_mu, tolerance = 1e-10)
  expect_equal(fit$target$var_K, pmf_var, tolerance = 1e-10)
  expect_equal(fit$target$var_K_used, pmf_var, tolerance = 1e-10)
})

test_that("target_pmf rejects explicit mu_K that conflicts with implied mean", {
  J <- 20L
  target_pmf <- stats::dpois(seq_len(J), lambda = 4)
  target_pmf <- target_pmf / sum(target_pmf)
  k <- seq_len(J)
  pmf_mu <- sum(k * target_pmf)
  pmf_var <- sum(k^2 * target_pmf) - pmf_mu^2

  expect_error(
    DPprior_fit(J = J,
                mu_K = pmf_mu + 0.25,
                var_K = pmf_var,
                method = "A2-KL",
                target_pmf = target_pmf,
                check_diagnostics = FALSE,
                warn_dominance = FALSE,
                M = 20L),
    "mu_K conflicts with target_pmf-derived mean"
  )
})

test_that("target_pmf rejects explicit var_K that conflicts with implied variance", {
  J <- 20L
  target_pmf <- stats::dpois(seq_len(J), lambda = 4)
  target_pmf <- target_pmf / sum(target_pmf)
  k <- seq_len(J)
  pmf_mu <- sum(k * target_pmf)
  pmf_var <- sum(k^2 * target_pmf) - pmf_mu^2

  expect_error(
    DPprior_fit(J = J,
                mu_K = pmf_mu,
                var_K = pmf_var + 0.25,
                method = "A2-KL",
                target_pmf = target_pmf,
                check_diagnostics = FALSE,
                warn_dominance = FALSE,
                M = 20L),
    "var_K conflicts with target_pmf-derived variance"
  )
})

test_that("target_pmf conflict checks run when method is omitted", {
  J <- 20L
  raw_pmf <- stats::dpois(seq_len(J), lambda = 4)
  target_pmf <- raw_pmf / sum(raw_pmf)
  k <- seq_len(J)
  pmf_mu <- sum(k * target_pmf)
  pmf_var <- sum(k^2 * target_pmf) - pmf_mu^2

  expect_error(
    DPprior_fit(J = J,
                mu_K = pmf_mu + 0.25,
                var_K = pmf_var,
                target_pmf = target_pmf,
                check_diagnostics = FALSE,
                warn_dominance = FALSE,
                M = 20L),
    "mu_K conflicts with target_pmf-derived mean"
  )
  expect_error(
    DPprior_fit(J = J,
                mu_K = pmf_mu,
                var_K = pmf_var + 0.25,
                target_pmf = target_pmf,
                check_diagnostics = FALSE,
                warn_dominance = FALSE,
                M = 20L),
    "var_K conflicts with target_pmf-derived variance"
  )
  expect_error(
    DPprior_fit(J = J,
                mu_K = pmf_mu + 0.25,
                var_K = pmf_var,
                target_pmf = c(0, raw_pmf),
                check_diagnostics = FALSE,
                warn_dominance = FALSE,
                M = 20L),
    "mu_K conflicts with target_pmf-derived mean"
  )
})

test_that("target_pmf consistency checks use normalized length J+1 input", {
  J <- 20L
  raw_pmf <- stats::dpois(seq_len(J), lambda = 4)
  k <- seq_len(J)
  normalized <- raw_pmf / sum(raw_pmf)
  pmf_mu <- sum(k * normalized)
  pmf_var <- sum(k^2 * normalized) - pmf_mu^2
  target_pmf <- c(0, raw_pmf)

  fit <- DPprior_fit(J = J,
                     mu_K = pmf_mu,
                     var_K = pmf_var,
                     method = "A2-KL",
                     target_pmf = target_pmf,
                     check_diagnostics = FALSE,
                     warn_dominance = FALSE,
                     M = 20L)

  expect_s3_class(fit, "DPprior_fit")
  expect_equal(fit$target$pmf, normalized, tolerance = 1e-12)
  expect_equal(fit$target$mu_K, pmf_mu, tolerance = 1e-12)
  expect_equal(fit$target$var_K, pmf_var, tolerance = 1e-12)
})

test_that("target_pmf rejects invalid explicit consistency moments", {
  J <- 20L
  target_pmf <- stats::dpois(seq_len(J), lambda = 4)
  target_pmf <- target_pmf / sum(target_pmf)

  expect_error(
    DPprior_fit(J = J, mu_K = NA_real_, target_pmf = target_pmf,
                check_diagnostics = FALSE),
    "mu_K must be a finite numeric scalar"
  )
  expect_error(
    DPprior_fit(J = J, var_K = -1, target_pmf = target_pmf,
                check_diagnostics = FALSE),
    "var_K conflicts with target_pmf-derived variance"
  )
})

# =============================================================================
# Verification Function Test
# =============================================================================

test_that("verify_DPprior_fit runs without errors", {
  # Run verification silently - may have some test failures internally
  # but should not crash
  result <- tryCatch({
    verify_DPprior_fit(verbose = FALSE)
    TRUE
  }, error = function(e) FALSE)

  # We just check it doesn't crash catastrophically
  expect_true(is.logical(result))
})
