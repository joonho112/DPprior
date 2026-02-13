# test-dual-anchor.R
# =============================================================================
# Module 15 (Dual-Anchor) Unit Tests - FINAL VERSION
# =============================================================================
#
# Comprehensive tests for the CORRECTED dual-anchor implementation.
#
# Note: Uses testthat 3rd edition conventions (no context(), no label argument)
#
# Author: JoonHo Lee
# Date: December 2025
# =============================================================================


# =============================================================================
# Helper Functions
# =============================================================================

.create_K_fit <- function(J, mu_K, var_K) {
  fit <- DPprior_a2_newton(J = J, mu_K = mu_K, var_K = var_K, verbose = FALSE)
  if (is.null(fit$J)) fit$J <- J
  if (is.null(fit$target)) {
    fit$target <- list(mu_K = mu_K, var_K = var_K)
  }
  if (!inherits(fit, "DPprior_fit")) {
    class(fit) <- "DPprior_fit"
  }
  fit
}


# =============================================================================
# Test Group 1: Loss Type Comparison (Core Fix Verification)
# =============================================================================

test_that("ABSOLUTE loss barely reduces dominance (confirms original bug)", {
  fit_K <- .create_K_fit(J = 50, mu_K = 5, var_K = 8)

  fit_dual <- DPprior_dual(
    fit_K,
    w1_target = list(prob = list(threshold = 0.5, value = 0.30)),
    lambda = 0.5,
    loss_type = "absolute"
  )

  p_before <- prob_w1_exceeds(0.5, fit_K$a, fit_K$b)
  p_after <- fit_dual$dual_anchor$w1_achieved$prob_gt_50
  reduction <- (p_before - p_after) / p_before

  # ABSOLUTE loss should give < 2% reduction (confirms bug exists)
  expect_lt(reduction, 0.02)
})


test_that("RELATIVE loss reduces dominance significantly", {
  fit_K <- .create_K_fit(J = 50, mu_K = 5, var_K = 8)

  fit_dual <- DPprior_dual(
    fit_K,
    w1_target = list(prob = list(threshold = 0.5, value = 0.30)),
    lambda = 0.5,
    loss_type = "relative"
  )

  p_before <- prob_w1_exceeds(0.5, fit_K$a, fit_K$b)
  p_after <- fit_dual$dual_anchor$w1_achieved$prob_gt_50
  reduction <- (p_before - p_after) / p_before
  gap_before <- abs(p_before - 0.30)
  gap_after <- abs(p_after - 0.30)

  # RELATIVE loss should give > 5% reduction
  expect_gt(reduction, 0.05)
  expect_lt(gap_after, gap_before)
})


test_that("ADAPTIVE loss is more aggressive than RELATIVE", {
  fit_K <- .create_K_fit(J = 50, mu_K = 5, var_K = 8)

  fit_rel <- DPprior_dual(
    fit_K,
    w1_target = list(prob = list(threshold = 0.5, value = 0.30)),
    lambda = 0.5,
    loss_type = "relative"
  )

  fit_adp <- DPprior_dual(
    fit_K,
    w1_target = list(prob = list(threshold = 0.5, value = 0.30)),
    lambda = 0.5,
    loss_type = "adaptive"
  )

  p_rel <- fit_rel$dual_anchor$w1_achieved$prob_gt_50
  p_adp <- fit_adp$dual_anchor$w1_achieved$prob_gt_50

  # ADAPTIVE should give lower (better) P(w1 > 0.5) than RELATIVE
  expect_lt(p_adp, p_rel)
})


test_that("ADAPTIVE scaling factors are computed from both extremes", {
  fit_K <- .create_K_fit(J = 50, mu_K = 5, var_K = 8)

  output <- capture.output({
    fit_adp <- DPprior_dual(
      fit_K,
      w1_target = list(prob = list(threshold = 0.5, value = 0.30)),
      lambda = 0.5,
      loss_type = "adaptive",
      verbose = TRUE
    )
  })

  # Check that scaling info was printed and stored
  scaling_line <- grep("L_K_scale", output, value = TRUE)
  expect_true(length(scaling_line) >= 1)
  expect_true(!is.null(fit_adp$dual_anchor$scaling))
})


# =============================================================================
# Test Group 2: Lambda Boundary Conditions
# =============================================================================

test_that("lambda = 1 recovers K-only solution exactly", {
  fit_K <- .create_K_fit(J = 50, mu_K = 5, var_K = 8)

  fit_dual <- DPprior_dual(
    fit_K,
    w1_target = list(prob = list(threshold = 0.5, value = 0.2)),
    lambda = 1.0
  )

  expect_equal(fit_dual$a, fit_K$a, tolerance = 0.001)
  expect_equal(fit_dual$b, fit_K$b, tolerance = 0.001)
})


test_that("lambda = 0 achieves weight target closely", {
  fit_K <- .create_K_fit(J = 50, mu_K = 5, var_K = 8)
  target_prob <- 0.30

  fit_dual <- DPprior_dual(
    fit_K,
    w1_target = list(prob = list(threshold = 0.5, value = target_prob)),
    lambda = 0.0,
    loss_type = "relative"
  )

  p_achieved <- fit_dual$dual_anchor$w1_achieved$prob_gt_50
  expect_equal(p_achieved, target_prob, tolerance = 0.02)
})


test_that("Intermediate lambda values produce intermediate results", {
  fit_K <- .create_K_fit(J = 50, mu_K = 5, var_K = 8)
  w1_target <- list(prob = list(threshold = 0.5, value = 0.30))

  fit_0 <- DPprior_dual(fit_K, w1_target, lambda = 0.0, loss_type = "relative")
  fit_05 <- DPprior_dual(fit_K, w1_target, lambda = 0.5, loss_type = "relative")
  fit_1 <- DPprior_dual(fit_K, w1_target, lambda = 1.0, loss_type = "relative")

  p_0 <- fit_0$dual_anchor$w1_achieved$prob_gt_50
  p_05 <- fit_05$dual_anchor$w1_achieved$prob_gt_50
  p_1 <- fit_1$dual_anchor$w1_achieved$prob_gt_50

  # P(w1>0.5) should be monotonically increasing with lambda
  expect_lt(p_0, p_05)
  expect_lt(p_05, p_1)
})


# =============================================================================
# Test Group 3: Trade-off Curve
# =============================================================================

test_that("Trade-off curve is computed correctly", {
  curve <- compute_tradeoff_curve(
    J = 50,
    K_target = list(mu_K = 5, var_K = 8),
    w1_target = list(prob = list(threshold = 0.5, value = 0.25)),
    lambda_seq = c(0, 0.5, 1),
    loss_type = "relative"
  )

  expect_equal(nrow(curve), 3)
  expect_true(all(c("lambda", "a", "b", "mu_K", "var_K",
                    "K_loss", "w_loss", "w1_prob_gt_50") %in% names(curve)))
})


test_that("Trade-off curve has expected monotonicity", {
  curve <- compute_tradeoff_curve(
    J = 50,
    K_target = list(mu_K = 5, var_K = 8),
    w1_target = list(prob = list(threshold = 0.5, value = 0.25)),
    lambda_seq = seq(0, 1, by = 0.2),
    loss_type = "relative"
  )

  # P(w1 > 0.5) should increase with lambda, K_loss should decrease
  expect_lt(curve$w1_prob_gt_50[1], curve$w1_prob_gt_50[6])
  expect_gt(curve$K_loss[1], curve$K_loss[6])
})


# =============================================================================
# Test Group 4: Different Target Types
# =============================================================================

test_that("Probability target works correctly", {
  fit_K <- .create_K_fit(J = 50, mu_K = 5, var_K = 8)

  fit_prob <- DPprior_dual(
    fit_K,
    w1_target = list(prob = list(threshold = 0.5, value = 0.30)),
    lambda = 0.5,
    loss_type = "relative"
  )

  expect_true(!is.null(fit_prob$a))
  expect_true(fit_prob$converged || fit_prob$dual_anchor$optim$convergence == 0)
  expect_true(!is.null(fit_prob$dual_anchor$w1_achieved$prob_gt_50))
})


test_that("Quantile target works correctly", {
  fit_K <- .create_K_fit(J = 50, mu_K = 5, var_K = 8)

  fit_q <- DPprior_dual(
    fit_K,
    w1_target = list(quantile = list(prob = 0.9, value = 0.40)),
    lambda = 0.5,
    loss_type = "relative"
  )

  expect_true(!is.null(fit_q$a))

  # Check that 90th percentile moved toward target
  q90_before <- quantile_w1(0.9, fit_K$a, fit_K$b)
  q90_after <- quantile_w1(0.9, fit_q$a, fit_q$b)
  expect_lt(abs(q90_after - 0.40), abs(q90_before - 0.40))
})


test_that("Mean target works correctly", {
  fit_K <- .create_K_fit(J = 50, mu_K = 5, var_K = 8)

  fit_m <- DPprior_dual(
    fit_K,
    w1_target = list(mean = 0.30),
    lambda = 0.5,
    loss_type = "relative"
  )

  expect_true(!is.null(fit_m$a))

  # Check that E[w1] moved toward target
  mean_before <- mean_w1(fit_K$a, fit_K$b)
  mean_after <- mean_w1(fit_m$a, fit_m$b)
  expect_lt(abs(mean_after - 0.30), abs(mean_before - 0.30))
})


# =============================================================================
# Test Group 5: K Moments Quality
# =============================================================================

test_that("K moments remain reasonable at balanced lambda", {
  fit_K <- .create_K_fit(J = 50, mu_K = 5, var_K = 8)

  fit_dual <- DPprior_dual(
    fit_K,
    w1_target = list(prob = list(threshold = 0.5, value = 0.30)),
    lambda = 0.5,
    loss_type = "relative"
  )

  # K moments should not deviate too much at lambda = 0.5
  expect_equal(fit_dual$fit$mu_K, 5, tolerance = 1.5)
  expect_equal(fit_dual$fit$var_K, 8, tolerance = 3.0)
})


test_that("K moments degrade gracefully as lambda decreases", {
  fit_K <- .create_K_fit(J = 50, mu_K = 5, var_K = 8)
  w1_target <- list(prob = list(threshold = 0.5, value = 0.30))

  fit_08 <- DPprior_dual(fit_K, w1_target, lambda = 0.8, loss_type = "relative")
  fit_02 <- DPprior_dual(fit_K, w1_target, lambda = 0.2, loss_type = "relative")

  # K_loss should be smaller at higher lambda
  K_loss_08 <- (fit_08$fit$mu_K - 5)^2 + (fit_08$fit$var_K - 8)^2
  K_loss_02 <- (fit_02$fit$mu_K - 5)^2 + (fit_02$fit$var_K - 8)^2

  expect_lt(K_loss_08, K_loss_02)
})


# =============================================================================
# Test Group 6: Input Validation
# =============================================================================

test_that("Invalid lambda values are rejected", {
  fit_K <- .create_K_fit(J = 50, mu_K = 5, var_K = 8)
  w1_target <- list(prob = list(threshold = 0.5, value = 0.3))

  expect_error(DPprior_dual(fit_K, w1_target, lambda = 1.5))
  expect_error(DPprior_dual(fit_K, w1_target, lambda = -0.1))
  expect_error(DPprior_dual(fit_K, w1_target, lambda = NA))
  expect_error(DPprior_dual(fit_K, w1_target, lambda = c(0.3, 0.5)))
})


test_that("Invalid w1_target specifications are rejected", {
  fit_K <- .create_K_fit(J = 50, mu_K = 5, var_K = 8)

  # Multiple targets
  expect_error(DPprior_dual(fit_K, w1_target = list(
    prob = list(threshold = 0.5, value = 0.3), mean = 0.3)))

  # No target
  expect_error(DPprior_dual(fit_K, w1_target = list()))

  # Malformed targets
  expect_error(DPprior_dual(fit_K, w1_target = list(prob = list(threshold = 0.5))))
  expect_error(DPprior_dual(fit_K, w1_target = list(quantile = list(prob = 0.9))))
})


test_that("Non-DPprior_fit input is rejected", {
  fake_fit <- list(a = 2, b = 1.5, J = 50)
  expect_error(DPprior_dual(fake_fit,
                            w1_target = list(prob = list(threshold = 0.5, value = 0.3))))
})


# =============================================================================
# Test Group 7: Output Structure
# =============================================================================

test_that("Result has correct class and structure", {
  fit_K <- .create_K_fit(J = 50, mu_K = 5, var_K = 8)

  fit_dual <- DPprior_dual(
    fit_K,
    w1_target = list(prob = list(threshold = 0.5, value = 0.30)),
    lambda = 0.5,
    loss_type = "relative"
  )

  expect_s3_class(fit_dual, "DPprior_fit")
  expect_true(all(c("a", "b", "J", "target", "method", "converged",
                    "fit", "dual_anchor") %in% names(fit_dual)))
  expect_equal(fit_dual$method, "dual-anchor")
  expect_true(all(c("w1_target", "lambda", "loss_type", "w1_achieved",
                    "K_loss", "w_loss", "init", "optim") %in%
                    names(fit_dual$dual_anchor)))
  expect_true(all(c("mean", "prob_gt_50", "prob_gt_90") %in%
                    names(fit_dual$dual_anchor$w1_achieved)))
  expect_equal(fit_dual$dual_anchor$loss_type, "relative")
})


test_that("Adaptive result includes scaling information", {
  fit_K <- .create_K_fit(J = 50, mu_K = 5, var_K = 8)

  fit_adp <- DPprior_dual(
    fit_K,
    w1_target = list(prob = list(threshold = 0.5, value = 0.30)),
    lambda = 0.5,
    loss_type = "adaptive"
  )

  expect_true(!is.null(fit_adp$dual_anchor$scaling))
  expect_true(all(c("L_K_scale", "L_w_scale") %in%
                    names(fit_adp$dual_anchor$scaling)))
  expect_gt(fit_adp$dual_anchor$scaling$L_K_scale, 0)
  expect_gt(fit_adp$dual_anchor$scaling$L_w_scale, 0)
})


test_that("Init parameters are stored correctly", {
  fit_K <- .create_K_fit(J = 50, mu_K = 5, var_K = 8)

  fit_dual <- DPprior_dual(
    fit_K,
    w1_target = list(prob = list(threshold = 0.5, value = 0.30)),
    lambda = 0.5,
    loss_type = "relative"
  )

  expect_equal(fit_dual$dual_anchor$init$a, fit_K$a, tolerance = 0.001)
  expect_equal(fit_dual$dual_anchor$init$b, fit_K$b, tolerance = 0.001)
})


# =============================================================================
# Test Group 8: Verification Function
# =============================================================================

test_that("verify_dual_anchor runs successfully", {
  skip_if_not(exists("verify_dual_anchor", mode = "function"),
              "verify_dual_anchor function not available")
  result <- verify_dual_anchor(verbose = FALSE)
  expect_true(result)
})


# =============================================================================
# Test Group 9: Edge Cases and Robustness
# =============================================================================

test_that("Works with different J values", {
  for (J in c(20, 50, 100)) {
    fit_K <- .create_K_fit(J = J, mu_K = 4, var_K = 6)
    fit_dual <- DPprior_dual(
      fit_K,
      w1_target = list(prob = list(threshold = 0.5, value = 0.30)),
      lambda = 0.5,
      loss_type = "relative"
    )
    expect_true(!is.null(fit_dual$a))
  }
})


test_that("Works with extreme lambda values near boundaries", {
  fit_K <- .create_K_fit(J = 50, mu_K = 5, var_K = 8)
  w1_target <- list(prob = list(threshold = 0.5, value = 0.30))

  fit_001 <- DPprior_dual(fit_K, w1_target, lambda = 0.001, loss_type = "relative")
  fit_999 <- DPprior_dual(fit_K, w1_target, lambda = 0.999, loss_type = "relative")

  expect_true(!is.null(fit_001$a))
  expect_true(!is.null(fit_999$a))
})


test_that("Handles challenging targets gracefully", {
  fit_K <- .create_K_fit(J = 50, mu_K = 5, var_K = 8)

  # Very aggressive target
  fit_aggressive <- DPprior_dual(
    fit_K,
    w1_target = list(prob = list(threshold = 0.5, value = 0.10)),
    lambda = 0.3,
    loss_type = "relative"
  )

  expect_true(!is.null(fit_aggressive$a))
  expect_true(is.finite(fit_aggressive$a))
  expect_true(is.finite(fit_aggressive$b))

  # P(w1>0.5) should at least decrease
  p_before <- prob_w1_exceeds(0.5, fit_K$a, fit_K$b)
  p_after <- fit_aggressive$dual_anchor$w1_achieved$prob_gt_50
  expect_lt(p_after, p_before)
})


# =============================================================================
# Test Group 10: Diagnostics Helper
# =============================================================================

test_that("dual_anchor_diagnostics works correctly", {
  skip_if_not(exists("dual_anchor_diagnostics", mode = "function"),
              "dual_anchor_diagnostics function not available")

  fit_K <- .create_K_fit(J = 50, mu_K = 5, var_K = 8)
  fit_dual <- DPprior_dual(
    fit_K,
    w1_target = list(prob = list(threshold = 0.5, value = 0.30)),
    lambda = 0.5,
    loss_type = "relative"
  )

  diag <- dual_anchor_diagnostics(fit_dual, fit_K)

  expect_true(is.data.frame(diag))
  expect_equal(nrow(diag), 7)
  expect_true(all(c("Metric", "K_only", "Dual_anchor") %in% names(diag)))
})
