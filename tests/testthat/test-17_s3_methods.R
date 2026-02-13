# =============================================================================
# testthat: Module 17 - S3 Methods for DPprior_fit Class
# =============================================================================
#
# Comprehensive test suite for print, summary, plot, and as.data.frame S3
# methods defined in R/17_s3_methods.R. Uses both a controlled mock_fit
# fixture (for deterministic output testing) and real fits from DPprior_a1,
# DPprior_a2_newton, and DPprior_dual (for integration testing).
#
# Author: JoonHo Lee (test suite)
# Date: February 2026
# =============================================================================


# =============================================================================
# Fixtures (shared across tests)
# =============================================================================

# --- Mock fit for controlled print/summary testing ---
mock_fit <- list(
  J = 50L,
  a = 1.597,
  b = 1.222,
  target = list(mu_K = 5.0, var_K = 8.0, var_K_used = 8.0,
                confidence = NULL, type = "moments"),
  fit = list(mu_K = 5.000001, var_K = 7.999998, residual = 1.5e-10),
  method = "A2-MN",
  converged = TRUE,
  iterations = 5L,
  status = "success",
  diagnostics = list(
    alpha = list(
      mean = 1.307, sd = 1.034, cv = 0.791,
      quantiles = c(q5 = 0.15, q25 = 0.55, q50 = 1.05,
                    q75 = 1.75, q95 = 3.25)
    ),
    K = list(mean = 5.0, sd = 2.83, mode = 3L),
    weights = list(
      mean = 0.509,
      prob_exceeds = c("prob_gt_0.3" = 0.72, "prob_gt_0.5" = 0.49,
                       "prob_gt_0.7" = 0.25, "prob_gt_0.9" = 0.08),
      dominance_risk = "high"
    )
  ),
  trace = data.frame(
    iter = 1:5,
    a = seq(1.5, 1.597, length.out = 5),
    b = seq(1.1, 1.222, length.out = 5),
    residual = c(0.1, 0.01, 1e-4, 1e-7, 1.5e-10),
    stringsAsFactors = FALSE
  )
)
class(mock_fit) <- "DPprior_fit"

# --- Real fits for integration testing ---
fit_a1 <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
fit_a2 <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8)
fit_dual <- DPprior_dual(
  fit_a2,
  w1_target = list(prob = list(threshold = 0.5, value = 0.3))
)


# =============================================================================
# Section 1: print.DPprior_fit
# =============================================================================

test_that("print.DPprior_fit displays header line", {
  out <- capture.output(print(mock_fit))
  expect_true(any(grepl("DPprior Prior Elicitation Result", out)))
})

test_that("print.DPprior_fit displays Gamma hyperprior parameters", {
  out <- capture.output(print(mock_fit))
  expect_true(any(grepl("Gamma Hyperprior", out)))
  # Should contain the 'a' and 'b' values
  expect_true(any(grepl("1\\.597", out)))
  expect_true(any(grepl("1\\.222", out)))
})

test_that("print.DPprior_fit displays method name", {
  out <- capture.output(print(mock_fit))
  expect_true(any(grepl("Method.*A2-MN", out)))
})

test_that("print.DPprior_fit returns x invisibly", {
  result <- withVisible(print(mock_fit))
  expect_false(result$visible)
  expect_identical(result$value, mock_fit)
})

test_that("print.DPprior_fit respects digits argument", {
  out_2 <- capture.output(print(mock_fit, digits = 2))
  out_6 <- capture.output(print(mock_fit, digits = 6))
  # With digits=2, Gamma params printed to 2 decimal places
  expect_true(any(grepl("1\\.60", out_2)))
  # With digits=6, more precision shown
  expect_true(any(grepl("1\\.597000", out_6)))
})

test_that("print.DPprior_fit shows dominance risk", {
  out <- capture.output(print(mock_fit))
  expect_true(any(grepl("Dominance Risk", out)))
  expect_true(any(grepl("HIGH", out)))
  expect_true(any(grepl("49%", out)))
})

test_that("print.DPprior_fit shows var_K projection when var_K differs from var_K_used", {
  projected_fit <- mock_fit
  projected_fit$target$var_K <- 3.0
  projected_fit$target$var_K_used <- 4.0

  out <- capture.output(print(projected_fit))
  expect_true(any(grepl("used after projection", out)))
  expect_true(any(grepl("3\\.00", out)))
  expect_true(any(grepl("4\\.0000", out)))
})


# =============================================================================
# Section 2: summary.DPprior_fit
# =============================================================================

test_that("summary.DPprior_fit returns summary.DPprior_fit class", {
  summ <- summary(mock_fit, print_output = FALSE)
  expect_s3_class(summ, "summary.DPprior_fit")
})

test_that("summary.DPprior_fit contains all expected components", {
  summ <- summary(mock_fit, print_output = FALSE)

  expect_true("method" %in% names(summ))
  expect_true("status" %in% names(summ))
  expect_true("gamma_prior" %in% names(summ))
  expect_true("alpha_summary" %in% names(summ))
  expect_true("target" %in% names(summ))
  expect_true("achieved" %in% names(summ))
  expect_true("errors" %in% names(summ))
  expect_true("scaling" %in% names(summ))
  expect_true("converged" %in% names(summ))
  expect_true("iterations" %in% names(summ))
  expect_true("diagnostics" %in% names(summ))
})

test_that("summary.DPprior_fit gamma_prior matches original a and b", {
  summ <- summary(mock_fit, print_output = FALSE)

  expect_equal(summ$gamma_prior$a, mock_fit$a)
  expect_equal(summ$gamma_prior$b, mock_fit$b)
})

test_that("summary.DPprior_fit alpha_summary is derived correctly", {
  summ <- summary(mock_fit, print_output = FALSE)

  a <- mock_fit$a
  b <- mock_fit$b
  expect_equal(summ$alpha_summary$E_alpha, a / b)
  expect_equal(summ$alpha_summary$Var_alpha, a / b^2)
  expect_equal(summ$alpha_summary$CV_alpha, 1 / sqrt(a))
  expect_equal(summ$alpha_summary$SD_alpha, sqrt(a / b^2))
})

test_that("summary.DPprior_fit computes errors correctly", {
  summ <- summary(mock_fit, print_output = FALSE)

  expect_equal(summ$errors$mu_K_abs, abs(5.0 - 5.000001), tolerance = 1e-10)
  expect_equal(summ$errors$var_K_abs, abs(8.0 - 7.999998), tolerance = 1e-10)
  expect_true(is.numeric(summ$errors$mu_K_rel_pct))
  expect_true(is.numeric(summ$errors$var_K_rel_pct))
})

test_that("summary.DPprior_fit print_output=FALSE suppresses console output", {
  out <- capture.output(summ <- summary(mock_fit, print_output = FALSE))
  # When print_output=FALSE, no output should be produced
  expect_equal(length(out), 0)
})

test_that("summary.DPprior_fit print_output=TRUE produces console output", {
  out <- capture.output(summ <- summary(mock_fit, print_output = TRUE))
  expect_true(length(out) > 0)
  expect_true(any(grepl("DPprior Prior Elicitation Summary", out)))
})

test_that("summary.DPprior_fit works for all fit types", {
  summ_a1 <- summary(fit_a1, print_output = FALSE)
  expect_s3_class(summ_a1, "summary.DPprior_fit")
  expect_equal(summ_a1$method, "A1")

  summ_a2 <- summary(fit_a2, print_output = FALSE)
  expect_s3_class(summ_a2, "summary.DPprior_fit")
  expect_equal(summ_a2$method, "A2-MN")

  summ_dual <- summary(fit_dual, print_output = FALSE)
  expect_s3_class(summ_dual, "summary.DPprior_fit")
  expect_equal(summ_dual$method, "dual-anchor")
})


# =============================================================================
# Section 3: print.summary.DPprior_fit
# =============================================================================

test_that("print.summary.DPprior_fit displays header", {
  summ <- summary(mock_fit, print_output = FALSE)
  out <- capture.output(print(summ))
  expect_true(any(grepl("DPprior Prior Elicitation Summary", out)))
})

test_that("print.summary.DPprior_fit displays target vs achieved", {
  summ <- summary(mock_fit, print_output = FALSE)
  out <- capture.output(print(summ))
  expect_true(any(grepl("Target vs Achieved", out)))
  expect_true(any(grepl("E\\[K_J\\]", out)))
  expect_true(any(grepl("Var\\(K_J\\)", out)))
})

test_that("print.summary.DPprior_fit diagnostics=TRUE shows full diagnostics", {
  summ <- summary(mock_fit, print_output = FALSE)
  out <- capture.output(print(summ, diagnostics = TRUE))
  expect_true(any(grepl("Full Diagnostics", out)))
})

test_that("print.summary.DPprior_fit diagnostics=FALSE shows reminder", {
  summ <- summary(mock_fit, print_output = FALSE)
  out <- capture.output(print(summ, diagnostics = FALSE))
  expect_true(any(grepl("diagnostics=TRUE", out)))
})

test_that("print.summary.DPprior_fit max_trace limits rows", {
  summ <- summary(mock_fit, print_output = FALSE)
  # mock_fit trace has 5 rows; max_trace=3 should truncate
  out <- capture.output(print(summ, max_trace = 3L))
  expect_true(any(grepl("Iteration Trace", out)))
  expect_true(any(grepl("more rows", out)))
})

test_that("print.summary.DPprior_fit returns x invisibly", {
  summ <- summary(mock_fit, print_output = FALSE)
  result <- withVisible(print(summ))
  expect_false(result$visible)
  expect_s3_class(result$value, "summary.DPprior_fit")
})


# =============================================================================
# Section 4: plot.DPprior_fit
# =============================================================================

test_that("plot.DPprior_fit type='alpha' produces a ggplot", {
  skip_if_not_installed("ggplot2")
  p <- plot(fit_a2, type = "alpha")
  expect_s3_class(p, "ggplot")
})

test_that("plot.DPprior_fit type='K' produces a ggplot", {
  skip_if_not_installed("ggplot2")
  p <- plot(fit_a2, type = "K")
  expect_s3_class(p, "ggplot")
})

test_that("plot.DPprior_fit type='w1' produces a ggplot", {
  skip_if_not_installed("ggplot2")
  p <- plot(fit_a2, type = "w1")
  expect_s3_class(p, "ggplot")
})

test_that("plot.DPprior_fit type='dashboard' produces output", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gridExtra")
  p <- plot(fit_a2, type = "dashboard")
  # dashboard returns a gtable or arrangeGrob, not necessarily ggplot
  expect_true(!is.null(p))
})

test_that("plot.DPprior_fit auto-detects dual for dual-anchor fits", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gridExtra")
  # For dual fits, auto should select "dual" type
  # This should not error
  p <- plot(fit_dual, type = "auto")
  expect_true(!is.null(p))
})

test_that("plot.DPprior_fit engine='base' runs without error", {
  # Base engine should not require ggplot2
  expect_no_error(plot(fit_a2, type = "alpha", engine = "base"))
})


# =============================================================================
# Section 5: as.data.frame.DPprior_fit
# =============================================================================

test_that("as.data.frame.DPprior_fit returns a data.frame", {
  df <- as.data.frame(fit_a2)
  expect_s3_class(df, "data.frame")
})

test_that("as.data.frame.DPprior_fit has expected columns", {
  df <- as.data.frame(fit_a2)
  expected_cols <- c("method", "status", "a", "b", "J",
                     "mu_K", "var_K", "mean_alpha", "cv_alpha",
                     "converged", "iterations")
  for (col in expected_cols) {
    expect_true(col %in% names(df),
                info = paste("missing column:", col))
  }
})

test_that("as.data.frame.DPprior_fit column values match fit", {
  df <- as.data.frame(fit_a2)
  expect_equal(df$a, fit_a2$a)
  expect_equal(df$b, fit_a2$b)
  expect_equal(df$J, fit_a2$J)
  expect_equal(df$mu_K, fit_a2$target$mu_K)
  expect_equal(df$var_K, fit_a2$target$var_K)
  expect_equal(df$mean_alpha, fit_a2$a / fit_a2$b)
  expect_equal(df$cv_alpha, 1 / sqrt(fit_a2$a))
  expect_equal(df$method, fit_a2$method)
})

test_that("as.data.frame.DPprior_fit works for A1 and A2 fit types", {
  df_a1 <- as.data.frame(fit_a1)
  df_a2 <- as.data.frame(fit_a2)

  expect_s3_class(df_a1, "data.frame")
  expect_s3_class(df_a2, "data.frame")

  expect_equal(nrow(df_a1), 1)
  expect_equal(nrow(df_a2), 1)

  # Methods should differ
  expect_equal(df_a1$method, "A1")
  expect_equal(df_a2$method, "A2-MN")
})

test_that("as.data.frame.DPprior_fit errors on dual-anchor (missing status field)", {
  # DPprior_dual does not set $status, causing data.frame() to fail
  # with differing number of rows. This documents a known limitation.
  expect_error(as.data.frame(fit_dual), "differing number of rows")
})


# =============================================================================
# Section 6: .dpprior_is_dual
# =============================================================================

test_that(".dpprior_is_dual returns FALSE for K-only fits", {
  expect_false(.dpprior_is_dual(fit_a1))
  expect_false(.dpprior_is_dual(fit_a2))
  expect_false(.dpprior_is_dual(mock_fit))
})

test_that(".dpprior_is_dual returns TRUE for dual-anchor fits", {
  expect_true(.dpprior_is_dual(fit_dual))
})

test_that(".dpprior_is_dual returns FALSE for partial dual_anchor without init", {
  partial <- mock_fit
  partial$dual_anchor <- list(lambda = 0.5)
  expect_false(.dpprior_is_dual(partial))
})


# =============================================================================
# Section 7: %||% (null-coalescing operator)
# =============================================================================

test_that("%||% returns x when x is not NULL", {
  expect_equal(5 %||% 10, 5)
  expect_equal("hello" %||% "world", "hello")
  expect_equal(FALSE %||% TRUE, FALSE)
})

test_that("%||% returns y when x is NULL", {
  expect_equal(NULL %||% 10, 10)
  expect_equal(NULL %||% "default", "default")
  expect_equal(NULL %||% NA, NA)
})


# =============================================================================
# Section 8: Edge Cases & Integration Tests
# =============================================================================

test_that("print.DPprior_fit works without diagnostics", {
  no_diag_fit <- mock_fit
  no_diag_fit$diagnostics <- NULL
  out <- capture.output(print(no_diag_fit))
  expect_true(any(grepl("DPprior Prior Elicitation Result", out)))
  expect_true(any(grepl("Gamma Hyperprior", out)))
  # Should NOT contain dominance risk line

  expect_false(any(grepl("Dominance Risk", out)))
})

test_that("print.DPprior_fit works without residual in fit", {
  no_resid_fit <- mock_fit
  no_resid_fit$fit$residual <- NULL
  out <- capture.output(print(no_resid_fit))
  expect_true(any(grepl("DPprior Prior Elicitation Result", out)))
  # Should still show achieved section
  expect_true(any(grepl("Achieved", out)))
})

test_that("print.DPprior_fit displays confidence when present", {
  conf_fit <- mock_fit
  conf_fit$target$confidence <- "low"
  out <- capture.output(print(conf_fit))
  expect_true(any(grepl("confidence.*low", out, ignore.case = TRUE)))
})

test_that("summary of dual-anchor fit includes dual_anchor component", {
  summ <- summary(fit_dual, print_output = FALSE)
  expect_true(!is.null(summ$dual_anchor))
  expect_true(!is.null(summ$dual_anchor$lambda))
  expect_true(!is.null(summ$dual_anchor$init))
})

test_that("print.summary.DPprior_fit shows dual-anchor section", {
  summ <- summary(fit_dual, print_output = FALSE)
  out <- capture.output(print(summ))
  expect_true(any(grepl("Dual-Anchor", out)))
  expect_true(any(grepl("Lambda", out, ignore.case = TRUE)))
})

test_that("summary without diagnostics still works", {
  no_diag_fit <- mock_fit
  no_diag_fit$diagnostics <- NULL
  summ <- summary(no_diag_fit, print_output = FALSE)
  expect_s3_class(summ, "summary.DPprior_fit")
  expect_null(summ$diagnostics)
})

test_that("print.DPprior_fit works on real A1 fit", {
  out <- capture.output(print(fit_a1))
  expect_true(length(out) > 0)
  # A1 may use the module-16 or module-17 print depending on load order,
  # but should produce some output about the elicitation
  expect_true(any(grepl("DPprior", out)))
})

test_that("print.DPprior_fit works on real A2 fit", {
  out <- capture.output(print(fit_a2))
  expect_true(length(out) > 0)
  expect_true(any(grepl("DPprior", out)))
})

test_that("print.DPprior_fit works on real dual-anchor fit", {
  out <- capture.output(print(fit_dual))
  expect_true(length(out) > 0)
  expect_true(any(grepl("DPprior", out)))
})

test_that("summary.DPprior_fit preserves trace when present", {
  summ <- summary(mock_fit, print_output = FALSE)
  expect_true(!is.null(summ$trace))
  expect_s3_class(summ$trace, "data.frame")
  expect_equal(nrow(summ$trace), 5)
})

test_that("summary.DPprior_fit handles missing iterations gracefully", {
  no_iter_fit <- mock_fit
  no_iter_fit$iterations <- NA_integer_
  summ <- summary(no_iter_fit, print_output = FALSE)
  expect_true(is.na(summ$iterations))
})

test_that("print.summary.DPprior_fit shows sample size and method", {
  summ <- summary(mock_fit, print_output = FALSE)
  out <- capture.output(print(summ))
  expect_true(any(grepl("J = 50", out)))
  expect_true(any(grepl("A2-MN", out)))
})
