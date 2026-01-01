# =============================================================================
# Tests for Module 09: rho (Co-Clustering) Distribution
# =============================================================================
#
# Author: JoonHo Lee
# Date: December 2025
# Part of: DPprior R Package
# =============================================================================

# Note: context() is deprecated in testthat 3rd edition


# =============================================================================
# Test: Conditional Mean E[rho|alpha] = 1/(1+alpha)
# =============================================================================

test_that("E[rho|alpha] = 1/(1+alpha) for various alpha", {
  for (alpha in c(0.5, 1, 2, 5, 10)) {
    expect_equal(mean_rho_given_alpha(alpha), 1 / (1 + alpha))
  }
})

test_that("E[rho|alpha] is vectorized", {
  alphas <- c(0.5, 1, 2, 5, 10)
  result <- mean_rho_given_alpha(alphas)
  expected <- 1 / (1 + alphas)

  expect_length(result, 5)
  expect_equal(result, expected, tolerance = 1e-12)
})


# =============================================================================
# Test: Conditional Variance Formula
# =============================================================================

test_that("Var(rho|alpha) formula is correct for alpha = 2", {
  alpha <- 2
  # Var(rho|alpha) = 2*alpha / [(1+alpha)^2(2+alpha)(3+alpha)]
  # For alpha=2: 2*2 / [(3)^2(4)(5)] = 4 / [9*4*5] = 4/180
  expected <- 2 * 2 / (9 * 4 * 5)
  expect_equal(var_rho_given_alpha(alpha), expected)
})

test_that("Var(rho|alpha) is vectorized and matches formula", {
  alphas <- c(0.5, 1, 2, 5, 10)
  computed <- var_rho_given_alpha(alphas)
  expected <- 2 * alphas / ((1 + alphas)^2 * (2 + alphas) * (3 + alphas))

  expect_length(computed, 5)
  expect_equal(computed, expected, tolerance = 1e-12)
})


# =============================================================================
# Test: Conditional Second Moment
# =============================================================================

test_that("E[rho^2|alpha] formula is correct", {
  alpha <- 2
  # E[rho^2|alpha] = (alpha+6) / [(alpha+1)(alpha+2)(alpha+3)]
  # For alpha=2: (2+6) / [(3)(4)(5)] = 8/60 = 2/15
  expected <- (2 + 6) / (3 * 4 * 5)
  computed <- DPprior:::mean_rho_sq_given_alpha(alpha)

  expect_equal(computed, expected, tolerance = 1e-12)
})

test_that("Variance formula consistent with second moment: Var = E[rho^2] - E[rho]^2", {
  for (alpha in c(0.5, 1, 2, 5, 10)) {
    E_rho <- mean_rho_given_alpha(alpha)
    E_rho_sq <- DPprior:::mean_rho_sq_given_alpha(alpha)
    var_from_moments <- E_rho_sq - E_rho^2
    var_direct <- var_rho_given_alpha(alpha)

    expect_equal(var_from_moments, var_direct, tolerance = 1e-12,
                 info = sprintf("alpha = %.1f", alpha))
  }
})


# =============================================================================
# Test: Monotonicity (GPT addition)
# =============================================================================

test_that("E[rho|alpha] decreases with alpha", {
  alpha_seq <- seq(0.1, 10, by = 0.5)
  means <- mean_rho_given_alpha(alpha_seq)
  expect_true(all(diff(means) < 0))
})


# =============================================================================
# Test: E[w1] = E[rho] Identity
# =============================================================================

test_that("E[w1] = E[rho] identity holds", {
  for (a in c(0.5, 1, 2)) {
    for (b in c(0.5, 1, 2)) {
      expect_true(DPprior:::verify_w1_rho_identity(a, b),
                  info = sprintf("a=%.1f, b=%.1f", a, b))
    }
  }
})


# =============================================================================
# Test: Variance is Positive
# =============================================================================

test_that("Var(rho) is positive for valid parameters", {
  for (a in c(0.5, 1, 2, 5)) {
    for (b in c(0.5, 1, 2)) {
      expect_gt(var_rho(a, b), 0)
    }
  }
})


# =============================================================================
# Test: Law of Total Variance
# =============================================================================

test_that("Law of total variance holds", {
  a <- 2.0
  b <- 1.0
  M <- 100

  var_direct <- var_rho(a, b, M)

  E_var_cond <- integrate_gamma(var_rho_given_alpha, a, b, M)
  E_mean_sq <- integrate_gamma(function(alpha) (1 / (1 + alpha))^2, a, b, M)
  E_mean <- mean_rho(a, b, M)
  var_components <- E_var_cond + (E_mean_sq - E_mean^2)

  expect_equal(var_direct, var_components, tolerance = 1e-10)
})


# =============================================================================
# Test: cv_rho Function
# =============================================================================

test_that("cv_rho equals sd/mean", {
  a <- 2
  b <- 1

  cv_val <- cv_rho(a, b)
  expected <- sqrt(var_rho(a, b)) / mean_rho(a, b)

  expect_equal(cv_val, expected, tolerance = 1e-12)
})


# =============================================================================
# Test: Summary Function
# =============================================================================

test_that("summary_rho returns correct structure", {
  result <- summary_rho(a = 2, b = 1)

  expect_s3_class(result, "rho_summary")
  expect_true(all(c("mean", "var", "sd", "cv", "interpretation",
                    "params", "alpha_prior", "conditional_at_alpha_mean") %in% names(result)))
})

test_that("summary_rho returns consistent results", {
  summary_stats <- summary_rho(2.0, 1.0)
  expect_equal(summary_stats$sd, sqrt(summary_stats$var))
  expect_equal(summary_stats$cv, summary_stats$sd / summary_stats$mean)
})

test_that("summary_rho conditional_at_alpha_mean is correct", {
  a <- 2
  b <- 1
  result <- summary_rho(a, b)

  alpha_mean <- a / b
  expect_equal(result$conditional_at_alpha_mean$alpha, alpha_mean)
  expect_equal(result$conditional_at_alpha_mean$mean, mean_rho_given_alpha(alpha_mean))
  expect_equal(result$conditional_at_alpha_mean$var, var_rho_given_alpha(alpha_mean))
})

test_that("summary_rho interpretation is correct", {
  # High co-clustering (low E[alpha])
  result_high <- summary_rho(a = 0.5, b = 1)
  expect_match(result_high$interpretation, "High co-clustering", ignore.case = TRUE)

  # Low co-clustering (high E[alpha])
  result_low <- summary_rho(a = 10, b = 1)
  expect_match(result_low$interpretation, "Low co-clustering|fragmented", ignore.case = TRUE)
})


# =============================================================================
# Test: Compare rho and w1
# =============================================================================

test_that("compare_rho_w1 confirms mean equality", {
  result <- compare_rho_w1(a = 2, b = 1)

  expect_true(result$mean_equal)
  expect_equal(result$mean_rho, result$mean_w1, tolerance = 1e-10)
})

test_that("compare_rho_w1 shows variance difference", {
  result <- compare_rho_w1(a = 2, b = 1)

  # Variances should be different
  expect_false(isTRUE(all.equal(result$var_rho, result$var_w1)))
  # Typically Var(rho) < Var(w1)
  expect_lt(result$var_ratio, 1)
})


# =============================================================================
# Test: Grid Function
# =============================================================================

test_that("rho_conditional_grid returns correct structure", {
  df <- rho_conditional_grid()

  expect_true(is.data.frame(df))
  expect_true(all(c("alpha", "mean", "var", "sd") %in% names(df)))
  expect_equal(nrow(df), 100)
})

test_that("rho_conditional_grid values are consistent", {
  df <- rho_conditional_grid(alpha_grid = c(1, 2, 5))

  expect_equal(df$mean, 1 / (1 + c(1, 2, 5)), tolerance = 1e-12)
  expect_equal(df$sd, sqrt(df$var), tolerance = 1e-12)
})


# =============================================================================
# Test: Random Generation
# =============================================================================

test_that("rrho generates correct number of samples", {
  set.seed(42)
  samples <- rrho(100, a = 2, b = 1)

  expect_length(samples, 100)
  expect_true(all(samples > 0))
  expect_true(all(samples < 1))
})

test_that("rrho mean approximates theoretical mean", {
  skip_on_cran()  # Skip expensive test on CRAN

  set.seed(123)
  samples <- rrho(5000, a = 2, b = 1)
  theoretical_mean <- mean_rho(a = 2, b = 1)

  # Allow 5% relative error for MC estimate
  expect_equal(mean(samples), theoretical_mean, tolerance = 0.05)
})

test_that("rrho rejects invalid input", {
  expect_error(rrho(-1, a = 2, b = 1), "must be a positive integer")
  expect_error(rrho(100, a = -1, b = 1), "must be finite and positive")
  expect_error(rrho(100, a = 2, b = 0), "must be finite and positive")
})


# =============================================================================
# Test: Input Validation
# =============================================================================

test_that("mean_rho_given_alpha rejects invalid input", {
  expect_error(mean_rho_given_alpha(-1), "must be finite and positive")
  expect_error(mean_rho_given_alpha(0), "must be finite and positive")
  expect_error(mean_rho_given_alpha(NA), "must be finite and positive")
  expect_error(mean_rho_given_alpha(Inf), "must be finite and positive")
})

test_that("var_rho_given_alpha rejects invalid input", {
  expect_error(var_rho_given_alpha(-1), "must be finite and positive")
  expect_error(var_rho_given_alpha(0), "must be finite and positive")
})

test_that("mean_rho rejects invalid a, b", {
  expect_error(mean_rho(0, 1), "must be finite and positive")
  expect_error(mean_rho(1, 0), "must be finite and positive")
  expect_error(mean_rho(-1, 1), "must be finite and positive")
})


# =============================================================================
# Test: Numerical Stability
# =============================================================================

test_that("Functions are stable for extreme parameters", {
  # Very small alpha
  expect_true(is.finite(mean_rho_given_alpha(0.01)))
  expect_true(is.finite(var_rho_given_alpha(0.01)))

  # Very large alpha
  expect_true(is.finite(mean_rho_given_alpha(100)))
  expect_true(is.finite(var_rho_given_alpha(100)))
})

test_that("Functions are stable for challenging (a, b) combinations", {
  # Small a (high variance in alpha)
  result_small_a <- mean_rho(a = 0.5, b = 0.5)
  expect_true(is.finite(result_small_a))
  expect_gt(result_small_a, 0)
  expect_lt(result_small_a, 1)

  # Large b (small mean alpha)
  result_large_b <- mean_rho(a = 1, b = 5)
  expect_true(is.finite(result_large_b))
  expect_gt(result_large_b, 0)
  expect_lt(result_large_b, 1)
})


# =============================================================================
# Test: Golden Reference Values (GPT addition)
# =============================================================================

test_that("Golden values for mean/variance are stable (if available)", {
  golden_path <- system.file("extdata/golden/golden_rho.csv", package = "DPprior")
  if (golden_path == "") {
    skip("golden_rho.csv not installed")
  }

  g <- read.csv(golden_path)
  for (i in seq_len(nrow(g))) {
    a <- g$a[i]
    b <- g$b[i]
    expect_equal(mean_rho(a, b), g$mean_rho[i], tolerance = 5e-4,
                 info = sprintf("mean_rho: a=%.1f, b=%.1f", a, b))
    expect_equal(var_rho(a, b), g$var_rho[i], tolerance = 5e-4,
                 info = sprintf("var_rho: a=%.1f, b=%.1f", a, b))
  }
})
