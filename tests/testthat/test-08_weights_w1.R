# =============================================================================
# Unit Tests: Module 08 - First Stick-Breaking Weight (w₁) Distribution
# =============================================================================
#
# This file contains comprehensive tests for the w₁ distribution functions.
# Golden values are computed directly from the closed-form formulas and
# verified by the inverse identity CDF(Q(u)) = u.
#
# Author: JoonHo Lee (jlee296@ua.edu)
# Date: December 2025
# Part of: DPprior R Package
# =============================================================================

library(testthat)

# =============================================================================
# Test 1: CDF and Quantile are Inverses
# =============================================================================

test_that("CDF and quantile are inverses", {
  for (a in c(0.5, 1, 2, 5)) {
    for (b in c(0.5, 1, 2)) {
      for (u in c(0.1, 0.5, 0.9)) {
        q <- quantile_w1(u, a, b)

        # For some (a, b, u) combinations, the true quantile is so close to 1
        # that it rounds to 1.0 in double precision. Skip strict inversion
        # checks in that numerically-degenerate regime.
        if (q >= 1 - 1e-15) next

        p <- cdf_w1(q, a, b)
        expect_equal(p, u, tolerance = 1e-10,
                     info = sprintf("a=%.1f, b=%.1f, u=%.1f", a, b, u))
      }
    }
  }
})


# =============================================================================
# Test 2: CDF is Monotonically Increasing
# =============================================================================

test_that("CDF is monotonically increasing", {
  a <- 2; b <- 1
  x_seq <- seq(0.01, 0.99, by = 0.01)
  cdf_vals <- cdf_w1(x_seq, a, b)
  expect_true(all(diff(cdf_vals) > 0))
})


# =============================================================================
# Test 3: Density Integrates to 1
# =============================================================================

test_that("Density integrates to 1", {
  a <- 2; b <- 1
  # Use upper limit very close to 1 for better coverage of the tail
  integral <- integrate(function(x) density_w1(x, a, b), 1e-6, 1 - 1e-16)$value
  expect_equal(integral, 1, tolerance = 1e-3)
})


# =============================================================================
# Test 4: Survival Function Equals 1 - CDF
# =============================================================================

test_that("prob_w1_exceeds matches 1 - CDF", {
  a <- 1.6; b <- 1.22
  t <- 0.5
  expect_equal(prob_w1_exceeds(t, a, b), 1 - cdf_w1(t, a, b), tolerance = 1e-12)
})


# =============================================================================
# Test 5: Conditional Expectation E[w₁|α] = 1/(1+α)
# =============================================================================

test_that("Conditional expectation E[w1|alpha] = 1/(1+alpha)", {
  alpha <- 2
  cond_mean <- integrate(function(x) x * alpha * (1 - x)^(alpha - 1), 0, 1)$value
  expect_equal(cond_mean, 1 / (1 + alpha), tolerance = 1e-12)
})


# =============================================================================
# Test 6: Golden Values (High-Precision, Verified by Inverse Identity)
# =============================================================================

test_that("w1 quantiles match golden values (closed form)", {

  # These values are computed from the closed-form quantile function and
  # verified by the inverse identity F(Q(u)) = u.
  # Note: Using tolerance = 1e-9 to account for floating point variations
  # across different platforms and R versions.
  golden <- data.frame(
    a = c(1.6, 2.0, 0.5),
    b = c(1.22, 1.0, 1.0),
    q50 = c(0.4839219199376489, 0.3391401985931721, 0.9502129316321360),
    prob_gt_50 = c(0.4868311039310716, 0.3488273884247458, 0.7685155230091048)
  )

  for (i in seq_len(nrow(golden))) {
    actual_q50 <- quantile_w1(0.5, golden$a[i], golden$b[i])
    actual_prob <- prob_w1_exceeds(0.5, golden$a[i], golden$b[i])

    expect_equal(actual_q50, golden$q50[i], tolerance = 1e-9,
                 info = sprintf("q50 for a=%.1f, b=%.2f", golden$a[i], golden$b[i]))
    expect_equal(actual_prob, golden$prob_gt_50[i], tolerance = 1e-9,
                 info = sprintf("P(w1>0.5) for a=%.1f, b=%.2f",
                                golden$a[i], golden$b[i]))
  }
})


# =============================================================================
# Test 7: Boundary Handling
# =============================================================================

test_that("CDF handles boundaries correctly", {
  a <- 2; b <- 1

  # At or below 0
  expect_equal(cdf_w1(0, a, b), 0)
  expect_equal(cdf_w1(-1, a, b), 0)
  expect_equal(cdf_w1(-Inf, a, b), 0)

  # At or above 1
  expect_equal(cdf_w1(1, a, b), 1)

  expect_equal(cdf_w1(2, a, b), 1)
  expect_equal(cdf_w1(Inf, a, b), 1)

  # NA handling
  expect_true(is.na(cdf_w1(NA, a, b)))
  expect_true(is.na(cdf_w1(NaN, a, b)))
})


test_that("prob_w1_exceeds handles boundaries correctly", {
  a <- 2; b <- 1

  # At or below 0
  expect_equal(prob_w1_exceeds(0, a, b), 1)
  expect_equal(prob_w1_exceeds(-1, a, b), 1)

  # At or above 1
  expect_equal(prob_w1_exceeds(1, a, b), 0)
  expect_equal(prob_w1_exceeds(2, a, b), 0)
})


test_that("quantile_w1 handles boundaries correctly", {
  a <- 2; b <- 1

  expect_equal(quantile_w1(0, a, b), 0)
  expect_equal(quantile_w1(1, a, b), 1)
})


test_that("quantile_w1 keeps a representable tiny-rate log product", {
  u <- 0.5
  a <- -log1p(-u) / 710
  b <- 1e-308
  reference <- 0.89290026396189548

  value <- quantile_w1(u, a, b)
  expect_true(is.finite(value))
  expect_lt(value, 1)
  expect_equal(value, reference, tolerance = 1e-13)
  expect_equal(cdf_w1(value, a, b), u, tolerance = 2e-14)
})


test_that("density_w1 handles boundaries correctly", {
  a <- 2; b <- 1

  # Outside (0, 1) -> 0 (or -Inf on log scale)
  expect_equal(density_w1(0, a, b), 0)
  expect_equal(density_w1(1, a, b), 0)
  expect_equal(density_w1(-1, a, b), 0)
  expect_equal(density_w1(2, a, b), 0)

  expect_equal(density_w1(0, a, b, log = TRUE), -Inf)
  expect_equal(density_w1(1, a, b, log = TRUE), -Inf)
})


# =============================================================================
# Test 8: Input Validation
# =============================================================================

test_that("input validation works correctly", {
  # Negative parameters
  expect_error(cdf_w1(0.5, a = -1, b = 1))
  expect_error(cdf_w1(0.5, a = 1, b = -1))
  expect_error(cdf_w1(0.5, a = 0, b = 1))
  expect_error(cdf_w1(0.5, a = 1, b = 0))

  # Invalid probability
  expect_error(quantile_w1(-0.1, a = 1, b = 1))
  expect_error(quantile_w1(1.1, a = 1, b = 1))
})


# =============================================================================
# Test 9: Vectorization
# =============================================================================

test_that("functions are properly vectorized", {
  a <- 2; b <- 1

  # CDF
  x_vec <- c(0.1, 0.3, 0.5, 0.7)
  cdf_vals <- cdf_w1(x_vec, a, b)
  expect_length(cdf_vals, 4)
  expect_true(all(cdf_vals >= 0 & cdf_vals <= 1))

  # Quantile
  u_vec <- c(0.1, 0.5, 0.9)
  q_vals <- quantile_w1(u_vec, a, b)
  expect_length(q_vals, 3)
  expect_true(all(q_vals >= 0 & q_vals <= 1))

  # Density
  dens_vals <- density_w1(x_vec, a, b)
  expect_length(dens_vals, 4)
  expect_true(all(dens_vals > 0))
})


# =============================================================================
# Test 10: Mean and Variance via Monte Carlo
# =============================================================================

test_that("mean_w1 matches Monte Carlo estimate", {
  skip_on_cran()  # Skip on CRAN due to Monte Carlo variance

  set.seed(42)
  a <- 2; b <- 1
  n_samples <- 50000

  samples <- rw1(n_samples, a, b)
  mc_mean <- mean(samples)
  quad_mean <- mean_w1(a, b)

  # Allow 5% relative error for MC
  expect_equal(quad_mean, mc_mean, tolerance = 0.05)
})


test_that("var_w1 matches Monte Carlo estimate", {
  skip_on_cran()

  set.seed(42)
  a <- 2; b <- 1
  n_samples <- 50000

  samples <- rw1(n_samples, a, b)
  mc_var <- var(samples)
  quad_var <- var_w1(a, b)

  # Allow 10% relative error for MC variance estimate
  expect_equal(quad_var, mc_var, tolerance = 0.1)
})


# =============================================================================
# Test 11: summary_w1 Output Structure
# =============================================================================

test_that("summary_w1 returns correct structure", {
  s <- summary_w1(a = 2, b = 1)

  expect_s3_class(s, "w1_summary")
  expect_named(s, c("mean", "var", "sd", "median", "quantiles",
                    "prob_gt_50", "prob_gt_90", "params", "estimand",
                    "label", "conditioning", "provenance"))

  expect_true(s$mean > 0 && s$mean < 1)
  expect_true(s$var > 0)
  expect_equal(s$sd, sqrt(s$var))
  expect_true(s$median > 0 && s$median < 1)
  expect_length(s$quantiles, 5)  # Default 5 quantiles
  expect_true(s$prob_gt_50 >= 0 && s$prob_gt_50 <= 1)
  expect_true(s$prob_gt_90 >= 0 && s$prob_gt_90 <= 1)
  expect_equal(s$params, list(a = 2, b = 1))
})


# =============================================================================
# Test 12: rw1 Random Generation
# =============================================================================

test_that("rw1 generates valid samples", {
  set.seed(42)
  samples <- rw1(1000, a = 2, b = 1)

  expect_length(samples, 1000)
  expect_true(all(samples > 0))
  expect_true(all(samples < 1))
})


test_that("rw1 validates input", {
  expect_error(rw1(-1, a = 2, b = 1))
  expect_error(rw1(0, a = 2, b = 1))
  expect_error(rw1(1.5, a = 2, b = 1))
  expect_error(rw1(100, a = -1, b = 1))
  expect_error(
    rw1(10, a = 1e-300, b = 1),
    class = "dpprior_weight_rng_error"
  )
})


# =============================================================================
# Test 13: Log-Density Consistency
# =============================================================================

test_that("log-density is consistent with density", {
  a <- 2; b <- 1
  x <- c(0.1, 0.3, 0.5, 0.7)

  log_dens <- density_w1(x, a, b, log = TRUE)
  dens <- density_w1(x, a, b, log = FALSE)

  expect_equal(log_dens, log(dens), tolerance = 1e-10)
})


# =============================================================================
# Test 14: Numerical Stability - expm1 Test
# =============================================================================

test_that("CDF is stable for very small x (expm1 test)", {
  a <- 2; b <- 1

  # For very small x, CDF should be close to 0
  # Using expm1 is critical here
  x_small <- c(1e-10, 1e-8, 1e-6, 1e-4)
  cdf_vals <- cdf_w1(x_small, a, b)

  # All should be positive and increasing
  expect_true(all(cdf_vals > 0))
  expect_true(all(diff(cdf_vals) > 0))

  # Check against naive calculation for moderate x
  # where both methods should agree
  x_moderate <- 0.1
  cdf_stable <- cdf_w1(x_moderate, a, b)

  # Naive calculation (less stable)
  denom <- b - log1p(-x_moderate)
  cdf_naive <- 1 - (b / denom)^a

  expect_equal(cdf_stable, cdf_naive, tolerance = 1e-10)
})


# =============================================================================
# Test 15: Consistency Check - CDF Derivation
# =============================================================================

test_that("CDF matches numerical integration of conditional", {
  # Verify: F(x) = ∫ [1 - (1-x)^α] × p(α) dα
  test_cases <- list(
    c(a = 2.0, b = 1.0, x = 0.3),
    c(a = 1.6, b = 1.22, x = 0.5),
    c(a = 0.5, b = 1.0, x = 0.3)
  )

  for (case in test_cases) {
    a <- case["a"]
    b <- case["b"]
    x <- case["x"]

    # Numerical integration
    cdf_numerical <- integrate(
      function(alpha) (1 - (1 - x)^alpha) * dgamma(alpha, shape = a, rate = b),
      lower = 0, upper = 100
    )$value

    # Closed form
    cdf_closed <- cdf_w1(x, a, b)

    expect_equal(cdf_numerical, cdf_closed, tolerance = 1e-6,
                 info = sprintf("Derivation check: a=%.1f, b=%.2f, x=%.1f",
                                a, b, x))
  }
})


# =============================================================================
# Phase 5 contract tests: W_SB semantics and W_max diagnostics
# =============================================================================

test_that("W_SB public metadata and print labels name the estimand exactly", {
  expect_equal(
    prob_wsb_exceeds(0.5, a = 2, b = 1),
    prob_w1_exceeds(0.5, a = 2, b = 1),
    tolerance = 0
  )

  s <- summary_w1(a = 2, b = 1)
  expect_identical(s$estimand, "W_SB")
  expect_identical(s$label, "First size-biased DP weight")
  expect_identical(s$conditioning, "gamma_mixed")
  expect_identical(s$provenance$gamma_parameterization, "shape_rate")

  printed <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(printed, "First size-biased DP weight")
  expect_match(printed, "P\\(W_SB > 0.5\\)")
  expect_false(any(grepl("dominance|largest", printed, ignore.case = TRUE)))
})


test_that("W_SB closed forms remain stable on endpoint and adversarial grids", {
  ab_grid <- expand.grid(
    a = c(0.1, 0.5, 2.5, 10),
    b = c(0.2, 1, 5)
  )
  thresholds <- c(0, 1e-12, 0.5, 0.9, 1 - 1e-12, 1)

  for (i in seq_len(nrow(ab_grid))) {
    a <- ab_grid$a[i]
    b <- ab_grid$b[i]
    actual <- prob_wsb_exceeds(thresholds, a, b)
    c_t <- -log1p(-thresholds[thresholds > 0 & thresholds < 1])
    expected <- c(
      1,
      exp(-a * log1p(c_t / b)),
      0
    )
    expect_equal(actual, expected, tolerance = 2e-14)
    expect_true(all(diff(actual) <= 0))
  }

  expect_equal(cdf_w1(c(-Inf, 0, 1, Inf), 2, 1), c(0, 0, 1, 1))
  expect_equal(prob_wsb_exceeds(c(-Inf, 0, 1, Inf), 2, 1), c(1, 1, 0, 0))
  expect_error(cdf_w1(0.5, c(1, 2), 1), class = "dpprior_length_error")
  expect_error(mean_w1(c(1, 2), 1), class = "dpprior_length_error")
})


test_that("subnormal Gamma parameters do not overflow an intermediate ratio", {
  tiny <- .Machine$double.xmin * .Machine$double.eps
  c_t <- -log1p(-0.5)
  log_ratio <- log(c_t) - log(tiny)
  log1p_ratio <- log_ratio + log1p(exp(-log_ratio))
  reference_log_tail <- -tiny * log1p_ratio

  tail <- prob_wsb_exceeds(0.5, a = tiny, b = tiny)
  cdf <- cdf_w1(0.5, a = tiny, b = tiny)
  density <- density_w1(0.5, a = tiny, b = tiny)
  bounds <- wmax_tail_bounds(0.5, a = tiny, b = tiny)

  expect_equal(tail, exp(reference_log_tail), tolerance = 0)
  expect_gt(tail, 0.999999999999)
  expect_true(is.finite(cdf))
  expect_gt(cdf, 0)
  expect_true(is.finite(density))
  expect_gt(density, 0)
  expect_gt(bounds$lower_bound, 0.999999999999)
  expect_equal(bounds$log_lower_bound, reference_log_tail, tolerance = 0)
})


test_that("W_SB variance remains nonnegative for concentrated Gamma priors", {
  a <- 1e-8
  b <- 1e8
  variance <- var_w1(a = a, b = b, M = 256L)
  summary <- summary_w1(a = a, b = b, M = 256L)

  # For alpha near zero, E[Var(W_SB | alpha)] is asymptotic to E[alpha]/2.
  leading_term <- (a / b) / 2
  expect_true(is.finite(variance))
  expect_gte(variance, 0)
  expect_lt(abs(variance - leading_term), 1e-24)
  expect_true(is.finite(summary$sd))
  expect_lt(abs(summary$sd^2 - summary$var), 1e-30)

  # Avoid forming the cubic-size conditional-variance denominator. With a
  # concentrated alpha near 1e110, the true variance is still representable.
  large_a <- 100
  large_b <- 1e-108
  large_variance <- var_w1(a = large_a, b = large_b, M = 256L)
  reciprocal_alpha_limit <-
    large_b^2 * large_a / ((large_a - 1)^2 * (large_a - 2))
  expect_true(is.finite(large_variance))
  expect_gt(large_variance, 0)
  expect_lt(
    abs(large_variance - reciprocal_alpha_limit) / reciprocal_alpha_limit,
    1e-10
  )
})


test_that("certified W_max bounds are distinct from the W_SB tail", {
  fixed <- wmax_tail_bounds(0.5, alpha = 1)
  expect_s3_class(fixed, "wmax_tail_bounds")
  expect_identical(fixed$estimand, "W_max")
  expect_identical(fixed$conditioning, "fixed_alpha")
  expect_true(fixed$certified)
  expect_lte(fixed$lower_bound, exp(fixed$log_lower_bound))
  expect_gte(fixed$upper_bound, exp(fixed$log_upper_bound))
  expect_equal(fixed$size_biased_tail, 0.5, tolerance = 1e-15)
  expect_equal(fixed$lower_bound, 0.5, tolerance = 1e-15)
  expect_equal(fixed$upper_bound, 1, tolerance = 1e-15)

  mixed <- wmax_tail_bounds(0.4, a = 2, b = 1)
  q <- prob_wsb_exceeds(0.4, 2, 1)
  expect_lte(mixed$lower_bound, q)
  expect_gte(mixed$upper_bound, min(q / 0.4, 1))
  expect_lte(q - mixed$lower_bound, 2 * .Machine$double.eps)
  expect_identical(mixed$conditioning, "gamma_mixed")
  expect_identical(mixed$method, "size_biased_mass_identity")
})


test_that("fixed-alpha W_max tail passes exact half-threshold identities", {
  result <- prob_wmax_exceeds(0.5, alpha = 1)

  expect_s3_class(result, "wmax_tail_result")
  expect_identical(result$estimand, "W_max")
  expect_identical(result$conditioning, "fixed_alpha")
  expect_identical(result$status, "converged")
  expect_true(result$usable)
  expect_true(result$bounds_usable)
  expect_true(result$verified)
  expect_true(result$estimate_usable)
  expect_equal(result$estimate, log(2), tolerance = 2e-12)
  expect_lte(result$lower, log(2))
  expect_gte(result$upper, log(2))
  expect_lte(result$lower_bound, result$lower)
  expect_gte(result$upper_bound, result$upper)
  expect_lte(result$abs_error_bound, 1e-8)
  expect_gt(result$estimate, result$lower_bound)
  expect_lte(result$estimate, result$upper_bound)
  expect_identical(
    result$method,
    "one_dimensional_quadrature_verified_by_positive_series"
  )
})


test_that("Gamma-mixed W_max tail is independently verified", {
  cases <- data.frame(
    a = c(0.5, 1, 2, 5),
    b = c(1, 1, 1, 2),
    threshold = c(0.5, 0.9, 0.5, 0.9),
    reference = c(
      0.8693316004509434,
      0.3092400892506535,
      0.4813041443928161,
      0.0228958329067216
    )
  )

  for (i in seq_len(nrow(cases))) {
    z <- prob_wmax_exceeds(
      cases$threshold[i], a = cases$a[i], b = cases$b[i]
    )
    expect_identical(z$conditioning, "gamma_mixed")
    expect_identical(z$status, "converged")
    expect_true(z$verified)
    expect_equal(z$estimate, cases$reference[i], tolerance = 2e-9)
    expect_lte(z$lower_bound, z$estimate + 1e-10)
    expect_gte(z$upper_bound, z$estimate - 1e-10)
    expect_lte(z$lower_bound, z$lower)
    expect_lte(z$lower, z$estimate)
    expect_lte(z$estimate, z$upper)
    expect_lte(z$upper, z$upper_bound)
    expect_lte(z$abs_error_bound, 1e-8)
    expect_identical(z$numerical$independent_method,
                     "gamma_quadrature_of_conditional_positive_series")
    expect_identical(z$numerical$M_selected, 256L)
    expect_identical(z$numerical$M_verification, 512L)
  }
})


test_that("thresholds below one half never receive a silent direct estimate", {
  auto <- prob_wmax_exceeds(0.4, alpha = 1)
  expect_identical(auto$status, "approximate")
  expect_false(auto$usable)
  expect_true(auto$bounds_usable)
  expect_false(auto$verified)
  expect_false(auto$estimate_usable)
  expect_true(is.na(auto$estimate))
  expect_identical(auto$method, "certified_bounds_only")
  expect_identical(auto$reason, "deterministic_method_unsupported_below_half")
  expect_true(auto$provenance$general_threshold_exact_method_deferred)

  explicit <- prob_wmax_exceeds(0.4, alpha = 1, method = "deterministic")
  expect_identical(explicit$status, "approximate")
  expect_false(explicit$usable)
  expect_true(explicit$bounds_usable)
  expect_true(is.na(explicit$estimate))
  expect_equal(explicit$lower_bound, auto$lower_bound, tolerance = 0)
  expect_equal(explicit$upper_bound, auto$upper_bound, tolerance = 0)

  requested_bounds <- prob_wmax_exceeds(0.4, alpha = 1, method = "bounds")
  expect_identical(requested_bounds$status, "approximate")
  expect_false(requested_bounds$usable)
  expect_true(requested_bounds$bounds_usable)
  expect_identical(requested_bounds$reason, "bounds_requested")
})


test_that("W_max probability boundaries short-circuit every requested method", {
  for (method in c("auto", "deterministic", "bounds", "monte_carlo")) {
    lower <- expect_no_warning(
      prob_wmax_exceeds(0, alpha = 1, method = method)
    )
    upper <- expect_no_warning(
      prob_wmax_exceeds(1, a = 2, b = 1, method = method)
    )
    expect_identical(lower$status, "converged")
    expect_identical(upper$status, "converged")
    expect_equal(lower$estimate, 1)
    expect_equal(upper$estimate, 0)
    expect_identical(lower$method, "probability_boundary_identity")
    expect_identical(upper$method, "probability_boundary_identity")
    expect_identical(lower$provenance$requested_method, method)
    expect_identical(upper$provenance$requested_method, method)
  }
})


test_that("large shape-rate Gamma mixing avoids logarithmic cancellation", {
  large <- exp(15)
  mixed <- prob_wmax_exceeds(0.5, a = large, b = large)
  fixed <- prob_wmax_exceeds(0.5, alpha = 1)

  expect_identical(mixed$status, "converged")
  expect_true(mixed$verified)
  expect_true(is.finite(mixed$estimate))
  expect_equal(mixed$estimate, fixed$estimate, tolerance = 2e-7)
  expect_lte(mixed$lower_bound, mixed$lower)
  expect_lte(mixed$lower, mixed$estimate)
  expect_lte(mixed$estimate, mixed$upper)
  expect_lte(mixed$upper, mixed$upper_bound)
})


test_that("strictly positive underflow tails retain conservative log bounds", {
  bounds <- wmax_tail_bounds(0.5, alpha = 2000)
  expect_true(is.finite(bounds$log_lower_bound))
  expect_true(is.finite(bounds$log_upper_bound))
  expect_true(bounds$probability_scale_underflow)
  expect_true(bounds$upper_scale_underflow)
  expect_equal(bounds$lower_bound, 0)
  expect_gt(bounds$upper_bound, 0)
  expect_identical(
    bounds$provenance$probability_representation,
    "log_bounds_with_conservative_smallest_positive_real_ceiling"
  )

  direct <- prob_wmax_exceeds(0.5, alpha = 2000)
  expect_identical(direct$status, "approximate")
  expect_false(direct$usable)
  expect_true(direct$bounds_usable)
  expect_false(direct$verified)
  expect_identical(direct$reason, "probability_scale_underflow")
  expect_true(all(is.na(c(
    direct$estimate, direct$lower, direct$upper, direct$abs_error_bound
  ))))
  expect_true("candidate" %in% names(direct$numerical))
  expect_false(direct$numerical$candidate_published)
})


test_that("certified ordinary bounds are rounded outward near probability one", {
  bounds <- wmax_tail_bounds(0.5, a = 1, b = 1e308)

  expect_lt(bounds$log_lower_bound, 0)
  expect_equal(bounds$size_biased_tail, 1)
  expect_lt(bounds$lower_bound, 1)
  expect_equal(bounds$upper_bound, 1)
  expect_true(bounds$lower_outward_rounded)
  expect_lte(log(bounds$lower_bound), bounds$log_lower_bound)
  expect_identical(
    bounds$provenance$probability_representation,
    "outward_rounded_ordinary_bounds_with_exact_log_bounds"
  )
})


test_that("failed deterministic verification withholds its raw candidate", {
  result <- prob_wmax_exceeds(0.5, a = 1, b = 1e308)

  expect_identical(result$status, "approximate")
  expect_false(result$verified)
  expect_false(result$usable)
  expect_false(result$estimate_usable)
  expect_true(result$bounds_usable)
  expect_identical(result$reason, "independent_verification_disagreement")
  expect_true(all(is.na(c(
    result$estimate, result$lower, result$upper, result$abs_error_bound
  ))))
  expect_lt(result$lower_bound, 1)
  expect_equal(result$upper_bound, 1)
  expect_true(is.finite(result$numerical$candidate))
  expect_false(result$numerical$candidate_published)

  printed <- paste(capture.output(print(result)), collapse = "\n")
  expect_match(printed, "Estimate: unavailable", fixed = TRUE)
  expect_match(printed, "independent_verification_disagreement", fixed = TRUE)
  expect_false(grepl(format(result$numerical$candidate, digits = 6),
                     printed, fixed = TRUE))
})


test_that("seeded W_max Monte Carlo is explicit and exactly stopped", {
  z1 <- prob_wmax_exceeds(
    0.4, alpha = 1, method = "monte_carlo",
    n = 5000, seed = 20260819, warn_low_successes = FALSE
  )
  z2 <- prob_wmax_exceeds(
    0.4, alpha = 1, method = "monte_carlo",
    n = 5000, seed = 20260819, warn_low_successes = FALSE
  )

  expect_identical(z1$estimate, z2$estimate)
  expect_identical(z1$sampling$successes, z2$sampling$successes)
  expect_identical(z1$status, "approximate")
  expect_false(z1$usable)
  expect_true(z1$bounds_usable)
  expect_true(z1$verified)
  expect_false(z1$estimate_usable)
  expect_identical(z1$method, "seeded_gem_monte_carlo_exact_stopping")
  expect_identical(z1$sampling$seed, 20260819L)
  expect_identical(z1$sampling$n, 5000L)
  expect_identical(z1$sampling$unresolved, 0L)
  expect_identical(z1$sampling$stopping_rule,
                   "remainder_not_greater_than_current_maximum")
  expect_true(z1$sampling$wilson_interval[1] <= z1$estimate)
  expect_true(z1$sampling$wilson_interval[2] >= z1$estimate)
  expect_true(z1$sampling$interval_intersection_nonempty)
  expect_true(z1$sampling$interval_contains_estimate)
  expect_lte(z1$lower_bound, z1$lower)
  expect_lte(z1$lower, z1$estimate)
  expect_lte(z1$estimate, z1$upper)
  expect_lte(z1$upper, z1$upper_bound)
  expect_lte(
    z1$sampling$hoeffding_interval_unintersected[["lower"]], z1$lower
  )
  expect_gte(
    z1$sampling$hoeffding_interval_unintersected[["upper"]], z1$upper
  )
  expect_lte(z1$lower_bound, z1$estimate + z1$abs_error_bound)
  expect_gte(z1$upper_bound, z1$estimate - z1$abs_error_bound)
})


test_that("failed exact-stopping MC remains unusable while bounds survive", {
  failed <- prob_wmax_exceeds(
    0.4, alpha = 100, method = "monte_carlo",
    n = 50, seed = 20260819, max_sticks = 1,
    warn_low_successes = FALSE
  )
  expect_identical(failed$status, "failed")
  expect_false(failed$usable)
  expect_true(failed$bounds_usable)
  expect_identical(failed$reason, "exact_stopping_limit_reached")
  expect_gt(failed$sampling$unresolved, 0L)
  expect_true(is.na(failed$estimate))
})


test_that("MC estimates outside certified bounds are explicit and unusable", {
  outside <- prob_wmax_exceeds(
    0.9, alpha = 1, method = "monte_carlo",
    n = 1, seed = 20260819, warn_low_successes = FALSE
  )
  expect_identical(outside$status, "approximate")
  expect_false(outside$usable)
  expect_true(outside$bounds_usable)
  expect_false(outside$verified)
  expect_false(outside$estimate_usable)
  expect_identical(
    outside$reason, "sampling_estimate_outside_certified_interval"
  )
  expect_true(outside$sampling$interval_intersection_nonempty)
  expect_false(outside$sampling$interval_contains_estimate)
  expect_true(is.na(outside$estimate))
  expect_true(is.na(outside$lower))
  expect_true(is.na(outside$upper))
  expect_true(is.na(outside$abs_error_bound))
  expect_equal(outside$sampling$raw_estimate, 0)
  outside_text <- paste(capture.output(print(outside)), collapse = "\n")
  expect_match(outside_text, "Estimate: unavailable", fixed = TRUE)
  expect_match(outside_text, "sampling_estimate_outside_certified_interval",
               fixed = TRUE)
})


test_that("W_max API rejects ambiguous specifications with typed conditions", {
  expect_error(
    wmax_tail_bounds(0.5, alpha = 1, a = 2, b = 1),
    class = "dpprior_weight_specification_error"
  )
  expect_error(
    wmax_tail_bounds(0.5, a = 2),
    class = "dpprior_weight_specification_error"
  )
  expect_error(
    wmax_tail_bounds(-0.1, alpha = 1),
    class = "dpprior_bounds_error"
  )
  expect_error(
    prob_wmax_exceeds(0.4, alpha = 1, method = "monte_carlo", n = 1000),
    class = "dpprior_weight_mc_error"
  )
  expect_error(
    prob_wmax_exceeds(0.5, alpha = c(1, 2)),
    class = "dpprior_length_error"
  )
})


test_that("W_max print methods expose status, method, and conditioning", {
  direct <- prob_wmax_exceeds(0.5, alpha = 1)
  direct_text <- paste(capture.output(print(direct)), collapse = "\n")
  expect_match(direct_text, "Status: converged")
  expect_match(direct_text, "Verified: yes")
  expect_match(direct_text, "W_max")
  expect_match(direct_text, "fixed alpha")
  expect_false(grepl("W_SB.*largest|dominance", direct_text,
                     ignore.case = TRUE))

  bounds <- wmax_tail_bounds(0.4, a = 2, b = 1)
  bounds_text <- paste(capture.output(print(bounds)), collapse = "\n")
  expect_match(bounds_text, "Certified bounds")
  expect_match(bounds_text, "Gamma")
  expect_match(bounds_text, "W_SB")
})
