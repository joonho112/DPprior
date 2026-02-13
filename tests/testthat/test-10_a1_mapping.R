# =============================================================================
# Tests for Module 10: A1 Closed-Form Prior Elicitation
# =============================================================================
#
# Author: JoonHo Lee
# Date: December 2025
# Part of: DPprior R Package
# =============================================================================

# Note: context() is deprecated in testthat 3rd edition


# =============================================================================
# Test: compute_scaling_constant
# =============================================================================

test_that("log scaling returns log(J)", {
  expect_equal(compute_scaling_constant(50, "log"), log(50))
  expect_equal(compute_scaling_constant(100, "log"), log(100))
  expect_equal(compute_scaling_constant(2, "log"), log(2))
})

test_that("harmonic scaling returns H_{J-1} = digamma(J) + gamma", {
  euler_gamma <- 0.5772156649015329
  expect_equal(compute_scaling_constant(50, "harmonic"),
               digamma(50) + euler_gamma)
  expect_equal(compute_scaling_constant(10, "harmonic"),
               digamma(10) + euler_gamma)
})

test_that("digamma scaling requires mu_K argument", {
  expect_error(compute_scaling_constant(50, "digamma"),
               "mu_K required for digamma scaling")
  # Should succeed when mu_K is provided
  result <- compute_scaling_constant(50, "digamma", mu_K = 5)
  expect_true(is.finite(result))
  expect_true(result > 0)
})

test_that("all scalings converge toward log(J) for large J", {
  J_large <- 500
  cJ_log <- compute_scaling_constant(J_large, "log")
  cJ_harm <- compute_scaling_constant(J_large, "harmonic")
  cJ_dig <- compute_scaling_constant(J_large, "digamma", mu_K = 5)

  # Harmonic and digamma should be within ~15% of log(J) for large J
  # The harmonic sum H_{J-1} converges to log(J) + gamma, so the relative

  # difference shrinks with J but does not vanish entirely
  expect_equal(cJ_harm / cJ_log, 1, tolerance = 0.15)
  expect_equal(cJ_dig / cJ_log, 1, tolerance = 0.25)
})

test_that("digamma scaling uses correct formula", {
  J <- 50
  mu_K <- 5
  alpha_tilde <- (mu_K - 1) / log(J)
  expected <- digamma(alpha_tilde + J) - digamma(alpha_tilde)
  expect_equal(compute_scaling_constant(J, "digamma", mu_K = mu_K), expected)
})


# =============================================================================
# Test: DPprior_a1 basic functionality
# =============================================================================

test_that("DPprior_a1 returns valid positive parameters", {
  fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
  expect_true(fit$a > 0)
  expect_true(fit$b > 0)
})

test_that("DPprior_a1 returns correct structure", {
  fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)

  # Check all required components are present
  expect_true("a" %in% names(fit))
  expect_true("b" %in% names(fit))
  expect_true("J" %in% names(fit))
  expect_true("target" %in% names(fit))
  expect_true("method" %in% names(fit))
  expect_true("status" %in% names(fit))
  expect_true("scaling" %in% names(fit))
  expect_true("cJ" %in% names(fit))
  expect_true("var_K_used" %in% names(fit))
  expect_true("converged" %in% names(fit))
  expect_true("iterations" %in% names(fit))

  # Check values
  expect_equal(fit$method, "A1")
  expect_equal(fit$J, 50)
  expect_equal(fit$target$mu_K, 5)
  expect_equal(fit$target$var_K, 8)
  expect_equal(fit$target$type, "moments")
  expect_true(fit$converged)
  expect_identical(fit$iterations, 0L)
})

test_that("DPprior_a1 output has class DPprior_fit", {
  fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
  expect_s3_class(fit, "DPprior_fit")
})

test_that("DPprior_a1 round-trip: NegBin moments recover targets", {
  fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
  a <- fit$a
  b <- fit$b
  cJ <- fit$cJ

  # Forward NegBin model
  p <- b / (b + cJ)
  mu_S <- a * (1 - p) / p
  var_S <- a * (1 - p) / p^2

  mu_K_recovered <- mu_S + 1
  var_K_recovered <- var_S

  expect_equal(mu_K_recovered, 5, tolerance = 1e-10)
  expect_equal(var_K_recovered, 8, tolerance = 1e-10)
})

test_that("DPprior_a1 works with all scaling methods", {
  for (scaling in c("log", "harmonic", "digamma")) {
    fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8, scaling = scaling)
    expect_true(fit$a > 0)
    expect_true(fit$b > 0)
    expect_equal(fit$scaling, scaling)
    expect_equal(fit$status, "success")
  }
})

test_that("DPprior_a1 different scalings give different parameters", {
  fit_log <- DPprior_a1(J = 50, mu_K = 5, var_K = 8, scaling = "log")
  fit_harm <- DPprior_a1(J = 50, mu_K = 5, var_K = 8, scaling = "harmonic")

  # Different cJ values should lead to different b (but same a)
  expect_false(identical(fit_log$cJ, fit_harm$cJ))
  expect_false(identical(fit_log$b, fit_harm$b))
  # Shape a should be identical since it depends only on mu_K and var_K
  expect_equal(fit_log$a, fit_harm$a, tolerance = 1e-10)
})


# =============================================================================
# Test: DPprior_a1 input validation
# =============================================================================

test_that("DPprior_a1 errors on invalid J", {
  expect_error(DPprior_a1(J = 1, mu_K = 5, var_K = 8),
               "J must be an integer >= 2")
  expect_error(DPprior_a1(J = 0, mu_K = 5, var_K = 8),
               "J must be an integer >= 2")
  expect_error(DPprior_a1(J = -5, mu_K = 5, var_K = 8),
               "J must be an integer >= 2")
  expect_error(DPprior_a1(J = 3.5, mu_K = 2, var_K = 3),
               "J must be an integer >= 2")
})

test_that("DPprior_a1 errors on invalid mu_K", {
  expect_error(DPprior_a1(J = 50, mu_K = 0.5, var_K = 8),
               "mu_K must be > 1")
  expect_error(DPprior_a1(J = 50, mu_K = 1, var_K = 8),
               "mu_K must be > 1")
  expect_error(DPprior_a1(J = 50, mu_K = 51, var_K = 8),
               "mu_K must be <= J")
})

test_that("DPprior_a1 errors on invalid var_K", {
  expect_error(DPprior_a1(J = 50, mu_K = 5, var_K = 0),
               "var_K must be a positive finite numeric scalar")
  expect_error(DPprior_a1(J = 50, mu_K = 5, var_K = -1),
               "var_K must be a positive finite numeric scalar")
})


# =============================================================================
# Test: DPprior_a1 variance feasibility (boundary projection)
# =============================================================================

test_that("infeasible variance triggers projection with warning", {
  # mu_K = 5 means mu_S = 4, so var_K must be > 4 for feasibility
  # var_K = 3 < 4 is infeasible
  expect_warning(
    fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 3),
    "var_K <= mu_K - 1: projected to feasible boundary"
  )
  expect_true(grepl("projected", fit$status))
  expect_true(fit$var_K_used > fit$target$var_K)
  expect_true(fit$a > 0)
  expect_true(fit$b > 0)
})

test_that("feasible variance does not trigger projection", {
  fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
  expect_equal(fit$status, "success")
  expect_equal(fit$var_K_used, 8)
})

test_that("var_K exactly at boundary triggers projection", {
  # mu_K = 5 means mu_S = 4, boundary is var_K = mu_S = 4
  # Due to epsilon buffering, var_K = 4 should be projected
  expect_warning(
    fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 4),
    "var_K <= mu_K - 1: projected to feasible boundary"
  )
  expect_true(grepl("projected", fit$status))
})


# =============================================================================
# Test: DPprior_a1 closed-form algebra
# =============================================================================

test_that("A1 closed-form formulas are algebraically correct", {
  J <- 50
  mu_K <- 5
  var_K <- 8

  cJ <- log(J)
  m <- mu_K - 1       # shifted mean
  D <- var_K - m       # denominator

  a_expected <- m^2 / D
  b_expected <- m * cJ / D

  fit <- DPprior_a1(J = J, mu_K = mu_K, var_K = var_K, scaling = "log")
  expect_equal(fit$a, a_expected, tolerance = 1e-12)
  expect_equal(fit$b, b_expected, tolerance = 1e-12)
})


# =============================================================================
# Test: vif_to_variance
# =============================================================================

test_that("VIF = 1 gives Poisson baseline variance (mu_K - 1)", {
  expect_equal(vif_to_variance(mu_K = 5, vif = 1), 4)
  expect_equal(vif_to_variance(mu_K = 10, vif = 1), 9)
})

test_that("VIF = 2 gives twice the Poisson baseline variance", {
  expect_equal(vif_to_variance(mu_K = 5, vif = 2), 8)
  expect_equal(vif_to_variance(mu_K = 10, vif = 2), 18)
})

test_that("VIF < 1 errors", {
  expect_error(vif_to_variance(mu_K = 5, vif = 0.5),
               "vif must be >= 1")
  expect_error(vif_to_variance(mu_K = 5, vif = 0),
               "vif must be >= 1")
})

test_that("vif_to_variance errors on invalid mu_K", {
  expect_error(vif_to_variance(mu_K = 0.5, vif = 2),
               "mu_K must be > 1")
})


# =============================================================================
# Test: confidence_to_vif
# =============================================================================

test_that("confidence levels map to decreasing VIF (low > medium > high)", {
  vif_low <- confidence_to_vif("low")
  vif_med <- confidence_to_vif("medium")
  vif_high <- confidence_to_vif("high")

  expect_true(vif_low > vif_med)
  expect_true(vif_med > vif_high)
})

test_that("confidence_to_vif returns known values", {
  expect_equal(confidence_to_vif("low"), 5.0)
  expect_equal(confidence_to_vif("medium"), 2.5)
  expect_equal(confidence_to_vif("high"), 1.5)
})

test_that("confidence_to_vif returns scalar numeric", {
  result <- confidence_to_vif("medium")
  expect_true(is.numeric(result))
  expect_length(result, 1)
})

test_that("confidence_to_vif rejects invalid input", {
  expect_error(confidence_to_vif("extreme"))
})


# =============================================================================
# Test: cv_alpha_to_variance
# =============================================================================

test_that("cv_alpha_to_variance follows correct algebra", {
  mu_K <- 5
  cv_alpha <- 0.5
  m <- mu_K - 1
  expected <- m + (cv_alpha * m)^2   # m(1 + cv^2 * m) = m + cv^2 * m^2
  expect_equal(cv_alpha_to_variance(mu_K, cv_alpha), expected)
})

test_that("larger cv_alpha gives larger variance", {
  var1 <- cv_alpha_to_variance(mu_K = 5, cv_alpha = 0.5)
  var2 <- cv_alpha_to_variance(mu_K = 5, cv_alpha = 1.0)
  var3 <- cv_alpha_to_variance(mu_K = 5, cv_alpha = 2.0)

  expect_true(var1 < var2)
  expect_true(var2 < var3)
})

test_that("cv_alpha round-trip: recovered CV matches target", {
  mu_K <- 5
  cv_target <- 0.5

  var_K <- cv_alpha_to_variance(mu_K, cv_target)
  fit <- DPprior_a1(J = 50, mu_K = mu_K, var_K = var_K)
  cv_recovered <- 1 / sqrt(fit$a)

  expect_equal(cv_recovered, cv_target, tolerance = 1e-8)
})

test_that("cv_alpha_to_variance errors on invalid inputs", {
  expect_error(cv_alpha_to_variance(mu_K = 0.5, cv_alpha = 1),
               "mu_K must be > 1")
  expect_error(cv_alpha_to_variance(mu_K = 5, cv_alpha = -0.5),
               "cv_alpha must be positive")
  expect_error(cv_alpha_to_variance(mu_K = 5, cv_alpha = 0),
               "cv_alpha must be positive")
})


# =============================================================================
# Test: compare_a1_a2
# =============================================================================

test_that("compare_a1_a2 returns list with a1 and a2 components", {
  result <- compare_a1_a2(J = 50, mu_K = 5, var_K = 8, verbose = FALSE)

  expect_true(is.list(result))
  expect_true("a1" %in% names(result))
  expect_true("a2" %in% names(result))
  expect_true("improvement_ratio" %in% names(result))

  # Both components have expected fields
  expect_true("a" %in% names(result$a1))
  expect_true("b" %in% names(result$a1))
  expect_true("residual" %in% names(result$a1))
  expect_true("a" %in% names(result$a2))
  expect_true("b" %in% names(result$a2))
  expect_true("residual" %in% names(result$a2))
})

test_that("A2 achieves smaller or equal residual than A1", {
  result <- compare_a1_a2(J = 50, mu_K = 5, var_K = 8, verbose = FALSE)
  expect_true(result$a2$residual <= result$a1$residual)
  expect_true(result$improvement_ratio >= 1)
})

test_that("compare_a1_a2 works for various parameter combos", {
  # Moderate case

  result <- compare_a1_a2(J = 100, mu_K = 10, var_K = 20, verbose = FALSE)
  expect_true(is.list(result))
  expect_true(result$a1$a > 0)
  expect_true(result$a2$a > 0)
})


# =============================================================================
# Test: S3 methods (print, summary, as.data.frame)
# =============================================================================

test_that("print.DPprior_fit runs without error", {
  fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
  expect_output(print(fit), "DPprior")
  expect_output(print(fit), "Method: A1")
})

test_that("summary.DPprior_fit returns correct structure", {
  fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
  s <- summary(fit, print_output = FALSE)

  expect_true(is.list(s))
  expect_equal(s$method, "A1")
  expect_equal(s$gamma_prior$a, fit$a)
  expect_equal(s$gamma_prior$b, fit$b)
  # Module 17 summary uses E_alpha and CV_alpha field names
  expect_equal(s$alpha_summary$E_alpha, fit$a / fit$b)
  expect_equal(s$alpha_summary$CV_alpha, 1 / sqrt(fit$a))
  expect_true(s$converged)
  expect_equal(s$iterations, 0L)
})

test_that("as.data.frame.DPprior_fit returns one-row data.frame", {
  fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
  df <- as.data.frame(fit)

  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 1)
  expect_equal(df$method, "A1")
  expect_equal(df$a, fit$a)
  expect_equal(df$b, fit$b)
  expect_equal(df$J, 50)
  expect_equal(df$mu_K, 5)
  expect_equal(df$var_K, 8)
})


# =============================================================================
# Test: verify_a1_roundtrip (internal verification function)
# =============================================================================

test_that("verify_a1_roundtrip passes for feasible cases", {
  fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
  expect_true(verify_a1_roundtrip(fit, verbose = FALSE))

  fit2 <- DPprior_a1(J = 100, mu_K = 10, var_K = 20)
  expect_true(verify_a1_roundtrip(fit2, verbose = FALSE))
})

test_that("verify_a1_roundtrip passes for projected cases", {
  fit <- suppressWarnings(DPprior_a1(J = 50, mu_K = 5, var_K = 3))
  # Round-trip should still pass using var_K_used (projected value)
  expect_true(verify_a1_roundtrip(fit, verbose = FALSE))
})
