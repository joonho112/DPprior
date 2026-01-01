# =============================================================================
# Test Suite: Module 02 - Gauss-Laguerre Quadrature
# =============================================================================
#
# Tests for numerical integration against Gamma distributions using
# generalized Gauss-Laguerre quadrature.
#
# Note on accuracy:
# - Polynomial integrands: machine precision (< 1e-10)
# - Non-polynomial smooth functions: 1e-4 to 1e-6 range
# - Functions with singularities near 0: depends on shape parameter
#
# Author: JoonHo Lee
# Date: December 2025
# Part of: DPprior R Package
# =============================================================================

test_that("gauss_laguerre_nodes returns correct structure", {
  quad <- gauss_laguerre_nodes(40, alpha_param = 0)

  expect_type(quad, "list")
  expect_named(quad, c("nodes", "weights", "weights_log", "weights_norm"))
  expect_length(quad$nodes, 40)
  expect_length(quad$weights, 40)
  expect_length(quad$weights_log, 40)
  expect_length(quad$weights_norm, 40)

  # weights_norm should sum to 1
  expect_equal(sum(quad$weights_norm), 1, tolerance = 1e-14)
})


test_that("gauss_laguerre_nodes handles M=1 correctly", {
  quad <- gauss_laguerre_nodes(1, alpha_param = 0)

  expect_length(quad$nodes, 1)
  expect_equal(quad$nodes, 1)
  expect_equal(quad$weights, 1, tolerance = 1e-15)
})


test_that("gauss_laguerre_nodes produces positive nodes", {
  test_cases <- list(
    list(M = 80, alpha_param = 0),
    list(M = 80, alpha_param = 1.5),
    list(M = 80, alpha_param = -0.5),
    list(M = 100, alpha_param = 4.0)
  )

  for (tc in test_cases) {
    quad <- gauss_laguerre_nodes(tc$M, tc$alpha_param)
    expect_true(all(quad$nodes > 0),
                info = sprintf("M=%d, alpha_param=%.1f", tc$M, tc$alpha_param))
  }
})


test_that("gauss_laguerre_nodes weights match mu_0 = Gamma(alpha + 1)", {
  test_params <- c(0, 0.5, 1.5, 3.0)

  for (alpha_param in test_params) {
    quad <- gauss_laguerre_nodes(80, alpha_param)
    weight_sum <- sum(quad$weights)
    mu0_expected <- exp(lgamma(alpha_param + 1))

    expect_equal(weight_sum, mu0_expected, tolerance = 1e-10,
                 info = sprintf("alpha_param = %.1f", alpha_param))
  }
})


test_that("build_gamma_quadrature returns normalized weights", {
  test_cases <- list(
    list(a = 1.0, b = 1.0),
    list(a = 2.5, b = 1.5),
    list(a = 0.5, b = 2.0),
    list(a = 5.0, b = 0.5)
  )

  for (tc in test_cases) {
    quad <- build_gamma_quadrature(tc$a, tc$b, M = 80)

    expect_equal(sum(quad$weights_normalized), 1, tolerance = 1e-15,
                 info = sprintf("a=%.1f, b=%.1f", tc$a, tc$b))
  }
})


test_that("Gamma moments via quadrature match closed-form", {
  test_cases <- list(
    list(a = 1.0, b = 1.0),
    list(a = 2.5, b = 1.5),
    list(a = 0.5, b = 2.0),
    list(a = 5.0, b = 0.5)
  )

  for (tc in test_cases) {
    a <- tc$a
    b <- tc$b

    # E[alpha] = a/b
    E_alpha <- integrate_gamma(identity, a, b, M = 80)
    expect_equal(E_alpha, a/b, tolerance = 1e-10,
                 info = sprintf("E[alpha]: a=%.1f, b=%.1f", a, b))

    # E[alpha^2] = a(a+1)/b^2
    E_alpha_sq <- integrate_gamma(function(x) x^2, a, b, M = 80)
    expect_equal(E_alpha_sq, a * (a + 1) / b^2, tolerance = 1e-10,
                 info = sprintf("E[alpha^2]: a=%.1f, b=%.1f", a, b))

    # Var(alpha) = a/b^2
    Var_alpha <- E_alpha_sq - E_alpha^2
    expect_equal(Var_alpha, a/b^2, tolerance = 1e-10,
                 info = sprintf("Var(alpha): a=%.1f, b=%.1f", a, b))
  }
})


test_that("Quadrature converges with increasing M", {
  a <- 2.5
  b <- 1.5
  true_mean <- a / b

  M_values <- c(10, 20, 50, 100)
  errors <- sapply(M_values, function(M) {
    abs(integrate_gamma(identity, a, b, M) - true_mean)
  })

  # All errors should be very small (near machine precision)
  expect_true(all(errors < 1e-12),
              info = "All errors should be near machine precision")

  # Final error should not be worse than first (allowing for machine precision noise)
  expect_true(errors[length(errors)] <= errors[1] + 1e-14,
              info = "Error at high M should not be worse than at low M")
})


test_that("Higher moments are accurate", {
  a <- 3.0
  b <- 1.5
  M <- 100

  # E[alpha^n] = a(a+1)...(a+n-1) / b^n
  for (n in 1:4) {
    true_moment <- prod(a + 0:(n-1)) / b^n
    quad_moment <- integrate_gamma(function(x) x^n, a, b, M)

    expect_equal(quad_moment, true_moment, tolerance = 1e-9,
                 info = sprintf("E[alpha^%d]", n))
  }
})


test_that("verify_quadrature returns TRUE for standard cases", {
  expect_true(verify_quadrature(1.0, 1.0, M = 80, verbose = FALSE))
  expect_true(verify_quadrature(2.5, 1.5, M = 80, verbose = FALSE))
  expect_true(verify_quadrature(5.0, 0.5, M = 80, verbose = FALSE))
})


test_that("verify_quadrature handles challenging parameters", {
  expect_true(verify_quadrature(0.5, 2.0, M = 80, tol = 1e-10, verbose = FALSE))
  expect_true(verify_quadrature(0.1, 1.0, M = 100, tol = 1e-8, verbose = FALSE))
  expect_true(verify_quadrature(10.0, 0.1, M = 100, tol = 1e-8, verbose = FALSE))
})


test_that("Input validation works for gauss_laguerre_nodes", {
  expect_error(gauss_laguerre_nodes(0), "positive integer")
  expect_error(gauss_laguerre_nodes(-5), "positive integer")
  expect_error(gauss_laguerre_nodes(1.5), "positive integer")
  expect_error(gauss_laguerre_nodes(40, alpha_param = -1.5), "must be.*> -1")
  expect_error(gauss_laguerre_nodes(40, alpha_param = -1), "must be.*> -1")
})


test_that("Input validation works for build_gamma_quadrature", {
  expect_error(build_gamma_quadrature(0, 1), "must be.*positive")
  expect_error(build_gamma_quadrature(-1, 1), "must be.*positive")
  expect_error(build_gamma_quadrature(1, 0), "must be.*positive")
  expect_error(build_gamma_quadrature(1, -1), "must be.*positive")
})


test_that("summary_quadrature returns expected structure", {
  summ <- summary_quadrature(2.5, 1.5, M = 80)

  expect_type(summ, "list")
  expect_equal(summ$n_nodes, 80)
  expect_equal(summ$gamma_mean, 2.5/1.5, tolerance = 1e-15)
  expect_equal(summ$gamma_sd, sqrt(2.5)/1.5, tolerance = 1e-15)
})


test_that("convergence_quadrature provides diagnostic information", {
  result <- convergence_quadrature(identity, 2.5, 1.5,
                                   M_values = c(10, 20, 50),
                                   true_value = 2.5/1.5)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3)
  expect_true("error" %in% names(result))

  # Errors should be very small at high M
  expect_lt(result$error[3], 1e-12)
})


test_that("Quadrature handles smooth non-polynomial functions", {
  # E[log(alpha)] for Gamma(a, b) = digamma(a) - log(b)
  # Use larger shape parameter for better accuracy away from 0
  a <- 5.0
  b <- 1.0

  E_log_alpha_true <- digamma(a) - log(b)
  E_log_alpha_quad <- integrate_gamma(log, a, b, M = 100)

  # Expect reasonable accuracy for non-polynomial
  expect_equal(E_log_alpha_quad, E_log_alpha_true, tolerance = 1e-4)
})


test_that("Quadrature handles reciprocal function with large shape", {
  # E[1/alpha] for Gamma(a, b) = b/(a-1) when a > 1
  # Use large shape parameter to avoid singularity at 0
  a <- 5.0
  b <- 1.0

  E_recip_true <- b / (a - 1)
  E_recip_quad <- integrate_gamma(function(x) 1/x, a, b, M = 100)

  # Expect good accuracy with large shape
  expect_equal(E_recip_quad, E_recip_true, tolerance = 1e-6)
})


# =============================================================================
# Golden Data Tests
# =============================================================================

test_that("Quadrature matches golden test data", {
  golden_data <- data.frame(
    a = c(1.0, 2.5, 0.5, 5.0),
    b = c(1.0, 1.5, 2.0, 0.5),
    E_alpha = c(1.0, 1.6666666666666667, 0.25, 10.0),
    Var_alpha = c(1.0, 1.1111111111111112, 0.125, 20.0),
    E_alpha_sq = c(2.0, 3.888888888888889, 0.1875, 120.0)
  )

  for (i in seq_len(nrow(golden_data))) {
    a <- golden_data$a[i]
    b <- golden_data$b[i]

    E_alpha_quad <- integrate_gamma(identity, a, b, M = 80)
    E_alpha_sq_quad <- integrate_gamma(function(x) x^2, a, b, M = 80)
    Var_alpha_quad <- E_alpha_sq_quad - E_alpha_quad^2

    expect_equal(E_alpha_quad, golden_data$E_alpha[i], tolerance = 1e-10,
                 info = sprintf("E[alpha] golden test: a=%.1f, b=%.1f", a, b))
    expect_equal(Var_alpha_quad, golden_data$Var_alpha[i], tolerance = 1e-10,
                 info = sprintf("Var(alpha) golden test: a=%.1f, b=%.1f", a, b))
    expect_equal(E_alpha_sq_quad, golden_data$E_alpha_sq[i], tolerance = 1e-10,
                 info = sprintf("E[alpha^2] golden test: a=%.1f, b=%.1f", a, b))
  }
})
