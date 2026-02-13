# =============================================================================
# Test Suite for Module 07: Score-Based Jacobian
# =============================================================================
#
# This file contains unit tests for the score-based Jacobian computation
# functions in R/07_jacobian.R.
#
# Author: JoonHo Lee
# Date: December 2025
# Part of: DPprior R Package
# =============================================================================

# =============================================================================
# Score Function Tests
# =============================================================================

test_that("score_a computes correctly", {
  # s_a(α) = log(b) - ψ(a) + log(α)
  a <- 2.0
  b <- 1.0
  alpha <- 2.0

  # Manual calculation
  expected <- log(b) - digamma(a) + log(alpha)
  expect_equal(score_a(alpha, a, b), expected, tolerance = 1e-12)

  # Vectorized
  alphas <- c(0.5, 1.0, 2.0, 5.0)
  expected_vec <- log(b) - digamma(a) + log(alphas)
  expect_equal(score_a(alphas, a, b), expected_vec, tolerance = 1e-12)
})


test_that("score_b computes correctly", {
  # s_b(α) = a/b - α
  a <- 2.0
  b <- 1.0
  alpha <- 2.0

  # Manual calculation
  expected <- a / b - alpha
  expect_equal(score_b(alpha, a, b), expected, tolerance = 1e-12)

  # Vectorized
  alphas <- c(0.5, 1.0, 2.0, 5.0)
  expected_vec <- a / b - alphas
  expect_equal(score_b(alphas, a, b), expected_vec, tolerance = 1e-12)
})


test_that("score functions validate input", {
  expect_error(score_a(-1, 2, 1))
  expect_error(score_a(1, -2, 1))
  expect_error(score_a(1, 2, -1))
  expect_error(score_b(-1, 2, 1))
  expect_error(score_b(1, -2, 1))
  expect_error(score_b(1, 2, -1))
})


test_that("score_b expectation is approximately zero", {
  # E[s_b] = E[a/b - α] = a/b - E[α] = a/b - a/b = 0
  # This should be very close to zero with quadrature
  test_cases <- list(
    c(a = 1.0, b = 0.5),
    c(a = 2.0, b = 1.0),
    c(a = 1.5, b = 2.0),
    c(a = 3.0, b = 1.5)
  )

  for (tc in test_cases) {
    a <- tc["a"]
    b <- tc["b"]
    quad <- build_gamma_quadrature(a, b, 80)
    E_sb <- sum(quad$weights_normalized * score_b(quad$alpha_nodes, a, b))

    # Should be very close to zero (s_b is linear in α)
    expect_equal(E_sb, 0, tolerance = 1e-10)
  }
})


test_that("Score functions have zero expectation (adaptive check)", {
  # Use R's integrate() for independent verification
  a <- 2.0
  b <- 1.0

  sa <- function(x) score_a(x, a, b) * dgamma(x, shape = a, rate = b)
  sb <- function(x) score_b(x, a, b) * dgamma(x, shape = a, rate = b)

  E_sa <- integrate(sa, lower = 0, upper = Inf, subdivisions = 2000)$value
  E_sb <- integrate(sb, lower = 0, upper = Inf, subdivisions = 2000)$value

  # E_sa has larger error due to log(alpha) term near zero
  # E_sb is linear in alpha, so converges faster
  expect_equal(E_sa, 0, tolerance = 1e-6)
  expect_equal(E_sb, 0, tolerance = 1e-8)
})


# =============================================================================
# Jacobian Computation Tests
# =============================================================================

test_that("moments_with_jacobian returns correct structure", {
  result <- moments_with_jacobian(J = 50, a = 2.0, b = 1.0)

  expect_true(is.list(result))
  expect_true("mean" %in% names(result))
  expect_true("var" %in% names(result))
  expect_true("jacobian" %in% names(result))

  expect_true(is.numeric(result$mean))
  expect_true(is.numeric(result$var))
  expect_true(is.matrix(result$jacobian))

  expect_equal(nrow(result$jacobian), 2L)
  expect_equal(ncol(result$jacobian), 2L)
})


test_that("moments_with_jacobian matches exact_K_moments", {
  test_cases <- list(
    list(J = 30, a = 1.5, b = 0.5),
    list(J = 50, a = 2.0, b = 1.0),
    list(J = 100, a = 2.0, b = 1.0),
    list(J = 50, a = 3.0, b = 1.5)
  )

  for (tc in test_cases) {
    result_jac <- moments_with_jacobian(tc$J, tc$a, tc$b)
    result_mom <- exact_K_moments(tc$J, tc$a, tc$b)

    expect_equal(result_jac$mean, result_mom$mean, tolerance = 1e-10)
    expect_equal(result_jac$var, result_mom$var, tolerance = 1e-10)
  }
})


test_that("Jacobian (mean row) matches finite differences", {
  # Use higher M for derivatives involving log(alpha)
  J <- 50
  a <- 2.0
  b <- 1.0
  eps <- 1e-6
  M_use <- 200

  result <- moments_with_jacobian(J, a, b, M = M_use)

  m_pp <- exact_K_moments(J, a + eps, b, M = M_use)$mean
  m_pm <- exact_K_moments(J, a - eps, b, M = M_use)$mean
  dM1_da_num <- (m_pp - m_pm) / (2 * eps)

  m_mp <- exact_K_moments(J, a, b + eps, M = M_use)$mean
  m_mm <- exact_K_moments(J, a, b - eps, M = M_use)$mean
  dM1_db_num <- (m_mp - m_mm) / (2 * eps)

  expect_equal(result$jacobian[1, 1], dM1_da_num, tolerance = 1e-3)
  expect_equal(result$jacobian[1, 2], dM1_db_num, tolerance = 1e-6)
})


test_that("Jacobian matches finite differences (full matrix)", {
  eps <- 1e-6
  tol <- 0.01  # 1% relative tolerance

  test_cases <- list(
    list(J = 50, a = 2.0, b = 1.0),
    list(J = 50, a = 3.0, b = 1.5),
    list(J = 100, a = 2.0, b = 1.0)
  )

  for (tc in test_cases) {
    J <- tc$J
    a <- tc$a
    b <- tc$b

    # Use higher M for more accurate comparison
    M <- 200L
    result <- moments_with_jacobian(J, a, b, M = M)
    J_analytic <- result$jacobian

    # Numerical Jacobian with even higher M
    M_fd <- 250L
    m_a_plus <- exact_K_moments(J, a + eps, b, M_fd)
    m_a_minus <- exact_K_moments(J, a - eps, b, M_fd)
    m_b_plus <- exact_K_moments(J, a, b + eps, M_fd)
    m_b_minus <- exact_K_moments(J, a, b - eps, M_fd)

    dM1_da_num <- (m_a_plus$mean - m_a_minus$mean) / (2 * eps)
    dM1_db_num <- (m_b_plus$mean - m_b_minus$mean) / (2 * eps)
    dV_da_num <- (m_a_plus$var - m_a_minus$var) / (2 * eps)
    dV_db_num <- (m_b_plus$var - m_b_minus$var) / (2 * eps)

    J_numeric <- matrix(c(dM1_da_num, dM1_db_num, dV_da_num, dV_db_num),
                        nrow = 2, byrow = TRUE)

    # Compare with relative tolerance
    rel_err <- abs(J_analytic - J_numeric) / (abs(J_numeric) + 1e-10)
    max_rel_err <- max(rel_err)

    # Use expect_true with label (not info - deprecated in testthat 3.0)
    expect_true(
      max_rel_err < tol,
      label = sprintf("J=%d, a=%.1f, b=%.1f: max_rel_err=%.4f < %.4f",
                      J, a, b, max_rel_err, tol)
    )
  }
})


test_that("Jacobian is non-singular for typical parameters", {
  test_cases <- list(
    list(J = 30, a = 1.5, b = 0.5),
    list(J = 50, a = 2.0, b = 1.0),
    list(J = 100, a = 2.0, b = 1.0),
    list(J = 50, a = 3.0, b = 1.5)
  )

  for (tc in test_cases) {
    result <- moments_with_jacobian(tc$J, tc$a, tc$b)
    det_val <- det(result$jacobian)
    cond_val <- kappa(result$jacobian)

    # Use expect_true with label (not expect_lt with info)
    expect_true(
      abs(det_val) > 1e-10,
      label = sprintf("J=%d, a=%.1f, b=%.1f: det=%.2e > 1e-10",
                      tc$J, tc$a, tc$b, det_val)
    )

    expect_true(
      cond_val < 1e10,
      label = sprintf("J=%d, a=%.1f, b=%.1f: cond=%.2e < 1e10",
                      tc$J, tc$a, tc$b, cond_val)
    )
  }
})


test_that("Jacobian allows Newton convergence", {
  # Simple Newton iteration to verify the Jacobian works for moment matching
  J <- 50
  mu_target <- 5.0
  var_target <- 8.0

  # Start from rough initial guess
  a <- 2.0
  b <- 1.0

  for (iter in 1:15) {
    result <- moments_with_jacobian(J, a, b, M = 120L)
    F_vec <- c(result$mean - mu_target, result$var - var_target)
    norm_F <- sqrt(sum(F_vec^2))

    if (norm_F < 1e-8) break

    delta <- solve(result$jacobian, -F_vec)
    a <- a + delta[1]
    b <- b + delta[2]
  }

  # Should converge to target moments
  final <- exact_K_moments(J, a, b)
  expect_equal(final$mean, mu_target, tolerance = 1e-6)
  expect_equal(final$var, var_target, tolerance = 1e-6)
})


test_that("moments_with_jacobian handles different M values", {
  J <- 50
  a <- 2.0
  b <- 1.0

  result_40 <- moments_with_jacobian(J, a, b, M = 40)
  result_80 <- moments_with_jacobian(J, a, b, M = 80)
  result_120 <- moments_with_jacobian(J, a, b, M = 120)

  # Moments should converge with increasing M (tight tolerance)
  expect_equal(result_40$mean, result_80$mean, tolerance = 1e-6)
  expect_equal(result_80$mean, result_120$mean, tolerance = 1e-8)

  # =====================================================================
  # IMPORTANT: Jacobians converge MORE SLOWLY than moments!
  # =====================================================================
  # The score function s_a contains log(alpha) which is not smooth at alpha->0
  # This causes slower quadrature convergence for score-weighted integrands.
  # Use LOOSER tolerance (1e-2 instead of 1e-4) - this is EXPECTED BEHAVIOR
  # =====================================================================

  expect_equal(result_40$jacobian, result_80$jacobian, tolerance = 0.01)
  expect_equal(result_80$jacobian, result_120$jacobian, tolerance = 0.005)
})


test_that("moments_with_jacobian validates input", {
  expect_error(moments_with_jacobian(J = -1, a = 2, b = 1))
  expect_error(moments_with_jacobian(J = 50, a = -2, b = 1))
  expect_error(moments_with_jacobian(J = 50, a = 2, b = -1))
  expect_error(moments_with_jacobian(J = 50, a = 2, b = 1, M = 5))
})


# =============================================================================
# Golden Test Data (from adaptive integration reference)
# =============================================================================

test_that("Jacobian matches golden test values", {
  # Golden test data from Python adaptive integration (scipy.integrate.quad)
  # This is the GROUND TRUTH independent of quadrature
  golden <- list(
    list(J = 50, a = 2.0, b = 1.0,
         mean = 6.6396928911, var = 12.9545022869,
         dM1_da = 2.2455202979, dM1_db = -4.1355851664,
         dV_da = 2.9446327031, dV_db = -13.0382284995),

    list(J = 100, a = 2.0, b = 1.0,
         mean = 7.9785512638, var = 20.4068433344,
         dM1_da = 2.8965938003, dM1_db = -5.4201373588,
         dV_da = 5.5851245747, dV_db = -23.8421067085),

    list(J = 100, a = 3.0, b = 1.5,
         mean = 8.1073904979, var = 15.9908796535,
         dM1_da = 1.9408087472, dM1_db = -3.7036385429,
         dV_da = 2.8881654955, dV_db = -11.9313891758)
  )

  for (g in golden) {
    # Use high M for comparison with adaptive reference
    result <- moments_with_jacobian(g$J, g$a, g$b, M = 200)

    expect_equal(result$mean, g$mean, tolerance = 1e-6)
    expect_equal(result$var, g$var, tolerance = 1e-6)

    # Jacobian elements - looser tolerance due to quadrature approximation
    # of log(alpha) weighted integrands
    expect_equal(result$jacobian[1, 1], g$dM1_da, tolerance = 1e-3)
    expect_equal(result$jacobian[1, 2], g$dM1_db, tolerance = 1e-3)
    expect_equal(result$jacobian[2, 1], g$dV_da, tolerance = 1e-3)
    expect_equal(result$jacobian[2, 2], g$dV_db, tolerance = 1e-3)
  }
})


# =============================================================================
# Verification Function Tests
# =============================================================================

test_that("verify_jacobian runs without error", {
  result <- verify_jacobian(J = 50, a = 2.0, b = 1.0, verbose = FALSE)

  expect_true(is.list(result))
  expect_true("analytic" %in% names(result))
  expect_true("numeric" %in% names(result))
  expect_true("pass" %in% names(result))
})


test_that("verify_score_expectation runs without error", {
  result <- verify_score_expectation(a = 2.0, b = 1.0, verbose = FALSE)

  expect_true(is.list(result))
  expect_true("E_score_a" %in% names(result))
  expect_true("E_score_b" %in% names(result))
})


test_that("verify_jacobian_all completes", {
  result <- verify_jacobian_all(verbose = FALSE)
  expect_true(is.logical(result))
})
