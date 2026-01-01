# =============================================================================
# testthat Tests: Module 12 A2-KL
# =============================================================================

test_that("KL divergence is non-negative", {
  J <- 50
  target <- discretize_chisq(J, df = 10, scale = 1)

  for (a in c(1, 2, 5)) {
    for (b in c(0.5, 1, 2)) {
      kl <- kl_divergence_K(target, a = a, b = b, J = J)
      expect_gte(kl, 0)
    }
  }
})


test_that("KL divergence is zero for identical distributions", {
  J <- 30
  logS <- compute_log_stirling(J)

  # Generate induced PMF for known (a, b)
  a_true <- 2.0
  b_true <- 1.5
  induced <- pmf_K_marginal(J, a_true, b_true, logS)[-1L]
  induced <- induced / sum(induced)

  # KL should be essentially zero
  kl <- kl_divergence_K(induced, a = a_true, b = b_true, J = J)
  expect_lt(kl, 1e-10)
})


test_that("A2-KL converges with chisq method", {
  J <- 50
  target <- list(mu_K = 5, var_K = 8)

  fit <- DPprior_a2_kl(J, target, method = "chisq")

  expect_true(fit$converged || fit$status == "fallback_used")
  expect_lt(fit$fit$kl, 0.1)
  expect_equal(fit$method, "A2-KL")
  expect_equal(fit$target$type, "chisq")
})


test_that("A2-KL converges with pmf method", {
  J <- 50

  # Binomial-shaped target
  target_pmf <- dbinom(1:J, size = J, prob = 0.1)
  target_pmf <- target_pmf / sum(target_pmf)

  fit <- DPprior_a2_kl(J, target_pmf, method = "pmf")

  expect_true(fit$converged || fit$status == "fallback_used")
  expect_lt(fit$fit$kl, 0.5)
  expect_equal(fit$method, "A2-KL")
  expect_equal(fit$target$type, "pmf")
})


test_that("A2-KL stores chi-square parameters for chisq method", {
  J <- 40
  mu_K <- 6
  var_K <- 10

  fit <- DPprior_a2_kl(J, list(mu_K = mu_K, var_K = var_K), method = "chisq")

  # Check that df and scale are stored

  expect_true(!is.null(fit$target$df))
  expect_true(!is.null(fit$target$scale))
  expect_true(!is.null(fit$target$mu_K_discrete))
  expect_true(!is.null(fit$target$var_K_discrete))

  # Verify df/scale match the expected formulas
  expected_df <- 2 * mu_K^2 / var_K
  expected_scale <- var_K / (2 * mu_K)
  expect_equal(fit$target$df, expected_df)
  expect_equal(fit$target$scale, expected_scale)
})


test_that("A2-KL provides comprehensive diagnostics", {
  J <- 50
  target <- list(mu_K = 5, var_K = 8)

  fit <- DPprior_a2_kl(J, target, method = "chisq")

  # Check diagnostics structure
  expect_true(!is.null(fit$diagnostics$init$a0))
  expect_true(!is.null(fit$diagnostics$init$b0))
  expect_true(!is.null(fit$diagnostics$init$init_method))
  expect_true(!is.null(fit$diagnostics$optim$method))
  expect_true(!is.null(fit$diagnostics$fallback_used))
  expect_true(!is.null(fit$diagnostics$kl_init))
  expect_true(!is.null(fit$diagnostics$kl_final))
})


test_that("A2-KL provides trace information", {
  J <- 50
  target <- list(mu_K = 5, var_K = 8)

  fit <- DPprior_a2_kl(J, target, method = "chisq")

  # Check trace structure
  expect_true(!is.null(fit$trace))
  expect_s3_class(fit$trace, "data.frame")
  expect_true(all(c("eval", "a", "b", "kl") %in% names(fit$trace)))
  expect_gt(nrow(fit$trace), 0)
})


test_that("A2-KL handles edge cases gracefully", {
  # Small J
  fit_small <- DPprior_a2_kl(J = 10, target = list(mu_K = 3, var_K = 2),
                             method = "chisq")
  expect_true(is.finite(fit_small$fit$kl))

  # Target with low variance (may trigger A1 projection warning)
  fit_low_var <- suppressWarnings(
    DPprior_a2_kl(J = 50, target = list(mu_K = 10, var_K = 3),
                  method = "chisq")
  )
  expect_true(is.finite(fit_low_var$fit$kl))
})


test_that("A2-KL fallback mechanism works", {
  # This test checks that fallback doesn't crash
  # and returns finite results even for challenging targets
  # Note: suppressWarnings because challenging targets may trigger
  # projection warnings in A1 initialization

  J <- 30

  # Use a challenging but valid target
  fit <- suppressWarnings(
    DPprior_a2_kl(J = J, target = list(mu_K = 3, var_K = 1.5),
                  method = "chisq")
  )

  expect_true(is.finite(fit$a))
  expect_true(is.finite(fit$b))
  expect_true(is.finite(fit$fit$kl))
  expect_true(is.logical(fit$diagnostics$fallback_used))
})


test_that("discretize_chisq produces valid PMF", {
  J <- 50

  for (df in c(5, 10, 20)) {
    for (scale in c(0.5, 1, 2)) {
      pmf <- discretize_chisq(J, df = df, scale = scale)

      expect_length(pmf, J)
      expect_true(all(pmf >= 0))
      expect_equal(sum(pmf), 1, tolerance = 1e-10)
    }
  }
})


test_that("kl_divergence_pmf has correct properties", {
  # Same distribution: KL = 0
  p <- c(0.2, 0.5, 0.3)
  expect_equal(kl_divergence_pmf(p, p), 0, tolerance = 1e-10)

  # Non-negative
  q <- c(0.3, 0.4, 0.3)
  expect_gte(kl_divergence_pmf(p, q), 0)

  # Asymmetric
  kl_pq <- kl_divergence_pmf(p, q)
  kl_qp <- kl_divergence_pmf(q, p)
  expect_false(isTRUE(all.equal(kl_pq, kl_qp)))
})


test_that("construct_target_pmf works with PMF input", {
  J <- 50

  # Uniform PMF
  target <- rep(1, J)
  result <- construct_target_pmf(J, target)

  expect_length(result$pmf, J)
  expect_equal(sum(result$pmf), 1, tolerance = 1e-10)
  expect_true(!is.null(result$mu_K))
  expect_true(!is.null(result$var_K))

  # Verify moments
  k_vals <- 1:J
  expect_equal(result$mu_K, sum(k_vals * result$pmf), tolerance = 1e-10)
})


test_that("construct_target_pmf works with moment input", {
  J <- 50
  mu_K <- 5
  var_K <- 8

  result <- construct_target_pmf(J, list(mu_K = mu_K, var_K = var_K))

  expect_length(result$pmf, J)
  expect_equal(sum(result$pmf), 1, tolerance = 1e-10)
  expect_equal(result$mu_K, mu_K)
  expect_equal(result$var_K, var_K)

  # Check chi-square parameters
  expect_true(!is.null(result$df))
  expect_true(!is.null(result$scale))
  expect_equal(result$df, 2 * mu_K^2 / var_K)
  expect_equal(result$scale, var_K / (2 * mu_K))

  # Check discretized moments
  expect_true(!is.null(result$mu_K_discrete))
  expect_true(!is.null(result$var_K_discrete))
})
