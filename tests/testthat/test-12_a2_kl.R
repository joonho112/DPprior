# =============================================================================
# testthat Tests: Module 12 A2-KL
# =============================================================================

.a2_kl_expect_valid_result <- function(fit, status = NULL) {
  expect_s3_class(fit, "DPprior_fit")
  expect_s3_class(fit, "dpprior_result")
  expect_identical(fit$schema, list(name = "dpprior.result", version = 1L))
  expect_identical(fit$object_type, "fit")
  expect_identical(fit$mode, "a2_kl")
  expect_identical(fit$method, "A2-KL")
  if (!is.null(status)) expect_identical(fit$status, status)
  report <- .dpprior_validate_result_v1(fit, collect = TRUE)
  details <- if (length(report$errors)) {
    paste(vapply(report$errors, conditionMessage, character(1)), collapse = "\n")
  } else {
    ""
  }
  expect_true(report$valid, info = details)
  expect_length(report$errors, 0L)
  invisible(fit)
}


.a2_kl_expect_schema_rejection <- function(object) {
  condition <- tryCatch(
    {
      .dpprior_validate_result_v1(object)
      NULL
    },
    error = identity
  )
  expect_s3_class(condition, "dpprior_schema_error")
  expect_true(
    is.character(condition$code) && length(condition$code) == 1L &&
      !is.na(condition$code) && nzchar(condition$code)
  )
  expect_true(
    is.character(condition$path) && length(condition$path) == 1L &&
      !is.na(condition$path) && nzchar(condition$path)
  )
  invisible(condition)
}

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

  expect_true(fit$status %in% c("converged", "boundary"))
  expect_true(fit$usable)
  expect_true(fit$verified)
  expect_lt(fit$fit$kl, 0.1)
  expect_equal(fit$method, "A2-KL")
  expect_identical(fit$schema$name, "dpprior.result")
  expect_identical(fit$schema$version, 1L)
  expect_identical(fit$mode, "a2_kl")
  expect_identical(fit$target$K$kind, "moments")
  expect_identical(
    fit$target$K$request,
    list(J = 50L, mu_K = 5, var_K = 8)
  )
})


test_that("A2-KL converges with pmf method", {
  J <- 50

  # Binomial-shaped target
  target_pmf <- dbinom(1:J, size = J, prob = 0.1)
  target_pmf <- target_pmf / sum(target_pmf)

  fit <- DPprior_a2_kl(J, target_pmf, method = "pmf")

  expect_true(fit$status %in% c("converged", "boundary"))
  expect_true(fit$usable)
  expect_true(fit$verified)
  expect_lt(fit$fit$kl, 0.5)
  expect_equal(fit$method, "A2-KL")
  expect_identical(fit$target$K$kind, "pmf")
  expect_identical(fit$target$K$request$pmf, target_pmf)
  expect_identical(fit$target$K$normalized$pmf, target_pmf)
  expect_identical(fit$target$K$used$pmf, target_pmf)
})


test_that("A2-KL stores chi-square parameters for chisq method", {
  J <- 40
  mu_K <- 6
  var_K <- 10

  fit <- DPprior_a2_kl(J, list(mu_K = mu_K, var_K = var_K), method = "chisq")

  # The requested moments remain authoritative. The conditioned chi-square
  # PMF is objective evidence, not a promoted user request.
  objective <- fit$target$K$derivation$request_to_normalized$
    evidence$A2_KL_objective
  expect_identical(fit$target$K$kind, "moments")
  expect_identical(
    fit$target$K$request,
    list(J = 40L, mu_K = mu_K, var_K = var_K)
  )
  expect_null(fit$target$K$pmf)
  expect_identical(objective$method, "chisq")
  expect_identical(objective$source, "A2_KL_backend_target")
  expect_length(objective$pmf, J)
  expect_equal(sum(objective$pmf), 1, tolerance = 1e-12)
  expect_true(is.finite(objective$mu_K_discrete))
  expect_true(is.finite(objective$var_K_discrete))

  # Verify df/scale match the expected formulas
  expected_df <- 2 * mu_K^2 / var_K
  expected_scale <- var_K / (2 * mu_K)
  expect_equal(objective$df, expected_df)
  expect_equal(objective$scale, expected_scale)
})


test_that("A2-KL canonicalization preserves exact selected-order identity", {
  chisq <- DPprior_a2_kl(
    50L, list(mu_K = 5, var_K = 8), method = "chisq"
  )
  target_pmf <- dbinom(1:50, size = 50, prob = 0.1)
  target_pmf <- target_pmf / sum(target_pmf)
  explicit <- DPprior_a2_kl(50L, target_pmf, method = "pmf")

  .a2_kl_expect_valid_result(chisq, "converged")
  .a2_kl_expect_valid_result(explicit, "converged")
  expect_identical(
    sprintf("%.17g", c(
      a = chisq$a, b = chisq$b, mean = chisq$achieved$K$mean,
      variance = chisq$achieved$K$variance,
      kl = chisq$residuals$distribution$kl,
      l1 = chisq$residuals$distribution$l1
    )),
    c(
      "2.2363025420148892", "1.7626540125938372",
      "5.0199840940404021", "7.6403352607781523",
      "0.0050296580334789698", "0.069116930667064411"
    )
  )
  expect_identical(
    sprintf("%.17g", c(
      a = explicit$a, b = explicit$b, mean = explicit$achieved$K$mean,
      variance = explicit$achieved$K$variance,
      kl = explicit$residuals$distribution$kl,
      l1 = explicit$residuals$distribution$l1
    )),
    c(
      "8.8137350305399131", "7.2918393773854957",
      "5.0237318007775178", "4.3913618182239045",
      "0.0018771831867521467", "0.037877755676274862"
    )
  )

  for (fit in list(chisq, explicit)) {
    legacy <- fit$compatibility$views$a2_kl_v0
    expect_identical(fit$a, legacy$a)
    expect_identical(fit$b, legacy$b)
    expect_identical(fit$achieved$K$mean, legacy$achieved$mu_K)
    expect_identical(fit$achieved$K$variance, legacy$achieved$var_K)
    expect_identical(fit$residuals$distribution$kl, legacy$achieved$kl)
    expect_identical(fit$residuals$distribution$l1, legacy$achieved$l1)
  }
})


test_that("A2-KL target authority distinguishes PMF input from chisq objective", {
  target_pmf <- rep(1 / 12, 12L)
  explicit <- DPprior_a2_kl(12L, target_pmf, method = "pmf", M = 40L)
  chisq <- DPprior_a2_kl(
    12L, list(mu_K = 4, var_K = 5), method = "chisq", M = 40L
  )

  expect_identical(explicit$target$K$request, list(
    J = 12L, pmf = target_pmf
  ))
  expect_identical(
    explicit$target$K$derivation$request_to_normalized$rule,
    "validate_strict_pmf"
  )
  expect_identical(
    explicit$target$K$derivation$request_to_normalized$before,
    explicit$target$K$request
  )
  expect_identical(
    explicit$target$K$derivation$request_to_normalized$after,
    explicit$target$K$normalized
  )
  expect_null(explicit$target$K$derivation$normalized_to_used)
  expect_null(explicit$target$K$interval)
  expect_null(explicit$target$K$family)

  expect_identical(
    chisq$target$K$request,
    list(J = 12L, mu_K = 4, var_K = 5)
  )
  expect_identical(chisq$target$K$normalized, chisq$target$K$used)
  expect_null(chisq$target$K$pmf)
  expect_null(chisq$target$K$interval)
  expect_null(chisq$target$K$family)
  objective <- chisq$target$K$derivation$request_to_normalized$
    evidence$A2_KL_objective
  expect_identical(objective$support, c(lower = 1L, upper = 12L))
  expect_length(objective$pmf, 12L)
  expect_identical(
    objective$pmf,
    chisq$compatibility$views$a2_kl_v0$target$pmf
  )
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


test_that("A2-KL records heterogeneous attempts and fresh KL objectives", {
  fit <- DPprior_a2_kl(
    30L, list(mu_K = 5, var_K = 8), method = "chisq"
  )
  attempts <- fit$computation$attempts
  candidates <- fit$computation$candidate_evaluations

  expect_identical(
    vapply(attempts, `[[`, character(1), "method"),
    c("A2-MN", "L-BFGS-B")
  )
  expect_identical(
    vapply(attempts, `[[`, character(1), "stage"),
    c("initialization", "primary")
  )
  expect_identical(
    vapply(candidates, `[[`, character(1), "objective_kind"),
    c("kl", "kl")
  )
  expect_null(candidates[[1L]]$recorded_objective)
  expect_false(candidates[[1L]]$recorded_objective_available)
  expect_identical(candidates[[1L]]$generator, "initialization")
  expect_true(is.finite(candidates[[1L]]$fresh_objective))
  expect_true(candidates[[2L]]$recorded_objective_available)
  expect_identical(candidates[[2L]]$recorded_objective_kind, "kl")
  expect_identical(candidates[[2L]]$selection_objective_kind, "kl")
  expect_identical(
    fit$computation$selected_candidate_id,
    candidates[[which(vapply(candidates, `[[`, logical(1), "selected"))]]$id
  )
})


test_that("A2-KL records A1 and heuristic initializer fallback routes", {
  local_mocked_bindings(
    DPprior_a2_newton = function(...) stop("injected A2-MN failure"),
    .package = "DPprior"
  )
  target <- dbinom(1:20, size = 20, prob = 0.2)
  target <- target / sum(target)
  fit <- DPprior_a2_kl(20L, target, method = "pmf", M = 40L)

  .a2_kl_expect_valid_result(fit)
  expect_identical(
    vapply(fit$computation$attempts[1:2], `[[`, character(1), "method"),
    c("A2-MN", "A1")
  )
  expect_identical(
    fit$computation$attempts[[1L]]$reason_code,
    "optimizer_error"
  )
  expect_identical(
    fit$computation$resources$initialization_method, "A1"
  )
})


test_that("A2-KL retains deterministic heuristic initialization evidence", {
  local_mocked_bindings(
    DPprior_a2_newton = function(...) stop("injected A2-MN failure"),
    DPprior_a1 = function(...) stop("injected A1 failure"),
    .package = "DPprior"
  )
  target <- dbinom(1:20, size = 20, prob = 0.2)
  target <- target / sum(target)
  fit <- DPprior_a2_kl(20L, target, method = "pmf", M = 40L)

  .a2_kl_expect_valid_result(fit)
  expect_identical(
    vapply(fit$computation$attempts[1:3], `[[`, character(1), "method"),
    c("A2-MN", "A1", "heuristic")
  )
  expect_identical(
    fit$computation$resources$initialization_method, "heuristic"
  )
  heuristic <- fit$computation$candidate_evaluations[[which(vapply(
    fit$computation$candidate_evaluations,
    function(candidate) identical(candidate$method, "heuristic"),
    logical(1)
  ))]]
  expect_null(heuristic$recorded_objective)
  expect_true(is.finite(heuristic$fresh_objective))
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

  expect_error(
    DPprior_a2_kl(J = 10, target = list(mu_K = 9, var_K = 9),
                  method = "chisq"),
    "var_K = 9.*maximum possible variance 8.*K in \\{1,...,10\\}.*mu_K = 9"
  )
})


test_that("A2-KL preserves fallback provenance for challenging targets", {

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
  expect_true(is.list(fit$attempts))
  expect_identical(
    fit$converged,
    fit$status %in% c("converged", "boundary") && fit$usable && fit$verified
  )
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


test_that("discretize_chisq rejects recyclable vector parameters", {
  bad_df <- tryCatch(
    discretize_chisq(10, df = c(2, 5), scale = 1),
    error = identity
  )
  expect_s3_class(bad_df, "dpprior_invalid_input")
  expect_s3_class(bad_df, "dpprior_chisq_parameter_error")
  expect_s3_class(bad_df, "dpprior_length_error")
  expect_identical(bad_df$argument, "df")
  expect_identical(bad_df$value, c(2, 5))
  expect_identical(bad_df$code, "length")

  bad_scale <- tryCatch(
    discretize_chisq(10, df = 2, scale = c(1, 2)),
    error = identity
  )
  expect_s3_class(bad_scale, "dpprior_invalid_input")
  expect_s3_class(bad_scale, "dpprior_chisq_parameter_error")
  expect_s3_class(bad_scale, "dpprior_length_error")
  expect_identical(bad_scale$argument, "scale")
  expect_identical(bad_scale$value, c(1, 2))
  expect_identical(bad_scale$code, "length")
})


test_that("PMF inputs must be ordinary vectors and never flattened arrays", {
  matrix_target <- matrix(rep(0.1, 10L), nrow = 2L)

  direct <- tryCatch(
    DPprior_a2_kl(10L, matrix_target, method = "pmf", M = 20L),
    error = identity
  )
  expect_s3_class(direct, "dpprior_pmf_error")
  expect_s3_class(direct, "dpprior_type_error")
  expect_identical(direct$code, "vector_required")
  expect_identical(direct$value, matrix_target)

  constructed <- tryCatch(
    construct_target_pmf(10L, matrix_target),
    error = identity
  )
  expect_s3_class(constructed, "dpprior_pmf_error")
  expect_s3_class(constructed, "dpprior_type_error")
  expect_identical(constructed$code, "vector_required")
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

  # Mathematical zero-support conventions are exact, not epsilon-smoothed.
  expect_equal(kl_divergence_pmf(c(0, 1), c(0.5, 0.5)), log(2))
  expect_identical(kl_divergence_pmf(c(1, 0), c(0, 1)), Inf)
})


test_that("construct_target_pmf works with PMF input", {
  J <- 50

  # Uniform PMF
  target <- rep(1 / J, J)
  result <- construct_target_pmf(J, target)

  expect_length(result$pmf, J)
  expect_equal(sum(result$pmf), 1, tolerance = 1e-10)
  expect_true(!is.null(result$mu_K))
  expect_true(!is.null(result$var_K))

  # Verify moments
  k_vals <- 1:J
  expect_equal(result$mu_K, sum(k_vals * result$pmf), tolerance = 1e-10)
})


test_that("constructor and direct solver share stable centered PMF moments", {
  J <- 500L
  epsilon <- 1e-12
  target <- numeric(J)
  target[J - 1L] <- epsilon
  target[J] <- 1 - epsilon

  constructed <- construct_target_pmf(J, target)
  direct <- suppressWarnings(DPprior_a2_kl(
    J, target, method = "pmf", M = 20L,
    max_iter = 1L, fallback_max_iter = 1L
  ))
  naive_variance <- sum(seq_len(J)^2 * target) -
    sum(seq_len(J) * target)^2

  expect_lt(naive_variance, 0)
  expect_gt(constructed$var_K, 0)
  expect_lt(
    abs(constructed$var_K - epsilon * (1 - epsilon)), 1e-24
  )
  expect_identical(direct$target$K$implied$mean, constructed$mu_K)
  expect_identical(direct$target$K$implied$variance, constructed$var_K)
})

test_that("length J+1 target PMF validates K=0 entry before dropping", {
  J <- 20L
  expected <- rep(1 / J, J)

  result <- DPprior:::.a2_kl_normalize_pmf(c(0, expected), J)
  expect_equal(result, expected, tolerance = 1e-12)

  expect_error(
    DPprior:::.a2_kl_normalize_pmf(c(1e-13, expected), J),
    class = "dpprior_pmf_support_error"
  )

  expect_error(
    DPprior:::.a2_kl_normalize_pmf(c(0.1, rep(1, J)), J),
    class = "dpprior_pmf_support_error"
  )
  expect_error(
    DPprior:::.a2_kl_normalize_pmf(c(NA_real_, rep(1, J)), J),
    "finite and non-missing"
  )
  expect_error(
    DPprior:::.a2_kl_normalize_pmf(c(Inf, rep(1, J)), J),
    "finite and non-missing"
  )
  expect_error(
    DPprior:::.a2_kl_normalize_pmf(c(-0.1, rep(1, J)), J),
    "non-negative"
  )
  expect_error(
    DPprior:::.a2_kl_normalize_pmf(c(0, rep(0, J)), J),
    "positive total mass on K=1:J"
  )
})

test_that("construct_target_pmf rejects positive K=0 mass", {
  J <- 20L

  expect_error(
    construct_target_pmf(J, c(0.25, rep(0.75 / J, J))),
    class = "dpprior_pmf_support_error"
  )
})

test_that("exported A2-KL PMF paths reject positive K=0 mass", {
  J <- 20L
  target_pmf <- c(0.25, rep(0.75 / J, J))

  expect_error(
    DPprior_a2_kl(J, target_pmf, method = "pmf", M = 20L),
    class = "dpprior_pmf_support_error"
  )
  expect_error(
    kl_divergence_K(target_pmf, a = 2, b = 1, J = J, M = 20L),
    class = "dpprior_pmf_support_error"
  )
})

test_that("A2-KL moment workflows reject mu_K at upper support boundary", {
  expect_error(
    construct_target_pmf(20, list(mu_K = 20, var_K = 1)),
    "mu_K must be < J"
  )
  expect_error(
    DPprior_a2_kl(20, list(mu_K = 20, var_K = 1), method = "chisq"),
    "mu_K must be < J"
  )
})


test_that("A2-KL target validation is typed and preserves rejected values", {
  fractional_J <- tryCatch(
    construct_target_pmf(2.7, c(0.5, 0.5)),
    error = identity
  )
  expect_s3_class(fractional_J, "dpprior_invalid_input")
  expect_s3_class(fractional_J, "dpprior_integer_error")
  expect_identical(fractional_J$argument, "J")
  expect_identical(fractional_J$value, 2.7)
  expect_identical(fractional_J$code, "integer")

  lower_support <- tryCatch(
    DPprior_a2_kl(
      10, list(mu_K = 1, var_K = 0.1), method = "chisq"
    ),
    error = identity
  )
  expect_s3_class(lower_support, "dpprior_invalid_input")
  expect_s3_class(lower_support, "dpprior_target_moment_support_error")
  expect_s3_class(lower_support, "dpprior_bounds_error")
  expect_identical(lower_support$argument, "mu_K")
  expect_identical(lower_support$value, 1)
  expect_identical(lower_support$code, "support_lower")

  upper_support <- tryCatch(
    construct_target_pmf(10, list(mu_K = 10, var_K = 0.1)),
    error = identity
  )
  expect_s3_class(upper_support, "dpprior_invalid_input")
  expect_s3_class(upper_support, "dpprior_target_moment_support_error")
  expect_identical(upper_support$argument, "mu_K")
  expect_identical(upper_support$value, 10)
  expect_identical(upper_support$code, "support_upper")

  bad_mu_type <- tryCatch(
    construct_target_pmf(10, list(mu_K = "2", var_K = 1)),
    error = identity
  )
  expect_s3_class(bad_mu_type, "dpprior_invalid_input")
  expect_s3_class(bad_mu_type, "dpprior_target_moment_error")
  expect_s3_class(bad_mu_type, "dpprior_type_error")
  expect_identical(bad_mu_type$argument, "mu_K")
  expect_identical(bad_mu_type$value, "2")
  expect_identical(bad_mu_type$code, "type")

  bad_var_nonfinite <- tryCatch(
    DPprior_a2_kl(
      10, list(mu_K = 2, var_K = Inf), method = "chisq"
    ),
    error = identity
  )
  expect_s3_class(bad_var_nonfinite, "dpprior_invalid_input")
  expect_s3_class(bad_var_nonfinite, "dpprior_target_moment_error")
  expect_s3_class(bad_var_nonfinite, "dpprior_nonfinite_error")
  expect_identical(bad_var_nonfinite$argument, "var_K")
  expect_identical(bad_var_nonfinite$value, Inf)
  expect_identical(bad_var_nonfinite$code, "nonfinite")

  nonpositive_var <- tryCatch(
    construct_target_pmf(10, list(mu_K = 2, var_K = 0)),
    error = identity
  )
  expect_s3_class(nonpositive_var, "dpprior_invalid_input")
  expect_s3_class(nonpositive_var, "dpprior_target_moment_error")
  expect_s3_class(nonpositive_var, "dpprior_bounds_error")
  expect_identical(nonpositive_var$argument, "var_K")
  expect_identical(nonpositive_var$value, 0)
  expect_identical(nonpositive_var$code, "positive")

  bad_structure <- tryCatch(
    construct_target_pmf(10, list(mu_K = 2)),
    error = identity
  )
  expect_s3_class(bad_structure, "dpprior_invalid_input")
  expect_s3_class(bad_structure, "dpprior_target_structure_error")
  expect_identical(bad_structure$argument, "target")
  expect_identical(bad_structure$value, list(mu_K = 2))
  expect_identical(bad_structure$code, "target_structure")

  fractional_M <- tryCatch(
    kl_divergence_K(rep(0.1, 10), a = 2, b = 1, J = 10, M = 20.5),
    error = identity
  )
  expect_s3_class(fractional_M, "dpprior_invalid_input")
  expect_s3_class(fractional_M, "dpprior_a2_kl_control_error")
  expect_s3_class(fractional_M, "dpprior_integer_error")
  expect_identical(fractional_M$argument, "M")
  expect_identical(fractional_M$value, 20.5)
  expect_identical(fractional_M$code, "integer")
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


test_that("exact public KL rejects malformed or silently repairable inputs", {
  valid <- c(0.25, 0.75)

  expect_error(
    kl_divergence_pmf(c(-0.1, 1.1), valid),
    class = "dpprior_pmf_negative_error"
  )
  expect_error(
    kl_divergence_pmf(c(NA_real_, 1), valid),
    class = "dpprior_pmf_nonfinite_error"
  )
  expect_error(
    kl_divergence_pmf(c(Inf, 0), valid),
    class = "dpprior_pmf_nonfinite_error"
  )
  expect_error(
    kl_divergence_pmf(c(0, 0), valid),
    class = "dpprior_pmf_zero_mass_error"
  )
  expect_error(
    kl_divergence_pmf(c(0.2, 0.8), c(0.1, 0.2, 0.7)),
    class = "dpprior_pmf_length_error"
  )
  expect_error(
    kl_divergence_pmf(c(2, 8), c(3, 7)),
    class = "dpprior_pmf_normalization_error"
  )
  expect_error(
    kl_divergence_pmf(valid, valid, eps = 1e-15),
    class = "dpprior_kl_smoothing_error"
  )
})


test_that("custom target PMFs are validated but never normalized", {
  J <- 10L
  normalized <- rep(1 / J, J)
  result <- construct_target_pmf(J, normalized)

  expect_identical(result$pmf, normalized)
  expect_identical(
    result$input_provenance$normalization,
    "validated_not_modified"
  )
  expect_error(
    construct_target_pmf(J, rep(1, J)),
    class = "dpprior_pmf_normalization_error"
  )
  expect_error(
    DPprior_a2_kl(J, c(0, rep(1, J)), method = "pmf", M = 20L),
    class = "dpprior_pmf_normalization_error"
  )
  expect_error(
    kl_divergence_K(rep(1, J), a = 2, b = 1, J = J, M = 20L),
    class = "dpprior_pmf_normalization_error"
  )
})


test_that("A2-KL result separates optimizer exit from verified adequacy", {
  fit <- DPprior_a2_kl(
    J = 50, target = list(mu_K = 5, var_K = 8), method = "chisq"
  )

  expect_true(fit$status %in% c(
    "converged", "boundary", "approximate", "infeasible", "failed"
  ))
  expect_identical(
    fit$converged,
    fit$status %in% c("converged", "boundary") && fit$usable && fit$verified
  )
  expect_true(fit$diagnostics$optim$optimizer_converged)
  expect_true(fit$verification$components$pmf_adequacy$passed)
  expect_true(fit$verification$components$order_stability$passed)
  expect_true(fit$verification$components$target_identity$passed)
  expect_true(fit$verification$components$candidate_selection$passed)
  expect_named(
    fit$verification$components$pmf_adequacy$value,
    c(
      "selected.kl", "selected.l1", "selected.mean_scaled",
      "selected.variance_scaled", "refined.kl", "refined.l1",
      "refined.mean_scaled", "refined.variance_scaled"
    )
  )
  expect_named(fit$verification$stability$delta, "pmf.l1")
  expect_false(
    "pmf.l1" %in% names(fit$verification$components$pmf_adequacy$value)
  )
  expect_equal(fit$verification$settings$M_selected, 80L)
  expect_equal(fit$verification$settings$M_verification, 160L)
  expect_identical(
    unname(unlist(fit$tolerances$distribution$adequacy[1:4])),
    c(0.015, 0.11, 0.01, 0.065)
  )
  expect_true(all(c("kl", "l1", "mean", "variance") %in%
                    names(fit$residuals$distribution)))
  expect_true(all(c("raw", "scaled", "ratios") %in%
                    names(fit$diagnostics$adequacy)))
  expect_identical(names(fit)[1:18], c(
    "schema", "object_type", "mode", "method", "J", "status", "usable",
    "verified", "message", "parameters", "target", "achieved",
    "residuals", "tolerances", "computation", "verification",
    "provenance", "compatibility"
  ))
  expect_invisible(DPprior:::.dpprior_validate_result_v1(fit))
})


test_that("higher-order verification never replaces selected-M public values", {
  fit <- DPprior_a2_kl(
    J = 50, target = list(mu_K = 5, var_K = 8), method = "chisq"
  )

  expect_identical(fit$achieved$K$M, 80L)
  expect_identical(
    fit$achieved$K$pmf,
    fit$verification$selected_snapshot$achieved$K$pmf
  )
  expect_identical(
    fit$residuals$distribution$kl,
    fit$verification$selected_snapshot$residuals$distribution$kl
  )
  expect_identical(
    fit$residuals$distribution$l1,
    fit$verification$selected_snapshot$residuals$distribution$l1
  )
  expect_false(identical(
    fit$achieved$K$pmf,
    fit$verification$verifier_snapshot$achieved$K$pmf
  ))
  expect_identical(fit$verification$verifier_snapshot$M, 160L)
  expect_equal(
    unname(unlist(fit$verification$verifier_snapshot$residuals$
                    distribution[c("kl", "l1")])),
    c(
      fit$compatibility$views$a2_kl_v0$verification$
        verification_metrics$kl,
      fit$compatibility$views$a2_kl_v0$verification$
        verification_metrics$l1
    )
  )
  expect_true(fit$verification$passed)
})


test_that("optimizer exit zero cannot hide target-family inadequacy", {
  target <- numeric(30L)
  target[c(2L, 25L)] <- 0.5

  fit <- DPprior_a2_kl(30L, target, method = "pmf")

  expect_true(fit$diagnostics$optim$optimizer_converged)
  expect_false(fit$verification$components$pmf_adequacy$passed)
  expect_identical(fit$status, "approximate")
  expect_false(fit$verified)
  expect_false(fit$usable)
  expect_match(fit$message, "adequacy failed")
})


test_that("fail-first optimizer fallback is explicit and independently verified", {
  local_mocked_bindings(
    .a2_kl_run_lbfgsb = function(par, fn, lower, upper, control) {
      list(
        par = par, value = fn(par),
        counts = setNames(c(1L, 1L), c("function", "gradient")),
        convergence = 52L, message = "injected primary failure"
      )
    },
    .package = "DPprior"
  )
  target <- dbinom(1:30, size = 30, prob = 0.15)
  target <- target / sum(target)

  fit <- DPprior_a2_kl(30L, target, method = "pmf")

  expect_true(fit$computation$fallback$attempted)
  expect_true(fit$computation$fallback$used)
  expect_true(fit$provenance$is_fallback)
  expect_identical(fit$provenance$selected_method, "A2-KL")
  expect_identical(
    fit$computation$attempts[[match(
      fit$computation$selected_attempt_id,
      vapply(fit$computation$attempts, `[[`, character(1), "id")
    )]]$method,
    "nlminb"
  )
  expect_identical(fit$status, "converged")
  expect_true(fit$verified)
  expect_identical(
    vapply(fit$attempts, `[[`, character(1), "method"),
    c("A2-MN", "L-BFGS-B", "nlminb")
  )
  expect_identical(fit$attempts[[2L]]$exit_code, 52L)
  expect_identical(fit$attempts[[3L]]$exit_code, 0L)
})


test_that("fallback-to-initialization is never labelled ordinary success", {
  local_mocked_bindings(
    .a2_kl_run_lbfgsb = function(...) stop("injected primary failure"),
    .a2_kl_run_nlminb = function(...) stop("injected fallback failure"),
    .package = "DPprior"
  )
  target <- dbinom(1:30, size = 30, prob = 0.15)
  target <- target / sum(target)

  fit <- DPprior_a2_kl(30L, target, method = "pmf")

  expect_identical(fit$status, "approximate")
  expect_false(fit$converged)
  expect_false(fit$verified)
  expect_false(fit$usable)
  expect_identical(fit$provenance$selected_method, "A2-KL")
  selected <- fit$computation$candidate_evaluations[[match(
    fit$computation$selected_candidate_id,
    vapply(
      fit$computation$candidate_evaluations, `[[`, character(1), "id"
    )
  )]]
  expect_identical(selected$generator, "initialization")
  expect_false(selected$optimizer_supported)
  expect_identical(fit$computation$termination$source, "candidate_evaluation")
  expect_match(fit$message, "initialization")
  expect_true(fit$computation$fallback$attempted)
  expect_true(all(vapply(
    fit$attempts[2:3], function(x) !is.null(x$error), logical(1)
  )))
})


test_that("verified solver-bound candidates receive boundary status", {
  J <- 20L
  logS <- compute_log_stirling(J)
  target <- pmf_K_marginal(
    J, exp(-1), exp(-0.95), logS, M = 160L
  )[-1L]
  target <- target / sum(target)

  fit <- DPprior_a2_kl(
    J, target, method = "pmf", M = 120L,
    log_bounds = c(-1, -0.9), boundary_tol = 1e-6
  )

  expect_identical(fit$status, "boundary")
  expect_true(fit$usable)
  expect_true(fit$verified)
  expect_true(fit$converged)
  expect_identical(fit$computation$termination$code, "boundary")
  expect_identical(fit$verification$selected_snapshot$M, 120L)
  expect_identical(fit$verification$verifier_snapshot$M, 240L)
  expect_false(fit$compatibility$views$a2_kl_v0$converged)
  expect_true(fit$diagnostics$optim$boundary_hit)
  expect_true(length(fit$diagnostics$optim$boundary_sides) >= 1L)
})


test_that("unavailable higher-order verification returns an honest status", {
  J <- 5L
  target <- rep(1 / J, J)

  fit <- DPprior_a2_kl(
    J, target, method = "pmf", M = 257L, max_iter = 1L
  )

  expect_identical(fit$status, "approximate")
  expect_false(fit$usable)
  expect_false(fit$verified)
  expect_true(fit$verification$performed)
  expect_false(fit$verification$passed)
  expect_identical(
    fit$verification$settings$M_verification_required, 257L
  )
  expect_identical(fit$verification$settings$M_verification, 257L)
  expect_true(fit$computation$resources$same_order_verifier_only)
  expect_identical(
    fit$computation$resources$source_required_M_verification, 514L
  )
  expect_false(
    fit$compatibility$views$a2_kl_v0$verification$performed
  )
  expect_match(fit$message, "independent PMF verification did not pass")
})


test_that("unknown dots and malformed chi-square targets fail before fitting", {
  side_effect <- 0L
  unknown <- tryCatch(
    DPprior_a2_kl(
      10L, rep(0.1, 10L), method = "pmf",
      typo = {
        side_effect <<- side_effect + 1L
        stop("unknown promise was forced")
      }
    ),
    error = identity
  )
  expect_s3_class(unknown, "dpprior_a2_kl_control_error")
  expect_identical(unknown$code, "unknown_control")
  expect_identical(side_effect, 0L)

  extra <- tryCatch(
    DPprior_a2_kl(
      10L, list(mu_K = 3, var_K = 2, var_k = 999), method = "chisq"
    ),
    error = identity
  )
  expect_s3_class(extra, "dpprior_target_structure_error")
  expect_identical(extra$code, "target_structure")

  duplicated <- structure(
    list(3, 2, 999), names = c("mu_K", "var_K", "var_K")
  )
  duplicate_error <- tryCatch(
    DPprior_a2_kl(10L, duplicated, method = "chisq"),
    error = identity
  )
  expect_s3_class(duplicate_error, "dpprior_target_structure_error")
  expect_identical(duplicate_error$code, "target_structure")
})


test_that("chi-square target verification uses a separate reconstruction", {
  fit <- DPprior_a2_kl(
    30L, list(mu_K = 5, var_K = 8), method = "chisq"
  )

  expect_true(fit$verification$components$target_identity$passed)
  objective <- fit$target$K$derivation$request_to_normalized$
    evidence$A2_KL_objective
  expect_identical(objective$method, "chisq")
  expect_identical(objective$source, "A2_KL_backend_target")
  legacy_verification <- fit$compatibility$views$a2_kl_v0$verification$target
  expect_true(legacy_verification$independently_reconstructed)
  expect_identical(
    legacy_verification$method,
    "separate_continuity_edge_cdf_difference_reconstruction"
  )
  expect_true(legacy_verification$passed)
  expect_lte(legacy_verification$l1_difference, .TOL_PMF_SUM)
})


test_that("A2-KL canonical result rejects one-field authority mutations", {
  target <- dbinom(1:30, size = 30, prob = 0.15)
  target <- target / sum(target)
  fit <- DPprior_a2_kl(30L, target, method = "pmf")
  mutations <- list(
    schema_version = function(x) {
      x$schema$version <- 1
      x
    },
    target_used = function(x) {
      x$target$K$used$pmf <- rev(x$target$K$used$pmf)
      x
    },
    achieved_pmf = function(x) {
      x$achieved$K$pmf <- rev(x$achieved$K$pmf)
      x
    },
    public_kl = function(x) {
      x$residuals$distribution$kl <-
        x$residuals$distribution$kl + 1e-4
      x
    },
    selected_snapshot = function(x) {
      x$verification$selected_snapshot$achieved$K$mean <-
        x$verification$selected_snapshot$achieved$K$mean + 1e-3
      x
    },
    verifier_order = function(x) {
      x$verification$verifier_snapshot$M <-
        x$verification$verifier_snapshot$M + 1L
      x
    },
    adequacy_check = function(x) {
      x$verification$components$pmf_adequacy$value[[1L]] <-
        x$verification$components$pmf_adequacy$value[[1L]] + 1e-3
      x
    },
    candidate_objective = function(x) {
      selected <- match(
        x$computation$selected_candidate_id,
        vapply(
          x$computation$candidate_evaluations, `[[`, character(1), "id"
        )
      )
      x$computation$candidate_evaluations[[selected]]$fresh_objective <-
        x$computation$candidate_evaluations[[selected]]$fresh_objective +
        1e-3
      x
    },
    selected_id = function(x) {
      x$computation$selected_candidate_id <- "candidate-forged"
      x
    },
    alias = function(x) {
      x$a <- x$a + 1
      x
    }
  )

  for (name in names(mutations)) {
    condition <- .a2_kl_expect_schema_rejection(mutations[[name]](fit))
    expect_false(
      inherits(condition, "simpleError") &&
        !inherits(condition, "dpprior_schema_error"),
      label = name
    )
  }

  chisq <- DPprior_a2_kl(
    30L, list(mu_K = 5, var_K = 8), method = "chisq"
  )
  chisq$target$K$derivation$request_to_normalized$evidence$
    A2_KL_objective$pmf <- rev(
      chisq$target$K$derivation$request_to_normalized$evidence$
        A2_KL_objective$pmf
    )
  .a2_kl_expect_schema_rejection(chisq)
})


test_that("A2-KL canonical result round-trips through serialization", {
  target <- dbinom(1:25, size = 25, prob = 0.16)
  target <- target / sum(target)
  fit <- DPprior_a2_kl(25L, target, method = "pmf", M = 40L)
  round_trip <- unserialize(serialize(fit, NULL, version = 3L))

  expect_identical(round_trip, fit)
  .a2_kl_expect_valid_result(round_trip, fit$status)
  expect_identical(
    round_trip$target$K$request,
    list(J = 25L, pmf = target)
  )
  expect_identical(
    round_trip$verification$verifier_snapshot$parameters,
    round_trip$parameters
  )
})


test_that("A2-KL S3 consumers use canonical fields and ignore legacy science", {
  fit <- DPprior_a2_kl(
    20L, list(mu_K = 4, var_K = 5), method = "chisq", M = 40L
  )
  legacy <- fit$compatibility$views$a2_kl_v0
  expect_identical(legacy$authority, "non_authoritative")
  expect_true(legacy$lossy)
  expect_identical(
    legacy$consumer_policy,
    "ignored_by_scientific_and_decision_consumers"
  )
  expect_identical(
    fit$compatibility$deprecations$a2_kl_v0[c(
      "authority", "lossy", "consumer_policy"
    )],
    legacy[c("authority", "lossy", "consumer_policy")]
  )

  output <- capture.output(print(fit))
  compact <- summary(fit, print_output = FALSE)
  expect_true(any(grepl("A2-KL", output, fixed = TRUE)))
  expect_s3_class(compact, "summary.DPprior_fit")
  expect_identical(compact$gamma_prior$a, fit$parameters$a)
  expect_identical(compact$gamma_prior$b, fit$parameters$b)

  forged <- fit
  forged$compatibility$views$a2_kl_v0$a <- 999
  forged$compatibility$views$a2_kl_v0$b <- 999
  forged$compatibility$views$a2_kl_v0$status <- "failed"
  forged$compatibility$views$a2_kl_v0$achieved$mu_K <- 999
  expect_invisible(.dpprior_validate_result_v1(forged))
  forged_summary <- summary(forged, print_output = FALSE)
  expect_identical(forged_summary$gamma_prior, compact$gamma_prior)
  expect_identical(forged_summary$achieved, compact$achieved)
})


test_that("A2-KL public controls cannot relax canonical truth caps", {
  target <- rep(0.1, 10L)
  relaxed <- list(
    tol = 1e-5,
    kl_tol = 0.02,
    l1_tol = 0.2,
    mean_scaled_tol = 0.02,
    var_scaled_tol = 0.1,
    pmf_abs_tol = 1e-9,
    pmf_rel_tol = 1e-7,
    boundary_tol = 1e-5
  )
  for (field in names(relaxed)) {
    args <- c(
      list(J = 10L, target = target, method = "pmf", M = 20L),
      setNames(list(relaxed[[field]]), field)
    )
    condition <- tryCatch(do.call(DPprior_a2_kl, args), error = identity)
    expect_s3_class(condition, "dpprior_a2_kl_control_error")
    expect_identical(condition$code, "canonical_control_cap", info = field)
    expect_identical(condition$argument, field, info = field)
  }
})


test_that("A2-KL exact objective does not accept epsilon smoothing", {
  target <- rep(0.1, 10L)
  expect_error(
    DPprior_a2_kl(10L, target, method = "pmf", M = 20L, eps = 1e-12),
    class = "dpprior_kl_smoothing_error"
  )
})

test_that("A2-KL method requires an ordinary character scalar", {
  target <- rep(0.1, 10L)
  malformed <- list(
    matrix("pmf", nrow = 1L),
    array("pmf", dim = c(1L, 1L, 1L)),
    structure("pmf", class = "dpprior_test_character")
  )
  for (value in malformed) {
    expect_error(
      DPprior_a2_kl(10L, target, method = value),
      class = "dpprior_type_error"
    )
  }
  expect_error(
    DPprior_a2_kl(
      10L, target, method = "pmf",
      log_bounds = matrix(c(-15, 15), nrow = 1L)
    ),
    class = "dpprior_a2_kl_control_error"
  )
})
