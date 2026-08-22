# =============================================================================
# Tests for Module 11: canonical A2-MN calibrated-moment solver
# =============================================================================

.a2_expect_valid_result <- function(fit, status = NULL) {
  expect_s3_class(fit, "DPprior_fit")
  expect_s3_class(fit, "dpprior_result")
  expect_identical(fit$schema, list(name = "dpprior.result", version = 1L))
  expect_identical(fit$object_type, "fit")
  expect_identical(fit$mode, "a2_moment")
  if (!is.null(status)) expect_identical(fit$status, status)
  report <- .dpprior_validate_result_v1(fit, collect = TRUE)
  details <- if (length(report$errors)) {
    paste(vapply(report$errors, conditionMessage, character(1)), collapse = "\n")
  } else {
    ""
  }
  expect_true(report$valid, info = details)
  expect_length(report$errors, 0L)
  expect_identical(
    fit$converged,
    fit$status %in% c("converged", "boundary") &&
      fit$usable && fit$verified
  )
  invisible(fit)
}


.a2_expect_schema_rejection <- function(object) {
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
  expect_false(is.null(condition$expected))
  invisible(condition)
}


test_that("A2-MN emits a valid canonical result with exact baseline identity", {
  fit <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8)

  .a2_expect_valid_result(fit, "converged")
  expect_true(fit$usable)
  expect_true(fit$verified)
  expect_identical(fit$method, "A2-MN")
  expect_identical(
    sprintf("%.17g", c(
      a = fit$a, b = fit$b, mean = fit$achieved$K$mean,
      variance = fit$achieved$K$variance,
      residual = fit$fit$residual,
      scaled = fit$fit$scaled_residual
    )),
    c(
      "2.0360925614372438", "1.6050540689499431",
      "5.0000000000011573", "8.0000000000018794",
      "2.2071305212604966e-12", "2.0882061512060722e-05"
    )
  )
  expect_identical(fit$achieved$K$M, 80L)
  expect_identical(fit$fit$mu_K, fit$achieved$K$mean)
  expect_identical(fit$fit$var_K, fit$achieved$K$variance)
})


test_that("canonical target, selected order, and verifier order are distinct", {
  fit <- DPprior_a2_newton(50, 5, 8)
  orders <- fit$computation$orders
  selected <- fit$verification$selected_snapshot
  verifier <- fit$verification$verifier_snapshot

  expect_identical(fit$target$K$request,
                   list(J = 50L, mu_K = 5, var_K = 8))
  expect_identical(fit$target$K$normalized$J, 50L)
  expect_identical(fit$target$K$used$J, 50L)
  expect_identical(fit$target$K$used$mu_K, 5)
  expect_identical(fit$target$K$used$var_K, 8)
  expect_identical(orders$M_requested, 80L)
  expect_identical(orders$M_selected, 80L)
  expect_identical(orders$M_verification_required, 160L)
  expect_identical(orders$M_verification_used, 160L)
  expect_identical(fit$verification$settings, list(
    M_selected = 80L,
    M_verification_required = 160L,
    M_verification = 160L
  ))
  expect_identical(selected$M, 80L)
  expect_identical(verifier$M, 160L)
  expect_identical(fit$achieved, selected$achieved)
  expect_identical(fit$residuals, selected$residuals)
  expect_identical(selected$parameters, fit$parameters)
  expect_identical(verifier$parameters, fit$parameters)
  expect_false(identical(selected$achieved$K$M, verifier$achieved$K$M))
})


test_that("adequacy and stability are separate componentwise truth controls", {
  target <- c(mean = 5, variance = 100)
  baseline <- .a2_residual_contract(
    target, target, abs_tol = 1e-8, rel_tol = 1e-8
  )
  observed <- target
  observed[["mean"]] <- target[["mean"]] +
    1.1 * baseline$tolerance[["mean"]]
  result <- .a2_residual_contract(
    observed, target, abs_tol = 1e-8, rel_tol = 1e-8
  )

  expect_lt(result$standardized_norm, 1)
  expect_false(result$component_pass[["mean"]])
  expect_true(result$component_pass[["variance"]])
  expect_false(result$passed)

  fit <- DPprior_a2_newton(50, 5, 8)
  expect_identical(fit$tolerances$K_adequacy, list(
    absolute = 1e-8, relative = 1e-8,
    scale_formula = "max(abs(target),1)"
  ))
  expect_identical(fit$tolerances$K_stability, list(
    absolute = 1e-10, relative = 1e-8, scale_floor = 1
  ))
  expect_true(fit$verification$components$residual_adequacy$passed)
  expect_true(fit$verification$components$order_stability$passed)
  expect_true(fit$verification$stability$passed)
  expect_named(
    fit$verification$components$residual_adequacy$value,
    c(
      "selected.mean", "selected.variance", "refined.mean",
      "refined.variance", "selected.standardized_norm",
      "refined.standardized_norm"
    )
  )
})


test_that("canonical evidence is fresh and ignores legacy verifier self-report", {
  fake_verifier <- function(...) {
    list(
      status = "approximate", passed = FALSE, available = TRUE,
      performed = TRUE, reason = "injected_legacy_failure",
      message = "injected legacy verifier", M_selected = 80L,
      M_required = 160L, M_verification = 160L,
      selected = c(mean = 500, variance = 800), recomputed = NULL,
      residuals = NULL, stability = NULL
    )
  }
  testthat::local_mocked_bindings(
    .a2_verify_candidate = fake_verifier,
    .package = "DPprior"
  )

  fit <- DPprior_a2_newton(50, 5, 8)
  .a2_expect_valid_result(fit, "converged")
  expect_true(fit$verified)
  expect_identical(
    fit$verification$reason, "independent_verification_passed"
  )
  expect_identical(
    fit$compatibility$views$legacy_v2$source_status, "approximate"
  )
  expect_false(fit$compatibility$views$legacy_v2$source_verified)
  expect_identical(
    fit$compatibility$views$legacy_v2$authority, "non_authoritative"
  )
})


test_that("verification order follows max(2M, M+40)", {
  fit80 <- DPprior_a2_newton(50, 5, 8, M = 80)
  fit81 <- DPprior_a2_newton(50, 5, 8, M = 81)

  expect_identical(fit80$computation$orders$M_verification_required, 160L)
  expect_identical(fit80$computation$orders$M_verification_used, 160L)
  expect_identical(fit81$computation$orders$M_verification_required, 162L)
  expect_identical(fit81$computation$orders$M_verification_used, 162L)
  expect_identical(fit80$diagnostics$M_verify_required, 160L)
  expect_identical(fit81$diagnostics$M_verify_required, 162L)
  expect_true(fit80$verified)
  expect_true(fit81$verified)
})


test_that("invalid quadrature orders are rejected with typed conditions", {
  expect_error(
    DPprior_a2_newton(50, 5, 8, M = 80.5),
    class = "dpprior_integer_error"
  )
  expect_error(
    DPprior_a2_newton(50, 5, 8, M = 80, M_verify = 159),
    class = "dpprior_a2_verification_error"
  )
  expect_error(
    DPprior_a2_newton(50, 5, 8, M = 81, M_verify = 161),
    class = "dpprior_a2_verification_error"
  )
  expect_error(
    DPprior_a2_newton(50, 5, 8, M = 513),
    class = "dpprior_bounds_error"
  )
})


test_that("M above the verifier ceiling fails closed without a public candidate", {
  reference <- DPprior_a2_newton(50, 5, 8)
  condition <- tryCatch(
    DPprior_a2_newton(
      50, 5, 8, a0 = reference$a, b0 = reference$b,
      M = 257, max_iter = 1, tol_F = 1
    ),
    error = identity
  )

  expect_s3_class(condition, "dpprior_a2_no_candidate")
  expect_s3_class(condition, "dpprior_a2_verification_error")
  expect_identical(
    condition$code, "a2_verification_order_exceeds_ceiling"
  )
  expect_identical(
    condition$action, "refit_with_smaller_M_at_or_below_256"
  )
  expect_match(conditionMessage(condition), "M <= 256")
  expect_match(conditionMessage(condition), "non-authoritative")

  failed <- condition$result
  .a2_expect_valid_result(failed, "failed")
  expect_false(failed$usable)
  expect_false(failed$verified)
  expect_null(failed$parameters)
  expect_identical(failed$achieved, list())
  expect_identical(failed$residuals, list())
  expect_false("a" %in% names(failed))
  expect_false("b" %in% names(failed))
  expect_identical(
    failed$computation$orders$M_verification_required, 514L
  )
  expect_null(failed$computation$orders$M_verification_used)
  expect_identical(failed$computation$termination$code, "no_candidate")
  expect_identical(failed$computation$termination$source, "no_candidate")
  expect_identical(
    failed$computation$attempts[[1L]]$error$class,
    "dpprior_a2_verification_error"
  )
  expect_identical(
    failed$computation$attempts[[1L]]$error$code,
    "a2_verification_order_exceeds_ceiling"
  )
  legacy <- failed$compatibility$views$legacy_v2
  expect_identical(legacy$authority, "non_authoritative")
  expect_true(legacy$lossy)
  expect_identical(legacy$numerical_candidate$a, reference$a)
  expect_identical(legacy$numerical_candidate$b, reference$b)
  expect_identical(legacy$required_action, condition$action)

  expect_error(
    DPprior_a2_newton(
      50, 5, 8, a0 = reference$a, b0 = reference$b,
      M = 257, M_verify = 512, max_iter = 1
    ),
    "verification is unavailable",
    class = "dpprior_a2_verification_error"
  )
})


test_that("typical and high-dispersion scenarios retain valid verification", {
  scenarios <- list(
    list(J = 25, mu = 3, variance = 4),
    list(J = 50, mu = 3, variance = 10),
    list(J = 100, mu = 10, variance = 15)
  )
  for (scenario in scenarios) {
    fit <- DPprior_a2_newton(
      scenario$J, scenario$mu, scenario$variance
    )
    label <- sprintf(
      "J=%d, mean=%g, variance=%g",
      scenario$J, scenario$mu, scenario$variance
    )
    expect_identical(fit$status, "converged", label = label)
    expect_true(fit$verified, label = label)
    expect_true(
      .dpprior_validate_result_v1(fit, collect = TRUE)$valid,
      label = label
    )
  }
})


test_that("A1 projection is initialization-only and target identity is fixed", {
  fit <- DPprior_a2_newton(50, 5, 3)
  initialization <- fit$compatibility$views$legacy_v2$initialization

  .a2_expect_valid_result(fit)
  expect_identical(fit$target$K$request$mu_K, 5)
  expect_identical(fit$target$K$request$var_K, 3)
  expect_identical(fit$target$K$normalized$var_K, 3)
  expect_identical(fit$target$K$used$var_K, 3)
  expect_identical(initialization$method, "A1")
  expect_identical(
    initialization$projection_policy, "nearest_for_start_only"
  )
  expect_true(initialization$target_was_not_projected_for_A2)
  expect_identical(initialization$target_projection$policy, "nearest")
  expect_true(initialization$target_projection$applied)
  expect_identical(
    initialization$target_projection$original_target$var_K, 3
  )
  expect_gt(initialization$target_projection$projected_target$var_K, 3)
})


test_that("A1 failure uses deterministic grid without changing A2 target", {
  expect_error(
    DPprior_a1(10, 9.5, 4, projection = "nearest"),
    class = "dpprior_a1_projection_impossible"
  )
  fit <- DPprior_a2_newton(10, 9.5, 4)
  initialization <- fit$compatibility$views$legacy_v2$initialization

  .a2_expect_valid_result(fit)
  expect_identical(initialization$method, "fixed_log_parameter_grid")
  expect_identical(
    initialization$projection_policy, "A1_failed_no_A2_projection"
  )
  expect_identical(
    initialization$attempt_methods,
    c("A1_nearest_start", "fixed_log_parameter_grid")
  )
  expect_identical(initialization$attempt_statuses, c("failed", "selected"))
  expect_identical(initialization$attempt_reason_codes, c(
    "a1_projection_exceeds_support", "fixed_grid_candidate_selected"
  ))
  expect_identical(fit$target$K$used$mu_K, 9.5)
  expect_identical(fit$target$K$used$var_K, 4)
})


test_that("fallback is explicit and retains exact baseline numerical identity", {
  fit <- DPprior_a2_newton(
    50, 3, 10, max_iter = 2, use_fallback = TRUE
  )
  attempts <- fit$computation$attempts

  .a2_expect_valid_result(fit, "converged")
  expect_identical(fit$method, "A2-MN+NM")
  expect_identical(
    sprintf("%.17g", c(
      a = fit$a, b = fit$b, mean = fit$achieved$K$mean,
      variance = fit$achieved$K$variance
    )),
    c(
      "0.2934715158344367", "0.43384348527125577",
      "2.9999999946841776", "10.000000014665607"
    )
  )
  expect_identical(fit$provenance$selected_method, "A2-MN+NM")
  expect_true(fit$provenance$is_fallback)
  expect_length(attempts, 2L)
  expect_identical(attempts[[1L]]$method, "scaled_log_newton")
  expect_identical(attempts[[1L]]$reason_code, "eligible_not_selected")
  expect_identical(attempts[[2L]]$method, "nelder_mead_log")
  expect_identical(attempts[[2L]]$reason_code, "selected")
  expect_identical(fit$computation$selected_attempt_id, "attempt-2")
  expect_identical(fit$computation$selected_candidate_id, "candidate-2")
  expect_true(fit$computation$fallback$attempted)
  expect_true(fit$computation$fallback$used)
  expect_identical(fit$computation$fallback$trigger_attempt_id, "attempt-1")
  expect_identical(fit$computation$fallback$selected_attempt_id, "attempt-2")
  expect_identical(
    fit$computation$termination$source, "fallback_optimizer"
  )
  expect_identical(
    fit$compatibility$views$legacy_v2$source_attempt_methods,
    c("scaled_log_newton", "nelder_mead_log")
  )
})


test_that("fallback rejection and error paths remain honest canonical attempts", {
  testthat::local_mocked_bindings(
    .a2_nelder_mead = function(start, objective, maxit = 1000L) {
      list(
        par = c(log(100), log(100)), value = 0, convergence = 0L,
        counts = structure(
          c(1L, NA_integer_), names = c("function", "gradient")
        ),
        message = "injected exit zero"
      )
    },
    .package = "DPprior"
  )
  rejected <- DPprior_a2_newton(
    50, 3, 10, max_iter = 1, use_fallback = TRUE
  )
  .a2_expect_valid_result(rejected, "approximate")
  expect_false(rejected$usable)
  expect_false(rejected$verified)
  expect_false(rejected$computation$fallback$used)
  expect_identical(rejected$computation$attempts[[2L]]$exit_code, 0L)
  expect_identical(
    rejected$computation$attempts[[2L]]$reason_code,
    "eligible_not_selected"
  )

  testthat::local_mocked_bindings(
    .a2_nelder_mead = function(...) stop("injected fallback failure"),
    .package = "DPprior"
  )
  failed <- DPprior_a2_newton(
    50, 3, 10, max_iter = 1, use_fallback = TRUE
  )
  .a2_expect_valid_result(failed, "approximate")
  attempt <- failed$computation$attempts[[2L]]
  expect_identical(attempt$reason_code, "optimizer_error")
  expect_identical(attempt$error$class, "simpleError")
  expect_match(attempt$error$message, "injected fallback failure")
  expect_false(failed$computation$fallback$used)
})


test_that("nonconvergence retains a finite but unusable diagnostic candidate", {
  fit <- DPprior_a2_newton(
    50, 5, 8, max_iter = 1, use_fallback = FALSE,
    tol_F = 1e-15, tol_rel = 0
  )
  selected <- fit$computation$candidate_evaluations[[1L]]

  .a2_expect_valid_result(fit, "approximate")
  expect_false(fit$usable)
  expect_false(fit$verified)
  expect_identical(
    sprintf("%.17g", c(
      a = fit$a, b = fit$b, mean = fit$achieved$K$mean,
      variance = fit$achieved$K$variance
    )),
    c(
      "1.1786481619925739", "0.91196785529358704",
      "4.9090450247150317", "10.854546576808852"
    )
  )
  expect_true(fit$verification$performed)
  expect_false(fit$verification$passed)
  expect_identical(fit$computation$termination$code, "approximate")
  expect_identical(
    fit$computation$termination$source, "candidate_evaluation"
  )
  expect_true(selected$selected)
  expect_identical(selected$outcome, "selected_diagnostic")
  expect_false(selected$execution_success)
  expect_false(selected$optimizer_supported)
  expect_true("execution_failed" %in% selected$rejection_codes)
  expect_true("optimizer_unsupported" %in% selected$rejection_codes)
})


test_that("boundary status is usable only with full fresh verification", {
  parameters <- c(a = exp(-15), b = 1)
  target <- exact_K_moments(50, parameters[["a"]], parameters[["b"]], M = 80)
  fit <- DPprior_a2_newton(
    50, target$mean, target$var,
    a0 = parameters[["a"]], b0 = parameters[["b"]],
    max_iter = 1, use_fallback = FALSE
  )

  .a2_expect_valid_result(fit, "boundary")
  expect_true(fit$usable)
  expect_true(fit$verified)
  expect_identical(fit$a, parameters[["a"]])
  expect_identical(fit$b, parameters[["b"]])
  expect_identical(fit$computation$termination$code, "boundary")
  expect_identical(fit$computation$termination$source, "optimizer")
  expect_match(
    fit$computation$termination$boundary_reason, "boundary tolerance"
  )
  expect_true(fit$verification$components$residual_adequacy$passed)
  expect_true(fit$verification$components$order_stability$passed)
})


test_that("relaxed solver controls cannot relax canonical truth controls", {
  fit <- DPprior_a2_newton(
    50, 5, 8, max_iter = 1, use_fallback = FALSE,
    tol_F = 1, tol_rel = 1,
    verification_abs_tol = 1, verification_rel_tol = 1
  )

  .a2_expect_valid_result(fit, "approximate")
  expect_identical(fit$computation$attempts[[1L]]$exit_code, 0L)
  expect_identical(fit$tolerances$K_adequacy$absolute, 1e-8)
  expect_identical(fit$tolerances$K_adequacy$relative, 1e-8)
  expect_identical(fit$tolerances$K_stability$absolute, 1e-10)
  expect_identical(fit$tolerances$K_stability$relative, 1e-8)
  expect_identical(
    fit$verification$reason, "selected_residual_tolerance_not_met"
  )
  expect_false(fit$usable)
  expect_false(fit$verified)
})


test_that("trace, attempts, scaling, and compatibility remain explicit", {
  fit1 <- DPprior_a2_newton(50, 5, 8)
  fit2 <- DPprior_a2_newton(50, 5, 8)
  required_trace <- c(
    "iter", "a", "b", "M1", "V", "residual", "residual_mean",
    "residual_variance", "scaled_mean", "scaled_variance",
    "standardized_norm", "max_budget_ratio", "step", "step_norm",
    "line_search_iterations", "accepted", "det_Jlog",
    "reciprocal_condition", "jacobian_status", "derivative_status",
    "reason_code"
  )

  expect_true(all(required_trace %in% names(fit1$computation$trace)))
  expect_identical(fit1$trace, fit1$computation$trace)
  expect_identical(fit1$attempts, fit1$computation$attempts)
  expect_equal(fit1$trace, fit2$trace)
  expect_identical(fit1$computation$scaling$requested,
                   fit1$computation$scaling$used)
  expect_true(fit1$computation$scaling$fixed_from_input)
  expect_identical(
    fit1$computation$scaling$formula, "max(abs(target),1)"
  )
  legacy <- fit1$compatibility$views$legacy_v2
  expect_identical(legacy$authority, "non_authoritative")
  expect_true(legacy$lossy)
  expect_identical(
    legacy$consumer_policy,
    "ignored_by_scientific_and_decision_consumers"
  )
  expect_identical(legacy$numerical_candidate$a, fit1$a)
  expect_identical(legacy$selected_fit, fit1$fit)
})


test_that("custom initialization is validated and audit-recorded", {
  start <- DPprior_a1(50, 5, 8)
  fit <- DPprior_a2_newton(
    50, 5, 8, a0 = start$a, b0 = start$b
  )
  initialization <- fit$compatibility$views$legacy_v2$initialization

  .a2_expect_valid_result(fit, "converged")
  expect_identical(initialization$method, "user")
  expect_identical(initialization$projection_policy, "not_applicable")
  expect_identical(initialization$attempt_reason_codes, "user_initialization")
  expect_identical(fit$diagnostics$a0, start$a)
  expect_identical(fit$diagnostics$b0, start$b)
  expect_error(
    DPprior_a2_newton(50, 5, 8, a0 = 1, b0 = NULL),
    class = "dpprior_initialization_error"
  )
  expect_error(
    DPprior_a2_newton(50, 5, 8, a0 = -1, b0 = 1),
    class = "dpprior_initialization_error"
  )
})


test_that("canonical result rejects one-field scientific and alias mutations", {
  fit <- DPprior_a2_newton(50, 5, 8)
  mutations <- list(
    schema_version = function(x) {
      x$schema$version <- 1
      x
    },
    target_used = function(x) {
      x$target$K$used$mu_K <- x$target$K$used$mu_K + 0.1
      x
    },
    achieved_order = function(x) {
      x$achieved$K$M <- 81L
      x
    },
    public_residual = function(x) {
      x$residuals$K$mean <- x$residuals$K$mean + 1e-4
      x
    },
    tolerance_authority = function(x) {
      x$tolerances$K_adequacy$absolute <- 1e-3
      x
    },
    selected_snapshot = function(x) {
      x$verification$selected_snapshot$achieved$K$mean <-
        x$verification$selected_snapshot$achieved$K$mean + 1e-3
      x
    },
    verifier_order = function(x) {
      x$verification$verifier_snapshot$M <- 161L
      x
    },
    stability = function(x) {
      x$verification$stability$delta[[1L]] <-
        x$verification$stability$delta[[1L]] + 1e-3
      x
    },
    central_order = function(x) {
      x$computation$orders$M_verification_required <- 159L
      x
    },
    selected_id = function(x) {
      x$computation$selected_candidate_id <- "candidate-forged"
      x
    },
    candidate_objective = function(x) {
      x$computation$candidate_evaluations[[1L]]$fresh_objective <- 0
      x
    },
    method = function(x) {
      x$method <- "A2-MN+NM"
      x
    },
    provenance = function(x) {
      x$provenance$selected_method <- "A2-MN+NM"
      x
    },
    alias = function(x) {
      x$a <- x$a + 1
      x
    },
    custom_class = function(x) {
      x$parameters$a <- structure(x$parameters$a, class = "forged")
      x
    },
    duplicate_name = function(x) {
      names(x)[match("message", names(x))] <- "status"
      x
    }
  )

  for (name in names(mutations)) {
    condition <- .a2_expect_schema_rejection(mutations[[name]](fit))
    expect_false(
      inherits(condition, "simpleError") &&
        !inherits(condition, "dpprior_schema_error"),
      label = name
    )
  }
})


test_that("canonical finite and failed objects round-trip serialization", {
  fit <- DPprior_a2_newton(50, 5, 8)
  round_trip <- unserialize(serialize(fit, NULL, version = 3L))
  expect_identical(round_trip, fit)
  .a2_expect_valid_result(round_trip, "converged")

  condition <- tryCatch(
    DPprior_a2_newton(
      50, 5, 8, a0 = fit$a, b0 = fit$b,
      M = 257, max_iter = 1, tol_F = 1
    ),
    error = identity
  )
  failed_round_trip <- unserialize(
    serialize(condition$result, NULL, version = 3L)
  )
  expect_identical(failed_round_trip, condition$result)
  .a2_expect_valid_result(failed_round_trip, "failed")
  expect_null(failed_round_trip$parameters)
  expect_identical(
    failed_round_trip$compatibility$views$legacy_v2$required_action,
    "refit_with_smaller_M_at_or_below_256"
  )
})


test_that("A2-MN validates targets and controls with typed conditions", {
  expect_error(DPprior_a2_newton(1, 1.1, 0.1))
  expect_error(
    DPprior_a2_newton(50, 0.5, 8), class = "dpprior_moment_target_error"
  )
  expect_error(
    DPprior_a2_newton(50, 50, 1), class = "dpprior_moment_target_error"
  )
  expect_error(
    DPprior_a2_newton(50, 5, -1), class = "dpprior_moment_target_error"
  )
  expect_error(DPprior_a2_newton(10, 9, 9), "maximum possible variance")
  expect_error(
    DPprior_a2_newton(50, 5, 8, tol_F = -1),
    class = "dpprior_control_error"
  )
  expect_error(
    DPprior_a2_newton(50, 5, 8, tol_rel = -1),
    class = "dpprior_control_error"
  )
  expect_error(
    DPprior_a2_newton(50, 5, 8, tol_step = 0),
    class = "dpprior_control_error"
  )
  expect_error(
    DPprior_a2_newton(50, 5, 8, max_iter = "x"),
    class = "dpprior_control_error"
  )
  expect_error(
    DPprior_a2_newton(50, 5, 8, damping = NA),
    class = "dpprior_control_error"
  )
  expect_error(
    DPprior_a2_newton(
      50, 5, 8,
      verification_abs_tol = 0, verification_rel_tol = 0
    ),
    class = "dpprior_control_error"
  )
})


test_that("verification helper and S3 print consume the canonical result", {
  expect_true(
    verify_a2_moment_matching(50, 5, 8, tol = 1e-6, verbose = FALSE)
  )
  fit <- DPprior_a2_newton(50, 5, 8)
  expect_output(print(fit), "DPprior.*Elicitation Result")
  expect_output(print(fit), "A2-MN")
  expect_output(print(fit), "Gamma")
})


test_that("A2 consumers ignore forged non-authoritative legacy views", {
  producer <- DPprior_a2_newton
  forge_legacy_view <- function(fit) {
    forged <- fit
    legacy <- forged$compatibility$views$legacy_v2
    legacy$source_status <- "failed"
    legacy$source_usable <- FALSE
    legacy$source_verified <- FALSE
    legacy$numerical_candidate <- list(a = 999, b = 999, J = 50L)
    legacy$selected_fit <- list(
      mu_K = 500, var_K = 800, residual = 0, scaled_residual = 0,
      residual_components = c(mean = 0, variance = 0)
    )
    legacy$iterations <- 999L
    legacy$termination <- "max_iter"
    legacy$solver_diagnostics$a0 <- 999
    forged$compatibility$views$legacy_v2 <- legacy

    # These deprecated top-level fields remain valid aliases of the forged
    # legacy view.  Scientific consumers must nevertheless ignore them.
    forged$fit <- legacy$selected_fit
    forged$iterations <- legacy$iterations
    forged$termination <- legacy$termination
    forged$diagnostics <- legacy$solver_diagnostics
    .dpprior_validate_result_v1(forged)
  }

  baseline_comparison <- compare_a1_a2(50, 5, 8, verbose = FALSE)
  baseline_single <- verify_a2_moment_matching(
    50, 5, 8, verbose = FALSE
  )
  baseline_all <- verify_a2_all(verbose = FALSE)
  testthat::local_mocked_bindings(
    DPprior_a2_newton = function(...) forge_legacy_view(producer(...)),
    .package = "DPprior"
  )

  expect_identical(
    compare_a1_a2(50, 5, 8, verbose = FALSE), baseline_comparison
  )
  expect_identical(
    verify_a2_moment_matching(50, 5, 8, verbose = FALSE),
    baseline_single
  )
  expect_identical(verify_a2_all(verbose = FALSE), baseline_all)
  all_output <- capture.output(verify_a2_all(verbose = TRUE))
  expect_false(any(grepl("999|max_iter", all_output, fixed = FALSE)))
  expect_true(any(grepl("term=converged", all_output, fixed = TRUE)))
  expect_output(
    verify_a2_moment_matching(50, 5, 8, verbose = TRUE),
    "Termination: converged \\(optimizer\\)"
  )
})


test_that("A2 consumers reject malformed canonical evidence with typed errors", {
  malformed <- DPprior_a2_newton(50, 5, 8)
  malformed$achieved$K$mean <- malformed$achieved$K$mean + 1e-3
  testthat::local_mocked_bindings(
    DPprior_a2_newton = function(...) malformed,
    .package = "DPprior"
  )

  for (consumer in list(
    function() verify_a2_moment_matching(50, 5, 8, verbose = FALSE),
    function() compare_a1_a2(50, 5, 8, verbose = FALSE),
    function() verify_a2_all(verbose = FALSE)
  )) {
    condition <- tryCatch(consumer(), error = identity)
    expect_s3_class(condition, "dpprior_schema_error")
    expect_true(nzchar(condition$code))
    expect_true(nzchar(condition$path))
    expect_false(is.null(condition$expected))
    expect_false(
      inherits(condition, "simpleError") &&
        !inherits(condition, "dpprior_schema_error")
    )
  }
})
