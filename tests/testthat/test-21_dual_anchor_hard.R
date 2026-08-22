.hard21_fast_control <- function() {
  list(
    scan_points = 7L, scan_keep = 5L, profile_starts = 2L, maxit = 40L
  )
}


.hard21_osm_control <- function() {
  list(
    scan_points = 9L, scan_keep = 5L, profile_starts = 2L, maxit = 100L
  )
}


.hard21_constraint <- function(bound = 0.25, relation = "<=") {
  list(
    metric = "wsb_tail", threshold = 0.5,
    relation = relation, bound = bound
  )
}


.hard21_catch <- function(expr) {
  tryCatch(expr, error = function(condition) condition)
}


.hard21_result <- function(x) {
  if (inherits(x, "condition") && !is.null(x$result)) x$result else x
}


.hard21_clone <- function(x) {
  unserialize(serialize(x, NULL, version = 3L))
}


.hard21_expect_schema_error <- function(x) {
  condition <- .hard21_catch(.dpprior_validate_result_v1(x))
  expect_s3_class(condition, "dpprior_schema_error")
  expect_true(is.character(condition$code) && length(condition$code) == 1L &&
                nzchar(condition$code))
  expect_true(is.character(condition$path) && length(condition$path) == 1L &&
                nzchar(condition$path))
  expect_false(
    inherits(condition, "simpleError") &&
      !inherits(condition, "dpprior_schema_error")
  )
  condition
}


.hard21_expect_canonical <- function(x, status = NULL) {
  expect_s3_class(x, "DPprior_dual_hard")
  expect_s3_class(x, "DPprior_fit")
  expect_s3_class(x, "dpprior_result")
  expect_identical(names(x)[seq_len(18L)], c(
    "schema", "object_type", "mode", "method", "J", "status", "usable",
    "verified", "message", "parameters", "target", "achieved",
    "residuals", "tolerances", "computation", "verification",
    "provenance", "compatibility"
  ))
  expect_identical(names(x)[[19L]], "constraint")
  expect_identical(x$schema$name, "dpprior.result")
  expect_identical(x$schema$version, 1L)
  expect_identical(x$object_type, "fit")
  expect_identical(x$mode, "dual_hard")
  expect_identical(x$method, "dual_anchor_hard_inequality")
  if (!is.null(status)) expect_identical(x$status, status)
  expect_invisible(.dpprior_validate_result_v1(x))
  x
}


.hard21_fit <- function(J = 20L, mu_K = 5, var_K = 8) {
  DPprior_fit(
    J, mu_K, var_K, method = "A2-MN", M = 80L,
    check_diagnostics = FALSE
  )
}


.hard21_osm_fit <- local({
  value <- NULL
  function() {
    if (is.null(value)) {
      value <<- DPprior_fit(
        50L, 3, 10, method = "A2-MN", M = 80L,
        check_diagnostics = FALSE
      )
    }
    .hard21_clone(value)
  }
})


.hard21_osm_result <- local({
  value <- NULL
  function() {
    if (is.null(value)) {
      value <<- DPprior_dual_hard(
        .hard21_osm_fit(), .hard21_constraint(),
        M = 20L, M_verify = 60L,
        control = .hard21_osm_control()
      )
    }
    .hard21_clone(value)
  }
})


test_that("hard inputs are exact, typed, and mode separated", {
  fit <- .hard21_fit()
  constraint <- .hard21_constraint(0.4)

  lambda_error <- .hard21_catch(DPprior_dual_hard(
    fit, constraint, lambda = stop("must not evaluate")
  ))
  expect_s3_class(lambda_error, "dpprior_dual_lambda_error")
  expect_s3_class(lambda_error, "dpprior_invalid_input")
  expect_identical(lambda_error$code, "lambda_in_hard_mode")

  invalid <- list(
    .hard21_catch(DPprior_dual_hard(
      fit, c(constraint, value = 0.4)
    )),
    .hard21_catch(DPprior_dual_hard(
      fit, constraint, M = 80L, M_verify = 159L
    )),
    .hard21_catch(DPprior_dual_hard(
      fit, constraint, control = list(mxit = 10L)
    )),
    .hard21_catch(DPprior_dual_hard(
      fit, constraint, log_bounds = matrix(c(-2, 2), nrow = 1L)
    )),
    .hard21_catch(DPprior_dual_hard(
      fit,
      list(metric = "unknown", relation = "<=", bound = 0.5)
    ))
  )
  for (condition in invalid) {
    expect_s3_class(condition, "dpprior_invalid_input")
    expect_true(is.character(condition$code) && nzchar(condition$code))
    expect_false(
      inherits(condition, "simpleError") &&
        !inherits(condition, "dpprior_invalid_input")
    )
  }

  a1 <- DPprior_fit(
    20L, 4, 8, method = "A1", check_diagnostics = FALSE
  )
  a1_error <- .hard21_catch(DPprior_dual_hard(a1, constraint))
  expect_s3_class(a1_error, "dpprior_dual_fit_error")
  expect_identical(a1_error$code, "input_fit_not_decision_ready")

  approximate_A2 <- .hard21_result(.hard21_catch(DPprior_fit(
    20L, 4, 8, method = "A2-KL", M = 80L,
    check_diagnostics = FALSE
  )))
  expect_identical(approximate_A2$status, "approximate")
  expect_false(approximate_A2$usable)
  expect_false(approximate_A2$verified)
  A2_error <- .hard21_catch(DPprior_dual_hard(
    approximate_A2, constraint
  ))
  expect_s3_class(A2_error, "dpprior_dual_fit_error")
  expect_identical(A2_error$code, "input_fit_not_decision_ready")
})


test_that("canonical input normalization binds orders, targets, and A2-KL evidence", {
  moment <- DPprior_a2_newton(
    50L, 5, 8, M = 80L, M_verify = 200L
  )
  moment_normalized <- .dpprior_v2_normalize_fit(moment)
  moment_reference <- moment_normalized$source$canonical_reference
  expect_identical(names(moment_reference), c(
    "schema", "mode", "method", "J", "status", "usable", "verified",
    "parameters", "target", "decision_evidence", "selected_snapshot",
    "verifier_snapshot"
  ))
  expect_null(moment_reference$decision_evidence)
  expect_identical(moment_reference$selected_snapshot$M, 80L)
  expect_identical(moment_reference$verifier_snapshot$M, 200L)
  expect_true(moment_normalized$source$input_fit_evidence$passed)
  expect_identical(
    moment_normalized$target_K$canonical, moment$target$K
  )

  chisq <- DPprior_a2_kl(
    50L, list(mu_K = 5, var_K = 8), method = "chisq",
    M = 80L, M_verify = 200L
  )
  expect_true(chisq$usable)
  expect_true(chisq$verified)
  kl_normalized <- .dpprior_v2_normalize_fit(chisq)
  decision <- kl_normalized$source$canonical_reference$decision_evidence
  expect_identical(names(decision), c(
    "target_K", "distribution_tolerances"
  ))
  expect_identical(decision$target_K, chisq$target$K)
  expect_identical(
    decision$distribution_tolerances, chisq$tolerances$distribution
  )
  expect_identical(
    kl_normalized$source$canonical_reference$verifier_snapshot$M, 200L
  )

  custom_bounds <- c(-5, 5)
  for (input in list(moment, chisq)) {
    result <- DPprior_dual_hard(
      input, .hard21_constraint(0.3), M = 20L, M_verify = 80L,
      log_bounds = custom_bounds, control = .hard21_fast_control(),
      allow_approximate = TRUE
    )
    .hard21_expect_canonical(result)
    expect_identical(result$target$K, input$target$K)
    expect_identical(
      result$computation$used$controls$log_bounds, custom_bounds
    )
    expect_identical(result$verification$settings$log_bounds, custom_bounds)
    expect_identical(
      result$computation$orders$M_verification_required, 60L
    )
    expect_identical(result$computation$orders$M_verification_used, 80L)
    expect_identical(
      result$provenance$input_fit$verifier_snapshot$M, 200L
    )
  }
})


test_that("all five K target routes retain byte-identical authority", {
  J_pmf <- 8L
  M_pmf <- 40L
  target_pmf <- pmf_K_marginal(
    J_pmf, 2, 1, compute_log_stirling(J_pmf),
    M = M_pmf, M_verify = 80L, strict = FALSE
  )[-1L]
  routes <- list(
    direct = DPprior_fit(
      20L, 5, 8, method = "A2-MN", M = 80L,
      check_diagnostics = FALSE
    ),
    confidence = DPprior_fit(
      20L, 5, method = "A2-MN", M = 80L,
      check_diagnostics = FALSE
    ),
    cv = DPprior_fit(
      20L, 5, cv_K = sqrt(8) / 5, method = "A2-MN", M = 80L,
      check_diagnostics = FALSE
    ),
    pmf = DPprior_fit(
      J_pmf, target_pmf = target_pmf, method = "A2-KL", M = M_pmf,
      check_diagnostics = FALSE
    ),
    interval = DPprior_fit(
      10L,
      K_interval = list(
        lower = 1L, upper = 1L, type = "hard_bounds", family = "maxent"
      ),
      M = 80L, check_diagnostics = FALSE
    )
  )
  expect_identical(
    vapply(routes, function(fit) fit$target$K$kind, character(1)),
    c(
      direct = "moments", confidence = "moments", cv = "cv",
      pmf = "pmf", interval = "interval"
    )
  )

  for (name in names(routes)) {
    fit <- routes[[name]]
    result <- DPprior_dual_hard(
      fit, .hard21_constraint(0.99), M = 80L, M_verify = 160L,
      control = .hard21_fast_control()
    )
    .hard21_expect_canonical(result)
    expect_true(result$status %in% c("converged", "boundary"), label = name)
    expect_true(result$usable, label = name)
    expect_true(result$verified, label = name)
    expect_identical(result$target$K, fit$target$K, label = name)
    expect_identical(
      result$provenance$input_fit$target,
      list(
        schema = fit$target$K$schema,
        kind = fit$target$K$kind,
        J = fit$target$K$J,
        used = fit$target$K$used,
        implied = fit$target$K$implied
      ),
      label = name
    )
  }

  poisoned <- .hard21_clone(routes$confidence)
  poisoned$target$K$compatibility$views$target_v0$request$mu_K <- -999
  poisoned$target$K$compatibility$views$target_v0$request$confidence <-
    "not_authoritative"
  clean_normalized <- .dpprior_v2_normalize_fit(routes$confidence)
  poison_normalized <- .dpprior_v2_normalize_fit(poisoned)
  expect_identical(
    poison_normalized$target_K[c("mu_K", "var_K", "units", "source")],
    clean_normalized$target_K[c("mu_K", "var_K", "units", "source")]
  )
  expect_identical(
    poison_normalized$source$canonical_reference,
    clean_normalized$source$canonical_reference
  )
})


test_that("OSM hard result preserves exact numerical identity in canonical fields", {
  fit <- .hard21_osm_fit()
  result <- .hard21_expect_canonical(.hard21_osm_result(), "boundary")

  expect_true(result$usable)
  expect_true(result$verified)
  expect_identical(result$target$K, fit$target$K)
  expect_identical(result$target$weight$metric, "wsb_tail")
  expect_identical(result$target$weight$relation, "at_most")
  expect_identical(result$target$weight$operator, "<=")
  expect_identical(result$constraint$feasibility$classification,
                   "feasible_candidate")
  expect_true(result$constraint$satisfied)
  expect_true(result$constraint$active)
  expect_lte(result$constraint$residual,
             result$constraint$tolerance$effective)
  expect_identical(result$constraint$slack, -result$constraint$residual)
  expect_equal(
    result$achieved$weight$value,
    prob_wsb_exceeds(
      0.5, result$parameters$a, result$parameters$b
    ),
    tolerance = 1e-12
  )
  expect_identical(
    sprintf("%.17g", c(
      a = result$parameters$a, b = result$parameters$b,
      K_loss = result$constraint$optimality$selected_K_loss,
      weight = result$achieved$weight$value,
      residual = result$constraint$residual
    )),
    c(
      "28.6944707486544", "14.003452300316328",
      "2.1040157620362652", "0.24999999999999939",
      "-6.106226635438361e-16"
    )
  )

  attempts <- result$computation$attempts
  evaluations <- result$computation$candidate_evaluations
  expect_true(length(attempts) > 0L)
  expect_true(length(evaluations) > 0L)
  expect_identical(
    result$constraint$feasibility$candidate_count,
    as.integer(length(evaluations))
  )
  expect_identical(
    result$constraint$feasibility$verified_candidate_count,
    as.integer(sum(vapply(
      evaluations, function(x) !is.null(x$verifier_snapshot), logical(1)
    )))
  )
  expect_identical(
    result$constraint$feasibility$feasible_candidate_count,
    as.integer(sum(vapply(
      evaluations, function(x) x$selection_eligible, logical(1)
    )))
  )
  expect_identical(
    vapply(evaluations, function(x) x$objective_kind, character(1)),
    rep("K_loss", length(evaluations))
  )
  expect_identical(
    names(result$verification$components),
    c(
      "constraint_selected", "constraint_refined", "order_stability",
      "metric_certification", "candidate_selection", "perturbation"
    )
  )
  expect_true(all(vapply(
    result$verification$components, function(x) x$passed, logical(1)
  )))
  expect_identical(
    result$parameters, result$verification$selected_snapshot$parameters
  )
  expect_identical(
    result$parameters, result$verification$verifier_snapshot$parameters
  )
  expect_false(any(c(
    "lambda", "weight_loss", "total_loss", "target_v0"
  ) %in% names(result)))
  expect_identical(result$compatibility$top_level_aliases, character())
  expect_identical(result$compatibility$views, list())
})


test_that("certificate, signed diagnostic, and empty failure tiers are honest", {
  fit <- .hard21_fit()
  common <- list(
    fit = fit, M = 20L, M_verify = 60L, log_bounds = c(-2, 2),
    control = .hard21_fast_control()
  )

  infeasible <- .hard21_catch(do.call(
    DPprior_dual_hard,
    c(common, list(constraint = .hard21_constraint(1e-8)))
  ))
  expect_s3_class(infeasible, "dpprior_dual_infeasible")
  expect_identical(infeasible$code, "dual_hard_infeasible")
  certified <- .hard21_expect_canonical(infeasible$result, "infeasible")
  expect_false(certified$usable)
  expect_true(certified$verified)
  expect_null(certified$parameters)
  expect_identical(certified$achieved, list())
  expect_identical(certified$residuals, list())
  expect_identical(
    certified$constraint$feasibility$classification,
    "certified_infeasible"
  )
  expect_true(certified$constraint$feasibility$certified_infeasible)
  expect_true(certified$constraint$feasibility$certificate$certified)
  expect_length(certified$computation$attempts, 1L)
  expect_length(certified$computation$candidate_evaluations, 0L)
  certificate <- certified$constraint$feasibility$certificate
  probe <- certified$computation$attempts[[1L]]
  expect_identical(
    probe$method, "analytic_monotonicity_feasibility_probe"
  )
  expect_identical(probe$stage, "feasibility")
  expect_identical(probe$reason_code,
                   "globally_infeasible_by_certificate")
  expect_identical(
    probe$bounds,
    list(log_a = certificate$domain$log_a,
         log_b = certificate$domain$log_b)
  )
  expect_identical(
    probe$candidate_objective, certificate$minimum$refined$value
  )
  expect_identical(
    certified$computation$termination[c("code", "source")],
    list(code = "certified_infeasible", source = "analytic_certificate")
  )

  fit_info <- .dpprior_v2_normalize_fit(fit)
  spec <- .dpprior_v2_normalize_weight_spec(
    .hard21_constraint(1e-6), "hard"
  )
  settings <- .dpprior_v2_normalize_hard_controls(
    list(abs = 1e-6, rel = 1e-6), 20L, 60L, c(-2, 2),
    .hard21_fast_control(), FALSE
  )
  fixed_diagnostic <- .dpprior_v2_hard_backend(fit_info, spec, settings)
  expect_identical(fixed_diagnostic$status, "approximate")
  local_mocked_bindings(
    .dpprior_v2_hard_backend = function(...) fixed_diagnostic,
    .package = "DPprior"
  )
  diagnostic_condition <- .hard21_catch(DPprior_dual_hard(
    fit, .hard21_constraint(1e-6), M = 20L, M_verify = 60L,
    log_bounds = c(-2, 2), control = .hard21_fast_control()
  ))
  expect_s3_class(diagnostic_condition, "dpprior_dual_approximation_error")
  expect_identical(diagnostic_condition$code, "dual_hard_approximate")
  diagnostic <- .hard21_expect_canonical(
    diagnostic_condition$result, "approximate"
  )
  expect_false(diagnostic$usable)
  expect_false(diagnostic$verified)
  expect_false(diagnostic$constraint$satisfied)
  expect_identical(
    diagnostic$constraint$feasibility$classification, "unknown"
  )
  expect_identical(
    diagnostic$computation$termination[c("code", "source")],
    list(code = "approximate", source = "candidate_evaluation")
  )
  returned_diagnostic <- DPprior_dual_hard(
    fit, .hard21_constraint(1e-6), M = 20L, M_verify = 60L,
    log_bounds = c(-2, 2), control = .hard21_fast_control(),
    allow_approximate = TRUE
  )
  expect_identical(returned_diagnostic, diagnostic_condition$result)

  empty <- .dpprior_v2_empty_hard_result(
    fit_info, spec, settings,
    list(
      classification = "unknown", certified_infeasible = FALSE,
      feasibility_unknown = TRUE, certificate = list(certified = FALSE)
    ),
    attempts = list(), status = "failed", message = "no candidate",
    candidates = list(), verified_candidates = list()
  )
  local_mocked_bindings(
    .dpprior_v2_hard_backend = function(...) empty,
    .package = "DPprior"
  )
  failed_condition <- .hard21_catch(DPprior_dual_hard(
    fit, .hard21_constraint(1e-6), M = 20L, M_verify = 60L,
    log_bounds = c(-2, 2), control = .hard21_fast_control()
  ))
  expect_s3_class(failed_condition, "dpprior_dual_hard_error")
  expect_identical(failed_condition$code, "dual_hard_failed")
  failed <- .hard21_expect_canonical(failed_condition$result, "failed")
  expect_false(failed$usable)
  expect_false(failed$verified)
  expect_null(failed$parameters)
  expect_length(failed$computation$attempts, 0L)
  expect_length(failed$computation$candidate_evaluations, 0L)
  expect_identical(
    failed$computation$termination[c("code", "source")],
    list(code = "no_candidate", source = "no_candidate")
  )

})


test_that("feasible tier outranks a lower-loss signed diagnostic", {
  result <- DPprior_dual_hard(
    .hard21_fit(), .hard21_constraint(0.001),
    M = 20L, M_verify = 60L, log_bounds = c(-2, 2),
    control = .hard21_fast_control(), allow_approximate = TRUE
  )
  .hard21_expect_canonical(result, "approximate")
  evaluations <- result$computation$candidate_evaluations
  feasible <- evaluations[vapply(
    evaluations, function(x) x$selection_eligible, logical(1)
  )]
  diagnostic <- evaluations[vapply(
    evaluations, function(x) x$diagnostic_eligible, logical(1)
  )]
  expect_true(length(feasible) > 0L)
  expect_true(length(diagnostic) > 0L)
  selected <- evaluations[[match(
    result$computation$selected_candidate_id,
    vapply(evaluations, function(x) x$id, character(1))
  )]]
  expect_true(selected$selection_eligible)
  expect_false(selected$diagnostic_eligible)
  expect_lt(
    min(vapply(diagnostic, function(x) x$selection_objective, numeric(1))),
    selected$selection_objective
  )
  expect_false(result$usable)
  expect_false(result$verified)
})


test_that("nonzero-exit finite science stays approximate and auto-caps are explicit", {
  fit <- .hard21_fit()
  fit_info <- .dpprior_v2_normalize_fit(fit)
  spec <- .dpprior_v2_normalize_weight_spec(
    .hard21_constraint(0.99), "hard"
  )
  settings <- .dpprior_v2_normalize_hard_controls(
    list(abs = 1e-6, rel = 1e-6), 80L, 160L, c(-15, 15),
    .hard21_fast_control(), FALSE
  )
  old <- .dpprior_v2_hard_backend(fit_info, spec, settings)
  retained <- old$.internal$verified_candidates[[which(vapply(
    old$.internal$verified_candidates,
    function(x) {
      identical(x$candidate$source, "K_only_L-BFGS-B") &&
        !.dpprior_v2_exit_zero(x$candidate$optimizer_exit_code) &&
        isTRUE(x$verification$passed)
    }, logical(1)
  ))[[1L]]]]
  retained$index <- 1L
  old$.internal$candidates <- list(retained$candidate)
  old$.internal$verified_candidates <- list(retained)
  old$.internal$selected <- retained
  old$.internal$diagnostic_selection <- FALSE
  old$status <- "approximate"
  old$message <- paste(
    "Finite independently verified candidate lacks successful optimizer",
    "exit evidence."
  )
  result <- .hard21_expect_canonical(
    .dpprior_v2_hard_canonical_result(old), "approximate"
  )
  selected <- result$computation$candidate_evaluations[[1L]]
  expect_false(selected$execution_success)
  expect_false(selected$optimizer_supported)
  expect_true(selected$selection_eligible)
  expect_false(selected$decision_eligible)
  expect_identical(selected$outcome, "selected_diagnostic")
  selected_attempt <- result$computation$attempts[[match(
    result$computation$selected_attempt_id,
    vapply(result$computation$attempts, function(x) x$id, character(1))
  )]]
  expect_false(.dpprior_v2_exit_zero(selected_attempt$exit_code))
  expect_identical(selected_attempt$reason_code, "selected")
  expect_false(result$usable)
  expect_false(result$verified)

  capped <- DPprior_dual_hard(
    .hard21_osm_fit(), .hard21_constraint(), M = 20L, M_verify = 60L,
    control = list(scan_points = 5L), allow_approximate = TRUE
  )
  .hard21_expect_canonical(capped)
  expect_identical(
    capped$computation$used$controls$auto_cap,
    list(scan_keep = TRUE, profile_starts = FALSE)
  )
  expect_identical(capped$computation$used$controls$scan_keep, 5L)
  expect_identical(capped$computation$used$controls$profile_starts, 4L)
})


test_that("failed unknown feasibility retains rejection evidence without a candidate", {
  condition <- .hard21_catch(DPprior_dual_hard(
    .hard21_fit(), .hard21_constraint(0.3),
    M = 20L, M_verify = 60L, log_bounds = c(-2, 2),
    control = .hard21_fast_control()
  ))
  expect_s3_class(condition, "dpprior_dual_hard_error")
  expect_identical(condition$code, "dual_hard_failed")
  result <- .hard21_expect_canonical(condition$result, "failed")
  expect_null(result$parameters)
  expect_identical(result$constraint$feasibility$classification, "unknown")
  expect_true(result$constraint$feasibility$feasibility_unknown)
  expect_true(length(result$computation$attempts) > 0L)
  expect_true(length(result$computation$candidate_evaluations) > 0L)
  expect_false(any(vapply(
    result$computation$attempts, function(x) x$selected, logical(1)
  )))
  expect_false(any(vapply(
    result$computation$candidate_evaluations,
    function(x) x$selected || x$selection_eligible || x$diagnostic_eligible,
    logical(1)
  )))
  expect_identical(
    result$computation$termination[c("code", "source")],
    list(code = "no_candidate", source = "no_candidate")
  )
})


test_that("canonical mutations fail closed with typed schema conditions", {
  base <- .hard21_osm_result()
  selected_index <- match(
    base$computation$selected_candidate_id,
    vapply(
      base$computation$candidate_evaluations,
      function(x) x$id, character(1)
    )
  )
  direct_index <- which(vapply(
    base$computation$candidate_evaluations,
    function(x) identical(x$generator, "direct_attempt") &&
      x$recorded_objective_available,
    logical(1)
  ))[[1L]]
  generated_index <- which(vapply(
    base$computation$candidate_evaluations,
    function(x) identical(x$generator, "deterministic_profile_scan"),
    logical(1)
  ))[[1L]]
  mutations <- list(
    candidate_count = function(x) {
      x$constraint$feasibility$candidate_count <-
        x$constraint$feasibility$candidate_count + 1L
      x
    },
    objective_kind = function(x) {
      x$computation$candidate_evaluations[[selected_index]]$objective_kind <-
        "penalized_diagnostic"
      x
    },
    recorded_kind = function(x) {
      x$computation$candidate_evaluations[[direct_index]]$
        recorded_objective_kind <- "penalized_diagnostic"
      x
    },
    rejection_codes = function(x) {
      x$computation$candidate_evaluations[[selected_index]]$
        rejection_codes <- "forged_rejection"
      x
    },
    generated_attempt_ownership = function(x) {
      x$computation$candidate_evaluations[[generated_index]]$attempt_id <-
        x$computation$attempts[[1L]]$id
      x
    },
    termination = function(x) {
      x$computation$termination$source <- "no_candidate"
      x
    },
    input_fit_J = function(x) {
      x$provenance$input_fit$J <- x$provenance$input_fit$J + 1L
      x
    },
    selected_id = function(x) {
      x$computation$selected_candidate_id <- "candidate-forged"
      x
    },
    verified_count = function(x) {
      x$constraint$feasibility$verified_candidate_count <-
        x$constraint$feasibility$verified_candidate_count - 1L
      x
    },
    feasible_count = function(x) {
      x$constraint$feasibility$feasible_candidate_count <-
        x$constraint$feasibility$feasible_candidate_count - 1L
      x
    },
    constraint_truth = function(x) {
      x$constraint$satisfied <- FALSE
      x
    },
    duplicate_top = function(x) {
      names(x)[[2L]] <- names(x)[[1L]]
      x
    }
  )
  for (name in names(mutations)) {
    condition <- .hard21_expect_schema_error(
      mutations[[name]](.hard21_clone(base))
    )
    expect_false(is.null(condition$expected), label = name)
  }
})


test_that("serialization and S3 access poisons cannot replace canonical truth", {
  fit <- .hard21_osm_fit()
  result <- .hard21_osm_result()
  round_trip <- unserialize(serialize(result, NULL, version = 3L))
  expect_identical(round_trip, result)
  .hard21_expect_canonical(round_trip, result$status)

  bindings <- c(
    "$.DPprior_fit", "[[.DPprior_fit",
    "$.dpprior_K_target", "[[.dpprior_K_target"
  )
  old <- lapply(bindings, get0, envir = .GlobalEnv, inherits = FALSE)
  names(old) <- bindings
  on.exit({
    for (name in bindings) {
      if (is.null(old[[name]])) {
        if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
          rm(list = name, envir = .GlobalEnv)
        }
      } else {
        assign(name, old[[name]], envir = .GlobalEnv)
      }
    }
  }, add = TRUE)
  poison <- function(...) stop("forged extraction dispatch", call. = FALSE)
  for (name in bindings) assign(name, poison, envir = .GlobalEnv)

  poisoned_result <- expect_no_error(DPprior_dual_hard(
    fit, .hard21_constraint(), M = 20L, M_verify = 60L,
    control = .hard21_osm_control()
  ))
  expect_no_error(.dpprior_validate_result_v1(poisoned_result))
  expect_no_error(capture.output(print(poisoned_result)))
  expect_no_error(summary(poisoned_result, print_output = FALSE))
  expect_no_error(as.data.frame(poisoned_result))
  expect_no_error(plot(
    poisoned_result, type = "alpha", engine = "base", show = FALSE
  ))
})
