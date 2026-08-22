.soft22_spine <- c(
  "schema", "object_type", "mode", "method", "J", "status", "usable",
  "verified", "message", "parameters", "target", "achieved", "residuals",
  "tolerances", "computation", "verification", "provenance", "compatibility"
)

.soft22_attempt_fields <- c(
  "id", "stage", "method", "start", "bounds", "control", "exit_code",
  "message", "iterations", "evaluations", "candidate_parameters",
  "candidate_objective", "elapsed_seconds", "warnings", "error", "selected",
  "reason_code", "unavailable"
)

.soft22_candidate_fields <- c(
  "id", "attempt_id", "method", "generator", "parameters",
  "objective_kind", "recorded_objective_kind", "selection_objective_kind",
  "recorded_objective", "recorded_objective_available",
  "recorded_objective_reason", "fresh_objective", "selection_objective",
  "objective_tolerance", "objective_passed", "selected_snapshot",
  "verifier_snapshot", "checks", "execution_success", "optimizer_supported",
  "scientifically_eligible", "selection_eligible", "diagnostic_eligible",
  "decision_eligible", "selected", "outcome", "rejection_codes", "source"
)

.soft22_cache <- new.env(parent = emptyenv())

.soft22_mn_fit <- function(M_verify = 160L) {
  key <- paste0("mn-", M_verify)
  if (!exists(key, envir = .soft22_cache, inherits = FALSE)) {
    assign(
      key,
      DPprior_a2_newton(
        20L, 4, 8, M = 80L, M_verify = M_verify, verbose = FALSE
      ),
      envir = .soft22_cache
    )
  }
  get(key, envir = .soft22_cache, inherits = FALSE)
}

.soft22_kl_fit <- function(M_verify = 160L) {
  key <- paste0("kl-", M_verify)
  if (!exists(key, envir = .soft22_cache, inherits = FALSE)) {
    target_pmf <- .get_K_pmf_support(
      20L, 2, 3, M = 80L, M_verify = M_verify
    )[["pmf", exact = TRUE]]
    assign(
      key,
      DPprior_a2_kl(
        20L, target_pmf, method = "pmf",
        M = 80L, M_verify = M_verify, verbose = FALSE
      ),
      envir = .soft22_cache
    )
  }
  get(key, envir = .soft22_cache, inherits = FALSE)
}

.soft22_tail_target <- function(value = 0.3, relation = "target",
                                threshold = 0.5) {
  list(
    metric = "wsb_tail", relation = relation,
    threshold = threshold, value = value
  )
}

.soft22_capture <- function(expr) {
  tryCatch(eval.parent(substitute(expr)), error = identity)
}

.soft22_retained_result <- function(x) {
  if (!inherits(x, "condition")) return(x)
  record <- if (is.list(x)) unclass(x) else NULL
  if (is.list(record)) record[["result", exact = TRUE]] else NULL
}

.soft22_expect_canonical <- function(x, status = NULL) {
  expect_identical(
    class(x),
    c("DPprior_dual_soft", "DPprior_fit", "dpprior_result", "list")
  )
  raw <- unclass(x)
  expect_identical(names(raw), c(.soft22_spine, "tradeoff"))
  expect_identical(raw[["object_type", exact = TRUE]], "fit")
  expect_identical(raw[["mode", exact = TRUE]], "dual_soft")
  expect_identical(raw[["method", exact = TRUE]], "dual-soft")
  if (!is.null(status)) {
    expect_identical(raw[["status", exact = TRUE]], status)
  }
  expect_silent(.dpprior_validate_result_v1(x))
  invisible(raw)
}

.soft22_scientific_core <- function(x) {
  raw <- unclass(x)
  raw[c(
    "parameters", "target", "achieved", "residuals", "tolerances",
    "verification", "tradeoff"
  )]
}


test_that("lambda one is an exact canonical no-optimizer endpoint", {
  never_call <- function(...) stop("optimizer must not run")
  for (fit in list(.soft22_mn_fit(200L), .soft22_kl_fit(200L))) {
    input_raw <- unclass(fit)
    input_info <- .dpprior_v2_normalize_fit(fit)
    result <- DPprior_dual_soft(
      fit, .soft22_tail_target(), lambda = 1,
      M_fit = 80L, M_verify = 200L, log_bounds = c(-10, 10),
      control = list(.optim_fun = never_call)
    )
    raw <- .soft22_expect_canonical(result, "converged")
    expect_true(raw[["usable", exact = TRUE]])
    expect_true(raw[["verified", exact = TRUE]])
    expect_identical(
      raw[["parameters", exact = TRUE]],
      input_raw[["parameters", exact = TRUE]]
    )
    target <- unclass(raw[["target", exact = TRUE]])
    expect_identical(
      target[["K", exact = TRUE]],
      unclass(input_raw[["target", exact = TRUE]])[["K", exact = TRUE]]
    )
    computation <- unclass(raw[["computation", exact = TRUE]])
    expect_length(computation[["attempts", exact = TRUE]], 0L)
    expect_length(computation[["candidate_evaluations", exact = TRUE]], 0L)
    expect_null(computation[["selected_attempt_id", exact = TRUE]])
    expect_null(computation[["selected_candidate_id", exact = TRUE]])
    termination <- unclass(computation[["termination", exact = TRUE]])
    expect_identical(termination[["code", exact = TRUE]], "endpoint")
    expect_identical(termination[["source", exact = TRUE]], "endpoint")
    verification <- unclass(raw[["verification", exact = TRUE]])
    expect_identical(
      names(verification[["components", exact = TRUE]]),
      c("endpoint_input_identity", "order_stability")
    )
    tradeoff <- unclass(raw[["tradeoff", exact = TRUE]])
    expect_true(tradeoff[["endpoint", exact = TRUE]])
    optimality <- unclass(tradeoff[["optimality", exact = TRUE]])
    expect_false(optimality[["performed", exact = TRUE]])
    expect_false(optimality[["passed", exact = TRUE]])
    provenance <- unclass(raw[["provenance", exact = TRUE]])
    expect_identical(
      provenance[["input_fit", exact = TRUE]],
      input_info[["source", exact = TRUE]][[
        "canonical_reference", exact = TRUE
      ]]
    )
  }
})


test_that("soft public inputs fail closed before scientific work", {
  fit <- .soft22_mn_fit()
  expect_error(
    DPprior_dual_soft(fit, .soft22_tail_target()),
    class = "dpprior_dual_soft_invalid_input"
  )
  for (bad in list(
    0, -0.1, 1.1, NA_real_, Inf, c(0.5, 0.7), matrix(0.5),
    structure(0.5, class = "evil")
  )) {
    expect_error(
      DPprior_dual_soft(fit, .soft22_tail_target(), lambda = bad),
      class = "dpprior_invalid_input"
    )
  }
  tiny <- .soft22_capture(DPprior_dual_soft(
    fit, .soft22_tail_target(), lambda = 1e-320
  ))
  expect_s3_class(tiny, "dpprior_invalid_input")
  expect_identical(
    unclass(tiny)[["code", exact = TRUE]],
    "dual_soft_lambda_not_representable"
  )
  forced <- FALSE
  unknown <- .soft22_capture(DPprior_dual_soft(
    fit, .soft22_tail_target(), lambda = 1,
    constraint_satisfied = {
      forced <<- TRUE
      stop("must not be forced")
    }
  ))
  expect_s3_class(unknown, "dpprior_invalid_input")
  expect_false(forced)
  expect_error(
    DPprior_dual_soft(
      fit, .soft22_tail_target(), lambda = 1, M_fit = 257L
    ),
    class = "dpprior_invalid_input"
  )
  expect_error(
    DPprior_dual_soft(
      fit, .soft22_tail_target(), lambda = 1,
      M_fit = 256L, M_verify = 513L
    ),
    class = "dpprior_invalid_input"
  )
  expect_error(
    DPprior_dual_soft(
      fit,
      list(
        metric = "wmax_tail", relation = "target",
        threshold = 0.5, value = 0.3
      ),
      lambda = 1
    ),
    class = "dpprior_invalid_input"
  )
})


test_that("interior result binds selected loss, ledger, and refined KKT evidence", {
  fit <- .soft22_mn_fit()
  result <- DPprior_dual_soft(
    fit, .soft22_tail_target(), lambda = 0.7,
    M_fit = 80L, M_verify = 160L
  )
  raw <- .soft22_expect_canonical(result, "converged")
  expect_true(raw[["usable", exact = TRUE]])
  expect_true(raw[["verified", exact = TRUE]])
  computation <- unclass(raw[["computation", exact = TRUE]])
  attempts <- computation[["attempts", exact = TRUE]]
  evaluations <- computation[["candidate_evaluations", exact = TRUE]]
  expect_true(length(attempts) >= 1L)
  expect_true(length(evaluations) >= 1L)
  for (attempt in attempts) {
    expect_identical(names(attempt), .soft22_attempt_fields)
  }
  for (evaluation in evaluations) {
    expect_identical(names(evaluation), .soft22_candidate_fields)
  }
  attempt_ids <- vapply(
    attempts,
    function(x) unclass(x)[["id", exact = TRUE]],
    character(1)
  )
  evaluation_ids <- vapply(
    evaluations,
    function(x) unclass(x)[["id", exact = TRUE]],
    character(1)
  )
  expect_identical(
    attempt_ids, sprintf("attempt-%03d", seq_along(attempts))
  )
  expect_identical(
    evaluation_ids, sprintf("candidate-%03d", seq_along(evaluations))
  )
  for (evaluation in evaluations) {
    evaluation <- unclass(evaluation)
    owner <- unclass(attempts[[match(
      evaluation[["attempt_id", exact = TRUE]], attempt_ids
    )]])
    expect_identical(evaluation[["generator", exact = TRUE]], "direct_attempt")
    expect_identical(
      evaluation[["objective_kind", exact = TRUE]], "soft_tradeoff"
    )
    expect_identical(
      evaluation[["recorded_objective_kind", exact = TRUE]],
      "soft_tradeoff"
    )
    expect_identical(
      evaluation[["selection_objective_kind", exact = TRUE]],
      "soft_tradeoff"
    )
    expect_identical(
      evaluation[["parameters", exact = TRUE]],
      owner[["candidate_parameters", exact = TRUE]]
    )
    expect_identical(
      evaluation[["recorded_objective", exact = TRUE]],
      owner[["candidate_objective", exact = TRUE]]
    )
  }
  selected_id <- computation[["selected_candidate_id", exact = TRUE]]
  selected_index <- match(
    selected_id,
    vapply(
      evaluations,
      function(x) unclass(x)[["id", exact = TRUE]],
      character(1)
    )
  )
  selected <- unclass(evaluations[[selected_index]])
  tradeoff <- unclass(raw[["tradeoff", exact = TRUE]])
  expect_identical(
    selected[["selection_objective", exact = TRUE]],
    tradeoff[["total_loss", exact = TRUE]]
  )
  expect_null(selected[["verifier_snapshot", exact = TRUE]])
  optimality <- unclass(tradeoff[["optimality", exact = TRUE]])
  expect_identical(
    optimality[["candidate_objective", exact = TRUE]],
    tradeoff[["total_loss", exact = TRUE]]
  )
  expect_length(optimality[["neighbor_objectives", exact = TRUE]], 8L)
  expect_identical(
    names(optimality[["neighbor_objectives", exact = TRUE]]),
    names(optimality[["neighbor_tolerances", exact = TRUE]])
  )
  expect_identical(
    names(optimality[["gradient", exact = TRUE]]),
    c("log_shape", "log_rate")
  )
  expect_identical(
    names(optimality[["component_pass", exact = TRUE]]),
    c("log_shape", "log_rate")
  )
  verification <- unclass(raw[["verification", exact = TRUE]])
  expect_identical(
    names(verification[["components", exact = TRUE]]),
    c(
      "objective_recomputation", "order_stability", "local_optimality",
      "candidate_selection"
    )
  )
  expect_identical(
    unclass(verification[["selected_snapshot", exact = TRUE]])[[
      "M", exact = TRUE
    ]],
    80L
  )
  expect_identical(
    unclass(verification[["verifier_snapshot", exact = TRUE]])[[
      "M", exact = TRUE
    ]],
    160L
  )
  fit_info <- .dpprior_v2_normalize_fit(fit)
  target_info <- .dpprior_soft_normalize_target(.soft22_tail_target())
  parameters <- unclass(raw[["parameters", exact = TRUE]])
  eta <- log(c(
    parameters[["a", exact = TRUE]], parameters[["b", exact = TRUE]]
  ))
  scales <- c(
    mean = max(abs(fit_info$target_K$mu_K), 1),
    variance = max(abs(fit_info$target_K$var_K), 1)
  )
  selected_objective <- .dpprior_soft_objective(
    fit_info, target_info, 0.7, 80L, scales, .LOG_BOUNDS_DEFAULT
  )
  verifier_objective <- .dpprior_soft_objective(
    fit_info, target_info, 0.7, 160L, scales, .LOG_BOUNDS_DEFAULT
  )
  expect_identical(
    selected[["selection_objective", exact = TRUE]], selected_objective(eta)
  )
  expect_identical(
    optimality[["recomputed_objective", exact = TRUE]],
    selected_objective(eta)
  )
  expect_identical(
    optimality[["local_base_objective", exact = TRUE]],
    verifier_objective(eta)
  )
  expect_identical(
    optimality[["start_objective", exact = TRUE]],
    verifier_objective(log(c(fit_info$a, fit_info$b)))
  )
  expect_identical(
    computation[["orders", exact = TRUE]][["M_selected", exact = TRUE]],
    80L
  )
  expect_identical(
    computation[["orders", exact = TRUE]][[
      "M_verification_used", exact = TRUE
    ]],
    160L
  )
})


test_that("actual A2-MN and strict-PMF A2-KL cover four weight metrics", {
  fits <- list(MN = .soft22_mn_fit(200L), KL = .soft22_kl_fit(200L))
  kl_raw <- unclass(fits[["KL"]])
  kl_target <- unclass(
    unclass(kl_raw[["target", exact = TRUE]])[["K", exact = TRUE]]
  )
  expect_identical(kl_target[["kind", exact = TRUE]], "pmf")
  targets <- list(
    tail = list(
      metric = "wsb_tail", relation = "at_most",
      threshold = 0.5, value = 0.5
    ),
    mean = list(
      metric = "wsb_mean", relation = "target", value = 0.4
    ),
    quantile = list(
      metric = "wsb_quantile", relation = "at_least",
      probability = 0.5, value = 0.2
    ),
    wmax = list(
      metric = "wmax_tail_upper", relation = "at_most",
      threshold = 0.5, value = 0.8
    )
  )
  for (fit_name in names(fits)) {
    for (target_name in names(targets)) {
      target <- targets[[target_name]]
      result <- DPprior_dual_soft(
        fits[[fit_name]], target, lambda = 0.5,
        M_fit = 80L, M_verify = 200L,
        log_bounds = c(-10, 10), allow_approximate = TRUE
      )
      raw <- .soft22_expect_canonical(result)
      expect_true(raw[["status", exact = TRUE]] %in%
                    c("converged", "boundary", "approximate"))
      if (identical(raw[["status", exact = TRUE]], "approximate")) {
        expect_false(raw[["usable", exact = TRUE]])
        expect_false(raw[["verified", exact = TRUE]])
      } else {
        expect_true(raw[["usable", exact = TRUE]])
        expect_true(raw[["verified", exact = TRUE]])
      }
      weight_target <- unclass(
        unclass(raw[["target", exact = TRUE]])[["weight", exact = TRUE]]
      )
      expect_identical(
        weight_target[["metric", exact = TRUE]],
        target[["metric", exact = TRUE]]
      )
      expect_identical(
        weight_target[["relation", exact = TRUE]],
        target[["relation", exact = TRUE]]
      )
      computation <- unclass(raw[["computation", exact = TRUE]])
      orders <- unclass(computation[["orders", exact = TRUE]])
      controls <- unclass(
        unclass(computation[["used", exact = TRUE]])[[
          "controls", exact = TRUE
        ]]
      )
      expect_identical(orders[["M_verification_used", exact = TRUE]], 200L)
      expect_identical(
        controls[["log_bounds", exact = TRUE]], c(-10, 10)
      )
    }
  }
})


test_that("boundary status is verified only with bound-aware KKT evidence", {
  fit <- .soft22_mn_fit()
  input_parameters <- unclass(unclass(fit)[["parameters", exact = TRUE]])
  eta <- log(c(
    input_parameters[["a", exact = TRUE]],
    input_parameters[["b", exact = TRUE]]
  ))
  log_bounds <- c(min(eta) - 0.02, max(eta) + 0.02)
  result <- DPprior_dual_soft(
    fit,
    list(metric = "wsb_mean", relation = "target", value = 0.9),
    lambda = 0.2, M_fit = 80L, M_verify = 160L,
    log_bounds = log_bounds, allow_approximate = TRUE
  )
  raw <- .soft22_expect_canonical(result, "boundary")
  expect_true(raw[["usable", exact = TRUE]])
  expect_true(raw[["verified", exact = TRUE]])
  computation <- unclass(raw[["computation", exact = TRUE]])
  termination <- unclass(computation[["termination", exact = TRUE]])
  expect_identical(termination[["code", exact = TRUE]], "boundary")
  expect_identical(
    termination[["boundary_reason", exact = TRUE]],
    "declared_log_parameter_boundary"
  )
  optimality <- unclass(
    unclass(raw[["tradeoff", exact = TRUE]])[["optimality", exact = TRUE]]
  )
  expect_true(any(
    optimality[["bound_state", exact = TRUE]] %in% c("lower", "upper")
  ))
  expect_true(optimality[["passed", exact = TRUE]])
})


test_that("fallback attempts and selection remain fully canonical", {
  fit <- .soft22_mn_fit()
  injected <- function(par, fn, method, lower = NULL, upper = NULL,
                       control = list(), ...) {
    if (identical(method, "L-BFGS-B")) return(list(partial = TRUE))
    stats::optim(par = par, fn = fn, method = method, control = control)
  }
  result <- DPprior_dual_soft(
    fit, .soft22_tail_target(), lambda = 0.7,
    M_fit = 80L, M_verify = 160L,
    control = list(.optim_fun = injected)
  )
  raw <- .soft22_expect_canonical(result, "converged")
  computation <- unclass(raw[["computation", exact = TRUE]])
  attempts <- computation[["attempts", exact = TRUE]]
  expect_length(attempts, 2L)
  expect_identical(
    unclass(attempts[[1L]])[["method", exact = TRUE]], "L-BFGS-B"
  )
  expect_identical(
    unclass(attempts[[2L]])[["method", exact = TRUE]], "Nelder-Mead"
  )
  expect_identical(
    computation[["selected_attempt_id", exact = TRUE]], "attempt-002"
  )
  fallback <- unclass(computation[["fallback", exact = TRUE]])
  expect_true(fallback[["attempted", exact = TRUE]])
  expect_true(fallback[["used", exact = TRUE]])
  expect_identical(
    fallback[["selected_attempt_id", exact = TRUE]], "attempt-002"
  )
})


test_that("a forged primary objective is rejected before fallback selection", {
  fit <- .soft22_mn_fit()
  counts <- stats::setNames(c(1, 0), c("function", "gradient"))
  injected <- function(par, fn, method, lower = NULL, upper = NULL,
                       control = list(), ...) {
    if (identical(method, "L-BFGS-B")) {
      return(list(
        par = par, value = 0, counts = counts,
        convergence = 0L, message = "forged primary objective"
      ))
    }
    stats::optim(par = par, fn = fn, method = method, control = control)
  }
  result <- DPprior_dual_soft(
    fit, .soft22_tail_target(), lambda = 0.7,
    M_fit = 80L, M_verify = 160L,
    control = list(.optim_fun = injected)
  )
  raw <- .soft22_expect_canonical(result, "converged")
  computation <- unclass(raw[["computation", exact = TRUE]])
  expect_identical(
    computation[["selected_attempt_id", exact = TRUE]], "attempt-002"
  )
  attempts <- computation[["attempts", exact = TRUE]]
  evaluations <- computation[["candidate_evaluations", exact = TRUE]]
  first_attempt <- unclass(attempts[[1L]])
  first_evaluation <- unclass(evaluations[[1L]])
  second_evaluation <- unclass(evaluations[[2L]])
  expect_identical(
    first_attempt[["reason_code", exact = TRUE]],
    "objective_recomputation_failed"
  )
  expect_false(first_evaluation[["objective_passed", exact = TRUE]])
  expect_false(first_evaluation[["selection_eligible", exact = TRUE]])
  expect_identical(
    first_evaluation[["rejection_codes", exact = TRUE]],
    "objective_mismatch"
  )
  expect_true(second_evaluation[["selected", exact = TRUE]])
  expect_true(second_evaluation[["decision_eligible", exact = TRUE]])
})


test_that("no candidate is a typed failed canonical result", {
  fit <- .soft22_mn_fit()
  injected <- function(...) stop("injected optimizer failure")
  caught <- .soft22_capture(DPprior_dual_soft(
    fit, .soft22_tail_target(), lambda = 0.7,
    M_fit = 80L, M_verify = 160L,
    control = list(.optim_fun = injected)
  ))
  expect_s3_class(caught, "dpprior_dual_soft_computation_error")
  expect_s3_class(caught, "dpprior_calibration_error")
  result <- .soft22_retained_result(caught)
  raw <- .soft22_expect_canonical(result, "failed")
  expect_false(raw[["usable", exact = TRUE]])
  expect_false(raw[["verified", exact = TRUE]])
  expect_null(raw[["parameters", exact = TRUE]])
  computation <- unclass(raw[["computation", exact = TRUE]])
  expect_length(computation[["attempts", exact = TRUE]], 2L)
  expect_length(computation[["candidate_evaluations", exact = TRUE]], 0L)
  expect_null(computation[["selected_attempt_id", exact = TRUE]])
  expect_null(computation[["selected_candidate_id", exact = TRUE]])
  termination <- unclass(computation[["termination", exact = TRUE]])
  expect_identical(termination[["code", exact = TRUE]], "no_candidate")
  expect_identical(termination[["source", exact = TRUE]], "no_candidate")
  verification <- unclass(raw[["verification", exact = TRUE]])
  expect_identical(verification[["method", exact = TRUE]], "no_candidate")
  expect_false(verification[["performed", exact = TRUE]])
  expect_false(verification[["passed", exact = TRUE]])
  tradeoff <- unclass(raw[["tradeoff", exact = TRUE]])
  expect_null(tradeoff[["total_loss", exact = TRUE]])
  expect_identical(
    unclass(caught)[["result", exact = TRUE]], result
  )
  caught_raw <- unclass(caught)
  expect_identical(caught_raw[["code", exact = TRUE]], "dual_soft_failed")
  expect_identical(
    caught_raw[["action", exact = TRUE]], "refit_or_change_controls"
  )
  expect_match(conditionMessage(caught), "allow_approximate does not return")

  caught_opt_in <- .soft22_capture(DPprior_dual_soft(
    fit, .soft22_tail_target(), lambda = 0.7,
    M_fit = 80L, M_verify = 160L,
    control = list(.optim_fun = injected), allow_approximate = TRUE
  ))
  expect_s3_class(caught_opt_in, "dpprior_dual_soft_computation_error")
  .soft22_expect_canonical(
    .soft22_retained_result(caught_opt_in), "failed"
  )
})


test_that("finite nonzero-exit candidates are approximate and opt-in only", {
  fit <- .soft22_mn_fit()
  counts <- stats::setNames(c(1, 0), c("function", "gradient"))
  injected <- function(par, fn, method, ...) {
    list(
      par = par, value = fn(par), counts = counts,
      convergence = 81L, message = paste("nonzero", method)
    )
  }
  caught <- .soft22_capture(DPprior_dual_soft(
    fit, .soft22_tail_target(), lambda = 0.7,
    M_fit = 80L, M_verify = 160L,
    control = list(.optim_fun = injected), allow_approximate = FALSE
  ))
  expect_s3_class(caught, "dpprior_dual_soft_approximation_error")
  expect_identical(
    unclass(caught)[["action", exact = TRUE]],
    "review_then_explicitly_allow_approximate"
  )
  retained <- .soft22_retained_result(caught)
  retained_raw <- .soft22_expect_canonical(retained, "approximate")
  expect_false(retained_raw[["usable", exact = TRUE]])
  expect_false(retained_raw[["verified", exact = TRUE]])
  returned <- DPprior_dual_soft(
    fit, .soft22_tail_target(), lambda = 0.7,
    M_fit = 80L, M_verify = 160L,
    control = list(.optim_fun = injected), allow_approximate = TRUE
  )
  returned_raw <- .soft22_expect_canonical(returned, "approximate")
  expect_false(returned_raw[["usable", exact = TRUE]])
  expect_false(returned_raw[["verified", exact = TRUE]])
  expect_identical(
    .soft22_scientific_core(returned),
    .soft22_scientific_core(retained)
  )
})


test_that("reported optimizer objective cannot manufacture a public candidate", {
  fit <- .soft22_mn_fit()
  counts <- stats::setNames(c(1, 0), c("function", "gradient"))
  forged <- function(par, fn, method, ...) {
    list(
      par = c(0, 0), value = 0, counts = counts,
      convergence = 0L, message = "forged objective"
    )
  }
  caught <- .soft22_capture(DPprior_dual_soft(
    fit, .soft22_tail_target(), lambda = 0.7,
    M_fit = 80L, M_verify = 160L,
    control = list(.optim_fun = forged)
  ))
  expect_s3_class(caught, "dpprior_dual_soft_computation_error")
  result <- .soft22_retained_result(caught)
  raw <- .soft22_expect_canonical(result, "failed")
  optimality <- unclass(
    unclass(raw[["tradeoff", exact = TRUE]])[["optimality", exact = TRUE]]
  )
  expect_false(optimality[["performed", exact = TRUE]])
  expect_false(optimality[["passed", exact = TRUE]])
  expect_match(
    optimality[["unavailable_reason", exact = TRUE]], "no finite optimizer"
  )
  expect_null(raw[["parameters", exact = TRUE]])
  expect_length(
    unclass(raw[["computation", exact = TRUE]])[[
      "candidate_evaluations", exact = TRUE
    ]],
    0L
  )
  expect_false(raw[["usable", exact = TRUE]])
  expect_false(raw[["verified", exact = TRUE]])
})


test_that("canonical input is accessor-safe and compatibility science is inert", {
  fit <- .soft22_mn_fit()
  target <- .soft22_tail_target()
  baseline <- DPprior_dual_soft(
    fit, target, 1, M_fit = 80L, M_verify = 160L
  )

  forged <- unserialize(serialize(fit, NULL))
  forged_raw <- unclass(forged)
  compatibility <- unclass(
    forged_raw[["compatibility", exact = TRUE]]
  )
  views <- unclass(compatibility[["views", exact = TRUE]])
  legacy <- unclass(views[["legacy_v2", exact = TRUE]])
  numerical <- unclass(legacy[["numerical_candidate", exact = TRUE]])
  numerical[["a"]] <- 999
  numerical[["b"]] <- 0.001
  legacy[["numerical_candidate"]] <- numerical
  views[["legacy_v2"]] <- legacy
  compatibility[["views"]] <- views
  forged_raw[["compatibility"]] <- compatibility
  class(forged_raw) <- class(forged)
  forged <- forged_raw
  expect_silent(.dpprior_validate_result_v1(forged))
  from_forged <- DPprior_dual_soft(
    forged, target, 1, M_fit = 80L, M_verify = 160L
  )
  expect_identical(
    .soft22_scientific_core(from_forged),
    .soft22_scientific_core(baseline)
  )

  local({
    calls <- 0L
    assign("$.DPprior_fit", function(...) {
      calls <<- calls + 1L
      stop("poison dollar dispatch")
    }, envir = environment())
    assign("[[.DPprior_fit", function(...) {
      calls <<- calls + 1L
      stop("poison bracket dispatch")
    }, envir = environment())
    assign("$.dpprior_K_target", function(...) {
      calls <<- calls + 1L
      stop("poison target dollar dispatch")
    }, envir = environment())
    assign("[[.dpprior_K_target", function(...) {
      calls <<- calls + 1L
      stop("poison target bracket dispatch")
    }, envir = environment())
    poisoned <- DPprior_dual_soft(
      fit, target, 1, M_fit = 80L, M_verify = 160L
    )
    .soft22_expect_canonical(poisoned, "converged")
    expect_identical(calls, 0L)
  })
})


test_that("serialization is lossless and canonical mutations fail typed", {
  result <- DPprior_dual_soft(
    .soft22_mn_fit(), .soft22_tail_target(), 0.7,
    M_fit = 80L, M_verify = 160L
  )
  roundtrip <- unserialize(serialize(result, NULL))
  expect_identical(roundtrip, result)
  .soft22_expect_canonical(roundtrip, "converged")

  loss_mutation <- unserialize(serialize(result, NULL))
  loss_raw <- unclass(loss_mutation)
  tradeoff <- unclass(loss_raw[["tradeoff", exact = TRUE]])
  tradeoff[["total_loss"]] <- tradeoff[["total_loss", exact = TRUE]] + 0.1
  loss_raw[["tradeoff"]] <- tradeoff
  class(loss_raw) <- class(result)
  expect_error(
    .dpprior_validate_result_v1(loss_raw), class = "dpprior_schema_error"
  )

  ledger_mutation <- unserialize(serialize(result, NULL))
  ledger_raw <- unclass(ledger_mutation)
  computation <- unclass(ledger_raw[["computation", exact = TRUE]])
  evaluations <- computation[["candidate_evaluations", exact = TRUE]]
  evaluation <- unclass(evaluations[[1L]])
  evaluation[["objective_kind"]] <- "K_loss"
  evaluations[[1L]] <- evaluation
  computation[["candidate_evaluations"]] <- evaluations
  ledger_raw[["computation"]] <- computation
  class(ledger_raw) <- class(result)
  expect_error(
    .dpprior_validate_result_v1(ledger_raw), class = "dpprior_schema_error"
  )

  duplicated <- unclass(result)
  names(duplicated)[[2L]] <- "schema"
  class(duplicated) <- class(result)
  expect_error(
    .dpprior_require_schema(
      duplicated, kind = "fit", allow_legacy = FALSE
    ),
    class = "dpprior_serialization_error"
  )

  extra_attribute <- result
  attr(extra_attribute, "evil") <- TRUE
  extra_condition <- .soft22_capture(
    .dpprior_validate_result_v1(extra_attribute)
  )
  expect_s3_class(extra_condition, "dpprior_schema_error")
  expect_identical(
    unclass(extra_condition)[["code", exact = TRUE]], "attributes"
  )

  owner_mutation <- unserialize(serialize(result, NULL))
  owner_raw <- unclass(owner_mutation)
  owner_computation <- unclass(owner_raw[["computation", exact = TRUE]])
  owner_evaluations <- owner_computation[[
    "candidate_evaluations", exact = TRUE
  ]]
  owner_evaluation <- unclass(owner_evaluations[[1L]])
  owner_evaluation[["attempt_id"]] <- "attempt-999"
  owner_evaluations[[1L]] <- owner_evaluation
  owner_computation[["candidate_evaluations"]] <- owner_evaluations
  owner_raw[["computation"]] <- owner_computation
  class(owner_raw) <- class(result)
  owner_condition <- .soft22_capture(
    .dpprior_validate_result_v1(owner_raw)
  )
  expect_s3_class(owner_condition, "dpprior_schema_error")
  expect_identical(
    unclass(owner_condition)[["code", exact = TRUE]], "candidate_attempt_id"
  )

  evil_fit <- .soft22_mn_fit()
  class(evil_fit) <- c("evil", class(evil_fit))
  expect_error(
    DPprior_dual_soft(
      evil_fit, .soft22_tail_target(), 1,
      M_fit = 80L, M_verify = 160L
    ),
    class = "dpprior_schema_error"
  )
  for (bad_fit in list(
    new.env(parent = emptyenv()),
    pairlist(schema = "dpprior.result/1")
  )) {
    expect_error(
      DPprior_dual_soft(
        bad_fit, .soft22_tail_target(), 1,
        M_fit = 80L, M_verify = 160L
      ),
      class = "dpprior_serialization_error"
    )
  }
})


test_that("unverified A1 and interval A2-KL inputs cannot enter soft decisions", {
  a1 <- DPprior_fit(
    20L, 4, 8, method = "A1", check_diagnostics = FALSE
  )
  a1_raw <- unclass(a1)
  expect_identical(a1_raw[["status", exact = TRUE]], "approximate")
  expect_false(a1_raw[["verified", exact = TRUE]])
  expect_error(
    DPprior_dual_soft(a1, .soft22_tail_target(), 1),
    class = "dpprior_invalid_input"
  )

  interval_condition <- .soft22_capture(DPprior_fit(
    20L, mu_K = 4,
    K_interval = list(
      lower = 2L, upper = 6L, type = "central_mass",
      coverage = 0.8, family = "maxent"
    ),
    M = 80L, check_diagnostics = FALSE
  ))
  interval_fit <- .soft22_retained_result(interval_condition)
  expect_false(is.null(interval_fit))
  interval_raw <- unclass(interval_fit)
  interval_target <- unclass(
    unclass(interval_raw[["target", exact = TRUE]])[["K", exact = TRUE]]
  )
  expect_identical(interval_target[["kind", exact = TRUE]], "interval")
  expect_false(interval_raw[["usable", exact = TRUE]])
  expect_false(interval_raw[["verified", exact = TRUE]])
  expect_error(
    DPprior_dual_soft(interval_fit, .soft22_tail_target(), 1),
    class = "dpprior_invalid_input"
  )
})


test_that("trade-off curve is canonical, order-independent, and warm-started", {
  target <- .soft22_tail_target()
  first <- compute_tradeoff_curve(
    20L, list(mu_K = 4, var_K = 8),
    lambda_seq = c(0.4, 1, 0.7), M = 80L, M_verify = 200L,
    log_bounds = c(-10, 10), target = target
  )
  second <- compute_tradeoff_curve(
    20L, list(mu_K = 4, var_K = 8),
    lambda_seq = c(0.7, 0.4, 1), M = 80L, M_verify = 200L,
    log_bounds = c(-10, 10), target = target
  )
  expect_s3_class(first, "dpprior_tradeoff_curve")
  expect_identical(first[["lambda"]], c(0.4, 0.7, 1))
  expect_identical(first[["mode"]], rep("dual_soft", 3L))
  expect_identical(first[["point_id"]], second[["point_id"]])
  expect_equal(
    first[["total_loss"]], second[["total_loss"]], tolerance = 1e-12
  )
  expect_identical(
    attr(first, "metadata")[["evaluation_order", exact = TRUE]],
    c(1, 0.7, 0.4)
  )
  expect_identical(
    first[["warm_start_from"]][first[["lambda"]] == 0.4],
    first[["point_id"]][first[["lambda"]] == 0.7]
  )
  fits <- attr(first, "fits")
  conditions <- attr(first, "conditions")
  expect_length(fits, nrow(first))
  expect_length(conditions, nrow(first))
  expect_identical(names(fits), first[["point_id"]])
  expect_identical(names(conditions), first[["point_id"]])
  expect_identical(names(fits), names(conditions))
  for (fit in fits) .soft22_expect_canonical(fit)
  expect_true(all(vapply(conditions, is.null, logical(1))))
})


test_that("curve condition slots preserve mixed success and failure identity", {
  injected <- function(fit, target, lambda, ...) {
    if (identical(lambda, 0.7)) stop("injected lambda 0.7 failure")
    DPprior_dual_soft(fit, target, lambda, ...)
  }
  curve <- compute_tradeoff_curve(
    20L, list(mu_K = 4, var_K = 8),
    lambda_seq = c(0.4, 1, 0.7), M = 80L, M_verify = 160L,
    target = .soft22_tail_target(),
    control = list(.fit_fun = injected)
  )
  fits <- attr(curve, "fits")
  conditions <- attr(curve, "conditions")
  expect_identical(curve[["lambda"]], c(0.4, 0.7, 1))
  expect_length(fits, nrow(curve))
  expect_length(conditions, nrow(curve))
  expect_identical(names(fits), curve[["point_id"]])
  expect_identical(names(conditions), curve[["point_id"]])
  expect_identical(names(fits), names(conditions))
  failure_index <- which(curve[["lambda"]] == 0.7)
  success_indices <- setdiff(seq_len(nrow(curve)), failure_index)
  expect_true(all(vapply(
    conditions[success_indices], is.null, logical(1)
  )))
  condition <- conditions[[failure_index]]
  expect_s3_class(condition, "dpprior_dual_soft_backend_contract_error")
  expect_identical(
    unclass(condition)[["result", exact = TRUE]], fits[[failure_index]]
  )
  .soft22_expect_canonical(fits[[failure_index]], "failed")
  for (index in success_indices) {
    .soft22_expect_canonical(fits[[index]])
  }
})


test_that("curve backend errors retain only validated canonical failed results", {
  curve <- compute_tradeoff_curve(
    20L, list(mu_K = 4, var_K = 8),
    lambda_seq = c(0.4, 1, 0.7), M = 80L, M_verify = 160L,
    target = .soft22_tail_target(),
    control = list(.fit_fun = function(...) stop("boom"))
  )
  expect_true(all(curve[["status"]] == "failed"))
  expect_false(any(curve[["usable"]]))
  expect_false(any(curve[["verified"]]))
  expect_true(all(
    curve[["condition_code"]] == "dual_soft_backend_contract"
  ))
  fits <- attr(curve, "fits")
  conditions <- attr(curve, "conditions")
  expect_length(fits, nrow(curve))
  expect_length(conditions, nrow(curve))
  expect_identical(names(fits), curve[["point_id"]])
  expect_identical(names(conditions), curve[["point_id"]])
  expect_identical(names(fits), names(conditions))
  for (index in seq_along(fits)) {
    .soft22_expect_canonical(fits[[index]], "failed")
    condition <- conditions[[index]]
    expect_s3_class(
      condition, "dpprior_dual_soft_backend_contract_error"
    )
    expect_identical(
      unclass(condition)[["result", exact = TRUE]], fits[[index]]
    )
  }
})


test_that("curve rejects forged canonical science and bypasses S3 accessors", {
  forged_fit <- function(fit, target, lambda, ...) {
    point <- DPprior_dual_soft(fit, target, lambda, ...)
    raw <- unclass(point)
    tradeoff <- unclass(raw[["tradeoff", exact = TRUE]])
    tradeoff[["total_loss"]] <- tradeoff[["total_loss", exact = TRUE]] + 1
    raw[["tradeoff"]] <- tradeoff
    class(raw) <- class(point)
    raw
  }
  forged_curve <- compute_tradeoff_curve(
    20L, list(mu_K = 4, var_K = 8), lambda_seq = 1,
    M = 80L, M_verify = 160L, target = .soft22_tail_target(),
    control = list(.fit_fun = forged_fit)
  )
  expect_identical(forged_curve[["status"]], "failed")
  expect_identical(
    forged_curve[["condition_code"]], "dual_soft_backend_contract"
  )
  .soft22_expect_canonical(attr(forged_curve, "fits")[[1L]], "failed")

  local({
    calls <- 0L
    assign("$.DPprior_dual_soft", function(...) {
      calls <<- calls + 1L
      stop("poison curve dollar dispatch")
    }, envir = environment())
    assign("[[.DPprior_dual_soft", function(...) {
      calls <<- calls + 1L
      stop("poison curve bracket dispatch")
    }, envir = environment())
    genuine <- compute_tradeoff_curve(
      20L, list(mu_K = 4, var_K = 8),
      lambda_seq = c(0.7, 1), M = 80L, M_verify = 160L,
      target = .soft22_tail_target()
    )
    expect_true(all(genuine[["status"]] == "converged"))
    expect_identical(calls, 0L)
  })
})


test_that("curve public controls and lambda policy remain fail closed", {
  common <- list(
    J = 20L, K_target = list(mu_K = 4, var_K = 8),
    target = .soft22_tail_target(), M = 80L
  )
  expect_error(
    do.call(
      compute_tradeoff_curve,
      c(common, list(lambda_seq = c(0, 1)))
    ),
    class = "dpprior_invalid_input"
  )
  expect_error(
    do.call(
      compute_tradeoff_curve,
      c(common, list(lambda_seq = c(0.5, 0.5)))
    ),
    class = "dpprior_invalid_input"
  )
  expect_error(
    do.call(
      compute_tradeoff_curve,
      c(common, list(lambda_seq = 0.5, loss_type = "adaptive"))
    ),
    class = "dpprior_invalid_input"
  )
  expect_error(
    compute_tradeoff_curve(
      20L, list(mu_Kish = 4, var_K = 8),
      lambda_seq = 1, M = 80L, target = .soft22_tail_target()
    ),
    class = "dpprior_invalid_input"
  )
  expect_error(
    compute_tradeoff_curve(
      20L, list(mu_K = 4, var_K = 8),
      lambda_seq = 1, M = 80L, M_verify = 513L,
      target = .soft22_tail_target()
    ),
    class = "dpprior_invalid_input"
  )
  legacy <- compute_tradeoff_curve(
    20L, list(mu_K = 4, var_K = 8),
    w1_target = list(
      prob = list(threshold = 0.5, value = 0.3)
    ),
    lambda_seq = c(0.7, 1), M = 80L
  )
  expect_identical(legacy[["mode"]], rep("dual_soft", 2L))
  expect_identical(legacy[["metric"]], rep("wsb_tail", 2L))
  expect_false(any(legacy[["lambda"]] == 0))
})
