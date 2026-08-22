.diag14_raw <- function(x) unclass(x)

.diag14_condition <- function(expr) {
  tryCatch(
    withCallingHandlers(
      force(expr),
      warning = function(warning) invokeRestart("muffleWarning")
    ),
    error = function(error) error
  )
}

.diag14_fit <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <<- DPprior_fit(
        20L, 4, 8, method = "A1", check_diagnostics = FALSE
      )
    }
    cached
  }
})

.diag14_result <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) cached <<- DPprior_diagnostics(.diag14_fit())
    cached
  }
})

.diag14_reclass <- function(raw, template) {
  structure(raw, class = class(template))
}

.diag14_recursive_names <- function(x) {
  if (!is.list(x)) return(character())
  raw <- if (is.object(x)) unclass(x) else x
  c(names(raw), unlist(lapply(unname(raw), .diag14_recursive_names),
                       use.names = FALSE))
}


test_that("component helpers retain alpha, K, W_SB, and rho methods", {
  alpha <- compute_alpha_diagnostics(4, 2)
  expect_equal(alpha$mean, 2)
  expect_equal(alpha$cv, 0.5)
  expect_identical(alpha$method, "closed_form_and_stats_qgamma")

  K <- compute_K_diagnostics(20L, 2, 1, M = 80L, M_verify = 160L)
  expect_identical(K$estimand, "K_J")
  expect_equal(sum(K$pmf), 1, tolerance = .TOL_PMF_SUM)
  expect_identical(K$method, "gauss-laguerre-marginal-pmf-and-moments")

  weights <- compute_weight_diagnostics(
    2, 1, thresholds = c(0.5, 0.9), M = 80L, M_verify = 160L
  )
  expect_identical(weights$size_biased$estimand, "W_SB")
  expect_equal(
    weights$size_biased$tail_probability[[1L]],
    prob_wsb_exceeds(0.5, 2, 1), tolerance = 1e-14
  )
  expect_false(isTRUE(weights$maximum$requested))
  expect_null(weights$maximum$values)

  rho <- compute_coclustering_diagnostics(
    2, 1, M = 80L, M_verify = 160L
  )
  expect_identical(rho$method, "gauss-laguerre-marginal-moments")
  expect_equal(rho$mean, mean_rho(2, 1, 80L), tolerance = 1e-14)
  expect_equal(rho$var, var_rho(2, 1, 80L), tolerance = 1e-14)
})


test_that("public diagnostics returns the exact canonical contract", {
  diagnostic <- .diag14_result()
  raw <- .diag14_raw(diagnostic)

  expect_identical(
    class(diagnostic), c("DPprior_diagnostics", "dpprior_result", "list")
  )
  expect_identical(
    names(raw), c(.DPPRIOR_RESULT_COMMON_FIELDS, "diagnostics")
  )
  expect_identical(raw$schema, list(name = "dpprior.result", version = 1L))
  expect_identical(raw$object_type, "diagnostics")
  expect_identical(raw$mode, "prior_diagnostics")
  expect_identical(raw$method, "canonical_prior_diagnostics")
  expect_identical(raw$status, "converged")
  expect_true(raw$usable)
  expect_true(raw$verified)
  expect_identical(
    names(raw$diagnostics),
    c("policy_results", "warnings", "alpha", "K", "weights",
      "coclustering")
  )
  expect_identical(raw$compatibility, .dpprior_new_compatibility())
  expect_no_error(.dpprior_validate_result_v1(diagnostic))
})


test_that("public evidence is freshly recomputed at selected and verifier orders", {
  raw <- .diag14_raw(.diag14_result())
  parameters <- raw$parameters
  orders <- raw$computation$orders
  selected_M <- orders$M_selected
  verifier_M <- orders$M_verification_used

  expect_identical(selected_M, 80L)
  expect_identical(orders$M_verification_required, 160L)
  expect_identical(verifier_M, 160L)
  expect_identical(raw$achieved$alpha, list(
    mean = parameters$a / parameters$b,
    CV = 1 / sqrt(parameters$a)
  ))

  K <- .get_K_pmf_support(
    raw$J, parameters$a, parameters$b,
    M = selected_M, M_verify = verifier_M,
    abs_tol = raw$tolerances$diagnostics$absolute,
    rel_tol = raw$tolerances$diagnostics$relative
  )
  selected_pmf <- unname(as.numeric(K$pmf))
  verifier_pmf <- unname(as.numeric(K$verification_pmf))
  selected_K <- .dpprior_target_pmf_moments(selected_pmf)
  verifier_K <- .dpprior_target_pmf_moments(verifier_pmf)

  expect_identical(raw$achieved$K$pmf, selected_pmf)
  expect_equal(raw$achieved$K$mean, selected_K[["mean"]], tolerance = 1e-14)
  expect_equal(
    raw$achieved$K$variance, selected_K[["variance"]], tolerance = 1e-14
  )
  expect_identical(raw$achieved$K$estimand, "K_J")
  expect_identical(raw$achieved$K$source, "fresh_diagnostics_selected_order")
  expect_equal(
    raw$achieved$weights$mean,
    mean_w1(parameters$a, parameters$b, selected_M), tolerance = 1e-14
  )
  expect_equal(
    raw$achieved$coclustering$mean,
    mean_rho(parameters$a, parameters$b, selected_M), tolerance = 1e-14
  )
  expect_equal(
    raw$achieved$coclustering$variance,
    var_rho(parameters$a, parameters$b, selected_M), tolerance = 1e-14
  )

  verifier <- raw$verification$verifier_snapshot
  expect_identical(verifier$achieved$K$pmf, verifier_pmf)
  expect_equal(
    verifier$achieved$K$mean, verifier_K[["mean"]], tolerance = 1e-14
  )
  expect_identical(verifier$parameters, raw$parameters)
  expect_identical(
    raw$verification$selected_snapshot$parameters, raw$parameters
  )
  expect_identical(
    raw$residuals$diagnostics[["K.pmf_l1"]],
    sum(abs(selected_pmf - verifier_pmf))
  )
})


test_that("computation and verification ledgers retain no optimizer claim", {
  raw <- .diag14_raw(.diag14_result())
  attempts <- raw$computation$attempts
  expect_identical(
    vapply(attempts, `[[`, character(1), "id"),
    paste0("diagnostic-", .DPPRIOR_DIAGNOSTIC_COMPONENTS)
  )
  expect_identical(
    unname(vapply(attempts, `[[`, character(1), "method")),
    unname(.DPPRIOR_DIAGNOSTIC_ATTEMPT_METHODS)
  )
  expect_true(all(vapply(attempts, function(attempt) {
    identical(attempt$stage, "diagnostic_component") &&
      identical(attempt$exit_code, 0L) &&
      identical(attempt$iterations, 0L) &&
      !attempt$selected && is.null(attempt$candidate_parameters) &&
      is.null(attempt$candidate_objective)
  }, logical(1))))
  expect_null(raw$computation$selected_candidate_id)
  expect_null(raw$computation$selected_attempt_id)
  expect_length(raw$computation$candidate_evaluations, 0L)
  expect_identical(raw$computation$termination$code, "diagnostics_recomputed")
  expect_identical(raw$computation$termination$source, "component_aggregation")
  expect_identical(
    names(raw$verification$components), "component_aggregation"
  )
  expect_identical(
    names(raw$verification$invariants),
    c("fixed_parameters", "dominance_category_removed")
  )
})


test_that("W_SB policy is exact and its warning carries the canonical result", {
  policy <- list(
    estimand = "W_SB", direction = "above",
    weight_threshold = 0.5, action_threshold = 0.1
  )
  captured <- NULL
  diagnostic <- withCallingHandlers(
    DPprior_diagnostics(.diag14_fit(), warning_policy = policy),
    warning = function(warning) {
      captured <<- warning
      invokeRestart("muffleWarning")
    }
  )
  raw <- .diag14_raw(diagnostic)
  record <- raw$diagnostics$policy_results[[1L]]
  expect_s3_class(captured, "dpprior_diagnostic_policy_warning")
  expect_identical(captured$code, "policy_triggered")
  expect_identical(captured$result, diagnostic)
  expect_identical(record$estimand, "W_SB")
  expect_identical(record$direction, "above")
  expect_identical(record$threshold, policy$action_threshold)
  expect_equal(
    record$value,
    prob_wsb_exceeds(policy$weight_threshold, raw$parameters$a,
                     raw$parameters$b),
    tolerance = 64 * .Machine$double.eps
  )
  expect_null(record$lower)
  expect_null(record$upper)
  expect_identical(record$outcome, "triggered")
  expect_identical(record$basis, "exact_tail_probability")
  expect_length(raw$diagnostics$warnings, 1L)
  expect_no_error(.dpprior_validate_result_v1(diagnostic))
})


test_that("W_max policy is explicitly indeterminate and makes no numeric claim", {
  policy <- list(
    estimand = "W_max", direction = "above",
    weight_threshold = 0.5, action_threshold = 0.3
  )
  expect_no_warning(
    diagnostic <- DPprior_diagnostics(.diag14_fit(), warning_policy = policy)
  )
  raw <- .diag14_raw(diagnostic)
  record <- raw$diagnostics$policy_results[[1L]]
  expect_identical(record, list(
    estimand = "W_max", direction = "above", threshold = 0.3,
    value = NULL, lower = NULL, upper = NULL,
    outcome = "indeterminate", basis = "backend_unavailable"
  ))
  expect_identical(names(raw$diagnostics$weights),
                   c("status", "usable", "verified", "mean"))
  expect_false(any(c(
    "maximum", "wmax", "W_max_point", "W_max_upper", "dominance_risk"
  ) %in% .diag14_recursive_names(raw)))
  expect_length(raw$diagnostics$warnings, 0L)
  expect_no_error(.dpprior_validate_result_v1(diagnostic))
})


test_that("direct W_max point backend is outside diagnostic construction", {
  calls <- 0L
  testthat::local_mocked_bindings(
    prob_wmax_exceeds = function(...) {
      calls <<- calls + 1L
      stop("direct W_max backend must not be called")
    },
    .package = "DPprior"
  )
  policy <- list(
    estimand = "W_max", direction = "below",
    weight_threshold = 0.9, action_threshold = 0.8
  )
  diagnostic <- DPprior_diagnostics(.diag14_fit(), warning_policy = policy)
  expect_identical(calls, 0L)
  expect_identical(
    .diag14_raw(diagnostic)$diagnostics$policy_results[[1L]]$basis,
    "backend_unavailable"
  )
  helper <- compute_weight_diagnostics(
    2, 1, thresholds = c(0.5, 0.9), M = 80L, M_verify = 160L
  )
  expect_identical(calls, 0L)
  expect_false(helper$maximum$requested)
})


test_that("approximation is fail-closed and condition$result is canonical", {
  condition <- .diag14_condition(DPprior_diagnostics(
    .diag14_fit(), abs_tol = 0, rel_tol = 0
  ))
  expect_s3_class(condition, "dpprior_diagnostics_approximation_error")
  expect_identical(condition$code, "approximation_not_accepted")
  result <- condition$result
  expect_s3_class(result, "DPprior_diagnostics")
  raw <- .diag14_raw(result)
  expect_identical(raw$status, "approximate")
  expect_false(raw$usable)
  expect_false(raw$verified)
  expect_true(raw$provenance$approximation$active)
  expect_false(raw$provenance$approximation$opt_in)
  expect_no_error(.dpprior_validate_result_v1(result))

  accepted <- DPprior_diagnostics(
    .diag14_fit(), abs_tol = 0, rel_tol = 0,
    allow_approximate = TRUE
  )
  accepted_raw <- .diag14_raw(accepted)
  expect_identical(accepted_raw$status, "approximate")
  expect_false(accepted_raw$usable)
  expect_false(accepted_raw$verified)
  expect_true(accepted_raw$provenance$approximation$opt_in)
  expect_no_error(.dpprior_validate_result_v1(accepted))
})


test_that("order controls are bound and unavailable paths retain the input fit", {
  explicit <- DPprior_diagnostics(.diag14_fit(), M_verify = 200L)
  orders <- .diag14_raw(explicit)$computation$orders
  expect_identical(orders$M_selected, 80L)
  expect_identical(orders$M_verification_required, 160L)
  expect_identical(orders$M_verification_used, 200L)
  expect_identical(
    orders$verification_used_reason, "explicit_diagnostics_verifier_order"
  )
  odd_ceiling_case <- DPprior_diagnostics(.diag14_fit(), M_verify = 257L)
  expect_identical(
    .diag14_raw(odd_ceiling_case)$computation$orders$M_verification_used,
    257L
  )

  condition <- .diag14_condition(
    DPprior_diagnostics(.diag14_fit(), M_verify = 100L)
  )
  expect_s3_class(condition, "dpprior_diagnostics_order_error")
  expect_identical(condition$code, "insufficient_diagnostic_verification_order")
  expect_identical(condition$result, .diag14_fit())
  expect_no_error(.dpprior_validate_result_v1(condition$result))

  ceiling <- .diag14_condition(
    DPprior_diagnostics(.diag14_fit(), M_verify = 513L)
  )
  expect_s3_class(ceiling, "dpprior_diagnostics_order_error")
  expect_identical(ceiling$code, "diagnostic_verification_order_bounds")
  expect_identical(ceiling$result, .diag14_fit())

  selected_257 <- .diag14_raw(.diag14_fit())
  selected_257$mode <- "a2_moment"
  selected_257$computation$orders$M_selected <- 257L
  unavailable <- .diag14_condition(.diagnostic_canonical_orders(
    .diag14_fit(), selected_257, NULL
  ))
  expect_s3_class(unavailable, "dpprior_diagnostics_order_error")
  expect_identical(
    unavailable$code, "diagnostic_verification_order_unavailable"
  )
  expect_identical(unavailable$result, .diag14_fit())
})


test_that("backend failure is typed and never fabricates a partial PMF", {
  testthat::local_mocked_bindings(
    .get_K_pmf_support = function(...) stop("forced PMF backend failure"),
    .package = "DPprior"
  )
  condition <- .diag14_condition(DPprior_diagnostics(.diag14_fit()))
  expect_s3_class(condition, "dpprior_diagnostics_computation_error")
  expect_identical(condition$code, "diagnostic_recomputation_failed")
  expect_match(conditionMessage(condition), "recomputation failed")
  expect_identical(condition$result, .diag14_fit())
  expect_s3_class(condition$cause, "simpleError")
  expect_no_error(.dpprior_validate_result_v1(condition$result))
})


test_that("flat, aliased, and hostile fit inputs fail before scientific access", {
  for (bad in list(
    list(a = 2, b = 1, J = 20L),
    list(a = 2, b = 1, J = 20L, M = 80L),
    NULL,
    structure(list(a = 2, b = 1, J = 20L), class = "DPprior_fit")
  )) {
    expect_s3_class(
      .diag14_condition(DPprior_diagnostics(bad)),
      "dpprior_error"
    )
  }

  hits <- 0L
  testthat::local_mocked_s3_method(
    "$", "evil_diag_fit",
    function(x, name) {
      hits <<- hits + 1L
      stop("hostile accessor dispatched")
    }
  )
  hostile <- .diag14_fit()
  class(hostile) <- c("evil_diag_fit", class(hostile))
  expect_s3_class(
    .diag14_condition(DPprior_diagnostics(hostile)),
    "dpprior_schema_error"
  )
  expect_identical(hits, 0L)

  nested_hits <- 0L
  testthat::local_mocked_s3_method(
    "$", "evil_diag_target",
    function(x, name) {
      nested_hits <<- nested_hits + 1L
      stop("hostile target accessor dispatched")
    }
  )
  testthat::local_mocked_s3_method(
    "[[", "evil_diag_target",
    function(x, i, ..., exact = TRUE) {
      nested_hits <<- nested_hits + 1L
      stop("hostile target accessor dispatched")
    }
  )
  nested_raw <- .diag14_raw(.diag14_fit())
  class(nested_raw$target) <- "evil_diag_target"
  nested <- .diag14_reclass(nested_raw, .diag14_fit())
  expect_s3_class(
    .diag14_condition(DPprior_diagnostics(nested)), "dpprior_error"
  )
  expect_identical(nested_hits, 0L)
})


test_that("warning policy and reserved thresholds reject compatibility aliases", {
  valid <- list(
    estimand = "W_SB", direction = "above",
    weight_threshold = 0.5, action_threshold = 0.3
  )
  invalid <- list(
    list(estimand = "W_SB", threshold = 0.5,
         direction = "above", action_threshold = 0.3),
    valid[c("direction", "estimand", "weight_threshold", "action_threshold")],
    c(valid, list(extra = TRUE)),
    structure(valid, class = "policy"),
    within(valid, estimand <- "unknown"),
    within(valid, direction <- "sideways"),
    within(valid, weight_threshold <- 1),
    within(valid, action_threshold <- -0.1)
  )
  for (policy in invalid) {
    expect_s3_class(
      .diag14_condition(DPprior_diagnostics(
        .diag14_fit(), warning_policy = policy
      )),
      "dpprior_warning_policy_error"
    )
  }
  for (thresholds in list(c(0.4, 0.8), 0.5, c(0.9, 0.5))) {
    condition <- .diag14_condition(DPprior_diagnostics(
      .diag14_fit(), thresholds = thresholds
    ))
    expect_s3_class(condition, "dpprior_diagnostics_input_error")
    expect_identical(condition$code, "reserved_thresholds")
  }
})


test_that("canonical mutation matrix is rejected by the R23 truth binder", {
  diagnostic <- .diag14_result()
  mutations <- list(
    status = function(x) { x$status <- "approximate"; x },
    usable = function(x) { x$usable <- FALSE; x },
    verified = function(x) { x$verified <- FALSE; x },
    method = function(x) { x$method <- "other"; x },
    parameters = function(x) { x$parameters$a <- x$parameters$a + 0.1; x },
    target_components = function(x) {
      x$target$requested_components <- rev(x$target$requested_components); x
    },
    achieved_alpha = function(x) { x$achieved$alpha$mean <- 999; x },
    achieved_K = function(x) { x$achieved$K$mean <- 999; x },
    achieved_weight = function(x) { x$achieved$weights$mean <- 0.01; x },
    achieved_rho = function(x) {
      x$achieved$coclustering$variance <- 0.5; x
    },
    residual = function(x) {
      x$residuals$diagnostics[["K.pmf_l1"]] <- 0.5; x
    },
    tolerance = function(x) {
      x$tolerances$diagnostics$refinement[["weights.mean"]] <- 1; x
    },
    request_control = function(x) {
      x$computation$request$controls$absolute_tolerance <- 0; x
    },
    used_control = function(x) {
      x$computation$used$controls$relative_tolerance <- 0; x
    },
    selected_order = function(x) { x$computation$orders$M_selected <- 81L; x },
    attempt_id = function(x) {
      x$computation$attempts[[1L]]$id <- "diagnostic-other"; x
    },
    attempt_method = function(x) {
      x$computation$attempts[[2L]]$method <- "other"; x
    },
    attempt_reason = function(x) {
      x$computation$attempts[[3L]]$reason_code <- "component_failed"; x
    },
    termination = function(x) {
      x$computation$termination$source <- "constructor"; x
    },
    verification_method = function(x) {
      x$verification$method <- "other"; x
    },
    verification_settings = function(x) {
      x$verification$settings$M_verification <- 200L; x
    },
    selected_snapshot = function(x) {
      x$verification$selected_snapshot$achieved$K$mean <- 999; x
    },
    verifier_snapshot = function(x) {
      x$verification$verifier_snapshot$source <- "other"; x
    },
    aggregation_check = function(x) {
      x$verification$components$component_aggregation$value[["K"]] <- FALSE; x
    },
    invariant = function(x) {
      x$verification$invariants$fixed_parameters$value <- FALSE; x
    },
    extension_alpha = function(x) { x$diagnostics$alpha$mean <- 999; x },
    extension_K = function(x) { x$diagnostics$K$pmf[[1L]] <- 0.5; x },
    extension_weight = function(x) { x$diagnostics$weights$mean <- 0.1; x },
    extension_rho = function(x) {
      x$diagnostics$coclustering$mean <- 0.1; x
    },
    forbidden_category = function(x) {
      x$diagnostics$dominance_risk <- "high"; x
    }
  )

  for (name in names(mutations)) {
    raw <- mutations[[name]](.diag14_raw(diagnostic))
    candidate <- .diag14_reclass(raw, diagnostic)
    condition <- .diag14_condition(.dpprior_validate_result_v1(candidate))
    expect_s3_class(condition, "dpprior_error")
  }
})


test_that("policy evidence mutations are rejected", {
  policy <- list(
    estimand = "W_max", direction = "above",
    weight_threshold = 0.5, action_threshold = 0.3
  )
  diagnostic <- DPprior_diagnostics(.diag14_fit(), warning_policy = policy)
  changes <- list(
    value = function(x) { x$diagnostics$policy_results[[1]]$value <- 0.2; x },
    lower = function(x) { x$diagnostics$policy_results[[1]]$lower <- 0.1; x },
    outcome = function(x) {
      x$diagnostics$policy_results[[1]]$outcome <- "triggered"; x
    },
    basis = function(x) {
      x$diagnostics$policy_results[[1]]$basis <- "certified_bounds"; x
    },
    threshold = function(x) {
      x$diagnostics$policy_results[[1]]$threshold <- 0.4; x
    }
  )
  for (name in names(changes)) {
    candidate <- .diag14_reclass(
      changes[[name]](.diag14_raw(diagnostic)), diagnostic
    )
    expect_s3_class(
      .diag14_condition(.dpprior_validate_result_v1(candidate)),
      "dpprior_schema_error"
    )
  }
})


test_that("S3 methods use only the canonical gate and reject stale science", {
  diagnostic <- .diag14_result()
  summary_row <- summary(diagnostic)
  expect_identical(summary_row$schema, "dpprior.result/1")
  expect_true(summary_row$W_max_available == FALSE)
  expect_true(is.na(summary_row$P_W_max_gt_50))
  expect_true(is.na(summary_row$P_W_max_gt_90))
  expect_identical(
    as.data.frame.DPprior_diagnostics(diagnostic), summary_row
  )
  output <- capture.output(returned <- print(diagnostic))
  expect_identical(returned, diagnostic)
  expect_true(any(grepl("W_max: unavailable", output, fixed = TRUE)))
  expect_false(any(grepl("W_max >", output, fixed = TRUE)))

  stale_raw <- .diag14_raw(diagnostic)
  stale_raw$achieved$K$mean <- 999
  stale <- .diag14_reclass(stale_raw, diagnostic)
  expect_s3_class(.diag14_condition(summary(stale)), "dpprior_schema_error")
  expect_s3_class(.diag14_condition(print(stale)), "dpprior_schema_error")
  expect_s3_class(
    .diag14_condition(as.data.frame.DPprior_diagnostics(stale)),
    "dpprior_schema_error"
  )

  flat <- structure(list(
    schema_version = 2L, status = "converged", usable = TRUE,
    verified = TRUE, message = "stale", method = "old", J = 20L,
    a = 2, b = 1, K = list(mean = 999)
  ), class = "DPprior_diagnostics")
  expect_s3_class(.diag14_condition(summary(flat)), "dpprior_error")
  expect_s3_class(
    .diag14_condition(as.data.frame.DPprior_diagnostics(flat)),
    "dpprior_error"
  )
})


test_that("nested hostile accessors do not dispatch during S3 validation", {
  hits <- 0L
  testthat::local_mocked_s3_method(
    "[[", "evil_diag_nested",
    function(x, i, ..., exact = TRUE) {
      hits <<- hits + 1L
      stop("hostile nested accessor dispatched")
    }
  )
  diagnostic <- .diag14_result()
  raw <- .diag14_raw(diagnostic)
  class(raw$achieved) <- "evil_diag_nested"
  hostile <- .diag14_reclass(raw, diagnostic)
  expect_s3_class(
    .diag14_condition(summary(hostile)), "dpprior_schema_error"
  )
  expect_identical(hits, 0L)
})


test_that("canonical diagnostics survive serialization", {
  diagnostic <- .diag14_result()
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(diagnostic, path, version = 3L)
  restored <- readRDS(path)
  expect_identical(restored, diagnostic)
  expect_no_error(.dpprior_validate_result_v1(restored))
  expect_identical(summary(restored), summary(diagnostic))
})


test_that("comparison consumes canonical fits without mutating them", {
  first <- .diag14_fit()
  before <- serialize(first, NULL, version = 3L)
  comparison <- compare_diagnostics(first = first, second = first)
  expect_identical(nrow(comparison), 2L)
  expect_identical(comparison$fit, c("first", "second"))
  expect_identical(serialize(first, NULL, version = 3L), before)
  condition <- .diag14_condition(compare_diagnostics(first = first, M = 100L))
  expect_s3_class(condition, "dpprior_diagnostics_input_error")
  expect_identical(condition$code, "reserved_order")
})


test_that("diagnostics self-check exercises the canonical producer", {
  expect_true(verify_diagnostics(verbose = FALSE))
})
