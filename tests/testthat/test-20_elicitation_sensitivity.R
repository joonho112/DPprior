.r20_cache <- new.env(parent = emptyenv())


.r20_cached <- function(name, factory) {
  if (!exists(name, envir = .r20_cache, inherits = FALSE)) {
    assign(name, factory(), envir = .r20_cache)
  }
  get(name, envir = .r20_cache, inherits = FALSE)
}


.r20_raw <- function(x) unclass(x)


.r20_key_for_route <- function(x, route) {
  raw <- .r20_raw(x)
  evidence <- raw$sensitivity$fit_evidence
  routes <- vapply(evidence, function(one) {
    one$input_provenance$target_route
  }, character(1L))
  names(evidence)[[match(route, routes)]]
}


.r20_actual_grid <- function() {
  .r20_cached("actual_grid", function() {
    generated <- pmf_K_marginal(
      20L, 2, 3, compute_log_stirling(20L),
      M = 80L, M_verify = 160L, strict = FALSE
    )
    target_pmf <- unname(as.numeric(generated[-1L]))
    .dpprior_run_elicitation_sensitivity(
      J = 20L,
      scenarios = list(
        list(
          mu_K = 5, var_K = 8, method = "A2-MN", M = 80L,
          check_diagnostics = TRUE, label = "direct"
        ),
        list(
          mu_K = 5, confidence = "high", method = "A2-MN", M = 80L,
          check_diagnostics = FALSE
        ),
        list(
          mu_K = 5, cv_K = 0.5, method = "A2-MN", M = 80L,
          check_diagnostics = TRUE
        ),
        list(
          K_interval = list(
            lower = 3L, upper = 10L, type = "equal_tail",
            coverage = 0.8, family = "maxent"
          ),
          method = "A2-KL", M = 80L, check_diagnostics = FALSE
        ),
        list(
          target_pmf = target_pmf, method = "A2-KL", M = 80L,
          check_diagnostics = TRUE
        )
      ),
      check_diagnostics = FALSE
    )
  })
}


.r20_direct <- function() {
  .r20_cached("direct", function() {
    .dpprior_run_elicitation_sensitivity(
      20L,
      list(list(mu_K = 5, var_K = 8, method = "A2-MN", M = 80L)),
      check_diagnostics = TRUE
    )
  })
}


.r20_reviewer_interval_grid <- function() {
  .r20_cached("reviewer_interval_grid", function() {
    .dpprior_run_elicitation_sensitivity(
      20L,
      lapply(c(0.7, 0.8, 0.9), function(coverage) list(
        K_interval = list(
          lower = 3L, upper = 10L, type = "equal_tail",
          coverage = coverage, family = "maxent"
        ),
        method = "A2-KL", M = 80L, check_diagnostics = FALSE
      ))
    )
  })
}


.r20_expect_schema_error <- function(x) {
  expect_error(
    .dpprior_validate_result_v1(x),
    class = "dpprior_schema_error"
  )
}


test_that("actual producers cover all five canonical elicitation routes", {
  fit <- .r20_actual_grid()
  expect_s3_class(fit, "dpprior_elicitation_sensitivity")
  expect_s3_class(fit, "dpprior_result")
  expect_invisible(.dpprior_validate_result_v1(fit))
  expect_invisible(.dpprior_validate_object(fit))

  raw <- .r20_raw(fit)
  extension <- raw$sensitivity
  expect_identical(
    names(extension),
    c(
      "scenarios", "scenario_results", "fit_evidence", "conditions",
      "interval_checks", "metrics_long", "local", "global", "metadata"
    )
  )
  expect_identical(nrow(extension$scenarios), 5L)
  expect_identical(
    sort(unname(vapply(extension$fit_evidence, function(one) {
      one$input_provenance$target_route
    }, character(1L)))),
    sort(c(
      "direct_variance", "qualitative_confidence",
      "coefficient_of_variation", "interval", "strict_pmf"
    ))
  )
  expect_true(all(is.na(extension$scenarios$scenario_label)))
  expect_identical(
    extension$scenarios$scenario_key,
    sort(extension$scenarios$scenario_key, method = "radix")
  )
})


test_that("fit evidence retains exact route and fallback provenance", {
  fit <- .r20_actual_grid()
  raw <- .r20_raw(fit)
  evidence <- raw$sensitivity$fit_evidence

  direct <- evidence[[.r20_key_for_route(fit, "direct_variance")]]
  confidence <- evidence[[.r20_key_for_route(
    fit, "qualitative_confidence"
  )]]
  cv <- evidence[[.r20_key_for_route(
    fit, "coefficient_of_variation"
  )]]
  interval <- evidence[[.r20_key_for_route(fit, "interval")]]
  pmf <- evidence[[.r20_key_for_route(fit, "strict_pmf")]]

  expected_fields <- c(
    "request", "input_provenance", "weight_target", "target", "method",
    "status", "usable", "verified", "parameters", "evaluator",
    "selected_snapshot", "verifier_snapshot", "condition_evidence", "source"
  )
  expect_true(all(vapply(evidence, function(one) {
    identical(names(one), expected_fields)
  }, logical(1L))))
  expect_identical(
    names(direct$input_provenance),
    c(
      "method_explicit", "confidence_explicit", "requested_method",
      "selected_method", "is_fallback", "target_route"
    )
  )
  expect_identical(direct$request[c("J", "mu_K", "var_K")],
                   list(J = 20L, mu_K = 5, var_K = 8))
  expect_true(direct$input_provenance$method_explicit)
  expect_false(direct$input_provenance$confidence_explicit)
  expect_identical(direct$input_provenance$requested_method, "A2-MN")
  expect_identical(direct$input_provenance$selected_method, "A2-MN")
  expect_false(direct$input_provenance$is_fallback)
  expect_true(confidence$input_provenance$method_explicit)
  expect_true(confidence$input_provenance$confidence_explicit)
  expect_identical(
    confidence$input_provenance$target_route,
    "qualitative_confidence"
  )
  expect_equal(confidence$target$used$variance, 1.5 * (5 - 1),
               tolerance = 0)
  expect_equal(cv$target$used$variance, (0.5 * 5)^2, tolerance = 0)
  expect_identical(names(interval$request),
                   c("J", "K_interval", "mu_K", "method", "M"))
  expect_identical(names(interval$request$K_interval),
                   c(
                     "lower", "upper", "type", "coverage", "family",
                     "mu_K", "support", "endpoints"
                   ))
  expect_identical(length(pmf$target$used$pmf), 20L)
  expect_equal(sum(pmf$target$used$pmf), 1, tolerance = 1e-12)
  expect_true(all(vapply(evidence, function(one) {
    identical(one$evaluator$M_selected, 80L) &&
      identical(one$evaluator$M_verification, 160L)
  }, logical(1L))))
})


test_that("a real Newton fallback is retained in input provenance", {
  backend <- DPprior_a2_newton(
    50L, 3, 10, max_iter = 2L, M = 80L, verbose = FALSE
  )
  expect_identical(.r20_raw(backend)$method, "A2-MN+NM")
  fit <- .dpprior_run_elicitation_sensitivity(
    50L,
    list(list(mu_K = 3, var_K = 10, method = "A2-MN", M = 80L)),
    check_diagnostics = FALSE,
    .fit_fun = function(...) backend
  )
  evidence <- .r20_raw(fit)$sensitivity$fit_evidence[[1L]]
  expect_identical(evidence$input_provenance$requested_method, "A2-MN")
  expect_identical(evidence$input_provenance$selected_method, "A2-MN+NM")
  expect_true(evidence$input_provenance$is_fallback)
  expect_identical(evidence$method, "A2-MN+NM")
  expect_invisible(.dpprior_validate_result_v1(fit))
})


test_that("preflight is transactional across the entire scenario set", {
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  adapter <- function(...) {
    calls$n <- calls$n + 1L
    DPprior_fit(...)
  }
  condition <- tryCatch(
    .dpprior_run_elicitation_sensitivity(
      20L,
      list(
        list(mu_K = 5, var_K = 8, method = "A2-MN", M = 80L),
        list(mu_K = 5, var_K = 8, unsupported = TRUE)
      ),
      .fit_fun = adapter
    ),
    error = identity
  )
  expect_s3_class(condition, "dpprior_sensitivity_scenarios_error")
  expect_identical(condition$code, "unknown_scenario_fields")
  expect_identical(calls$n, 0L)

  condition <- tryCatch(
    .dpprior_run_elicitation_sensitivity(
      20L,
      list(
        list(mu_K = 5, var_K = 8),
        list(mu_K = 5, confidence = "high", cv_K = 0.5)
      ),
      .fit_fun = adapter
    ),
    error = identity
  )
  expect_identical(condition$code, "multiple_uncertainty_sources")
  expect_identical(calls$n, 0L)
})


test_that("condition result dispatch distinguishes fit target and NULL", {
  fit_result <- .dpprior_run_elicitation_sensitivity(
    20L,
    list(list(mu_K = 5, var_K = 8, method = "A2-MN", M = 20L)),
    check_diagnostics = FALSE
  )
  fit_raw <- .r20_raw(fit_result)
  fit_evidence <- fit_raw$sensitivity$fit_evidence[[1L]]
  expect_identical(fit_evidence$status, "approximate")
  expect_false(fit_evidence$usable)
  expect_false(fit_evidence$verified)
  expect_false(is.null(fit_evidence$parameters))
  expect_identical(
    fit_evidence$condition_evidence$calibration$code,
    "calibration_unusable"
  )

  target_result <- .dpprior_run_elicitation_sensitivity(
    10L,
    list(list(
      K_interval = list(
        lower = 3L, upper = 10L, type = "equal_tail",
        coverage = 0.8, family = "maxent"
      ),
      M = 80L
    ))
  )
  target_raw <- .r20_raw(target_result)
  target_evidence <- target_raw$sensitivity$fit_evidence[[1L]]
  expect_identical(target_evidence$status, "infeasible")
  expect_null(target_evidence$parameters)
  expect_false(target_evidence$input_provenance$method_explicit)
  expect_false(target_evidence$input_provenance$confidence_explicit)
  expect_identical(target_evidence$input_provenance$requested_method,
                   "A2-KL")
  expect_identical(target_evidence$input_provenance$selected_method,
                   "A2-KL")
  expect_false(target_evidence$input_provenance$is_fallback)
  expect_identical(target_evidence$input_provenance$target_route, "interval")
  expect_identical(
    target_evidence$condition_evidence$calibration,
    target_evidence$condition_evidence$target
  )
  expect_identical(
    target_raw$sensitivity$interval_checks[[1L]]$source,
    "target_infeasibility_certificate"
  )

  null_result <- .dpprior_run_elicitation_sensitivity(
    20L,
    list(list(mu_K = 5, var_K = 8, method = "A2-MN", M = 80L)),
    check_diagnostics = FALSE,
    .fit_fun = function(...) stop("injected NULL-result failure", call. = FALSE)
  )
  null_raw <- .r20_raw(null_result)
  null_evidence <- null_raw$sensitivity$fit_evidence[[1L]]
  expect_identical(null_evidence$status, "failed")
  expect_null(null_evidence$parameters)
  expect_identical(
    null_evidence$condition_evidence$calibration$code,
    "sensitivity_backend_contract"
  )
  expect_invisible(.dpprior_validate_result_v1(null_result))
})


test_that("A1 diagnostics approximation stays in the diagnostics channel", {
  fit <- .dpprior_run_elicitation_sensitivity(
    20L,
    list(list(
      mu_K = 5, var_K = 8, method = "A1", M = 10L,
      check_diagnostics = TRUE
    ))
  )
  raw <- .r20_raw(fit)
  evidence <- raw$sensitivity$fit_evidence[[1L]]
  condition <- evidence$condition_evidence$diagnostics
  expect_identical(evidence$status, "approximate")
  expect_true(evidence$usable)
  expect_false(evidence$verified)
  expect_identical(condition$class,
                   "dpprior_diagnostics_approximation_error")
  expect_identical(
    condition$classes,
    c(
      "dpprior_diagnostics_approximation_error", "dpprior_fit_error",
      "dpprior_calibration_error", "dpprior_error", "error",
      "dpprior_condition", "condition"
    )
  )
  expect_identical(condition$code, "fit_diagnostics_approximate")
  expect_null(evidence$condition_evidence$calibration)
  expect_identical(evidence$source, "retained_canonical_fit_evidence")
  expect_invisible(.dpprior_validate_result_v1(fit))
})


test_that("J one strict PMF retains certified nonidentifiability", {
  fit <- .dpprior_run_elicitation_sensitivity(
    1L, list(list(target_pmf = 1, method = "A2-KL", M = 20L))
  )
  raw <- .r20_raw(fit)
  evidence <- raw$sensitivity$fit_evidence[[1L]]
  expect_identical(raw$status, "infeasible")
  expect_false(raw$usable)
  expect_true(raw$verified)
  expect_identical(
    evidence$condition_evidence$calibration$code,
    "calibration_nonidentifiable_j1"
  )
  expect_null(evidence$condition_evidence$target)
  expect_false(raw$sensitivity$scenarios$diagnostics_requested[[1L]])
})


test_that("only fit-attached diagnostics are accepted", {
  expect_false(".diagnostics_fun" %in%
                 names(formals(.dpprior_run_elicitation_sensitivity)))
  calls <- new.env(parent = emptyenv())
  calls$flags <- logical()
  adapter <- function(...) {
    args <- list(...)
    calls$flags <- c(calls$flags, args$check_diagnostics)
    do.call(DPprior_fit, args)
  }
  fit <- .dpprior_run_elicitation_sensitivity(
    20L,
    list(
      list(
        mu_K = 5, var_K = 8, method = "A2-MN", M = 80L,
        check_diagnostics = TRUE
      ),
      list(
        mu_K = 5, var_K = 8.5, method = "A2-MN", M = 80L,
        check_diagnostics = FALSE
      )
    ),
    .fit_fun = adapter
  )
  expect_setequal(calls$flags, c(TRUE, FALSE))
  raw <- .r20_raw(fit)
  for (index in seq_len(nrow(raw$sensitivity$scenarios))) {
    requested <- raw$sensitivity$scenarios$diagnostics_requested[[index]]
    value <- raw$sensitivity$scenario_results$E_W_SB[[index]]
    expect_identical(is.finite(value), requested)
  }

  unattached <- DPprior_fit(
    20L, 5, 8, method = "A2-MN", M = 80L,
    check_diagnostics = FALSE
  )
  quarantined <- .dpprior_run_elicitation_sensitivity(
    20L,
    list(list(mu_K = 5, var_K = 8, method = "A2-MN", M = 80L)),
    check_diagnostics = TRUE, .fit_fun = function(...) unattached
  )
  quarantine_raw <- .r20_raw(quarantined)
  expect_identical(quarantine_raw$status, "failed")
  expect_identical(
    quarantine_raw$sensitivity$conditions[[1L]]$diagnostics$code,
    "sensitivity_diagnostic_contract"
  )
})


test_that("fresh selected and verifier snapshots bind every finite fit", {
  fit <- .r20_actual_grid()
  evidence <- .r20_raw(fit)$sensitivity$fit_evidence
  finite <- evidence[!vapply(evidence, function(one) {
    is.null(one$parameters)
  }, logical(1L))]
  expect_true(length(finite) > 0L)
  for (one in finite) {
    expect_identical(one$selected_snapshot$M, one$request$M)
    expect_identical(
      one$verifier_snapshot$M,
      max(2L * one$request$M, one$request$M + 40L)
    )
    expect_equal(sum(one$selected_snapshot$K$pmf), 1, tolerance = 1e-10)
    expect_equal(sum(one$verifier_snapshot$K$pmf), 1, tolerance = 1e-10)
    support <- seq_len(one$request$J)
    selected_mean <- sum(support * one$selected_snapshot$K$pmf)
    selected_variance <- sum(
      (support - selected_mean)^2 * one$selected_snapshot$K$pmf
    )
    expect_equal(one$selected_snapshot$K$mean, selected_mean,
                 tolerance = 1e-13)
    expect_equal(one$selected_snapshot$K$variance, selected_variance,
                 tolerance = 1e-13)
  }
})


test_that("the exact 22 metrics and four-field global summary reconcile", {
  fit <- .r20_actual_grid()
  raw <- .r20_raw(fit)
  sensitivity <- raw$sensitivity
  expect_identical(
    names(sensitivity$scenario_results),
    c("scenario_key", "status", "usable", "verified",
      .DPPRIOR_SENSITIVITY_METRICS)
  )
  expect_identical(
    unique(sensitivity$metrics_long$metric),
    .DPPRIOR_SENSITIVITY_METRICS
  )
  expect_identical(
    nrow(sensitivity$metrics_long),
    5L * length(.DPPRIOR_SENSITIVITY_METRICS)
  )
  expect_identical(
    names(sensitivity$global),
    c("scenario_count", "converged_count", "failed_count", "metric_count")
  )
  expect_identical(sensitivity$global$scenario_count, 5L)
  expect_identical(
    sensitivity$global$metric_count,
    as.integer(nrow(sensitivity$metrics_long))
  )
  expect_identical(names(sensitivity$local), .DPPRIOR_SENSITIVITY_LOCAL_FIELDS)
  expect_identical(nrow(sensitivity$local), 0L)
  unavailable <- is.na(sensitivity$metrics_long$value)
  expect_true(all(!is.na(sensitivity$metrics_long$reason[unavailable])))
  expect_true(all(nzchar(sensitivity$metrics_long$reason[unavailable])))
  expect_true(all(is.na(
    sensitivity$metrics_long$reason[!unavailable]
  )))
})


test_that("reviewer [3,10] A2-KL cells preserve signed negative residuals", {
  fit <- .r20_reviewer_interval_grid()
  raw <- .r20_raw(fit)
  rows <- raw$sensitivity$scenario_results
  expect_true(all(rows$interval_residual < 0))
  expect_equal(
    rows$interval_residual,
    rows$interval_achieved - rows$interval_requested,
    tolerance = 0
  )
  expect_true(all(rows$status == "approximate"))
  for (index in seq_len(nrow(rows))) {
    key <- rows$scenario_key[[index]]
    audit <- raw$sensitivity$interval_checks[[key]]
    expect_equal(audit$selected$coverage_residual,
                 rows$interval_residual[[index]], tolerance = 0)
    expect_false(audit$verification$passed)
  }
  expect_invisible(.dpprior_validate_result_v1(fit))
})


test_that("scientific identity is order independent and labels are quarantined", {
  first <- list(
    list(
      mu_K = 4, var_K = 7, method = "A1", M = 80L,
      check_diagnostics = FALSE, label = "first label"
    ),
    list(
      M = 80L, method = "A1", var_K = 9, mu_K = 6,
      check_diagnostics = FALSE, scenario_label = "second label"
    )
  )
  second <- list(
    list(
      scenario_label = "changed", check_diagnostics = FALSE,
      mu_K = 6L, var_K = 9L, M = 80, method = "A1"
    ),
    list(
      label = "also changed", method = "A1", var_K = 7L,
      mu_K = 4L, check_diagnostics = FALSE, M = 80
    )
  )
  x <- .dpprior_run_elicitation_sensitivity(20L, first)
  y <- .dpprior_run_elicitation_sensitivity(20, second)
  expect_identical(x, y)
  expect_true(all(is.na(.r20_raw(x)$sensitivity$scenarios$scenario_label)))

  duplicate <- tryCatch(
    .dpprior_run_elicitation_sensitivity(
      20L,
      list(
        list(mu_K = 4, var_K = 7, method = "A1", label = "a"),
        list(mu_K = 4L, var_K = 7L, method = "A1", label = "b")
      )
    ),
    error = identity
  )
  expect_s3_class(duplicate, "dpprior_sensitivity_duplicate_scenario")
  expect_identical(duplicate$code, "duplicate_scenario_content")
})


test_that("malformed M and adapter contracts stop before calibration", {
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  adapter <- function(...) {
    calls$n <- calls$n + 1L
    DPprior_fit(...)
  }
  for (bad_M in list(9L, 257L, 20.5, NA_real_, structure(80, class = "bad"))) {
    condition <- tryCatch(
      .dpprior_run_elicitation_sensitivity(
        20L, list(list(mu_K = 5, var_K = 8)),
        M = bad_M, .fit_fun = adapter
      ),
      error = identity
    )
    expect_s3_class(condition, "dpprior_sensitivity_scenarios_error")
  }
  expect_identical(calls$n, 0L)

  condition <- tryCatch(
    .dpprior_run_elicitation_sensitivity(
      20L, list(list(mu_K = 5, var_K = 8)),
      .fit_fun = function(J) NULL
    ),
    error = identity
  )
  expect_s3_class(condition, "dpprior_sensitivity_adapter_error")
  expect_identical(condition$code, "fit_adapter_formals")

  condition <- tryCatch(
    .dpprior_run_elicitation_sensitivity(
      20L, list(list(mu_K = 5, var_K = 8)),
      .fit_fun = adapter,
      .key_fun = function(...) "scn_0000000000000000"
    ),
    error = identity
  )
  expect_s3_class(condition, "dpprior_sensitivity_adapter_error")
  expect_identical(condition$code, "scenario_key_noncanonical")
  expect_identical(calls$n, 0L)
})


test_that("the upper supported M uses an independent order of 512", {
  fit <- .dpprior_run_elicitation_sensitivity(
    20L,
    list(list(
      mu_K = 5, var_K = 8, method = "A1", M = 256L,
      check_diagnostics = FALSE
    ))
  )
  evidence <- .r20_raw(fit)$sensitivity$fit_evidence[[1L]]
  expect_identical(evidence$evaluator$M_selected, 256L)
  expect_identical(evidence$evaluator$M_verification, 512L)
  expect_identical(evidence$verifier_snapshot$M, 512L)
  expect_invisible(.dpprior_validate_result_v1(fit))
})


test_that("canonical fit mismatches fail closed without borrowing science", {
  wrong <- DPprior_fit(
    20L, 6, 9, method = "A2-MN", M = 80L,
    check_diagnostics = FALSE
  )
  fit <- .dpprior_run_elicitation_sensitivity(
    20L,
    list(list(mu_K = 5, var_K = 8, method = "A2-MN", M = 80L)),
    check_diagnostics = FALSE, .fit_fun = function(...) wrong
  )
  raw <- .r20_raw(fit)
  evidence <- raw$sensitivity$fit_evidence[[1L]]
  expect_identical(raw$status, "failed")
  expect_null(evidence$parameters)
  expect_identical(
    evidence$condition_evidence$calibration$code,
    "sensitivity_backend_contract"
  )
  expect_false(identical(evidence$target, .r20_raw(wrong)$target$K))
})


test_that("hostile fit, target, and condition accessors never dispatch", {
  fit_backend <- DPprior_fit(
    20L, 5, 8, method = "A2-MN", M = 80L,
    check_diagnostics = FALSE
  )
  counter <- new.env(parent = emptyenv())
  counter$dollar <- 0L
  counter$bracket <- 0L
  counter$target_dollar <- 0L
  counter$target_bracket <- 0L
  counter$condition_bracket <- 0L
  global <- .GlobalEnv
  method_names <- c(
    "$.DPprior_fit", "[[.DPprior_fit", "$.DPprior_target_K",
    "[[.DPprior_target_K", "[[.dpprior_diagnostics_approximation_error"
  )
  previous <- lapply(method_names, function(name) {
    get0(name, envir = global, inherits = FALSE)
  })
  on.exit({
    for (index in seq_along(method_names)) {
      name <- method_names[[index]]
      if (is.null(previous[[index]])) {
        if (exists(name, envir = global, inherits = FALSE)) {
          rm(list = name, envir = global)
        }
      } else {
        assign(name, previous[[index]], envir = global)
      }
    }
  }, add = TRUE)
  assign(method_names[[1L]], function(x, name) {
    counter$dollar <- counter$dollar + 1L
    stop("hostile $ dispatched", call. = FALSE)
  }, envir = global)
  assign(method_names[[2L]], function(x, i, ...) {
    counter$bracket <- counter$bracket + 1L
    stop("hostile [[ dispatched", call. = FALSE)
  }, envir = global)
  assign(method_names[[3L]], function(x, name) {
    counter$target_dollar <- counter$target_dollar + 1L
    stop("hostile target $ dispatched", call. = FALSE)
  }, envir = global)
  assign(method_names[[4L]], function(x, i, ...) {
    counter$target_bracket <- counter$target_bracket + 1L
    stop("hostile target [[ dispatched", call. = FALSE)
  }, envir = global)
  assign(method_names[[5L]], function(x, i, ...) {
    counter$condition_bracket <- counter$condition_bracket + 1L
    stop("hostile condition [[ dispatched", call. = FALSE)
  }, envir = global)

  fit <- .dpprior_run_elicitation_sensitivity(
    20L,
    list(list(mu_K = 5, var_K = 8, method = "A2-MN", M = 80L)),
    check_diagnostics = FALSE,
    .fit_fun = function(...) fit_backend
  )
  condition_fit <- .dpprior_run_elicitation_sensitivity(
    20L,
    list(list(
      mu_K = 5, var_K = 8, method = "A1", M = 10L,
      check_diagnostics = TRUE
    ))
  )
  expect_identical(counter$dollar, 0L)
  expect_identical(counter$bracket, 0L)
  expect_identical(counter$target_dollar, 0L)
  expect_identical(counter$target_bracket, 0L)
  expect_identical(counter$condition_bracket, 0L)
  expect_invisible(.dpprior_validate_result_v1(fit))
  expect_invisible(.dpprior_validate_result_v1(condition_fit))
})


test_that("serialization preserves canonical evidence byte for byte", {
  fit <- .r20_actual_grid()
  copy <- unserialize(serialize(fit, NULL, version = 3L))
  expect_identical(copy, fit)
  expect_invisible(.dpprior_validate_result_v1(copy))
  expect_identical(
    serialize(copy, NULL, version = 3L),
    serialize(fit, NULL, version = 3L)
  )
})


test_that("authority mutations are rejected independently", {
  fit <- .r20_direct()
  mutate <- function(path) {
    candidate <- unserialize(serialize(fit, NULL))
    path(candidate)
  }
  mutations <- list(
    request = function(x) {
      key <- names(x$sensitivity$fit_evidence)[[1L]]
      x$sensitivity$fit_evidence[[key]]$request$var_K <- 8.1
      x
    },
    input_provenance = function(x) {
      key <- names(x$sensitivity$fit_evidence)[[1L]]
      x$sensitivity$fit_evidence[[key]]$input_provenance$target_route <-
        "qualitative_confidence"
      x
    },
    target = function(x) {
      key <- names(x$sensitivity$fit_evidence)[[1L]]
      x$sensitivity$fit_evidence[[key]]$target$used$var_K <- 9
      x
    },
    parameter = function(x) {
      key <- names(x$sensitivity$fit_evidence)[[1L]]
      x$sensitivity$fit_evidence[[key]]$parameters$a <-
        x$sensitivity$fit_evidence[[key]]$parameters$a + 0.1
      x
    },
    selected_pmf = function(x) {
      key <- names(x$sensitivity$fit_evidence)[[1L]]
      pmf <- x$sensitivity$fit_evidence[[key]]$selected_snapshot$K$pmf
      pmf[1:2] <- rev(pmf[1:2])
      x$sensitivity$fit_evidence[[key]]$selected_snapshot$K$pmf <- pmf
      x
    },
    verifier_mean = function(x) {
      key <- names(x$sensitivity$fit_evidence)[[1L]]
      x$sensitivity$fit_evidence[[key]]$verifier_snapshot$K$mean <-
        x$sensitivity$fit_evidence[[key]]$verifier_snapshot$K$mean + 0.01
      x
    },
    scenario_flag = function(x) {
      x$sensitivity$scenarios$method_explicit[[1L]] <- FALSE
      x
    },
    wide_metric = function(x) {
      x$sensitivity$scenario_results$E_alpha[[1L]] <-
        x$sensitivity$scenario_results$E_alpha[[1L]] + 0.01
      x
    },
    long_metric = function(x) {
      row <- which(x$sensitivity$metrics_long$metric == "E_alpha")[[1L]]
      x$sensitivity$metrics_long$value[[row]] <-
        x$sensitivity$metrics_long$value[[row]] + 0.01
      x
    },
    conditions = function(x) {
      key <- names(x$sensitivity$conditions)[[1L]]
      x$sensitivity$conditions[[key]]$calibration <- list(
        class = "dpprior_calibration_unusable",
        classes = c(
          "dpprior_calibration_unusable", "dpprior_calibration_error",
          "dpprior_error", "error", "dpprior_condition", "condition"
        ),
        code = "calibration_unusable", message = "forged"
      )
      x
    },
    global = function(x) {
      x$sensitivity$global$metric_count <-
        x$sensitivity$global$metric_count - 1L
      x
    },
    reconciliation = function(x) {
      x$verification$settings$scenario_count <- 2L
      x
    }
  )
  for (name in names(mutations)) {
    candidate <- mutate(mutations[[name]])
    .r20_expect_schema_error(candidate)
  }
})


test_that("interval audit and signed metric mutations are rejected", {
  fit <- .r20_reviewer_interval_grid()
  raw <- .r20_raw(fit)
  key <- raw$sensitivity$scenarios$scenario_key[[1L]]

  changed_audit <- unserialize(serialize(fit, NULL))
  changed_audit$sensitivity$interval_checks[[key]]$selected$coverage <-
    changed_audit$sensitivity$interval_checks[[key]]$selected$coverage + 0.01
  .r20_expect_schema_error(changed_audit)

  changed_sign <- unserialize(serialize(fit, NULL))
  index <- match(key, changed_sign$sensitivity$scenario_results$scenario_key)
  changed_sign$sensitivity$scenario_results$interval_residual[[index]] <-
    abs(changed_sign$sensitivity$scenario_results$interval_residual[[index]])
  .r20_expect_schema_error(changed_sign)
})


test_that("R20 contains no legacy target or standalone diagnostics science", {
  source_candidates <- c(
    testthat::test_path("..", "..", "R",
                        "20_elicitation_sensitivity.R"),
    testthat::test_path("..", "..", "00_pkg_src", "DPprior", "R",
                        "20_elicitation_sensitivity.R")
  )
  source_path <- source_candidates[file.exists(source_candidates)][1L]
  testthat::skip_if(
    is.na(source_path),
    "Package source is unavailable for this static contract audit."
  )
  source <- paste(readLines(source_path, warn = FALSE), collapse = "\n")
  expect_false(grepl("target_v0", source, fixed = TRUE))
  expect_false(grepl("compatibility$views", source, fixed = TRUE))
  expect_false(grepl("DPprior_diagnostics(", source, fixed = TRUE))
})
