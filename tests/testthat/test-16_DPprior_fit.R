.fit16_raw <- function(x) unclass(x)

.fit16_condition <- function(expr) {
  tryCatch(
    withCallingHandlers(
      force(expr),
      warning = function(warning) invokeRestart("muffleWarning")
    ),
    error = function(error) error
  )
}

.fit16_expect_canonical <- function(x, mode = NULL) {
  expect_s3_class(x, "dpprior_result")
  expect_identical(.fit16_raw(x)[["schema", exact = TRUE]][[
    "name", exact = TRUE
  ]], "dpprior.result")
  expect_identical(.fit16_raw(x)[["object_type", exact = TRUE]], "fit")
  if (!is.null(mode)) {
    expect_identical(.fit16_raw(x)[["mode", exact = TRUE]], mode)
  }
  expect_no_error(.dpprior_validate_result_v1(x))
}

.fit16_expect_unusable_result <- function(expr, mode = NULL) {
  condition <- .fit16_condition(expr)
  expect_s3_class(condition, "dpprior_calibration_unusable")
  expect_identical(condition[["code", exact = TRUE]],
                   "calibration_unusable")
  result <- condition[["result", exact = TRUE]]
  .fit16_expect_canonical(result, mode)
  raw <- .fit16_raw(result)
  expect_false(raw[["usable", exact = TRUE]])
  expect_false(raw[["verified", exact = TRUE]])
  expect_identical(condition[["status", exact = TRUE]],
                   raw[["status", exact = TRUE]])
  expect_identical(condition[["method", exact = TRUE]],
                   switch(
                     raw[["mode", exact = TRUE]],
                     a2_moment = "A2-MN", a2_kl = "A2-KL",
                     raw[["method", exact = TRUE]]
                   ))
  expect_true(
    .dpprior_is_plain_numeric(condition[["residual", exact = TRUE]]) &&
      length(condition[["residual", exact = TRUE]]) == 1L
  )
  result
}

.fit16_scientific_core <- function(x) {
  raw <- .fit16_raw(x)
  raw[c(
    "mode", "method", "J", "status", "usable", "verified", "message",
    "parameters", "achieved", "residuals", "tolerances", "verification"
  )]
}


test_that("confidence utilities retain strict ordinary-scalar behavior", {
  expect_identical(confidence_to_vif_fit("low"), 5)
  expect_identical(confidence_to_vif_fit("medium"), 2.5)
  expect_identical(confidence_to_vif_fit("high"), 1.5)
  expect_identical(vif_to_variance_fit(5, 2.5), 10)

  for (bad in list("unknown", matrix("low"), structure("low", class = "x"))) {
    expect_s3_class(
      .fit16_condition(confidence_to_vif_fit(bad)),
      "dpprior_confidence_error"
    )
  }
})


test_that("A1 wrapper preserves the native scientific result", {
  direct <- DPprior_a1(20L, 4, 8)
  wrapped <- DPprior_fit(
    20L, 4, 8, method = "A1", check_diagnostics = FALSE
  )
  .fit16_expect_canonical(wrapped, "a1_proxy")
  expect_identical(.fit16_scientific_core(wrapped),
                   .fit16_scientific_core(direct))
  raw <- .fit16_raw(wrapped)
  expect_identical(raw[["computation", exact = TRUE]][[
    "resources", exact = TRUE
  ]][["wrapper", exact = TRUE]][["target_authority", exact = TRUE]],
  "result.target.K")
})


test_that("A2-MN wrapper preserves candidate, status, and verifier evidence", {
  direct <- DPprior_a2_newton(20L, 4, 8, M = 80L, verbose = FALSE)
  wrapped <- DPprior_fit(
    20L, 4, 8, method = "A2-MN", M = 80L,
    check_diagnostics = FALSE
  )
  .fit16_expect_canonical(wrapped, "a2_moment")
  expect_identical(.fit16_scientific_core(wrapped),
                   .fit16_scientific_core(direct))
  expect_identical(
    .fit16_raw(wrapped)[["computation", exact = TRUE]][[
      "candidate_evaluations", exact = TRUE
    ]],
    .fit16_raw(direct)[["computation", exact = TRUE]][[
      "candidate_evaluations", exact = TRUE
    ]]
  )
})


test_that("A2-KL wrapper preserves objective and verifier evidence", {
  direct <- DPprior_a2_kl(
    20L, list(mu_K = 4, var_K = 8), method = "chisq",
    M = 80L, verbose = FALSE
  )
  condition <- .fit16_condition(DPprior_fit(
    20L, 4, 8, method = "A2-KL", M = 80L,
    check_diagnostics = FALSE
  ))
  expect_s3_class(condition, "dpprior_calibration_unusable")
  expect_false(inherits(condition, "dpprior_fit_error"))
  expect_identical(condition[["code", exact = TRUE]],
                   "calibration_unusable")
  expect_identical(condition[["method", exact = TRUE]], "A2-KL")
  expect_identical(condition[["status", exact = TRUE]], "approximate")
  wrapped <- condition[["result", exact = TRUE]]
  .fit16_expect_canonical(wrapped, "a2_kl")
  expect_identical(.fit16_scientific_core(wrapped),
                   .fit16_scientific_core(direct))
  expect_identical(
    condition[["residual", exact = TRUE]],
    .fit16_raw(direct)[["residuals", exact = TRUE]][[
      "distribution", exact = TRUE
    ]][["kl", exact = TRUE]]
  )
  expect_identical(
    .fit16_raw(wrapped)[["target", exact = TRUE]][["K", exact = TRUE]],
    .fit16_raw(direct)[["target", exact = TRUE]][["K", exact = TRUE]]
  )
})


test_that("confidence and CV routes retain canonical request authority", {
  confidence <- DPprior_fit(
    20L, 4, confidence = "high", method = "A1",
    check_diagnostics = FALSE
  )
  confidence_raw <- .fit16_raw(confidence)
  confidence_target <- .fit16_raw(
    confidence_raw[["target", exact = TRUE]][["K", exact = TRUE]]
  )
  expect_identical(names(confidence_target[["request", exact = TRUE]]),
                   c("J", "mean", "confidence"))
  expect_identical(names(confidence_target[["used", exact = TRUE]]),
                   c("J", "mean", "variance", "interval", "pmf"))
  expect_false("var_K_used" %in% names(
    confidence_raw[["compatibility", exact = TRUE]][[
      "top_level_aliases", exact = TRUE
    ]]
  ))

  cv <- DPprior_fit(
    20L, 4, cv_K = sqrt(8) / 4, method = "A2-MN", M = 80L,
    check_diagnostics = FALSE
  )
  cv_target <- .fit16_raw(
    .fit16_raw(cv)[["target", exact = TRUE]][["K", exact = TRUE]]
  )
  expect_identical(cv_target[["kind", exact = TRUE]], "cv")
  expect_identical(names(cv_target[["request", exact = TRUE]]),
                   c("J", "mean", "cv"))
  expect_no_error(.dpprior_validate_result_v1(cv))
})


test_that("default confidence is explicit in the canonical target", {
  fit <- DPprior_fit(
    20L, mu_K = 4, method = "A2-MN", M = 80L,
    check_diagnostics = FALSE
  )
  target <- .fit16_raw(
    .fit16_raw(fit)[["target", exact = TRUE]][["K", exact = TRUE]]
  )
  expect_identical(target[["request", exact = TRUE]][[
    "confidence", exact = TRUE
  ]], "medium")
  expect_identical(target[["implied", exact = TRUE]][[
    "variance", exact = TRUE
  ]], 7.5)
})


test_that("strict PMF and interval routes dispatch to A2-KL", {
  pmf <- rep(1 / 12, 12)
  pmf_fit <- .fit16_expect_unusable_result(DPprior_fit(
    12L, target_pmf = pmf, M = 40L, check_diagnostics = FALSE
  ), "a2_kl")
  pmf_target <- .fit16_raw(
    .fit16_raw(pmf_fit)[["target", exact = TRUE]][["K", exact = TRUE]]
  )
  expect_identical(pmf_target[["kind", exact = TRUE]], "pmf")
  expect_identical(pmf_target[["pmf", exact = TRUE]], pmf)
  pmf_object <- DPprior_target_K(12L, target_pmf = pmf)
  pmf_object_fit <- .fit16_expect_unusable_result(DPprior_fit(
    12L, target_K = pmf_object, M = 40L,
    check_diagnostics = FALSE
  ), "a2_kl")
  expect_identical(
    .fit16_raw(pmf_object_fit)[["target", exact = TRUE]][["K", exact = TRUE]],
    pmf_object
  )

  interval <- list(
    lower = 2L, upper = 6L, type = "central_mass",
    coverage = 0.8, family = "maxent"
  )
  interval_fit <- .fit16_expect_unusable_result(DPprior_fit(
    20L, mu_K = 4, K_interval = interval, M = 40L,
    check_diagnostics = FALSE
  ), "a2_kl")
  interval_target <- .fit16_raw(
    .fit16_raw(interval_fit)[["target", exact = TRUE]][["K", exact = TRUE]]
  )
  expect_identical(interval_target[["kind", exact = TRUE]], "interval")
  expect_identical(interval_target[["request", exact = TRUE]][[
    "K_interval", exact = TRUE
  ]][["type", exact = TRUE]], "central_mass")
  expect_equal(sum(interval_target[["pmf", exact = TRUE]]), 1,
               tolerance = .TOL_PMF_SUM)
  interval_object <- DPprior_target_K(
    20L, mu_K = 4, K_interval = interval
  )
  interval_object_fit <- .fit16_expect_unusable_result(DPprior_fit(
    20L, target_K = interval_object, M = 40L,
    check_diagnostics = FALSE
  ), "a2_kl")
  expect_identical(
    .fit16_raw(interval_object_fit)[["target", exact = TRUE]][[
      "K", exact = TRUE
    ]],
    interval_object
  )
})


test_that("target_K input has exact scalar and J conflict rules", {
  target <- DPprior_target_K(20L, 4, confidence = "high")
  fit <- DPprior_fit(
    20L, target_K = target, method = "A2-MN", M = 80L,
    check_diagnostics = FALSE
  )
  .fit16_expect_canonical(fit, "a2_moment")
  expect_identical(
    .fit16_raw(fit)[["target", exact = TRUE]][["K", exact = TRUE]],
    target
  )
  resource <- .fit16_raw(fit)[["computation", exact = TRUE]][[
    "resources", exact = TRUE
  ]][["wrapper", exact = TRUE]]
  expect_identical(resource[["target_argument", exact = TRUE]], "target_K")

  conflict <- .fit16_condition(DPprior_fit(
    20L, mu_K = 4, target_K = target, method = "A2-MN",
    check_diagnostics = FALSE
  ))
  expect_s3_class(conflict, "dpprior_target_method_conflict")
  expect_identical(conflict[["code", exact = TRUE]],
                   "target_K_scalar_conflict")
  expect_s3_class(conflict[["result", exact = TRUE]], "dpprior_target")

  wrong_J <- .fit16_condition(DPprior_fit(
    21L, target_K = target, method = "A2-MN",
    check_diagnostics = FALSE
  ))
  expect_s3_class(wrong_J, "dpprior_target_integrity_error")
  expect_identical(wrong_J[["code", exact = TRUE]], "target_J_mismatch")
})


test_that("target_K routes match scalar normalized and used authority", {
  routes <- list(
    direct = list(
      target = DPprior_target_K(20L, 4, var_K = 8),
      scalar = function() DPprior_fit(
        20L, 4, 8, method = "A2-MN", M = 80L,
        check_diagnostics = FALSE
      )
    ),
    confidence = list(
      target = DPprior_target_K(20L, 4, confidence = "high"),
      scalar = function() DPprior_fit(
        20L, 4, confidence = "high", method = "A2-MN", M = 80L,
        check_diagnostics = FALSE
      )
    ),
    cv = list(
      target = DPprior_target_K(20L, 4, cv_K = sqrt(8) / 4),
      scalar = function() DPprior_fit(
        20L, 4, cv_K = sqrt(8) / 4, method = "A2-MN", M = 80L,
        check_diagnostics = FALSE
      )
    )
  )
  for (route in routes) {
    scalar <- route[["scalar"]]()
    object <- DPprior_fit(
      20L, target_K = route[["target"]], method = "A2-MN", M = 80L,
      check_diagnostics = FALSE
    )
    scalar_target <- .fit16_raw(
      .fit16_raw(scalar)[["target", exact = TRUE]][["K", exact = TRUE]]
    )
    object_target <- .fit16_raw(
      .fit16_raw(object)[["target", exact = TRUE]][["K", exact = TRUE]]
    )
    expect_identical(object_target[["normalized", exact = TRUE]],
                     scalar_target[["normalized", exact = TRUE]])
    expect_identical(object_target[["used", exact = TRUE]],
                     scalar_target[["used", exact = TRUE]])
  }
})


test_that("method and target conflicts fail before numerical dispatch", {
  pmf <- DPprior_target_K(10L, target_pmf = rep(0.1, 10))
  wrong <- .fit16_condition(DPprior_fit(
    10L, target_K = pmf, method = "A2-MN", check_diagnostics = FALSE
  ))
  expect_s3_class(wrong, "dpprior_target_method_conflict")
  expect_identical(wrong[["code", exact = TRUE]], "target_pmf_wrong_method")

  cv <- DPprior_target_K(20L, 4, cv_K = 0.5)
  wrong_cv <- .fit16_condition(DPprior_fit(
    20L, target_K = cv, method = "A2-KL", check_diagnostics = FALSE
  ))
  expect_s3_class(wrong_cv, "dpprior_target_method_conflict")
  expect_identical(wrong_cv[["code", exact = TRUE]],
                   "a2_kl_target_kind_unsupported")
})


test_that("A1 projection is retained only for a direct target", {
  fit <- NULL
  expect_warning(
    fit <- DPprior_fit(
      20L, 4, 2, method = "A1", a1_projection = "nearest",
      check_diagnostics = FALSE
    ),
    class = "dpprior_a1_projection_warning"
  )
  .fit16_expect_canonical(fit, "a1_proxy")
  target <- .fit16_raw(
    .fit16_raw(fit)[["target", exact = TRUE]][["K", exact = TRUE]]
  )
  expect_true(target[["provenance", exact = TRUE]][[
    "projection", exact = TRUE
  ]][["applied", exact = TRUE]])
  expect_identical(
    .fit16_raw(fit)[["provenance", exact = TRUE]][[
      "projection", exact = TRUE
    ]],
    target[["provenance", exact = TRUE]][["projection", exact = TRUE]]
  )

  derived <- .fit16_condition(DPprior_fit(
    20L, 4, cv_K = 0.2, method = "A1", a1_projection = "nearest",
    check_diagnostics = FALSE
  ))
  expect_s3_class(derived, "dpprior_a1_projection_policy_error")
  expect_identical(derived[["code", exact = TRUE]],
                   "a1_projection_derived_route_unsupported")
})


test_that("A1 projection-required condition retains a canonical target", {
  condition <- .fit16_condition(DPprior_fit(
    20L, 4, 2, method = "A1", check_diagnostics = FALSE
  ))
  expect_s3_class(condition, "dpprior_a1_projection_required")
  expect_identical(condition[["code", exact = TRUE]],
                   "a1_projection_required")
  expect_s3_class(condition[["result", exact = TRUE]], "dpprior_target")
  expect_no_error(.dpprior_validate_target_v1(
    condition[["result", exact = TRUE]]
  ))
})


test_that("native A2-MN fallback evidence is unchanged", {
  fallback_backend <- DPprior_a2_newton(
    50L, 3, 10, max_iter = 2L, M = 80L, verbose = FALSE
  )
  testthat::local_mocked_bindings(
    DPprior_a2_newton = function(...) fallback_backend,
    .package = "DPprior"
  )
  fit <- DPprior_fit(
    50L, 3, 10, method = "A2-MN", M = 80L,
    check_diagnostics = FALSE
  )
  .fit16_expect_canonical(fit, "a2_moment")
  raw <- .fit16_raw(fit)
  expect_identical(raw[["method", exact = TRUE]], "A2-MN+NM")
  expect_true(raw[["computation", exact = TRUE]][[
    "fallback", exact = TRUE
  ]][["used", exact = TRUE]])
  expect_identical(raw[["computation", exact = TRUE]][[
    "attempts", exact = TRUE
  ]], .fit16_raw(fallback_backend)[["computation", exact = TRUE]][[
    "attempts", exact = TRUE
  ]])
})


test_that("A2-MN nonconvergence is retained only on a typed condition", {
  backend <- DPprior_a2_newton(
    50L, 5, 8, M = 80L, max_iter = 1L, use_fallback = FALSE,
    tol_F = 1e-15, tol_rel = 0, verbose = FALSE
  )
  backend_raw <- .fit16_raw(backend)
  expect_identical(backend_raw[["status", exact = TRUE]], "approximate")
  expect_false(backend_raw[["usable", exact = TRUE]])
  expect_false(backend_raw[["verified", exact = TRUE]])
  testthat::local_mocked_bindings(
    DPprior_a2_newton = function(...) backend,
    .package = "DPprior"
  )

  wrapped <- .fit16_expect_unusable_result(DPprior_fit(
    50L, 5, 8, method = "A2-MN", M = 80L,
    check_diagnostics = FALSE
  ), "a2_moment")
  expect_identical(.fit16_scientific_core(wrapped),
                   .fit16_scientific_core(backend))
})


test_that("verified A2-MN boundary evidence remains publicly usable", {
  parameters <- c(a = exp(-15), b = 1)
  target <- exact_K_moments(
    50L, parameters[["a"]], parameters[["b"]], M = 80L
  )
  backend <- DPprior_a2_newton(
    50L, target[["mean"]], target[["var"]],
    a0 = parameters[["a"]], b0 = parameters[["b"]],
    max_iter = 1L, use_fallback = FALSE, M = 80L, verbose = FALSE
  )
  backend_raw <- .fit16_raw(backend)
  expect_identical(backend_raw[["status", exact = TRUE]], "boundary")
  expect_true(backend_raw[["usable", exact = TRUE]])
  expect_true(backend_raw[["verified", exact = TRUE]])
  testthat::local_mocked_bindings(
    DPprior_a2_newton = function(...) backend,
    .package = "DPprior"
  )

  wrapped <- DPprior_fit(
    50L, target[["mean"]], target[["var"]], method = "A2-MN",
    M = 80L, check_diagnostics = FALSE
  )
  .fit16_expect_canonical(wrapped, "a2_moment")
  expect_identical(.fit16_scientific_core(wrapped),
                   .fit16_scientific_core(backend))
})


test_that("M above the verifier ceiling remains a canonical no-candidate", {
  condition <- .fit16_condition(DPprior_fit(
    20L, 4, 8, method = "A2-MN", M = 257L,
    check_diagnostics = FALSE
  ))
  expect_s3_class(condition, "dpprior_a2_no_candidate")
  expect_identical(condition[["code", exact = TRUE]],
                   "a2_verification_order_exceeds_ceiling")
  result <- condition[["result", exact = TRUE]]
  .fit16_expect_canonical(result, "a2_moment")
  raw <- .fit16_raw(result)
  expect_identical(raw[["status", exact = TRUE]], "failed")
  expect_false(raw[["usable", exact = TRUE]])
  expect_false(raw[["verified", exact = TRUE]])
  expect_null(raw[["parameters", exact = TRUE]])
  expect_identical(raw[["computation", exact = TRUE]][[
    "resources", exact = TRUE
  ]][["wrapper", exact = TRUE]][["M_requested", exact = TRUE]], 257L)
})


test_that("fit-attached diagnostics use the exact authority contract", {
  fit <- DPprior_fit(20L, 4, 8, method = "A1", M = 80L)
  .fit16_expect_canonical(fit, "a1_proxy")
  raw <- .fit16_raw(fit)
  diagnostics <- raw[["diagnostics", exact = TRUE]]
  expect_identical(names(diagnostics), c(
    "authority", "policy_results", "warnings", "alpha", "K",
    "weights", "coclustering"
  ))
  authority <- diagnostics[["authority", exact = TRUE]]
  expect_identical(authority[["method", exact = TRUE]],
                   "fresh_component_specific_diagnostics")
  expect_identical(authority[["M_selected", exact = TRUE]], 80L)
  expect_identical(authority[["M_verification_required", exact = TRUE]],
                   160L)
  expect_identical(authority[["M_verification_used", exact = TRUE]], 160L)
  expect_identical(authority[["absolute_tolerance", exact = TRUE]], 1e-10)
  expect_identical(authority[["relative_tolerance", exact = TRUE]], 1e-8)
})


test_that("A2 diagnostics bind byte-identically to central orders", {
  fit <- DPprior_fit(
    50L, 5, 8, method = "A2-KL", M = 80L,
    check_diagnostics = TRUE
  )
  .fit16_expect_canonical(fit, "a2_kl")
  raw <- .fit16_raw(fit)
  expect_identical(raw[["status", exact = TRUE]], "converged")
  expect_true(raw[["usable", exact = TRUE]])
  expect_true(raw[["verified", exact = TRUE]])
  authority <- raw[["diagnostics", exact = TRUE]][[
    "authority", exact = TRUE
  ]]
  orders <- raw[["computation", exact = TRUE]][["orders", exact = TRUE]]
  expect_identical(authority[["M_selected", exact = TRUE]],
                   orders[["M_selected", exact = TRUE]])
  expect_identical(authority[["M_verification_required", exact = TRUE]],
                   orders[["M_verification_required", exact = TRUE]])
  expect_identical(authority[["M_verification_used", exact = TRUE]],
                   orders[["M_verification_used", exact = TRUE]])
  expect_false("diagnostics" %in% names(
    raw[["compatibility", exact = TRUE]][[
      "top_level_aliases", exact = TRUE
    ]]
  ))
  expect_true("diagnostics" %in% names(
    raw[["compatibility", exact = TRUE]][["views", exact = TRUE]]
  ))
})


test_that("diagnostic truth mutations are rejected by the shared validator", {
  fit <- DPprior_fit(20L, 4, 8, method = "A1", M = 80L)
  mutations <- list(
    alpha = function(x) {
      x[["diagnostics"]][["alpha"]][["mean"]] <- 9
      x
    },
    K = function(x) {
      x[["diagnostics"]][["K"]][["mean"]] <- 9
      x
    },
    weights = function(x) {
      x[["diagnostics"]][["weights"]][["mean"]] <- 0.9
      x
    },
    coclustering = function(x) {
      x[["diagnostics"]][["coclustering"]][["mean"]] <- 0.9
      x
    }
  )
  for (mutate in mutations) {
    bad <- mutate(fit)
    condition <- .fit16_condition(.dpprior_validate_result_v1(bad))
    expect_s3_class(condition, "dpprior_schema_error")
    expect_identical(condition[["code", exact = TRUE]],
                     "fit_diagnostic_component_truth")
  }
})


test_that("approximate diagnostics stop with the canonical fit attached", {
  condition <- .fit16_condition(DPprior_fit(
    20L, 4, 8, method = "A1", M = 20L,
    check_diagnostics = TRUE
  ))
  expect_s3_class(condition, "dpprior_diagnostics_approximation_error")
  expect_identical(condition[["code", exact = TRUE]],
                   "fit_diagnostics_approximate")
  result <- condition[["result", exact = TRUE]]
  .fit16_expect_canonical(result, "a1_proxy")
  expect_true("approximate" %in% vapply(
    c("alpha", "K", "weights", "coclustering"),
    function(name) .fit16_raw(result)[["diagnostics", exact = TRUE]][[
      name, exact = TRUE
    ]][["status", exact = TRUE]],
    character(1)
  ))
})


test_that("diagnostic policy warning is emitted only after successful return", {
  warnings <- list()
  condition <- tryCatch(
    withCallingHandlers(
      DPprior_fit(
        20L, 4, 8, method = "A1", M = 20L,
        warning_policy = list(
          estimand = "W_SB", direction = "above",
          weight_threshold = 0.5, action_threshold = 0
        )
      ),
      warning = function(warning) {
        warnings[[length(warnings) + 1L]] <<- warning
        invokeRestart("muffleWarning")
      }
    ),
    error = function(error) error
  )
  expect_s3_class(condition, "dpprior_diagnostics_approximation_error")
  expect_length(warnings, 0L)
  .fit16_expect_canonical(condition[["result", exact = TRUE]])
})


test_that("W_SB warns after validation and W_max is explicit unavailable", {
  WSB_policy <- list(
    estimand = "W_SB", direction = "above",
    weight_threshold = 0.5, action_threshold = 0
  )
  fit <- NULL
  expect_warning(
    fit <- DPprior_fit(
      20L, 4, 8, method = "A1", M = 80L,
      warning_policy = WSB_policy
    ),
    class = "dpprior_diagnostic_policy_warning"
  )
  .fit16_expect_canonical(fit)
  record <- .fit16_raw(fit)[["diagnostics", exact = TRUE]][[
    "policy_results", exact = TRUE
  ]][[1L]]
  expect_identical(record[["estimand", exact = TRUE]], "W_SB")
  expect_identical(record[["basis", exact = TRUE]],
                   "exact_tail_probability")
  expect_identical(record[["outcome", exact = TRUE]], "triggered")

  Wmax_policy <- list(
    estimand = "W_max", direction = "above",
    weight_threshold = 0.4, action_threshold = 0.2
  )
  Wmax <- expect_no_warning(DPprior_fit(
    20L, 4, 8, method = "A1", M = 80L,
    warning_policy = Wmax_policy
  ))
  Wmax_record <- .fit16_raw(Wmax)[["diagnostics", exact = TRUE]][[
    "policy_results", exact = TRUE
  ]][[1L]]
  expect_identical(Wmax_record[["basis", exact = TRUE]],
                   "backend_unavailable")
  expect_identical(Wmax_record[["outcome", exact = TRUE]], "indeterminate")
  expect_null(Wmax_record[["value", exact = TRUE]])
})


test_that("fit warning policies share the canonical diagnostics grammar", {
  canonical <- list(
    estimand = "W_SB", direction = "above",
    weight_threshold = 0.5, action_threshold = 1
  )
  legacy <- list(
    estimand = "W_SB", threshold = 0.5,
    direction = "above", action_threshold = 1
  )

  canonical_fit <- expect_no_warning(DPprior_fit(
    20L, 4, 8, method = "A1", M = 80L,
    warning_policy = canonical
  ))
  legacy_fit <- expect_no_warning(DPprior_fit(
    20L, 4, 8, method = "A1", M = 80L,
    warning_policy = legacy
  ))
  for (fit in list(canonical_fit, legacy_fit)) {
    raw <- .fit16_raw(fit)
    expect_identical(
      raw[["diagnostics", exact = TRUE]][["authority", exact = TRUE]][[
        "warning_policy", exact = TRUE
      ]],
      canonical
    )
    expect_identical(
      raw[["computation", exact = TRUE]][["resources", exact = TRUE]][[
        "wrapper", exact = TRUE
      ]][["warning_policy", exact = TRUE]],
      canonical
    )
  }
  expect_identical(canonical_fit, legacy_fit)

  posthoc <- expect_no_warning(DPprior_diagnostics(
    canonical_fit, warning_policy = canonical
  ))
  expect_s3_class(posthoc, "DPprior_diagnostics")
  expect_no_error(.dpprior_validate_result_v1(posthoc))

  malformed <- list(
    c(canonical, list(threshold = 0.5)),
    c(canonical, list(extra = TRUE)),
    structure(canonical, class = "policy"),
    structure(canonical, source = "forged"),
    canonical[c("direction", "estimand", "weight_threshold",
                "action_threshold")]
  )
  duplicate <- canonical
  names(duplicate)[[4L]] <- "weight_threshold"
  malformed[[length(malformed) + 1L]] <- duplicate
  for (policy in malformed) {
    condition <- .fit16_condition(DPprior_fit(
      20L, 4, 8, method = "A1", M = 80L,
      warning_policy = policy
    ))
    expect_s3_class(condition, "dpprior_warning_policy_error")
    expect_identical(condition[["code", exact = TRUE]],
                     "invalid_policy_record")
  }

  semantic_malformed <- list(
    within(canonical, estimand <- "W_S"),
    within(canonical, estimand <- c("W_SB", "W_max")),
    within(canonical, direction <- "a"),
    within(canonical, direction <- c("above", "below")),
    within(canonical, direction <- structure("above", class = "policy_value"))
  )
  for (policy in semantic_malformed) {
    condition <- .fit16_condition(DPprior_fit(
      20L, 4, 8, method = "A1", M = 80L,
      warning_policy = policy
    ))
    expect_s3_class(condition, "dpprior_warning_policy_error")
    expect_identical(condition[["code", exact = TRUE]], "choice")
  }
})


test_that("malformed backends become typed conditions without raw errors", {
  testthat::local_mocked_bindings(
    DPprior_a1 = function(...) stop("raw backend failure", call. = FALSE),
    .package = "DPprior"
  )
  condition <- .fit16_condition(DPprior_fit(
    20L, 4, 8, method = "A1", check_diagnostics = FALSE
  ))
  expect_s3_class(condition, "dpprior_backend_contract_error")
  expect_false(inherits(condition, "simpleError"))
  expect_identical(condition[["code", exact = TRUE]], "backend_raw_error")
  expect_s3_class(condition[["result", exact = TRUE]], "dpprior_target")
})


test_that("backend mode and J mismatches fail closed", {
  wrong <- DPprior_a1(21L, 4, 8)
  testthat::local_mocked_bindings(
    DPprior_a1 = function(...) wrong,
    .package = "DPprior"
  )
  condition <- .fit16_condition(DPprior_fit(
    20L, 4, 8, method = "A1", check_diagnostics = FALSE
  ))
  expect_s3_class(condition, "dpprior_backend_contract_error")
  expect_identical(condition[["code", exact = TRUE]],
                   "backend_dispatch_mismatch")
})


test_that("canonical and target serialization preserve exact validity", {
  target <- DPprior_target_K(20L, 4, confidence = "high")
  fit <- DPprior_fit(
    20L, target_K = target, method = "A1", M = 80L,
    check_diagnostics = TRUE
  )
  target_copy <- unserialize(serialize(target, NULL))
  fit_copy <- unserialize(serialize(fit, NULL))
  expect_identical(target_copy, target)
  expect_identical(fit_copy, fit)
  expect_no_error(.dpprior_validate_target_v1(target_copy))
  expect_no_error(.dpprior_validate_result_v1(fit_copy))
})


test_that("duplicate target fields and malformed enum attributes are typed", {
  target <- DPprior_target_K(20L, 4, var_K = 8)
  duplicate <- unclass(target)
  duplicate <- c(duplicate, list(J = 20L))
  class(duplicate) <- class(target)
  duplicate_condition <- .fit16_condition(DPprior_fit(
    20L, target_K = duplicate, method = "A2-MN",
    check_diagnostics = FALSE
  ))
  expect_s3_class(duplicate_condition, "dpprior_serialization_error")

  for (bad in list(matrix("A1"), structure("A1", class = "evil"))) {
    condition <- .fit16_condition(DPprior_fit(
      20L, 4, 8, method = bad, check_diagnostics = FALSE
    ))
    expect_s3_class(condition, "dpprior_method_error")
  }
})


test_that("forged target accessors cannot redirect wrapper target reads", {
  target <- DPprior_target_K(20L, 4, confidence = "high")
  backend <- DPprior_a1(20L, 4, 4.5)
  testthat::local_mocked_bindings(
    DPprior_a1 = function(...) backend,
    .package = "DPprior"
  )
  assign("$.dpprior_K_target", function(...) stop("poison target dollar"),
         envir = .GlobalEnv)
  assign("[[.dpprior_K_target", function(...) stop("poison target bracket"),
         envir = .GlobalEnv)
  on.exit({
    rm("$.dpprior_K_target", envir = .GlobalEnv)
    rm("[[.dpprior_K_target", envir = .GlobalEnv)
  }, add = TRUE)

  fit <- DPprior_fit(
    20L, target_K = target, method = "A1",
    check_diagnostics = FALSE
  )
  .fit16_expect_canonical(fit, "a1_proxy")
  retained <- .fit16_raw(
    .fit16_raw(fit)[["target", exact = TRUE]][["K", exact = TRUE]]
  )
  expect_identical(retained[["request", exact = TRUE]][[
    "confidence", exact = TRUE
  ]], "high")
})


test_that("non-authoritative legacy-view poison cannot change science", {
  backend <- DPprior_a2_newton(20L, 4, 8, M = 80L, verbose = FALSE)
  poisoned <- .fit16_raw(backend)
  compatibility <- poisoned[["compatibility", exact = TRUE]]
  views <- compatibility[["views", exact = TRUE]]
  legacy <- views[["legacy_v2", exact = TRUE]]
  selected_fit <- legacy[["selected_fit", exact = TRUE]]
  selected_fit[["mu_K"]] <- 999
  legacy[["selected_fit"]] <- selected_fit
  views[["legacy_v2"]] <- legacy
  compatibility[["views"]] <- views
  poisoned[["compatibility"]] <- compatibility
  poisoned[["fit"]] <- selected_fit
  class(poisoned) <- class(backend)
  expect_no_error(.dpprior_validate_result_v1(poisoned))

  testthat::local_mocked_bindings(
    DPprior_a2_newton = function(...) poisoned,
    .package = "DPprior"
  )
  fit <- DPprior_fit(
    20L, 4, 8, method = "A2-MN", M = 80L,
    check_diagnostics = FALSE
  )
  .fit16_expect_canonical(fit, "a2_moment")
  expect_identical(.fit16_raw(fit)[["achieved", exact = TRUE]],
                   .fit16_raw(backend)[["achieved", exact = TRUE]])
  expect_false(identical(
    .fit16_raw(fit)[["achieved", exact = TRUE]][["K", exact = TRUE]][[
      "mean", exact = TRUE
    ]],
    .fit16_raw(fit)[["compatibility", exact = TRUE]][[
      "views", exact = TRUE
    ]][["legacy_v2", exact = TRUE]][["selected_fit", exact = TRUE]][[
      "mu_K", exact = TRUE
    ]]
  ))
})


test_that("forged dollar and bracket methods cannot redirect canonical reads", {
  assign("$.DPprior_fit", function(...) stop("poison dollar"),
         envir = .GlobalEnv)
  assign("[[.DPprior_fit", function(...) stop("poison bracket"),
         envir = .GlobalEnv)
  on.exit({
    rm("$.DPprior_fit", envir = .GlobalEnv)
    rm("[[.DPprior_fit", envir = .GlobalEnv)
  }, add = TRUE)

  fit <- DPprior_fit(
    20L, 4, 8, method = "A1", M = 80L,
    check_diagnostics = TRUE
  )
  .fit16_expect_canonical(fit, "a1_proxy")
  expect_no_error(capture.output(print(fit)))
  expect_no_error(capture.output(summary(fit)))
  expect_s3_class(as.data.frame(fit), "data.frame")
  expect_no_error(plot(
    fit, type = "alpha", engine = "base", show = FALSE
  ))
})


test_that("S3 consumers read canonical wrapper output", {
  fit <- DPprior_fit(
    20L, 4, 8, method = "A1", M = 80L,
    check_diagnostics = TRUE
  )
  printed <- capture.output(print(fit))
  expect_true(any(grepl("dpprior.result/1", printed, fixed = TRUE)))
  summary_fit <- NULL
  capture.output(summary_fit <- summary(fit))
  expect_true(is.list(summary_fit) || is.data.frame(summary_fit))
  expect_s3_class(as.data.frame(fit), "data.frame")
  expect_no_error(plot(
    fit, type = "K", engine = "base", show = FALSE
  ))
})


test_that("verify_DPprior_fit validates the canonical smoke route", {
  expect_true(verify_DPprior_fit(verbose = FALSE))
})
