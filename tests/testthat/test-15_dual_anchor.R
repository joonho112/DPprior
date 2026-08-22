.legacy_dual_fit <- function(J = 20, mu_K = 5, var_K = 8, M = 80) {
  DPprior_fit(
    J = J, mu_K = mu_K, var_K = var_K,
    method = "A2-MN", M = M, check_diagnostics = FALSE
  )
}

.legacy_dual_target <- function(value = 0.3, threshold = 0.5) {
  list(prob = list(threshold = threshold, value = value))
}

.capture_legacy_call <- function(code) {
  warnings <- list()
  value <- tryCatch(
    withCallingHandlers(
      code(),
      warning = function(condition) {
        warnings[[length(warnings) + 1L]] <<- condition
        invokeRestart("muffleWarning")
      }
    ),
    error = identity
  )
  list(value = value, warnings = warnings)
}

.legacy_weight_fresh <- function(raw) {
  target <- unclass(raw[["target", exact = TRUE]][["weight", exact = TRUE]])
  parameters <- raw[["parameters", exact = TRUE]]
  switch(
    target[["metric", exact = TRUE]],
    wsb_tail = prob_wsb_exceeds(
      target[["threshold", exact = TRUE]],
      parameters[["a", exact = TRUE]],
      parameters[["b", exact = TRUE]]
    ),
    wsb_mean = mean_w1(
      parameters[["a", exact = TRUE]],
      parameters[["b", exact = TRUE]],
      raw[["computation", exact = TRUE]][["orders", exact = TRUE]][[
        "M_selected", exact = TRUE
      ]]
    ),
    wsb_quantile = quantile_w1(
      target[["probability", exact = TRUE]],
      parameters[["a", exact = TRUE]],
      parameters[["b", exact = TRUE]]
    )
  )
}

.legacy_recursive_names <- function(value) {
  if (!is.list(value)) return(character())
  raw <- if (is.object(value)) unclass(value) else value
  c(
    names(raw),
    unlist(lapply(raw, .legacy_recursive_names), use.names = FALSE)
  )
}

.legacy_plain_tree <- function(value) {
  if (is.null(value)) return(TRUE)
  if (is.object(value) || !is.null(dim(value))) return(FALSE)
  if (!is.list(value)) return(TRUE)
  all(vapply(value, .legacy_plain_tree, logical(1)))
}

test_that("legacy warning is emitted exactly once and only after validation", {
  fit <- .legacy_dual_fit()
  success <- .capture_legacy_call(function() {
    DPprior_dual(fit, .legacy_dual_target(), lambda = 1, M = 80)
  })
  expect_length(success[["warnings"]], 1L)
  expect_s3_class(
    success[["warnings"]][[1L]], "dpprior_dual_legacy_warning"
  )
  expect_s3_class(
    success[["warnings"]][[1L]], "dpprior_deprecated_warning"
  )
  expect_identical(
    success[["warnings"]][[1L]][["code", exact = TRUE]],
    "deprecated_legacy_dual_anchor"
  )
  expect_s3_class(success[["value"]], "DPprior_fit")

  invalid_calls <- list(
    function() {
      DPprior_dual(fit, .legacy_dual_target(), lambda = matrix(0.5))
    },
    function() DPprior_dual(fit, list(), lambda = 0.5),
    function() {
      DPprior_dual(
        fit,
        list(prob = list(threshold = 0.5, value = 0.3), mean = 0.3),
        lambda = 0.5
      )
    },
    function() {
      DPprior_dual(
        list(a = 1, b = 1, J = 20), .legacy_dual_target(), lambda = 0.5
      )
    }
  )
  for (call in invalid_calls) {
    captured <- .capture_legacy_call(call)
    expect_s3_class(captured[["value"]], "error")
    expect_length(captured[["warnings"]], 0L)
  }
})

test_that("max_iter is bounded before warning and fallback arithmetic", {
  fit <- .legacy_dual_fit()
  safe_max_iter <- as.integer(floor(.Machine$integer.max / 2))
  success <- .capture_legacy_call(function() {
    DPprior_dual(
      fit, .legacy_dual_target(), lambda = 1,
      max_iter = safe_max_iter, M = 80
    )
  })
  expect_length(success[["warnings"]], 1L)
  expect_s3_class(
    success[["warnings"]][[1L]], "dpprior_dual_legacy_warning"
  )
  expect_s3_class(success[["value"]], "DPprior_fit")
  success_raw <- unclass(success[["value"]])
  expect_identical(
    success_raw[["computation", exact = TRUE]][["request", exact = TRUE]][[
      "controls", exact = TRUE
    ]][["max_iter", exact = TRUE]],
    safe_max_iter
  )
  expect_identical(
    success_raw[["computation", exact = TRUE]][["resources", exact = TRUE]][[
      "optimizer_controls", exact = TRUE
    ]][["fallback", exact = TRUE]][["maxit", exact = TRUE]],
    2L * safe_max_iter
  )
  .dpprior_validate_result_v1(success[["value"]])

  expected_class <- c(
    "dpprior_dual_legacy_control_error",
    "dpprior_bounds_error",
    "dpprior_invalid_input",
    "dpprior_error",
    "error",
    "dpprior_condition",
    "condition"
  )
  invalid_values <- list(
    as.numeric(safe_max_iter) + 1,
    .Machine$integer.max
  )
  for (bad in invalid_values) {
    captured <- .capture_legacy_call(function() {
      DPprior_dual(
        fit, .legacy_dual_target(), lambda = 1,
        max_iter = bad, M = 80
      )
    })
    condition <- captured[["value"]]
    expect_identical(class(condition), expected_class)
    expect_identical(condition[["code", exact = TRUE]], "bounds")
    expect_identical(condition[["argument", exact = TRUE]], "max_iter")
    expect_length(captured[["warnings"]], 0L)
  }
})

test_that("public loss_type failures are exact typed pre-warning errors", {
  fit <- .legacy_dual_fit()
  bad_values <- list(
    1,
    "bogus",
    c("relative", "adaptive"),
    structure("relative", class = "evil_loss_type")
  )
  expected_class <- c(
    "dpprior_dual_legacy_loss_type_error",
    "dpprior_dual_legacy_invalid_input",
    "dpprior_dual_anchor_error",
    "dpprior_invalid_input",
    "dpprior_error",
    "error",
    "dpprior_condition",
    "condition"
  )
  for (bad in bad_values) {
    captured <- .capture_legacy_call(function() {
      DPprior_dual(
        fit, .legacy_dual_target(), lambda = 0.5,
        loss_type = bad, M = 80
      )
    })
    condition <- captured[["value"]]
    expect_identical(class(condition), expected_class)
    expect_identical(
      condition[["code", exact = TRUE]], "legacy_dual_loss_type"
    )
    expect_identical(condition[["argument", exact = TRUE]], "loss_type")
    expect_identical(
      condition[["expected", exact = TRUE]],
      "one of relative, adaptive, or absolute"
    )
    expect_length(captured[["warnings"]], 0L)
  }
})

test_that("mutated canonical hard fields or mode are rejected, not stripped", {
  fit <- .legacy_dual_fit()
  raw <- unclass(fit)
  raw[["constraint_satisfied"]] <- TRUE
  raw[["constraint_residual"]] <- -1
  raw[["constraint_tolerance"]] <- 1
  raw[["mode"]] <- "dual_hard"
  raw[["provenance"]][["constraint_satisfied"]] <- TRUE
  class(raw) <- class(fit)

  captured <- .capture_legacy_call(function() {
    DPprior_dual(raw, .legacy_dual_target(), lambda = 1, M = 80)
  })
  expect_s3_class(captured[["value"]], "dpprior_schema_error")
  expect_true(nzchar(captured[["value"]][["code", exact = TRUE]]))
  expect_length(captured[["warnings"]], 0L)
})

test_that("legacy producer publishes the exact canonical spine and quartet", {
  fit <- .legacy_dual_fit()
  result <- suppressWarnings(DPprior_dual(
    fit, .legacy_dual_target(), lambda = 0.5, M = 80
  ))
  .dpprior_validate_result_v1(result)
  raw <- unclass(result)
  expected_common <- c(
    "schema", "object_type", "mode", "method", "J",
    "status", "usable", "verified", "message",
    "parameters", "target", "achieved", "residuals", "tolerances",
    "computation", "verification", "provenance", "compatibility"
  )
  expected_aliases <- c(
    "a", "b", "converged", "iterations", "fit", "attempts", "dual_anchor"
  )
  expect_identical(
    names(raw), c(expected_common, "legacy", expected_aliases)
  )
  expect_identical(class(result), c("DPprior_fit", "dpprior_result", "list"))
  expect_identical(raw[["object_type", exact = TRUE]], "fit")
  expect_identical(raw[["mode", exact = TRUE]], "dual_legacy")
  expect_identical(raw[["method", exact = TRUE]], "dual-anchor")
  expect_identical(raw[["status", exact = TRUE]], "approximate")
  expect_identical(raw[["usable", exact = TRUE]], TRUE)
  expect_identical(raw[["verified", exact = TRUE]], FALSE)
  expect_identical(raw[["converged", exact = TRUE]], FALSE)
  expect_true(raw[["provenance", exact = TRUE]][[
    "legacy", exact = TRUE
  ]][["active", exact = TRUE]])
  expect_true(raw[["provenance", exact = TRUE]][[
    "approximation", exact = TRUE
  ]][["opt_in", exact = TRUE]])
  expect_null(raw[["provenance", exact = TRUE]][[
    "input_fit", exact = TRUE
  ]])
  input_raw <- unclass(fit)
  expect_identical(
    raw[["target", exact = TRUE]][["K", exact = TRUE]],
    input_raw[["target", exact = TRUE]][["K", exact = TRUE]]
  )
  expect_identical(
    serialize(raw[["target", exact = TRUE]][["K", exact = TRUE]], NULL),
    serialize(input_raw[["target", exact = TRUE]][["K", exact = TRUE]], NULL)
  )
})

test_that("lambda one is deterministic with an empty canonical ledger", {
  fit <- .legacy_dual_fit()
  result <- suppressWarnings(DPprior_dual(
    fit, .legacy_dual_target(), lambda = 1L, M = 80
  ))
  raw <- unclass(result)
  input_raw <- unclass(fit)
  expect_identical(
    raw[["parameters", exact = TRUE]],
    input_raw[["parameters", exact = TRUE]]
  )
  expect_identical(raw[["computation", exact = TRUE]][[
    "attempts", exact = TRUE
  ]], list())
  expect_null(raw[["computation", exact = TRUE]][[
    "selected_attempt_id", exact = TRUE
  ]])
  expect_identical(
    raw[["computation", exact = TRUE]][["termination", exact = TRUE]][[
      "code", exact = TRUE
    ]],
    "deterministic"
  )
  expect_identical(
    raw[["computation", exact = TRUE]][["termination", exact = TRUE]][[
      "source", exact = TRUE
    ]],
    "legacy_adapter"
  )
  expect_identical(raw[["legacy", exact = TRUE]][[
    "lambda", exact = TRUE
  ]], 1)
  baseline <- raw[["computation", exact = TRUE]][[
    "resources", exact = TRUE
  ]][["K_only_baseline", exact = TRUE]]
  expect_true(.legacy_plain_tree(baseline))
  expect_identical(
    baseline[["parameters", exact = TRUE]],
    input_raw[["parameters", exact = TRUE]]
  )
  .dpprior_validate_result_v1(result)
})

test_that("decision-ready A2-KL input retains exact K authority", {
  fit <- DPprior_fit(
    50L, 5, 8, method = "A2-KL", M = 80L,
    check_diagnostics = TRUE
  )
  result <- suppressWarnings(DPprior_dual(
    fit, list(mean = 0.3), lambda = 1, M = 80
  ))
  raw <- unclass(result)
  fit_raw <- unclass(fit)
  expect_identical(
    raw[["target", exact = TRUE]][["K", exact = TRUE]],
    fit_raw[["target", exact = TRUE]][["K", exact = TRUE]]
  )
  expect_identical(
    raw[["computation", exact = TRUE]][["resources", exact = TRUE]][[
      "K_only_baseline", exact = TRUE
    ]][["mode", exact = TRUE]],
    "a2_kl"
  )
  expect_s3_class(dual_anchor_diagnostics(result, M = 80), "data.frame")
  .dpprior_validate_result_v1(result)
})

test_that("lambda zero remains a public legacy value with typed selection", {
  fit <- .legacy_dual_fit()
  result <- suppressWarnings(DPprior_dual(
    fit, .legacy_dual_target(0.3, 0.41),
    lambda = 0L, loss_type = "relative", M = 80
  ))
  raw <- unclass(result)
  attempts <- raw[["computation", exact = TRUE]][["attempts", exact = TRUE]]
  selected <- which(vapply(
    attempts,
    function(attempt) isTRUE(attempt[["selected", exact = TRUE]]),
    logical(1)
  ))
  expect_identical(raw[["legacy", exact = TRUE]][[
    "lambda", exact = TRUE
  ]], 0)
  expect_identical(
    raw[["computation", exact = TRUE]][["request", exact = TRUE]][[
      "controls", exact = TRUE
    ]][["lambda", exact = TRUE]],
    0
  )
  expect_length(selected, 1L)
  expect_identical(
    raw[["computation", exact = TRUE]][[
      "selected_attempt_id", exact = TRUE
    ]],
    attempts[[selected]][["id", exact = TRUE]]
  )
  expect_identical(
    raw[["computation", exact = TRUE]][["termination", exact = TRUE]][[
      "iterations", exact = TRUE
    ]],
    attempts[[selected]][["iterations", exact = TRUE]]
  )
  expect_identical(
    raw[["achieved", exact = TRUE]][["weight", exact = TRUE]][[
      "value", exact = TRUE
    ]],
    .legacy_weight_fresh(raw)
  )
  .dpprior_validate_result_v1(result)
})

test_that("all loss types retain historical scaling and numerical behavior", {
  fit <- .legacy_dual_fit()
  results <- list()
  for (loss_type in c("relative", "adaptive", "absolute")) {
    result <- suppressWarnings(DPprior_dual(
      fit, .legacy_dual_target(),
      lambda = 0.7, loss_type = loss_type, M = 80
    ))
    raw <- unclass(result)
    expect_identical(raw[["legacy", exact = TRUE]][[
      "losses", exact = TRUE
    ]][["loss_type", exact = TRUE]], loss_type)
    expect_true(is.finite(raw[["legacy", exact = TRUE]][[
      "losses", exact = TRUE
    ]][["total_loss", exact = TRUE]]))
    .dpprior_validate_result_v1(result)
    results[[loss_type]] <- result
  }
  adaptive_raw <- unclass(results[["adaptive"]])
  adaptive_attempts <- adaptive_raw[["computation", exact = TRUE]][[
    "attempts", exact = TRUE
  ]]
  expect_identical(
    adaptive_attempts[[1L]][["id", exact = TRUE]],
    "attempt-scaling-001"
  )
  expect_identical(
    adaptive_attempts[[1L]][["stage", exact = TRUE]], "scaling"
  )
  expect_false(adaptive_attempts[[1L]][["selected", exact = TRUE]])
  expect_identical(
    names(adaptive_raw[["computation", exact = TRUE]][[
      "scaling", exact = TRUE
    ]][["values", exact = TRUE]]),
    c("L_K_scale", "L_w_scale")
  )

  frozen_fit <- .legacy_dual_fit(J = 50, mu_K = 3, var_K = 10)
  frozen <- suppressWarnings(DPprior_dual(
    frozen_fit, .legacy_dual_target(0.25, 0.5),
    lambda = 0.7, loss_type = "adaptive", M = 80
  ))
  frozen_raw <- unclass(frozen)
  achieved <- frozen_raw[["compatibility", exact = TRUE]][[
    "views", exact = TRUE
  ]][["dual_anchor", exact = TRUE]][["w1_achieved", exact = TRUE]][[
    "prob_gt_50", exact = TRUE
  ]]
  expect_equal(achieved, 0.2650345, tolerance = 2e-5)
  expect_equal(achieved - 0.25, 0.0150345, tolerance = 2e-5)
})

test_that("standard relative-loss science values remain stable", {
  fit <- .legacy_dual_fit(J = 50, mu_K = 5, var_K = 8)
  result <- suppressWarnings(DPprior_dual(
    fit, .legacy_dual_target(0.3, 0.5),
    lambda = 0.5, loss_type = "relative", M = 80
  ))
  raw <- unclass(result)
  expect_equal(
    raw[["parameters", exact = TRUE]][["a", exact = TRUE]],
    2.575222603531012, tolerance = 1e-7
  )
  expect_equal(
    raw[["parameters", exact = TRUE]][["b", exact = TRUE]],
    1.833605239665517, tolerance = 1e-7
  )
  expect_equal(
    raw[["achieved", exact = TRUE]][["K", exact = TRUE]][[
      "mean", exact = TRUE
    ]],
    5.380060410495426, tolerance = 1e-7
  )
  expect_equal(
    raw[["achieved", exact = TRUE]][["K", exact = TRUE]][[
      "variance", exact = TRUE
    ]],
    7.930246930841893, tolerance = 1e-7
  )
  expect_equal(
    raw[["achieved", exact = TRUE]][["weight", exact = TRUE]][[
      "value", exact = TRUE
    ]],
    0.437907733205280, tolerance = 1e-7
  )
  expect_equal(
    raw[["legacy", exact = TRUE]][["losses", exact = TRUE]][[
      "total_loss", exact = TRUE
    ]],
    0.01243620139718597, tolerance = 1e-7
  )
  .dpprior_validate_result_v1(result)
})

test_that("all legacy target forms publish and freshly evaluate their metric", {
  fit <- .legacy_dual_fit()
  cases <- list(
    list(
      input = list(prob = list(threshold = 0.37, value = 0.3)),
      metric = "wsb_tail", threshold = 0.37, probability = NULL
    ),
    list(
      input = list(mean = 0.3),
      metric = "wsb_mean", threshold = NULL, probability = NULL
    ),
    list(
      input = list(quantile = list(prob = 0.73, value = 0.4)),
      metric = "wsb_quantile", threshold = NULL, probability = 0.73
    )
  )
  for (case in cases) {
    result <- suppressWarnings(DPprior_dual(
      fit, case[["input"]], lambda = 1, M = 80
    ))
    raw <- unclass(result)
    weight <- unclass(raw[["target", exact = TRUE]][[
      "weight", exact = TRUE
    ]])
    authority <- list(
      metric = weight[["metric", exact = TRUE]],
      relation = weight[["relation", exact = TRUE]],
      value = weight[["value", exact = TRUE]],
      threshold = weight[["threshold", exact = TRUE]],
      probability = weight[["probability", exact = TRUE]]
    )
    expect_identical(weight[["metric", exact = TRUE]], case[["metric"]])
    expect_identical(weight[["threshold", exact = TRUE]], case[["threshold"]])
    expect_identical(
      weight[["probability", exact = TRUE]], case[["probability"]]
    )
    expect_identical(weight[["request", exact = TRUE]], authority)
    expect_identical(weight[["normalized", exact = TRUE]], authority)
    expect_identical(weight[["used", exact = TRUE]], authority)
    expect_identical(
      raw[["achieved", exact = TRUE]][["weight", exact = TRUE]][[
        "value", exact = TRUE
      ]],
      .legacy_weight_fresh(raw)
    )
    .dpprior_validate_result_v1(result)
  }
})

test_that("natural nonzero primary exit selects and binds the typed fallback", {
  fit <- .legacy_dual_fit(J = 50, mu_K = 3, var_K = 10)
  result <- suppressWarnings(DPprior_dual(
    fit, .legacy_dual_target(0.25, 0.5), lambda = 0.5,
    loss_type = "adaptive", M = 80
  ))
  raw <- unclass(result)
  attempts <- raw[["computation", exact = TRUE]][["attempts", exact = TRUE]]
  fallback <- raw[["computation", exact = TRUE]][[
    "fallback", exact = TRUE
  ]]
  termination <- raw[["computation", exact = TRUE]][[
    "termination", exact = TRUE
  ]]
  expect_identical(
    vapply(attempts, function(x) x[["id", exact = TRUE]], character(1)),
    c(
      "attempt-scaling-001", "attempt-primary-001",
      "attempt-fallback-001"
    )
  )
  expect_identical(attempts[[1L]][["exit_code", exact = TRUE]], 0L)
  expect_identical(attempts[[2L]][["exit_code", exact = TRUE]], 1L)
  expect_identical(
    attempts[[2L]][["reason_code", exact = TRUE]],
    "optimizer_exit_nonzero"
  )
  expect_identical(attempts[[3L]][["exit_code", exact = TRUE]], 0L)
  expect_true(attempts[[3L]][["selected", exact = TRUE]])
  expect_identical(
    raw[["computation", exact = TRUE]][[
      "selected_attempt_id", exact = TRUE
    ]],
    "attempt-fallback-001"
  )
  expect_true(fallback[["attempted", exact = TRUE]])
  expect_true(fallback[["used", exact = TRUE]])
  expect_identical(
    fallback[["trigger_attempt_id", exact = TRUE]],
    "attempt-primary-001"
  )
  expect_identical(
    fallback[["selected_attempt_id", exact = TRUE]],
    "attempt-fallback-001"
  )
  expect_identical(termination[["source", exact = TRUE]], "fallback_optimizer")
  expect_identical(
    termination[["iterations", exact = TRUE]],
    attempts[[3L]][["iterations", exact = TRUE]]
  )
  expect_true(raw[["provenance", exact = TRUE]][[
    "is_fallback", exact = TRUE
  ]])
  .dpprior_validate_result_v1(result)
})

test_that("diagnostics consume canonical evidence and ignore compatibility", {
  fit <- .legacy_dual_fit()
  legacy <- suppressWarnings(DPprior_dual(
    fit, .legacy_dual_target(), lambda = 0.5, M = 80
  ))
  diagnostic <- dual_anchor_diagnostics(legacy, M = 80)
  expect_s3_class(diagnostic, "data.frame")
  expect_identical(nrow(diagnostic), 7L)
  expect_identical(
    names(diagnostic), c("Metric", "K_only", "Dual_anchor")
  )
  expect_true(all(c(
    "P(W_SB > 0.5)", "P(W_SB > 0.9)", "E[W_SB]"
  ) %in% diagnostic[["Metric"]]))
  expect_identical(
    dual_anchor_diagnostics(legacy, fit, M = 80), diagnostic
  )

  poisoned <- unclass(legacy)
  poisoned[["compatibility"]][["views"]][["dual_anchor"]][[
    "init"
  ]] <- list(a = 999, b = 999)
  poisoned[["dual_anchor"]] <- poisoned[["compatibility"]][[
    "views"
  ]][["dual_anchor"]]
  class(poisoned) <- class(legacy)
  .dpprior_validate_result_v1(poisoned)
  expect_identical(
    dual_anchor_diagnostics(poisoned, fit, M = 80), diagnostic
  )

  baseline_poison <- unclass(legacy)
  baseline_poison[["computation"]][["resources"]][[
    "K_only_baseline"
  ]][["parameters"]][["a"]] <- 999
  class(baseline_poison) <- class(legacy)
  expect_error(
    dual_anchor_diagnostics(baseline_poison, M = 80),
    class = "dpprior_dual_diagnostics_invalid_input"
  )

  unrelated <- .legacy_dual_fit(mu_K = 7, var_K = 10)
  expect_error(
    dual_anchor_diagnostics(legacy, unrelated, M = 80),
    class = "dpprior_dual_diagnostics_invalid_input"
  )
})

test_that("diagnostics accept the frozen canonical soft producer", {
  fit <- .legacy_dual_fit()
  soft <- DPprior_dual_soft(
    fit,
    list(
      metric = "wsb_tail", relation = "target",
      threshold = 0.37, value = 0.3
    ),
    lambda = 0.7, M_fit = 80
  )
  diagnostic <- dual_anchor_diagnostics(soft, M = 80)
  expect_s3_class(diagnostic, "data.frame")
  expect_identical(nrow(diagnostic), 7L)
  expect_identical(
    dual_anchor_diagnostics(soft, fit, M = 80), diagnostic
  )

  incomplete <- structure(
    list(schema = list(name = "dpprior.result", version = 1L)),
    class = c("DPprior_dual_soft", "DPprior_fit")
  )
  expect_error(
    dual_anchor_diagnostics(incomplete, M = 80),
    class = "dpprior_dual_diagnostics_invalid_input"
  )
})

test_that("hostile accessors, S3 grafts, and serialization fail closed", {
  fit <- .legacy_dual_fit()
  `[[.DPprior_fit` <- function(...) {
    stop("hostile bracket accessor ran")
  }
  `$.DPprior_fit` <- function(...) {
    stop("hostile dollar accessor ran")
  }
  result <- suppressWarnings(DPprior_dual(
    fit, .legacy_dual_target(), lambda = 1, M = 80
  ))
  result_raw <- unclass(result)
  expect_identical(result_raw[["mode", exact = TRUE]], "dual_legacy")
  expect_s3_class(
    dual_anchor_diagnostics(result, fit, M = 80), "data.frame"
  )

  wire <- serialize(result, NULL, xdr = TRUE)
  restored <- unserialize(wire)
  expect_identical(restored, result)
  expect_identical(class(restored), c(
    "DPprior_fit", "dpprior_result", "list"
  ))
  .dpprior_validate_result_v1(restored)

  grafted <- unclass(fit)
  class(grafted) <- c("evil_fit", class(fit))
  captured <- .capture_legacy_call(function() {
    DPprior_dual(grafted, .legacy_dual_target(), lambda = 1, M = 80)
  })
  expect_s3_class(captured[["value"]], "dpprior_schema_error")
  expect_length(captured[["warnings"]], 0L)

  classed_target <- list(
    prob = structure(
      list(threshold = 0.5, value = 0.3), class = "evil_target"
    )
  )
  captured_target <- .capture_legacy_call(function() {
    DPprior_dual(fit, classed_target, lambda = 1, M = 80)
  })
  expect_s3_class(captured_target[["value"]], "dpprior_invalid_input")
  expect_length(captured_target[["warnings"]], 0L)
})

test_that("compatibility is quarantined and canonical evidence has no hard keys", {
  fit <- .legacy_dual_fit()
  result <- suppressWarnings(DPprior_dual(
    fit, .legacy_dual_target(), lambda = 0.5, M = 80
  ))
  raw <- unclass(result)
  boundary <- .dpprior_compatibility_quarantine_boundary()
  views <- raw[["compatibility", exact = TRUE]][["views", exact = TRUE]]
  deprecation <- raw[["compatibility", exact = TRUE]][[
    "deprecations", exact = TRUE
  ]][["legacy_dual_v2", exact = TRUE]]
  expect_identical(
    views[["legacy_dual_v2", exact = TRUE]][names(boundary)], boundary
  )
  expect_identical(
    views[["dual_anchor", exact = TRUE]][names(boundary)], boundary
  )
  expect_identical(deprecation[names(boundary)], boundary)
  aliases <- names(raw[["compatibility", exact = TRUE]][[
    "top_level_aliases", exact = TRUE
  ]])
  canonical <- raw
  canonical[["compatibility"]] <- NULL
  canonical[aliases] <- NULL
  forbidden <- c(
    "constraint", "constraint_satisfied", "constraint_residual",
    "constraint_slack", "constraint_tolerance", "feasibility"
  )
  expect_length(
    intersect(.legacy_recursive_names(canonical), forbidden), 0L
  )
  .dpprior_validate_result_v1(result)
})

test_that("legacy verification helper validates the canonical migration", {
  expect_true(verify_dual_anchor(verbose = FALSE))
})
