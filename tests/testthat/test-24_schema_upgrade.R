make_v1_a2_mn_fit <- function() {
  structure(list(
    a = 1.27223464009856,
    b = 0.58075377814199,
    J = 20,
    target = list(mu_K = 5, var_K = 8, type = "moments"),
    method = "A2-MN",
    status = "success",
    converged = TRUE,
    iterations = 7L,
    termination = "residual",
    fit = list(
      mu_K = 5.00000000001361,
      var_K = 7.99999999987637,
      residual = 1.24381729905704e-10
    ),
    diagnostics = list(
      a0 = 4, b0 = 2.99573227355399, tol_F = 1e-8,
      tol_step = 1e-10, M = 80L, fallback_used = FALSE
    ),
    trace = data.frame(iteration = 1L, residual = 1.24381729905704e-10)
  ), class = "DPprior_fit")
}


make_v1_a1_fit <- function() {
  structure(list(
    a = 4,
    b = 2.99573227355399,
    J = 20,
    target = list(mu_K = 5, var_K = 8, type = "moments"),
    method = "A1",
    status = "success",
    scaling = list(cJ = log(20), rule = "legacy"),
    cJ = log(20),
    var_K_used = 8,
    converged = TRUE,
    iterations = 0L,
    fit = NULL,
    diagnostics = NULL,
    trace = NULL
  ), class = "DPprior_fit")
}


make_v1_a2_kl_fit <- function() {
  pmf <- rep(1 / 20, 20L)
  support <- seq_along(pmf)
  mean <- sum(support * pmf)
  variance <- sum((support - mean)^2 * pmf)
  structure(list(
    a = 1.5,
    b = 0.7,
    J = 20L,
    target = list(
      type = "chisq", pmf = pmf, mu_K = mean, var_K = variance,
      df = 6.25, scale = 0.8,
      mu_K_discrete = mean, var_K_discrete = variance
    ),
    method = "A2-KL",
    status = "success",
    converged = TRUE,
    iterations = 8L,
    termination = "optim_converged",
    fit = list(mu_K = 5, var_K = 8, kl = 0.02, residual = 0.02),
    diagnostics = list(M = 80L, fallback_used = FALSE),
    trace = data.frame(evaluation = 1L, kl = 0.02)
  ), class = "DPprior_fit")
}


make_v1_wrapper_fit <- function(method = "A2-MN",
                                include_diagnostics = FALSE) {
  raw <- unclass(make_v1_a2_mn_fit())
  raw$target <- list(
    mu_K = 5, var_K = 8, var_K_used = 8,
    confidence = NULL, type = "moments"
  )
  raw$method <- method
  raw$solver_diagnostics <- raw$diagnostics
  raw$diagnostics <- NULL
  raw <- raw[c(
    "a", "b", "J", "target", "method", "status", "converged",
    "iterations", "termination", "fit", "solver_diagnostics", "trace"
  )]
  if (include_diagnostics) {
    raw$diagnostics <- make_v1_diagnostics()
  }
  structure(raw, class = "DPprior_fit")
}


make_v1_dual_fit <- function() {
  fit <- make_v1_a2_mn_fit()
  raw <- unclass(fit)
  raw$target <- list(
    mu_K = 5, var_K = 8, var_K_used = 8,
    confidence = NULL, type = "moments"
  )
  raw$method <- "dual-anchor"
  raw$solver_diagnostics <- raw$diagnostics
  raw$diagnostics <- NULL
  raw$dual_anchor <- list(
    w1_target = list(prob = list(threshold = 0.5, value = 0.3)),
    lambda = 1,
    loss_type = "relative",
    w1_achieved = list(
      mean = 0.42, prob_gt_50 = 0.368, prob_gt_90 = 0.13
    ),
    K_loss = 0,
    init = list(a = raw$a, b = raw$b),
    note = "lambda = 1 returns the K-only solution"
  )
  raw <- raw[c(
    "a", "b", "J", "target", "method", "status", "converged",
    "iterations", "termination", "fit", "solver_diagnostics", "trace",
    "dual_anchor"
  )]
  structure(raw, class = "DPprior_fit")
}


make_v1_diagnostics <- function() {
  structure(list(
    J = 20L,
    a = 1.27223464009856,
    b = 0.58075377814199,
    alpha = list(
      mean = 2.19, sd = 1.94, cv = 0.887, median = 1.65,
      quantiles = c(`5%` = 0.19, `50%` = 1.65, `95%` = 6.03)
    ),
    K = list(
      mean = 5, var = 8, sd = sqrt(8), mode = 4L, median = 5L,
      quantiles = c(`5%` = 1L, `50%` = 5L, `95%` = 10L),
      pmf = rep(0.05, 20L)
    ),
    weights = list(
      mean = 0.42, median = 0.343,
      quantiles = c(`5%` = 0.024, `50%` = 0.343, `95%` = 0.996),
      prob_exceeds = c(`0.5` = 0.368, `0.9` = 0.13),
      dominance_risk = "moderate"
    ),
    coclustering = list(
      mean = 0.42, var = 0.073, sd = 0.269,
      interpretation = "Moderate prior co-clustering"
    ),
    warnings = character()
  ), class = "DPprior_diagnostics")
}


capture_upgrade_conditions <- function(expr) {
  warnings <- list()
  terminal <- NULL
  value <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        warnings[[length(warnings) + 1L]] <<- w
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      terminal <<- e
      NULL
    }
  )
  list(value = value, warnings = warnings, terminal = terminal)
}


has_recursive_name <- function(x, field) {
  if (!is.list(x)) {
    return(FALSE)
  }
  field %in% names(x) || any(vapply(
    unname(x), has_recursive_name, logical(1), field = field
  ))
}


test_that("schema detection recognizes only exact frozen v1.1 shapes", {
  detect <- DPprior:::.dpprior_detect_schema
  fits <- list(
    make_v1_a1_fit(),
    make_v1_a2_mn_fit(),
    make_v1_a2_kl_fit(),
    make_v1_wrapper_fit("A1"),
    make_v1_wrapper_fit("A2-MN"),
    make_v1_wrapper_fit("A2-MN", include_diagnostics = TRUE),
    make_v1_dual_fit()
  )
  diagnostics <- make_v1_diagnostics()

  for (fit in fits) {
    expect_identical(detect(fit), "DPprior/1.1/fit")
  }
  expect_identical(detect(diagnostics), "DPprior/1.1/diagnostics")

  fit <- fits[[2L]]
  extra <- fit
  extra$unexpected <- TRUE
  expect_identical(detect(extra), "unknown")

  partial <- fit
  names(partial)[names(partial) == "method"] <- "meth"
  expect_identical(detect(partial), "unknown")

  duplicate <- fit
  names(duplicate)[2L] <- "a"
  expect_identical(detect(duplicate), "malformed-schema")
})


test_that("schema records reject dimensions and custom scalar classes", {
  detect <- DPprior:::.dpprior_detect_schema
  canonical <- list(
    schema = list(name = "dpprior.result", version = 1L)
  )
  expect_identical(detect(canonical), "dpprior.result/1")

  matrix_name <- canonical
  matrix_name$schema$name <- matrix("dpprior.result", 1L, 1L)
  expect_identical(detect(matrix_name), "malformed-schema")

  array_version <- canonical
  array_version$schema$version <- array(1L, dim = 1L)
  expect_identical(detect(array_version), "malformed-schema")

  classed_version <- canonical
  classed_version$schema$version <- structure(1L, class = "forged_version")
  expect_identical(detect(classed_version), "malformed-schema")

  double_version <- canonical
  double_version$schema$version <- 1
  expect_identical(detect(double_version), "malformed-schema")

  attributed_schema <- canonical
  attr(attributed_schema$schema, "forged") <- TRUE
  expect_identical(detect(attributed_schema), "malformed-schema")
})


test_that("schema detection rejects pairlists and extra legacy attributes", {
  detect <- DPprior:::.dpprior_detect_schema
  outer_pairlist <- pairlist(
    schema = list(name = "dpprior.result", version = 1L)
  )
  nested_pairlist <- list(
    schema = pairlist(name = "dpprior.result", version = 1L)
  )
  attributed <- make_v1_a2_mn_fit()
  attr(attributed, "forged") <- TRUE
  nested_attributed <- make_v1_a2_mn_fit()
  attr(nested_attributed$target, "forged") <- TRUE

  expect_identical(detect(outer_pairlist), "unknown")
  expect_identical(detect(nested_pairlist), "malformed-schema")
  expect_identical(detect(attributed), "unknown")
  expect_identical(detect(nested_attributed), "unknown")
})


test_that("legacy detection bypasses forged extraction methods", {
  fit <- make_v1_a2_mn_fit()
  old_dollar <- get0("$.DPprior_fit", envir = .GlobalEnv, inherits = FALSE)
  old_bracket <- get0("[[.DPprior_fit", envir = .GlobalEnv, inherits = FALSE)
  old_names <- get0("names.DPprior_fit", envir = .GlobalEnv, inherits = FALSE)
  old_dim <- get0("dim.DPprior_fit", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (is.null(old_dollar)) {
      rm("$.DPprior_fit", envir = .GlobalEnv)
    } else {
      assign("$.DPprior_fit", old_dollar, envir = .GlobalEnv)
    }
    if (is.null(old_bracket)) {
      rm("[[.DPprior_fit", envir = .GlobalEnv)
    } else {
      assign("[[.DPprior_fit", old_bracket, envir = .GlobalEnv)
    }
    if (is.null(old_names)) {
      rm("names.DPprior_fit", envir = .GlobalEnv)
    } else {
      assign("names.DPprior_fit", old_names, envir = .GlobalEnv)
    }
    if (is.null(old_dim)) {
      rm("dim.DPprior_fit", envir = .GlobalEnv)
    } else {
      assign("dim.DPprior_fit", old_dim, envir = .GlobalEnv)
    }
  }, add = TRUE)
  assign("$.DPprior_fit", function(...) stop("forged dollar dispatch"),
         envir = .GlobalEnv)
  assign("[[.DPprior_fit", function(...) stop("forged bracket dispatch"),
         envir = .GlobalEnv)
  assign("names.DPprior_fit", function(...) stop("forged names dispatch"),
         envir = .GlobalEnv)
  assign("dim.DPprior_fit", function(...) stop("forged dim dispatch"),
         envir = .GlobalEnv)

  expect_identical(
    DPprior:::.dpprior_detect_schema(fit), "DPprior/1.1/fit"
  )

  guarded <- capture_upgrade_conditions(DPprior:::.dpprior_require_schema(
    fit, kind = "fit", allow_legacy = FALSE
  ))
  expect_s3_class(guarded$terminal, "dpprior_serialization_error")
  expect_identical(guarded$terminal$code, "legacy_upgrade_required")
})


test_that("require guard fails closed with typed migration guidance", {
  result <- capture_upgrade_conditions(
    DPprior:::.dpprior_require_schema(
      make_v1_a2_mn_fit(), kind = "fit", allow_legacy = FALSE
    )
  )
  expect_length(result$warnings, 0L)
  expect_s3_class(result$terminal, "dpprior_legacy_object_error")
  expect_s3_class(result$terminal, "dpprior_serialization_error")
  expect_s3_class(result$terminal, "dpprior_error")
  expect_s3_class(result$terminal, "dpprior_condition")
  expect_identical(result$terminal$code, "legacy_upgrade_required")
  expect_identical(result$terminal$detected_schema, "DPprior/1.1/fit")
  expect_match(result$terminal$upgrade_action, "upgrade_DPprior_object", fixed = TRUE)
})


test_that("require guard upgrades legacy objects only after explicit opt in", {
  result <- capture_upgrade_conditions(
    DPprior:::.dpprior_require_schema(
      make_v1_a2_mn_fit(), kind = "fit", allow_legacy = TRUE
    )
  )

  expect_null(result$terminal)
  expect_length(result$warnings, 1L)
  expect_s3_class(result$warnings[[1L]], "dpprior_legacy_object_warning")
  expect_s3_class(result$warnings[[1L]], "dpprior_condition")
  expect_identical(result$value$schema, list(
    name = "dpprior.result", version = 1L
  ))
  expect_identical(result$value$mode, "a2_moment")
  expect_identical(result$value$status, "approximate")
  expect_false(result$value$usable)
  expect_false(result$value$verified)
})


test_that("guard controls reject matrices arrays and custom classes", {
  bad <- list(
    matrix(FALSE, 1L, 1L),
    array(FALSE, dim = 1L),
    structure(FALSE, class = "forged_logical")
  )
  for (value in bad) {
    condition <- tryCatch(
      DPprior:::.dpprior_require_schema(
        make_v1_a2_mn_fit(), kind = "fit", allow_legacy = value
      ),
      error = identity
    )
    expect_s3_class(condition, "dpprior_upgrade_control_error")
    expect_identical(condition$code, "type")
  }
})


test_that("detection and guard require an exact ordinary kind scalar", {
  bad_kinds <- list(
    matrix("fit", 1L, 1L),
    array("fit", dim = 1L),
    structure("fit", class = "forged_kind"),
    "fi",
    c("fit", "result")
  )
  for (kind in bad_kinds) {
    detected <- tryCatch(
      DPprior:::.dpprior_detect_schema(make_v1_a2_mn_fit(), kind),
      error = identity
    )
    guarded <- tryCatch(
      DPprior:::.dpprior_require_schema(
        make_v1_a2_mn_fit(), kind = kind, allow_legacy = FALSE
      ),
      error = identity
    )
    expect_s3_class(detected, "dpprior_schema_guard_error")
    expect_s3_class(guarded, "dpprior_schema_guard_error")
    expect_true(detected$code %in% c("type", "choice"))
    expect_true(guarded$code %in% c("type", "choice"))
  }
})


test_that("all supported v1 fit shapes migrate conservatively", {
  cases <- list(
    direct_a1 = list(make_v1_a1_fit(), "a1_proxy", FALSE),
    direct_a2_mn = list(make_v1_a2_mn_fit(), "a2_moment", FALSE),
    direct_a2_kl = list(make_v1_a2_kl_fit(), "a2_kl", FALSE),
    wrapper_a1 = list(make_v1_wrapper_fit("A1"), "a1_proxy", TRUE),
    wrapper_a2_mn = list(
      make_v1_wrapper_fit("A2-MN"), "a2_moment", FALSE
    ),
    wrapper_with_diagnostics = list(
      make_v1_wrapper_fit("A2-MN", include_diagnostics = TRUE),
      "a2_moment", FALSE
    ),
    dual = list(make_v1_dual_fit(), "dual_legacy", TRUE)
  )

  for (case in cases) {
    result <- capture_upgrade_conditions(DPprior:::upgrade_DPprior_object(
      case[[1L]], verify = FALSE, allow_legacy = TRUE
    ))
    expect_null(result$terminal)
    expect_length(result$warnings, 1L)
    expect_s3_class(result$warnings[[1L]], "dpprior_legacy_object_warning")
    expect_identical(result$value$mode, case[[2L]])
    expect_identical(result$value$status, "approximate")
    expect_identical(result$value$usable, case[[3L]])
    expect_false(result$value$verified)
    expect_identical(
      DPprior:::.dpprior_validate_object(result$value), result$value
    )
  }
})


test_that("successful fit migration warns once and never trusts v1 status", {
  source <- make_v1_a2_mn_fit()
  before <- serialize(source, NULL, version = 3L)
  result <- capture_upgrade_conditions(
    DPprior:::upgrade_DPprior_object(source)
  )

  expect_null(result$terminal)
  expect_length(result$warnings, 1L)
  expect_s3_class(result$warnings[[1L]], "dpprior_legacy_object_warning")
  expect_s3_class(result$warnings[[1L]], "dpprior_deprecated_warning")
  expect_s3_class(result$warnings[[1L]], "dpprior_warning")
  expect_identical(result$warnings[[1L]]$code, "legacy_object_upgraded")
  expect_identical(result$warnings[[1L]]$source_schema, "DPprior/1.1/fit")
  expect_identical(result$warnings[[1L]]$target_schema, "dpprior.result/1")
  expect_true(
    "public_candidate_quarantined_pending_refit" %in%
      result$warnings[[1L]]$losses
  )
  expect_identical(result$value$schema, list(
    name = "dpprior.result", version = 1L
  ))
  expect_identical(result$value$mode, "a2_moment")
  expect_identical(result$value$status, "approximate")
  expect_false(result$value$usable)
  expect_false(result$value$verified)
  expect_false(result$value$verification$performed)
  expect_false(result$value$verification$passed)
  expect_length(result$value$computation$attempts, 0L)
  expect_true(all(vapply(
    result$value$computation$orders[1:4], is.null, logical(1)
  )))
  expect_false(
    result$value$compatibility$views$fixed_candidate_recomputation$decision_ready
  )
  expect_identical(serialize(source, NULL, version = 3L), before)
  expect_identical(
    DPprior:::.dpprior_validate_object(result$value), result$value
  )
})


test_that("A2 candidates are public-quarantined with explicit refit provenance", {
  cases <- list(
    make_v1_a2_mn_fit(),
    make_v1_a2_kl_fit(),
    make_v1_wrapper_fit("A2-MN", include_diagnostics = TRUE)
  )

  for (source in cases) {
    result <- suppressWarnings(DPprior:::upgrade_DPprior_object(
      source, allow_legacy = TRUE
    ))
    source_view <- result$compatibility$views$source
    audit <- result$compatibility$views$fixed_candidate_recomputation

    expect_null(result$parameters)
    expect_length(result$achieved, 0L)
    expect_length(result$residuals, 0L)
    expect_false(result$usable)
    expect_false(result$verified)
    expect_false(result$verification$performed)
    expect_true(source_view$public_candidate_quarantined)
    expect_identical(source_view$required_action, "refit_with_current_API")
    expect_identical(
      source_view$source_parameters[c("a", "b", "J")],
      unclass(source)[c("a", "b", "J")]
    )
    expect_identical(
      audit$candidate[c("a", "b")],
      unclass(source)[c("a", "b")]
    )
    expect_false(audit$decision_ready)
    expect_true(
      "public_candidate_quarantined_pending_refit" %in%
        result$provenance$migration$missing_evidence
    )
  }
})


test_that("nested legacy diagnostic bundles have an explicit loss disposition", {
  source <- make_v1_wrapper_fit(
    "A2-MN", include_diagnostics = TRUE
  )
  result <- suppressWarnings(DPprior:::upgrade_DPprior_object(source))

  expect_identical(
    result$compatibility$views$source$discarded_legacy_fields,
    c("diagnostics", "solver_diagnostics")
  )
  expect_true(result$compatibility$views$source$public_candidate_quarantined)
  expect_false("diagnostics" %in% names(result))
})


test_that("allow_legacy is limited to A1 and retained legacy Dual-Anchor", {
  a2 <- suppressWarnings(DPprior:::upgrade_DPprior_object(
    make_v1_a2_mn_fit(), allow_legacy = TRUE
  ))
  a1_default <- suppressWarnings(DPprior:::upgrade_DPprior_object(
    make_v1_a1_fit()
  ))
  a1_opt_in <- suppressWarnings(DPprior:::upgrade_DPprior_object(
    make_v1_a1_fit(), allow_legacy = TRUE
  ))
  dual <- suppressWarnings(DPprior:::upgrade_DPprior_object(
    make_v1_dual_fit(), allow_legacy = TRUE
  ))

  expect_false(a2$usable)
  expect_false(a1_default$usable)
  expect_true(a1_opt_in$usable)
  expect_false(a1_opt_in$verified)
  expect_identical(a1_opt_in$mode, "a1_proxy")
  expect_identical(
    a1_opt_in$verification$selected_snapshot$source,
    "fresh_migration_fixed_candidate_recomputation"
  )
  expect_identical(
    a1_opt_in$verification$selected_snapshot$M,
    a1_opt_in$compatibility$views$fixed_candidate_recomputation$M_selected
  )
  expect_identical(dual$mode, "dual_legacy")
  expect_true(dual$usable)
  expect_false(dual$verified)
  expect_true("legacy" %in% names(dual))
  expect_false("constraint" %in% names(dual))
  expect_false("tradeoff" %in% names(dual))
  dual_core <- unclass(dual)
  dual_core$compatibility <- NULL
  expect_false(has_recursive_name(
    dual_core, "constraint_satisfied"
  ))
})


test_that("verify=FALSE retains no fresh fit audit", {
  result <- suppressWarnings(DPprior:::upgrade_DPprior_object(
    make_v1_a2_mn_fit(), verify = FALSE
  ))
  expect_false(
    "fixed_candidate_recomputation" %in%
      names(result$compatibility$views)
  )
  expect_false(result$verification$performed)
  expect_false(result$verified)
})


test_that("diagnostics are freshly recomputed without dominance category", {
  source <- make_v1_diagnostics()
  before <- serialize(source, NULL, version = 3L)
  result <- capture_upgrade_conditions(
    DPprior:::upgrade_DPprior_object(source, verify = TRUE)
  )

  expect_null(result$terminal)
  expect_length(result$warnings, 1L)
  expect_s3_class(result$warnings[[1L]], "dpprior_legacy_object_warning")
  expect_s3_class(result$warnings[[1L]], "dpprior_deprecated_warning")
  expect_s3_class(result$warnings[[1L]], "dpprior_warning")
  expect_identical(
    result$warnings[[1L]]$source_schema,
    "DPprior/1.1/diagnostics"
  )
  expect_identical(result$value$mode, "prior_diagnostics")
  expect_identical(result$value$status, "converged")
  expect_true(result$value$usable)
  expect_true(result$value$verified)
  expect_true(result$value$verification$performed)
  expect_true(result$value$verification$passed)
  component_check <-
    result$value$verification$components$component_aggregation
  expect_true(component_check$passed)
  expect_identical(
    component_check$value,
    stats::setNames(rep(TRUE, 4L), c("alpha", "K", "weights", "coclustering"))
  )
  expect_identical(component_check$value, component_check$reference)
  expect_identical(component_check$source, "fresh_component_specific_checks")
  expect_false(has_recursive_name(result$value, "dominance_risk"))
  expect_identical(
    names(result$value$diagnostics$weights),
    c("status", "usable", "verified", "mean")
  )
  expect_identical(
    vapply(
      result$value$computation$attempts,
      function(attempt) attempt$id,
      character(1)
    ),
    paste0("diagnostic-", c("alpha", "K", "weights", "coclustering"))
  )
  expect_identical(
    result$value$compatibility$views$source$discarded_legacy_fields,
    "weights$dominance_risk"
  )
  expect_match(
    result$value$compatibility$views$source$source_digest,
    "^[0-9a-f]{32}$"
  )
  expect_identical(serialize(source, NULL, version = 3L), before)
  expect_identical(
    DPprior:::.dpprior_validate_object(result$value), result$value
  )
})


test_that("failed migration emits only a typed terminal condition", {
  malformed <- make_v1_a2_mn_fit()
  malformed$a <- structure(1, class = "forged_parameter")
  result <- capture_upgrade_conditions(
    DPprior:::upgrade_DPprior_object(malformed)
  )
  expect_length(result$warnings, 0L)
  expect_s3_class(result$terminal, "dpprior_legacy_object_error")
  expect_s3_class(result$terminal, "dpprior_serialization_error")
  expect_identical(result$terminal$code, "legacy_invalid_parameter")

  no_recompute <- capture_upgrade_conditions(
    DPprior:::upgrade_DPprior_object(
      make_v1_diagnostics(), verify = FALSE
    )
  )
  expect_length(no_recompute$warnings, 0L)
  expect_s3_class(no_recompute$terminal, "dpprior_legacy_object_error")
  expect_identical(
    no_recompute$terminal$code, "diagnostics_recomputation_required"
  )

  malformed_dual <- make_v1_dual_fit()
  malformed_dual$dual_anchor$w1_achieved <- 1
  raw_guard <- capture_upgrade_conditions(
    DPprior:::upgrade_DPprior_object(
      malformed_dual, allow_legacy = TRUE
    )
  )
  expect_length(raw_guard$warnings, 0L)
  expect_s3_class(raw_guard$terminal, "dpprior_legacy_object_error")
  expect_s3_class(raw_guard$terminal, "dpprior_serialization_error")
  expect_identical(
    raw_guard$terminal$code, "legacy_missing_weight_achieved"
  )

  normalization_failure <- make_v1_a2_mn_fit()
  normalization_failure$target$type <- environment()
  normalized <- capture_upgrade_conditions(
    DPprior:::upgrade_DPprior_object(normalization_failure)
  )
  expect_length(normalized$warnings, 0L)
  expect_s3_class(normalized$terminal, "dpprior_legacy_object_error")
  expect_s3_class(normalized$terminal, "dpprior_serialization_error")
  expect_identical(
    normalized$terminal$code, "legacy_normalization_failed"
  )
  expect_s3_class(normalized$terminal$cause, "dpprior_schema_error")
})


test_that("fit audit warnings become one terminal serialization condition", {
  testthat::local_mocked_bindings(
    .dpprior_upgrade_fixed_candidate_audit = function(...) {
      warning("injected recomputation warning")
      NULL
    },
    .package = "DPprior"
  )
  result <- capture_upgrade_conditions(
    DPprior:::upgrade_DPprior_object(make_v1_a2_mn_fit())
  )

  expect_length(result$warnings, 0L)
  expect_s3_class(result$terminal, "dpprior_legacy_object_error")
  expect_s3_class(result$terminal, "dpprior_serialization_error")
  expect_identical(result$terminal$code, "legacy_recomputation_warning")
  expect_s3_class(result$terminal$cause, "simpleWarning")
})


test_that("unknown malformed and future schemas fail closed", {
  objects <- list(
    structure(list(a = 1, b = 1, J = 10), class = "DPprior_fit"),
    list(schema = list(name = "dpprior.result", version = 99L)),
    list(schema = list(name = "future.result", version = 1L)),
    list(schema = list(name = "dpprior.result", version = matrix(1L))),
    pairlist(schema = list(name = "dpprior.result", version = 1L)),
    list(schema = pairlist(name = "dpprior.result", version = 1L))
  )
  for (object in objects) {
    result <- capture_upgrade_conditions(
      DPprior:::upgrade_DPprior_object(object)
    )
    expect_length(result$warnings, 0L)
    expect_s3_class(result$terminal, "dpprior_legacy_object_error")
    expect_s3_class(result$terminal, "dpprior_serialization_error")
    expect_identical(result$terminal$code, "legacy_schema_unrecognized")
  }
})


test_that("current canonical objects are validated idempotently and quietly", {
  migrated <- suppressWarnings(DPprior:::upgrade_DPprior_object(
    make_v1_a2_mn_fit()
  ))
  result <- capture_upgrade_conditions(
    DPprior:::upgrade_DPprior_object(migrated)
  )
  expect_null(result$terminal)
  expect_length(result$warnings, 0L)
  expect_identical(result$value, migrated)

  malformed <- migrated
  malformed$verified <- TRUE
  bad <- capture_upgrade_conditions(
    DPprior:::upgrade_DPprior_object(malformed)
  )
  expect_length(bad$warnings, 0L)
  expect_s3_class(bad$terminal, "dpprior_schema_error")
})


test_that("canonical upgrades preserve RDS round trips", {
  objects <- list(
    suppressWarnings(DPprior:::upgrade_DPprior_object(
      make_v1_a2_mn_fit()
    )),
    suppressWarnings(DPprior:::upgrade_DPprior_object(
      make_v1_diagnostics()
    ))
  )
  for (object in objects) {
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path), add = TRUE)
    saveRDS(object, path, version = 3L)
    restored <- readRDS(path)
    expect_identical(restored, object)
    expect_identical(DPprior:::.dpprior_validate_object(restored), restored)
  }
})


test_that("public upgrade controls reject shaped or classed scalars", {
  bad_flags <- list(
    matrix(TRUE, 1L, 1L),
    array(TRUE, dim = 1L),
    structure(TRUE, class = "forged_logical")
  )
  for (flag in bad_flags) {
    condition <- tryCatch(
      DPprior:::upgrade_DPprior_object(
        make_v1_a2_mn_fit(), verify = flag
      ),
      error = identity
    )
    expect_s3_class(condition, "dpprior_upgrade_control_error")
    expect_identical(condition$code, "type")
  }

  bad_orders <- list(
    matrix(160L, 1L, 1L),
    array(160L, dim = 1L),
    structure(160L, class = "forged_order")
  )
  for (order in bad_orders) {
    condition <- tryCatch(
      DPprior:::upgrade_DPprior_object(
        make_v1_a2_mn_fit(), M_verify = order
      ),
      error = identity
    )
    expect_s3_class(condition, "dpprior_upgrade_control_error")
  }
})
