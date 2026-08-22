# =============================================================================
# Canonical schema detection and conservative legacy-object migration
# =============================================================================

.DPPRIOR_RESULT_SCHEMA_V1 <- "dpprior.result/1"
.DPPRIOR_TARGET_SCHEMA_V1 <- "dpprior.target/1"
.DPPRIOR_WEIGHT_TARGET_SCHEMA_V1 <- "dpprior.weight-target/1"
.DPPRIOR_LEGACY_FIT_V1_1 <- "DPprior/1.1/fit"
.DPPRIOR_LEGACY_DIAGNOSTICS_V1_1 <- "DPprior/1.1/diagnostics"
.DPPRIOR_UNKNOWN_SCHEMA <- "unknown"


.dpprior_upgrade_exact_names <- function(x, expected) {
  typeof(x) == "list" && is.list(x) && !is.object(x) && is.null(dim(x)) &&
    identical(names(attributes(x)), "names") &&
    !anyDuplicated(names(x)) && identical(names(x), expected)
}


.dpprior_upgrade_scalar_character <- function(x) {
  is.character(x) && !is.object(x) && is.null(dim(x)) &&
    length(x) == 1L && !is.na(x) && nzchar(x)
}


.dpprior_upgrade_scalar_logical <- function(x) {
  is.logical(x) && !is.object(x) && is.null(dim(x)) &&
    length(x) == 1L && !is.na(x)
}


.dpprior_upgrade_schema_label <- function(schema) {
  if (!.dpprior_upgrade_exact_names(schema, c("name", "version")) ||
      !.dpprior_upgrade_scalar_character(schema[["name", exact = TRUE]]) ||
      !is.integer(schema[["version", exact = TRUE]]) ||
      is.object(schema[["version", exact = TRUE]]) ||
      !is.null(dim(schema[["version", exact = TRUE]])) ||
      length(schema[["version", exact = TRUE]]) != 1L ||
      is.na(schema[["version", exact = TRUE]]) ||
      !is.finite(schema[["version", exact = TRUE]])) {
    return("malformed-schema")
  }
  sprintf(
    "%s/%d", schema[["name", exact = TRUE]],
    as.integer(schema[["version", exact = TRUE]])
  )
}


.dpprior_upgrade_validate_kind <- function(kind) {
  allowed <- c(
    "auto", "result", "fit", "diagnostics", "target", "weight-target"
  )
  if (!.dpprior_upgrade_scalar_character(kind)) {
    .dpprior_abort_invalid(
      "kind must be one ordinary non-missing character scalar",
      c("dpprior_schema_guard_error", "dpprior_type_error"),
      "kind", kind, paste(allowed, collapse = ", "), "type"
    )
  }
  if (!kind %in% allowed) {
    .dpprior_abort_invalid(
      sprintf("Unsupported schema kind '%s'", kind),
      c("dpprior_schema_guard_error", "dpprior_choice_error"),
      "kind", kind, paste(allowed, collapse = ", "), "choice"
    )
  }
  kind
}


.dpprior_is_v1_1_fit_shape <- function(x) {
  attribute_names <- names(attributes(x))
  if (typeof(x) != "list" || !is.list(x) ||
      isS4(x) || !identical(class(x), "DPprior_fit") ||
      anyDuplicated(attribute_names) ||
      !identical(sort(attribute_names), c("class", "names"))) {
    return(FALSE)
  }
  raw <- unclass(x)
  if (!is.null(dim(raw)) || anyDuplicated(names(raw))) {
    return(FALSE)
  }

  direct_a1 <- c(
    "a", "b", "J", "target", "method", "status", "scaling", "cJ",
    "var_K_used", "converged", "iterations", "fit", "diagnostics", "trace"
  )
  direct_solver <- c(
    "a", "b", "J", "target", "method", "status", "converged",
    "iterations", "termination", "fit", "diagnostics", "trace"
  )
  wrapper <- c(
    "a", "b", "J", "target", "method", "status", "converged",
    "iterations", "termination", "fit", "solver_diagnostics", "trace"
  )
  wrapper_diagnostics <- c(wrapper, "diagnostics")
  dual <- c(wrapper, "dual_anchor")
  moment_target <- c("mu_K", "var_K", "type")
  wrapper_target <- c(
    "mu_K", "var_K", "var_K_used", "confidence", "type"
  )
  kl_target <- c(
    "type", "pmf", "mu_K", "var_K", "df", "scale",
    "mu_K_discrete", "var_K_discrete"
  )

  if (identical(names(raw), direct_a1)) {
    return(identical(raw[["method", exact = TRUE]], "A1") &&
             .dpprior_upgrade_exact_names(
               raw[["target", exact = TRUE]], moment_target
             ))
  }
  if (identical(names(raw), direct_solver)) {
    if (identical(raw[["method", exact = TRUE]], "A2-MN")) {
      return(.dpprior_upgrade_exact_names(
        raw[["target", exact = TRUE]], moment_target
      ))
    }
    if (identical(raw[["method", exact = TRUE]], "A2-KL")) {
      return(.dpprior_upgrade_exact_names(
        raw[["target", exact = TRUE]], kl_target
      ))
    }
    return(FALSE)
  }
  if (identical(names(raw), wrapper) ||
      identical(names(raw), wrapper_diagnostics)) {
    method <- raw[["method", exact = TRUE]]
    return((identical(method, "A1") || identical(method, "A2-MN")) &&
             .dpprior_upgrade_exact_names(
               raw[["target", exact = TRUE]], wrapper_target
             ))
  }
  if (identical(names(raw), dual)) {
    return(identical(raw[["method", exact = TRUE]], "dual-anchor") &&
             .dpprior_upgrade_exact_names(
               raw[["target", exact = TRUE]], wrapper_target
             ) &&
             .dpprior_upgrade_exact_names(
               raw[["dual_anchor", exact = TRUE]],
               c(
                 "w1_target", "lambda", "loss_type", "w1_achieved",
                 "K_loss", "init", "note"
               )
             ))
  }
  FALSE
}


.dpprior_is_v1_1_diagnostics_shape <- function(x) {
  attribute_names <- names(attributes(x))
  if (typeof(x) != "list" || !is.list(x) ||
      isS4(x) || !identical(class(x), "DPprior_diagnostics") ||
      anyDuplicated(attribute_names) ||
      !identical(sort(attribute_names), c("class", "names"))) {
    return(FALSE)
  }
  raw <- unclass(x)
  if (!is.null(dim(raw)) || anyDuplicated(names(raw))) {
    return(FALSE)
  }
  .dpprior_upgrade_exact_names(
      raw, c("J", "a", "b", "alpha", "K", "weights", "coclustering", "warnings")
    ) &&
    .dpprior_upgrade_exact_names(
      raw[["alpha", exact = TRUE]],
      c("mean", "sd", "cv", "median", "quantiles")
    ) &&
    .dpprior_upgrade_exact_names(
      raw[["K", exact = TRUE]],
      c("mean", "var", "sd", "mode", "median", "quantiles", "pmf")
    ) &&
    .dpprior_upgrade_exact_names(
      raw[["weights", exact = TRUE]],
      c("mean", "median", "quantiles", "prob_exceeds", "dominance_risk")
    ) &&
    .dpprior_upgrade_exact_names(
      raw[["coclustering", exact = TRUE]],
      c("mean", "var", "sd", "interpretation")
    )
}


# Detect only canonical schemas and frozen v1.1 shapes. Inheritance alone and
# absence of a schema marker are never sufficient evidence of a legacy shape.
.dpprior_detect_schema <- function(x, kind = "auto") {
  kind <- .dpprior_upgrade_validate_kind(kind)
  if (typeof(x) != "list" || !is.list(x)) {
    return(.DPPRIOR_UNKNOWN_SCHEMA)
  }
  if (isS4(x)) {
    return("malformed-schema")
  }
  raw <- if (is.object(x)) unclass(x) else x
  if (!is.null(dim(raw)) || anyDuplicated(names(raw))) {
    return("malformed-schema")
  }
  if ("schema" %in% names(raw)) {
    return(.dpprior_upgrade_schema_label(raw[["schema", exact = TRUE]]))
  }
  if ("schema_version" %in% names(raw)) {
    version <- raw[["schema_version", exact = TRUE]]
    value <- if (is.numeric(version) && !is.object(version) &&
                    is.null(dim(version)) && length(version) == 1L &&
                    !is.na(version) && is.finite(version)) {
      as.character(version)
    } else {
      "malformed"
    }
    return(sprintf("legacy-schema_version/%s", value))
  }
  if (kind %in% c("auto", "result", "fit") &&
      .dpprior_is_v1_1_fit_shape(x)) {
    return(.DPPRIOR_LEGACY_FIT_V1_1)
  }
  if (kind %in% c("auto", "result", "diagnostics") &&
      .dpprior_is_v1_1_diagnostics_shape(x)) {
    return(.DPPRIOR_LEGACY_DIAGNOSTICS_V1_1)
  }
  .DPPRIOR_UNKNOWN_SCHEMA
}


.dpprior_abort_serialization <- function(
    message,
    subclass = character(),
    code = "unsupported_schema",
    x = NULL,
    detected_schema = .DPPRIOR_UNKNOWN_SCHEMA,
    expected_schema = .DPPRIOR_RESULT_SCHEMA_V1,
    upgrade_action = NULL,
    cause = NULL,
    call = NULL) {
  stop(.dpprior_new_condition(
    message = message,
    classes = c(
      subclass, "dpprior_serialization_error", "dpprior_error", "error"
    ),
    call = call,
    code = code,
    object_class = if (is.null(x)) character() else class(x),
    detected_schema = detected_schema,
    expected_schema = expected_schema,
    upgrade_action = upgrade_action,
    cause = cause
  ))
}


.dpprior_abort_legacy_object <- function(
    x,
    kind,
    detected_schema,
    action = paste(
      "Call upgrade_DPprior_object(x, verify = TRUE), then saveRDS()",
      "the returned object."
    ),
    code = "legacy_upgrade_required",
    call = NULL) {
  .dpprior_abort_serialization(
    message = sprintf(
      paste(
        "%s uses unsupported serialized schema '%s'.",
        "%s"
      ),
      kind, detected_schema, action
    ),
    subclass = "dpprior_legacy_object_error",
    code = code,
    x = x,
    detected_schema = detected_schema,
    expected_schema = .DPPRIOR_RESULT_SCHEMA_V1,
    upgrade_action = action,
    call = call
  )
}


.dpprior_warn_legacy_upgrade <- function(
    source_schema,
    losses,
    replacement = paste(
      "Save the returned object and refit with the current API before any",
      "decision-readiness claim."
    ),
    call = NULL) {
  warning(.dpprior_new_condition(
    message = sprintf(
      "Upgraded '%s' conservatively; %s",
      source_schema, replacement
    ),
    classes = c(
      "dpprior_legacy_object_warning", "dpprior_deprecated_warning",
      "dpprior_warning", "warning"
    ),
    call = call,
    code = "legacy_object_upgraded",
    source_schema = source_schema,
    target_schema = .DPPRIOR_RESULT_SCHEMA_V1,
    losses = losses,
    replacement = replacement,
    since = "2.0.0",
    removal_earliest = "not_applicable"
  ))
  invisible(NULL)
}


.dpprior_upgrade_validate_flag <- function(x, argument) {
  if (!.dpprior_upgrade_scalar_logical(x)) {
    .dpprior_abort_invalid(
      sprintf("%s must be TRUE or FALSE", argument),
      c("dpprior_upgrade_control_error", "dpprior_type_error"),
      argument, x, "one non-missing logical value", "type"
    )
  }
  x
}


.dpprior_upgrade_kind_matches <- function(x, kind) {
  if (kind %in% c("auto", "result")) {
    return(TRUE)
  }
  if (kind == "fit") {
    return(identical(unclass(x)[["object_type", exact = TRUE]], "fit"))
  }
  if (kind == "diagnostics") {
    return(identical(
      unclass(x)[["object_type", exact = TRUE]], "diagnostics"
    ))
  }
  FALSE
}


.dpprior_upgrade_detected_kind_matches <- function(x, detected, kind) {
  if (identical(kind, "auto")) {
    return(TRUE)
  }
  if (identical(detected, .DPPRIOR_RESULT_SCHEMA_V1)) {
    return(.dpprior_upgrade_kind_matches(x, kind))
  }
  if (identical(detected, .DPPRIOR_TARGET_SCHEMA_V1)) {
    return(identical(kind, "target"))
  }
  if (identical(detected, .DPPRIOR_WEIGHT_TARGET_SCHEMA_V1)) {
    return(identical(kind, "weight-target"))
  }
  FALSE
}


# Require a canonical object at a consumer boundary. Legacy migration is never
# implicit unless the caller deliberately opts into the legacy adapter path.
.dpprior_require_schema <- function(x, kind = "result", allow_legacy = FALSE) {
  kind <- .dpprior_upgrade_validate_kind(kind)
  allow_legacy <- .dpprior_upgrade_validate_flag(
    allow_legacy, "allow_legacy"
  )
  detected <- .dpprior_detect_schema(x, kind)

  if (detected %in% c(
    .DPPRIOR_RESULT_SCHEMA_V1,
    .DPPRIOR_TARGET_SCHEMA_V1,
    .DPPRIOR_WEIGHT_TARGET_SCHEMA_V1
  )) {
    validated <- .dpprior_validate_object(x)
    if (!.dpprior_upgrade_detected_kind_matches(x, detected, kind)) {
      .dpprior_abort_serialization(
        sprintf("Expected a canonical %s object", kind),
        "dpprior_schema_kind_error", "schema_kind_mismatch", x,
        detected, .DPPRIOR_RESULT_SCHEMA_V1,
        "Pass an object whose canonical object_type matches the consumer."
      )
    }
    return(validated)
  }

  is_v1 <- detected %in% c(
    .DPPRIOR_LEGACY_FIT_V1_1,
    .DPPRIOR_LEGACY_DIAGNOSTICS_V1_1
  )
  if (is_v1 && allow_legacy) {
    return(upgrade_DPprior_object(
      x, verify = TRUE, M_verify = NULL, allow_legacy = TRUE
    ))
  }
  if (is_v1) {
    .dpprior_abort_legacy_object(x, kind, detected)
  }

  .dpprior_abort_serialization(
    message = sprintf(
      "Object schema '%s' is unknown, malformed, or unsupported", detected
    ),
    subclass = "dpprior_schema_unsupported_error",
    code = "unsupported_schema",
    x = x,
    detected_schema = detected,
    expected_schema = if (kind == "target") {
      .DPPRIOR_TARGET_SCHEMA_V1
    } else if (kind == "weight-target") {
      .DPPRIOR_WEIGHT_TARGET_SCHEMA_V1
    } else {
      .DPPRIOR_RESULT_SCHEMA_V1
    },
    upgrade_action = paste(
      "Use a supported canonical object or recreate the object with the",
      "current DPprior version."
    )
  )
}


.dpprior_upgrade_validate_positive <- function(
    x, path, source_schema = "DPprior/1.1") {
  ok <- is.numeric(x) && !is.object(x) && is.null(dim(x)) &&
    length(x) == 1L && !is.na(x) && is.finite(x) && x > 0
  if (!ok) {
    .dpprior_abort_serialization(
      sprintf("Legacy field '%s' must be one finite positive scalar", path),
      "dpprior_legacy_object_error", "legacy_invalid_parameter",
      detected_schema = source_schema,
      expected_schema = .DPPRIOR_RESULT_SCHEMA_V1,
      upgrade_action = "Recreate or refit the object from validated inputs."
    )
  }
  as.numeric(x)
}


.dpprior_upgrade_validate_J <- function(
    x, source_schema = "DPprior/1.1") {
  ok <- is.numeric(x) && !is.object(x) && is.null(dim(x)) &&
    length(x) == 1L && !is.na(x) && is.finite(x) &&
    x == floor(x) && x >= 1 && x <= .Machine$integer.max
  if (!ok) {
    .dpprior_abort_serialization(
      "Legacy field 'J' must be one positive integer",
      "dpprior_legacy_object_error", "legacy_invalid_sample_size",
      detected_schema = source_schema,
      expected_schema = .DPPRIOR_RESULT_SCHEMA_V1,
      upgrade_action = "Recreate or refit the object from validated inputs."
    )
  }
  as.integer(x)
}


.dpprior_upgrade_validate_M_verify <- function(M_verify, verify) {
  if (is.null(M_verify)) {
    return(NULL)
  }
  if (!isTRUE(verify)) {
    .dpprior_abort_invalid(
      "M_verify is meaningful only when verify=TRUE",
      "dpprior_upgrade_control_error", "M_verify", M_verify,
      "NULL when verify=FALSE", "unused_control"
    )
  }
  ok <- is.numeric(M_verify) && !is.object(M_verify) &&
    is.null(dim(M_verify)) && length(M_verify) == 1L &&
    !is.na(M_verify) && is.finite(M_verify) &&
    M_verify == floor(M_verify) && M_verify >= 1 &&
    M_verify <= .QUADRATURE_MAX_NODES
  if (!ok) {
    .dpprior_abort_invalid(
      sprintf(
        "M_verify must be one integer in [1, %d]",
        .QUADRATURE_MAX_NODES
      ),
      c("dpprior_upgrade_control_error", "dpprior_count_error"),
      "M_verify", M_verify,
      sprintf("integer in [1, %d]", .QUADRATURE_MAX_NODES), "bounds"
    )
  }
  as.integer(M_verify)
}


.dpprior_legacy_digest <- function(x) {
  path <- tempfile("dpprior-legacy-", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(x, path, version = 3L)
  unname(tools::md5sum(path))
}


.dpprior_upgrade_fit_mode <- function(method) {
  switch(
    method,
    A1 = "a1_proxy",
    `A2-MN` = "a2_moment",
    `A2-KL` = "a2_kl",
    `dual-anchor` = "dual_legacy",
    NULL
  )
}


.dpprior_upgrade_source_summary <- function(x, source_schema) {
  raw <- if (is.object(x)) unclass(x) else x
  discarded <- character()
  if (inherits(x, "DPprior_diagnostics") &&
      !is.null(
        raw[["weights", exact = TRUE]][["dominance_risk", exact = TRUE]]
      )) {
    discarded <- "weights$dominance_risk"
  }
  if (inherits(x, "DPprior_fit")) {
    for (field in c("diagnostics", "solver_diagnostics")) {
      if (field %in% names(raw) &&
          !is.null(raw[[field, exact = TRUE]])) {
        discarded <- c(discarded, field)
      }
    }
  }
  list(
    source_schema = source_schema,
    source_package_version = "1.1.0",
    source_class = class(x),
    source_fields = names(raw),
    source_status = raw[["status", exact = TRUE]],
    source_converged = raw[["converged", exact = TRUE]],
    source_parameters = if (all(c("a", "b", "J") %in% names(raw))) {
      list(
        a = raw[["a", exact = TRUE]],
        b = raw[["b", exact = TRUE]],
        J = raw[["J", exact = TRUE]],
        parameterization = "Gamma(shape=a, rate=b)"
      )
    } else {
      NULL
    },
    source_digest = .dpprior_legacy_digest(x),
    digest_method = "base-r-saveRDS-v3-md5",
    discarded_legacy_fields = if (length(discarded)) {
      unique(discarded)
    } else {
      NULL
    }
  )
}


.dpprior_upgrade_missing_evidence <- function(mode) {
  base <- c(
    "requested_controls", "M_requested", "M_selected",
    "M_verification_required", "M_verification_used", "scaling",
    "canonical_attempts", "candidate_selection", "fallback_lineage",
    "independent_verifier_snapshot", "source_commit"
  )
  if (mode == "a1_proxy") {
    return(c(base, "exact_prior_verification"))
  }
  if (mode == "a2_kl") {
    return(c(base, "optimizer_lineage", "target_family_verification"))
  }
  if (mode == "a2_moment") {
    return(c(base, "optimizer_lineage"))
  }
  c(base, "legacy_objective_scaling", "legacy_optimizer_lineage")
}


.dpprior_upgrade_recompute_orders <- function(M_verify) {
  M_selected <- as.integer(.QUAD_NODES_DEFAULT)
  required <- as.integer(.quadrature_verification_required_order(M_selected))
  if (is.null(M_verify)) {
    M_verify <- required
  }
  if (M_verify < required) {
    .dpprior_abort_invalid(
      sprintf(
        "M_verify must be at least %d for migration audit order M=%d",
        required, M_selected
      ),
      c("dpprior_upgrade_control_error", "dpprior_count_error"),
      "M_verify", M_verify, sprintf("integer at least %d", required),
      "insufficient_verification_order"
    )
  }
  list(M_selected = M_selected, M_verify = as.integer(M_verify))
}


.dpprior_upgrade_fit_truth <- function(mode) {
  if (identical(mode, "a2_moment")) {
    controls <- list(
      max_iter = 20L,
      damping = TRUE,
      use_fallback = TRUE,
      tol_step = 1e-10,
      log_bounds = c(-15, 15),
      boundary_tol = 1e-6,
      line_search_max = 20L,
      jacobian_rcond_singular = 1e-12,
      jacobian_rcond_ill = sqrt(.Machine$double.eps),
      fallback = list(
        maxit = 1000L,
        reltol = 1e-12,
        finite_penalty = .Machine$double.xmax / 1024
      ),
      selection_tolerance = 0
    )
    tolerances <- list(
      K_adequacy = list(
        absolute = 1e-8,
        relative = 1e-8,
        scale_formula = "max(abs(target),1)"
      ),
      K_stability = list(
        absolute = 1e-10,
        relative = 1e-8,
        scale_floor = 1
      ),
      step = 1e-10,
      boundary = 1e-6
    )
    return(list(controls = controls, tolerances = tolerances))
  }
  if (identical(mode, "a2_kl")) {
    controls <- list(
      max_iter = 100L,
      optimizer_tol = 1e-6,
      log_bounds = c(-15, 15),
      boundary_tol = 1e-6,
      fallback_max_iter = 100L,
      primary = list(
        maxit = 100L,
        factr = 1e-6 / .Machine$double.eps,
        pgtol = 1e-6
      ),
      fallback = list(
        iter.max = 100L,
        eval.max = 200L,
        rel.tol = 1e-6,
        x.tol = 1e-6
      ),
      fallback_trigger_worse_than_start = 1e-12,
      selection_tolerance = 0
    )
    order <- list(
      pmf_absolute = 1e-10,
      pmf_relative = 1e-8,
      pmf_l1 = 1e-10 + 1e-8,
      direct_moment_absolute = 1e-10,
      direct_moment_relative = 1e-8,
      target_identity_l1 = .TOL_PMF_SUM
    )
    tolerances <- list(
      distribution = list(
        adequacy = list(
          kl = 0.015,
          l1 = 0.11,
          mean_scaled = 0.01,
          variance_scaled = 0.065,
          mean_scale_formula = "max(1,sqrt(target_variance))",
          variance_scale_formula = "max(1,target_variance)"
        ),
        order = order
      ),
      boundary = 1e-6
    )
    return(list(controls = controls, tolerances = tolerances))
  }
  list(controls = list(), tolerances = list())
}


.dpprior_upgrade_setting <- function(method, controls = list()) {
  list(
    method = method,
    controls = controls,
    parameterization = "Gamma(shape=a, rate=b)"
  )
}


.dpprior_upgrade_unknown_orders <- function() {
  .dpprior_new_orders(
    M_requested = NULL,
    M_selected = NULL,
    M_verification_required = NULL,
    M_verification_used = NULL,
    requested_reason = "not_recorded_by_v1.1",
    selected_reason = "not_recorded_by_v1.1",
    verification_required_reason = "not_recorded_by_v1.1",
    verification_used_reason = "not_recorded_by_v1.1"
  )
}


.dpprior_upgrade_computation <- function(method,
                                         controls = list(),
                                         mode = NULL,
                                         M_selected = NULL) {
  setting <- .dpprior_upgrade_setting(method, controls)
  orders <- if (is.null(M_selected)) {
    .dpprior_upgrade_unknown_orders()
  } else {
    .dpprior_new_orders(
      M_requested = as.integer(M_selected),
      M_selected = as.integer(M_selected),
      M_verification_required = NULL,
      M_verification_used = NULL,
      requested_reason = "migration_fixed_candidate_recomputation",
      selected_reason = "fresh_fixed_candidate_recomputation",
      verification_required_reason = "not_decision_verification",
      verification_used_reason = "not_decision_verification"
    )
  }
  termination_code <- if (identical(mode, "a1_proxy")) {
    "deterministic"
  } else {
    "legacy_migration_no_selection"
  }
  .dpprior_new_computation(
    request = setting,
    used = setting,
    orders = orders,
    scaling = .dpprior_new_scaling(
      formula = "legacy_evidence_unavailable"
    ),
    attempts = list(),
    selected_attempt_id = NULL,
    fallback = .dpprior_new_fallback(
      message = "v1.1 fallback lineage was not reconstructed"
    ),
    termination = .dpprior_new_termination(
      code = termination_code,
      message = paste(
        "The v1.1 status and optimizer exit were not reused as candidate",
        "selection evidence."
      ),
      source = "constructor"
    ),
    trace = NULL,
    resources = list()
  )
}


.dpprior_upgrade_unverified <- function(reason) {
  .dpprior_new_verification(
    method = "legacy_evidence_quarantine",
    performed = FALSE,
    passed = FALSE,
    reason = reason,
    settings = list(),
    selected_snapshot = NULL,
    verifier_snapshot = NULL,
    stability = NULL,
    components = list(),
    invariants = list()
  )
}


.dpprior_upgrade_backend <- function() {
  package_version <- tryCatch(
    as.character(utils::packageVersion("DPprior")),
    error = function(e) "development"
  )
  list(
    package = "DPprior",
    package_version = package_version,
    implementation = "schema_upgrade_v1",
    source_commit = NULL
  )
}


.dpprior_upgrade_provenance <- function(
    method,
    source_schema,
    missing_evidence,
    approximation_opt_in = FALSE,
    legacy_active = FALSE,
    legacy_contract = NULL,
    approximation_active = TRUE) {
  .dpprior_new_provenance(
    requested_method = method,
    selected_method = method,
    is_fallback = FALSE,
    approximation = list(
      active = approximation_active,
      opt_in = approximation_opt_in,
      kind = if (approximation_active) "legacy_schema_migration" else NULL,
      warning_code = if (approximation_active) {
        "legacy_object_upgraded"
      } else {
        NULL
      }
    ),
    projection = list(
      applied = FALSE, opt_in = FALSE, policy = NULL, record = NULL
    ),
    parameterization = "Gamma(shape=a, rate=b)",
    backend = .dpprior_upgrade_backend(),
    input_fit = NULL,
    migration = list(
      source_schema = source_schema,
      adapter = "upgrade_DPprior_object",
      lossless = FALSE,
      missing_evidence = unique(missing_evidence),
      warnings = "legacy_object_upgraded"
    ),
    legacy = list(
      active = legacy_active,
      contract = legacy_contract,
      deprecation_stage = if (legacy_active) "v2_migration" else NULL
    )
  )
}


.dpprior_upgrade_compatibility <- function(source_summary,
                                           audit = NULL) {
  views <- list(source = source_summary)
  if (!is.null(audit)) {
    views$fixed_candidate_recomputation <- audit
  }
  .dpprior_new_compatibility(
    top_level_aliases = character(),
    views = views,
    deprecations = list(
      legacy_schema = list(
        code = "legacy_object_upgraded",
        first_deprecated_version = "2.0.0",
        removal_floor = "not_scheduled"
      )
    )
  )
}


.dpprior_upgrade_target <- function(raw, mode, source_schema) {
  target <- raw[["target", exact = TRUE]]
  J <- .dpprior_upgrade_validate_J(
    raw[["J", exact = TRUE]], source_schema
  )
  method <- raw[["method", exact = TRUE]]

  if (mode == "a2_kl") {
    pmf <- target[["pmf", exact = TRUE]]
    ok <- is.numeric(pmf) && !is.object(pmf) && is.null(dim(pmf)) &&
      length(pmf) == J && !anyNA(pmf) && all(is.finite(pmf)) &&
      all(pmf >= 0) && abs(sum(pmf) - 1) <= .TOL_PMF_SUM
    if (!ok) {
      .dpprior_abort_serialization(
        "Legacy A2-KL target PMF failed the strict current PMF contract",
        "dpprior_legacy_object_error", "legacy_invalid_target_pmf",
        detected_schema = source_schema,
        upgrade_action = "Refit from the original target specification."
      )
    }
    kind <- "pmf"
    family <- NULL
    implied <- list(
      mean = target[["mu_K_discrete", exact = TRUE]],
      variance = target[["var_K_discrete", exact = TRUE]]
    )
    request <- list(J = J, pmf = pmf)
    normalized <- list(J = J, pmf = pmf, interval = NULL)
    used <- normalized
    normalization_rule <- "validate_strict_pmf"
  } else {
    mu_K <- .dpprior_upgrade_validate_positive(
      target[["mu_K", exact = TRUE]], "target$mu_K", source_schema
    )
    var_K <- .dpprior_upgrade_validate_positive(
      target[["var_K", exact = TRUE]], "target$var_K", source_schema
    )
    if (mu_K < 1 || mu_K > J) {
      .dpprior_abort_serialization(
        "Legacy target mean lies outside the exact K_J support",
        "dpprior_legacy_object_error", "legacy_target_outside_support",
        detected_schema = source_schema,
        upgrade_action = "Refit from a target satisfying 1 <= mu_K <= J."
      )
    }
    variance_upper <- (mu_K - 1) * (J - mu_K)
    if (var_K > variance_upper + 1e-12) {
      .dpprior_abort_serialization(
        "Legacy target variance exceeds the exact finite-support bound",
        "dpprior_legacy_object_error", "legacy_target_outside_support",
        detected_schema = source_schema,
        upgrade_action = "Refit from a feasible K_J target."
      )
    }
    if (mode == "a1_proxy" && var_K <= mu_K - 1) {
      .dpprior_abort_serialization(
        paste(
          "Legacy A1 target requires a projection that was not recorded",
          "with the explicit v2 opt-in contract."
        ),
        "dpprior_legacy_object_error", "legacy_projection_not_explicit",
        detected_schema = source_schema,
        upgrade_action = paste(
          "Refit with projection='nearest' only after explicitly accepting",
          "the target change."
        )
      )
    }
    used_var <- target[["var_K_used", exact = TRUE]]
    if (is.null(used_var)) {
      used_var <- raw[["var_K_used", exact = TRUE]]
    }
    if (is.null(used_var)) {
      used_var <- var_K
    }
    used_var <- .dpprior_upgrade_validate_positive(
      used_var, "target$var_K_used", source_schema
    )
    if (!identical(as.numeric(var_K), as.numeric(used_var))) {
      .dpprior_abort_serialization(
        paste(
          "Legacy target records an automatic variance projection without",
          "the explicit v2 opt-in evidence required for migration."
        ),
        "dpprior_legacy_object_error", "legacy_projection_not_explicit",
        detected_schema = source_schema,
        upgrade_action = paste(
          "Refit and explicitly select the current projection policy, or",
          "supply a feasible unprojected target."
        )
      )
    }
    request <- list(J = J, mu_K = mu_K, var_K = var_K)
    normalized <- list(
      J = J, mu_K = mu_K, var_K = var_K,
      interval = NULL, pmf = NULL
    )
    used <- normalized
    kind <- "moments"
    pmf <- NULL
    family <- NULL
    implied <- list(mean = mu_K, variance = var_K)
    normalization_rule <- "canonicalize_direct_moments"
  }

  target_method <- paste0("legacy_target:", method)
  .dpprior_new_target_K(
    kind = kind,
    J = J,
    request = request,
    normalized = normalized,
    used = used,
    derivation = list(
      request_to_normalized = list(
        rule = normalization_rule,
        outcome = "canonicalized",
        opt_in = FALSE,
        before = request,
        after = normalized,
        evidence = list(
          source = "recognized_v1.1_shape",
          source_schema = source_schema,
          source_fields = names(target)
        )
      ),
      normalized_to_used = NULL
    ),
    interval = NULL,
    family = family,
    assumptions = list(source = "recognized_v1.1_shape"),
    pmf = pmf,
    implied = implied,
    achieved_interval = NULL,
    residuals = list(),
    tolerances = list(),
    status = "approximate",
    usable = FALSE,
    verified = FALSE,
    message = "Legacy target structure was preserved but lacks v2 verification lineage.",
    parameters = NULL,
    computation = .dpprior_upgrade_computation(target_method),
    verification = .dpprior_upgrade_unverified(
      "v1.1 did not preserve canonical target verification evidence"
    ),
    provenance = .dpprior_upgrade_provenance(
      target_method, source_schema,
      c("target_verification_lineage", "explicit_request_nulls")
    ),
    compatibility = .dpprior_new_compatibility(
      views = list(
        legacy_target = target,
        legacy_target_fields = names(target)
      ),
      deprecations = list()
    )
  )
}


.dpprior_upgrade_fixed_candidate_audit <- function(raw, mode, M_verify) {
  orders <- .dpprior_upgrade_recompute_orders(M_verify)
  a <- raw[["a", exact = TRUE]]
  b <- raw[["b", exact = TRUE]]
  J <- as.integer(raw[["J", exact = TRUE]])

  if (mode == "a2_kl") {
    target <- raw[["target", exact = TRUE]]
    target_pmf <- unname(target[["pmf", exact = TRUE]])
    logS <- compute_log_stirling(J)
    induced <- .a2_kl_induced_log_pmf(
      J, a, b, logS,
      M = orders$M_selected,
      M_verify = orders$M_verify,
      abs_tol = 1e-10,
      rel_tol = 1e-8
    )
    selected <- .a2_kl_pmf_metrics(target_pmf, induced$selected)
    verifier <- .a2_kl_pmf_metrics(target_pmf, induced$verification)
    return(list(
      method = "fresh_fixed_candidate_pmf_audit",
      candidate_unchanged = TRUE,
      candidate = list(
        a = a, b = b, parameterization = "Gamma(shape=a, rate=b)"
      ),
      M_selected = orders$M_selected,
      M_verification = orders$M_verify,
      selected = list(
        mean = selected$mean,
        variance = selected$variance,
        kl = selected$kl,
        l1 = selected$l1
      ),
      verifier = list(
        mean = verifier$mean,
        variance = verifier$variance,
        kl = verifier$kl,
        l1 = verifier$l1
      ),
      numerical_status = induced$metadata$status,
      numerical_reason = induced$metadata$reason,
      target_residuals = list(
        selected_mean = selected$mean_residual,
        selected_variance = selected$variance_residual,
        verifier_mean = verifier$mean_residual,
        verifier_variance = verifier$variance_residual
      ),
      decision_ready = FALSE,
      decision_ready_reason = paste(
        "A fresh fixed-candidate PMF audit cannot reconstruct v1.1",
        "optimizer selection or adequacy lineage; refit required."
      )
    ))
  }

  moments <- exact_K_moments(
    J, a, b,
    M = orders$M_selected,
    M_verify = orders$M_verify,
    abs_tol = 1e-10,
    rel_tol = 1e-8,
    strict = FALSE
  )
  target <- raw[["target", exact = TRUE]]
  target_mean <- target[["mu_K", exact = TRUE]]
  target_variance <- target[["var_K_used", exact = TRUE]]
  if (is.null(target_variance)) {
    target_variance <- raw[["var_K_used", exact = TRUE]]
  }
  if (is.null(target_variance)) {
    target_variance <- target[["var_K", exact = TRUE]]
  }
  list(
    method = "fresh_fixed_candidate_moment_audit",
    candidate_unchanged = TRUE,
    candidate = list(
      a = a, b = b, parameterization = "Gamma(shape=a, rate=b)"
    ),
    M_selected = orders$M_selected,
    M_verification = orders$M_verify,
    selected = list(mean = moments$mean, variance = moments$var),
    order_stability = list(
      status = moments$status,
      reason = moments$quadrature$reason,
      mean_difference = moments$quadrature$mean_difference,
      variance_difference = moments$quadrature$variance_difference,
      mean_tolerance = moments$quadrature$mean_tolerance,
      variance_tolerance = moments$quadrature$variance_tolerance
    ),
    target_residuals = list(
      mean = moments$mean - target_mean,
      variance = moments$var - target_variance
    ),
    decision_ready = FALSE,
    decision_ready_reason = paste(
      "A fresh fixed-candidate moment audit cannot reconstruct v1.1",
      "optimizer selection lineage; refit required."
    )
  )
}


.dpprior_upgrade_stored_achieved <- function(raw, mode, source_schema) {
  fit <- raw[["fit", exact = TRUE]]
  if (is.null(fit)) {
    return(list(achieved = list(), residuals = list()))
  }
  if (typeof(fit) != "list" || !is.list(fit) || is.object(fit) ||
      !is.null(dim(fit)) ||
      anyDuplicated(names(fit))) {
    .dpprior_abort_serialization(
      "Legacy fit evidence is not an ordinary named list",
      "dpprior_legacy_object_error", "legacy_invalid_achieved",
      detected_schema = source_schema,
      upgrade_action = "Refit from the original target specification."
    )
  }
  required <- c("mu_K", "var_K")
  if (!all(required %in% names(fit))) {
    .dpprior_abort_serialization(
      "Legacy fit evidence lacks stored K_J moments",
      "dpprior_legacy_object_error", "legacy_missing_achieved",
      detected_schema = source_schema,
      upgrade_action = "Refit from the original target specification."
    )
  }
  mean <- fit[["mu_K", exact = TRUE]]
  variance <- fit[["var_K", exact = TRUE]]
  valid <- is.numeric(mean) && !is.object(mean) && is.null(dim(mean)) &&
    length(mean) == 1L && !is.na(mean) && is.finite(mean) &&
    is.numeric(variance) && !is.object(variance) &&
    is.null(dim(variance)) && length(variance) == 1L &&
    !is.na(variance) && is.finite(variance) && variance >= 0
  if (!valid) {
    .dpprior_abort_serialization(
      "Legacy stored K_J moments are not finite scalar values",
      "dpprior_legacy_object_error", "legacy_invalid_achieved",
      detected_schema = source_schema,
      upgrade_action = "Refit from the original target specification."
    )
  }

  target <- raw[["target", exact = TRUE]]
  target_mean <- target[["mu_K", exact = TRUE]]
  target_variance <- target[["var_K_used", exact = TRUE]]
  if (is.null(target_variance)) {
    target_variance <- raw[["var_K_used", exact = TRUE]]
  }
  if (is.null(target_variance)) {
    target_variance <- target[["var_K", exact = TRUE]]
  }
  achieved <- list(K = list(
    mean = as.numeric(mean),
    variance = as.numeric(variance),
    estimand = "K_J",
    source = "legacy_serialized_selected_values",
    M = NULL
  ))
  if (mode == "a2_kl") {
    kl <- fit[["kl", exact = TRUE]]
    if (!is.numeric(kl) || is.object(kl) || !is.null(dim(kl)) ||
        length(kl) != 1L || is.na(kl) || !is.finite(kl) || kl < 0) {
      .dpprior_abort_serialization(
        "Legacy A2-KL fit lacks a finite non-negative KL value",
        "dpprior_legacy_object_error", "legacy_invalid_kl",
        detected_schema = source_schema,
        upgrade_action = "Refit from the original target PMF."
      )
    }
    achieved$distribution_fit <- list(
      kl = as.numeric(kl),
      l1 = NULL,
      source = "legacy_serialized_selected_values"
    )
  }
  if (mode == "dual_legacy") {
    dual <- raw[["dual_anchor", exact = TRUE]]
    weight_achieved <- dual[["w1_achieved", exact = TRUE]]
    if (!.dpprior_upgrade_exact_names(
      weight_achieved, c("mean", "prob_gt_50", "prob_gt_90")
    )) {
      .dpprior_abort_serialization(
        paste(
          "Legacy Dual-Anchor achieved weight evidence has an",
          "unsupported shape."
        ),
        "dpprior_legacy_object_error", "legacy_missing_weight_achieved",
        detected_schema = source_schema,
        upgrade_action = "Refit with the current explicit weight target API."
      )
    }
    threshold <- dual[["w1_target", exact = TRUE]][[
      "prob", exact = TRUE
    ]][["threshold", exact = TRUE]]
    weight_field <- if (identical(as.numeric(threshold), 0.5)) {
      "prob_gt_50"
    } else if (identical(as.numeric(threshold), 0.9)) {
      "prob_gt_90"
    } else {
      NULL
    }
    weight_value <- if (is.null(weight_field)) {
      NULL
    } else {
      weight_achieved[[weight_field, exact = TRUE]]
    }
    valid_weight <- is.numeric(weight_value) && !is.object(weight_value) &&
      is.null(dim(weight_value)) && length(weight_value) == 1L &&
      !is.na(weight_value) && is.finite(weight_value) &&
      weight_value >= 0 && weight_value <= 1
    if (!valid_weight) {
      .dpprior_abort_serialization(
        paste(
          "Legacy Dual-Anchor fit does not retain a supported achieved",
          "W_SB tail probability at its target threshold."
        ),
        "dpprior_legacy_object_error", "legacy_missing_weight_achieved",
        detected_schema = source_schema,
        upgrade_action = "Refit with the current explicit weight target API."
      )
    }
    achieved$weight <- list(
      metric = "wsb_tail",
      value = as.numeric(weight_value),
      source = "legacy_serialized_selected_values"
    )
  }
  list(
    achieved = achieved,
    residuals = list(K = list(
      mean = as.numeric(mean - target_mean),
      variance = as.numeric(variance - target_variance),
      source = "recomputed_from_legacy_stored_values"
    ))
  )
}


.dpprior_upgrade_fit_verification <- function(parameters,
                                               achieved,
                                               residuals,
                                               tolerances,
                                               M_selected = NULL) {
  snapshot_source <- achieved[["K", exact = TRUE]][["source", exact = TRUE]]
  selected_snapshot <- .dpprior_new_snapshot(
    parameters = parameters,
    M = M_selected,
    achieved = achieved,
    residuals = residuals,
    tolerances = tolerances,
    finite = TRUE,
    source = snapshot_source
  )
  .dpprior_new_verification(
    method = "legacy_evidence_quarantine",
    performed = FALSE,
    passed = FALSE,
    reason = paste(
      "Stored selected values were quarantined; v1.1 did not retain an",
      "independent verifier or candidate-selection lineage."
    ),
    settings = list(),
    selected_snapshot = selected_snapshot,
    verifier_snapshot = NULL,
    stability = NULL,
    components = list(),
    invariants = list()
  )
}


.dpprior_upgrade_fit_extension <- function(raw, mode, usable) {
  if (mode == "a1_proxy") {
    return(list(proxy = list(
      mapping = list(
        source = "v1.1_serialized_mapping",
        cJ = raw[["cJ", exact = TRUE]],
        scaling = raw[["scaling", exact = TRUE]]
      ),
      mapping_verification = list(
        performed = FALSE,
        passed = FALSE,
        reason = "v1.1 mapping verification lineage unavailable"
      ),
      projection = list(
        applied = FALSE,
        opt_in = FALSE,
        policy = NULL,
        record = NULL
      ),
      caveats = "A1 is a proxy and is not exact prior verification."
    )))
  }
  if (mode == "dual_legacy") {
    dual <- raw[["dual_anchor", exact = TRUE]]
    return(list(legacy = list(
      contract = "v1.1_path_scaled_soft_equality_loss",
      lambda = dual[["lambda", exact = TRUE]],
      losses = list(
        K_loss = dual[["K_loss", exact = TRUE]],
        loss_type = dual[["loss_type", exact = TRUE]]
      ),
      approximation_opt_in = usable,
      warning_code = "legacy_object_upgraded"
    )))
  }
  list()
}


.dpprior_upgrade_weight_target <- function(raw, source_schema) {
  dual <- raw[["dual_anchor", exact = TRUE]]
  weight_target <- dual[["w1_target", exact = TRUE]]
  if (!.dpprior_upgrade_exact_names(weight_target, "prob") ||
      !.dpprior_upgrade_exact_names(
        weight_target[["prob", exact = TRUE]], c("threshold", "value")
      )) {
    .dpprior_abort_serialization(
      "Legacy Dual-Anchor weight target has an unsupported shape",
      "dpprior_legacy_object_error", "legacy_weight_target_unrecognized",
      detected_schema = source_schema,
      upgrade_action = "Refit with DPprior_dual() for legacy reproduction."
    )
  }
  probability_target <- weight_target[["prob", exact = TRUE]]
  threshold <- probability_target[["threshold", exact = TRUE]]
  value <- probability_target[["value", exact = TRUE]]
  valid <- is.numeric(threshold) && !is.object(threshold) &&
    is.null(dim(threshold)) && length(threshold) == 1L &&
    !is.na(threshold) && is.finite(threshold) &&
    threshold > 0 && threshold < 1 &&
    is.numeric(value) && !is.object(value) && is.null(dim(value)) &&
    length(value) == 1L && !is.na(value) && is.finite(value) &&
    value >= 0 && value <= 1
  if (!valid) {
    .dpprior_abort_serialization(
      "Legacy Dual-Anchor probability target is outside [0,1]",
      "dpprior_legacy_object_error", "legacy_invalid_weight_target",
      detected_schema = source_schema,
      upgrade_action = "Refit with a valid explicit W_SB target."
    )
  }
  record <- list(
    metric = "wsb_tail",
    relation = "target",
    value = as.numeric(value),
    threshold = as.numeric(threshold),
    probability = NULL
  )
  .dpprior_new_weight_target(
    request = record,
    normalized = record,
    used = record,
    metric = "wsb_tail",
    relation = "target",
    operator = "target",
    value = as.numeric(value),
    threshold = as.numeric(threshold),
    probability = NULL,
    estimand = "P(W_SB > threshold)",
    units = "probability",
    certification = list(),
    provenance = list(
      source_schema = source_schema,
      transformation = NULL,
      selection = NULL,
      legacy_request = weight_target,
      legacy_semantics = "soft_equality_target_not_a_hard_certificate"
    )
  )
}


.dpprior_upgrade_v1_fit <- function(x, verify, M_verify, allow_legacy) {
  source_schema <- .DPPRIOR_LEGACY_FIT_V1_1
  raw <- unclass(x)
  a <- .dpprior_upgrade_validate_positive(
    raw[["a", exact = TRUE]], "a", source_schema
  )
  b <- .dpprior_upgrade_validate_positive(
    raw[["b", exact = TRUE]], "b", source_schema
  )
  J <- .dpprior_upgrade_validate_J(
    raw[["J", exact = TRUE]], source_schema
  )
  method <- raw[["method", exact = TRUE]]
  mode <- .dpprior_upgrade_fit_mode(method)
  if (is.null(mode)) {
    .dpprior_abort_serialization(
      sprintf("Legacy fit method '%s' has no safe migration adapter", method),
      "dpprior_legacy_object_error", "legacy_method_unrecognized", x,
      source_schema, .DPPRIOR_RESULT_SCHEMA_V1,
      "Refit with a current supported method."
    )
  }

  source_summary <- .dpprior_upgrade_source_summary(x, source_schema)
  target_K <- .dpprior_upgrade_target(raw, mode, source_schema)
  target <- if (mode == "dual_legacy") {
    list(
      K = target_K,
      weight = .dpprior_upgrade_weight_target(raw, source_schema)
    )
  } else {
    list(K = target_K)
  }
  stored <- .dpprior_upgrade_stored_achieved(raw, mode, source_schema)
  audit_warnings <- list()
  audit <- if (verify) {
    value <- tryCatch(
      withCallingHandlers(
        .dpprior_upgrade_fixed_candidate_audit(raw, mode, M_verify),
        warning = function(warning) {
          audit_warnings[[length(audit_warnings) + 1L]] <<- warning
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        .dpprior_abort_serialization(
          paste(
            "Fresh fixed-candidate recomputation failed during legacy",
            "migration."
          ),
          "dpprior_legacy_object_error", "legacy_recomputation_failed", x,
          source_schema, .DPPRIOR_RESULT_SCHEMA_V1,
          "Refit from the original target specification.", cause = e
        )
      }
    )
    if (length(audit_warnings)) {
      .dpprior_abort_serialization(
        "Fresh fixed-candidate recomputation emitted an unexpected warning",
        "dpprior_legacy_object_error", "legacy_recomputation_warning", x,
        source_schema, .DPPRIOR_RESULT_SCHEMA_V1,
        "Refit and inspect current numerical diagnostics directly.",
        cause = audit_warnings[[1L]]
      )
    }
    value
  } else {
    NULL
  }

  selected_M <- NULL
  if (identical(mode, "a1_proxy") &&
      !"K" %in% names(stored$achieved) && !is.null(audit)) {
    selected_M <- as.integer(audit$M_selected)
    stored$achieved <- list(K = list(
      mean = audit$selected$mean,
      variance = audit$selected$variance,
      estimand = "K_J",
      source = "fresh_migration_fixed_candidate_recomputation",
      M = selected_M
    ))
    stored$residuals <- list(K = list(
      mean = audit$target_residuals$mean,
      variance = audit$target_residuals$variance,
      source = "fresh_migration_fixed_candidate_recomputation"
    ))
  }
  candidate_complete <- "K" %in% names(stored$achieved) &&
    (!identical(mode, "dual_legacy") ||
       "weight" %in% names(stored$achieved))
  quarantine_candidate <- mode %in% c("a2_moment", "a2_kl") ||
    !candidate_complete
  approximate_usable <- allow_legacy && !quarantine_candidate &&
    mode %in% c("a1_proxy", "dual_legacy")
  missing_evidence <- .dpprior_upgrade_missing_evidence(mode)
  if (quarantine_candidate) {
    missing_evidence <- c(
      missing_evidence, "public_candidate_quarantined_pending_refit"
    )
  }
  stored_parameters <- .dpprior_new_parameters(
    a, b, "Gamma(shape=a, rate=b)"
  )
  truth <- .dpprior_upgrade_fit_truth(mode)
  public_parameters <- if (quarantine_candidate) NULL else stored_parameters
  public_achieved <- if (quarantine_candidate) list() else stored$achieved
  public_residuals <- if (quarantine_candidate) list() else stored$residuals
  verification <- if (quarantine_candidate) {
    .dpprior_upgrade_unverified(paste(
      "The v1.1 A2 parameter pair is retained only in compatibility",
      "evidence because optimizer candidate-selection lineage is absent."
    ))
  } else {
    .dpprior_upgrade_fit_verification(
      stored_parameters, public_achieved, public_residuals,
      truth$tolerances, M_selected = selected_M
    )
  }
  source_summary$public_candidate_quarantined <- quarantine_candidate
  source_summary$required_action <- "refit_with_current_API"
  provenance <- .dpprior_upgrade_provenance(
    method = method,
    source_schema = source_schema,
    missing_evidence = missing_evidence,
    approximation_opt_in = approximate_usable,
    legacy_active = identical(mode, "dual_legacy"),
    legacy_contract = if (identical(mode, "dual_legacy")) {
      "v1.1_path_scaled_soft_equality_loss"
    } else {
      NULL
    }
  )
  message <- paste(
    "Migrated v1.1 candidate is approximate and unverified because optimizer",
    "selection and independent verification lineage were not serialized;",
    "refit is required for decision readiness."
  )
  result <- .dpprior_new_fit(
    mode = mode,
    method = method,
    J = J,
    status = "approximate",
    usable = approximate_usable,
    verified = FALSE,
    message = message,
    parameters = public_parameters,
    target = target,
    achieved = public_achieved,
    residuals = public_residuals,
    tolerances = truth$tolerances,
    computation = .dpprior_upgrade_computation(
      method, truth$controls, mode = mode, M_selected = selected_M
    ),
    verification = verification,
    provenance = provenance,
    compatibility = .dpprior_upgrade_compatibility(source_summary, audit),
    extension = .dpprior_upgrade_fit_extension(raw, mode, approximate_usable)
  )
  .dpprior_validate_object(result)
  .dpprior_warn_legacy_upgrade(source_schema, missing_evidence)
  result
}


.dpprior_upgrade_diagnostics_truth <- function(parameters, J, orders) {
  absolute <- 1e-10
  relative <- 1e-8
  pmf_mass <- .TOL_PMF_SUM
  a <- parameters[["a", exact = TRUE]]
  b <- parameters[["b", exact = TRUE]]
  K <- .get_K_pmf_support(
    J, a, b,
    M = orders$M_selected,
    M_verify = orders$M_verify,
    abs_tol = absolute,
    rel_tol = relative
  )
  selected_pmf <- unname(as.numeric(K[["pmf", exact = TRUE]]))
  verifier_pmf <- unname(as.numeric(
    K[["verification_pmf", exact = TRUE]]
  ))
  selected_K <- .dpprior_target_pmf_moments(selected_pmf)
  verifier_K <- .dpprior_target_pmf_moments(verifier_pmf)
  selected_weight <- as.numeric(mean_w1(a, b, orders$M_selected))
  verifier_weight <- as.numeric(mean_w1(a, b, orders$M_verify))
  selected_rho <- c(
    mean = as.numeric(mean_rho(a, b, orders$M_selected)),
    variance = as.numeric(var_rho(a, b, orders$M_selected))
  )
  verifier_rho <- c(
    mean = as.numeric(mean_rho(a, b, orders$M_verify)),
    variance = as.numeric(var_rho(a, b, orders$M_verify))
  )
  delta <- c(
    K.mean = abs(selected_K[["mean"]] - verifier_K[["mean"]]),
    K.variance = abs(
      selected_K[["variance"]] - verifier_K[["variance"]]
    ),
    K.pmf_l1 = sum(abs(selected_pmf - verifier_pmf)),
    weights.mean = abs(selected_weight - verifier_weight),
    coclustering.mean = abs(selected_rho[["mean"]] - verifier_rho[["mean"]]),
    coclustering.variance = abs(
      selected_rho[["variance"]] - verifier_rho[["variance"]]
    )
  )
  scalar_tolerance <- function(selected, verifier) {
    absolute + relative * max(abs(selected), abs(verifier), 1)
  }
  refinement <- c(
    K.mean = scalar_tolerance(
      selected_K[["mean"]], verifier_K[["mean"]]
    ),
    K.variance = scalar_tolerance(
      selected_K[["variance"]], verifier_K[["variance"]]
    ),
    K.pmf_l1 = absolute + relative,
    weights.mean = scalar_tolerance(selected_weight, verifier_weight),
    coclustering.mean = scalar_tolerance(
      selected_rho[["mean"]], verifier_rho[["mean"]]
    ),
    coclustering.variance = scalar_tolerance(
      selected_rho[["variance"]], verifier_rho[["variance"]]
    )
  )
  component_pass <- c(
    alpha = TRUE,
    K = all(delta[c("K.mean", "K.variance", "K.pmf_l1")] <=
              refinement[c("K.mean", "K.variance", "K.pmf_l1")]) &&
      abs(sum(selected_pmf) - 1) <= pmf_mass &&
      abs(sum(verifier_pmf) - 1) <= pmf_mass,
    weights = delta[["weights.mean"]] <= refinement[["weights.mean"]],
    coclustering = all(delta[c(
      "coclustering.mean", "coclustering.variance"
    )] <= refinement[c(
      "coclustering.mean", "coclustering.variance"
    )])
  )
  list(
    controls = list(
      absolute_tolerance = absolute,
      relative_tolerance = relative,
      pmf_mass_tolerance = pmf_mass
    ),
    tolerances = list(diagnostics = list(
      absolute = absolute,
      relative = relative,
      pmf_mass = pmf_mass,
      refinement = refinement
    )),
    residuals = list(diagnostics = delta),
    selected = list(
      alpha = list(mean = a / b, CV = 1 / sqrt(a)),
      K = list(
        mean = unname(selected_K[["mean"]]),
        variance = unname(selected_K[["variance"]]),
        estimand = "K_J",
        source = "fresh_diagnostics_selected_order",
        M = orders$M_selected,
        pmf = selected_pmf
      ),
      weights = list(mean = selected_weight),
      coclustering = as.list(selected_rho)
    ),
    verifier = list(
      alpha = list(mean = a / b, CV = 1 / sqrt(a)),
      K = list(
        mean = unname(verifier_K[["mean"]]),
        variance = unname(verifier_K[["variance"]]),
        estimand = "K_J",
        source = "fresh_diagnostics_verifier_evidence",
        M = orders$M_verify,
        pmf = verifier_pmf
      ),
      weights = list(mean = verifier_weight),
      coclustering = as.list(verifier_rho)
    ),
    component_pass = component_pass
  )
}


.dpprior_upgrade_diagnostic_attempts <- function(component_pass) {
  lapply(seq_along(.DPPRIOR_DIAGNOSTIC_COMPONENTS), function(index) {
    name <- .DPPRIOR_DIAGNOSTIC_COMPONENTS[[index]]
    status <- if (component_pass[[name]]) "converged" else "approximate"
    unavailable <- c(
      start = "component diagnostic has no optimizer start",
      bounds = "component diagnostic has no optimizer bounds",
      candidate_parameters = "diagnostics do not select fit parameters",
      candidate_objective = "diagnostics do not optimize an objective"
    )
    .dpprior_new_attempt(
      id = paste0("diagnostic-", name),
      stage = "diagnostic_component",
      method = unname(.DPPRIOR_DIAGNOSTIC_ATTEMPT_METHODS[[name]]),
      start = NULL,
      bounds = NULL,
      control = list(component = name),
      exit_code = 0L,
      message = "Component diagnostic freshly reconstructed.",
      iterations = 0L,
      evaluations = list(function_count = 1L),
      candidate_parameters = NULL,
      candidate_objective = NULL,
      elapsed_seconds = 0,
      warnings = character(),
      error = NULL,
      selected = FALSE,
      reason_code = paste0("component_", status),
      unavailable = unavailable
    )
  })
}


.dpprior_upgrade_diagnostics_computation <- function(truth, orders) {
  method <- "canonical_prior_diagnostics"
  setting <- list(
    method = method,
    controls = truth$controls,
    parameterization = "Gamma(shape=a, rate=b)"
  )
  attempts <- .dpprior_upgrade_diagnostic_attempts(truth$component_pass)
  .dpprior_new_computation(
    request = setting,
    used = setting,
    orders = .dpprior_new_orders(
      M_requested = orders$M_selected,
      M_selected = orders$M_selected,
      M_verification_required = orders$M_verify,
      M_verification_used = orders$M_verify,
      requested_reason = "migration_recomputation_default",
      selected_reason = "fresh_diagnostics_selected_order",
      verification_required_reason = "current_quadrature_contract",
      verification_used_reason = "fresh_diagnostics_verifier_order"
    ),
    scaling = .dpprior_new_scaling(),
    attempts = attempts,
    selected_attempt_id = NULL,
    fallback = .dpprior_new_fallback(),
    termination = .dpprior_new_termination(
      code = "diagnostics_recomputed",
      message = "All diagnostic components were freshly reconstructed.",
      source = "component_aggregation"
    ),
    trace = NULL,
    resources = list(
      component_elapsed_seconds = stats::setNames(
        rep(0, length(.DPPRIOR_DIAGNOSTIC_COMPONENTS)),
        .DPPRIOR_DIAGNOSTIC_COMPONENTS
      )
    )
  )
}


.dpprior_upgrade_diagnostics_verification <- function(truth,
                                                       parameters,
                                                       orders) {
  selected <- .dpprior_new_snapshot(
    parameters = parameters,
    M = orders$M_selected,
    achieved = truth$selected,
    residuals = truth$residuals,
    tolerances = truth$tolerances,
    finite = TRUE,
    source = "fresh_diagnostics_selected_order"
  )
  verifier <- .dpprior_new_snapshot(
    parameters = parameters,
    M = orders$M_verify,
    achieved = truth$verifier,
    residuals = truth$residuals,
    tolerances = truth$tolerances,
    finite = TRUE,
    source = "fresh_diagnostics_verifier_evidence"
  )
  component_reference <- stats::setNames(
    rep(TRUE, length(truth$component_pass)), names(truth$component_pass)
  )
  check <- function(value, reference = TRUE) {
    .dpprior_new_check(
      value = value,
      reference = reference,
      tolerance = NULL,
      operator = "identical",
      source = "fresh_component_specific_checks"
    )
  }
  .dpprior_new_verification(
    method = "fresh_component_specific_diagnostics",
    performed = TRUE,
    passed = all(truth$component_pass),
    reason = "Fresh component-specific reconstruction completed.",
    settings = list(
      M_selected = orders$M_selected,
      M_verification = orders$M_verify
    ),
    selected_snapshot = selected,
    verifier_snapshot = verifier,
    stability = NULL,
    components = list(component_aggregation = check(
      truth$component_pass, component_reference
    )),
    invariants = list(
      fixed_parameters = check(TRUE),
      dominance_category_removed = check(TRUE)
    )
  )
}


.dpprior_upgrade_v1_diagnostics <- function(x, verify, M_verify) {
  source_schema <- .DPPRIOR_LEGACY_DIAGNOSTICS_V1_1
  if (!verify) {
    .dpprior_abort_serialization(
      paste(
        "A v1.1 diagnostics object cannot be migrated without fresh",
        "recomputation because its categorical dominance field is prohibited."
      ),
      "dpprior_legacy_object_error", "diagnostics_recomputation_required", x,
      source_schema, .DPPRIOR_RESULT_SCHEMA_V1,
      "Call upgrade_DPprior_object(x, verify = TRUE)."
    )
  }
  raw <- unclass(x)
  a <- .dpprior_upgrade_validate_positive(
    raw[["a", exact = TRUE]], "a", source_schema
  )
  b <- .dpprior_upgrade_validate_positive(
    raw[["b", exact = TRUE]], "b", source_schema
  )
  J <- .dpprior_upgrade_validate_J(
    raw[["J", exact = TRUE]], source_schema
  )
  orders <- .dpprior_upgrade_recompute_orders(M_verify)
  parameters <- .dpprior_new_parameters(a, b, "Gamma(shape=a, rate=b)")
  caught_warnings <- list()
  truth <- tryCatch(
    withCallingHandlers(
      .dpprior_upgrade_diagnostics_truth(parameters, J, orders),
      warning = function(w) {
        caught_warnings[[length(caught_warnings) + 1L]] <<- w
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      .dpprior_abort_serialization(
        "Fresh diagnostics recomputation failed during legacy migration",
        "dpprior_legacy_object_error", "legacy_recomputation_failed", x,
        source_schema, .DPPRIOR_RESULT_SCHEMA_V1,
        "Refit and recompute diagnostics with the current API.", cause = e
      )
    }
  )
  if (length(caught_warnings)) {
    .dpprior_abort_serialization(
      "Fresh diagnostics recomputation emitted an unexpected warning",
      "dpprior_legacy_object_error", "legacy_recomputation_warning", x,
      source_schema, .DPPRIOR_RESULT_SCHEMA_V1,
      "Recompute diagnostics directly and inspect the warning.",
      cause = caught_warnings[[1L]]
    )
  }

  source_summary <- .dpprior_upgrade_source_summary(x, source_schema)
  component_status <- ifelse(
    truth$component_pass, "converged", "approximate"
  )
  precedence <- c(
    converged = 1L, boundary = 2L, approximate = 3L,
    infeasible = 4L, failed = 5L
  )
  status <- component_status[[
    which.max(unname(precedence[component_status]))
  ]]
  usable <- status %in% c("converged", "boundary")
  verified <- usable && all(truth$component_pass)
  approximation_active <- identical(status, "approximate")
  diagnostics <- list(
    policy_results = list(),
    warnings = character(),
    alpha = c(
      list(
        status = component_status[["alpha"]],
        usable = truth$component_pass[["alpha"]],
        verified = truth$component_pass[["alpha"]]
      ),
      truth$selected$alpha
    ),
    K = c(
      list(
        status = component_status[["K"]],
        usable = truth$component_pass[["K"]],
        verified = truth$component_pass[["K"]]
      ),
      truth$selected$K[c("mean", "variance", "pmf", "M")]
    ),
    weights = c(
      list(
        status = component_status[["weights"]],
        usable = truth$component_pass[["weights"]],
        verified = truth$component_pass[["weights"]]
      ),
      truth$selected$weights
    ),
    coclustering = c(
      list(
        status = component_status[["coclustering"]],
        usable = truth$component_pass[["coclustering"]],
        verified = truth$component_pass[["coclustering"]]
      ),
      truth$selected$coclustering
    )
  )
  missing_evidence <- c(
    "original_diagnostic_orders", "original_component_attempts",
    "original_verifier_lineage", "categorical_dominance_semantics_discarded"
  )
  result <- .dpprior_new_diagnostics(
    method = "canonical_prior_diagnostics",
    J = J,
    status = status,
    usable = usable,
    verified = verified,
    message = paste(
      "All canonical diagnostic components were freshly reconstructed at",
      "selected and verifier orders; v1.1 dominance_risk was discarded."
    ),
    parameters = parameters,
    target = list(
      requested_components = .DPPRIOR_DIAGNOSTIC_COMPONENTS,
      warning_policy = NULL
    ),
    achieved = truth$selected,
    residuals = truth$residuals,
    tolerances = truth$tolerances,
    computation = .dpprior_upgrade_diagnostics_computation(truth, orders),
    verification = .dpprior_upgrade_diagnostics_verification(
      truth, parameters, orders
    ),
    provenance = .dpprior_upgrade_provenance(
      method = "canonical_prior_diagnostics",
      source_schema = source_schema,
      missing_evidence = missing_evidence,
      approximation_opt_in = FALSE,
      legacy_active = FALSE,
      approximation_active = approximation_active
    ),
    diagnostics = diagnostics,
    compatibility = .dpprior_upgrade_compatibility(source_summary)
  )
  if (.dpprior_has_recursive_field(result, "dominance_risk")) {
    .dpprior_abort_serialization(
      "Diagnostics migration retained the prohibited dominance_risk field",
      "dpprior_legacy_object_error", "prohibited_legacy_field_retained", x,
      source_schema, .DPPRIOR_RESULT_SCHEMA_V1,
      "Recompute diagnostics with the current API."
    )
  }
  .dpprior_validate_object(result)
  .dpprior_warn_legacy_upgrade(source_schema, missing_evidence)
  result
}


.dpprior_has_recursive_field <- function(x, field) {
  if (!is.list(x)) {
    return(FALSE)
  }
  raw <- if (is.object(x)) unclass(x) else x
  field %in% names(raw) || any(vapply(
    unname(raw), .dpprior_has_recursive_field, logical(1), field = field
  ))
}


.dpprior_upgrade_v1_object <- function(x,
                                       verify = TRUE,
                                       M_verify = NULL,
                                       allow_legacy = FALSE) {
  detected <- .dpprior_detect_schema(x, "auto")
  tryCatch(
    {
      if (identical(detected, .DPPRIOR_LEGACY_FIT_V1_1)) {
        return(.dpprior_upgrade_v1_fit(
          x, verify = verify, M_verify = M_verify,
          allow_legacy = allow_legacy
        ))
      }
      if (identical(detected, .DPPRIOR_LEGACY_DIAGNOSTICS_V1_1)) {
        return(.dpprior_upgrade_v1_diagnostics(
          x, verify = verify, M_verify = M_verify
        ))
      }
      .dpprior_abort_serialization(
        sprintf("No v1.1 migration adapter exists for schema '%s'", detected),
        c("dpprior_legacy_object_error", "dpprior_schema_unsupported_error"),
        "legacy_schema_unrecognized", x, detected,
        .DPPRIOR_RESULT_SCHEMA_V1,
        "Recreate or refit the object with the current DPprior API."
      )
    },
    error = function(error) {
      if (inherits(error, "dpprior_serialization_error") ||
          inherits(error, "dpprior_upgrade_control_error")) {
        stop(error)
      }
      .dpprior_abort_serialization(
        paste(
          "The recognized v1.1 object could not be normalized into the",
          "canonical schema."
        ),
        "dpprior_legacy_object_error", "legacy_normalization_failed", x,
        detected, .DPPRIOR_RESULT_SCHEMA_V1,
        "Refit or recompute the object with the current DPprior API.",
        cause = error
      )
    }
  )
}


#' Upgrade a serialized DPprior object conservatively
#'
#' Converts an exactly recognized DPprior 1.1 fit or diagnostics object into
#' the canonical result schema. A migrated fit never inherits the former
#' `success` or `converged` claim: it remains approximate and unverified because
#' v1.1 did not serialize candidate-selection and independent-verifier lineage.
#'
#' @param x A canonical object, or an exactly recognized DPprior 1.1 fit or
#'   diagnostics object.
#' @param verify Logical. For fits, record a fresh fixed-candidate numerical
#'   audit without changing the candidate or scientific status. Diagnostics
#'   migration requires `TRUE` and recomputes all metrics from `a`, `b`, and
#'   `J`.
#' @param M_verify Optional current verifier order used only for the fresh
#'   migration audit or diagnostics recomputation.
#' @param allow_legacy Logical. When `TRUE`, only migrated A1 and retained
#'   legacy Dual-Anchor fits may be usable approximations. It never makes an
#'   A2, hard, or soft result decision-ready.
#'
#' @return A validated canonical object. Current \code{dpprior.result/1},
#'   \code{dpprior.target/1}, and \code{dpprior.weight-target/1} inputs are
#'   returned after validation. An exactly recognized 1.1 object is converted
#'   conservatively and retains explicit migration and quarantine evidence;
#'   it is never promoted to a verified current calibration merely because a
#'   fixed-candidate audit succeeds. The input is not modified.
#'
#' @details Objects that are not exactly recognized fail closed with a typed
#'   serialization condition. Refit with the current API when decision-ready
#'   calibration is required; compatibility views are not scientific fields.
#'
#' @seealso \code{\link{DPprior_fit}}, \code{\link{DPprior_target_K}}
#' @export
upgrade_DPprior_object <- function(x,
                                   verify = TRUE,
                                   M_verify = NULL,
                                   allow_legacy = FALSE) {
  verify <- .dpprior_upgrade_validate_flag(verify, "verify")
  allow_legacy <- .dpprior_upgrade_validate_flag(
    allow_legacy, "allow_legacy"
  )
  M_verify <- .dpprior_upgrade_validate_M_verify(M_verify, verify)
  detected <- .dpprior_detect_schema(x, "auto")

  if (detected %in% c(
    .DPPRIOR_RESULT_SCHEMA_V1,
    .DPPRIOR_TARGET_SCHEMA_V1,
    .DPPRIOR_WEIGHT_TARGET_SCHEMA_V1
  )) {
    return(.dpprior_validate_object(x))
  }
  if (detected %in% c(
    .DPPRIOR_LEGACY_FIT_V1_1,
    .DPPRIOR_LEGACY_DIAGNOSTICS_V1_1
  )) {
    return(.dpprior_upgrade_v1_object(
      x,
      verify = verify,
      M_verify = M_verify,
      allow_legacy = allow_legacy
    ))
  }

  .dpprior_abort_serialization(
    message = sprintf(
      "Object schema '%s' is not recognized by the conservative upgrader",
      detected
    ),
    subclass = c(
      "dpprior_legacy_object_error", "dpprior_schema_unsupported_error"
    ),
    code = "legacy_schema_unrecognized",
    x = x,
    detected_schema = detected,
    expected_schema = .DPPRIOR_RESULT_SCHEMA_V1,
    upgrade_action = paste(
      "Use an unmodified DPprior 1.1 object with an exact supported shape,",
      "or recreate/refit it with the current API."
    )
  )
}
