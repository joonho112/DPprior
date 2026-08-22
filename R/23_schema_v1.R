# =============================================================================
# Module 23: Canonical schema v1 constructors and validators
# =============================================================================
#
# This module owns dependency-free construction and validation of the v1
# canonical object contract. It never selects candidates, repairs targets, or
# upgrades legacy objects. Validators may invoke frozen deterministic evaluators
# only to recompute and bind truth-making evidence already published in an
# object; that verification cannot create or improve a scientific claim.
# =============================================================================

.DPPRIOR_SCHEMA_VERSION <- 1L
.DPPRIOR_RESULT_SCHEMA_NAME <- "dpprior.result"
.DPPRIOR_TARGET_SCHEMA_NAME <- "dpprior.target"
.DPPRIOR_WEIGHT_TARGET_SCHEMA_NAME <- "dpprior.weight-target"

.DPPRIOR_RESULT_COMMON_FIELDS <- c(
  "schema", "object_type", "mode", "method", "J",
  "status", "usable", "verified", "message",
  "parameters", "target", "achieved", "residuals", "tolerances",
  "computation", "verification", "provenance", "compatibility"
)

.DPPRIOR_TARGET_FIELDS <- c(
  "schema", "kind", "J", "support", "request", "normalized", "used",
  "derivation", "interval", "family", "assumptions", "pmf", "implied",
  "achieved_interval", "residuals", "tolerances", "status", "usable",
  "verified", "message", "parameters", "computation", "verification",
  "provenance", "compatibility"
)

.DPPRIOR_WEIGHT_TARGET_FIELDS <- c(
  "schema", "kind", "request", "normalized", "used", "metric",
  "relation", "operator", "value", "threshold", "probability", "estimand",
  "units", "certification", "provenance"
)

.DPPRIOR_RESULT_MODES <- c(
  "a1_proxy", "a2_moment", "a2_kl", "dual_hard", "dual_soft",
  "dual_legacy", "prior_diagnostics", "elicitation_sensitivity"
)

.DPPRIOR_MODE_OBJECT_TYPE <- c(
  a1_proxy = "fit",
  a2_moment = "fit",
  a2_kl = "fit",
  dual_hard = "fit",
  dual_soft = "fit",
  dual_legacy = "fit",
  prior_diagnostics = "diagnostics",
  elicitation_sensitivity = "sensitivity"
)

.DPPRIOR_MODE_METHODS <- list(
  a1_proxy = "A1",
  a2_moment = c("A2-MN", "A2-MN+NM"),
  a2_kl = "A2-KL",
  dual_hard = "dual_anchor_hard_inequality",
  dual_soft = "dual-soft",
  dual_legacy = "dual-anchor",
  prior_diagnostics = "canonical_prior_diagnostics",
  elicitation_sensitivity = "elicitation_sensitivity"
)

.DPPRIOR_MODE_REQUIRED_EXTENSION <- c(
  a1_proxy = "proxy",
  a2_moment = NA_character_,
  a2_kl = NA_character_,
  dual_hard = "constraint",
  dual_soft = "tradeoff",
  dual_legacy = "legacy",
  prior_diagnostics = "diagnostics",
  elicitation_sensitivity = "sensitivity"
)

.DPPRIOR_RESULT_EXTENSIONS <- c(
  "proxy", "constraint", "tradeoff", "legacy", "diagnostics", "sensitivity"
)
.DPPRIOR_DIAGNOSTIC_COMPONENTS <- c(
  "alpha", "K", "weights", "coclustering"
)
.DPPRIOR_DIAGNOSTIC_ATTEMPT_METHODS <- c(
  alpha = "closed_form_and_stats_qgamma",
  K = "gauss-laguerre-marginal-pmf-and-moments",
  weights = "W_SB-closed-form-plus-W_max-contracted-backend",
  coclustering = "gauss-laguerre-rho-moments"
)
.DPPRIOR_SENSITIVITY_METRICS <- c(
  "a", "b", "E_alpha", "CV_alpha", "E_K_J", "Var_K_J", "CV_K_J",
  "interval_requested", "interval_achieved", "interval_residual",
  "interval_left_tail", "interval_right_tail", "E_W_SB",
  "P_W_SB_gt_50", "P_W_SB_gt_90", "P_W_max_gt_50",
  "P_W_max_gt_90", "P_W_max_gt_50_lower_bound",
  "P_W_max_gt_50_upper_bound", "P_W_max_gt_90_lower_bound",
  "P_W_max_gt_90_upper_bound", "E_rho"
)
.DPPRIOR_SENSITIVITY_METRIC_COMPONENT <- c(
  a = "solver", b = "solver", E_alpha = "alpha", CV_alpha = "alpha",
  E_K_J = "K", Var_K_J = "K", CV_K_J = "K",
  interval_requested = "interval", interval_achieved = "interval",
  interval_residual = "interval", interval_left_tail = "interval",
  interval_right_tail = "interval", E_W_SB = "weights",
  P_W_SB_gt_50 = "weights", P_W_SB_gt_90 = "weights",
  P_W_max_gt_50 = "weights", P_W_max_gt_90 = "weights",
  P_W_max_gt_50_lower_bound = "weights",
  P_W_max_gt_50_upper_bound = "weights",
  P_W_max_gt_90_lower_bound = "weights",
  P_W_max_gt_90_upper_bound = "weights", E_rho = "coclustering"
)
.DPPRIOR_SENSITIVITY_METRIC_SOURCE <- c(
  a = "calibration_fit", b = "calibration_fit",
  E_alpha = "analytic_parameters", CV_alpha = "analytic_parameters",
  E_K_J = "fit_achieved", Var_K_J = "fit_achieved",
  CV_K_J = "fit_achieved", interval_requested = "interval_backcheck",
  interval_achieved = "interval_backcheck",
  interval_residual = "interval_backcheck",
  interval_left_tail = "interval_backcheck",
  interval_right_tail = "interval_backcheck", E_W_SB = "diagnostics",
  P_W_SB_gt_50 = "diagnostics", P_W_SB_gt_90 = "diagnostics",
  P_W_max_gt_50 = "diagnostics", P_W_max_gt_90 = "diagnostics",
  P_W_max_gt_50_lower_bound = "diagnostics",
  P_W_max_gt_50_upper_bound = "diagnostics",
  P_W_max_gt_90_lower_bound = "diagnostics",
  P_W_max_gt_90_upper_bound = "diagnostics", E_rho = "diagnostics"
)
.DPPRIOR_SENSITIVITY_LOCAL_FIELDS <- c(
  "scenario_key", "axis", "axis_value", "settings_key",
  "lower_scenario_key", "upper_scenario_key", "lower_value", "upper_value",
  "metric", "component", "derivative", "method", "reason"
)

.DPPRIOR_ATTEMPT_FIELDS <- c(
  "id", "stage", "method", "start", "bounds", "control", "exit_code",
  "message", "iterations", "evaluations", "candidate_parameters",
  "candidate_objective", "elapsed_seconds", "warnings", "error", "selected",
  "reason_code", "unavailable"
)

.DPPRIOR_ATTEMPT_TERMINAL_FAILURE_REASONS <- c(
  "optimizer_error", "optimizer_exit_nonzero", "nonfinite_candidate",
  "nonfinite_objective", "candidate_evaluation_failed", "no_candidate"
)

.DPPRIOR_ATTEMPT_SCIENTIFIC_REJECTION_REASONS <- c(
  "residual_adequacy_failed", "pmf_adequacy_failed",
  "order_stability_failed", "constraint_verification_failed",
  "objective_recomputation_failed", "local_optimality_failed",
  "candidate_eligibility_failed", "diagnostic_only"
)

.DPPRIOR_COMPUTATION_FIELDS <- c(
  "request", "used", "orders", "scaling", "attempts",
  "candidate_evaluations", "selected_candidate_id", "selected_attempt_id",
  "fallback", "termination", "trace", "resources"
)

.DPPRIOR_CANDIDATE_EVALUATION_FIELDS <- c(
  "id", "attempt_id", "method", "generator", "parameters", "objective_kind",
  "recorded_objective_kind", "selection_objective_kind",
  "recorded_objective", "recorded_objective_available",
  "recorded_objective_reason", "fresh_objective", "selection_objective",
  "objective_tolerance", "objective_passed", "selected_snapshot",
  "verifier_snapshot", "checks", "execution_success", "optimizer_supported",
  "scientifically_eligible", "selection_eligible", "diagnostic_eligible",
  "decision_eligible",
  "selected", "outcome", "rejection_codes", "source"
)

.DPPRIOR_VERIFICATION_FIELDS <- c(
  "method", "performed", "passed", "reason", "settings",
  "selected_snapshot", "verifier_snapshot", "stability", "components",
  "invariants"
)

.DPPRIOR_DECISION_CHECK_FIELDS <- c(
  "passed", "value", "reference", "tolerance", "operator", "source"
)

.DPPRIOR_STABILITY_FIELDS <- c(
  "delta", "tolerance", "formula", "scale_floor", "passed", "source"
)

.DPPRIOR_TERMINATION_CODES <- list(
  converged = c(
    "selected", "converged", "deterministic", "closed_form", "endpoint",
    "diagnostics_recomputed"
  ),
  boundary = c("selected", "boundary", "deterministic", "closed_form",
               "endpoint"),
  approximate = c(
    "selected", "approximate", "deterministic", "closed_form",
    "legacy_migration_no_selection"
  ),
  infeasible = c("infeasible", "certified_infeasible"),
  failed = c("failed", "error", "no_candidate")
)

.DPPRIOR_PROVENANCE_FIELDS <- c(
  "requested_method", "selected_method", "is_fallback", "approximation",
  "projection", "parameterization", "backend", "input_fit", "migration",
  "legacy"
)

.DPPRIOR_COMPATIBILITY_FIELDS <- c(
  "top_level_aliases", "views", "deprecations"
)


# --- Typed schema conditions -------------------------------------------------

.dpprior_abort_schema <- function(message,
                                  code,
                                  path,
                                  expected = NULL,
                                  actual = NULL,
                                  call = NULL) {
  stop(.dpprior_new_condition(
    message = message,
    classes = c("dpprior_schema_error", "dpprior_error", "error"),
    call = call,
    code = as.character(code),
    path = as.character(path),
    expected = expected,
    actual = actual
  ))
}


.dpprior_schema_fail <- function(code, path, expected, actual = NULL) {
  .dpprior_abort_schema(
    sprintf("Invalid canonical schema at `%s`: expected %s.", path, expected),
    code = code,
    path = path,
    expected = expected,
    actual = actual
  )
}


.dpprior_schema_require <- function(ok, code, path, expected, actual = NULL) {
  if (!isTRUE(ok)) {
    .dpprior_schema_fail(code, path, expected, actual)
  }
  invisible(TRUE)
}


.dpprior_schema_collect_validation <- function(validator, x, collect = FALSE) {
  .dpprior_schema_require(
    .dpprior_schema_is_scalar_logical(collect), "type", "collect",
    "one non-missing logical value", collect
  )
  if (!collect) {
    return(validator(x))
  }
  tryCatch(
    {
      validator(x)
      list(valid = TRUE, errors = list())
    },
    dpprior_schema_error = function(condition) {
      list(valid = FALSE, errors = list(condition))
    }
  )
}


# --- Primitive structural checks -------------------------------------------

.dpprior_schema_has_only_attributes <- function(x,
                                                 allowed = character(),
                                                 required = character()) {
  attrs <- attributes(x)
  attr_names <- if (is.null(attrs)) character() else names(attrs)
  (!is.null(attr_names) || length(attrs) == 0L) &&
    !anyDuplicated(attr_names) &&
    all(attr_names %in% allowed) &&
    all(required %in% attr_names) &&
    (is.null(attrs) || all(vapply(
      attrs,
      function(value) {
        is.atomic(value) && !is.object(value) && is.null(attributes(value))
      },
      logical(1)
    )))
}


.dpprior_schema_validate_attributes <- function(x,
                                                path,
                                                allowed = character(),
                                                required = character()) {
  .dpprior_schema_require(
    .dpprior_schema_has_only_attributes(x, allowed, required),
    "attributes", path,
    if (length(allowed) == 0L) {
      "no attributes"
    } else {
      paste("only canonical attributes", paste(allowed, collapse = ", "))
    },
    names(attributes(x))
  )
  attrs <- attributes(x)
  if (!is.null(attrs)) {
    for (attribute_name in names(attrs)) {
      attribute_value <- attrs[[attribute_name, exact = TRUE]]
      .dpprior_schema_require(
        is.atomic(attribute_value) && !is.object(attribute_value) &&
          is.null(attributes(attribute_value)),
        "attribute_value", paste0(path, ".<", attribute_name, ">"),
        "an ordinary attribute value without nested attributes",
        attribute_value
      )
    }
  }
  invisible(TRUE)
}

.dpprior_schema_is_scalar_character <- function(x, allow_empty = FALSE) {
  is.character(x) && !is.object(x) && is.null(attributes(x)) &&
    length(x) == 1L && !is.na(x) && (allow_empty || nzchar(x))
}


.dpprior_schema_is_scalar_logical <- function(x) {
  is.logical(x) && !is.object(x) && is.null(attributes(x)) &&
    length(x) == 1L && !is.na(x)
}


.dpprior_schema_is_finite_scalar <- function(x) {
  is.numeric(x) && !is.object(x) && is.null(attributes(x)) && length(x) == 1L &&
    !is.na(x) && is.finite(x)
}


.dpprior_schema_is_count <- function(x, minimum = 0L) {
  .dpprior_schema_is_finite_scalar(x) && x == floor(x) &&
    x >= minimum && x <= .Machine$integer.max
}


.dpprior_schema_is_integer_scalar <- function(x) {
  .dpprior_schema_is_finite_scalar(x) && x == floor(x) &&
    abs(x) <= .Machine$integer.max
}


.dpprior_schema_validate_plain_record_value <- function(x, path,
                                                         numeric_only = FALSE) {
  if (is.null(x)) {
    return(invisible(TRUE))
  }
  if (typeof(x) == "list" && is.list(x)) {
    .dpprior_schema_validate_named_list(x, path)
    for (field in names(x)) {
      .dpprior_schema_validate_plain_record_value(
        x[[field, exact = TRUE]], paste0(path, ".", field), numeric_only
      )
    }
    return(invisible(TRUE))
  }
  allowed <- is.numeric(x) || (!numeric_only && (is.logical(x) ||
    is.character(x)))
  allowed_attributes <- if (is.null(names(x))) character() else "names"
  .dpprior_schema_require(
    allowed && !is.object(x) && is.null(dim(x)) &&
      .dpprior_schema_has_only_attributes(x, allowed_attributes) &&
      length(x) > 0L &&
      !anyNA(x) && (!is.numeric(x) || all(is.finite(x))),
    "record_value", path,
    if (numeric_only) "a finite ordinary numeric vector" else
      "a finite ordinary numeric/logical/character vector",
    x
  )
  if (!is.null(names(x))) {
    .dpprior_schema_require(
      !anyDuplicated(names(x)) && all(!is.na(names(x))) &&
        all(nzchar(names(x))),
      "record_names", path, "unique non-empty vector names", names(x)
    )
  }
  invisible(TRUE)
}


.dpprior_schema_validate_plain_data_frame <- function(x, path) {
  .dpprior_schema_require(
    is.data.frame(x) && identical(class(x), "data.frame"),
    "data_frame", path, "an ordinary base data.frame", class(x)
  )
  .dpprior_schema_validate_attributes(
    x, path, c("names", "class", "row.names"),
    c("names", "class", "row.names")
  )
  .dpprior_schema_require(
    !is.null(names(x)) && !anyDuplicated(names(x)) && all(nzchar(names(x))),
    "data_frame_names", path,
    "unique non-empty column names", names(x)
  )
  for (field in names(x)) {
    column <- x[[field, exact = TRUE]]
    allowed_attributes <- if (is.null(names(column))) character() else "names"
    .dpprior_schema_require(
      (is.numeric(column) || is.logical(column) || is.character(column)) &&
        !is.object(column) && is.null(dim(column)) &&
        .dpprior_schema_has_only_attributes(column, allowed_attributes),
      "data_frame_column", paste0(path, ".", field),
      "an ordinary atomic vector column without noncanonical attributes",
      column
    )
    .dpprior_schema_require(
      length(column) == nrow(x), "data_frame_rows",
      paste0(path, ".", field),
      "a column length equal to the data-frame row count", length(column)
    )
  }
  row_names <- attr(x, "row.names", exact = TRUE)
  .dpprior_schema_require(
    (is.integer(row_names) || is.character(row_names)) &&
      !is.object(row_names) && is.null(attributes(row_names)) &&
      length(row.names(x)) == nrow(x) && !anyNA(row.names(x)) &&
      !anyDuplicated(row.names(x)),
    "data_frame_rows", paste0(path, ".<row.names>"),
    "ordinary unique row names coherent with the row count", row_names
  )
  invisible(TRUE)
}


.dpprior_schema_exact_names <- function(x, expected, path) {
  .dpprior_schema_require(
    typeof(x) == "list" && is.list(x) && !is.object(x) && is.null(dim(x)),
    "type", path,
    "an ordinary unclassed list", class(x)
  )
  .dpprior_schema_validate_attributes(
    x, path, "names", if (length(expected) > 0L) "names" else character()
  )
  .dpprior_schema_require(
    !anyDuplicated(names(x)), "duplicate_names", path,
    "unique field names", names(x)
  )
  .dpprior_schema_require(
    identical(names(x), expected), "field_names", path,
    paste(expected, collapse = ", "), names(x)
  )
  invisible(TRUE)
}


.dpprior_schema_validate_named_list <- function(x, path, allow_empty = TRUE) {
  .dpprior_schema_require(
    typeof(x) == "list" && is.list(x) && !is.object(x) && is.null(dim(x)),
    "type", path,
    "an ordinary unclassed named list", class(x)
  )
  .dpprior_schema_validate_attributes(
    x, path, "names",
    if (length(x) > 0L || !allow_empty) "names" else character()
  )
  .dpprior_schema_require(
    !anyDuplicated(names(x)), "duplicate_names", path,
    "unique names", names(x)
  )
  if (length(x) > 0L || !allow_empty) {
    .dpprior_schema_require(
      !is.null(names(x)) && all(!is.na(names(x))) && all(nzchar(names(x))),
      "names", path, "non-empty names for every element", names(x)
    )
  }
  invisible(TRUE)
}


.dpprior_schema_validate_character_vector <- function(x, path,
                                                       named = FALSE) {
  allowed_attributes <- if (named && !is.null(names(x))) "names" else
    character()
  .dpprior_schema_require(
    is.character(x) && !is.object(x) && is.null(dim(x)) &&
      .dpprior_schema_has_only_attributes(x, allowed_attributes) && !anyNA(x),
    "type", path, "a non-missing character vector", x
  )
  if (length(x) > 0L) {
    .dpprior_schema_require(
      all(nzchar(x)), "empty_value", path, "non-empty values", x
    )
  }
  if (named && length(x) > 0L) {
    .dpprior_schema_require(
      !is.null(names(x)) && !anyDuplicated(names(x)) &&
        all(!is.na(names(x))) && all(nzchar(names(x))),
      "names", path, "unique non-empty names", names(x)
    )
  }
  invisible(TRUE)
}


.dpprior_schema_validate_scalar_character <- function(x, path,
                                                       allow_empty = FALSE) {
  .dpprior_schema_require(
    .dpprior_schema_is_scalar_character(x, allow_empty),
    "type", path,
    if (allow_empty) "one non-missing character value" else
      "one non-empty character value",
    x
  )
  invisible(TRUE)
}


.dpprior_schema_validate_scalar_logical <- function(x, path) {
  .dpprior_schema_require(
    .dpprior_schema_is_scalar_logical(x), "type", path,
    "one non-missing logical value", x
  )
  invisible(TRUE)
}


.dpprior_schema_validate_finite_scalar <- function(x, path,
                                                    lower = -Inf,
                                                    upper = Inf,
                                                    lower_open = FALSE,
                                                    upper_open = FALSE) {
  .dpprior_schema_require(
    .dpprior_schema_is_finite_scalar(x), "type", path,
    "one finite numeric value", x
  )
  lower_ok <- if (lower_open) x > lower else x >= lower
  upper_ok <- if (upper_open) x < upper else x <= upper
  .dpprior_schema_require(
    lower_ok && upper_ok, "bounds", path,
    sprintf(
      "%s%s, %s%s",
      if (lower_open) "(" else "[", format(lower), format(upper),
      if (upper_open) ")" else "]"
    ),
    x
  )
  invisible(TRUE)
}


.dpprior_schema_validate_nullable_character <- function(x, path) {
  if (!is.null(x)) {
    .dpprior_schema_validate_scalar_character(x, path)
  }
  invisible(TRUE)
}


.dpprior_schema_validate_count_or_null <- function(x, path,
                                                   minimum = 1L) {
  if (!is.null(x)) {
    .dpprior_schema_require(
      .dpprior_schema_is_count(x, minimum), "count", path,
      sprintf("NULL or an integer at least %d", minimum), x
    )
  }
  invisible(TRUE)
}


# --- Schema and status records ----------------------------------------------

.dpprior_schema <- function(name = "result",
                            version = .DPPRIOR_SCHEMA_VERSION) {
  .dpprior_schema_validate_scalar_character(name, "schema.name")
  .dpprior_schema_require(
    name %in% c("result", "target", "weight-target"),
    "schema_name", "schema.name", "result, target, or weight-target", name
  )
  .dpprior_schema_require(
    .dpprior_schema_is_count(version, 1L), "schema_version", "schema.version",
    "a positive integer", version
  )
  full_name <- switch(
    name,
    result = .DPPRIOR_RESULT_SCHEMA_NAME,
    target = .DPPRIOR_TARGET_SCHEMA_NAME,
    `weight-target` = .DPPRIOR_WEIGHT_TARGET_SCHEMA_NAME
  )
  list(name = full_name, version = as.integer(version))
}


.dpprior_validate_schema_record <- function(x, expected_name, path = "schema") {
  .dpprior_schema_exact_names(x, c("name", "version"), path)
  .dpprior_schema_validate_scalar_character(x$name, paste0(path, ".name"))
  .dpprior_schema_require(
    identical(x$name, expected_name), "schema_name", paste0(path, ".name"),
    expected_name, x$name
  )
  .dpprior_schema_require(
    identical(x$version, .DPPRIOR_SCHEMA_VERSION),
    "schema_version", paste0(path, ".version"),
    sprintf("integer %d", .DPPRIOR_SCHEMA_VERSION), x$version
  )
  invisible(TRUE)
}


.dpprior_validate_status_record <- function(x, path = "status_record") {
  .dpprior_schema_exact_names(
    x, c("status", "usable", "verified", "message"), path
  )
  .dpprior_schema_validate_scalar_character(x$status, paste0(path, ".status"))
  .dpprior_schema_require(
    x$status %in% c("converged", "boundary", "approximate", "infeasible",
                    "failed"),
    "status", paste0(path, ".status"),
    "converged, boundary, approximate, infeasible, or failed", x$status
  )
  .dpprior_schema_validate_scalar_logical(x$usable, paste0(path, ".usable"))
  .dpprior_schema_validate_scalar_logical(
    x$verified, paste0(path, ".verified")
  )
  .dpprior_schema_validate_scalar_character(
    x$message, paste0(path, ".message"), allow_empty = TRUE
  )

  if (x$status %in% c("converged", "boundary")) {
    .dpprior_schema_require(
      x$usable && x$verified, "status_quartet", path,
      "usable=TRUE and verified=TRUE for converged/boundary", x
    )
  }
  if (identical(x$status, "approximate")) {
    .dpprior_schema_require(
      !x$verified, "status_quartet", path,
      "verified=FALSE for approximate", x
    )
  }
  if (identical(x$status, "infeasible")) {
    .dpprior_schema_require(
      !x$usable && x$verified, "status_quartet", path,
      "usable=FALSE and verified=TRUE for certified infeasible", x
    )
  }
  if (identical(x$status, "failed")) {
    .dpprior_schema_require(
      !x$usable && !x$verified, "status_quartet", path,
      "usable=FALSE and verified=FALSE for failed", x
    )
  }
  invisible(TRUE)
}


.dpprior_new_status <- function(status, usable, verified, message = "") {
  out <- list(
    status = status,
    usable = usable,
    verified = verified,
    message = message
  )
  .dpprior_validate_status_record(out)
  out
}


.dpprior_new_parameters <- function(a, b, parameterization) {
  out <- list(a = a, b = b, parameterization = parameterization)
  .dpprior_validate_parameters(out)
  out
}


.dpprior_validate_parameters <- function(x, path = "parameters",
                                         nullable = FALSE) {
  if (is.null(x) && nullable) {
    return(invisible(TRUE))
  }
  .dpprior_schema_exact_names(x, c("a", "b", "parameterization"), path)
  .dpprior_schema_validate_finite_scalar(
    x$a, paste0(path, ".a"), lower = 0, lower_open = TRUE
  )
  .dpprior_schema_validate_finite_scalar(
    x$b, paste0(path, ".b"), lower = 0, lower_open = TRUE
  )
  .dpprior_schema_validate_scalar_character(
    x$parameterization, paste0(path, ".parameterization")
  )
  invisible(TRUE)
}


# --- Computation leaf records ------------------------------------------------

.dpprior_validate_attempt <- function(x, path = "attempt") {
  .dpprior_schema_exact_names(x, .DPPRIOR_ATTEMPT_FIELDS, path)
  for (field in c("id", "stage", "method", "reason_code")) {
    .dpprior_schema_validate_scalar_character(
      x[[field]], paste0(path, ".", field)
    )
  }
  .dpprior_schema_validate_scalar_character(
    x$message, paste0(path, ".message"), allow_empty = TRUE
  )
  .dpprior_schema_validate_scalar_logical(
    x$selected, paste0(path, ".selected")
  )
  .dpprior_schema_validate_character_vector(
    x$warnings, paste0(path, ".warnings")
  )

  for (field in c("start", "bounds")) {
    if (!is.null(x[[field, exact = TRUE]])) {
      .dpprior_schema_validate_plain_record_value(
        x[[field, exact = TRUE]], paste0(path, ".", field), numeric_only = TRUE
      )
    }
  }
  if (!is.null(x$control)) {
    .dpprior_schema_validate_plain_record_value(
      x$control, paste0(path, ".control"), numeric_only = FALSE
    )
  }
  if (!is.null(x$exit_code)) {
    .dpprior_schema_require(
      .dpprior_schema_is_integer_scalar(x$exit_code),
      "exit_code", paste0(path, ".exit_code"),
      "NULL or one finite whole-number scalar", x$exit_code
    )
  }
  if (!is.null(x$evaluations)) {
    .dpprior_schema_validate_named_list(
      x$evaluations, paste0(path, ".evaluations"), allow_empty = FALSE
    )
    for (field in names(x$evaluations)) {
      .dpprior_schema_require(
        .dpprior_schema_is_count(x$evaluations[[field, exact = TRUE]], 0L),
        "evaluations", paste0(path, ".evaluations.", field),
        "one non-negative integer count",
        x$evaluations[[field, exact = TRUE]]
      )
    }
  }
  .dpprior_schema_validate_character_vector(
    x$unavailable, paste0(path, ".unavailable"), named = TRUE
  )
  .dpprior_schema_require(
    all(names(x$unavailable) %in% .DPPRIOR_ATTEMPT_FIELDS),
    "unavailable_path", paste0(path, ".unavailable"),
    "names from the canonical attempt schema", names(x$unavailable)
  )

  if (!is.null(x$iterations)) {
    .dpprior_schema_require(
      .dpprior_schema_is_count(x$iterations, 0L), "count",
      paste0(path, ".iterations"), "NULL or a non-negative integer",
      x$iterations
    )
  }
  if (!is.null(x$elapsed_seconds)) {
    .dpprior_schema_validate_finite_scalar(
      x$elapsed_seconds, paste0(path, ".elapsed_seconds"), lower = 0
    )
  }
  if (!is.null(x$candidate_objective)) {
    .dpprior_schema_validate_finite_scalar(
      x$candidate_objective, paste0(path, ".candidate_objective")
    )
  }
  .dpprior_validate_parameters(
    x$candidate_parameters, paste0(path, ".candidate_parameters"),
    nullable = TRUE
  )
  if (!is.null(x$error)) {
    .dpprior_schema_exact_names(
      x$error, c("class", "code", "message"), paste0(path, ".error")
    )
    for (field in c("class", "code", "message")) {
      .dpprior_schema_validate_scalar_character(
        x$error[[field]], paste0(path, ".error.", field),
        allow_empty = identical(field, "message")
      )
    }
  }

  nullable_evidence <- c(
    "start", "bounds", "control", "exit_code", "iterations", "evaluations",
    "candidate_parameters", "candidate_objective", "elapsed_seconds"
  )
  null_fields <- nullable_evidence[vapply(
    nullable_evidence, function(field) is.null(x[[field]]), logical(1)
  )]
  unavailable_fields <- names(x$unavailable)
  if (is.null(unavailable_fields)) {
    unavailable_fields <- character()
  }
  .dpprior_schema_require(
    identical(sort(unavailable_fields), sort(null_fields)),
    "missing_unavailable_reason", paste0(path, ".unavailable"),
    paste(
      "exactly one named reason for every and only every NULL evidence field"
    ),
    list(NULL_fields = null_fields, declared = unavailable_fields)
  )
  if (x$selected) {
    .dpprior_schema_require(
      !is.null(x$candidate_parameters), "selected_candidate", path,
      "selected attempt with candidate parameters", NULL
    )
  }
  invisible(TRUE)
}


.dpprior_new_attempt <- function(id,
                                 stage,
                                 method,
                                 start = NULL,
                                 bounds = NULL,
                                 control = NULL,
                                 exit_code = NULL,
                                 message = "",
                                 iterations = NULL,
                                 evaluations = NULL,
                                 candidate_parameters = NULL,
                                 candidate_objective = NULL,
                                 elapsed_seconds = NULL,
                                 warnings = character(),
                                 error = NULL,
                                 selected = FALSE,
                                 reason_code,
                                 unavailable = character()) {
  out <- list(
    id = id,
    stage = stage,
    method = method,
    start = start,
    bounds = bounds,
    control = control,
    exit_code = exit_code,
    message = message,
    iterations = iterations,
    evaluations = evaluations,
    candidate_parameters = candidate_parameters,
    candidate_objective = candidate_objective,
    elapsed_seconds = elapsed_seconds,
    warnings = warnings,
    error = error,
    selected = selected,
    reason_code = reason_code,
    unavailable = unavailable
  )
  .dpprior_validate_attempt(out)
  out
}


.dpprior_validate_setting_record <- function(x, path) {
  .dpprior_schema_exact_names(
    x, c("method", "controls", "parameterization"), path
  )
  .dpprior_schema_validate_scalar_character(x$method, paste0(path, ".method"))
  .dpprior_schema_validate_named_list(x$controls, paste0(path, ".controls"))
  .dpprior_schema_validate_plain_record_value(
    x$controls, paste0(path, ".controls")
  )
  .dpprior_schema_validate_scalar_character(
    x$parameterization, paste0(path, ".parameterization")
  )
  invisible(TRUE)
}


.dpprior_validate_orders <- function(x, path = "computation.orders") {
  fields <- c(
    "M_requested", "M_selected", "M_verification_required",
    "M_verification_used", "requested_reason", "selected_reason",
    "verification_required_reason", "verification_used_reason"
  )
  .dpprior_schema_exact_names(x, fields, path)
  values <- fields[seq_len(4L)]
  reasons <- fields[5:8]
  for (i in seq_along(values)) {
    .dpprior_schema_validate_count_or_null(
      x[[values[[i]]]], paste0(path, ".", values[[i]])
    )
    .dpprior_schema_validate_scalar_character(
      x[[reasons[[i]]]], paste0(path, ".", reasons[[i]])
    )
  }
  if (!is.null(x$M_verification_required) &&
      !is.null(x$M_verification_used)) {
    .dpprior_schema_require(
      x$M_verification_used >= x$M_verification_required,
      "verification_order", paste0(path, ".M_verification_used"),
      "an order at least M_verification_required", x$M_verification_used
    )
  }
  invisible(TRUE)
}


.dpprior_new_orders <- function(M_requested = NULL,
                                M_selected = NULL,
                                M_verification_required = NULL,
                                M_verification_used = NULL,
                                requested_reason,
                                selected_reason,
                                verification_required_reason,
                                verification_used_reason) {
  out <- list(
    M_requested = M_requested,
    M_selected = M_selected,
    M_verification_required = M_verification_required,
    M_verification_used = M_verification_used,
    requested_reason = requested_reason,
    selected_reason = selected_reason,
    verification_required_reason = verification_required_reason,
    verification_used_reason = verification_used_reason
  )
  .dpprior_validate_orders(out)
  out
}


.dpprior_validate_scaling <- function(x, path = "computation.scaling") {
  .dpprior_schema_exact_names(
    x,
    c("requested", "used", "formula", "values", "fixed_from_input",
      "change_reason"),
    path
  )
  for (field in c("requested", "used")) {
    if (!is.null(x[[field]])) {
      .dpprior_schema_validate_named_list(
        x[[field]], paste0(path, ".", field)
      )
      .dpprior_schema_validate_plain_record_value(
        x[[field]], paste0(path, ".", field)
      )
    }
  }
  .dpprior_schema_validate_scalar_character(
    x$formula, paste0(path, ".formula")
  )
  .dpprior_schema_validate_named_list(x$values, paste0(path, ".values"))
  .dpprior_schema_validate_plain_record_value(x$values, paste0(path, ".values"))
  .dpprior_schema_validate_scalar_logical(
    x$fixed_from_input, paste0(path, ".fixed_from_input")
  )
  .dpprior_schema_validate_scalar_character(
    x$change_reason, paste0(path, ".change_reason"), allow_empty = TRUE
  )
  if (!identical(x$requested, x$used)) {
    .dpprior_schema_require(
      nzchar(x$change_reason), "scaling_change", paste0(path, ".change_reason"),
      "a non-empty reason when requested and used scaling differ",
      x$change_reason
    )
  }
  invisible(TRUE)
}


.dpprior_new_scaling <- function(requested = NULL,
                                 used = NULL,
                                 formula = "not_applicable",
                                 values = list(),
                                 fixed_from_input = FALSE,
                                 change_reason = "") {
  out <- list(
    requested = requested,
    used = used,
    formula = formula,
    values = values,
    fixed_from_input = fixed_from_input,
    change_reason = change_reason
  )
  .dpprior_validate_scaling(out)
  out
}


.dpprior_validate_fallback <- function(x, path = "computation.fallback") {
  .dpprior_schema_exact_names(
    x,
    c("attempted", "used", "trigger_attempt_id", "selected_attempt_id",
      "reason_code", "message", "outcome"),
    path
  )
  .dpprior_schema_validate_scalar_logical(x$attempted, paste0(path, ".attempted"))
  .dpprior_schema_validate_scalar_logical(x$used, paste0(path, ".used"))
  for (field in c("trigger_attempt_id", "selected_attempt_id", "reason_code")) {
    .dpprior_schema_validate_nullable_character(
      x[[field]], paste0(path, ".", field)
    )
  }
  .dpprior_schema_validate_scalar_character(
    x$message, paste0(path, ".message"), allow_empty = TRUE
  )
  .dpprior_schema_validate_scalar_character(x$outcome, paste0(path, ".outcome"))
  if (x$used) {
    .dpprior_schema_require(
      x$attempted && !is.null(x$trigger_attempt_id) &&
        !is.null(x$selected_attempt_id) && !is.null(x$reason_code) &&
        identical(x$outcome, "selected"),
      "fallback", path,
      paste(
        "used fallback with attempted=TRUE, trigger/selected/reason IDs,",
        "and outcome=selected"
      ), x
    )
  }
  if (x$attempted && !x$used) {
    .dpprior_schema_require(
      !is.null(x$trigger_attempt_id) && is.null(x$selected_attempt_id) &&
        !is.null(x$reason_code) &&
        x$outcome %in% c("failed", "rejected", "attempted_not_selected"),
      "fallback", path,
      paste(
        "an unused attempted fallback with trigger/reason, no selected ID,",
        "and a closed failed/rejected outcome"
      ), x
    )
  }
  if (!x$attempted) {
    .dpprior_schema_require(
      !x$used && is.null(x$trigger_attempt_id) &&
        is.null(x$selected_attempt_id) && is.null(x$reason_code) &&
        identical(x$outcome, "not_attempted"),
      "fallback", path,
      "no fallback IDs/reason and outcome=not_attempted when not attempted", x
    )
  }
  invisible(TRUE)
}


.dpprior_new_fallback <- function(attempted = FALSE,
                                  used = FALSE,
                                  trigger_attempt_id = NULL,
                                  selected_attempt_id = NULL,
                                  reason_code = NULL,
                                  message = "",
                                  outcome = "not_attempted") {
  out <- list(
    attempted = attempted,
    used = used,
    trigger_attempt_id = trigger_attempt_id,
    selected_attempt_id = selected_attempt_id,
    reason_code = reason_code,
    message = message,
    outcome = outcome
  )
  .dpprior_validate_fallback(out)
  out
}


.dpprior_validate_termination <- function(x, path = "computation.termination") {
  .dpprior_schema_exact_names(
    x, c("code", "message", "source", "iterations", "boundary_reason"), path
  )
  .dpprior_schema_validate_scalar_character(x$code, paste0(path, ".code"))
  .dpprior_schema_validate_scalar_character(
    x$message, paste0(path, ".message"), allow_empty = TRUE
  )
  .dpprior_schema_validate_scalar_character(x$source, paste0(path, ".source"))
  if (!is.null(x$iterations)) {
    .dpprior_schema_require(
      .dpprior_schema_is_count(x$iterations, 0L), "count",
      paste0(path, ".iterations"), "NULL or a non-negative integer",
      x$iterations
    )
  }
  .dpprior_schema_validate_nullable_character(
    x$boundary_reason, paste0(path, ".boundary_reason")
  )
  invisible(TRUE)
}


.dpprior_new_termination <- function(code,
                                     message = "",
                                     source,
                                     iterations = NULL,
                                     boundary_reason = NULL) {
  out <- list(
    code = code,
    message = message,
    source = source,
    iterations = iterations,
    boundary_reason = boundary_reason
  )
  .dpprior_validate_termination(out)
  out
}


.dpprior_validate_candidate_evaluation <- function(
    x, path = "candidate_evaluation") {
  .dpprior_schema_exact_names(
    x, .DPPRIOR_CANDIDATE_EVALUATION_FIELDS, path
  )
  for (field in c(
    "id", "method", "generator", "objective_kind", "outcome", "source"
  )) {
    .dpprior_schema_validate_scalar_character(
      x[[field, exact = TRUE]], paste0(path, ".", field)
    )
  }
  .dpprior_schema_validate_nullable_character(
    x[["attempt_id", exact = TRUE]], paste0(path, ".attempt_id")
  )
  .dpprior_schema_require(
    x[["generator", exact = TRUE]] %in% c(
      "direct_attempt", "deterministic_profile_scan",
      "input_fit", "feasibility_extreme", "initialization",
      "derived_diagnostic_attempt"
    ),
    "candidate_generator", paste0(path, ".generator"),
    paste(
      "direct_attempt, deterministic_profile_scan,",
      "input_fit, feasibility_extreme, initialization, or",
      "derived_diagnostic_attempt"
    ),
    x[["generator", exact = TRUE]]
  )
  .dpprior_schema_require(
    !is.null(x[["attempt_id", exact = TRUE]]) ||
      x[["generator", exact = TRUE]] %in%
        c("deterministic_profile_scan", "input_fit", "feasibility_extreme"),
    "candidate_generator", paste0(path, ".attempt_id"),
    "a parent attempt ID except for authorized attempt-free generators",
    x[["attempt_id", exact = TRUE]]
  )
  .dpprior_schema_require(
    x[["objective_kind", exact = TRUE]] %in% c(
      "standardized_residual", "kl", "K_loss", "weight_metric_probe",
      "penalized_diagnostic", "soft_tradeoff"
    ),
    "candidate_objective_kind", paste0(path, ".objective_kind"),
    paste(
      "standardized_residual, kl, K_loss, weight_metric_probe,",
      "penalized_diagnostic, or soft_tradeoff"
    ),
    x[["objective_kind", exact = TRUE]]
  )
  for (field in c("recorded_objective_kind", "selection_objective_kind")) {
    .dpprior_schema_validate_nullable_character(
      x[[field, exact = TRUE]], paste0(path, ".", field)
    )
    if (!is.null(x[[field, exact = TRUE]])) {
      .dpprior_schema_require(
        x[[field, exact = TRUE]] %in% c(
          "standardized_residual", "kl", "K_loss", "weight_metric_probe",
          "penalized_diagnostic", "soft_tradeoff"
        ),
        "candidate_objective_kind", paste0(path, ".", field),
        "a closed canonical objective kind", x[[field, exact = TRUE]]
      )
    }
  }
  .dpprior_validate_parameters(
    x[["parameters", exact = TRUE]], paste0(path, ".parameters")
  )
  .dpprior_schema_validate_scalar_logical(
    x[["recorded_objective_available", exact = TRUE]],
    paste0(path, ".recorded_objective_available")
  )
  if (x[["recorded_objective_available", exact = TRUE]]) {
    .dpprior_schema_validate_finite_scalar(
      x[["recorded_objective", exact = TRUE]],
      paste0(path, ".recorded_objective"), lower = 0
    )
    same_objective_kind <- identical(
      x[["recorded_objective_kind", exact = TRUE]],
      x[["objective_kind", exact = TRUE]]
    )
    expected_reason <- if (same_objective_kind) NULL else paste0(
      "source_objective_kind_mismatch:",
      x[["recorded_objective_kind", exact = TRUE]], "->",
      x[["objective_kind", exact = TRUE]]
    )
    .dpprior_schema_require(
      !is.null(x[["recorded_objective_kind", exact = TRUE]]) &&
        identical(x[["recorded_objective_reason", exact = TRUE]],
                  expected_reason),
      "candidate_recorded_objective", paste0(path, ".recorded_objective_reason"),
      "NULL for comparable kinds or the exact closed kind-mismatch reason",
      x[["recorded_objective_reason", exact = TRUE]]
    )
    if (!is.null(x[["objective_passed", exact = TRUE]])) {
      .dpprior_schema_validate_scalar_logical(
        x[["objective_passed", exact = TRUE]], paste0(path, ".objective_passed")
      )
    }
  } else {
    .dpprior_schema_require(
      is.null(x[["recorded_objective", exact = TRUE]]) &&
        is.null(x[["recorded_objective_kind", exact = TRUE]]) &&
        is.null(x[["objective_passed", exact = TRUE]]),
      "candidate_recorded_objective", path,
      "NULL recorded objective and objective_passed when unavailable",
      list(
        objective = x[["recorded_objective", exact = TRUE]],
        passed = x[["objective_passed", exact = TRUE]]
      )
    )
    .dpprior_schema_validate_scalar_character(
      x[["recorded_objective_reason", exact = TRUE]],
      paste0(path, ".recorded_objective_reason")
    )
  }
  .dpprior_schema_validate_finite_scalar(
    x[["fresh_objective", exact = TRUE]], paste0(path, ".fresh_objective"),
    lower = 0
  )
  if (!is.null(x[["selection_objective", exact = TRUE]])) {
    .dpprior_schema_validate_finite_scalar(
      x[["selection_objective", exact = TRUE]],
      paste0(path, ".selection_objective"), lower = 0
    )
    .dpprior_schema_require(
      identical(x[["selection_objective_kind", exact = TRUE]],
                x[["objective_kind", exact = TRUE]]),
      "candidate_selection_objective_kind",
      paste0(path, ".selection_objective_kind"),
      "identity with the freshly recomputed objective kind",
      x[["selection_objective_kind", exact = TRUE]]
    )
  } else {
    .dpprior_schema_require(
      is.null(x[["selection_objective_kind", exact = TRUE]]),
      "candidate_selection_objective_kind",
      paste0(path, ".selection_objective_kind"),
      "NULL when selection_objective is unavailable",
      x[["selection_objective_kind", exact = TRUE]]
    )
  }
  .dpprior_schema_validate_finite_scalar(
    x[["objective_tolerance", exact = TRUE]],
    paste0(path, ".objective_tolerance"), lower = 0
  )
  for (field in c(
    "execution_success", "optimizer_supported", "scientifically_eligible",
    "selection_eligible", "diagnostic_eligible", "decision_eligible",
    "selected"
  )) {
    .dpprior_schema_validate_scalar_logical(
      x[[field, exact = TRUE]], paste0(path, ".", field)
    )
  }
  .dpprior_validate_snapshot(
    x[["selected_snapshot", exact = TRUE]],
    paste0(path, ".selected_snapshot"), nullable = FALSE
  )
  .dpprior_validate_snapshot(
    x[["verifier_snapshot", exact = TRUE]],
    paste0(path, ".verifier_snapshot")
  )
  .dpprior_schema_require(
    x[["selected_snapshot", exact = TRUE]][["finite", exact = TRUE]] &&
      identical(
        x[["selected_snapshot", exact = TRUE]][["parameters", exact = TRUE]],
        x[["parameters", exact = TRUE]]
      ) &&
      (is.null(x[["verifier_snapshot", exact = TRUE]]) || (
        x[["verifier_snapshot", exact = TRUE]][["finite", exact = TRUE]] &&
          identical(
            x[["verifier_snapshot", exact = TRUE]][[
              "parameters", exact = TRUE
            ]],
            x[["parameters", exact = TRUE]]
          )
      )),
    "candidate_snapshot", path,
    paste(
      "a finite selected snapshot and, when retained, a finite fixed-",
      "candidate verifier snapshot"
    ),
    list(
      selected = x[["selected_snapshot", exact = TRUE]],
      verifier = x[["verifier_snapshot", exact = TRUE]]
    )
  )
  .dpprior_schema_validate_named_list(
    x[["checks", exact = TRUE]], paste0(path, ".checks"), allow_empty = FALSE
  )
  .dpprior_schema_require(
    length(x[["checks", exact = TRUE]]) > 0L,
    "candidate_checks", paste0(path, ".checks"),
    "at least one substantive scientific eligibility check", NULL
  )
  for (check_name in names(x[["checks", exact = TRUE]])) {
    check <- x[["checks", exact = TRUE]][[check_name, exact = TRUE]]
    .dpprior_validate_decision_check(
      check, paste0(path, ".checks.", check_name)
    )
    .dpprior_schema_require(
      identical(check[["source", exact = TRUE]],
                paste0("candidate:", x[["id", exact = TRUE]])),
      "candidate_check_source", paste0(path, ".checks.", check_name, ".source"),
      "the exact owning candidate-evaluation ID", check[["source", exact = TRUE]]
    )
  }
  .dpprior_schema_validate_character_vector(
    x[["rejection_codes", exact = TRUE]], paste0(path, ".rejection_codes")
  )
  .dpprior_schema_require(
    !anyDuplicated(x[["rejection_codes", exact = TRUE]]),
    "candidate_rejection_codes", paste0(path, ".rejection_codes"),
    "unique rejection codes", x[["rejection_codes", exact = TRUE]]
  )
  objective_kinds_comparable <-
    x[["recorded_objective_available", exact = TRUE]] && identical(
      x[["recorded_objective_kind", exact = TRUE]],
      x[["objective_kind", exact = TRUE]]
    )
  expected_objective_pass <- if (objective_kinds_comparable
  ) {
    abs(x[["recorded_objective", exact = TRUE]] -
          x[["fresh_objective", exact = TRUE]]) <=
      x[["objective_tolerance", exact = TRUE]]
  } else {
    NULL
  }
  objective_consistent <- is.null(expected_objective_pass) ||
    expected_objective_pass
  check_passes <- vapply(
    x[["checks", exact = TRUE]],
    function(check) isTRUE(check[["passed", exact = TRUE]]),
    logical(1)
  )
  expected_scientific <- all(check_passes)
  expected_selection <- objective_consistent && expected_scientific &&
    !is.null(x[["selection_objective", exact = TRUE]])
  expected_decision <- expected_selection &&
    x[["optimizer_supported", exact = TRUE]]
  expected_rejections <- c(
    if (!x[["execution_success", exact = TRUE]]) "execution_failed",
    if (!is.null(expected_objective_pass) && !expected_objective_pass)
      "objective_mismatch",
    if (x[["recorded_objective_available", exact = TRUE]] &&
        !objective_kinds_comparable) "recorded_objective_noncomparable",
    if (any(!check_passes)) {
      paste0("check_failed:", names(check_passes)[!check_passes])
    } else {
      character()
    },
    if (is.null(x[["selection_objective", exact = TRUE]]))
      "noncomparable_objective",
    if (!x[["optimizer_supported", exact = TRUE]]) "optimizer_unsupported"
  )
  expected_outcome <- if (x[["selected", exact = TRUE]]) {
    if (expected_decision) "selected" else "selected_diagnostic"
  } else if (expected_selection) {
    "eligible_not_selected"
  } else {
    "rejected"
  }
  .dpprior_schema_require(
    identical(x[["objective_passed", exact = TRUE]], expected_objective_pass) &&
      (!x[["optimizer_supported", exact = TRUE]] ||
         x[["execution_success", exact = TRUE]]) &&
      identical(x[["scientifically_eligible", exact = TRUE]],
                expected_scientific) &&
      identical(x[["selection_eligible", exact = TRUE]],
                expected_selection) &&
      (!x[["diagnostic_eligible", exact = TRUE]] || (
        !expected_selection &&
          !is.null(x[["selection_objective", exact = TRUE]]) &&
          !is.null(x[["verifier_snapshot", exact = TRUE]])
      )) &&
      (!x[["selected", exact = TRUE]] ||
         expected_selection || x[["diagnostic_eligible", exact = TRUE]]) &&
      identical(x[["decision_eligible", exact = TRUE]], expected_decision) &&
      identical(x[["outcome", exact = TRUE]], expected_outcome) &&
      identical(x[["rejection_codes", exact = TRUE]], expected_rejections),
    "candidate_evaluation_truth", path,
    paste(
      "objective/scientific/decision eligibility, outcome, and ordered",
      "rejection codes recomputed from retained evidence"
    ),
    x
  )
  invisible(TRUE)
}


.dpprior_new_candidate_evaluation <- function(
    id,
    attempt_id,
    method,
    generator = if (is.null(attempt_id)) "input_fit" else "direct_attempt",
    parameters,
    objective_kind,
    recorded_objective_kind = if (is.null(recorded_objective)) {
      NULL
    } else {
      objective_kind
    },
    selection_objective_kind = if (is.null(selection_objective)) {
      NULL
    } else {
      objective_kind
    },
    recorded_objective = NULL,
    recorded_objective_reason = NULL,
    fresh_objective,
    selection_objective = fresh_objective,
    objective_tolerance,
    selected_snapshot,
    verifier_snapshot = NULL,
    checks,
    execution_success,
    optimizer_supported = execution_success,
    diagnostic_eligible = FALSE,
    selected,
    source = "canonical_candidate_evaluation") {
  recorded_objective_available <- !is.null(recorded_objective)
  if (is.null(recorded_objective_reason)) {
    recorded_objective_reason <- if (!recorded_objective_available) {
      "not_recorded_by_source_attempt"
    } else if (!identical(recorded_objective_kind, objective_kind)) {
      paste0(
        "source_objective_kind_mismatch:", recorded_objective_kind, "->",
        objective_kind
      )
    } else {
      NULL
    }
  }
  objective_passed <- if (recorded_objective_available &&
                          identical(recorded_objective_kind, objective_kind)) {
    abs(recorded_objective - fresh_objective) <= objective_tolerance
  } else {
    NULL
  }
  objective_consistent <- is.null(objective_passed) || objective_passed
  check_passes <- vapply(
    checks, function(check) isTRUE(check[["passed", exact = TRUE]]), logical(1)
  )
  scientifically_eligible <- length(check_passes) > 0L && all(check_passes)
  selection_eligible <- objective_consistent && scientifically_eligible &&
    !is.null(selection_objective)
  decision_eligible <- selection_eligible && optimizer_supported
  rejection_codes <- c(
    if (!execution_success) "execution_failed",
    if (!is.null(objective_passed) && !objective_passed) "objective_mismatch",
    if (recorded_objective_available &&
        !identical(recorded_objective_kind, objective_kind)) {
      "recorded_objective_noncomparable"
    },
    if (any(!check_passes)) {
      paste0("check_failed:", names(check_passes)[!check_passes])
    } else {
      character()
    },
    if (is.null(selection_objective)) "noncomparable_objective",
    if (!optimizer_supported) "optimizer_unsupported"
  )
  outcome <- if (selected) {
    if (decision_eligible) "selected" else "selected_diagnostic"
  } else if (selection_eligible) {
    "eligible_not_selected"
  } else {
    "rejected"
  }
  out <- list(
    id = id,
    attempt_id = attempt_id,
    method = method,
    generator = generator,
    parameters = parameters,
    objective_kind = objective_kind,
    recorded_objective_kind = recorded_objective_kind,
    selection_objective_kind = selection_objective_kind,
    recorded_objective = recorded_objective,
    recorded_objective_available = recorded_objective_available,
    recorded_objective_reason = recorded_objective_reason,
    fresh_objective = fresh_objective,
    selection_objective = selection_objective,
    objective_tolerance = objective_tolerance,
    objective_passed = objective_passed,
    selected_snapshot = selected_snapshot,
    verifier_snapshot = verifier_snapshot,
    checks = checks,
    execution_success = execution_success,
    optimizer_supported = optimizer_supported,
    scientifically_eligible = scientifically_eligible,
    selection_eligible = selection_eligible,
    diagnostic_eligible = diagnostic_eligible,
    decision_eligible = decision_eligible,
    selected = selected,
    outcome = outcome,
    rejection_codes = rejection_codes,
    source = source
  )
  .dpprior_validate_candidate_evaluation(out)
  out
}


.dpprior_validate_computation <- function(x, path = "computation") {
  .dpprior_schema_exact_names(x, .DPPRIOR_COMPUTATION_FIELDS, path)
  .dpprior_validate_setting_record(x$request, paste0(path, ".request"))
  .dpprior_validate_setting_record(x$used, paste0(path, ".used"))
  .dpprior_validate_orders(x$orders, paste0(path, ".orders"))
  .dpprior_validate_scaling(x$scaling, paste0(path, ".scaling"))
  .dpprior_schema_require(
    typeof(x$attempts) == "list" && is.list(x$attempts) &&
      !is.object(x$attempts) &&
      is.null(dim(x$attempts)) && .dpprior_schema_has_only_attributes(
        x$attempts,
        if (is.null(names(x$attempts))) character() else "names"
      ),
    "type", paste0(path, ".attempts"),
    "an ordinary unclassed ordered list", class(x$attempts)
  )
  for (i in seq_along(x$attempts)) {
    .dpprior_validate_attempt(
      x$attempts[[i]], sprintf("%s.attempts[[%d]]", path, i)
    )
  }
  .dpprior_schema_require(
    typeof(x$candidate_evaluations) == "list" &&
      is.list(x$candidate_evaluations) &&
      !is.object(x$candidate_evaluations) &&
      is.null(dim(x$candidate_evaluations)) &&
      .dpprior_schema_has_only_attributes(
        x$candidate_evaluations,
        if (is.null(names(x$candidate_evaluations))) character() else "names"
      ),
    "type", paste0(path, ".candidate_evaluations"),
    "an ordinary unclassed ordered candidate-evaluation ledger",
    class(x$candidate_evaluations)
  )
  for (i in seq_along(x$candidate_evaluations)) {
    .dpprior_validate_candidate_evaluation(
      x$candidate_evaluations[[i]],
      sprintf("%s.candidate_evaluations[[%d]]", path, i)
    )
  }
  evaluation_ids <- vapply(
    x$candidate_evaluations,
    function(evaluation) evaluation[["id", exact = TRUE]], character(1)
  )
  evaluation_attempt_ids <- vapply(
    x$candidate_evaluations,
    function(evaluation) evaluation[["attempt_id", exact = TRUE]] %||% "",
    character(1)
  )
  .dpprior_schema_require(
    !anyDuplicated(evaluation_ids),
    "candidate_evaluation_ids", paste0(path, ".candidate_evaluations"),
    "unique evaluation IDs (multiple candidates may originate in one attempt)",
    list(ids = evaluation_ids, attempt_ids = evaluation_attempt_ids)
  )
  selected_evaluations <- vapply(
    x$candidate_evaluations,
    function(evaluation) isTRUE(evaluation[["selected", exact = TRUE]]),
    logical(1)
  )
  .dpprior_schema_require(
    sum(selected_evaluations) <= 1L,
    "selected_candidate", paste0(path, ".candidate_evaluations"),
    "at most one selected candidate evaluation",
    evaluation_ids[selected_evaluations]
  )
  .dpprior_schema_validate_nullable_character(
    x$selected_candidate_id, paste0(path, ".selected_candidate_id")
  )
  .dpprior_schema_require(
    if (sum(selected_evaluations) == 1L) {
      identical(x$selected_candidate_id, evaluation_ids[selected_evaluations])
    } else {
      is.null(x$selected_candidate_id)
    },
    "selected_candidate", paste0(path, ".selected_candidate_id"),
    "the sole selected candidate ID, or NULL when none is selected",
    x$selected_candidate_id
  )
  ids <- vapply(x$attempts, function(attempt) attempt$id, character(1))
  .dpprior_schema_require(
    !anyDuplicated(ids), "attempt_ids", paste0(path, ".attempts"),
    "unique attempt IDs", ids
  )
  selected <- vapply(
    x$attempts, function(attempt) isTRUE(attempt$selected), logical(1)
  )
  .dpprior_schema_require(
    sum(selected) <= 1L, "selected_attempt", paste0(path, ".attempts"),
    "at most one selected attempt", ids[selected]
  )
  .dpprior_schema_validate_nullable_character(
    x$selected_attempt_id, paste0(path, ".selected_attempt_id")
  )
  if (sum(selected) == 1L) {
    .dpprior_schema_require(
      identical(x$selected_attempt_id, ids[selected]),
      "selected_attempt", paste0(path, ".selected_attempt_id"),
      "the ID of the sole selected attempt", x$selected_attempt_id
    )
  } else {
    .dpprior_schema_require(
      is.null(x$selected_attempt_id),
      "selected_attempt", paste0(path, ".selected_attempt_id"),
      "NULL when no attempt is selected", x$selected_attempt_id
    )
  }

  .dpprior_validate_fallback(x$fallback, paste0(path, ".fallback"))
  if (x$fallback$attempted) {
    trigger_index <- match(x$fallback$trigger_attempt_id, ids)
    later_indices <- if (!is.na(trigger_index) && trigger_index < length(ids)) {
      seq.int(trigger_index + 1L, length(ids))
    } else {
      integer()
    }
    fallback_stage_indices <- later_indices[vapply(
      x$attempts[later_indices],
      function(attempt) attempt[["stage", exact = TRUE]] %in%
        c("fallback", "recovery"),
      logical(1)
    )]
    .dpprior_schema_require(
      !is.na(trigger_index) && length(fallback_stage_indices) > 0L,
      "fallback_trigger", paste0(path, ".fallback.trigger_attempt_id"),
      paste(
        "an existing trigger followed by at least one ordered",
        "fallback/recovery-stage attempt"
      ), x$fallback$trigger_attempt_id
    )
  }
  if (x$fallback$used) {
    trigger_index <- match(x$fallback$trigger_attempt_id, ids)
    selected_index <- match(x$fallback$selected_attempt_id, ids)
    .dpprior_schema_require(
      !is.na(selected_index) && selected_index > trigger_index &&
        selected_index %in% fallback_stage_indices &&
        identical(x$fallback$selected_attempt_id, x$selected_attempt_id),
      "fallback_order", paste0(path, ".fallback"),
      "a later selected fallback attempt matching selected_attempt_id", x$fallback
    )
    .dpprior_schema_require(
      !isTRUE(x$attempts[[trigger_index]]$selected),
      "fallback_trigger", paste0(path, ".fallback.trigger_attempt_id"),
      "a non-selected failed or rejected primary attempt",
      x$fallback$trigger_attempt_id
    )
  }

  settings_changed <- !identical(x$request$method, x$used$method) ||
    !identical(x$request$controls, x$used$controls) ||
    !identical(x$request$parameterization, x$used$parameterization)
  if (settings_changed) {
    .dpprior_schema_require(
      x$fallback$attempted && !is.null(x$fallback$reason_code),
      "requested_used", path,
      "explicit fallback evidence for changed requested/used settings",
      list(request = x$request, used = x$used)
    )
  }

  .dpprior_validate_termination(x$termination, paste0(path, ".termination"))
  if (!is.null(x$trace)) {
    .dpprior_schema_validate_plain_data_frame(
      x[["trace", exact = TRUE]], paste0(path, ".trace")
    )
  }
  .dpprior_schema_validate_named_list(
    x$resources, paste0(path, ".resources")
  )
  .dpprior_schema_validate_plain_record_value(
    x[["resources", exact = TRUE]], paste0(path, ".resources")
  )
  invisible(TRUE)
}


.dpprior_new_computation <- function(request,
                                     used,
                                     orders,
                                     scaling,
                                     attempts = list(),
                                     candidate_evaluations = list(),
                                     selected_candidate_id = NULL,
                                     selected_attempt_id = NULL,
                                     fallback = .dpprior_new_fallback(),
                                     termination,
                                     trace = NULL,
                                     resources = list()) {
  out <- list(
    request = request,
    used = used,
    orders = orders,
    scaling = scaling,
    attempts = attempts,
    candidate_evaluations = candidate_evaluations,
    selected_candidate_id = selected_candidate_id,
    selected_attempt_id = selected_attempt_id,
    fallback = fallback,
    termination = termination,
    trace = trace,
    resources = resources
  )
  .dpprior_validate_computation(out)
  out
}


# --- Verification -----------------------------------------------------------

.dpprior_validate_snapshot <- function(x, path = "snapshot", nullable = TRUE) {
  if (is.null(x) && nullable) {
    return(invisible(TRUE))
  }
  .dpprior_schema_exact_names(
    x,
    c("parameters", "M", "achieved", "residuals", "tolerances", "finite",
      "source"),
    path
  )
  .dpprior_validate_parameters(
    x$parameters, paste0(path, ".parameters"), nullable = TRUE
  )
  .dpprior_schema_validate_count_or_null(x$M, paste0(path, ".M"))
  for (field in c("achieved", "residuals", "tolerances")) {
    .dpprior_schema_validate_named_list(
      x[[field]], paste0(path, ".", field)
    )
    .dpprior_schema_validate_plain_record_value(
      x[[field]], paste0(path, ".", field)
    )
  }
  .dpprior_schema_validate_scalar_logical(x$finite, paste0(path, ".finite"))
  .dpprior_schema_validate_scalar_character(x$source, paste0(path, ".source"))
  invisible(TRUE)
}


.dpprior_new_snapshot <- function(parameters,
                                  M,
                                  achieved,
                                  residuals,
                                  tolerances,
                                  finite,
                                  source) {
  out <- list(
    parameters = parameters,
    M = M,
    achieved = achieved,
    residuals = residuals,
    tolerances = tolerances,
    finite = finite,
    source = source
  )
  .dpprior_validate_snapshot(out, nullable = FALSE)
  out
}


.dpprior_schema_check_pass <- function(x) {
  if (.dpprior_schema_is_scalar_logical(x)) {
    return(isTRUE(x))
  }
  if (typeof(x) == "list" && is.list(x) && !is.object(x) && is.null(dim(x)) &&
      !anyDuplicated(names(x)) &&
      sum(names(x) == "passed") == 1L) {
    passed <- x[["passed", exact = TRUE]]
    if (.dpprior_schema_is_scalar_logical(passed)) {
      return(isTRUE(passed))
    }
  }
  NA
}


.dpprior_validate_check_list <- function(x, path) {
  .dpprior_schema_validate_named_list(x, path)
  if (length(x) > 0L) {
    passes <- vapply(x, .dpprior_schema_check_pass, logical(1))
    .dpprior_schema_require(
      !anyNA(passes), "check_shape", path,
      "logical checks or records with one logical `passed` field", x
    )
  }
  invisible(TRUE)
}


.dpprior_validate_decision_check <- function(x, path = "decision_check") {
  .dpprior_schema_exact_names(x, .DPPRIOR_DECISION_CHECK_FIELDS, path)
  .dpprior_schema_validate_scalar_logical(
    x[["passed", exact = TRUE]], paste0(path, ".passed")
  )
  .dpprior_schema_require(
    !is.null(x[["value", exact = TRUE]]), "decision_evidence",
    paste0(path, ".value"), "retained substantive check evidence", NULL
  )
  .dpprior_schema_validate_plain_record_value(
    x[["value", exact = TRUE]], paste0(path, ".value")
  )
  .dpprior_schema_validate_plain_record_value(
    x[["reference", exact = TRUE]], paste0(path, ".reference")
  )
  .dpprior_schema_validate_plain_record_value(
    x[["tolerance", exact = TRUE]], paste0(path, ".tolerance"),
    numeric_only = TRUE
  )
  .dpprior_schema_validate_scalar_character(
    x[["operator", exact = TRUE]], paste0(path, ".operator")
  )
  .dpprior_schema_require(
    x[["operator", exact = TRUE]] %in%
      c("abs_lte", "lte", "gte", "lt", "gt", "identical"),
    "decision_operator", paste0(path, ".operator"),
    "abs_lte, lte, gte, lt, gt, or identical",
    x[["operator", exact = TRUE]]
  )
  .dpprior_schema_validate_scalar_character(
    x[["source", exact = TRUE]], paste0(path, ".source")
  )
  operator <- x[["operator", exact = TRUE]]
  expected_pass <- if (identical(operator, "identical")) {
    is.null(x[["tolerance", exact = TRUE]]) &&
      identical(x[["value", exact = TRUE]], x[["reference", exact = TRUE]])
  } else {
    value <- x[["value", exact = TRUE]]
    reference <- x[["reference", exact = TRUE]]
    tolerance <- x[["tolerance", exact = TRUE]]
    numeric_comparison <- is.numeric(value) && !is.object(value) &&
      is.null(dim(value)) && is.numeric(reference) && !is.object(reference) &&
      is.null(dim(reference)) && is.numeric(tolerance) &&
      !is.object(tolerance) && is.null(dim(tolerance)) &&
      length(value) > 0L && length(reference) %in% c(1L, length(value)) &&
      length(tolerance) %in% c(1L, length(value)) &&
      !anyNA(value) && !anyNA(reference) && !anyNA(tolerance) &&
      all(is.finite(value)) && all(is.finite(reference)) &&
      all(is.finite(tolerance)) && all(tolerance >= 0)
    if (!numeric_comparison) {
      FALSE
    } else {
      switch(
        operator,
        abs_lte = all(abs(value - reference) <= tolerance),
        lte = all(value <= reference + tolerance),
        gte = all(value >= reference - tolerance),
        lt = all(value < reference - tolerance),
        gt = all(value > reference + tolerance)
      )
    }
  }
  .dpprior_schema_require(
    identical(x[["passed", exact = TRUE]], expected_pass),
    "decision_pass", paste0(path, ".passed"),
    "the result recomputed from value/reference/tolerance/operator",
    x[["passed", exact = TRUE]]
  )
  invisible(TRUE)
}


.dpprior_new_check <- function(passed = NULL,
                               value,
                               reference = 0,
                               tolerance = NULL,
                               operator = if (is.null(tolerance)) {
                                 "identical"
                               } else {
                                 "abs_lte"
                               },
                               source) {
  if (is.null(passed)) {
    passed <- if (identical(operator, "identical")) {
      identical(value, reference) && is.null(tolerance)
    } else if (is.numeric(value) && is.numeric(reference) &&
               is.numeric(tolerance)) {
      switch(
        operator,
        abs_lte = all(abs(value - reference) <= tolerance),
        lte = all(value <= reference + tolerance),
        gte = all(value >= reference - tolerance),
        lt = all(value < reference - tolerance),
        gt = all(value > reference + tolerance),
        FALSE
      )
    } else {
      FALSE
    }
  }
  out <- list(
    passed = passed,
    value = value,
    reference = reference,
    tolerance = tolerance,
    operator = operator,
    source = source
  )
  .dpprior_validate_decision_check(out)
  out
}


.dpprior_bind_decision_check <- function(x,
                                         value,
                                         reference,
                                         tolerance,
                                         operator,
                                         path) {
  .dpprior_validate_decision_check(x, path)
  expected <- list(
    value = value, reference = reference, tolerance = tolerance,
    operator = operator
  )
  for (field in names(expected)) {
    .dpprior_schema_require(
      identical(x[[field, exact = TRUE]], expected[[field]]),
      "decision_identity", paste0(path, ".", field),
      "identity with independently recomputed canonical decision evidence",
      x[[field, exact = TRUE]]
    )
  }
  invisible(TRUE)
}


.dpprior_required_checks_pass <- function(x, required, path,
                                          substantive = FALSE) {
  missing <- setdiff(required, names(x))
  .dpprior_schema_require(
    length(missing) == 0L, "required_checks", path,
    paste("required checks", paste(required, collapse = ", ")), missing
  )
  passes <- vapply(x[required], .dpprior_schema_check_pass, logical(1))
  .dpprior_schema_require(
    all(passes), "failed_checks", path,
    "all required checks passing", names(passes)[!passes]
  )
  if (substantive) {
    for (field in required) {
      .dpprior_validate_decision_check(
        x[[field, exact = TRUE]], paste0(path, ".", field)
      )
    }
  }
  invisible(TRUE)
}


.dpprior_validate_stability <- function(x, path = "verification.stability") {
  .dpprior_schema_exact_names(x, .DPPRIOR_STABILITY_FIELDS, path)
  for (field in c("delta", "tolerance")) {
    value <- x[[field, exact = TRUE]]
    .dpprior_schema_require(
      is.numeric(value) && !is.object(value) && is.null(dim(value)) &&
        .dpprior_schema_has_only_attributes(value, "names", "names") &&
        length(value) > 0L && !anyNA(value) && all(is.finite(value)) &&
        !is.null(names(value)) && !anyDuplicated(names(value)) &&
        all(nzchar(names(value))),
      "stability", paste0(path, ".", field),
      "a finite uniquely named numeric vector", value
    )
  }
  .dpprior_schema_require(
    identical(names(x[["delta", exact = TRUE]]),
              names(x[["tolerance", exact = TRUE]])) &&
      all(x[["delta", exact = TRUE]] >= 0) &&
      all(x[["tolerance", exact = TRUE]] >= 0),
    "stability", path,
    "matching non-negative delta and tolerance components", x
  )
  formula <- x[["formula", exact = TRUE]]
  .dpprior_schema_require(
    is.character(formula) && !is.object(formula) && is.null(dim(formula)) &&
      .dpprior_schema_has_only_attributes(formula, "names", "names") &&
      length(formula) == length(x[["delta", exact = TRUE]]) &&
      !anyNA(formula) && all(nzchar(formula)) &&
      identical(names(formula), names(x[["delta", exact = TRUE]])) &&
      all(formula %in% c(
        "absolute_plus_relative_max", "direct_pmf_l1_tolerance",
        "fixed_constraint_tolerance", "fixed_precomputed_tolerance"
      )),
    "stability_formula", paste0(path, ".formula"),
    "a named closed formula for every stability component", formula
  )
  scale_floor <- x[["scale_floor", exact = TRUE]]
  .dpprior_schema_require(
    is.numeric(scale_floor) && !is.object(scale_floor) &&
      is.null(dim(scale_floor)) &&
      .dpprior_schema_has_only_attributes(scale_floor, "names", "names") &&
      length(scale_floor) == length(x[["delta", exact = TRUE]]) &&
      !anyNA(scale_floor) && all(is.finite(scale_floor)) &&
      all(scale_floor >= 0) &&
      identical(names(scale_floor), names(x[["delta", exact = TRUE]])),
    "stability_scale_floor", paste0(path, ".scale_floor"),
    "a named nonnegative scale floor for every component", scale_floor
  )
  .dpprior_schema_validate_scalar_logical(
    x[["passed", exact = TRUE]], paste0(path, ".passed")
  )
  .dpprior_schema_validate_scalar_character(
    x[["source", exact = TRUE]], paste0(path, ".source")
  )
  expected_pass <- all(
    x[["delta", exact = TRUE]] <= x[["tolerance", exact = TRUE]]
  )
  .dpprior_schema_require(
    identical(x[["passed", exact = TRUE]], expected_pass),
    "stability_pass", paste0(path, ".passed"),
    "all recorded deltas within their fixed tolerances",
    x[["passed", exact = TRUE]]
  )
  invisible(TRUE)
}


.dpprior_new_stability <- function(
    delta,
    tolerance,
    formula = setNames(
      rep("fixed_precomputed_tolerance", length(delta)), names(delta)
    ),
    scale_floor = setNames(rep(0, length(delta)), names(delta)),
    passed = NULL,
    source) {
  if (is.null(passed)) {
    passed <- is.numeric(delta) && is.numeric(tolerance) &&
      length(delta) == length(tolerance) &&
      all(is.finite(delta)) && all(is.finite(tolerance)) &&
      all(delta <= tolerance)
  }
  out <- list(
    delta = delta,
    tolerance = tolerance,
    formula = formula,
    scale_floor = scale_floor,
    passed = passed,
    source = source
  )
  .dpprior_validate_stability(out)
  out
}


.dpprior_validate_verification <- function(x, path = "verification") {
  .dpprior_schema_exact_names(x, .DPPRIOR_VERIFICATION_FIELDS, path)
  .dpprior_schema_validate_scalar_character(x$method, paste0(path, ".method"))
  .dpprior_schema_validate_scalar_logical(
    x$performed, paste0(path, ".performed")
  )
  .dpprior_schema_validate_scalar_logical(x$passed, paste0(path, ".passed"))
  .dpprior_schema_validate_scalar_character(
    x$reason, paste0(path, ".reason"), allow_empty = TRUE
  )
  .dpprior_schema_validate_named_list(x$settings, paste0(path, ".settings"))
  .dpprior_schema_validate_plain_record_value(
    x$settings, paste0(path, ".settings")
  )
  .dpprior_validate_snapshot(
    x$selected_snapshot, paste0(path, ".selected_snapshot")
  )
  .dpprior_validate_snapshot(
    x$verifier_snapshot, paste0(path, ".verifier_snapshot")
  )
  if (!is.null(x$stability)) {
    .dpprior_validate_stability(x$stability, paste0(path, ".stability"))
  }
  .dpprior_validate_check_list(x$components, paste0(path, ".components"))
  .dpprior_validate_check_list(x$invariants, paste0(path, ".invariants"))

  if (x$passed) {
    snapshot_pair <- !is.null(x$selected_snapshot) &&
      !is.null(x$verifier_snapshot)
    certificate_only <- is.null(x$selected_snapshot) &&
      is.null(x$verifier_snapshot) && is.null(x$stability) &&
      identical(names(x$components), "infeasibility_certificate")
    .dpprior_schema_require(
      x$performed && (snapshot_pair || certificate_only),
      "verification_pass", path,
      paste(
        "performed verification with selected/verifier snapshots or a",
        "certificate-only analytic infeasibility proof"
      ), x
    )
    if (snapshot_pair) {
      .dpprior_schema_require(
        isTRUE(x$selected_snapshot[["finite", exact = TRUE]]) &&
          isTRUE(x$verifier_snapshot[["finite", exact = TRUE]]),
        "verification_finite", path,
        "finite selected and verifier snapshots for passed verification",
        list(
          selected = x$selected_snapshot[["finite", exact = TRUE]],
          verifier = x$verifier_snapshot[["finite", exact = TRUE]]
        )
      )
    }
    if (!is.null(x$selected_snapshot[["parameters", exact = TRUE]]) ||
        !is.null(x$verifier_snapshot[["parameters", exact = TRUE]])) {
      .dpprior_schema_require(
        identical(
          x$selected_snapshot[["parameters", exact = TRUE]],
          x$verifier_snapshot[["parameters", exact = TRUE]]
        ),
        "verifier_candidate", paste0(path, ".verifier_snapshot.parameters"),
        "identity with selected snapshot parameters",
        x$verifier_snapshot[["parameters", exact = TRUE]]
      )
    }
    if (length(x$components) > 0L) {
      .dpprior_required_checks_pass(
        x$components, names(x$components), paste0(path, ".components")
      )
    }
    if (length(x$invariants) > 0L) {
      .dpprior_required_checks_pass(
        x$invariants, names(x$invariants), paste0(path, ".invariants")
      )
    }
  }
  if (!x$performed) {
    .dpprior_schema_require(
      !x$passed && is.null(x$verifier_snapshot),
      "verification_not_performed", path,
      "passed=FALSE and verifier_snapshot=NULL when not performed", x
    )
  }
  invisible(TRUE)
}


.dpprior_new_verification <- function(method,
                                      performed,
                                      passed,
                                      reason = "",
                                      settings = list(),
                                      selected_snapshot = NULL,
                                      verifier_snapshot = NULL,
                                      stability = NULL,
                                      components = list(),
                                      invariants = list()) {
  out <- list(
    method = method,
    performed = performed,
    passed = passed,
    reason = reason,
    settings = settings,
    selected_snapshot = selected_snapshot,
    verifier_snapshot = verifier_snapshot,
    stability = stability,
    components = components,
    invariants = invariants
  )
  .dpprior_validate_verification(out)
  out
}


# --- Provenance and compatibility -------------------------------------------

.dpprior_validate_input_fit_target_reference <- function(x, J, path) {
  .dpprior_schema_exact_names(
    x, c("schema", "kind", "J", "used", "implied"), path
  )
  .dpprior_validate_schema_record(
    x[["schema", exact = TRUE]], .DPPRIOR_TARGET_SCHEMA_NAME,
    paste0(path, ".schema")
  )
  .dpprior_schema_validate_scalar_character(
    x[["kind", exact = TRUE]], paste0(path, ".kind")
  )
  .dpprior_schema_require(
    x[["kind", exact = TRUE]] %in% c("moments", "cv", "pmf", "interval", "family") &&
      identical(x[["J", exact = TRUE]], J),
    "input_fit_target", path,
    "a supported target kind and J identical to the input-fit reference", x
  )
  for (field in c("used", "implied")) {
    .dpprior_schema_validate_named_list(
      x[[field, exact = TRUE]], paste0(path, ".", field)
    )
    .dpprior_schema_validate_plain_record_value(
      x[[field, exact = TRUE]], paste0(path, ".", field)
    )
  }
  invisible(TRUE)
}


.dpprior_validate_input_fit_snapshot_reference <- function(
    x, J, expected_source, path) {
  .dpprior_schema_exact_names(
    x, c("parameters", "M", "achieved_K", "finite", "source"), path
  )
  .dpprior_validate_parameters(
    x[["parameters", exact = TRUE]], paste0(path, ".parameters")
  )
  .dpprior_schema_require(
    .dpprior_schema_is_count(x[["M", exact = TRUE]], 10L),
    "input_fit_snapshot_order", paste0(path, ".M"),
    "a retained quadrature order of at least 10", x[["M", exact = TRUE]]
  )
  .dpprior_validate_achieved_K(
    x[["achieved_K", exact = TRUE]], J, paste0(path, ".achieved_K")
  )
  .dpprior_schema_validate_scalar_logical(
    x[["finite", exact = TRUE]], paste0(path, ".finite")
  )
  .dpprior_schema_validate_scalar_character(
    x[["source", exact = TRUE]], paste0(path, ".source")
  )
  .dpprior_schema_require(
    x[["finite", exact = TRUE]] &&
      identical(x[["source", exact = TRUE]], expected_source) &&
      identical(
        x[["achieved_K", exact = TRUE]][["M", exact = TRUE]],
        x[["M", exact = TRUE]]
      ),
    "input_fit_snapshot", path,
    paste(
      "finite", expected_source,
      "input-fit evidence with snapshot M identical to achieved_K.M"
    ), x
  )
  invisible(TRUE)
}


.dpprior_validate_input_fit_A2_KL_tolerances <- function(x, path) {
  .dpprior_schema_exact_names(x, c("adequacy", "order"), path)
  adequacy <- x[["adequacy", exact = TRUE]]
  .dpprior_schema_exact_names(
    adequacy,
    c(
      "kl", "l1", "mean_scaled", "variance_scaled",
      "mean_scale_formula", "variance_scale_formula"
    ), paste0(path, ".adequacy")
  )
  adequacy_caps <- c(
    kl = 0.015, l1 = 0.11, mean_scaled = 0.01, variance_scaled = 0.065
  )
  for (field in names(adequacy_caps)) {
    .dpprior_schema_validate_finite_scalar(
      adequacy[[field, exact = TRUE]], paste0(path, ".adequacy.", field),
      lower = 0
    )
    .dpprior_schema_require(
      adequacy[[field, exact = TRUE]] <= adequacy_caps[[field]],
      "input_fit_A2_KL_tolerance", paste0(path, ".adequacy.", field),
      paste("no greater than", format(adequacy_caps[[field]])),
      adequacy[[field, exact = TRUE]]
    )
  }
  .dpprior_schema_require(
    identical(
      adequacy[["mean_scale_formula", exact = TRUE]],
      "max(1,sqrt(target_variance))"
    ) && identical(
      adequacy[["variance_scale_formula", exact = TRUE]],
      "max(1,target_variance)"
    ),
    "input_fit_A2_KL_formula", paste0(path, ".adequacy"),
    "the exact four-gate A2-KL scale formulas", adequacy
  )

  order <- x[["order", exact = TRUE]]
  .dpprior_schema_exact_names(
    order,
    c(
      "pmf_absolute", "pmf_relative", "pmf_l1",
      "direct_moment_absolute", "direct_moment_relative",
      "target_identity_l1"
    ), paste0(path, ".order")
  )
  for (field in names(order)) {
    .dpprior_schema_validate_finite_scalar(
      order[[field, exact = TRUE]], paste0(path, ".order.", field), lower = 0
    )
  }
  .dpprior_schema_require(
    order[["pmf_absolute", exact = TRUE]] <= 1e-10 &&
      order[["pmf_relative", exact = TRUE]] <= 1e-8 &&
      identical(
        order[["pmf_l1", exact = TRUE]],
        order[["pmf_absolute", exact = TRUE]] +
          order[["pmf_relative", exact = TRUE]]
      ) && identical(
        order[["direct_moment_absolute", exact = TRUE]],
        order[["pmf_absolute", exact = TRUE]]
      ) && identical(
        order[["direct_moment_relative", exact = TRUE]],
        order[["pmf_relative", exact = TRUE]]
      ) && identical(order[["target_identity_l1", exact = TRUE]],
                     .TOL_PMF_SUM),
    "input_fit_A2_KL_tolerance", paste0(path, ".order"),
    paste(
      "capped PMF-order and direct-moment tolerances with exact derived",
      "PMF L1 and target-identity gates"
    ), order
  )
  invisible(TRUE)
}


.dpprior_validate_input_fit_A2_KL_decision_evidence <- function(
    x, input_fit, path) {
  .dpprior_schema_exact_names(
    x, c("target_K", "distribution_tolerances"), path
  )
  target_K <- x[["target_K", exact = TRUE]]
  .dpprior_validate_target_v1_impl(target_K)
  target_raw <- unclass(target_K)
  expected_compact_target <- list(
    schema = target_raw[["schema", exact = TRUE]],
    kind = target_raw[["kind", exact = TRUE]],
    J = target_raw[["J", exact = TRUE]],
    used = target_raw[["used", exact = TRUE]],
    implied = target_raw[["implied", exact = TRUE]]
  )
  .dpprior_schema_require(
    identical(input_fit[["target", exact = TRUE]], expected_compact_target) &&
      identical(target_raw[["J", exact = TRUE]],
                input_fit[["J", exact = TRUE]]),
    "input_fit_A2_KL_target_identity", paste0(path, ".target_K"),
    "a full canonical target whose compact identity is retained by input_fit",
    list(full = target_K, compact = input_fit[["target", exact = TRUE]])
  )

  distribution <- x[["distribution_tolerances", exact = TRUE]]
  .dpprior_validate_input_fit_A2_KL_tolerances(
    distribution, paste0(path, ".distribution_tolerances")
  )
  adequacy <- distribution[["adequacy", exact = TRUE]]
  order <- distribution[["order", exact = TRUE]]
  J <- input_fit[["J", exact = TRUE]]
  parameters <- input_fit[["parameters", exact = TRUE]]
  selected_M <- input_fit[["selected_snapshot", exact = TRUE]][[
    "M", exact = TRUE
  ]]
  verifier_M <- input_fit[["verifier_snapshot", exact = TRUE]][[
    "M", exact = TRUE
  ]]
  fresh <- tryCatch(
    .get_K_pmf_support(
      J, parameters[["a", exact = TRUE]], parameters[["b", exact = TRUE]],
      M = selected_M, M_verify = verifier_M,
      abs_tol = order[["pmf_absolute", exact = TRUE]],
      rel_tol = order[["pmf_relative", exact = TRUE]]
    ),
    error = function(error) error
  )
  .dpprior_schema_require(
    is.list(fresh) && !inherits(fresh, "condition") &&
      is.numeric(fresh[["pmf", exact = TRUE]]) &&
      is.numeric(fresh[["verification_pmf", exact = TRUE]]) &&
      length(fresh[["pmf", exact = TRUE]]) == J &&
      length(fresh[["verification_pmf", exact = TRUE]]) == J,
    "input_fit_A2_KL_recomputation", path,
    "successful selected and actual-used-order verifier PMF recomputation",
    fresh
  )
  fresh_pmfs <- list(
    selected_snapshot = unname(as.numeric(fresh[["pmf", exact = TRUE]])),
    verifier_snapshot = unname(as.numeric(
      fresh[["verification_pmf", exact = TRUE]]
    ))
  )
  fresh_moments <- lapply(fresh_pmfs, .dpprior_target_pmf_moments)
  for (snapshot_name in names(fresh_pmfs)) {
    achieved_K <- input_fit[[snapshot_name, exact = TRUE]][[
      "achieved_K", exact = TRUE
    ]]
    recorded_pmf <- achieved_K[["pmf", exact = TRUE]]
    expected_pmf <- fresh_pmfs[[snapshot_name, exact = TRUE]]
    pmf_tolerance <- 64 * .Machine$double.eps * pmax(
      1, abs(recorded_pmf), abs(expected_pmf)
    )
    expected_moments <- fresh_moments[[snapshot_name, exact = TRUE]]
    .dpprior_schema_require(
      !is.null(recorded_pmf) && length(recorded_pmf) == J &&
        all(abs(recorded_pmf - expected_pmf) <= pmf_tolerance) &&
        .dpprior_sensitivity_close_numeric(
          achieved_K[["mean", exact = TRUE]],
          expected_moments[["mean", exact = TRUE]]
        ) && .dpprior_sensitivity_close_numeric(
          achieved_K[["variance", exact = TRUE]],
          expected_moments[["variance", exact = TRUE]]
        ),
      "input_fit_A2_KL_snapshot_truth",
      paste0(path, ".", snapshot_name, ".achieved_K"),
      "the freshly recomputed PMF and its moments at the retained actual order",
      list(recorded = achieved_K, recomputed = list(
        pmf = expected_pmf, moments = expected_moments
      ))
    )
  }

  target_pmf <- .dpprior_result_A2_KL_objective_pmf(
    target_raw, J, paste0(path, ".target_K")
  )
  support <- seq_len(J)
  target_mean <- sum(support * target_pmf)
  target_moments <- c(
    mean = target_mean,
    variance = sum((support - target_mean)^2 * target_pmf)
  )
  adequacy_limits <- c(
    kl = adequacy[["kl", exact = TRUE]],
    l1 = adequacy[["l1", exact = TRUE]],
    mean_scaled = adequacy[["mean_scaled", exact = TRUE]],
    variance_scaled = adequacy[["variance_scaled", exact = TRUE]]
  )
  adequacy_values <- lapply(fresh_pmfs, function(pmf) {
    .dpprior_schema_require(
      all(pmf[target_pmf > 0] > 0),
      "input_fit_A2_KL_support", path,
      "positive achieved mass wherever the objective target has mass", pmf
    )
    moments <- .dpprior_target_pmf_moments(pmf)
    positive <- target_pmf > 0
    c(
      kl = sum(target_pmf[positive] * log(
        target_pmf[positive] / pmf[positive]
      )),
      l1 = sum(abs(target_pmf - pmf)),
      mean_scaled = abs(moments[["mean"]] - target_moments[["mean"]]) /
        max(1, sqrt(target_moments[["variance"]])),
      variance_scaled = abs(
        moments[["variance"]] - target_moments[["variance"]]
      ) / max(1, target_moments[["variance"]])
    )
  })
  .dpprior_schema_require(
    all(vapply(adequacy_values, function(values) {
      all(is.finite(values)) && all(values <= adequacy_limits)
    }, logical(1))),
    "input_fit_A2_KL_adequacy", path,
    "all four A2-KL adequacy gates at selected and verifier orders",
    list(values = adequacy_values, tolerances = adequacy_limits)
  )
  pmf_order_l1 <- sum(abs(
    fresh_pmfs[["selected_snapshot", exact = TRUE]] -
      fresh_pmfs[["verifier_snapshot", exact = TRUE]]
  ))
  moment_order_delta <- abs(
    fresh_moments[["selected_snapshot", exact = TRUE]] -
      fresh_moments[["verifier_snapshot", exact = TRUE]]
  )
  moment_order_tolerance <- order[["direct_moment_absolute", exact = TRUE]] +
    order[["direct_moment_relative", exact = TRUE]] * pmax(
      abs(fresh_moments[["selected_snapshot", exact = TRUE]]),
      abs(fresh_moments[["verifier_snapshot", exact = TRUE]]), 1
    )
  .dpprior_schema_require(
    is.finite(pmf_order_l1) &&
      pmf_order_l1 <= order[["pmf_l1", exact = TRUE]] &&
      all(moment_order_delta <= moment_order_tolerance),
    "input_fit_A2_KL_order_stability", path,
    "PMF L1 and direct moments stable at the retained verifier order",
    list(
      pmf_l1 = pmf_order_l1,
      pmf_tolerance = order[["pmf_l1", exact = TRUE]],
      moment_delta = moment_order_delta,
      moment_tolerance = moment_order_tolerance
    )
  )
  invisible(TRUE)
}


.dpprior_validate_provenance <- function(x, path = "provenance") {
  .dpprior_schema_exact_names(x, .DPPRIOR_PROVENANCE_FIELDS, path)
  for (field in c("requested_method", "selected_method", "parameterization")) {
    .dpprior_schema_validate_scalar_character(
      x[[field]], paste0(path, ".", field)
    )
  }
  .dpprior_schema_validate_scalar_logical(
    x$is_fallback, paste0(path, ".is_fallback")
  )

  .dpprior_schema_exact_names(
    x$approximation, c("active", "opt_in", "kind", "warning_code"),
    paste0(path, ".approximation")
  )
  .dpprior_schema_validate_scalar_logical(
    x$approximation$active, paste0(path, ".approximation.active")
  )
  .dpprior_schema_validate_scalar_logical(
    x$approximation$opt_in, paste0(path, ".approximation.opt_in")
  )
  for (field in c("kind", "warning_code")) {
    .dpprior_schema_validate_nullable_character(
      x$approximation[[field]], paste0(path, ".approximation.", field)
    )
  }
  if (x$approximation$active) {
    .dpprior_schema_require(
      !is.null(x$approximation$kind), "approximation",
      paste0(path, ".approximation.kind"),
      "a named approximation kind when active", NULL
    )
  } else {
    .dpprior_schema_require(
      !x$approximation$opt_in && is.null(x$approximation$kind) &&
        is.null(x$approximation$warning_code),
      "approximation", paste0(path, ".approximation"),
      "no opt-in, kind, or warning code when approximation is inactive",
      x$approximation
    )
  }

  .dpprior_schema_exact_names(
    x$projection, c("applied", "opt_in", "policy", "record"),
    paste0(path, ".projection")
  )
  .dpprior_schema_validate_scalar_logical(
    x$projection$applied, paste0(path, ".projection.applied")
  )
  .dpprior_schema_validate_scalar_logical(
    x$projection$opt_in, paste0(path, ".projection.opt_in")
  )
  .dpprior_schema_validate_nullable_character(
    x$projection$policy, paste0(path, ".projection.policy")
  )
  if (!is.null(x$projection$record)) {
    .dpprior_schema_validate_named_list(
      x$projection$record, paste0(path, ".projection.record")
    )
    .dpprior_schema_validate_plain_record_value(
      x$projection$record, paste0(path, ".projection.record")
    )
  }
  if (x$projection$applied) {
    .dpprior_schema_require(
      x$projection$opt_in && !is.null(x$projection$policy) &&
        !is.null(x$projection$record),
      "projection", paste0(path, ".projection"),
      "opt-in policy and before/after record for applied projection",
      x$projection
    )
  }

  .dpprior_schema_exact_names(
    x$backend,
    c("package", "package_version", "implementation", "source_commit"),
    paste0(path, ".backend")
  )
  for (field in c("package", "package_version", "implementation")) {
    .dpprior_schema_validate_scalar_character(
      x$backend[[field]], paste0(path, ".backend.", field)
    )
  }
  .dpprior_schema_validate_nullable_character(
    x$backend$source_commit, paste0(path, ".backend.source_commit")
  )

  if (!is.null(x$input_fit)) {
    .dpprior_schema_exact_names(
      x$input_fit,
      c(
        "schema", "mode", "method", "J", "status", "usable", "verified",
        "parameters", "target", "decision_evidence", "selected_snapshot",
        "verifier_snapshot"
      ),
      paste0(path, ".input_fit")
    )
    .dpprior_schema_validate_scalar_character(
      x$input_fit$schema, paste0(path, ".input_fit.schema")
    )
    .dpprior_schema_validate_scalar_character(
      x$input_fit$mode, paste0(path, ".input_fit.mode")
    )
    .dpprior_schema_validate_scalar_character(
      x$input_fit$method, paste0(path, ".input_fit.method")
    )
    .dpprior_schema_require(
      .dpprior_schema_is_count(x$input_fit$J, 1L),
      "input_fit_J", paste0(path, ".input_fit.J"),
      "a positive integer J", x$input_fit$J
    )
    .dpprior_schema_validate_scalar_character(
      x$input_fit$status, paste0(path, ".input_fit.status")
    )
    .dpprior_schema_validate_scalar_logical(
      x$input_fit$usable, paste0(path, ".input_fit.usable")
    )
    .dpprior_schema_validate_scalar_logical(
      x$input_fit$verified, paste0(path, ".input_fit.verified")
    )
    .dpprior_validate_parameters(
      x$input_fit$parameters, paste0(path, ".input_fit.parameters")
    )
    .dpprior_validate_input_fit_target_reference(
      x$input_fit$target, x$input_fit$J, paste0(path, ".input_fit.target")
    )
    .dpprior_validate_input_fit_snapshot_reference(
      x$input_fit$selected_snapshot, x$input_fit$J,
      "selected_order",
      paste0(path, ".input_fit.selected_snapshot")
    )
    .dpprior_validate_input_fit_snapshot_reference(
      x$input_fit$verifier_snapshot, x$input_fit$J,
      "independent_verifier",
      paste0(path, ".input_fit.verifier_snapshot")
    )
    .dpprior_schema_require(
      identical(x$input_fit$schema, "dpprior.result/1"),
      "input_fit_schema", paste0(path, ".input_fit.schema"),
      "dpprior.result/1", x$input_fit$schema
    )
    .dpprior_schema_require(
      x$input_fit$mode %in%
        c("a2_moment", "a2_kl", "dual_hard", "dual_soft"),
      "input_fit_mode", paste0(path, ".input_fit.mode"),
      "a supported decision-ready fit mode", x$input_fit$mode
    )
    .dpprior_schema_require(
      x$input_fit$method %in% .DPPRIOR_MODE_METHODS[[x$input_fit$mode]],
      "input_fit_method", paste0(path, ".input_fit.method"),
      "a method authorized for input_fit.mode", x$input_fit$method
    )
    .dpprior_schema_require(
      x$input_fit$status %in% c("converged", "boundary") &&
        x$input_fit$usable && x$input_fit$verified &&
        identical(x$input_fit$selected_snapshot$parameters,
                  x$input_fit$parameters) &&
        identical(x$input_fit$verifier_snapshot$parameters,
                  x$input_fit$parameters),
      "input_fit_status", paste0(path, ".input_fit"),
      paste(
        "converged/boundary with usable=verified=TRUE and exact fixed-",
        "candidate selected/verifier parameter identity"
      ), x$input_fit
    )
    selected_input_M <- x$input_fit$selected_snapshot$M
    verifier_input_M <- x$input_fit$verifier_snapshot$M
    required_input_M <- max(
      2 * as.numeric(selected_input_M), as.numeric(selected_input_M) + 40
    )
    .dpprior_schema_require(
      is.finite(required_input_M) && required_input_M <=
        .QUADRATURE_MAX_NODES &&
        verifier_input_M >= required_input_M &&
        verifier_input_M <= .QUADRATURE_MAX_NODES,
      "input_fit_verifier_order", paste0(path, ".input_fit.verifier_snapshot.M"),
      paste(
        "an actual used verifier order at least",
        "max(2*M_selected, M_selected+40) and no greater than 512"
      ),
      list(
        selected = selected_input_M, required = required_input_M,
        used = verifier_input_M
      )
    )
    expected_input_moments <- lapply(
      c(selected_input_M, verifier_input_M),
      function(M) tryCatch(
        exact_K_moments(
          x$input_fit$J,
          x$input_fit$parameters$a,
          x$input_fit$parameters$b,
          M = M
        ),
        error = function(error) error
      )
    )
    .dpprior_schema_require(
      all(vapply(expected_input_moments, is.list, logical(1))),
      "input_fit_recomputation", paste0(path, ".input_fit"),
      "successful fixed-order recomputation of both retained K snapshots",
      expected_input_moments
    )
    for (index in seq_along(expected_input_moments)) {
      snapshot_name <- c("selected_snapshot", "verifier_snapshot")[[index]]
      achieved_K <- x$input_fit[[snapshot_name]]$achieved_K
      recomputed <- expected_input_moments[[index]]
      .dpprior_schema_require(
        .dpprior_sensitivity_close_numeric(achieved_K$mean, recomputed$mean) &&
          .dpprior_sensitivity_close_numeric(achieved_K$variance, recomputed$var),
        "input_fit_snapshot_truth",
        paste0(path, ".input_fit.", snapshot_name, ".achieved_K"),
        "mean/variance freshly recomputed from input-fit J/a/b/M",
        list(recorded = achieved_K, recomputed = recomputed[c("mean", "var")])
      )
    }
    if (identical(x$input_fit$mode, "a2_kl")) {
      .dpprior_schema_require(
        !is.null(x$input_fit$decision_evidence),
        "input_fit_A2_KL_decision_evidence",
        paste0(path, ".input_fit.decision_evidence"),
        "full canonical target and distribution tolerances for A2-KL", NULL
      )
      .dpprior_validate_input_fit_A2_KL_decision_evidence(
        x$input_fit$decision_evidence, x$input_fit,
        paste0(path, ".input_fit.decision_evidence")
      )
    } else {
      .dpprior_schema_require(
        is.null(x$input_fit$decision_evidence),
        "input_fit_decision_evidence",
        paste0(path, ".input_fit.decision_evidence"),
        "NULL outside A2-KL input-fit references",
        x$input_fit$decision_evidence
      )
      input_target_moments <- x$input_fit$target$implied
      input_selected_K <- x$input_fit$selected_snapshot$achieved_K
      .dpprior_schema_require(
        .dpprior_sensitivity_close_numeric(
          input_selected_K$mean, input_target_moments$mean,
          1e-8 + 1e-8 * max(abs(input_target_moments$mean), 1)
        ) && .dpprior_sensitivity_close_numeric(
          input_selected_K$variance, input_target_moments$variance,
          1e-8 + 1e-8 * max(abs(input_target_moments$variance), 1)
        ),
        "input_fit_target_truth", paste0(path, ".input_fit.target"),
        "decision-ready selected K moments meeting the retained input target",
        list(selected = input_selected_K, target = input_target_moments)
      )
    }
  }

  .dpprior_schema_exact_names(
    x$migration,
    c("source_schema", "adapter", "lossless", "missing_evidence", "warnings"),
    paste0(path, ".migration")
  )
  for (field in c("source_schema", "adapter")) {
    .dpprior_schema_validate_scalar_character(
      x$migration[[field]], paste0(path, ".migration.", field)
    )
  }
  .dpprior_schema_validate_scalar_logical(
    x$migration$lossless, paste0(path, ".migration.lossless")
  )
  .dpprior_schema_validate_character_vector(
    x$migration$missing_evidence,
    paste0(path, ".migration.missing_evidence")
  )
  .dpprior_schema_validate_character_vector(
    x$migration$warnings, paste0(path, ".migration.warnings")
  )
  if (x$migration$lossless) {
    .dpprior_schema_require(
      length(x$migration$missing_evidence) == 0L &&
        length(x$migration$warnings) == 0L,
      "migration_lossless", paste0(path, ".migration"),
      "no missing evidence or warnings when lossless=TRUE", x$migration
    )
  }

  .dpprior_schema_exact_names(
    x$legacy, c("active", "contract", "deprecation_stage"),
    paste0(path, ".legacy")
  )
  .dpprior_schema_validate_scalar_logical(
    x$legacy$active, paste0(path, ".legacy.active")
  )
  for (field in c("contract", "deprecation_stage")) {
    .dpprior_schema_validate_nullable_character(
      x$legacy[[field]], paste0(path, ".legacy.", field)
    )
  }
  if (x$legacy$active) {
    .dpprior_schema_require(
      !is.null(x$legacy$contract), "legacy", paste0(path, ".legacy.contract"),
      "a legacy contract name when active", NULL
    )
    .dpprior_schema_require(
      !is.null(x$legacy$deprecation_stage), "legacy",
      paste0(path, ".legacy.deprecation_stage"),
      "a deprecation stage when legacy is active", NULL
    )
  } else {
    .dpprior_schema_require(
      is.null(x$legacy$contract) && is.null(x$legacy$deprecation_stage),
      "legacy", paste0(path, ".legacy"),
      "no contract or deprecation stage when legacy is inactive", x$legacy
    )
  }
  invisible(TRUE)
}


.dpprior_new_provenance <- function(requested_method,
                                    selected_method,
                                    is_fallback,
                                    approximation,
                                    projection,
                                    parameterization,
                                    backend,
                                    input_fit = NULL,
                                    migration,
                                    legacy) {
  out <- list(
    requested_method = requested_method,
    selected_method = selected_method,
    is_fallback = is_fallback,
    approximation = approximation,
    projection = projection,
    parameterization = parameterization,
    backend = backend,
    input_fit = input_fit,
    migration = migration,
    legacy = legacy
  )
  .dpprior_validate_provenance(out)
  out
}


.dpprior_compatibility_quarantine_boundary <- function() {
  list(
    authority = "non_authoritative",
    lossy = TRUE,
    consumer_policy = "ignored_by_scientific_and_decision_consumers"
  )
}


.dpprior_native_compatibility_specs <- function() {
  boundary <- .dpprior_compatibility_quarantine_boundary()
  list(
    a1_v0 = list(
      mode = "a1_proxy", method = "A1",
      implementation = "R/10_a1_mapping.R:DPprior_a1",
      view_names = c(
        "a1_v0", "mapping_verified", "at_target_boundary",
        "boundary_reason", "scaling", "cJ", "converged", "iterations",
        "fit", "diagnostics", "trace"
      ),
      deprecation = c(list(code = "a1_flat_v0_view_quarantined"), boundary)
    ),
    legacy_v2 = list(
      mode = "a2_moment", method = NULL,
      implementation = "R/11_a2_newton.R:DPprior_a2_newton",
      view_names = c("converged", "legacy_v2"),
      deprecation = c(list(
        code = "a2_moment_legacy_view",
        first_deprecated_version = "2.0.0",
        removal_floor = "not_scheduled"
      ), boundary)
    ),
    a2_kl_v0 = list(
      mode = "a2_kl", method = "A2-KL",
      implementation = "R/12_a2_kl.R:DPprior_a2_kl",
      view_names = c(
        "a2_kl_v0", "converged", "iterations", "termination", "fit",
        "diagnostics"
      ),
      deprecation = c(list(code = "a2_kl_flat_v0_view_quarantined"),
                      boundary)
    ),
    legacy_dual_v2 = list(
      mode = "dual_legacy", method = "dual-anchor",
      implementation = "R/15_dual_anchor.R:DPprior_dual",
      view_names = c(
        "legacy_dual_v2", "converged", "iterations", "fit", "attempts",
        "dual_anchor"
      ),
      deprecation = c(list(
        code = "legacy_dual_flat_view_quarantined",
        first_deprecated_version = "2.0.0",
        removal_floor = "not_scheduled"
      ), boundary)
    )
  )
}


.dpprior_new_compatibility <- function(top_level_aliases = character(),
                                       views = list(),
                                       deprecations = list()) {
  # A v1.1 source summary is retained for auditability only.  Its digest and
  # descriptive leaves cannot be re-authenticated after the raw source bytes
  # have deliberately been discarded, so make the trust boundary explicit at
  # construction time.  Canonical result/target/verification fields never read
  # these leaves; an optional fixed-candidate audit is separately recomputed.
  source <- views[["source", exact = TRUE]]
  quarantined_source <- typeof(source) == "list" && is.list(source) &&
    !is.object(source) && is.null(dim(source)) &&
    isTRUE(source[["public_candidate_quarantined", exact = TRUE]]) &&
    identical(source[["required_action", exact = TRUE]],
              "refit_with_current_API")
  if (quarantined_source) {
    boundary <- .dpprior_compatibility_quarantine_boundary()
    boundary_names <- names(boundary)
    present <- boundary_names %in% names(source)
    .dpprior_schema_require(
      !any(present) || (all(present) &&
        identical(source[boundary_names], boundary)),
      "compatibility_quarantine_boundary", "compatibility.views.source",
      paste(
        "either no pre-normalized quarantine boundary or the exact",
        "non-authoritative lossy consumer boundary"
      ), source[intersect(boundary_names, names(source))]
    )
    if (!any(present)) {
      views[["source"]] <- c(source, boundary)
    }
    legacy_schema <- deprecations[["legacy_schema", exact = TRUE]]
    .dpprior_schema_require(
      typeof(legacy_schema) == "list" && is.list(legacy_schema) &&
        !is.object(legacy_schema) && is.null(dim(legacy_schema)),
      "compatibility_quarantine_boundary",
      "compatibility.deprecations.legacy_schema",
      "one ordinary legacy-schema deprecation record", legacy_schema
    )
    deprecation_present <- boundary_names %in% names(legacy_schema)
    .dpprior_schema_require(
      !any(deprecation_present) || (all(deprecation_present) &&
        identical(legacy_schema[boundary_names], boundary)),
      "compatibility_quarantine_boundary",
      "compatibility.deprecations.legacy_schema",
      paste(
        "either no pre-normalized quarantine boundary or the exact",
        "non-authoritative lossy consumer boundary"
      ), legacy_schema[intersect(boundary_names, names(legacy_schema))]
    )
    if (!any(deprecation_present)) {
      deprecations[["legacy_schema"]] <- c(legacy_schema, boundary)
    }
  }

  # Migrated target views mix canonical mirrors with additional historical
  # descriptors.  The canonical target fields above the compatibility layer
  # are authoritative; the retained legacy view itself is not.
  legacy_target_view <- all(c("legacy_target", "legacy_target_fields") %in%
                              names(views))
  if (legacy_target_view) {
    boundary <- .dpprior_compatibility_quarantine_boundary()
    expected <- c(list(
      code = "legacy_target_view_quarantined",
      first_deprecated_version = "2.0.0",
      removal_floor = "not_scheduled"
    ), boundary)
    recorded <- deprecations[["legacy_target_view", exact = TRUE]]
    .dpprior_schema_require(
      is.null(recorded) || identical(recorded, expected),
      "compatibility_quarantine_boundary",
      "compatibility.deprecations.legacy_target_view",
      "the exact non-authoritative lossy legacy-target-view boundary",
      recorded
    )
    deprecations[["legacy_target_view"]] <- expected
  }

  # Phase-8 native producers carry large historical records solely so old
  # callers can inspect them.  Normalize the frozen producers' marker records
  # at the constructor boundary (A1 predates the common spelling), while the
  # result validator below rejects any later marker or route mutation.
  native_specs <- .dpprior_native_compatibility_specs()
  native_markers <- intersect(names(native_specs), names(views))
  for (marker in native_markers) {
    spec <- native_specs[[marker, exact = TRUE]]
    view <- views[[marker, exact = TRUE]]
    .dpprior_schema_require(
      typeof(view) == "list" && is.list(view) && !is.object(view) &&
        is.null(dim(view)),
      "compatibility_quarantine_boundary",
      paste0("compatibility.views.", marker),
      "one ordinary non-authoritative legacy view", view
    )
    boundary <- .dpprior_compatibility_quarantine_boundary()
    boundary_names <- names(boundary)
    recorded_boundary <- view[intersect(boundary_names, names(view))]
    .dpprior_schema_require(
      length(recorded_boundary) == 0L ||
        identical(recorded_boundary, boundary[names(recorded_boundary)]),
      "compatibility_quarantine_boundary",
      paste0("compatibility.views.", marker),
      "only the exact non-authoritative lossy consumer boundary",
      recorded_boundary
    )
    for (field in boundary_names) {
      view[[field]] <- boundary[[field, exact = TRUE]]
    }
    views[[marker]] <- view

    recorded_deprecation <- deprecations[[marker, exact = TRUE]]
    legacy_a1_deprecation <- list(
      code = "a1_flat_v0_view_quarantined",
      authority = "non_authoritative",
      consumer_policy = "canonical_fields_only"
    )
    .dpprior_schema_require(
      identical(recorded_deprecation, spec$deprecation) ||
        (identical(marker, "a1_v0") &&
           identical(recorded_deprecation, legacy_a1_deprecation)),
      "compatibility_quarantine_boundary",
      paste0("compatibility.deprecations.", marker),
      "the exact route-specific native compatibility deprecation",
      recorded_deprecation
    )
    deprecations[[marker]] <- spec$deprecation
  }

  out <- list(
    top_level_aliases = top_level_aliases,
    views = views,
    deprecations = deprecations
  )
  .dpprior_validate_compatibility(out)
  out
}


.dpprior_validate_compatibility <- function(x, path = "compatibility") {
  .dpprior_schema_exact_names(x, .DPPRIOR_COMPATIBILITY_FIELDS, path)
  .dpprior_schema_validate_character_vector(
    x$top_level_aliases, paste0(path, ".top_level_aliases"), named = TRUE
  )
  .dpprior_schema_validate_named_list(x$views, paste0(path, ".views"))
  .dpprior_schema_validate_named_list(
    x$deprecations, paste0(path, ".deprecations")
  )
  .dpprior_schema_validate_plain_record_value(
    x[["views", exact = TRUE]], paste0(path, ".views")
  )
  .dpprior_schema_validate_plain_record_value(
    x[["deprecations", exact = TRUE]], paste0(path, ".deprecations")
  )
  invisible(TRUE)
}


.dpprior_validate_native_compatibility_boundary <- function(
    raw, allow_constructor_pending = FALSE) {
  compatibility <- raw[["compatibility", exact = TRUE]]
  views <- compatibility[["views", exact = TRUE]]
  deprecations <- compatibility[["deprecations", exact = TRUE]]
  aliases <- compatibility[["top_level_aliases", exact = TRUE]]
  specs <- .dpprior_native_compatibility_specs()
  markers <- names(specs)
  recorded_markers <- intersect(
    markers, unique(c(names(views), names(deprecations)))
  )
  implementation <- raw[["provenance", exact = TRUE]][[
    "backend", exact = TRUE
  ]][["implementation", exact = TRUE]]
  implementation_markers <- markers[vapply(
    specs,
    function(spec) identical(implementation, spec$implementation),
    logical(1)
  )]
  pending_empty <- isTRUE(allow_constructor_pending) &&
    length(recorded_markers) == 0L &&
    length(implementation_markers) == 1L &&
    length(aliases) == 0L && length(views) == 0L &&
    length(deprecations) == 0L
  if (pending_empty) {
    # R10/R12 build the canonical spine first and immediately pass it through
    # .dpprior_append_compatibility_v2().  Only that private constructor call
    # may observe the empty intermediate; every public validation below uses
    # allow_constructor_pending=FALSE and therefore rejects marker removal.
    return(invisible(TRUE))
  }
  if (length(recorded_markers) == 0L &&
      length(implementation_markers) == 0L) {
    return(invisible(TRUE))
  }

  route_ok <- length(recorded_markers) == 1L &&
    length(implementation_markers) == 1L &&
    identical(recorded_markers, implementation_markers)
  .dpprior_schema_require(
    route_ok,
    "native_compatibility_route", "result.compatibility",
    paste(
      "exactly one producer marker matching the canonical backend",
      "implementation; native markers cannot be grafted onto another mode"
    ),
    list(
      recorded_markers = recorded_markers,
      implementation = implementation,
      implementation_markers = implementation_markers
    )
  )
  marker <- recorded_markers[[1L]]
  spec <- specs[[marker, exact = TRUE]]
  .dpprior_schema_require(
    identical(raw[["object_type", exact = TRUE]], "fit") &&
      identical(raw[["mode", exact = TRUE]], spec$mode) &&
      (is.null(spec$method) ||
         identical(raw[["method", exact = TRUE]], spec$method)) &&
      identical(names(views), spec$view_names) &&
      identical(names(deprecations), marker),
    "native_compatibility_route", "result.compatibility",
    paste(
      "the exact mode, method, ordered view vocabulary, and sole",
      "route-specific deprecation emitted by the frozen native producer"
    ),
    list(
      object_type = raw[["object_type", exact = TRUE]],
      mode = raw[["mode", exact = TRUE]],
      method = raw[["method", exact = TRUE]],
      views = names(views), deprecations = names(deprecations)
    )
  )
  boundary <- .dpprior_compatibility_quarantine_boundary()
  marker_view <- views[[marker, exact = TRUE]]
  .dpprior_schema_require(
    identical(marker_view[names(boundary)], boundary) &&
      identical(deprecations[[marker, exact = TRUE]], spec$deprecation),
    "native_compatibility_boundary", "result.compatibility",
    paste(
      "the exact non-authoritative lossy consumer boundary on both the",
      "legacy view and its route-specific deprecation"
    ),
    list(
      view = marker_view[names(boundary)],
      deprecation = deprecations[[marker, exact = TRUE]]
    )
  )
  if (identical(marker, "legacy_dual_v2")) {
    .dpprior_schema_require(
      identical(
        views[["dual_anchor", exact = TRUE]][names(boundary)], boundary
      ),
      "native_compatibility_boundary",
      "result.compatibility.views.dual_anchor",
      paste(
        "the exact non-authoritative lossy consumer boundary on the",
        "second retained legacy dual view"
      ),
      views[["dual_anchor", exact = TRUE]][names(boundary)]
    )
  }

  expected_converged <- raw[["status", exact = TRUE]] %in%
    c("converged", "boundary") &&
    isTRUE(raw[["usable", exact = TRUE]]) &&
    isTRUE(raw[["verified", exact = TRUE]])
  .dpprior_schema_require(
    identical(views[["converged", exact = TRUE]], expected_converged),
    "native_compatibility_decision_mirror",
    "result.compatibility.views.converged",
    "the canonical status/usable/verified convergence decision",
    views[["converged", exact = TRUE]]
  )
  allowed_aliases <- if (identical(marker, "a1_v0")) {
    value <- c(
      mapping_verified = "compatibility.views.mapping_verified",
      at_target_boundary = "compatibility.views.at_target_boundary",
      scaling = "compatibility.views.scaling",
      cJ = "compatibility.views.cJ",
      var_K_used = "target.K.used.var_K",
      projection = "proxy.projection", caveats = "proxy.caveats",
      attempts = "computation.attempts",
      mapping_verification = "proxy.mapping_verification",
      converged = "compatibility.views.converged",
      iterations = "compatibility.views.iterations"
    )
    if (!is.null(raw[["parameters", exact = TRUE]])) {
      value <- c(a = "parameters.a", b = "parameters.b", value)
    }
    value
  } else if (identical(marker, "legacy_v2")) {
    if (is.null(raw[["parameters", exact = TRUE]])) {
      c(converged = "compatibility.views.converged")
    } else {
      c(
        a = "parameters.a", b = "parameters.b",
        converged = "compatibility.views.converged",
        iterations = "compatibility.views.legacy_v2.iterations",
        termination = "compatibility.views.legacy_v2.termination",
        fit = "compatibility.views.legacy_v2.selected_fit",
        diagnostics = "compatibility.views.legacy_v2.solver_diagnostics",
        trace = "computation.trace", attempts = "computation.attempts"
      )
    }
  } else if (identical(marker, "legacy_dual_v2")) {
    c(
      a = "parameters.a", b = "parameters.b",
      converged = "compatibility.views.converged",
      iterations = "compatibility.views.iterations",
      fit = "compatibility.views.fit",
      attempts = "compatibility.views.attempts",
      dual_anchor = "compatibility.views.dual_anchor"
    )
  } else {
    c(
      a = "parameters.a", b = "parameters.b",
      attempts = "computation.attempts",
      converged = "compatibility.views.converged",
      iterations = "compatibility.views.iterations",
      termination = "compatibility.views.termination",
      fit = "compatibility.views.fit",
      diagnostics = "compatibility.views.diagnostics",
      trace = "computation.trace"
    )
  }
  alias_subset <- allowed_aliases[
    names(allowed_aliases) %in% names(aliases)
  ]
  alias_ok <- length(aliases) == 0L || identical(aliases, alias_subset)
  .dpprior_schema_require(
    alias_ok,
    "native_compatibility_aliases",
    "result.compatibility.top_level_aliases",
    paste(
      "an ordered subset of the frozen route's derived aliases; no legacy",
      "view may mint a new scientific or decision alias"
    ), aliases
  )
  invisible(TRUE)
}


# --- Canonical targets -------------------------------------------------------

.dpprior_validate_weight_target_v1 <- function(x, collect = FALSE) {
  validate <- function(object) {
    .dpprior_schema_require(
      typeof(object) == "list" && is.list(object),
      "type", "weight_target", "an ordinary list", typeof(object)
    )
    raw <- unclass(object)
    .dpprior_schema_require(
      identical(
        class(object), c("dpprior_weight_target", "dpprior_target", "list")
      ),
      "class", "weight_target",
      "the exact canonical weight-target class vector", class(object)
    )
    .dpprior_schema_exact_names(
      raw, .DPPRIOR_WEIGHT_TARGET_FIELDS, "weight_target"
    )
    .dpprior_validate_schema_record(
      raw$schema, .DPPRIOR_WEIGHT_TARGET_SCHEMA_NAME, "weight_target.schema"
    )
    .dpprior_schema_require(
      identical(raw$kind, "weight"), "kind", "weight_target.kind",
      "weight", raw$kind
    )
    for (field in c("request", "normalized", "used", "certification",
                    "provenance")) {
      .dpprior_schema_validate_named_list(
        raw[[field, exact = TRUE]], paste0("weight_target.", field)
      )
      .dpprior_schema_validate_plain_record_value(
        raw[[field, exact = TRUE]], paste0("weight_target.", field)
      )
    }
    for (field in c("metric", "relation", "operator", "estimand", "units")) {
      .dpprior_schema_validate_scalar_character(
        raw[[field]], paste0("weight_target.", field)
      )
    }
    for (field in c("value", "threshold", "probability")) {
      if (!is.null(raw[[field]])) {
        .dpprior_schema_validate_finite_scalar(
          raw[[field]], paste0("weight_target.", field)
        )
      }
    }
    .dpprior_schema_require(
      raw$metric %in% c(
        "wsb_mean", "wsb_tail", "wsb_quantile", "wmax_tail_upper"
      ),
      "weight_metric", "weight_target.metric",
      "an approved Phase 8 weight metric", raw$metric
    )
    .dpprior_schema_require(
      raw$relation %in% c("at_most", "at_least", "target"),
      "weight_relation", "weight_target.relation",
      "at_most, at_least, or target", raw$relation
    )
    expected_operator <- switch(
      raw$relation, at_most = "<=", at_least = ">=", target = "target"
    )
    .dpprior_schema_require(
      identical(raw$operator, expected_operator),
      "weight_operator", "weight_target.operator",
      sprintf("%s for relation %s", expected_operator, raw$relation),
      raw$operator
    )
    .dpprior_schema_validate_finite_scalar(
      raw$value, "weight_target.value", lower = 0, upper = 1
    )
    needs_threshold <- raw$metric %in% c("wsb_tail", "wmax_tail_upper")
    needs_probability <- identical(raw$metric, "wsb_quantile")
    if (needs_threshold) {
      .dpprior_schema_validate_finite_scalar(
        raw$threshold, "weight_target.threshold", lower = 0, upper = 1
      )
    } else {
      .dpprior_schema_require(
        is.null(raw$threshold), "irrelevant_component",
        "weight_target.threshold", "NULL for this metric", raw$threshold
      )
    }
    if (needs_probability) {
      .dpprior_schema_validate_finite_scalar(
        raw$probability, "weight_target.probability", lower = 0, upper = 1
      )
    } else {
      .dpprior_schema_require(
        is.null(raw$probability), "irrelevant_component",
        "weight_target.probability", "NULL for this metric", raw$probability
      )
    }
    estimand <- switch(
      raw$metric,
      wsb_tail = "P(W_SB > threshold)",
      wsb_mean = "E(W_SB)",
      wsb_quantile = "Q_probability(W_SB)",
      wmax_tail_upper = "certified upper bound for P(W_max > threshold)"
    )
    .dpprior_schema_require(
      identical(raw$estimand, estimand), "weight_estimand",
      "weight_target.estimand", estimand, raw$estimand
    )
    .dpprior_schema_require(
      identical(raw$units, "probability"), "weight_units",
      "weight_target.units", "probability", raw$units
    )
    authority_fields <- c(
      "metric", "relation", "value", "threshold", "probability"
    )
    for (record_name in c("normalized", "used")) {
      record <- raw[[record_name, exact = TRUE]]
      .dpprior_schema_exact_names(
        record, authority_fields, paste0("weight_target.", record_name)
      )
      for (field in authority_fields) {
        .dpprior_schema_require(
          identical(record[[field, exact = TRUE]], raw[[field, exact = TRUE]]),
          "weight_authority",
          paste0("weight_target.", record_name, ".", field),
          paste0("identity with weight_target.", field),
          record[[field, exact = TRUE]]
        )
      }
    }
    if (!identical(raw[["request", exact = TRUE]],
                   raw[["normalized", exact = TRUE]])) {
      transformation <- raw[["provenance", exact = TRUE]][[
        "transformation", exact = TRUE
      ]]
      .dpprior_schema_require(
        !is.null(transformation), "weight_transformation",
        "weight_target.provenance.transformation",
        "explicit request-to-normalized transformation evidence", NULL
      )
      .dpprior_schema_exact_names(
        transformation,
        c("rule", "opt_in", "before", "after", "evidence"),
        "weight_target.provenance.transformation"
      )
      request <- raw[["request", exact = TRUE]]
      normalized <- raw[["normalized", exact = TRUE]]
      evidence <- transformation[["evidence", exact = TRUE]]
      .dpprior_schema_exact_names(
        evidence,
        c("mode", "value_field", "relation_from", "relation_to",
          "probability_field"),
        "weight_target.provenance.transformation.evidence"
      )
      for (field in c(
        "mode", "value_field", "relation_from", "relation_to",
        "probability_field"
      )) {
        .dpprior_schema_validate_scalar_character(
          evidence[[field, exact = TRUE]],
          paste0("weight_target.provenance.transformation.evidence.", field)
        )
      }
      hard_route <- "bound" %in% names(request) &&
        !("value" %in% names(request))
      soft_route <- "value" %in% names(request) &&
        !("bound" %in% names(request))
      expected_rule <- if (hard_route) {
        "canonicalize_hard_weight_target"
      } else if (soft_route) {
        "canonicalize_soft_weight_target"
      } else {
        NA_character_
      }
      expected_mode <- if (hard_route) "hard" else if (soft_route) "soft" else
        NA_character_
      expected_value_field <- if (hard_route) "bound" else if (soft_route) {
        "value"
      } else {
        NA_character_
      }
      expected_request_fields <- if (hard_route) {
        c("metric", "relation", "bound", "threshold", "probability")
      } else if (soft_route) {
        c("metric", "relation", "value", "threshold", "probability")
      } else {
        character()
      }
      .dpprior_schema_require(
        length(expected_request_fields) > 0L &&
          identical(names(request), expected_request_fields),
        "weight_request_fields", "weight_target.request",
        paste(
          "the exact mode-specific metric/relation/value/threshold/",
          "probability request fields"
        ), names(request)
      )
      relation_map <- c(
        "<=" = "at_most", "at_most" = "at_most",
        ">=" = "at_least", "at_least" = "at_least",
        "target" = "target"
      )
      request_relation <- request[["relation", exact = TRUE]]
      expected_relation <- if (
        .dpprior_schema_is_scalar_character(request_relation) &&
          request_relation %in% names(relation_map)
      ) {
        unname(relation_map[[request_relation]])
      } else {
        NA_character_
      }
      expected_probability_field <- if (
        is.null(request[["probability", exact = TRUE]])
      ) {
        "none"
      } else {
        "probability"
      }
      request_probability <- switch(
        expected_probability_field,
        probability = request[["probability", exact = TRUE]],
        none = NULL,
        NULL
      )
      request_threshold <- request[["threshold", exact = TRUE]]
      request_value <- if (hard_route) {
        request[["bound", exact = TRUE]]
      } else if (soft_route) {
        request[["value", exact = TRUE]]
      } else {
        NULL
      }
      .dpprior_schema_require(
        transformation[["rule", exact = TRUE]] %in% c(
          "canonicalize_hard_weight_target",
          "canonicalize_soft_weight_target"
        ) &&
          identical(transformation[["rule", exact = TRUE]], expected_rule) &&
          identical(transformation[["opt_in", exact = TRUE]], FALSE) &&
          identical(transformation[["before", exact = TRUE]],
                    raw[["request", exact = TRUE]]) &&
          identical(transformation[["after", exact = TRUE]],
                    normalized) &&
          identical(evidence[["mode", exact = TRUE]], expected_mode) &&
          identical(evidence[["value_field", exact = TRUE]],
                    expected_value_field) &&
          identical(evidence[["relation_from", exact = TRUE]],
                    request_relation) &&
          identical(evidence[["relation_to", exact = TRUE]],
                    expected_relation) &&
          identical(evidence[["probability_field", exact = TRUE]],
                    expected_probability_field) &&
          identical(request[["metric", exact = TRUE]],
                    normalized[["metric", exact = TRUE]]) &&
          identical(request_value, normalized[["value", exact = TRUE]]) &&
          identical(expected_relation,
                    normalized[["relation", exact = TRUE]]) &&
          identical(request_threshold,
                    normalized[["threshold", exact = TRUE]]) &&
          identical(request_probability,
                    normalized[["probability", exact = TRUE]]) &&
          (!hard_route || expected_relation %in% c("at_most", "at_least")) &&
          (!soft_route || expected_relation %in%
             c("target", "at_most", "at_least")),
        "weight_transformation_identity",
        "weight_target.provenance.transformation",
        paste(
          "an exact hard-bound or soft-value canonicalization with frozen",
          "relation mapping and explicit probability/prob alias evidence"
        ),
        transformation
      )
      .dpprior_schema_validate_plain_record_value(
        transformation, "weight_target.provenance.transformation"
      )
    } else {
      .dpprior_schema_exact_names(
        raw[["request", exact = TRUE]], authority_fields,
        "weight_target.request"
      )
      .dpprior_schema_require(
        is.null(raw[["provenance", exact = TRUE]][[
          "transformation", exact = TRUE
        ]]),
        "weight_transformation_identity",
        "weight_target.provenance.transformation",
        "NULL when request equals normalized", raw$provenance
      )
    }
    if (!identical(raw[["normalized", exact = TRUE]],
                   raw[["used", exact = TRUE]])) {
      selection <- raw[["provenance", exact = TRUE]][[
        "selection", exact = TRUE
      ]]
      .dpprior_schema_require(
        !is.null(selection), "weight_selection",
        "weight_target.provenance.selection",
        "explicit normalized-to-used selection evidence", NULL
      )
      .dpprior_schema_exact_names(
        selection, c("rule", "opt_in", "before", "after", "evidence"),
        "weight_target.provenance.selection"
      )
      .dpprior_schema_require(
        identical(selection[["rule", exact = TRUE]],
                  "explicit_weight_target_selection") &&
          identical(selection[["opt_in", exact = TRUE]], TRUE) &&
          identical(selection[["before", exact = TRUE]],
                    raw[["normalized", exact = TRUE]]) &&
          identical(selection[["after", exact = TRUE]],
                    raw[["used", exact = TRUE]]),
        "weight_selection_identity", "weight_target.provenance.selection",
        "the exact approved normalized-before/used-after opt-in selection",
        selection
      )
      .dpprior_schema_validate_plain_record_value(
        selection, "weight_target.provenance.selection"
      )
    } else {
      .dpprior_schema_require(
        is.null(raw[["provenance", exact = TRUE]][[
          "selection", exact = TRUE
        ]]),
        "weight_selection_identity", "weight_target.provenance.selection",
        "NULL when normalized equals used", raw$provenance
      )
    }
    if (identical(raw$metric, "wmax_tail_upper")) {
      .dpprior_schema_exact_names(
        raw[["certification", exact = TRUE]],
        c("kind", "certified", "passed", "method", "source"),
        "weight_target.certification"
      )
      certification_kind <- raw[["certification", exact = TRUE]][[
        "kind", exact = TRUE
      ]]
      .dpprior_schema_require(
        identical(certification_kind, "upper_bound") &&
          isTRUE(raw[["certification", exact = TRUE]][[
            "certified", exact = TRUE
          ]]) &&
          isTRUE(raw[["certification", exact = TRUE]][[
            "passed", exact = TRUE
          ]]) &&
          identical(raw[["certification", exact = TRUE]][[
            "method", exact = TRUE
          ]], "certified_size_biased_mass_upper_bound") &&
          identical(raw[["certification", exact = TRUE]][[
            "source", exact = TRUE
          ]], "wmax_tail_bounds"),
        "weight_certification", "weight_target.certification",
        paste(
          "a passing certified size-biased-mass upper bound from",
          "wmax_tail_bounds"
        ), raw$certification
      )
      for (field in c("method", "source")) {
        .dpprior_schema_validate_scalar_character(
          raw[["certification", exact = TRUE]][[field, exact = TRUE]],
          paste0("weight_target.certification.", field)
        )
      }
    } else {
      .dpprior_schema_require(
        identical(raw[["certification", exact = TRUE]], list()),
        "weight_certification", "weight_target.certification",
        "an exact empty record for metrics without W_max certification",
        raw[["certification", exact = TRUE]]
      )
    }
    invisible(object)
  }
  .dpprior_schema_collect_validation(validate, x, collect)
}


.dpprior_new_weight_target <- function(request,
                                       normalized,
                                       used,
                                       metric,
                                       relation,
                                       operator,
                                       value = NULL,
                                       threshold = NULL,
                                       probability = NULL,
                                       estimand,
                                       units,
                                       certification = list(),
                                       provenance = list()) {
  out <- list(
    schema = .dpprior_schema("weight-target"),
    kind = "weight",
    request = request,
    normalized = normalized,
    used = used,
    metric = metric,
    relation = relation,
    operator = operator,
    value = value,
    threshold = threshold,
    probability = probability,
    estimand = estimand,
    units = units,
    certification = certification,
    provenance = provenance
  )
  class(out) <- c("dpprior_weight_target", "dpprior_target", "list")
  .dpprior_validate_weight_target_v1(out)
  out
}


.DPPRIOR_TARGET_DERIVATION_FIELDS <- c(
  "request_to_normalized", "normalized_to_used"
)

.DPPRIOR_TARGET_DERIVATION_ENTRY_FIELDS <- c(
  "rule", "outcome", "opt_in", "before", "after", "evidence"
)

.DPPRIOR_TARGET_NORMALIZATION_RULES <- c(
  "canonicalize_direct_moments", "canonicalize_confidence_target",
  "canonicalize_cv_target", "validate_strict_pmf",
  "drop_structural_k0_zero", "canonicalize_interval_request",
  "canonicalize_scaled_chisq_family_request"
)

.DPPRIOR_TARGET_USE_RULES <- c(
  "derive_variance_from_confidence_vif", "derive_variance_from_cv",
  "construct_maxent_hard_bounds_pmf",
  "construct_maxent_equal_tail_pmf",
  "construct_maxent_central_mass_pmf",
  "construct_scaled_chisq_conditioned_pmf",
  "project_a1_variance_to_nearest_interior"
)


.dpprior_target_pmf_moments <- function(pmf) {
  support <- seq_along(pmf)
  mean_value <- sum(support * pmf)
  c(
    mean = mean_value,
    variance = sum((support - mean_value)^2 * pmf)
  )
}


.dpprior_target_interval_masses <- function(pmf, interval) {
  support <- seq_along(pmf)
  lower <- interval[["lower", exact = TRUE]]
  upper <- interval[["upper", exact = TRUE]]
  c(
    left_mass = sum(pmf[support < lower]),
    inside_mass = sum(pmf[support >= lower & support <= upper]),
    right_mass = sum(pmf[support > upper])
  )
}


.dpprior_expected_target_stability <- function(raw) {
  selected <- raw[["verification", exact = TRUE]][[
    "selected_snapshot", exact = TRUE
  ]][["achieved", exact = TRUE]]
  verifier <- raw[["verification", exact = TRUE]][[
    "verifier_snapshot", exact = TRUE
  ]][["achieved", exact = TRUE]]
  selected_pmf <- selected[["pmf", exact = TRUE]]
  verifier_pmf <- verifier[["pmf", exact = TRUE]]
  selected_moments <- .dpprior_target_pmf_moments(selected_pmf)
  verifier_moments <- .dpprior_target_pmf_moments(verifier_pmf)
  delta <- c(
    K.mean = abs(selected_moments[["mean"]] - verifier_moments[["mean"]]),
    K.variance = abs(
      selected_moments[["variance"]] - verifier_moments[["variance"]]
    ),
    pmf.l1 = sum(abs(selected_pmf - verifier_pmf))
  )
  truth <- raw[["tolerances", exact = TRUE]]
  tolerance <- c(
    K.mean = truth[["moment_relative", exact = TRUE]] * max(
      truth[["moment_scale_floor", exact = TRUE]],
      abs(selected_moments[["mean"]]), abs(verifier_moments[["mean"]])
    ),
    K.variance = truth[["moment_relative", exact = TRUE]] * max(
      truth[["moment_scale_floor", exact = TRUE]],
      abs(selected_moments[["variance"]]),
      abs(verifier_moments[["variance"]])
    ),
    pmf.l1 = truth[["pmf_l1", exact = TRUE]]
  )
  formula <- c(
    K.mean = "absolute_plus_relative_max",
    K.variance = "absolute_plus_relative_max",
    pmf.l1 = "direct_pmf_l1_tolerance"
  )
  scale_floor <- c(
    K.mean = truth[["moment_scale_floor", exact = TRUE]],
    K.variance = truth[["moment_scale_floor", exact = TRUE]],
    pmf.l1 = 0
  )
  if (identical(raw[["kind", exact = TRUE]], "interval")) {
    selected_mass <- .dpprior_target_interval_masses(
      selected_pmf, raw[["interval", exact = TRUE]]
    )
    verifier_mass <- .dpprior_target_interval_masses(
      verifier_pmf, raw[["interval", exact = TRUE]]
    )
    interval_delta <- abs(selected_mass - verifier_mass)
    names(interval_delta) <- paste0("interval.", names(interval_delta))
    interval_tolerance <- setNames(
      rep(truth[["constraint", exact = TRUE]], length(interval_delta)),
      names(interval_delta)
    )
    interval_formula <- setNames(
      rep("fixed_constraint_tolerance", length(interval_delta)),
      names(interval_delta)
    )
    interval_floor <- setNames(rep(0, length(interval_delta)), names(interval_delta))
    delta <- c(delta, interval_delta)
    tolerance <- c(tolerance, interval_tolerance)
    formula <- c(formula, interval_formula)
    scale_floor <- c(scale_floor, interval_floor)
  }
  .dpprior_new_stability(
    delta = delta, tolerance = tolerance, formula = formula,
    scale_floor = scale_floor, source = "independent_target_reconstruction"
  )
}


.dpprior_validate_target_truth_authority <- function(raw) {
  if (!(raw[["kind", exact = TRUE]] %in% c("interval", "family")) ||
      !(raw[["status", exact = TRUE]] %in% c("converged", "boundary"))) {
    return(invisible(TRUE))
  }
  truth <- raw[["tolerances", exact = TRUE]]
  .dpprior_schema_exact_names(
    truth, c("constraint", "pmf_l1", "moment_relative", "moment_scale_floor"),
    "target.tolerances"
  )
  for (field in c("constraint", "pmf_l1", "moment_relative",
                  "moment_scale_floor")) {
    .dpprior_schema_validate_finite_scalar(
      truth[[field, exact = TRUE]], paste0("target.tolerances.", field),
      lower = 0, lower_open = field %in% c("constraint", "pmf_l1",
                                           "moment_relative")
    )
  }
  .dpprior_schema_require(
    truth[["constraint", exact = TRUE]] >= 64 * .Machine$double.eps &&
      truth[["constraint", exact = TRUE]] <= 1e-9 &&
      identical(truth[["pmf_l1", exact = TRUE]], .TOL_PMF_SUM) &&
      identical(truth[["moment_relative", exact = TRUE]], 1e-8) &&
      identical(truth[["moment_scale_floor", exact = TRUE]], 1),
    "target_truth_tolerances", "target.tolerances",
    paste(
      "constraint tolerance between the numerical floor and the fixed 1e-9",
      "verified-target ceiling, plus fixed PMF/moment tolerances"
    ), truth
  )
  controls <- raw[["computation", exact = TRUE]][["used", exact = TRUE]][[
    "controls", exact = TRUE
  ]]
  .dpprior_schema_exact_names(
    controls,
    c(
      "constraint_tolerance", "root_tolerance", "max_iterations",
      "pmf_l1_tolerance", "moment_relative_tolerance", "moment_scale_floor"
    ),
    "target.computation.used.controls"
  )
  if (identical(raw[["kind", exact = TRUE]], "interval")) {
    .dpprior_schema_validate_finite_scalar(
      controls[["root_tolerance", exact = TRUE]],
      "target.computation.used.controls.root_tolerance",
      lower = .Machine$double.eps
    )
    .dpprior_schema_require(
      controls[["root_tolerance", exact = TRUE]] <= max(
        .Machine$double.eps,
        min(1e-12, truth[["constraint", exact = TRUE]] / 10)
      ) && .dpprior_schema_is_count(
        controls[["max_iterations", exact = TRUE]], 1L
      ) && controls[["max_iterations", exact = TRUE]] <= 1000L,
      "target_root_policy", "target.computation.used.controls",
      paste(
        "a root iteration cap in 1:1000 and root tolerance no looser than",
        "max(machine epsilon, min(1e-12, constraint_tolerance/10))"
      ), controls
    )
  } else {
    .dpprior_schema_require(
      is.null(controls[["root_tolerance", exact = TRUE]]) &&
        is.null(controls[["max_iterations", exact = TRUE]]),
      "target_family_controls", "target.computation.used.controls",
      "explicit NULL root controls for an analytic family construction",
      controls[c("root_tolerance", "max_iterations")]
    )
  }
  .dpprior_schema_require(
    identical(controls[["constraint_tolerance", exact = TRUE]],
              truth[["constraint", exact = TRUE]]) &&
      identical(controls[["pmf_l1_tolerance", exact = TRUE]],
                truth[["pmf_l1", exact = TRUE]]) &&
      identical(controls[["moment_relative_tolerance", exact = TRUE]],
                truth[["moment_relative", exact = TRUE]]) &&
      identical(controls[["moment_scale_floor", exact = TRUE]],
                truth[["moment_scale_floor", exact = TRUE]]),
    "target_truth_controls", "target.computation.used.controls",
    "exact identity with the canonical target truth-control spine", controls
  )
  request_controls <- raw[["computation", exact = TRUE]][[
    "request", exact = TRUE
  ]][["controls", exact = TRUE]]
  .dpprior_schema_require(
    identical(request_controls, controls) &&
      !raw[["computation", exact = TRUE]][["fallback", exact = TRUE]][[
        "attempted", exact = TRUE
      ]] &&
      !raw[["computation", exact = TRUE]][["fallback", exact = TRUE]][[
        "used", exact = TRUE
      ]],
    "target_requested_controls", "target.computation",
    "exact requested/used truth controls with no hidden fallback", controls
  )
  invisible(TRUE)
}


.dpprior_validate_target_derivation_entry <- function(entry, path) {
  .dpprior_schema_exact_names(
    entry, .DPPRIOR_TARGET_DERIVATION_ENTRY_FIELDS, path
  )
  for (field in c("rule", "outcome")) {
    .dpprior_schema_validate_scalar_character(
      entry[[field, exact = TRUE]], paste0(path, ".", field)
    )
  }
  .dpprior_schema_validate_scalar_logical(
    entry[["opt_in", exact = TRUE]], paste0(path, ".opt_in")
  )
  for (field in c("before", "after", "evidence")) {
    .dpprior_schema_validate_named_list(
      entry[[field, exact = TRUE]], paste0(path, ".", field)
    )
    .dpprior_schema_validate_plain_record_value(
      entry[[field, exact = TRUE]], paste0(path, ".", field)
    )
  }
  invisible(TRUE)
}


.dpprior_validate_target_derivation <- function(raw) {
  derivation <- raw[["derivation", exact = TRUE]]
  .dpprior_schema_exact_names(
    derivation, .DPPRIOR_TARGET_DERIVATION_FIELDS, "target.derivation"
  )
  first <- derivation[["request_to_normalized", exact = TRUE]]
  .dpprior_validate_target_derivation_entry(
    first, "target.derivation.request_to_normalized"
  )
  .dpprior_schema_require(
    first[["rule", exact = TRUE]] %in% .DPPRIOR_TARGET_NORMALIZATION_RULES &&
      identical(first[["outcome", exact = TRUE]], "canonicalized") &&
      !first[["opt_in", exact = TRUE]] &&
      identical(first[["before", exact = TRUE]],
                raw[["request", exact = TRUE]]) &&
      identical(first[["after", exact = TRUE]],
                raw[["normalized", exact = TRUE]]),
    "target_derivation_identity",
    "target.derivation.request_to_normalized",
    paste(
      "an approved non-opt-in canonicalization with exact request-before",
      "and normalized-after identities"
    ), first
  )
  expected_first_rules <- switch(
    raw[["kind", exact = TRUE]],
    moments = c("canonicalize_direct_moments",
                "canonicalize_confidence_target"),
    cv = "canonicalize_cv_target",
    pmf = c("validate_strict_pmf", "drop_structural_k0_zero"),
    interval = "canonicalize_interval_request",
    family = "canonicalize_scaled_chisq_family_request"
  )
  .dpprior_schema_require(
    first[["rule", exact = TRUE]] %in% expected_first_rules,
    "target_derivation_kind",
    "target.derivation.request_to_normalized.rule",
    "a normalization rule authorized for target.kind",
    first[["rule", exact = TRUE]]
  )
  first_rule <- first[["rule", exact = TRUE]]
  request_record <- first[["before", exact = TRUE]]
  normalized_record <- first[["after", exact = TRUE]]
  if (identical(first_rule, "canonicalize_direct_moments")) {
    .dpprior_schema_exact_names(
      request_record, c("J", "mu_K", "var_K"),
      "target.derivation.request_to_normalized.before"
    )
    .dpprior_schema_exact_names(
      normalized_record, c("J", "mu_K", "var_K", "interval", "pmf"),
      "target.derivation.request_to_normalized.after"
    )
  } else if (identical(first_rule, "canonicalize_confidence_target")) {
    .dpprior_schema_exact_names(
      request_record, c("J", "mean", "confidence"),
      "target.derivation.request_to_normalized.before"
    )
    .dpprior_schema_exact_names(
      normalized_record, c("J", "mean", "confidence", "interval", "pmf"),
      "target.derivation.request_to_normalized.after"
    )
  } else if (identical(first_rule, "canonicalize_cv_target")) {
    .dpprior_schema_exact_names(
      request_record, c("J", "mean", "cv"),
      "target.derivation.request_to_normalized.before"
    )
    .dpprior_schema_exact_names(
      normalized_record, c("J", "mean", "cv", "interval", "pmf"),
      "target.derivation.request_to_normalized.after"
    )
  } else if (first_rule %in%
             c("validate_strict_pmf", "drop_structural_k0_zero")) {
    .dpprior_schema_exact_names(
      request_record, c("J", "pmf"),
      "target.derivation.request_to_normalized.before"
    )
    .dpprior_schema_exact_names(
      normalized_record, c("J", "pmf", "interval"),
      "target.derivation.request_to_normalized.after"
    )
  } else if (identical(first_rule, "canonicalize_interval_request")) {
    .dpprior_schema_exact_names(
      request_record, c("J", "K_interval", "mu_K"),
      "target.derivation.request_to_normalized.before"
    )
    .dpprior_schema_exact_names(
      normalized_record, c("J", "interval", "family", "pmf"),
      "target.derivation.request_to_normalized.after"
    )
    request_mean <- request_record[["mu_K", exact = TRUE]]
    if (!is.null(request_mean)) {
      .dpprior_schema_validate_finite_scalar(
        request_mean,
        "target.derivation.request_to_normalized.before.mu_K",
        lower = 1, upper = raw[["J", exact = TRUE]]
      )
    }
    .dpprior_schema_require(
      identical(request_record[["J", exact = TRUE]],
                normalized_record[["J", exact = TRUE]]) &&
        identical(request_record[["K_interval", exact = TRUE]],
                  normalized_record[["interval", exact = TRUE]]) &&
        identical(
          request_mean,
          normalized_record[["interval", exact = TRUE]][["mu_K", exact = TRUE]]
        ) && is.null(normalized_record[["pmf", exact = TRUE]]),
      "interval_request_authority",
      "target.derivation.request_to_normalized",
      paste(
        "exact J, requested K_interval, and requested mu_K identities in",
        "the pre-construction normalized interval record"
      ),
      list(request = request_record, normalized = normalized_record)
    )
  } else if (identical(
    first_rule, "canonicalize_scaled_chisq_family_request"
  )) {
    .dpprior_schema_exact_names(
      request_record, c("J", "mean", "variance", "family"),
      "target.derivation.request_to_normalized.before"
    )
    .dpprior_schema_exact_names(
      normalized_record,
      c("J", "mean", "variance", "interval", "family", "pmf"),
      "target.derivation.request_to_normalized.after"
    )
    expected_family <- list(
      name = "scaled_chisq",
      parameterization = "df_and_scale_from_requested_moments",
      explicit = TRUE,
      reference_measure = paste(
        "continuity-corrected chi-square bins conditioned on 1:J"
      )
    )
    .dpprior_schema_require(
      identical(request_record[["J", exact = TRUE]],
                normalized_record[["J", exact = TRUE]]) &&
        identical(request_record[["mean", exact = TRUE]],
                  normalized_record[["mean", exact = TRUE]]) &&
        identical(request_record[["variance", exact = TRUE]],
                  normalized_record[["variance", exact = TRUE]]) &&
        identical(request_record[["family", exact = TRUE]], expected_family) &&
        identical(normalized_record[["family", exact = TRUE]],
                  expected_family) &&
        identical(raw[["family", exact = TRUE]], expected_family) &&
        is.null(normalized_record[["interval", exact = TRUE]]) &&
        is.null(normalized_record[["pmf", exact = TRUE]]),
      "scaled_chisq_request_authority",
      "target.derivation.request_to_normalized",
      paste(
        "exact J/moments and the closed scaled-chi-square family metadata",
        "from request through normalized authority"
      ),
      list(request = request_record, normalized = normalized_record)
    )
  }
  moment_pair <- function(record) {
    fields <- if (all(c("mean", "variance") %in% names(record))) {
      c("mean", "variance")
    } else if (all(c("mu_K", "var_K") %in% names(record))) {
      c("mu_K", "var_K")
    } else {
      return(NULL)
    }
    unname(unlist(record[fields], use.names = FALSE))
  }
  if (identical(first[["rule", exact = TRUE]],
                "canonicalize_direct_moments")) {
    .dpprior_schema_require(
      identical(moment_pair(first[["before", exact = TRUE]]),
                moment_pair(first[["after", exact = TRUE]])),
      "target_moment_identity",
      "target.derivation.request_to_normalized",
      "unchanged requested moments under direct canonicalization", first
    )
  }
  if (first[["rule", exact = TRUE]] %in%
      c("validate_strict_pmf", "drop_structural_k0_zero")) {
    input_pmf <- first[["before", exact = TRUE]][["pmf", exact = TRUE]]
    normalized_pmf <- first[["after", exact = TRUE]][["pmf", exact = TRUE]]
    pmf_identity <- if (identical(first[["rule", exact = TRUE]],
                                  "validate_strict_pmf")) {
      identical(input_pmf, normalized_pmf)
    } else {
      is.numeric(input_pmf) && length(input_pmf) == raw$J + 1L &&
        identical(input_pmf[[1L]], 0) &&
        identical(input_pmf[-1L], normalized_pmf)
    }
    .dpprior_schema_require(
      pmf_identity, "target_pmf_derivation",
      "target.derivation.request_to_normalized",
      "strict PMF identity or documented removal of one exact structural K=0 zero",
      first
    )
  }
  if (identical(first[["rule", exact = TRUE]],
                "canonicalize_confidence_target")) {
    .dpprior_schema_require(
      identical(first[["after", exact = TRUE]][["J", exact = TRUE]],
                first[["before", exact = TRUE]][["J", exact = TRUE]]) &&
        identical(first[["after", exact = TRUE]][["mean", exact = TRUE]],
                  first[["before", exact = TRUE]][["mean", exact = TRUE]]) &&
        identical(first[["after", exact = TRUE]][[
          "confidence", exact = TRUE
        ]], first[["before", exact = TRUE]][[
          "confidence", exact = TRUE
        ]]) &&
        is.null(first[["after", exact = TRUE]][["interval", exact = TRUE]]) &&
        is.null(first[["after", exact = TRUE]][["pmf", exact = TRUE]]),
      "confidence_canonicalization",
      "target.derivation.request_to_normalized",
      "unchanged J/mean/confidence plus explicit NULL interval/pmf", first
    )
  }
  if (identical(first[["rule", exact = TRUE]],
                "canonicalize_cv_target")) {
    .dpprior_schema_require(
      identical(first[["after", exact = TRUE]][["J", exact = TRUE]],
                first[["before", exact = TRUE]][["J", exact = TRUE]]) &&
        identical(first[["after", exact = TRUE]][["mean", exact = TRUE]],
                  first[["before", exact = TRUE]][["mean", exact = TRUE]]) &&
        identical(first[["after", exact = TRUE]][["cv", exact = TRUE]],
                  first[["before", exact = TRUE]][["cv", exact = TRUE]]) &&
        is.null(first[["after", exact = TRUE]][["interval", exact = TRUE]]) &&
        is.null(first[["after", exact = TRUE]][["pmf", exact = TRUE]]),
      "cv_canonicalization", "target.derivation.request_to_normalized",
      "unchanged J/mean/CV plus explicit NULL interval/pmf", first
    )
  }

  second <- derivation[["normalized_to_used", exact = TRUE]]
  changed <- !identical(
    raw[["normalized", exact = TRUE]], raw[["used", exact = TRUE]]
  )
  if (!changed) {
    .dpprior_schema_require(
      is.null(second), "target_derivation_identity",
      "target.derivation.normalized_to_used",
      "NULL exactly when normalized equals used", second
    )
    return(invisible(TRUE))
  }
  .dpprior_validate_target_derivation_entry(
    second, "target.derivation.normalized_to_used"
  )
  projection_rule <- identical(
    second[["rule", exact = TRUE]],
    "project_a1_variance_to_nearest_interior"
  )
  .dpprior_schema_require(
    second[["rule", exact = TRUE]] %in% .DPPRIOR_TARGET_USE_RULES &&
      second[["outcome", exact = TRUE]] %in% c("derived", "projected") &&
      identical(second[["outcome", exact = TRUE]],
                if (projection_rule) "projected" else "derived") &&
      identical(second[["opt_in", exact = TRUE]], projection_rule) &&
      identical(second[["before", exact = TRUE]],
                raw[["normalized", exact = TRUE]]) &&
      identical(second[["after", exact = TRUE]],
                raw[["used", exact = TRUE]]),
    "target_derivation_identity", "target.derivation.normalized_to_used",
    paste(
      "an approved derivation/projection with exact normalized-before and",
      "used-after identities; opt-in only for A1 projection"
    ), second
  )
  rule <- second[["rule", exact = TRUE]]
  before <- second[["before", exact = TRUE]]
  after <- second[["after", exact = TRUE]]
  evidence <- second[["evidence", exact = TRUE]]
  if (identical(rule, "derive_variance_from_confidence_vif")) {
    .dpprior_schema_exact_names(
      after, c("J", "mean", "variance", "interval", "pmf"),
      "target.derivation.normalized_to_used.after"
    )
    .dpprior_schema_exact_names(
      evidence, c("confidence", "vif", "formula"),
      "target.derivation.normalized_to_used.evidence"
    )
    .dpprior_schema_validate_scalar_character(
      evidence[["confidence", exact = TRUE]],
      "target.derivation.normalized_to_used.evidence.confidence"
    )
    confidence_vif <- c(low = 5, medium = 2.5, high = 1.5)
    .dpprior_schema_require(
      evidence[["confidence", exact = TRUE]] %in% names(confidence_vif),
      "confidence_level",
      "target.derivation.normalized_to_used.evidence.confidence",
      "low, medium, or high", evidence[["confidence", exact = TRUE]]
    )
    .dpprior_schema_validate_finite_scalar(
      evidence[["vif", exact = TRUE]],
      "target.derivation.normalized_to_used.evidence.vif",
      lower = 0, lower_open = TRUE
    )
    expected_variance <- evidence[["vif", exact = TRUE]] *
      (before[["mean", exact = TRUE]] - 1)
    .dpprior_schema_require(
      identical(first[["rule", exact = TRUE]],
                "canonicalize_confidence_target") &&
        identical(evidence[["confidence", exact = TRUE]],
                  before[["confidence", exact = TRUE]]) &&
        identical(evidence[["vif", exact = TRUE]],
                  unname(confidence_vif[[
                    evidence[["confidence", exact = TRUE]]
                  ]])) &&
        identical(evidence[["formula", exact = TRUE]],
                  "variance = vif * (mean - 1)") &&
        identical(after[["mean", exact = TRUE]],
                  before[["mean", exact = TRUE]]) &&
        identical(after[["variance", exact = TRUE]], expected_variance),
      "confidence_derivation_math",
      "target.derivation.normalized_to_used",
      "the retained confidence/VIF variance formula and exact result", second
    )
  }
  if (identical(rule, "derive_variance_from_cv")) {
    .dpprior_schema_exact_names(
      after, c("J", "mean", "variance", "interval", "pmf"),
      "target.derivation.normalized_to_used.after"
    )
    .dpprior_schema_exact_names(
      evidence, c("cv", "definition", "formula"),
      "target.derivation.normalized_to_used.evidence"
    )
    .dpprior_schema_validate_finite_scalar(
      evidence[["cv", exact = TRUE]],
      "target.derivation.normalized_to_used.evidence.cv", lower = 0
    )
    expected_variance <- (evidence[["cv", exact = TRUE]] *
                            before[["mean", exact = TRUE]])^2
    .dpprior_schema_require(
      identical(first[["rule", exact = TRUE]], "canonicalize_cv_target") &&
        identical(evidence[["cv", exact = TRUE]],
                  before[["cv", exact = TRUE]]) &&
        identical(evidence[["definition", exact = TRUE]],
                  "SD(K_J) / E(K_J)") &&
        identical(evidence[["formula", exact = TRUE]],
                  "variance = (cv * mean)^2") &&
        identical(after[["mean", exact = TRUE]],
                  before[["mean", exact = TRUE]]) &&
        identical(after[["variance", exact = TRUE]], expected_variance),
      "cv_derivation_math", "target.derivation.normalized_to_used",
      "the retained CV variance formula and exact result", second
    )
  }
  if (rule %in% c(
    "construct_maxent_hard_bounds_pmf",
    "construct_maxent_equal_tail_pmf",
    "construct_maxent_central_mass_pmf"
  )) {
    .dpprior_schema_exact_names(
      after, c("J", "interval", "family", "pmf"),
      "target.derivation.normalized_to_used.after"
    )
    maxent_fields <- c(
      "constructor_method", "group_masses", "common_tilt", "boundary",
      "mean_constraint", "feasible_mean_lower", "feasible_mean_upper",
      "feasibility_tolerance", "root_tolerance"
    )
    .dpprior_schema_exact_names(
      evidence, maxent_fields,
      "target.derivation.normalized_to_used.evidence"
    )
    .dpprior_schema_validate_scalar_character(
      evidence[["constructor_method", exact = TRUE]],
      "target.derivation.normalized_to_used.evidence.constructor_method"
    )
    for (field in c("feasible_mean_lower", "feasible_mean_upper",
                    "feasibility_tolerance", "root_tolerance")) {
      .dpprior_schema_validate_finite_scalar(
        evidence[[field, exact = TRUE]],
        paste0("target.derivation.normalized_to_used.evidence.", field),
        lower = if (field %in% c("feasibility_tolerance", "root_tolerance"))
          0 else -Inf
      )
    }
    interval <- raw[["interval", exact = TRUE]]
    pmf <- raw[["pmf", exact = TRUE]]
    family <- raw[["family", exact = TRUE]]
    .dpprior_schema_exact_names(
      interval,
      c(
        "lower", "upper", "type", "coverage", "family", "mu_K",
        "support", "endpoints"
      ),
      "target.interval"
    )
    .dpprior_schema_exact_names(
      family,
      c("name", "parameterization", "explicit", "reference_measure"),
      "target.family"
    )
    .dpprior_schema_require(
      identical(family[["name", exact = TRUE]], "maxent") &&
        identical(family[["parameterization", exact = TRUE]],
                  "common exponential tilt within fixed-mass groups") &&
        identical(family[["explicit", exact = TRUE]], TRUE) &&
        identical(family[["reference_measure", exact = TRUE]],
                  "counting measure on 1:J") &&
        identical(interval[["family", exact = TRUE]], "maxent") &&
        identical(interval[["support", exact = TRUE]],
                  c(lower = 1L, upper = as.integer(raw[["J", exact = TRUE]]))) &&
        identical(interval[["endpoints", exact = TRUE]], "inclusive"),
      "maxent_family", "target.family",
      "the closed canonical MaxEnt construction family and support", family
    )
    lower <- interval[["lower", exact = TRUE]]
    upper <- interval[["upper", exact = TRUE]]
    .dpprior_schema_validate_finite_scalar(
      interval[["coverage", exact = TRUE]], "target.interval.coverage",
      lower = 0, upper = 1
    )
    .dpprior_schema_require(
      .dpprior_schema_is_count(lower, 1L) &&
        .dpprior_schema_is_count(upper, lower) && upper <= raw$J,
      "maxent_interval", "target.interval",
      "integer inclusive lower/upper endpoints within 1:J", interval
    )
    support <- seq_len(raw$J)
    expected_type <- switch(
      rule,
      construct_maxent_hard_bounds_pmf = "hard_bounds",
      construct_maxent_equal_tail_pmf = "equal_tail",
      construct_maxent_central_mass_pmf = "central_mass"
    )
    coverage <- interval[["coverage", exact = TRUE]]
    groups <- switch(
      expected_type,
      hard_bounds = list(inside = support[support >= lower & support <= upper]),
      equal_tail = list(
        left = support[support < lower],
        inside = support[support >= lower & support <= upper],
        right = support[support > upper]
      ),
      central_mass = list(
        inside = support[support >= lower & support <= upper],
        outside = support[support < lower | support > upper]
      )
    )
    expected_group_masses <- switch(
      expected_type,
      hard_bounds = c(inside = 1),
      equal_tail = c(
        left = (1 - coverage) / 2,
        inside = coverage,
        right = (1 - coverage) / 2
      ),
      central_mass = c(inside = coverage, outside = 1 - coverage)
    )
    .dpprior_schema_exact_names(
      evidence[["group_masses", exact = TRUE]], names(expected_group_masses),
      "target.derivation.normalized_to_used.evidence.group_masses"
    )
    group_masses <- vapply(
      names(expected_group_masses),
      function(group) {
        value <- evidence[["group_masses", exact = TRUE]][[group, exact = TRUE]]
        .dpprior_schema_validate_finite_scalar(
          value,
          paste0(
            "target.derivation.normalized_to_used.evidence.group_masses.",
            group
          ),
          lower = 0, upper = 1
        )
        value
      },
      numeric(1)
    )
    .dpprior_schema_require(
      identical(interval[["type", exact = TRUE]], expected_type) &&
        (if (identical(expected_type, "hard_bounds")) {
          identical(coverage, 1)
        } else {
          TRUE
        }) &&
        all(abs(group_masses - expected_group_masses) <= .TOL_PMF_SUM) &&
        abs(sum(group_masses) - 1) <= .TOL_PMF_SUM,
      "maxent_group_masses",
      "target.derivation.normalized_to_used.evidence.group_masses",
      "the exact interval-type constraint masses", group_masses
    )
    for (group in names(groups)) {
      .dpprior_schema_require(
        group_masses[[group]] == 0 || length(groups[[group]]) > 0L,
        "maxent_positive_group",
        "target.derivation.normalized_to_used.evidence.group_masses",
        "positive mass only for a non-empty canonical support group",
        list(group = group, mass = group_masses[[group]])
      )
    }
    positive_groups <- names(groups)[group_masses > 0]
    expected_lower <- sum(vapply(
      positive_groups,
      function(group) group_masses[[group]] * min(groups[[group]]), numeric(1)
    ))
    expected_upper <- sum(vapply(
      positive_groups,
      function(group) group_masses[[group]] * max(groups[[group]]), numeric(1)
    ))
    interval_mean <- interval[["mu_K", exact = TRUE]]
    if (!is.null(interval_mean)) {
      .dpprior_schema_validate_finite_scalar(
        interval_mean, "target.interval.mu_K", lower = 1, upper = raw$J
      )
    }
    before_means <- c()
    for (field in c("mean", "mu_K")) {
      if (field %in% names(before) && !is.null(before[[field, exact = TRUE]])) {
        .dpprior_schema_validate_finite_scalar(
          before[[field, exact = TRUE]],
          paste0("target.derivation.normalized_to_used.before.", field),
          lower = 1, upper = raw$J
        )
        before_means <- c(before_means, before[[field, exact = TRUE]])
      }
    }
    mean_values <- c(before_means, interval_mean)
    expected_mean <- if (length(mean_values)) mean_values[[1L]] else NULL
    .dpprior_schema_require(
      !length(mean_values) || all(mean_values == expected_mean),
      "maxent_mean_authority", "target.interval.mu_K",
      "one consistent normalized/interval mean constraint", mean_values
    )
    .dpprior_schema_require(
      identical(evidence[["mean_constraint", exact = TRUE]], expected_mean),
      "maxent_mean_authority",
      "target.derivation.normalized_to_used.evidence.mean_constraint",
      "identity with the authoritative interval mean constraint",
      evidence[["mean_constraint", exact = TRUE]]
    )
    expected_feasibility_tolerance <- 64 * .Machine$double.eps * max(
      1, abs(expected_lower), abs(expected_upper),
      if (is.null(expected_mean)) 1 else abs(expected_mean)
    )
    .dpprior_schema_require(
      abs(evidence[["feasible_mean_lower", exact = TRUE]] - expected_lower) <=
        64 * .Machine$double.eps * max(1, abs(expected_lower)) &&
        abs(evidence[["feasible_mean_upper", exact = TRUE]] - expected_upper) <=
          64 * .Machine$double.eps * max(1, abs(expected_upper)) &&
        identical(evidence[["feasibility_tolerance", exact = TRUE]],
                  expected_feasibility_tolerance) &&
        (is.null(expected_mean) ||
          (expected_mean >= expected_lower - expected_feasibility_tolerance &&
             expected_mean <= expected_upper + expected_feasibility_tolerance)),
      "maxent_feasibility", "target.derivation.normalized_to_used.evidence",
      "the analytic support-group mean hull and fixed rounding tolerance",
      evidence[c(
        "feasible_mean_lower", "feasible_mean_upper", "feasibility_tolerance"
      )]
    )
    expected_side <- if (!is.null(expected_mean) &&
                         abs(expected_mean - expected_lower) <=
                           expected_feasibility_tolerance) {
      "minimum"
    } else if (!is.null(expected_mean) &&
               abs(expected_mean - expected_upper) <=
                 expected_feasibility_tolerance) {
      "maximum"
    } else if (is.null(expected_mean) &&
               all(vapply(groups[group_masses > 0], length, integer(1)) == 1L)) {
      "minimum"
    } else {
      NULL
    }
    boundary <- evidence[["boundary", exact = TRUE]]
    .dpprior_schema_exact_names(
      boundary, c("active", "side"),
      "target.derivation.normalized_to_used.evidence.boundary"
    )
    .dpprior_schema_validate_scalar_logical(
      boundary[["active", exact = TRUE]],
      "target.derivation.normalized_to_used.evidence.boundary.active"
    )
    if (!is.null(boundary[["side", exact = TRUE]])) {
      .dpprior_schema_validate_scalar_character(
        boundary[["side", exact = TRUE]],
        "target.derivation.normalized_to_used.evidence.boundary.side"
      )
    }
    .dpprior_schema_require(
      identical(boundary[["active", exact = TRUE]], !is.null(expected_side)) &&
        identical(boundary[["side", exact = TRUE]], expected_side),
      "maxent_boundary", "target.derivation.normalized_to_used.evidence.boundary",
      "the boundary state recomputed from the analytic mean hull", boundary
    )
    common_tilt <- evidence[["common_tilt", exact = TRUE]]
    if (is.null(expected_side)) {
      .dpprior_schema_validate_finite_scalar(
        common_tilt,
        "target.derivation.normalized_to_used.evidence.common_tilt"
      )
    } else {
      .dpprior_schema_require(
        is.null(common_tilt), "maxent_boundary_tilt",
        "target.derivation.normalized_to_used.evidence.common_tilt",
        "NULL for an analytic boundary construction", common_tilt
      )
    }
    reconstructed_pmf <- numeric(raw$J)
    for (group in names(groups)) {
      group_support <- groups[[group]]
      mass <- group_masses[[group]]
      if (mass > 0) {
        if (!is.null(expected_side)) {
          point <- if (identical(expected_side, "minimum")) {
            min(group_support)
          } else {
            max(group_support)
          }
          reconstructed_pmf[[point]] <- reconstructed_pmf[[point]] + mass
        } else {
          logits <- common_tilt * group_support
          weights <- exp(logits - max(logits))
          reconstructed_pmf[group_support] <- mass * weights / sum(weights)
        }
      }
    }
    expected_constructor <- if (!is.null(expected_side)) {
      "boundary"
    } else if (is.null(expected_mean)) {
      "analytic"
    } else {
      "uniroot"
    }
    expected_root_tolerance <- raw[["computation", exact = TRUE]][[
      "used", exact = TRUE
    ]][["controls", exact = TRUE]][["root_tolerance", exact = TRUE]]
    mean_value <- sum(support * pmf)
    mean_ok <- is.null(expected_mean) || abs(mean_value - expected_mean) <=
      raw[["tolerances", exact = TRUE]][["constraint", exact = TRUE]] *
        max(1, abs(expected_mean))
    .dpprior_schema_require(
      identical(first[["rule", exact = TRUE]],
                "canonicalize_interval_request") &&
        identical(raw[["kind", exact = TRUE]], "interval") &&
        mean_ok &&
        identical(evidence[["constructor_method", exact = TRUE]],
                  expected_constructor) &&
        identical(evidence[["root_tolerance", exact = TRUE]],
                  expected_root_tolerance) &&
        (!is.null(expected_mean) || identical(common_tilt, 0)) &&
        max(abs(reconstructed_pmf - pmf)) <= .TOL_PMF_SUM &&
        all(abs(group_masses - expected_group_masses) <= .TOL_PMF_SUM),
      "maxent_derivation_math", "target.derivation.normalized_to_used",
      paste(
        "PMF reconstructed from canonical support groups/common tilt, fixed",
        "truth controls, analytic feasibility bounds, and mean authority"
      ),
      list(evidence = evidence, reconstructed_pmf = reconstructed_pmf)
    )
  }
  if (identical(rule, "construct_scaled_chisq_conditioned_pmf")) {
    .dpprior_schema_exact_names(
      after, c("J", "interval", "family", "pmf"),
      "target.derivation.normalized_to_used.after"
    )
    chisq_fields <- c(
      "df", "scale", "binning", "retained_mass_before_normalization",
      "omitted_mass", "normalization", "support"
    )
    .dpprior_schema_exact_names(
      evidence, chisq_fields,
      "target.derivation.normalized_to_used.evidence"
    )
    for (field in c("df", "scale", "retained_mass_before_normalization",
                    "omitted_mass")) {
      .dpprior_schema_validate_finite_scalar(
        evidence[[field, exact = TRUE]],
        paste0("target.derivation.normalized_to_used.evidence.", field),
        lower = 0, lower_open = field %in% c("df", "scale")
      )
    }
    requested_moments <- moment_pair(before)
    .dpprior_schema_require(
      length(requested_moments) == 2L && requested_moments[[2L]] > 0,
      "chisq_requested_moments",
      "target.derivation.normalized_to_used.before",
      "a finite positive requested variance for scaled-chi-square construction",
      requested_moments
    )
    expected_df <- 2 * requested_moments[[1L]]^2 /
      requested_moments[[2L]]
    expected_scale <- requested_moments[[2L]] /
      (2 * requested_moments[[1L]])
    numerical_tolerance <- 1e-12 * max(
      1, abs(expected_df), abs(expected_scale),
      abs(evidence[["df", exact = TRUE]]),
      abs(evidence[["scale", exact = TRUE]])
    )
    edges_lower <- (seq_len(raw$J) - 0.5) /
      evidence[["scale", exact = TRUE]]
    edges_upper <- (seq_len(raw$J) + 0.5) /
      evidence[["scale", exact = TRUE]]
    raw_mass <- stats::pchisq(
      edges_upper, df = evidence[["df", exact = TRUE]]
    ) - stats::pchisq(
      edges_lower, df = evidence[["df", exact = TRUE]]
    )
    raw_mass <- pmax(raw_mass, 0)
    retained_mass <- sum(raw_mass)
    reconstructed_pmf <- raw_mass / retained_mass
    .dpprior_schema_require(
      identical(first[["rule", exact = TRUE]],
                "canonicalize_scaled_chisq_family_request") &&
        identical(raw[["kind", exact = TRUE]], "family") &&
        identical(evidence[["binning", exact = TRUE]],
                  "continuity_corrected_half_integer_bins") &&
        identical(evidence[["normalization", exact = TRUE]],
                  "explicit_support_conditioning") &&
        identical(evidence[["support", exact = TRUE]],
                  c(lower = 1L, upper = as.integer(raw$J))) &&
        abs(evidence[["retained_mass_before_normalization", exact = TRUE]] +
              evidence[["omitted_mass", exact = TRUE]] - 1) <=
          .TOL_PMF_SUM &&
        abs(evidence[["df", exact = TRUE]] - expected_df) <=
          numerical_tolerance &&
        abs(evidence[["scale", exact = TRUE]] - expected_scale) <=
          numerical_tolerance &&
        abs(evidence[["retained_mass_before_normalization", exact = TRUE]] -
              retained_mass) <= 1e-12 * max(1, abs(retained_mass)) &&
        max(abs(reconstructed_pmf - raw[["pmf", exact = TRUE]])) <=
          .TOL_PMF_SUM,
      "chisq_derivation_math", "target.derivation.normalized_to_used",
      paste(
        "df/scale from requested moments plus exact continuity-corrected",
        "binning and explicit support conditioning"
      ),
      evidence
    )
  }
  if (projection_rule) {
    projection_fields <- c(
      "policy", "opt_in", "applied", "reason", "original_target",
      "projected_target", "distance", "a1_lower_bound",
      "numerical_interior_lower", "support_upper_bound", "epsilon", "buffer",
      "requested_buffer", "preferred_buffer", "representability_floor",
      "effective_buffer", "available_gap", "buffer_was_capped",
      "buffer_was_floored", "representable_interior"
    )
    .dpprior_schema_exact_names(
      evidence, projection_fields,
      "target.derivation.normalized_to_used.evidence"
    )
    for (field in c(
      "distance", "a1_lower_bound", "numerical_interior_lower",
      "support_upper_bound", "epsilon", "buffer", "requested_buffer",
      "preferred_buffer", "representability_floor", "effective_buffer",
      "available_gap"
    )) {
      .dpprior_schema_validate_finite_scalar(
        evidence[[field, exact = TRUE]],
        paste0("target.derivation.normalized_to_used.evidence.", field),
        lower = if (field %in% c(
          "distance", "epsilon", "buffer", "requested_buffer", "preferred_buffer",
          "representability_floor", "effective_buffer", "available_gap"
        )) 0 else -Inf
      )
    }
    for (field in c(
      "opt_in", "applied", "buffer_was_capped", "buffer_was_floored",
      "representable_interior"
    )) {
      .dpprior_schema_validate_scalar_logical(
        evidence[[field, exact = TRUE]],
        paste0("target.derivation.normalized_to_used.evidence.", field)
      )
    }
    .dpprior_schema_validate_scalar_character(
      evidence[["policy", exact = TRUE]],
      "target.derivation.normalized_to_used.evidence.policy"
    )
    .dpprior_schema_validate_scalar_character(
      evidence[["reason", exact = TRUE]],
      "target.derivation.normalized_to_used.evidence.reason"
    )
    for (field in c("original_target", "projected_target")) {
      .dpprior_schema_exact_names(
        evidence[[field, exact = TRUE]], c("mu_K", "var_K"),
        paste0("target.derivation.normalized_to_used.evidence.", field)
      )
      for (moment in c("mu_K", "var_K")) {
        .dpprior_schema_validate_finite_scalar(
          evidence[[field, exact = TRUE]][[moment, exact = TRUE]],
          paste0(
            "target.derivation.normalized_to_used.evidence.", field, ".",
            moment
          )
        )
      }
    }
    before_moments <- moment_pair(before)
    after_moments <- moment_pair(after)
    .dpprior_schema_require(
      !is.null(before_moments) && length(before_moments) == 2L &&
        !is.null(after_moments) && length(after_moments) == 2L,
      "projection_moments", "target.derivation.normalized_to_used",
      "authoritative before/after mean and variance pairs", second
    )
    expected_lower <- before_moments[[1L]] - 1
    expected_upper <- expected_lower * (raw[["J", exact = TRUE]] -
      before_moments[[1L]])
    expected_requested_buffer <- evidence[["epsilon", exact = TRUE]] *
      (1 + expected_lower^2)
    expected_available_gap <- expected_upper - expected_lower
    expected_preferred_buffer <- if (expected_available_gap > 0) {
      min(expected_requested_buffer, expected_available_gap / 2)
    } else {
      expected_requested_buffer
    }
    expected_floor <- .Machine$double.eps * max(1, abs(expected_lower))
    expected_selected_buffer <- if (expected_available_gap > 0) {
      min(
        max(expected_preferred_buffer, expected_floor), expected_available_gap
      )
    } else {
      expected_preferred_buffer
    }
    expected_numerical_lower <- expected_lower + expected_selected_buffer
    if (expected_available_gap > 0 &&
        !(is.finite(expected_numerical_lower) &&
          expected_numerical_lower > expected_lower &&
          expected_numerical_lower <= expected_upper)) {
      expected_numerical_lower <- expected_upper
    }
    expected_effective_buffer <- expected_numerical_lower - expected_lower
    expected_representable <- expected_available_gap > 0 &&
      is.finite(expected_numerical_lower) &&
      expected_numerical_lower > expected_lower &&
      expected_numerical_lower <= expected_upper
    expected_capped <- expected_available_gap > 0 &&
      expected_requested_buffer > expected_available_gap / 2
    expected_floored <- expected_selected_buffer > expected_preferred_buffer
    expected_original <- list(
      mu_K = before_moments[[1L]], var_K = before_moments[[2L]]
    )
    expected_projected <- list(
      mu_K = after_moments[[1L]], var_K = expected_numerical_lower
    )
    .dpprior_schema_require(
      !is.null(before_moments) && !is.null(after_moments) &&
        identical(before_moments[[1L]], after_moments[[1L]]) &&
        identical(evidence[["policy", exact = TRUE]], "nearest") &&
        evidence[["opt_in", exact = TRUE]] &&
        evidence[["applied", exact = TRUE]] &&
        identical(evidence[["reason", exact = TRUE]],
                  "a1_strict_lower_bound") &&
        identical(evidence[["original_target", exact = TRUE]],
                  expected_original) &&
        identical(evidence[["projected_target", exact = TRUE]],
                  expected_projected) &&
        identical(evidence[["a1_lower_bound", exact = TRUE]],
                  expected_lower) &&
        identical(evidence[["support_upper_bound", exact = TRUE]],
                  expected_upper) &&
        identical(evidence[["requested_buffer", exact = TRUE]],
                  expected_requested_buffer) &&
        identical(evidence[["preferred_buffer", exact = TRUE]],
                  expected_preferred_buffer) &&
        identical(evidence[["representability_floor", exact = TRUE]],
                  expected_floor) &&
        identical(evidence[["numerical_interior_lower", exact = TRUE]],
                  expected_numerical_lower) &&
        identical(evidence[["effective_buffer", exact = TRUE]],
                  expected_effective_buffer) &&
        identical(evidence[["buffer", exact = TRUE]],
                  expected_effective_buffer) &&
        identical(evidence[["available_gap", exact = TRUE]],
                  expected_available_gap) &&
        identical(evidence[["buffer_was_capped", exact = TRUE]],
                  expected_capped) &&
        identical(evidence[["buffer_was_floored", exact = TRUE]],
                  expected_floored) &&
        identical(evidence[["representable_interior", exact = TRUE]],
                  expected_representable) &&
        identical(after_moments[[2L]], expected_numerical_lower) &&
        identical(evidence[["distance", exact = TRUE]],
                  expected_numerical_lower - before_moments[[2L]]) &&
        before_moments[[2L]] <= expected_lower && expected_representable,
      "projection_derivation_math",
      "target.derivation.normalized_to_used.evidence",
      paste(
        "unchanged mean, exact projection distance/gap, and projected",
        "variance inside the retained representable interior"
      ), evidence
    )
    projection <- raw[["provenance", exact = TRUE]][[
      "projection", exact = TRUE
    ]]
    .dpprior_schema_exact_names(
      projection[["record", exact = TRUE]], c("before", "after", "authority"),
      "target.provenance.projection.record"
    )
    .dpprior_schema_require(
      projection[["applied", exact = TRUE]] &&
        projection[["opt_in", exact = TRUE]] &&
        identical(projection[["record", exact = TRUE]][[
          "before", exact = TRUE
        ]], second[["before", exact = TRUE]]) &&
        identical(projection[["record", exact = TRUE]][[
          "after", exact = TRUE
        ]], second[["after", exact = TRUE]]) &&
        identical(projection[["record", exact = TRUE]][[
          "authority", exact = TRUE
        ]], evidence) &&
        identical(projection[["policy", exact = TRUE]],
                  evidence[["policy", exact = TRUE]]),
      "target_projection_authority", "target.provenance.projection",
      "matching opt-in projection provenance with exact before/after",
      projection
    )
  }
  invisible(TRUE)
}


.dpprior_validate_constructed_target_snapshot <- function(raw, snapshot, path) {
  .dpprior_schema_require(
    !is.null(snapshot), "target_verifier_snapshot", path,
    "a retained finite target-construction snapshot", snapshot
  )
  achieved <- snapshot[["achieved", exact = TRUE]]
  .dpprior_schema_exact_names(
    achieved, c("implied", "interval", "pmf"), paste0(path, ".achieved")
  )
  pmf <- achieved[["pmf", exact = TRUE]]
  .dpprior_schema_require(
    is.numeric(pmf) && !is.object(pmf) && is.null(dim(pmf)) &&
      .dpprior_schema_has_only_attributes(pmf) &&
      length(pmf) == raw[["J", exact = TRUE]] && !anyNA(pmf) &&
      all(is.finite(pmf)) && all(pmf >= 0) &&
      abs(sum(pmf) - 1) <= raw[["tolerances", exact = TRUE]][[
        "pmf_l1", exact = TRUE
      ]],
    "target_snapshot_pmf", paste0(path, ".achieved.pmf"),
    "a finite nonnegative normalized PMF of length J", pmf
  )
  implied <- achieved[["implied", exact = TRUE]]
  .dpprior_schema_exact_names(
    implied, c("mean", "variance"), paste0(path, ".achieved.implied")
  )
  for (field in c("mean", "variance")) {
    .dpprior_schema_validate_finite_scalar(
      implied[[field, exact = TRUE]],
      paste0(path, ".achieved.implied.", field),
      lower = if (identical(field, "mean")) 1 else 0,
      upper = if (identical(field, "mean")) raw[["J", exact = TRUE]] else Inf
    )
  }
  recomputed <- .dpprior_target_pmf_moments(pmf)
  moment_tolerance <- raw[["tolerances", exact = TRUE]][[
    "moment_relative", exact = TRUE
  ]] * pmax(
    raw[["tolerances", exact = TRUE]][["moment_scale_floor", exact = TRUE]],
    abs(unlist(implied, use.names = FALSE)), abs(recomputed)
  )
  .dpprior_schema_require(
    all(abs(unlist(implied, use.names = FALSE) - recomputed) <=
          moment_tolerance),
    "target_snapshot_moments", paste0(path, ".achieved.implied"),
    "mean and variance independently recomputed from the snapshot PMF",
    list(recorded = implied, recomputed = recomputed)
  )
  if (identical(raw[["kind", exact = TRUE]], "interval")) {
    interval_achieved <- achieved[["interval", exact = TRUE]]
    .dpprior_schema_exact_names(
      interval_achieved, c("left_mass", "inside_mass", "right_mass"),
      paste0(path, ".achieved.interval")
    )
    for (field in c("left_mass", "inside_mass", "right_mass")) {
      .dpprior_schema_validate_finite_scalar(
        interval_achieved[[field, exact = TRUE]],
        paste0(path, ".achieved.interval.", field), lower = 0, upper = 1
      )
    }
    recomputed_mass <- .dpprior_target_interval_masses(
      pmf, raw[["interval", exact = TRUE]]
    )
    .dpprior_schema_require(
      all(abs(unlist(interval_achieved, use.names = FALSE) - recomputed_mass) <=
            raw[["tolerances", exact = TRUE]][["constraint", exact = TRUE]]),
      "target_snapshot_interval", paste0(path, ".achieved.interval"),
      "left/inside/right masses recomputed from the snapshot PMF",
      list(recorded = interval_achieved, recomputed = recomputed_mass)
    )
  } else {
    .dpprior_schema_require(
      is.null(achieved[["interval", exact = TRUE]]),
      "target_snapshot_interval", paste0(path, ".achieved.interval"),
      "NULL outside interval target construction",
      achieved[["interval", exact = TRUE]]
    )
  }
  invisible(TRUE)
}


.dpprior_validate_constructed_target_verification <- function(raw) {
  verification <- raw[["verification", exact = TRUE]]
  selected <- verification[["selected_snapshot", exact = TRUE]]
  verifier <- verification[["verifier_snapshot", exact = TRUE]]
  stability <- verification[["stability", exact = TRUE]]
  truth <- raw[["tolerances", exact = TRUE]]
  controls <- raw[["computation", exact = TRUE]][["used", exact = TRUE]][[
    "controls", exact = TRUE
  ]]
  expected_settings <- list(
    pmf_l1_tolerance = truth[["pmf_l1", exact = TRUE]],
    moment_relative_tolerance = truth[["moment_relative", exact = TRUE]],
    moment_scale_floor = truth[["moment_scale_floor", exact = TRUE]],
    constraint_tolerance = truth[["constraint", exact = TRUE]],
    root_tolerance = controls[["root_tolerance", exact = TRUE]],
    max_iterations = controls[["max_iterations", exact = TRUE]]
  )
  .dpprior_schema_require(
    identical(verification[["method", exact = TRUE]],
              "independent_target_reconstruction") &&
      verification[["performed", exact = TRUE]] &&
      verification[["passed", exact = TRUE]] &&
      identical(verification[["reason", exact = TRUE]], "verified") &&
      identical(verification[["settings", exact = TRUE]], expected_settings) &&
      identical(selected[["source", exact = TRUE]], "target_constructor") &&
      identical(verifier[["source", exact = TRUE]],
                "independent_target_reconstruction") &&
      identical(stability[["source", exact = TRUE]],
                "independent_target_reconstruction") &&
      identical(names(verification[["components", exact = TRUE]]),
                c("target_reconstruction", "order_stability")) &&
      identical(names(verification[["invariants", exact = TRUE]]),
                c("support_identity", "authority_identity")),
    "target_verification_contract", "target.verification",
    paste(
      "the closed independent-reconstruction method, sources, settings,",
      "component vocabulary, and invariant vocabulary"
    ), verification
  )
  selected_pmf <- selected[["achieved", exact = TRUE]][["pmf", exact = TRUE]]
  verifier_pmf <- verifier[["achieved", exact = TRUE]][["pmf", exact = TRUE]]
  reconstruction_value <- sum(abs(
    verifier_pmf - raw[["pmf", exact = TRUE]]
  ))
  .dpprior_bind_decision_check(
    verification[["components", exact = TRUE]][[
      "target_reconstruction", exact = TRUE
    ]], reconstruction_value, 0, truth[["pmf_l1", exact = TRUE]], "lte",
    "target.verification.components.target_reconstruction"
  )
  .dpprior_bind_decision_check(
    verification[["components", exact = TRUE]][[
      "order_stability", exact = TRUE
    ]], stability[["delta", exact = TRUE]],
    setNames(
      rep(0, length(stability[["delta", exact = TRUE]])),
      names(stability[["delta", exact = TRUE]])
    ), stability[["tolerance", exact = TRUE]], "lte",
    "target.verification.components.order_stability"
  )
  support_identity <- c(
    target_support = identical(
      raw[["support", exact = TRUE]], seq_len(raw[["J", exact = TRUE]])
    ),
    selected_length = length(selected_pmf) == raw[["J", exact = TRUE]],
    verifier_length = length(verifier_pmf) == raw[["J", exact = TRUE]]
  )
  authority_identity <- c(
    selected_pmf = identical(selected_pmf, raw[["pmf", exact = TRUE]]),
    used_pmf = identical(
      raw[["used", exact = TRUE]][["pmf", exact = TRUE]],
      raw[["pmf", exact = TRUE]]
    ),
    normalized_pmf_unidentified = is.null(
      raw[["normalized", exact = TRUE]][["pmf", exact = TRUE]]
    ),
    request_J = identical(
      raw[["request", exact = TRUE]][["J", exact = TRUE]],
      raw[["J", exact = TRUE]]
    )
  )
  .dpprior_bind_decision_check(
    verification[["invariants", exact = TRUE]][[
      "support_identity", exact = TRUE
    ]], support_identity,
    setNames(rep(TRUE, length(support_identity)), names(support_identity)),
    NULL, "identical", "target.verification.invariants.support_identity"
  )
  .dpprior_bind_decision_check(
    verification[["invariants", exact = TRUE]][[
      "authority_identity", exact = TRUE
    ]], authority_identity,
    setNames(rep(TRUE, length(authority_identity)), names(authority_identity)),
    NULL, "identical", "target.verification.invariants.authority_identity"
  )
  all_checks <- c(
    verification[["components", exact = TRUE]],
    verification[["invariants", exact = TRUE]]
  )
  .dpprior_schema_require(
    all(vapply(
      all_checks,
      function(check) identical(
        check[["source", exact = TRUE]], "independent_target_reconstruction"
      ), logical(1)
    )),
    "target_verification_source", "target.verification",
    "every truth-making check sourced to independent target reconstruction",
    vapply(all_checks, function(check) check[["source", exact = TRUE]],
           character(1))
  )
  invisible(TRUE)
}


.dpprior_validate_target_infeasibility_certificate <- function(raw) {
  .dpprior_schema_require(
    identical(raw[["kind", exact = TRUE]], "interval") &&
      identical(raw[["status", exact = TRUE]], "infeasible") &&
      !raw[["usable", exact = TRUE]] && raw[["verified", exact = TRUE]] &&
      is.null(raw[["pmf", exact = TRUE]]) &&
      is.null(raw[["implied", exact = TRUE]]) &&
      is.null(raw[["achieved_interval", exact = TRUE]]) &&
      is.null(raw[["parameters", exact = TRUE]]) &&
      identical(raw[["normalized", exact = TRUE]], raw[["used", exact = TRUE]]),
    "target_infeasibility_shape", "target",
    paste(
      "an interval-only certified-infeasible target with no PMF, achieved",
      "claim, finite parameters, or fictional normalized-to-used change"
    ), raw[c(
      "kind", "status", "usable", "verified", "pmf", "implied",
      "achieved_interval", "parameters"
    )]
  )
  interval <- raw[["interval", exact = TRUE]]
  family <- raw[["family", exact = TRUE]]
  .dpprior_schema_exact_names(
    interval,
    c(
      "lower", "upper", "type", "coverage", "family", "mu_K",
      "support", "endpoints"
    ), "target.interval"
  )
  .dpprior_schema_exact_names(
    family,
    c("name", "parameterization", "explicit", "reference_measure"),
    "target.family"
  )
  .dpprior_schema_require(
    identical(family, list(
      name = "maxent",
      parameterization = "common exponential tilt within fixed-mass groups",
      explicit = TRUE,
      reference_measure = "counting measure on 1:J"
    )) && identical(interval[["family", exact = TRUE]], "maxent") &&
      identical(raw[["assumptions", exact = TRUE]], list(
        estimand = "K_J", construction = "maxent"
      )) &&
      identical(interval[["support", exact = TRUE]],
                c(lower = 1L, upper = as.integer(raw[["J", exact = TRUE]]))) &&
      identical(interval[["endpoints", exact = TRUE]], "inclusive"),
    "target_infeasibility_family", "target.family",
    "the closed MaxEnt interval family and exact 1:J support", family
  )
  support <- seq_len(raw[["J", exact = TRUE]])
  lower <- interval[["lower", exact = TRUE]]
  upper <- interval[["upper", exact = TRUE]]
  coverage <- interval[["coverage", exact = TRUE]]
  .dpprior_schema_require(
    .dpprior_schema_is_count(lower, 1L) &&
      .dpprior_schema_is_count(upper, lower) && upper <= raw[["J", exact = TRUE]],
    "target_infeasibility_interval", "target.interval",
    "integer inclusive endpoints within 1:J", interval
  )
  .dpprior_schema_validate_finite_scalar(
    coverage, "target.interval.coverage", lower = 0, upper = 1
  )
  groups <- switch(
    interval[["type", exact = TRUE]],
    hard_bounds = list(inside = support[support >= lower & support <= upper]),
    equal_tail = list(
      left = support[support < lower],
      inside = support[support >= lower & support <= upper],
      right = support[support > upper]
    ),
    central_mass = list(
      inside = support[support >= lower & support <= upper],
      outside = support[support < lower | support > upper]
    ),
    NULL
  )
  .dpprior_schema_require(
    !is.null(groups), "target_infeasibility_interval", "target.interval.type",
    "hard_bounds, equal_tail, or central_mass", interval[["type", exact = TRUE]]
  )
  group_masses <- switch(
    interval[["type", exact = TRUE]],
    hard_bounds = c(inside = 1),
    equal_tail = c(
      left = (1 - coverage) / 2,
      inside = coverage,
      right = (1 - coverage) / 2
    ),
    central_mass = c(inside = coverage, outside = 1 - coverage)
  )
  group_counts <- vapply(groups, length, integer(1))
  empty_groups <- names(groups)[group_masses > 0 & group_counts == 0L]
  positive_groups <- names(groups)[group_masses > 0 & group_counts > 0L]
  requested_mean <- interval[["mu_K", exact = TRUE]]
  if (!is.null(requested_mean)) {
    .dpprior_schema_validate_finite_scalar(
      requested_mean, "target.interval.mu_K", lower = 1,
      upper = raw[["J", exact = TRUE]
      ]
    )
  }
  lower_bound <- if (!length(empty_groups)) sum(vapply(
    positive_groups,
    function(group) group_masses[[group]] * min(groups[[group]]),
    numeric(1)
  )) else NULL
  upper_bound <- if (!length(empty_groups)) sum(vapply(
    positive_groups,
    function(group) group_masses[[group]] * max(groups[[group]]),
    numeric(1)
  )) else NULL
  feasibility_tolerance <- 64 * .Machine$double.eps * max(
    1,
    if (is.null(lower_bound)) 0 else abs(lower_bound),
    if (is.null(upper_bound)) 0 else abs(upper_bound),
    if (is.null(requested_mean)) 0 else abs(requested_mean)
  )
  outside_side <- if (!is.null(requested_mean) && !is.null(lower_bound) &&
                      requested_mean < lower_bound - feasibility_tolerance) {
    "below"
  } else if (!is.null(requested_mean) && !is.null(upper_bound) &&
             requested_mean > upper_bound + feasibility_tolerance) {
    "above"
  } else {
    NULL
  }
  outside_distance <- if (identical(outside_side, "below")) {
    lower_bound - requested_mean
  } else if (identical(outside_side, "above")) {
    requested_mean - upper_bound
  } else {
    NULL
  }
  expected_kind <- if (length(empty_groups)) {
    "positive_mass_group_has_empty_support"
  } else {
    "mean_outside_group_mass_hull"
  }
  certificate <- raw[["verification", exact = TRUE]][[
    "settings", exact = TRUE
  ]][["certificate", exact = TRUE]]
  certificate_fields <- c(
    "version", "method", "kind", "assumptions", "request", "J", "support",
    "interval", "family", "group_masses", "group_support_counts",
    "empty_groups", "requested_mean", "feasible_mean_lower",
    "feasible_mean_upper", "feasibility_tolerance", "side",
    "outside_distance", "certified", "source"
  )
  .dpprior_schema_exact_names(
    certificate, certificate_fields,
    "target.verification.settings.certificate"
  )
  if (is.null(certificate[["empty_groups", exact = TRUE]])) {
    .dpprior_schema_require(
      length(empty_groups) == 0L, "target_infeasibility_empty_groups",
      "target.verification.settings.certificate.empty_groups",
      "NULL exactly when no positive-mass support group is empty",
      certificate[["empty_groups", exact = TRUE]]
    )
  } else {
    .dpprior_schema_validate_character_vector(
      certificate[["empty_groups", exact = TRUE]],
      "target.verification.settings.certificate.empty_groups"
    )
  }
  .dpprior_schema_require(
    identical(certificate[["version", exact = TRUE]], "1") &&
      identical(certificate[["method", exact = TRUE]],
                "analytic_interval_group_support_feasibility") &&
      identical(certificate[["kind", exact = TRUE]], expected_kind) &&
      identical(certificate[["assumptions", exact = TRUE]],
                raw[["assumptions", exact = TRUE]]) &&
      identical(certificate[["request", exact = TRUE]],
                raw[["request", exact = TRUE]]) &&
      identical(certificate[["J", exact = TRUE]], raw[["J", exact = TRUE]]) &&
      identical(certificate[["support", exact = TRUE]],
                raw[["support", exact = TRUE]]) &&
      identical(certificate[["interval", exact = TRUE]], interval) &&
      identical(certificate[["family", exact = TRUE]], family) &&
      identical(certificate[["group_masses", exact = TRUE]],
                as.list(group_masses)) &&
      identical(certificate[["group_support_counts", exact = TRUE]],
                as.list(group_counts)) &&
      identical(
        certificate[["empty_groups", exact = TRUE]],
        if (length(empty_groups)) empty_groups else NULL
      ) &&
      identical(certificate[["requested_mean", exact = TRUE]],
                requested_mean) &&
      identical(certificate[["feasible_mean_lower", exact = TRUE]],
                lower_bound) &&
      identical(certificate[["feasible_mean_upper", exact = TRUE]],
                upper_bound) &&
      identical(certificate[["feasibility_tolerance", exact = TRUE]],
                feasibility_tolerance) &&
      identical(certificate[["side", exact = TRUE]], outside_side) &&
      identical(certificate[["outside_distance", exact = TRUE]],
                outside_distance) &&
      identical(certificate[["certified", exact = TRUE]], TRUE) &&
      identical(certificate[["source", exact = TRUE]],
                "analytic_maxent_interval_feasibility") &&
      (if (identical(expected_kind,
                     "positive_mass_group_has_empty_support")) {
        length(empty_groups) > 0L
      } else {
        !is.null(outside_distance) &&
          outside_distance > feasibility_tolerance
      }),
    "target_infeasibility_certificate",
    "target.verification.settings.certificate",
    paste(
      "a request-bound analytic group-support certificate whose empty-group",
      "or mean-hull contradiction is independently recomputed"
    ), certificate
  )
  verification <- raw[["verification", exact = TRUE]]
  .dpprior_schema_require(
    identical(verification[["method", exact = TRUE]],
              "analytic_interval_infeasibility_certificate") &&
      verification[["performed", exact = TRUE]] &&
      verification[["passed", exact = TRUE]] &&
      identical(verification[["reason", exact = TRUE]],
                "certified_infeasible") &&
      identical(names(verification[["settings", exact = TRUE]]),
                "certificate") &&
      is.null(verification[["selected_snapshot", exact = TRUE]]) &&
      is.null(verification[["verifier_snapshot", exact = TRUE]]) &&
      is.null(verification[["stability", exact = TRUE]]) &&
      identical(names(verification[["components", exact = TRUE]]),
                "infeasibility_certificate") &&
      identical(names(verification[["invariants", exact = TRUE]]),
                c("request_identity", "support_identity")),
    "target_infeasibility_verification", "target.verification",
    "the exact analytic certificate-only target verification contract",
    verification
  )
  component_value <- if (length(empty_groups)) {
    as.integer(length(empty_groups))
  } else {
    outside_distance
  }
  component_tolerance <- if (length(empty_groups)) 0 else
    feasibility_tolerance
  .dpprior_bind_decision_check(
    verification[["components", exact = TRUE]][[
      "infeasibility_certificate", exact = TRUE
    ]], component_value, 0, component_tolerance, "gt",
    "target.verification.components.infeasibility_certificate"
  )
  request_identity <- c(
    J = identical(raw[["request", exact = TRUE]][["J", exact = TRUE]],
                  raw[["J", exact = TRUE]]),
    interval = identical(
      raw[["request", exact = TRUE]][["K_interval", exact = TRUE]], interval
    ),
    mean = identical(
      raw[["request", exact = TRUE]][["mu_K", exact = TRUE]], requested_mean
    )
  )
  support_identity <- c(
    top = identical(raw[["support", exact = TRUE]], support),
    interval = identical(interval[["support", exact = TRUE]],
                         c(lower = 1L, upper = as.integer(raw[["J", exact = TRUE]])))
  )
  .dpprior_bind_decision_check(
    verification[["invariants", exact = TRUE]][[
      "request_identity", exact = TRUE
    ]], request_identity,
    setNames(rep(TRUE, length(request_identity)), names(request_identity)),
    NULL, "identical", "target.verification.invariants.request_identity"
  )
  .dpprior_bind_decision_check(
    verification[["invariants", exact = TRUE]][[
      "support_identity", exact = TRUE
    ]], support_identity,
    setNames(rep(TRUE, length(support_identity)), names(support_identity)),
    NULL, "identical", "target.verification.invariants.support_identity"
  )
  expected_source <- "analytic_maxent_interval_feasibility"
  all_checks <- c(
    verification[["components", exact = TRUE]],
    verification[["invariants", exact = TRUE]]
  )
  .dpprior_schema_require(
    all(vapply(
      all_checks,
      function(check) identical(check[["source", exact = TRUE]],
                                  expected_source), logical(1)
    )),
    "target_infeasibility_source", "target.verification",
    "the fixed analytic MaxEnt feasibility source on every proof check",
    vapply(all_checks, function(check) check[["source", exact = TRUE]],
           character(1))
  )
  computation <- raw[["computation", exact = TRUE]]
  attempts <- computation[["attempts", exact = TRUE]]
  .dpprior_schema_require(
    identical(computation[["request", exact = TRUE]][["method", exact = TRUE]],
              "target_interval_maxent") &&
      identical(computation[["used", exact = TRUE]][["method", exact = TRUE]],
                "target_interval_maxent") &&
      length(attempts) == 1L &&
      is.null(computation[["selected_attempt_id", exact = TRUE]]) &&
      is.null(computation[["selected_candidate_id", exact = TRUE]]) &&
      length(computation[["candidate_evaluations", exact = TRUE]]) == 0L &&
      !computation[["fallback", exact = TRUE]][["attempted", exact = TRUE]] &&
      !computation[["fallback", exact = TRUE]][["used", exact = TRUE]] &&
      identical(computation[["termination", exact = TRUE]][[
        "code", exact = TRUE
      ]], "certified_infeasible") &&
      identical(computation[["termination", exact = TRUE]][[
        "source", exact = TRUE
      ]], "analytic_certificate") &&
      is.null(computation[["termination", exact = TRUE]][[
        "iterations", exact = TRUE
      ]]),
    "target_infeasibility_computation", "target.computation",
    paste(
      "one typed analytic feasibility attempt, no selected candidate or",
      "fallback, and certified_infeasible/analytic_certificate termination"
    ), computation
  )
  attempt <- attempts[[1L]]
  .dpprior_schema_require(
      identical(attempt[["stage", exact = TRUE]], "feasibility") &&
      identical(attempt[["method", exact = TRUE]],
                "analytic_interval_group_support_feasibility") &&
      identical(attempt[["start", exact = TRUE]], list(
        J = raw[["J", exact = TRUE]], requested_mean = requested_mean,
        interval_lower = lower, interval_upper = upper,
        coverage = coverage
      )) &&
      identical(attempt[["bounds", exact = TRUE]], list(
        feasible_mean_lower = lower_bound,
        feasible_mean_upper = upper_bound
      )) &&
      identical(attempt[["control", exact = TRUE]], list(
        feasibility_tolerance = feasibility_tolerance
      )) && identical(attempt[["exit_code", exact = TRUE]], 0L) &&
      identical(attempt[["iterations", exact = TRUE]], 0L) &&
      identical(attempt[["evaluations", exact = TRUE]],
                list(function_count = 1L)) &&
      identical(attempt[["message", exact = TRUE]],
                "analytic interval infeasibility certified") &&
      is.null(attempt[["error", exact = TRUE]]) &&
      identical(attempt[["warnings", exact = TRUE]], character()) &&
      is.null(attempt[["candidate_parameters", exact = TRUE]]) &&
      is.null(attempt[["candidate_objective", exact = TRUE]]) &&
      !attempt[["selected", exact = TRUE]] &&
      identical(attempt[["reason_code", exact = TRUE]],
                "globally_infeasible_by_certificate") &&
      identical(names(attempt[["unavailable", exact = TRUE]]),
                c("candidate_parameters", "candidate_objective")),
    "target_infeasibility_attempt", "target.computation.attempts[[1]]",
    "the exact request/bounds/tolerance-bound analytic failure attempt",
    attempt
  )
  invisible(TRUE)
}


.dpprior_validate_target_v1_impl <- function(x) {
  .dpprior_schema_require(
    typeof(x) == "list" && is.list(x),
    "type", "target", "an ordinary list", typeof(x)
  )
  raw <- unclass(x)
  .dpprior_schema_validate_attributes(raw, "target", "names", "names")
  .dpprior_schema_require(
    identical(class(x), c("dpprior_K_target", "dpprior_target", "list")),
    "class", "target", "the exact canonical K-target class vector", class(x)
  )
  .dpprior_schema_require(
    !anyDuplicated(names(raw)), "duplicate_names", "target",
    "unique top-level names", names(raw)
  )
  .dpprior_schema_require(
    length(raw) >= length(.DPPRIOR_TARGET_FIELDS) &&
      identical(
        names(raw)[seq_along(.DPPRIOR_TARGET_FIELDS)],
        .DPPRIOR_TARGET_FIELDS
      ),
    "field_names", "target",
    "the exact 25-field target spine in construction order", names(raw)
  )
  .dpprior_validate_compatibility(
    raw[["compatibility", exact = TRUE]], "target.compatibility"
  )
  target_aliases <- raw[["compatibility", exact = TRUE]][[
    "top_level_aliases", exact = TRUE
  ]]
  target_alias_names <- names(target_aliases)
  if (is.null(target_alias_names)) {
    target_alias_names <- character()
  }
  target_tail <- names(raw)[-(seq_along(.DPPRIOR_TARGET_FIELDS))]
  .dpprior_schema_require(
    identical(target_tail, target_alias_names), "alias_order", "target",
    "only registered compatibility aliases after the target spine", target_tail
  )
  .dpprior_validate_schema_record(
    raw$schema, .DPPRIOR_TARGET_SCHEMA_NAME, "target.schema"
  )
  .dpprior_schema_validate_scalar_character(raw$kind, "target.kind")
  .dpprior_schema_require(
    raw$kind %in% c("moments", "cv", "interval", "pmf", "family"),
    "target_kind", "target.kind",
    "moments, cv, interval, pmf, or family", raw$kind
  )
  .dpprior_schema_require(
    .dpprior_schema_is_count(raw$J, 1L), "count", "target.J",
    "an integer at least 1", raw$J
  )
  .dpprior_schema_require(
    is.integer(raw$support) && !is.object(raw$support) &&
      is.null(dim(raw$support)) && identical(raw$support, seq_len(raw$J)),
    "support", "target.support", "the exact integer support 1:J", raw$support
  )
  for (field in c("request", "normalized", "used", "derivation",
                  "assumptions", "residuals", "tolerances")) {
    .dpprior_schema_validate_named_list(
      raw[[field, exact = TRUE]], paste0("target.", field)
    )
    .dpprior_schema_validate_plain_record_value(
      raw[[field, exact = TRUE]], paste0("target.", field)
    )
  }
  for (field in c("interval", "family", "implied", "achieved_interval")) {
    if (!is.null(raw[[field]])) {
      .dpprior_schema_validate_named_list(
        raw[[field, exact = TRUE]], paste0("target.", field)
      )
      .dpprior_schema_validate_plain_record_value(
        raw[[field, exact = TRUE]], paste0("target.", field)
      )
    }
  }
  has_interval <- !is.null(raw[["interval", exact = TRUE]])
  has_family <- !is.null(raw[["family", exact = TRUE]])
  kind_fields_ok <- switch(
    raw[["kind", exact = TRUE]],
    moments = !has_interval && !has_family,
    cv = !has_interval && !has_family,
    pmf = !has_interval && !has_family,
    interval = has_interval,
    family = has_family && !has_interval
  )
  .dpprior_schema_require(
    kind_fields_ok,
    "target_kind_fields", "target",
    "kind-specific interval/family presence and absence", list(
      kind = raw[["kind", exact = TRUE]],
      interval = has_interval,
      family = has_family
    )
  )
  for (record_name in c("request", "normalized", "used")) {
    record <- raw[[record_name]]
    .dpprior_schema_require(
      "J" %in% names(record), "target_J",
      paste0("target.", record_name), "an exact J field", names(record)
    )
    record_J <- record[["J", exact = TRUE]]
    .dpprior_schema_require(
      .dpprior_schema_is_count(record_J, 1L),
      "target_J", paste0("target.", record_name, ".J"),
      "one finite integer-valued scalar", record_J
    )
    .dpprior_schema_require(
      identical(as.integer(record_J), as.integer(raw$J)),
      "target_J", paste0("target.", record_name, ".J"),
      "the same J as target.J", record_J
    )
  }

  target_moments <- function(record, record_path, required = FALSE) {
    has_canonical <- all(c("mean", "variance") %in% names(record))
    has_legacy_spelling <- all(c("mu_K", "var_K") %in% names(record))
    .dpprior_schema_require(
      !(has_canonical && has_legacy_spelling),
      "target_moments", record_path,
      "only one of mean/variance or mu_K/var_K", names(record)
    )
    if (!has_canonical && !has_legacy_spelling) {
      .dpprior_schema_require(
        !required, "target_moments", record_path,
        "an authoritative mean/variance pair", names(record)
      )
      return(NULL)
    }
    fields <- if (has_canonical) c("mean", "variance") else c("mu_K", "var_K")
    mean_value <- record[[fields[[1L]], exact = TRUE]]
    variance_value <- record[[fields[[2L]], exact = TRUE]]
    .dpprior_schema_validate_finite_scalar(
      mean_value, paste0(record_path, ".", fields[[1L]]),
      lower = 1, upper = raw$J
    )
    .dpprior_schema_validate_finite_scalar(
      variance_value, paste0(record_path, ".", fields[[2L]]), lower = 0
    )
    bound <- (mean_value - 1) * (raw$J - mean_value)
    tolerance <- 1e-8 * max(1, abs(bound), abs(variance_value))
    .dpprior_schema_require(
      variance_value <= bound + tolerance,
      "target_variance_bound", paste0(record_path, ".", fields[[2L]]),
      "variance <= (mean-1)*(J-mean)", variance_value
    )
    c(mean = mean_value, variance = variance_value)
  }

  used_moments <- target_moments(
    raw$used, "target.used", required = raw$kind %in% c("moments", "cv")
  )
  implied_moments <- if (is.null(raw$implied)) {
    NULL
  } else {
    target_moments(raw$implied, "target.implied", required = TRUE)
  }
  normalized <- raw[["normalized", exact = TRUE]]
  used <- raw[["used", exact = TRUE]]
  for (field in c("interval", "pmf")) {
    .dpprior_schema_require(
      field %in% names(used) &&
        identical(used[[field, exact = TRUE]], raw[[field, exact = TRUE]]),
      "target_authority", paste0("target.used.", field),
      paste0("identity with target.", field), used[[field, exact = TRUE]]
    )
  }
  .dpprior_schema_require(
    "interval" %in% names(normalized) &&
      identical(normalized[["interval", exact = TRUE]],
                raw[["interval", exact = TRUE]]),
    "target_authority", "target.normalized.interval",
    "identity with target.interval", normalized[["interval", exact = TRUE]]
  )
  if (!is.null(raw[["family", exact = TRUE]])) {
    for (record_name in c("normalized", "used")) {
      record <- raw[[record_name, exact = TRUE]]
      .dpprior_schema_require(
        "family" %in% names(record) &&
          identical(record[["family", exact = TRUE]],
                    raw[["family", exact = TRUE]]),
        "target_authority", paste0("target.", record_name, ".family"),
        "identity with target.family", record[["family", exact = TRUE]]
      )
    }
  }
  normalized_pmf <- normalized[["pmf", exact = TRUE]]
  .dpprior_schema_require(
    if (identical(raw[["kind", exact = TRUE]], "pmf")) {
      identical(normalized_pmf, raw[["pmf", exact = TRUE]])
    } else {
      is.null(normalized_pmf)
    },
    "target_normalized_pmf", "target.normalized.pmf",
    paste(
      "the strict input PMF for kind=pmf, otherwise NULL before any",
      "interval/family construction"
    ), normalized_pmf
  )

  if (!is.null(raw$pmf)) {
    .dpprior_schema_require(
      is.numeric(raw$pmf) && !is.object(raw$pmf) && is.null(dim(raw$pmf)) &&
        .dpprior_schema_has_only_attributes(raw$pmf) &&
        length(raw$pmf) == raw$J && !anyNA(raw$pmf) &&
        all(is.finite(raw$pmf)) && all(raw$pmf >= 0),
      "pmf", "target.pmf",
      "a finite nonnegative numeric vector of length J", raw$pmf
    )
    .dpprior_schema_require(
      abs(sum(raw$pmf) - 1) <= .TOL_PMF_SUM,
      "pmf_normalization", "target.pmf",
      sprintf("mass 1 within %.3g without repair", .TOL_PMF_SUM), sum(raw$pmf)
    )
    .dpprior_schema_require(
      identical(used[["pmf", exact = TRUE]], raw$pmf),
      "target_pmf_identity", "target.used.pmf",
      "identity with target.pmf", used[["pmf", exact = TRUE]]
    )
    support <- seq_len(raw$J)
    pmf_moments <- c(
      mean = sum(support * raw$pmf),
      variance = {
        pmf_mean <- sum(support * raw$pmf)
        sum((support - pmf_mean)^2 * raw$pmf)
      }
    )
    .dpprior_schema_require(
      !is.null(implied_moments), "target_implied", "target.implied",
      "implied mean/variance for an identified PMF", raw$implied
    )
    moment_tolerance <- 1e-8 * pmax(1, abs(pmf_moments), abs(implied_moments))
    .dpprior_schema_require(
      all(abs(implied_moments - pmf_moments) <= moment_tolerance),
      "target_pmf_moments", "target.implied",
      "moments recomputed from target.pmf",
      list(implied = implied_moments, recomputed = pmf_moments)
    )
  }
  if (raw$kind %in% c("moments", "cv")) {
    .dpprior_schema_require(
      is.null(raw$pmf), "unidentified_pmf", "target.pmf",
      "NULL when moments/CV do not identify a distribution", raw$pmf
    )
    if (!is.null(implied_moments)) {
      moment_tolerance <- 1e-8 * pmax(
        1, abs(used_moments), abs(implied_moments)
      )
      .dpprior_schema_require(
        all(abs(used_moments - implied_moments) <= moment_tolerance),
        "target_moment_identity", "target.implied",
        "agreement with target.used moments",
        list(used = used_moments, implied = implied_moments)
      )
    }
  }
  if (identical(raw$kind, "pmf")) {
    .dpprior_schema_require(
      !is.null(raw$pmf), "target_pmf", "target.pmf",
      "an identified PMF for kind=pmf", NULL
    )
  }
  if (raw$kind %in% c("family", "interval") &&
      raw$status %in% c("converged", "boundary")) {
    .dpprior_schema_require(
      !is.null(raw$pmf), "target_pmf", "target.pmf",
      "an identified PMF for a successful family/interval target", NULL
    )
  }

  status_record <- list(
    status = raw$status,
    usable = raw$usable,
    verified = raw$verified,
    message = raw$message
  )
  .dpprior_validate_status_record(status_record, "target.status_record")
  .dpprior_validate_parameters(raw$parameters, "target.parameters", TRUE)
  .dpprior_validate_computation(raw$computation, "target.computation")
  .dpprior_validate_verification(raw$verification, "target.verification")
  .dpprior_validate_provenance(raw$provenance, "target.provenance")
  .dpprior_validate_compatibility(raw$compatibility, "target.compatibility")
  .dpprior_validate_target_truth_authority(raw)
  .dpprior_validate_target_derivation(raw)
  projection_entry <- raw[["derivation", exact = TRUE]][[
    "normalized_to_used", exact = TRUE
  ]]
  projection_changed <- !is.null(projection_entry) && identical(
    projection_entry[["rule", exact = TRUE]],
    "project_a1_variance_to_nearest_interior"
  ) && !identical(
    raw[["normalized", exact = TRUE]], raw[["used", exact = TRUE]]
  ) && isTRUE(projection_entry[["opt_in", exact = TRUE]])
  target_projection <- raw[["provenance", exact = TRUE]][[
    "projection", exact = TRUE
  ]]
  .dpprior_schema_require(
    identical(target_projection[["applied", exact = TRUE]],
              projection_changed) && if (!projection_changed) {
      identical(target_projection, list(
        applied = FALSE, opt_in = FALSE, policy = NULL, record = NULL
      ))
    } else {
      TRUE
    },
    "target_projection_iff", "target.provenance.projection",
    paste(
      "applied exactly for a changed, opt-in A1 nearest projection and the",
      "exact inactive record for every non-projection target"
    ), target_projection
  )

  .dpprior_schema_require(
    identical(raw$computation$request$method,
              raw$provenance$requested_method),
    "requested_method", "target.provenance.requested_method",
    "identity with target.computation.request.method",
    raw$provenance$requested_method
  )
  .dpprior_schema_require(
    identical(raw$computation$used$method, raw$provenance$selected_method),
    "selected_method", "target.provenance.selected_method",
    "identity with target.computation.used.method",
    raw$provenance$selected_method
  )
  .dpprior_schema_require(
    identical(raw$computation$fallback$used, raw$provenance$is_fallback),
    "fallback", "target.provenance.is_fallback",
    "identity with target.computation.fallback.used",
    raw$provenance$is_fallback
  )
  .dpprior_schema_require(
    identical(
      raw[["verified", exact = TRUE]],
      raw[["verification", exact = TRUE]][["passed", exact = TRUE]]
    ),
    "verification", "target.verified",
    "identity with target.verification.passed",
    raw[["verified", exact = TRUE]]
  )
  selected_snapshot <- raw[["verification", exact = TRUE]][[
    "selected_snapshot", exact = TRUE
  ]]
  verifier_snapshot <- raw[["verification", exact = TRUE]][[
    "verifier_snapshot", exact = TRUE
  ]]
  constructed_target <- raw[["kind", exact = TRUE]] %in% c(
    "interval", "family"
  ) && raw[["status", exact = TRUE]] %in% c("converged", "boundary")
  if (constructed_target) {
    .dpprior_validate_constructed_target_snapshot(
      raw, selected_snapshot, "target.verification.selected_snapshot"
    )
    .dpprior_validate_constructed_target_snapshot(
      raw, verifier_snapshot, "target.verification.verifier_snapshot"
    )
    expected_stability <- .dpprior_expected_target_stability(raw)
    .dpprior_schema_require(
      identical(raw[["verification", exact = TRUE]][[
        "stability", exact = TRUE
      ]], expected_stability) && expected_stability[["passed", exact = TRUE]],
      "target_stability_identity", "target.verification.stability",
      paste(
        "the independently recomputed selected-versus-verifier PMF, moment,",
        "and interval-mass stability record, with every gate passing"
      ),
      raw[["verification", exact = TRUE]][["stability", exact = TRUE]]
    )
    .dpprior_validate_constructed_target_verification(raw)
  }
  if (!is.null(selected_snapshot)) {
    expected_achieved <- list(
      implied = raw[["implied", exact = TRUE]],
      interval = raw[["achieved_interval", exact = TRUE]],
      pmf = raw[["pmf", exact = TRUE]]
    )
    .dpprior_schema_require(
      identical(selected_snapshot[["parameters", exact = TRUE]],
                raw[["parameters", exact = TRUE]]) &&
        identical(selected_snapshot[["M", exact = TRUE]],
                  raw[["computation", exact = TRUE]][[
                    "orders", exact = TRUE
                  ]][["M_selected", exact = TRUE]]) &&
        identical(selected_snapshot[["achieved", exact = TRUE]],
                  expected_achieved) &&
        identical(selected_snapshot[["residuals", exact = TRUE]],
                  raw[["residuals", exact = TRUE]]) &&
        identical(selected_snapshot[["tolerances", exact = TRUE]],
                  raw[["tolerances", exact = TRUE]]),
      "target_snapshot_identity", "target.verification.selected_snapshot",
      paste(
        "identity with target parameters, selected order, achieved target,",
        "residuals, and tolerances"
      ),
      selected_snapshot
    )
  }
  if (!is.null(raw[["parameters", exact = TRUE]])) {
    attempts <- raw[["computation", exact = TRUE]][[
      "attempts", exact = TRUE
    ]]
    selected_id <- raw[["computation", exact = TRUE]][[
      "selected_attempt_id", exact = TRUE
    ]]
    .dpprior_schema_require(
      !is.null(selected_snapshot) && !is.null(selected_id),
      "target_candidate", "target.computation.selected_attempt_id",
      "a selected snapshot and selected attempt for finite parameters", NULL
    )
    ids <- vapply(
      attempts, function(attempt) attempt[["id", exact = TRUE]], character(1)
    )
    selected_attempt <- attempts[[match(selected_id, ids)]]
    .dpprior_schema_require(
      identical(
        selected_attempt[["candidate_parameters", exact = TRUE]],
        raw[["parameters", exact = TRUE]]
      ),
      "target_candidate", "target.parameters",
      "identity with the selected attempt candidate parameters",
      raw[["parameters", exact = TRUE]]
    )
  }
  if (raw[["verified", exact = TRUE]]) {
    .dpprior_schema_require(
      length(raw[["verification", exact = TRUE]][[
        "components", exact = TRUE
      ]]) > 0L &&
        length(raw[["verification", exact = TRUE]][[
          "invariants", exact = TRUE
        ]]) > 0L,
      "verification_evidence", "target.verification",
      "non-empty components and invariants for a verified target",
      raw[["verification", exact = TRUE]]
    )
    if (raw[["kind", exact = TRUE]] %in% c("moments", "cv", "pmf")) {
      .dpprior_schema_require(
        !is.null(verifier_snapshot) &&
          identical(
            verifier_snapshot[["achieved", exact = TRUE]][[
              "implied", exact = TRUE
            ]],
            raw[["implied", exact = TRUE]]
          ),
        "target_verifier_identity",
        "target.verification.verifier_snapshot.achieved.implied",
        "identity with the authoritative target moments",
        if (is.null(verifier_snapshot)) NULL else
          verifier_snapshot[["achieved", exact = TRUE]][[
            "implied", exact = TRUE
          ]]
      )
    }
    if (!is.null(raw[["pmf", exact = TRUE]])) {
      verifier_pmf <- if (is.null(verifier_snapshot)) NULL else
        verifier_snapshot[["achieved", exact = TRUE]][["pmf", exact = TRUE]]
      if (identical(raw[["kind", exact = TRUE]], "pmf")) {
        .dpprior_schema_require(
          identical(verifier_pmf, raw[["pmf", exact = TRUE]]),
          "target_verifier_identity",
          "target.verification.verifier_snapshot.achieved.pmf",
          "identity with the strict authoritative input PMF", verifier_pmf
        )
      } else {
        .dpprior_schema_require(
          is.numeric(verifier_pmf) && !is.object(verifier_pmf) &&
            is.null(dim(verifier_pmf)) &&
            .dpprior_schema_has_only_attributes(verifier_pmf) &&
            length(verifier_pmf) == raw$J &&
            !anyNA(verifier_pmf) && all(is.finite(verifier_pmf)) &&
            all(verifier_pmf >= 0) &&
            abs(sum(verifier_pmf) - 1) <= .TOL_PMF_SUM &&
            !is.null(raw[["verification", exact = TRUE]][[
              "stability", exact = TRUE
            ]]),
          "target_verifier_distribution",
          "target.verification.verifier_snapshot.achieved.pmf",
          paste(
            "an independently reconstructed valid PMF with a typed",
            "selected-versus-verifier stability record"
          ), verifier_pmf
        )
      }
    }
    if (!is.null(raw[["achieved_interval", exact = TRUE]]) &&
        !constructed_target) {
      .dpprior_schema_require(
        !is.null(verifier_snapshot) &&
          identical(
            verifier_snapshot[["achieved", exact = TRUE]][[
              "interval", exact = TRUE
            ]],
            raw[["achieved_interval", exact = TRUE]]
          ),
        "target_verifier_identity",
        "target.verification.verifier_snapshot.achieved.interval",
        "identity with the authoritative achieved interval", NULL
      )
    }
  }
  if (identical(raw[["status", exact = TRUE]], "infeasible")) {
    .dpprior_validate_target_infeasibility_certificate(raw)
  }
  target_orders <- raw[["computation", exact = TRUE]][[
    "orders", exact = TRUE
  ]]
  if (raw[["verified", exact = TRUE]] &&
      !is.null(target_orders[["M_verification_required", exact = TRUE]])) {
    .dpprior_schema_require(
      !is.null(target_orders[["M_verification_used", exact = TRUE]]) &&
        !is.null(verifier_snapshot) &&
        identical(
          verifier_snapshot[["M", exact = TRUE]],
          target_orders[["M_verification_used", exact = TRUE]]
        ),
      "verification_order", "target.computation.orders",
      "a retained verifier snapshot at the used required order",
      target_orders
    )
  }
  target_termination <- raw[["computation", exact = TRUE]][[
    "termination", exact = TRUE
  ]]
  .dpprior_schema_require(
    target_termination[["code", exact = TRUE]] %in%
      .DPPRIOR_TERMINATION_CODES[[raw[["status", exact = TRUE]]]] &&
      if (identical(raw[["status", exact = TRUE]], "boundary")) {
        !is.null(target_termination[["boundary_reason", exact = TRUE]])
      } else {
        is.null(target_termination[["boundary_reason", exact = TRUE]])
      },
    "termination_status", "target.computation.termination",
    "a termination code/reason coherent with target status",
    target_termination
  )
  if (identical(raw$status, "approximate") && raw$usable) {
    .dpprior_schema_require(
      raw$provenance$approximation$active &&
        raw$provenance$approximation$opt_in,
      "approximation_opt_in", "target.usable",
      "explicit approximation opt-in", raw$provenance$approximation
    )
  }
  .dpprior_validate_alias_identity(raw, target_alias_names, path = "target")
  invisible(x)
}


.dpprior_validate_target_v1 <- function(x, collect = FALSE) {
  .dpprior_schema_collect_validation(.dpprior_validate_target_v1_impl, x, collect)
}


.dpprior_new_target_K <- function(kind,
                                  J,
                                  support = NULL,
                                  request,
                                  normalized,
                                  used,
                                  derivation = list(),
                                  interval = NULL,
                                  family = NULL,
                                  assumptions = list(),
                                  pmf = NULL,
                                  implied = NULL,
                                  achieved_interval = NULL,
                                  residuals = list(),
                                  tolerances = list(),
                                  status,
                                  usable,
                                  verified,
                                  message = "",
                                  parameters = NULL,
                                  computation,
                                  verification,
                                  provenance,
                                  compatibility = .dpprior_new_compatibility()) {
  .dpprior_schema_require(
    .dpprior_schema_is_count(J, 1L), "count", "target.J",
    "one ordinary integer-valued scalar at least 1", J
  )
  J <- as.integer(J)
  if (is.null(support)) {
    support <- seq_len(J)
  }
  out <- list(
    schema = .dpprior_schema("target"),
    kind = kind,
    J = J,
    support = support,
    request = request,
    normalized = normalized,
    used = used,
    derivation = derivation,
    interval = interval,
    family = family,
    assumptions = assumptions,
    pmf = pmf,
    implied = implied,
    achieved_interval = achieved_interval,
    residuals = residuals,
    tolerances = tolerances,
    status = status,
    usable = usable,
    verified = verified,
    message = message,
    parameters = parameters,
    computation = computation,
    verification = verification,
    provenance = provenance,
    compatibility = compatibility
  )
  class(out) <- c("dpprior_K_target", "dpprior_target", "list")
  .dpprior_validate_target_v1(out)
  out
}


# --- Mode-specific result extensions ----------------------------------------

.dpprior_validate_proxy_extension <- function(x, path = "proxy") {
  .dpprior_schema_exact_names(
    x, c("mapping", "mapping_verification", "projection", "caveats"), path
  )
  for (field in c("mapping", "mapping_verification", "projection")) {
    .dpprior_schema_validate_named_list(
      x[[field]], paste0(path, ".", field)
    )
    .dpprior_schema_validate_plain_record_value(
      x[[field, exact = TRUE]], paste0(path, ".", field)
    )
  }
  .dpprior_schema_validate_character_vector(x$caveats, paste0(path, ".caveats"))
  invisible(TRUE)
}


.dpprior_new_proxy <- function(mapping,
                               mapping_verification,
                               projection = list(),
                               caveats = character()) {
  out <- list(
    mapping = mapping,
    mapping_verification = mapping_verification,
    projection = projection,
    caveats = caveats
  )
  .dpprior_validate_proxy_extension(out)
  out
}


.DPPRIOR_HARD_OPTIMALITY_FIELDS <- c(
  "performed", "passed", "selection_rule", "K_scales",
  "selected_K_loss", "minimum_K_loss", "tie_tolerance",
  "perturbation_passed", "source", "unavailable_reason"
)

.DPPRIOR_SOFT_OPTIMALITY_FIELDS <- c(
  "performed", "passed", "recorded_objective", "recomputed_objective",
  "objective_difference", "objective_tolerance", "objective_passed",
  "gradient", "gradient_method", "bound_state", "boundary_tolerance",
  "stationarity_operator", "stationarity_tolerance", "component_pass",
  "stationarity_passed",
  "local_base_objective", "neighbor_objectives", "neighbor_tolerances",
  "local_minimum_passed", "start_objective", "candidate_objective",
  "start_tolerance", "no_worse_start", "selection_tolerance", "source",
  "unavailable_reason"
)

.DPPRIOR_INFEASIBILITY_CERTIFICATE_FIELDS <- c(
  "method", "version", "J", "metric", "relation", "target_value",
  "threshold", "probability", "support",
  "domain", "monotonicity", "M_selected", "M_verification",
  "effective_tolerance", "minimum", "maximum", "lower_bound",
  "upper_bound", "tolerance", "certified", "source"
)

.DPPRIOR_CERTIFICATE_CORNER_FIELDS <- c(
  "eta", "parameters", "selected", "refined", "finite", "uncertainty",
  "lower", "upper"
)

.DPPRIOR_CERTIFICATE_EVALUATION_FIELDS <- c(
  "metric", "value", "method", "error_bound", "certified", "M", "source"
)


.dpprior_certificate_expected_evaluation <- function(
    metric, threshold, probability, J, parameters, M, path) {
  estimand <- switch(
    metric,
    wsb_tail = "P(W_SB > threshold)",
    wsb_mean = "E(W_SB)",
    wsb_quantile = "Q_probability(W_SB)",
    wmax_tail_upper = "certified upper bound for P(W_max > threshold)",
    NULL
  )
  .dpprior_schema_require(
    !is.null(estimand) &&
      exists(".dpprior_v2_eval_metric", mode = "function", inherits = TRUE),
    "certificate_evaluator", path,
    "the authoritative Phase 8 metric evaluator for a supported metric",
    metric
  )
  spec <- list(
    metric = metric, estimand = estimand, units = "probability",
    threshold = threshold, probability = probability
  )
  evaluated <- tryCatch(
    .dpprior_v2_eval_metric(
      spec,
      parameters[["a", exact = TRUE]],
      parameters[["b", exact = TRUE]],
      J = J, M = M
    ),
    error = identity
  )
  .dpprior_schema_require(
    !inherits(evaluated, "condition") &&
      typeof(evaluated) == "list" && is.list(evaluated) &&
      .dpprior_schema_is_finite_scalar(evaluated[["value", exact = TRUE]]) &&
      .dpprior_schema_is_finite_scalar(
        evaluated[["error_bound", exact = TRUE]]
      ) &&
      typeof(evaluated[["certification", exact = TRUE]]) == "list" &&
      isTRUE(evaluated[["certification", exact = TRUE]][[
        "certified", exact = TRUE
      ]]),
    "certificate_evaluator", path,
    paste(
      "a finite independently certified value/error result from the",
      "authoritative Phase 8 metric evaluator"
    ),
    if (inherits(evaluated, "condition")) conditionMessage(evaluated) else
      evaluated
  )
  evaluated
}


.dpprior_validate_certificate_evaluation <- function(
    x, metric, threshold, probability, J, parameters, M, path) {
  .dpprior_schema_exact_names(
    x, .DPPRIOR_CERTIFICATE_EVALUATION_FIELDS, path
  )
  for (field in c("metric", "method", "source")) {
    .dpprior_schema_validate_scalar_character(
      x[[field, exact = TRUE]], paste0(path, ".", field)
    )
  }
  .dpprior_schema_require(
    identical(x[["metric", exact = TRUE]], metric) &&
      identical(x[["M", exact = TRUE]], M),
    "certificate_evaluation_identity", path,
    "metric and quadrature order identical to the certificate authority",
    x[c("metric", "M")]
  )
  .dpprior_schema_validate_finite_scalar(
    x[["value", exact = TRUE]], paste0(path, ".value"), lower = 0, upper = 1
  )
  .dpprior_schema_validate_finite_scalar(
    x[["error_bound", exact = TRUE]], paste0(path, ".error_bound"), lower = 0
  )
  .dpprior_schema_validate_scalar_logical(
    x[["certified", exact = TRUE]], paste0(path, ".certified")
  )
  expected_method <- switch(
    metric,
    wsb_tail = "closed_form_beta_gamma_survival",
    wsb_quantile = "closed_form_beta_gamma_quantile",
    wmax_tail_upper = "certified_size_biased_mass_upper_bound",
    wsb_mean = "generalized_gauss_laguerre",
    NA_character_
  )
  .dpprior_schema_require(
    !is.na(expected_method) &&
      identical(x[["method", exact = TRUE]], expected_method) &&
      identical(x[["source", exact = TRUE]],
                "independent_metric_evaluator") &&
      identical(x[["certified", exact = TRUE]],
                !identical(metric, "wsb_mean")),
    "certificate_evaluation_method", path,
    paste(
      "the closed metric evaluator/method and independent certification",
      "policy"
    ),
    x[c("metric", "method", "certified", "source")]
  )
  expected <- .dpprior_certificate_expected_evaluation(
    metric, threshold, probability, J, parameters, M, path
  )
  numeric_tolerance <- 64 * .Machine$double.eps * max(
    1,
    abs(x[["value", exact = TRUE]]),
    abs(expected[["value", exact = TRUE]]),
    abs(x[["error_bound", exact = TRUE]]),
    abs(expected[["error_bound", exact = TRUE]])
  )
  .dpprior_schema_require(
    abs(x[["value", exact = TRUE]] -
          expected[["value", exact = TRUE]]) <= numeric_tolerance &&
      identical(x[["method", exact = TRUE]],
                expected[["method", exact = TRUE]]) &&
      abs(x[["error_bound", exact = TRUE]] -
            expected[["error_bound", exact = TRUE]]) <= numeric_tolerance &&
      identical(
        x[["certified", exact = TRUE]],
        expected[["certification", exact = TRUE]][["certified", exact = TRUE]]
      ),
    "certificate_evaluation_truth", path,
    paste(
      "value, method, error bound, and certification freshly recomputed at",
      "the retained J/a/b/M/threshold/probability"
    ),
    list(recorded = x, expected = expected)
  )
  invisible(TRUE)
}


.dpprior_validate_certificate_corner <- function(
    x, metric, threshold, probability, J, M_selected, M_verification,
    expected_eta, path) {
  .dpprior_schema_exact_names(x, .DPPRIOR_CERTIFICATE_CORNER_FIELDS, path)
  eta <- x[["eta", exact = TRUE]]
  .dpprior_schema_require(
    is.numeric(eta) && !is.object(eta) && is.null(dim(eta)) &&
      .dpprior_schema_has_only_attributes(eta, "names", "names") &&
      identical(names(eta), c("log_a", "log_b")) &&
      length(eta) == 2L && !anyNA(eta) && all(is.finite(eta)) &&
      identical(eta, expected_eta),
    "certificate_corner_eta", paste0(path, ".eta"),
    "the exact analytic monotonicity corner", eta
  )
  .dpprior_schema_exact_names(
    x[["parameters", exact = TRUE]], c("a", "b"),
    paste0(path, ".parameters")
  )
  for (field in c("a", "b")) {
    .dpprior_schema_validate_finite_scalar(
      x[["parameters", exact = TRUE]][[field, exact = TRUE]],
      paste0(path, ".parameters.", field), lower = 0, lower_open = TRUE
    )
  }
  .dpprior_schema_require(
    all(abs(
      unlist(x[["parameters", exact = TRUE]], use.names = FALSE) - exp(eta)
    ) <= 1e-12 * pmax(
      1, unlist(x[["parameters", exact = TRUE]], use.names = FALSE)
    )),
    "certificate_corner_parameters", paste0(path, ".parameters"),
    "a/b exactly derived from exp(log_a/log_b)", x[["parameters", exact = TRUE]]
  )
  .dpprior_validate_certificate_evaluation(
    x[["selected", exact = TRUE]], metric, threshold, probability, J,
    x[["parameters", exact = TRUE]], M_selected,
    paste0(path, ".selected")
  )
  .dpprior_validate_certificate_evaluation(
    x[["refined", exact = TRUE]], metric, threshold, probability, J,
    x[["parameters", exact = TRUE]], M_verification,
    paste0(path, ".refined")
  )
  .dpprior_schema_validate_scalar_logical(
    x[["finite", exact = TRUE]], paste0(path, ".finite")
  )
  for (field in c("uncertainty", "lower", "upper")) {
    .dpprior_schema_validate_finite_scalar(
      x[[field, exact = TRUE]], paste0(path, ".", field),
      lower = 0, upper = if (identical(field, "uncertainty")) Inf else 1
    )
  }
  expected_uncertainty <- max(
    abs(x[["selected", exact = TRUE]][["value", exact = TRUE]] -
          x[["refined", exact = TRUE]][["value", exact = TRUE]]),
    x[["selected", exact = TRUE]][["error_bound", exact = TRUE]],
    x[["refined", exact = TRUE]][["error_bound", exact = TRUE]],
    64 * .Machine$double.eps
  )
  expected_lower <- max(
    0, x[["refined", exact = TRUE]][["value", exact = TRUE]] -
      expected_uncertainty
  )
  expected_upper <- min(
    1, x[["refined", exact = TRUE]][["value", exact = TRUE]] +
      expected_uncertainty
  )
  tolerance <- 1e-12 * max(
    1, expected_uncertainty, expected_lower, expected_upper,
    x[["uncertainty", exact = TRUE]], x[["lower", exact = TRUE]],
    x[["upper", exact = TRUE]]
  )
  .dpprior_schema_require(
    x[["finite", exact = TRUE]] &&
      x[["selected", exact = TRUE]][["certified", exact = TRUE]] &&
      x[["refined", exact = TRUE]][["certified", exact = TRUE]] &&
      abs(x[["uncertainty", exact = TRUE]] - expected_uncertainty) <= tolerance &&
      abs(x[["lower", exact = TRUE]] - expected_lower) <= tolerance &&
      abs(x[["upper", exact = TRUE]] - expected_upper) <= tolerance,
    "certificate_corner_truth", path,
    paste(
      "two independently certified orders with uncertainty=max(order delta,",
      "error bounds, rounding floor) and outward lower/upper bounds"
    ),
    x
  )
  invisible(TRUE)
}


.dpprior_validate_infeasibility_certificate <- function(
    x, path = "constraint.feasibility.certificate") {
  .dpprior_schema_exact_names(
    x, .DPPRIOR_INFEASIBILITY_CERTIFICATE_FIELDS, path
  )
  for (field in c("method", "version", "metric", "relation", "source")) {
    .dpprior_schema_validate_scalar_character(
      x[[field, exact = TRUE]], paste0(path, ".", field)
    )
  }
  .dpprior_schema_require(
    identical(
      x[["method", exact = TRUE]],
      "analytic_global_monotonicity_with_refined_order_corner_enclosures"
    ) && identical(x[["version", exact = TRUE]], "1"),
    "certificate_method", path,
    "the versioned global analytic monotonicity certificate", x
  )
  .dpprior_schema_require(
    x[["relation", exact = TRUE]] %in% c("at_most", "at_least"),
    "certificate_relation", paste0(path, ".relation"),
    "at_most or at_least", x[["relation", exact = TRUE]]
  )
  .dpprior_schema_require(
    .dpprior_schema_is_count(x[["J", exact = TRUE]], 1L),
    "certificate_J", paste0(path, ".J"),
    "a positive integer cluster-support size", x[["J", exact = TRUE]]
  )
  needs_threshold <- x[["metric", exact = TRUE]] %in%
    c("wsb_tail", "wmax_tail_upper")
  needs_probability <- identical(
    x[["metric", exact = TRUE]], "wsb_quantile"
  )
  if (needs_threshold) {
    .dpprior_schema_validate_finite_scalar(
      x[["threshold", exact = TRUE]], paste0(path, ".threshold"),
      lower = 0, upper = 1
    )
  } else {
    .dpprior_schema_require(
      is.null(x[["threshold", exact = TRUE]]), "certificate_threshold",
      paste0(path, ".threshold"), "NULL for this metric",
      x[["threshold", exact = TRUE]]
    )
  }
  if (needs_probability) {
    .dpprior_schema_validate_finite_scalar(
      x[["probability", exact = TRUE]], paste0(path, ".probability"),
      lower = 0, upper = 1
    )
  } else {
    .dpprior_schema_require(
      is.null(x[["probability", exact = TRUE]]), "certificate_probability",
      paste0(path, ".probability"), "NULL for this metric",
      x[["probability", exact = TRUE]]
    )
  }
  .dpprior_schema_require(
    is.numeric(x[["support", exact = TRUE]]) &&
      !is.object(x[["support", exact = TRUE]]) &&
      is.null(dim(x[["support", exact = TRUE]])) &&
      .dpprior_schema_has_only_attributes(x[["support", exact = TRUE]]) &&
      identical(as.numeric(x[["support", exact = TRUE]]), c(0, 1)),
    "certificate_support", paste0(path, ".support"),
    "the exact probability support c(0, 1)",
    x[["support", exact = TRUE]]
  )
  .dpprior_schema_exact_names(
    x[["domain", exact = TRUE]], c("log_a", "log_b", "a", "b"),
    paste0(path, ".domain")
  )
  for (field in c("log_a", "log_b", "a", "b")) {
    bounds <- x[["domain", exact = TRUE]][[field, exact = TRUE]]
    .dpprior_schema_require(
      is.numeric(bounds) && !is.object(bounds) && is.null(dim(bounds)) &&
        .dpprior_schema_has_only_attributes(bounds) &&
        length(bounds) == 2L && !anyNA(bounds) && all(is.finite(bounds)) &&
        bounds[[1L]] < bounds[[2L]],
      "certificate_domain", paste0(path, ".domain.", field),
      "two finite increasing bounds", bounds
    )
  }
  .dpprior_schema_require(
    all(abs(exp(x[["domain", exact = TRUE]][["log_a", exact = TRUE]]) -
              x[["domain", exact = TRUE]][["a", exact = TRUE]]) <=
          1e-12 * pmax(1, x[["domain", exact = TRUE]][["a", exact = TRUE]])) &&
      all(abs(exp(x[["domain", exact = TRUE]][["log_b", exact = TRUE]]) -
                x[["domain", exact = TRUE]][["b", exact = TRUE]]) <=
            1e-12 * pmax(1, x[["domain", exact = TRUE]][["b", exact = TRUE]])),
    "certificate_domain", paste0(path, ".domain"),
    "a/b bounds equal exp(log_a/log_b) bounds", x[["domain", exact = TRUE]]
  )
  .dpprior_schema_exact_names(
    x[["monotonicity", exact = TRUE]], c("a", "b"),
    paste0(path, ".monotonicity")
  )
  .dpprior_schema_require(
    identical(
      x[["monotonicity", exact = TRUE]],
      list(a = "non_increasing", b = "non_decreasing")
    ),
    "certificate_monotonicity", paste0(path, ".monotonicity"),
    "the proved global monotonicity directions", x[["monotonicity", exact = TRUE]]
  )
  for (field in c("M_selected", "M_verification")) {
    .dpprior_schema_require(
      .dpprior_schema_is_count(x[[field, exact = TRUE]], 1L),
      "certificate_order", paste0(path, ".", field),
      "a positive integer quadrature order", x[[field, exact = TRUE]]
    )
  }
  .dpprior_schema_require(
    x[["M_verification", exact = TRUE]] >=
      x[["M_selected", exact = TRUE]],
    "certificate_order", paste0(path, ".M_verification"),
    "an order at least M_selected", x[["M_verification", exact = TRUE]]
  )
  .dpprior_schema_validate_finite_scalar(
    x[["effective_tolerance", exact = TRUE]],
    paste0(path, ".effective_tolerance"), lower = 0
  )
  minimum_eta <- c(
    log_a = x[["domain", exact = TRUE]][["log_a", exact = TRUE]][[2L]],
    log_b = x[["domain", exact = TRUE]][["log_b", exact = TRUE]][[1L]]
  )
  maximum_eta <- c(
    log_a = x[["domain", exact = TRUE]][["log_a", exact = TRUE]][[1L]],
    log_b = x[["domain", exact = TRUE]][["log_b", exact = TRUE]][[2L]]
  )
  .dpprior_validate_certificate_corner(
    x[["minimum", exact = TRUE]], x[["metric", exact = TRUE]],
    x[["threshold", exact = TRUE]], x[["probability", exact = TRUE]],
    x[["J", exact = TRUE]],
    x[["M_selected", exact = TRUE]],
    x[["M_verification", exact = TRUE]], minimum_eta,
    paste0(path, ".minimum")
  )
  .dpprior_validate_certificate_corner(
    x[["maximum", exact = TRUE]], x[["metric", exact = TRUE]],
    x[["threshold", exact = TRUE]], x[["probability", exact = TRUE]],
    x[["J", exact = TRUE]],
    x[["M_selected", exact = TRUE]],
    x[["M_verification", exact = TRUE]], maximum_eta,
    paste0(path, ".maximum")
  )
  for (field in c("target_value", "lower_bound", "upper_bound")) {
    .dpprior_schema_validate_finite_scalar(
      x[[field, exact = TRUE]], paste0(path, ".", field), lower = 0, upper = 1
    )
  }
  .dpprior_schema_validate_finite_scalar(
    x[["tolerance", exact = TRUE]], paste0(path, ".tolerance"), lower = 0
  )
  bound_tolerance <- 1e-12 * max(
    1, x[["minimum", exact = TRUE]][["lower", exact = TRUE]],
    x[["maximum", exact = TRUE]][["upper", exact = TRUE]],
    x[["lower_bound", exact = TRUE]], x[["upper_bound", exact = TRUE]]
  )
  .dpprior_schema_require(
    abs(x[["lower_bound", exact = TRUE]] -
          x[["minimum", exact = TRUE]][["lower", exact = TRUE]]) <=
        bound_tolerance &&
      abs(x[["upper_bound", exact = TRUE]] -
            x[["maximum", exact = TRUE]][["upper", exact = TRUE]]) <=
        bound_tolerance &&
      identical(x[["tolerance", exact = TRUE]],
                x[["effective_tolerance", exact = TRUE]]) &&
      x[["minimum", exact = TRUE]][["uncertainty", exact = TRUE]] <= 1e-6 &&
      x[["maximum", exact = TRUE]][["uncertainty", exact = TRUE]] <= 1e-6 &&
      x[["lower_bound", exact = TRUE]] <= x[["upper_bound", exact = TRUE]] &&
      identical(x[["source", exact = TRUE]],
                "phase8_metric_extrema_corner_enclosures"),
    "certificate_bounds", path,
    paste(
      "global bounds from the exact minimum/maximum refined corner",
      "enclosures, fixed effective tolerance, stable uncertainty, and closed",
      "source"
    ),
    x
  )
  .dpprior_schema_validate_scalar_logical(
    x[["certified", exact = TRUE]], paste0(path, ".certified")
  )
  expected <- if (identical(x[["relation", exact = TRUE]], "at_most")) {
    x[["lower_bound", exact = TRUE]] >
      x[["target_value", exact = TRUE]] + x[["tolerance", exact = TRUE]]
  } else {
    x[["upper_bound", exact = TRUE]] <
      x[["target_value", exact = TRUE]] - x[["tolerance", exact = TRUE]]
  }
  .dpprior_schema_require(
    identical(x[["certified", exact = TRUE]], expected) && expected,
    "certificate_truth", paste0(path, ".certified"),
    "TRUE only when a global bound proves the declared inequality infeasible",
    x[["certified", exact = TRUE]]
  )
  invisible(TRUE)
}


.dpprior_validate_hard_optimality <- function(x, available,
                                              path = "constraint.optimality") {
  .dpprior_schema_exact_names(x, .DPPRIOR_HARD_OPTIMALITY_FIELDS, path)
  .dpprior_schema_validate_scalar_logical(
    x[["performed", exact = TRUE]], paste0(path, ".performed")
  )
  .dpprior_schema_validate_scalar_logical(
    x[["passed", exact = TRUE]], paste0(path, ".passed")
  )
  .dpprior_schema_validate_scalar_character(
    x[["source", exact = TRUE]], paste0(path, ".source")
  )
  if (available) {
    .dpprior_schema_require(
      x[["performed", exact = TRUE]] &&
        is.null(x[["unavailable_reason", exact = TRUE]]),
      "hard_optimality", path,
      "performed evidence and no unavailable reason for a finite candidate",
      x
    )
    .dpprior_schema_validate_scalar_character(
      x[["selection_rule", exact = TRUE]], paste0(path, ".selection_rule")
    )
    .dpprior_schema_exact_names(
      x[["K_scales", exact = TRUE]], c("mean", "variance"),
      paste0(path, ".K_scales")
    )
    for (field in c("mean", "variance")) {
      .dpprior_schema_validate_finite_scalar(
        x[["K_scales", exact = TRUE]][[field, exact = TRUE]],
        paste0(path, ".K_scales.", field), lower = 0, lower_open = TRUE
      )
    }
    for (field in c("selected_K_loss", "minimum_K_loss", "tie_tolerance")) {
      .dpprior_schema_validate_finite_scalar(
        x[[field, exact = TRUE]], paste0(path, ".", field), lower = 0
      )
    }
    .dpprior_schema_validate_scalar_logical(
      x[["perturbation_passed", exact = TRUE]],
      paste0(path, ".perturbation_passed")
    )
    expected_pass <-
      x[["selected_K_loss", exact = TRUE]] <=
        x[["minimum_K_loss", exact = TRUE]] +
        x[["tie_tolerance", exact = TRUE]] &&
      x[["perturbation_passed", exact = TRUE]]
    .dpprior_schema_require(
      identical(x[["passed", exact = TRUE]], expected_pass),
      "hard_optimality", paste0(path, ".passed"),
      "selection within tie tolerance and passing perturbation evidence",
      x[["passed", exact = TRUE]]
    )
  } else {
    unavailable <- c(
      "selection_rule", "K_scales", "selected_K_loss", "minimum_K_loss",
      "tie_tolerance", "perturbation_passed"
    )
    .dpprior_schema_require(
      !x[["performed", exact = TRUE]] && !x[["passed", exact = TRUE]] &&
        all(vapply(x[unavailable], is.null, logical(1))),
      "hard_optimality", path,
      "no active optimality claims without a finite candidate", x
    )
    .dpprior_schema_validate_scalar_character(
      x[["unavailable_reason", exact = TRUE]],
      paste0(path, ".unavailable_reason")
    )
  }
  invisible(TRUE)
}


.dpprior_validate_soft_optimality <- function(x, available,
                                              path = "tradeoff.optimality") {
  .dpprior_schema_exact_names(x, .DPPRIOR_SOFT_OPTIMALITY_FIELDS, path)
  for (field in c("performed", "passed")) {
    .dpprior_schema_validate_scalar_logical(
      x[[field, exact = TRUE]], paste0(path, ".", field)
    )
  }
  .dpprior_schema_validate_scalar_character(
    x[["source", exact = TRUE]], paste0(path, ".source")
  )
  if (available) {
    .dpprior_schema_require(
      x[["performed", exact = TRUE]] &&
        is.null(x[["unavailable_reason", exact = TRUE]]),
      "soft_optimality", path,
      "performed evidence and no unavailable reason for a finite candidate",
      x
    )
    gradient <- x[["gradient", exact = TRUE]]
    .dpprior_schema_require(
      is.numeric(gradient) && !is.object(gradient) && is.null(dim(gradient)) &&
        .dpprior_schema_has_only_attributes(gradient, "names", "names") &&
        length(gradient) == 2L && !anyNA(gradient) && all(is.finite(gradient)) &&
        identical(names(gradient), c("log_shape", "log_rate")),
      "soft_gradient", paste0(path, ".gradient"),
      "finite named c(log_shape, log_rate) projected-gradient evidence",
      gradient
    )
    for (field in c("gradient_method", "bound_state", "stationarity_operator")) {
      values <- x[[field, exact = TRUE]]
      .dpprior_schema_require(
        is.character(values) && !is.object(values) && is.null(dim(values)) &&
          .dpprior_schema_has_only_attributes(values, "names", "names") &&
          identical(names(values), c("log_shape", "log_rate")) &&
          length(values) == 2L && !anyNA(values) && all(nzchar(values)),
        "soft_stationarity", paste0(path, ".", field),
        "an exact named two-component stationarity record", values
      )
    }
    .dpprior_schema_require(
      all(x[["bound_state", exact = TRUE]] %in%
            c("interior", "lower", "upper")),
      "soft_bound_state", paste0(path, ".bound_state"),
      "interior, lower, or upper for each log parameter",
      x[["bound_state", exact = TRUE]]
    )
    expected_stationarity_operator <- setNames(vapply(
      x[["bound_state", exact = TRUE]],
      function(state) switch(
        state, interior = "abs_lte", lower = "gte", upper = "lte"
      ),
      character(1)
    ), names(x[["bound_state", exact = TRUE]]))
    expected_gradient_method <- setNames(vapply(
      x[["bound_state", exact = TRUE]],
      function(state) switch(
        state,
        interior = "central finite-difference interior gradient",
        lower = "forward feasible-direction KKT gradient",
        upper = "backward feasible-direction KKT gradient"
      ),
      character(1)
    ), names(x[["bound_state", exact = TRUE]]))
    .dpprior_schema_require(
      identical(x[["stationarity_operator", exact = TRUE]],
                expected_stationarity_operator) &&
        identical(x[["gradient_method", exact = TRUE]],
                  expected_gradient_method),
      "soft_stationarity_method", path,
      "closed bound-state/operator/finite-difference-method identities",
      list(
        state = x[["bound_state", exact = TRUE]],
        operator = x[["stationarity_operator", exact = TRUE]],
        method = x[["gradient_method", exact = TRUE]]
      )
    )
    for (field in c(
      "recorded_objective", "recomputed_objective", "objective_difference",
      "objective_tolerance", "boundary_tolerance", "stationarity_tolerance",
      "local_base_objective", "start_objective", "candidate_objective",
      "start_tolerance", "selection_tolerance"
    )) {
      .dpprior_schema_validate_finite_scalar(
        x[[field, exact = TRUE]], paste0(path, ".", field),
        lower = 0
      )
    }
    for (field in c(
      "objective_passed", "stationarity_passed", "local_minimum_passed",
      "no_worse_start"
    )) {
      .dpprior_schema_validate_scalar_logical(
        x[[field, exact = TRUE]], paste0(path, ".", field)
      )
    }
    component_pass <- x[["component_pass", exact = TRUE]]
    .dpprior_schema_require(
      is.logical(component_pass) && !is.object(component_pass) &&
        is.null(dim(component_pass)) &&
        .dpprior_schema_has_only_attributes(component_pass, "names", "names") &&
        identical(names(component_pass), c("log_shape", "log_rate")) &&
        length(component_pass) == 2L && !anyNA(component_pass),
      "soft_stationarity_components", paste0(path, ".component_pass"),
      "an exact named two-component logical result", component_pass
    )
    neighbor_values <- x[["neighbor_objectives", exact = TRUE]]
    neighbor_tolerances <- x[["neighbor_tolerances", exact = TRUE]]
    for (field in c("neighbor_objectives", "neighbor_tolerances")) {
      values <- x[[field, exact = TRUE]]
      .dpprior_schema_require(
        is.numeric(values) && !is.object(values) && is.null(dim(values)) &&
          .dpprior_schema_has_only_attributes(values, "names", "names") &&
          length(values) >= 1L && !anyNA(values) && all(is.finite(values)) &&
          !is.null(names(values)) && !anyDuplicated(names(values)) &&
          all(nzchar(names(values))),
        "soft_local_optimality", paste0(path, ".", field),
        "a finite uniquely named vector of retained feasible neighbors",
        values
      )
    }
    .dpprior_schema_require(
      identical(names(neighbor_values), names(neighbor_tolerances)) &&
        all(names(neighbor_values) %in% c(
          "log_shape_plus", "log_shape_minus",
          "log_rate_plus", "log_rate_minus",
          "diagonal_pp", "diagonal_pm", "diagonal_mp", "diagonal_mm"
        )) &&
        all(neighbor_values >= 0) && all(neighbor_tolerances >= 0),
      "soft_local_optimality", path,
      "matching nonnegative objectives/tolerances for every retained neighbor",
      list(values = neighbor_values, tolerances = neighbor_tolerances)
    )
    objective_difference <- abs(
      x[["recorded_objective", exact = TRUE]] -
        x[["recomputed_objective", exact = TRUE]]
    )
    objective_pass <- objective_difference <=
      x[["objective_tolerance", exact = TRUE]]
    stationarity_components <- vapply(seq_along(gradient), function(index) {
      switch(
        x[["stationarity_operator", exact = TRUE]][[index]],
        abs_lte = abs(gradient[[index]]) <=
          x[["stationarity_tolerance", exact = TRUE]],
        gte = gradient[[index]] >=
          -x[["stationarity_tolerance", exact = TRUE]],
        lte = gradient[[index]] <=
          x[["stationarity_tolerance", exact = TRUE]]
      )
    }, logical(1))
    names(stationarity_components) <- names(gradient)
    stationarity <- all(stationarity_components)
    local_minimum <- all(
      x[["local_base_objective", exact = TRUE]] <=
        neighbor_values + neighbor_tolerances
    )
    no_worse_start <- x[["local_base_objective", exact = TRUE]] <=
      x[["start_objective", exact = TRUE]] +
      x[["start_tolerance", exact = TRUE]]
    expected_pass <- stationarity && objective_pass && local_minimum &&
      no_worse_start
    .dpprior_schema_require(
      identical(x[["objective_difference", exact = TRUE]],
                objective_difference) &&
        identical(x[["objective_passed", exact = TRUE]], objective_pass) &&
        identical(x[["component_pass", exact = TRUE]],
                  stationarity_components) &&
        identical(x[["stationarity_passed", exact = TRUE]], stationarity) &&
        identical(x[["local_minimum_passed", exact = TRUE]], local_minimum) &&
        identical(x[["no_worse_start", exact = TRUE]], no_worse_start) &&
        identical(x[["passed", exact = TRUE]], expected_pass),
      "soft_optimality", path,
      paste(
        "recomputed objective identity, projected stationarity, retained",
        "coordinate/diagonal neighbors, and no-worse-start evidence"
      ),
      x
    )
  } else {
    unavailable <- setdiff(
      .DPPRIOR_SOFT_OPTIMALITY_FIELDS,
      c("performed", "passed", "source", "unavailable_reason")
    )
    .dpprior_schema_require(
      !x[["performed", exact = TRUE]] && !x[["passed", exact = TRUE]] &&
        all(vapply(x[unavailable], is.null, logical(1))),
      "soft_optimality", path,
      "no active optimality claims without a finite candidate", x
    )
    .dpprior_schema_validate_scalar_character(
      x[["unavailable_reason", exact = TRUE]],
      paste0(path, ".unavailable_reason")
    )
  }
  invisible(TRUE)
}


.dpprior_expected_soft_optimality <- function(raw, lower_bounds,
                                               upper_bounds) {
  .dpprior_schema_require(
    exists(".dpprior_soft_objective", mode = "function", inherits = TRUE) &&
      length(lower_bounds) == 2L && length(upper_bounds) == 2L &&
      identical(lower_bounds[[1L]], lower_bounds[[2L]]) &&
      identical(upper_bounds[[1L]], upper_bounds[[2L]]),
    "soft_objective_authority", "result.computation.attempts.selected.bounds",
    "the authoritative soft objective and one common Phase 8 log-bound pair",
    list(lower = lower_bounds, upper = upper_bounds)
  )
  target_K <- raw[["target", exact = TRUE]][["K", exact = TRUE]]
  target_moments <- target_K[["implied", exact = TRUE]]
  input_fit <- raw[["provenance", exact = TRUE]][["input_fit", exact = TRUE]]
  fit_info <- list(
    J = raw[["J", exact = TRUE]],
    a = input_fit[["parameters", exact = TRUE]][["a", exact = TRUE]],
    b = input_fit[["parameters", exact = TRUE]][["b", exact = TRUE]],
    target_K = list(
      mu_K = target_moments[["mean", exact = TRUE]],
      var_K = target_moments[["variance", exact = TRUE]]
    )
  )
  weight_target <- raw[["target", exact = TRUE]][["weight", exact = TRUE]]
  tradeoff <- raw[["tradeoff", exact = TRUE]]
  K_scales <- unlist(
    tradeoff[["scales", exact = TRUE]][["K", exact = TRUE]],
    use.names = TRUE
  )
  orders <- raw[["computation", exact = TRUE]][["orders", exact = TRUE]]
  M_selected <- orders[["M_selected", exact = TRUE]]
  M_verification <- orders[["M_verification_used", exact = TRUE]]
  log_bounds <- c(lower_bounds[[1L]], upper_bounds[[1L]])
  selected_objective <- .dpprior_soft_objective(
    fit_info, weight_target, tradeoff[["lambda", exact = TRUE]], M_selected,
    K_scales, log_bounds
  )
  refined_objective <- .dpprior_soft_objective(
    fit_info, weight_target, tradeoff[["lambda", exact = TRUE]], M_verification,
    K_scales, log_bounds
  )
  eta <- log(c(
    raw[["parameters", exact = TRUE]][["a", exact = TRUE]],
    raw[["parameters", exact = TRUE]][["b", exact = TRUE]]
  ))
  names(eta) <- c("log_shape", "log_rate")
  selected_value <- selected_objective(eta)
  refined_value <- refined_objective(eta)
  stationarity_step <- raw[["tolerances", exact = TRUE]][[
    "stationarity", exact = TRUE
  ]][["step", exact = TRUE]]
  boundary_tolerance <- raw[["tolerances", exact = TRUE]][[
    "boundary", exact = TRUE
  ]]
  gradient <- setNames(numeric(2L), names(eta))
  bound_state <- gradient_method <- stationarity_operator <-
    setNames(character(2L), names(eta))
  for (index in seq_along(eta)) {
    lower_distance <- eta[[index]] - lower_bounds[[index]]
    upper_distance <- upper_bounds[[index]] - eta[[index]]
    on_lower <- lower_distance <= boundary_tolerance
    on_upper <- upper_distance <= boundary_tolerance
    if (!on_lower && !on_upper) {
      step <- min(
        stationarity_step, max(lower_distance, 0) / 2,
        max(upper_distance, 0) / 2
      )
      plus <- minus <- eta
      plus[[index]] <- plus[[index]] + step
      minus[[index]] <- minus[[index]] - step
      gradient[[index]] <-
        (refined_objective(plus) - refined_objective(minus)) / (2 * step)
      bound_state[[index]] <- "interior"
      gradient_method[[index]] <-
        "central finite-difference interior gradient"
      stationarity_operator[[index]] <- "abs_lte"
    } else if (on_lower) {
      step <- min(stationarity_step, max(upper_distance, 0) / 2)
      plus <- eta
      plus[[index]] <- plus[[index]] + step
      gradient[[index]] <- (refined_objective(plus) - refined_value) / step
      bound_state[[index]] <- "lower"
      gradient_method[[index]] <- "forward feasible-direction KKT gradient"
      stationarity_operator[[index]] <- "gte"
    } else {
      step <- min(stationarity_step, max(lower_distance, 0) / 2)
      minus <- eta
      minus[[index]] <- minus[[index]] - step
      gradient[[index]] <- (refined_value - refined_objective(minus)) / step
      bound_state[[index]] <- "upper"
      gradient_method[[index]] <- "backward feasible-direction KKT gradient"
      stationarity_operator[[index]] <- "lte"
    }
  }
  directions <- rbind(
    log_shape_plus = c(1, 0), log_shape_minus = c(-1, 0),
    log_rate_plus = c(0, 1), log_rate_minus = c(0, -1),
    diagonal_pp = c(1, 1) / sqrt(2),
    diagonal_pm = c(1, -1) / sqrt(2),
    diagonal_mp = c(-1, 1) / sqrt(2),
    diagonal_mm = c(-1, -1) / sqrt(2)
  )
  local_step <- raw[["tolerances", exact = TRUE]][[
    "neighborhood", exact = TRUE
  ]][["step", exact = TRUE]]
  proposals <- t(vapply(seq_len(nrow(directions)), function(index) {
    eta + local_step * directions[index, ]
  }, numeric(2L)))
  feasible <- apply(proposals, 1L, function(proposal) {
    all(proposal >= lower_bounds) && all(proposal <= upper_bounds)
  })
  neighbor_values <- vapply(
    which(feasible), function(index) refined_objective(proposals[index, ]),
    numeric(1)
  )
  names(neighbor_values) <- rownames(directions)[feasible]
  input_eta <- log(c(
    input_fit[["parameters", exact = TRUE]][["a", exact = TRUE]],
    input_fit[["parameters", exact = TRUE]][["b", exact = TRUE]]
  ))
  names(input_eta) <- names(eta)
  expected <- list(
    eta = eta, selected_objective = selected_value,
    refined_objective = refined_value, gradient = gradient,
    gradient_method = gradient_method, bound_state = bound_state,
    stationarity_operator = stationarity_operator,
    neighbor_objectives = neighbor_values,
    start_objective = refined_objective(input_eta), input_eta = input_eta
  )
  objective_values <- c(
    selected_value, refined_value, neighbor_values, expected$start_objective
  )
  .dpprior_schema_require(
    length(neighbor_values) >= 1L && !anyNA(c(objective_values, gradient)) &&
      all(is.finite(c(objective_values, gradient))) &&
      all(objective_values >= 0),
    "soft_objective_authority", "result.tradeoff.optimality",
    "finite nonnegative selected/refined/start/neighbor objectives and gradient",
    expected
  )
  expected
}


.dpprior_validate_constraint_extension <- function(x, path = "constraint") {
  .dpprior_schema_exact_names(
    x,
    c("relation", "operator", "residual", "slack", "tolerance", "satisfied",
      "active", "feasibility", "optimality"),
    path
  )
  .dpprior_schema_validate_scalar_character(
    x$relation, paste0(path, ".relation")
  )
  .dpprior_schema_validate_scalar_character(
    x$operator, paste0(path, ".operator")
  )
  .dpprior_schema_require(
    x$relation %in% c("at_most", "at_least"),
    "constraint_relation", paste0(path, ".relation"),
    "at_most or at_least", x$relation
  )
  expected_operator <- if (identical(x$relation, "at_most")) "<=" else ">="
  .dpprior_schema_require(
    identical(x$operator, expected_operator),
    "constraint_operator", paste0(path, ".operator"), expected_operator,
    x$operator
  )
  candidate_available <- !is.null(x[["residual", exact = TRUE]])
  candidate_fields <- c("residual", "slack", "tolerance", "satisfied", "active")
  if (candidate_available) {
    .dpprior_schema_validate_finite_scalar(
      x$residual, paste0(path, ".residual")
    )
    .dpprior_schema_validate_finite_scalar(x$slack, paste0(path, ".slack"))
    .dpprior_schema_validate_named_list(
      x$tolerance, paste0(path, ".tolerance")
    )
    .dpprior_schema_require(
      "effective" %in% names(x$tolerance), "constraint_tolerance",
      paste0(path, ".tolerance"), "an exact effective tolerance field",
      names(x$tolerance)
    )
    .dpprior_schema_validate_finite_scalar(
      x$tolerance[["effective", exact = TRUE]],
      paste0(path, ".tolerance.effective"), lower = 0
    )
    coherence_tolerance <- 1e-12 * max(1, abs(x$residual), abs(x$slack))
    .dpprior_schema_require(
      abs(x$slack + x$residual) <= coherence_tolerance,
      "constraint_slack", paste0(path, ".slack"),
      "-residual within central tolerance", x$slack
    )
    .dpprior_schema_validate_scalar_logical(
      x$satisfied, paste0(path, ".satisfied")
    )
    .dpprior_schema_validate_scalar_logical(x$active, paste0(path, ".active"))
  } else {
    .dpprior_schema_require(
      all(vapply(
        x[candidate_fields], is.null, logical(1)
      )),
      "constraint_availability", path,
      paste(
        "all candidate constraint fields NULL when no candidate residual",
        "is available"
      ),
      x[candidate_fields]
    )
  }
  .dpprior_schema_exact_names(
    x$feasibility,
    c("classification", "certified_infeasible", "feasibility_unknown",
      "certificate", "candidate_count", "verified_candidate_count",
      "feasible_candidate_count"),
    paste0(path, ".feasibility")
  )
  .dpprior_schema_validate_scalar_character(
    x$feasibility$classification, paste0(path, ".feasibility.classification")
  )
  .dpprior_schema_require(
    x$feasibility$classification %in%
      c("feasible_candidate", "certified_infeasible", "unknown"),
    "feasibility", paste0(path, ".feasibility.classification"),
    "feasible_candidate, certified_infeasible, or unknown",
    x$feasibility$classification
  )
  .dpprior_schema_validate_scalar_logical(
    x$feasibility$certified_infeasible,
    paste0(path, ".feasibility.certified_infeasible")
  )
  .dpprior_schema_validate_scalar_logical(
    x$feasibility$feasibility_unknown,
    paste0(path, ".feasibility.feasibility_unknown")
  )
  for (field in c(
    "candidate_count", "verified_candidate_count", "feasible_candidate_count"
  )) {
    .dpprior_schema_require(
      .dpprior_schema_is_count(x$feasibility[[field]], 0L), "count",
      paste0(path, ".feasibility.", field), "a non-negative integer",
      x$feasibility[[field]]
    )
  }
  if (!is.null(x$feasibility$certificate)) {
    .dpprior_validate_infeasibility_certificate(
      x$feasibility$certificate, paste0(path, ".feasibility.certificate")
    )
  }
  .dpprior_schema_require(
    !(x$feasibility$certified_infeasible &&
        x$feasibility$feasibility_unknown),
    "feasibility", paste0(path, ".feasibility"),
    "certified_infeasible and feasibility_unknown not both TRUE",
    x$feasibility
  )
  .dpprior_schema_require(
    identical(
      x$feasibility$certified_infeasible,
      identical(x$feasibility$classification, "certified_infeasible")
    ) && identical(
      x$feasibility$feasibility_unknown,
      identical(x$feasibility$classification, "unknown")
    ),
    "feasibility", paste0(path, ".feasibility"),
    "classification consistent with the two feasibility flags", x$feasibility
  )
  .dpprior_schema_require(
    x$feasibility$verified_candidate_count <= x$feasibility$candidate_count &&
      x$feasibility$feasible_candidate_count <=
        x$feasibility$verified_candidate_count,
    "feasibility_counts", paste0(path, ".feasibility"),
    paste(
      "feasible_candidate_count <= verified_candidate_count <=",
      "candidate_count"
    ), x$feasibility
  )
  if (identical(x$feasibility$classification, "feasible_candidate")) {
    .dpprior_schema_require(
      x$feasibility$candidate_count >= 1L,
      "feasibility_counts", paste0(path, ".feasibility.candidate_count"),
      "at least one candidate for feasible_candidate", x$feasibility$candidate_count
    )
  }
  .dpprior_schema_require(
    identical(x$feasibility$classification, "certified_infeasible") ==
      !is.null(x$feasibility$certificate),
    "certificate_presence", paste0(path, ".feasibility.certificate"),
    "a substantive certificate exactly for certified_infeasible",
    x$feasibility$certificate
  )
  if (candidate_available) {
    residual_satisfied <- x$residual <=
      x$tolerance[["effective", exact = TRUE]]
    .dpprior_schema_require(
      !x$satisfied || residual_satisfied,
      "constraint_satisfaction", paste0(path, ".satisfied"),
      paste(
        "TRUE only when the selected residual is within tolerance; the final",
        "FALSE/TRUE two-order identity is enforced by the hard result validator"
      ), x$satisfied
    )
  }
  .dpprior_validate_hard_optimality(
    x[["optimality", exact = TRUE]], candidate_available,
    paste0(path, ".optimality")
  )
  invisible(TRUE)
}


.dpprior_new_constraint <- function(relation,
                                    operator,
                                    residual,
                                    slack,
                                    tolerance,
                                    satisfied,
                                    active,
                                    feasibility,
                                    optimality = list()) {
  out <- list(
    relation = relation,
    operator = operator,
    residual = residual,
    slack = slack,
    tolerance = tolerance,
    satisfied = satisfied,
    active = active,
    feasibility = feasibility,
    optimality = optimality
  )
  .dpprior_validate_constraint_extension(out)
  out
}


.dpprior_validate_tradeoff_extension <- function(x, path = "tradeoff") {
  .dpprior_schema_exact_names(
    x,
    c("lambda", "K_loss", "weight_loss", "total_loss", "target_residual",
      "directed_residual", "scales", "endpoint", "optimality"),
    path
  )
  .dpprior_schema_validate_finite_scalar(
    x$lambda, paste0(path, ".lambda"), lower = 0, upper = 1
  )
  decision_fields <- c(
    "K_loss", "weight_loss", "total_loss", "target_residual",
    "directed_residual"
  )
  candidate_available <- !is.null(x[["K_loss", exact = TRUE]])
  if (candidate_available) {
    for (field in c("K_loss", "weight_loss", "total_loss")) {
      .dpprior_schema_validate_finite_scalar(
        x[[field]], paste0(path, ".", field), lower = 0
      )
    }
    for (field in c("target_residual", "directed_residual")) {
      .dpprior_schema_validate_finite_scalar(
        x[[field]], paste0(path, ".", field)
      )
    }
  } else {
    .dpprior_schema_require(
      all(vapply(x[decision_fields], is.null, logical(1))),
      "tradeoff_availability", path,
      "all soft loss/residual fields NULL when no candidate is available",
      x[decision_fields]
    )
  }
  .dpprior_schema_validate_named_list(x$scales, paste0(path, ".scales"))
  .dpprior_schema_require(
    identical(names(x$scales), c("K", "weight")),
    "soft_scales", paste0(path, ".scales"),
    "exact K and weight scales", names(x$scales)
  )
  .dpprior_schema_exact_names(
    x$scales[["K", exact = TRUE]], c("mean", "variance"),
    paste0(path, ".scales.K")
  )
  for (field in c("mean", "variance")) {
    .dpprior_schema_validate_finite_scalar(
      x$scales[["K", exact = TRUE]][[field, exact = TRUE]],
      paste0(path, ".scales.K.", field),
      lower = 0, lower_open = TRUE
    )
  }
  .dpprior_schema_validate_finite_scalar(
    x$scales[["weight", exact = TRUE]], paste0(path, ".scales.weight"),
    lower = 0, lower_open = TRUE
  )
  .dpprior_schema_validate_scalar_logical(x$endpoint, paste0(path, ".endpoint"))
  if (candidate_available) {
    expected <- x$lambda * x$K_loss + (1 - x$lambda) * x$weight_loss
    tolerance <- 1e-12 * max(1, abs(expected), abs(x$total_loss))
    .dpprior_schema_require(
      abs(x$total_loss - expected) <= tolerance,
      "soft_total_loss", paste0(path, ".total_loss"),
      "lambda*K_loss + (1-lambda)*weight_loss", x$total_loss
    )
  }
  .dpprior_validate_soft_optimality(
    x[["optimality", exact = TRUE]],
    candidate_available && !x[["endpoint", exact = TRUE]],
    paste0(path, ".optimality")
  )
  .dpprior_schema_require(
    identical(x$endpoint, identical(x$lambda, 1)),
    "soft_endpoint", paste0(path, ".endpoint"),
    "TRUE exactly when lambda equals 1", x$endpoint
  )
  invisible(TRUE)
}


.dpprior_new_tradeoff <- function(lambda,
                                   K_loss,
                                   weight_loss,
                                   total_loss,
                                   target_residual,
                                   directed_residual,
                                   scales,
                                   endpoint = identical(lambda, 1),
                                   optimality = list()) {
  out <- list(
    lambda = lambda,
    K_loss = K_loss,
    weight_loss = weight_loss,
    total_loss = total_loss,
    target_residual = target_residual,
    directed_residual = directed_residual,
    scales = scales,
    endpoint = endpoint,
    optimality = optimality
  )
  .dpprior_validate_tradeoff_extension(out)
  out
}


.dpprior_validate_legacy_extension <- function(x, path = "legacy") {
  .dpprior_schema_exact_names(
    x,
    c("contract", "lambda", "losses", "approximation_opt_in", "warning_code"),
    path
  )
  .dpprior_schema_validate_scalar_character(x$contract, paste0(path, ".contract"))
  .dpprior_schema_validate_finite_scalar(
    x$lambda, paste0(path, ".lambda"), lower = 0, upper = 1
  )
  .dpprior_schema_validate_named_list(x$losses, paste0(path, ".losses"))
  .dpprior_schema_validate_plain_record_value(
    x[["losses", exact = TRUE]], paste0(path, ".losses")
  )
  .dpprior_schema_validate_scalar_logical(
    x$approximation_opt_in, paste0(path, ".approximation_opt_in")
  )
  .dpprior_schema_validate_scalar_character(
    x$warning_code, paste0(path, ".warning_code")
  )
  invisible(TRUE)
}


.dpprior_legacy_numeric_close <- function(x, y) {
  .dpprior_schema_is_finite_scalar(x) &&
    .dpprior_schema_is_finite_scalar(y) &&
    abs(x - y) <= 64 * .Machine$double.eps * max(1, abs(x), abs(y))
}


.dpprior_legacy_weight_truth <- function(target, parameters, J, M, path) {
  .dpprior_schema_require(
    exists(".dpprior_v2_eval_metric", mode = "function", inherits = TRUE),
    "legacy_weight_evaluator", path,
    "the authoritative named Phase-8 metric evaluator", NULL
  )
  spec <- list(
    metric = target[["metric", exact = TRUE]],
    estimand = target[["estimand", exact = TRUE]],
    units = target[["units", exact = TRUE]],
    threshold = target[["threshold", exact = TRUE]],
    probability = target[["probability", exact = TRUE]]
  )
  evaluated <- tryCatch(
    .dpprior_v2_eval_metric(
      spec,
      parameters[["a", exact = TRUE]],
      parameters[["b", exact = TRUE]],
      J = J, M = M
    ),
    error = identity
  )
  .dpprior_schema_require(
    !inherits(evaluated, "condition") &&
      typeof(evaluated) == "list" && is.list(evaluated) &&
      .dpprior_schema_is_finite_scalar(evaluated[["value", exact = TRUE]]),
    "legacy_weight_evaluator", path,
    "one finite value from the authoritative named metric evaluator",
    if (inherits(evaluated, "condition")) conditionMessage(evaluated) else
      evaluated
  )
  evaluated[["value", exact = TRUE]]
}


.dpprior_legacy_loss_truth <- function(parameters, J, M, target_K,
                                       target_weight, loss_type, scales,
                                       path) {
  moments <- tryCatch(
    exact_K_moments(
      J,
      parameters[["a", exact = TRUE]],
      parameters[["b", exact = TRUE]],
      M = M
    ),
    error = identity
  )
  .dpprior_schema_require(
    !inherits(moments, "condition") &&
      typeof(moments) == "list" && is.list(moments) &&
      .dpprior_schema_is_finite_scalar(moments[["mean", exact = TRUE]]) &&
      .dpprior_schema_is_finite_scalar(moments[["var", exact = TRUE]]),
    "legacy_K_recomputation", path,
    "finite exact K moments at the retained order",
    if (inherits(moments, "condition")) conditionMessage(moments) else moments
  )
  weight_value <- .dpprior_legacy_weight_truth(
    target_weight, parameters, J, M, paste0(path, ".weight")
  )
  implied <- target_K[["implied", exact = TRUE]]
  K_residual <- c(
    mean = moments[["mean", exact = TRUE]] -
      implied[["mean", exact = TRUE]],
    variance = moments[["var", exact = TRUE]] -
      implied[["variance", exact = TRUE]]
  )
  weight_residual <- weight_value - target_weight[["value", exact = TRUE]]
  if (identical(loss_type, "absolute")) {
    raw_K_loss <- sum(K_residual^2)
  } else {
    .dpprior_schema_require(
      implied[["mean", exact = TRUE]] > 0 &&
        implied[["variance", exact = TRUE]] > 0,
      "legacy_relative_target", paste0(path, ".target_K.implied"),
      "strictly positive target mean and variance for relative loss",
      implied
    )
    raw_K_loss <-
      (K_residual[["mean"]] / implied[["mean", exact = TRUE]])^2 +
      (K_residual[["variance"]] /
         implied[["variance", exact = TRUE]])^2
  }
  raw_weight_loss <- weight_residual^2
  K_loss <- raw_K_loss
  weight_loss <- raw_weight_loss
  if (identical(loss_type, "adaptive") &&
      !is.null(scales[["L_K_scale", exact = TRUE]]) &&
      !is.null(scales[["L_w_scale", exact = TRUE]])) {
    K_loss <- raw_K_loss / scales[["L_K_scale", exact = TRUE]]
    weight_loss <- raw_weight_loss / scales[["L_w_scale", exact = TRUE]]
  }
  list(
    moments = moments,
    weight_value = weight_value,
    K_residual = K_residual,
    weight_residual = weight_residual,
    raw_K_loss = raw_K_loss,
    raw_weight_loss = raw_weight_loss,
    K_loss = K_loss,
    weight_loss = weight_loss
  )
}


.dpprior_schema_is_retained_v2_package_version <- function(x) {
  .dpprior_schema_is_scalar_character(x) &&
    grepl("^2(?:[.-][0-9]+)+\\z", x, perl = TRUE)
}


.dpprior_is_dual_legacy_migration_route <- function(raw) {
  provenance <- raw[["provenance", exact = TRUE]]
  computation <- raw[["computation", exact = TRUE]]
  verification <- raw[["verification", exact = TRUE]]
  compatibility <- raw[["compatibility", exact = TRUE]]
  legacy <- raw[["legacy", exact = TRUE]]
  migration <- provenance[["migration", exact = TRUE]]
  backend <- provenance[["backend", exact = TRUE]]
  termination <- computation[["termination", exact = TRUE]]
  deprecation <- compatibility[["deprecations", exact = TRUE]][[
    "legacy_schema", exact = TRUE
  ]]

  identical(migration[["source_schema", exact = TRUE]], "DPprior/1.1/fit") &&
    identical(migration[["adapter", exact = TRUE]],
              "upgrade_DPprior_object") &&
    identical(migration[["lossless", exact = TRUE]], FALSE) &&
    identical(backend[["package", exact = TRUE]], "DPprior") &&
    .dpprior_schema_is_retained_v2_package_version(
      backend[["package_version", exact = TRUE]]
    ) &&
    identical(backend[["implementation", exact = TRUE]],
              "schema_upgrade_v1") &&
    is.null(backend[["source_commit", exact = TRUE]]) &&
    identical(provenance[["legacy", exact = TRUE]], list(
      active = TRUE,
      contract = "v1.1_path_scaled_soft_equality_loss",
      deprecation_stage = "v2_migration"
    )) &&
    identical(legacy[["contract", exact = TRUE]],
              "v1.1_path_scaled_soft_equality_loss") &&
    identical(legacy[["warning_code", exact = TRUE]],
              "legacy_object_upgraded") &&
    identical(verification[["method", exact = TRUE]],
              "legacy_evidence_quarantine") &&
    identical(termination[["code", exact = TRUE]],
              "legacy_migration_no_selection") &&
    identical(termination[["source", exact = TRUE]], "constructor") &&
    length(computation[["attempts", exact = TRUE]]) == 0L &&
    length(computation[["candidate_evaluations", exact = TRUE]]) == 0L &&
    is.null(computation[["selected_candidate_id", exact = TRUE]]) &&
    is.null(computation[["selected_attempt_id", exact = TRUE]]) &&
    is.null(provenance[["input_fit", exact = TRUE]]) &&
    identical(deprecation, list(
      code = "legacy_object_upgraded",
      first_deprecated_version = "2.0.0",
      removal_floor = "not_scheduled"
    ))
}


.dpprior_validate_dual_legacy_migration_boundary <- function(raw) {
  provenance <- raw[["provenance", exact = TRUE]]
  computation <- raw[["computation", exact = TRUE]]
  verification <- raw[["verification", exact = TRUE]]
  compatibility <- raw[["compatibility", exact = TRUE]]
  legacy <- raw[["legacy", exact = TRUE]]
  parameters <- raw[["parameters", exact = TRUE]]
  achieved <- raw[["achieved", exact = TRUE]]
  residuals <- raw[["residuals", exact = TRUE]]
  target_K <- raw[["target", exact = TRUE]][["K", exact = TRUE]]
  target_weight <- raw[["target", exact = TRUE]][["weight", exact = TRUE]]
  if (is.object(target_K)) target_K <- unclass(target_K)
  if (is.object(target_weight)) target_weight <- unclass(target_weight)
  opt_in <- legacy[["approximation_opt_in", exact = TRUE]]
  backend <- provenance[["backend", exact = TRUE]]
  target_backend <- target_K[["provenance", exact = TRUE]][[
    "backend", exact = TRUE
  ]]
  target_source_schema <- target_K[["provenance", exact = TRUE]][[
    "migration", exact = TRUE
  ]][["source_schema", exact = TRUE]]
  target_derivation_source_schema <- target_K[["derivation", exact = TRUE]][[
    "request_to_normalized", exact = TRUE
  ]][["evidence", exact = TRUE]][["source_schema", exact = TRUE]]
  weight_source_schema <- target_weight[["provenance", exact = TRUE]][[
    "source_schema", exact = TRUE
  ]]
  migration_source_schema <- provenance[["migration", exact = TRUE]][[
    "source_schema", exact = TRUE
  ]]

  expected_message <- paste(
    "Migrated v1.1 candidate is approximate and unverified because optimizer",
    "selection and independent verification lineage were not serialized;",
    "refit is required for decision readiness."
  )
  .dpprior_schema_require(
    .dpprior_is_dual_legacy_migration_route(raw) &&
      identical(raw[["object_type", exact = TRUE]], "fit") &&
      identical(raw[["mode", exact = TRUE]], "dual_legacy") &&
      identical(raw[["method", exact = TRUE]], "dual-anchor") &&
      identical(raw[["status", exact = TRUE]], "approximate") &&
      identical(raw[["usable", exact = TRUE]], opt_in) &&
      identical(raw[["verified", exact = TRUE]], FALSE) &&
      identical(raw[["message", exact = TRUE]], expected_message) &&
      identical(raw[["tolerances", exact = TRUE]], list()) &&
      identical(target_backend, backend) &&
      identical(target_source_schema, migration_source_schema) &&
      identical(target_derivation_source_schema, migration_source_schema) &&
      identical(weight_source_schema, migration_source_schema) &&
      identical(parameters[["parameterization", exact = TRUE]],
                "Gamma(shape=a, rate=b)"),
    "legacy_migration_route", "result",
    paste(
      "the exact finite v1.1 dual migration route, opt-in usability,",
      "message, empty tolerances, nested v1.1 source identity, and Gamma",
      "parameterization"
    ), raw[c("mode", "method", "status", "usable", "verified")]
  )

  expected_setting <- list(
    method = "dual-anchor", controls = list(),
    parameterization = "Gamma(shape=a, rate=b)"
  )
  expected_computation <- list(
    request = expected_setting, used = expected_setting,
    orders = list(
      M_requested = NULL, M_selected = NULL,
      M_verification_required = NULL, M_verification_used = NULL,
      requested_reason = "not_recorded_by_v1.1",
      selected_reason = "not_recorded_by_v1.1",
      verification_required_reason = "not_recorded_by_v1.1",
      verification_used_reason = "not_recorded_by_v1.1"
    ),
    scaling = list(
      requested = NULL, used = NULL,
      formula = "legacy_evidence_unavailable", values = list(),
      fixed_from_input = FALSE, change_reason = ""
    ),
    attempts = list(), candidate_evaluations = list(),
    selected_candidate_id = NULL, selected_attempt_id = NULL,
    fallback = list(
      attempted = FALSE, used = FALSE, trigger_attempt_id = NULL,
      selected_attempt_id = NULL, reason_code = NULL,
      message = "v1.1 fallback lineage was not reconstructed",
      outcome = "not_attempted"
    ),
    termination = list(
      code = "legacy_migration_no_selection",
      message = paste(
        "The v1.1 status and optimizer exit were not reused as candidate",
        "selection evidence."
      ), source = "constructor", iterations = NULL, boundary_reason = NULL
    ),
    trace = NULL, resources = list()
  )
  .dpprior_schema_require(
    identical(computation, expected_computation),
    "legacy_migration_computation", "result.computation",
    "the exact control-free v1.1 migration computation quarantine",
    computation
  )

  selected_snapshot <- verification[["selected_snapshot", exact = TRUE]]
  expected_reason <- paste(
    "Stored selected values were quarantined; v1.1 did not retain an",
    "independent verifier or candidate-selection lineage."
  )
  .dpprior_schema_require(
    identical(verification[["method", exact = TRUE]],
              "legacy_evidence_quarantine") &&
      identical(verification[["performed", exact = TRUE]], FALSE) &&
      identical(verification[["passed", exact = TRUE]], FALSE) &&
      identical(verification[["reason", exact = TRUE]], expected_reason) &&
      identical(verification[["settings", exact = TRUE]], list()) &&
      identical(verification[["verifier_snapshot", exact = TRUE]], NULL) &&
      identical(verification[["stability", exact = TRUE]], NULL) &&
      identical(verification[["components", exact = TRUE]], list()) &&
      identical(verification[["invariants", exact = TRUE]], list()) &&
      identical(selected_snapshot[["parameters", exact = TRUE]], parameters) &&
      identical(selected_snapshot[["M", exact = TRUE]], NULL) &&
      identical(selected_snapshot[["achieved", exact = TRUE]], achieved) &&
      identical(selected_snapshot[["residuals", exact = TRUE]], residuals) &&
      identical(selected_snapshot[["tolerances", exact = TRUE]], list()) &&
      identical(selected_snapshot[["finite", exact = TRUE]], TRUE) &&
      identical(selected_snapshot[["source", exact = TRUE]],
                "legacy_serialized_selected_values"),
    "legacy_migration_verification", "result.verification",
    "the exact unperformed quarantine byte-bound to the public candidate",
    verification
  )

  expected_missing <- c(
    "requested_controls", "M_requested", "M_selected",
    "M_verification_required", "M_verification_used", "scaling",
    "canonical_attempts", "candidate_selection", "fallback_lineage",
    "independent_verifier_snapshot", "source_commit",
    "legacy_objective_scaling", "legacy_optimizer_lineage"
  )
  .dpprior_schema_require(
    identical(provenance[["requested_method", exact = TRUE]],
              "dual-anchor") &&
      identical(provenance[["selected_method", exact = TRUE]],
                "dual-anchor") &&
      identical(provenance[["is_fallback", exact = TRUE]], FALSE) &&
      identical(provenance[["approximation", exact = TRUE]], list(
        active = TRUE, opt_in = opt_in, kind = "legacy_schema_migration",
        warning_code = "legacy_object_upgraded"
      )) &&
      identical(provenance[["parameterization", exact = TRUE]],
                "Gamma(shape=a, rate=b)") &&
      identical(provenance[["migration", exact = TRUE]], list(
        source_schema = "DPprior/1.1/fit",
        adapter = "upgrade_DPprior_object", lossless = FALSE,
        missing_evidence = expected_missing,
        warnings = "legacy_object_upgraded"
      )),
    "legacy_migration_provenance", "result.provenance",
    "the exact v1.1 schema-upgrade provenance and missing-evidence inventory",
    provenance
  )

  .dpprior_schema_require(
    identical(names(legacy), c(
      "contract", "lambda", "losses", "approximation_opt_in", "warning_code"
    )) &&
      identical(names(legacy[["losses", exact = TRUE]]),
                c("K_loss", "loss_type")) &&
      .dpprior_schema_is_finite_scalar(legacy[["lambda", exact = TRUE]]) &&
      legacy[["lambda", exact = TRUE]] >= 0 &&
      legacy[["lambda", exact = TRUE]] <= 1 &&
      .dpprior_schema_is_finite_scalar(legacy[["losses", exact = TRUE]][[
        "K_loss", exact = TRUE
      ]]) && legacy[["losses", exact = TRUE]][[
        "K_loss", exact = TRUE
      ]] >= 0 &&
      legacy[["losses", exact = TRUE]][["loss_type", exact = TRUE]] %in%
        c("relative", "adaptive", "absolute") &&
      identical(opt_in, provenance[["approximation", exact = TRUE]][[
        "opt_in", exact = TRUE
      ]]),
    "legacy_migration_extension", "result.legacy",
    "the exact typed retained legacy loss and opt-in extension", legacy
  )

  .dpprior_schema_require(
    identical(names(achieved), c("K", "weight")) &&
      identical(names(residuals), "K") &&
      identical(achieved[["K", exact = TRUE]][["M", exact = TRUE]], NULL) &&
      identical(achieved[["K", exact = TRUE]][["source", exact = TRUE]],
                "legacy_serialized_selected_values") &&
      identical(achieved[["weight", exact = TRUE]][["metric", exact = TRUE]],
                target_weight[["metric", exact = TRUE]]) &&
      identical(achieved[["weight", exact = TRUE]][["source", exact = TRUE]],
                "legacy_serialized_selected_values") &&
      identical(residuals[["K", exact = TRUE]][["source", exact = TRUE]],
                "recomputed_from_legacy_stored_values") &&
      .dpprior_legacy_numeric_close(
        residuals[["K", exact = TRUE]][["mean", exact = TRUE]],
        achieved[["K", exact = TRUE]][["mean", exact = TRUE]] -
          target_K[["implied", exact = TRUE]][["mean", exact = TRUE]]
      ) &&
      .dpprior_legacy_numeric_close(
        residuals[["K", exact = TRUE]][["variance", exact = TRUE]],
        achieved[["K", exact = TRUE]][["variance", exact = TRUE]] -
          target_K[["implied", exact = TRUE]][["variance", exact = TRUE]]
      ),
    "legacy_migration_retained_truth", "result.achieved",
    "stored selected evidence with exact source markers and residual arithmetic",
    list(achieved = achieved, residuals = residuals)
  )

  views <- compatibility[["views", exact = TRUE]]
  source <- views[["source", exact = TRUE]]
  expected_view_names <- if (
    "fixed_candidate_recomputation" %in% names(views)
  ) c("source", "fixed_candidate_recomputation") else "source"
  .dpprior_schema_require(
    identical(compatibility[["top_level_aliases", exact = TRUE]],
              character()) &&
      identical(names(views), expected_view_names) &&
      identical(names(compatibility[["deprecations", exact = TRUE]]),
                "legacy_schema") &&
      identical(source[["source_schema", exact = TRUE]],
                "DPprior/1.1/fit") &&
      identical(source[["source_package_version", exact = TRUE]], "1.1.0") &&
      identical(source[["source_class", exact = TRUE]], "DPprior_fit") &&
      identical(source[["public_candidate_quarantined", exact = TRUE]],
                FALSE) &&
      identical(source[["required_action", exact = TRUE]],
                "refit_with_current_API"),
    "legacy_migration_compatibility", "result.compatibility",
    "the exact source/audit view boundary and public-candidate migration policy",
    compatibility
  )
  invisible(TRUE)
}


.dpprior_validate_dual_legacy_invariants <- function(raw) {
  path <- "result"
  .dpprior_schema_require(
    identical(raw[["object_type", exact = TRUE]], "fit") &&
      identical(raw[["mode", exact = TRUE]], "dual_legacy") &&
      !is.null(raw[["parameters", exact = TRUE]]),
    "legacy_route", path,
    "one finite canonical dual_legacy fit", raw[c("object_type", "mode")]
  )

  contract <- "DPprior_dual_v1_path_scaled_soft_equality_loss"
  warning_code <- "deprecated_legacy_dual_anchor"
  legacy <- raw[["legacy", exact = TRUE]]
  provenance <- raw[["provenance", exact = TRUE]]
  computation <- raw[["computation", exact = TRUE]]
  verification <- raw[["verification", exact = TRUE]]
  parameters <- raw[["parameters", exact = TRUE]]
  parameterization <- parameters[["parameterization", exact = TRUE]]
  target_K <- raw[["target", exact = TRUE]][["K", exact = TRUE]]
  target_weight <- raw[["target", exact = TRUE]][["weight", exact = TRUE]]
  # The enclosing result validator authenticates both canonical target class
  # vectors before reaching this route.  Strip them again defensively here so
  # no scientific read below can dispatch a hostile registered `$`/`[[` method.
  if (is.object(target_K)) target_K <- unclass(target_K)
  if (is.object(target_weight)) target_weight <- unclass(target_weight)
  J <- raw[["J", exact = TRUE]]

  .dpprior_schema_require(
    identical(raw[["method", exact = TRUE]], "dual-anchor") &&
      identical(raw[["status", exact = TRUE]], "approximate") &&
      isTRUE(raw[["usable", exact = TRUE]]) &&
      identical(raw[["verified", exact = TRUE]], FALSE) &&
      parameterization %in% c("log(shape), log(rate)", "log_ab"),
    "legacy_literal_contract", path,
    paste(
      "dual-anchor/approximate/usable/unverified with one retained canonical",
      "log shape-rate parameterization"
    ),
    raw[c("method", "status", "usable", "verified", "parameters")]
  )
  .dpprior_schema_require(
    identical(legacy[["contract", exact = TRUE]], contract) &&
      identical(legacy[["approximation_opt_in", exact = TRUE]], TRUE) &&
      identical(legacy[["warning_code", exact = TRUE]], warning_code) &&
      identical(provenance[["requested_method", exact = TRUE]],
                "dual-anchor") &&
      identical(provenance[["selected_method", exact = TRUE]],
                "dual-anchor") &&
      identical(provenance[["parameterization", exact = TRUE]],
                parameterization) &&
      identical(provenance[["input_fit", exact = TRUE]], NULL) &&
      identical(provenance[["approximation", exact = TRUE]], list(
        active = TRUE, opt_in = TRUE, kind = "legacy_soft_equality_loss",
        warning_code = warning_code
      )) &&
      identical(provenance[["legacy", exact = TRUE]], list(
        active = TRUE, contract = contract,
        deprecation_stage = "compatibility_window"
      )),
    "legacy_provenance", "result.provenance",
    "the exact native legacy approximation and retention contract",
    list(legacy = legacy, provenance = provenance)
  )
  .dpprior_schema_require(
    identical(provenance[["backend", exact = TRUE]][[
      "package", exact = TRUE
    ]], "DPprior") &&
      identical(provenance[["backend", exact = TRUE]][[
        "implementation", exact = TRUE
      ]], "R/15_dual_anchor.R:DPprior_dual") &&
      is.null(provenance[["backend", exact = TRUE]][[
        "source_commit", exact = TRUE
      ]]) &&
      identical(provenance[["migration", exact = TRUE]], list(
      source_schema = "native", adapter = "none", lossless = TRUE,
      missing_evidence = character(), warnings = character()
    )),
    "legacy_provenance", "result.provenance",
    paste(
      "the exact native package/implementation and lossless non-migration",
      "provenance; package_version remains retained release metadata"
    ),
    list(backend = provenance$backend, migration = provenance$migration)
  )

  settings <- computation[["used", exact = TRUE]]
  .dpprior_schema_require(
    identical(computation[["request", exact = TRUE]], settings) &&
      identical(names(settings[["controls", exact = TRUE]]),
                c("lambda", "loss_type", "max_iter", "M")) &&
      identical(settings[["method", exact = TRUE]], "dual-anchor") &&
      identical(settings[["parameterization", exact = TRUE]],
                parameterization),
    "legacy_settings", "result.computation",
    "byte-identical requested/used dual-anchor settings in canonical order",
    list(request = computation$request, used = settings)
  )
  controls <- settings[["controls", exact = TRUE]]
  lambda <- controls[["lambda", exact = TRUE]]
  loss_type <- controls[["loss_type", exact = TRUE]]
  max_iter <- controls[["max_iter", exact = TRUE]]
  M <- controls[["M", exact = TRUE]]
  .dpprior_schema_require(
    .dpprior_schema_is_finite_scalar(lambda) && lambda >= 0 && lambda <= 1 &&
      loss_type %in% c("relative", "adaptive", "absolute") &&
      .dpprior_schema_is_count(max_iter, 1L) &&
      max_iter <= floor(.Machine$integer.max / 2) &&
      .dpprior_schema_is_count(M, 1L) && M <= .QUADRATURE_MAX_NODES &&
      identical(legacy[["lambda", exact = TRUE]], lambda),
    "legacy_controls", "result.computation.used.controls",
    paste(
      "the exact typed lambda/loss_type/max_iter/M control record with",
      "overflow-safe max_iter"
    ),
    controls
  )
  orders <- computation[["orders", exact = TRUE]]
  .dpprior_schema_require(
    identical(orders, list(
      M_requested = M, M_selected = M,
      M_verification_required = NULL, M_verification_used = NULL,
      requested_reason = "public_argument",
      selected_reason = "legacy_selected_order",
      verification_required_reason =
        "legacy_adapter_did_not_request_independent_order",
      verification_used_reason =
        "legacy_adapter_did_not_perform_independent_order"
    )),
    "legacy_orders", "result.computation.orders",
    "the exact selected-only legacy quadrature-order record", orders
  )

  .dpprior_schema_require(
    identical(names(computation[["resources", exact = TRUE]]), c(
      "K_only_baseline", "optimizer_controls", "legacy_weight_request",
      "objective"
    )),
    "legacy_resources", "result.computation.resources",
    "the exact ordered legacy resource vocabulary",
    names(computation[["resources", exact = TRUE]])
  )
  resources <- computation[["resources", exact = TRUE]]
  .dpprior_schema_require(
    identical(resources[["optimizer_controls", exact = TRUE]], list(
      adaptive = list(maxit = 50L),
      primary = list(maxit = max_iter),
      fallback = list(maxit = 2L * max_iter)
    )),
    "legacy_optimizer_controls",
    "result.computation.resources.optimizer_controls",
    "the exact fixed scaling and normalized primary/fallback maxit controls",
    resources[["optimizer_controls", exact = TRUE]]
  )

  authority <- list(
    metric = target_weight[["metric", exact = TRUE]],
    relation = target_weight[["relation", exact = TRUE]],
    value = target_weight[["value", exact = TRUE]],
    threshold = target_weight[["threshold", exact = TRUE]],
    probability = target_weight[["probability", exact = TRUE]]
  )
  expected_estimand <- switch(
    authority$metric,
    wsb_tail = "P(W_SB > threshold)",
    wsb_mean = "E(W_SB)",
    wsb_quantile = "Q_probability(W_SB)",
    NULL
  )
  expected_legacy_request <- switch(
    authority$metric,
    wsb_tail = list(prob = list(
      threshold = authority$threshold, value = authority$value
    )),
    wsb_mean = list(mean = authority$value),
    wsb_quantile = list(quantile = list(
      prob = authority$probability, value = authority$value
    )),
    NULL
  )
  .dpprior_schema_require(
    !is.null(expected_estimand) &&
      identical(target_weight[["request", exact = TRUE]], authority) &&
      identical(target_weight[["normalized", exact = TRUE]], authority) &&
      identical(target_weight[["used", exact = TRUE]], authority) &&
      identical(target_weight[["relation", exact = TRUE]], "target") &&
      identical(target_weight[["operator", exact = TRUE]], "target") &&
      identical(target_weight[["estimand", exact = TRUE]],
                expected_estimand) &&
      identical(target_weight[["units", exact = TRUE]], "probability") &&
      identical(target_weight[["certification", exact = TRUE]], list()) &&
      identical(target_weight[["provenance", exact = TRUE]], list(
        source = "R/15_dual_anchor.R:legacy_weight_target",
        transformation = NULL, selection = NULL,
        legacy_request = expected_legacy_request
      )) &&
      identical(resources[["legacy_weight_request", exact = TRUE]],
                expected_legacy_request),
    "legacy_weight_target", "result.target.weight",
    paste(
      "one exact direct canonical weight target and its retained safe legacy",
      "request"
    ),
    list(target = target_weight,
         resource = resources[["legacy_weight_request", exact = TRUE]])
  )

  baseline <- resources[["K_only_baseline", exact = TRUE]]
  .dpprior_schema_exact_names(
    baseline,
    c(
      "schema", "mode", "method", "J", "status", "usable", "verified",
      "parameters", "target", "selected_snapshot"
    ),
    "result.computation.resources.K_only_baseline"
  )
  expected_target_reference <- list(
    schema = target_K[["schema", exact = TRUE]],
    kind = target_K[["kind", exact = TRUE]],
    J = target_K[["J", exact = TRUE]],
    used = target_K[["used", exact = TRUE]],
    implied = target_K[["implied", exact = TRUE]]
  )
  baseline_mode <- baseline[["mode", exact = TRUE]]
  expected_parameterization <- switch(
    baseline_mode,
    a2_moment = "log(shape), log(rate)",
    a2_kl = "log_ab",
    NULL
  )
  baseline_method_ok <- switch(
    baseline_mode,
    a2_moment = baseline[["method", exact = TRUE]] %in%
      c("A2-MN", "A2-MN+NM"),
    a2_kl = identical(baseline[["method", exact = TRUE]], "A2-KL"),
    FALSE
  )
  .dpprior_schema_require(
    identical(baseline[["schema", exact = TRUE]], "dpprior.result/1") &&
      baseline_method_ok && identical(baseline[["J", exact = TRUE]], J) &&
      !is.null(expected_parameterization) &&
      identical(parameterization, expected_parameterization) &&
      baseline[["status", exact = TRUE]] %in% c("converged", "boundary") &&
      isTRUE(baseline[["usable", exact = TRUE]]) &&
      isTRUE(baseline[["verified", exact = TRUE]]) &&
      identical(baseline[["parameters", exact = TRUE]][[
        "parameterization", exact = TRUE
      ]], parameterization) &&
      identical(baseline[["target", exact = TRUE]],
                expected_target_reference),
    "legacy_baseline", "result.computation.resources.K_only_baseline",
    "a decision-ready A2 baseline bound to the exact public K target",
    baseline
  )
  .dpprior_validate_parameters(
    baseline[["parameters", exact = TRUE]],
    "result.computation.resources.K_only_baseline.parameters"
  )
  baseline_snapshot <- baseline[["selected_snapshot", exact = TRUE]]
  .dpprior_schema_exact_names(
    baseline_snapshot,
    c("parameters", "M", "achieved_K", "finite", "source"),
    "result.computation.resources.K_only_baseline.selected_snapshot"
  )
  .dpprior_validate_achieved_K(
    baseline_snapshot[["achieved_K", exact = TRUE]], J,
    paste0(
      "result.computation.resources.K_only_baseline.",
      "selected_snapshot.achieved_K"
    )
  )
  baseline_M <- baseline_snapshot[["M", exact = TRUE]]
  baseline_moments <- tryCatch(
    exact_K_moments(
      J, baseline$parameters$a, baseline$parameters$b, M = baseline_M
    ),
    error = identity
  )
  baseline_K <- baseline_snapshot[["achieved_K", exact = TRUE]]
  .dpprior_schema_require(
    identical(baseline_snapshot[["parameters", exact = TRUE]],
              baseline[["parameters", exact = TRUE]]) &&
      .dpprior_schema_is_count(baseline_M, 1L) &&
      isTRUE(baseline_snapshot[["finite", exact = TRUE]]) &&
      identical(baseline_snapshot[["source", exact = TRUE]],
                "selected_order") &&
      identical(baseline_K[["M", exact = TRUE]], baseline_M) &&
      identical(baseline_K[["source", exact = TRUE]], "selected_order") &&
      !inherits(baseline_moments, "condition") &&
      .dpprior_legacy_numeric_close(
        baseline_K[["mean", exact = TRUE]],
        baseline_moments[["mean", exact = TRUE]]
      ) && .dpprior_legacy_numeric_close(
        baseline_K[["variance", exact = TRUE]],
        baseline_moments[["var", exact = TRUE]]
      ),
    "legacy_baseline_truth",
    "result.computation.resources.K_only_baseline.selected_snapshot",
    "the exact retained baseline parameters/order and freshly recomputed K",
    list(snapshot = baseline_snapshot, recomputed = baseline_moments)
  )
  if ("pmf" %in% names(baseline_K)) {
    baseline_verifier_M <- min(
      .QUADRATURE_MAX_NODES,
      max(2L * baseline_M, baseline_M + 40L)
    )
    fresh_baseline_pmf <- tryCatch(
      .get_K_pmf_support(
        J, baseline$parameters$a, baseline$parameters$b,
        M = baseline_M, M_verify = baseline_verifier_M,
        abs_tol = 1e-10, rel_tol = 1e-8
      )[["pmf", exact = TRUE]],
      error = identity
    )
    .dpprior_schema_require(
      !inherits(fresh_baseline_pmf, "condition") &&
        is.numeric(fresh_baseline_pmf) &&
        length(fresh_baseline_pmf) == J &&
        max(abs(
          unname(as.numeric(baseline_K[["pmf", exact = TRUE]])) -
            unname(as.numeric(fresh_baseline_pmf))
        )) <= .TOL_PMF_SUM,
      "legacy_baseline_PMF_truth",
      paste0(
        "result.computation.resources.K_only_baseline.",
        "selected_snapshot.achieved_K.pmf"
      ),
      "the freshly recomputed selected-order marginal K PMF",
      list(recorded = baseline_K[["pmf", exact = TRUE]],
           recomputed = fresh_baseline_pmf)
    )
  }
  if (identical(baseline_mode, "a2_moment")) {
    target_implied <- target_K[["implied", exact = TRUE]]
    selected_tolerance <- c(
      mean = 1e-8 + 1e-8 * max(abs(target_implied$mean), 1),
      variance = 1e-8 + 1e-8 * max(abs(target_implied$variance), 1)
    )
    .dpprior_schema_require(
      abs(baseline_K$mean - target_implied$mean) <=
        selected_tolerance[["mean"]] &&
        abs(baseline_K$variance - target_implied$variance) <=
          selected_tolerance[["variance"]],
      "legacy_baseline_target_truth",
      "result.computation.resources.K_only_baseline.target",
      "the canonical A2-MN selected-order target adequacy gates",
      list(
        selected = baseline_K[c("mean", "variance")],
        target = target_implied, tolerance = selected_tolerance
      )
    )
  } else {
    objective_pmf <- .dpprior_result_A2_KL_objective_pmf(
      target_K, J,
      "result.computation.resources.K_only_baseline.target"
    )
    selected_pmf <- baseline_K[["pmf", exact = TRUE]]
    .dpprior_schema_require(
      !is.null(selected_pmf) && length(selected_pmf) == J &&
        all(selected_pmf[objective_pmf > 0] > 0),
      "legacy_baseline_target_truth",
      paste0(
        "result.computation.resources.K_only_baseline.",
        "selected_snapshot.achieved_K.pmf"
      ),
      "positive selected mass wherever the canonical A2-KL objective has mass",
      selected_pmf
    )
    objective_moments <- .dpprior_target_pmf_moments(objective_pmf)
    positive <- objective_pmf > 0
    adequacy <- c(
      kl = sum(objective_pmf[positive] * log(
        objective_pmf[positive] / selected_pmf[positive]
      )),
      l1 = sum(abs(objective_pmf - selected_pmf)),
      mean_scaled = abs(
        baseline_K$mean - objective_moments[["mean"]]
      ) / max(1, sqrt(objective_moments[["variance"]])),
      variance_scaled = abs(
        baseline_K$variance - objective_moments[["variance"]]
      ) / max(1, objective_moments[["variance"]])
    )
    adequacy_limits <- c(
      kl = 0.015, l1 = 0.11, mean_scaled = 0.01,
      variance_scaled = 0.065
    )
    .dpprior_schema_require(
      all(is.finite(adequacy)) && all(adequacy <= adequacy_limits),
      "legacy_baseline_target_truth",
      "result.computation.resources.K_only_baseline.target",
      "all four canonical A2-KL selected-order adequacy gates",
      list(values = adequacy, tolerances = adequacy_limits)
    )
  }

  scaling <- computation[["scaling", exact = TRUE]]
  expected_formula <- switch(
    loss_type,
    relative = "legacy_relative_squared_loss",
    adaptive = "legacy_path_derived_adaptive_loss",
    absolute = "legacy_absolute_squared_loss"
  )
  .dpprior_schema_require(
    identical(scaling[["requested", exact = TRUE]],
              list(loss_type = loss_type)) &&
      identical(scaling[["used", exact = TRUE]],
                list(loss_type = loss_type)) &&
      identical(scaling[["formula", exact = TRUE]], expected_formula) &&
      identical(scaling[["fixed_from_input", exact = TRUE]], FALSE) &&
      identical(scaling[["change_reason", exact = TRUE]], ""),
    "legacy_scaling", "result.computation.scaling",
    "the exact loss-type-specific legacy scaling policy", scaling
  )
  scale_values <- scaling[["values", exact = TRUE]]
  if (identical(loss_type, "adaptive")) {
    .dpprior_schema_require(
      identical(names(scale_values), c("L_K_scale", "L_w_scale")) &&
        (if (identical(lambda, 1)) {
          is.null(scale_values[["L_K_scale", exact = TRUE]]) &&
            is.null(scale_values[["L_w_scale", exact = TRUE]])
        } else {
          .dpprior_schema_is_finite_scalar(
            scale_values[["L_K_scale", exact = TRUE]]
          ) && scale_values[["L_K_scale", exact = TRUE]] > 0 &&
            .dpprior_schema_is_finite_scalar(
              scale_values[["L_w_scale", exact = TRUE]]
            ) && scale_values[["L_w_scale", exact = TRUE]] > 0
        }),
      "legacy_scaling_values", "result.computation.scaling.values",
      paste(
        "positive adaptive scales except the exact named NULL endpoint",
        "short-circuit"
      ), scale_values
    )
  } else {
    .dpprior_schema_require(
      identical(scale_values, list()),
      "legacy_scaling_values", "result.computation.scaling.values",
      "an empty scaling record for non-adaptive legacy losses", scale_values
    )
  }

  .dpprior_schema_exact_names(
    legacy[["losses", exact = TRUE]],
    c("loss_type", "K_loss", "weight_loss", "total_loss", "scaling"),
    "result.legacy.losses"
  )
  objective_record <- resources[["objective", exact = TRUE]]
  .dpprior_schema_exact_names(
    objective_record,
    c("lambda", "loss_type", "K_loss", "weight_loss", "total_loss"),
    "result.computation.resources.objective"
  )
  .dpprior_schema_require(
    identical(legacy$losses$loss_type, loss_type) &&
      identical(legacy$losses$scaling, scale_values) &&
      identical(objective_record$lambda, lambda) &&
      identical(objective_record$loss_type, loss_type) &&
      identical(legacy$losses[c("K_loss", "weight_loss", "total_loss")],
                objective_record[c("K_loss", "weight_loss", "total_loss")]),
    "legacy_loss_mirrors", "result.legacy.losses",
    "exact identity across controls, scaling, legacy losses, and resources",
    list(legacy = legacy$losses, objective = objective_record)
  )

  truth <- .dpprior_legacy_loss_truth(
    parameters, J, M, target_K, target_weight, loss_type, scale_values,
    "result"
  )
  achieved_K <- raw[["achieved", exact = TRUE]][["K", exact = TRUE]]
  achieved_weight <- raw[["achieved", exact = TRUE]][[
    "weight", exact = TRUE
  ]]
  K_residual <- raw[["residuals", exact = TRUE]][["K", exact = TRUE]]
  weight_residual <- raw[["residuals", exact = TRUE]][[
    "weight", exact = TRUE
  ]]
  .dpprior_schema_require(
    identical(names(raw[["achieved", exact = TRUE]]), c("K", "weight")) &&
      identical(names(raw[["residuals", exact = TRUE]]), c("K", "weight")) &&
      identical(raw[["tolerances", exact = TRUE]], list()) &&
      identical(names(achieved_K),
                c("mean", "variance", "estimand", "source", "M")) &&
      identical(achieved_K[["estimand", exact = TRUE]], "K_J") &&
      identical(achieved_K[["source", exact = TRUE]],
                "legacy_selected_order") &&
      identical(achieved_K[["M", exact = TRUE]], M) &&
      identical(achieved_weight[["metric", exact = TRUE]],
                target_weight[["metric", exact = TRUE]]) &&
      identical(achieved_weight[["source", exact = TRUE]],
                "legacy_selected_order") &&
      identical(names(K_residual), c("mean", "variance")) &&
      identical(names(weight_residual), c("raw", "directed")) &&
      .dpprior_legacy_numeric_close(achieved_K$mean, truth$moments$mean) &&
      .dpprior_legacy_numeric_close(achieved_K$variance,
                                    truth$moments$var) &&
      .dpprior_legacy_numeric_close(achieved_weight$value,
                                    truth$weight_value) &&
      .dpprior_legacy_numeric_close(K_residual$mean,
                                    truth$K_residual[["mean"]]) &&
      .dpprior_legacy_numeric_close(K_residual$variance,
                                    truth$K_residual[["variance"]]) &&
      .dpprior_legacy_numeric_close(weight_residual$raw,
                                    truth$weight_residual) &&
      .dpprior_legacy_numeric_close(weight_residual$directed,
                                    truth$weight_residual),
    "legacy_fresh_truth", "result.achieved",
    paste(
      "public/snapshot K and weight values plus residuals freshly",
      "recomputed from J/a/b/M and canonical targets"
    ),
    list(
      achieved = raw$achieved, residuals = raw$residuals,
      recomputed = truth
    )
  )
  expected_total <- lambda * truth$K_loss +
    (1 - lambda) * truth$weight_loss
  .dpprior_schema_require(
    .dpprior_legacy_numeric_close(objective_record$K_loss, truth$K_loss) &&
      .dpprior_legacy_numeric_close(
        objective_record$weight_loss, truth$weight_loss
      ) && .dpprior_legacy_numeric_close(
        objective_record$total_loss, expected_total
      ),
    "legacy_objective_truth", "result.computation.resources.objective",
    "fresh exact legacy K/weight/total loss arithmetic", list(
      recorded = objective_record,
      recomputed = list(
        K_loss = truth$K_loss, weight_loss = truth$weight_loss,
        total_loss = expected_total
      )
    )
  )

  .dpprior_schema_require(
    identical(verification[["method", exact = TRUE]], "legacy_unverified") &&
      identical(verification[["performed", exact = TRUE]], FALSE) &&
      identical(verification[["passed", exact = TRUE]], FALSE) &&
      identical(verification[["reason", exact = TRUE]],
                "legacy adapter retains no independent verification") &&
      identical(verification[["selected_snapshot", exact = TRUE]], list(
        parameters = parameters, M = M, achieved = raw$achieved,
        residuals = raw$residuals, tolerances = list(), finite = TRUE,
        source = "legacy_selected_order"
      )) &&
      is.null(verification[["verifier_snapshot", exact = TRUE]]) &&
      is.null(verification[["stability", exact = TRUE]]) &&
      identical(verification[["settings", exact = TRUE]], list()) &&
      identical(verification[["components", exact = TRUE]], list()) &&
      identical(verification[["invariants", exact = TRUE]], list()),
    "legacy_verification", "result.verification",
    "the exact selected-only unverified legacy evidence contract",
    verification
  )

  attempts <- computation[["attempts", exact = TRUE]]
  .dpprior_schema_require(
    identical(computation[["candidate_evaluations", exact = TRUE]], list()) &&
      is.null(computation[["selected_candidate_id", exact = TRUE]]) &&
      is.null(computation[["trace", exact = TRUE]]),
    "legacy_computation", "result.computation",
    "no candidate-evaluation ledger or trace for the legacy adapter",
    computation[c("candidate_evaluations", "selected_candidate_id", "trace")]
  )

  baseline_relative <- .dpprior_legacy_loss_truth(
    baseline$parameters, J, M, target_K, target_weight, "relative", list(),
    "result.computation.resources.K_only_baseline"
  )
  if (identical(lambda, 1)) {
    .dpprior_schema_require(
      identical(attempts, list()) &&
        is.null(computation[["selected_attempt_id", exact = TRUE]]) &&
        identical(computation[["fallback", exact = TRUE]], list(
          attempted = FALSE, used = FALSE, trigger_attempt_id = NULL,
          selected_attempt_id = NULL, reason_code = NULL, message = "",
          outcome = "not_attempted"
        )) && identical(computation[["termination", exact = TRUE]], list(
          code = "deterministic",
          message = "lambda=1 retained the K-only parameters",
          source = "legacy_adapter", iterations = 0L,
          boundary_reason = NULL
        )) && identical(parameters, baseline$parameters) &&
        identical(provenance[["is_fallback", exact = TRUE]], FALSE),
      "legacy_endpoint", "result.computation",
      paste(
        "the exact attempt-free lambda=1 K-only endpoint with baseline",
        "parameter identity"
      ), computation
    )
  } else {
    adaptive <- identical(loss_type, "adaptive")
    stages <- vapply(attempts, `[[`, character(1), "stage")
    methods <- vapply(attempts, `[[`, character(1), "method")
    expected_prefix <- if (adaptive) c("scaling", "primary") else "primary"
    expected_methods <- if (adaptive) c("BFGS-weight-scaling", "BFGS") else
      "BFGS"
    fallback_record <- computation[["fallback", exact = TRUE]]
    fallback_attempted <- fallback_record[["attempted", exact = TRUE]]
    expected_stages <- c(expected_prefix,
                         if (fallback_attempted) "fallback" else character())
    expected_methods <- c(expected_methods,
                          if (fallback_attempted) "Nelder-Mead" else
                            character())
    expected_ids <- c(
      if (adaptive) "attempt-scaling-001" else character(),
      "attempt-primary-001",
      if (fallback_attempted) "attempt-fallback-001" else character()
    )
    ids <- vapply(attempts, `[[`, character(1), "id")
    .dpprior_schema_require(
      identical(stages, expected_stages) &&
        identical(methods, expected_methods) &&
        identical(ids, expected_ids),
      "legacy_attempt_route", "result.computation.attempts",
      "the exact ordered adaptive?/primary/fallback producer attempt route",
      list(ids = ids, stages = stages, methods = methods)
    )
    expected_start <- c(
      log_shape = log(baseline$parameters$a),
      log_rate = log(baseline$parameters$b)
    )
    # Optimizer message/warning text is retained descriptive evidence only.
    # Selection is derived below from typed availability, objective, exit,
    # stage, and fallback fields; the scaling-stage exit is likewise not a
    # convergence claim and cannot become the selected public execution.
    for (i in seq_along(attempts)) {
      attempt <- attempts[[i]]
      stage <- attempt[["stage", exact = TRUE]]
      expected_control <- list(maxit = if (identical(stage, "scaling")) {
        50L
      } else if (identical(stage, "fallback")) {
        2L * max_iter
      } else {
        max_iter
      })
      .dpprior_schema_require(
        identical(attempt[["start", exact = TRUE]], expected_start) &&
          is.null(attempt[["bounds", exact = TRUE]]) &&
          identical(attempt[["control", exact = TRUE]], expected_control) &&
          identical(
            attempt[["iterations", exact = TRUE]],
            if (is.null(attempt[["evaluations", exact = TRUE]]) ||
                !"function_count" %in%
                  names(attempt[["evaluations", exact = TRUE]])) {
              NULL
            } else {
              attempt[["evaluations", exact = TRUE]][[
                "function_count", exact = TRUE
              ]]
            }
          ) &&
          is.null(attempt[["elapsed_seconds", exact = TRUE]]) &&
          identical(
            attempt[["unavailable", exact = TRUE]][[
              "elapsed_seconds", exact = TRUE
            ]],
            "legacy optimizer did not retain elapsed seconds"
          ),
        "legacy_attempt_control",
        sprintf("result.computation.attempts[[%d]]", i),
        paste(
          "the exact baseline start, absent bounds, stage-specific maxit,",
          paste(
            "iterations derived from retained function evaluations, and an",
            "explicitly unavailable elapsed time"
          )
        ),
        attempt[c(
          "start", "bounds", "control", "iterations", "evaluations",
          "elapsed_seconds", "unavailable"
        )]
      )
      candidate <- attempt[["candidate_parameters", exact = TRUE]]
      candidate_objective <- attempt[["candidate_objective", exact = TRUE]]
      if (!is.null(candidate)) {
        .dpprior_schema_require(
          identical(candidate[["parameterization", exact = TRUE]],
                    parameterization) && !is.null(candidate_objective),
          "legacy_attempt_candidate",
          sprintf("result.computation.attempts[[%d]]", i),
          "a complete parameter/objective pair in the fixed parameterization",
          attempt[c("candidate_parameters", "candidate_objective")]
        )
        candidate_truth <- .dpprior_legacy_loss_truth(
          candidate, J, M, target_K, target_weight,
          if (identical(stage, "scaling")) "relative" else loss_type,
          if (identical(stage, "scaling")) list() else scale_values,
          sprintf("result.computation.attempts[[%d]]", i)
        )
        expected_attempt_objective <- if (identical(stage, "scaling")) {
          candidate_truth$weight_loss
        } else {
          lambda * candidate_truth$K_loss +
            (1 - lambda) * candidate_truth$weight_loss
        }
        .dpprior_schema_require(
          .dpprior_legacy_numeric_close(
            candidate_objective, expected_attempt_objective
          ),
          "legacy_attempt_objective",
          sprintf("result.computation.attempts[[%d]].candidate_objective", i),
          paste(
            "fresh relative weight-only loss for scaling or the declared",
            "fresh total objective for primary/fallback"
          ),
          list(recorded = candidate_objective,
               recomputed = expected_attempt_objective)
        )
      }
    }

    scaling_attempt <- if (adaptive) attempts[[1L]] else NULL
    primary_index <- if (adaptive) 2L else 1L
    primary <- attempts[[primary_index]]
    primary_valid <- is.null(primary[["error", exact = TRUE]]) &&
      !is.null(primary[["candidate_parameters", exact = TRUE]]) &&
      !is.null(primary[["candidate_objective", exact = TRUE]])
    expected_fallback_attempted <- !primary_valid ||
      !identical(primary[["exit_code", exact = TRUE]], 0L)
    fallback_attempt <- if (fallback_attempted) attempts[[length(attempts)]] else
      NULL
    fallback_valid <- !is.null(fallback_attempt) &&
      is.null(fallback_attempt[["error", exact = TRUE]]) &&
      !is.null(fallback_attempt[["candidate_parameters", exact = TRUE]]) &&
      !is.null(fallback_attempt[["candidate_objective", exact = TRUE]])
    expected_fallback_used <- fallback_attempted && fallback_valid &&
      (!primary_valid || fallback_attempt$candidate_objective <
         primary$candidate_objective)
    expected_selected_id <- if (expected_fallback_used) {
      "attempt-fallback-001"
    } else {
      "attempt-primary-001"
    }
    .dpprior_schema_require(
      identical(fallback_attempted, expected_fallback_attempted) &&
        identical(fallback_record[["used", exact = TRUE]],
                  expected_fallback_used) &&
        identical(computation[["selected_attempt_id", exact = TRUE]],
                  expected_selected_id) &&
        identical(provenance[["is_fallback", exact = TRUE]],
                  expected_fallback_used),
      "legacy_fallback_truth", "result.computation.fallback",
      "fallback attempt/use and selected ID derived from retained attempts",
      list(fallback = fallback_record, selected = computation$selected_attempt_id)
    )
    if (fallback_attempted) {
      .dpprior_schema_require(
        identical(fallback_record[["trigger_attempt_id", exact = TRUE]],
                  "attempt-primary-001") &&
          identical(fallback_record[["selected_attempt_id", exact = TRUE]],
                    if (expected_fallback_used) expected_selected_id else NULL) &&
          identical(fallback_record[["reason_code", exact = TRUE]],
                    "primary_optimizer_exit_nonzero") &&
          identical(fallback_record[["message", exact = TRUE]],
                    "fallback was attempted after the primary optimizer") &&
          identical(fallback_record[["outcome", exact = TRUE]],
                    if (expected_fallback_used) "selected" else
                      "attempted_not_selected"),
        "legacy_fallback_truth", "result.computation.fallback",
        "the exact primary-triggered fallback provenance record",
        fallback_record
      )
    } else {
      .dpprior_schema_require(
        identical(fallback_record, list(
          attempted = FALSE, used = FALSE, trigger_attempt_id = NULL,
          selected_attempt_id = NULL, reason_code = NULL, message = "",
          outcome = "not_attempted"
        )),
        "legacy_fallback_truth", "result.computation.fallback",
        "the exact not-attempted fallback record", fallback_record
      )
    }
    selected_attempt <- attempts[[match(expected_selected_id, ids)]]
    expected_termination_code <- if (identical(
      selected_attempt[["exit_code", exact = TRUE]], 0L
    )) "selected" else "approximate"
    .dpprior_schema_require(
      identical(selected_attempt[["candidate_parameters", exact = TRUE]],
                parameters) &&
        .dpprior_legacy_numeric_close(
          selected_attempt[["candidate_objective", exact = TRUE]],
          objective_record[["total_loss", exact = TRUE]]
        ) && identical(computation[["termination", exact = TRUE]], list(
          code = expected_termination_code,
          message = paste(
            "legacy equality-loss candidate retained; optimizer exit is",
            "descriptive evidence only"
          ),
          source = if (expected_fallback_used) "fallback_optimizer" else
            "optimizer",
          iterations = selected_attempt[["iterations", exact = TRUE]],
          boundary_reason = NULL
        )),
      "legacy_selected_attempt", "result.computation",
      paste(
        "public parameters/objective and termination derived from the exact",
        "selected producer attempt"
      ),
      list(selected = selected_attempt, termination = computation$termination)
    )

    if (adaptive) {
      scaling_valid <- is.null(scaling_attempt[["error", exact = TRUE]]) &&
        !is.null(scaling_attempt[["candidate_parameters", exact = TRUE]]) &&
        !is.null(scaling_attempt[["candidate_objective", exact = TRUE]])
      scale_parameters <- if (scaling_valid) {
        scaling_attempt[["candidate_parameters", exact = TRUE]]
      } else {
        baseline[["parameters", exact = TRUE]]
      }
      scale_truth <- .dpprior_legacy_loss_truth(
        scale_parameters, J, M, target_K, target_weight, "relative", list(),
        "result.computation.scaling.values"
      )
      expected_scales <- list(
        L_K_scale = max(scale_truth$K_loss, 0.01),
        L_w_scale = max(baseline_relative$weight_loss, 0.001)
      )
      .dpprior_schema_require(
        .dpprior_legacy_numeric_close(
          scale_values$L_K_scale, expected_scales$L_K_scale
        ) && .dpprior_legacy_numeric_close(
          scale_values$L_w_scale, expected_scales$L_w_scale
        ),
        "legacy_scaling_truth", "result.computation.scaling.values",
        "fresh path-derived adaptive K/weight scales", list(
          recorded = scale_values, recomputed = expected_scales
        )
      )
    }
  }

  expected_message <- if (identical(lambda, 1)) {
    paste(
      "legacy lambda=1 returned the K-only parameters; no independent",
      "verification was performed"
    )
  } else {
    paste(
      "legacy equality-loss candidate retained for compatibility; optimizer",
      "exit was not treated as independent verification"
    )
  }
  .dpprior_schema_require(
    identical(raw[["message", exact = TRUE]], expected_message),
    "legacy_message", "result.message",
    "the exact route-specific legacy approximation message",
    raw[["message", exact = TRUE]]
  )
  invisible(TRUE)
}


.dpprior_new_legacy_details <- function(contract,
                                        lambda,
                                        losses = list(),
                                        approximation_opt_in,
                                        warning_code) {
  out <- list(
    contract = contract,
    lambda = lambda,
    losses = losses,
    approximation_opt_in = approximation_opt_in,
    warning_code = warning_code
  )
  .dpprior_validate_legacy_extension(out)
  out
}


.dpprior_validate_diagnostics_extension <- function(
    x, path = "diagnostics", fit_raw = NULL) {
  fit_attached <- !is.null(fit_raw)
  .dpprior_schema_exact_names(
    x,
    c(
      if (fit_attached) "authority" else character(),
      "policy_results", "warnings", "alpha", "K", "weights",
      "coclustering"
    ),
    path
  )
  policy_results <- x[["policy_results", exact = TRUE]]
  .dpprior_schema_require(
    typeof(policy_results) == "list" && is.list(policy_results) &&
      !is.object(policy_results) && is.null(dim(policy_results)) &&
      .dpprior_schema_has_only_attributes(policy_results) &&
      is.null(names(policy_results)),
    "diagnostic_policy_results", paste0(path, ".policy_results"),
    "an ordinary unclassed unnamed ordered policy-result list",
    policy_results
  )
  .dpprior_schema_validate_character_vector(
    x[["warnings", exact = TRUE]], paste0(path, ".warnings")
  )
  component_fields <- list(
    alpha = c("status", "usable", "verified", "mean", "CV"),
    K = c("status", "usable", "verified", "mean", "variance", "pmf", "M"),
    weights = c("status", "usable", "verified", "mean"),
    coclustering = c(
      "status", "usable", "verified", "mean", "variance"
    )
  )
  for (field in .DPPRIOR_DIAGNOSTIC_COMPONENTS) {
    component <- x[[field, exact = TRUE]]
    if (is.null(component)) {
      next
    }
    .dpprior_schema_exact_names(
      component, component_fields[[field]], paste0(path, ".", field)
    )
    .dpprior_schema_validate_scalar_character(
      component[["status", exact = TRUE]], paste0(path, ".", field, ".status")
    )
    .dpprior_schema_require(
      component[["status", exact = TRUE]] %in% c(
        "converged", "boundary", "approximate", "infeasible", "failed"
      ),
      "diagnostic_component_status", paste0(path, ".", field, ".status"),
      "a closed canonical status", component[["status", exact = TRUE]]
    )
    for (flag in c("usable", "verified")) {
      .dpprior_schema_validate_scalar_logical(
        component[[flag, exact = TRUE]], paste0(path, ".", field, ".", flag)
      )
    }
    component_status <- component[["status", exact = TRUE]]
    .dpprior_schema_require(
      if (component_status %in% c("converged", "boundary")) {
        component[["usable", exact = TRUE]] &&
          component[["verified", exact = TRUE]]
      } else if (identical(component_status, "infeasible")) {
        !component[["usable", exact = TRUE]] &&
          component[["verified", exact = TRUE]]
      } else {
        !component[["verified", exact = TRUE]] &&
          (!identical(component_status, "failed") ||
             !component[["usable", exact = TRUE]])
      },
      "diagnostic_component_status", paste0(path, ".", field),
      "status/usable/verified coherence", component
    )
    numeric_fields <- setdiff(component_fields[[field]], c(
      "status", "usable", "verified", "pmf", "M"
    ))
    for (numeric_field in numeric_fields) {
      .dpprior_schema_validate_finite_scalar(
        component[[numeric_field, exact = TRUE]],
        paste0(path, ".", field, ".", numeric_field),
        lower = if (numeric_field %in% c("CV", "variance")) 0 else -Inf
      )
    }
    if (identical(field, "K")) {
      .dpprior_schema_require(
        .dpprior_schema_is_count(component[["M", exact = TRUE]], 1L),
        "diagnostic_K_order", paste0(path, ".K.M"),
        "a positive integer quadrature order", component[["M", exact = TRUE]]
      )
      .dpprior_schema_require(
        is.numeric(component[["pmf", exact = TRUE]]) &&
          !is.object(component[["pmf", exact = TRUE]]) &&
          is.null(dim(component[["pmf", exact = TRUE]])) &&
          .dpprior_schema_has_only_attributes(component[["pmf", exact = TRUE]]) &&
          !anyNA(component[["pmf", exact = TRUE]]) &&
          all(is.finite(component[["pmf", exact = TRUE]])) &&
          all(component[["pmf", exact = TRUE]] >= 0) &&
          abs(sum(component[["pmf", exact = TRUE]]) - 1) <= .TOL_PMF_SUM,
        "diagnostic_K_pmf", paste0(path, ".K.pmf"),
        "an ordinary normalized nonnegative PMF", component[["pmf", exact = TRUE]]
      )
    }
  }
  for (index in seq_along(x[["policy_results", exact = TRUE]])) {
    record <- x[["policy_results", exact = TRUE]][[index]]
    record_path <- sprintf("%s.policy_results[[%d]]", path, index)
    .dpprior_schema_exact_names(
      record,
      c("estimand", "direction", "threshold", "value", "lower", "upper",
        "outcome", "basis"),
      record_path
    )
    for (field in c("estimand", "direction", "outcome", "basis")) {
      .dpprior_schema_validate_scalar_character(
        record[[field, exact = TRUE]], paste0(record_path, ".", field)
      )
    }
    .dpprior_schema_require(
      record[["estimand", exact = TRUE]] %in% c("W_SB", "W_max") &&
        record[["direction", exact = TRUE]] %in% c("above", "below") &&
        record[["outcome", exact = TRUE]] %in%
          c("triggered", "not_triggered", "indeterminate") &&
        record[["basis", exact = TRUE]] %in% c(
          "exact_tail_probability", "verified_numerical_interval",
          "certified_bounds", "backend_unavailable"
        ),
      "diagnostic_policy_vocabulary", record_path,
      "the closed policy estimand/direction/outcome/basis vocabulary", record
    )
    .dpprior_schema_validate_finite_scalar(
      record[["threshold", exact = TRUE]], paste0(record_path, ".threshold"),
      lower = 0, upper = 1
    )
    for (field in c("value", "lower", "upper")) {
      if (!is.null(record[[field, exact = TRUE]])) {
        .dpprior_schema_validate_finite_scalar(
          record[[field, exact = TRUE]], paste0(record_path, ".", field),
          lower = 0, upper = 1
        )
      }
    }
    expected_outcome <- if (identical(record[["basis", exact = TRUE]],
                                      "backend_unavailable")) {
      .dpprior_schema_require(
        all(vapply(record[c("value", "lower", "upper")], is.null, logical(1))),
        "diagnostic_policy_evidence", record_path,
        "no numeric claims for backend_unavailable", record
      )
      "indeterminate"
    } else if (identical(record[["basis", exact = TRUE]],
                         "exact_tail_probability")) {
      .dpprior_schema_require(
        !is.null(record[["value", exact = TRUE]]) &&
          is.null(record[["lower", exact = TRUE]]) &&
          is.null(record[["upper", exact = TRUE]]),
        "diagnostic_policy_evidence", record_path,
        "one exact value and no interval bounds", record
      )
      if (identical(record[["direction", exact = TRUE]], "above")) {
        if (record[["value", exact = TRUE]] > record[["threshold", exact = TRUE]])
          "triggered" else "not_triggered"
      } else if (record[["value", exact = TRUE]] <
                 record[["threshold", exact = TRUE]]) {
        "triggered"
      } else "not_triggered"
    } else {
      .dpprior_schema_require(
        is.null(record[["value", exact = TRUE]]) &&
          !is.null(record[["lower", exact = TRUE]]) &&
          !is.null(record[["upper", exact = TRUE]]) &&
          record[["lower", exact = TRUE]] <= record[["upper", exact = TRUE]],
        "diagnostic_policy_evidence", record_path,
        "ordered lower/upper evidence and no point value", record
      )
      if (identical(record[["direction", exact = TRUE]], "above")) {
        if (record[["lower", exact = TRUE]] > record[["threshold", exact = TRUE]])
          "triggered" else if (record[["upper", exact = TRUE]] <=
                               record[["threshold", exact = TRUE]])
          "not_triggered" else "indeterminate"
      } else if (record[["upper", exact = TRUE]] <
                 record[["threshold", exact = TRUE]]) {
        "triggered"
      } else if (record[["lower", exact = TRUE]] >=
                 record[["threshold", exact = TRUE]]) {
        "not_triggered"
      } else "indeterminate"
    }
    .dpprior_schema_require(
      identical(record[["outcome", exact = TRUE]], expected_outcome),
      "diagnostic_policy_outcome", paste0(record_path, ".outcome"),
      "the outcome recomputed from retained exact/interval evidence", record
    )
  }
  triggered <- sum(vapply(
    x[["policy_results", exact = TRUE]],
    function(record) identical(record[["outcome", exact = TRUE]], "triggered"),
    logical(1)
  ))
  .dpprior_schema_require(
    length(x[["warnings", exact = TRUE]]) == triggered,
    "diagnostic_policy_warnings", paste0(path, ".warnings"),
    "exactly one warning per triggered policy and none otherwise",
    x[["warnings", exact = TRUE]]
  )
  if (fit_attached) {
    .dpprior_validate_fit_diagnostics_truth(x, fit_raw, path)
  }
  invisible(TRUE)
}


.dpprior_validate_fit_diagnostics_truth <- function(
    x, raw, path = "diagnostics") {
  authority <- x[["authority", exact = TRUE]]
  authority_path <- paste0(path, ".authority")
  .dpprior_schema_exact_names(
    authority,
    c(
      "method", "M_selected", "M_verification_required",
      "M_verification_used", "absolute_tolerance", "relative_tolerance",
      "pmf_mass_tolerance", "warning_policy", "allow_approximate"
    ),
    authority_path
  )
  .dpprior_schema_validate_scalar_character(
    authority[["method", exact = TRUE]], paste0(authority_path, ".method")
  )
  .dpprior_schema_require(
    identical(
      authority[["method", exact = TRUE]],
      "fresh_component_specific_diagnostics"
    ),
    "fit_diagnostic_authority", paste0(authority_path, ".method"),
    "fresh_component_specific_diagnostics",
    authority[["method", exact = TRUE]]
  )
  for (field in c(
    "M_selected", "M_verification_required", "M_verification_used"
  )) {
    value <- authority[[field, exact = TRUE]]
    .dpprior_schema_require(
      is.integer(value) && .dpprior_schema_is_count(value, 10L) &&
        value <= .QUADRATURE_MAX_NODES,
      "fit_diagnostic_order", paste0(authority_path, ".", field),
      paste("an ordinary integer quadrature order in [10,",
            .QUADRATURE_MAX_NODES, "]"), value
    )
  }
  for (field in c(
    "absolute_tolerance", "relative_tolerance", "pmf_mass_tolerance"
  )) {
    .dpprior_schema_validate_finite_scalar(
      authority[[field, exact = TRUE]], paste0(authority_path, ".", field),
      lower = 0
    )
  }
  .dpprior_schema_require(
    identical(authority[["absolute_tolerance", exact = TRUE]], 1e-10) &&
      identical(authority[["relative_tolerance", exact = TRUE]], 1e-8) &&
      identical(
        authority[["pmf_mass_tolerance", exact = TRUE]], .TOL_PMF_SUM
      ),
    "fit_diagnostic_tolerance_authority", authority_path,
    paste(
      "the fixed diagnostic absolute=1e-10, relative=1e-8, and",
      "canonical PMF-mass tolerances"
    ),
    authority[c(
      "absolute_tolerance", "relative_tolerance", "pmf_mass_tolerance"
    )]
  )
  .dpprior_schema_validate_scalar_logical(
    authority[["allow_approximate", exact = TRUE]],
    paste0(authority_path, ".allow_approximate")
  )

  warning_policy <- authority[["warning_policy", exact = TRUE]]
  if (!is.null(warning_policy)) {
    .dpprior_schema_exact_names(
      warning_policy,
      c("estimand", "direction", "weight_threshold", "action_threshold"),
      paste0(authority_path, ".warning_policy")
    )
    for (field in c("estimand", "direction")) {
      .dpprior_schema_validate_scalar_character(
        warning_policy[[field, exact = TRUE]],
        paste0(authority_path, ".warning_policy.", field)
      )
    }
    .dpprior_schema_require(
      warning_policy[["estimand", exact = TRUE]] %in% c("W_SB", "W_max") &&
        warning_policy[["direction", exact = TRUE]] %in% c("above", "below"),
      "fit_diagnostic_policy_authority",
      paste0(authority_path, ".warning_policy"),
      "a W_SB/W_max estimand and above/below direction", warning_policy
    )
    .dpprior_schema_validate_finite_scalar(
      warning_policy[["weight_threshold", exact = TRUE]],
      paste0(authority_path, ".warning_policy.weight_threshold"),
      lower = 0, upper = 1, lower_open = TRUE, upper_open = TRUE
    )
    .dpprior_schema_validate_finite_scalar(
      warning_policy[["action_threshold", exact = TRUE]],
      paste0(authority_path, ".warning_policy.action_threshold"),
      lower = 0, upper = 1
    )
  }

  parameters <- raw[["parameters", exact = TRUE]]
  .dpprior_schema_require(
    !is.null(parameters), "fit_diagnostic_parameters", "result.parameters",
    "finite canonical parameters for attached diagnostics", parameters
  )
  .dpprior_validate_parameters(parameters, "result.parameters")
  J <- raw[["J", exact = TRUE]]
  .dpprior_schema_require(
    .dpprior_schema_is_count(J, 1L), "fit_diagnostic_J", "result.J",
    "a positive integer J", J
  )

  M_selected <- authority[["M_selected", exact = TRUE]]
  M_required <- authority[["M_verification_required", exact = TRUE]]
  M_verification <- authority[["M_verification_used", exact = TRUE]]
  expected_required <- max(2L * M_selected, M_selected + 40L)
  .dpprior_schema_require(
    identical(M_required, as.integer(expected_required)) &&
      M_verification >= M_required && M_verification > M_selected,
    "fit_diagnostic_order_policy", authority_path,
    paste(
      "M_verification_required=max(2*M_selected,M_selected+40) and an",
      "independently higher used order"
    ),
    authority[c(
      "M_selected", "M_verification_required", "M_verification_used"
    )]
  )
  mode <- raw[["mode", exact = TRUE]]
  central_orders <- raw[["computation", exact = TRUE]][[
    "orders", exact = TRUE
  ]]
  if (identical(mode, "a1_proxy")) {
    .dpprior_schema_require(
      all(vapply(
        central_orders[c(
          "M_requested", "M_selected", "M_verification_required",
          "M_verification_used"
        )],
        is.null, logical(1)
      )) && identical(M_verification, M_required),
      "fit_diagnostic_A1_orders", authority_path,
      paste(
        "A1 calibration orders unavailable and a diagnostic verifier exactly",
        "at the required higher order"
      ),
      list(central = central_orders, diagnostic = authority)
    )
  } else {
    .dpprior_schema_require(
      identical(M_selected, central_orders[["M_selected", exact = TRUE]]) &&
        identical(
          M_required,
          central_orders[["M_verification_required", exact = TRUE]]
        ) && identical(
          M_verification,
          central_orders[["M_verification_used", exact = TRUE]]
        ),
      "fit_diagnostic_central_orders", authority_path,
      "orders byte-identical to the fit's central selected/required/used orders",
      list(central = central_orders, diagnostic = authority)
    )
  }

  a <- parameters[["a", exact = TRUE]]
  b <- parameters[["b", exact = TRUE]]
  K_evidence <- tryCatch(
    .get_K_pmf_support(
      J, a, b, M = M_selected, M_verify = M_verification,
      abs_tol = authority[["absolute_tolerance", exact = TRUE]],
      rel_tol = authority[["relative_tolerance", exact = TRUE]]
    ),
    error = function(error) error
  )
  .dpprior_schema_require(
    is.list(K_evidence) && !inherits(K_evidence, "condition") &&
      is.numeric(K_evidence[["pmf", exact = TRUE]]) &&
      is.numeric(K_evidence[["verification_pmf", exact = TRUE]]),
    "fit_diagnostic_recomputation", path,
    "successful selected/verifier diagnostic PMF recomputation", K_evidence
  )
  selected_pmf <- unname(as.numeric(K_evidence[["pmf", exact = TRUE]]))
  verifier_pmf <- unname(as.numeric(
    K_evidence[["verification_pmf", exact = TRUE]]
  ))
  .dpprior_schema_require(
    length(selected_pmf) == J && length(verifier_pmf) == J &&
      all(is.finite(selected_pmf)) && all(is.finite(verifier_pmf)) &&
      all(selected_pmf >= 0) && all(verifier_pmf >= 0),
    "fit_diagnostic_recomputation", paste0(path, ".K.pmf"),
    "finite nonnegative selected/verifier PMFs on support 1:J",
    list(selected = selected_pmf, verifier = verifier_pmf)
  )
  selected_K <- .dpprior_target_pmf_moments(selected_pmf)
  verifier_K <- .dpprior_target_pmf_moments(verifier_pmf)
  selected_weight <- as.numeric(mean_w1(a, b, M_selected))
  verifier_weight <- as.numeric(mean_w1(a, b, M_verification))
  selected_rho <- c(
    mean = as.numeric(mean_rho(a, b, M_selected)),
    variance = as.numeric(var_rho(a, b, M_selected))
  )
  verifier_rho <- c(
    mean = as.numeric(mean_rho(a, b, M_verification)),
    variance = as.numeric(var_rho(a, b, M_verification))
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
  absolute_tolerance <- authority[["absolute_tolerance", exact = TRUE]]
  relative_tolerance <- authority[["relative_tolerance", exact = TRUE]]
  scalar_tolerance <- function(selected, verifier) {
    absolute_tolerance + relative_tolerance *
      max(abs(selected), abs(verifier), 1)
  }
  refinement_tolerance <- c(
    K.mean = scalar_tolerance(selected_K[["mean"]], verifier_K[["mean"]]),
    K.variance = scalar_tolerance(
      selected_K[["variance"]], verifier_K[["variance"]]
    ),
    K.pmf_l1 = absolute_tolerance + relative_tolerance,
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
              refinement_tolerance[c("K.mean", "K.variance", "K.pmf_l1")]) &&
      abs(sum(selected_pmf) - 1) <=
        authority[["pmf_mass_tolerance", exact = TRUE]] &&
      abs(sum(verifier_pmf) - 1) <=
        authority[["pmf_mass_tolerance", exact = TRUE]],
    weights = delta[["weights.mean"]] <=
      refinement_tolerance[["weights.mean"]],
    coclustering = all(delta[c(
      "coclustering.mean", "coclustering.variance"
    )] <= refinement_tolerance[c(
      "coclustering.mean", "coclustering.variance"
    )])
  )
  expected_status <- ifelse(component_pass, "converged", "approximate")
  expected_usable <- component_pass |
    (!component_pass & authority[["allow_approximate", exact = TRUE]])
  expected_components <- list(
    alpha = list(
      status = unname(expected_status[["alpha"]]),
      usable = unname(expected_usable[["alpha"]]),
      verified = unname(component_pass[["alpha"]]),
      mean = a / b, CV = 1 / sqrt(a)
    ),
    K = list(
      status = unname(expected_status[["K"]]),
      usable = unname(expected_usable[["K"]]),
      verified = unname(component_pass[["K"]]),
      mean = unname(selected_K[["mean"]]),
      variance = unname(selected_K[["variance"]]),
      pmf = selected_pmf, M = M_selected
    ),
    weights = list(
      status = unname(expected_status[["weights"]]),
      usable = unname(expected_usable[["weights"]]),
      verified = unname(component_pass[["weights"]]),
      mean = selected_weight
    ),
    coclustering = list(
      status = unname(expected_status[["coclustering"]]),
      usable = unname(expected_usable[["coclustering"]]),
      verified = unname(component_pass[["coclustering"]]),
      mean = unname(selected_rho[["mean"]]),
      variance = unname(selected_rho[["variance"]])
    )
  )
  .dpprior_schema_require(
    identical(x[.DPPRIOR_DIAGNOSTIC_COMPONENTS], expected_components),
    "fit_diagnostic_component_truth", path,
    paste(
      "component values and status flags freshly recomputed from canonical",
      "J, parameters, orders, and fixed tolerances"
    ),
    list(
      recorded = x[.DPPRIOR_DIAGNOSTIC_COMPONENTS],
      expected = expected_components
    )
  )

  expected_policy_results <- if (is.null(warning_policy)) {
    list()
  } else if (identical(warning_policy[["estimand", exact = TRUE]], "W_SB")) {
    value <- as.numeric(.diagnostic_wsb_tail(
      warning_policy[["weight_threshold", exact = TRUE]], a, b
    ))
    outcome <- if (identical(
      warning_policy[["direction", exact = TRUE]], "above"
    )) {
      if (value > warning_policy[["action_threshold", exact = TRUE]]) {
        "triggered"
      } else "not_triggered"
    } else if (value < warning_policy[["action_threshold", exact = TRUE]]) {
      "triggered"
    } else "not_triggered"
    list(list(
      estimand = "W_SB",
      direction = warning_policy[["direction", exact = TRUE]],
      threshold = warning_policy[["action_threshold", exact = TRUE]],
      value = value, lower = NULL, upper = NULL,
      outcome = outcome, basis = "exact_tail_probability"
    ))
  } else {
    list(list(
      estimand = "W_max",
      direction = warning_policy[["direction", exact = TRUE]],
      threshold = warning_policy[["action_threshold", exact = TRUE]],
      value = NULL, lower = NULL, upper = NULL,
      outcome = "indeterminate", basis = "backend_unavailable"
    ))
  }
  .dpprior_schema_require(
    identical(x[["policy_results", exact = TRUE]], expected_policy_results),
    "fit_diagnostic_policy_truth", paste0(path, ".policy_results"),
    paste(
      "the normalized policy bound to a fresh W_SB value, or an explicit",
      "backend-unavailable W_max result"
    ),
    list(
      authority = warning_policy,
      recorded = x[["policy_results", exact = TRUE]],
      expected = expected_policy_results
    )
  )
  expected_warning_count <- sum(vapply(
    expected_policy_results,
    function(record) identical(record[["outcome", exact = TRUE]], "triggered"),
    logical(1)
  ))
  .dpprior_schema_require(
    length(x[["warnings", exact = TRUE]]) == expected_warning_count,
    "fit_diagnostic_warning_truth", paste0(path, ".warnings"),
    "exactly one warning for a freshly recomputed triggered policy",
    x[["warnings", exact = TRUE]]
  )
  invisible(TRUE)
}


.dpprior_sensitivity_close_numeric <- function(x, y, tolerance = NULL) {
  if (!.dpprior_schema_is_finite_scalar(x) ||
      !.dpprior_schema_is_finite_scalar(y)) return(FALSE)
  if (is.null(tolerance)) {
    tolerance <- 64 * .Machine$double.eps * max(1, abs(x), abs(y))
  }
  .dpprior_schema_is_finite_scalar(tolerance) && tolerance >= 0 &&
    abs(x - y) <= tolerance
}


.dpprior_sensitivity_validate_condition_summary <- function(
    x, path, nullable = TRUE) {
  if (is.null(x)) {
    .dpprior_schema_require(
      nullable, "sensitivity_condition", path,
      "a typed normalized condition summary", x
    )
    return(invisible(TRUE))
  }
  .dpprior_schema_exact_names(
    x, c("class", "classes", "code", "message"), path
  )
  for (field in c("class", "code", "message")) {
    .dpprior_schema_validate_scalar_character(
      x[[field, exact = TRUE]], paste0(path, ".", field),
      allow_empty = identical(field, "message")
    )
  }
  .dpprior_schema_validate_character_vector(
    x[["classes", exact = TRUE]], paste0(path, ".classes")
  )
  .dpprior_schema_require(
    length(x[["classes", exact = TRUE]]) > 0L &&
      !anyDuplicated(x[["classes", exact = TRUE]]) &&
      identical(x[["classes", exact = TRUE]][[1L]],
                x[["class", exact = TRUE]]) &&
      x[["class", exact = TRUE]] %in% x[["classes", exact = TRUE]],
    "sensitivity_condition_classes", paste0(path, ".classes"),
    "unique ordered classes led by and containing the primary class",
    x[["classes", exact = TRUE]]
  )
  invisible(TRUE)
}


.dpprior_sensitivity_validate_condition_evidence <- function(x, path) {
  .dpprior_schema_exact_names(
    x,
    c(
      "calibration", "calibration_warnings", "diagnostics",
      "diagnostic_warnings", "target", "interval"
    ),
    path
  )
  for (slot in c("calibration", "diagnostics", "target", "interval")) {
    .dpprior_sensitivity_validate_condition_summary(
      x[[slot, exact = TRUE]], paste0(path, ".", slot)
    )
  }
  for (slot in c("calibration_warnings", "diagnostic_warnings")) {
    warnings <- x[[slot, exact = TRUE]]
    .dpprior_schema_require(
      typeof(warnings) == "list" && is.list(warnings) &&
        !is.object(warnings) && is.null(dim(warnings)) &&
        .dpprior_schema_has_only_attributes(warnings),
      "sensitivity_warning_records", paste0(path, ".", slot),
      "an ordered ordinary list of typed normalized condition summaries",
      warnings
    )
    for (index in seq_along(warnings)) {
      .dpprior_sensitivity_validate_condition_summary(
        warnings[[index]], sprintf("%s.%s[[%d]]", path, slot, index),
        nullable = FALSE
      )
      .dpprior_schema_require(
        "warning" %in% warnings[[index]][["classes", exact = TRUE]],
        "sensitivity_warning_class",
        sprintf("%s.%s[[%d]].classes", path, slot, index),
        "a typed producer warning class chain containing warning",
        warnings[[index]][["classes", exact = TRUE]]
      )
    }
  }
  invisible(TRUE)
}


.dpprior_sensitivity_condition_contract <- function(x) {
  if (is.null(x)) return(NULL)
  contracts <- list(
    calibration_unusable = list(
      primary = "dpprior_calibration_unusable",
      classes = c(
        "dpprior_calibration_unusable", "dpprior_calibration_error",
        "dpprior_error", "error", "dpprior_condition", "condition"
      )
    ),
    fit_diagnostics_approximate = list(
      primary = "dpprior_diagnostics_approximation_error",
      classes = c(
        "dpprior_diagnostics_approximation_error", "dpprior_fit_error",
        "dpprior_calibration_error", "dpprior_error", "error",
        "dpprior_condition", "condition"
      )
    ),
    calibration_nonidentifiable_j1 = list(
      primary = "dpprior_calibration_nonidentifiable",
      classes = c(
        "dpprior_calibration_nonidentifiable", "dpprior_fit_error",
        "dpprior_calibration_error", "dpprior_error", "error",
        "dpprior_condition", "condition"
      )
    ),
    empty_interval_group = list(
      primary = "dpprior_interval_empty_tail_error",
      classes = c(
        "dpprior_interval_empty_tail_error", "dpprior_interval_infeasible",
        "dpprior_calibration_error", "dpprior_error", "error",
        "dpprior_condition", "condition"
      )
    ),
    mean_outside_group_mass_hull = list(
      primary = "dpprior_interval_mean_infeasible",
      classes = c(
        "dpprior_interval_mean_infeasible", "dpprior_interval_infeasible",
        "dpprior_calibration_error", "dpprior_error", "error",
        "dpprior_condition", "condition"
      )
    ),
    target_certified_infeasible = list(
      primary = "dpprior_interval_infeasible",
      classes = c(
        "dpprior_interval_infeasible", "dpprior_calibration_error",
        "dpprior_error", "error", "dpprior_condition", "condition"
      )
    ),
    maxent_root_not_bracketed = list(
      primary = "dpprior_interval_solver_failed",
      classes = c(
        "dpprior_interval_solver_failed", "dpprior_numerical_error",
        "dpprior_calibration_error", "dpprior_error", "error",
        "dpprior_condition", "condition"
      )
    ),
    maxent_root_solver_error = list(
      primary = "dpprior_interval_solver_failed",
      classes = c(
        "dpprior_interval_solver_failed", "dpprior_numerical_error",
        "dpprior_calibration_error", "dpprior_error", "error",
        "dpprior_condition", "condition"
      )
    ),
    target_failed = list(
      primary = "dpprior_interval_solver_failed",
      classes = c(
        "dpprior_interval_solver_failed", "dpprior_numerical_error",
        "dpprior_calibration_error", "dpprior_error", "error",
        "dpprior_condition", "condition"
      )
    ),
    sensitivity_backend_contract = list(
      primary = "dpprior_sensitivity_backend_contract_error",
      classes = c(
        "dpprior_sensitivity_backend_contract_error",
        "dpprior_backend_contract_error", "dpprior_calibration_error",
        "dpprior_sensitivity_error", "dpprior_error", "error",
        "dpprior_condition", "condition"
      )
    ),
    sensitivity_diagnostic_contract = list(
      primary = "dpprior_sensitivity_diagnostic_contract_error",
      classes = c(
        "dpprior_sensitivity_diagnostic_contract_error",
        "dpprior_diagnostics_error", "dpprior_sensitivity_error",
        "dpprior_error", "error", "dpprior_condition", "condition"
      )
    )
  )
  code <- x[["code", exact = TRUE]]
  contract <- contracts[[code, exact = TRUE]]
  if (is.null(contract)) return(NULL)
  if (!identical(x[["class", exact = TRUE]], contract$primary) ||
      !identical(x[["classes", exact = TRUE]], contract$classes)) {
    return(NULL)
  }
  code
}


.dpprior_sensitivity_validate_input_provenance <- function(
    x, request, method, path) {
  .dpprior_schema_exact_names(
    x,
    c(
      "method_explicit", "confidence_explicit", "requested_method",
      "selected_method", "is_fallback", "target_route"
    ),
    path
  )
  for (field in c("method_explicit", "confidence_explicit", "is_fallback")) {
    .dpprior_schema_validate_scalar_logical(
      x[[field, exact = TRUE]], paste0(path, ".", field)
    )
  }
  for (field in c("requested_method", "selected_method", "target_route")) {
    .dpprior_schema_validate_scalar_character(
      x[[field, exact = TRUE]], paste0(path, ".", field)
    )
  }
  pair <- x[c("requested_method", "selected_method", "is_fallback")]
  approved_pair <- identical(pair, list(
    requested_method = "A1", selected_method = "A1", is_fallback = FALSE
  )) || identical(pair, list(
    requested_method = "A2-MN", selected_method = "A2-MN",
    is_fallback = FALSE
  )) || identical(pair, list(
    requested_method = "A2-MN", selected_method = "A2-MN+NM",
    is_fallback = TRUE
  )) || identical(pair, list(
    requested_method = "A2-KL", selected_method = "A2-KL",
    is_fallback = FALSE
  ))
  route <- x[["target_route", exact = TRUE]]
  requested_method <- x[["requested_method", exact = TRUE]]
  route_method_ok <- if (route %in% c("interval", "strict_pmf")) {
    identical(requested_method, "A2-KL")
  } else if (identical(route, "coefficient_of_variation")) {
    requested_method %in% c("A1", "A2-MN")
  } else {
    requested_method %in% c("A1", "A2-MN", "A2-KL")
  }
  default_method <- if (route %in% c("interval", "strict_pmf")) {
    "A2-KL"
  } else {
    "A2-MN"
  }
  .dpprior_schema_require(
    approved_pair &&
      identical(request[["method", exact = TRUE]], requested_method) &&
      identical(method, x[["selected_method", exact = TRUE]]) &&
      route %in% c(
        "direct_variance", "qualitative_confidence",
        "coefficient_of_variation", "interval", "strict_pmf"
      ) && route_method_ok &&
      (x[["method_explicit", exact = TRUE]] ||
         identical(requested_method, default_method)) &&
      if (identical(route, "qualitative_confidence")) {
        request[["confidence", exact = TRUE]] %in% c("low", "medium", "high") &&
          (x[["confidence_explicit", exact = TRUE]] ||
             identical(request[["confidence", exact = TRUE]], "medium"))
      } else {
        !x[["confidence_explicit", exact = TRUE]] &&
          !("confidence" %in% names(request))
      },
    "sensitivity_input_provenance", path,
    paste(
      "a closed requested/selected/fallback method tuple, exact route,",
      "retained explicitness provenance, and request-method identity"
    ),
    x
  )
  invisible(TRUE)
}


.dpprior_sensitivity_validate_interval_authority <- function(x, J, path) {
  .dpprior_schema_exact_names(
    x,
    c(
      "lower", "upper", "type", "coverage", "family", "mu_K",
      "support", "endpoints"
    ),
    path
  )
  .dpprior_schema_require(
    .dpprior_schema_is_count(x[["lower", exact = TRUE]], 1L) &&
      .dpprior_schema_is_count(
        x[["upper", exact = TRUE]], x[["lower", exact = TRUE]]
      ) && x[["upper", exact = TRUE]] <= J &&
      .dpprior_schema_is_finite_scalar(x[["coverage", exact = TRUE]]) &&
      x[["coverage", exact = TRUE]] > 0 &&
      x[["coverage", exact = TRUE]] <= 1 &&
      x[["type", exact = TRUE]] %in%
        c("hard_bounds", "central_mass", "equal_tail") &&
      identical(x[["family", exact = TRUE]], "maxent") &&
      identical(
        x[["support", exact = TRUE]],
        c(lower = 1L, upper = as.integer(J))
      ) && identical(x[["endpoints", exact = TRUE]], "inclusive") &&
      (is.null(x[["mu_K", exact = TRUE]]) ||
         (.dpprior_schema_is_finite_scalar(x[["mu_K", exact = TRUE]]) &&
            x[["mu_K", exact = TRUE]] >= 1 &&
            x[["mu_K", exact = TRUE]] <= J)) &&
      if (identical(x[["type", exact = TRUE]], "hard_bounds")) {
        identical(x[["coverage", exact = TRUE]], 1) &&
          if (identical(x[["lower", exact = TRUE]],
                        x[["upper", exact = TRUE]])) {
            identical(
              x[["mu_K", exact = TRUE]],
              as.numeric(x[["lower", exact = TRUE]])
            )
          } else TRUE
      } else if (identical(x[["type", exact = TRUE]], "central_mass")) {
        !is.null(x[["mu_K", exact = TRUE]]) &&
          x[["coverage", exact = TRUE]] < 1
      } else {
        x[["coverage", exact = TRUE]] < 1
      },
    "sensitivity_interval_authority", path,
    paste(
      "the exact normalized inclusive interval authority, fixed support,",
      "closed type/family, and required/inferred mean semantics"
    ),
    x
  )
  invisible(TRUE)
}


.dpprior_sensitivity_validate_plain_pmf <- function(x, length, path) {
  .dpprior_schema_require(
    is.numeric(x) && !is.object(x) && is.null(dim(x)) && is.null(names(x)) &&
      length(x) == length && !anyNA(x) && all(is.finite(x)) &&
      all(x >= 0) && abs(sum(x) - 1) <= .TOL_PMF_SUM,
    "sensitivity_target_pmf", path,
    sprintf(
      "an unnamed finite nonnegative PMF of length %d with unit mass", length
    ),
    x
  )
  invisible(TRUE)
}


.dpprior_sensitivity_A2_KL_target_pmf <- function(target, J, path) {
  used <- target[["used", exact = TRUE]]
  pmf <- used[["pmf", exact = TRUE]]
  if (!is.null(pmf)) {
    pmf <- unname(as.numeric(pmf))
    .dpprior_sensitivity_validate_plain_pmf(pmf, J, path)
    return(pmf)
  }
  mean_name <- if ("mu_K" %in% names(used)) "mu_K" else "mean"
  variance_name <- if ("var_K" %in% names(used)) "var_K" else "variance"
  .dpprior_schema_require(
    identical(target[["kind", exact = TRUE]], "moments") &&
      .dpprior_schema_is_finite_scalar(used[[mean_name, exact = TRUE]]) &&
      .dpprior_schema_is_finite_scalar(used[[variance_name, exact = TRUE]]) &&
      used[[mean_name, exact = TRUE]] > 0 &&
      used[[variance_name, exact = TRUE]] > 0,
    "sensitivity_A2_KL_target_authority", path,
    paste(
      "a canonical target PMF or positive finite moment authority for the",
      "fixed continuity-corrected scaled-chi-square objective"
    ),
    target
  )
  mean <- used[[mean_name, exact = TRUE]]
  variance <- used[[variance_name, exact = TRUE]]
  reconstructed <- tryCatch(
    .a2_kl_verify_chisq_pmf(
      J, df = 2 * mean^2 / variance, scale = variance / (2 * mean)
    ),
    error = function(error) error
  )
  .dpprior_schema_require(
    is.numeric(reconstructed) && !is.object(reconstructed) &&
      is.null(dim(reconstructed)) && length(reconstructed) == J &&
      !anyNA(reconstructed) && all(is.finite(reconstructed)) &&
      all(reconstructed >= 0) && abs(sum(reconstructed) - 1) <= .TOL_PMF_SUM,
    "sensitivity_A2_KL_target_reconstruction", path,
    "successful independent scaled-chi-square objective-PMF reconstruction",
    reconstructed
  )
  unname(as.numeric(reconstructed))
}


.dpprior_sensitivity_canonical_target <- function(request, route, path) {
  .dpprior_schema_require(
    typeof(request) == "list" && is.list(request) && !is.object(request) &&
      is.null(dim(request)) && !is.null(names(request)) &&
      !anyNA(names(request)) && all(nzchar(names(request))) &&
      !anyDuplicated(names(request)),
    "sensitivity_request_authority", path,
    "one ordinary uniquely named normalized scientific request", request
  )
  required_names <- switch(
    route,
    direct_variance = c("J", "mu_K", "var_K", "method", "M"),
    qualitative_confidence = c("J", "mu_K", "confidence", "method", "M"),
    coefficient_of_variation = c("J", "mu_K", "cv_K", "method", "M"),
    interval = c("J", "K_interval", "mu_K", "method", "M"),
    strict_pmf = c("J", "target_pmf", "method", "M"),
    character()
  )
  .dpprior_schema_require(
    identical(names(request), required_names),
    "sensitivity_request_fields", path,
    paste(required_names, collapse = ", "), names(request)
  )
  J <- request[["J", exact = TRUE]]
  .dpprior_schema_require(
    .dpprior_schema_is_count(J, 1L) &&
      .dpprior_schema_is_count(request[["M", exact = TRUE]], 10L) &&
      request[["M", exact = TRUE]] <= 256L,
    "sensitivity_request_authority", path,
    "positive integer J and selected quadrature order M in [10, 256]", request
  )
  if (identical(route, "direct_variance")) {
    .dpprior_schema_require(
      .dpprior_schema_is_finite_scalar(request[["mu_K", exact = TRUE]]) &&
        .dpprior_schema_is_finite_scalar(request[["var_K", exact = TRUE]]),
      "sensitivity_direct_target", path,
      "finite direct mean and variance", request
    )
  } else if (identical(route, "qualitative_confidence")) {
    .dpprior_schema_require(
      .dpprior_schema_is_finite_scalar(request[["mu_K", exact = TRUE]]) &&
        request[["confidence", exact = TRUE]] %in% c("low", "medium", "high"),
      "sensitivity_confidence_target", path,
      "finite mean and low/medium/high confidence", request
    )
  } else if (identical(route, "coefficient_of_variation")) {
    .dpprior_schema_require(
      .dpprior_schema_is_finite_scalar(request[["mu_K", exact = TRUE]]) &&
        .dpprior_schema_is_finite_scalar(request[["cv_K", exact = TRUE]]) &&
        request[["cv_K", exact = TRUE]] > 0,
      "sensitivity_cv_target", path,
      "finite mean and positive finite coefficient of variation", request
    )
  } else if (identical(route, "interval")) {
    .dpprior_sensitivity_validate_interval_authority(
      request[["K_interval", exact = TRUE]], J,
      paste0(path, ".K_interval")
    )
    .dpprior_schema_require(
      identical(
        request[["mu_K", exact = TRUE]],
        request[["K_interval", exact = TRUE]][["mu_K", exact = TRUE]]
      ),
      "sensitivity_interval_mean_authority", paste0(path, ".mu_K"),
      "top-level mu_K exactly identical to normalized K_interval$mu_K",
      request[["mu_K", exact = TRUE]]
    )
  } else if (identical(route, "strict_pmf")) {
    pmf <- request[["target_pmf", exact = TRUE]]
    .dpprior_schema_require(
      is.numeric(pmf) && !is.object(pmf) && is.null(dim(pmf)) &&
        is.null(names(pmf)) && length(pmf) %in% c(J, J + 1L) &&
        !anyNA(pmf) && all(is.finite(pmf)) && all(pmf >= 0) &&
        if (length(pmf) == J + 1L) identical(as.numeric(pmf[[1L]]), 0) else TRUE,
      "sensitivity_strict_pmf_request", paste0(path, ".target_pmf"),
      paste(
        "an unnamed finite nonnegative J or J+1 PMF, with an exact",
        "structural K=0 zero for J+1"
      ),
      pmf
    )
  }
  constructor_args <- switch(
    route,
    direct_variance = list(
      J = J, mu_K = request[["mu_K", exact = TRUE]],
      var_K = request[["var_K", exact = TRUE]]
    ),
    qualitative_confidence = list(
      J = J, mu_K = request[["mu_K", exact = TRUE]],
      confidence = request[["confidence", exact = TRUE]]
    ),
    coefficient_of_variation = list(
      J = J, mu_K = request[["mu_K", exact = TRUE]],
      cv_K = request[["cv_K", exact = TRUE]]
    ),
    interval = list(
      J = J, mu_K = request[["mu_K", exact = TRUE]],
      K_interval = request[["K_interval", exact = TRUE]][if (identical(
        request[["K_interval", exact = TRUE]][["type", exact = TRUE]],
        "hard_bounds"
      )) {
        c("lower", "upper", "type", "family")
      } else {
        c("lower", "upper", "type", "coverage", "family")
      }]
    ),
    strict_pmf = list(
      J = J, target_pmf = request[["target_pmf", exact = TRUE]]
    )
  )
  canonical <- tryCatch(
    list(target = do.call(.dp_target_K, constructor_args), condition = NULL),
    error = function(condition) list(
      target = condition[["result", exact = TRUE]], condition = condition
    )
  )
  validation <- if (inherits(canonical[["target", exact = TRUE]],
                             "dpprior_K_target")) {
    tryCatch(
      .dpprior_validate_target_v1(
        canonical[["target", exact = TRUE]], collect = TRUE
      ),
      error = function(condition) NULL
    )
  } else {
    NULL
  }
  .dpprior_schema_require(
    is.list(validation) && isTRUE(validation[["valid", exact = TRUE]]),
    "sensitivity_target_canonicalization", path,
    "successful producer-canonical target reconstruction", canonical
  )
  canonical[["target", exact = TRUE]]
}


.dpprior_sensitivity_validate_fit_evidence <- function(
    x, key, path = "sensitivity.fit_evidence") {
  .dpprior_schema_exact_names(
    x,
    c(
      "request", "input_provenance", "weight_target", "target", "method",
      "status", "usable", "verified", "parameters", "evaluator",
      "selected_snapshot", "verifier_snapshot", "condition_evidence", "source"
    ),
    path
  )
  .dpprior_schema_validate_named_list(
    x[["request", exact = TRUE]], paste0(path, ".request"),
    allow_empty = FALSE
  )
  .dpprior_schema_validate_plain_record_value(
    x[["request", exact = TRUE]], paste0(path, ".request")
  )
  .dpprior_schema_require(
    is.null(x[["weight_target", exact = TRUE]]),
    "sensitivity_weight_target_authority", paste0(path, ".weight_target"),
    paste(
      "NULL until a typed Phase-9 weight-target sensitivity authority is",
      "retained; opaque producer metadata cannot authorize scientific claims"
    ),
    x[["weight_target", exact = TRUE]]
  )
  .dpprior_schema_validate_scalar_character(
    x[["method", exact = TRUE]], paste0(path, ".method")
  )
  .dpprior_schema_require(
    x[["method", exact = TRUE]] %in%
      unique(unlist(.DPPRIOR_MODE_METHODS[c(
        "a1_proxy", "a2_moment", "a2_kl"
      )], use.names = FALSE)),
    "sensitivity_fit_method", paste0(path, ".method"),
    "a closed elicitation-fit method", x[["method", exact = TRUE]]
  )
  .dpprior_sensitivity_validate_input_provenance(
    x[["input_provenance", exact = TRUE]],
    x[["request", exact = TRUE]], x[["method", exact = TRUE]],
    paste0(path, ".input_provenance")
  )
  .dpprior_validate_status_record(
    list(
      status = x[["status", exact = TRUE]],
      usable = x[["usable", exact = TRUE]],
      verified = x[["verified", exact = TRUE]], message = ""
    ),
    paste0(path, ".status_record")
  )
  .dpprior_schema_validate_scalar_character(
    x[["source", exact = TRUE]], paste0(path, ".source")
  )
  .dpprior_schema_require(
    identical(
      x[["source", exact = TRUE]], "retained_canonical_fit_evidence"
    ),
    "sensitivity_fit_source", paste0(path, ".source"),
    "retained_canonical_fit_evidence", x[["source", exact = TRUE]]
  )
  .dpprior_sensitivity_validate_condition_evidence(
    x[["condition_evidence", exact = TRUE]],
    paste0(path, ".condition_evidence")
  )

  evaluator <- x[["evaluator", exact = TRUE]]
  .dpprior_schema_exact_names(
    evaluator,
    c(
      "method", "M_selected", "M_verification", "absolute_tolerance",
      "relative_tolerance", "W_max_point_policy"
    ),
    paste0(path, ".evaluator")
  )
  .dpprior_schema_validate_scalar_character(
    evaluator[["method", exact = TRUE]], paste0(path, ".evaluator.method")
  )
  .dpprior_schema_require(
    identical(
      evaluator[["method", exact = TRUE]],
      "gauss_laguerre_marginal_pmf_and_moments"
    ) &&
      .dpprior_schema_is_count(evaluator[["M_selected", exact = TRUE]], 10L) &&
      .dpprior_schema_is_count(
        evaluator[["M_verification", exact = TRUE]],
        max(
          2L * evaluator[["M_selected", exact = TRUE]],
          evaluator[["M_selected", exact = TRUE]] + 40L
        )
      ) && identical(
        evaluator[["M_verification", exact = TRUE]],
        max(
          2L * evaluator[["M_selected", exact = TRUE]],
          evaluator[["M_selected", exact = TRUE]] + 40L
        )
      ) && identical(
        evaluator[["absolute_tolerance", exact = TRUE]], 1e-10
      ) && identical(
        evaluator[["relative_tolerance", exact = TRUE]], 1e-8
      ) && identical(
        evaluator[["W_max_point_policy", exact = TRUE]],
        "unavailable_without_retained_typed_backend_evidence"
      ),
    "sensitivity_evaluator_authority", paste0(path, ".evaluator"),
    paste(
      "the retained selected order, required independent order, fixed",
      "numerical tolerances, and fail-closed W_max point policy"
    ),
    evaluator
  )

  request <- x[["request", exact = TRUE]]
  input_provenance <- x[["input_provenance", exact = TRUE]]
  .dpprior_schema_require(
    .dpprior_schema_is_count(request[["J", exact = TRUE]], 1L) &&
      identical(
        request[["M", exact = TRUE]], evaluator[["M_selected", exact = TRUE]]
      ),
    "sensitivity_request_authority", paste0(path, ".request"),
    "request J/M identical to the retained evaluator authority",
    request
  )
  J <- request[["J", exact = TRUE]]
  canonical_target <- .dpprior_sensitivity_canonical_target(
    request, input_provenance[["target_route", exact = TRUE]],
    paste0(path, ".request")
  )
  canonical_target_raw <- unclass(canonical_target)
  expected_target <- list(
    kind = canonical_target_raw[["kind", exact = TRUE]],
    J = canonical_target_raw[["J", exact = TRUE]],
    request = canonical_target_raw[["request", exact = TRUE]],
    used = canonical_target_raw[["used", exact = TRUE]]
  )
  target <- x[["target", exact = TRUE]]
  .dpprior_schema_exact_names(
    target, c("kind", "J", "request", "used"),
    paste0(path, ".target")
  )
  .dpprior_schema_require(
    identical(target, expected_target),
    "sensitivity_target_authority", paste0(path, ".target"),
    paste(
      "the exact producer-reconstructed target kind/J plus canonical",
      "request and used scientific authority"
    ),
    list(recorded = target, expected = expected_target)
  )
  if (identical(input_provenance[["target_route", exact = TRUE]],
                "strict_pmf")) {
    .dpprior_sensitivity_validate_plain_pmf(
      target[["used", exact = TRUE]][["pmf", exact = TRUE]], J,
      paste0(path, ".target.used.pmf")
    )
  }
  target_condition <- if (canonical_target_raw[["status", exact = TRUE]] %in%
                          c("infeasible", "failed")) {
    tryCatch(
      .dp_target_K_stop_unusable(canonical_target),
      error = function(condition) condition
    )
  } else {
    NULL
  }
  target_condition_summary <- if (inherits(target_condition, "condition")) {
    list(
      class = class(target_condition)[[1L]],
      classes = class(target_condition),
      code = target_condition[["code", exact = TRUE]],
      message = conditionMessage(target_condition)
    )
  } else {
    NULL
  }

  parameters <- x[["parameters", exact = TRUE]]
  nonidentifiable_J1 <- identical(J, 1L) &&
    canonical_target_raw[["status", exact = TRUE]] %in%
      c("converged", "boundary")
  if (nonidentifiable_J1) {
    .dpprior_schema_require(
      is.null(parameters) &&
        identical(x[["status", exact = TRUE]], "infeasible") &&
        !x[["usable", exact = TRUE]] && x[["verified", exact = TRUE]] &&
        is.null(x[["selected_snapshot", exact = TRUE]]) &&
        is.null(x[["verifier_snapshot", exact = TRUE]]),
      "sensitivity_J1_nonidentifiable", path,
      paste(
        "the producer-canonical nonidentifiable J=1 strict-PMF outcome:",
        "infeasible/FALSE/TRUE with no parameters or numerical snapshots"
      ),
      x
    )
  }
  if (is.null(parameters)) {
    target_status <- canonical_target_raw[["status", exact = TRUE]]
    target_outcome_ok <- if (identical(
      x[["status", exact = TRUE]], "infeasible"
    )) {
      nonidentifiable_J1 || identical(target_status, "infeasible")
    } else {
      !identical(target_status, "infeasible")
    }
    .dpprior_schema_require(
      x[["status", exact = TRUE]] %in% c("failed", "infeasible") &&
        !x[["usable", exact = TRUE]] &&
        target_outcome_ok &&
        is.null(x[["selected_snapshot", exact = TRUE]]) &&
        is.null(x[["verifier_snapshot", exact = TRUE]]),
      "sensitivity_fit_unavailable", path,
      paste(
        "failed/infeasible non-usable fit evidence bound to the canonical",
        "target outcome, with no fabricated numerical snapshots"
      ),
      x
    )
    return(list(
      key = key, J = J, request = request, target = target,
      input_provenance = input_provenance,
      target_status = canonical_target_raw[["status", exact = TRUE]],
      target_usable = canonical_target_raw[["usable", exact = TRUE]],
      target_verified = canonical_target_raw[["verified", exact = TRUE]],
      target_condition = target_condition_summary,
      method = x[["method", exact = TRUE]], status = x[["status", exact = TRUE]],
      usable = x[["usable", exact = TRUE]], verified = x[["verified", exact = TRUE]],
      parameters = NULL, evaluator = evaluator, selected = NULL,
      verifier = NULL
    ))
  }

  .dpprior_validate_parameters(parameters, paste0(path, ".parameters"))
  .dpprior_schema_require(
    canonical_target_raw[["status", exact = TRUE]] %in%
      c("converged", "boundary") &&
      canonical_target_raw[["usable", exact = TRUE]] &&
      canonical_target_raw[["verified", exact = TRUE]] &&
      !x[["status", exact = TRUE]] %in% c("failed", "infeasible") &&
      !is.null(x[["selected_snapshot", exact = TRUE]]) &&
      !is.null(x[["verifier_snapshot", exact = TRUE]]) &&
      if (identical(x[["method", exact = TRUE]], "A1")) {
        identical(x[["status", exact = TRUE]], "approximate") &&
          x[["usable", exact = TRUE]] && !x[["verified", exact = TRUE]]
      } else if (identical(x[["status", exact = TRUE]], "approximate")) {
        !x[["usable", exact = TRUE]] && !x[["verified", exact = TRUE]]
      } else {
        x[["status", exact = TRUE]] %in% c("converged", "boundary") &&
          x[["usable", exact = TRUE]] && x[["verified", exact = TRUE]]
      },
    "sensitivity_fit_available", path,
    paste(
      "a usable/verified canonical target, finite parameters, both retained",
      "evaluator snapshots, and the exact A1 versus A2 status quartet"
    ),
    x
  )
  numerical <- tryCatch({
    logS <- compute_log_stirling(J)
    pmf <- pmf_K_marginal(
      J, parameters[["a", exact = TRUE]], parameters[["b", exact = TRUE]],
      logS, M = evaluator[["M_selected", exact = TRUE]],
      M_verify = evaluator[["M_verification", exact = TRUE]],
      abs_tol = evaluator[["absolute_tolerance", exact = TRUE]],
      rel_tol = evaluator[["relative_tolerance", exact = TRUE]], strict = FALSE
    )
    verification_pmf <- attr(pmf, ".marginal_verification_pmf", exact = TRUE)
    list(
      selected_pmf = unname(as.numeric(pmf[-1L])),
      verifier_pmf = unname(as.numeric(verification_pmf[-1L]))
    )
  }, error = function(error) error)
  .dpprior_schema_require(
    is.list(numerical) &&
      identical(names(numerical), c("selected_pmf", "verifier_pmf")) &&
      length(numerical[["selected_pmf", exact = TRUE]]) == J &&
      length(numerical[["verifier_pmf", exact = TRUE]]) == J,
    "sensitivity_evaluator_failure", path,
    "successful selected and independent marginal-PMF recomputation", numerical
  )
  expected <- list(
    selected_snapshot = list(
      M = evaluator[["M_selected", exact = TRUE]],
      K = c(
        as.list(.dpprior_target_pmf_moments(numerical$selected_pmf)),
        list(pmf = numerical$selected_pmf)
      ),
      finite = TRUE, source = "selected_order"
    ),
    verifier_snapshot = list(
      M = evaluator[["M_verification", exact = TRUE]],
      K = c(
        as.list(.dpprior_target_pmf_moments(numerical$verifier_pmf)),
        list(pmf = numerical$verifier_pmf)
      ),
      finite = TRUE, source = "independent_verifier"
    )
  )
  for (snapshot_name in names(expected)) {
    recorded <- x[[snapshot_name, exact = TRUE]]
    expected_snapshot <- expected[[snapshot_name, exact = TRUE]]
    .dpprior_schema_exact_names(
      recorded, c("M", "K", "finite", "source"),
      paste0(path, ".", snapshot_name)
    )
    .dpprior_schema_exact_names(
      recorded[["K", exact = TRUE]], c("mean", "variance", "pmf"),
      paste0(path, ".", snapshot_name, ".K")
    )
    .dpprior_schema_require(
      identical(recorded[["M", exact = TRUE]], expected_snapshot$M) &&
        identical(recorded[["finite", exact = TRUE]], TRUE) &&
        identical(recorded[["source", exact = TRUE]], expected_snapshot$source) &&
        .dpprior_sensitivity_close_numeric(
          recorded[["K", exact = TRUE]][["mean", exact = TRUE]],
          expected_snapshot$K$mean
        ) && .dpprior_sensitivity_close_numeric(
          recorded[["K", exact = TRUE]][["variance", exact = TRUE]],
          expected_snapshot$K$variance
        ) && is.numeric(recorded[["K", exact = TRUE]][["pmf", exact = TRUE]]) &&
        !is.object(recorded[["K", exact = TRUE]][["pmf", exact = TRUE]]) &&
        is.null(dim(recorded[["K", exact = TRUE]][["pmf", exact = TRUE]])) &&
        length(recorded[["K", exact = TRUE]][["pmf", exact = TRUE]]) == J &&
        sum(abs(
          recorded[["K", exact = TRUE]][["pmf", exact = TRUE]] -
            expected_snapshot$K$pmf
        )) <= .TOL_PMF_SUM,
      "sensitivity_fit_snapshot_truth", paste0(path, ".", snapshot_name),
      "the freshly recomputed fixed-order K PMF/moments and exact source/M",
      recorded
    )
  }
  selected <- expected[["selected_snapshot", exact = TRUE]]
  verifier <- expected[["verifier_snapshot", exact = TRUE]]
  stability_tolerance <- evaluator[["absolute_tolerance", exact = TRUE]] +
    evaluator[["relative_tolerance", exact = TRUE]] * c(
      mean = max(abs(selected$K$mean), abs(verifier$K$mean), 1),
      variance = max(abs(selected$K$variance), abs(verifier$K$variance), 1)
    )
  stability_delta <- abs(c(
    mean = selected$K$mean - verifier$K$mean,
    variance = selected$K$variance - verifier$K$variance
  ))
  if (x[["verified", exact = TRUE]]) {
    .dpprior_schema_require(
      all(stability_delta <= stability_tolerance),
      "sensitivity_fit_order_stability", path,
      "selected/verifier K-moment agreement at fixed central tolerances",
      list(delta = stability_delta, tolerance = stability_tolerance)
    )
  }
  if (identical(x[["method", exact = TRUE]], "A2-KL") &&
      x[["verified", exact = TRUE]]) {
    target_pmf <- .dpprior_sensitivity_A2_KL_target_pmf(
      target, J, paste0(path, ".target.used")
    )
    target_moments <- .dpprior_target_pmf_moments(target_pmf)
    adequacy_tolerance <- c(
      kl = 0.015, l1 = 0.11, mean_scaled = 0.01,
      variance_scaled = 0.065
    )
    distribution_metrics <- function(snapshot) {
      pmf <- snapshot[["K", exact = TRUE]][["pmf", exact = TRUE]]
      positive <- target_pmf > 0
      if (any(pmf[positive] <= 0)) {
        return(setNames(rep(Inf, 4L), names(adequacy_tolerance)))
      }
      moments <- .dpprior_target_pmf_moments(pmf)
      c(
        kl = sum(target_pmf[positive] * log(
          target_pmf[positive] / pmf[positive]
        )),
        l1 = sum(abs(target_pmf - pmf)),
        mean_scaled = abs(
          moments[["mean"]] - target_moments[["mean"]]
        ) / max(1, sqrt(target_moments[["variance"]])),
        variance_scaled = abs(
          moments[["variance"]] - target_moments[["variance"]]
        ) / max(1, target_moments[["variance"]])
      )
    }
    selected_adequacy <- distribution_metrics(selected)
    verifier_adequacy <- distribution_metrics(verifier)
    pmf_order_l1 <- sum(abs(
      selected[["K", exact = TRUE]][["pmf", exact = TRUE]] -
        verifier[["K", exact = TRUE]][["pmf", exact = TRUE]]
    ))
    .dpprior_schema_require(
      all(is.finite(selected_adequacy)) &&
        all(is.finite(verifier_adequacy)) &&
        all(selected_adequacy <= adequacy_tolerance) &&
        all(verifier_adequacy <= adequacy_tolerance) &&
        pmf_order_l1 <= 1e-10 + 1e-8,
      "sensitivity_fit_pmf_target_truth", paste0(path, ".target"),
      paste(
        "the fixed four-gate A2-KL adequacy contract at selected and",
        "independent orders plus fixed PMF-order stability"
      ),
      list(
        selected = selected_adequacy, verifier = verifier_adequacy,
        tolerance = adequacy_tolerance, pmf_order_l1 = pmf_order_l1,
        pmf_order_tolerance = 1e-10 + 1e-8
      )
    )
  }
  if (target[["kind", exact = TRUE]] %in% c("moments", "cv") &&
      x[["verified", exact = TRUE]] &&
      x[["method", exact = TRUE]] %in% c("A2-MN", "A2-MN+NM")) {
    target_implied <- canonical_target_raw[["implied", exact = TRUE]]
    .dpprior_schema_require(
      !is.null(target_implied) &&
      .dpprior_sensitivity_close_numeric(
        selected$K$mean, target_implied[["mean", exact = TRUE]],
        evaluator$absolute_tolerance + evaluator$relative_tolerance *
          max(abs(target_implied[["mean", exact = TRUE]]), 1)
      ) && .dpprior_sensitivity_close_numeric(
        selected$K$variance, target_implied[["variance", exact = TRUE]],
        evaluator$absolute_tolerance + evaluator$relative_tolerance *
          max(abs(target_implied[["variance", exact = TRUE]]), 1)
      ),
      "sensitivity_fit_target_truth", paste0(path, ".target"),
      "verified selected K moments meeting the retained request target",
      list(selected = selected$K, target = target)
    )
  }
  list(
    key = key, J = J, request = request, target = target,
    input_provenance = input_provenance,
    target_status = canonical_target_raw[["status", exact = TRUE]],
    target_usable = canonical_target_raw[["usable", exact = TRUE]],
    target_verified = canonical_target_raw[["verified", exact = TRUE]],
    target_condition = target_condition_summary,
    method = x[["method", exact = TRUE]], status = x[["status", exact = TRUE]],
    usable = x[["usable", exact = TRUE]], verified = x[["verified", exact = TRUE]],
    parameters = parameters, evaluator = evaluator,
    selected = selected, verifier = verifier
  )
}


.dpprior_validate_sensitivity_extension <- function(x, path = "sensitivity") {
  .dpprior_schema_exact_names(
    x,
    c(
      "scenarios", "scenario_results", "fit_evidence", "conditions",
      "interval_checks", "metrics_long", "local", "global", "metadata"
    ),
    path
  )
  .dpprior_schema_validate_plain_data_frame(
    x[["scenarios", exact = TRUE]], paste0(path, ".scenarios")
  )
  .dpprior_schema_validate_plain_data_frame(
    x[["scenario_results", exact = TRUE]],
    paste0(path, ".scenario_results")
  )
  .dpprior_schema_require(
    identical(names(x[["scenarios", exact = TRUE]]),
              c(
                "scenario_key", "scenario_label", "canonical_content",
                "base_key", "diagnostics_requested", "method_explicit",
                "confidence_explicit", "effective_method",
                "effective_confidence"
              )) &&
      identical(names(x[["scenario_results", exact = TRUE]]),
                c(
                  "scenario_key", "status", "usable", "verified",
                  .DPPRIOR_SENSITIVITY_METRICS
                )),
    "sensitivity_table_fields", path,
    "the exact canonical scenario and result table columns",
    list(
      scenarios = names(x[["scenarios", exact = TRUE]]),
      results = names(x[["scenario_results", exact = TRUE]])
    )
  )
  metadata <- x[["metadata", exact = TRUE]]
  .dpprior_schema_exact_names(metadata, "J", paste0(path, ".metadata"))
  .dpprior_schema_require(
    .dpprior_schema_is_count(metadata[["J", exact = TRUE]], 1L),
    "sensitivity_J", paste0(path, ".metadata.J"),
    "a positive integer J", metadata[["J", exact = TRUE]]
  )
  sensitivity_J <- metadata[["J", exact = TRUE]]
  keys <- x[["scenarios", exact = TRUE]][["scenario_key", exact = TRUE]]
  scenarios <- x[["scenarios", exact = TRUE]]
  canonical_content <- scenarios[["canonical_content", exact = TRUE]]
  expected_base_keys <- if (is.character(canonical_content) &&
                            !anyNA(canonical_content)) {
    unname(vapply(
      canonical_content, .dp_sensitivity_content_key, character(1)
    ))
  } else {
    rep(NA_character_, length(keys))
  }
  expected_final_keys <- expected_base_keys
  colliding_base_keys <- unique(expected_base_keys[
    duplicated(expected_base_keys) |
      duplicated(expected_base_keys, fromLast = TRUE)
  ])
  for (base_key in colliding_base_keys) {
    members <- which(expected_base_keys == base_key)
    rank <- match(
      canonical_content[members],
      sort(canonical_content[members], method = "radix")
    )
    assigned <- character()
    for (member in members[order(rank)]) {
      suffix <- rank[[match(member, members)]]
      repeat {
        candidate <- paste0(base_key, "_c", sprintf("%03d", suffix))
        if (!(candidate %in% expected_base_keys) &&
            !(candidate %in% expected_final_keys[
              setdiff(seq_along(expected_final_keys), members)
            ]) && !(candidate %in% assigned)) {
          break
        }
        suffix <- suffix + 1L
      }
      expected_final_keys[[member]] <- candidate
      assigned <- c(assigned, candidate)
    }
  }
  .dpprior_schema_require(
    is.character(keys) && length(keys) > 0L && !anyNA(keys) &&
      all(nzchar(keys)) &&
      !anyDuplicated(keys) && identical(keys, sort(keys, method = "radix")) &&
      all(grepl("^scn_[0-9a-f]{16}(_c[0-9]{3})?$", keys)) &&
      identical(
        scenarios[["base_key", exact = TRUE]], expected_base_keys
      ) && identical(keys, expected_final_keys) &&
      identical(
        x[["scenario_results", exact = TRUE]][["scenario_key", exact = TRUE]],
        keys
      ),
    "sensitivity_keys", paste0(path, ".scenarios.scenario_key"),
    "unique nonempty lexicographically ordered keys shared with results", keys
  )
  .dpprior_schema_require(
    is.character(scenarios[["scenario_label", exact = TRUE]]) &&
      all(is.na(scenarios[["scenario_label", exact = TRUE]])) &&
      is.character(scenarios[["canonical_content", exact = TRUE]]) &&
      !anyNA(scenarios[["canonical_content", exact = TRUE]]) &&
      all(nzchar(scenarios[["canonical_content", exact = TRUE]])) &&
      !anyDuplicated(scenarios[["canonical_content", exact = TRUE]]) &&
      is.character(scenarios[["base_key", exact = TRUE]]) &&
      !anyNA(scenarios[["base_key", exact = TRUE]]) &&
      all(grepl(
        "^scn_[0-9a-f]{16}$", scenarios[["base_key", exact = TRUE]]
      )) && all(vapply(
        c("diagnostics_requested", "method_explicit", "confidence_explicit"),
        function(field) {
          is.logical(scenarios[[field, exact = TRUE]]) &&
            !anyNA(scenarios[[field, exact = TRUE]])
        }, logical(1)
      )) && is.character(scenarios[["effective_method", exact = TRUE]]) &&
      !anyNA(scenarios[["effective_method", exact = TRUE]]) &&
      all(nzchar(scenarios[["effective_method", exact = TRUE]])) &&
      all(scenarios[["effective_method", exact = TRUE]] %in%
            unique(unlist(.DPPRIOR_MODE_METHODS[c(
              "a1_proxy", "a2_moment", "a2_kl"
            )], use.names = FALSE))) &&
      is.character(scenarios[["effective_confidence", exact = TRUE]]) &&
      all(
        is.na(scenarios[["effective_confidence", exact = TRUE]]) |
          scenarios[["effective_confidence", exact = TRUE]] %in%
            c("low", "medium", "high")
      ) && all(
        !scenarios[["confidence_explicit", exact = TRUE]] |
          !is.na(scenarios[["effective_confidence", exact = TRUE]])
      ),
    "sensitivity_scenario_identity", paste0(path, ".scenarios"),
    paste(
      "unique canonical scenario content, frozen key syntax, typed request",
      "flags, closed effective method/confidence vocabulary, and labels",
      "quarantined outside scientific identity"
    ), scenarios
  )
  .dpprior_schema_validate_named_list(
    x[["fit_evidence", exact = TRUE]], paste0(path, ".fit_evidence")
  )
  without_warning_sequences <- function(record, nested = FALSE) {
    if (typeof(record) != "list" || !is.list(record) || is.object(record) ||
        is.null(names(record))) {
      return(record)
    }
    condition_record <- if (nested &&
                            "condition_evidence" %in% names(record)) {
      record[["condition_evidence", exact = TRUE]]
    } else {
      record
    }
    if (typeof(condition_record) == "list" && is.list(condition_record) &&
        !is.object(condition_record) && !is.null(names(condition_record))) {
      warning_slots <- intersect(
        c("calibration_warnings", "diagnostic_warnings"),
        names(condition_record)
      )
      condition_record[warning_slots] <- rep(list(NULL), length(warning_slots))
      if (nested) {
        record[["condition_evidence"]] <- condition_record
      } else {
        record <- condition_record
      }
    }
    record
  }
  serializable_fit_evidence <- lapply(
    x[["fit_evidence", exact = TRUE]], without_warning_sequences,
    nested = TRUE
  )
  .dpprior_schema_validate_plain_record_value(
    serializable_fit_evidence, paste0(path, ".fit_evidence")
  )
  .dpprior_schema_require(
    identical(names(x[["fit_evidence", exact = TRUE]]), keys),
    "sensitivity_fit_evidence_keys", paste0(path, ".fit_evidence"),
    "one ordered retained fit-evidence record for every scenario key",
    names(x[["fit_evidence", exact = TRUE]])
  )
  fit_truth <- lapply(seq_along(keys), function(index) {
    key <- keys[[index]]
    .dpprior_sensitivity_validate_fit_evidence(
      x[["fit_evidence", exact = TRUE]][[key, exact = TRUE]], key,
      paste0(path, ".fit_evidence.", key)
    )
  })
  expected_canonical_content <- vapply(seq_along(keys), function(index) {
    .dp_sensitivity_canonical_string(list(
      request = .dp_sensitivity_identity_request(
        fit_truth[[index]][["request", exact = TRUE]]
      ),
      diagnostics_requested = scenarios[[
        "diagnostics_requested", exact = TRUE
      ]][[index]],
      weight_target = x[["fit_evidence", exact = TRUE]][[
        keys[[index]], exact = TRUE
      ]][["weight_target", exact = TRUE]]
    ))
  }, character(1))
  expected_confidence <- vapply(seq_along(keys), function(index) {
    request <- fit_truth[[index]][["request", exact = TRUE]]
    if (identical(
      fit_truth[[index]][["input_provenance", exact = TRUE]][[
        "target_route", exact = TRUE
      ]],
      "qualitative_confidence"
    )) {
      request[["confidence", exact = TRUE]]
    } else {
      NA_character_
    }
  }, character(1))
  .dpprior_schema_require(
    identical(canonical_content, unname(expected_canonical_content)) &&
      identical(
        scenarios[["effective_method", exact = TRUE]],
        vapply(fit_truth, function(one) {
          one[["input_provenance", exact = TRUE]][[
            "requested_method", exact = TRUE
          ]]
        }, character(1))
      ) && identical(
        scenarios[["method_explicit", exact = TRUE]],
        vapply(fit_truth, function(one) {
          one[["input_provenance", exact = TRUE]][[
            "method_explicit", exact = TRUE
          ]]
        }, logical(1))
      ) && identical(
        scenarios[["confidence_explicit", exact = TRUE]],
        vapply(fit_truth, function(one) {
          one[["input_provenance", exact = TRUE]][[
            "confidence_explicit", exact = TRUE
          ]]
        }, logical(1))
      ) && identical(
        scenarios[["effective_confidence", exact = TRUE]],
        unname(expected_confidence)
      ),
    "sensitivity_canonical_content", paste0(path, ".scenarios"),
    paste(
      "the producer canonical serialization of the full normalized request,",
      "diagnostics flag, and weight target plus exact effective provenance"
    ),
    scenarios
  )
  result_status <- x[["scenario_results", exact = TRUE]][["status", exact = TRUE]]
  result_usable <- x[["scenario_results", exact = TRUE]][["usable", exact = TRUE]]
  result_verified <- x[["scenario_results", exact = TRUE]][["verified", exact = TRUE]]
  results <- x[["scenario_results", exact = TRUE]]
  .dpprior_schema_require(
    is.character(result_status) && !anyNA(result_status) &&
      all(result_status %in% c(
        "converged", "boundary", "approximate", "infeasible", "failed"
      )) && is.logical(result_usable) && !anyNA(result_usable) &&
      is.logical(result_verified) && !anyNA(result_verified),
    "sensitivity_result_status", paste0(path, ".scenario_results"),
    "closed row statuses with ordinary usable/verified flags",
    x[["scenario_results", exact = TRUE]]
  )
  sensitivity_status_rank <- c(
    converged = 1L, boundary = 2L, approximate = 3L,
    infeasible = 4L, failed = 5L
  )
  for (index in seq_along(keys)) {
    .dpprior_validate_status_record(
      list(
        status = result_status[[index]], usable = result_usable[[index]],
        verified = result_verified[[index]], message = ""
      ), sprintf("%s.scenario_results[%d,]", path, index)
    )
    authority <- fit_truth[[index]]
    .dpprior_schema_require(
      identical(authority[["J", exact = TRUE]], sensitivity_J) &&
        sensitivity_status_rank[[result_status[[index]]]] >=
          sensitivity_status_rank[[authority[["status", exact = TRUE]]]] &&
        (!result_usable[[index]] || authority[["usable", exact = TRUE]]) &&
        (!result_verified[[index]] || authority[["verified", exact = TRUE]]),
      "sensitivity_fit_status_authority",
      sprintf("%s.scenario_results[%d,]", path, index),
      paste(
        "J identity and a scenario row no less conservative than its retained",
        "fit status/usable/verified evidence"
      ),
      list(
        row = c(
          status = result_status[[index]], usable = result_usable[[index]],
          verified = result_verified[[index]]
        ),
        fit = authority[c("status", "usable", "verified")]
      )
    )
    if (result_status[[index]] %in% c("converged", "boundary") &&
        (result_usable[[index]] || result_verified[[index]])) {
      essential <- vapply(
        c("a", "b", "E_alpha", "CV_alpha", "E_K_J", "Var_K_J", "CV_K_J"),
        function(field) results[[field, exact = TRUE]][[index]], numeric(1)
      )
      .dpprior_schema_require(
        all(is.finite(essential)), "sensitivity_usable_evidence",
        sprintf("%s.scenario_results[%d,]", path, index),
        "finite calibration and K evidence for every usable/verified row",
        essential
      )
    }
  }
  for (metric in .DPPRIOR_SENSITIVITY_METRICS) {
    metric_values <- x[["scenario_results", exact = TRUE]][[
      metric, exact = TRUE
    ]]
    .dpprior_schema_require(
      is.numeric(metric_values) &&
        all(is.finite(metric_values) |
              (is.na(metric_values) & !is.nan(metric_values))),
      "sensitivity_wide_metric", paste0(path, ".scenario_results.", metric),
      "finite numeric values or canonical NA-unavailable values", metric_values
    )
  }
  probability_metrics <- c(
    "interval_requested", "interval_achieved", "interval_left_tail",
    "interval_right_tail", "E_W_SB", "P_W_SB_gt_50", "P_W_SB_gt_90",
    "P_W_max_gt_50", "P_W_max_gt_90",
    "P_W_max_gt_50_lower_bound", "P_W_max_gt_50_upper_bound",
    "P_W_max_gt_90_lower_bound", "P_W_max_gt_90_upper_bound", "E_rho"
  )
  for (metric in probability_metrics) {
    finite <- is.finite(results[[metric, exact = TRUE]])
    .dpprior_schema_require(
      all(!finite | (results[[metric, exact = TRUE]] >= 0 &
                       results[[metric, exact = TRUE]] <= 1)),
      "sensitivity_probability_metric",
      paste0(path, ".scenario_results.", metric),
      "finite probability/expectation values in [0,1]", results[[metric]]
    )
  }
  finite_interval_residual <- is.finite(
    results[["interval_residual", exact = TRUE]]
  )
  .dpprior_schema_require(
    all(!finite_interval_residual |
          abs(results[["interval_residual", exact = TRUE]]) <= 1),
    "sensitivity_interval_metric",
    paste0(path, ".scenario_results.interval_residual"),
    "a finite interval coverage residual in [-1,1] or canonical NA",
    results[["interval_residual", exact = TRUE]]
  )
  close_numeric <- function(x_value, y_value) {
    is.finite(x_value) && is.finite(y_value) &&
      abs(x_value - y_value) <=
        64 * .Machine$double.eps * max(1, abs(x_value), abs(y_value))
  }
  for (index in seq_along(keys)) {
    a <- results[["a", exact = TRUE]][[index]]
    b <- results[["b", exact = TRUE]][[index]]
    authority <- fit_truth[[index]]
    .dpprior_schema_require(
      (!is.finite(a) || a > 0) && (!is.finite(b) || b > 0),
      "sensitivity_parameter_domain",
      sprintf("%s.scenario_results[%d,]", path, index),
      "positive finite a/b values or canonical NA", c(a = a, b = b)
    )
    if (is.null(authority[["parameters", exact = TRUE]])) {
      unavailable_metrics <- c(
        "a", "b", "E_alpha", "CV_alpha", "E_K_J", "Var_K_J", "CV_K_J",
        "E_W_SB", "P_W_SB_gt_50", "P_W_SB_gt_90", "P_W_max_gt_50",
        "P_W_max_gt_90", "P_W_max_gt_50_lower_bound",
        "P_W_max_gt_50_upper_bound", "P_W_max_gt_90_lower_bound",
        "P_W_max_gt_90_upper_bound", "E_rho"
      )
      .dpprior_schema_require(
        all(vapply(unavailable_metrics, function(metric) {
          value <- results[[metric, exact = TRUE]][[index]]
          is.na(value) && !is.nan(value)
        }, logical(1))),
        "sensitivity_fit_metric_availability",
        sprintf("%s.scenario_results[%d,]", path, index),
        paste(
          "canonical NA for every fit/diagnostic metric when retained fit",
          "parameters and evaluator snapshots are unavailable"
        ),
        results[index, unavailable_metrics, drop = FALSE]
      )
    } else {
      expected_a <- authority[["parameters", exact = TRUE]][["a", exact = TRUE]]
      expected_b <- authority[["parameters", exact = TRUE]][["b", exact = TRUE]]
      selected_K <- authority[["selected", exact = TRUE]][["K", exact = TRUE]]
      M_selected <- authority[["evaluator", exact = TRUE]][[
        "M_selected", exact = TRUE
      ]]
      expected_direct <- c(
        a = expected_a,
        b = expected_b,
        E_alpha = expected_a / expected_b,
        CV_alpha = 1 / sqrt(expected_a),
        E_K_J = selected_K[["mean", exact = TRUE]],
        Var_K_J = selected_K[["variance", exact = TRUE]],
        CV_K_J = sqrt(selected_K[["variance", exact = TRUE]]) /
          selected_K[["mean", exact = TRUE]]
      )
      direct_identity <- vapply(names(expected_direct), function(metric) {
        close_numeric(
          results[[metric, exact = TRUE]][[index]], expected_direct[[metric]]
        )
      }, logical(1))
      .dpprior_schema_require(
        all(direct_identity), "sensitivity_fit_metric_truth",
        sprintf("%s.scenario_results[%d,]", path, index),
        paste(
          "a/b identical to retained fit evidence and alpha/K metrics freshly",
          "recomputed from its selected fixed-order marginal PMF"
        ),
        list(recorded = results[index, names(expected_direct), drop = FALSE],
             expected = expected_direct)
      )
      diagnostic_names <- c(
        "E_W_SB", "P_W_SB_gt_50", "P_W_SB_gt_90", "P_W_max_gt_50",
        "P_W_max_gt_90", "P_W_max_gt_50_lower_bound",
        "P_W_max_gt_50_upper_bound", "P_W_max_gt_90_lower_bound",
        "P_W_max_gt_90_upper_bound", "E_rho"
      )
      if (scenarios[["diagnostics_requested", exact = TRUE]][[index]]) {
        wmax_50 <- unclass(wmax_tail_bounds(0.5, a = expected_a, b = expected_b))
        wmax_90 <- unclass(wmax_tail_bounds(0.9, a = expected_a, b = expected_b))
        expected_diagnostics <- c(
          E_W_SB = as.numeric(mean_w1(expected_a, expected_b, M_selected)),
          P_W_SB_gt_50 = as.numeric(.diagnostic_wsb_tail(
            0.5, expected_a, expected_b
          )),
          P_W_SB_gt_90 = as.numeric(.diagnostic_wsb_tail(
            0.9, expected_a, expected_b
          )),
          P_W_max_gt_50 = NA_real_, P_W_max_gt_90 = NA_real_,
          P_W_max_gt_50_lower_bound = wmax_50[["lower_bound", exact = TRUE]],
          P_W_max_gt_50_upper_bound = wmax_50[["upper_bound", exact = TRUE]],
          P_W_max_gt_90_lower_bound = wmax_90[["lower_bound", exact = TRUE]],
          P_W_max_gt_90_upper_bound = wmax_90[["upper_bound", exact = TRUE]],
          E_rho = as.numeric(mean_rho(expected_a, expected_b, M_selected))
        )
        diagnostic_identity <- vapply(names(expected_diagnostics), function(metric) {
          expected_value <- expected_diagnostics[[metric]]
          recorded_value <- results[[metric, exact = TRUE]][[index]]
          if (is.na(expected_value) && !is.nan(expected_value)) {
            is.na(recorded_value) && !is.nan(recorded_value)
          } else {
            close_numeric(recorded_value, expected_value)
          }
        }, logical(1))
        .dpprior_schema_require(
          all(diagnostic_identity), "sensitivity_diagnostic_metric_truth",
          sprintf("%s.scenario_results[%d,]", path, index),
          paste(
            "fresh fixed-order W_SB/rho values, analytic certified W_max",
            "bounds, and unavailable W_max points without typed backend evidence"
          ),
          list(
            recorded = results[index, names(expected_diagnostics), drop = FALSE],
            expected = expected_diagnostics
          )
        )
      } else {
        .dpprior_schema_require(
          all(vapply(diagnostic_names, function(metric) {
            value <- results[[metric, exact = TRUE]][[index]]
            is.na(value) && !is.nan(value)
          }, logical(1))),
          "sensitivity_diagnostic_metric_truth",
          sprintf("%s.scenario_results[%d,]", path, index),
          "canonical NA diagnostic metrics when diagnostics were not requested",
          results[index, diagnostic_names, drop = FALSE]
        )
      }
    }
    .dpprior_schema_require(
      !(is.finite(results[["E_alpha", exact = TRUE]][[index]]) ||
          is.finite(results[["CV_alpha", exact = TRUE]][[index]])) ||
        (is.finite(a) && is.finite(b)),
      "sensitivity_alpha_availability",
      sprintf("%s.scenario_results[%d,]", path, index),
      "finite E_alpha/CV_alpha claims only with finite positive a/b",
      results[index, c("a", "b", "E_alpha", "CV_alpha"), drop = FALSE]
    )
    if (is.finite(a) && is.finite(b)) {
      .dpprior_schema_require(
        a > 0 && b > 0 &&
          close_numeric(
            results[["E_alpha", exact = TRUE]][[index]], a / b
          ) && close_numeric(
            results[["CV_alpha", exact = TRUE]][[index]], 1 / sqrt(a)
          ),
        "sensitivity_alpha_identity",
        sprintf("%s.scenario_results[%d,]", path, index),
        "positive a/b and exact E_alpha=a/b, CV_alpha=1/sqrt(a)",
        results[index, c("a", "b", "E_alpha", "CV_alpha"), drop = FALSE]
      )
      for (entry in list(
        list(metric = "P_W_SB_gt_50", threshold = 0.5),
        list(metric = "P_W_SB_gt_90", threshold = 0.9)
      )) {
        recorded_tail <- results[[entry$metric, exact = TRUE]][[index]]
        if (is.finite(recorded_tail)) {
          expected_tail <- as.numeric(.diagnostic_wsb_tail(
            entry$threshold, a, b
          ))
          .dpprior_schema_require(
            close_numeric(recorded_tail, expected_tail),
            "sensitivity_W_SB_tail_identity",
            sprintf("%s.scenario_results[%d,].%s", path, index, entry$metric),
            "the exact analytic W_SB tail probability from retained a/b",
            recorded_tail
          )
        }
      }
    }
    diagnostics_metrics <- c(
      "E_W_SB", "P_W_SB_gt_50", "P_W_SB_gt_90", "P_W_max_gt_50",
      "P_W_max_gt_90", "P_W_max_gt_50_lower_bound",
      "P_W_max_gt_50_upper_bound", "P_W_max_gt_90_lower_bound",
      "P_W_max_gt_90_upper_bound", "E_rho"
    )
    if (!scenarios[["diagnostics_requested", exact = TRUE]][[index]]) {
      .dpprior_schema_require(
        all(vapply(
          diagnostics_metrics,
          function(metric) is.na(results[[metric, exact = TRUE]][[index]]) &&
            !is.nan(results[[metric, exact = TRUE]][[index]]),
          logical(1)
        )),
        "sensitivity_diagnostics_request",
        sprintf("%s.scenario_results[%d,]", path, index),
        "no diagnostic metric claims when diagnostics were not requested",
        results[index, diagnostics_metrics, drop = FALSE]
      )
    } else if (result_usable[[index]] || result_verified[[index]]) {
      required_diagnostics <- c(
        "E_W_SB", "P_W_SB_gt_50", "P_W_SB_gt_90", "E_rho"
      )
      .dpprior_schema_require(
        all(vapply(
          required_diagnostics,
          function(metric) is.finite(results[[metric, exact = TRUE]][[index]]),
          logical(1)
        )),
        "sensitivity_diagnostics_request",
        sprintf("%s.scenario_results[%d,]", path, index),
        "finite core diagnostic metrics for a usable requested diagnostic row",
        results[index, required_diagnostics, drop = FALSE]
      )
    }
    K_values <- vapply(
      c("E_K_J", "Var_K_J", "CV_K_J"),
      function(field) results[[field, exact = TRUE]][[index]], numeric(1)
    )
    .dpprior_schema_require(
      (!is.finite(K_values[["E_K_J"]]) ||
         (K_values[["E_K_J"]] >= 1 &&
            K_values[["E_K_J"]] <= sensitivity_J)) &&
        (!is.finite(K_values[["Var_K_J"]]) ||
           (K_values[["Var_K_J"]] >= 0 &&
              K_values[["Var_K_J"]] <= (sensitivity_J - 1)^2 / 4)) &&
        (!is.finite(K_values[["CV_K_J"]]) ||
           K_values[["CV_K_J"]] >= 0) &&
        (!is.finite(K_values[["CV_K_J"]]) ||
           (is.finite(K_values[["E_K_J"]]) &&
              is.finite(K_values[["Var_K_J"]]))),
      "sensitivity_K_domain",
      sprintf("%s.scenario_results[%d,]", path, index),
      paste(
        "each available K moment in its finite-support domain and a finite",
        "CV only when its mean and variance are available"
      ),
      K_values
    )
    if (all(is.finite(K_values))) {
      .dpprior_schema_require(
        K_values[["E_K_J"]] >= 1 &&
          K_values[["E_K_J"]] <= sensitivity_J &&
          K_values[["Var_K_J"]] >= 0 &&
          K_values[["Var_K_J"]] <=
            (K_values[["E_K_J"]] - 1) *
              (sensitivity_J - K_values[["E_K_J"]]) &&
          close_numeric(
            K_values[["CV_K_J"]],
            sqrt(K_values[["Var_K_J"]]) / K_values[["E_K_J"]]
          ),
        "sensitivity_K_identity",
        sprintf("%s.scenario_results[%d,]", path, index),
        paste(
          "finite-support K moments and",
          "CV_K_J=sqrt(Var_K_J)/E_K_J on the same finite source"
        ), K_values
      )
    }
    interval_values <- vapply(
      c("interval_requested", "interval_achieved", "interval_residual"),
      function(field) results[[field, exact = TRUE]][[index]], numeric(1)
    )
    if (all(is.finite(interval_values))) {
      .dpprior_schema_require(
        close_numeric(
          interval_values[["interval_residual"]],
          interval_values[["interval_achieved"]] -
            interval_values[["interval_requested"]]
        ),
        "sensitivity_interval_identity",
        sprintf("%s.scenario_results[%d,]", path, index),
        "interval_residual=interval_achieved-interval_requested",
        interval_values
      )
    }
    for (threshold in c("50", "90")) {
      point <- results[[paste0("P_W_max_gt_", threshold), exact = TRUE]][[index]]
      lower <- results[[paste0(
        "P_W_max_gt_", threshold, "_lower_bound"
      ), exact = TRUE]][[index]]
      upper <- results[[paste0(
        "P_W_max_gt_", threshold, "_upper_bound"
      ), exact = TRUE]][[index]]
      .dpprior_schema_require(
        identical(is.finite(lower), is.finite(upper)),
        "sensitivity_W_max_bounds",
        sprintf("%s.scenario_results[%d,].P_W_max_gt_%s", path,
                index, threshold),
        "paired finite lower/upper bounds or paired canonical NA values",
        c(lower = lower, upper = upper)
      )
      if (is.finite(lower) && is.finite(upper)) {
        .dpprior_schema_require(
          lower <= upper &&
            (!is.finite(point) || (point >= lower && point <= upper)),
          "sensitivity_W_max_bounds",
          sprintf("%s.scenario_results[%d,].P_W_max_gt_%s", path,
                  index, threshold),
          "ordered certified bounds containing any retained point estimate",
          c(point = point, lower = lower, upper = upper)
        )
      }
    }
  }
  .dpprior_schema_validate_named_list(
    x[["conditions", exact = TRUE]], paste0(path, ".conditions")
  )
  .dpprior_schema_validate_named_list(
    x[["interval_checks", exact = TRUE]], paste0(path, ".interval_checks")
  )
  serializable_conditions <- lapply(
    x[["conditions", exact = TRUE]], without_warning_sequences
  )
  .dpprior_schema_validate_plain_record_value(
    serializable_conditions, paste0(path, ".conditions")
  )
  .dpprior_schema_validate_plain_record_value(
    x[["interval_checks", exact = TRUE]], paste0(path, ".interval_checks")
  )
  .dpprior_schema_require(
    identical(names(x[["conditions", exact = TRUE]]), keys) &&
      identical(names(x[["interval_checks", exact = TRUE]]), keys),
    "sensitivity_key_identity", path,
    "condition and interval ledgers keyed exactly like scenario tables",
    list(
      conditions = names(x[["conditions", exact = TRUE]]),
      intervals = names(x[["interval_checks", exact = TRUE]])
    )
  )
  for (key in keys) {
    scenario_index <- match(key, keys)
    condition_record <- x[["conditions", exact = TRUE]][[key, exact = TRUE]]
    .dpprior_sensitivity_validate_condition_evidence(
      condition_record, paste0(path, ".conditions.", key)
    )
    .dpprior_schema_require(
      length(condition_record[["diagnostic_warnings", exact = TRUE]]) == 0L,
      "sensitivity_diagnostic_warning_authority",
      paste0(path, ".conditions.", key, ".diagnostic_warnings"),
      paste(
        "the producer-canonical empty diagnostic-warning slot; captured",
        "calibration-call warnings belong only to calibration_warnings"
      ),
      condition_record[["diagnostic_warnings", exact = TRUE]]
    )
    retained_conditions <- x[["fit_evidence", exact = TRUE]][[
      key, exact = TRUE
    ]][["condition_evidence", exact = TRUE]]
    .dpprior_schema_require(
      identical(condition_record, retained_conditions),
      "sensitivity_condition_authority",
      paste0(path, ".conditions.", key),
      paste(
        "the ordered typed conditions/warnings retained with the canonical",
        "per-scenario fit evidence"
      ),
      condition_record
    )
    condition_slots <- c("calibration", "diagnostics", "target", "interval")
    has_condition <- !vapply(
      condition_record[condition_slots], is.null, logical(1)
    )
    condition_codes <- vapply(condition_slots, function(slot) {
      code <- .dpprior_sensitivity_condition_contract(
        condition_record[[slot, exact = TRUE]]
      )
      if (is.null(code)) NA_character_ else code
    }, character(1))
    .dpprior_schema_require(
      all(!has_condition | !is.na(condition_codes)),
      "sensitivity_condition_contract",
      paste0(path, ".conditions.", key),
      paste(
        "closed producer class/code/class-chain summaries, with no raw",
        "condition objects or open-ended class/code vocabulary"
      ),
      condition_record
    )
    authority <- fit_truth[[scenario_index]]
    fit_status <- authority[["status", exact = TRUE]]
    fit_method <- authority[["method", exact = TRUE]]
    calibration_code <- unname(condition_codes[["calibration"]])
    diagnostics_code <- unname(condition_codes[["diagnostics"]])
    warnings_empty <-
      length(condition_record[["calibration_warnings", exact = TRUE]]) == 0L &&
      length(condition_record[["diagnostic_warnings", exact = TRUE]]) == 0L
    expected_target_condition <- authority[["target_condition", exact = TRUE]]
    target_chain_ok <- if (is.null(expected_target_condition)) {
      !has_condition[["target"]] && !has_condition[["interval"]]
    } else {
      identical(
        condition_record[["target", exact = TRUE]],
        expected_target_condition
      ) && !has_condition[["interval"]]
    }
    native_condition_ok <- if (fit_status %in% c("converged", "boundary")) {
      !has_condition[["calibration"]] &&
        (!has_condition[["diagnostics"]] ||
           (identical(result_status[[scenario_index]], "failed") &&
              identical(diagnostics_code, "sensitivity_diagnostic_contract")))
    } else if (identical(fit_status, "approximate") &&
               identical(fit_method, "A1")) {
      !has_condition[["calibration"]] &&
        (!has_condition[["diagnostics"]] ||
           (identical(result_status[[scenario_index]], "approximate") &&
              identical(diagnostics_code, "fit_diagnostics_approximate")) ||
           (identical(result_status[[scenario_index]], "failed") &&
              identical(diagnostics_code, "sensitivity_diagnostic_contract")))
    } else if (identical(fit_status, "approximate")) {
      identical(calibration_code, "calibration_unusable") &&
        (!has_condition[["diagnostics"]] ||
           (identical(result_status[[scenario_index]], "failed") &&
              identical(diagnostics_code, "sensitivity_diagnostic_contract")))
    } else if (identical(fit_status, "infeasible")) {
      if (identical(sensitivity_J, 1L) && is.null(
        authority[["target_condition", exact = TRUE]]
      )) {
        identical(calibration_code, "calibration_nonidentifiable_j1") &&
          target_chain_ok
      } else {
        !is.null(expected_target_condition) && identical(
          condition_record[["calibration", exact = TRUE]],
          expected_target_condition
        ) && target_chain_ok
      }
    } else {
      (calibration_code %in% c(
        "calibration_unusable", "sensitivity_backend_contract"
      ) || (
        calibration_code %in% c(
          "maxent_root_not_bracketed", "maxent_root_solver_error",
          "target_failed"
        ) && identical(
          authority[["input_provenance", exact = TRUE]][[
            "target_route", exact = TRUE
          ]],
          "interval"
        ) && !is.null(expected_target_condition) && identical(
          condition_record[["calibration", exact = TRUE]],
          expected_target_condition
        )
      ) || identical(
        diagnostics_code, "sensitivity_diagnostic_contract"
      )) && target_chain_ok
    }
    row_condition_ok <- if (identical(
      result_status[[scenario_index]], "failed"
    )) {
      any(has_condition)
    } else if (identical(result_status[[scenario_index]], "infeasible")) {
      identical(fit_status, "infeasible") && has_condition[["calibration"]]
    } else {
      TRUE
    }
    warning_condition_ok <- if (!warnings_empty) {
      identical(result_status[[scenario_index]], "failed") &&
        has_condition[["diagnostics"]] &&
        identical(diagnostics_code, "sensitivity_diagnostic_contract")
    } else {
      TRUE
    }
    .dpprior_schema_require(
      native_condition_ok && target_chain_ok && row_condition_ok &&
        warning_condition_ok,
      "sensitivity_condition_status",
      paste0(path, ".conditions.", key),
      paste(
        "the exact method/status condition matrix: condition-free successful",
        "fits, typed A1/A2 approximations, certified infeasibility, or typed",
        "failed evidence; warning-only success is forbidden"
      ),
      list(
        fit_status = fit_status, fit_method = fit_method,
        row_status = result_status[[scenario_index]],
        conditions = condition_record
      )
    )
    interval_record <- x[["interval_checks", exact = TRUE]][[
      key, exact = TRUE
    ]]
    .dpprior_schema_exact_names(
      interval_record,
      c(
        "requested", "selected", "verification", "status", "source",
        "usable", "verified", "reason"
      ), paste0(path, ".interval_checks.", key)
    )
    for (field in c("usable", "verified")) {
      .dpprior_schema_validate_scalar_logical(
        interval_record[[field, exact = TRUE]],
        paste0(path, ".interval_checks.", key, ".", field)
      )
    }
    if (is.null(interval_record[["status", exact = TRUE]])) {
      .dpprior_schema_require(
        all(vapply(
          interval_record[c(
            "requested", "selected", "verification", "status", "source"
          )], is.null, logical(1)
        )) && !interval_record[["usable", exact = TRUE]] &&
          !interval_record[["verified", exact = TRUE]] &&
          !identical(
            fit_truth[[scenario_index]][["target", exact = TRUE]][[
              "kind", exact = TRUE
            ]],
            "interval"
          ) &&
          identical(
            interval_record[["reason", exact = TRUE]],
            "not_interval_scenario"
          ),
        "sensitivity_interval_na", paste0(path, ".interval_checks.", key),
        "the exact non-interval unavailable record", interval_record
      )
      .dpprior_schema_require(
        all(vapply(
          c(
            "interval_requested", "interval_achieved", "interval_residual",
            "interval_left_tail", "interval_right_tail"
          ),
          function(metric) {
            value <- results[[metric, exact = TRUE]][[scenario_index]]
            is.na(value) && !is.nan(value)
          }, logical(1)
        )),
        "sensitivity_interval_wide_identity",
        paste0(path, ".scenario_results.", key),
        "canonical NA interval metrics for a non-interval scenario",
        results[scenario_index, c(
          "interval_requested", "interval_achieved", "interval_residual",
          "interval_left_tail", "interval_right_tail"
        ), drop = FALSE]
      )
    } else {
      .dpprior_schema_validate_scalar_character(
        interval_record[["status", exact = TRUE]],
        paste0(path, ".interval_checks.", key, ".status")
      )
      .dpprior_schema_validate_scalar_character(
        interval_record[["source", exact = TRUE]],
        paste0(path, ".interval_checks.", key, ".source")
      )
      .dpprior_schema_require(
        interval_record[["status", exact = TRUE]] %in% c(
          "converged", "boundary", "approximate", "infeasible", "failed"
        ) && interval_record[["source", exact = TRUE]] %in% c(
          "wrapper_backcheck_selected",
          "target_infeasibility_certificate", "target_construction_failure"
        ),
        "sensitivity_interval_status",
        paste0(path, ".interval_checks.", key, ".status"),
        "a closed canonical interval status and producer-backed source",
        interval_record[c("status", "source")]
      )
      .dpprior_schema_validate_scalar_character(
        interval_record[["reason", exact = TRUE]],
        paste0(path, ".interval_checks.", key, ".reason"), allow_empty = TRUE
      )
      .dpprior_validate_status_record(
        list(
          status = interval_record[["status", exact = TRUE]],
          usable = interval_record[["usable", exact = TRUE]],
          verified = interval_record[["verified", exact = TRUE]],
          message = interval_record[["reason", exact = TRUE]]
        ),
        paste0(path, ".interval_checks.", key, ".status_record")
      )
      certificate_route <- interval_record[["source", exact = TRUE]] %in%
        c("target_infeasibility_certificate", "target_construction_failure")
      if (certificate_route) {
        requested <- interval_record[["requested", exact = TRUE]]
        .dpprior_sensitivity_validate_interval_authority(
          requested, sensitivity_J,
          paste0(path, ".interval_checks.", key, ".requested")
        )
        authority <- fit_truth[[scenario_index]]
        expected_status <- if (identical(
          interval_record[["source", exact = TRUE]],
          "target_infeasibility_certificate"
        )) "infeasible" else "failed"
        expected_condition <- authority[["target_condition", exact = TRUE]]
        expected_interval_wide <- c(
          interval_requested = requested[["coverage", exact = TRUE]],
          interval_achieved = NA_real_, interval_residual = NA_real_,
          interval_left_tail = NA_real_, interval_right_tail = NA_real_
        )
        recorded_interval_wide <- vapply(
          names(expected_interval_wide),
          function(metric) results[[metric, exact = TRUE]][[scenario_index]],
          numeric(1)
        )
        .dpprior_schema_require(
          identical(
            authority[["target", exact = TRUE]][["kind", exact = TRUE]],
            "interval"
          ) && identical(
            requested,
            authority[["target", exact = TRUE]][["used", exact = TRUE]][[
              "interval", exact = TRUE
            ]]
          ) && identical(authority[["status", exact = TRUE]], expected_status) &&
            identical(authority[["target_status", exact = TRUE]], expected_status) &&
            identical(interval_record[["status", exact = TRUE]], expected_status) &&
            is.null(interval_record[["selected", exact = TRUE]]) &&
            is.null(interval_record[["verification", exact = TRUE]]) &&
            !interval_record[["usable", exact = TRUE]] &&
            identical(
              interval_record[["verified", exact = TRUE]],
              identical(expected_status, "infeasible")
            ) && !is.null(expected_condition) && identical(
              interval_record[["reason", exact = TRUE]],
              expected_condition[["message", exact = TRUE]]
            ) && identical(recorded_interval_wide, expected_interval_wide),
          "sensitivity_interval_target_outcome",
          paste0(path, ".interval_checks.", key),
          paste(
            "an exact PMF-free certified-infeasible or construction-failed",
            "interval route bound to the reconstructed target condition"
          ),
          interval_record
        )
        next
      }
      .dpprior_schema_require(
        sensitivity_status_rank[[result_status[[scenario_index]]]] >=
          sensitivity_status_rank[[interval_record[["status", exact = TRUE]]]],
        "sensitivity_interval_status_aggregation",
        paste0(path, ".scenario_results.", key, ".status"),
        "a scenario status at least as conservative as its interval audit",
        list(
          scenario = result_status[[scenario_index]],
          interval = interval_record[["status", exact = TRUE]]
        )
      )
      requested <- interval_record[["requested", exact = TRUE]]
      selected <- interval_record[["selected", exact = TRUE]]
      .dpprior_sensitivity_validate_interval_authority(
        requested, sensitivity_J,
        paste0(path, ".interval_checks.", key, ".requested")
      )
      .dpprior_schema_require(
        identical(
            fit_truth[[scenario_index]][["target", exact = TRUE]][[
              "kind", exact = TRUE
            ]],
            "interval"
          ) && identical(
            requested,
            fit_truth[[scenario_index]][["target", exact = TRUE]][[
              "used", exact = TRUE
            ]][[
              "interval", exact = TRUE
            ]]
          ),
        "sensitivity_interval_request",
        paste0(path, ".interval_checks.", key, ".requested"),
        paste(
          "a closed interval type/family exactly identical to the retained",
          "normalized per-scenario request authority"
        ),
        requested
      )
      .dpprior_schema_exact_names(
        selected,
        c("coverage", "lower_tail", "upper_tail", "coverage_residual"),
        paste0(path, ".interval_checks.", key, ".selected")
      )
      for (field in names(selected)) {
        .dpprior_schema_validate_finite_scalar(
          selected[[field, exact = TRUE]],
          paste0(path, ".interval_checks.", key, ".selected.", field)
        )
      }
      interval_rounding <- 64 * .Machine$double.eps
      .dpprior_schema_require(
        all(unlist(selected[c(
          "coverage", "lower_tail", "upper_tail"
        )], use.names = FALSE) >= 0) &&
          all(unlist(selected[c(
            "coverage", "lower_tail", "upper_tail"
          )], use.names = FALSE) <= 1) &&
          abs(
            selected[["coverage", exact = TRUE]] +
              selected[["lower_tail", exact = TRUE]] +
              selected[["upper_tail", exact = TRUE]] - 1
          ) <= interval_rounding && close_numeric(
            selected[["coverage_residual", exact = TRUE]],
            selected[["coverage", exact = TRUE]] -
              requested[["coverage", exact = TRUE]]
          ),
        "sensitivity_interval_selected",
        paste0(path, ".interval_checks.", key, ".selected"),
        "probability partition and exact selected coverage residual", selected
      )
      fit_selected <- fit_truth[[scenario_index]][["selected", exact = TRUE]]
      fit_verifier <- fit_truth[[scenario_index]][["verifier", exact = TRUE]]
      .dpprior_schema_require(
        !is.null(fit_selected) && !is.null(fit_verifier),
        "sensitivity_interval_pmf_authority",
        paste0(path, ".fit_evidence.", key),
        "retained selected and independent verifier K PMFs for interval claims",
        fit_truth[[scenario_index]]
      )
      selected_masses <- .dpprior_target_interval_masses(
        fit_selected[["K", exact = TRUE]][["pmf", exact = TRUE]], requested
      )
      verifier_masses <- .dpprior_target_interval_masses(
        fit_verifier[["K", exact = TRUE]][["pmf", exact = TRUE]], requested
      )
      expected_selected_interval <- c(
        coverage = selected_masses[["inside_mass"]],
        lower_tail = selected_masses[["left_mass"]],
        upper_tail = selected_masses[["right_mass"]],
        coverage_residual = selected_masses[["inside_mass"]] -
          requested[["coverage", exact = TRUE]]
      )
      .dpprior_schema_require(
        all(vapply(names(expected_selected_interval), function(field) {
          .dpprior_sensitivity_close_numeric(
            selected[[field, exact = TRUE]],
            expected_selected_interval[[field]]
          )
        }, logical(1))),
        "sensitivity_interval_selected_truth",
        paste0(path, ".interval_checks.", key, ".selected"),
        "the selected interval masses freshly summed from retained selected PMF",
        list(recorded = selected, expected = expected_selected_interval)
      )
      if (identical(
        interval_record[["source", exact = TRUE]],
        "wrapper_backcheck_selected"
      )) {
        interval_verification <- interval_record[["verification", exact = TRUE]]
        .dpprior_schema_exact_names(
          interval_verification,
          c(
            "coverage", "lower_tail", "upper_tail", "coverage_residual",
            "tolerance", "passed", "source"
          ),
          paste0(path, ".interval_checks.", key, ".verification")
        )
        for (field in c(
          "coverage", "lower_tail", "upper_tail", "coverage_residual",
          "tolerance"
        )) {
          .dpprior_schema_validate_finite_scalar(
            interval_verification[[field, exact = TRUE]],
            paste0(
              path, ".interval_checks.", key, ".verification.", field
            ),
            lower = if (identical(field, "coverage_residual")) -1 else 0,
            upper = if (identical(field, "coverage_residual")) 1 else 1
          )
        }
        .dpprior_schema_validate_scalar_logical(
          interval_verification[["passed", exact = TRUE]],
          paste0(path, ".interval_checks.", key, ".verification.passed")
        )
        .dpprior_schema_validate_scalar_character(
          interval_verification[["source", exact = TRUE]],
          paste0(path, ".interval_checks.", key, ".verification.source")
        )
        interval_tolerance <- interval_verification[["tolerance", exact = TRUE]]
        expected_verifier_interval <- c(
          coverage = verifier_masses[["inside_mass"]],
          lower_tail = verifier_masses[["left_mass"]],
          upper_tail = verifier_masses[["right_mass"]],
          coverage_residual = verifier_masses[["inside_mass"]] -
            requested[["coverage", exact = TRUE]]
        )
        verifier_identity <- all(vapply(
          names(expected_verifier_interval), function(field) {
            .dpprior_sensitivity_close_numeric(
              interval_verification[[field, exact = TRUE]],
              expected_verifier_interval[[field]]
            )
          }, logical(1)
        ))
        constraint_pass <- function(masses) {
          if (identical(requested[["type", exact = TRUE]], "hard_bounds")) {
            abs(masses[["left_mass"]]) <= interval_tolerance &&
              abs(masses[["inside_mass"]] - 1) <= interval_tolerance &&
              abs(masses[["right_mass"]]) <= interval_tolerance
          } else if (identical(
            requested[["type", exact = TRUE]], "equal_tail"
          )) {
            tail_target <- (1 - requested[["coverage", exact = TRUE]]) / 2
            abs(masses[["left_mass"]] - tail_target) <= interval_tolerance &&
              abs(masses[["inside_mass"]] -
                    requested[["coverage", exact = TRUE]]) <=
                interval_tolerance &&
              abs(masses[["right_mass"]] - tail_target) <= interval_tolerance
          } else {
            abs(masses[["inside_mass"]] -
                  requested[["coverage", exact = TRUE]]) <= interval_tolerance
          }
        }
        expected_interval_pass <-
          verifier_identity && constraint_pass(selected_masses) &&
          constraint_pass(verifier_masses) &&
          all(abs(c(
            coverage = selected[["coverage", exact = TRUE]] -
              interval_verification[["coverage", exact = TRUE]],
            lower_tail = selected[["lower_tail", exact = TRUE]] -
              interval_verification[["lower_tail", exact = TRUE]],
            upper_tail = selected[["upper_tail", exact = TRUE]] -
              interval_verification[["upper_tail", exact = TRUE]]
          )) <= interval_tolerance)
        .dpprior_schema_require(
          identical(interval_tolerance, 1e-8) &&
            identical(
            interval_verification[["source", exact = TRUE]],
              "independent_interval_backcheck"
            ) && verifier_identity && abs(
              interval_verification[["coverage", exact = TRUE]] +
                interval_verification[["lower_tail", exact = TRUE]] +
                interval_verification[["upper_tail", exact = TRUE]] - 1
            ) <= interval_rounding && close_numeric(
              interval_verification[["coverage_residual", exact = TRUE]],
              interval_verification[["coverage", exact = TRUE]] -
                requested[["coverage", exact = TRUE]]
            ) && identical(
              interval_verification[["passed", exact = TRUE]],
              expected_interval_pass
            ) && if (expected_interval_pass) {
              identical(
                interval_record[["status", exact = TRUE]], "converged"
              ) &&
                interval_record[["usable", exact = TRUE]] &&
                interval_record[["verified", exact = TRUE]] &&
                identical(interval_record[["reason", exact = TRUE]], "")
            } else {
              identical(
                interval_record[["status", exact = TRUE]], "approximate"
              ) &&
                !interval_record[["usable", exact = TRUE]] &&
                !interval_record[["verified", exact = TRUE]] &&
                identical(
                  interval_record[["reason", exact = TRUE]],
                  paste(
                    "Selected and verifier interval constraints did not pass",
                    "fixed tolerance"
                  )
                )
            },
          "sensitivity_interval_verification",
          paste0(path, ".interval_checks.", key, ".verification"),
          paste(
            "an independent fixed-tolerance interval backcheck whose pass",
            "and status quartet are freshly recomputed"
          ),
          interval_verification
        )
      }
      expected_interval_wide <- c(
        interval_requested = requested[["coverage", exact = TRUE]],
        interval_achieved = selected[["coverage", exact = TRUE]],
        interval_residual = selected[["coverage_residual", exact = TRUE]],
        interval_left_tail = selected[["lower_tail", exact = TRUE]],
        interval_right_tail = selected[["upper_tail", exact = TRUE]]
      )
      recorded_interval_wide <- vapply(
        names(expected_interval_wide),
        function(metric) results[[metric, exact = TRUE]][[scenario_index]],
        numeric(1)
      )
      .dpprior_schema_require(
        all(vapply(names(expected_interval_wide), function(metric) {
          close_numeric(
            recorded_interval_wide[[metric]], expected_interval_wide[[metric]]
          )
        }, logical(1))),
        "sensitivity_interval_wide_identity",
        paste0(path, ".scenario_results.", key),
        "the five wide interval metrics bound to the retained interval audit",
        recorded_interval_wide
      )
    }
  }
  .dpprior_schema_validate_plain_data_frame(
    x[["metrics_long", exact = TRUE]], paste0(path, ".metrics_long")
  )
  .dpprior_schema_require(
    identical(
      names(x[["metrics_long", exact = TRUE]]),
      c(
        "scenario_key", "metric", "value", "reason", "source", "component",
        "status", "usable", "verified"
      )
    ),
    "sensitivity_metric_fields", paste0(path, ".metrics_long"),
    "the exact canonical long-metric columns",
    names(x[["metrics_long", exact = TRUE]])
  )
  expected_metric_keys <- rep(keys, each = length(.DPPRIOR_SENSITIVITY_METRICS))
  expected_metrics <- rep(
    .DPPRIOR_SENSITIVITY_METRICS, times = length(keys)
  )
  metrics <- x[["metrics_long", exact = TRUE]]
  values <- metrics[["value", exact = TRUE]]
  reasons <- metrics[["reason", exact = TRUE]]
  finite_values <- is.finite(values)
  unavailable_values <- is.na(values) & !is.nan(values)
  expected_wide_values <- unname(unlist(lapply(seq_along(keys), function(index) {
    vapply(
      .DPPRIOR_SENSITIVITY_METRICS,
      function(metric) results[[metric, exact = TRUE]][[index]], numeric(1)
    )
  }), use.names = FALSE))
  expected_sources <- rep(
    unname(.DPPRIOR_SENSITIVITY_METRIC_SOURCE), times = length(keys)
  )
  expected_components <- rep(
    unname(.DPPRIOR_SENSITIVITY_METRIC_COMPONENT), times = length(keys)
  )
  .dpprior_schema_require(
    identical(metrics[["scenario_key", exact = TRUE]], expected_metric_keys) &&
      identical(metrics[["metric", exact = TRUE]], expected_metrics) &&
      identical(values, expected_wide_values) &&
      is.numeric(values) && all(finite_values | unavailable_values) &&
      is.character(reasons) && length(reasons) == length(values) &&
      all(is.na(reasons[finite_values])) &&
      all(!is.na(reasons[unavailable_values]) &
            nzchar(reasons[unavailable_values])) &&
      is.character(metrics[["source", exact = TRUE]]) &&
      !anyNA(metrics[["source", exact = TRUE]]) &&
      identical(metrics[["source", exact = TRUE]], expected_sources) &&
      is.character(metrics[["component", exact = TRUE]]) &&
      !anyNA(metrics[["component", exact = TRUE]]) &&
      identical(metrics[["component", exact = TRUE]], expected_components) &&
      is.character(metrics[["status", exact = TRUE]]) &&
      !anyNA(metrics[["status", exact = TRUE]]) &&
      is.logical(metrics[["usable", exact = TRUE]]) &&
      !anyNA(metrics[["usable", exact = TRUE]]) &&
      is.logical(metrics[["verified", exact = TRUE]]) &&
      !anyNA(metrics[["verified", exact = TRUE]]),
    "sensitivity_metric_grid", paste0(path, ".metrics_long"),
    paste(
      "the exact scenario-by-22-metric grid, with finite values carrying NA",
      "reasons and unavailable NA values carrying a nonempty reason"
    ), metrics
  )
  .dpprior_schema_validate_plain_data_frame(
    x[["local", exact = TRUE]], paste0(path, ".local")
  )
  .dpprior_schema_require(
    identical(names(x[["local", exact = TRUE]]),
              .DPPRIOR_SENSITIVITY_LOCAL_FIELDS),
    "sensitivity_local_fields", paste0(path, ".local"),
    "the exact canonical local-derivative columns",
    names(x[["local", exact = TRUE]])
  )
  local <- x[["local", exact = TRUE]]
  .dpprior_schema_require(
    nrow(local) == 0L,
    "sensitivity_local_authority", paste0(path, ".local"),
    paste(
      "an empty canonical local table until per-axis request trees and",
      "axis-removed settings identities are retained for recomputation"
    ),
    local
  )
  if (nrow(local) > 0L) {
    character_fields <- c(
      "scenario_key", "axis", "settings_key", "lower_scenario_key",
      "upper_scenario_key", "metric", "component", "method", "reason"
    )
    numeric_fields <- c(
      "axis_value", "lower_value", "upper_value", "derivative"
    )
    .dpprior_schema_require(
      all(vapply(character_fields, function(field) {
        is.character(local[[field, exact = TRUE]])
      }, logical(1))) && all(vapply(numeric_fields, function(field) {
        is.numeric(local[[field, exact = TRUE]])
      }, logical(1))),
      "sensitivity_local_types", paste0(path, ".local"),
      "the exact character/numeric local-derivative column types", local
    )
    finite_derivative <- is.finite(local[["derivative", exact = TRUE]])
    unavailable_derivative <- is.na(local[["derivative", exact = TRUE]]) &
      !is.nan(local[["derivative", exact = TRUE]])
    lower_key <- local[["lower_scenario_key", exact = TRUE]]
    upper_key <- local[["upper_scenario_key", exact = TRUE]]
    expected_local_component <- unname(
      .DPPRIOR_SENSITIVITY_METRIC_COMPONENT[
        local[["metric", exact = TRUE]]
      ]
    )
    .dpprior_schema_require(
      all(local[["scenario_key", exact = TRUE]] %in% keys) &&
        !anyNA(local[["axis", exact = TRUE]]) &&
        all(local[["axis", exact = TRUE]] %in%
              c("mu_K", "var_K", "cv_K", "interval_coverage")) &&
        !anyNA(local[["settings_key", exact = TRUE]]) &&
        all(grepl(
          "^scn_[0-9a-f]{16}$", local[["settings_key", exact = TRUE]]
        )) &&
        all(local[["metric", exact = TRUE]] %in%
              .DPPRIOR_SENSITIVITY_METRICS) &&
        identical(local[["component", exact = TRUE]],
                  expected_local_component) &&
        all(local[["method", exact = TRUE]] ==
              "bracketed_secant_across_nearest_same-setting_neighbors") &&
        all(is.finite(local[["axis_value", exact = TRUE]])) &&
        all(is.finite(local[["lower_value", exact = TRUE]]) |
              (is.na(local[["lower_value", exact = TRUE]]) &
                 !is.nan(local[["lower_value", exact = TRUE]]))) &&
        all(is.finite(local[["upper_value", exact = TRUE]]) |
              (is.na(local[["upper_value", exact = TRUE]]) &
                 !is.nan(local[["upper_value", exact = TRUE]]))) &&
        all(finite_derivative | unavailable_derivative) &&
        all(!anyDuplicated(paste(
          local[["scenario_key", exact = TRUE]],
          local[["axis", exact = TRUE]],
          format(local[["axis_value", exact = TRUE]], digits = 17),
          local[["settings_key", exact = TRUE]],
          local[["metric", exact = TRUE]], sep = "\r"
        ))),
      "sensitivity_local_identity", paste0(path, ".local"),
      "unique producer-shaped local derivative rows on canonical keys/metrics",
      local
    )
    no_bracket <- is.na(lower_key) & is.na(upper_key) &
      is.na(local[["lower_value", exact = TRUE]]) &
      is.na(local[["upper_value", exact = TRUE]])
    valid_bracket <- !is.na(lower_key) & lower_key %in% keys &
      !is.na(upper_key) & upper_key %in% keys & lower_key != upper_key &
      is.finite(local[["lower_value", exact = TRUE]]) &
      is.finite(local[["upper_value", exact = TRUE]]) &
      local[["lower_value", exact = TRUE]] <
        local[["upper_value", exact = TRUE]]
    .dpprior_schema_require(
      all(
        !finite_derivative |
          (valid_bracket & is.na(local[["reason", exact = TRUE]]))
      ) && all(
        !unavailable_derivative |
          ((no_bracket | valid_bracket) &
             !is.na(local[["reason", exact = TRUE]]) &
             nzchar(local[["reason", exact = TRUE]]))
      ),
      "sensitivity_local_availability", paste0(path, ".local"),
      paste(
        "finite derivatives require an ordered canonical bracket and no",
        "reason; unavailable derivatives require either no bracket or a valid",
        "canonical bracket plus a nonempty reason"
      ), local
    )
    if (any(finite_derivative)) {
      derivative_truth <- vapply(which(finite_derivative), function(index) {
        lower_index <- match(lower_key[[index]], keys)
        upper_index <- match(upper_key[[index]], keys)
        metric <- local[["metric", exact = TRUE]][[index]]
        lower_metric <- results[[metric, exact = TRUE]][[lower_index]]
        upper_metric <- results[[metric, exact = TRUE]][[upper_index]]
        expected_derivative <- (upper_metric - lower_metric) /
          (local[["upper_value", exact = TRUE]][[index]] -
             local[["lower_value", exact = TRUE]][[index]])
        is.finite(expected_derivative) && close_numeric(
          local[["derivative", exact = TRUE]][[index]], expected_derivative
        )
      }, logical(1))
      .dpprior_schema_require(
        all(derivative_truth), "sensitivity_local_derivative",
        paste0(path, ".local.derivative"),
        "the bracketed secant freshly recomputed from canonical wide metrics",
        local[["derivative", exact = TRUE]]
      )
    }
  }
  for (field in c("global", "metadata")) {
    .dpprior_schema_validate_named_list(
      x[[field]], paste0(path, ".", field)
    )
    .dpprior_schema_validate_plain_record_value(
      x[[field, exact = TRUE]], paste0(path, ".", field)
    )
  }
  .dpprior_schema_exact_names(
    x[["global", exact = TRUE]],
    c("scenario_count", "converged_count", "failed_count", "metric_count"),
    paste0(path, ".global")
  )
  .dpprior_schema_exact_names(
    x[["metadata", exact = TRUE]], "J", paste0(path, ".metadata")
  )
  .dpprior_schema_require(
    identical(x[["global", exact = TRUE]][["scenario_count", exact = TRUE]],
              as.integer(length(keys))) &&
      identical(x[["global", exact = TRUE]][["converged_count", exact = TRUE]],
                as.integer(sum(result_status %in% c("converged", "boundary")))) &&
      identical(x[["global", exact = TRUE]][["failed_count", exact = TRUE]],
                as.integer(sum(result_status %in% c("failed", "infeasible")))) &&
      identical(x[["global", exact = TRUE]][["metric_count", exact = TRUE]],
                as.integer(nrow(metrics))),
    "sensitivity_global", paste0(path, ".global"),
    "global counts recomputed from canonical rows and metric grid",
    x[["global", exact = TRUE]]
  )
  invisible(TRUE)
}


.dpprior_expected_diagnostics_evidence <- function(raw) {
  parameters <- raw[["parameters", exact = TRUE]]
  J <- raw[["J", exact = TRUE]]
  orders <- raw[["computation", exact = TRUE]][["orders", exact = TRUE]]
  M_selected <- orders[["M_selected", exact = TRUE]]
  M_verification <- orders[["M_verification_used", exact = TRUE]]
  .dpprior_schema_require(
    !is.null(parameters), "diagnostic_parameters", "result.parameters",
    "a finite positive canonical a/b parameter record", parameters
  )
  .dpprior_validate_parameters(parameters, "result.parameters")
  .dpprior_schema_require(
    .dpprior_schema_is_count(J, 1L), "diagnostic_J", "result.J",
    "a positive integer J before any numerical diagnostic call", J
  )
  .dpprior_schema_require(
    .dpprior_schema_is_count(M_selected, 1L) &&
      .dpprior_schema_is_count(M_verification, 1L),
    "diagnostic_orders", "result.computation.orders",
    paste(
      "positive integer M_selected and M_verification_used before any",
      "numerical diagnostic call"
    ),
    list(M_selected = M_selected, M_verification_used = M_verification)
  )
  tolerances <- raw[["tolerances", exact = TRUE]][[
    "diagnostics", exact = TRUE
  ]]
  .dpprior_schema_exact_names(
    raw[["tolerances", exact = TRUE]], "diagnostics", "result.tolerances"
  )
  .dpprior_schema_exact_names(
    raw[["residuals", exact = TRUE]], "diagnostics", "result.residuals"
  )
  .dpprior_schema_exact_names(
    tolerances, c("absolute", "relative", "pmf_mass", "refinement"),
    "result.tolerances.diagnostics"
  )
  for (field in c("absolute", "relative", "pmf_mass")) {
    .dpprior_schema_validate_finite_scalar(
      tolerances[[field, exact = TRUE]],
      paste0("result.tolerances.diagnostics.", field), lower = 0
    )
  }
  .dpprior_schema_require(
    tolerances[["absolute", exact = TRUE]] <= 1e-10 &&
      tolerances[["relative", exact = TRUE]] <= 1e-8 &&
      identical(tolerances[["pmf_mass", exact = TRUE]], .TOL_PMF_SUM),
    "diagnostic_truth_tolerances", "result.tolerances.diagnostics",
    paste(
      "absolute <= 1e-10, relative <= 1e-8, and the fixed canonical",
      "PMF-mass tolerance"
    ), tolerances
  )
  controls <- raw[["computation", exact = TRUE]][["used", exact = TRUE]][[
    "controls", exact = TRUE
  ]]
  .dpprior_schema_exact_names(
    controls, c("absolute_tolerance", "relative_tolerance", "pmf_mass_tolerance"),
    "result.computation.used.controls"
  )
  .dpprior_schema_require(
    identical(controls[["absolute_tolerance", exact = TRUE]],
              tolerances[["absolute", exact = TRUE]]) &&
      identical(controls[["relative_tolerance", exact = TRUE]],
                tolerances[["relative", exact = TRUE]]) &&
      identical(controls[["pmf_mass_tolerance", exact = TRUE]],
                tolerances[["pmf_mass", exact = TRUE]]),
    "diagnostic_truth_tolerances", "result.computation.used.controls",
    "exact identity with the central diagnostic tolerance authority", controls
  )
  .dpprior_schema_require(
    identical(
      raw[["computation", exact = TRUE]][["request", exact = TRUE]][[
        "controls", exact = TRUE
      ]], controls
    ),
    "diagnostic_truth_tolerances", "result.computation.request.controls",
    "the retained requested controls identical to the fixed used controls",
    raw[["computation", exact = TRUE]][["request", exact = TRUE]][[
      "controls", exact = TRUE
    ]]
  )
  a <- parameters[["a", exact = TRUE]]
  b <- parameters[["b", exact = TRUE]]
  K <- .get_K_pmf_support(
    J, a, b, M = M_selected, M_verify = M_verification,
    abs_tol = tolerances[["absolute", exact = TRUE]],
    rel_tol = tolerances[["relative", exact = TRUE]]
  )
  selected_pmf <- unname(as.numeric(K[["pmf", exact = TRUE]]))
  verifier_pmf <- unname(as.numeric(K[["verification_pmf", exact = TRUE]]))
  selected_K <- .dpprior_target_pmf_moments(selected_pmf)
  verifier_K <- .dpprior_target_pmf_moments(verifier_pmf)
  selected_weight <- as.numeric(mean_w1(a, b, M_selected))
  verifier_weight <- as.numeric(mean_w1(a, b, M_verification))
  selected_rho <- c(
    mean = as.numeric(mean_rho(a, b, M_selected)),
    variance = as.numeric(var_rho(a, b, M_selected))
  )
  verifier_rho <- c(
    mean = as.numeric(mean_rho(a, b, M_verification)),
    variance = as.numeric(var_rho(a, b, M_verification))
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
    tolerances[["absolute", exact = TRUE]] +
      tolerances[["relative", exact = TRUE]] *
        max(abs(selected), abs(verifier), 1)
  }
  refinement_tolerance <- c(
    K.mean = scalar_tolerance(
      selected_K[["mean"]], verifier_K[["mean"]]
    ),
    K.variance = scalar_tolerance(
      selected_K[["variance"]], verifier_K[["variance"]]
    ),
    K.pmf_l1 = tolerances[["absolute", exact = TRUE]] +
      tolerances[["relative", exact = TRUE]],
    weights.mean = scalar_tolerance(selected_weight, verifier_weight),
    coclustering.mean = scalar_tolerance(
      selected_rho[["mean"]], verifier_rho[["mean"]]
    ),
    coclustering.variance = scalar_tolerance(
      selected_rho[["variance"]], verifier_rho[["variance"]]
    )
  )
  .dpprior_schema_require(
    is.numeric(tolerances[["refinement", exact = TRUE]]) &&
      !is.object(tolerances[["refinement", exact = TRUE]]) &&
      is.null(dim(tolerances[["refinement", exact = TRUE]])) &&
      .dpprior_schema_has_only_attributes(
        tolerances[["refinement", exact = TRUE]], "names"
      ) && identical(
        tolerances[["refinement", exact = TRUE]], refinement_tolerance
      ),
    "diagnostic_refinement_tolerances",
    "result.tolerances.diagnostics.refinement",
    "the named refinement tolerances recomputed from the central formula",
    tolerances[["refinement", exact = TRUE]]
  )
  component_pass <- c(
    alpha = TRUE,
    K = all(delta[c("K.mean", "K.variance", "K.pmf_l1")] <=
              refinement_tolerance[c("K.mean", "K.variance", "K.pmf_l1")]) &&
      abs(sum(selected_pmf) - 1) <= tolerances[["pmf_mass", exact = TRUE]] &&
      abs(sum(verifier_pmf) - 1) <= tolerances[["pmf_mass", exact = TRUE]],
    weights = delta[["weights.mean"]] <=
      refinement_tolerance[["weights.mean"]],
    coclustering = all(delta[c(
      "coclustering.mean", "coclustering.variance"
    )] <= refinement_tolerance[c(
      "coclustering.mean", "coclustering.variance"
    )])
  )
  list(
    selected = list(
      alpha = list(mean = a / b, CV = 1 / sqrt(a)),
      K = list(
        mean = unname(selected_K[["mean"]]),
        variance = unname(selected_K[["variance"]]), estimand = "K_J",
        source = "fresh_diagnostics_selected_order", M = M_selected,
        pmf = selected_pmf
      ),
      weights = list(mean = selected_weight),
      coclustering = as.list(selected_rho)
    ),
    verifier = list(
      alpha = list(mean = a / b, CV = 1 / sqrt(a)),
      K = list(
        mean = unname(verifier_K[["mean"]]),
        variance = unname(verifier_K[["variance"]]), estimand = "K_J",
        source = "fresh_diagnostics_verifier_evidence", M = M_verification,
        pmf = verifier_pmf
      ),
      weights = list(mean = verifier_weight),
      coclustering = as.list(verifier_rho)
    ),
    delta = delta, refinement_tolerance = refinement_tolerance,
    component_pass = component_pass
  )
}


.dpprior_A2_migration_truth <- function(mode) {
  if (identical(mode, "a2_moment")) {
    return(list(
      controls = list(
        max_iter = 20L, damping = TRUE, use_fallback = TRUE,
        tol_step = 1e-10, log_bounds = c(-15, 15), boundary_tol = 1e-6,
        line_search_max = 20L, jacobian_rcond_singular = 1e-12,
        jacobian_rcond_ill = sqrt(.Machine$double.eps),
        fallback = list(
          maxit = 1000L, reltol = 1e-12,
          finite_penalty = .Machine$double.xmax / 1024
        ),
        selection_tolerance = 0
      ),
      tolerances = list(
        K_adequacy = list(
          absolute = 1e-8, relative = 1e-8,
          scale_formula = "max(abs(target),1)"
        ),
        K_stability = list(
          absolute = 1e-10, relative = 1e-8, scale_floor = 1
        ),
        step = 1e-10, boundary = 1e-6
      )
    ))
  }
  list(
    controls = list(
      max_iter = 100L, optimizer_tol = 1e-6,
      log_bounds = c(-15, 15), boundary_tol = 1e-6,
      fallback_max_iter = 100L,
      primary = list(
        maxit = 100L, factr = 1e-6 / .Machine$double.eps, pgtol = 1e-6
      ),
      fallback = list(
        iter.max = 100L, eval.max = 200L, rel.tol = 1e-6, x.tol = 1e-6
      ),
      fallback_trigger_worse_than_start = 1e-12,
      selection_tolerance = 0
    ),
    tolerances = list(
      distribution = list(
        adequacy = list(
          kl = 0.015, l1 = 0.11, mean_scaled = 0.01,
          variance_scaled = 0.065,
          mean_scale_formula = "max(1,sqrt(target_variance))",
          variance_scale_formula = "max(1,target_variance)"
        ),
        order = list(
          pmf_absolute = 1e-10, pmf_relative = 1e-8,
          pmf_l1 = 1e-10 + 1e-8,
          direct_moment_absolute = 1e-10,
          direct_moment_relative = 1e-8,
          target_identity_l1 = .TOL_PMF_SUM
        )
      ),
      boundary = 1e-6
    )
  )
}


.dpprior_A2_migration_computation <- function(method, controls) {
  setting <- list(
    method = method, controls = controls,
    parameterization = "Gamma(shape=a, rate=b)"
  )
  list(
    request = setting, used = setting,
    orders = list(
      M_requested = NULL, M_selected = NULL,
      M_verification_required = NULL, M_verification_used = NULL,
      requested_reason = "not_recorded_by_v1.1",
      selected_reason = "not_recorded_by_v1.1",
      verification_required_reason = "not_recorded_by_v1.1",
      verification_used_reason = "not_recorded_by_v1.1"
    ),
    scaling = list(
      requested = NULL, used = NULL,
      formula = "legacy_evidence_unavailable", values = list(),
      fixed_from_input = FALSE, change_reason = ""
    ),
    attempts = list(), candidate_evaluations = list(),
    selected_candidate_id = NULL, selected_attempt_id = NULL,
    fallback = list(
      attempted = FALSE, used = FALSE, trigger_attempt_id = NULL,
      selected_attempt_id = NULL, reason_code = NULL,
      message = "v1.1 fallback lineage was not reconstructed",
      outcome = "not_attempted"
    ),
    termination = list(
      code = "legacy_migration_no_selection",
      message = paste(
        "The v1.1 status and optimizer exit were not reused as candidate",
        "selection evidence."
      ),
      source = "constructor", iterations = NULL, boundary_reason = NULL
    ),
    trace = NULL, resources = list()
  )
}


.dpprior_A2_migration_verification <- function(target = FALSE) {
  reason <- if (target) {
    "v1.1 did not preserve canonical target verification evidence"
  } else {
    paste(
      "The v1.1 A2 parameter pair is retained only in compatibility",
      "evidence because optimizer candidate-selection lineage is absent."
    )
  }
  list(
    method = "legacy_evidence_quarantine", performed = FALSE,
    passed = FALSE, reason = reason, settings = list(),
    selected_snapshot = NULL, verifier_snapshot = NULL, stability = NULL,
    components = list(), invariants = list()
  )
}


.dpprior_A2_migration_provenance <- function(
    method, missing_evidence, backend_package_version) {
  list(
    requested_method = method, selected_method = method,
    is_fallback = FALSE,
    approximation = list(
      active = TRUE, opt_in = FALSE, kind = "legacy_schema_migration",
      warning_code = "legacy_object_upgraded"
    ),
    projection = list(
      applied = FALSE, opt_in = FALSE, policy = NULL, record = NULL
    ),
    parameterization = "Gamma(shape=a, rate=b)",
    backend = list(
      package = "DPprior", package_version = backend_package_version,
      implementation = "schema_upgrade_v1", source_commit = NULL
    ),
    input_fit = NULL,
    migration = list(
      source_schema = "DPprior/1.1/fit",
      adapter = "upgrade_DPprior_object", lossless = FALSE,
      missing_evidence = missing_evidence,
      warnings = "legacy_object_upgraded"
    ),
    legacy = list(active = FALSE, contract = NULL, deprecation_stage = NULL)
  )
}


.dpprior_A2_migration_expected_audit <- function(
    mode, J, source_parameters, target_K, M_verification) {
  M_selected <- as.integer(.QUAD_NODES_DEFAULT)
  a <- source_parameters[["a", exact = TRUE]]
  b <- source_parameters[["b", exact = TRUE]]
  candidate <- list(
    a = a, b = b, parameterization = "Gamma(shape=a, rate=b)"
  )
  if (identical(mode, "a2_kl")) {
    numerical <- tryCatch(
      withCallingHandlers({
        logS <- compute_log_stirling(J)
        induced <- .a2_kl_induced_log_pmf(
          J, a, b, logS, M = M_selected, M_verify = M_verification,
          abs_tol = 1e-10, rel_tol = 1e-8
        )
        target_pmf <- unname(target_K[["pmf", exact = TRUE]])
        list(
          induced = induced,
          selected = .a2_kl_pmf_metrics(target_pmf, induced$selected),
          verifier = .a2_kl_pmf_metrics(target_pmf, induced$verification)
        )
      }, warning = function(warning) stop(warning)),
      error = function(error) error
    )
    if (inherits(numerical, "condition")) return(numerical)
    selected <- numerical$selected
    verifier <- numerical$verifier
    induced <- numerical$induced
    return(list(
      method = "fresh_fixed_candidate_pmf_audit",
      candidate_unchanged = TRUE, candidate = candidate,
      M_selected = M_selected, M_verification = M_verification,
      selected = list(
        mean = selected$mean, variance = selected$variance,
        kl = selected$kl, l1 = selected$l1
      ),
      verifier = list(
        mean = verifier$mean, variance = verifier$variance,
        kl = verifier$kl, l1 = verifier$l1
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
  moments <- tryCatch(
    withCallingHandlers(
      exact_K_moments(
        J, a, b, M = M_selected, M_verify = M_verification,
        abs_tol = 1e-10, rel_tol = 1e-8, strict = FALSE
      ),
      warning = function(warning) stop(warning)
    ),
    error = function(error) error
  )
  if (inherits(moments, "condition")) return(moments)
  list(
    method = "fresh_fixed_candidate_moment_audit",
    candidate_unchanged = TRUE, candidate = candidate,
    M_selected = M_selected, M_verification = M_verification,
    selected = list(mean = moments$mean, variance = moments$var),
    order_stability = list(
      status = moments$status, reason = moments$quadrature$reason,
      mean_difference = moments$quadrature$mean_difference,
      variance_difference = moments$quadrature$variance_difference,
      mean_tolerance = moments$quadrature$mean_tolerance,
      variance_tolerance = moments$quadrature$variance_tolerance
    ),
    target_residuals = list(
      mean = moments$mean - target_K[["used", exact = TRUE]][[
        "mu_K", exact = TRUE
      ]],
      variance = moments$var - target_K[["used", exact = TRUE]][[
        "var_K", exact = TRUE
      ]]
    ),
    decision_ready = FALSE,
    decision_ready_reason = paste(
      "A fresh fixed-candidate moment audit cannot reconstruct v1.1",
      "optimizer selection lineage; refit required."
    )
  )
}


.dpprior_validate_A2_parameterless_contract <- function(raw) {
  mode <- raw[["mode", exact = TRUE]]
  if (!identical(raw[["object_type", exact = TRUE]], "fit") ||
      !mode %in% c("a2_moment", "a2_kl")) {
    return(invisible(TRUE))
  }

  computation <- raw[["computation", exact = TRUE]]
  verification <- raw[["verification", exact = TRUE]]
  provenance <- raw[["provenance", exact = TRUE]]
  migration <- provenance[["migration", exact = TRUE]]
  approximation <- provenance[["approximation", exact = TRUE]]
  backend <- provenance[["backend", exact = TRUE]]
  target_K <- raw[["target", exact = TRUE]][["K", exact = TRUE]]
  target_provenance <- target_K[["provenance", exact = TRUE]]
  compatibility_views <- raw[["compatibility", exact = TRUE]][[
    "views", exact = TRUE
  ]]
  source_view <- if ("source" %in% names(compatibility_views)) {
    compatibility_views[["source", exact = TRUE]]
  } else {
    NULL
  }
  token_is <- function(x, token) {
    is.character(x) && !is.object(x) && is.null(dim(x)) &&
      length(x) == 1L && !is.na(x) && identical(x, token)
  }
  contains_token <- function(x, token) {
    is.character(x) && !is.object(x) && is.null(dim(x)) &&
      !anyNA(x) && token %in% x
  }
  source_is_list <- typeof(source_view) == "list" && is.list(source_view) &&
    !is.object(source_view) && is.null(dim(source_view))
  target_migration <- target_provenance[["migration", exact = TRUE]]
  migration_marker <- any(c(
    token_is(migration[["source_schema", exact = TRUE]], "DPprior/1.1/fit"),
    token_is(migration[["adapter", exact = TRUE]], "upgrade_DPprior_object"),
    token_is(backend[["implementation", exact = TRUE]], "schema_upgrade_v1"),
    token_is(approximation[["kind", exact = TRUE]], "legacy_schema_migration"),
    token_is(approximation[["warning_code", exact = TRUE]],
             "legacy_object_upgraded"),
    contains_token(migration[["warnings", exact = TRUE]],
                   "legacy_object_upgraded"),
    contains_token(migration[["missing_evidence", exact = TRUE]],
                   "public_candidate_quarantined_pending_refit"),
    token_is(verification[["method", exact = TRUE]],
             "legacy_evidence_quarantine"),
    token_is(computation[["termination", exact = TRUE]][[
      "code", exact = TRUE
    ]], "legacy_migration_no_selection"),
    token_is(target_migration[["source_schema", exact = TRUE]],
             "DPprior/1.1/fit"),
    token_is(target_migration[["adapter", exact = TRUE]],
             "upgrade_DPprior_object"),
    source_is_list && isTRUE(source_view[[
      "public_candidate_quarantined", exact = TRUE
    ]])
  ))
  parameterless_approximate <-
    is.null(raw[["parameters", exact = TRUE]]) &&
      identical(raw[["status", exact = TRUE]], "approximate")

  if (migration_marker || parameterless_approximate) {
    expected_method <- if (identical(mode, "a2_moment")) "A2-MN" else "A2-KL"
    truth <- .dpprior_A2_migration_truth(mode)
    expected_message <- paste(
      "Migrated v1.1 candidate is approximate and unverified because optimizer",
      "selection and independent verification lineage were not serialized;",
      "refit is required for decision readiness."
    )
    no_extensions <- !any(c(
      "proxy", "constraint", "tradeoff", "legacy", "diagnostics",
      "sensitivity"
    ) %in% names(raw))
    .dpprior_schema_require(
      identical(raw[["object_type", exact = TRUE]], "fit") &&
        identical(raw[["method", exact = TRUE]], expected_method) &&
        identical(raw[["status", exact = TRUE]], "approximate") &&
        !raw[["usable", exact = TRUE]] && !raw[["verified", exact = TRUE]] &&
        is.null(raw[["parameters", exact = TRUE]]) &&
        length(raw[["achieved", exact = TRUE]]) == 0L &&
        length(raw[["residuals", exact = TRUE]]) == 0L && no_extensions &&
        identical(raw[["message", exact = TRUE]], expected_message) &&
        identical(raw[["tolerances", exact = TRUE]], truth$tolerances) &&
        identical(names(raw[["target", exact = TRUE]]), "K"),
      "migration_quarantine_status", "result",
      paste(
        "the exact parameterless, non-usable, unverified approximate A2",
        "migration quarantine, fixed message/tolerances, one K target, and",
        "no public scientific extension"
      ),
      list(
        mode = mode, method = raw[["method", exact = TRUE]],
        status = raw[["status", exact = TRUE]],
        usable = raw[["usable", exact = TRUE]],
        verified = raw[["verified", exact = TRUE]],
        parameters = raw[["parameters", exact = TRUE]],
        achieved = raw[["achieved", exact = TRUE]],
        residuals = raw[["residuals", exact = TRUE]]
      )
    )

    expected_computation <- .dpprior_A2_migration_computation(
      expected_method, truth$controls
    )
    .dpprior_schema_require(
      identical(computation, expected_computation),
      "migration_quarantine_computation", "result.computation",
      paste(
        "the exact migrated method/control setting, unknown orders, empty",
        "candidate/attempt ledgers, fixed quarantine scaling/fallback, and",
        "legacy_migration_no_selection/constructor termination"
      ), computation
    )

    expected_verification <- .dpprior_A2_migration_verification()
    .dpprior_schema_require(
      identical(verification, expected_verification),
      "migration_quarantine_verification", "result.verification",
      paste(
        "the exact unperformed/unpassed legacy evidence quarantine with no",
        "decision snapshots, settings, components, stability, or invariants"
      ), verification
    )

    expected_missing <- c(
      "requested_controls", "M_requested", "M_selected",
      "M_verification_required", "M_verification_used", "scaling",
      "canonical_attempts", "candidate_selection", "fallback_lineage",
      "independent_verifier_snapshot", "source_commit", "optimizer_lineage"
    )
    if (identical(mode, "a2_kl")) {
      expected_missing <- c(expected_missing, "target_family_verification")
    }
    expected_missing <- c(
      expected_missing, "public_candidate_quarantined_pending_refit"
    )
    destination_package_version <- backend[["package_version", exact = TRUE]]
    .dpprior_schema_require(
      identical(backend[["package", exact = TRUE]], "DPprior") &&
        .dpprior_schema_is_retained_v2_package_version(
          destination_package_version
        ) && identical(
          backend[["implementation", exact = TRUE]], "schema_upgrade_v1"
        ) && is.null(backend[["source_commit", exact = TRUE]]),
      "migration_quarantine_provenance", "result.provenance.backend",
      paste(
        "a retained plain v2.x destination package version with the exact",
        "DPprior/schema_upgrade_v1 backend identity"
      ),
      backend
    )
    expected_provenance <- .dpprior_A2_migration_provenance(
      expected_method, expected_missing, destination_package_version
    )
    .dpprior_schema_require(
      identical(provenance, expected_provenance),
      "migration_quarantine_provenance", "result.provenance",
      paste(
        "the exact v1.1 schema-upgrade provenance, parameterization/backend,",
        "complete ordered missing-evidence inventory, and inactive legacy",
        "semantics"
      ), provenance
    )

    source_fields <- c(
      "source_schema", "source_package_version", "source_class",
      "source_fields", "source_status", "source_converged",
      "source_parameters", "source_digest", "digest_method",
      "discarded_legacy_fields", "public_candidate_quarantined",
      "required_action", "authority", "lossy", "consumer_policy"
    )
    source_shape <- source_is_list &&
      identical(names(source_view), source_fields)
    .dpprior_schema_require(
      source_shape, "migration_quarantine_compatibility",
      "result.compatibility.views.source",
      "the exact retained v1.1 source-summary field vocabulary", source_view
    )
    source_parameters <- source_view[["source_parameters", exact = TRUE]]
    source_parameter_shape <- typeof(source_parameters) == "list" &&
      is.list(source_parameters) && !is.object(source_parameters) &&
      is.null(dim(source_parameters)) &&
      identical(names(source_parameters), c("a", "b", "J", "parameterization"))
    .dpprior_schema_require(
      source_parameter_shape &&
        .dpprior_schema_is_finite_scalar(source_parameters[["a", exact = TRUE]]) &&
        source_parameters[["a", exact = TRUE]] > 0 &&
        .dpprior_schema_is_finite_scalar(source_parameters[["b", exact = TRUE]]) &&
        source_parameters[["b", exact = TRUE]] > 0 &&
        .dpprior_schema_is_count(source_parameters[["J", exact = TRUE]], 1L) &&
        source_parameters[["J", exact = TRUE]] == raw[["J", exact = TRUE]] &&
        identical(source_parameters[["parameterization", exact = TRUE]],
                  "Gamma(shape=a, rate=b)") &&
        identical(source_view[["source_schema", exact = TRUE]],
                  "DPprior/1.1/fit") &&
        identical(source_view[["source_package_version", exact = TRUE]],
                  "1.1.0") &&
        identical(source_view[["source_class", exact = TRUE]], "DPprior_fit") &&
        is.character(source_view[["source_fields", exact = TRUE]]) &&
        length(source_view[["source_fields", exact = TRUE]]) > 0L &&
        !anyNA(source_view[["source_fields", exact = TRUE]]) &&
        all(nzchar(source_view[["source_fields", exact = TRUE]])) &&
        .dpprior_schema_is_scalar_character(
          source_view[["source_status", exact = TRUE]], allow_empty = TRUE
        ) &&
        .dpprior_schema_is_scalar_logical(
          source_view[["source_converged", exact = TRUE]]
        ) &&
        is.character(source_view[["source_digest", exact = TRUE]]) &&
        length(source_view[["source_digest", exact = TRUE]]) == 1L &&
        grepl("^[0-9a-f]{32}$", source_view[["source_digest", exact = TRUE]]) &&
        identical(source_view[["digest_method", exact = TRUE]],
                  "base-r-saveRDS-v3-md5") &&
        isTRUE(source_view[["public_candidate_quarantined", exact = TRUE]]) &&
        identical(source_view[["required_action", exact = TRUE]],
                  "refit_with_current_API") &&
        identical(
          source_view[c("authority", "lossy", "consumer_policy")],
          .dpprior_compatibility_quarantine_boundary()
        ) &&
        identical(raw[["compatibility", exact = TRUE]][[
          "deprecations", exact = TRUE
        ]], list(legacy_schema = c(list(
          code = "legacy_object_upgraded",
          first_deprecated_version = "2.0.0",
          removal_floor = "not_scheduled"
        ), .dpprior_compatibility_quarantine_boundary()))),
      "migration_quarantine_compatibility",
      "result.compatibility.views.source",
      paste(
        "typed retained legacy source parameters/digest plus explicit public",
        "candidate quarantine, refit action, and deprecation record"
      ), source_view
    )

    direct_source_fields <- c(
      "a", "b", "J", "target", "method", "status", "converged",
      "iterations", "termination", "fit", "diagnostics", "trace"
    )
    wrapper_source_fields <- c(
      "a", "b", "J", "target", "method", "status", "converged",
      "iterations", "termination", "fit", "solver_diagnostics", "trace"
    )
    wrapper_diagnostic_fields <- c(wrapper_source_fields, "diagnostics")
    source_route <- if (identical(
      source_view[["source_fields", exact = TRUE]], direct_source_fields
    )) {
      "direct"
    } else if (identical(mode, "a2_moment") && identical(
      source_view[["source_fields", exact = TRUE]], wrapper_source_fields
    )) {
      "wrapper"
    } else if (identical(mode, "a2_moment") && identical(
      source_view[["source_fields", exact = TRUE]], wrapper_diagnostic_fields
    )) {
      "wrapper_diagnostics"
    } else {
      NULL
    }
    discarded <- source_view[["discarded_legacy_fields", exact = TRUE]]
    allowed_discarded <- intersect(
      c("diagnostics", "solver_diagnostics"),
      source_view[["source_fields", exact = TRUE]]
    )
    discarded_ok <- is.null(discarded) || (
      is.character(discarded) && !is.object(discarded) &&
        is.null(attributes(discarded)) && length(discarded) > 0L &&
        !anyNA(discarded) && !anyDuplicated(discarded) &&
        identical(discarded, allowed_discarded[allowed_discarded %in% discarded])
    )
    .dpprior_schema_require(
      !is.null(source_route) && discarded_ok,
      "migration_quarantine_compatibility",
      "result.compatibility.views.source",
      paste(
        "one exact frozen direct/wrapper A2 source field route and an ordered",
        "subset of genuinely discardable diagnostic fields"
      ), source_view[c("source_fields", "discarded_legacy_fields")]
    )

    target_compatibility <- target_K[["compatibility", exact = TRUE]]
    target_views <- target_compatibility[["views", exact = TRUE]]
    legacy_target <- target_views[["legacy_target", exact = TRUE]]
    target_method <- paste0("legacy_target:", expected_method)
    target_J <- as.integer(raw[["J", exact = TRUE]])
    legacy_target_ok <- typeof(legacy_target) == "list" &&
      is.list(legacy_target) && !is.object(legacy_target) &&
      is.null(dim(legacy_target))
    .dpprior_schema_require(
      legacy_target_ok,
      "migration_target_compatibility",
      "result.target.K.compatibility.views.legacy_target",
      "one ordinary retained frozen v1.1 target record", legacy_target
    )

    if (identical(mode, "a2_moment")) {
      legacy_fields <- if (identical(source_route, "direct")) {
        c("mu_K", "var_K", "type")
      } else {
        c("mu_K", "var_K", "var_K_used", "confidence", "type")
      }
      scalar_target <- function(field) {
        .dpprior_schema_is_finite_scalar(legacy_target[[field, exact = TRUE]]) &&
          legacy_target[[field, exact = TRUE]] > 0
      }
      wrapper_target_ok <- identical(source_route, "direct") || (
        scalar_target("var_K_used") &&
          identical(
            as.numeric(legacy_target[["var_K_used", exact = TRUE]]),
            as.numeric(legacy_target[["var_K", exact = TRUE]])
          ) &&
          (is.null(legacy_target[["confidence", exact = TRUE]]) ||
             token_is(legacy_target[["confidence", exact = TRUE]], "low") ||
             token_is(legacy_target[["confidence", exact = TRUE]], "medium") ||
             token_is(legacy_target[["confidence", exact = TRUE]], "high"))
      )
      .dpprior_schema_require(
        identical(names(legacy_target), legacy_fields) &&
          scalar_target("mu_K") && scalar_target("var_K") &&
          identical(legacy_target[["type", exact = TRUE]], "moments") &&
          wrapper_target_ok,
        "migration_target_compatibility",
        "result.target.K.compatibility.views.legacy_target",
        paste(
          "the exact direct or wrapper moment-target vocabulary, positive",
          "moments, unchanged used variance, and closed confidence label"
        ), legacy_target
      )
      request <- list(
        J = target_J,
        mu_K = as.numeric(legacy_target[["mu_K", exact = TRUE]]),
        var_K = as.numeric(legacy_target[["var_K", exact = TRUE]])
      )
      normalized <- c(request, list(interval = NULL, pmf = NULL))
      target_kind <- "moments"
      target_pmf <- NULL
      target_implied <- list(mean = request$mu_K, variance = request$var_K)
      target_rule <- "canonicalize_direct_moments"
    } else {
      legacy_fields <- c(
        "type", "pmf", "mu_K", "var_K", "df", "scale",
        "mu_K_discrete", "var_K_discrete"
      )
      legacy_pmf <- legacy_target[["pmf", exact = TRUE]]
      legacy_pmf_ok <- is.numeric(legacy_pmf) && !is.object(legacy_pmf) &&
        is.null(attributes(legacy_pmf)) && length(legacy_pmf) == target_J &&
        !anyNA(legacy_pmf) && all(is.finite(legacy_pmf)) &&
        all(legacy_pmf >= 0) && abs(sum(legacy_pmf) - 1) <= .TOL_PMF_SUM
      positive_fields <- c("mu_K", "var_K", "df", "scale", "mu_K_discrete")
      positive_ok <- all(vapply(
        positive_fields,
        function(field) {
          .dpprior_schema_is_finite_scalar(
            legacy_target[[field, exact = TRUE]]
          ) && legacy_target[[field, exact = TRUE]] > 0
        },
        logical(1)
      ))
      variance_discrete <- legacy_target[["var_K_discrete", exact = TRUE]]
      .dpprior_schema_require(
        identical(names(legacy_target), legacy_fields) &&
          identical(legacy_target[["type", exact = TRUE]], "chisq") &&
          legacy_pmf_ok && positive_ok &&
          .dpprior_schema_is_finite_scalar(variance_discrete) &&
          variance_discrete >= 0,
        "migration_target_compatibility",
        "result.target.K.compatibility.views.legacy_target",
        paste(
          "the exact frozen chi-square target vocabulary with a strict PMF",
          "and finite positive retained source parameters"
        ), legacy_target
      )
      request <- list(J = target_J, pmf = legacy_pmf)
      normalized <- c(request, list(interval = NULL))
      target_kind <- "pmf"
      target_pmf <- legacy_pmf
      target_implied <- list(
        mean = legacy_target[["mu_K_discrete", exact = TRUE]],
        variance = variance_discrete
      )
      target_rule <- "validate_strict_pmf"
    }

    target_evidence <- list(
      source = "recognized_v1.1_shape",
      source_schema = "DPprior/1.1/fit",
      source_fields = legacy_fields
    )
    expected_target <- list(
      schema = list(name = "dpprior.target", version = 1L),
      kind = target_kind, J = target_J, support = seq_len(target_J),
      request = request, normalized = normalized, used = normalized,
      derivation = list(
        request_to_normalized = list(
          rule = target_rule, outcome = "canonicalized", opt_in = FALSE,
          before = request, after = normalized, evidence = target_evidence
        ),
        normalized_to_used = NULL
      ),
      interval = NULL, family = NULL,
      assumptions = list(source = "recognized_v1.1_shape"),
      pmf = target_pmf, implied = target_implied,
      achieved_interval = NULL, residuals = list(), tolerances = list(),
      status = "approximate", usable = FALSE, verified = FALSE,
      message = paste(
        "Legacy target structure was preserved but lacks v2 verification",
        "lineage."
      ),
      parameters = NULL,
      computation = .dpprior_A2_migration_computation(target_method, list()),
      verification = .dpprior_A2_migration_verification(target = TRUE),
      provenance = .dpprior_A2_migration_provenance(
        target_method,
        c("target_verification_lineage", "explicit_request_nulls"),
        destination_package_version
      ),
      compatibility = list(
        top_level_aliases = character(),
        views = list(
          legacy_target = legacy_target,
          legacy_target_fields = legacy_fields
        ),
        deprecations = list(legacy_target_view = c(list(
          code = "legacy_target_view_quarantined",
          first_deprecated_version = "2.0.0",
          removal_floor = "not_scheduled"
        ), .dpprior_compatibility_quarantine_boundary()))
      )
    )
    .dpprior_schema_require(
      identical(target_K[["J", exact = TRUE]], target_J) &&
        source_parameters[["J", exact = TRUE]] == target_J &&
        identical(unclass(target_K), expected_target),
      "migration_target_contract", "result.target.K",
      paste(
        "the complete canonical migrated target derived from the retained",
        "frozen v1.1 target, including exact computation, verification,",
        "provenance, compatibility, message, and scientific authority"
      ), target_K
    )

    compatibility <- raw[["compatibility", exact = TRUE]]
    expected_deprecations <- list(legacy_schema = c(list(
      code = "legacy_object_upgraded",
      first_deprecated_version = "2.0.0",
      removal_floor = "not_scheduled"
    ), .dpprior_compatibility_quarantine_boundary()))
    audit_present <- "fixed_candidate_recomputation" %in%
      names(compatibility_views)
    expected_view_names <- if (audit_present) {
      c("source", "fixed_candidate_recomputation")
    } else {
      "source"
    }
    .dpprior_schema_require(
      identical(compatibility[["top_level_aliases", exact = TRUE]],
                character()) &&
        identical(names(compatibility_views), expected_view_names) &&
        identical(compatibility[["deprecations", exact = TRUE]],
                  expected_deprecations),
      "migration_quarantine_compatibility", "result.compatibility",
      paste(
        "no aliases, exactly source plus an optional fixed-candidate audit,",
        "no extra views, and the one canonical legacy deprecation"
      ), compatibility
    )

    if (audit_present) {
      audit <- compatibility_views[[
        "fixed_candidate_recomputation", exact = TRUE
      ]]
      expected_audit_fields <- if (identical(mode, "a2_kl")) {
        c(
          "method", "candidate_unchanged", "candidate", "M_selected",
          "M_verification", "selected", "verifier", "numerical_status",
          "numerical_reason", "target_residuals", "decision_ready",
          "decision_ready_reason"
        )
      } else {
        c(
          "method", "candidate_unchanged", "candidate", "M_selected",
          "M_verification", "selected", "order_stability",
          "target_residuals", "decision_ready", "decision_ready_reason"
        )
      }
      .dpprior_schema_require(
        typeof(audit) == "list" && is.list(audit) && !is.object(audit) &&
          is.null(dim(audit)) &&
          identical(names(audit), expected_audit_fields),
        "migration_fixed_candidate_audit",
        "result.compatibility.views.fixed_candidate_recomputation",
        "the exact ordered mode-specific fixed-candidate audit vocabulary",
        audit
      )
      M_selected <- audit[["M_selected", exact = TRUE]]
      M_verification <- audit[["M_verification", exact = TRUE]]
      required_order <- as.integer(
        .quadrature_verification_required_order(.QUAD_NODES_DEFAULT)
      )
      order_ok <- identical(M_selected, as.integer(.QUAD_NODES_DEFAULT)) &&
        is.integer(M_verification) && !is.object(M_verification) &&
        is.null(attributes(M_verification)) && length(M_verification) == 1L &&
        !is.na(M_verification) && M_verification >= required_order &&
        M_verification <= .QUADRATURE_MAX_NODES
      .dpprior_schema_require(
        order_ok,
        "migration_fixed_candidate_audit",
        "result.compatibility.views.fixed_candidate_recomputation",
        paste(
          "M_selected=80 and an integer verifier order from the canonical",
          "required order through the quadrature ceiling"
        ), audit[c("M_selected", "M_verification")]
      )
      expected_audit <- .dpprior_A2_migration_expected_audit(
        mode, target_J, source_parameters, target_K, M_verification
      )
      .dpprior_schema_require(
        !inherits(expected_audit, "condition") &&
          identical(audit, expected_audit),
        "migration_fixed_candidate_audit",
        "result.compatibility.views.fixed_candidate_recomputation",
        paste(
          "the exact fresh fixed-candidate audit recomputed from retained",
          "source a/b/J, canonical target, M=80, and the recorded verifier",
          "order, with decision_ready=FALSE"
        ),
        if (inherits(expected_audit, "condition")) {
          conditionMessage(expected_audit)
        } else {
          list(recorded = audit, expected = expected_audit)
        }
      )
    }
    # R24 retains no independent verify-request marker. Consequently, absence
    # of the entire optional audit is the canonical verify=FALSE shape and is
    # intentionally accepted as a conservative evidence downgrade. Any
    # present audit is closed and freshly recomputed above.
    return(invisible(TRUE))
  }

  if (is.null(raw[["parameters", exact = TRUE]]) &&
      identical(raw[["status", exact = TRUE]], "failed")) {
    attempts <- computation[["attempts", exact = TRUE]]
    retained_execution_failure <- length(attempts) > 0L && any(vapply(
      attempts,
      function(attempt) {
        !is.null(attempt[["error", exact = TRUE]]) ||
          (!is.null(attempt[["exit_code", exact = TRUE]]) &&
             !identical(attempt[["exit_code", exact = TRUE]], 0L) &&
             !identical(attempt[["exit_code", exact = TRUE]], 0))
      },
      logical(1)
    ))
    .dpprior_schema_require(
      retained_execution_failure,
      "parameterless_A2_failure_evidence", "result.computation.attempts",
      paste(
        "at least one retained typed attempt error or non-zero optimizer exit",
        "for a parameterless failed A2 fit"
      ), attempts
    )
  }
  invisible(TRUE)
}


.dpprior_validate_result_verification_contract <- function(raw) {
  mode <- raw[["mode", exact = TRUE]]
  status <- raw[["status", exact = TRUE]]
  verification <- raw[["verification", exact = TRUE]]
  computation <- raw[["computation", exact = TRUE]]
  orders <- computation[["orders", exact = TRUE]]
  controls <- computation[["used", exact = TRUE]][["controls", exact = TRUE]]
  components <- verification[["components", exact = TRUE]]
  invariants <- verification[["invariants", exact = TRUE]]
  finite_fit <- identical(raw[["object_type", exact = TRUE]], "fit") &&
    !is.null(raw[["parameters", exact = TRUE]])

  validate_checks <- function(records, expected_names, source, path) {
    .dpprior_schema_require(
      identical(names(records), expected_names),
      "verification_vocabulary", path,
      paste("the exact ordered check vocabulary:",
            paste(expected_names, collapse = ", ")),
      names(records)
    )
    for (name in expected_names) {
      check <- records[[name, exact = TRUE]]
      .dpprior_validate_decision_check(check, paste0(path, ".", name))
      .dpprior_schema_require(
        identical(check[["source", exact = TRUE]], source),
        "verification_source", paste0(path, ".", name, ".source"),
        source, check[["source", exact = TRUE]]
      )
    }
  }

  if (identical(mode, "prior_diagnostics")) {
    diagnostics <- raw[["diagnostics", exact = TRUE]]
    expected_evidence <- .dpprior_expected_diagnostics_evidence(raw)
    diagnostic_target <- raw[["target", exact = TRUE]]
    .dpprior_schema_exact_names(
      diagnostic_target, c("requested_components", "warning_policy"),
      "result.target"
    )
    .dpprior_schema_require(
      identical(
        diagnostic_target[["requested_components", exact = TRUE]],
        .DPPRIOR_DIAGNOSTIC_COMPONENTS
      ),
      "diagnostic_target", "result.target.requested_components",
      "the four ordered canonical diagnostic components",
      diagnostic_target[["requested_components", exact = TRUE]]
    )
    warning_policy <- diagnostic_target[["warning_policy", exact = TRUE]]
    if (is.null(warning_policy)) {
      .dpprior_schema_require(
        length(diagnostics[["policy_results", exact = TRUE]]) == 0L &&
          length(diagnostics[["warnings", exact = TRUE]]) == 0L,
        "diagnostic_policy_authority", "result.diagnostics.policy_results",
        "no policy result or warning when warning_policy is NULL",
        diagnostics[c("policy_results", "warnings")]
      )
    } else {
      .dpprior_schema_exact_names(
        warning_policy,
        c("estimand", "direction", "weight_threshold", "action_threshold"),
        "result.target.warning_policy"
      )
      for (field in c("estimand", "direction")) {
        .dpprior_schema_validate_scalar_character(
          warning_policy[[field, exact = TRUE]],
          paste0("result.target.warning_policy.", field)
        )
      }
      .dpprior_schema_require(
        warning_policy[["estimand", exact = TRUE]] %in% c("W_SB", "W_max") &&
          warning_policy[["direction", exact = TRUE]] %in% c("above", "below"),
        "diagnostic_policy_authority", "result.target.warning_policy",
        "the closed W_SB/W_max and above/below policy vocabulary",
        warning_policy
      )
      .dpprior_schema_validate_finite_scalar(
        warning_policy[["weight_threshold", exact = TRUE]],
        "result.target.warning_policy.weight_threshold",
        lower = 0, upper = 1, lower_open = TRUE, upper_open = TRUE
      )
      .dpprior_schema_validate_finite_scalar(
        warning_policy[["action_threshold", exact = TRUE]],
        "result.target.warning_policy.action_threshold", lower = 0, upper = 1
      )
      .dpprior_schema_require(
        length(diagnostics[["policy_results", exact = TRUE]]) == 1L,
        "diagnostic_policy_authority", "result.diagnostics.policy_results",
        "exactly one result for the one retained warning policy",
        diagnostics[["policy_results", exact = TRUE]]
      )
      policy_result <- diagnostics[["policy_results", exact = TRUE]][[1L]]
      .dpprior_schema_require(
        identical(policy_result[["estimand", exact = TRUE]],
                  warning_policy[["estimand", exact = TRUE]]) &&
          identical(policy_result[["direction", exact = TRUE]],
                    warning_policy[["direction", exact = TRUE]]) &&
          identical(policy_result[["threshold", exact = TRUE]],
                    warning_policy[["action_threshold", exact = TRUE]]),
        "diagnostic_policy_authority", "result.diagnostics.policy_results[[1]]",
        paste(
          "estimand/direction identity and the canonical action threshold",
          "from result.target.warning_policy"
        ),
        policy_result
      )
      if (identical(warning_policy[["estimand", exact = TRUE]], "W_SB")) {
        expected_policy_value <- as.numeric(.diagnostic_wsb_tail(
          warning_policy[["weight_threshold", exact = TRUE]],
          raw[["parameters", exact = TRUE]][["a", exact = TRUE]],
          raw[["parameters", exact = TRUE]][["b", exact = TRUE]]
        ))
        .dpprior_schema_require(
          identical(policy_result[["basis", exact = TRUE]],
                    "exact_tail_probability") &&
            !is.null(policy_result[["value", exact = TRUE]]) &&
            is.null(policy_result[["lower", exact = TRUE]]) &&
            is.null(policy_result[["upper", exact = TRUE]]) &&
            abs(policy_result[["value", exact = TRUE]] -
                  expected_policy_value) <=
              64 * .Machine$double.eps * max(
                1, abs(policy_result[["value", exact = TRUE]]),
                abs(expected_policy_value)
              ),
          "diagnostic_policy_truth", "result.diagnostics.policy_results[[1]]",
          paste(
            "the exact W_SB tail probability freshly recomputed from",
            "parameters and the retained weight threshold"
          ),
          policy_result
        )
      } else {
        .dpprior_schema_require(
          identical(policy_result[["basis", exact = TRUE]],
                    "backend_unavailable") &&
            all(vapply(
              policy_result[c("value", "lower", "upper")],
              is.null, logical(1)
            )),
          "diagnostic_policy_truth", "result.diagnostics.policy_results[[1]]",
          paste(
            "backend_unavailable for W_max unless independently retained",
            "maximum-weight evidence is added to the canonical extension"
          ),
          policy_result
        )
      }
    }
    expected_component_status <- ifelse(
      expected_evidence[["component_pass", exact = TRUE]],
      "converged", "approximate"
    )
    expected_component_records <- list(
      alpha = c(
        list(
          status = expected_component_status[["alpha"]],
          usable = expected_evidence[["component_pass", exact = TRUE]][["alpha"]],
          verified = expected_evidence[["component_pass", exact = TRUE]][["alpha"]]
        ),
        expected_evidence[["selected", exact = TRUE]][["alpha", exact = TRUE]]
      ),
      K = c(
        list(
          status = expected_component_status[["K"]],
          usable = expected_evidence[["component_pass", exact = TRUE]][["K"]],
          verified = expected_evidence[["component_pass", exact = TRUE]][["K"]]
        ),
        expected_evidence[["selected", exact = TRUE]][["K", exact = TRUE]][
          c("mean", "variance", "pmf", "M")
        ]
      ),
      weights = c(
        list(
          status = expected_component_status[["weights"]],
          usable = expected_evidence[["component_pass", exact = TRUE]][["weights"]],
          verified = expected_evidence[["component_pass", exact = TRUE]][["weights"]]
        ),
        expected_evidence[["selected", exact = TRUE]][["weights", exact = TRUE]]
      ),
      coclustering = c(
        list(
          status = expected_component_status[["coclustering"]],
          usable = expected_evidence[["component_pass", exact = TRUE]][["coclustering"]],
          verified = expected_evidence[["component_pass", exact = TRUE]][["coclustering"]]
        ),
        expected_evidence[["selected", exact = TRUE]][[
          "coclustering", exact = TRUE
        ]]
      )
    )
    .dpprior_schema_require(
      identical(
        diagnostics[.DPPRIOR_DIAGNOSTIC_COMPONENTS],
        expected_component_records
      ),
      "diagnostic_component_truth", "result.diagnostics",
      paste(
        "all selected-order diagnostic values and component status flags",
        "recomputed from J, parameters, orders, and fixed tolerances"
      ), diagnostics[.DPPRIOR_DIAGNOSTIC_COMPONENTS]
    )
    component_status <- vapply(
      .DPPRIOR_DIAGNOSTIC_COMPONENTS,
      function(name) {
        component <- diagnostics[[name, exact = TRUE]]
        if (is.null(component)) "failed" else
          component[["status", exact = TRUE]]
      }, character(1)
    )
    precedence <- c(
      converged = 1L, boundary = 2L, approximate = 3L,
      infeasible = 4L, failed = 5L
    )
    expected_status <- component_status[[
      which.max(unname(precedence[component_status]))
    ]]
    component_verified <- setNames(vapply(
      .DPPRIOR_DIAGNOSTIC_COMPONENTS,
      function(name) {
        component <- diagnostics[[name, exact = TRUE]]
        !is.null(component) && isTRUE(component[["verified", exact = TRUE]])
      }, logical(1)
    ), .DPPRIOR_DIAGNOSTIC_COMPONENTS)
    expected_verified <- expected_status %in% c("converged", "boundary") &&
      all(component_verified)
    expected_usable <- expected_status %in% c("converged", "boundary")
    .dpprior_schema_require(
      identical(status, expected_status) &&
        identical(raw[["verified", exact = TRUE]], expected_verified) &&
        identical(raw[["usable", exact = TRUE]], expected_usable),
      "diagnostic_status_aggregation", "result.status",
      paste(
        "component precedence failed>infeasible>approximate>boundary>",
        "converged with recomputed usable/verified flags"
      ),
      list(
        recorded = raw[c("status", "usable", "verified")],
        component_status = component_status
      )
    )
    attempts <- computation[["attempts", exact = TRUE]]
    .dpprior_schema_require(
      length(attempts) == length(.DPPRIOR_DIAGNOSTIC_COMPONENTS) &&
        identical(
          unname(vapply(
            attempts, function(attempt) attempt[["id", exact = TRUE]],
            character(1)
          )),
          paste0("diagnostic-", .DPPRIOR_DIAGNOSTIC_COMPONENTS)
        ) && identical(
          unname(vapply(
            attempts, function(attempt) attempt[["method", exact = TRUE]],
            character(1)
          )),
          unname(.DPPRIOR_DIAGNOSTIC_ATTEMPT_METHODS)
        ) && all(vapply(attempts, function(attempt) {
          identical(attempt[["stage", exact = TRUE]], "diagnostic_component") &&
            !attempt[["selected", exact = TRUE]]
        }, logical(1))) &&
        is.null(computation[["selected_attempt_id", exact = TRUE]]) &&
        is.null(computation[["selected_candidate_id", exact = TRUE]]) &&
        length(computation[["candidate_evaluations", exact = TRUE]]) == 0L,
      "diagnostic_attempts", "result.computation.attempts",
      "the four ordered component-specific diagnostic attempts and no selection",
      attempts
    )
    for (index in seq_along(attempts)) {
      expected_reason <- paste0("component_", component_status[[index]])
      expected_unavailable <- c(
        start = "component diagnostic has no optimizer start",
        bounds = "component diagnostic has no optimizer bounds",
        candidate_parameters = "diagnostics do not select fit parameters",
        candidate_objective = "diagnostics do not optimize an objective"
      )
      .dpprior_schema_require(
        identical(attempts[[index]][["reason_code", exact = TRUE]],
                  expected_reason) &&
          identical(attempts[[index]][["exit_code", exact = TRUE]], 0L) &&
          identical(attempts[[index]][["iterations", exact = TRUE]], 0L) &&
          identical(
            attempts[[index]][["control", exact = TRUE]],
            list(component = .DPPRIOR_DIAGNOSTIC_COMPONENTS[[index]])
          ) && identical(
            attempts[[index]][["evaluations", exact = TRUE]],
            list(function_count = 1L)
          ) && is.null(attempts[[index]][["start", exact = TRUE]]) &&
          is.null(attempts[[index]][["bounds", exact = TRUE]]) &&
          length(attempts[[index]][["warnings", exact = TRUE]]) == 0L &&
          is.null(attempts[[index]][["error", exact = TRUE]]) &&
          is.null(attempts[[index]][["candidate_parameters", exact = TRUE]]) &&
          is.null(attempts[[index]][["candidate_objective", exact = TRUE]]) &&
          identical(attempts[[index]][["unavailable", exact = TRUE]],
                    expected_unavailable),
        "diagnostic_attempt_status",
        sprintf("result.computation.attempts[[%d]]", index),
        paste(
          "the exact normalized component execution ledger with no optimizer",
          "candidate claims"
        ), attempts[[index]]
      )
    }
    .dpprior_schema_require(
      identical(
        computation[["termination", exact = TRUE]][["code", exact = TRUE]],
        "diagnostics_recomputed"
      ) && identical(
        computation[["termination", exact = TRUE]][["source", exact = TRUE]],
        "component_aggregation"
      ) && is.null(computation[["termination", exact = TRUE]][[
        "iterations", exact = TRUE
      ]]),
      "diagnostic_termination", "result.computation.termination",
      "diagnostics_recomputed/component_aggregation with no selected iterations",
      computation[["termination", exact = TRUE]]
    )
    expected_settings <- list(
      M_selected = orders[["M_selected", exact = TRUE]],
      M_verification = orders[["M_verification_used", exact = TRUE]]
    )
    .dpprior_schema_require(
      identical(verification[["method", exact = TRUE]],
                "fresh_component_specific_diagnostics") &&
        verification[["performed", exact = TRUE]] &&
        identical(verification[["passed", exact = TRUE]], expected_verified) &&
        identical(verification[["settings", exact = TRUE]], expected_settings) &&
        is.null(verification[["stability", exact = TRUE]]) &&
        !is.null(verification[["selected_snapshot", exact = TRUE]]) &&
        !is.null(verification[["verifier_snapshot", exact = TRUE]]) &&
        identical(
          verification[["selected_snapshot", exact = TRUE]][["source", exact = TRUE]],
          "fresh_diagnostics_selected_order"
        ) && identical(
          verification[["verifier_snapshot", exact = TRUE]][["source", exact = TRUE]],
          "fresh_diagnostics_verifier_evidence"
        ),
      "diagnostic_verification_contract", "result.verification",
      "the exact fresh component-specific diagnostic verifier contract",
      verification
    )
    validate_checks(
      components, "component_aggregation", "fresh_component_specific_checks",
      "result.verification.components"
    )
    validate_checks(
      invariants, c("fixed_parameters", "dominance_category_removed"),
      "fresh_component_specific_checks", "result.verification.invariants"
    )
    .dpprior_bind_decision_check(
      components[["component_aggregation", exact = TRUE]],
      component_verified, setNames(rep(TRUE, 4L), names(component_verified)),
      NULL, "identical", "result.verification.components.component_aggregation"
    )
    fixed_parameters <- identical(
      verification[["selected_snapshot", exact = TRUE]][[
        "parameters", exact = TRUE
      ]], raw[["parameters", exact = TRUE]]
    ) && identical(
      verification[["verifier_snapshot", exact = TRUE]][[
        "parameters", exact = TRUE
      ]], raw[["parameters", exact = TRUE]]
    )
    dominance_removed <- length(.dpprior_schema_find_forbidden_names(
      diagnostics, "dominance_risk", "result.diagnostics"
    )) == 0L
    .dpprior_bind_decision_check(
      invariants[["fixed_parameters", exact = TRUE]], fixed_parameters, TRUE,
      NULL, "identical", "result.verification.invariants.fixed_parameters"
    )
    .dpprior_bind_decision_check(
      invariants[["dominance_category_removed", exact = TRUE]],
      dominance_removed, TRUE, NULL, "identical",
      "result.verification.invariants.dominance_category_removed"
    )
    selected_evidence <- verification[["selected_snapshot", exact = TRUE]]
    verifier_evidence <- verification[["verifier_snapshot", exact = TRUE]]
    expected_verifier_achieved <- raw[["achieved", exact = TRUE]]
    expected_verifier_achieved <- expected_evidence[["verifier", exact = TRUE]]
    .dpprior_schema_require(
      identical(raw[["achieved", exact = TRUE]],
                expected_evidence[["selected", exact = TRUE]]) &&
        identical(selected_evidence[["achieved", exact = TRUE]],
                  expected_evidence[["selected", exact = TRUE]]) &&
        identical(selected_evidence[["residuals", exact = TRUE]],
                  raw[["residuals", exact = TRUE]]) &&
        identical(selected_evidence[["tolerances", exact = TRUE]],
                  raw[["tolerances", exact = TRUE]]) &&
        identical(selected_evidence[["M", exact = TRUE]],
                  orders[["M_selected", exact = TRUE]]) &&
        identical(verifier_evidence[["achieved", exact = TRUE]],
                  expected_verifier_achieved) &&
        identical(verifier_evidence[["residuals", exact = TRUE]],
                  raw[["residuals", exact = TRUE]]) &&
        identical(verifier_evidence[["tolerances", exact = TRUE]],
                  raw[["tolerances", exact = TRUE]]) &&
        identical(verifier_evidence[["M", exact = TRUE]],
                  orders[["M_verification_used", exact = TRUE]]),
      "diagnostic_snapshot_identity", "result.verification",
      paste(
        "selected public evidence at M_selected and fixed-parameter",
        "component verifier evidence at M_verification"
      ), verification
    )
    expected_residuals <- list(
      diagnostics = expected_evidence[["delta", exact = TRUE]]
    )
    .dpprior_schema_require(
      identical(raw[["residuals", exact = TRUE]], expected_residuals),
      "diagnostic_component_identity", "result.achieved",
      paste(
        "public achieved values and retained refinement deltas identical",
        "to independently recomputed component diagnostics"
      ), list(
        achieved = raw[["achieved", exact = TRUE]],
        residuals = raw[["residuals", exact = TRUE]]
      )
    )
    return(invisible(TRUE))
  }

  if (identical(mode, "elicitation_sensitivity")) {
    sensitivity <- raw[["sensitivity", exact = TRUE]]
    .dpprior_schema_require(
      is.null(raw[["parameters", exact = TRUE]]) &&
        all(vapply(
          orders[c(
            "M_requested", "M_selected", "M_verification_required",
            "M_verification_used"
          )],
          is.null, logical(1)
        )),
      "sensitivity_parameter_order_authority", "result",
      paste(
        "no public fit parameters and no quadrature-order claims for a",
        "table-reconciliation object"
      ),
      list(parameters = raw[["parameters", exact = TRUE]], orders = orders)
    )
    keys <- sensitivity[["scenarios", exact = TRUE]][[
      "scenario_key", exact = TRUE
    ]]
    metrics <- sensitivity[["metrics_long", exact = TRUE]]
    row_status <- sensitivity[["scenario_results", exact = TRUE]][[
      "status", exact = TRUE
    ]]
    row_usable <- sensitivity[["scenario_results", exact = TRUE]][[
      "usable", exact = TRUE
    ]]
    row_verified <- sensitivity[["scenario_results", exact = TRUE]][[
      "verified", exact = TRUE
    ]]
    precedence <- c(
      converged = 1L, boundary = 2L, approximate = 3L,
      infeasible = 4L, failed = 5L
    )
    expected_status <- row_status[[
      which.max(unname(precedence[row_status]))
    ]]
    expected_verified <- expected_status %in%
      c("converged", "boundary", "infeasible") &&
      all(row_verified)
    expected_usable <- expected_status %in% c("converged", "boundary") &&
      all(row_usable)
    .dpprior_schema_require(
      identical(status, expected_status) &&
        identical(raw[["verified", exact = TRUE]], expected_verified) &&
        identical(raw[["usable", exact = TRUE]], expected_usable),
      "sensitivity_status_aggregation", "result.status",
      "top status/flags recomputed from every canonical scenario row",
      raw[c("status", "usable", "verified")]
    )
    expected_settings <- list(
      scenario_count = as.integer(length(keys)), scenario_keys = keys,
      metric_names = .DPPRIOR_SENSITIVITY_METRICS
    )
    .dpprior_schema_require(
      identical(verification[["method", exact = TRUE]], "reconciliation") &&
        verification[["performed", exact = TRUE]] &&
        identical(verification[["passed", exact = TRUE]], expected_verified) &&
        identical(verification[["settings", exact = TRUE]], expected_settings) &&
        is.null(verification[["stability", exact = TRUE]]) &&
        !is.null(verification[["selected_snapshot", exact = TRUE]]) &&
        !is.null(verification[["verifier_snapshot", exact = TRUE]]) &&
        identical(
          verification[["selected_snapshot", exact = TRUE]][["source", exact = TRUE]],
          "sensitivity_tables"
        ) && identical(
          verification[["verifier_snapshot", exact = TRUE]][["source", exact = TRUE]],
          "sensitivity_reconciliation"
        ),
      "sensitivity_verification_contract", "result.verification",
      "the exact independent table-reconciliation contract", verification
    )
    validate_checks(
      components, "reconciliation", "sensitivity_reconciliation",
      "result.verification.components"
    )
    invariant_names <- c(
      "unique_keys", "lexicographic_order", "finite_or_reason",
      "failure_preservation"
    )
    validate_checks(
      invariants, invariant_names, "sensitivity_reconciliation",
      "result.verification.invariants"
    )
    expected_metric_keys <- rep(
      keys, each = length(.DPPRIOR_SENSITIVITY_METRICS)
    )
    truth <- c(
      scenario_key_identity = identical(
        sensitivity[["scenario_results", exact = TRUE]][[
          "scenario_key", exact = TRUE
        ]], keys
      ),
      condition_key_identity = identical(
        names(sensitivity[["conditions", exact = TRUE]]), keys
      ),
      interval_key_identity = identical(
        names(sensitivity[["interval_checks", exact = TRUE]]), keys
      ),
      metric_grid_identity = identical(
        metrics[["scenario_key", exact = TRUE]], expected_metric_keys
      ) && identical(
        metrics[["metric", exact = TRUE]],
        rep(.DPPRIOR_SENSITIVITY_METRICS, times = length(keys))
      ),
      row_status_identity = all(vapply(seq_along(keys), function(index) {
        rows <- metrics[["scenario_key", exact = TRUE]] == keys[[index]]
        all(metrics[["status", exact = TRUE]][rows] == row_status[[index]]) &&
          all(metrics[["usable", exact = TRUE]][rows] == row_usable[[index]]) &&
          all(metrics[["verified", exact = TRUE]][rows] == row_verified[[index]])
      }, logical(1))),
      top_status_identity = identical(status, expected_status),
      global_summary_identity = identical(
        sensitivity[["global", exact = TRUE]][[
          "scenario_count", exact = TRUE
        ]], as.integer(length(keys))
      ) && identical(
        sensitivity[["global", exact = TRUE]][[
          "metric_count", exact = TRUE
        ]], as.integer(nrow(metrics))
      )
    )
    .dpprior_bind_decision_check(
      components[["reconciliation", exact = TRUE]], truth,
      setNames(rep(TRUE, length(truth)), names(truth)), NULL, "identical",
      "result.verification.components.reconciliation"
    )
    finite_or_reason <- all(
      is.finite(metrics[["value", exact = TRUE]]) |
        (is.na(metrics[["value", exact = TRUE]]) &
           !is.nan(metrics[["value", exact = TRUE]]) &
           !is.na(metrics[["reason", exact = TRUE]]) &
           nzchar(metrics[["reason", exact = TRUE]]))
    )
    failed_keys <- keys[row_status %in% c("failed", "infeasible")]
    failure_preservation <- all(vapply(failed_keys, function(key) {
      rows <- metrics[["scenario_key", exact = TRUE]] == key
      key %in% names(sensitivity[["conditions", exact = TRUE]]) &&
        any(
          !is.na(metrics[["reason", exact = TRUE]][rows]) &
            nzchar(metrics[["reason", exact = TRUE]][rows])
        )
    }, logical(1)))
    invariant_truth <- c(
      unique_keys = !anyDuplicated(keys),
      lexicographic_order = identical(keys, sort(keys, method = "radix")),
      finite_or_reason = finite_or_reason,
      failure_preservation = failure_preservation
    )
    for (name in invariant_names) {
      .dpprior_bind_decision_check(
        invariants[[name, exact = TRUE]], invariant_truth[[name]], TRUE, NULL,
        "identical", paste0("result.verification.invariants.", name)
      )
    }
    sensitivity_target <- raw[["target", exact = TRUE]]
    .dpprior_schema_exact_names(
      sensitivity_target, "defaults", "result.target"
    )
    .dpprior_schema_exact_names(
      sensitivity_target[["defaults", exact = TRUE]], "J",
      "result.target.defaults"
    )
    .dpprior_schema_require(
      identical(
        sensitivity_target[["defaults", exact = TRUE]][["J", exact = TRUE]],
        raw[["J", exact = TRUE]]
      ) && identical(
        sensitivity[["metadata", exact = TRUE]][["J", exact = TRUE]],
        raw[["J", exact = TRUE]]
      ),
      "sensitivity_J_identity", "result.target.defaults.J",
      "one exact J across result, target defaults, metadata, and content",
      list(
        result = raw[["J", exact = TRUE]],
        target = sensitivity_target[["defaults", exact = TRUE]][[
          "J", exact = TRUE
        ]],
        metadata = sensitivity[["metadata", exact = TRUE]][["J", exact = TRUE]]
      )
    )
    for (snapshot_name in c("selected_snapshot", "verifier_snapshot")) {
      snapshot <- verification[[snapshot_name, exact = TRUE]]
      .dpprior_schema_require(
        is.null(snapshot[["parameters", exact = TRUE]]) &&
          is.null(snapshot[["M", exact = TRUE]]) &&
          snapshot[["finite", exact = TRUE]] &&
          identical(snapshot[["achieved", exact = TRUE]],
                  raw[["achieved", exact = TRUE]]) &&
          identical(snapshot[["residuals", exact = TRUE]],
                    raw[["residuals", exact = TRUE]]) &&
          identical(snapshot[["tolerances", exact = TRUE]],
                    raw[["tolerances", exact = TRUE]]),
        "sensitivity_snapshot_identity",
        paste0("result.verification.", snapshot_name),
        paste(
          "parameter-free, order-free finite table reconciliation summaries",
          "identical in both evidence snapshots"
        ),
        snapshot
      )
    }
    .dpprior_schema_require(
      identical(raw[["achieved", exact = TRUE]][[
        "scenario_count", exact = TRUE
      ]], as.integer(length(keys))) &&
        identical(raw[["residuals", exact = TRUE]][[
          "missing_metric_count", exact = TRUE
        ]], as.integer(sum(is.na(metrics[["value", exact = TRUE]])))) &&
        identical(raw[["tolerances", exact = TRUE]][[
          "expected_metric_count", exact = TRUE
        ]], as.integer(length(keys) * length(.DPPRIOR_SENSITIVITY_METRICS))) &&
        identical(sensitivity[["metadata", exact = TRUE]][["J", exact = TRUE]],
                  raw[["J", exact = TRUE]]) &&
        length(computation[["attempts", exact = TRUE]]) == 0L &&
        length(computation[["candidate_evaluations", exact = TRUE]]) == 0L &&
        is.null(computation[["selected_attempt_id", exact = TRUE]]) &&
        is.null(computation[["selected_candidate_id", exact = TRUE]]) &&
        identical(
          computation[["termination", exact = TRUE]][["code", exact = TRUE]],
          "deterministic"
        ) && identical(
          computation[["termination", exact = TRUE]][["source", exact = TRUE]],
          "sensitivity_reconciliation"
        ),
      "sensitivity_reconciliation", "result.sensitivity",
      "public summaries, J, empty execution ledger, and deterministic termination",
      sensitivity
    )
    return(invisible(TRUE))
  }

  if (identical(mode, "dual_hard") && identical(status, "infeasible")) {
    .dpprior_schema_require(
      identical(
        verification[["method", exact = TRUE]],
        "analytic_global_monotonicity_corner_certificate"
      ) && verification[["performed", exact = TRUE]] &&
        verification[["passed", exact = TRUE]] &&
        identical(verification[["settings", exact = TRUE]], list()) &&
        is.null(verification[["selected_snapshot", exact = TRUE]]) &&
        is.null(verification[["verifier_snapshot", exact = TRUE]]) &&
        is.null(verification[["stability", exact = TRUE]]),
      "verification_contract", "result.verification",
      paste(
        "the certificate-only hard verification method with empty settings",
        "and no fabricated numerical snapshots or stability"
      ), verification
    )
    validate_checks(
      components, "infeasibility_certificate", "analytic_certificate",
      "result.verification.components"
    )
    validate_checks(
      invariants, "domain_monotonicity", "analytic_certificate",
      "result.verification.invariants"
    )
    .dpprior_bind_decision_check(
      invariants[["domain_monotonicity", exact = TRUE]], TRUE, TRUE, NULL,
      "identical", "result.verification.invariants.domain_monotonicity"
    )
    return(invisible(TRUE))
  }

  if (raw[["object_type", exact = TRUE]] == "fit" && !finite_fit &&
      identical(status, "failed")) {
    .dpprior_schema_require(
      identical(verification[["method", exact = TRUE]], "no_candidate") &&
        !verification[["performed", exact = TRUE]] &&
        !verification[["passed", exact = TRUE]] &&
        identical(verification[["settings", exact = TRUE]], list()) &&
        is.null(verification[["selected_snapshot", exact = TRUE]]) &&
        is.null(verification[["verifier_snapshot", exact = TRUE]]) &&
        is.null(verification[["stability", exact = TRUE]]) &&
        length(components) == 0L,
      "verification_contract", "result.verification",
      "the exact no-candidate verification record without decision claims",
      verification
    )
    validate_checks(
      invariants, "no_public_candidate", "independent_verifier",
      "result.verification.invariants"
    )
    no_public_candidate <- is.null(raw[["parameters", exact = TRUE]]) &&
      is.null(computation[["selected_candidate_id", exact = TRUE]]) &&
      is.null(computation[["selected_attempt_id", exact = TRUE]])
    .dpprior_bind_decision_check(
      invariants[["no_public_candidate", exact = TRUE]],
      no_public_candidate, TRUE, NULL, "identical",
      "result.verification.invariants.no_public_candidate"
    )
    return(invisible(TRUE))
  }

  if (!finite_fit || !mode %in%
      c("a2_moment", "a2_kl", "dual_hard", "dual_soft")) {
    return(invisible(TRUE))
  }

  contract <- switch(
    mode,
    a2_moment = list(
      method = "independent_higher_order_moment_recomputation",
      settings = list(
        M_selected = orders[["M_selected", exact = TRUE]],
        M_verification_required = orders[[
          "M_verification_required", exact = TRUE
        ]],
        M_verification = orders[["M_verification_used", exact = TRUE]]
      ),
      components = c(
        "residual_adequacy", "order_stability", "parameter_identity"
      ),
      invariants = c("finite_parameters", "K_support")
    ),
    a2_kl = list(
      method = "fresh higher-order marginal PMF and direct-moment audit",
      settings = list(
        M_selected = orders[["M_selected", exact = TRUE]],
        M_verification = orders[["M_verification_used", exact = TRUE]],
        M_verification_required = orders[[
          "M_verification_required", exact = TRUE
        ]],
        pmf_abs_tol = raw[["tolerances", exact = TRUE]][[
          "distribution", exact = TRUE
        ]][["order", exact = TRUE]][["pmf_absolute", exact = TRUE]],
        pmf_rel_tol = raw[["tolerances", exact = TRUE]][[
          "distribution", exact = TRUE
        ]][["order", exact = TRUE]][["pmf_relative", exact = TRUE]]
      ),
      components = c(
        "target_identity", "pmf_adequacy", "order_stability",
        "candidate_selection"
      ),
      invariants = c(
        "finite_parameters", "parameter_identity", "pmf_probability",
        "support_identity"
      )
    ),
    dual_hard = list(
      method = "fresh_higher_order_recomputation_and_local_perturbation",
      settings = list(
        M_fit = orders[["M_selected", exact = TRUE]],
        M_verify = orders[["M_verification_used", exact = TRUE]],
        M_verify_required = orders[["M_verification_required", exact = TRUE]],
        log_bounds = controls[["log_bounds", exact = TRUE]],
        verification_abs_tol = controls[["verification_abs_tol", exact = TRUE]],
        verification_rel_tol = controls[["verification_rel_tol", exact = TRUE]],
        perturbation_step = controls[["perturbation_step", exact = TRUE]],
        perturbation_abs_tol = controls[["perturbation_abs_tol", exact = TRUE]]
      ),
      components = c(
        "constraint_selected", "constraint_refined", "order_stability",
        "metric_certification", "candidate_selection", "perturbation"
      ),
      invariants = c(
        "probability", "K_support", "finite_parameters_inside_domain"
      )
    ),
    dual_soft = list(
      method = "fresh higher-order K moments and named weight metric",
      settings = list(
        M_fit = orders[["M_selected", exact = TRUE]],
        M_verify = orders[["M_verification_used", exact = TRUE]],
        log_bounds = controls[["log_bounds", exact = TRUE]],
        verification_abs_tol = controls[["verification_abs_tol", exact = TRUE]],
        verification_rel_tol = controls[["verification_rel_tol", exact = TRUE]],
        objective_abs_tol = controls[["objective_abs_tol", exact = TRUE]],
        objective_rel_tol = controls[["objective_rel_tol", exact = TRUE]],
        boundary_tol = controls[["boundary_tol", exact = TRUE]],
        stationarity_step = controls[["stationarity_step", exact = TRUE]],
        stationarity_tol = controls[["stationarity_tol", exact = TRUE]],
        neighborhood_step = controls[["local_neighbor_step", exact = TRUE]]
      ),
      components = if (isTRUE(raw[["tradeoff", exact = TRUE]][[
        "endpoint", exact = TRUE
      ]])) c("endpoint_input_identity", "order_stability") else c(
        "objective_recomputation", "order_stability", "local_optimality",
        "candidate_selection"
      ),
      invariants = c(
        "probability", "K_support", "finite_parameters_inside_domain",
        "fixed_input_scales"
      )
    )
  )
  .dpprior_schema_require(
    verification[["performed", exact = TRUE]] &&
      identical(verification[["method", exact = TRUE]], contract$method) &&
      identical(verification[["settings", exact = TRUE]], contract$settings),
    "verification_contract", "result.verification",
    "the exact mode-specific independent verification method and settings",
    list(
      method = verification[["method", exact = TRUE]],
      settings = verification[["settings", exact = TRUE]]
    )
  )
  .dpprior_schema_require(
    !is.null(verification[["selected_snapshot", exact = TRUE]]) &&
      !is.null(verification[["verifier_snapshot", exact = TRUE]]) &&
      !is.null(verification[["stability", exact = TRUE]]) &&
      identical(
        verification[["selected_snapshot", exact = TRUE]][["source", exact = TRUE]],
        "selected_order"
      ) && identical(
        verification[["verifier_snapshot", exact = TRUE]][["source", exact = TRUE]],
        "independent_verifier"
      ) && identical(
        verification[["stability", exact = TRUE]][["source", exact = TRUE]],
        "independent_verifier"
      ),
    "verification_independence", "result.verification",
    paste(
      "selected-order evidence plus an independent verifier snapshot and",
      "independent stability source, never optimizer self-report"
    ), verification
  )
  validate_checks(
    components, contract$components, "independent_verifier",
    "result.verification.components"
  )
  validate_checks(
    invariants, contract$invariants, "independent_verifier",
    "result.verification.invariants"
  )
  invariant_truth <- switch(
    mode,
    a2_moment = c(
      finite_parameters = all(is.finite(unlist(
        raw[["parameters", exact = TRUE]][c("a", "b")], use.names = FALSE
      ))),
      K_support = TRUE
    ),
    a2_kl = c(
      finite_parameters = all(is.finite(unlist(
        raw[["parameters", exact = TRUE]][c("a", "b")], use.names = FALSE
      ))),
      parameter_identity = identical(
        verification[["selected_snapshot", exact = TRUE]][[
          "parameters", exact = TRUE
        ]],
        verification[["verifier_snapshot", exact = TRUE]][[
          "parameters", exact = TRUE
        ]]
      ),
      pmf_probability = TRUE,
      support_identity = TRUE
    ),
    dual_hard = c(
      probability = raw[["achieved", exact = TRUE]][[
        "weight", exact = TRUE
      ]][["value", exact = TRUE]] >= 0 &&
        raw[["achieved", exact = TRUE]][["weight", exact = TRUE]][[
          "value", exact = TRUE
        ]] <= 1,
      K_support = TRUE,
      finite_parameters_inside_domain = all(
        log(unlist(raw[["parameters", exact = TRUE]][c("a", "b")])) >=
          controls[["log_bounds", exact = TRUE]][[1L]] &
          log(unlist(raw[["parameters", exact = TRUE]][c("a", "b")])) <=
          controls[["log_bounds", exact = TRUE]][[2L]]
      )
    ),
    dual_soft = c(
      probability = raw[["achieved", exact = TRUE]][[
        "weight", exact = TRUE
      ]][["value", exact = TRUE]] >= 0 &&
        raw[["achieved", exact = TRUE]][["weight", exact = TRUE]][[
          "value", exact = TRUE
        ]] <= 1,
      K_support = TRUE,
      finite_parameters_inside_domain = all(
        log(unlist(raw[["parameters", exact = TRUE]][c("a", "b")])) >=
          controls[["log_bounds", exact = TRUE]][[1L]] &
          log(unlist(raw[["parameters", exact = TRUE]][c("a", "b")])) <=
          controls[["log_bounds", exact = TRUE]][[2L]]
      ),
      fixed_input_scales = isTRUE(computation[["scaling", exact = TRUE]][[
        "fixed_from_input", exact = TRUE
      ]])
    )
  )
  for (name in contract$invariants) {
    .dpprior_bind_decision_check(
      invariants[[name, exact = TRUE]], invariant_truth[[name]], TRUE, NULL,
      "identical", paste0("result.verification.invariants.", name)
    )
  }
  invisible(TRUE)
}


.dpprior_schema_find_forbidden_names <- function(x, forbidden, path = "result") {
  if (typeof(x) != "list" || !is.list(x)) {
    return(character())
  }
  raw <- if (is.object(x)) unclass(x) else x
  object_names <- names(raw)
  hits <- character()
  if (!is.null(object_names)) {
    matched <- which(object_names %in% forbidden)
    if (length(matched) > 0L) {
      hits <- paste0(path, ".", object_names[matched])
    }
  }
  for (i in seq_along(raw)) {
    child_name <- if (!is.null(object_names) && nzchar(object_names[[i]])) {
      object_names[[i]]
    } else {
      sprintf("[[%d]]", i)
    }
    hits <- c(
      hits,
      .dpprior_schema_find_forbidden_names(
        raw[[i]], forbidden, paste0(path, ".", child_name)
      )
    )
  }
  unique(hits)
}


# --- Canonical result validation --------------------------------------------

.dpprior_validate_achieved_K <- function(x, J, path = "achieved.K") {
  .dpprior_schema_validate_named_list(x, path, allow_empty = FALSE)
  required <- c("mean", "variance", "estimand", "source", "M")
  .dpprior_schema_require(
    all(required %in% names(x)), "achieved_K", path,
    paste("required fields", paste(required, collapse = ", ")), names(x)
  )
  .dpprior_schema_require(
    all(names(x) %in% c(required, "pmf")), "achieved_K", path,
    "only canonical K fields (plus optional pmf)", names(x)
  )
  .dpprior_schema_require(
    .dpprior_schema_is_count(J, 1L), "count", paste0(path, ".J"),
    "an authoritative integer J at least 1", J
  )
  .dpprior_schema_validate_finite_scalar(
    x$mean, paste0(path, ".mean"), lower = 1, upper = J
  )
  .dpprior_schema_validate_finite_scalar(
    x$variance, paste0(path, ".variance"), lower = 0
  )
  variance_bound <- (x$mean - 1) * (J - x$mean)
  moment_tolerance <- 1e-8 * max(1, abs(variance_bound), abs(x$variance))
  .dpprior_schema_require(
    x$variance <= variance_bound + moment_tolerance,
    "K_variance_bound", paste0(path, ".variance"),
    "variance <= (mean-1)*(J-mean) for support 1:J", x$variance
  )
  .dpprior_schema_validate_scalar_character(
    x$estimand, paste0(path, ".estimand")
  )
  .dpprior_schema_validate_scalar_character(x$source, paste0(path, ".source"))
  .dpprior_schema_validate_count_or_null(x$M, paste0(path, ".M"))
  if ("pmf" %in% names(x)) {
    .dpprior_schema_require(
      is.numeric(x$pmf) && !is.object(x$pmf) && is.null(dim(x$pmf)) &&
        .dpprior_schema_has_only_attributes(x$pmf) &&
        length(x$pmf) == J && !anyNA(x$pmf) && all(is.finite(x$pmf)) &&
        all(x$pmf >= 0),
      "pmf", paste0(path, ".pmf"),
      "a finite nonnegative numeric vector of length J", x$pmf
    )
    .dpprior_schema_require(
      abs(sum(x$pmf) - 1) <= .TOL_PMF_SUM,
      "pmf_normalization", paste0(path, ".pmf"),
      sprintf("mass 1 within %.3g", .TOL_PMF_SUM), sum(x$pmf)
    )
    support <- seq_len(J)
    pmf_mean <- sum(support * x$pmf)
    pmf_variance <- sum((support - pmf_mean)^2 * x$pmf)
    mean_tolerance <- 1e-8 * max(1, abs(pmf_mean), abs(x$mean))
    variance_tolerance <- 1e-8 * max(
      1, abs(pmf_variance), abs(x$variance)
    )
    .dpprior_schema_require(
      abs(x$mean - pmf_mean) <= mean_tolerance &&
        abs(x$variance - pmf_variance) <= variance_tolerance,
      "pmf_moments", path,
      "mean/variance agreeing with optional PMF moments",
      list(recorded = c(x$mean, x$variance),
           recomputed = c(pmf_mean, pmf_variance))
    )
  }
  invisible(TRUE)
}


.dpprior_validate_achieved_weight <- function(x, expected_metric,
                                              path = "achieved.weight") {
  .dpprior_schema_exact_names(x, c("metric", "value", "source"), path)
  .dpprior_schema_validate_scalar_character(
    x[["metric", exact = TRUE]], paste0(path, ".metric")
  )
  .dpprior_schema_require(
    identical(x[["metric", exact = TRUE]], expected_metric),
    "weight_metric", paste0(path, ".metric"),
    "identity with the canonical weight-target metric",
    x[["metric", exact = TRUE]]
  )
  .dpprior_schema_validate_finite_scalar(
    x[["value", exact = TRUE]], paste0(path, ".value"), lower = 0, upper = 1
  )
  .dpprior_schema_validate_scalar_character(
    x[["source", exact = TRUE]], paste0(path, ".source")
  )
  invisible(TRUE)
}


.dpprior_validate_result_target <- function(x, mode, object_type,
                                            path = "target") {
  .dpprior_schema_validate_named_list(x, path)
  if (!identical(object_type, "fit")) {
    return(invisible(TRUE))
  }
  expected <- if (mode %in% c("dual_hard", "dual_soft", "dual_legacy")) {
    c("K", "weight")
  } else {
    "K"
  }
  .dpprior_schema_require(
    identical(names(x), expected), "target_bundle", path,
    paste(expected, collapse = ", "), names(x)
  )
  .dpprior_validate_target_v1(x[["K", exact = TRUE]])
  if ("weight" %in% expected) {
    .dpprior_validate_weight_target_v1(x[["weight", exact = TRUE]])
  }
  invisible(TRUE)
}


.dpprior_validate_truth_numeric_record <- function(x, fields, path,
                                                    positive = character()) {
  .dpprior_schema_exact_names(x, fields, path)
  for (field in fields) {
    .dpprior_schema_validate_finite_scalar(
      x[[field, exact = TRUE]], paste0(path, ".", field), lower = 0,
      lower_open = field %in% positive
    )
  }
  invisible(TRUE)
}


.dpprior_validate_result_truth_authority <- function(raw) {
  mode <- raw[["mode", exact = TRUE]]
  controls <- raw[["computation", exact = TRUE]][[
    "used", exact = TRUE
  ]][["controls", exact = TRUE]]
  tolerances <- raw[["tolerances", exact = TRUE]]
  control_path <- "result.computation.used.controls"
  tolerance_path <- "result.tolerances"

  validate_bounds <- function(x, path) {
    .dpprior_schema_require(
      is.numeric(x) && !is.object(x) && is.null(dim(x)) &&
        .dpprior_schema_has_only_attributes(x) && length(x) == 2L &&
        !anyNA(x) && all(is.finite(x)) && x[[1L]] < x[[2L]],
      "truth_bounds", path, "two finite increasing unclassed bounds", x
    )
  }
  validate_cap <- function(value, maximum, path, exact = FALSE) {
    .dpprior_schema_validate_finite_scalar(value, path, lower = 0)
    .dpprior_schema_require(
      if (exact) identical(value, maximum) else value <= maximum,
      "truth_control_policy", path,
      if (exact) paste("exactly", format(maximum)) else
        paste("no greater than", format(maximum)),
      value
    )
  }

  if (identical(mode, "a2_moment")) {
    .dpprior_schema_exact_names(
      controls,
      c(
        "max_iter", "damping", "use_fallback", "tol_step", "log_bounds",
        "boundary_tol", "line_search_max", "jacobian_rcond_singular",
        "jacobian_rcond_ill", "fallback", "selection_tolerance"
      ), control_path
    )
    for (field in c("max_iter", "line_search_max")) {
      .dpprior_schema_require(
        .dpprior_schema_is_count(controls[[field, exact = TRUE]], 1L),
        "truth_control", paste0(control_path, ".", field),
        "a positive integer", controls[[field, exact = TRUE]]
      )
    }
    for (field in c("damping", "use_fallback")) {
      .dpprior_schema_validate_scalar_logical(
        controls[[field, exact = TRUE]], paste0(control_path, ".", field)
      )
    }
    validate_bounds(controls[["log_bounds", exact = TRUE]],
                    paste0(control_path, ".log_bounds"))
    validate_cap(controls[["tol_step", exact = TRUE]], 1e-10,
                 paste0(control_path, ".tol_step"))
    validate_cap(controls[["boundary_tol", exact = TRUE]], 1e-6,
                 paste0(control_path, ".boundary_tol"), exact = TRUE)
    validate_cap(
      controls[["jacobian_rcond_singular", exact = TRUE]], 1e-12,
      paste0(control_path, ".jacobian_rcond_singular"), exact = TRUE
    )
    validate_cap(
      controls[["jacobian_rcond_ill", exact = TRUE]],
      sqrt(.Machine$double.eps),
      paste0(control_path, ".jacobian_rcond_ill"), exact = TRUE
    )
    .dpprior_schema_exact_names(
      controls[["fallback", exact = TRUE]],
      c("maxit", "reltol", "finite_penalty"),
      paste0(control_path, ".fallback")
    )
    .dpprior_schema_require(
      .dpprior_schema_is_count(
        controls[["fallback", exact = TRUE]][["maxit", exact = TRUE]], 1L
      ),
      "truth_control", paste0(control_path, ".fallback.maxit"),
      "a positive integer", controls[["fallback", exact = TRUE]]
    )
    validate_cap(
      controls[["fallback", exact = TRUE]][["reltol", exact = TRUE]],
      1e-12, paste0(control_path, ".fallback.reltol")
    )
    validate_cap(
      controls[["fallback", exact = TRUE]][["finite_penalty", exact = TRUE]],
      .Machine$double.xmax / 1024,
      paste0(control_path, ".fallback.finite_penalty"), exact = TRUE
    )
    validate_cap(
      controls[["selection_tolerance", exact = TRUE]], 0,
      paste0(control_path, ".selection_tolerance"), exact = TRUE
    )
    .dpprior_schema_exact_names(
      tolerances, c("K_adequacy", "K_stability", "step", "boundary"),
      tolerance_path
    )
    .dpprior_schema_exact_names(
      tolerances[["K_adequacy", exact = TRUE]],
      c("absolute", "relative", "scale_formula"),
      paste0(tolerance_path, ".K_adequacy")
    )
    validate_cap(
      tolerances[["K_adequacy", exact = TRUE]][["absolute", exact = TRUE]],
      1e-8, paste0(tolerance_path, ".K_adequacy.absolute")
    )
    validate_cap(
      tolerances[["K_adequacy", exact = TRUE]][["relative", exact = TRUE]],
      1e-8, paste0(tolerance_path, ".K_adequacy.relative")
    )
    .dpprior_schema_require(
      identical(
        tolerances[["K_adequacy", exact = TRUE]][[
          "scale_formula", exact = TRUE
        ]], "max(abs(target),1)"
      ),
      "truth_formula", paste0(tolerance_path, ".K_adequacy.scale_formula"),
      "max(abs(target),1)",
      tolerances[["K_adequacy", exact = TRUE]][["scale_formula", exact = TRUE]]
    )
    .dpprior_validate_truth_numeric_record(
      tolerances[["K_stability", exact = TRUE]],
      c("absolute", "relative", "scale_floor"),
      paste0(tolerance_path, ".K_stability")
    )
    .dpprior_schema_require(
      tolerances[["K_stability", exact = TRUE]][["absolute", exact = TRUE]] <=
        1e-10 &&
        tolerances[["K_stability", exact = TRUE]][["relative", exact = TRUE]] <=
          1e-8 &&
        identical(
          tolerances[["K_stability", exact = TRUE]][[
            "scale_floor", exact = TRUE
          ]], 1
        ) && identical(tolerances[["step", exact = TRUE]],
                       controls[["tol_step", exact = TRUE]]) &&
        identical(tolerances[["boundary", exact = TRUE]],
                  controls[["boundary_tol", exact = TRUE]]),
      "truth_tolerance_identity", tolerance_path,
      "separate capped A2-MN adequacy/stability, step, and boundary controls",
      tolerances
    )
  }

  if (identical(mode, "a2_kl")) {
    .dpprior_schema_exact_names(
      controls,
      c(
        "max_iter", "optimizer_tol", "log_bounds", "boundary_tol",
        "fallback_max_iter", "primary", "fallback",
        "fallback_trigger_worse_than_start", "selection_tolerance"
      ), control_path
    )
    for (field in c("max_iter", "fallback_max_iter")) {
      .dpprior_schema_require(
        .dpprior_schema_is_count(controls[[field, exact = TRUE]], 1L),
        "truth_control", paste0(control_path, ".", field),
        "a positive integer", controls[[field, exact = TRUE]]
      )
    }
    validate_bounds(controls[["log_bounds", exact = TRUE]],
                    paste0(control_path, ".log_bounds"))
    for (field in c("optimizer_tol", "boundary_tol")) {
      validate_cap(controls[[field, exact = TRUE]], 1e-6,
                   paste0(control_path, ".", field))
    }
    validate_cap(
      controls[["fallback_trigger_worse_than_start", exact = TRUE]], 1e-12,
      paste0(control_path, ".fallback_trigger_worse_than_start"), exact = TRUE
    )
    validate_cap(
      controls[["selection_tolerance", exact = TRUE]], 0,
      paste0(control_path, ".selection_tolerance"), exact = TRUE
    )
    .dpprior_schema_exact_names(
      controls[["primary", exact = TRUE]], c("maxit", "factr", "pgtol"),
      paste0(control_path, ".primary")
    )
    .dpprior_schema_exact_names(
      controls[["fallback", exact = TRUE]],
      c("iter.max", "eval.max", "rel.tol", "x.tol"),
      paste0(control_path, ".fallback")
    )
    .dpprior_schema_require(
      identical(
        controls[["primary", exact = TRUE]][["maxit", exact = TRUE]],
        controls[["max_iter", exact = TRUE]]
      ) && identical(
        controls[["primary", exact = TRUE]][["factr", exact = TRUE]],
        controls[["optimizer_tol", exact = TRUE]] / .Machine$double.eps
      ) && identical(
        controls[["primary", exact = TRUE]][["pgtol", exact = TRUE]],
        controls[["optimizer_tol", exact = TRUE]]
      ) && identical(
        controls[["fallback", exact = TRUE]][["iter.max", exact = TRUE]],
        controls[["fallback_max_iter", exact = TRUE]]
      ) && identical(
        controls[["fallback", exact = TRUE]][["eval.max", exact = TRUE]],
        as.integer(max(200L, 2L * controls[["fallback_max_iter", exact = TRUE]]))
      ) && identical(
        controls[["fallback", exact = TRUE]][["rel.tol", exact = TRUE]],
        controls[["optimizer_tol", exact = TRUE]]
      ) && identical(
        controls[["fallback", exact = TRUE]][["x.tol", exact = TRUE]],
        controls[["optimizer_tol", exact = TRUE]]
      ),
      "truth_optimizer_controls", control_path,
      "normalized A2-KL primary/fallback controls derived from fixed inputs",
      controls
    )
    .dpprior_schema_exact_names(
      tolerances, c("distribution", "boundary"), tolerance_path
    )
    distribution <- tolerances[["distribution", exact = TRUE]]
    .dpprior_schema_exact_names(
      distribution, c("adequacy", "order"),
      paste0(tolerance_path, ".distribution")
    )
    adequacy <- distribution[["adequacy", exact = TRUE]]
    .dpprior_schema_exact_names(
      adequacy,
      c(
        "kl", "l1", "mean_scaled", "variance_scaled",
        "mean_scale_formula", "variance_scale_formula"
      ), paste0(tolerance_path, ".distribution.adequacy")
    )
    for (pair in list(
      c("kl", 0.015), c("l1", 0.11), c("mean_scaled", 0.01),
      c("variance_scaled", 0.065)
    )) {
      validate_cap(
        adequacy[[pair[[1L]], exact = TRUE]], as.numeric(pair[[2L]]),
        paste0(tolerance_path, ".distribution.adequacy.", pair[[1L]])
      )
    }
    .dpprior_schema_require(
      identical(adequacy[["mean_scale_formula", exact = TRUE]],
                "max(1,sqrt(target_variance))") &&
        identical(adequacy[["variance_scale_formula", exact = TRUE]],
                  "max(1,target_variance)"),
      "truth_formula", paste0(tolerance_path, ".distribution.adequacy"),
      "the exact Phase 8 mean/variance scale formulas", adequacy
    )
    order <- distribution[["order", exact = TRUE]]
    .dpprior_validate_truth_numeric_record(
      order,
      c(
        "pmf_absolute", "pmf_relative", "pmf_l1",
        "direct_moment_absolute", "direct_moment_relative",
        "target_identity_l1"
      ), paste0(tolerance_path, ".distribution.order")
    )
    .dpprior_schema_require(
      order[["pmf_absolute", exact = TRUE]] <= 1e-10 &&
        order[["pmf_relative", exact = TRUE]] <= 1e-8 &&
        identical(
          order[["pmf_l1", exact = TRUE]],
          order[["pmf_absolute", exact = TRUE]] +
            order[["pmf_relative", exact = TRUE]]
        ) && identical(
          order[["direct_moment_absolute", exact = TRUE]],
          order[["pmf_absolute", exact = TRUE]]
        ) && identical(
          order[["direct_moment_relative", exact = TRUE]],
          order[["pmf_relative", exact = TRUE]]
        ) && identical(order[["target_identity_l1", exact = TRUE]],
                       .TOL_PMF_SUM) &&
        identical(tolerances[["boundary", exact = TRUE]],
                  controls[["boundary_tol", exact = TRUE]]),
      "truth_tolerance_identity", paste0(tolerance_path, ".distribution"),
      "four adequacy gates plus a separate fixed PMF order/identity contract",
      distribution
    )
  }

  if (identical(mode, "dual_hard")) {
    .dpprior_schema_exact_names(
      controls,
      c(
        "maxit", "scan_points", "scan_keep", "profile_starts", "log_bounds",
        "root_tol", "optim_reltol", "penalty", "boundary_tol", "constraint_abs_tol",
        "constraint_rel_tol", "verification_abs_tol", "verification_rel_tol",
        "perturbation_step", "perturbation_abs_tol", "auto_cap"
      ), control_path
    )
    for (field in c("maxit", "scan_points", "scan_keep", "profile_starts")) {
      .dpprior_schema_require(
        .dpprior_schema_is_count(controls[[field, exact = TRUE]], 1L),
        "truth_control", paste0(control_path, ".", field),
        "a positive integer", controls[[field, exact = TRUE]]
      )
    }
    validate_bounds(controls[["log_bounds", exact = TRUE]],
                    paste0(control_path, ".log_bounds"))
    .dpprior_schema_require(
      controls[["log_bounds", exact = TRUE]][[1L]] >= -.EXP_MAX &&
        controls[["log_bounds", exact = TRUE]][[2L]] <= .EXP_MAX,
      "truth_bounds", paste0(control_path, ".log_bounds"),
      "bounds within the safe exponential range", controls[[
        "log_bounds", exact = TRUE
      ]]
    )
    .dpprior_schema_require(
      identical(
        raw[["computation", exact = TRUE]][["request", exact = TRUE]][[
          "controls", exact = TRUE
        ]][["log_bounds", exact = TRUE]],
        controls[["log_bounds", exact = TRUE]]
      ),
      "truth_control_identity", paste0(control_path, ".log_bounds"),
      "request and used hard log bounds identical", controls[[
        "log_bounds", exact = TRUE
      ]]
    )
    .dpprior_schema_require(
      controls[["scan_keep", exact = TRUE]] <=
        controls[["scan_points", exact = TRUE]] &&
        controls[["profile_starts", exact = TRUE]] <=
          controls[["scan_points", exact = TRUE]],
      "truth_control", control_path,
      "scan_keep/profile_starts no greater than scan_points", controls
    )
    .dpprior_schema_exact_names(
      controls[["auto_cap", exact = TRUE]], c("scan_keep", "profile_starts"),
      paste0(control_path, ".auto_cap")
    )
    for (field in c("scan_keep", "profile_starts")) {
      .dpprior_schema_validate_scalar_logical(
        controls[["auto_cap", exact = TRUE]][[field, exact = TRUE]],
        paste0(control_path, ".auto_cap.", field)
      )
    }
    for (pair in list(
      c("root_tol", 1e-10), c("optim_reltol", 1e-10),
      c("constraint_abs_tol", 1e-6), c("constraint_rel_tol", 1e-6),
      c("verification_abs_tol", 1e-8),
      c("verification_rel_tol", 1e-6),
      c("perturbation_abs_tol", 1e-4), c("penalty", 1e12)
    )) {
      validate_cap(
        controls[[pair[[1L]], exact = TRUE]], as.numeric(pair[[2L]]),
        paste0(control_path, ".", pair[[1L]])
      )
    }
    validate_cap(controls[["boundary_tol", exact = TRUE]], 1e-5,
                 paste0(control_path, ".boundary_tol"), exact = TRUE)
    validate_cap(controls[["perturbation_step", exact = TRUE]], 1e-6,
                 paste0(control_path, ".perturbation_step"), exact = TRUE)
    .dpprior_schema_exact_names(
      tolerances,
      c("constraint", "K", "weight", "perturbation", "boundary", "certificate"),
      tolerance_path
    )
    .dpprior_validate_truth_numeric_record(
      tolerances[["constraint", exact = TRUE]],
      c("absolute", "relative", "effective"),
      paste0(tolerance_path, ".constraint")
    )
    for (field in c("K", "weight")) {
      .dpprior_validate_truth_numeric_record(
        tolerances[[field, exact = TRUE]],
        c("absolute", "relative", "scale_floor"),
        paste0(tolerance_path, ".", field)
      )
    }
    .dpprior_validate_truth_numeric_record(
      tolerances[["perturbation", exact = TRUE]],
      c("absolute", "relative", "scale_floor", "step", "min_evaluations"),
      paste0(tolerance_path, ".perturbation")
    )
    .dpprior_validate_truth_numeric_record(
      tolerances[["certificate", exact = TRUE]],
      c("corner_uncertainty", "rounding_floor"),
      paste0(tolerance_path, ".certificate")
    )
    .dpprior_schema_require(
      identical(
        tolerances[["constraint", exact = TRUE]][["absolute", exact = TRUE]],
        controls[["constraint_abs_tol", exact = TRUE]]
      ) && identical(
        tolerances[["constraint", exact = TRUE]][["relative", exact = TRUE]],
        controls[["constraint_rel_tol", exact = TRUE]]
      ) && identical(
        tolerances[["constraint", exact = TRUE]][["effective", exact = TRUE]],
        controls[["constraint_abs_tol", exact = TRUE]] +
          controls[["constraint_rel_tol", exact = TRUE]] * max(
            abs(raw[["target", exact = TRUE]][["weight", exact = TRUE]][[
              "value", exact = TRUE
            ]]), 1e-8
          )
      ) && identical(tolerances[["K", exact = TRUE]],
                tolerances[["weight", exact = TRUE]]) &&
        identical(
          tolerances[["K", exact = TRUE]][["absolute", exact = TRUE]],
          controls[["verification_abs_tol", exact = TRUE]]
        ) && identical(
          tolerances[["K", exact = TRUE]][["relative", exact = TRUE]],
          controls[["verification_rel_tol", exact = TRUE]]
        ) && identical(
          tolerances[["K", exact = TRUE]][["scale_floor", exact = TRUE]], 1e-8
        ) && identical(
          tolerances[["perturbation", exact = TRUE]][["absolute", exact = TRUE]],
          controls[["perturbation_abs_tol", exact = TRUE]]
        ) && identical(
          tolerances[["perturbation", exact = TRUE]][["relative", exact = TRUE]],
          controls[["verification_rel_tol", exact = TRUE]]
        ) && identical(
          tolerances[["perturbation", exact = TRUE]][["scale_floor", exact = TRUE]],
          1e-8
        ) && identical(
          tolerances[["perturbation", exact = TRUE]][["step", exact = TRUE]],
          controls[["perturbation_step", exact = TRUE]]
        ) && identical(
          tolerances[["perturbation", exact = TRUE]][["min_evaluations", exact = TRUE]],
          2L
        ) && identical(tolerances[["boundary", exact = TRUE]],
                       controls[["boundary_tol", exact = TRUE]]) &&
        identical(
          tolerances[["certificate", exact = TRUE]][[
            "corner_uncertainty", exact = TRUE
          ]], 1e-6
        ) && identical(
          tolerances[["certificate", exact = TRUE]][[
            "rounding_floor", exact = TRUE
          ]], 64 * .Machine$double.eps
        ),
      "truth_tolerance_identity", tolerance_path,
      "hard verification, perturbation, boundary, and certificate policies",
      tolerances
    )
  }

  if (identical(mode, "dual_soft")) {
    .dpprior_schema_exact_names(
      controls,
      c(
        "max_iter", "log_bounds", "primary", "fallback", "boundary_tol",
        "verification_abs_tol", "verification_rel_tol", "objective_abs_tol",
        "objective_rel_tol", "stationarity_step", "stationarity_tol",
        "local_neighbor_step", "selection_tolerance", "optimizer_adapter"
      ), control_path
    )
    .dpprior_schema_require(
      .dpprior_schema_is_count(controls[["max_iter", exact = TRUE]], 1L),
      "truth_control", paste0(control_path, ".max_iter"),
      "a positive integer", controls[["max_iter", exact = TRUE]]
    )
    validate_bounds(controls[["log_bounds", exact = TRUE]],
                    paste0(control_path, ".log_bounds"))
    for (record_name in c("primary", "fallback")) {
      record <- controls[[record_name, exact = TRUE]]
      expected <- if (identical(record_name, "primary")) {
        c("maxit", "fnscale", "parscale", "ndeps")
      } else {
        c("maxit", "reltol", "fnscale", "parscale", "ndeps")
      }
      .dpprior_schema_exact_names(
        record, expected, paste0(control_path, ".", record_name)
      )
      .dpprior_schema_require(
        .dpprior_schema_is_count(record[["maxit", exact = TRUE]], 1L) &&
          identical(record[["fnscale", exact = TRUE]], 1),
        "truth_optimizer_controls", paste0(control_path, ".", record_name),
        "positive maxit and fnscale exactly one", record
      )
      for (field in c("parscale", "ndeps")) {
        values <- record[[field, exact = TRUE]]
        .dpprior_schema_require(
          is.numeric(values) && !is.object(values) && is.null(dim(values)) &&
            .dpprior_schema_has_only_attributes(values, "names", "names") &&
            identical(names(values), c("log_shape", "log_rate")) &&
            length(values) == 2L && !anyNA(values) && all(is.finite(values)) &&
            all(values > 0),
          "truth_optimizer_controls",
          paste0(control_path, ".", record_name, ".", field),
          "two positive named log-parameter controls", values
        )
      }
    }
    validate_cap(controls[["boundary_tol", exact = TRUE]], 1e-4,
                 paste0(control_path, ".boundary_tol"))
    validate_cap(controls[["verification_abs_tol", exact = TRUE]], 1e-8,
                 paste0(control_path, ".verification_abs_tol"))
    validate_cap(controls[["verification_rel_tol", exact = TRUE]], 1e-6,
                 paste0(control_path, ".verification_rel_tol"))
    validate_cap(controls[["objective_abs_tol", exact = TRUE]], 1e-10,
                 paste0(control_path, ".objective_abs_tol"))
    validate_cap(controls[["objective_rel_tol", exact = TRUE]], 1e-8,
                 paste0(control_path, ".objective_rel_tol"))
    validate_cap(controls[["stationarity_step", exact = TRUE]], 1e-5,
                 paste0(control_path, ".stationarity_step"), exact = TRUE)
    validate_cap(controls[["stationarity_tol", exact = TRUE]], 1e-4,
                 paste0(control_path, ".stationarity_tol"))
    validate_cap(controls[["local_neighbor_step", exact = TRUE]], 1e-3,
                 paste0(control_path, ".local_neighbor_step"), exact = TRUE)
    validate_cap(controls[["selection_tolerance", exact = TRUE]], 0,
                 paste0(control_path, ".selection_tolerance"), exact = TRUE)
    .dpprior_schema_validate_scalar_character(
      controls[["optimizer_adapter", exact = TRUE]],
      paste0(control_path, ".optimizer_adapter")
    )
    .dpprior_schema_require(
      controls[["optimizer_adapter", exact = TRUE]] %in%
        c("stats::optim", "injected_optimizer_adapter"),
      "truth_optimizer_controls", paste0(control_path, ".optimizer_adapter"),
      "stats::optim or injected_optimizer_adapter",
      controls[["optimizer_adapter", exact = TRUE]]
    )
    .dpprior_schema_exact_names(
      tolerances,
      c("K", "weight", "objective", "boundary", "stationarity",
        "neighborhood", "selection"), tolerance_path
    )
    for (field in c("K", "weight", "objective")) {
      .dpprior_validate_truth_numeric_record(
        tolerances[[field, exact = TRUE]],
        c("absolute", "relative", "scale_floor"),
        paste0(tolerance_path, ".", field)
      )
    }
    .dpprior_validate_truth_numeric_record(
      tolerances[["stationarity", exact = TRUE]], c("tolerance", "step"),
      paste0(tolerance_path, ".stationarity")
    )
    .dpprior_validate_truth_numeric_record(
      tolerances[["neighborhood", exact = TRUE]], "step",
      paste0(tolerance_path, ".neighborhood")
    )
    .dpprior_schema_require(
      identical(tolerances[["K", exact = TRUE]],
                tolerances[["weight", exact = TRUE]]) &&
        identical(
          tolerances[["K", exact = TRUE]][["absolute", exact = TRUE]],
          controls[["verification_abs_tol", exact = TRUE]]
        ) && identical(
          tolerances[["K", exact = TRUE]][["relative", exact = TRUE]],
          controls[["verification_rel_tol", exact = TRUE]]
        ) && identical(
          tolerances[["K", exact = TRUE]][["scale_floor", exact = TRUE]], 1
        ) && identical(
          tolerances[["objective", exact = TRUE]][["absolute", exact = TRUE]],
          controls[["objective_abs_tol", exact = TRUE]]
        ) && identical(
          tolerances[["objective", exact = TRUE]][["relative", exact = TRUE]],
          controls[["objective_rel_tol", exact = TRUE]]
        ) && identical(
          tolerances[["objective", exact = TRUE]][["scale_floor", exact = TRUE]],
          1
        ) && identical(tolerances[["boundary", exact = TRUE]],
                       controls[["boundary_tol", exact = TRUE]]) &&
        identical(tolerances[["stationarity", exact = TRUE]], list(
          tolerance = controls[["stationarity_tol", exact = TRUE]],
          step = controls[["stationarity_step", exact = TRUE]]
        )) && identical(tolerances[["neighborhood", exact = TRUE]], list(
          step = controls[["local_neighbor_step", exact = TRUE]]
        )) && identical(tolerances[["selection", exact = TRUE]],
                        controls[["selection_tolerance", exact = TRUE]]),
      "truth_tolerance_identity", tolerance_path,
      "soft stability/objective/KKT/neighborhood/exact-selection policies",
      tolerances
    )
  }
  invisible(TRUE)
}


.dpprior_expected_result_stability <- function(raw) {
  selected <- raw[["verification", exact = TRUE]][[
    "selected_snapshot", exact = TRUE
  ]]
  verifier <- raw[["verification", exact = TRUE]][[
    "verifier_snapshot", exact = TRUE
  ]]
  selected_K <- selected[["achieved", exact = TRUE]][["K", exact = TRUE]]
  verifier_K <- verifier[["achieved", exact = TRUE]][["K", exact = TRUE]]
  mode <- raw[["mode", exact = TRUE]]
  delta <- tolerance <- formula <- scale_floor <- numeric()
  if (identical(mode, "a2_kl")) {
    selected_pmf <- selected_K[["pmf", exact = TRUE]]
    verifier_pmf <- verifier_K[["pmf", exact = TRUE]]
    distribution_tolerance <- raw[["tolerances", exact = TRUE]][[
      "distribution", exact = TRUE
    ]][["order", exact = TRUE]][["pmf_l1", exact = TRUE]]
    .dpprior_schema_validate_finite_scalar(
      distribution_tolerance,
      "result.tolerances.distribution.order.pmf_l1", lower = 0
    )
    delta <- c(pmf.l1 = sum(abs(selected_pmf - verifier_pmf)))
    tolerance <- c(pmf.l1 = distribution_tolerance)
    formula <- c(pmf.l1 = "direct_pmf_l1_tolerance")
    scale_floor <- c(pmf.l1 = 0)
  } else {
    K_tolerance_name <- if (identical(mode, "a2_moment")) {
      "K_stability"
    } else {
      "K"
    }
    K_tolerance <- raw[["tolerances", exact = TRUE]][[
      K_tolerance_name, exact = TRUE
    ]]
    K_path <- paste0("result.tolerances.", K_tolerance_name)
    .dpprior_schema_exact_names(
      K_tolerance, c("absolute", "relative", "scale_floor"), K_path
    )
    for (field in c("absolute", "relative", "scale_floor")) {
      .dpprior_schema_validate_finite_scalar(
        K_tolerance[[field, exact = TRUE]], paste0(K_path, ".", field),
        lower = 0
      )
    }
    selected_values <- c(
      K.mean = selected_K[["mean", exact = TRUE]],
      K.variance = selected_K[["variance", exact = TRUE]]
    )
    verifier_values <- c(
      K.mean = verifier_K[["mean", exact = TRUE]],
      K.variance = verifier_K[["variance", exact = TRUE]]
    )
    delta <- abs(selected_values - verifier_values)
    K_scale_floor <- K_tolerance[["scale_floor", exact = TRUE]]
    tolerance <- K_tolerance[["absolute", exact = TRUE]] +
      K_tolerance[["relative", exact = TRUE]] * pmax(
        abs(selected_values), abs(verifier_values), K_scale_floor
      )
    formula <- setNames(
      rep("absolute_plus_relative_max", length(delta)), names(delta)
    )
    scale_floor <- setNames(rep(K_scale_floor, length(delta)), names(delta))
  }
  if (mode %in% c("dual_hard", "dual_soft")) {
    selected_weight <- selected[["achieved", exact = TRUE]][[
      "weight", exact = TRUE
    ]][["value", exact = TRUE]]
    verifier_weight <- verifier[["achieved", exact = TRUE]][[
      "weight", exact = TRUE
    ]][["value", exact = TRUE]]
    weight_tolerance <- raw[["tolerances", exact = TRUE]][[
      "weight", exact = TRUE
    ]]
    .dpprior_schema_validate_named_list(
      weight_tolerance, "result.tolerances.weight"
    )
    .dpprior_schema_require(
      "absolute" %in% names(weight_tolerance), "stability_tolerance",
      "result.tolerances.weight", "an absolute tolerance",
      names(weight_tolerance)
    )
    .dpprior_schema_validate_finite_scalar(
      weight_tolerance[["absolute", exact = TRUE]],
      "result.tolerances.weight.absolute", lower = 0
    )
    .dpprior_schema_exact_names(
      weight_tolerance, c("absolute", "relative", "scale_floor"),
      "result.tolerances.weight"
    )
    relative <- weight_tolerance[["relative", exact = TRUE]]
    .dpprior_schema_validate_finite_scalar(
      relative, "result.tolerances.weight.relative", lower = 0
    )
    weight_scale_floor <- weight_tolerance[["scale_floor", exact = TRUE]]
    .dpprior_schema_validate_finite_scalar(
      weight_scale_floor, "result.tolerances.weight.scale_floor", lower = 0
    )
    delta <- c(delta, weight.value = abs(selected_weight - verifier_weight))
    tolerance <- c(
      tolerance,
      weight.value = weight_tolerance[["absolute", exact = TRUE]] +
        relative * max(
          abs(selected_weight), abs(verifier_weight),
          weight_scale_floor
        )
    )
    formula <- c(formula, weight.value = "absolute_plus_relative_max")
    scale_floor <- c(
      scale_floor,
      weight.value = weight_scale_floor
    )
  }
  list(
    delta = delta, tolerance = tolerance, formula = formula,
    scale_floor = scale_floor
  )
}


.dpprior_result_target_moments <- function(target,
                                           path = "result.target.K") {
  implied <- target[["implied", exact = TRUE]]
  .dpprior_schema_exact_names(implied, c("mean", "variance"),
                              paste0(path, ".implied"))
  for (field in c("mean", "variance")) {
    .dpprior_schema_validate_finite_scalar(
      implied[[field, exact = TRUE]], paste0(path, ".implied.", field),
      lower = if (identical(field, "variance")) 0 else 1
    )
  }
  c(
    mean = implied[["mean", exact = TRUE]],
    variance = implied[["variance", exact = TRUE]]
  )
}


.dpprior_result_K_residuals <- function(snapshot, target_moments, path) {
  achieved <- snapshot[["achieved", exact = TRUE]][["K", exact = TRUE]]
  residuals <- snapshot[["residuals", exact = TRUE]][["K", exact = TRUE]]
  .dpprior_schema_exact_names(residuals, c("mean", "variance"), path)
  expected <- c(
    mean = achieved[["mean", exact = TRUE]] - target_moments[["mean"]],
    variance = achieved[["variance", exact = TRUE]] -
      target_moments[["variance"]]
  )
  recorded <- unlist(residuals, use.names = TRUE)
  .dpprior_schema_require(
    identical(recorded, expected), "K_residual_identity", path,
    "residuals freshly recomputed from achieved K minus target moments",
    residuals
  )
  expected
}


.dpprior_result_K_tolerance <- function(tolerances, achieved,
                                        target_moments, path) {
  .dpprior_schema_exact_names(
    tolerances, c("absolute", "relative", "scale_formula"), path
  )
  for (field in c("absolute", "relative")) {
    .dpprior_schema_validate_finite_scalar(
      tolerances[[field, exact = TRUE]], paste0(path, ".", field), lower = 0
    )
  }
  .dpprior_schema_require(
    identical(tolerances[["scale_formula", exact = TRUE]],
              "max(abs(target),1)"),
    "K_adequacy_formula", paste0(path, ".scale_formula"),
    "max(abs(target),1)", tolerances[["scale_formula", exact = TRUE]]
  )
  c(
    mean = tolerances[["absolute", exact = TRUE]] +
      tolerances[["relative", exact = TRUE]] * max(
        abs(target_moments[["mean"]]), 1
      ),
    variance = tolerances[["absolute", exact = TRUE]] +
      tolerances[["relative", exact = TRUE]] * max(
        abs(target_moments[["variance"]]), 1
      )
  )
}


.dpprior_result_distribution_metrics <- function(snapshot,
                                                  target_pmf,
                                                  target_moments,
                                                  path) {
  achieved_K <- snapshot[["achieved", exact = TRUE]][["K", exact = TRUE]]
  achieved_pmf <- achieved_K[["pmf", exact = TRUE]]
  .dpprior_schema_require(
    !is.null(achieved_pmf), "a2_kl_pmf", paste0(path, ".achieved.K.pmf"),
    "a retained achieved PMF", NULL
  )
  .dpprior_schema_require(
    all(achieved_pmf[target_pmf > 0] > 0),
    "a2_kl_support", paste0(path, ".achieved.K.pmf"),
    "positive achieved mass wherever the target PMF is positive",
    achieved_pmf
  )
  tolerances <- snapshot[["tolerances", exact = TRUE]][[
    "distribution", exact = TRUE
  ]][["adequacy", exact = TRUE]]
  .dpprior_schema_exact_names(
    tolerances,
    c(
      "kl", "l1", "mean_scaled", "variance_scaled",
      "mean_scale_formula", "variance_scale_formula"
    ),
    paste0(path, ".tolerances.distribution.adequacy")
  )
  for (field in c("kl", "l1", "mean_scaled", "variance_scaled")) {
    .dpprior_schema_validate_finite_scalar(
      tolerances[[field, exact = TRUE]],
      paste0(path, ".tolerances.distribution.adequacy.", field), lower = 0
    )
  }
  mean_residual <- achieved_K[["mean", exact = TRUE]] -
    target_moments[["mean"]]
  variance_residual <- achieved_K[["variance", exact = TRUE]] -
    target_moments[["variance"]]
  positive <- target_pmf > 0
  expected_residuals <- c(
    kl = sum(target_pmf[positive] * log(
      target_pmf[positive] / achieved_pmf[positive]
    )),
    l1 = sum(abs(target_pmf - achieved_pmf)),
    mean = mean_residual,
    variance = variance_residual
  )
  residuals <- snapshot[["residuals", exact = TRUE]][[
    "distribution", exact = TRUE
  ]]
  .dpprior_schema_exact_names(
    residuals, names(expected_residuals),
    paste0(path, ".residuals.distribution")
  )
  recorded <- unlist(residuals, use.names = TRUE)
  numerical_tolerance <- 1e-12 * pmax(
    1, abs(recorded), abs(expected_residuals)
  )
  .dpprior_schema_require(
    identical(names(recorded), names(expected_residuals)) &&
      all(abs(recorded - expected_residuals) <= numerical_tolerance),
    "a2_kl_metric_identity", paste0(path, ".residuals.distribution"),
    "KL/L1 and moment residuals freshly recomputed from the PMFs",
    list(recorded = recorded, recomputed = expected_residuals)
  )
  mean_scale <- max(1, sqrt(target_moments[["variance"]]))
  variance_scale <- max(1, target_moments[["variance"]])
  values <- c(
    kl = expected_residuals[["kl"]],
    l1 = expected_residuals[["l1"]],
    mean_scaled = abs(mean_residual) / mean_scale,
    variance_scaled = abs(variance_residual) / variance_scale
  )
  list(values = values, tolerances = c(
    kl = tolerances[["kl", exact = TRUE]],
    l1 = tolerances[["l1", exact = TRUE]],
    mean_scaled = tolerances[["mean_scaled", exact = TRUE]],
    variance_scaled = tolerances[["variance_scaled", exact = TRUE]]
  ))
}


.dpprior_result_A2_KL_objective_pmf <- function(target, J,
                                                path = "result.target.K") {
  if (!is.null(target[["pmf", exact = TRUE]])) {
    return(target[["pmf", exact = TRUE]])
  }
  evidence <- target[["derivation", exact = TRUE]][[
    "request_to_normalized", exact = TRUE
  ]][["evidence", exact = TRUE]][["A2_KL_objective", exact = TRUE]]
  .dpprior_schema_require(
    identical(target[["kind", exact = TRUE]], "moments") &&
      !is.null(evidence),
    "a2_kl_objective_target", path,
    paste(
      "an authoritative target PMF or a retained scaled-chi-square",
      "objective-distribution derivation for a moment target"
    ), evidence
  )
  .dpprior_schema_exact_names(
    evidence,
    c("method", "df", "scale", "binning", "normalization", "support",
      "retained_mass_before_normalization", "omitted_mass", "pmf",
      "mu_K_discrete", "var_K_discrete", "source"),
    paste0(path, ".derivation.request_to_normalized.evidence.A2_KL_objective")
  )
  .dpprior_schema_require(
    identical(evidence[["method", exact = TRUE]], "chisq") &&
      identical(evidence[["binning", exact = TRUE]],
                "continuity_corrected_half_integer_bins") &&
      identical(evidence[["normalization", exact = TRUE]],
                "explicit_support_conditioning") &&
      identical(evidence[["support", exact = TRUE]],
                c(lower = 1L, upper = as.integer(J))) &&
      identical(evidence[["source", exact = TRUE]],
                "A2_KL_backend_target"),
    "a2_kl_objective_method", paste0(path, ".derivation"),
    paste(
      "the closed Phase 8 scaled-chi-square method, binning, support,",
      "normalization, and source"
    ), evidence
  )
  for (field in c("df", "scale")) {
    .dpprior_schema_validate_finite_scalar(
      evidence[[field, exact = TRUE]], paste0(path, ".derivation.", field),
      lower = 0, lower_open = TRUE
    )
  }
  requested_moments <- .dpprior_result_target_moments(target, path)
  .dpprior_schema_require(
    requested_moments[["variance"]] > 0,
    "a2_kl_objective_variance", paste0(path, ".implied.variance"),
    "positive variance for a scaled-chi-square objective",
    requested_moments[["variance"]]
  )
  expected_df <- 2 * requested_moments[["mean"]]^2 /
    requested_moments[["variance"]]
  expected_scale <- requested_moments[["variance"]] /
    (2 * requested_moments[["mean"]])
  parameter_tolerance <- 1e-12 * max(
    1, abs(expected_df), abs(expected_scale),
    abs(evidence[["df", exact = TRUE]]),
    abs(evidence[["scale", exact = TRUE]])
  )
  .dpprior_schema_require(
    abs(evidence[["df", exact = TRUE]] - expected_df) <=
      parameter_tolerance &&
      abs(evidence[["scale", exact = TRUE]] - expected_scale) <=
        parameter_tolerance,
    "a2_kl_objective_parameters", paste0(path, ".derivation"),
    "scaled-chi-square df and scale derived from requested moments",
    evidence
  )
  .dpprior_schema_validate_scalar_character(
    evidence[["source", exact = TRUE]], paste0(path, ".derivation.source")
  )
  for (field in c("retained_mass_before_normalization", "omitted_mass")) {
    .dpprior_schema_validate_finite_scalar(
      evidence[[field, exact = TRUE]], paste0(path, ".derivation.", field),
      lower = 0, upper = 1
    )
  }
  pmf <- evidence[["pmf", exact = TRUE]]
  .dpprior_schema_require(
    is.numeric(pmf) && !is.object(pmf) && is.null(dim(pmf)) &&
      .dpprior_schema_has_only_attributes(pmf) &&
      length(pmf) == J && !anyNA(pmf) && all(is.finite(pmf)) &&
      all(pmf >= 0) && abs(sum(pmf) - 1) <= .TOL_PMF_SUM,
    "a2_kl_objective_pmf", paste0(path, ".derivation.pmf"),
    "a normalized nonnegative objective PMF of length J", pmf
  )
  lower_edges <- (seq_len(J) - 0.5) / evidence[["scale", exact = TRUE]]
  upper_edges <- (seq_len(J) + 0.5) / evidence[["scale", exact = TRUE]]
  raw_mass <- stats::pchisq(
    upper_edges, df = evidence[["df", exact = TRUE]]
  ) - stats::pchisq(
    lower_edges, df = evidence[["df", exact = TRUE]]
  )
  raw_mass <- pmax(raw_mass, 0)
  retained_mass <- sum(raw_mass)
  reconstructed_pmf <- raw_mass / retained_mass
  .dpprior_schema_require(
    is.finite(retained_mass) && retained_mass > 0 &&
      abs(evidence[["retained_mass_before_normalization", exact = TRUE]] -
            retained_mass) <= 1e-12 * max(1, abs(retained_mass)) &&
      abs(evidence[["omitted_mass", exact = TRUE]] -
            max(0, 1 - retained_mass)) <=
        1e-12 * max(1, abs(1 - retained_mass)) &&
      max(abs(pmf - reconstructed_pmf)) <= .TOL_PMF_SUM,
    "a2_kl_objective_reconstruction", paste0(path, ".derivation.pmf"),
    paste(
      "the independently reconstructed continuity-corrected conditioned",
      "scaled-chi-square PMF and retained/omitted masses"
    ), evidence
  )
  support <- seq_len(J)
  expected_mean <- sum(support * pmf)
  expected_variance <- sum((support - expected_mean)^2 * pmf)
  for (field in c("mu_K_discrete", "var_K_discrete")) {
    .dpprior_schema_validate_finite_scalar(
      evidence[[field, exact = TRUE]], paste0(path, ".derivation.", field),
      lower = if (identical(field, "var_K_discrete")) 0 else 1
    )
  }
  tolerance <- 1e-12 * pmax(
    1, abs(c(expected_mean, expected_variance)),
    abs(c(evidence[["mu_K_discrete", exact = TRUE]],
          evidence[["var_K_discrete", exact = TRUE]]))
  )
  .dpprior_schema_require(
    all(abs(c(evidence[["mu_K_discrete", exact = TRUE]],
                  evidence[["var_K_discrete", exact = TRUE]]) -
              c(expected_mean, expected_variance)) <= tolerance),
    "a2_kl_objective_moments", paste0(path, ".derivation"),
    "objective PMF moments recomputed from the retained PMF", evidence
  )
  pmf
}


.dpprior_candidate_attempt_stage_ok <- function(mode, method, stage) {
  allowed_stages <- switch(
    mode,
    a2_moment = switch(
      method,
      `A2-MN` = c("primary", "optimizer"),
      `A2-MN+NM` = c("fallback", "recovery"),
      scaled_log_newton = c("primary", "optimizer"),
      nelder_mead_log = c("fallback", "recovery"),
      fixed_log_parameter_grid = "initialization",
      character()
    ),
    a2_kl = switch(
      method,
      `A2-KL` = c("primary", "optimizer"),
      `L-BFGS-B` = c("primary", "optimizer"),
      nlminb = c("fallback", "optimizer", "recovery"),
      `A2-MN` = "initialization",
      A1 = "initialization",
      heuristic = "initialization",
      character()
    ),
    dual_hard = switch(
      method,
      dual_anchor_hard_inequality = c("primary", "optimizer"),
      `L-BFGS-B` = c("primary", "optimizer"),
      `K_only_L-BFGS-B` = c("primary", "optimizer"),
      constrained_profile_optimize = c("profile", "optimizer"),
      verified_candidate = "verification",
      analytic_monotonicity_feasibility_probe = "feasibility",
      deterministic_feasible_profile_scan = "profile_scan",
      `penalty_L-BFGS-B_diagnostic` = "diagnostic",
      character()
    ),
    dual_soft = switch(
      method,
      `dual-soft` = c("primary", "optimizer"),
      `L-BFGS-B` = c("primary", "optimizer"),
      `Nelder-Mead` = c("fallback", "recovery"),
      nlminb = c("fallback", "recovery"),
      character()
    ),
    character()
  )
  length(allowed_stages) > 0L && stage %in% allowed_stages
}


.dpprior_candidate_optimizer_method <- function(mode, method) {
  method %in% switch(
    mode,
    a2_moment = c("A2-MN", "A2-MN+NM", "scaled_log_newton",
                  "nelder_mead_log"),
    a2_kl = c("A2-KL", "L-BFGS-B", "nlminb"),
    dual_hard = c("dual_anchor_hard_inequality", "L-BFGS-B",
                  "K_only_L-BFGS-B", "constrained_profile_optimize"),
    dual_soft = c("dual-soft", "L-BFGS-B", "Nelder-Mead", "nlminb"),
    character()
  )
}


.dpprior_candidate_expected_recorded_kind <- function(mode, attempt,
                                                       generator) {
  method <- attempt[["method", exact = TRUE]]
  stage <- attempt[["stage", exact = TRUE]]
  switch(
    mode,
    a2_moment = "standardized_residual",
    a2_kl = if (identical(generator, "initialization")) NULL else "kl",
    dual_hard = switch(
      method,
      analytic_monotonicity_feasibility_probe = "weight_metric_probe",
      deterministic_feasible_profile_scan = "K_loss",
      `penalty_L-BFGS-B_diagnostic` = "penalized_diagnostic",
      `K_only_L-BFGS-B` = "K_loss",
      constrained_profile_optimize = "K_loss",
      `L-BFGS-B` = "K_loss",
      dual_anchor_hard_inequality = "K_loss",
      NULL
    ),
    dual_soft = if (
      method %in% c("dual-soft", "L-BFGS-B", "Nelder-Mead", "nlminb") &&
        stage %in% c("primary", "optimizer", "fallback", "recovery")
    ) "soft_tradeoff" else NULL,
    NULL
  )
}


.dpprior_candidate_generator_allowed <- function(mode, generator,
                                                  has_parent_attempt) {
  allowed <- switch(
    mode,
    a2_moment = c("direct_attempt", "initialization"),
    a2_kl = c("direct_attempt", "initialization"),
    dual_hard = c(
      "direct_attempt", "deterministic_profile_scan", "input_fit",
      "feasibility_extreme", "derived_diagnostic_attempt"
    ),
    dual_soft = "direct_attempt",
    character()
  )
  if (!generator %in% allowed) {
    return(FALSE)
  }
  if (!has_parent_attempt) {
    return(
      identical(mode, "dual_hard") &&
        generator %in% c(
          "deterministic_profile_scan", "input_fit", "feasibility_extreme"
        )
    )
  }
  !identical(generator, "input_fit")
}


.dpprior_candidate_check <- function(evaluation, name, path) {
  checks <- evaluation[["checks", exact = TRUE]]
  .dpprior_schema_require(
    name %in% names(checks), "candidate_checks", paste0(path, ".checks"),
    paste("a substantive", name, "candidate check"), names(checks)
  )
  checks[[name, exact = TRUE]]
}


.dpprior_candidate_expected_stability <- function(raw, evaluation) {
  candidate_raw <- raw
  candidate_verification <- candidate_raw[["verification", exact = TRUE]]
  candidate_verification[["selected_snapshot"]] <- evaluation[[
    "selected_snapshot", exact = TRUE
  ]]
  candidate_verification[["verifier_snapshot"]] <- evaluation[[
    "verifier_snapshot", exact = TRUE
  ]]
  candidate_raw[["verification"]] <- candidate_verification
  .dpprior_expected_result_stability(candidate_raw)
}


.dpprior_validate_result_candidate_science <- function(raw,
                                                       evaluation,
                                                       attempt = NULL,
                                                       path) {
  mode <- raw[["mode", exact = TRUE]]
  selected <- evaluation[["selected_snapshot", exact = TRUE]]
  verifier <- evaluation[["verifier_snapshot", exact = TRUE]]
  J <- raw[["J", exact = TRUE]]
  target_K <- raw[["target", exact = TRUE]][["K", exact = TRUE]]
  target_moments <- .dpprior_result_target_moments(target_K)
  selected_K <- selected[["achieved", exact = TRUE]][["K", exact = TRUE]]
  orders <- raw[["computation", exact = TRUE]][["orders", exact = TRUE]]
  selected_order <- orders[["M_selected", exact = TRUE]]
  verifier_order <- orders[["M_verification_used", exact = TRUE]]
  .dpprior_validate_achieved_K(
    selected_K, J, paste0(path, ".selected_snapshot.achieved.K")
  )
  .dpprior_schema_require(
    !is.null(selected_order) &&
      identical(selected[["M", exact = TRUE]], selected_order) &&
      identical(selected_K[["M", exact = TRUE]], selected_order),
    "candidate_selected_order", paste0(path, ".selected_snapshot"),
    paste(
      "snapshot M and achieved.K.M identical to the central",
      "computation.orders.M_selected authority"
    ),
    list(
      snapshot_M = selected[["M", exact = TRUE]],
      achieved_M = selected_K[["M", exact = TRUE]],
      authority = selected_order
    )
  )
  if (!is.null(verifier)) {
    verifier_K_for_order <- verifier[["achieved", exact = TRUE]][[
      "K", exact = TRUE
    ]]
    .dpprior_schema_require(
      !is.null(verifier_order) &&
        identical(verifier[["M", exact = TRUE]], verifier_order) &&
        identical(verifier_K_for_order[["M", exact = TRUE]], verifier_order),
      "candidate_verifier_order", paste0(path, ".verifier_snapshot"),
      paste(
        "snapshot M and achieved.K.M identical to the central",
        "computation.orders.M_verification_used authority"
      ),
      list(
        snapshot_M = verifier[["M", exact = TRUE]],
        achieved_M = verifier_K_for_order[["M", exact = TRUE]],
        authority = verifier_order
      )
    )
  }
  .dpprior_schema_require(
    identical(selected[["tolerances", exact = TRUE]],
              raw[["tolerances", exact = TRUE]]),
    "candidate_tolerances", paste0(path, ".selected_snapshot.tolerances"),
    "the central canonical result tolerances", selected[["tolerances", exact = TRUE]]
  )
  if (!is.null(verifier)) {
    verifier_K <- verifier[["achieved", exact = TRUE]][["K", exact = TRUE]]
    .dpprior_validate_achieved_K(
      verifier_K, J, paste0(path, ".verifier_snapshot.achieved.K")
    )
    .dpprior_schema_require(
      identical(verifier[["tolerances", exact = TRUE]],
                raw[["tolerances", exact = TRUE]]),
      "candidate_tolerances", paste0(path, ".verifier_snapshot.tolerances"),
      "the central canonical result tolerances",
      verifier[["tolerances", exact = TRUE]]
    )
  }

  expected_kind <- NULL
  expected_fresh <- NULL
  expected_selection <- NULL

  if (identical(mode, "a2_moment")) {
    selected_residual <- .dpprior_result_K_residuals(
      selected, target_moments, paste0(path, ".selected_snapshot.residuals.K")
    )
    selected_tolerance <- .dpprior_result_K_tolerance(
      raw[["tolerances", exact = TRUE]][["K_adequacy", exact = TRUE]],
      selected_K,
      target_moments, paste0(path, ".selected_snapshot.tolerances.K")
    )
    expected_kind <- "standardized_residual"
    expected_fresh <- sum((selected_residual / selected_tolerance)^2)
    expected_selection <- expected_fresh
    .dpprior_bind_decision_check(
      .dpprior_candidate_check(evaluation, "candidate_finite", path),
      c(parameters = TRUE, snapshot = TRUE, objective = TRUE),
      c(parameters = TRUE, snapshot = TRUE, objective = TRUE),
      NULL, "identical", paste0(path, ".checks.candidate_finite")
    )
    .dpprior_schema_require(
      identical(names(evaluation[["checks", exact = TRUE]]),
                "candidate_finite"),
      "candidate_check_set", paste0(path, ".checks"),
      "only the finite selected-order gate for A2-MN candidate selection",
      names(evaluation[["checks", exact = TRUE]])
    )
  }

  if (identical(mode, "a2_kl")) {
    target_pmf <- .dpprior_result_A2_KL_objective_pmf(target_K, J)
    target_support <- seq_len(J)
    target_objective_mean <- sum(target_support * target_pmf)
    target_objective_moments <- c(
      mean = target_objective_mean,
      variance = sum(
        (target_support - target_objective_mean)^2 * target_pmf
      )
    )
    selected_metrics <- .dpprior_result_distribution_metrics(
      selected, target_pmf, target_objective_moments,
      paste0(path, ".selected_snapshot")
    )
    expected_kind <- "kl"
    expected_fresh <- selected_metrics[["values"]][["kl"]]
    expected_selection <- expected_fresh
    .dpprior_bind_decision_check(
      .dpprior_candidate_check(evaluation, "candidate_distribution", path),
      c(
        pmf_mass_error = abs(sum(selected_K[["pmf", exact = TRUE]]) - 1),
        pmf_minimum_violation = max(
          0, -min(selected_K[["pmf", exact = TRUE]])
        )
      ),
      c(pmf_mass_error = 0, pmf_minimum_violation = 0),
      c(pmf_mass_error = .TOL_PMF_SUM, pmf_minimum_violation = 0),
      "lte", paste0(path, ".checks.candidate_distribution")
    )
    .dpprior_schema_require(
      identical(names(evaluation[["checks", exact = TRUE]]),
                "candidate_distribution"),
      "candidate_check_set", paste0(path, ".checks"),
      paste(
        "only the finite selected-order distribution gate for A2-KL",
        "candidate selection; selected-result adequacy and order checks are",
        "validated separately"
      ),
      names(evaluation[["checks", exact = TRUE]])
    )
  }

  if (identical(mode, "dual_hard")) {
    weight_target <- raw[["target", exact = TRUE]][["weight", exact = TRUE]]
    selected_weight_record <- selected[["achieved", exact = TRUE]][[
      "weight", exact = TRUE
    ]]
    .dpprior_validate_achieved_weight(
      selected_weight_record, weight_target[["metric", exact = TRUE]],
      paste0(path, ".selected_snapshot.achieved.weight")
    )
    K_scales <- raw[["constraint", exact = TRUE]][[
      "optimality", exact = TRUE
    ]][["K_scales", exact = TRUE]]
    if (is.null(K_scales)) {
      K_scales <- raw[["computation", exact = TRUE]][[
        "scaling", exact = TRUE
      ]][["values", exact = TRUE]][["K", exact = TRUE]]
    }
    .dpprior_schema_exact_names(
      K_scales, c("mean", "variance"), paste0(path, ".K_scales")
    )
    for (field in c("mean", "variance")) {
      .dpprior_schema_validate_finite_scalar(
        K_scales[[field, exact = TRUE]], paste0(path, ".K_scales.", field),
        lower = 0, lower_open = TRUE
      )
    }
    K_loss <-
      ((selected_K[["mean", exact = TRUE]] - target_moments[["mean"]]) /
         K_scales[["mean", exact = TRUE]])^2 +
      ((selected_K[["variance", exact = TRUE]] -
          target_moments[["variance"]]) /
         K_scales[["variance", exact = TRUE]])^2
    expected_kind <- "K_loss"
    expected_fresh <- K_loss
    if (is.null(verifier)) {
      expected_selection <- NULL
      .dpprior_bind_decision_check(
        .dpprior_candidate_check(evaluation, "candidate_domain", path),
        c(
          K_support = selected_K[["mean", exact = TRUE]] >= 1 &&
            selected_K[["mean", exact = TRUE]] <= J,
          weight_support = selected_weight_record[["value", exact = TRUE]] >= 0 &&
            selected_weight_record[["value", exact = TRUE]] <= 1
        ),
        c(K_support = TRUE, weight_support = TRUE), NULL, "identical",
        paste0(path, ".checks.candidate_domain")
      )
      .dpprior_schema_require(
        identical(names(evaluation[["checks", exact = TRUE]]),
                  "candidate_domain"),
        "candidate_check_set", paste0(path, ".checks"),
        paste(
          "only a finite-domain gate when no independent hard verifier was",
          "performed; the fresh K-loss is retained but not selection-eligible"
        ),
        names(evaluation[["checks", exact = TRUE]])
      )
    } else {
      expected_selection <- K_loss
      verifier_weight_record <- verifier[["achieved", exact = TRUE]][[
        "weight", exact = TRUE
      ]]
      .dpprior_validate_achieved_weight(
        verifier_weight_record, weight_target[["metric", exact = TRUE]],
        paste0(path, ".verifier_snapshot.achieved.weight")
      )
      residual_for <- function(value) {
        if (identical(weight_target[["relation", exact = TRUE]], "at_most")) {
          value - weight_target[["value", exact = TRUE]]
        } else {
          weight_target[["value", exact = TRUE]] - value
        }
      }
      constraint_tolerance <- raw[["constraint", exact = TRUE]][[
        "tolerance", exact = TRUE
      ]][["effective", exact = TRUE]]
      if (is.null(constraint_tolerance)) {
        constraint_tolerance <- raw[["tolerances", exact = TRUE]][[
          "constraint", exact = TRUE
        ]][["effective", exact = TRUE]]
      }
      .dpprior_schema_validate_finite_scalar(
        constraint_tolerance, paste0(path, ".constraint_tolerance"),
        lower = 0
      )
      .dpprior_bind_decision_check(
        .dpprior_candidate_check(evaluation, "constraint_selected", path),
        residual_for(selected_weight_record[["value", exact = TRUE]]), 0,
        constraint_tolerance, "lte",
        paste0(path, ".checks.constraint_selected")
      )
      .dpprior_bind_decision_check(
        .dpprior_candidate_check(evaluation, "constraint_refined", path),
        residual_for(verifier_weight_record[["value", exact = TRUE]]), 0,
        constraint_tolerance, "lte",
        paste0(path, ".checks.constraint_refined")
      )
      stability <- .dpprior_candidate_expected_stability(raw, evaluation)
      .dpprior_bind_decision_check(
        .dpprior_candidate_check(evaluation, "order_stability", path),
        stability[["delta"]],
        setNames(rep(0, length(stability[["delta"]])),
                 names(stability[["delta"]])),
        stability[["tolerance"]], "lte",
        paste0(path, ".checks.order_stability")
      )
      certification_identity <- c(
        selected_metric = identical(
          selected_weight_record[["metric", exact = TRUE]],
          weight_target[["metric", exact = TRUE]]
        ),
        refined_metric = identical(
          verifier_weight_record[["metric", exact = TRUE]],
          weight_target[["metric", exact = TRUE]]
        ),
        selected_finite = is.finite(
          selected_weight_record[["value", exact = TRUE]]
        ),
        refined_finite = is.finite(
          verifier_weight_record[["value", exact = TRUE]]
        )
      )
      .dpprior_bind_decision_check(
        .dpprior_candidate_check(evaluation, "metric_certification", path),
        certification_identity,
        setNames(rep(TRUE, length(certification_identity)),
                 names(certification_identity)),
        NULL, "identical", paste0(path, ".checks.metric_certification")
      )
      perturbation <- .dpprior_candidate_check(
        evaluation, "perturbation", path
      )
      perturbation_policy <- raw[["tolerances", exact = TRUE]][[
        "perturbation", exact = TRUE
      ]]
      expected_perturbation_tolerance <-
        perturbation_policy[["absolute", exact = TRUE]] +
        perturbation_policy[["relative", exact = TRUE]] * max(
          abs(selected_weight_record[["value", exact = TRUE]]),
          perturbation_policy[["scale_floor", exact = TRUE]]
        )
      .dpprior_schema_require(
        identical(perturbation[["operator", exact = TRUE]], "lte") &&
          is.numeric(perturbation[["value", exact = TRUE]]) &&
          identical(names(perturbation[["value", exact = TRUE]]),
                    "maximum_metric_delta") &&
          identical(perturbation[["reference", exact = TRUE]], 0) &&
          identical(perturbation[["tolerance", exact = TRUE]],
                    c(maximum_metric_delta =
                        expected_perturbation_tolerance)),
        "candidate_perturbation", paste0(path, ".checks.perturbation"),
        paste(
          "the retained maximum perturbation delta versus zero using the",
          "central fixed perturbation tolerance"
        ),
        perturbation
      )
      invariant_identity <- c(
        selected_finite = selected[["finite", exact = TRUE]],
        refined_finite = verifier[["finite", exact = TRUE]],
        parameter_identity = identical(
          selected[["parameters", exact = TRUE]],
          verifier[["parameters", exact = TRUE]]
        ),
        support_valid = selected_K[["mean", exact = TRUE]] >= 1 &&
          selected_K[["mean", exact = TRUE]] <= J &&
          verifier_K[["mean", exact = TRUE]] >= 1 &&
          verifier_K[["mean", exact = TRUE]] <= J
      )
      .dpprior_bind_decision_check(
        .dpprior_candidate_check(evaluation, "invariants", path),
        invariant_identity,
        setNames(rep(TRUE, length(invariant_identity)),
                 names(invariant_identity)), NULL,
        "identical", paste0(path, ".checks.invariants")
      )
      .dpprior_schema_require(
        identical(
          names(evaluation[["checks", exact = TRUE]]),
          c(
            "constraint_selected", "constraint_refined", "order_stability",
            "metric_certification", "perturbation", "invariants"
          )
        ),
        "candidate_check_set", paste0(path, ".checks"),
        paste(
          "the complete hard-candidate selected/refined constraint, order,",
          "metric, perturbation, and invariant gate set"
        ),
        names(evaluation[["checks", exact = TRUE]])
      )
      check_passes <- vapply(
        evaluation[["checks", exact = TRUE]],
        function(check) isTRUE(check[["passed", exact = TRUE]]), logical(1)
      )
      expected_diagnostic <-
        evaluation[["execution_success", exact = TRUE]] &&
        all(check_passes[c(
          "order_stability", "metric_certification", "perturbation",
          "invariants"
        )]) &&
        !all(check_passes[c("constraint_selected", "constraint_refined")])
      .dpprior_schema_require(
        identical(evaluation[["diagnostic_eligible", exact = TRUE]],
                  expected_diagnostic),
        "hard_diagnostic_eligibility", paste0(path, ".diagnostic_eligible"),
        paste(
          "TRUE exactly for verifier-complete finite K-loss candidates whose",
          "nonconstraint science passes but selected/refined constraint fails"
        ), evaluation[["diagnostic_eligible", exact = TRUE]]
      )
    }
  }

  if (identical(mode, "dual_soft")) {
    weight_target <- raw[["target", exact = TRUE]][["weight", exact = TRUE]]
    selected_weight_record <- selected[["achieved", exact = TRUE]][[
      "weight", exact = TRUE
    ]]
    .dpprior_validate_achieved_weight(
      selected_weight_record, weight_target[["metric", exact = TRUE]],
      paste0(path, ".selected_snapshot.achieved.weight")
    )
    tradeoff <- raw[["tradeoff", exact = TRUE]]
    K_scales <- tradeoff[["scales", exact = TRUE]][["K", exact = TRUE]]
    K_loss <-
      ((selected_K[["mean", exact = TRUE]] - target_moments[["mean"]]) /
         K_scales[["mean", exact = TRUE]])^2 +
      ((selected_K[["variance", exact = TRUE]] -
          target_moments[["variance"]]) /
         K_scales[["variance", exact = TRUE]])^2
    raw_residual <- selected_weight_record[["value", exact = TRUE]] -
      weight_target[["value", exact = TRUE]]
    directed_residual <- switch(
      weight_target[["relation", exact = TRUE]],
      target = raw_residual,
      at_most = max(0, raw_residual),
      at_least = max(0, -raw_residual)
    )
    weight_loss <- (directed_residual /
      tradeoff[["scales", exact = TRUE]][["weight", exact = TRUE]])^2
    expected_kind <- "soft_tradeoff"
    expected_fresh <- tradeoff[["lambda", exact = TRUE]] * K_loss +
      (1 - tradeoff[["lambda", exact = TRUE]]) * weight_loss
    expected_selection <- expected_fresh
    .dpprior_bind_decision_check(
      .dpprior_candidate_check(evaluation, "candidate_domain", path),
      c(
        K_support = selected_K[["mean", exact = TRUE]] >= 1 &&
          selected_K[["mean", exact = TRUE]] <= J,
        K_variance = selected_K[["variance", exact = TRUE]] >= 0,
        weight_support = selected_weight_record[["value", exact = TRUE]] >= 0 &&
          selected_weight_record[["value", exact = TRUE]] <= 1
      ),
      c(K_support = TRUE, K_variance = TRUE, weight_support = TRUE),
      NULL, "identical",
      paste0(path, ".checks.candidate_domain")
    )
    .dpprior_schema_require(
      identical(names(evaluation[["checks", exact = TRUE]]),
                "candidate_domain"),
      "candidate_check_set", paste0(path, ".checks"),
      paste(
        "only the finite selected-order domain gate for soft candidate",
        "selection; selected-result optimality and order checks are",
        "validated separately"
      ),
      names(evaluation[["checks", exact = TRUE]])
    )
  }

  if (!identical(mode, "dual_hard")) {
    .dpprior_schema_require(
      !evaluation[["diagnostic_eligible", exact = TRUE]],
      "candidate_diagnostic_mode", paste0(path, ".diagnostic_eligible"),
      "FALSE outside the signed hard diagnostic tier",
      evaluation[["diagnostic_eligible", exact = TRUE]]
    )
  } else if (is.null(verifier)) {
    .dpprior_schema_require(
      !evaluation[["diagnostic_eligible", exact = TRUE]],
      "hard_diagnostic_eligibility", paste0(path, ".diagnostic_eligible"),
      "FALSE without an independent hard verifier snapshot",
      evaluation[["diagnostic_eligible", exact = TRUE]]
    )
  }

  numerical_tolerance <- 1e-12 * max(
    1, abs(expected_fresh), abs(evaluation[["fresh_objective", exact = TRUE]])
  )
  selection_identity <- if (is.null(expected_selection)) {
    is.null(evaluation[["selection_objective", exact = TRUE]])
  } else {
    !is.null(evaluation[["selection_objective", exact = TRUE]]) &&
      abs(evaluation[["selection_objective", exact = TRUE]] -
            expected_selection) <= numerical_tolerance
  }
  .dpprior_schema_require(
    identical(evaluation[["objective_kind", exact = TRUE]], expected_kind) &&
      abs(evaluation[["fresh_objective", exact = TRUE]] - expected_fresh) <=
        numerical_tolerance && selection_identity,
    "candidate_objective_identity", path,
    paste(
      "mode-specific fresh objective, objective kind, and comparable",
      "selection objective recomputed from candidate snapshots"
    ),
    list(
      recorded_kind = evaluation[["objective_kind", exact = TRUE]],
      expected_kind = expected_kind,
      recorded_fresh = evaluation[["fresh_objective", exact = TRUE]],
      expected_fresh = expected_fresh,
      recorded_selection = evaluation[["selection_objective", exact = TRUE]],
      expected_selection = expected_selection
    )
  )
  invisible(TRUE)
}


.dpprior_attempt_exit_zero <- function(attempt) {
  identical(attempt[["exit_code", exact = TRUE]], 0L) ||
    identical(attempt[["exit_code", exact = TRUE]], 0)
}


.dpprior_result_selection_tolerance <- function(raw) {
  mode <- raw[["mode", exact = TRUE]]
  value <- switch(
    mode,
    a2_moment = raw[["computation", exact = TRUE]][[
      "used", exact = TRUE
    ]][["controls", exact = TRUE]][["selection_tolerance", exact = TRUE]],
    a2_kl = raw[["computation", exact = TRUE]][[
      "used", exact = TRUE
    ]][["controls", exact = TRUE]][["selection_tolerance", exact = TRUE]],
    dual_hard = raw[["constraint", exact = TRUE]][[
      "optimality", exact = TRUE
    ]][["tie_tolerance", exact = TRUE]],
    dual_soft = raw[["tradeoff", exact = TRUE]][[
      "optimality", exact = TRUE
    ]][["selection_tolerance", exact = TRUE]],
    NULL
  )
  if (!is.null(value)) {
    .dpprior_schema_validate_finite_scalar(
      value, "result.candidate_selection_tolerance", lower = 0
    )
  }
  value
}


.dpprior_validate_result_candidate_ledger <- function(raw) {
  mode <- raw[["mode", exact = TRUE]]
  status <- raw[["status", exact = TRUE]]
  computation <- raw[["computation", exact = TRUE]]
  verification <- raw[["verification", exact = TRUE]]
  attempts <- computation[["attempts", exact = TRUE]]
  evaluations <- computation[["candidate_evaluations", exact = TRUE]]
  attempt_ids <- vapply(
    attempts, function(attempt) attempt[["id", exact = TRUE]], character(1)
  )
  selected_attempt_flags <- vapply(
    attempts, function(attempt) isTRUE(attempt[["selected", exact = TRUE]]),
    logical(1)
  )
  evaluation_ids <- vapply(
    evaluations,
    function(evaluation) evaluation[["id", exact = TRUE]], character(1)
  )
  selected_evaluation_flags <- vapply(
    evaluations,
    function(evaluation) isTRUE(evaluation[["selected", exact = TRUE]]),
    logical(1)
  )
  candidate_modes <- c("a2_moment", "a2_kl", "dual_hard", "dual_soft")
  finite_fit <- identical(raw[["object_type", exact = TRUE]], "fit") &&
    !is.null(raw[["parameters", exact = TRUE]])
  soft_endpoint <- finite_fit && identical(mode, "dual_soft") &&
    isTRUE(raw[["tradeoff", exact = TRUE]][["endpoint", exact = TRUE]])

  if (!(identical(raw[["object_type", exact = TRUE]], "fit") &&
        mode %in% candidate_modes)) {
    .dpprior_schema_require(
      length(evaluations) == 0L &&
        is.null(computation[["selected_candidate_id", exact = TRUE]]),
      "candidate_ledger_mode", "result.computation.candidate_evaluations",
      "an empty candidate ledger outside optimizer/dual fit modes",
      evaluation_ids
    )
    return(invisible(TRUE))
  }

  for (i in seq_along(attempts)) {
    attempt <- attempts[[i]]
    .dpprior_schema_require(
      .dpprior_candidate_attempt_stage_ok(
        mode, attempt[["method", exact = TRUE]],
        attempt[["stage", exact = TRUE]]
      ),
      "attempt_stage_method",
      sprintf("result.computation.attempts[[%d]]", i),
      "a producer-authorized stage/method pair for the declared mode",
      attempt[c("stage", "method")]
    )
    if (!selected_attempt_flags[[i]]) {
      .dpprior_schema_require(
        !identical(attempt[["reason_code", exact = TRUE]], "selected"),
        "attempt_reason",
        sprintf("result.computation.attempts[[%d]].reason_code", i),
        "a non-selected reason code when selected=FALSE",
        attempt[["reason_code", exact = TRUE]]
      )
    }
    reason <- attempt[["reason_code", exact = TRUE]]
    has_parameters <- !is.null(attempt[["candidate_parameters", exact = TRUE]])
    has_objective <- !is.null(attempt[["candidate_objective", exact = TRUE]])
    has_error <- !is.null(attempt[["error", exact = TRUE]])
    zero_exit <- .dpprior_attempt_exit_zero(attempt)
    certificate_probe <- identical(mode, "dual_hard") &&
      identical(status, "infeasible") &&
      identical(attempt[["stage", exact = TRUE]], "feasibility") &&
      identical(
        attempt[["method", exact = TRUE]],
        "analytic_monotonicity_feasibility_probe"
      ) && identical(reason, "globally_infeasible_by_certificate")
    aggregate_profile_scan <- identical(mode, "dual_hard") &&
      identical(attempt[["stage", exact = TRUE]], "profile_scan") &&
      identical(
        attempt[["method", exact = TRUE]],
        "deterministic_feasible_profile_scan"
      ) && identical(reason, "diagnostic_only") && !has_parameters &&
      !has_objective
    .dpprior_schema_require(
      !has_objective || has_parameters || certificate_probe,
      "attempt_candidate_coherence",
      sprintf("result.computation.attempts[[%d]]", i),
      paste(
        "a finite objective only with retained candidate parameters, except",
        "for the exact analytic certificate probe whose objective is the",
        "active refined corner"
      ), attempt
    )
    if (identical(reason, "optimizer_error")) {
      .dpprior_schema_require(
        has_error, "attempt_reason_evidence",
        sprintf("result.computation.attempts[[%d]].reason_code", i),
        "optimizer_error only with a retained typed error", reason
      )
    }
    if (identical(reason, "optimizer_exit_nonzero")) {
      .dpprior_schema_require(
        !is.null(attempt[["exit_code", exact = TRUE]]) && !zero_exit,
        "attempt_reason_evidence",
        sprintf("result.computation.attempts[[%d]].reason_code", i),
        "optimizer_exit_nonzero only with a retained nonzero exit", reason
      )
    }
    if (identical(reason, "nonfinite_candidate")) {
      .dpprior_schema_require(
        !has_parameters, "attempt_reason_evidence",
        sprintf("result.computation.attempts[[%d]].reason_code", i),
        "nonfinite_candidate only when canonical parameters are unavailable",
        reason
      )
    }
    if (identical(reason, "nonfinite_objective")) {
      .dpprior_schema_require(
        !has_objective, "attempt_reason_evidence",
        sprintf("result.computation.attempts[[%d]].reason_code", i),
        "nonfinite_objective only when canonical objective is unavailable",
        reason
      )
    }
    if (!certificate_probe && !aggregate_profile_scan && zero_exit &&
        !has_error && !has_parameters) {
      .dpprior_schema_require(
        !has_objective && reason %in% c(
          "nonfinite_candidate", "candidate_evaluation_failed", "no_candidate"
        ),
        "attempt_reason_evidence",
        sprintf("result.computation.attempts[[%d]]", i),
        paste(
          "an exit-zero attempt without parameters must retain no objective",
          "and a closed candidate-unavailability reason"
        ), attempt
      )
    }
  }

  for (i in seq_along(evaluations)) {
    evaluation <- evaluations[[i]]
    path <- sprintf("result.computation.candidate_evaluations[[%d]]", i)
    attempt_id <- evaluation[["attempt_id", exact = TRUE]]
    attempt_index <- if (is.null(attempt_id)) NA_integer_ else
      match(attempt_id, attempt_ids)
    .dpprior_schema_require(
      is.null(attempt_id) || !is.na(attempt_index),
      "candidate_attempt_id", paste0(path, ".attempt_id"),
      "NULL or an exact ID from computation.attempts", attempt_id
    )
    attempt <- if (is.na(attempt_index)) NULL else attempts[[attempt_index]]
    generator <- evaluation[["generator", exact = TRUE]]
    .dpprior_schema_require(
      .dpprior_candidate_generator_allowed(
        mode, generator, !is.null(attempt)
      ),
      "candidate_generator_mode", paste0(path, ".generator"),
      paste(
        "a mode-authorized generator and parent-attempt presence: A2 modes",
        "and soft interior require attempts; hard alone may retain input-fit,",
        "feasibility-extreme, or deterministic generated candidates"
      ),
      list(mode = mode, generator = generator, attempt_id = attempt_id)
    )

    if (!is.null(attempt)) {
      expected_generator <- switch(
        mode,
        a2_moment = if (identical(attempt[["stage", exact = TRUE]],
                                  "initialization")) {
          "initialization"
        } else {
          "direct_attempt"
        },
        a2_kl = if (identical(attempt[["stage", exact = TRUE]],
                              "initialization")) {
          "initialization"
        } else {
          "direct_attempt"
        },
        dual_hard = switch(
          attempt[["method", exact = TRUE]],
          analytic_monotonicity_feasibility_probe = "feasibility_extreme",
          deterministic_feasible_profile_scan = "deterministic_profile_scan",
          `penalty_L-BFGS-B_diagnostic` = "derived_diagnostic_attempt",
          "direct_attempt"
        ),
        dual_soft = "direct_attempt"
      )
      .dpprior_schema_require(
        identical(generator, expected_generator),
        "candidate_generator_parent", paste0(path, ".generator"),
        "the exact mode/stage/method-specific parent-attempt generator",
        list(generator = generator, expected = expected_generator,
             attempt = attempt[c("stage", "method")])
      )
      if (generator %in% c("direct_attempt", "initialization")) {
        .dpprior_schema_require(
          identical(evaluation[["method", exact = TRUE]],
                    attempt[["method", exact = TRUE]]) &&
            identical(evaluation[["parameters", exact = TRUE]],
                      attempt[["candidate_parameters", exact = TRUE]]),
          "candidate_attempt_identity", path,
          "method and parameters identical to the owning direct attempt",
          list(evaluation = evaluation[c("method", "parameters")],
               attempt = attempt[c("method", "candidate_parameters")])
        )
        recorded_available <- !is.null(
          attempt[["candidate_objective", exact = TRUE]]
        )
        expected_recorded_kind <- .dpprior_candidate_expected_recorded_kind(
          mode, attempt, generator
        )
        .dpprior_schema_require(
          if (recorded_available) {
            !is.null(expected_recorded_kind) && identical(
              evaluation[["recorded_objective_kind", exact = TRUE]],
              expected_recorded_kind
            )
          } else {
            is.null(evaluation[["recorded_objective_kind", exact = TRUE]])
          },
          "candidate_recorded_objective_kind",
          paste0(path, ".recorded_objective_kind"),
          paste(
            "the exact mode/stage/method recorded objective kind, with",
            "A2-KL initialization objectives explicitly unavailable"
          ),
          list(
            recorded = evaluation[["recorded_objective_kind", exact = TRUE]],
            expected = expected_recorded_kind,
            attempt = attempt[c("stage", "method")]
          )
        )
        expected_recorded_reason <- if (recorded_available && !identical(
          evaluation[["recorded_objective_kind", exact = TRUE]],
          evaluation[["objective_kind", exact = TRUE]]
        )) {
          paste0(
            "source_objective_kind_mismatch:",
            evaluation[["recorded_objective_kind", exact = TRUE]], "->",
            evaluation[["objective_kind", exact = TRUE]]
          )
        } else if (recorded_available) {
          NULL
        } else {
          evaluation[["recorded_objective_reason", exact = TRUE]]
        }
        .dpprior_schema_require(
          identical(
            evaluation[["recorded_objective_available", exact = TRUE]],
            recorded_available
          ) && if (recorded_available) {
            identical(evaluation[["recorded_objective", exact = TRUE]],
                      attempt[["candidate_objective", exact = TRUE]]) &&
              identical(
                evaluation[["recorded_objective_reason", exact = TRUE]],
                expected_recorded_reason
              )
          } else {
            is.null(evaluation[["recorded_objective", exact = TRUE]]) &&
              !is.null(evaluation[["recorded_objective_reason", exact = TRUE]])
          },
          "candidate_recorded_objective", path,
          "recorded-objective availability/value identical to the attempt",
          list(evaluation = evaluation[c(
            "recorded_objective", "recorded_objective_available",
            "recorded_objective_reason"
          )], attempt = attempt[["candidate_objective", exact = TRUE]])
        )
      } else if (identical(generator, "deterministic_profile_scan")) {
        .dpprior_schema_require(
          identical(attempt[["method", exact = TRUE]],
                    "deterministic_feasible_profile_scan") &&
            identical(evaluation[["method", exact = TRUE]],
                      "deterministic_profile_scan") &&
            !evaluation[["recorded_objective_available", exact = TRUE]],
          "candidate_generator_identity", path,
          paste(
            "an aggregate deterministic-profile parent attempt and an",
            "explicitly unavailable per-candidate recorded objective"
          ),
          list(evaluation = evaluation, attempt = attempt)
        )
      } else if (identical(generator, "feasibility_extreme")) {
        recorded_available <- !is.null(
          attempt[["candidate_objective", exact = TRUE]]
        )
        .dpprior_schema_require(
          identical(mode, "dual_hard") &&
            identical(
              attempt[["method", exact = TRUE]],
              "analytic_monotonicity_feasibility_probe"
            ) && identical(attempt[["stage", exact = TRUE]], "feasibility") &&
            evaluation[["method", exact = TRUE]] %in% c(
              "analytic_monotonicity_feasibility_probe",
              "analytic_feasibility_extreme"
            ) && identical(
              evaluation[["parameters", exact = TRUE]],
              attempt[["candidate_parameters", exact = TRUE]]
            ) && identical(
              evaluation[["recorded_objective_available", exact = TRUE]],
              recorded_available
            ) && if (recorded_available) {
              identical(
                evaluation[["recorded_objective", exact = TRUE]],
                attempt[["candidate_objective", exact = TRUE]]
              ) && identical(
                evaluation[["recorded_objective_kind", exact = TRUE]],
                "weight_metric_probe"
              ) && identical(
                evaluation[["objective_kind", exact = TRUE]], "K_loss"
              )
            } else {
              is.null(evaluation[["recorded_objective", exact = TRUE]]) &&
                !is.null(
                  evaluation[["recorded_objective_reason", exact = TRUE]]
                )
            },
          "candidate_generator_identity", path,
          paste(
            "a hard analytic-feasibility parent with exact retained",
            "parameters and recorded-objective availability/value"
          ),
          list(evaluation = evaluation, attempt = attempt)
        )
      } else if (identical(generator, "derived_diagnostic_attempt")) {
        .dpprior_schema_require(
          identical(mode, "dual_hard") &&
            identical(attempt[["method", exact = TRUE]],
                      "penalty_L-BFGS-B_diagnostic") &&
            identical(attempt[["stage", exact = TRUE]], "diagnostic") &&
            identical(evaluation[["method", exact = TRUE]],
                      "penalty_L-BFGS-B_diagnostic") &&
            identical(evaluation[["parameters", exact = TRUE]],
                      attempt[["candidate_parameters", exact = TRUE]]) &&
            identical(evaluation[["recorded_objective", exact = TRUE]],
                      attempt[["candidate_objective", exact = TRUE]]) &&
            identical(evaluation[["recorded_objective_kind", exact = TRUE]],
                      "penalized_diagnostic") &&
            identical(evaluation[["objective_kind", exact = TRUE]], "K_loss"),
          "candidate_generator_identity", path,
          paste(
            "a hard penalty diagnostic parent whose raw penalized objective",
            "is retained separately from the fresh comparable K-loss"
          ), list(evaluation = evaluation, attempt = attempt)
        )
      }
    } else {
      expected_attempt_free_method <- switch(
        generator,
        deterministic_profile_scan = "deterministic_profile_scan",
        input_fit = "input_fit",
        feasibility_extreme = "analytic_feasibility_extreme",
        NA_character_
      )
      .dpprior_schema_require(
        !is.na(expected_attempt_free_method) &&
          identical(evaluation[["method", exact = TRUE]],
                    expected_attempt_free_method) &&
          !evaluation[["recorded_objective_available", exact = TRUE]],
        "candidate_generator_identity", path,
        "an authorized attempt-free generator/method with no recorded objective",
        evaluation[c("generator", "method", "recorded_objective_available")]
      )
    }

    expected_execution <- if (is.null(attempt)) {
      TRUE
    } else if (identical(generator, "initialization")) {
      is.null(attempt[["error", exact = TRUE]]) &&
        !is.null(attempt[["candidate_parameters", exact = TRUE]])
    } else {
      .dpprior_attempt_exit_zero(attempt) &&
        is.null(attempt[["error", exact = TRUE]])
    }
    expected_optimizer <- !is.null(attempt) &&
      identical(generator, "direct_attempt") && expected_execution &&
      .dpprior_candidate_optimizer_method(
        mode, evaluation[["method", exact = TRUE]]
      ) && .dpprior_candidate_attempt_stage_ok(
        mode, attempt[["method", exact = TRUE]],
        attempt[["stage", exact = TRUE]]
      ) && attempt[["stage", exact = TRUE]] %in%
        c("primary", "optimizer", "fallback", "recovery")
    .dpprior_schema_require(
      identical(evaluation[["execution_success", exact = TRUE]],
                expected_execution) &&
        identical(evaluation[["optimizer_supported", exact = TRUE]],
                  expected_optimizer),
      "candidate_execution_truth", path,
      paste(
        "execution success and ordinary optimizer support recomputed from",
        "generator, exact attempt exit/error, and approved stage/method"
      ),
      evaluation[c("execution_success", "optimizer_supported")]
    )

    .dpprior_validate_result_candidate_science(
      raw, evaluation, attempt, path
    )

    if (!is.null(attempt) &&
        generator %in% c("direct_attempt", "initialization")) {
      expected_reason <- if (evaluation[["selected", exact = TRUE]]) {
        "selected"
      } else if (evaluation[["selection_eligible", exact = TRUE]]) {
        "eligible_not_selected"
      } else {
        NULL
      }
      .dpprior_schema_require(
        if (is.null(expected_reason)) {
          attempt[["reason_code", exact = TRUE]] %in% c(
            .DPPRIOR_ATTEMPT_TERMINAL_FAILURE_REASONS,
            .DPPRIOR_ATTEMPT_SCIENTIFIC_REJECTION_REASONS
          )
        } else {
          identical(attempt[["reason_code", exact = TRUE]], expected_reason)
        },
        "candidate_attempt_outcome",
        sprintf("result.computation.attempts[[%d]].reason_code", attempt_index),
        "a reason code derived from candidate selection eligibility/outcome",
        attempt[["reason_code", exact = TRUE]]
      )
    }
  }

  attempts_with_candidates <- which(vapply(
    attempts,
    function(attempt) !is.null(attempt[["candidate_parameters", exact = TRUE]]),
    logical(1)
  ))
  for (attempt_index in attempts_with_candidates) {
    .dpprior_schema_require(
      any(vapply(
        evaluations,
        function(evaluation) identical(
          evaluation[["attempt_id", exact = TRUE]], attempt_ids[[attempt_index]]
        ),
        logical(1)
      )),
      "candidate_attempt_coverage",
      sprintf("result.computation.attempts[[%d]]", attempt_index),
      "at least one candidate-evaluation ledger entry for retained parameters",
      attempt_ids[[attempt_index]]
    )
  }

  if (finite_fit && !soft_endpoint) {
    .dpprior_schema_require(
      sum(selected_evaluation_flags) == 1L,
      "finite_fit_candidate", "result.computation.candidate_evaluations",
      "exactly one selected candidate evaluation for a finite fit",
      evaluation_ids[selected_evaluation_flags]
    )
    selected_index <- which(selected_evaluation_flags)
    selected_evaluation <- evaluations[[selected_index]]
    selected_attempt_id <- selected_evaluation[["attempt_id", exact = TRUE]]
    selected_generator <- selected_evaluation[["generator", exact = TRUE]]
    owns_selected_attempt <- !is.null(selected_attempt_id) &&
      selected_generator %in% c("direct_attempt", "initialization")
    .dpprior_schema_require(
      identical(selected_evaluation[["parameters", exact = TRUE]],
                raw[["parameters", exact = TRUE]]) &&
        identical(selected_evaluation[["selected_snapshot", exact = TRUE]],
                  verification[["selected_snapshot", exact = TRUE]]) &&
        (is.null(selected_evaluation[["verifier_snapshot", exact = TRUE]]) ||
           identical(
             selected_evaluation[["verifier_snapshot", exact = TRUE]],
             verification[["verifier_snapshot", exact = TRUE]]
           )),
      "selected_candidate_identity",
      "result.computation.candidate_evaluations.selected",
      paste(
        "public parameters/selected snapshot and any retained verifier",
        "snapshot identical to the selected ledger candidate"
      ),
      selected_evaluation
    )
    .dpprior_schema_require(
      if (owns_selected_attempt) {
        identical(computation[["selected_attempt_id", exact = TRUE]],
                  selected_attempt_id) &&
          sum(selected_attempt_flags) == 1L &&
          identical(attempt_ids[selected_attempt_flags], selected_attempt_id)
      } else {
        is.null(computation[["selected_attempt_id", exact = TRUE]]) &&
          sum(selected_attempt_flags) == 0L
      },
      "selected_candidate_attempt",
      "result.computation.selected_attempt_id",
      paste(
        "the owning direct/initialization attempt, or NULL for an",
        "attempt-free/generated selected candidate"
      ),
      computation[["selected_attempt_id", exact = TRUE]]
    )
    if (owns_selected_attempt) {
      selected_attempt <- attempts[[match(selected_attempt_id, attempt_ids)]]
      fallback_stage <- selected_attempt[["stage", exact = TRUE]] %in%
        c("fallback", "recovery")
      .dpprior_schema_require(
        identical(fallback_stage,
                  computation[["fallback", exact = TRUE]][["used", exact = TRUE]]),
        "selected_fallback_identity",
        "result.computation.candidate_evaluations.selected",
        paste(
          "selected fallback/recovery stage exactly when fallback.used=TRUE;",
          "primary/optimizer stage otherwise"
        ),
        list(
          stage = selected_attempt[["stage", exact = TRUE]],
          fallback = computation[["fallback", exact = TRUE]]
        )
      )
    } else {
      .dpprior_schema_require(
        !computation[["fallback", exact = TRUE]][["used", exact = TRUE]],
        "selected_fallback_identity",
        "result.computation.candidate_evaluations.selected",
        "fallback.used=FALSE for an attempt-free/generated selected candidate",
        computation[["fallback", exact = TRUE]]
      )
    }

    selection_pool <- which(vapply(
      evaluations,
      function(evaluation) {
        isTRUE(evaluation[["selection_eligible", exact = TRUE]]) &&
          identical(evaluation[["selection_objective_kind", exact = TRUE]],
                    selected_evaluation[[
                      "selection_objective_kind", exact = TRUE
                    ]])
      },
      logical(1)
    ))
    diagnostic_pool <- if (identical(mode, "dual_hard")) {
      which(vapply(
        evaluations,
        function(evaluation) {
          isTRUE(evaluation[["diagnostic_eligible", exact = TRUE]]) &&
            identical(evaluation[["selection_objective_kind", exact = TRUE]],
                      "K_loss")
        }, logical(1)
      ))
    } else {
      integer()
    }
    pool <- if (identical(mode, "dual_hard") &&
                length(selection_pool) == 0L) {
      diagnostic_pool
    } else {
      selection_pool
    }
    .dpprior_schema_require(
      selected_index %in% pool, "candidate_selection",
      "result.computation.selected_candidate_id",
      paste(
        "a selected feasible candidate, or only when that tier is empty,",
        "a signed hard unsatisfied diagnostic candidate"
      ),
      selected_evaluation[["id", exact = TRUE]]
    )
    selection_tolerance <- .dpprior_result_selection_tolerance(raw)
    .dpprior_schema_require(
      !is.null(selection_tolerance), "candidate_selection_tolerance",
      "result.computation.selected_candidate_id",
      "a fixed central selection tolerance", NULL
    )
    objectives <- vapply(
      evaluations[pool],
      function(evaluation) evaluation[["selection_objective", exact = TRUE]],
      numeric(1)
    )
    minimum_objective <- min(objectives)
    selected_objective <- selected_evaluation[[
      "selection_objective", exact = TRUE
    ]]
    if (identical(mode, "dual_hard")) {
      expected_tie_tolerance <- max(
        raw[["tolerances", exact = TRUE]][["K", exact = TRUE]][[
          "absolute", exact = TRUE
        ]],
        64 * .Machine$double.eps * max(1, abs(minimum_objective))
      )
      .dpprior_schema_require(
        identical(selection_tolerance, expected_tie_tolerance),
        "hard_tie_tolerance", "result.constraint.optimality.tie_tolerance",
        paste(
          "max(K verification absolute tolerance, 64 machine eps times",
          "the minimum K-loss scale)"
        ), selection_tolerance
      )
    }
    .dpprior_schema_require(
      selected_objective <= minimum_objective + selection_tolerance,
      "candidate_selection", "result.computation.selected_candidate_id",
      "the minimum comparable eligible objective within fixed tolerance",
      c(selected = selected_objective, minimum = minimum_objective,
        tolerance = selection_tolerance)
    )
    tied_pool <- pool[objectives <= minimum_objective + selection_tolerance]
    if (identical(mode, "dual_hard")) {
      supported_ties <- tied_pool[vapply(
        evaluations[tied_pool],
        function(evaluation) isTRUE(
          evaluation[["optimizer_supported", exact = TRUE]]
        ), logical(1)
      )]
      if (length(supported_ties) > 0L) {
        tied_pool <- supported_ties
      }
      .dpprior_schema_require(
        identical(selected_index, tied_pool[[1L]]),
        "hard_candidate_tie_break",
        "result.computation.selected_candidate_id",
        paste(
          "the first generated minimum-K-loss candidate after preferring",
          "ordinary optimizer support within the fixed numerical tie"
        ), selected_evaluation[["id", exact = TRUE]]
      )
    }
    if ("candidate_selection" %in%
        names(verification[["components", exact = TRUE]])) {
      .dpprior_bind_decision_check(
        verification[["components", exact = TRUE]][[
          "candidate_selection", exact = TRUE
        ]],
        selected_objective - minimum_objective, 0, selection_tolerance,
        "lte", "result.verification.components.candidate_selection"
      )
    }
    if (identical(mode, "dual_hard")) {
      optimality <- raw[["constraint", exact = TRUE]][[
        "optimality", exact = TRUE
      ]]
      .dpprior_schema_require(
        identical(selected_evaluation[["objective_kind", exact = TRUE]],
                  "K_loss") &&
          identical(optimality[["selection_rule", exact = TRUE]],
                    "minimum_K_loss") &&
          identical(optimality[["selected_K_loss", exact = TRUE]],
                    selected_objective) &&
          identical(optimality[["minimum_K_loss", exact = TRUE]],
                    minimum_objective) &&
          identical(optimality[["tie_tolerance", exact = TRUE]],
                    selection_tolerance),
        "hard_candidate_selection", "result.constraint.optimality",
        "the selected/minimum K-loss ledger proof and fixed tie tolerance",
        optimality
      )
    }
    if (status %in% c("converged", "boundary")) {
      .dpprior_schema_require(
        selected_evaluation[["optimizer_supported", exact = TRUE]] &&
          selected_evaluation[["decision_eligible", exact = TRUE]],
        "ordinary_convergence_evidence",
        "result.computation.candidate_evaluations.selected",
        paste(
          "ordinary optimizer support plus the mode-specific candidate",
          "selection gates for a converged/boundary result"
        ),
        selected_evaluation
      )
    }
  } else if (soft_endpoint) {
    .dpprior_schema_require(
      length(attempts) == 0L && length(evaluations) == 0L &&
        is.null(computation[["selected_attempt_id", exact = TRUE]]) &&
        is.null(computation[["selected_candidate_id", exact = TRUE]]),
      "soft_endpoint_selection", "result.computation",
      "no optimizer attempts/candidate ledger for a fixed input-fit endpoint",
      computation[c("attempts", "candidate_evaluations",
                    "selected_attempt_id", "selected_candidate_id")]
    )
  } else {
    .dpprior_schema_require(
      sum(selected_evaluation_flags) == 0L &&
        !any(vapply(
          evaluations,
          function(evaluation) isTRUE(
            evaluation[["selection_eligible", exact = TRUE]]
          ),
          logical(1)
        )) && !any(vapply(
          evaluations,
          function(evaluation) isTRUE(
            evaluation[["diagnostic_eligible", exact = TRUE]]
          ), logical(1)
        )) && sum(selected_attempt_flags) == 0L &&
        is.null(computation[["selected_attempt_id", exact = TRUE]]) &&
        is.null(computation[["selected_candidate_id", exact = TRUE]]),
      "no_public_candidate", "result.computation",
      paste(
        "no selected, selection-eligible, or signed diagnostic-eligible",
        "candidate when no finite public fit candidate is retained"
      ),
      list(attempts = attempt_ids[selected_attempt_flags],
           candidates = evaluation_ids[selected_evaluation_flags])
    )
  }

  if (identical(mode, "dual_hard")) {
    hard_candidates <- vapply(
      evaluations,
      function(evaluation) identical(
        evaluation[["objective_kind", exact = TRUE]], "K_loss"
      ),
      logical(1)
    )
    verified_candidates <- hard_candidates & vapply(
      evaluations,
      function(evaluation) !is.null(
        evaluation[["verifier_snapshot", exact = TRUE]]
      ),
      logical(1)
    )
    feasible_candidates <- hard_candidates & vapply(
      evaluations,
      function(evaluation) isTRUE(
        evaluation[["selection_eligible", exact = TRUE]]
      ),
      logical(1)
    )
    feasibility <- raw[["constraint", exact = TRUE]][[
      "feasibility", exact = TRUE
    ]]
    .dpprior_schema_require(
      identical(feasibility[["candidate_count", exact = TRUE]],
                as.integer(sum(hard_candidates))) &&
        identical(feasibility[["verified_candidate_count", exact = TRUE]],
                  as.integer(sum(verified_candidates))) &&
        identical(feasibility[["feasible_candidate_count", exact = TRUE]],
                  as.integer(sum(feasible_candidates))),
      "hard_candidate_counts", "result.constraint.feasibility",
      paste(
        "finite K-loss ledger rows and the subset on which the independent",
        "verifier was actually invoked (regardless of pass/fail), and the",
        "verified scientifically feasible selection subset"
      ),
      feasibility[c(
        "candidate_count", "verified_candidate_count", "feasible_candidate_count"
      )]
    )
  }

  if (length(evaluations) > 1L) {
    for (left in seq_len(length(evaluations) - 1L)) {
      for (right in seq.int(left + 1L, length(evaluations))) {
        lhs <- evaluations[[left]]
        rhs <- evaluations[[right]]
        same_parameters <- identical(
          lhs[["parameters", exact = TRUE]], rhs[["parameters", exact = TRUE]]
        )
        same_selection_kind <-
          !is.null(lhs[["selection_objective_kind", exact = TRUE]]) &&
          identical(
            lhs[["selection_objective_kind", exact = TRUE]],
            rhs[["selection_objective_kind", exact = TRUE]]
          )
        if (same_parameters && same_selection_kind) {
          tolerance <- max(
            lhs[["objective_tolerance", exact = TRUE]],
            rhs[["objective_tolerance", exact = TRUE]]
          )
          .dpprior_schema_require(
            abs(lhs[["fresh_objective", exact = TRUE]] -
                  rhs[["fresh_objective", exact = TRUE]]) <= tolerance,
            "candidate_objective_determinism",
            sprintf("result.computation.candidate_evaluations[[%d]]", right),
            paste(
              "one deterministic fresh selection objective for identical",
              "parameters and selection-objective kind"
            ),
            c(left = lhs[["fresh_objective", exact = TRUE]],
              right = rhs[["fresh_objective", exact = TRUE]])
          )
        }
        same_recorded_kind <-
          lhs[["recorded_objective_available", exact = TRUE]] &&
          rhs[["recorded_objective_available", exact = TRUE]] &&
          identical(
            lhs[["recorded_objective_kind", exact = TRUE]],
            rhs[["recorded_objective_kind", exact = TRUE]]
          )
        if (same_parameters && same_recorded_kind) {
          recorded_tolerance <- 1e-12 * max(
            1, abs(lhs[["recorded_objective", exact = TRUE]]),
            abs(rhs[["recorded_objective", exact = TRUE]])
          )
          .dpprior_schema_require(
            abs(lhs[["recorded_objective", exact = TRUE]] -
                  rhs[["recorded_objective", exact = TRUE]]) <=
              recorded_tolerance,
            "candidate_recorded_objective_determinism",
            sprintf(
              "result.computation.candidate_evaluations[[%d]]", right
            ),
            paste(
              "one recorded deterministic objective for identical",
              "parameters and recorded-objective kind"
            ),
            c(
              left = lhs[["recorded_objective", exact = TRUE]],
              right = rhs[["recorded_objective", exact = TRUE]]
            )
          )
        }
      }
    }
  }

  termination <- computation[["termination", exact = TRUE]]
  boundary_reason <- termination[["boundary_reason", exact = TRUE]]
  .dpprior_schema_require(
    if (identical(status, "boundary")) !is.null(boundary_reason) else
      is.null(boundary_reason),
    "termination_status", "result.computation.termination.boundary_reason",
    "a reason exactly when status is boundary", boundary_reason
  )
  .dpprior_schema_require(
    termination[["code", exact = TRUE]] %in%
      .DPPRIOR_TERMINATION_CODES[[status]],
    "termination_status", "result.computation.termination.code",
    paste("a closed code for status", status), termination[["code", exact = TRUE]]
  )
  fallback <- computation[["fallback", exact = TRUE]]
  if (identical(status, "infeasible")) {
    .dpprior_schema_require(
      identical(mode, "dual_hard") && raw[["verified", exact = TRUE]] &&
        identical(termination[["code", exact = TRUE]],
                  "certified_infeasible") &&
        identical(termination[["source", exact = TRUE]],
                  "analytic_certificate"),
      "termination_certificate", "result.computation.termination",
      "a verified hard analytic certificate termination", termination
    )
  }
  if (identical(status, "failed")) {
    failure_pair <- paste(
      termination[["code", exact = TRUE]],
      termination[["source", exact = TRUE]], sep = "/"
    )
    has_attempt_error <- any(vapply(
      attempts,
      function(attempt) !is.null(attempt[["error", exact = TRUE]]),
      logical(1)
    ))
    expected_pair <- failure_pair %in% c(
      "no_candidate/no_candidate", "failed/optimizer",
      "failed/fallback_optimizer", "error/error"
    ) &&
      (!identical(failure_pair, "failed/fallback_optimizer") ||
         fallback[["attempted", exact = TRUE]]) &&
      (!identical(failure_pair, "failed/optimizer") ||
         (length(attempts) > 0L && !fallback[["attempted", exact = TRUE]])) &&
      (!identical(failure_pair, "error/error") || has_attempt_error) &&
      (!(length(evaluations) > 0L || length(attempts) == 0L) ||
         identical(failure_pair, "no_candidate/no_candidate"))
    .dpprior_schema_require(
      expected_pair,
      "termination_source", "result.computation.termination",
      paste(
        "an exact code/source pair: no_candidate/no_candidate for a checked",
        "or attempt-free empty result, failed/optimizer or fallback_optimizer",
        "for retained execution failures, or error/error with retained error"
      ), termination
    )
  }
  if (finite_fit && !soft_endpoint) {
    selected_evaluation <- evaluations[[which(selected_evaluation_flags)]]
    if (identical(status, "approximate")) {
      .dpprior_schema_require(
        identical(termination[["code", exact = TRUE]], "approximate") &&
          identical(termination[["source", exact = TRUE]],
                    "candidate_evaluation"),
        "termination_source", "result.computation.termination",
        paste(
          "candidate_evaluation/approximate for every retained diagnostic",
          "or optimizer-unsupported hard/soft/A2 candidate"
        ), termination
      )
    } else if (selected_evaluation[["optimizer_supported", exact = TRUE]]) {
      expected_source <- if (fallback[["used", exact = TRUE]]) {
        "fallback_optimizer"
      } else {
        "optimizer"
      }
      .dpprior_schema_require(
        identical(termination[["source", exact = TRUE]], expected_source) &&
          termination[["code", exact = TRUE]] %in%
            if (identical(status, "boundary")) {
              c("selected", "boundary")
            } else {
              c("selected", "converged")
            },
        "termination_source", "result.computation.termination.source",
        paste(
          "optimizer/fallback_optimizer matching the selected attempt stage",
          "plus selected/converged-or-boundary code"
        ), termination
      )
    } else {
      .dpprior_schema_require(
        identical(status, "approximate") &&
          identical(termination[["code", exact = TRUE]], "approximate") &&
          identical(termination[["source", exact = TRUE]],
                    "candidate_evaluation"),
        "termination_source", "result.computation.termination",
        paste(
          "candidate_evaluation/approximate for a retained selected",
          "candidate without ordinary optimizer support"
        ),
        termination
      )
    }
  }
  if (soft_endpoint) {
    .dpprior_schema_require(
      identical(termination[["code", exact = TRUE]], "endpoint") &&
        identical(termination[["source", exact = TRUE]], "endpoint"),
      "termination_endpoint", "result.computation.termination",
      "the exact endpoint/endpoint code and source for lambda=1", termination
    )
  }
  if (sum(selected_attempt_flags) == 1L) {
    selected_attempt <- attempts[[which(selected_attempt_flags)]]
    .dpprior_schema_require(
      identical(termination[["iterations", exact = TRUE]],
                selected_attempt[["iterations", exact = TRUE]]) &&
        (!status %in% c("converged", "boundary") ||
           .dpprior_attempt_exit_zero(selected_attempt)),
      "termination_selected_attempt", "result.computation.termination",
      paste(
        "iterations identical to the selected attempt and zero exit for",
        "converged/boundary status"
      ),
      termination
    )
  } else {
    .dpprior_schema_require(
      is.null(termination[["iterations", exact = TRUE]]),
      "termination_iterations", "result.computation.termination.iterations",
      "NULL when no execution attempt is selected",
      termination[["iterations", exact = TRUE]]
    )
  }
  invisible(TRUE)
}


.dpprior_validate_result_top_names <- function(raw) {
  .dpprior_schema_require(
    typeof(raw) == "list" && is.list(raw) && !is.object(raw), "type", "result",
    "an ordinary unclassed list", class(raw)
  )
  .dpprior_schema_validate_attributes(raw, "result", "names", "names")
  .dpprior_schema_require(
    !anyDuplicated(names(raw)), "duplicate_names", "result",
    "unique top-level names", names(raw)
  )
  .dpprior_schema_require(
    length(raw) >= length(.DPPRIOR_RESULT_COMMON_FIELDS) &&
      identical(
        names(raw)[seq_along(.DPPRIOR_RESULT_COMMON_FIELDS)],
        .DPPRIOR_RESULT_COMMON_FIELDS
      ),
    "field_names", "result",
    "the exact 18-field canonical spine in construction order", names(raw)
  )
  .dpprior_validate_compatibility(
    raw[["compatibility", exact = TRUE]], "result.compatibility"
  )
  alias_names <- names(
    raw[["compatibility", exact = TRUE]][["top_level_aliases", exact = TRUE]]
  )
  if (is.null(alias_names)) {
    alias_names <- character()
  }
  tail_names <- names(raw)[-(seq_along(.DPPRIOR_RESULT_COMMON_FIELDS))]
  .dpprior_schema_require(
    length(alias_names) <= length(tail_names) &&
      (length(alias_names) == 0L ||
         identical(tail(tail_names, length(alias_names)), alias_names)),
    "alias_order", "result",
    "registered compatibility aliases appended after canonical extensions",
    tail_names
  )
  extension_names <- if (length(alias_names) == 0L) {
    tail_names
  } else {
    head(tail_names, -length(alias_names))
  }
  .dpprior_schema_require(
    all(extension_names %in% .DPPRIOR_RESULT_EXTENSIONS) &&
      identical(
        extension_names,
        .DPPRIOR_RESULT_EXTENSIONS[.DPPRIOR_RESULT_EXTENSIONS %in%
                                     extension_names]
      ),
    "extension_names", "result",
    "known canonical extensions in fixed order", extension_names
  )
  list(extensions = extension_names, aliases = alias_names)
}


.dpprior_validate_mode_extensions <- function(raw, extension_names) {
  mode <- raw[["mode", exact = TRUE]]
  required <- unname(.DPPRIOR_MODE_REQUIRED_EXTENSION[[mode]])
  if (!is.na(required)) {
    .dpprior_schema_require(
      required %in% extension_names, "required_extension", "result",
      paste("required extension", required), extension_names
    )
  }
  allowed <- switch(
    mode,
    a1_proxy = c("proxy", "diagnostics"),
    a2_moment = "diagnostics",
    a2_kl = "diagnostics",
    dual_hard = c("constraint", "diagnostics"),
    dual_soft = c("tradeoff", "diagnostics"),
    dual_legacy = c("legacy", "diagnostics"),
    prior_diagnostics = "diagnostics",
    elicitation_sensitivity = "sensitivity"
  )
  .dpprior_schema_require(
    all(extension_names %in% allowed), "forbidden_extension", "result",
    paste("only extensions", paste(allowed, collapse = ", ")), extension_names
  )
  if ("proxy" %in% extension_names) {
    .dpprior_validate_proxy_extension(raw[["proxy", exact = TRUE]])
  }
  if ("constraint" %in% extension_names) {
    .dpprior_validate_constraint_extension(raw[["constraint", exact = TRUE]])
  }
  if ("tradeoff" %in% extension_names) {
    .dpprior_validate_tradeoff_extension(raw[["tradeoff", exact = TRUE]])
  }
  if ("legacy" %in% extension_names) {
    .dpprior_validate_legacy_extension(raw[["legacy", exact = TRUE]])
  }
  if ("diagnostics" %in% extension_names) {
    .dpprior_validate_diagnostics_extension(
      raw[["diagnostics", exact = TRUE]],
      fit_raw = if (identical(raw[["object_type", exact = TRUE]], "fit")) {
        raw
      } else NULL
    )
  }
  if ("sensitivity" %in% extension_names) {
    .dpprior_validate_sensitivity_extension(raw[["sensitivity", exact = TRUE]])
  }
  invisible(TRUE)
}


.dpprior_validate_result_v1_impl <- function(
    x, native_compatibility_constructor_pending = FALSE) {
  .dpprior_schema_require(
    typeof(x) == "list" && is.list(x),
    "type", "result", "an ordinary list", typeof(x)
  )
  .dpprior_schema_require(
    inherits(x, "dpprior_result"), "class", "result",
    "an object inheriting from dpprior_result", class(x)
  )
  raw <- unclass(x)
  top <- .dpprior_validate_result_top_names(raw)
  .dpprior_validate_schema_record(
    raw[["schema", exact = TRUE]], .DPPRIOR_RESULT_SCHEMA_NAME,
    "result.schema"
  )
  .dpprior_schema_validate_scalar_character(
    raw[["object_type", exact = TRUE]], "result.object_type"
  )
  .dpprior_schema_require(
    raw[["object_type", exact = TRUE]] %in%
      c("fit", "diagnostics", "sensitivity"),
    "object_type", "result.object_type", "fit, diagnostics, or sensitivity",
    raw[["object_type", exact = TRUE]]
  )
  .dpprior_schema_validate_scalar_character(
    raw[["mode", exact = TRUE]], "result.mode"
  )
  .dpprior_schema_require(
    raw[["mode", exact = TRUE]] %in% .DPPRIOR_RESULT_MODES,
    "mode", "result.mode", paste(.DPPRIOR_RESULT_MODES, collapse = ", "),
    raw[["mode", exact = TRUE]]
  )
  .dpprior_schema_require(
    identical(
      unname(.DPPRIOR_MODE_OBJECT_TYPE[[raw[["mode", exact = TRUE]]]]),
      raw[["object_type", exact = TRUE]]
    ),
    "mode_object_type", "result.object_type",
    "the object_type assigned to result.mode",
    raw[["object_type", exact = TRUE]]
  )
  .dpprior_schema_require(
    identical(
      class(x),
      .dpprior_result_classes(
        raw[["mode", exact = TRUE]], raw[["object_type", exact = TRUE]]
      )
    ),
    "mode_class", "result",
    "the exact canonical class vector for mode and object_type", class(x)
  )
  .dpprior_schema_validate_scalar_character(
    raw[["method", exact = TRUE]], "result.method"
  )
  allowed_methods <- .DPPRIOR_MODE_METHODS[[raw[["mode", exact = TRUE]]]]
  .dpprior_schema_require(
    raw[["method", exact = TRUE]] %in% allowed_methods,
    "mode_method", "result.method",
    paste("one of", paste(allowed_methods, collapse = ", "),
          "for the declared mode"),
    raw[["method", exact = TRUE]]
  )
  .dpprior_schema_require(
    .dpprior_schema_is_count(raw[["J", exact = TRUE]], 1L),
    "count", "result.J", "an integer-valued scalar at least 1",
    raw[["J", exact = TRUE]]
  )
  status_record <- list(
    status = raw[["status", exact = TRUE]],
    usable = raw[["usable", exact = TRUE]],
    verified = raw[["verified", exact = TRUE]],
    message = raw[["message", exact = TRUE]]
  )
  .dpprior_validate_status_record(status_record, "result.status_record")
  .dpprior_validate_parameters(
    raw[["parameters", exact = TRUE]], "result.parameters", nullable = TRUE
  )
  if (raw[["status", exact = TRUE]] %in% c("converged", "boundary")) {
    .dpprior_schema_require(
      !is.null(raw[["parameters", exact = TRUE]]) ||
        !identical(raw[["object_type", exact = TRUE]], "fit"),
      "parameters", "result.parameters",
      "finite parameters for a converged/boundary fit", NULL
    )
  }
  if (raw[["status", exact = TRUE]] %in% c("failed", "infeasible") &&
      identical(raw[["object_type", exact = TRUE]], "fit")) {
    .dpprior_schema_require(
      is.null(raw[["parameters", exact = TRUE]]),
      "parameters", "result.parameters",
      "NULL for a failed/infeasible fit", raw[["parameters", exact = TRUE]]
    )
  }
  .dpprior_validate_result_target(
    raw[["target", exact = TRUE]], raw[["mode", exact = TRUE]],
    raw[["object_type", exact = TRUE]], "result.target"
  )
  if (identical(raw[["object_type", exact = TRUE]], "fit")) {
    # The target validators above authenticate the exact canonical class
    # vectors.  Strip those classes from this local validation payload before
    # any nested scientific reads so a registered `$`/`[[` S3 method cannot
    # redirect target evidence after the class boundary has been checked.
    target_payload <- raw[["target", exact = TRUE]]
    target_payload[["K"]] <- unclass(target_payload[["K", exact = TRUE]])
    if ("weight" %in% names(target_payload)) {
      target_payload[["weight"]] <- unclass(
        target_payload[["weight", exact = TRUE]]
      )
    }
    raw[["target"]] <- target_payload
  }
  for (field in c("achieved", "residuals", "tolerances")) {
    .dpprior_schema_validate_named_list(
      raw[[field, exact = TRUE]], paste0("result.", field)
    )
    .dpprior_schema_validate_plain_record_value(
      raw[[field, exact = TRUE]], paste0("result.", field)
    )
  }
  if ("K" %in% names(raw[["achieved", exact = TRUE]])) {
    .dpprior_validate_achieved_K(
      raw[["achieved", exact = TRUE]][["K", exact = TRUE]],
      raw[["J", exact = TRUE]],
      "result.achieved.K"
    )
  }
  if (raw[["mode", exact = TRUE]] %in%
      c("dual_hard", "dual_soft", "dual_legacy") &&
      "weight" %in% names(raw[["achieved", exact = TRUE]])) {
    .dpprior_validate_achieved_weight(
      raw[["achieved", exact = TRUE]][["weight", exact = TRUE]],
      raw[["target", exact = TRUE]][["weight", exact = TRUE]][[
        "metric", exact = TRUE
      ]],
      "result.achieved.weight"
    )
  }
  .dpprior_validate_computation(
    raw[["computation", exact = TRUE]], "result.computation"
  )
  .dpprior_validate_verification(
    raw[["verification", exact = TRUE]], "result.verification"
  )
  .dpprior_validate_provenance(
    raw[["provenance", exact = TRUE]], "result.provenance"
  )
  .dpprior_validate_native_compatibility_boundary(
    raw,
    allow_constructor_pending = native_compatibility_constructor_pending
  )
  if (identical(raw[["object_type", exact = TRUE]], "fit")) {
    .dpprior_schema_require(
      identical(
        raw[["provenance", exact = TRUE]][["projection", exact = TRUE]],
        raw[["target", exact = TRUE]][["K", exact = TRUE]][[
          "provenance", exact = TRUE
        ]][["projection", exact = TRUE]]
      ),
      "result_target_projection", "result.provenance.projection",
      "exact identity with the canonical K-target projection provenance",
      raw[["provenance", exact = TRUE]][["projection", exact = TRUE]]
    )
  }
  .dpprior_validate_mode_extensions(raw, top$extensions)
  .dpprior_validate_A2_parameterless_contract(raw)
  .dpprior_validate_result_verification_contract(raw)

  computation <- raw[["computation", exact = TRUE]]
  provenance <- raw[["provenance", exact = TRUE]]
  .dpprior_schema_require(
    identical(raw[["method", exact = TRUE]],
              computation[["used", exact = TRUE]][["method", exact = TRUE]]),
    "used_method", "result.method",
    "identity with computation.used.method", raw[["method", exact = TRUE]]
  )
  .dpprior_schema_require(
    identical(
      provenance[["requested_method", exact = TRUE]],
      computation[["request", exact = TRUE]][["method", exact = TRUE]]
    ),
    "requested_method", "result.provenance.requested_method",
    "identity with computation.request.method",
    provenance[["requested_method", exact = TRUE]]
  )
  .dpprior_schema_require(
    identical(
      provenance[["selected_method", exact = TRUE]],
      computation[["used", exact = TRUE]][["method", exact = TRUE]]
    ),
    "selected_method", "result.provenance.selected_method",
    "identity with computation.used.method",
    provenance[["selected_method", exact = TRUE]]
  )
  .dpprior_schema_require(
    identical(
      provenance[["is_fallback", exact = TRUE]],
      computation[["fallback", exact = TRUE]][["used", exact = TRUE]]
    ),
    "fallback", "result.provenance.is_fallback",
    "identity with computation.fallback.used",
    provenance[["is_fallback", exact = TRUE]]
  )
  if (identical(raw[["method", exact = TRUE]], "A2-MN+NM")) {
    .dpprior_schema_require(
      computation[["fallback", exact = TRUE]][["used", exact = TRUE]],
      "method_fallback", "result.method",
      "A2-MN+NM only with an explicitly selected fallback", FALSE
    )
  }
  if (!is.null(raw[["parameters", exact = TRUE]])) {
    .dpprior_schema_require(
      identical(
        raw[["parameters", exact = TRUE]][["parameterization", exact = TRUE]],
        provenance[["parameterization", exact = TRUE]]
      ) && identical(
        raw[["parameters", exact = TRUE]][["parameterization", exact = TRUE]],
        computation[["used", exact = TRUE]][["parameterization", exact = TRUE]]
      ),
      "parameterization", "result.parameters.parameterization",
      "identity across parameters, computation, and provenance",
      raw[["parameters", exact = TRUE]][["parameterization", exact = TRUE]]
    )
  }

  verification <- raw[["verification", exact = TRUE]]
  .dpprior_schema_require(
    identical(raw[["verified", exact = TRUE]],
              verification[["passed", exact = TRUE]]),
    "verification", "result.verified",
    "identity with verification.passed", raw[["verified", exact = TRUE]]
  )
  selected_snapshot <- verification[["selected_snapshot", exact = TRUE]]
  if (identical(raw[["object_type", exact = TRUE]], "fit")) {
    for (snapshot_name in c("selected_snapshot", "verifier_snapshot")) {
      snapshot <- verification[[snapshot_name, exact = TRUE]]
      if (!is.null(snapshot)) {
        if ("K" %in% names(snapshot[["achieved", exact = TRUE]])) {
          snapshot_K <- snapshot[["achieved", exact = TRUE]][[
            "K", exact = TRUE
          ]]
          .dpprior_validate_achieved_K(
            snapshot_K,
            raw[["J", exact = TRUE]],
            paste0("result.verification.", snapshot_name, ".achieved.K")
          )
          if (raw[["mode", exact = TRUE]] %in%
              c("a2_moment", "a2_kl", "dual_hard", "dual_soft")) {
            order_name <- if (identical(snapshot_name, "selected_snapshot")) {
              "M_selected"
            } else {
              "M_verification_used"
            }
            authoritative_M <- computation[["orders", exact = TRUE]][[
              order_name, exact = TRUE
            ]]
            .dpprior_schema_require(
              !is.null(authoritative_M) &&
                identical(snapshot[["M", exact = TRUE]], authoritative_M) &&
                identical(snapshot_K[["M", exact = TRUE]], authoritative_M),
              "verification_snapshot_order",
              paste0("result.verification.", snapshot_name),
              paste(
                "snapshot M and achieved.K.M identical to computation.orders.",
                order_name, sep = ""
              ),
              list(
                snapshot_M = snapshot[["M", exact = TRUE]],
                achieved_K_M = snapshot_K[["M", exact = TRUE]],
                authority = authoritative_M
              )
            )
          }
        }
        if (raw[["mode", exact = TRUE]] %in%
            c("dual_hard", "dual_soft", "dual_legacy") &&
            "weight" %in% names(snapshot[["achieved", exact = TRUE]])) {
          .dpprior_validate_achieved_weight(
            snapshot[["achieved", exact = TRUE]][[
              "weight", exact = TRUE
            ]],
            raw[["target", exact = TRUE]][["weight", exact = TRUE]][[
              "metric", exact = TRUE
            ]],
            paste0(
              "result.verification.", snapshot_name, ".achieved.weight"
            )
          )
        }
      }
    }
  }
  if (identical(raw[["object_type", exact = TRUE]], "fit") &&
      !is.null(raw[["parameters", exact = TRUE]])) {
    required_achieved <- if (raw[["mode", exact = TRUE]] %in%
        c("dual_hard", "dual_soft", "dual_legacy")) {
      c("K", "weight")
    } else {
      "K"
    }
    .dpprior_schema_require(
      all(required_achieved %in% names(raw[["achieved", exact = TRUE]])),
      "candidate_achieved", "result.achieved",
      paste("finite candidate fields", paste(required_achieved, collapse = ", ")),
      names(raw[["achieved", exact = TRUE]])
    )
    .dpprior_schema_require(
      !is.null(selected_snapshot), "selected_snapshot",
      "result.verification.selected_snapshot",
      "a selected snapshot for every finite fit candidate", NULL
    )
    for (field in c("parameters", "achieved", "residuals", "tolerances")) {
      .dpprior_schema_require(
        identical(raw[[field, exact = TRUE]],
                  selected_snapshot[[field, exact = TRUE]]),
        "selected_snapshot_identity", paste0("result.", field),
        paste("identity with verification.selected_snapshot.", field,
              sep = ""),
        raw[[field, exact = TRUE]]
      )
    }
    .dpprior_schema_require(
      identical(
        selected_snapshot[["M", exact = TRUE]],
        computation[["orders", exact = TRUE]][["M_selected", exact = TRUE]]
      ),
      "selected_order", "result.verification.selected_snapshot.M",
      "identity with computation.orders.M_selected",
      selected_snapshot[["M", exact = TRUE]]
    )
    attempts <- computation[["attempts", exact = TRUE]]
    selected_id <- computation[["selected_attempt_id", exact = TRUE]]
    if (!is.null(selected_id)) {
      ids <- vapply(attempts, function(attempt) attempt$id, character(1))
      selected_attempt <- attempts[[match(selected_id, ids)]]
      .dpprior_schema_require(
        identical(
          selected_attempt[["candidate_parameters", exact = TRUE]],
          raw[["parameters", exact = TRUE]]
        ),
        "selected_attempt_identity", "result.parameters",
        "identity with the selected attempt candidate parameters",
        raw[["parameters", exact = TRUE]]
      )
    }
  }
  verifier_snapshot <- verification[["verifier_snapshot", exact = TRUE]]
  if (!is.null(verifier_snapshot)) {
    if (!is.null(selected_snapshot)) {
      .dpprior_schema_require(
        identical(
          verifier_snapshot[["parameters", exact = TRUE]],
          selected_snapshot[["parameters", exact = TRUE]]
        ),
        "verifier_candidate", "result.verification.verifier_snapshot.parameters",
        "the fixed selected candidate parameters",
        verifier_snapshot[["parameters", exact = TRUE]]
      )
    }
    .dpprior_schema_require(
      identical(
        verifier_snapshot[["M", exact = TRUE]],
        computation[["orders", exact = TRUE]][["M_verification_used",
                                                exact = TRUE]]
      ),
      "verification_order", "result.verification.verifier_snapshot.M",
      "identity with computation.orders.M_verification_used",
      verifier_snapshot[["M", exact = TRUE]]
    )
  }
  orders <- computation[["orders", exact = TRUE]]
  if (raw[["verified", exact = TRUE]] &&
      !is.null(orders[["M_verification_required", exact = TRUE]])) {
    .dpprior_schema_require(
      !is.null(orders[["M_verification_used", exact = TRUE]]) &&
        !is.null(verifier_snapshot) &&
        !is.null(verifier_snapshot[["M", exact = TRUE]]),
      "verification_order", "result.computation.orders",
      "a retained verifier order and verifier snapshot M for a passed claim",
      orders
    )
  }
  if ("K" %in% names(raw[["achieved", exact = TRUE]])) {
    .dpprior_schema_require(
      identical(
        raw[["achieved", exact = TRUE]][["K", exact = TRUE]][["M",
                                                               exact = TRUE]],
        computation[["orders", exact = TRUE]][["M_selected", exact = TRUE]]
      ),
      "selected_order", "result.achieved.K.M",
      "identity with computation.orders.M_selected",
      raw[["achieved", exact = TRUE]][["K", exact = TRUE]][["M",
                                                              exact = TRUE]]
    )
  }

  mode <- raw[["mode", exact = TRUE]]
  status <- raw[["status", exact = TRUE]]
  usable <- raw[["usable", exact = TRUE]]
  finite_fit <- identical(raw[["object_type", exact = TRUE]], "fit") &&
    !is.null(raw[["parameters", exact = TRUE]])
  attempts <- computation[["attempts", exact = TRUE]]
  selected_flags <- vapply(
    attempts,
    function(attempt) isTRUE(attempt[["selected", exact = TRUE]]),
    logical(1)
  )
  candidate_ledger_mode <- identical(
    raw[["object_type", exact = TRUE]], "fit"
  ) && mode %in% c("a2_moment", "a2_kl", "dual_hard", "dual_soft")
  if (candidate_ledger_mode) {
    .dpprior_validate_result_truth_authority(raw)
    .dpprior_validate_result_candidate_ledger(raw)
  }
  if (!candidate_ledger_mode) {
  for (i in seq_along(attempts)) {
    if (identical(raw[["object_type", exact = TRUE]], "fit") &&
        !selected_flags[[i]]) {
      .dpprior_schema_require(
        !identical(attempts[[i]][["reason_code", exact = TRUE]], "selected"),
        "attempt_reason",
        sprintf("result.computation.attempts[[%d]].reason_code", i),
        "a non-selected reason code when selected=FALSE",
        attempts[[i]][["reason_code", exact = TRUE]]
      )
    }
  }
  soft_endpoint <- finite_fit && identical(mode, "dual_soft") &&
    isTRUE(raw[["tradeoff", exact = TRUE]][["endpoint", exact = TRUE]])
  optimizer_derived <- finite_fit && (
    mode %in% c("a2_moment", "a2_kl", "dual_hard") ||
      (identical(mode, "dual_soft") && !soft_endpoint)
  )
  if (finite_fit) {
    if (optimizer_derived) {
      .dpprior_schema_require(
        length(attempts) > 0L && sum(selected_flags) == 1L &&
          !is.null(computation[["selected_attempt_id", exact = TRUE]]),
        "finite_fit_selection", "result.computation.attempts",
        paste(
          "at least one attempt and exactly one selected attempt for an",
          "optimizer-derived finite candidate"
        ),
        list(count = length(attempts), selected = sum(selected_flags))
      )
    } else if (length(attempts) > 0L) {
      .dpprior_schema_require(
        sum(selected_flags) == 1L &&
          !is.null(computation[["selected_attempt_id", exact = TRUE]]),
        "finite_fit_selection", "result.computation.attempts",
        "exactly one selected attempt when finite-fit attempts are retained",
        sum(selected_flags)
      )
    } else {
      attempt_free_route <- switch(
        mode,
        a1_proxy =
          computation[["termination", exact = TRUE]][[
            "source", exact = TRUE
          ]] %in% c("closed_form", "constructor") &&
          computation[["termination", exact = TRUE]][[
            "code", exact = TRUE
          ]] %in% c("closed_form", "deterministic"),
        dual_soft = soft_endpoint &&
          identical(
            computation[["termination", exact = TRUE]][[
              "source", exact = TRUE
            ]],
            "endpoint"
          ) && identical(
            computation[["termination", exact = TRUE]][[
              "code", exact = TRUE
            ]],
            "endpoint"
          ),
        dual_legacy =
          computation[["termination", exact = TRUE]][[
            "source", exact = TRUE
          ]] %in% c("constructor", "legacy_adapter") &&
          computation[["termination", exact = TRUE]][[
            "code", exact = TRUE
          ]] %in% c("deterministic", "legacy_migration_no_selection"),
        FALSE
      )
      .dpprior_schema_require(
        attempt_free_route,
        "finite_fit_selection", "result.computation.termination.source",
        "a mode-authorized substantive closed-form/endpoint attempt-free route",
        computation[["termination", exact = TRUE]]
      )
    }
  } else if (identical(raw[["object_type", exact = TRUE]], "fit")) {
    snapshot_candidate_claims <- vapply(
      c("selected_snapshot", "verifier_snapshot"),
      function(snapshot_name) {
        snapshot <- verification[[snapshot_name, exact = TRUE]]
        !is.null(snapshot) && (
          !is.null(snapshot[["parameters", exact = TRUE]]) ||
            any(c("K", "weight") %in%
                names(snapshot[["achieved", exact = TRUE]]))
        )
      },
      logical(1)
    )
    .dpprior_schema_require(
      sum(selected_flags) == 0L &&
        is.null(computation[["selected_attempt_id", exact = TRUE]]) &&
        !any(c("K", "weight") %in%
             names(raw[["achieved", exact = TRUE]])) &&
        !any(snapshot_candidate_claims),
      "finite_fit_selection", "result.computation",
      paste(
        "no selected attempt or candidate achieved values when no finite",
        "public fit candidate is retained"
      ),
      list(
        selected_attempt_id = computation[[
          "selected_attempt_id", exact = TRUE
        ]],
        achieved = names(raw[["achieved", exact = TRUE]]),
        snapshot_candidate_claims = snapshot_candidate_claims
      )
    )
    rejection_components <- c(
      residual_adequacy_failed = "residual_adequacy",
      pmf_adequacy_failed = "pmf_adequacy",
      order_stability_failed = "order_stability",
      constraint_verification_failed = "constraint_refined",
      local_optimality_failed = "local_optimality",
      candidate_eligibility_failed = "candidate_eligibility"
    )
    terminal_failure_reasons <- c(
      "optimizer_error", "optimizer_exit_nonzero", "nonfinite_candidate",
      "nonfinite_objective", "candidate_evaluation_failed", "no_candidate"
    )
    for (i in seq_along(attempts)) {
      attempt <- attempts[[i]]
      zero_exit <- identical(attempt[["exit_code", exact = TRUE]], 0L) ||
        identical(attempt[["exit_code", exact = TRUE]], 0)
      finite_success_like <- zero_exit &&
        is.null(attempt[["error", exact = TRUE]]) &&
        !is.null(attempt[["candidate_parameters", exact = TRUE]]) &&
        !is.null(attempt[["candidate_objective", exact = TRUE]])
      reason <- attempt[["reason_code", exact = TRUE]]
      if (finite_success_like && identical(status, "infeasible")) {
        .dpprior_schema_require(
          identical(mode, "dual_hard") &&
            identical(reason, "globally_infeasible_by_certificate") &&
            "infeasibility_certificate" %in%
              names(verification[["components", exact = TRUE]]) &&
            isTRUE(.dpprior_schema_check_pass(
              verification[["components", exact = TRUE]][[
                "infeasibility_certificate", exact = TRUE
              ]]
            )),
          "no_candidate_attempt",
          sprintf("result.computation.attempts[[%d]].reason_code", i),
          paste(
            "global-certificate ineligibility for every numerically finite",
            "attempt retained under certified infeasibility"
          ),
          reason
        )
      } else if (finite_success_like) {
        component_name <- unname(rejection_components[[reason]])
        component <- if (length(component_name) == 1L &&
                         component_name %in%
                           names(verification[["components", exact = TRUE]])) {
          verification[["components", exact = TRUE]][[
            component_name, exact = TRUE
          ]]
        } else {
          NULL
        }
        .dpprior_schema_require(
          length(component_name) == 1L && !is.null(component),
          "no_candidate_attempt",
          sprintf("result.computation.attempts[[%d]].reason_code", i),
          paste(
            "a closed scientific rejection reason linked to a retained",
            "failed verification decision"
          ),
          reason
        )
        .dpprior_validate_decision_check(
          component,
          paste0("result.verification.components.", component_name)
        )
        .dpprior_schema_require(
          !component[["passed", exact = TRUE]] &&
            identical(component[["source", exact = TRUE]],
                      paste0("attempt:", attempt[["id", exact = TRUE]])),
          "no_candidate_attempt",
          sprintf("result.computation.attempts[[%d]]", i),
          paste(
            "a failed substantive decision check whose source identifies",
            "this exact numerical attempt"
          ),
          component
        )
      } else {
        .dpprior_schema_require(
          reason %in% c(terminal_failure_reasons, names(rejection_components),
                        "globally_infeasible_by_certificate"),
          "no_candidate_attempt",
          sprintf("result.computation.attempts[[%d]].reason_code", i),
          "a closed terminal-failure or scientific-rejection reason", reason
        )
      }
    }
  }
  if (finite_fit && status %in% c("converged", "boundary") &&
      length(attempts) > 0L) {
    selected_attempt <- attempts[[which(selected_flags)]]
    .dpprior_schema_require(
      identical(selected_attempt[["exit_code", exact = TRUE]], 0L) ||
        identical(selected_attempt[["exit_code", exact = TRUE]], 0),
      "termination_exit", "result.computation.attempts.selected.exit_code",
      "zero for a converged/boundary selected optimizer attempt",
      selected_attempt[["exit_code", exact = TRUE]]
    )
  }
  boundary_reason <- computation[["termination", exact = TRUE]][[
    "boundary_reason", exact = TRUE
  ]]
  .dpprior_schema_require(
    if (identical(status, "boundary")) !is.null(boundary_reason) else
      is.null(boundary_reason),
    "termination_status", "result.computation.termination.boundary_reason",
    "a reason exactly when status is boundary", boundary_reason
  )
  termination_code <- computation[["termination", exact = TRUE]][[
    "code", exact = TRUE
  ]]
  allowed_termination_codes <- .DPPRIOR_TERMINATION_CODES[[status]]
  if (identical(mode, "prior_diagnostics") &&
      status %in% c("approximate", "failed")) {
    allowed_termination_codes <- unique(c(
      allowed_termination_codes, "diagnostics_recomputed"
    ))
  }
  if (identical(mode, "elicitation_sensitivity") &&
      status %in% c("infeasible", "failed")) {
    allowed_termination_codes <- unique(c(
      allowed_termination_codes, "deterministic"
    ))
  }
  .dpprior_schema_require(
    termination_code %in% allowed_termination_codes,
    "termination_status", "result.computation.termination.code",
    paste(
      "one of", paste(allowed_termination_codes, collapse = ", "),
      "for status", status
    ),
    termination_code
  )
  termination_source <- computation[["termination", exact = TRUE]][[
    "source", exact = TRUE
  ]]
  fallback <- computation[["fallback", exact = TRUE]]
  if (identical(raw[["object_type", exact = TRUE]], "fit") &&
      identical(status, "infeasible")) {
    .dpprior_schema_require(
      identical(mode, "dual_hard") &&
        identical(termination_code, "certified_infeasible") &&
        identical(termination_source, "analytic_certificate") &&
        raw[["verified", exact = TRUE]],
      "termination_certificate", "result.computation.termination",
      paste(
        "dual_hard certified_infeasible from analytic_certificate with a",
        "verified proof"
      ),
      list(mode = mode, termination = computation[["termination", exact = TRUE]])
    )
  }
  if (identical(raw[["object_type", exact = TRUE]], "fit") &&
      identical(status, "failed")) {
    allowed_failure_sources <- switch(
      mode,
      a1_proxy = c("constructor", "mapping", "no_candidate", "error"),
      a2_moment = c("optimizer", "fallback_optimizer", "no_candidate", "error"),
      a2_kl = c("optimizer", "fallback_optimizer", "no_candidate", "error"),
      dual_hard = c("optimizer", "fallback_optimizer", "no_candidate", "error"),
      dual_soft = c("optimizer", "fallback_optimizer", "no_candidate", "error"),
      dual_legacy = c("constructor", "legacy_adapter", "error"),
      character()
    )
    .dpprior_schema_require(
      termination_source %in% allowed_failure_sources,
      "termination_source", "result.computation.termination.source",
      paste(
        "one of", paste(allowed_failure_sources, collapse = ", "),
        "for a failed", mode, "fit"
      ),
      termination_source
    )
    if (identical(termination_source, "fallback_optimizer")) {
      .dpprior_schema_require(
        fallback[["attempted", exact = TRUE]],
        "termination_fallback", "result.computation.termination.source",
        "fallback_optimizer only after a recorded fallback attempt", fallback
      )
    }
    if (identical(termination_source, "optimizer")) {
      .dpprior_schema_require(
        !fallback[["attempted", exact = TRUE]],
        "termination_fallback", "result.computation.termination.source",
        "optimizer when no fallback attempt was recorded", fallback
      )
    }
  }
  if (finite_fit && sum(selected_flags) == 1L) {
    selected_attempt <- attempts[[which(selected_flags)]]
    attempt_iterations <- selected_attempt[["iterations", exact = TRUE]]
    termination_iterations <- computation[["termination", exact = TRUE]][[
      "iterations", exact = TRUE
    ]]
    if (!is.null(attempt_iterations) || !is.null(termination_iterations)) {
      .dpprior_schema_require(
        identical(attempt_iterations, termination_iterations),
        "termination_iterations", "result.computation.termination.iterations",
        "identity with the selected attempt iterations",
        termination_iterations
      )
    }
  }
  selection_modes <- c("a2_kl", "dual_hard", "dual_soft")
  if (!candidate_ledger_mode && finite_fit && identical(mode, "a2_moment") &&
      length(attempts) > 0L) {
    allowed_attempt_methods <- c(
      "A2-MN", "A2-MN+NM", "scaled_log_newton", "nelder_mead_log",
      "fixed_log_parameter_grid"
    )
    selected_attempt <- attempts[[which(selected_flags)]]
    .dpprior_schema_require(
      selected_attempt[["method", exact = TRUE]] %in%
        allowed_attempt_methods &&
        identical(selected_attempt[["reason_code", exact = TRUE]],
                  "selected") &&
        is.null(selected_attempt[["error", exact = TRUE]]),
      "selected_attempt_evidence", "result.computation.attempts.selected",
      paste(
        "an approved A2-MN attempt method with reason=selected and no error"
      ), selected_attempt
    )
    for (i in which(!selected_flags)) {
      .dpprior_schema_require(
        !identical(attempts[[i]][["reason_code", exact = TRUE]], "selected"),
        "attempt_reason",
        sprintf("result.computation.attempts[[%d]].reason_code", i),
        "a non-selected reason code", attempts[[i]][[
          "reason_code", exact = TRUE
        ]]
      )
    }
  }
  if (!candidate_ledger_mode && finite_fit && mode %in% selection_modes &&
      length(attempts) > 0L) {
    allowed_attempt_methods <- switch(
      mode,
      a2_kl = c(
        "A2-KL", "L-BFGS-B", "nlminb", "A2-MN", "A1", "heuristic"
      ),
      dual_hard = c(
        "dual_anchor_hard_inequality", "L-BFGS-B",
        "constrained_profile_optimize", "K_only_L-BFGS-B",
        "verified_candidate", "analytic_monotonicity_feasibility_probe",
        "deterministic_feasible_profile_scan",
        "penalty_L-BFGS-B_diagnostic"
      ),
      dual_soft = c("dual-soft", "L-BFGS-B", "Nelder-Mead", "nlminb")
    )
    for (i in seq_along(attempts)) {
      attempt <- attempts[[i]]
      .dpprior_schema_require(
        attempt[["method", exact = TRUE]] %in% allowed_attempt_methods,
        "attempt_method", sprintf("result.computation.attempts[[%d]].method", i),
        paste("one of", paste(allowed_attempt_methods, collapse = ", ")),
        attempt[["method", exact = TRUE]]
      )
      if (identical(mode, "a2_kl") &&
          attempt[["method", exact = TRUE]] %in%
            c("A2-MN", "A1", "heuristic")) {
        .dpprior_schema_require(
          identical(attempt[["stage", exact = TRUE]], "initialization"),
          "attempt_stage", sprintf(
            "result.computation.attempts[[%d]].stage", i
          ),
          "stage=initialization for retained A2-MN/A1/heuristic starts",
          attempt[["stage", exact = TRUE]]
        )
      }
      if (!is.null(attempt[["candidate_objective", exact = TRUE]])) {
        .dpprior_schema_require(
          attempt[["candidate_objective", exact = TRUE]] >= 0,
          "candidate_objective",
          sprintf("result.computation.attempts[[%d]].candidate_objective", i),
          "a non-negative KL/loss objective", attempt[["candidate_objective",
                                                        exact = TRUE]]
        )
      }
      if (!attempt[["selected", exact = TRUE]]) {
        .dpprior_schema_require(
          !identical(attempt[["reason_code", exact = TRUE]], "selected"),
          "attempt_reason",
          sprintf("result.computation.attempts[[%d]].reason_code", i),
          "a non-selected reason code", attempt[["reason_code", exact = TRUE]]
        )
      }
    }
    selected_attempt <- attempts[[which(selected_flags)]]
    .dpprior_schema_require(
      identical(selected_attempt[["reason_code", exact = TRUE]], "selected") &&
        is.null(selected_attempt[["error", exact = TRUE]]) &&
        (identical(selected_attempt[["exit_code", exact = TRUE]], 0L) ||
           identical(selected_attempt[["exit_code", exact = TRUE]], 0)),
      "selected_attempt_evidence", "result.computation.attempts.selected",
      "reason=selected, no error, and zero exit code", selected_attempt
    )
    eligible <- vapply(attempts, function(attempt) {
      !is.null(attempt[["candidate_parameters", exact = TRUE]]) &&
        !is.null(attempt[["candidate_objective", exact = TRUE]]) &&
        is.null(attempt[["error", exact = TRUE]]) &&
        (identical(attempt[["exit_code", exact = TRUE]], 0L) ||
           identical(attempt[["exit_code", exact = TRUE]], 0)) &&
        (attempt[["selected", exact = TRUE]] ||
           identical(attempt[["reason_code", exact = TRUE]],
                     "eligible_not_selected"))
    }, logical(1))
    .dpprior_schema_require(
      eligible[[which(selected_flags)]], "candidate_selection",
      "result.computation.attempts.selected",
      "the selected attempt among eligible finite successful candidates",
      selected_attempt
    )
    eligible_indices <- which(eligible)
    if (length(eligible_indices) > 1L) {
      for (left_position in seq_len(length(eligible_indices) - 1L)) {
        left_index <- eligible_indices[[left_position]]
        for (right_index in eligible_indices[seq.int(
          left_position + 1L, length(eligible_indices)
        )]) {
          left <- attempts[[left_index]]
          right <- attempts[[right_index]]
          if (identical(left[["method", exact = TRUE]],
                        right[["method", exact = TRUE]]) &&
              identical(left[["candidate_parameters", exact = TRUE]],
                        right[["candidate_parameters", exact = TRUE]])) {
            .dpprior_schema_require(
              identical(left[["candidate_objective", exact = TRUE]],
                        right[["candidate_objective", exact = TRUE]]),
              "candidate_objective_identity",
              sprintf(
                "result.computation.attempts[[%d]].candidate_objective",
                right_index
              ),
              paste(
                "the same deterministic objective for identical method and",
                "candidate parameters"
              ),
              c(
                first = left[["candidate_objective", exact = TRUE]],
                second = right[["candidate_objective", exact = TRUE]]
              )
            )
          }
        }
      }
    }
    objectives <- vapply(
      attempts[eligible],
      function(attempt) attempt[["candidate_objective", exact = TRUE]],
      numeric(1)
    )
    selected_objective <- selected_attempt[["candidate_objective", exact = TRUE]]
    minimum_objective <- min(objectives)
    selection_tolerance <- switch(
      mode,
      a2_kl = computation[["used", exact = TRUE]][[
        "controls", exact = TRUE
      ]][["selection_tolerance", exact = TRUE]],
      dual_hard = raw[["constraint", exact = TRUE]][[
        "optimality", exact = TRUE
      ]][["tie_tolerance", exact = TRUE]],
      dual_soft = raw[["tradeoff", exact = TRUE]][[
        "optimality", exact = TRUE
      ]][["selection_tolerance", exact = TRUE]]
    )
    .dpprior_schema_validate_finite_scalar(
      selection_tolerance, "result.computation.used.controls.selection_tolerance",
      lower = 0
    )
    .dpprior_schema_require(
      selected_objective <= minimum_objective + selection_tolerance,
      "candidate_selection", "result.computation.selected_attempt_id",
      "the minimum eligible objective within the fixed tie tolerance",
      c(selected = selected_objective, minimum = minimum_objective,
        tolerance = selection_tolerance)
    )
    if (raw[["verified", exact = TRUE]]) {
      .dpprior_bind_decision_check(
        verification[["components", exact = TRUE]][[
          "candidate_selection", exact = TRUE
        ]],
        selected_objective - minimum_objective, 0, selection_tolerance, "lte",
        "result.verification.components.candidate_selection"
      )
    }
    if (identical(mode, "dual_hard")) {
      optimality <- raw[["constraint", exact = TRUE]][[
        "optimality", exact = TRUE
      ]]
      target_moments <- .dpprior_result_target_moments(
        raw[["target", exact = TRUE]][["K", exact = TRUE]]
      )
      achieved_K <- raw[["achieved", exact = TRUE]][["K", exact = TRUE]]
      scales <- optimality[["K_scales", exact = TRUE]]
      expected_K_loss <-
        ((achieved_K[["mean", exact = TRUE]] - target_moments[["mean"]]) /
           scales[["mean", exact = TRUE]])^2 +
        ((achieved_K[["variance", exact = TRUE]] -
            target_moments[["variance"]]) /
           scales[["variance", exact = TRUE]])^2
      .dpprior_schema_require(
        identical(optimality[["selection_rule", exact = TRUE]],
                  "minimum_K_loss") &&
          identical(optimality[["selected_K_loss", exact = TRUE]],
                    expected_K_loss) &&
          identical(optimality[["minimum_K_loss", exact = TRUE]],
                    minimum_objective) &&
          identical(optimality[["tie_tolerance", exact = TRUE]],
                    selection_tolerance) &&
          identical(selected_objective, expected_K_loss),
        "hard_optimality_identity", "result.constraint.optimality",
        paste(
          "fresh fixed-scale K loss, eligible minimum, tie tolerance, and",
          "selected attempt objective identities"
        ),
        optimality
      )
      .dpprior_bind_decision_check(
        verification[["components", exact = TRUE]][[
          "perturbation", exact = TRUE
        ]],
        optimality[["perturbation_passed", exact = TRUE]], TRUE, NULL,
        "identical", "result.verification.components.perturbation"
      )
    }
    if (identical(mode, "a2_kl")) {
      .dpprior_schema_require(
        identical(
          selected_objective,
          raw[["residuals", exact = TRUE]][[
            "distribution", exact = TRUE
          ]][["kl", exact = TRUE]]
        ),
        "a2_kl_objective_identity",
        "result.computation.attempts.selected.candidate_objective",
        "identity with freshly recomputed selected-order KL",
        selected_objective
      )
    }
  }
  if (!candidate_ledger_mode && finite_fit && length(attempts) > 0L &&
      sum(selected_flags) == 1L &&
      status %in% c("converged", "boundary", "approximate")) {
    expected_termination_source <- if (fallback[["used", exact = TRUE]]) {
      "fallback_optimizer"
    } else {
      "optimizer"
    }
    .dpprior_schema_require(
      identical(termination_source, expected_termination_source) &&
        computation[["termination", exact = TRUE]][["code", exact = TRUE]] %in%
        c("selected", "converged", "boundary", "approximate"),
      "termination_source", "result.computation.termination",
      paste(
        "optimizer/fallback_optimizer exactly matching fallback.used for a",
        "selected optimizer attempt"
      ),
      computation[["termination", exact = TRUE]]
    )
  }
  }
  if (raw[["verified", exact = TRUE]]) {
    .dpprior_schema_require(
      length(verification[["components", exact = TRUE]]) > 0L &&
        length(verification[["invariants", exact = TRUE]]) > 0L,
      "verification_evidence", "result.verification",
      "non-empty passing components and invariants for a verified result",
      verification
    )
    required_components <- if (
      identical(status, "infeasible") &&
        !identical(mode, "elicitation_sensitivity")
    ) {
      "infeasibility_certificate"
    } else switch(
      mode,
      a2_moment = c(
        "residual_adequacy", "order_stability", "parameter_identity"
      ),
      a2_kl = c(
        "target_identity", "pmf_adequacy", "order_stability",
        "candidate_selection"
      ),
      dual_hard = c(
        "constraint_selected", "constraint_refined", "order_stability",
        "metric_certification", "candidate_selection", "perturbation"
      ),
      dual_soft = if (raw[["tradeoff", exact = TRUE]][[
        "endpoint", exact = TRUE
      ]]) {
        c("endpoint_input_identity", "order_stability")
      } else {
        c(
          "objective_recomputation", "order_stability", "local_optimality",
          "candidate_selection"
        )
      },
      prior_diagnostics = "component_aggregation",
      elicitation_sensitivity = "reconciliation",
      character()
    )
    .dpprior_required_checks_pass(
      verification[["components", exact = TRUE]], required_components,
      "result.verification.components", substantive = TRUE
    )
    if (finite_fit &&
        mode %in% c("a2_moment", "a2_kl", "dual_hard", "dual_soft")) {
      stability <- verification[["stability", exact = TRUE]]
      .dpprior_schema_require(
        !is.null(stability), "stability", "result.verification.stability",
        "a selected-versus-verifier stability record", NULL
      )
      expected_stability <- .dpprior_expected_result_stability(raw)
      recorded_delta <- stability[["delta", exact = TRUE]]
      recorded_tolerance <- stability[["tolerance", exact = TRUE]]
      comparison_scale <- max(
        1, abs(expected_stability$delta),
        abs(expected_stability$tolerance), abs(recorded_delta),
        abs(recorded_tolerance)
      )
      .dpprior_schema_require(
        identical(names(recorded_delta), names(expected_stability$delta)) &&
          identical(
            names(recorded_tolerance), names(expected_stability$tolerance)
          ) &&
          max(abs(recorded_delta - expected_stability$delta)) <=
            1e-12 * comparison_scale &&
          max(abs(recorded_tolerance - expected_stability$tolerance)) <=
            1e-12 * comparison_scale &&
          identical(stability[["formula", exact = TRUE]],
                    expected_stability$formula) &&
          identical(stability[["scale_floor", exact = TRUE]],
                    expected_stability$scale_floor) &&
          stability[["passed", exact = TRUE]],
        "stability_identity", "result.verification.stability",
        paste(
          "deltas recomputed from selected/verifier snapshots and",
          "tolerances recomputed from canonical tolerances"
        ),
        stability
      )
    }
  }
  if (identical(mode, "a2_moment") && finite_fit) {
    .dpprior_schema_require(
      !is.null(selected_snapshot) && !is.null(verifier_snapshot),
      "a2_moment_evidence", "result.verification",
      "selected and independent verifier snapshots", verification
    )
    target_moments <- .dpprior_result_target_moments(
      raw[["target", exact = TRUE]][["K", exact = TRUE]]
    )
    selected_residual <- .dpprior_result_K_residuals(
      selected_snapshot, target_moments,
      "result.verification.selected_snapshot.residuals.K"
    )
    refined_residual <- .dpprior_result_K_residuals(
      verifier_snapshot, target_moments,
      "result.verification.verifier_snapshot.residuals.K"
    )
    selected_tolerance <- .dpprior_result_K_tolerance(
      selected_snapshot[["tolerances", exact = TRUE]][[
        "K_adequacy", exact = TRUE
      ]],
      selected_snapshot[["achieved", exact = TRUE]][["K", exact = TRUE]],
      target_moments,
      "result.verification.selected_snapshot.tolerances.K"
    )
    refined_tolerance <- .dpprior_result_K_tolerance(
      verifier_snapshot[["tolerances", exact = TRUE]][[
        "K_adequacy", exact = TRUE
      ]],
      verifier_snapshot[["achieved", exact = TRUE]][["K", exact = TRUE]],
      target_moments,
      "result.verification.verifier_snapshot.tolerances.K"
    )
    adequacy_value <- c(
      selected.mean = abs(selected_residual[["mean"]]),
      selected.variance = abs(selected_residual[["variance"]]),
      refined.mean = abs(refined_residual[["mean"]]),
      refined.variance = abs(refined_residual[["variance"]]),
      selected.standardized_norm = sqrt(mean(
        (selected_residual / selected_tolerance)^2
      )),
      refined.standardized_norm = sqrt(mean(
        (refined_residual / refined_tolerance)^2
      ))
    )
    adequacy_tolerance <- c(
      selected.mean = selected_tolerance[["mean"]],
      selected.variance = selected_tolerance[["variance"]],
      refined.mean = refined_tolerance[["mean"]],
      refined.variance = refined_tolerance[["variance"]],
      selected.standardized_norm = 1,
      refined.standardized_norm = 1
    )
    components <- verification[["components", exact = TRUE]]
    .dpprior_bind_decision_check(
      components[["residual_adequacy", exact = TRUE]], adequacy_value,
      setNames(rep(0, length(adequacy_value)), names(adequacy_value)),
      adequacy_tolerance, "lte",
      "result.verification.components.residual_adequacy"
    )
    stability <- verification[["stability", exact = TRUE]]
    .dpprior_bind_decision_check(
      components[["order_stability", exact = TRUE]],
      stability[["delta", exact = TRUE]],
      setNames(rep(0, length(stability[["delta", exact = TRUE]])),
               names(stability[["delta", exact = TRUE]])),
      stability[["tolerance", exact = TRUE]], "lte",
      "result.verification.components.order_stability"
    )
    .dpprior_bind_decision_check(
      components[["parameter_identity", exact = TRUE]],
      unlist(selected_snapshot[["parameters", exact = TRUE]][c("a", "b")]),
      unlist(verifier_snapshot[["parameters", exact = TRUE]][c("a", "b")]),
      NULL, "identical",
      "result.verification.components.parameter_identity"
    )
  }
  if (identical(mode, "a2_kl") && finite_fit) {
    target_K <- raw[["target", exact = TRUE]][["K", exact = TRUE]]
    target_pmf <- .dpprior_result_A2_KL_objective_pmf(
      target_K, raw[["J", exact = TRUE]]
    )
    .dpprior_schema_require(
      !is.null(selected_snapshot) && !is.null(verifier_snapshot),
      "a2_kl_evidence", "result.verification",
      "selected and independent verifier PMF snapshots", verification
    )
    target_support <- seq_len(raw[["J", exact = TRUE]])
    target_mean <- sum(target_support * target_pmf)
    target_moments <- c(
      mean = target_mean,
      variance = sum((target_support - target_mean)^2 * target_pmf)
    )
    selected_metrics <- .dpprior_result_distribution_metrics(
      selected_snapshot, target_pmf, target_moments,
      "result.verification.selected_snapshot"
    )
    refined_metrics <- .dpprior_result_distribution_metrics(
      verifier_snapshot, target_pmf, target_moments,
      "result.verification.verifier_snapshot"
    )
    adequacy_value <- c(
      selected.kl = selected_metrics$values[["kl"]],
      selected.l1 = selected_metrics$values[["l1"]],
      selected.mean_scaled = selected_metrics$values[["mean_scaled"]],
      selected.variance_scaled =
        selected_metrics$values[["variance_scaled"]],
      refined.kl = refined_metrics$values[["kl"]],
      refined.l1 = refined_metrics$values[["l1"]],
      refined.mean_scaled = refined_metrics$values[["mean_scaled"]],
      refined.variance_scaled =
        refined_metrics$values[["variance_scaled"]]
    )
    adequacy_tolerance <- c(
      selected.kl = selected_metrics$tolerances[["kl"]],
      selected.l1 = selected_metrics$tolerances[["l1"]],
      selected.mean_scaled = selected_metrics$tolerances[["mean_scaled"]],
      selected.variance_scaled =
        selected_metrics$tolerances[["variance_scaled"]],
      refined.kl = refined_metrics$tolerances[["kl"]],
      refined.l1 = refined_metrics$tolerances[["l1"]],
      refined.mean_scaled = refined_metrics$tolerances[["mean_scaled"]],
      refined.variance_scaled =
        refined_metrics$tolerances[["variance_scaled"]]
    )
    components <- verification[["components", exact = TRUE]]
    .dpprior_bind_decision_check(
      components[["target_identity", exact = TRUE]], target_pmf,
      target_pmf, NULL,
      "identical", "result.verification.components.target_identity"
    )
    .dpprior_bind_decision_check(
      components[["pmf_adequacy", exact = TRUE]], adequacy_value,
      setNames(rep(0, length(adequacy_value)), names(adequacy_value)),
      adequacy_tolerance, "lte",
      "result.verification.components.pmf_adequacy"
    )
    stability <- verification[["stability", exact = TRUE]]
    .dpprior_bind_decision_check(
      components[["order_stability", exact = TRUE]],
      stability[["delta", exact = TRUE]],
      setNames(rep(0, length(stability[["delta", exact = TRUE]])),
               names(stability[["delta", exact = TRUE]])),
      stability[["tolerance", exact = TRUE]], "lte",
      "result.verification.components.order_stability"
    )
  }
  if (identical(status, "approximate") && usable) {
    .dpprior_schema_require(
      mode %in% c("a1_proxy", "dual_legacy", "prior_diagnostics") &&
        provenance[["approximation", exact = TRUE]][["active", exact = TRUE]] &&
        provenance[["approximation", exact = TRUE]][["opt_in", exact = TRUE]],
      "approximation_opt_in", "result.usable",
      "explicit mode-allowed approximation opt-in", provenance$approximation
    )
  }
  if (mode %in% c("dual_hard", "dual_soft") &&
      identical(status, "approximate")) {
    .dpprior_schema_require(
      !usable, "approximation_policy", "result.usable",
      "FALSE for approximate hard/soft results", usable
    )
  }
  if (identical(mode, "a1_proxy")) {
    .dpprior_schema_require(
      identical(status, "approximate") && !raw[["verified", exact = TRUE]],
      "a1_status", "result.status",
      "approximate and unverified for A1 proxy", status
    )
  }
  if (identical(mode, "dual_legacy")) {
    .dpprior_schema_require(
      identical(status, "approximate") &&
        provenance[["legacy", exact = TRUE]][["active", exact = TRUE]],
      "legacy_status", "result.status",
      "approximate with active legacy provenance", status
    )
  } else {
    .dpprior_schema_require(
      !provenance[["legacy", exact = TRUE]][["active", exact = TRUE]],
      "legacy_mode", "result.provenance.legacy.active",
      "FALSE outside dual_legacy", TRUE
    )
  }

  canonical_without_compatibility <- raw
  canonical_without_compatibility[["compatibility"]] <- NULL
  canonical_without_compatibility[top$aliases] <- NULL
  if (mode %in% c("dual_hard", "dual_soft")) {
    input_fit <- provenance[["input_fit", exact = TRUE]]
    target_K <- raw[["target", exact = TRUE]][["K", exact = TRUE]]
    expected_input_target <- list(
      schema = target_K[["schema", exact = TRUE]],
      kind = target_K[["kind", exact = TRUE]],
      J = target_K[["J", exact = TRUE]],
      used = target_K[["used", exact = TRUE]],
      implied = target_K[["implied", exact = TRUE]]
    )
    input_decision_evidence <- if (is.null(input_fit)) NULL else input_fit[[
      "decision_evidence", exact = TRUE
    ]]
    input_decision_target_bound <- if (is.null(input_fit)) {
      FALSE
    } else if (identical(input_fit[["mode", exact = TRUE]], "a2_kl")) {
      !is.null(input_decision_evidence) && identical(
        unclass(input_decision_evidence[["target_K", exact = TRUE]]), target_K
      )
    } else {
      is.null(input_decision_evidence)
    }
    input_authority <- c(
      present = !is.null(input_fit),
      J = !is.null(input_fit) && identical(
        input_fit[["J", exact = TRUE]], raw[["J", exact = TRUE]]
      ),
      compact_target = !is.null(input_fit) && identical(
        input_fit[["target", exact = TRUE]], expected_input_target
      ),
      decision_target = input_decision_target_bound,
      selected_parameters = !is.null(input_fit) && identical(
        input_fit[["selected_snapshot", exact = TRUE]][[
          "parameters", exact = TRUE
        ]], input_fit[["parameters", exact = TRUE]]
      ),
      verifier_parameters = !is.null(input_fit) && identical(
        input_fit[["verifier_snapshot", exact = TRUE]][[
          "parameters", exact = TRUE
        ]], input_fit[["parameters", exact = TRUE]]
      )
    )
    .dpprior_schema_require(
      all(input_authority),
      "input_fit_current_authority", "result.provenance.input_fit",
      paste(
        "a decision-ready input fit with J/target identical to the current",
        "dual result, full A2-KL decision target when applicable, and both",
        "retained snapshots fixed at its parameters"
      ), input_authority
    )
  }
  if (identical(mode, "dual_hard")) {
    hits <- .dpprior_schema_find_forbidden_names(
      canonical_without_compatibility,
      c("tradeoff", "lambda", "weight_loss", "total_loss", "target_residual")
    )
    .dpprior_schema_require(
      length(hits) == 0L, "hard_soft_separation", "result",
      "no soft fields in hard canonical evidence", hits
    )
    constraint <- raw[["constraint", exact = TRUE]]
    weight_target <- raw[["target", exact = TRUE]][["weight", exact = TRUE]]
    .dpprior_schema_require(
      !is.null(provenance[["input_fit", exact = TRUE]]),
      "input_fit", "result.provenance.input_fit",
      "a decision-ready canonical input fit for every hard run", NULL
    )
    .dpprior_schema_require(
      weight_target[["relation", exact = TRUE]] %in% c("at_most", "at_least"),
      "hard_relation", "result.target.weight.relation",
      "at_most or at_least", weight_target[["relation", exact = TRUE]]
    )
    .dpprior_schema_require(
      identical(constraint[["relation", exact = TRUE]],
                weight_target[["relation", exact = TRUE]]) &&
        identical(constraint[["operator", exact = TRUE]],
                  weight_target[["operator", exact = TRUE]]),
      "hard_target_identity", "result.constraint",
      "relation/operator identical to the canonical weight target",
      constraint[c("relation", "operator")]
    )
    if (finite_fit) {
      target_K_moments <- raw[["target", exact = TRUE]][["K", exact = TRUE]][[
        "implied", exact = TRUE
      ]]
      expected_hard_scales <- list(K = list(
        mean = max(abs(target_K_moments[["mean", exact = TRUE]]), 1),
        variance = max(abs(target_K_moments[["variance", exact = TRUE]]), 1)
      ))
      hard_scaling <- computation[["scaling", exact = TRUE]]
      .dpprior_schema_require(
        hard_scaling[["fixed_from_input", exact = TRUE]] &&
          identical(hard_scaling[["requested", exact = TRUE]],
                    expected_hard_scales) &&
          identical(hard_scaling[["used", exact = TRUE]],
                    expected_hard_scales) &&
          identical(hard_scaling[["values", exact = TRUE]],
                    expected_hard_scales) &&
          identical(hard_scaling[["formula", exact = TRUE]],
                    "fixed_from_input_target_max_abs_one") &&
          identical(
            constraint[["optimality", exact = TRUE]][["K_scales", exact = TRUE]],
            expected_hard_scales[["K", exact = TRUE]]
          ),
        "hard_fixed_scales", "result.computation.scaling",
        "K scales fixed once as max(abs(input target moment), 1)",
        list(
          scaling = hard_scaling,
          optimality = constraint[["optimality", exact = TRUE]]
        )
      )
      achieved_weight <- raw[["achieved", exact = TRUE]][[
        "weight", exact = TRUE
      ]][["value", exact = TRUE]]
      .dpprior_schema_validate_finite_scalar(
        achieved_weight, "result.achieved.weight.value", lower = 0, upper = 1
      )
      constraint_tolerances <- raw[["tolerances", exact = TRUE]][[
        "constraint", exact = TRUE
      ]]
      .dpprior_schema_exact_names(
        constraint_tolerances, c("absolute", "relative", "effective"),
        "result.tolerances.constraint"
      )
      expected_effective <-
        constraint_tolerances[["absolute", exact = TRUE]] +
        constraint_tolerances[["relative", exact = TRUE]] * max(
          abs(weight_target[["value", exact = TRUE]]), 1e-8
        )
      constraint_tolerance <- constraint[["tolerance", exact = TRUE]]
      .dpprior_schema_validate_plain_record_value(
        constraint_tolerance, "result.constraint.tolerance",
        numeric_only = TRUE
      )
      .dpprior_schema_require(
        identical(
          constraint_tolerance[["absolute", exact = TRUE]],
          constraint_tolerances[["absolute", exact = TRUE]]
        ) && identical(
          constraint_tolerance[["relative", exact = TRUE]],
          constraint_tolerances[["relative", exact = TRUE]]
        ) && identical(
          constraint_tolerance[["effective", exact = TRUE]],
          expected_effective
        ) && identical(
          constraint_tolerance, constraint_tolerances
        ),
        "constraint_tolerance_identity", "result.constraint.tolerance",
        paste(
          "exact identity with result.tolerances.constraint and",
          "effective = absolute + relative*max(abs(target), 1e-8)"
        ),
        constraint_tolerance
      )
      expected_constraint_residual <- if (identical(
        weight_target[["relation", exact = TRUE]], "at_most"
      )) {
        achieved_weight - weight_target[["value", exact = TRUE]]
      } else {
        weight_target[["value", exact = TRUE]] - achieved_weight
      }
      residual_tolerance <- 1e-12 * max(
        1, abs(expected_constraint_residual),
        abs(constraint[["residual", exact = TRUE]])
      )
      expected_active <- abs(expected_constraint_residual) <= expected_effective
      .dpprior_schema_require(
        abs(constraint[["residual", exact = TRUE]] -
              expected_constraint_residual) <= residual_tolerance &&
          identical(constraint[["active", exact = TRUE]], expected_active),
        "hard_residual_identity", "result.constraint",
        paste(
          "residual recomputed from selected achieved weight and target,",
          "with active determined by effective tolerance"
        ),
        constraint
      )
      selected_weight <- verification[["selected_snapshot", exact = TRUE]][[
        "achieved", exact = TRUE
      ]][["weight", exact = TRUE]][["value", exact = TRUE]]
      refined_weight <- verification[["verifier_snapshot", exact = TRUE]][[
        "achieved", exact = TRUE
      ]][["weight", exact = TRUE]][["value", exact = TRUE]]
      residual_for <- function(value) {
        if (identical(weight_target[["relation", exact = TRUE]], "at_most")) {
          value - weight_target[["value", exact = TRUE]]
        } else {
          weight_target[["value", exact = TRUE]] - value
        }
      }
      selected_residual <- residual_for(selected_weight)
      refined_residual <- residual_for(refined_weight)
      components <- verification[["components", exact = TRUE]]
      .dpprior_bind_decision_check(
        components[["constraint_selected", exact = TRUE]],
        selected_residual, 0, expected_effective, "lte",
        "result.verification.components.constraint_selected"
      )
      .dpprior_bind_decision_check(
        components[["constraint_refined", exact = TRUE]],
        refined_residual, 0, expected_effective, "lte",
        "result.verification.components.constraint_refined"
      )
      stability <- verification[["stability", exact = TRUE]]
      .dpprior_bind_decision_check(
        components[["order_stability", exact = TRUE]],
        stability[["delta", exact = TRUE]][["weight.value"]], 0,
        stability[["tolerance", exact = TRUE]][["weight.value"]], "lte",
        "result.verification.components.order_stability"
      )
      expected_two_order_satisfaction <-
        selected_residual <= expected_effective &&
        refined_residual <= expected_effective
      .dpprior_schema_require(
        identical(constraint[["satisfied", exact = TRUE]],
                  expected_two_order_satisfaction),
        "hard_two_order_satisfaction", "result.constraint.satisfied",
        paste(
          "the exact conjunction of selected and independently refined",
          "constraint residuals against the central constraint tolerance"
        ),
        list(
          recorded = constraint[["satisfied", exact = TRUE]],
          selected = selected_residual, refined = refined_residual,
          tolerance = expected_effective
        )
      )
    }
    if (identical(weight_target[["metric", exact = TRUE]], "wmax_tail_upper")) {
      .dpprior_schema_require(
        identical(weight_target[["relation", exact = TRUE]], "at_most"),
        "unsafe_relation", "result.target.weight.relation",
        "at_most for wmax_tail_upper hard constraints",
        weight_target[["relation", exact = TRUE]]
      )
    }
    if (isTRUE(constraint[["satisfied", exact = TRUE]])) {
      .dpprior_required_checks_pass(
        verification[["components", exact = TRUE]],
        c("constraint_selected", "constraint_refined", "order_stability",
          "metric_certification", "candidate_selection", "perturbation"),
        "result.verification.components"
      )
      selected_candidate_id <- computation[[
        "selected_candidate_id", exact = TRUE
      ]]
      candidate_ids <- vapply(
        computation[["candidate_evaluations", exact = TRUE]],
        function(candidate) candidate[["id", exact = TRUE]], character(1)
      )
      selected_candidate_index <- match(selected_candidate_id, candidate_ids)
      selected_candidate <- if (is.na(selected_candidate_index)) NULL else
        computation[["candidate_evaluations", exact = TRUE]][[
          selected_candidate_index
        ]]
      required_candidate_checks <- c(
        "constraint_selected", "constraint_refined", "order_stability",
        "metric_certification", "perturbation", "invariants"
      )
      candidate_checks_pass <- !is.null(selected_candidate) && identical(
        names(selected_candidate[["checks", exact = TRUE]]),
        required_candidate_checks
      ) && all(vapply(
        selected_candidate[["checks", exact = TRUE]],
        function(check) isTRUE(check[["passed", exact = TRUE]]), logical(1)
      ))
      central_invariants_pass <- all(vapply(
        verification[["invariants", exact = TRUE]],
        function(check) isTRUE(check[["passed", exact = TRUE]]), logical(1)
      ))
      .dpprior_schema_require(
        candidate_checks_pass && central_invariants_pass,
        "hard_satisfied_evidence", "result.constraint.satisfied",
        paste(
          "a satisfied hard claim only when the selected candidate passes",
          "both constraint orders, stability, metric certification,",
          "perturbation, candidate invariants, and central invariants"
        ),
        list(candidate = selected_candidate, invariants = verification$invariants)
      )
    }
    feasibility <- constraint[["feasibility", exact = TRUE]]
    if (status %in% c("converged", "boundary")) {
      .dpprior_schema_require(
        constraint[["satisfied", exact = TRUE]] &&
          identical(
            feasibility[["classification", exact = TRUE]],
            "feasible_candidate"
          ) &&
          feasibility[["verified_candidate_count", exact = TRUE]] >= 1L,
        "hard_status", "result.constraint",
        paste(
          "a satisfied verified feasible candidate for a",
          "converged/boundary hard result"
        ),
        list(status = status, constraint = constraint)
      )
    } else if (status %in% c("failed", "infeasible")) {
      .dpprior_schema_require(
        is.null(constraint[["satisfied", exact = TRUE]]) &&
          !identical(
            feasibility[["classification", exact = TRUE]],
            "feasible_candidate"
          ),
        "hard_status", "result.constraint",
        paste(
          "unavailable constraint decision and no feasible claim for a",
          "failed/infeasible hard status"
        ),
        list(status = status, constraint = constraint)
      )
    } else if (identical(status, "approximate")) {
      .dpprior_schema_require(
        !usable && (
          !isTRUE(constraint[["satisfied", exact = TRUE]]) ||
            identical(
              feasibility[["classification", exact = TRUE]],
              "feasible_candidate"
            )
        ),
        "hard_status", "result.constraint",
        paste(
          "unusable top status while independently verified nested",
          "constraint satisfaction may be retained for an approximate hard",
          "candidate"
        ),
        list(status = status, constraint = constraint)
      )
    } else {
      .dpprior_schema_require(
        !isTRUE(constraint[["satisfied", exact = TRUE]]),
        "hard_status", "result.constraint",
        "no satisfied claim for an unsupported hard status",
        list(status = status, constraint = constraint)
      )
    }
    if (identical(status, "infeasible") && raw[["verified", exact = TRUE]]) {
      certificate <- constraint[["feasibility", exact = TRUE]][[
        "certificate", exact = TRUE
      ]]
      .dpprior_schema_require(
        constraint[["feasibility", exact = TRUE]][["certified_infeasible",
                                                   exact = TRUE]] &&
          !constraint[["feasibility", exact = TRUE]][["feasibility_unknown",
                                                       exact = TRUE]] &&
          !is.null(constraint[["feasibility", exact = TRUE]][["certificate",
                                                              exact = TRUE]]),
        "infeasibility_certificate", "result.constraint.feasibility",
        "a retained certified infeasibility record", constraint$feasibility
      )
      .dpprior_schema_require(
        identical(certificate[["J", exact = TRUE]],
                  raw[["J", exact = TRUE]]) &&
          identical(certificate[["metric", exact = TRUE]],
                  weight_target[["metric", exact = TRUE]]) &&
          identical(certificate[["relation", exact = TRUE]],
                    weight_target[["relation", exact = TRUE]]) &&
          identical(certificate[["target_value", exact = TRUE]],
                    weight_target[["value", exact = TRUE]]) &&
          identical(certificate[["threshold", exact = TRUE]],
                    weight_target[["threshold", exact = TRUE]]) &&
          identical(certificate[["probability", exact = TRUE]],
                    weight_target[["probability", exact = TRUE]]) &&
          identical(
            certificate[["domain", exact = TRUE]][["log_a", exact = TRUE]],
            computation[["used", exact = TRUE]][["controls", exact = TRUE]][[
              "log_bounds", exact = TRUE
            ]]
          ) && identical(
            certificate[["domain", exact = TRUE]][["log_b", exact = TRUE]],
            computation[["used", exact = TRUE]][["controls", exact = TRUE]][[
              "log_bounds", exact = TRUE
            ]]
          ) &&
          identical(certificate[["effective_tolerance", exact = TRUE]],
                    raw[["tolerances", exact = TRUE]][[
                      "constraint", exact = TRUE
                    ]][["effective", exact = TRUE]]) &&
          identical(certificate[["tolerance", exact = TRUE]],
                    raw[["tolerances", exact = TRUE]][[
                      "constraint", exact = TRUE
                    ]][["effective", exact = TRUE]]) &&
          certificate[["minimum", exact = TRUE]][[
            "uncertainty", exact = TRUE
          ]] <= raw[["tolerances", exact = TRUE]][[
            "certificate", exact = TRUE
          ]][["corner_uncertainty", exact = TRUE]] &&
          certificate[["maximum", exact = TRUE]][[
            "uncertainty", exact = TRUE
          ]] <= raw[["tolerances", exact = TRUE]][[
            "certificate", exact = TRUE
          ]][["corner_uncertainty", exact = TRUE]],
        "certificate_target_identity",
        "result.constraint.feasibility.certificate",
        paste(
          "metric/relation/target identity plus the central effective",
          "constraint and corner-uncertainty authorities"
        ),
        certificate
      )
      probe_indices <- which(vapply(attempts, function(attempt) {
        identical(attempt[["stage", exact = TRUE]], "feasibility") &&
          identical(
            attempt[["method", exact = TRUE]],
            "analytic_monotonicity_feasibility_probe"
          )
      }, logical(1)))
      .dpprior_schema_require(
        length(probe_indices) == 1L,
        "certificate_probe_attempt", "result.computation.attempts",
        "exactly one analytic monotonicity feasibility-probe attempt",
        probe_indices
      )
      probe <- attempts[[probe_indices[[1L]]]]
      certificate_corner <- if (identical(
        certificate[["relation", exact = TRUE]], "at_most"
      )) certificate[["minimum", exact = TRUE]] else
        certificate[["maximum", exact = TRUE]]
      expected_probe_bounds <- list(
        log_a = certificate[["domain", exact = TRUE]][["log_a", exact = TRUE]],
        log_b = certificate[["domain", exact = TRUE]][["log_b", exact = TRUE]]
      )
      expected_probe_control <- list(
        M = certificate[["M_selected", exact = TRUE]],
        M_verify = certificate[["M_verification", exact = TRUE]]
      )
      .dpprior_schema_require(
        is.null(probe[["start", exact = TRUE]]) &&
          identical(probe[["bounds", exact = TRUE]], expected_probe_bounds) &&
          identical(probe[["control", exact = TRUE]], expected_probe_control) &&
          identical(probe[["exit_code", exact = TRUE]], 0L) &&
          identical(probe[["iterations", exact = TRUE]], 0L) &&
          is.null(probe[["error", exact = TRUE]]) &&
          is.null(probe[["candidate_parameters", exact = TRUE]]) &&
          identical(
            probe[["candidate_objective", exact = TRUE]],
            certificate_corner[["refined", exact = TRUE]][["value", exact = TRUE]]
          ) &&
          !probe[["selected", exact = TRUE]] &&
          identical(probe[["reason_code", exact = TRUE]],
                    "globally_infeasible_by_certificate") &&
          identical(names(probe[["unavailable", exact = TRUE]]),
                    c("start", "candidate_parameters")),
        "certificate_probe_attempt",
        sprintf("result.computation.attempts[[%d]]", probe_indices[[1L]]),
        paste(
          "the exact certificate domain/orders, zero-exit analytic probe,",
          "refined active-corner objective, and unavailable public parameters"
        ), probe
      )
      certificate_value <- if (identical(
        certificate[["relation", exact = TRUE]], "at_most"
      )) certificate[["lower_bound", exact = TRUE]] else
        certificate[["upper_bound", exact = TRUE]]
      certificate_operator <- if (identical(
        certificate[["relation", exact = TRUE]], "at_most"
      )) "gt" else "lt"
      .dpprior_bind_decision_check(
        verification[["components", exact = TRUE]][[
          "infeasibility_certificate", exact = TRUE
        ]],
        certificate_value, certificate[["target_value", exact = TRUE]],
        certificate[["tolerance", exact = TRUE]], certificate_operator,
        "result.verification.components.infeasibility_certificate"
      )
    }
  }
  if (identical(mode, "dual_soft")) {
    hits <- .dpprior_schema_find_forbidden_names(
      canonical_without_compatibility,
      c("constraint", "constraint_satisfied", "constraint_residual",
        "constraint_slack", "constraint_tolerance", "feasibility",
        "feasibility_unknown")
    )
    .dpprior_schema_require(
      length(hits) == 0L, "hard_soft_separation", "result",
      "no hard fields in soft canonical evidence", hits
    )
    tradeoff <- raw[["tradeoff", exact = TRUE]]
    scaling <- computation[["scaling", exact = TRUE]]
    .dpprior_schema_require(
      !is.null(provenance[["input_fit", exact = TRUE]]),
      "input_fit", "result.provenance.input_fit",
      "a decision-ready canonical input fit for every soft run", NULL
    )
    .dpprior_schema_require(
      scaling[["fixed_from_input", exact = TRUE]] &&
        identical(tradeoff[["scales", exact = TRUE]],
                  scaling[["values", exact = TRUE]]),
      "soft_scaling", "result.tradeoff.scales",
      "fixed scales identical to computation.scaling.values",
      tradeoff[["scales", exact = TRUE]]
    )
    failed_endpoint <- !finite_fit &&
      isTRUE(tradeoff[["endpoint", exact = TRUE]])
    if (failed_endpoint) {
      .dpprior_schema_require(
        identical(status, "failed") &&
          length(computation[["attempts", exact = TRUE]]) == 0L &&
          length(computation[["candidate_evaluations", exact = TRUE]]) == 0L &&
          is.null(computation[["selected_attempt_id", exact = TRUE]]) &&
          is.null(computation[["selected_candidate_id", exact = TRUE]]) &&
          identical(computation[["termination", exact = TRUE]][[
            "code", exact = TRUE
          ]], "no_candidate") &&
          identical(computation[["termination", exact = TRUE]][[
            "source", exact = TRUE
          ]], "no_candidate") &&
          is.null(computation[["termination", exact = TRUE]][[
            "iterations", exact = TRUE
          ]]) &&
          identical(verification[["method", exact = TRUE]], "no_candidate") &&
          !verification[["performed", exact = TRUE]] &&
          !verification[["passed", exact = TRUE]] &&
          identical(verification[["settings", exact = TRUE]], list()) &&
          is.null(verification[["selected_snapshot", exact = TRUE]]) &&
          is.null(verification[["verifier_snapshot", exact = TRUE]]) &&
          is.null(verification[["stability", exact = TRUE]]) &&
          length(verification[["components", exact = TRUE]]) == 0L,
        "soft_failed_endpoint", "result.tradeoff.endpoint",
        paste(
          "lambda=1 without a finite public candidate retains only the exact",
          "attempt-free failed/no-candidate route and no endpoint science claims"
        ),
        list(
          status = status,
          computation = computation[c(
            "attempts", "candidate_evaluations", "selected_attempt_id",
            "selected_candidate_id", "termination"
          )],
          verification = verification
        )
      )
    }
    if (finite_fit) {
      K_target <- raw[["target", exact = TRUE]][["K", exact = TRUE]]
      target_moments <- K_target[["implied", exact = TRUE]]
      .dpprior_schema_require(
        !is.null(target_moments) &&
          all(c("mean", "variance") %in% names(target_moments)),
        "soft_target_moments", "result.target.K.implied",
        "authoritative implied mean and variance", target_moments
      )
      expected_soft_scales <- list(
        K = list(
          mean = max(abs(target_moments[["mean", exact = TRUE]]), 1),
          variance = max(abs(target_moments[["variance", exact = TRUE]]), 1)
        ),
        weight = 1
      )
      .dpprior_schema_require(
        identical(scaling[["requested", exact = TRUE]],
                  expected_soft_scales) &&
          identical(scaling[["used", exact = TRUE]], expected_soft_scales) &&
          identical(scaling[["values", exact = TRUE]], expected_soft_scales) &&
          identical(scaling[["formula", exact = TRUE]],
                    "fixed_scaled_squared_loss") &&
          identical(tradeoff[["scales", exact = TRUE]],
                    expected_soft_scales),
        "soft_fixed_scales", "result.computation.scaling",
        "K scales fixed from the input target and weight scale fixed to one",
        list(scaling = scaling, tradeoff_scales = tradeoff[["scales"]])
      )
      achieved_K <- raw[["achieved", exact = TRUE]][["K", exact = TRUE]]
      K_scales <- tradeoff[["scales", exact = TRUE]][["K", exact = TRUE]]
      expected_K_loss <-
        ((achieved_K[["mean", exact = TRUE]] -
            target_moments[["mean", exact = TRUE]]) /
           K_scales[["mean", exact = TRUE]])^2 +
        ((achieved_K[["variance", exact = TRUE]] -
            target_moments[["variance", exact = TRUE]]) /
           K_scales[["variance", exact = TRUE]])^2
      weight_target <- raw[["target", exact = TRUE]][["weight", exact = TRUE]]
      achieved_weight <- raw[["achieved", exact = TRUE]][[
        "weight", exact = TRUE
      ]][["value", exact = TRUE]]
      expected_raw_residual <- achieved_weight -
        weight_target[["value", exact = TRUE]]
      expected_directed_residual <- switch(
        weight_target[["relation", exact = TRUE]],
        target = expected_raw_residual,
        at_most = max(0, expected_raw_residual),
        at_least = max(0, -expected_raw_residual)
      )
      expected_weight_loss <- (
        expected_directed_residual /
          tradeoff[["scales", exact = TRUE]][["weight", exact = TRUE]]
      )^2
      expected_total_loss <- tradeoff[["lambda", exact = TRUE]] *
        expected_K_loss + (1 - tradeoff[["lambda", exact = TRUE]]) *
        expected_weight_loss
      expected_losses <- c(
        K_loss = expected_K_loss,
        weight_loss = expected_weight_loss,
        total_loss = expected_total_loss,
        target_residual = expected_raw_residual,
        directed_residual = expected_directed_residual
      )
      recorded_losses <- unlist(
        tradeoff[names(expected_losses)], use.names = TRUE
      )
      loss_tolerance <- 1e-12 * max(
        1, abs(expected_losses), abs(recorded_losses)
      )
      .dpprior_schema_require(
        identical(names(recorded_losses), names(expected_losses)) &&
          all(abs(recorded_losses - expected_losses) <= loss_tolerance),
        "soft_loss_identity", "result.tradeoff",
        paste(
          "K/weight/total losses and residuals recomputed from canonical",
          "selected achieved values, targets, lambda, and fixed scales"
        ),
        tradeoff
      )
      optimality <- tradeoff[["optimality", exact = TRUE]]
      selected_attempt <- if (length(attempts) > 0L &&
                              sum(selected_flags) == 1L) {
        attempts[[which(selected_flags)]]
      } else {
        NULL
      }
      if (!tradeoff[["endpoint", exact = TRUE]]) {
        .dpprior_schema_exact_names(
          selected_attempt[["bounds", exact = TRUE]], c("lower", "upper"),
          "result.computation.attempts.selected.bounds"
        )
        lower_bounds <- selected_attempt[["bounds", exact = TRUE]][[
          "lower", exact = TRUE
        ]]
        upper_bounds <- selected_attempt[["bounds", exact = TRUE]][[
          "upper", exact = TRUE
        ]]
        .dpprior_schema_require(
          is.numeric(lower_bounds) && is.numeric(upper_bounds) &&
            !is.object(lower_bounds) && !is.object(upper_bounds) &&
            is.null(dim(lower_bounds)) && is.null(dim(upper_bounds)) &&
            .dpprior_schema_has_only_attributes(lower_bounds) &&
            .dpprior_schema_has_only_attributes(upper_bounds) &&
            length(lower_bounds) == 2L && length(upper_bounds) == 2L &&
            !anyNA(lower_bounds) && !anyNA(upper_bounds) &&
            all(is.finite(lower_bounds)) && all(is.finite(upper_bounds)) &&
            all(lower_bounds < upper_bounds),
          "soft_parameter_bounds",
          "result.computation.attempts.selected.bounds",
          "two finite increasing log-parameter bound pairs",
          selected_attempt[["bounds", exact = TRUE]]
        )
        expected_optimality <- .dpprior_expected_soft_optimality(
          raw, lower_bounds, upper_bounds
        )
        eta <- expected_optimality[["eta", exact = TRUE]]
        boundary_tolerance <- optimality[["boundary_tolerance", exact = TRUE]]
        expected_bound_state <- expected_optimality[[
          "bound_state", exact = TRUE
        ]]
        .dpprior_schema_require(
          all(eta >= lower_bounds - boundary_tolerance) &&
            all(eta <= upper_bounds + boundary_tolerance) &&
            identical(optimality[["bound_state", exact = TRUE]],
                      expected_bound_state),
          "soft_bound_state_identity", "result.tradeoff.optimality.bound_state",
          paste(
            "componentwise lower/interior/upper state recomputed from selected",
            "log parameters, retained bounds, and boundary tolerance"
          ),
          optimality[["bound_state", exact = TRUE]]
        )
        .dpprior_schema_require(
          !is.null(selected_attempt) &&
            (identical(status, "approximate") ||
               optimality[["passed", exact = TRUE]]),
          "soft_optimality", "result.tradeoff.optimality",
          paste(
            "retained optimality evidence, passing for converged/boundary",
            "and allowed to document failure for approximate status"
          ),
          optimality
        )
        selected_objective <- selected_attempt[[
          "candidate_objective", exact = TRUE
        ]]
        soft_numeric_tolerance <- function(recorded, expected) {
          128 * .Machine$double.eps * max(
            1, abs(recorded), abs(expected)
          )
        }
        .dpprior_schema_require(
          abs(optimality[["recorded_objective", exact = TRUE]] -
                selected_objective) <= soft_numeric_tolerance(
                  optimality[["recorded_objective", exact = TRUE]],
                  selected_objective
                ) &&
            abs(optimality[["recomputed_objective", exact = TRUE]] -
                  expected_optimality[["selected_objective", exact = TRUE]]) <=
              soft_numeric_tolerance(
                optimality[["recomputed_objective", exact = TRUE]],
                expected_optimality[["selected_objective", exact = TRUE]]
              ) &&
            abs(optimality[["candidate_objective", exact = TRUE]] -
                  expected_optimality[["selected_objective", exact = TRUE]]) <=
              soft_numeric_tolerance(
                optimality[["candidate_objective", exact = TRUE]],
                expected_optimality[["selected_objective", exact = TRUE]]
              ) &&
            abs(expected_total_loss -
                  expected_optimality[["selected_objective", exact = TRUE]]) <=
              soft_numeric_tolerance(
                expected_total_loss,
                expected_optimality[["selected_objective", exact = TRUE]]
              ),
          "soft_objective_identity", "result.tradeoff.optimality",
          paste(
            "recorded/selected candidate objective at M_selected and a",
            "fresh authoritative soft-backend objective matching public loss"
          ),
          optimality
        )
        objective_policy <- raw[["tolerances", exact = TRUE]][[
          "objective", exact = TRUE
        ]]
        objective_tolerance_for <- function(left, right) {
          objective_policy[["absolute", exact = TRUE]] +
            objective_policy[["relative", exact = TRUE]] * max(
              abs(left), abs(right),
              objective_policy[["scale_floor", exact = TRUE]]
            )
        }
        expected_objective_tolerance <- objective_tolerance_for(
          optimality[["recorded_objective", exact = TRUE]],
          optimality[["recomputed_objective", exact = TRUE]]
        )
        expected_start_tolerance <- objective_tolerance_for(
          optimality[["local_base_objective", exact = TRUE]],
          optimality[["start_objective", exact = TRUE]]
        )
        expected_neighbor_tolerances <- vapply(
          optimality[["neighbor_objectives", exact = TRUE]],
          function(neighbor) objective_tolerance_for(
            optimality[["local_base_objective", exact = TRUE]], neighbor
          ), numeric(1)
        )
        .dpprior_schema_require(
          identical(optimality[["boundary_tolerance", exact = TRUE]],
                    raw[["tolerances", exact = TRUE]][["boundary", exact = TRUE]]) &&
            identical(optimality[["stationarity_tolerance", exact = TRUE]],
                      raw[["tolerances", exact = TRUE]][[
                        "stationarity", exact = TRUE
                      ]][["tolerance", exact = TRUE]]) &&
            identical(optimality[["objective_tolerance", exact = TRUE]],
                      expected_objective_tolerance) &&
            identical(optimality[["start_tolerance", exact = TRUE]],
                      expected_start_tolerance) &&
            identical(optimality[["neighbor_tolerances", exact = TRUE]],
                      expected_neighbor_tolerances) &&
            identical(optimality[["selection_tolerance", exact = TRUE]],
                      raw[["tolerances", exact = TRUE]][[
                        "selection", exact = TRUE
                      ]]) &&
            abs(optimality[["local_base_objective", exact = TRUE]] -
                  expected_optimality[["refined_objective", exact = TRUE]]) <=
              soft_numeric_tolerance(
                optimality[["local_base_objective", exact = TRUE]],
                expected_optimality[["refined_objective", exact = TRUE]]
              ),
          "soft_truth_tolerances", "result.tradeoff.optimality",
          paste(
            "objective/start/neighborhood/boundary/stationarity/selection",
            "tolerances derived from the central fixed policy"
          ), optimality
        )
        .dpprior_schema_require(
          identical(names(optimality[["neighbor_objectives", exact = TRUE]]),
                    names(expected_optimality[[
                      "neighbor_objectives", exact = TRUE
                    ]])) &&
            all(abs(
              optimality[["neighbor_objectives", exact = TRUE]] -
                expected_optimality[["neighbor_objectives", exact = TRUE]]
            ) <= vapply(seq_along(optimality[[
              "neighbor_objectives", exact = TRUE
            ]]), function(index) soft_numeric_tolerance(
              optimality[["neighbor_objectives", exact = TRUE]][[index]],
              expected_optimality[["neighbor_objectives", exact = TRUE]][[index]]
            ), numeric(1))) &&
            all(abs(
              optimality[["gradient", exact = TRUE]] -
                expected_optimality[["gradient", exact = TRUE]]
            ) <= vapply(seq_along(optimality[["gradient", exact = TRUE]]),
                       function(index) soft_numeric_tolerance(
                         optimality[["gradient", exact = TRUE]][[index]],
                         expected_optimality[["gradient", exact = TRUE]][[index]]
                       ), numeric(1))) &&
            abs(optimality[["start_objective", exact = TRUE]] -
                  expected_optimality[["start_objective", exact = TRUE]]) <=
              soft_numeric_tolerance(
                optimality[["start_objective", exact = TRUE]],
                expected_optimality[["start_objective", exact = TRUE]]
              ) &&
            identical(optimality[["gradient_method", exact = TRUE]],
                      expected_optimality[["gradient_method", exact = TRUE]]) &&
            identical(optimality[["stationarity_operator", exact = TRUE]],
                      expected_optimality[[
                        "stationarity_operator", exact = TRUE
                      ]]) &&
            identical(optimality[["source", exact = TRUE]],
                      "independent_refined_objective_verification"),
          "soft_optimality_truth", "result.tradeoff.optimality",
          paste(
            "the exact feasible fixed-neighborhood direction set plus",
            "independently recomputed refined-order neighbors, finite-",
            "difference gradient, KKT methods/operators, and input-fit start"
          ),
          list(recorded = optimality, expected = expected_optimality)
        )
        components <- verification[["components", exact = TRUE]]
        .dpprior_bind_decision_check(
          components[["objective_recomputation", exact = TRUE]],
          optimality[["objective_difference", exact = TRUE]], 0,
          optimality[["objective_tolerance", exact = TRUE]], "abs_lte",
          "result.verification.components.objective_recomputation"
        )
        .dpprior_bind_decision_check(
          components[["local_optimality", exact = TRUE]],
          optimality[["component_pass", exact = TRUE]],
          setNames(
            rep(TRUE, length(optimality[["component_pass", exact = TRUE]])),
            names(optimality[["component_pass", exact = TRUE]])
          ),
          NULL, "identical",
          "result.verification.components.local_optimality"
        )
      }
    }
    if (finite_fit && tradeoff[["endpoint", exact = TRUE]]) {
      input_fit <- provenance[["input_fit", exact = TRUE]]
      target_K <- raw[["target", exact = TRUE]][["K", exact = TRUE]]
      expected_input_target <- list(
        schema = target_K[["schema", exact = TRUE]],
        kind = target_K[["kind", exact = TRUE]],
        J = target_K[["J", exact = TRUE]],
        used = target_K[["used", exact = TRUE]],
        implied = target_K[["implied", exact = TRUE]]
      )
      endpoint_snapshot_reference <- function(snapshot) {
        list(
          parameters = snapshot[["parameters", exact = TRUE]],
          M = snapshot[["M", exact = TRUE]],
          achieved_K = snapshot[["achieved", exact = TRUE]][[
            "K", exact = TRUE
          ]],
          finite = snapshot[["finite", exact = TRUE]],
          source = snapshot[["source", exact = TRUE]]
        )
      }
      endpoint_identity <- c(
        schema = identical(input_fit[["schema", exact = TRUE]],
                           "dpprior.result/1"),
        mode = input_fit[["mode", exact = TRUE]] %in%
          c("a2_moment", "a2_kl", "dual_hard", "dual_soft"),
        method = input_fit[["method", exact = TRUE]] %in%
          .DPPRIOR_MODE_METHODS[[input_fit[["mode", exact = TRUE]]]],
        J = identical(input_fit[["J", exact = TRUE]],
                      raw[["J", exact = TRUE]]),
        status = input_fit[["status", exact = TRUE]] %in%
          c("converged", "boundary"),
        usable = isTRUE(input_fit[["usable", exact = TRUE]]),
        verified = isTRUE(input_fit[["verified", exact = TRUE]]),
        parameters = identical(
          input_fit[["parameters", exact = TRUE]],
          raw[["parameters", exact = TRUE]]
        ),
        target = identical(
          input_fit[["target", exact = TRUE]], expected_input_target
        ),
        decision_evidence = if (identical(
          input_fit[["mode", exact = TRUE]], "a2_kl"
        )) {
          !is.null(input_fit[["decision_evidence", exact = TRUE]]) &&
            identical(
              unclass(input_fit[["decision_evidence", exact = TRUE]][[
                "target_K", exact = TRUE
              ]]), target_K
            )
        } else {
          is.null(input_fit[["decision_evidence", exact = TRUE]])
        },
        selected_snapshot = identical(
          input_fit[["selected_snapshot", exact = TRUE]],
          endpoint_snapshot_reference(
            verification[["selected_snapshot", exact = TRUE]]
          )
        ),
        verifier_snapshot = identical(
          input_fit[["verifier_snapshot", exact = TRUE]],
          endpoint_snapshot_reference(
            verification[["verifier_snapshot", exact = TRUE]]
          )
        )
      )
      .dpprior_bind_decision_check(
        verification[["components", exact = TRUE]][[
          "endpoint_input_identity", exact = TRUE
        ]], endpoint_identity,
        setNames(rep(TRUE, length(endpoint_identity)), names(endpoint_identity)),
        NULL, "identical",
        "result.verification.components.endpoint_input_identity"
      )
      endpoint_stability <- verification[["stability", exact = TRUE]]
      .dpprior_bind_decision_check(
        verification[["components", exact = TRUE]][[
          "order_stability", exact = TRUE
        ]], endpoint_stability[["delta", exact = TRUE]],
        setNames(
          rep(0, length(endpoint_stability[["delta", exact = TRUE]])),
          names(endpoint_stability[["delta", exact = TRUE]])
        ), endpoint_stability[["tolerance", exact = TRUE]], "lte",
        "result.verification.components.order_stability"
      )
      .dpprior_schema_require(
        length(computation[["attempts", exact = TRUE]]) == 0L &&
          identical(computation[["termination", exact = TRUE]][[
            "code", exact = TRUE]], "endpoint") &&
          identical(computation[["termination", exact = TRUE]][[
            "source", exact = TRUE]], "endpoint") &&
          all(endpoint_identity),
        "soft_endpoint", "result.tradeoff.endpoint",
        paste(
          "empty attempts, exact endpoint termination, and public parameters",
          "identical to a decision-ready input fit"
        ), endpoint_identity
      )
    }
  }
  if (identical(mode, "dual_legacy")) {
    if (.dpprior_is_dual_legacy_migration_route(raw)) {
      .dpprior_validate_dual_legacy_migration_boundary(raw)
    } else {
      .dpprior_validate_dual_legacy_invariants(raw)
    }
    hits <- .dpprior_schema_find_forbidden_names(
      canonical_without_compatibility,
      c("constraint", "constraint_satisfied", "constraint_residual",
        "constraint_slack", "constraint_tolerance", "feasibility")
    )
    .dpprior_schema_require(
      length(hits) == 0L, "legacy_hard_separation", "result",
      "no hard feasibility claims in legacy canonical evidence", hits
    )
  }

  .dpprior_validate_alias_identity(raw, top$aliases)
  invisible(x)
}


.dpprior_validate_result_v1 <- function(x, collect = FALSE) {
  .dpprior_schema_collect_validation(.dpprior_validate_result_v1_impl, x, collect)
}


# --- Result constructors -----------------------------------------------------

.dpprior_result_classes <- function(mode, object_type) {
  leading <- switch(
    mode,
    dual_hard = "DPprior_dual_hard",
    dual_soft = "DPprior_dual_soft",
    prior_diagnostics = "DPprior_diagnostics",
    elicitation_sensitivity = "dpprior_elicitation_sensitivity",
    "DPprior_fit"
  )
  unique(c(
    leading,
    if (identical(object_type, "fit")) "DPprior_fit" else character(),
    "dpprior_result", "list"
  ))
}


.dpprior_new_result <- function(object_type,
                                mode,
                                method,
                                J,
                                status,
                                usable,
                                verified,
                                message,
                                parameters,
                                target,
                                achieved,
                                residuals,
                                tolerances,
                                computation,
                                verification,
                                provenance,
                                compatibility = .dpprior_new_compatibility(),
                                extension = list()) {
  .dpprior_schema_validate_scalar_character(object_type, "object_type")
  .dpprior_schema_validate_scalar_character(mode, "mode")
  .dpprior_schema_require(
    mode %in% .DPPRIOR_RESULT_MODES, "mode", "mode",
    paste(.DPPRIOR_RESULT_MODES, collapse = ", "), mode
  )
  .dpprior_schema_require(
    .dpprior_schema_is_count(J, 1L), "count", "J",
    "one ordinary integer-valued scalar at least 1", J
  )
  J <- as.integer(J)
  .dpprior_schema_validate_named_list(extension, "extension")
  .dpprior_schema_require(
    all(names(extension) %in% .DPPRIOR_RESULT_EXTENSIONS),
    "extension_names", "extension", "only canonical extension names",
    names(extension)
  )
  extension <- extension[.DPPRIOR_RESULT_EXTENSIONS[
    .DPPRIOR_RESULT_EXTENSIONS %in% names(extension)
  ]]
  out <- c(
    list(
      schema = .dpprior_schema("result"),
      object_type = object_type,
      mode = mode,
      method = method,
      J = J,
      status = status,
      usable = usable,
      verified = verified,
      message = message,
      parameters = parameters,
      target = target,
      achieved = achieved,
      residuals = residuals,
      tolerances = tolerances,
      computation = computation,
      verification = verification,
      provenance = provenance,
      compatibility = compatibility
    ),
    extension
  )
  class(out) <- .dpprior_result_classes(mode, object_type)
  provenance_list <- typeof(provenance) == "list" && is.list(provenance)
  backend_record <- if (provenance_list) {
    provenance[["backend", exact = TRUE]]
  } else {
    NULL
  }
  native_implementation <- if (typeof(backend_record) == "list" &&
      is.list(backend_record)) {
    backend_record[["implementation", exact = TRUE]]
  } else {
    NULL
  }
  native_specs <- .dpprior_native_compatibility_specs()
  compatibility_list <- typeof(compatibility) == "list" &&
    is.list(compatibility)
  native_pending <- compatibility_list &&
    length(compatibility[["top_level_aliases", exact = TRUE]]) == 0L &&
    any(vapply(
      native_specs,
      function(spec) identical(native_implementation, spec$implementation),
      logical(1)
    ))
  .dpprior_validate_result_v1_impl(
    out, native_compatibility_constructor_pending = native_pending
  )
  out
}


.dpprior_new_fit <- function(mode,
                             method,
                             J,
                             status,
                             usable,
                             verified,
                             message = "",
                             parameters,
                             target,
                             achieved,
                             residuals,
                             tolerances,
                             computation,
                             verification,
                             provenance,
                             compatibility = .dpprior_new_compatibility(),
                             extension = list()) {
  .dpprior_new_result(
    object_type = "fit", mode = mode, method = method, J = J,
    status = status, usable = usable, verified = verified, message = message,
    parameters = parameters, target = target, achieved = achieved,
    residuals = residuals, tolerances = tolerances,
    computation = computation, verification = verification,
    provenance = provenance, compatibility = compatibility,
    extension = extension
  )
}


.dpprior_new_diagnostics <- function(method,
                                     J,
                                     status,
                                     usable,
                                     verified,
                                     message = "",
                                     parameters = NULL,
                                     target = list(),
                                     achieved = list(),
                                     residuals = list(),
                                     tolerances = list(),
                                     computation,
                                     verification,
                                     provenance,
                                     diagnostics,
                                     compatibility = .dpprior_new_compatibility()) {
  .dpprior_new_result(
    object_type = "diagnostics", mode = "prior_diagnostics", method = method,
    J = J, status = status, usable = usable, verified = verified,
    message = message, parameters = parameters, target = target,
    achieved = achieved, residuals = residuals, tolerances = tolerances,
    computation = computation, verification = verification,
    provenance = provenance, compatibility = compatibility,
    extension = list(diagnostics = diagnostics)
  )
}


.dpprior_new_sensitivity <- function(method,
                                     J,
                                     status,
                                     usable,
                                     verified,
                                     message = "",
                                     target = list(),
                                     achieved = list(),
                                     residuals = list(),
                                     tolerances = list(),
                                     computation,
                                     verification,
                                     provenance,
                                     sensitivity,
                                     compatibility = .dpprior_new_compatibility()) {
  .dpprior_new_result(
    object_type = "sensitivity", mode = "elicitation_sensitivity",
    method = method, J = J, status = status, usable = usable,
    verified = verified, message = message, parameters = NULL, target = target,
    achieved = achieved, residuals = residuals, tolerances = tolerances,
    computation = computation, verification = verification,
    provenance = provenance, compatibility = compatibility,
    extension = list(sensitivity = sensitivity)
  )
}


# --- Compatibility aliases and dispatch -------------------------------------

.dpprior_get_exact_path <- function(x, path) {
  .dpprior_schema_validate_scalar_character(path, "alias_path")
  parts <- strsplit(path, ".", fixed = TRUE)[[1L]]
  value <- x
  traversed <- character()
  for (part in parts) {
    traversed <- c(traversed, part)
    .dpprior_schema_require(
      typeof(value) == "list" && is.list(value),
      "alias_path", paste(traversed, collapse = "."),
      "a list record", typeof(value)
    )
    raw <- if (is.object(value)) unclass(value) else value
    .dpprior_schema_require(
      part %in% names(raw), "alias_path", paste(traversed, collapse = "."),
      "an existing exact canonical path", names(raw)
    )
    value <- raw[[part, exact = TRUE]]
  }
  value
}


.dpprior_validate_alias_identity <- function(raw, alias_names,
                                             path = "result") {
  aliases <- raw[["compatibility", exact = TRUE]][[
    "top_level_aliases", exact = TRUE
  ]]
  registry_names <- names(aliases)
  if (is.null(registry_names)) {
    registry_names <- character()
  }
  .dpprior_schema_require(
    identical(registry_names, alias_names), "aliases",
    paste0(path, ".compatibility"),
    "alias registry matching appended aliases", registry_names
  )
  for (alias in alias_names) {
    canonical <- .dpprior_get_exact_path(raw, aliases[[alias]])
    .dpprior_schema_require(
      identical(raw[[alias, exact = TRUE]], canonical),
      "alias_identity", paste0(path, ".", alias),
      paste("identity with", aliases[[alias]]), raw[[alias, exact = TRUE]]
    )
    if (identical(alias, "converged")) {
      expected_converged <- raw[["status", exact = TRUE]] %in%
        c("converged", "boundary") && raw[["usable", exact = TRUE]] &&
        raw[["verified", exact = TRUE]]
      .dpprior_schema_require(
        identical(raw[[alias, exact = TRUE]], expected_converged) &&
          identical(aliases[[alias]], "compatibility.views.converged"),
        "converged_alias", paste0(path, ".converged"),
        "the canonical derived convergence formula via compatibility view",
        raw[[alias, exact = TRUE]]
      )
    }
  }
  invisible(TRUE)
}


.dpprior_append_compatibility_v2 <- function(x,
                                             aliases = character(),
                                             views = list(),
                                             deprecations = list()) {
  .dpprior_schema_validate_character_vector(aliases, "aliases", named = TRUE)
  .dpprior_schema_validate_named_list(views, "views")
  .dpprior_schema_validate_named_list(deprecations, "deprecations")
  raw <- unclass(x)
  .dpprior_schema_require(
    "compatibility" %in% names(raw), "compatibility", "compatibility",
    "a canonical compatibility field", names(raw)
  )
  existing <- raw[["compatibility", exact = TRUE]][[
    "top_level_aliases", exact = TRUE
  ]]
  .dpprior_schema_require(
    length(existing) == 0L, "compatibility", "compatibility.top_level_aliases",
    "an empty alias registry before alias derivation", existing
  )
  .dpprior_schema_require(
    !any(names(aliases) %in% names(raw)), "alias_collision", "aliases",
    "alias names that do not collide with canonical fields",
    intersect(names(aliases), names(raw))
  )
  compatibility <- .dpprior_new_compatibility(aliases, views, deprecations)
  raw[["compatibility"]] <- compatibility
  for (alias in names(aliases)) {
    raw[[alias]] <- .dpprior_get_exact_path(raw, aliases[[alias]])
  }
  class(raw) <- class(x)
  .dpprior_validate_object(raw)
  raw
}


.dpprior_validate_object_impl <- function(x) {
  .dpprior_schema_require(
    typeof(x) == "list" && is.list(x),
    "type", "object", "an ordinary list", typeof(x)
  )
  raw <- unclass(x)
  .dpprior_schema_require(
    "schema" %in% names(raw), "schema", "object.schema",
    "an exact schema field", names(raw)
  )
  schema <- raw[["schema", exact = TRUE]]
  .dpprior_schema_require(
    typeof(schema) == "list" && is.list(schema) && !is.object(schema) &&
      "name" %in% names(schema),
    "schema", "object.schema", "an ordinary schema record", class(schema)
  )
  name <- schema[["name", exact = TRUE]]
  if (identical(name, .DPPRIOR_RESULT_SCHEMA_NAME)) {
    return(.dpprior_validate_result_v1_impl(x))
  }
  if (identical(name, .DPPRIOR_TARGET_SCHEMA_NAME)) {
    return(.dpprior_validate_target_v1_impl(x))
  }
  if (identical(name, .DPPRIOR_WEIGHT_TARGET_SCHEMA_NAME)) {
    return(.dpprior_validate_weight_target_v1(x))
  }
  .dpprior_schema_fail(
    "schema_name", "object.schema.name",
    "a supported canonical schema name", name
  )
}


.dpprior_validate_object <- function(x, collect = FALSE) {
  .dpprior_schema_collect_validation(.dpprior_validate_object_impl, x, collect)
}
