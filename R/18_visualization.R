# =============================================================================
# Module 18: Visualization (Final Version)
# =============================================================================
#
# Publication-ready visualization functions for DPprior.
#
# Key features:
# - Direct parameter API: plot_alpha_prior(a = 1.6, b = 1.2)
# - Canonical status-aware K_J PMF data with no invented distribution fallback
# - No annotation boxes (stats in subtitle only)
# - gtable dashboard (no patchwork dependency)
# - Estimand-labelled W_SB threshold diagnostics
#
# Author: JoonHo Lee (jlee296@ua.edu)
# Date: December 2025
# Part of: DPprior R Package
# =============================================================================

# =============================================================================
# Internal Utilities
# =============================================================================

.dpprior_has_ggplot2 <- function() {
  requireNamespace("ggplot2", quietly = TRUE)
}

.dpprior_require_ggplot2 <- function() {
  if (!.dpprior_has_ggplot2()) {
    stop("ggplot2 is required for engine = 'ggplot2'. Install ggplot2 or use engine = 'base'.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.dpprior_coalesce <- function(x, y) {
  if (is.null(x)) y else x
}

.dpprior_get_quad_nodes_default <- function() {
  if (exists(".QUAD_NODES_DEFAULT", inherits = TRUE)) {
    get(".QUAD_NODES_DEFAULT", inherits = TRUE)
  } else {
    80L
  }
}

.dpprior_visualization_abort <- function(message, classes, code, ...) {
  stop(.dpprior_new_condition(
    message = message,
    classes = unique(c(classes, "dpprior_error", "error")),
    code = code,
    ...
  ))
}

.dpprior_visualization_plain_record <- function(x) {
  if (is.null(x)) return(NULL)
  if (typeof(x) == "list") {
    record <- if (is.object(x)) unclass(x) else x
    output <- lapply(seq_along(record), function(index) {
      .dpprior_visualization_plain_record(record[[index]])
    })
    names(output) <- attr(record, "names", exact = TRUE)
    return(output)
  }
  if (is.atomic(x) && is.object(x)) return(unclass(x))
  x
}

.dpprior_visualization_candidate_unavailable <- function(view, operation) {
  migrated <- isTRUE(view[["core", exact = TRUE]][["migrated", exact = TRUE]])
  code <- if (migrated) {
    "migrated_fit_candidate_unavailable"
  } else {
    "canonical_fit_candidate_unavailable"
  }
  .dpprior_visualization_abort(
    paste(
      "The canonical fit has no public parameter and achieved-K candidate for",
      operation,
      "; no compatibility value or recomputed substitute was used."
    ),
    classes = c(
      "dpprior_s3_unavailable_error", "dpprior_visualization_data_error"
    ),
    code = code,
    operation = operation,
    schema = "dpprior.result/1",
    mode = view[["core", exact = TRUE]][["mode", exact = TRUE]],
    status = view[["core", exact = TRUE]][["status", exact = TRUE]],
    usable = view[["core", exact = TRUE]][["usable", exact = TRUE]],
    verified = view[["core", exact = TRUE]][["verified", exact = TRUE]],
    feasibility = if (identical(
      view[["core", exact = TRUE]][["mode", exact = TRUE]], "dual_hard"
    )) {
      constraint <- view[["constraint", exact = TRUE]]
      if (is.null(constraint)) NULL else {
        constraint[["feasibility", exact = TRUE]]
      }
    } else {
      NULL
    },
    action = "inspect_status_and_refit_if_appropriate"
  )
}

.dpprior_visualization_fit_view <- function(
    fit, operation = "visualization", require_candidate = TRUE) {
  validated <- .dpprior_require_schema(
    fit, kind = "fit", allow_legacy = FALSE
  )
  raw <- .dpprior_visualization_plain_record(unclass(validated))
  target <- raw[["target", exact = TRUE]]
  achieved <- raw[["achieved", exact = TRUE]]
  parameters <- raw[["parameters", exact = TRUE]]
  target_K <- target[["K", exact = TRUE]]
  achieved_K <- achieved[["K", exact = TRUE]]
  provenance <- raw[["provenance", exact = TRUE]]
  migration <- provenance[["migration", exact = TRUE]]
  mode <- raw[["mode", exact = TRUE]]
  view <- list(
    raw = raw,
    core = list(
      schema = "dpprior.result/1",
      mode = mode,
      method = raw[["method", exact = TRUE]],
      J = raw[["J", exact = TRUE]],
      status = raw[["status", exact = TRUE]],
      usable = raw[["usable", exact = TRUE]],
      verified = raw[["verified", exact = TRUE]],
      message = raw[["message", exact = TRUE]],
      migrated = !is.null(migration) && identical(
        migration[["adapter", exact = TRUE]], "upgrade_DPprior_object"
      )
    ),
    parameters = parameters,
    target = list(
      K = target_K,
      weight = target[["weight", exact = TRUE]]
    ),
    achieved = list(
      K = achieved_K,
      weight = achieved[["weight", exact = TRUE]]
    ),
    M = if (is.null(achieved_K)) NULL else {
      achieved_K[["M", exact = TRUE]]
    },
    constraint = if (identical(mode, "dual_hard")) {
      raw[["constraint", exact = TRUE]]
    } else {
      NULL
    },
    tradeoff = if (identical(mode, "dual_soft")) {
      raw[["tradeoff", exact = TRUE]]
    } else {
      NULL
    },
    legacy = if (identical(mode, "dual_legacy")) {
      raw[["legacy", exact = TRUE]]
    } else {
      NULL
    },
    input_fit = provenance[["input_fit", exact = TRUE]]
  )
  candidate_available <- !is.null(parameters) && !is.null(achieved_K)
  view[["candidate_available"]] <- candidate_available
  if (isTRUE(require_candidate) && !candidate_available) {
    .dpprior_visualization_candidate_unavailable(view, operation)
  }
  view
}

.dpprior_visualization_input_fit_view <- function(input_fit) {
  if (is.null(input_fit)) return(NULL)
  list(
    core = list(
      schema = input_fit[["schema", exact = TRUE]],
      mode = input_fit[["mode", exact = TRUE]],
      method = input_fit[["method", exact = TRUE]],
      J = input_fit[["J", exact = TRUE]],
      status = input_fit[["status", exact = TRUE]],
      usable = input_fit[["usable", exact = TRUE]],
      verified = input_fit[["verified", exact = TRUE]]
    ),
    parameters = input_fit[["parameters", exact = TRUE]],
    target = input_fit[["target", exact = TRUE]],
    achieved = list(
      K = input_fit[["selected_snapshot", exact = TRUE]][[
        "achieved_K", exact = TRUE
      ]]
    ),
    M = input_fit[["selected_snapshot", exact = TRUE]][["M", exact = TRUE]],
    verifier = input_fit[["verifier_snapshot", exact = TRUE]]
  )
}

.dpprior_visualization_compact_target_K <- function(target_K) {
  list(
    schema = target_K[["schema", exact = TRUE]],
    kind = target_K[["kind", exact = TRUE]],
    J = target_K[["J", exact = TRUE]],
    used = target_K[["used", exact = TRUE]],
    implied = target_K[["implied", exact = TRUE]]
  )
}

.dpprior_visualization_dual_view <- function(
    fit_dual, fit_K_only = NULL, operation = "dual visualization") {
  current <- .dpprior_visualization_fit_view(
    fit_dual, operation = operation, require_candidate = TRUE
  )
  mode <- current[["core", exact = TRUE]][["mode", exact = TRUE]]
  if (!mode %in% c("dual_hard", "dual_soft", "dual_legacy")) {
    .dpprior_visualization_abort(
      paste(operation, "requires a canonical hard, soft, or legacy dual fit."),
      classes = c("dpprior_visualization_input_error", "dpprior_invalid_input"),
      code = "visualization_not_dual_fit",
      operation = operation,
      mode = mode
    )
  }
  input_fit <- current[["input_fit", exact = TRUE]]
  if (is.null(input_fit)) {
    migrated <- isTRUE(current[["core", exact = TRUE]][["migrated", exact = TRUE]])
    .dpprior_visualization_abort(
      paste(
        "The canonical dual fit has no authoritative provenance.input_fit",
        "lineage for", operation,
        "; no compatibility initialization or explicit substitute was used."
      ),
      classes = c(
        "dpprior_s3_unavailable_error", "dpprior_visualization_data_error"
      ),
      code = if (migrated) {
        "migrated_comparison_lineage_unavailable"
      } else {
        "canonical_comparison_lineage_unavailable"
      },
      operation = operation,
      mode = mode,
      status = current[["core", exact = TRUE]][["status", exact = TRUE]],
      action = "refit_with_current_API"
    )
  }
  baseline <- .dpprior_visualization_input_fit_view(input_fit)
  if (!is.null(fit_K_only)) {
    supplied <- .dpprior_visualization_fit_view(
      fit_K_only, operation = paste(operation, "explicit baseline"),
      require_candidate = TRUE
    )
    reference_error <- NULL
    supplied_reference <- tryCatch({
      if (!exists(
        ".dpprior_v2_canonical_input_fit_reference",
        mode = "function", inherits = TRUE
      )) {
        stop("canonical input-fit reference builder is unavailable")
      }
      supplied_validated <- .dpprior_require_schema(
        fit_K_only, kind = "fit", allow_legacy = FALSE
      )
      supplied_raw <- unclass(supplied_validated)
      supplied_target <- unclass(supplied_raw[["target", exact = TRUE]])
      .dpprior_visualization_plain_record(
        .dpprior_v2_canonical_input_fit_reference(
          supplied_raw, supplied_target[["K", exact = TRUE]]
        )
      )
    }, error = function(condition) {
      reference_error <<- condition
      NULL
    })
    if (is.null(supplied_reference) ||
        !identical(supplied_reference, input_fit)) {
      .dpprior_visualization_abort(
        paste(
          "fit_K_only does not exactly match the authoritative",
          "provenance.input_fit lineage."
        ),
        classes = c(
          "dpprior_visualization_input_error", "dpprior_invalid_input"
        ),
        code = "comparison_input_fit_mismatch",
        operation = operation,
        mismatch = "canonical_input_fit_reference",
        parent = reference_error
      )
    }
  }
  list(
    mode = mode,
    current = current,
    baseline = baseline,
    target_weight = current[["target", exact = TRUE]][["weight", exact = TRUE]],
    achieved_weight = current[["achieved", exact = TRUE]][["weight", exact = TRUE]],
    constraint = current[["constraint", exact = TRUE]],
    tradeoff = current[["tradeoff", exact = TRUE]],
    legacy = current[["legacy", exact = TRUE]]
  )
}

.dpprior_visualization_fit_caption <- function(view, estimand) {
  if (is.null(view)) return(NULL)
  core <- view[["core", exact = TRUE]]
  caption <- sprintf(
    "fit mode=%s; status=%s; verified=%s; estimand=%s",
    core[["mode", exact = TRUE]],
    core[["status", exact = TRUE]],
    if (isTRUE(core[["verified", exact = TRUE]])) "yes" else "no",
    estimand
  )
  if (identical(core[["mode", exact = TRUE]], "dual_legacy")) {
    caption <- paste0(
      caption,
      "; deprecated legacy dual fit (descriptive single-fit view only)"
    )
  }
  caption
}

.dpprior_visualization_weight_target_text <- function(target_weight) {
  if (is.null(target_weight)) return("weight target unavailable")
  estimand <- target_weight[["estimand", exact = TRUE]]
  operator <- target_weight[["operator", exact = TRUE]]
  value <- target_weight[["value", exact = TRUE]]
  metric <- target_weight[["metric", exact = TRUE]]
  sprintf("%s [%s] %s %.4g", estimand, metric, operator, value)
}

.dpprior_visualization_dual_title <- function(dual) {
  if (identical(dual[["mode", exact = TRUE]], "dual_hard")) {
    return("Hard-Constraint Dual Comparison")
  }
  tradeoff <- dual[["tradeoff", exact = TRUE]]
  sprintf(
    "Soft-Trade-off Dual Comparison (lambda=%.3g)",
    tradeoff[["lambda", exact = TRUE]]
  )
}

.dpprior_visualization_curve_view <- function(tradeoff_data, metric) {
  expected_columns <- c(
    "point_id", "lambda", "mode", "metric", "relation", "target_value",
    "status", "usable", "verified", "converged", "outcome",
    "condition_class", "condition_code", "condition_message", "a", "b",
    "mu_K", "var_K", "achieved_weight", "target_residual", "K_loss",
    "weight_loss", "total_loss", "attempt_count", "selected_method",
    "warm_start_from", "w_loss", "w1_prob_gt_50", "E_w1"
  )
  expected_attributes <- c(
    "names", "row.names", "class", "fits", "conditions", "metadata"
  )
  if (!identical(class(tradeoff_data), c(
    "dpprior_tradeoff_curve", "data.frame"
  ))) {
    .dpprior_visualization_abort(
      paste(
        "tradeoff_data must be a dpprior_tradeoff_curve returned by",
        "compute_tradeoff_curve()."
      ),
      classes = c(
        "dpprior_visualization_input_error", "dpprior_invalid_input"
      ),
      code = "tradeoff_curve_contract",
      operation = "trade-off curve visualization"
    )
  }

  observed_columns <- attr(tradeoff_data, "names", exact = TRUE)
  ordinary_column_names <- identical(typeof(observed_columns), "character") &&
    !is.object(observed_columns) && is.null(attributes(observed_columns))
  if (!ordinary_column_names ||
      !identical(observed_columns, expected_columns)) {
    .dpprior_visualization_abort(
      "tradeoff_data does not have the exact canonical curve columns.",
      classes = c(
        "dpprior_visualization_input_error", "dpprior_invalid_input"
      ),
      code = "tradeoff_curve_columns",
      operation = "trade-off curve visualization"
    )
  }
  observed_attributes <- names(attributes(tradeoff_data))
  if (!setequal(observed_attributes, expected_attributes)) {
    .dpprior_visualization_abort(
      "tradeoff_data has missing or unexpected retained-evidence attributes.",
      classes = c(
        "dpprior_visualization_input_error", "dpprior_invalid_input"
      ),
      code = "tradeoff_curve_retained_evidence",
      operation = "trade-off curve visualization",
      missing = setdiff(expected_attributes, observed_attributes),
      extra = setdiff(observed_attributes, expected_attributes)
    )
  }

  raw <- unclass(tradeoff_data)
  point_id_column <- raw[["point_id", exact = TRUE]]
  point_id_is_ordinary <- identical(typeof(point_id_column), "character") &&
    !is.object(point_id_column) && is.null(attributes(point_id_column))
  if (!point_id_is_ordinary) {
    .dpprior_visualization_abort(
      "tradeoff_data point_id must be one ordinary character vector.",
      classes = c(
        "dpprior_visualization_input_error", "dpprior_invalid_input"
      ),
      code = "tradeoff_curve_values",
      operation = "trade-off curve visualization",
      mismatched_fields = "point_id"
    )
  }
  row_count <- length(point_id_column)
  expected_types <- c(
    point_id = "character", lambda = "double", mode = "character",
    metric = "character", relation = "character", target_value = "double",
    status = "character", usable = "logical", verified = "logical",
    converged = "logical", outcome = "character",
    condition_class = "character", condition_code = "character",
    condition_message = "character", a = "double", b = "double",
    mu_K = "double", var_K = "double", achieved_weight = "double",
    target_residual = "double", K_loss = "double", weight_loss = "double",
    total_loss = "double", attempt_count = "integer",
    selected_method = "character", warm_start_from = "character",
    w_loss = "double", w1_prob_gt_50 = "double", E_w1 = "double"
  )
  ordinary_columns <- vapply(expected_columns, function(field) {
    column <- raw[[field, exact = TRUE]]
    identical(typeof(column), expected_types[[field]]) &&
      !is.object(column) && is.null(attributes(column)) &&
      identical(length(column), row_count)
  }, logical(1L))
  if (row_count < 1L || !all(ordinary_columns) ||
      !identical(attr(tradeoff_data, "row.names", exact = TRUE),
                 seq_len(row_count))) {
    .dpprior_visualization_abort(
      paste(
        "tradeoff_data columns must be ordinary unclassed atomic vectors of",
        "the exact producer types and row count."
      ),
      classes = c(
        "dpprior_visualization_input_error", "dpprior_invalid_input"
      ),
      code = "tradeoff_curve_values",
      operation = "trade-off curve visualization",
      mismatched_fields = names(ordinary_columns)[!ordinary_columns]
    )
  }

  lambda <- raw[["lambda", exact = TRUE]]
  value <- raw[[metric, exact = TRUE]]
  point_id <- raw[["point_id", exact = TRUE]]
  required_finite <- c("lambda", "target_value")
  nullable_numeric <- c(
    "a", "b", "mu_K", "var_K", "achieved_weight", "target_residual",
    "K_loss", "weight_loss", "total_loss", "w_loss",
    "w1_prob_gt_50", "E_w1"
  )
  nonnullable_character <- c(
    "point_id", "mode", "metric", "relation", "status", "outcome",
    "selected_method", "warm_start_from"
  )
  nullable_character <- c(
    "condition_class", "condition_code", "condition_message"
  )
  numeric_values_ok <- all(vapply(required_finite, function(field) {
    all(is.finite(raw[[field, exact = TRUE]]))
  }, logical(1L))) && all(vapply(nullable_numeric, function(field) {
    column <- raw[[field, exact = TRUE]]
    !any(is.nan(column)) && all(is.na(column) | is.finite(column))
  }, logical(1L)))
  character_values_ok <- all(vapply(nonnullable_character, function(field) {
    column <- raw[[field, exact = TRUE]]
    !anyNA(column) && all(nzchar(column))
  }, logical(1L))) && all(vapply(nullable_character, function(field) {
    column <- raw[[field, exact = TRUE]]
    all(is.na(column) | nzchar(column))
  }, logical(1L)))
  logical_values_ok <- all(vapply(
    c("usable", "verified", "converged"),
    function(field) !anyNA(raw[[field, exact = TRUE]]),
    logical(1L)
  ))
  attempt_values <- raw[["attempt_count", exact = TRUE]]
  attempt_values_ok <- !anyNA(attempt_values) && all(attempt_values >= 0L)
  ordering_ok <- all(lambda > 0 & lambda <= 1) &&
    !anyDuplicated(lambda) && !is.unsorted(lambda, strictly = TRUE) &&
    identical(order(lambda, point_id, method = "radix"), seq_len(row_count))
  if (!numeric_values_ok || !character_values_ok || !logical_values_ok ||
      !attempt_values_ok || !ordering_ok || anyDuplicated(point_id)) {
    .dpprior_visualization_abort(
      "tradeoff_data contains invalid canonical curve values or ordering.",
      classes = c(
        "dpprior_visualization_input_error", "dpprior_invalid_input"
      ),
      code = "tradeoff_curve_values",
      operation = "trade-off curve visualization",
      metric = metric
    )
  }

  ordinary_named_list <- function(x, expected_names = NULL,
                                  expected_length = NULL) {
    if (typeof(x) != "list" || is.object(x)) return(FALSE)
    attribute_names <- names(attributes(x))
    if (!all(attribute_names %in% "names")) return(FALSE)
    container_names <- attr(x, "names", exact = TRUE)
    names_are_ordinary <- is.null(container_names) || (
      identical(typeof(container_names), "character") &&
        !is.object(container_names) && is.null(attributes(container_names)) &&
        identical(length(container_names), length(x))
    )
    if (!names_are_ordinary) return(FALSE)
    if (!is.null(expected_length) &&
        !identical(length(x), expected_length)) return(FALSE)
    if (!is.null(expected_names) &&
        !identical(container_names, expected_names)) {
      return(FALSE)
    }
    TRUE
  }
  ordinary_metadata_tree <- function(x) {
    if (is.null(x)) return(TRUE)
    if (identical(typeof(x), "list")) {
      if (!ordinary_named_list(x)) return(FALSE)
      return(all(vapply(seq_along(x), function(index) {
        ordinary_metadata_tree(x[[index]])
      }, logical(1L))))
    }
    is.atomic(x) && !is.object(x) && is.null(attributes(x))
  }
  fits <- attr(tradeoff_data, "fits", exact = TRUE)
  conditions <- attr(tradeoff_data, "conditions", exact = TRUE)
  metadata <- attr(tradeoff_data, "metadata", exact = TRUE)
  metadata_names <- c(
    "mode", "evaluation_order", "output_order", "warm_start_policy",
    "failure_policy", "selection_policy", "target_K", "target_weight",
    "fixed_scales"
  )
  retained_ok <- ordinary_named_list(
    fits, point_id, row_count
  ) && ordinary_named_list(
    conditions, point_id, row_count
  ) && ordinary_named_list(
    metadata, metadata_names, length(metadata_names)
  ) && ordinary_metadata_tree(metadata)
  if (!retained_ok) {
    .dpprior_visualization_abort(
      paste(
        "tradeoff_data does not retain ordinary row-aligned canonical fits,",
        "conditions, and exact metadata."
      ),
      classes = c(
        "dpprior_visualization_input_error", "dpprior_invalid_input"
      ),
      code = "tradeoff_curve_retained_evidence",
      operation = "trade-off curve visualization"
    )
  }

  abort_mismatch <- function(point, fields) {
    .dpprior_visualization_abort(
      sprintf(
        "tradeoff_data row '%s' differs from its retained canonical evidence.",
        point
      ),
      classes = c(
        "dpprior_visualization_input_error", "dpprior_invalid_input"
      ),
      code = "tradeoff_curve_science_mismatch",
      operation = "trade-off curve visualization",
      point_id = point,
      mismatched_fields = fields
    )
  }
  scalar_matches <- function(observed, expected) {
    if (is.null(expected)) return(identical(observed, NA_real_))
    identical(observed, expected)
  }
  fit_records <- vector("list", row_count)
  path_reference <- NULL

  for (index in seq_len(row_count)) {
    fit <- .dpprior_require_schema(
      fits[[index]], kind = "fit", allow_legacy = FALSE
    )
    fit_raw <- .dpprior_visualization_plain_record(unclass(fit))
    parameters <- fit_raw[["parameters", exact = TRUE]]
    target_bundle <- fit_raw[["target", exact = TRUE]]
    target_K <- target_bundle[["K", exact = TRUE]]
    target_K_used <- target_K[["used", exact = TRUE]]
    target_K_science <- list(
      mu_K = target_K_used[["mu_K", exact = TRUE]],
      var_K = target_K_used[["var_K", exact = TRUE]]
    )
    target_weight <- target_bundle[["weight", exact = TRUE]]
    target_weight_science <- list(
      metric = target_weight[["metric", exact = TRUE]],
      relation = target_weight[["relation", exact = TRUE]],
      operator = target_weight[["operator", exact = TRUE]],
      value = target_weight[["value", exact = TRUE]],
      threshold = target_weight[["threshold", exact = TRUE]],
      probability = target_weight[["probability", exact = TRUE]],
      estimand = target_weight[["estimand", exact = TRUE]],
      units = target_weight[["units", exact = TRUE]]
    )
    point_target <- target_weight_science[
      c("metric", "relation", "value", "threshold", "probability")
    ]
    achieved <- fit_raw[["achieved", exact = TRUE]]
    achieved_K <- achieved[["K", exact = TRUE]]
    achieved_weight <- achieved[["weight", exact = TRUE]]
    tradeoff <- fit_raw[["tradeoff", exact = TRUE]]
    computation <- fit_raw[["computation", exact = TRUE]]
    orders <- computation[["orders", exact = TRUE]]
    attempts <- computation[["attempts", exact = TRUE]]
    resources <- computation[["resources", exact = TRUE]]
    provenance <- fit_raw[["provenance", exact = TRUE]]
    selected_M <- orders[["M_selected", exact = TRUE]]
    effective_orders <- if (!is.null(selected_M)) {
      list(
        M_selected = selected_M,
        M_verification_used =
          orders[["M_verification_used", exact = TRUE]]
      )
    } else {
      list(
        M_selected = resources[["requested_M", exact = TRUE]],
        M_verification_used =
          resources[["requested_M_verify", exact = TRUE]]
      )
    }
    expected_point_id <- .dpprior_soft_point_id(
      fit_raw[["J", exact = TRUE]], target_K_science, point_target,
      tradeoff[["lambda", exact = TRUE]]
    )

    condition <- conditions[[index]]
    condition_classes <- if (is.null(condition)) {
      character()
    } else {
      attr(condition, "class", exact = TRUE)
    }
    condition_raw <- if (typeof(condition) == "list") {
      unclass(condition)
    } else {
      NULL
    }
    plain_condition_character <- function(value, allow_null = FALSE) {
      if (is.null(value)) return(isTRUE(allow_null))
      identical(typeof(value), "character") && !is.object(value) &&
        is.null(attributes(value)) && identical(length(value), 1L) &&
        !is.na(value) && nzchar(value)
    }
    condition_evidence_plain <- if (is.null(condition)) {
      TRUE
    } else {
      typeof(condition_raw) == "list" &&
        ordinary_named_list(condition_raw) &&
        identical(typeof(condition_classes), "character") &&
        !is.object(condition_classes) &&
        is.null(attributes(condition_classes)) &&
        length(condition_classes) >= 1L &&
        !anyNA(condition_classes) && all(nzchar(condition_classes)) &&
        plain_condition_character(
          condition_raw[["code", exact = TRUE]]
        ) &&
        plain_condition_character(
          condition_raw[["message", exact = TRUE]]
        ) &&
        plain_condition_character(
          condition_raw[["status", exact = TRUE]], allow_null = TRUE
        ) &&
        plain_condition_character(
          condition_raw[["action", exact = TRUE]], allow_null = TRUE
        )
    }
    if (!condition_evidence_plain) {
      abort_mismatch(expected_point_id, "condition_retained_evidence")
    }
    condition_fields <- if (is.null(condition)) {
      list(
        class = NA_character_, code = NA_character_, message = NA_character_
      )
    } else {
      list(
        class = condition_classes[[1L]],
        code = condition_raw[["code", exact = TRUE]],
        message = condition_raw[["message", exact = TRUE]]
      )
    }
    status <- fit_raw[["status", exact = TRUE]]
    allow_approximate <- resources[["allow_approximate_return", exact = TRUE]]
    expected_condition <- identical(status, "failed") ||
      identical(status, "infeasible") ||
      (identical(status, "approximate") && !isTRUE(allow_approximate))
    expected_outcome <- if (expected_condition) {
      "condition_retained"
    } else {
      "fit_returned"
    }
    condition_code <- condition_fields[["code", exact = TRUE]]
    allowed_codes <- switch(
      status,
      approximate = "dual_soft_approximate",
      infeasible = "dual_soft_infeasible",
      failed = c("dual_soft_failed", "dual_soft_backend_contract"),
      character()
    )
    first_class <- switch(
      condition_code,
      dual_soft_approximate = "dpprior_dual_soft_approximation_error",
      dual_soft_infeasible = "dpprior_dual_soft_infeasible_error",
      dual_soft_failed = "dpprior_dual_soft_computation_error",
      dual_soft_backend_contract =
        "dpprior_dual_soft_backend_contract_error",
      NA_character_
    )
    expected_classes <- if (is.na(first_class)) {
      character()
    } else {
      c(
        first_class, "dpprior_dual_soft_error",
        "dpprior_calibration_error", "dpprior_error", "error",
        "dpprior_condition", "condition"
      )
    }
    standard_code <- condition_code %in% c(
      "dual_soft_approximate", "dual_soft_infeasible", "dual_soft_failed"
    )
    expected_action <- if (identical(status, "approximate")) {
      "review_then_explicitly_allow_approximate"
    } else {
      "refit_or_change_controls"
    }
    condition_contract <- c(
      condition_presence = identical(!is.null(condition), expected_condition),
      condition_outcome = identical(
        raw[["outcome", exact = TRUE]][[index]], expected_outcome
      ),
      condition_class = identical(
        raw[["condition_class", exact = TRUE]][[index]],
        condition_fields[["class", exact = TRUE]]
      ),
      condition_code = identical(
        raw[["condition_code", exact = TRUE]][[index]], condition_code
      ),
      condition_message = identical(
        raw[["condition_message", exact = TRUE]][[index]],
        condition_fields[["message", exact = TRUE]]
      ),
      condition_allowed = if (expected_condition) {
        condition_code %in% allowed_codes &&
          identical(condition_classes, expected_classes)
      } else {
        is.null(condition) &&
          identical(condition_fields, list(
            class = NA_character_, code = NA_character_,
            message = NA_character_
          ))
      },
      condition_result = if (expected_condition) {
        typeof(condition_raw) == "list" &&
          identical(condition_raw[["result", exact = TRUE]], fits[[index]])
      } else {
        TRUE
      },
      condition_status = if (expected_condition && standard_code) {
        identical(condition_raw[["status", exact = TRUE]], status) &&
          identical(
            condition_raw[["action", exact = TRUE]], expected_action
          )
      } else {
        TRUE
      }
    )

    finite_parameters <- !is.null(parameters)
    fresh_wsb_tail <- if (finite_parameters && !is.null(selected_M)) {
      prob_wsb_exceeds(
        0.5, parameters[["a", exact = TRUE]],
        parameters[["b", exact = TRUE]]
      )
    } else {
      NULL
    }
    fresh_wsb_mean <- if (finite_parameters && !is.null(selected_M)) {
      mean_w1(
        parameters[["a", exact = TRUE]],
        parameters[["b", exact = TRUE]], selected_M
      )
    } else {
      NULL
    }
    expected_scalars <- list(
      a = if (is.null(parameters)) NULL else
        parameters[["a", exact = TRUE]],
      b = if (is.null(parameters)) NULL else
        parameters[["b", exact = TRUE]],
      mu_K = if (is.null(achieved_K)) NULL else
        achieved_K[["mean", exact = TRUE]],
      var_K = if (is.null(achieved_K)) NULL else
        achieved_K[["variance", exact = TRUE]],
      achieved_weight = if (is.null(achieved_weight)) NULL else
        achieved_weight[["value", exact = TRUE]],
      target_residual = tradeoff[["target_residual", exact = TRUE]],
      K_loss = tradeoff[["K_loss", exact = TRUE]],
      weight_loss = tradeoff[["weight_loss", exact = TRUE]],
      total_loss = tradeoff[["total_loss", exact = TRUE]],
      w_loss = tradeoff[["weight_loss", exact = TRUE]],
      w1_prob_gt_50 = fresh_wsb_tail,
      E_w1 = fresh_wsb_mean
    )
    scalar_ok <- vapply(names(expected_scalars), function(field) {
      scalar_matches(
        raw[[field, exact = TRUE]][[index]],
        expected_scalars[[field]]
      )
    }, logical(1L))

    current_path <- list(
      J = fit_raw[["J", exact = TRUE]],
      target_K = target_K,
      target_K_science = target_K_science,
      target_weight = target_weight,
      target_weight_science = target_weight_science,
      input_fit = provenance[["input_fit", exact = TRUE]],
      request = computation[["request", exact = TRUE]],
      used = computation[["used", exact = TRUE]],
      effective_orders = effective_orders,
      scales = tradeoff[["scales", exact = TRUE]]
    )
    if (is.null(path_reference)) path_reference <- current_path
    path_ok <- c(
      path_J = identical(current_path[["J"]], path_reference[["J"]]),
      path_target_K = identical(
        current_path[["target_K"]], path_reference[["target_K"]]
      ),
      path_target_weight = identical(
        current_path[["target_weight"]],
        path_reference[["target_weight"]]
      ),
      path_input_fit = identical(
        current_path[["input_fit"]], path_reference[["input_fit"]]
      ),
      path_request = identical(
        current_path[["request"]], path_reference[["request"]]
      ),
      path_used = identical(
        current_path[["used"]], path_reference[["used"]]
      ),
      path_orders = identical(
        current_path[["effective_orders"]],
        path_reference[["effective_orders"]]
      ),
      path_scales = identical(
        current_path[["scales"]], path_reference[["scales"]]
      )
    )
    row_ok <- c(
      point_id = identical(
        raw[["point_id", exact = TRUE]][[index]], expected_point_id
      ),
      lambda = identical(
        raw[["lambda", exact = TRUE]][[index]],
        tradeoff[["lambda", exact = TRUE]]
      ),
      mode = identical(
        raw[["mode", exact = TRUE]][[index]],
        fit_raw[["mode", exact = TRUE]]
      ) && identical(fit_raw[["mode", exact = TRUE]], "dual_soft"),
      metric = identical(
        raw[["metric", exact = TRUE]][[index]],
        target_weight_science[["metric", exact = TRUE]]
      ),
      relation = identical(
        raw[["relation", exact = TRUE]][[index]],
        target_weight_science[["relation", exact = TRUE]]
      ),
      target_value = identical(
        raw[["target_value", exact = TRUE]][[index]],
        target_weight_science[["value", exact = TRUE]]
      ),
      status = identical(
        raw[["status", exact = TRUE]][[index]], status
      ),
      usable = identical(
        raw[["usable", exact = TRUE]][[index]],
        fit_raw[["usable", exact = TRUE]]
      ),
      verified = identical(
        raw[["verified", exact = TRUE]][[index]],
        fit_raw[["verified", exact = TRUE]]
      ),
      converged = identical(
        raw[["converged", exact = TRUE]][[index]],
        identical(status, "converged")
      ),
      attempt_count = identical(
        raw[["attempt_count", exact = TRUE]][[index]],
        as.integer(length(attempts))
      ),
      selected_method = identical(
        raw[["selected_method", exact = TRUE]][[index]],
        provenance[["selected_method", exact = TRUE]]
      ),
      fit_method = identical(
        fit_raw[["method", exact = TRUE]], "dual-soft"
      ),
      condition_contract,
      scalar_ok,
      path_ok
    )
    if (!all(row_ok)) {
      abort_mismatch(expected_point_id, names(row_ok)[!row_ok])
    }

    fit_records[[index]] <- list(
      fit = fit,
      raw = fit_raw,
      parameters = parameters,
      resources = resources,
      expected_point_id = expected_point_id
    )
  }

  metadata_weight <- metadata[["target_weight", exact = TRUE]]
  weight_names <- c(
    "metric", "relation", "operator", "value", "threshold", "probability",
    "estimand", "units", "raw_input"
  )
  weight_structure_ok <- ordinary_named_list(
    metadata_weight, weight_names, length(weight_names)
  )
  normalized_weight <- if (weight_structure_ok &&
                           ordinary_named_list(
                             metadata_weight[["raw_input", exact = TRUE]]
                           )) {
    tryCatch(
      .dpprior_soft_normalize_target(
        metadata_weight[["raw_input", exact = TRUE]]
      ),
      error = function(condition) condition
    )
  } else {
    NULL
  }
  weight_science_ok <- weight_structure_ok && all(vapply(
    names(path_reference[["target_weight_science"]]), function(field) {
      identical(
        metadata_weight[[field, exact = TRUE]],
        path_reference[["target_weight_science"]][[field, exact = TRUE]]
      )
    }, logical(1L)
  ))
  expected_fixed_scales <- list(
    K_mean = max(abs(
      path_reference[["target_K_science"]][["mu_K", exact = TRUE]]
    ), 1),
    K_variance = max(abs(
      path_reference[["target_K_science"]][["var_K", exact = TRUE]]
    ), 1),
    weight = 1
  )
  metadata_ok <- c(
    metadata_mode = identical(
      metadata[["mode", exact = TRUE]], "soft_tradeoff"
    ),
    metadata_evaluation_order = identical(
      metadata[["evaluation_order", exact = TRUE]],
      sort(lambda, decreasing = TRUE)
    ),
    metadata_output_order = identical(
      metadata[["output_order", exact = TRUE]],
      "ascending lambda then point_id"
    ),
    metadata_warm_start_policy = identical(
      metadata[["warm_start_policy", exact = TRUE]],
      "descending lambda; verified usable predecessor only"
    ),
    metadata_failure_policy = identical(
      metadata[["failure_policy", exact = TRUE]],
      "all requested points retained"
    ),
    metadata_selection_policy = identical(
      metadata[["selection_policy", exact = TRUE]],
      "descriptive path; no best lambda is selected"
    ),
    metadata_target_K = identical(
      metadata[["target_K", exact = TRUE]],
      path_reference[["target_K_science"]]
    ),
    metadata_target_weight = weight_science_ok &&
      !inherits(normalized_weight, "condition") &&
      identical(normalized_weight, metadata_weight),
    metadata_fixed_scales = identical(
      metadata[["fixed_scales", exact = TRUE]], expected_fixed_scales
    ) && identical(
      path_reference[["scales"]],
      list(
        K = list(
          mean = expected_fixed_scales[["K_mean", exact = TRUE]],
          variance = expected_fixed_scales[["K_variance", exact = TRUE]]
        ),
        weight = expected_fixed_scales[["weight", exact = TRUE]]
      )
    )
  )
  if (!all(metadata_ok)) {
    abort_mismatch("<metadata>", names(metadata_ok)[!metadata_ok])
  }

  evaluation_index <- order(lambda, decreasing = TRUE)
  warm_source <- "K-only input fit"
  source_parameters <- path_reference[["input_fit", exact = TRUE]][[
    "parameters", exact = TRUE
  ]]
  for (index in evaluation_index) {
    record <- fit_records[[index]]
    parameters <- record[["parameters", exact = TRUE]]
    resources <- record[["resources", exact = TRUE]]
    fit_raw <- record[["raw", exact = TRUE]]
    expected_resource <- if (identical(
      raw[["lambda", exact = TRUE]][[index]], 1
    ) || is.null(parameters)) {
      NULL
    } else {
      source_parameters[c("a", "b")]
    }
    warm_ok <- c(
      warm_start_from = identical(
        raw[["warm_start_from", exact = TRUE]][[index]], warm_source
      ),
      warm_start_parameters = identical(
        resources[["warm_start", exact = TRUE]], expected_resource
      )
    )
    if (!all(warm_ok)) {
      abort_mismatch(
        record[["expected_point_id", exact = TRUE]],
        names(warm_ok)[!warm_ok]
      )
    }
    if (isTRUE(fit_raw[["usable", exact = TRUE]]) &&
        isTRUE(fit_raw[["verified", exact = TRUE]])) {
      warm_source <- record[["expected_point_id", exact = TRUE]]
      source_parameters <- parameters
    }
  }

  available <- is.finite(value)
  if (!any(available)) {
    .dpprior_visualization_abort(
      sprintf("metric '%s' is unavailable at every retained curve point.", metric),
      classes = c(
        "dpprior_s3_unavailable_error", "dpprior_visualization_data_error"
      ),
      code = "tradeoff_metric_unavailable",
      operation = "trade-off curve visualization",
      metric = metric
    )
  }
  plot_data <- data.frame(
    point_id = raw[["point_id", exact = TRUE]],
    lambda = lambda,
    value = value,
    status = raw[["status", exact = TRUE]],
    usable = raw[["usable", exact = TRUE]],
    verified = raw[["verified", exact = TRUE]],
    outcome = raw[["outcome", exact = TRUE]],
    condition_code = raw[["condition_code", exact = TRUE]],
    available = available,
    stringsAsFactors = FALSE
  )
  available_index <- which(available)
  line_data <- plot_data[available_index, , drop = FALSE]
  line_data[["segment"]] <- cumsum(c(
    TRUE, diff(available_index) != 1L
  ))
  list(
    data = plot_data,
    line_data = line_data,
    requested_points = nrow(plot_data),
    available_points = nrow(line_data),
    unavailable_points = sum(!available)
  )
}

.dpprior_assert_fit <- function(fit) {
  .dpprior_visualization_fit_view(
    fit, operation = "fit validation", require_candidate = TRUE
  )
  invisible(TRUE)
}

.dpprior_extract_abJ <- function(fit) {
  view <- .dpprior_visualization_fit_view(
    fit, operation = "parameter extraction", require_candidate = TRUE
  )
  list(
    a = as.numeric(view[["parameters", exact = TRUE]][["a", exact = TRUE]]),
    b = as.numeric(view[["parameters", exact = TRUE]][["b", exact = TRUE]]),
    J = as.integer(view[["core", exact = TRUE]][["J", exact = TRUE]])
  )
}

# =============================================================================
# Theme and Colors
# =============================================================================

#' Publication-Quality Theme for DPprior Plots
#'
#' @param base_size Numeric; base font size (default: 11).
#' @param base_family Character; base font family (default: "").
#' @return A ggplot2 theme object.
#'
#' @family visualization
#'
#' @export
theme_DPprior <- function(base_size = 11, base_family = "") {
  .dpprior_require_ggplot2()
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey40",
                                            margin = ggplot2::margin(b = 8)),
      plot.caption = ggplot2::element_text(color = "grey30", hjust = 0),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(face = "plain"),
      axis.text = ggplot2::element_text(color = "grey20"),
      legend.position = "none"
    )
}

#' DPprior Color Palette
#' @return A named list with color values.
#'
#' @family visualization
#'
#' @export
DPprior_colors <- function() {
  list(
    primary = "#4682B4",      # steelblue
    secondary = "#41B6C4",    # teal
    accent = "#2C7FB8",       # blue
    ink = "#1F2D3A",          # dark
    warning = "#E34A33",      # red-orange
    shade = "#D9D9D9",        # light gray

    # Comparison colors (K-only vs Dual)
    k_only = "#4682B4",       # steelblue (baseline)
    dual = "#E67E22"          # carrot orange (new solution)
  )
}

# =============================================================================
# Distribution Helpers
# =============================================================================

.dpprior_density_w1 <- function(x, a, b) {
  if (exists("density_w1", mode = "function")) {
    return(density_w1(x, a, b))
  }
  # Closed-form fallback: p(x|a,b) = a*b^a / ((1-x)[b - log(1-x)]^(a+1))
  x <- pmin(pmax(x, .Machine$double.eps), 1 - .Machine$double.eps)
  denom <- (1 - x) * (b - log(1 - x))^(a + 1)
  (a * b^a) / denom
}

.dpprior_prob_w1_exceeds <- function(t, a, b) {
  if (exists("prob_w1_exceeds", mode = "function")) {
    return(prob_w1_exceeds(t, a, b))
  }
  # Closed-form W_SB tail: (b / (b - log(1-t)))^a
  t <- pmin(pmax(t, .Machine$double.eps), 1 - .Machine$double.eps)
  (b / (b - log(1 - t)))^a
}

.dpprior_quantile_w1 <- function(u, a, b) {
  if (exists("quantile_w1", mode = "function")) {
    return(quantile_w1(u, a, b))
  }
  u <- pmin(pmax(u, .Machine$double.eps), 1 - .Machine$double.eps)
  exponent <- b * (1 - (1 - u)^(-1 / a))
  1 - exp(exponent)
}

.dpprior_mean_w1 <- function(a, b) {
  if (exists("mean_w1", mode = "function")) {
    return(mean_w1(a, b))
  }
  f <- function(w) w * .dpprior_density_w1(w, a, b)
  tryCatch(
    stats::integrate(f, 0, 1, rel.tol = 1e-6)$value,
    error = function(e) NA_real_
  )
}

#' Canonical K PMF Computation for Visualization
#' @keywords internal
.dpprior_get_K_pmf <- function(J, a, b, M = NULL) {
  M <- .dpprior_coalesce(M, .dpprior_get_quad_nodes_default())

  if (!exists(".get_K_pmf_support", mode = "function")) {
    stop(.dpprior_new_condition(
      paste(
        "canonical K_J diagnostic PMF backend is unavailable; no",
        "approximate distribution was substituted"
      ),
      classes = c(
        "dpprior_visualization_data_error", "dpprior_numerical_error",
        "dpprior_error", "error"
      ),
      method = "canonical_diagnostics_pmf",
      code = "canonical_pmf_backend_unavailable"
    ))
  }

  obj <- tryCatch(
    .get_K_pmf_support(J, a, b, M),
    error = function(error) error
  )
  if (inherits(obj, "error")) {
    stop(.dpprior_new_condition(
      paste(
        "canonical K_J PMF computation failed; visualization data are",
        "unavailable and no uniform or negative-binomial PMF was substituted: ",
        conditionMessage(obj)
      ),
      classes = c(
        "dpprior_visualization_data_error", "dpprior_numerical_error",
        "dpprior_error", "error"
      ),
      method = "canonical_diagnostics_pmf",
      parent = obj,
      result = obj$result %||% NULL,
      code = "canonical_pmf_failed"
    ))
  }

  list(
    k = obj$support,
    pmf = obj$pmf,
    method = obj$metadata$engine %||% "canonical_diagnostics_pmf",
    status = obj$metadata$status %||% obj$status %||% "approximate",
    verified = identical(
      obj$metadata$status %||% obj$status, "converged"
    ) && isTRUE(obj$metadata$verification$passed),
    verification = obj$metadata$verification %||% NULL,
    provenance = list(
      requested_method = "canonical_diagnostic_pmf",
      selected_method =
        obj$metadata$engine %||% "canonical_diagnostics_pmf",
      is_fallback = FALSE,
      marginal_metadata = obj$metadata
    )
  )
}

# =============================================================================
# Summary Statistics Computation
# =============================================================================

#' Compute Summary Statistics
#' @keywords internal
.dpprior_compute_summary <- function(fit = NULL, a = NULL, b = NULL, J = NULL,
                                     ci_level = 0.95) {
  view <- NULL
  if (!is.null(fit)) {
    view <- .dpprior_visualization_fit_view(
      fit, operation = "summary visualization", require_candidate = TRUE
    )
    parameters <- view[["parameters", exact = TRUE]]
    a <- as.numeric(parameters[["a", exact = TRUE]])
    b <- as.numeric(parameters[["b", exact = TRUE]])
    J <- as.integer(view[["core", exact = TRUE]][["J", exact = TRUE]])
  } else {
    if (is.null(a) || is.null(b)) {
      stop("Either 'fit' or both 'a' and 'b' must be provided", call. = FALSE)
    }
    J <- .dpprior_coalesce(J, 50L)
  }

  # Alpha statistics
  alpha_mean <- a / b
  alpha_cv <- 1 / sqrt(a)
  alpha_ci <- stats::qgamma(
    c((1 - ci_level) / 2, (1 + ci_level) / 2),
    shape = a, rate = b
  )

  # K statistics from the canonical status-aware PMF backend
  M <- if (is.null(view)) NULL else view[["M", exact = TRUE]]
  Kp <- .dpprior_get_K_pmf(J, a, b, M = M)
  k <- Kp$k; pmf <- Kp$pmf
  cdf <- cumsum(pmf)

  K_mean <- sum(k * pmf)
  K_var <- max(0, sum((k^2) * pmf) - K_mean^2)
  K_mode <- k[which.max(pmf)]

  # w1 statistics
  w1_mean <- .dpprior_mean_w1(a, b)
  w1_median <- .dpprior_quantile_w1(0.5, a, b)
  w1_p50 <- .dpprior_prob_w1_exceeds(0.5, a, b)
  w1_p90 <- .dpprior_prob_w1_exceeds(0.9, a, b)

  if (is.null(view)) {
    target_mu <- target_var <- NA_real_
    achieved_mu <- K_mean
    achieved_var <- K_var
    target_estimand <- "not_supplied"
    achieved_estimand <- "fresh_visualization_PMF"
    core <- list(
      mode = "direct_parameters", method = NA_character_,
      status = Kp$status %||% "approximate", usable = TRUE,
      verified = Kp$verified %||% FALSE, message = "direct parameter input"
    )
  } else {
    target_K <- view[["target", exact = TRUE]][["K", exact = TRUE]]
    implied <- target_K[["implied", exact = TRUE]]
    achieved_K <- view[["achieved", exact = TRUE]][["K", exact = TRUE]]
    target_mu <- implied[["mean", exact = TRUE]]
    target_var <- implied[["variance", exact = TRUE]]
    achieved_mu <- achieved_K[["mean", exact = TRUE]]
    achieved_var <- achieved_K[["variance", exact = TRUE]]
    target_estimand <- "canonical_target_K_implied_moments"
    achieved_estimand <- achieved_K[["estimand", exact = TRUE]]
    core <- view[["core", exact = TRUE]]
  }

  list(
    a = a, b = b, J = J,
    alpha = list(mean = alpha_mean, cv = alpha_cv, ci = alpha_ci),
    K = list(k = k, pmf = pmf, cdf = cdf, mean = K_mean, var = K_var,
             mode = K_mode, pmf_method = Kp$method,
             status = Kp$status %||% "approximate",
             verified = Kp$verified %||% FALSE,
             verification = Kp$verification %||% NULL,
             provenance = Kp$provenance %||% NULL),
    w1 = list(
      estimand = "W_SB",
      label = "First size-biased DP weight",
      mean = w1_mean,
      median = w1_median,
      p_gt_50 = w1_p50,
      p_gt_90 = w1_p90,
      method = "closed_form_beta_gamma_mixture"
    ),
    target = list(
      mu_K = target_mu, var_K = target_var, estimand = target_estimand
    ),
    achieved = list(
      mu_K = achieved_mu, var_K = achieved_var,
      estimand = achieved_estimand
    ),
    method = core[["method", exact = TRUE]],
    mode = core[["mode", exact = TRUE]],
    status = core[["status", exact = TRUE]],
    usable = core[["usable", exact = TRUE]],
    verified = core[["verified", exact = TRUE]],
    message = core[["message", exact = TRUE]],
    constraint = if (is.null(view)) NULL else {
      view[["constraint", exact = TRUE]]
    },
    tradeoff = if (is.null(view)) NULL else view[["tradeoff", exact = TRUE]],
    legacy = if (is.null(view)) NULL else view[["legacy", exact = TRUE]]
  )
}

.dpprior_K_provenance_text <- function(summary) {
  sprintf(
    "status=%s; verified=%s; method=%s",
    summary$K$status %||% "unavailable",
    if (isTRUE(summary$K$verified)) "yes" else "no",
    summary$K$pmf_method %||% "unavailable"
  )
}

# =============================================================================
# Individual Plot Functions
# =============================================================================

#' Plot Prior Density of Alpha
#'
#' @param fit A DPprior_fit object, or NULL if a and b provided directly.
#' @param a Numeric; shape parameter (used if fit is NULL).
#' @param b Numeric; rate parameter (used if fit is NULL).
#' @param engine "ggplot2" (default) or "base".
#' @param base_size Base font size.
#' @param ci_level Credible interval level (default 0.95).
#' @param n_grid Number of grid points.
#' @param show If TRUE, print the plot.
#' @return A ggplot object or invisible(NULL) for base.
#'
#' @examples
#' # From fit object
#' fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8)
#' plot_alpha_prior(fit)
#'
#' # Direct parameter specification
#' plot_alpha_prior(a = 1.6, b = 1.2)
#'
#' @seealso \code{\link{DPprior_fit}} for fitting, \code{\link{plot.DPprior_fit}} for S3 plot method
#'
#' @family visualization
#'
#' @export
plot_alpha_prior <- function(fit = NULL, a = NULL, b = NULL,
                             engine = c("ggplot2", "base"),
                             base_size = 11,
                             ci_level = 0.95,
                             n_grid = 500,
                             show = TRUE) {
  engine <- match.arg(engine)
  view <- NULL

  # Extract parameters from fit or use direct values
  if (!is.null(fit)) {
    view <- .dpprior_visualization_fit_view(
      fit, operation = "alpha-prior visualization", require_candidate = TRUE
    )
    parameters <- view[["parameters", exact = TRUE]]
    a <- parameters[["a", exact = TRUE]]
    b <- parameters[["b", exact = TRUE]]
  } else {
    if (is.null(a) || is.null(b)) {
      stop("Either 'fit' or both 'a' and 'b' must be provided", call. = FALSE)
    }
  }

  # Compute statistics
  alpha_mean <- a / b
  alpha_cv <- 1 / sqrt(a)
  ci <- stats::qgamma(c((1 - ci_level) / 2, (1 + ci_level) / 2), shape = a, rate = b)

  # Grid
  x_max <- max(stats::qgamma(0.999, shape = a, rate = b), alpha_mean * 3)
  x <- seq(0.001, x_max, length.out = n_grid)
  dens <- stats::dgamma(x, shape = a, rate = b)
  df <- data.frame(x = x, density = dens)

  # Base R fallback
  if (engine == "base" || !.dpprior_has_ggplot2()) {
    if (!isTRUE(show)) {
      return(invisible(NULL))
    }
    .dpprior_base_plot_alpha(
      df, alpha_mean, ci, alpha_cv, a, b,
      context = .dpprior_visualization_fit_caption(
        view, "Gamma concentration hyperprior"
      )
    )
    return(invisible(NULL))
  }

  .dpprior_require_ggplot2()
  colors <- DPprior_colors()

  # Build subtitle with stats + Gamma parameters (two lines)
  subtitle_text <- sprintf(
    "E[alpha]=%.2f,  CV(alpha)=%.2f,  %d%% CI=[%.2f, %.2f]\nGamma(a=%.4f, b=%.3f)",
    alpha_mean, alpha_cv, round(100 * ci_level), ci[1], ci[2], a, b
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = density)) +
    ggplot2::geom_area(fill = colors$primary, alpha = 0.3) +
    ggplot2::geom_line(color = colors$primary, linewidth = 1) +
    ggplot2::annotate("rect", xmin = ci[1], xmax = ci[2],
                      ymin = -Inf, ymax = Inf, fill = colors$shade, alpha = 0.3) +
    ggplot2::geom_vline(xintercept = alpha_mean, linetype = "dashed",
                        color = colors$accent, linewidth = 0.9) +
    ggplot2::geom_vline(xintercept = ci, linetype = "dotted",
                        color = "gray50", linewidth = 0.5) +
    theme_DPprior(base_size = base_size) +
    ggplot2::labs(
      x = expression(Concentration~parameter~alpha),
      y = "Density",
      title = expression(bold("(A)")~Prior~on~alpha),
      subtitle = subtitle_text,
      caption = .dpprior_visualization_fit_caption(
        view, "Gamma concentration hyperprior"
      )
    )

  if (isTRUE(show)) print(p)
  p
}


#' Plot Prior PMF of K_J
#'
#' @param fit A DPprior_fit object, or NULL if J, a, b provided directly.
#' @param J Integer; sample size (used if fit is NULL).
#' @param a Numeric; shape parameter (used if fit is NULL).
#' @param b Numeric; rate parameter (used if fit is NULL).
#' @param engine "ggplot2" (default) or "base".
#' @param base_size Base font size.
#' @param max_k Maximum k to display. If NULL, auto-determined by CDF >= 0.999.
#' @param show_cdf If TRUE, overlay CDF line. Default is FALSE.
#' @param show If TRUE, print the plot.
#' @return A ggplot object or invisible(NULL) for base.
#'
#' @examples
#' fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8)
#' plot_K_prior(fit)
#'
#' plot_K_prior(J = 50, a = 1.6, b = 1.2)
#'
#' @seealso \code{\link{DPprior_fit}} for fitting, \code{\link{plot.DPprior_fit}} for S3 plot method
#'
#' @family visualization
#'
#' @export
plot_K_prior <- function(fit = NULL, J = NULL, a = NULL, b = NULL,
                         engine = c("ggplot2", "base"),
                         base_size = 11,
                         max_k = NULL,
                         show_cdf = FALSE,
                         show = TRUE) {
  engine <- match.arg(engine)

  # Extract parameters
  if (!is.null(fit)) {
    .dpprior_assert_fit(fit)
    abJ <- .dpprior_extract_abJ(fit)
    J <- abJ$J; a <- abJ$a; b <- abJ$b
  } else {
    if (is.null(J) || is.null(a) || is.null(b)) {
      stop("Either 'fit' or all of 'J', 'a', 'b' must be provided", call. = FALSE)
    }
  }

  # Get a canonical status-aware summary
  summ <- .dpprior_compute_summary(fit = fit, a = a, b = b, J = J)
  k <- summ$K$k
  pmf <- summ$K$pmf
  cdf <- summ$K$cdf
  K_mean <- summ$K$mean
  K_var <- summ$K$var
  K_mode <- summ$K$mode
  target_mu <- summ$target$mu_K
  achieved_mu <- summ$achieved$mu_K

  # Auto max_k: smallest k with CDF >= 0.999
  if (is.null(max_k)) {
    cut_idx <- which(cdf >= 0.999)[1]
    if (is.na(cut_idx)) cut_idx <- length(k)
    max_k <- min(k[cut_idx], J, 100)
    max_k <- max(max_k, K_mode + 10, 15)
  }
  max_k <- min(max_k, J)

  # Truncate
  idx <- k <= max_k
  df <- data.frame(k = k[idx], pmf = pmf[idx], cdf = cdf[idx])

  # Base R fallback
  if (engine == "base" || !.dpprior_has_ggplot2()) {
    if (!isTRUE(show)) {
      return(invisible(NULL))
    }
    .dpprior_base_plot_K(
      df, target_mu, achieved_mu, K_mean, K_var, K_mode, a, b,
      status = summ$K$status, verified = summ$K$verified,
      method = summ$K$pmf_method,
      context = if (is.null(fit)) NULL else sprintf(
        paste0(
          "fit mode=%s; status=%s; verified=%s; ",
          "target estimand=%s; achieved estimand=%s%s"
        ),
        summ$mode, summ$status,
        if (isTRUE(summ$verified)) "yes" else "no",
        summ$target$estimand, summ$achieved$estimand,
        if (identical(summ$mode, "dual_legacy")) {
          "; deprecated legacy dual fit"
        } else {
          ""
        }
      )
    )
    return(invisible(NULL))
  }

  .dpprior_require_ggplot2()
  colors <- DPprior_colors()
  y_max <- max(df$pmf)

  # Build subtitle with stats + Gamma parameters (two lines)
  subtitle_text <- sprintf("E[K]=%.2f,  Var(K)=%.2f,  Mode=%d\nGamma(a=%.4f, b=%.3f)",
                           K_mean, K_var, K_mode, a, b)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = k, y = pmf)) +
    ggplot2::geom_col(fill = colors$primary, alpha = 0.7, width = 0.8) +
    theme_DPprior(base_size = base_size)

  # CDF overlay (scaled to PMF axis)
  if (show_cdf) {
    p <- p +
      ggplot2::geom_line(ggplot2::aes(y = cdf * y_max),
                         color = colors$ink, linewidth = 0.8) +
      ggplot2::geom_point(ggplot2::aes(y = cdf * y_max),
                          color = colors$ink, size = 1)
  }

  # Target mean line (if available)
  if (is.finite(target_mu)) {
    p <- p + ggplot2::geom_vline(xintercept = target_mu, linetype = "dashed",
                                 color = colors$warning, linewidth = 0.8)
  }

  # Achieved mean line
  p <- p + ggplot2::geom_vline(xintercept = achieved_mu, linetype = "solid",
                               color = colors$accent, linewidth = 0.8)

  caption_text <- sprintf(
    "PMF status: %s; verified: %s; method: %s",
    summ$K$status,
    if (isTRUE(summ$K$verified)) "yes" else "no",
    summ$K$pmf_method
  )
  if (!is.null(fit)) {
    fit_caption <- sprintf(
      paste0(
        "fit mode=%s; status=%s; verified=%s; ",
        "target estimand=%s; achieved estimand=%s"
      ),
      summ$mode, summ$status,
      if (isTRUE(summ$verified)) "yes" else "no",
      summ$target$estimand, summ$achieved$estimand
    )
    if (identical(summ$mode, "dual_legacy")) {
      fit_caption <- paste0(
        fit_caption,
        "; deprecated legacy dual fit (descriptive single-fit view only)"
      )
    }
    caption_text <- paste(fit_caption, caption_text, sep = "\n")
  }
  if (show_cdf) {
    caption_text <- paste(
      caption_text, "Line: CDF (scaled to PMF axis)", sep = "\n"
    )
  }
  p <- p + ggplot2::labs(
    x = expression(Number~of~clusters~K[J]),
    y = "Probability mass",
    title = expression(bold("(B)")~Prior~PMF~of~K[J]),
    subtitle = subtitle_text,
    caption = caption_text
  )

  if (isTRUE(show)) print(p)
  p
}


#' Plot Prior Density of the First Size-Biased DP Weight
#'
#' @param fit A DPprior_fit object, or NULL if a, b provided directly.
#' @param a Numeric; shape parameter (used if fit is NULL).
#' @param b Numeric; rate parameter (used if fit is NULL).
#' @param engine "ggplot2" (default) or "base".
#' @param base_size Base font size.
#' @param thresholds Strictly increasing numeric vector of exactly two tail
#'   thresholds in \eqn{(0,1)} for the first size-biased DP weight
#'   (default: \code{c(0.5, 0.9)}).
#' @param n_grid Number of grid points.
#' @param show If TRUE, print the plot.
#' @return A ggplot object or invisible(NULL) for base.
#'
#' @examples
#' fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8)
#' plot_w1_prior(fit)
#'
#' plot_w1_prior(a = 1.6, b = 1.2)
#'
#' @seealso \code{\link{DPprior_fit}} for fitting, \code{\link{plot.DPprior_fit}} for S3 plot method
#'
#' @family visualization
#'
#' @export
plot_w1_prior <- function(fit = NULL, a = NULL, b = NULL,
                          engine = c("ggplot2", "base"),
                          base_size = 11,
                          thresholds = c(0.5, 0.9),
                          n_grid = 500,
                          show = TRUE) {
  engine <- match.arg(engine)
  view <- NULL
  thresholds <- .dpprior_validate_probability(
    thresholds, "thresholds", scalar = FALSE, open = TRUE
  )
  if (length(thresholds) != 2L) {
    .dpprior_abort_invalid(
      "thresholds must contain exactly two probabilities",
      c("dpprior_visualization_input_error", "dpprior_length_error"),
      "thresholds", thresholds, "numeric vector of length 2", "length"
    )
  }
  if (thresholds[[1L]] >= thresholds[[2L]]) {
    .dpprior_abort_invalid(
      "thresholds must be strictly increasing",
      c("dpprior_visualization_input_error", "dpprior_bounds_error"),
      "thresholds", thresholds, "thresholds[1] < thresholds[2]",
      "not_increasing"
    )
  }

  # Extract parameters
  if (!is.null(fit)) {
    view <- .dpprior_visualization_fit_view(
      fit, operation = "size-biased-weight visualization",
      require_candidate = TRUE
    )
    parameters <- view[["parameters", exact = TRUE]]
    a <- parameters[["a", exact = TRUE]]
    b <- parameters[["b", exact = TRUE]]
  } else {
    if (is.null(a) || is.null(b)) {
      stop("Either 'fit' or both 'a' and 'b' must be provided", call. = FALSE)
    }
  }

  # Grid
  x <- seq(1e-6, 1 - 1e-6, length.out = n_grid)
  dens <- .dpprior_density_w1(x, a, b)
  dens[!is.finite(dens)] <- NA

  df <- data.frame(x = x, density = dens)
  df_shade <- df[df$x >= thresholds[1], , drop = FALSE]

  # Statistics
  mean_w <- .dpprior_mean_w1(a, b)
  median_w <- .dpprior_quantile_w1(0.5, a, b)
  p_gt_first <- .dpprior_prob_w1_exceeds(thresholds[1], a, b)
  p_gt_second <- .dpprior_prob_w1_exceeds(thresholds[2], a, b)

  # Base R fallback
  if (engine == "base" || !.dpprior_has_ggplot2()) {
    if (!isTRUE(show)) {
      return(invisible(NULL))
    }
    .dpprior_base_plot_w1(
      df, thresholds, mean_w, median_w, p_gt_first, p_gt_second, a, b,
      context = .dpprior_visualization_fit_caption(
        view, "W_SB (first size-biased DP weight)"
      )
    )
    return(invisible(NULL))
  }

  .dpprior_require_ggplot2()
  colors <- DPprior_colors()

  y_max <- max(dens, na.rm = TRUE)
  if (!is.finite(y_max) || y_max > 20) y_max <- 20

  # Build subtitle with stats + Gamma parameters (two lines)
  subtitle_text <- sprintf(
    paste0(
      "E[W_SB]=%.3f,  Median=%.3f,  ",
      "P(W_SB>%.3g)=%.1f%%,  P(W_SB>%.3g)=%.1f%%\n",
      "Gamma(a=%.4f, b=%.3f)"
    ),
    mean_w, median_w,
    thresholds[[1L]], 100 * p_gt_first,
    thresholds[[2L]], 100 * p_gt_second,
    a, b
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = density)) +
    ggplot2::geom_area(data = df_shade, fill = colors$secondary, alpha = 0.3) +
    ggplot2::geom_line(color = colors$primary, linewidth = 1, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = thresholds, linetype = "dashed",
                        color = colors$warning, linewidth = 0.7) +
    theme_DPprior(base_size = base_size) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, y_max * 1.15)) +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.25)) +
    ggplot2::labs(
      x = expression(First~size-biased~DP~weight~W[SB]),
      y = "Density",
      title = expression(bold("(C)")~Prior~Density~of~W[SB]),
      subtitle = subtitle_text,
      caption = .dpprior_visualization_fit_caption(
        view, "W_SB (first size-biased DP weight)"
      )
    )

  if (isTRUE(show)) print(p)
  p
}


#' Summary Table Plot
#' @keywords internal
.dpprior_summary_table_plot <- function(fit, base_size = 11, ci_level = 0.95) {
  .dpprior_require_ggplot2()

  summ <- .dpprior_compute_summary(fit, ci_level = ci_level)
  colors <- DPprior_colors()

  # Build a descriptive table from canonical fields.  Legacy results remain
  # plotable as single fits, but are labelled deprecated and are never used as
  # a source of comparison lineage.
  metrics <- c(
    "J", "Mode", "Method", "Status", "Verified", "Gamma prior",
    "E[alpha]", "CV(alpha)", sprintf("%d%% CI(alpha)", round(100*ci_level)),
    "Target E[K]", "Target estimand", "Achieved E[K]", "Achieved Var(K)",
    "Achieved estimand", "E[W_SB]", "P(W_SB > 0.5)", "P(W_SB > 0.9)"
  )
  values <- c(
    sprintf("%d", summ$J),
    as.character(summ$mode),
    as.character(summ$method),
    as.character(summ$status),
    if (isTRUE(summ$verified)) "yes" else "no",
    sprintf("Gamma(%.4f, %.3f)", summ$a, summ$b),
    sprintf("%.3f", summ$alpha$mean),
    sprintf("%.3f", summ$alpha$cv),
    sprintf("[%.2f, %.2f]", summ$alpha$ci[1], summ$alpha$ci[2]),
    if (is.finite(summ$target$mu_K)) sprintf("%.2f", summ$target$mu_K) else "NA",
    as.character(summ$target$estimand),
    sprintf("%.3f", summ$achieved$mu_K),
    sprintf("%.3f", summ$achieved$var_K),
    as.character(summ$achieved$estimand),
    sprintf("%.3f", summ$w1$mean),
    sprintf("%.1f%%", 100 * summ$w1$p_gt_50),
    sprintf("%.1f%%", 100 * summ$w1$p_gt_90)
  )
  if (identical(summ$mode, "dual_legacy")) {
    metrics <- c(metrics, "Legacy contract")
    values <- c(values, "deprecated; descriptive single-fit view only")
  }
  tbl <- data.frame(
    Metric = metrics, Value = values, stringsAsFactors = FALSE
  )

  tbl$row <- seq_len(nrow(tbl))

  p <- ggplot2::ggplot(tbl, ggplot2::aes(y = -row)) +
    ggplot2::geom_text(ggplot2::aes(x = 0, label = Metric),
                       hjust = 0, fontface = "bold", color = colors$ink, size = 3.2) +
    ggplot2::geom_text(ggplot2::aes(x = 1, label = Value),
                       hjust = 1, color = colors$ink, size = 3.2) +
    ggplot2::scale_x_continuous(limits = c(-0.05, 1.05)) +
    ggplot2::scale_y_continuous(limits = c(-nrow(tbl) - 2, 0)) +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    ) +
    ggplot2::labs(title = expression(bold("(D)")~Summary~Statistics))

  p
}


#' Create Dashboard using gtable (no patchwork needed)
#' @keywords internal
.dpprior_dashboard_gtable <- function(p1, p2, p3, p4, title = NULL) {
  if (!requireNamespace("gtable", quietly = TRUE) ||
      !requireNamespace("grid", quietly = TRUE)) {
    return(NULL)
  }

  grobs <- matrix(list(
    ggplot2::ggplotGrob(p1), ggplot2::ggplotGrob(p2),
    ggplot2::ggplotGrob(p3), ggplot2::ggplotGrob(p4)
  ), nrow = 2, byrow = TRUE)

  g <- gtable::gtable_matrix(
    name = "dpprior_dashboard",
    grobs = grobs,
    widths = grid::unit(c(1, 1), "null"),
    heights = grid::unit(c(1, 1), "null")
  )

  # Add title if provided
  if (!is.null(title) && nzchar(title)) {
    title_grob <- grid::textGrob(
      title,
      gp = grid::gpar(fontsize = 14, fontface = "bold"),
      vjust = 0.5
    )
    g <- gtable::gtable_add_rows(g, heights = grid::unit(1.5, "lines"), pos = 0)
    g <- gtable::gtable_add_grob(g, title_grob, t = 1, l = 1, r = 2)
  }

  g
}


#' 4-Panel Prior Dashboard
#'
#' @param fit A DPprior_fit object.
#' @param engine "ggplot2" (default) or "base".
#' @param base_size Base font size.
#' @param ci_level Credible interval level.
#' @param title Optional overall title for the dashboard.
#' @param show If TRUE, draw the dashboard.
#' @return A gtable grob (for ggplot2) or invisible(NULL) for base.
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso \code{\link{DPprior_fit}} for fitting, \code{\link{plot.DPprior_fit}} for S3 plot method
#'
#' @family visualization
#'
#' @export
plot_prior_dashboard <- function(fit,
                                 engine = c("ggplot2", "base"),
                                 base_size = 11,
                                 ci_level = 0.95,
                                 title = NULL,
                                 show = TRUE) {
  engine <- match.arg(engine)
  .dpprior_assert_fit(fit)

  # Base R fallback
  if (engine == "base" || !.dpprior_has_ggplot2()) {
    if (!isTRUE(show)) {
      return(invisible(NULL))
    }
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)

    if (!is.null(title) && nzchar(title)) {
      graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
    } else {
      graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    }

    plot_alpha_prior(fit, engine = "base", ci_level = ci_level)
    plot_K_prior(fit, engine = "base")
    plot_w1_prior(fit, engine = "base")
    .dpprior_base_summary_panel(fit, ci_level = ci_level)

    if (!is.null(title) && nzchar(title)) {
      graphics::mtext(title, outer = TRUE, cex = 1.2, font = 2)
    }

    return(invisible(NULL))
  }

  .dpprior_require_ggplot2()

  # Create individual plots (don't show yet)
  p_alpha <- plot_alpha_prior(fit, engine = "ggplot2", base_size = base_size,
                              ci_level = ci_level, show = FALSE)
  p_K <- plot_K_prior(fit, engine = "ggplot2", base_size = base_size, show = FALSE)
  p_w1 <- plot_w1_prior(fit, engine = "ggplot2", base_size = base_size, show = FALSE)
  p_tbl <- .dpprior_summary_table_plot(fit, base_size = base_size, ci_level = ci_level)

  # Try gtable dashboard
  g <- .dpprior_dashboard_gtable(p_alpha, p_K, p_w1, p_tbl, title = title)

  if (is.null(g)) {
    # Fallback: return list
    if (isTRUE(show)) {
      print(p_alpha); print(p_K); print(p_w1); print(p_tbl)
    }
    return(invisible(list(alpha = p_alpha, K = p_K, w1 = p_w1, summary = p_tbl)))
  }

  if (isTRUE(show)) {
    grid::grid.newpage()
    grid::grid.draw(g)
  }

  g
}


# =============================================================================
# Dual-Anchor Visualization Functions
# =============================================================================

#' Plot Dual-Anchor Comparison Dashboard
#'
#' Creates a comparison dashboard showing K-only vs dual-anchor solutions.
#' Displays changes in alpha, K, and the first size-biased weight
#' \eqn{W_{SB}} side-by-side.
#'
#' @param fit_dual A canonical hard or soft dual fit from
#'   \code{DPprior_dual_hard()} or \code{DPprior_dual_soft()}. Retained legacy
#'   fits have no authoritative comparison lineage and raise a typed
#'   unavailable condition.
#' @param fit_K_only Optional canonical K-only fit. If supplied, it must exactly
#'   match the authoritative \code{provenance.input_fit} reference retained by
#'   \code{fit_dual}.
#' @param engine "ggplot2" (default) or "base".
#' @param base_size Base font size.
#' @param title Optional title for the dashboard.
#' @param show If TRUE, draw the plot.
#'
#' @return A gtable grob or list of ggplot objects.
#'
#' @examples
#' fit_K <- DPprior_fit(J = 50, mu_K = 5, var_K = 8)
#' fit_dual <- DPprior_dual_soft(
#'   fit_K,
#'   target = list(
#'     metric = "wsb_tail", relation = "target",
#'     threshold = 0.5, value = 0.3
#'   ),
#'   lambda = 0.5
#' )
#' plot_dual_comparison(fit_dual)
#'
#' @seealso \code{\link{DPprior_dual_hard}},
#'   \code{\link{DPprior_dual_soft}}, and \code{\link{plot.DPprior_fit}}
#'
#' @family visualization
#'
#' @export
plot_dual_comparison <- function(fit_dual,
                                 fit_K_only = NULL,
                                 engine = c("ggplot2", "base"),
                                 base_size = 10,
                                 title = NULL,
                                 show = TRUE) {
  engine <- match.arg(engine)
  dual <- .dpprior_visualization_dual_view(
    fit_dual, fit_K_only,
    operation = "dual comparison visualization"
  )
  baseline <- dual[["baseline", exact = TRUE]]
  current <- dual[["current", exact = TRUE]]
  baseline_parameters <- baseline[["parameters", exact = TRUE]]
  current_parameters <- current[["parameters", exact = TRUE]]
  a_K <- baseline_parameters[["a", exact = TRUE]]
  b_K <- baseline_parameters[["b", exact = TRUE]]
  a_dual <- current_parameters[["a", exact = TRUE]]
  b_dual <- current_parameters[["b", exact = TRUE]]
  J <- current[["core", exact = TRUE]][["J", exact = TRUE]]
  mode_title <- .dpprior_visualization_dual_title(dual)

  if (engine == "base" || !.dpprior_has_ggplot2()) {
    if (!isTRUE(show)) {
      return(invisible(NULL))
    }
    .dpprior_base_dual_comparison(
      dual, if (is.null(title)) mode_title else title
    )
    return(invisible(NULL))
  }

  .dpprior_require_ggplot2()
  colors <- DPprior_colors()

  # ---- Panel 1: Alpha density comparison ----
  x_max <- max(
    stats::qgamma(0.999, shape = a_K, rate = b_K),
    stats::qgamma(0.999, shape = a_dual, rate = b_dual)
  ) * 1.1
  x <- seq(0.001, x_max, length.out = 400)

  df_alpha <- data.frame(
    x = rep(x, 2),
    density = c(
      stats::dgamma(x, shape = a_K, rate = b_K),
      stats::dgamma(x, shape = a_dual, rate = b_dual)
    ),
    Method = rep(c("K-only", "Dual-anchor"), each = length(x))
  )

  p_alpha <- ggplot2::ggplot(df_alpha, ggplot2::aes(x = x, y = density,
                                                    color = Method, linetype = Method)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_color_manual(values = c("K-only" = colors$k_only, "Dual-anchor" = colors$dual)) +
    ggplot2::scale_linetype_manual(values = c("K-only" = "solid", "Dual-anchor" = "solid")) +
    theme_DPprior(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
      x = expression(alpha), y = "Density",
      title = "Alpha Prior Comparison",
      subtitle = sprintf("K-only: Gamma(%.3f, %.3f)  |  Dual: Gamma(%.3f, %.3f)",
                         a_K, b_K, a_dual, b_dual)
    )

  # ---- Panel 2: K PMF comparison ----
  summ_K <- .dpprior_compute_summary(a = a_K, b = b_K, J = J)
  summ_dual <- .dpprior_compute_summary(a = a_dual, b = b_dual, J = J)

  max_k <- max(which(summ_K$K$cdf >= 0.999)[1], which(summ_dual$K$cdf >= 0.999)[1], na.rm = TRUE)
  max_k <- min(max_k, J, 50)

  df_K <- data.frame(
    k = rep(summ_K$K$k[1:max_k], 2),
    pmf = c(summ_K$K$pmf[1:max_k], summ_dual$K$pmf[1:max_k]),
    Method = rep(c("K-only", "Dual-anchor"), each = max_k)
  )

  p_K <- ggplot2::ggplot(df_K, ggplot2::aes(x = k, y = pmf, fill = Method)) +
    ggplot2::geom_col(position = "dodge", alpha = 0.8, width = 0.8) +
    ggplot2::scale_fill_manual(values = c("K-only" = colors$k_only, "Dual-anchor" = colors$dual)) +
    theme_DPprior(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
      x = expression(K[J]), y = "PMF",
      title = "K Distribution Comparison",
      subtitle = sprintf(
        "Canonical achieved E[K]: K-only=%.2f | dual=%.2f",
        baseline[["achieved", exact = TRUE]][["K", exact = TRUE]][[
          "mean", exact = TRUE
        ]],
        current[["achieved", exact = TRUE]][["K", exact = TRUE]][[
          "mean", exact = TRUE
        ]]
      ),
      caption = sprintf(
        "K-only PMF: %s | Dual PMF: %s",
        .dpprior_K_provenance_text(summ_K),
        .dpprior_K_provenance_text(summ_dual)
      )
    )

  # ---- Panel 3: first size-biased weight density comparison ----
  x_w <- seq(1e-6, 1 - 1e-6, length.out = 400)
  dens_K <- .dpprior_density_w1(x_w, a_K, b_K)
  dens_dual <- .dpprior_density_w1(x_w, a_dual, b_dual)
  dens_K[!is.finite(dens_K)] <- NA
  dens_dual[!is.finite(dens_dual)] <- NA

  df_w1 <- data.frame(
    x = rep(x_w, 2),
    density = c(dens_K, dens_dual),
    Method = rep(c("K-only", "Dual-anchor"), each = length(x_w))
  )

  y_max_w <- min(max(c(dens_K, dens_dual), na.rm = TRUE), 20)

  p_gt_50_K <- .dpprior_prob_w1_exceeds(0.5, a_K, b_K)
  p_gt_50_dual <- .dpprior_prob_w1_exceeds(0.5, a_dual, b_dual)

  p_w1 <- ggplot2::ggplot(df_w1, ggplot2::aes(x = x, y = density,
                                              color = Method, linetype = Method)) +
    ggplot2::geom_line(linewidth = 1, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed", color = colors$warning, linewidth = 0.5) +
    ggplot2::scale_color_manual(values = c("K-only" = colors$k_only, "Dual-anchor" = colors$dual)) +
    ggplot2::scale_linetype_manual(values = c("K-only" = "solid", "Dual-anchor" = "solid")) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, y_max_w * 1.1)) +
    theme_DPprior(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
      x = expression(W[SB]), y = "Density",
      title = "First Size-Biased Weight Distribution Comparison",
      subtitle = sprintf(
        "P(W_SB>0.5): K-only=%.1f%% | dual=%.1f%%",
        100 * p_gt_50_K, 100 * p_gt_50_dual
      ),
      caption = .dpprior_visualization_weight_target_text(
        dual[["target_weight", exact = TRUE]]
      )
    )

  # ---- Panel 4: Summary comparison table ----
  p_tbl <- .dpprior_dual_comparison_table(
    dual, summ_K, summ_dual, base_size
  )

  # Combine into dashboard
  if (is.null(title)) title <- mode_title

  g <- .dpprior_dashboard_gtable(p_alpha, p_K, p_w1, p_tbl, title = title)

  if (is.null(g)) {
    if (isTRUE(show)) {
      print(p_alpha); print(p_K); print(p_w1); print(p_tbl)
    }
    return(invisible(list(alpha = p_alpha, K = p_K, w1 = p_w1, summary = p_tbl)))
  }

  if (isTRUE(show)) {
    grid::grid.newpage()
    grid::grid.draw(g)
  }

  g
}


#' Dual comparison summary table
#' @keywords internal
.dpprior_dual_comparison_table <- function(dual, summ_K, summ_dual,
                                           base_size = 11) {
  .dpprior_require_ggplot2()
  colors <- DPprior_colors()
  baseline <- dual[["baseline", exact = TRUE]]
  current <- dual[["current", exact = TRUE]]
  baseline_parameters <- baseline[["parameters", exact = TRUE]]
  current_parameters <- current[["parameters", exact = TRUE]]
  baseline_K <- baseline[["achieved", exact = TRUE]][["K", exact = TRUE]]
  current_K <- current[["achieved", exact = TRUE]][["K", exact = TRUE]]
  metrics <- c(
    "Gamma(a, b)", "E[alpha]", "Canonical achieved E[K]",
    "Canonical achieved Var(K)", "K PMF provenance", "E[W_SB]",
    "P(W_SB > 0.5)", "P(W_SB > 0.9)", "Status", "Verified"
  )
  K_only <- c(
    sprintf(
      "(%.3f, %.3f)", baseline_parameters[["a", exact = TRUE]],
      baseline_parameters[["b", exact = TRUE]]
    ),
    sprintf("%.3f", summ_K$alpha$mean),
    sprintf("%.2f", baseline_K[["mean", exact = TRUE]]),
    sprintf("%.2f", baseline_K[["variance", exact = TRUE]]),
    .dpprior_K_provenance_text(summ_K),
    sprintf("%.3f", summ_K$w1$mean),
    sprintf("%.1f%%", 100 * summ_K$w1$p_gt_50),
    sprintf("%.1f%%", 100 * summ_K$w1$p_gt_90),
    baseline[["core", exact = TRUE]][["status", exact = TRUE]],
    if (isTRUE(baseline[["core", exact = TRUE]][[
      "verified", exact = TRUE
    ]])) "yes" else "no"
  )
  Dual <- c(
    sprintf(
      "(%.3f, %.3f)", current_parameters[["a", exact = TRUE]],
      current_parameters[["b", exact = TRUE]]
    ),
    sprintf("%.3f", summ_dual$alpha$mean),
    sprintf("%.2f", current_K[["mean", exact = TRUE]]),
    sprintf("%.2f", current_K[["variance", exact = TRUE]]),
    .dpprior_K_provenance_text(summ_dual),
    sprintf("%.3f", summ_dual$w1$mean),
    sprintf("%.1f%%", 100 * summ_dual$w1$p_gt_50),
    sprintf("%.1f%%", 100 * summ_dual$w1$p_gt_90),
    current[["core", exact = TRUE]][["status", exact = TRUE]],
    if (isTRUE(current[["core", exact = TRUE]][[
      "verified", exact = TRUE
    ]])) "yes" else "no"
  )
  target_text <- .dpprior_visualization_weight_target_text(
    dual[["target_weight", exact = TRUE]]
  )
  achieved_weight <- dual[["achieved_weight", exact = TRUE]][[
    "value", exact = TRUE
  ]]
  if (identical(dual[["mode", exact = TRUE]], "dual_hard")) {
    constraint <- dual[["constraint", exact = TRUE]]
    feasibility <- constraint[["feasibility", exact = TRUE]]
    metrics <- c(
      metrics, "Constraint", "Achieved weight metric", "Residual",
      "Satisfied", "Feasibility"
    )
    K_only <- c(K_only, rep("", 5L))
    Dual <- c(
      Dual, target_text, sprintf("%.4g", achieved_weight),
      sprintf("%.4g", constraint[["residual", exact = TRUE]]),
      if (isTRUE(constraint[["satisfied", exact = TRUE]])) "yes" else "no",
      feasibility[["classification", exact = TRUE]]
    )
  } else {
    tradeoff <- dual[["tradeoff", exact = TRUE]]
    metrics <- c(
      metrics, "Weight target", "Achieved weight metric", "Lambda",
      "K loss", "Weight loss", "Total loss", "Endpoint"
    )
    K_only <- c(K_only, rep("", 7L))
    Dual <- c(
      Dual, target_text, sprintf("%.4g", achieved_weight),
      sprintf("%.4g", tradeoff[["lambda", exact = TRUE]]),
      sprintf("%.4g", tradeoff[["K_loss", exact = TRUE]]),
      sprintf("%.4g", tradeoff[["weight_loss", exact = TRUE]]),
      sprintf("%.4g", tradeoff[["total_loss", exact = TRUE]]),
      if (isTRUE(tradeoff[["endpoint", exact = TRUE]])) "yes" else "no"
    )
  }
  tbl <- data.frame(
    Metric = metrics, K_only = K_only, Dual = Dual,
    stringsAsFactors = FALSE
  )

  tbl$row <- seq_len(nrow(tbl))

  p <- ggplot2::ggplot(tbl, ggplot2::aes(y = -row)) +
    ggplot2::geom_text(ggplot2::aes(x = 0, label = Metric),
                       hjust = 0, fontface = "bold", color = colors$ink, size = 2.8) +
    ggplot2::geom_text(ggplot2::aes(x = 0.45, label = K_only),
                       hjust = 0.5, color = colors$k_only, size = 2.8, fontface = "bold") +
    ggplot2::geom_text(ggplot2::aes(x = 0.85, label = Dual),
                       hjust = 0.5, color = colors$dual, size = 2.8, fontface = "bold") +
    # Column headers
    ggplot2::annotate("text", x = 0.45, y = 0.5, label = "K-only",
                      fontface = "bold", color = colors$k_only, size = 3) +
    ggplot2::annotate("text", x = 0.85, y = 0.5, label = "Dual",
                      fontface = "bold", color = colors$dual, size = 3) +
    ggplot2::scale_x_continuous(limits = c(-0.05, 1.05)) +
    ggplot2::scale_y_continuous(limits = c(-nrow(tbl) - 1, 1)) +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = base_size),
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    ) +
    ggplot2::labs(title = "Comparison Summary")

  p
}


#' Plot Trade-off Curve
#'
#' Visualizes the Pareto trade-off between K_J fit and weight constraint
#' across different lambda values.
#'
#' @param tradeoff_data Data frame from compute_tradeoff_curve().
#' @param metric Which metric to plot on y-axis: \code{"w1_prob_gt_50"}
#'   (default), \code{"E_w1"}, \code{"K_loss"}, \code{"mu_K"}, or
#'   \code{"var_K"}. The two legacy column names refer specifically to
#'   \eqn{W_{SB}}, the first size-biased DP weight.
#' @param target_value Optional target value to mark with horizontal line.
#' @param engine "ggplot2" (default) or "base".
#' @param base_size Base font size.
#' @param title Optional title.
#' @param show If TRUE, print the plot.
#'
#' @return A ggplot object or invisible(NULL).
#'
#' @examples
#' curve <- compute_tradeoff_curve(
#'   J = 50,
#'   K_target = list(mu_K = 5, var_K = 8),
#'   w1_target = list(prob = list(threshold = 0.5, value = 0.25)),
#'   lambda_seq = seq(0.1, 1, by = 0.1)
#' )
#' plot_tradeoff_curve(curve, target_value = 0.25)
#'
#' @seealso \code{\link{DPprior_fit}} for fitting, \code{\link{plot.DPprior_fit}} for S3 plot method
#'
#' @family visualization
#'
#' @export
plot_tradeoff_curve <- function(tradeoff_data,
                                metric = c("w1_prob_gt_50", "E_w1", "K_loss",
                                           "mu_K", "var_K"),
                                target_value = NULL,
                                engine = c("ggplot2", "base"),
                                base_size = 11,
                                title = NULL,
                                show = TRUE) {
  engine <- match.arg(engine)
  metric <- match.arg(metric)
  curve <- .dpprior_visualization_curve_view(tradeoff_data, metric)
  plot_data <- curve[["data", exact = TRUE]]
  line_data <- curve[["line_data", exact = TRUE]]

  # Labels
  metric_labels <- list(
    w1_prob_gt_50 = "P(W_SB > 0.5)",
    E_w1 = "E[W_SB]",
    K_loss = "Scaled K loss (fixed input-derived scale)",
    mu_K = "E[K]",
    var_K = "Var(K)"
  )
  y_label <- metric_labels[[metric]]

  if (engine == "base" || !.dpprior_has_ggplot2()) {
    if (!isTRUE(show)) {
      return(invisible(NULL))
    }
    graphics::plot(plot_data[["lambda", exact = TRUE]],
                   plot_data[["value", exact = TRUE]],
                   type = "p", pch = 19, col = "steelblue4",
                   xlab = "lambda (weight on K anchor)",
                   ylab = y_label,
                   main = if (is.null(title)) "Trade-off Curve" else title)
    for (segment in unique(line_data[["segment", exact = TRUE]])) {
      selected <- line_data[["segment", exact = TRUE]] == segment
      if (sum(selected) > 1L) {
        graphics::lines(
          line_data[["lambda", exact = TRUE]][selected],
          line_data[["value", exact = TRUE]][selected],
          lwd = 2, col = "steelblue4"
        )
      }
    }
    if (!is.null(target_value)) {
      graphics::abline(h = target_value, lty = 2, col = "firebrick3")
    }
    graphics::grid()
    return(invisible(NULL))
  }

  .dpprior_require_ggplot2()
  colors <- DPprior_colors()

  p <- ggplot2::ggplot(
    plot_data, ggplot2::aes(x = lambda, y = value)
  ) +
    ggplot2::geom_line(
      data = line_data, ggplot2::aes(group = segment),
      color = colors$primary, linewidth = 1
    ) +
    ggplot2::geom_point(color = colors$primary, size = 3) +
    theme_DPprior(base_size = base_size) +
    ggplot2::labs(
      x = "lambda (weight on K anchor)",
      y = y_label,
      title = if (is.null(title)) "Dual-Anchor Trade-off Curve" else title,
      subtitle = sprintf(
        paste0(
          "0 < lambda < 1: soft trade-off | lambda=1: K-only | ",
          "%d/%d metric values available; gaps are not imputed"
        ),
        curve[["available_points", exact = TRUE]],
        curve[["requested_points", exact = TRUE]]
      )
    )

  if (!is.null(target_value)) {
    p <- p + ggplot2::geom_hline(yintercept = target_value,
                                 linetype = "dashed", color = colors$warning, linewidth = 0.8) +
      ggplot2::annotate("text", x = 0.02, y = target_value,
                        label = sprintf("Target: %.2f", target_value),
                        hjust = 0, vjust = -0.5, color = colors$warning, size = 3)
  }

  # Mark the lowest requested trade-off, midpoint, and K-only endpoint.
  key_lambdas <- unique(c(
    min(plot_data[["lambda", exact = TRUE]], na.rm = TRUE), 0.5, 1
  ))
  key_data <- plot_data[
    plot_data[["available", exact = TRUE]] &
      plot_data[["lambda", exact = TRUE]] %in% key_lambdas,
    , drop = FALSE
  ]
  if (nrow(key_data) > 0) {
    p <- p + ggplot2::geom_point(data = key_data, color = colors$accent, size = 4, shape = 18)
  }

  if (isTRUE(show)) print(p)
  p
}


#' Plot Trade-off Multi-Panel
#'
#' Creates a multi-panel view of the trade-off curve showing multiple metrics.
#'
#' @param tradeoff_data Data frame from compute_tradeoff_curve().
#' @param w1_target_prob Optional target probability for
#'   \eqn{P(W_{SB} > 0.5)}.
#' @param engine "ggplot2" (default) or "base".
#' @param base_size Base font size.
#' @param title Optional title.
#' @param show If TRUE, draw the plot.
#'
#' @return A gtable grob or list of ggplot objects.
#'
#' @seealso \code{\link{DPprior_fit}} for fitting, \code{\link{plot.DPprior_fit}} for S3 plot method
#'
#' @family visualization
#'
#' @export
plot_tradeoff_dashboard <- function(tradeoff_data,
                                    w1_target_prob = NULL,
                                    engine = c("ggplot2", "base"),
                                    base_size = 10,
                                    title = NULL,
                                    show = TRUE) {
  engine <- match.arg(engine)
  # Validate the complete retained path before either engine can return early.
  .dpprior_visualization_curve_view(tradeoff_data, "K_loss")

  if (engine == "base" || !.dpprior_has_ggplot2()) {
    if (!isTRUE(show)) {
      return(invisible(NULL))
    }
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)
    graphics::par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

    plot_tradeoff_curve(tradeoff_data, "w1_prob_gt_50", w1_target_prob,
                        "base", show = show)
    plot_tradeoff_curve(tradeoff_data, "E_w1", NULL, "base", show = show)
    plot_tradeoff_curve(tradeoff_data, "mu_K", NULL, "base", show = show)
    plot_tradeoff_curve(tradeoff_data, "var_K", NULL, "base", show = show)

    return(invisible(NULL))
  }

  .dpprior_require_ggplot2()

  # Panel 1: P(W_SB > 0.5) vs lambda
  p1 <- plot_tradeoff_curve(tradeoff_data, "w1_prob_gt_50", w1_target_prob,
                            "ggplot2", base_size, "P(W_SB > 0.5) vs lambda", FALSE)

  # Panel 2: E[W_SB] vs lambda
  p2 <- plot_tradeoff_curve(tradeoff_data, "E_w1", NULL,
                            "ggplot2", base_size, "E[W_SB] vs lambda", FALSE)

  # Panels 3 and 4 use the same canonical status-aware path validation and do
  # not connect across unavailable retained points.
  p3 <- plot_tradeoff_curve(
    tradeoff_data, "mu_K", NULL, "ggplot2", base_size,
    "E[K] vs lambda", FALSE
  )
  p4 <- plot_tradeoff_curve(
    tradeoff_data, "var_K", NULL, "ggplot2", base_size,
    "Var(K) vs lambda", FALSE
  )

  # Combine
  if (is.null(title)) {
    title <- "Dual-Anchor Trade-off Analysis"
  }

  g <- .dpprior_dashboard_gtable(p1, p2, p3, p4, title = title)

  if (is.null(g)) {
    if (isTRUE(show)) {
      print(p1); print(p2); print(p3); print(p4)
    }
    return(invisible(list(p_w1 = p1, p_Ew1 = p2, p_EK = p3, p_VarK = p4)))
  }

  if (isTRUE(show)) {
    grid::grid.newpage()
    grid::grid.draw(g)
  }

  g
}


#' Plot Canonical Dual Comparison Dashboard
#'
#' Creates the canonical four-panel comparison view: alpha, \eqn{K_J}, and
#' \eqn{W_{SB}} distributions plus a mode-specific evidence table. Hard fits
#' show constraint evidence; soft fits show the recorded fixed-scale trade-off.
#' Retained legacy fits raise a typed comparison-lineage unavailable condition.
#'
#' @param fit_dual A canonical hard or soft dual fit from
#'   \code{DPprior_dual_hard()} or \code{DPprior_dual_soft()}.
#' @param tradeoff_data Optional canonical \code{dpprior_tradeoff_curve}. When
#'   supplied, its retained evidence is validated before the fit comparison is
#'   drawn; it is not substituted for the fit's authoritative trade-off record.
#' @param engine "ggplot2" (default) or "base".
#' @param base_size Base font size.
#' @param title Optional title.
#' @param show If TRUE, draw the plot.
#'
#' @return A gtable grob, a list of plots when gtable assembly is unavailable,
#'   or invisible(NULL) for the base engine.
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso \code{\link{DPprior_fit}} for fitting, \code{\link{plot.DPprior_fit}} for S3 plot method
#'
#' @family visualization
#'
#' @export
plot_dual_dashboard <- function(fit_dual,
                                tradeoff_data = NULL,
                                engine = c("ggplot2", "base"),
                                base_size = 10,
                                title = NULL,
                                show = TRUE) {
  engine <- match.arg(engine)
  # The comparison dashboard is the sole canonical dual view. Its final panel
  # is mode-specific: a hard fit exposes constraint evidence, while a soft fit
  # exposes the recorded trade-off. A legacy fit has no authoritative input-fit
  # lineage and therefore fails with a typed unavailable condition.
  if (!is.null(tradeoff_data)) {
    .dpprior_visualization_curve_view(tradeoff_data, "K_loss")
  }
  plot_dual_comparison(
    fit_dual = fit_dual, engine = engine, base_size = base_size,
    title = title, show = show
  )
}


#' Base R dual comparison fallback
#' @keywords internal
.dpprior_base_dual_comparison <- function(dual, title) {
  baseline <- dual[["baseline", exact = TRUE]]
  current <- dual[["current", exact = TRUE]]
  baseline_parameters <- baseline[["parameters", exact = TRUE]]
  current_parameters <- current[["parameters", exact = TRUE]]
  a_K <- baseline_parameters[["a", exact = TRUE]]
  b_K <- baseline_parameters[["b", exact = TRUE]]
  a_dual <- current_parameters[["a", exact = TRUE]]
  b_dual <- current_parameters[["b", exact = TRUE]]
  J <- current[["core", exact = TRUE]][["J", exact = TRUE]]
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)

  if (!is.null(title) && nzchar(title)) {
    graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  } else {
    graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  }

  # Alpha comparison
  x_max <- max(stats::qgamma(0.999, a_K, b_K), stats::qgamma(0.999, a_dual, b_dual)) * 1.1
  x <- seq(0.001, x_max, length.out = 200)
  graphics::plot(x, stats::dgamma(x, a_K, b_K), type = "l", col = "#4682B4", lwd = 2,
                 xlab = expression(alpha), ylab = "Density", main = "(A) Alpha Prior")
  graphics::lines(x, stats::dgamma(x, a_dual, b_dual), col = "#E67E22", lwd = 2)
  graphics::legend("topright", c("K-only", "Dual"), lty = c(1, 1),
                   col = c("#4682B4", "#E67E22"), lwd = 2, cex = 0.8)

  # First size-biased weight comparison
  x_w <- seq(1e-4, 1 - 1e-4, length.out = 200)
  dens_K <- .dpprior_density_w1(x_w, a_K, b_K)
  dens_dual <- .dpprior_density_w1(x_w, a_dual, b_dual)
  y_max <- min(max(c(dens_K, dens_dual), na.rm = TRUE), 15)

  graphics::plot(x_w, dens_K, type = "l", col = "#4682B4", lwd = 2,
                 xlab = expression(W[SB]), ylab = "Density",
                 main = "(B) First Size-Biased Weight Distribution",
                 ylim = c(0, y_max))
  graphics::lines(x_w, dens_dual, col = "#E67E22", lwd = 2)
  graphics::abline(v = 0.5, lty = 2, col = "firebrick3")

  # Summary text
  summ_K <- .dpprior_compute_summary(a = a_K, b = b_K, J = J)
  summ_dual <- .dpprior_compute_summary(a = a_dual, b = b_dual, J = J)

  graphics::plot.new()
  graphics::title(main = "(C) Comparison Summary")
  mode_lines <- if (identical(dual[["mode", exact = TRUE]], "dual_hard")) {
    constraint <- dual[["constraint", exact = TRUE]]
    c(
      "Mode: hard constraint",
      .dpprior_visualization_weight_target_text(
        dual[["target_weight", exact = TRUE]]
      ),
      sprintf(
        "Residual=%.4g; satisfied=%s",
        constraint[["residual", exact = TRUE]],
        if (isTRUE(constraint[["satisfied", exact = TRUE]])) "yes" else "no"
      )
    )
  } else {
    tradeoff <- dual[["tradeoff", exact = TRUE]]
    c(
      "Mode: soft trade-off",
      sprintf("Lambda=%.3g", tradeoff[["lambda", exact = TRUE]]),
      sprintf(
        "K loss=%.4g; weight loss=%.4g; total=%.4g",
        tradeoff[["K_loss", exact = TRUE]],
        tradeoff[["weight_loss", exact = TRUE]],
        tradeoff[["total_loss", exact = TRUE]]
      )
    )
  }
  lines <- c(
    mode_lines,
    "",
    "K-only:",
    sprintf("  Gamma(%.3f, %.3f)", a_K, b_K),
    sprintf("  K PMF %s", .dpprior_K_provenance_text(summ_K)),
    sprintf("  P(W_SB>0.5) = %.1f%%", 100 * summ_K$w1$p_gt_50),
    "",
    "Dual-anchor:",
    sprintf("  Gamma(%.3f, %.3f)", a_dual, b_dual),
    sprintf("  K PMF %s", .dpprior_K_provenance_text(summ_dual)),
    sprintf("  P(W_SB>0.5) = %.1f%%", 100 * summ_dual$w1$p_gt_50)
  )
  y_pos <- seq(0.9, 0.1, length.out = length(lines))
  for (i in seq_along(lines)) {
    graphics::text(0.1, y_pos[i], lines[i], adj = 0, cex = 0.9)
  }

  # Reduction summary
  graphics::plot.new()
  graphics::title(main = "(D) Reduction")
  red_p50 <- (summ_K$w1$p_gt_50 - summ_dual$w1$p_gt_50) / summ_K$w1$p_gt_50
  red_mean <- (summ_K$w1$mean - summ_dual$w1$mean) / summ_K$w1$mean

  lines2 <- c(
    sprintf("P(W_SB>0.5) reduction: %.1f%%", 100 * red_p50),
    sprintf("E[W_SB] reduction: %.1f%%", 100 * red_mean),
    "",
    sprintf("E[K] change: %.2f -> %.2f", summ_K$K$mean, summ_dual$K$mean),
    sprintf("Var(K) change: %.2f -> %.2f", summ_K$K$var, summ_dual$K$var)
  )
  y_pos2 <- seq(0.8, 0.2, length.out = length(lines2))
  for (i in seq_along(lines2)) {
    graphics::text(0.1, y_pos2[i], lines2[i], adj = 0, cex = 0.9)
  }

  if (!is.null(title) && nzchar(title)) {
    graphics::mtext(title, outer = TRUE, cex = 1.2, font = 2)
  }
}


# =============================================================================
# Base R Fallbacks
# =============================================================================

.dpprior_base_plot_alpha <- function(df, alpha_mean, ci, alpha_cv, a, b,
                                     context = NULL) {
  graphics::plot(df$x, df$density, type = "l", lwd = 2,
                 xlab = expression(alpha), ylab = "Density",
                 main = expression(paste("(A) Prior on ", alpha)))
  graphics::polygon(df$x, df$density, col = "grey90", border = NA)
  graphics::lines(df$x, df$density, lwd = 2, col = "steelblue4")
  graphics::abline(v = alpha_mean, col = "darkred", lwd = 2, lty = 2)
  graphics::abline(v = ci, col = "gray50", lwd = 1, lty = 3)
  graphics::mtext(sprintf("Mean=%.3f, CV=%.2f, CI=[%.2f, %.2f], Gamma(%.4f, %.3f)",
                          alpha_mean, alpha_cv, ci[1], ci[2], a, b),
                  side = 3, line = 0.2, adj = 0, cex = 0.8)
  if (!is.null(context)) {
    graphics::mtext(context, side = 1, line = 3.1, adj = 0, cex = 0.65)
  }
}

.dpprior_base_plot_K <- function(
    df, target_mu, achieved_mu, K_mean, K_var, K_mode, a, b,
    status, verified, method, context = NULL) {
  graphics::barplot(height = df$pmf, names.arg = df$k, border = NA,
                    col = "steelblue3", xlab = "k", ylab = "PMF",
                    main = expression(paste("(B) Prior PMF of ", K[J])))
  mids <- seq_along(df$k)
  cdf_scaled <- df$cdf * max(df$pmf)
  graphics::lines(mids, cdf_scaled, lwd = 2, col = "grey30")
  graphics::mtext(sprintf("E[K]=%.2f, Var=%.2f, Mode=%d, Gamma(%.4f, %.3f)",
                          K_mean, K_var, K_mode, a, b),
                  side = 3, line = 0.2, adj = 0, cex = 0.8)
  graphics::mtext(
    paste(c(
      if (is.null(context)) character() else context,
      sprintf("PMF status: %s; verified: %s; method: %s",
              status, if (isTRUE(verified)) "yes" else "no", method)
    ), collapse = "; "),
    side = 1, line = 3.1, adj = 0, cex = 0.65
  )
}

.dpprior_base_plot_w1 <- function(df, thresholds, mean_w, median_w,
                                  p_gt_first, p_gt_second, a, b,
                                  context = NULL) {
  graphics::plot(df$x, df$density, type = "l", lwd = 2, col = "steelblue4",
                 xlab = expression(W[SB]), ylab = "Density",
                 main = expression(paste("(C) Prior Density of ", W[SB])))
  shade <- df[df$x >= thresholds[1], , drop = FALSE]
  graphics::polygon(c(shade$x, rev(shade$x)),
                    c(shade$density, rep(0, nrow(shade))),
                    col = "grey85", border = NA)
  graphics::lines(df$x, df$density, lwd = 2, col = "steelblue4")
  graphics::abline(v = thresholds, col = "firebrick3", lwd = 1.5, lty = 2)
  graphics::mtext(sprintf(
    paste0(
      "E[W_SB]=%.3f, P(W_SB>%.3g)=%.1f%%, ",
      "P(W_SB>%.3g)=%.1f%%, Gamma(%.4f, %.3f)"
    ),
    mean_w, thresholds[[1L]], 100 * p_gt_first,
    thresholds[[2L]], 100 * p_gt_second, a, b
  ),
                  side = 3, line = 0.2, adj = 0, cex = 0.8)
  if (!is.null(context)) {
    graphics::mtext(context, side = 1, line = 3.1, adj = 0, cex = 0.65)
  }
}

.dpprior_base_summary_panel <- function(fit, ci_level = 0.95) {
  summ <- .dpprior_compute_summary(fit, ci_level = ci_level)

  graphics::plot.new()
  graphics::title(main = "(D) Summary Statistics")

  lines <- c(
    sprintf("J = %d", summ$J),
    sprintf("Mode = %s", summ$mode),
    sprintf(
      "Status = %s; verified = %s", summ$status,
      if (isTRUE(summ$verified)) "yes" else "no"
    ),
    sprintf("Gamma(a=%.4f, b=%.3f)", summ$a, summ$b),
    sprintf("E[alpha] = %.3f, CV = %.2f", summ$alpha$mean, summ$alpha$cv),
    sprintf("E[K] = %.2f, Var(K) = %.2f", summ$achieved$mu_K, summ$achieved$var_K),
    sprintf("Achieved estimand = %s", summ$achieved$estimand),
    sprintf("E[W_SB] = %.3f", summ$w1$mean),
    sprintf("P(W_SB>0.5) = %.1f%%", 100 * summ$w1$p_gt_50),
    if (identical(summ$mode, "dual_legacy")) {
      "DEPRECATED legacy dual fit; descriptive single-fit view only"
    } else {
      character()
    }
  )

  y_pos <- seq(0.85, 0.15, length.out = length(lines))
  for (i in seq_along(lines)) {
    graphics::text(0.1, y_pos[i], lines[i], adj = 0, cex = 0.9)
  }
}


# =============================================================================
# Module Verification
# =============================================================================

#' Verify Visualization Module
#'
#' @param verbose Logical; if TRUE, print progress messages. Default is TRUE.
#'
#' @keywords internal
verify_visualization <- function(verbose = TRUE) {
  if (isTRUE(verbose)) cat("Verifying canonical visualization consumers...\n")
  tryCatch({
    fit <- DPprior_fit(
      20L, 4, 8, method = "A2-MN", M = 80L,
      check_diagnostics = FALSE
    )
    soft <- DPprior_dual_soft(
      fit,
      list(
        metric = "wsb_tail", relation = "target",
        threshold = 0.5, value = 0.3
      ),
      lambda = 1, M_fit = 80L, M_verify = 160L
    )
    stopifnot(
      is.null(plot_alpha_prior(fit, engine = "base", show = FALSE)),
      is.null(plot_K_prior(fit, engine = "base", show = FALSE)),
      is.null(plot_w1_prior(fit, engine = "base", show = FALSE)),
      is.null(plot_prior_dashboard(fit, engine = "base", show = FALSE)),
      is.null(plot_dual_comparison(soft, engine = "base", show = FALSE))
    )
    if (.dpprior_has_ggplot2()) {
      stopifnot(
        inherits(plot_alpha_prior(fit, show = FALSE), "ggplot"),
        inherits(plot_K_prior(fit, show = FALSE), "ggplot"),
        inherits(plot_w1_prior(fit, show = FALSE), "ggplot")
      )
    }
    if (isTRUE(verbose)) cat("Canonical visualization verification passed\n")
    TRUE
  }, error = function(condition) {
    if (isTRUE(verbose)) cat("FAILED:", conditionMessage(condition), "\n")
    FALSE
  })
}
