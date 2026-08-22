# =============================================================================
# Module 19: bounded-discrete elicitation targets
# =============================================================================
# Internal target machinery plus the Phase 9 DPprior_target_K() entry point.

.DP_TARGET_STATUS_CODES <- c(
  "converged", "boundary", "approximate", "infeasible", "failed"
)

.DP_INTERVAL_TYPES <- c("hard_bounds", "central_mass", "equal_tail")
.DP_INTERVAL_FAMILIES <- c("maxent", "discrete_normal", "beta_binomial")
.DP_ELICITATION_TOL_MIN <- 64 * .Machine$double.eps

.dp_interval_abort <- function(message, classes, argument = "K_interval",
                               value = NULL, expected = NULL,
                               code = "invalid_interval", ...) {
  stop(.dpprior_new_condition(
    message = message,
    classes = c(classes, "dpprior_invalid_input", "dpprior_error", "error"),
    argument = argument,
    value = value,
    expected = expected,
    code = code,
    ...
  ))
}

.dp_interval_ambiguous <- function(value, J, code = "bare_range") {
  .dp_interval_abort(
    paste0(
      "K_interval is ambiguous: it could mean hard support, central ",
      "coverage, equal tails, or an informal plausible range. Supply ",
      "inclusive integer endpoints plus type and an explicit family; ",
      "probabilistic intervals also require explicit coverage. No variance ",
      "was inferred. For example: list(lower=3L, upper=10L, ",
      "type='equal_tail', coverage=.80, family='maxent')."
    ),
    c("dpprior_ambiguous_interval_error", "dpprior_interval_error"),
    value = value,
    expected = "named interval list with explicit type and family",
    code = code,
    received = value,
    missing = c("type", "family", "coverage for probabilistic intervals"),
    allowed_types = .DP_INTERVAL_TYPES,
    support = c(1L, as.integer(J)),
    inferred_variance = NULL
  )
}

.dp_interval_scalar_character <- function(x, name, choices = NULL) {
  if (!is.character(x) || !is.null(dim(x)) || is.object(x) ||
      length(x) != 1L || is.na(x) || !nzchar(x)) {
    .dp_interval_abort(
      sprintf("K_interval$%s must be one non-missing character value.", name),
      c("dpprior_interval_structure_error", "dpprior_type_error"),
      paste0("K_interval$", name), x, "character scalar", "type"
    )
  }
  if (!is.null(choices) && !(x %in% choices)) {
    .dp_interval_abort(
      sprintf(
        "K_interval$%s must be one of: %s.",
        name, paste(choices, collapse = ", ")
      ),
      c("dpprior_interval_structure_error", "dpprior_choice_error"),
      paste0("K_interval$", name), x, paste(choices, collapse = ", "),
      "choice"
    )
  }
  x
}

# Validate and canonicalize interval syntax without constructing a PMF.
.dp_validate_K_interval <- function(K_interval, J) {
  assert_valid_J(J)
  J <- as.integer(J)

  if (is.numeric(K_interval) && length(K_interval) == 2L) {
    .dp_interval_ambiguous(K_interval, J)
  }
  if (!is.list(K_interval) || is.null(names(K_interval)) ||
      any(!nzchar(names(K_interval))) || anyDuplicated(names(K_interval))) {
    .dp_interval_abort(
      "K_interval must be a uniquely named list; a bare range is not enough.",
      c("dpprior_interval_structure_error", "dpprior_type_error"),
      value = K_interval,
      expected = paste0(
        "list(lower=..., upper=..., type=..., coverage=..., family=...)"
      ),
      code = "structure"
    )
  }

  allowed_fields <- c("lower", "upper", "type", "coverage", "family", "mu_K")
  unknown <- setdiff(names(K_interval), allowed_fields)
  if (length(unknown)) {
    .dp_interval_abort(
      sprintf("Unknown K_interval field(s): %s.", paste(unknown, collapse = ", ")),
      c("dpprior_interval_structure_error", "dpprior_unknown_argument_error"),
      value = K_interval,
      expected = paste(allowed_fields, collapse = ", "),
      code = "unknown_fields",
      unknown = unknown
    )
  }

  # A "plausible" range is ambiguous by definition. Detect that semantic
  # request before demanding an otherwise irrelevant family or coverage so
  # the same ambiguity condition is raised for every plausible-range spelling.
  type <- NULL
  if ("type" %in% names(K_interval)) {
    type <- .dp_interval_scalar_character(K_interval$type, "type")
    if (type %in% c("plausible", "plausible_range")) {
      .dp_interval_ambiguous(K_interval, J, code = "plausible_range")
    }
  }

  required <- c("lower", "upper", "type", "family")
  missing <- required[!required %in% names(K_interval)]
  if (length(missing)) {
    .dp_interval_abort(
      sprintf("K_interval is missing required field(s): %s.",
              paste(missing, collapse = ", ")),
      c("dpprior_interval_under_specified", "dpprior_interval_structure_error"),
      value = K_interval,
      expected = paste(required, collapse = ", "),
      code = "missing_fields",
      missing = missing
    )
  }

  type <- .dp_interval_scalar_character(type, "type", .DP_INTERVAL_TYPES)
  family <- .dp_interval_scalar_character(K_interval$family, "family")
  if (!(family %in% .DP_INTERVAL_FAMILIES)) {
    .dp_interval_abort(
      sprintf(
        "K_interval$family='%s' is unsupported; choose maxent, or provide a strict target_pmf.",
        family
      ),
      c("dpprior_interval_family_error", "dpprior_choice_error"),
      "K_interval$family", family,
      paste(.DP_INTERVAL_FAMILIES, collapse = ", "), "family"
    )
  }

  validate_endpoint <- function(x, name) {
    .dpprior_validate_count(
      x, paste0("K_interval$", name), minimum = 1L, maximum = J,
      .subclass = "dpprior_interval_support_error"
    )
  }
  lower <- validate_endpoint(K_interval$lower, "lower")
  upper <- validate_endpoint(K_interval$upper, "upper")
  if (lower > upper) {
    .dp_interval_abort(
      "K_interval endpoints are reversed: lower must be <= upper.",
      c("dpprior_interval_support_error", "dpprior_bounds_error"),
      value = c(lower = lower, upper = upper),
      expected = sprintf("1 <= lower <= upper <= %d", J),
      code = "reversed_endpoints"
    )
  }

  has_coverage <- "coverage" %in% names(K_interval)
  if (identical(type, "hard_bounds")) {
    if (has_coverage) {
      .dp_interval_abort(
        paste0(
          "coverage is forbidden for type='hard_bounds'; hard bounds already ",
          "mean zero mass outside the inclusive interval."
        ),
        c("dpprior_interval_semantic_conflict_error",
          "dpprior_conflicting_input"),
        "K_interval$coverage", K_interval$coverage,
        "omit coverage for hard_bounds", "hard_bounds_coverage"
      )
    }
    coverage <- 1
  } else {
    if (!has_coverage || is.null(K_interval$coverage)) {
      .dp_interval_abort(
        sprintf("coverage is required for type='%s'; there is no default.", type),
        c("dpprior_interval_missing_coverage_error",
          "dpprior_interval_under_specified"),
        "K_interval$coverage", NULL, "one scalar in (0, 1)",
        "missing_coverage"
      )
    }
    if (.dpprior_is_plain_numeric(K_interval$coverage) &&
        length(K_interval$coverage) == 1L &&
        !is.na(K_interval$coverage) && is.finite(K_interval$coverage) &&
        K_interval$coverage == 1) {
      .dp_interval_abort(
        "coverage=1 changes the assertion to hard bounds; use type='hard_bounds' and omit coverage.",
        c("dpprior_interval_semantic_conflict_error",
          "dpprior_conflicting_input"),
        "K_interval$coverage", K_interval$coverage,
        "one scalar in (0, 1)", "probability_one"
      )
    }
    coverage <- tryCatch(
      .dpprior_validate_probability(
        K_interval$coverage, "K_interval$coverage", open = TRUE,
        .subclass = "dpprior_interval_probability_error"
      ),
      dpprior_invalid_input = function(condition) stop(condition)
    )
  }

  mu <- NULL
  if ("mu_K" %in% names(K_interval) && !is.null(K_interval$mu_K)) {
    mu <- .dpprior_validate_scalar(
      K_interval$mu_K, "K_interval$mu_K", lower = 1, upper = J,
      .subclass = "dpprior_interval_mean_error"
    )
  }

  list(
    lower = lower,
    upper = upper,
    type = type,
    coverage = coverage,
    family = family,
    mu_K = mu,
    support = c(lower = 1L, upper = J),
    endpoints = "inclusive"
  )
}

.dp_K_pmf_metrics <- function(pmf) {
  support <- seq_along(pmf)
  mean_K <- sum(support * pmf)
  variance_K <- sum((support - mean_K)^2 * pmf)
  probs <- c(`2.5%` = 0.025, `25%` = 0.25, `50%` = 0.5,
             `75%` = 0.75, `97.5%` = 0.975)
  cdf <- cumsum(pmf)
  quantiles <- vapply(probs, function(prob) which(cdf >= prob)[1L], integer(1L))
  list(
    mean = mean_K,
    variance = variance_K,
    sd = sqrt(variance_K),
    cv = if (mean_K > 0) sqrt(variance_K) / mean_K else NA_real_,
    quantiles = quantiles
  )
}

# Separate direct postcondition implementation used only for verification.
# It intentionally avoids .dp_K_pmf_metrics() so a helper-level regression
# cannot make selected and verification metrics agree tautologically.
.dp_K_pmf_metrics_verify <- function(pmf) {
  support <- seq_along(pmf)
  mean_K <- drop(crossprod(pmf, support))
  centered <- support - mean_K
  variance_K <- drop(crossprod(pmf, centered * centered))
  probs <- c(`2.5%` = 0.025, `25%` = 0.25, `50%` = 0.5,
             `75%` = 0.75, `97.5%` = 0.975)
  cumulative <- cumsum(pmf)
  quantiles <- vapply(
    probs,
    function(probability) min(support[cumulative >= probability]),
    integer(1L)
  )
  list(
    mean = mean_K,
    variance = variance_K,
    sd = sqrt(variance_K),
    cv = if (mean_K > 0) sqrt(variance_K) / mean_K else NA_real_,
    quantiles = quantiles
  )
}

.dp_target_plain_view_value <- function(value) {
  if (is.null(value)) return(NULL)
  if (inherits(value, "condition")) {
    return(list(
      classes = as.character(class(value)),
      fields = .dp_target_plain_view_value(unclass(value))
    ))
  }
  if (is.language(value)) {
    return(paste(deparse(value, width.cutoff = 500L), collapse = " "))
  }
  if (is.environment(value) || is.function(value) ||
      typeof(value) %in% c("externalptr", "weakref")) {
    return(sprintf("<%s omitted from plain compatibility view>",
                   typeof(value)))
  }
  if (is.list(value)) {
    raw <- unclass(value)
    if (!length(raw)) return(list())
    record_names <- names(raw)
    if (is.null(record_names)) {
      record_names <- sprintf("item_%03d", seq_along(raw))
    } else {
      invalid <- is.na(record_names) | !nzchar(record_names)
      record_names[invalid] <- sprintf("item_%03d", which(invalid))
      record_names <- make.unique(record_names, sep = "_")
    }
    out <- lapply(raw, .dp_target_plain_view_value)
    names(out) <- record_names
    return(out)
  }
  if (!is.atomic(value) || !length(value)) return(NULL)
  original_names <- names(value)
  if (is.factor(value)) {
    value <- as.character(value)
  } else if (is.object(value) || !is.null(dim(value))) {
    value <- as.vector(value)
  }
  if (!(is.numeric(value) || is.logical(value) || is.character(value))) {
    value <- as.character(value)
  }
  invalid <- is.na(value) | (is.numeric(value) & !is.finite(value))
  if (any(invalid)) {
    type <- typeof(value)
    value <- vapply(seq_along(value), function(index) {
      item <- value[[index]]
      if (is.na(item)) return(paste0("NA_", type, "_"))
      if (is.numeric(item) && is.infinite(item)) {
        return(if (item > 0) "Inf" else "-Inf")
      }
      format(item, digits = 17L, scientific = TRUE, trim = TRUE)
    }, character(1L))
  }
  attributes(value) <- NULL
  if (!is.null(original_names) && length(original_names) == length(value) &&
      !anyNA(original_names) && all(nzchar(original_names)) &&
      !anyDuplicated(original_names)) {
    names(value) <- original_names
  }
  value
}

.dp_target_condition_view <- function(condition) {
  if (is.null(condition)) return(NULL)
  .dp_target_plain_view_value(condition)
}

.dp_target_attempt_v1 <- function(attempt, index) {
  record_or_null <- function(value, numeric_only = FALSE) {
    if (is.null(value)) return(NULL)
    tryCatch(
      {
        .dpprior_schema_validate_plain_record_value(
          value, "target_attempt_evidence", numeric_only = numeric_only
        )
        value
      },
      dpprior_schema_error = function(condition) NULL
    )
  }
  available <- function(value) {
    !is.null(value) && length(value) > 0L && !anyNA(value) &&
      (!is.numeric(value) || all(is.finite(value)))
  }
  count_or_null <- function(value) {
    if (!available(value) || !is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value < 0 || value != floor(value)) {
      return(NULL)
    }
    as.integer(value)
  }
  numeric_or_null <- function(value) {
    if (!available(value) || !is.numeric(value) || length(value) != 1L ||
        !is.finite(value)) {
      return(NULL)
    }
    unname(value)
  }
  start <- record_or_null(attempt$start, numeric_only = TRUE)
  bounds <- record_or_null(attempt$bounds, numeric_only = TRUE)
  control <- record_or_null(attempt$control)
  exit_code <- count_or_null(attempt$exit_code)
  iterations <- count_or_null(attempt$iterations)
  evaluation_count <- count_or_null(attempt$evaluations)
  evaluations <- if (is.null(evaluation_count)) NULL else
    list(total = evaluation_count)
  candidate_objective <- numeric_or_null(attempt$candidate_objective)
  elapsed_seconds <- numeric_or_null(attempt$elapsed)
  unavailable_values <- list(
    start = start, bounds = bounds, control = control, exit_code = exit_code,
    iterations = iterations, evaluations = evaluations,
    candidate_parameters = NULL, candidate_objective = candidate_objective,
    elapsed_seconds = elapsed_seconds
  )
  unavailable <- vapply(
    names(unavailable_values)[vapply(unavailable_values, is.null, logical(1))],
    function(field) paste0(field, "_not_recorded_by_target_constructor"),
    character(1L)
  )
  error_message <- if (is.character(attempt$error) &&
      length(attempt$error) == 1L && !is.na(attempt$error) &&
      nzchar(attempt$error)) attempt$error else NULL
  error_record <- if (is.null(error_message)) NULL else list(
    class = if (is.character(attempt$error_class) &&
        length(attempt$error_class) > 0L && nzchar(attempt$error_class[[1L]])) {
      attempt$error_class[[1L]]
    } else {
      "dpprior_target_attempt_error"
    },
    code = "target_attempt_error",
    message = error_message
  )
  warnings <- if (is.character(attempt$warnings)) {
    unique(attempt$warnings[!is.na(attempt$warnings) & nzchar(attempt$warnings)])
  } else {
    character()
  }
  exit_success <- !is.null(exit_code) && identical(exit_code, 0L)
  .dpprior_new_attempt(
    id = sprintf("target-%03d", index),
    stage = "target_construction",
    method = if (is.character(attempt$method) &&
        length(attempt$method) == 1L && !is.na(attempt$method) &&
        nzchar(attempt$method)) attempt$method else "target_constructor",
    start = start,
    bounds = bounds,
    control = control,
    exit_code = exit_code,
    message = if (is.character(attempt$message) &&
        length(attempt$message) == 1L && !is.na(attempt$message)) {
      attempt$message
    } else {
      ""
    },
    iterations = iterations,
    evaluations = evaluations,
    candidate_parameters = NULL,
    candidate_objective = candidate_objective,
    elapsed_seconds = elapsed_seconds,
    warnings = warnings,
    error = error_record,
    selected = FALSE,
    reason_code = if (exit_success) "candidate_evaluated" else
      "candidate_not_selected",
    unavailable = unavailable
  )
}

.dp_target_snapshot_v1 <- function(implied, interval, pmf, residuals,
                                   tolerances, source) {
  .dpprior_new_snapshot(
    parameters = NULL,
    M = NULL,
    achieved = list(implied = implied, interval = interval, pmf = pmf),
    residuals = residuals,
    tolerances = tolerances,
    finite = if (is.null(pmf)) TRUE else all(is.finite(pmf)),
    source = source
  )
}

.dp_target_derivation_entry <- function(rule, before, after, evidence,
                                        outcome = "canonicalized",
                                        opt_in = FALSE) {
  list(
    rule = rule, outcome = outcome, opt_in = opt_in,
    before = before, after = after, evidence = evidence
  )
}

.dp_target_orders_v1 <- function() {
  .dpprior_new_orders(
    requested_reason = "not_applicable_to_finite_support_target",
    selected_reason = "not_applicable_to_finite_support_target",
    verification_required_reason =
      "not_applicable_to_finite_support_target",
    verification_used_reason =
      "not_applicable_to_finite_support_target"
  )
}

.dp_target_settings_v1 <- function(method, controls) {
  list(
    method = method,
    controls = controls,
    parameterization = "bounded_discrete_K_target"
  )
}

.dp_target_interval_truth_v1 <- function(tolerance, root_control) {
  tolerances <- list(
    constraint = tolerance,
    pmf_l1 = .TOL_PMF_SUM,
    moment_relative = 1e-8,
    moment_scale_floor = 1
  )
  controls <- list(
    constraint_tolerance = tolerances$constraint,
    root_tolerance = root_control$root_tol,
    max_iterations = as.integer(root_control$max_iterations),
    pmf_l1_tolerance = tolerances$pmf_l1,
    moment_relative_tolerance = tolerances$moment_relative,
    moment_scale_floor = tolerances$moment_scale_floor
  )
  list(tolerances = tolerances, controls = controls)
}

.dp_target_interval_masses_v1 <- function(pmf, interval) {
  masses <- .dpprior_target_interval_masses(pmf, interval)
  masses[abs(masses) <= .TOL_PMF_SUM] <- 0
  masses[abs(masses - 1) <= .TOL_PMF_SUM] <- 1
  masses
}

.dp_target_interval_residuals_v1 <- function(pmf, interval) {
  achieved <- .dp_target_interval_masses_v1(pmf, interval)
  requested <- .dp_interval_requested_masses(interval)
  constrained <- is.finite(requested)
  interval_residual <- achieved[constrained] - requested[constrained]
  list(
    mean = if (is.null(interval$mu_K)) NULL else
      sum(seq_along(pmf) * pmf) - interval$mu_K,
    interval = as.list(interval_residual)
  )
}

.dp_target_interval_evidence_v1 <- function(x, controls) {
  groups <- .dp_interval_groups(x$interval, x$J)
  masses <- groups$masses
  positive <- names(masses)[masses > 0]
  lower <- sum(vapply(
    positive,
    function(group) masses[[group]] * min(groups$values[[group]]),
    numeric(1L)
  ))
  upper <- sum(vapply(
    positive,
    function(group) masses[[group]] * max(groups$values[[group]]),
    numeric(1L)
  ))
  requested_mean <- x$interval$mu_K
  feasibility_tolerance <- 64 * .Machine$double.eps * max(
    1, abs(lower), abs(upper),
    if (is.null(requested_mean)) 1 else abs(requested_mean)
  )
  boundary <- x$parameters$boundary
  list(
    constructor_method = if (!is.null(boundary)) {
      "boundary"
    } else if (is.null(requested_mean)) {
      "analytic"
    } else {
      "uniroot"
    },
    group_masses = as.list(masses),
    common_tilt = if (is.null(boundary)) x$parameters$common_tilt else NULL,
    boundary = list(active = !is.null(boundary), side = boundary),
    mean_constraint = requested_mean,
    feasible_mean_lower = lower,
    feasible_mean_upper = upper,
    feasibility_tolerance = feasibility_tolerance,
    root_tolerance = controls$root_tolerance
  )
}

.dp_target_interval_verification_v1 <- function(
    x, implied, achieved_interval, residuals, tolerances, controls) {
  groups <- .dp_interval_groups(x$interval, x$J)
  verifier_pmf <- .dp_interval_group_pmf_verify(
    x$J, groups,
    tilt = if (is.null(x$parameters$common_tilt)) 0 else
      x$parameters$common_tilt,
    boundary = x$parameters$boundary
  )
  verifier_pmf <- unname(as.numeric(verifier_pmf))
  verifier_implied_raw <- .dp_K_pmf_metrics_verify(verifier_pmf)
  verifier_implied <- list(
    mean = verifier_implied_raw$mean,
    variance = verifier_implied_raw$variance
  )
  verifier_interval <- as.list(
    .dp_target_interval_masses_v1(verifier_pmf, x$interval)
  )
  verifier_residuals <- .dp_target_interval_residuals_v1(
    verifier_pmf, x$interval
  )
  selected <- .dp_target_snapshot_v1(
    implied, achieved_interval, x$pmf, residuals, tolerances,
    "target_constructor"
  )
  verifier <- .dp_target_snapshot_v1(
    verifier_implied, verifier_interval, verifier_pmf,
    verifier_residuals, tolerances, "independent_target_reconstruction"
  )
  provisional <- list(
    kind = "interval", interval = x$interval, tolerances = tolerances,
    verification = list(
      selected_snapshot = selected, verifier_snapshot = verifier
    )
  )
  stability <- .dpprior_expected_target_stability(provisional)
  check_source <- "independent_target_reconstruction"
  components <- list(
    target_reconstruction = .dpprior_new_check(
      value = sum(abs(verifier_pmf - x$pmf)), reference = 0,
      tolerance = tolerances$pmf_l1, operator = "lte",
      source = check_source
    ),
    order_stability = .dpprior_new_check(
      value = stability$delta,
      reference = setNames(rep(0, length(stability$delta)),
                           names(stability$delta)),
      tolerance = stability$tolerance, operator = "lte",
      source = check_source
    )
  )
  if (!isTRUE(x$verified)) {
    components$constructor_postcondition <- .dpprior_new_check(
      value = FALSE, reference = TRUE, tolerance = NULL,
      operator = "identical", source = check_source
    )
  }
  invariants <- list(
    support_identity = .dpprior_new_check(
      value = c(
        target_support = identical(as.integer(x$support), seq_len(x$J)),
        selected_length = length(x$pmf) == x$J,
        verifier_length = length(verifier_pmf) == x$J
      ),
      reference = c(
        target_support = TRUE, selected_length = TRUE,
        verifier_length = TRUE
      ),
      tolerance = NULL, operator = "identical", source = check_source
    ),
    authority_identity = .dpprior_new_check(
      value = c(
        selected_pmf = TRUE, used_pmf = TRUE,
        normalized_pmf_unidentified = TRUE, request_J = TRUE
      ),
      reference = c(
        selected_pmf = TRUE, used_pmf = TRUE,
        normalized_pmf_unidentified = TRUE, request_J = TRUE
      ),
      tolerance = NULL, operator = "identical", source = check_source
    )
  )
  .dpprior_new_verification(
    method = "independent_target_reconstruction",
    performed = TRUE,
    passed = isTRUE(x$verified),
    reason = if (isTRUE(x$verified)) "verified" else x$message,
    settings = list(
      pmf_l1_tolerance = tolerances$pmf_l1,
      moment_relative_tolerance = tolerances$moment_relative,
      moment_scale_floor = tolerances$moment_scale_floor,
      constraint_tolerance = tolerances$constraint,
      root_tolerance = controls$root_tolerance,
      max_iterations = controls$max_iterations
    ),
    selected_snapshot = selected,
    verifier_snapshot = verifier,
    stability = stability,
    components = components,
    invariants = invariants
  )
}

.dp_target_simple_verification_v1 <- function(
    x, kind, implied, residuals, tolerances, request, normalized, used) {
  selected <- .dp_target_snapshot_v1(
    implied, NULL, x$pmf, residuals, tolerances, "target_constructor"
  )
  independent_implied <- implied
  if (identical(kind, "pmf")) {
    recomputed <- .dp_K_pmf_metrics_verify(x$pmf)
    moment_delta <- c(
      mean = recomputed$mean - implied$mean,
      variance = recomputed$variance - implied$variance
    )
    moment_tolerance <- x$tolerances$moment_check * pmax(
      1, abs(unlist(implied, use.names = FALSE)),
      abs(c(recomputed$mean, recomputed$variance))
    )
  } else {
    reference_variance <- if (identical(x$provenance$source,
                                        "qualitative_confidence")) {
      x$provenance$vif * (implied$mean - 1)
    } else if (identical(x$provenance$source,
                         "coefficient_of_variation")) {
      (x$provenance$raw_cv_K * implied$mean)^2
    } else {
      request$var_K
    }
    moment_delta <- c(
      mean = implied$mean - if ("mu_K" %in% names(request)) {
        request$mu_K
      } else {
        request$mean
      },
      variance = implied$variance - reference_variance
    )
    moment_tolerance <- c(mean = 0, variance = 0)
  }
  verifier <- .dp_target_snapshot_v1(
    independent_implied, NULL, x$pmf,
    list(independent_moment_delta = moment_delta), tolerances,
    "independent_target_verification"
  )
  upper <- (implied$mean - 1) * (x$J - implied$mean)
  upper_tolerance <- 1e-8 * max(1, abs(upper), abs(implied$variance))
  components <- list(
    moment_reconstruction = .dpprior_new_check(
      value = moment_delta,
      reference = setNames(rep(0, length(moment_delta)), names(moment_delta)),
      tolerance = moment_tolerance,
      operator = "abs_lte",
      source = "independent_target_verification"
    ),
    finite_support_feasibility = .dpprior_new_check(
      value = implied$variance, reference = upper,
      tolerance = upper_tolerance, operator = "lte",
      source = "analytic_finite_support_moment_bound"
    )
  )
  if (identical(kind, "pmf")) {
    components$pmf_mass <- .dpprior_new_check(
      value = sum(x$pmf), reference = 1, tolerance = .TOL_PMF_SUM,
      operator = "abs_lte", source = "independent_target_verification"
    )
  }
  invariants <- list(
    support_identity = .dpprior_new_check(
      value = identical(as.integer(x$support), seq_len(x$J)),
      reference = TRUE, tolerance = NULL, operator = "identical",
      source = "independent_target_verification"
    ),
    authority_identity = .dpprior_new_check(
      value = c(
        request_J = identical(request$J, as.integer(x$J)),
        normalized_J = identical(normalized$J, as.integer(x$J)),
        used_J = identical(used$J, as.integer(x$J)),
        pmf_identity = is.null(x$pmf) ||
          identical(used$pmf, x$pmf)
      ),
      reference = c(
        request_J = TRUE, normalized_J = TRUE, used_J = TRUE,
        pmf_identity = TRUE
      ),
      tolerance = NULL, operator = "identical",
      source = "independent_target_verification"
    )
  )
  .dpprior_new_verification(
    method = "independent_finite_support_target_verification",
    performed = TRUE, passed = isTRUE(x$verified),
    reason = if (isTRUE(x$verified)) "verified" else x$message,
    settings = list(
      support = c(lower = 1L, upper = as.integer(x$J)),
      pmf_l1_tolerance = .TOL_PMF_SUM,
      moment_relative_tolerance = if (identical(kind, "pmf")) {
        x$tolerances$moment_check
      } else {
        0
      }
    ),
    selected_snapshot = selected,
    verifier_snapshot = verifier,
    stability = NULL,
    components = components,
    invariants = invariants
  )
}

.dp_target_infeasible_v1 <- function(x, request, normalized, family,
                                     tolerances, controls, settings) {
  support <- seq_len(x$J)
  groups_record <- .dp_interval_groups(x$interval, x$J)
  groups <- groups_record$values
  group_masses <- groups_record$masses
  positive_groups <- names(group_masses)[group_masses > 0]
  group_counts <- vapply(groups, length, integer(1L))
  empty_groups <- positive_groups[group_counts[positive_groups] == 0L]
  requested_mean <- x$interval$mu_K
  lower_bound <- if (!length(empty_groups)) sum(vapply(
    positive_groups,
    function(group) group_masses[[group]] * min(groups[[group]]),
    numeric(1L)
  )) else NULL
  upper_bound <- if (!length(empty_groups)) sum(vapply(
    positive_groups,
    function(group) group_masses[[group]] * max(groups[[group]]),
    numeric(1L)
  )) else NULL
  feasibility_tolerance <- 64 * .Machine$double.eps * max(
    1,
    if (is.null(lower_bound)) 0 else abs(lower_bound),
    if (is.null(upper_bound)) 0 else abs(upper_bound),
    if (is.null(requested_mean)) 0 else abs(requested_mean)
  )
  side <- if (!is.null(requested_mean) && !is.null(lower_bound) &&
              requested_mean < lower_bound - feasibility_tolerance) {
    "below"
  } else if (!is.null(requested_mean) && !is.null(upper_bound) &&
             requested_mean > upper_bound + feasibility_tolerance) {
    "above"
  } else {
    NULL
  }
  outside_distance <- if (identical(side, "below")) {
    lower_bound - requested_mean
  } else if (identical(side, "above")) {
    requested_mean - upper_bound
  } else {
    NULL
  }
  assumptions <- list(estimand = "K_J", construction = "maxent")
  certificate <- list(
    version = "1",
    method = "analytic_interval_group_support_feasibility",
    kind = if (length(empty_groups)) {
      "positive_mass_group_has_empty_support"
    } else {
      "mean_outside_group_mass_hull"
    },
    assumptions = assumptions,
    request = request,
    J = as.integer(x$J),
    support = support,
    interval = x$interval,
    family = family,
    group_masses = as.list(group_masses),
    group_support_counts = as.list(group_counts),
    empty_groups = if (length(empty_groups)) empty_groups else NULL,
    requested_mean = requested_mean,
    feasible_mean_lower = lower_bound,
    feasible_mean_upper = upper_bound,
    feasibility_tolerance = feasibility_tolerance,
    side = side,
    outside_distance = outside_distance,
    certified = TRUE,
    source = "analytic_maxent_interval_feasibility"
  )
  check_source <- "analytic_maxent_interval_feasibility"
  verification <- .dpprior_new_verification(
    method = "analytic_interval_infeasibility_certificate",
    performed = TRUE, passed = TRUE, reason = "certified_infeasible",
    settings = list(certificate = certificate),
    selected_snapshot = NULL, verifier_snapshot = NULL, stability = NULL,
    components = list(
      infeasibility_certificate = .dpprior_new_check(
        value = if (length(empty_groups)) {
          as.integer(length(empty_groups))
        } else {
          outside_distance
        },
        reference = 0,
        tolerance = if (length(empty_groups)) 0 else feasibility_tolerance,
        operator = "gt", source = check_source
      )
    ),
    invariants = list(
      request_identity = .dpprior_new_check(
        value = c(J = TRUE, interval = TRUE, mean = TRUE),
        reference = c(J = TRUE, interval = TRUE, mean = TRUE),
        tolerance = NULL, operator = "identical", source = check_source
      ),
      support_identity = .dpprior_new_check(
        value = c(top = TRUE, interval = TRUE),
        reference = c(top = TRUE, interval = TRUE),
        tolerance = NULL, operator = "identical", source = check_source
      )
    )
  )
  attempt <- .dpprior_new_attempt(
    id = "target-interval-feasibility-001",
    stage = "feasibility",
    method = "analytic_interval_group_support_feasibility",
    start = list(
      J = as.integer(x$J), requested_mean = requested_mean,
      interval_lower = x$interval$lower,
      interval_upper = x$interval$upper,
      coverage = x$interval$coverage
    ),
    bounds = list(
      feasible_mean_lower = lower_bound,
      feasible_mean_upper = upper_bound
    ),
    control = list(feasibility_tolerance = feasibility_tolerance),
    exit_code = 0L,
    message = "analytic interval infeasibility certified",
    iterations = 0L,
    evaluations = list(function_count = 1L),
    candidate_parameters = NULL,
    candidate_objective = NULL,
    elapsed_seconds = 0,
    warnings = character(), error = NULL, selected = FALSE,
    reason_code = "globally_infeasible_by_certificate",
    unavailable = c(
      candidate_parameters =
        "analytic certificate route has no candidate parameters",
      candidate_objective =
        "analytic certificate route has no candidate objective"
    )
  )
  computation <- .dpprior_new_computation(
    request = settings, used = settings,
    orders = .dp_target_orders_v1(), scaling = .dpprior_new_scaling(),
    attempts = list(attempt), candidate_evaluations = list(),
    selected_candidate_id = NULL, selected_attempt_id = NULL,
    fallback = .dpprior_new_fallback(),
    termination = .dpprior_new_termination(
      code = "certified_infeasible", message = x$message,
      source = "analytic_certificate", iterations = NULL
    ),
    resources = list()
  )
  list(
    assumptions = assumptions,
    residuals = list(
      unavailable_reason = "analytic interval construction infeasible"
    ),
    verification = verification,
    computation = computation
  )
}

.dp_K_target_class <- function(x) {
  stopifnot(x$status %in% .DP_TARGET_STATUS_CODES)
  legacy_x <- x
  J <- as.integer(x$J)
  source <- x$provenance$source
  target_kind <- if (identical(source, "coefficient_of_variation")) {
    "cv"
  } else if (identical(source, "target_pmf")) {
    "pmf"
  } else if (identical(source, "K_interval")) {
    "interval"
  } else {
    "moments"
  }
  if (identical(target_kind, "interval") &&
      identical(x$status, "approximate")) {
    # A retained candidate that failed its declared construction
    # postcondition is evidence, not an authoritative canonical target PMF.
    # The full candidate remains in compatibility.views.target_v0.
    x$pmf <- NULL
    x$implied <- NULL
    x$achieved_interval <- NULL
    x$parameters <- NULL
  }
  implied <- if (is.null(x$implied)) NULL else list(
    mean = x$implied$mean, variance = x$implied$variance
  )
  interval <- if (identical(target_kind, "interval")) list(
    lower = x$interval$lower,
    upper = x$interval$upper,
    type = x$interval$type,
    coverage = x$interval$coverage,
    family = x$interval$family,
    mu_K = x$interval$mu_K,
    support = x$interval$support,
    endpoints = x$interval$endpoints
  ) else NULL
  if (!is.null(interval) && is.null(interval$mu_K) &&
      identical(interval$type, "hard_bounds") &&
      identical(interval$lower, interval$upper)) {
    interval$mu_K <- as.numeric(interval$lower)
  }
  if (identical(target_kind, "interval")) x$interval <- interval
  family <- if (identical(target_kind, "interval")) x$family else NULL

  if (identical(source, "direct_variance")) {
    request <- list(J = J, mu_K = implied$mean, var_K = implied$variance)
    normalized <- c(request, list(interval = NULL, pmf = NULL))
    used <- normalized
    first_rule <- "canonicalize_direct_moments"
    second <- NULL
    method <- "target_moments"
  } else if (identical(source, "qualitative_confidence")) {
    request <- list(
      J = J, mean = implied$mean, confidence = x$provenance$confidence
    )
    normalized <- c(request, list(interval = NULL, pmf = NULL))
    used <- list(
      J = J, mean = implied$mean, variance = implied$variance,
      interval = NULL, pmf = NULL
    )
    first_rule <- "canonicalize_confidence_target"
    second <- .dp_target_derivation_entry(
      "derive_variance_from_confidence_vif", normalized, used,
      list(
        confidence = x$provenance$confidence,
        vif = x$provenance$vif,
        formula = "variance = vif * (mean - 1)"
      ),
      outcome = "derived"
    )
    method <- "target_confidence_vif"
  } else if (identical(source, "coefficient_of_variation")) {
    request <- list(
      J = J, mean = implied$mean, cv = x$provenance$raw_cv_K
    )
    normalized <- c(request, list(interval = NULL, pmf = NULL))
    used <- list(
      J = J, mean = implied$mean, variance = implied$variance,
      interval = NULL, pmf = NULL
    )
    first_rule <- "canonicalize_cv_target"
    second <- .dp_target_derivation_entry(
      "derive_variance_from_cv", normalized, used,
      list(
        cv = x$provenance$raw_cv_K,
        definition = "SD(K_J) / E(K_J)",
        formula = "variance = (cv * mean)^2"
      ),
      outcome = "derived"
    )
    method <- "target_cv"
  } else if (identical(source, "target_pmf")) {
    input_pmf <- unname(as.numeric(x$request$target_pmf))
    request <- list(J = J, pmf = input_pmf)
    normalized <- list(J = J, pmf = unname(x$pmf), interval = NULL)
    used <- normalized
    first_rule <- if (isTRUE(x$provenance$k0_entry_removed)) {
      "drop_structural_k0_zero"
    } else {
      "validate_strict_pmf"
    }
    second <- NULL
    method <- "target_pmf"
  } else {
    request <- list(J = J, K_interval = interval, mu_K = interval$mu_K)
    normalized <- list(
      J = J, interval = interval, family = family, pmf = NULL
    )
    used <- if (is.null(x$pmf)) normalized else list(
      J = J, interval = interval, family = family, pmf = unname(x$pmf)
    )
    first_rule <- "canonicalize_interval_request"
    method <- "target_interval_maxent"
  }
  first_evidence <- switch(
    first_rule,
    canonicalize_direct_moments = list(
      source = "validated_public_arguments",
      original_request = "compatibility.views.target_v0.request"
    ),
    canonicalize_confidence_target = list(
      source = "validated_public_arguments",
      original_request = "compatibility.views.target_v0.request"
    ),
    canonicalize_cv_target = list(
      source = "validated_public_arguments",
      original_request = "compatibility.views.target_v0.request"
    ),
    validate_strict_pmf = list(
      source = "strict_pmf_validation", input_length = as.integer(length(input_pmf)),
      normalization = "none"
    ),
    drop_structural_k0_zero = list(
      source = "strict_pmf_validation", input_length = as.integer(length(input_pmf)),
      removed_entry = "K=0 exact structural zero"
    ),
    canonicalize_interval_request = list(
      source = "validated_interval_request",
      original_request = "compatibility.views.target_v0.request"
    )
  )
  if (identical(first_rule, "canonicalize_interval_request") &&
      is.null(legacy_x$interval$mu_K) && !is.null(interval$mu_K)) {
    first_evidence$singleton_support_implied_mean <- TRUE
    first_evidence$implied_mean <- interval$mu_K
  }

  root_control <- x$.canonical_controls
  if (identical(target_kind, "interval")) {
    if (is.null(root_control)) {
      root_control <- list(
        root_tol = max(
          .Machine$double.eps,
          min(1e-12, x$tolerances$constraint / 10)
        ),
        max_iterations = 1000L
      )
    }
    truth <- .dp_target_interval_truth_v1(
      x$tolerances$constraint, root_control
    )
    tolerances <- truth$tolerances
    controls <- truth$controls
    if (!is.null(x$pmf)) {
      second <- .dp_target_derivation_entry(
        switch(
          interval$type,
          hard_bounds = "construct_maxent_hard_bounds_pmf",
          equal_tail = "construct_maxent_equal_tail_pmf",
          central_mass = "construct_maxent_central_mass_pmf"
        ),
        normalized, used,
        .dp_target_interval_evidence_v1(x, controls),
        outcome = "derived"
      )
    } else {
      second <- NULL
    }
  } else {
    controls <- if (identical(target_kind, "pmf")) {
      list(
        pmf_l1_tolerance = .TOL_PMF_SUM,
        moment_check_tolerance = x$tolerances$moment_check
      )
    } else {
      list()
    }
    tolerances <- if (identical(target_kind, "pmf")) {
      list(
        pmf_sum = x$tolerances$pmf_sum,
        moment_check = x$tolerances$moment_check
      )
    } else {
      list(moment = x$tolerances$moment)
    }
  }
  derivation <- list(
    request_to_normalized = .dp_target_derivation_entry(
      first_rule, request, normalized, first_evidence
    ),
    normalized_to_used = second
  )
  settings <- .dp_target_settings_v1(method, controls)

  if (identical(x$status, "infeasible")) {
    infeasible <- .dp_target_infeasible_v1(
      x, request, normalized, family, tolerances, controls, settings
    )
    assumptions <- infeasible$assumptions
    residuals <- infeasible$residuals
    computation <- infeasible$computation
    verification <- infeasible$verification
  } else {
    assumptions <- if (identical(target_kind, "interval")) {
      list(estimand = "K_J", construction = "maxent")
    } else if (identical(target_kind, "pmf")) {
      list(estimand = "K_J", pmf_authoritative = TRUE)
    } else {
      list(estimand = "K_J", pmf_identified = FALSE)
    }
    residuals <- if (identical(target_kind, "interval")) {
      if (is.null(x$pmf)) {
        list(unavailable_reason = "target construction produced no PMF")
      } else {
        .dp_target_interval_residuals_v1(x$pmf, interval)
      }
    } else if (identical(target_kind, "pmf")) {
      list(pmf_sum = sum(x$pmf) - 1)
    } else {
      list(mean = 0, variance = 0)
    }
    attempts <- lapply(seq_along(x$attempts), function(index) {
      .dp_target_attempt_v1(x$attempts[[index]], index)
    })
    iterations <- if (length(attempts) == 1L) attempts[[1L]]$iterations else
      NULL
    computation <- .dpprior_new_computation(
      request = settings, used = settings,
      orders = .dp_target_orders_v1(), scaling = .dpprior_new_scaling(),
      attempts = attempts, candidate_evaluations = list(),
      selected_candidate_id = NULL, selected_attempt_id = NULL,
      fallback = .dpprior_new_fallback(),
      termination = .dpprior_new_termination(
        code = x$status, message = x$message,
        source = "bounded_discrete_target_constructor",
        iterations = iterations,
        boundary_reason = if (identical(x$status, "boundary")) {
          "analytic_feasibility_boundary"
        } else {
          NULL
        }
      ),
      resources = list()
    )
    achieved_interval <- if (identical(target_kind, "interval") &&
        !is.null(x$pmf)) {
      as.list(.dp_target_interval_masses_v1(x$pmf, interval))
    } else {
      NULL
    }
    if (identical(target_kind, "interval") && !is.null(x$pmf)) {
      verification <- .dp_target_interval_verification_v1(
        x, implied, achieved_interval, residuals, tolerances, controls
      )
    } else if (target_kind %in% c("moments", "cv", "pmf") &&
               !is.null(implied)) {
      verification <- .dp_target_simple_verification_v1(
        x, target_kind, implied, residuals, tolerances,
        request, normalized, used
      )
    } else {
      check_source <- "independent_target_failure_audit"
      verification <- .dpprior_new_verification(
        method = check_source,
        performed = TRUE, passed = FALSE, reason = x$message,
        settings = list(
          method = method,
          feasibility_source = "analytic_interval_group_support"
        ),
        selected_snapshot = NULL, verifier_snapshot = NULL,
        stability = NULL,
        components = list(
          solver_postcondition = .dpprior_new_check(
            value = FALSE, reference = TRUE, tolerance = NULL,
            operator = "identical", source = check_source
          )
        ),
        invariants = list(
          request_identity = .dpprior_new_check(
            value = identical(request$J, J), reference = TRUE,
            tolerance = NULL, operator = "identical", source = check_source
          ),
          support_identity = .dpprior_new_check(
            value = identical(as.integer(x$support), seq_len(J)),
            reference = TRUE, tolerance = NULL, operator = "identical",
            source = check_source
          )
        )
      )
    }
  }
  achieved_interval <- if (identical(target_kind, "interval") &&
      !is.null(x$pmf)) {
    as.list(.dp_target_interval_masses_v1(x$pmf, interval))
  } else {
    NULL
  }

  legacy_view <- legacy_x
  legacy_view$.canonical_controls <- NULL
  legacy_view$condition <- .dp_target_condition_view(legacy_x$condition)
  legacy_view <- .dp_target_plain_view_value(legacy_view)
  compatibility <- .dpprior_new_compatibility(
    top_level_aliases = stats::setNames(character(), character()),
    views = list(target_v0 = legacy_view),
    deprecations = list(
      target_v0 = paste(
        "Use canonical target fields; the plain migration view retains",
        "legacy scientific values and text representations of non-plain evidence."
      )
    )
  )
  source_commit <- getOption("DPprior.source_commit", NULL)
  if (!is.character(source_commit) || length(source_commit) != 1L ||
      is.na(source_commit) || !nzchar(source_commit)) {
    source_commit <- NULL
  }
  provenance <- .dpprior_new_provenance(
    requested_method = method,
    selected_method = method,
    is_fallback = FALSE,
    approximation = list(
      active = identical(x$status, "approximate"),
      opt_in = identical(x$status, "approximate") && isTRUE(x$usable),
      kind = if (identical(x$status, "approximate")) {
        "target_postcondition_not_met"
      } else {
        NULL
      },
      warning_code = if (identical(x$status, "approximate")) {
        "target_approximate"
      } else {
        NULL
      }
    ),
    projection = list(
      applied = FALSE, opt_in = FALSE, policy = NULL, record = NULL
    ),
    parameterization = "bounded_discrete_K_target",
    backend = list(
      package = "DPprior",
      package_version = tryCatch(
        as.character(utils::packageVersion("DPprior")),
        error = function(condition) "development"
      ),
      implementation = paste0(
        "R/19_elicitation_targets.R:.dp_K_target_class/", method
      ),
      source_commit = source_commit
    ),
    migration = list(
      source_schema = "phase7_target_v0",
      adapter = ".dp_K_target_class",
      lossless = FALSE,
      missing_evidence = c(
        "legacy_condition_object_identity_not_retained",
        "legacy_nonfinite_or_nonplain_evidence_encoded_as_plain_text",
        "legacy_wide_request_relocated_to_compatibility_view",
        "legacy_distribution_diagnostics_relocated_to_compatibility_view",
        if (is.null(source_commit)) "source_commit_not_embedded" else character()
      ),
      warnings = "canonical_target_is_producer_native_not_a_v0_object_clone"
    ),
    legacy = list(active = FALSE, contract = NULL, deprecation_stage = NULL)
  )

  target <- .dpprior_new_target_K(
    kind = target_kind,
    J = J,
    support = as.integer(x$support),
    request = request,
    normalized = normalized,
    used = used,
    derivation = derivation,
    interval = interval,
    family = family,
    assumptions = assumptions,
    pmf = if (is.null(x$pmf)) NULL else unname(x$pmf),
    implied = implied,
    achieved_interval = achieved_interval,
    residuals = residuals,
    tolerances = tolerances,
    status = x$status,
    usable = x$usable,
    verified = x$verified,
    message = x$message,
    parameters = NULL,
    computation = computation,
    verification = verification,
    provenance = provenance,
    compatibility = compatibility
  )
  .dpprior_append_compatibility_v2(
    target,
    aliases = c(
      attempts = "computation.attempts",
      constraint_residuals = "residuals"
    ),
    views = compatibility$views,
    deprecations = compatibility$deprecations
  )
}

.dp_K_target_moment <- function(J, mu_K, var_K, source, request,
                                derivation = NULL) {
  moments <- .a2_kl_validate_target_moments(J, mu_K, var_K)
  implied <- list(
    mean = moments$mu_K,
    variance = moments$var_K,
    sd = sqrt(moments$var_K),
    cv = sqrt(moments$var_K) / moments$mu_K,
    quantiles = NULL
  )
  .dp_K_target_class(list(
    kind = "K_target",
    J = J,
    support = seq_len(J),
    request = request,
    interval = NULL,
    family = list(
      name = NULL,
      parameterization = "moments_only_no_distribution_identified",
      explicit = FALSE
    ),
    assumptions = list(
      estimand = "finite-design occupied-cluster count K_J",
      pmf_identified = FALSE
    ),
    pmf = NULL,
    implied = implied,
    achieved_interval = NULL,
    constraint_residuals = list(
      raw = c(mean = 0, variance = 0),
      scaled = c(mean = 0, variance = 0),
      scale = c(mean = max(1, abs(moments$mu_K)),
                variance = max(1, abs(moments$var_K)))
    ),
    tolerances = list(moment = 0),
    verification = list(
      method = "independent finite-support moment feasibility bound",
      performed = TRUE,
      passed = TRUE,
      support_variance_upper = .max_var_K_fixed_mean(J, moments$mu_K)
    ),
    status = "converged",
    usable = TRUE,
    verified = TRUE,
    message = "Exact moment target validated; no target PMF was inferred.",
    parameters = list(mu_K = moments$mu_K, var_K = moments$var_K),
    attempts = list(),
    provenance = c(list(
      source = source,
      selected_method = "exact_moment_adapter",
      pmf_identified = FALSE,
      hidden_family_used = FALSE,
      normalization = "not_applicable"
    ), derivation)
  ))
}

.dp_K_target_custom_pmf <- function(J, target_pmf, mu_K, var_K, request,
                                    tolerance) {
  target_pmf <- .dpprior_validate_plain_vector(
    target_pmf, "target_pmf", "dpprior_pmf_error"
  )
  input_length <- length(target_pmf)
  if (is.numeric(target_pmf) && input_length == J + 1L) {
    if (!identical(as.numeric(target_pmf[1L]), 0)) {
      .dp_interval_abort(
        "target_pmf[1] represents K=0 and must be exactly zero.",
        c("dpprior_pmf_support_error", "dpprior_pmf_error"),
        "target_pmf[1]", target_pmf[1L], "exactly zero", "support"
      )
    }
    pmf <- target_pmf[-1L]
    k0_entry <- TRUE
  } else {
    pmf <- target_pmf
    k0_entry <- FALSE
  }
  pmf <- .dpprior_validate_pmf(
    pmf, "target_pmf", expected_length = J,
    require_sum_one = TRUE, tolerance = .TOL_PMF_SUM
  )
  pmf <- unname(as.numeric(pmf))
  implied <- .dp_K_pmf_metrics(pmf)

  check_moment <- function(supplied, achieved, name) {
    if (is.null(supplied)) return(invisible(TRUE))
    supplied <- .dpprior_validate_scalar(
      supplied, name, .subclass = "dpprior_pmf_moment_error"
    )
    scale <- max(1, abs(supplied), abs(achieved))
    if (abs(supplied - achieved) > tolerance * scale) {
      .dp_interval_abort(
        sprintf(
          "%s conflicts with the strict target_pmf (got %.12g; PMF implies %.12g).",
          name, supplied, achieved
        ),
        c("dpprior_pmf_moment_conflict", "dpprior_conflicting_input"),
        name, supplied, sprintf("%.12g within tolerance %.3g", achieved,
                                tolerance),
        "target_pmf_moment_conflict",
        implied = achieved,
        tolerance = tolerance
      )
    }
    invisible(TRUE)
  }
  check_moment(mu_K, implied$mean, "mu_K")
  check_moment(var_K, implied$variance, "var_K")

  verify <- .dpprior_validate_pmf(
    pmf, "verification_pmf", expected_length = J,
    require_sum_one = TRUE, tolerance = .TOL_PMF_SUM
  )
  verification_implied <- .dp_K_pmf_metrics_verify(verify)
  verification_passed <- isTRUE(all.equal(
    implied, verification_implied, tolerance = tolerance
  ))
  .dp_K_target_class(list(
    kind = "K_target",
    J = J,
    support = seq_len(J),
    request = request,
    interval = NULL,
    family = list(
      name = "custom_pmf",
      parameterization = "probability vector on inclusive support 1:J",
      explicit = TRUE
    ),
    assumptions = list(
      estimand = "finite-design occupied-cluster count K_J",
      pmf_authoritative = TRUE
    ),
    pmf = pmf,
    implied = implied,
    achieved_interval = NULL,
    constraint_residuals = list(
      raw = c(pmf_sum = sum(pmf) - 1),
      scaled = c(pmf_sum = sum(pmf) - 1),
      scale = c(pmf_sum = 1)
    ),
    tolerances = list(pmf_sum = .TOL_PMF_SUM, moment_check = tolerance),
    verification = list(
      method = paste0(
        "strict second-pass PMF validation plus separately coded direct ",
        "support-sum, centered-moment, and quantile postcondition audit"
      ),
      performed = TRUE,
      passed = verification_passed,
      pmf_sum = sum(verify),
      implied = verification_implied
    ),
    status = if (verification_passed) "converged" else "approximate",
    usable = verification_passed,
    verified = verification_passed,
    message = if (verification_passed) {
      "Strict custom target PMF validated without normalization."
    } else {
      "Strict target PMF passed input validation but not the independent postcondition audit."
    },
    parameters = NULL,
    attempts = list(),
    provenance = list(
      source = "target_pmf",
      selected_method = "strict_pmf_adapter",
      pmf_identified = TRUE,
      normalization = "validated_not_modified",
      input_length = input_length,
      input_sum = sum(target_pmf),
      k0_entry_removed = k0_entry
    )
  ))
}

# Canonical internal target constructor. Moment-only routes intentionally leave
# pmf=NULL because moments do not identify a distribution on 1:J.
.dp_target_K <- function(J, mu_K = NULL, var_K = NULL, confidence = NULL,
                         cv_K = NULL, K_interval = NULL, target_pmf = NULL,
                         tolerance = 1e-9, root_control = NULL) {
  assert_valid_J(J)
  J <- as.integer(J)
  tolerance <- .dpprior_validate_control(
    tolerance, "tolerance", type = "numeric",
    lower = .DP_ELICITATION_TOL_MIN
  )
  if (!is.null(root_control) && !is.list(root_control)) {
    .dpprior_abort_invalid(
      "root_control must be NULL or a named list.",
      c("dpprior_elicitation_control_error", "dpprior_type_error"),
      "root_control", root_control, "NULL or named list", "type"
    )
  }
  if (is.null(root_control)) root_control <- list()
  if (length(root_control) &&
      (is.null(names(root_control)) || any(!nzchar(names(root_control))) ||
       anyDuplicated(names(root_control)))) {
    .dpprior_abort_invalid(
      "root_control must be a uniquely named list.",
      c("dpprior_elicitation_control_error", "dpprior_type_error"),
      "root_control", root_control, "uniquely named list", "names"
    )
  }
  allowed_controls <- c("root_tol", "max_iterations")
  unknown_controls <- setdiff(names(root_control), allowed_controls)
  if (length(unknown_controls)) {
    .dpprior_abort_invalid(
      sprintf("Unknown root_control field(s): %s.",
              paste(unknown_controls, collapse = ", ")),
      c("dpprior_elicitation_control_error",
        "dpprior_unknown_argument_error"),
      "root_control", root_control, paste(allowed_controls, collapse = ", "),
      "unknown_fields"
    )
  }
  root_control <- list(
    root_tol = if (is.null(root_control$root_tol)) {
      max(.Machine$double.eps, min(1e-12, tolerance / 10))
    } else {
      .dpprior_validate_control(
        root_control$root_tol, "root_control$root_tol", "numeric",
        lower = .Machine$double.eps,
        upper = if (is.null(K_interval)) NULL else max(
          .Machine$double.eps, min(1e-12, tolerance / 10)
        )
      )
    },
    max_iterations = if (is.null(root_control$max_iterations)) {
      1000L
    } else {
      .dpprior_validate_control(
        root_control$max_iterations, "root_control$max_iterations", "count",
        lower = 1L,
        upper = if (is.null(K_interval)) .Machine$integer.max else 1000L
      )
    }
  )

  request <- list(
    J = J, mu_K = mu_K, var_K = var_K, confidence = confidence,
    cv_K = cv_K, K_interval = K_interval, target_pmf = target_pmf
  )

  # A direct variance accompanying target_pmf is a backward-compatible
  # consistency annotation, not a second source of uncertainty.
  source_flags <- c(
    var_K = !is.null(var_K) && is.null(target_pmf),
    confidence = !is.null(confidence),
    cv_K = !is.null(cv_K),
    K_interval = !is.null(K_interval),
    target_pmf = !is.null(target_pmf)
  )
  if (sum(source_flags) != 1L) {
    .dpprior_abort_invalid(
      paste0(
        "Specify exactly one uncertainty source: var_K, confidence, cv_K, ",
        "K_interval, or target_pmf."
      ),
      c("dpprior_elicitation_source_error", "dpprior_conflicting_input"),
      "uncertainty source", names(source_flags)[source_flags],
      "exactly one source", if (sum(source_flags) == 0L) "missing" else "conflict"
    )
  }

  source <- names(source_flags)[source_flags]
  if (identical(source, "target_pmf")) {
    return(.dp_K_target_custom_pmf(
      J, target_pmf, mu_K, var_K, request, tolerance
    ))
  }
  if (identical(source, "K_interval")) {
    if (tolerance > 1e-9) {
      .dpprior_abort_invalid(
        paste0(
          "tolerance must be <= 1e-9 for a verified interval target; ",
          "this is the canonical truth-tolerance ceiling."
        ),
        c("dpprior_elicitation_control_error", "dpprior_bounds_error"),
        "tolerance", tolerance,
        sprintf("in [%.17g, 1e-9]", .DP_ELICITATION_TOL_MIN),
        "interval_tolerance_ceiling"
      )
    }
    return(.dp_target_K_interval(
      J, K_interval, mu_K, request, tolerance, root_control
    ))
  }

  if (is.null(mu_K)) {
    .dpprior_abort_invalid(
      sprintf("mu_K is required when uncertainty is supplied through %s.", source),
      c("dpprior_target_moment_error", "dpprior_elicitation_source_error"),
      "mu_K", NULL, "one scalar in (1, J)", "missing_mean"
    )
  }

  if (identical(source, "var_K")) {
    return(.dp_K_target_moment(
      J, mu_K, var_K, "direct_variance", request,
      list(raw_var_K = var_K, derived_variance = var_K)
    ))
  }
  if (identical(source, "confidence")) {
    vif <- confidence_to_vif_fit(confidence)
    derived <- vif * (mu_K - 1)
    return(.dp_K_target_moment(
      J, mu_K, derived, "qualitative_confidence", request,
      list(confidence = confidence, vif = vif, derived_variance = derived)
    ))
  }

  cv_K <- .dpprior_validate_scalar(
    cv_K, "cv_K", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_cv_K_error"
  )
  mu_valid <- .dpprior_validate_scalar(
    mu_K, "mu_K", lower = 1, upper = J,
    lower_open = TRUE, upper_open = TRUE,
    .subclass = "dpprior_target_moment_error"
  )
  derived <- (cv_K * mu_valid)^2
  .dp_K_target_moment(
    J, mu_valid, derived, "coefficient_of_variation", request,
    list(
      raw_cv_K = cv_K,
      definition = "SD(K_J) / E(K_J)",
      mapping = "Var(K_J) = (cv_K * mu_K)^2",
      derived_variance = derived
    )
  )
}

#' Construct a Canonical Bounded-Discrete Cluster-Count Target
#'
#' Constructs and validates a canonical elicitation target for the occupied
#' cluster count \eqn{K_J} on the exact support \eqn{1,\ldots,J}. Supply exactly
#' one uncertainty source: \code{var_K}, \code{confidence}, \code{cv_K},
#' \code{K_interval}, or \code{target_pmf}.
#'
#' @param J Integer design size and upper support point for \eqn{K_J}.
#' @param mu_K Optional target mean for \eqn{K_J}.
#' @param var_K Optional target variance for \eqn{K_J}.
#' @param confidence Optional qualitative confidence label.
#' @param cv_K Optional coefficient of variation, defined as
#'   \eqn{SD(K_J)/E(K_J)}.
#' @param K_interval Optional explicit interval-target specification.
#' @param target_pmf Optional strict PMF on the support \eqn{1,\ldots,J}; a
#'   structural zero for \eqn{K=0} may be supplied as the first of \code{J+1}
#'   entries.
#' @param tolerance Positive numerical postcondition tolerance.
#' @param root_control Optional named controls for interval-target root solving.
#'
#' @return A validated \code{dpprior.target/1} \code{dpprior_K_target}
#'   object. The exact support and original request are retained alongside
#'   normalized and used values, derivation evidence, implied moments,
#'   tolerances, verification, and provenance. Pass the object unchanged as
#'   \code{target_K} to \code{\link{DPprior_fit}}.
#'
#' @examples
#' target <- DPprior_target_K(J = 50, mu_K = 5, var_K = 8)
#' target$used[c("mu_K", "var_K")]
#' fit <- DPprior_fit(
#'   J = 50, target_K = target, check_diagnostics = FALSE
#' )
#'
#' @family elicitation
#' @export
DPprior_target_K <- function(J, mu_K = NULL, var_K = NULL,
                             confidence = NULL, cv_K = NULL,
                             K_interval = NULL, target_pmf = NULL,
                             tolerance = 1e-9, root_control = NULL) {
  .dp_target_K(
    J = J, mu_K = mu_K, var_K = var_K, confidence = confidence,
    cv_K = cv_K, K_interval = K_interval, target_pmf = target_pmf,
    tolerance = tolerance, root_control = root_control
  )
}

.dp_interval_groups <- function(interval, J) {
  support <- seq_len(J)
  inside <- support >= interval$lower & support <= interval$upper
  if (identical(interval$type, "hard_bounds")) {
    return(list(
      values = list(inside = support[inside]),
      masses = c(inside = 1)
    ))
  }
  if (identical(interval$type, "central_mass")) {
    return(list(
      values = list(inside = support[inside], outside = support[!inside]),
      masses = c(inside = interval$coverage, outside = 1 - interval$coverage)
    ))
  }
  tail_mass <- (1 - interval$coverage) / 2
  list(
    values = list(
      left = support[support < interval$lower],
      inside = support[inside],
      right = support[support > interval$upper]
    ),
    masses = c(left = tail_mass, inside = interval$coverage,
               right = tail_mass)
  )
}

.dp_interval_group_feasibility <- function(groups) {
  empty <- names(groups$values)[lengths(groups$values) == 0L & groups$masses > 0]
  if (length(empty)) {
    return(list(
      feasible = FALSE,
      reason = "positive_mass_group_has_empty_support",
      empty_groups = empty,
      minimum = NA_real_,
      maximum = NA_real_
    ))
  }
  list(
    feasible = TRUE,
    reason = "analytic_group_mass_convex_hull",
    empty_groups = character(),
    minimum = sum(vapply(
      names(groups$values),
      function(name) groups$masses[[name]] * min(groups$values[[name]]),
      numeric(1L)
    )),
    maximum = sum(vapply(
      names(groups$values),
      function(name) groups$masses[[name]] * max(groups$values[[name]]),
      numeric(1L)
    ))
  )
}

.dp_interval_group_pmf <- function(J, groups, tilt = 0,
                                   boundary = NULL) {
  pmf <- numeric(J)
  for (name in names(groups$values)) {
    values <- groups$values[[name]]
    mass <- groups$masses[[name]]
    if (mass == 0) next
    if (!is.null(boundary)) {
      point <- if (identical(boundary, "minimum")) min(values) else max(values)
      pmf[point] <- pmf[point] + mass
    } else {
      log_weights <- tilt * values
      probabilities <- exp(log_weights - logsumexp_vec(log_weights))
      pmf[values] <- pmf[values] + mass * probabilities
    }
  }
  unname(pmf)
}

# Verification implementation deliberately does not call
# .dp_interval_group_pmf() or logsumexp_vec(). It reconstructs the same
# finite-support exponential-family PMF with an independently coded max-shift
# normalization, then the caller audits the requested group masses and mean.
.dp_interval_group_pmf_verify <- function(J, groups, tilt = 0,
                                          boundary = NULL) {
  pmf <- numeric(J)
  for (name in names(groups$values)) {
    values <- groups$values[[name]]
    mass <- groups$masses[[name]]
    if (mass == 0) next
    if (!is.null(boundary)) {
      point <- if (identical(boundary, "minimum")) min(values) else max(values)
      pmf[point] <- pmf[point] + mass
      next
    }
    log_weights <- tilt * values
    shift <- max(log_weights)
    weights <- exp(log_weights - shift)
    probabilities <- weights / sum(weights)
    pmf[values] <- pmf[values] + mass * probabilities
  }
  unname(pmf)
}

.dp_interval_metrics <- function(pmf, interval) {
  support <- seq_along(pmf)
  left <- sum(pmf[support < interval$lower])
  inside <- sum(pmf[
    support >= interval$lower & support <= interval$upper
  ])
  right <- sum(pmf[support > interval$upper])
  list(left_mass = left, inside_mass = inside, right_mass = right)
}

.dp_interval_requested_masses <- function(interval) {
  if (identical(interval$type, "hard_bounds")) {
    return(c(left_mass = 0, inside_mass = 1, right_mass = 0))
  }
  if (identical(interval$type, "central_mass")) {
    return(c(left_mass = NA_real_, inside_mass = interval$coverage,
             right_mass = NA_real_))
  }
  tail <- (1 - interval$coverage) / 2
  c(left_mass = tail, inside_mass = interval$coverage, right_mass = tail)
}

.dp_K_target_infeasible <- function(J, request, interval, family, message,
                                    code, condition_class, feasibility,
                                    tolerance, root_control) {
  condition <- .dpprior_new_condition(
    message = message,
    classes = c(condition_class, "dpprior_interval_infeasible",
                "dpprior_calibration_error", "dpprior_error", "error"),
    code = code,
    interval = interval,
    feasibility = feasibility
  )
  .dp_K_target_class(list(
    kind = "K_target",
    J = J,
    support = seq_len(J),
    request = request,
    interval = interval,
    family = family,
    assumptions = list(
      estimand = "finite-design occupied-cluster count K_J",
      endpoints = "inclusive",
      reference_measure = "counting measure on 1:J"
    ),
    pmf = NULL,
    implied = NULL,
    achieved_interval = NULL,
    constraint_residuals = NULL,
    tolerances = list(constraint = tolerance),
    verification = list(
      method = "analytic finite-support group-mass feasibility",
      performed = TRUE,
      passed = TRUE,
      constraint_feasible = FALSE,
      infeasibility_certified = TRUE,
      feasibility = feasibility
    ),
    status = "infeasible",
    usable = FALSE,
    verified = TRUE,
    message = message,
    parameters = NULL,
    attempts = list(),
    .canonical_controls = root_control,
    condition = condition,
    provenance = list(
      source = "K_interval",
      selected_method = "analytic_feasibility",
      feasibility_unknown = FALSE,
      fallback_used = FALSE,
      projection_used = FALSE
    )
  ))
}

.dp_K_target_failed <- function(J, request, interval, family, message,
                                code, source_condition, feasibility,
                                tolerance, root_control,
                                attempts = list()) {
  condition <- .dpprior_new_condition(
    message = message,
    classes = c("dpprior_interval_solver_failed", "dpprior_numerical_error",
                "dpprior_calibration_error", "dpprior_error", "error"),
    code = code,
    interval = interval,
    feasibility = feasibility,
    source_condition = source_condition
  )
  .dp_K_target_class(list(
    kind = "K_target",
    J = J,
    support = seq_len(J),
    request = request,
    interval = interval,
    family = family,
    assumptions = list(
      estimand = "finite-design occupied-cluster count K_J",
      endpoints = "inclusive",
      reference_measure = "counting measure on 1:J"
    ),
    pmf = NULL,
    implied = NULL,
    achieved_interval = NULL,
    constraint_residuals = NULL,
    tolerances = list(constraint = tolerance),
    verification = list(
      method = "analytic feasibility plus numerical-solver postcondition audit",
      performed = TRUE,
      passed = FALSE,
      constraint_feasible = isTRUE(feasibility$feasible),
      infeasibility_certified = FALSE,
      numerical_solver_passed = FALSE,
      feasibility = feasibility
    ),
    status = "failed",
    usable = FALSE,
    verified = FALSE,
    message = message,
    parameters = NULL,
    attempts = attempts,
    .canonical_controls = root_control,
    condition = condition,
    provenance = list(
      source = "K_interval",
      selected_method = "numerical_solver_failure",
      feasibility_unknown = FALSE,
      fallback_used = FALSE,
      projection_used = FALSE
    )
  ))
}

.dp_target_K_unusable_condition <- function(target) {
  status <- target[["status", exact = TRUE]]
  message <- target[["message", exact = TRUE]]
  interval <- target[["interval", exact = TRUE]]
  verification <- target[["verification", exact = TRUE]]
  computation <- target[["computation", exact = TRUE]]

  if (identical(status, "infeasible")) {
    certificate <- verification[["settings", exact = TRUE]][[
      "certificate", exact = TRUE
    ]]
    certificate_kind <- certificate[["kind", exact = TRUE]]
    if (identical(certificate_kind,
                  "positive_mass_group_has_empty_support")) {
      primary_class <- "dpprior_interval_empty_tail_error"
      code <- "empty_interval_group"
    } else if (identical(certificate_kind,
                         "mean_outside_group_mass_hull")) {
      primary_class <- "dpprior_interval_mean_infeasible"
      code <- "mean_outside_group_mass_hull"
    } else {
      primary_class <- "dpprior_interval_infeasible"
      code <- "target_certified_infeasible"
    }
    return(.dpprior_new_condition(
      message = message,
      classes = unique(c(
        primary_class, "dpprior_interval_infeasible",
        "dpprior_calibration_error", "dpprior_error", "error"
      )),
      code = code,
      interval = interval,
      certificate = certificate,
      result = target
    ))
  }

  if (identical(status, "failed")) {
    attempts <- computation[["attempts", exact = TRUE]]
    methods <- if (length(attempts)) {
      vapply(
        attempts,
        function(attempt) attempt[["method", exact = TRUE]],
        character(1)
      )
    } else {
      character()
    }
    code <- if ("expanding_root_bracket" %in% methods) {
      "maxent_root_not_bracketed"
    } else if ("uniroot_common_exponential_tilt" %in% methods) {
      "maxent_root_solver_error"
    } else {
      "target_failed"
    }
    return(.dpprior_new_condition(
      message = message,
      classes = c(
        "dpprior_interval_solver_failed", "dpprior_numerical_error",
        "dpprior_calibration_error", "dpprior_error", "error"
      ),
      code = code,
      interval = interval,
      result = target
    ))
  }

  .dpprior_new_condition(
    message = message,
    classes = c(
      "dpprior_elicitation_unusable", "dpprior_calibration_error",
      "dpprior_error", "error"
    ),
    code = paste0("target_", status),
    result = target
  )
}

.dp_target_K_stop_unusable <- function(target) {
  if (!inherits(target, "dpprior_K_target")) {
    .dpprior_abort_invalid(
      "target must inherit from dpprior_K_target.",
      c("dpprior_target_integrity_error", "dpprior_type_error"),
      "target", target, "dpprior_K_target", "type"
    )
  }
  .dpprior_validate_target_v1(target)
  if (isTRUE(target$usable)) return(target)
  stop(.dp_target_K_unusable_condition(target))
}

.dp_target_from_interval_maxent <- function(J, interval, mu_K, request,
                                            tolerance, root_control) {
  groups <- .dp_interval_groups(interval, J)
  feasibility <- .dp_interval_group_feasibility(groups)
  family <- list(
    name = "maxent",
    parameterization = "common exponential tilt within fixed-mass groups",
    explicit = TRUE,
    reference_measure = "counting measure on 1:J"
  )
  if (!isTRUE(feasibility$feasible)) {
    return(.dp_K_target_infeasible(
      J, request, interval, family,
      sprintf(
        "The %s request is infeasible because group(s) %s have empty support.",
        interval$type, paste(feasibility$empty_groups, collapse = ", ")
      ),
      "empty_interval_group", "dpprior_interval_empty_tail_error",
      feasibility, tolerance, root_control
    ))
  }

  needs_mean <- !is.null(mu_K)
  if (identical(interval$type, "central_mass") && !needs_mean) {
    .dp_interval_abort(
      paste0(
        "type='central_mass' requires an actual discrete mu_K location ",
        "anchor; coverage alone does not identify the target."
      ),
      c("dpprior_interval_missing_mean_error",
        "dpprior_interval_under_specified"),
      "mu_K", NULL, sprintf("one scalar in [%d, %d]", interval$lower,
                             interval$upper),
      "missing_mean"
    )
  }

  feasibility_epsilon <- 64 * .Machine$double.eps * max(
    1, abs(feasibility$minimum), abs(feasibility$maximum),
    if (needs_mean) abs(mu_K) else 1
  )
  if (needs_mean &&
      (mu_K < feasibility$minimum - feasibility_epsilon ||
       mu_K > feasibility$maximum + feasibility_epsilon)) {
    return(.dp_K_target_infeasible(
      J, request, interval, family,
      sprintf(
        paste0(
          "Requested mu_K=%.12g is outside the analytic feasible mean ",
          "interval [%.12g, %.12g] for this interval constraint."
        ),
        mu_K, feasibility$minimum, feasibility$maximum
      ),
      "mean_outside_group_mass_hull", "dpprior_interval_mean_infeasible",
      feasibility, tolerance, root_control
    ))
  }

  boundary <- NULL
  tilt <- 0
  solver_passed <- TRUE
  attempt <- list(
    method = "analytic_uniform_within_groups",
    start = 0,
    bounds = c(-Inf, Inf),
    exit_code = 0L,
    message = "No root required.",
    iterations = 0L
  )
  if (!needs_mean && identical(interval$type, "hard_bounds") &&
      interval$lower == interval$upper) {
    boundary <- "minimum"
    attempt$method <- "analytic_singleton_boundary"
    attempt$message <- "Singleton hard support fixes a sparse boundary PMF."
  }
  if (needs_mean) {
    boundary_tol <- feasibility_epsilon
    if (abs(mu_K - feasibility$minimum) <= boundary_tol) {
      boundary <- "minimum"
    } else if (abs(mu_K - feasibility$maximum) <= boundary_tol) {
      boundary <- "maximum"
    } else {
      mean_at <- function(lambda) {
        pmf <- .dp_interval_group_pmf(J, groups, tilt = lambda)
        sum(seq_len(J) * pmf)
      }
      objective <- function(lambda) mean_at(lambda) - mu_K
      span <- 1
      lower_value <- objective(-span)
      upper_value <- objective(span)
      expansions <- 0L
      while ((lower_value >= 0 || upper_value <= 0) && expansions < 60L) {
        span <- span * 2
        lower_value <- objective(-span)
        upper_value <- objective(span)
        expansions <- expansions + 1L
      }
      if (lower_value >= 0 || upper_value <= 0) {
        failed_feasibility <- c(
          feasibility,
          list(bracket = c(-span, span),
               values = c(lower_value, upper_value))
        )
        return(.dp_K_target_failed(
          J, request, interval, family,
          "The maximum-entropy tilt could not bracket the requested interior mean.",
          "maxent_root_not_bracketed", NULL, failed_feasibility, tolerance,
          root_control,
          attempts = list(list(
            method = "expanding_root_bracket",
            bounds = c(-span, span),
            values = c(lower_value, upper_value),
            expansions = expansions,
            exit_code = 1L,
            message = "Analytically feasible request was not numerically bracketed."
          ))
        ))
      }
      root_tol <- root_control$root_tol
      root_warnings <- character()
      root <- tryCatch(
        withCallingHandlers(
          stats::uniroot(
            objective, interval = c(-span, span), tol = root_tol,
            maxiter = root_control$max_iterations
          ),
          warning = function(condition) {
            root_warnings <<- c(root_warnings, conditionMessage(condition))
            invokeRestart("muffleWarning")
          }
        ),
        error = function(condition) condition
      )
      if (inherits(root, "error")) {
        return(.dp_K_target_failed(
          J, request, interval, family,
          paste("The maximum-entropy root solver failed:",
                conditionMessage(root)),
          "maxent_root_solver_error", root, feasibility, tolerance,
          root_control,
          attempts = list(list(
            method = "uniroot_common_exponential_tilt",
            bounds = c(-span, span),
            exit_code = 1L,
            message = conditionMessage(root),
            warnings = root_warnings
          ))
        ))
      }
      tilt <- root$root
      root_residual <- objective(tilt)
      root_passed <- is.finite(root_residual) &&
        abs(root_residual) <= tolerance * max(1, abs(mu_K)) &&
        !length(root_warnings)
      solver_passed <- root_passed
      attempt <- list(
        method = "uniroot_common_exponential_tilt",
        start = 0,
        bounds = c(-span, span),
        exit_code = if (root_passed) 0L else 1L,
        message = if (root_passed) {
          "Root bracketed and solved; postcondition passed."
        } else {
          "Root returned a candidate but did not pass the convergence postcondition."
        },
        iterations = root$iter,
        estim_prec = root$estim.prec,
        residual = root_residual,
        warnings = root_warnings,
        converged = root_passed
      )
    }
  }

  pmf <- .dp_interval_group_pmf(
    J, groups, tilt = tilt, boundary = boundary
  )
  pmf <- .dpprior_validate_pmf(
    pmf, "maxent_target_pmf", expected_length = J,
    tolerance = .TOL_PMF_SUM
  )
  verification_pmf <- .dp_interval_group_pmf_verify(
    J, groups, tilt = tilt, boundary = boundary
  )
  verification_pmf <- .dpprior_validate_pmf(
    verification_pmf, "maxent_verification_pmf", expected_length = J,
    tolerance = .TOL_PMF_SUM
  )
  backcheck <- .dp_backcheck_K_interval(
    pmf, interval, verification_pmf = verification_pmf,
    J = J, tolerance = tolerance, source = "target_constructor"
  )
  implied <- .dp_K_pmf_metrics(pmf)
  verification_implied <- .dp_K_pmf_metrics_verify(verification_pmf)
  mean_residual <- if (needs_mean) implied$mean - mu_K else NA_real_
  mean_verified_residual <- if (needs_mean) {
    verification_implied$mean - mu_K
  } else {
    NA_real_
  }
  interval_passed <- identical(backcheck$status, "converged")
  mean_passed <- !needs_mean || abs(mean_verified_residual) <=
    tolerance * max(1, abs(mu_K))
  requested_masses <- .dp_interval_requested_masses(interval)
  selected_achieved <- unlist(
    .dp_interval_metrics(pmf, interval), use.names = TRUE
  )
  verification_achieved <- unlist(
    .dp_interval_metrics(verification_pmf, interval), use.names = TRUE
  )
  positive_mass_names <- names(requested_masses)[
    !is.na(requested_masses) & requested_masses > 0
  ]
  positive_mass_representable <- all(
    requested_masses[positive_mass_names] >= .Machine$double.xmin
  )
  structural_mass_passed <- positive_mass_representable &&
    all(selected_achieved[positive_mass_names] > 0) &&
    all(verification_achieved[positive_mass_names] > 0)
  verified <- interval_passed && mean_passed && solver_passed &&
    structural_mass_passed
  status <- if (!is.null(boundary)) {
    if (verified) "boundary" else "approximate"
  } else if (verified) {
    "converged"
  } else {
    "approximate"
  }
  achieved <- selected_achieved
  interval_residual <- achieved - requested_masses
  raw <- c(mean = mean_residual, interval_residual)
  scale <- c(
    mean = if (needs_mean) max(1, abs(mu_K)) else NA_real_,
    left_mass = 1, inside_mass = 1, right_mass = 1
  )

  .dp_K_target_class(list(
    kind = "K_target",
    J = J,
    support = seq_len(J),
    request = request,
    interval = interval,
    family = family,
    assumptions = list(
      estimand = "finite-design occupied-cluster count K_J",
      endpoints = "inclusive",
      maximum_entropy_reference = "counting measure on 1:J",
      group_constraints = groups$masses,
      mean_constraint = needs_mean
    ),
    pmf = unname(pmf),
    implied = implied,
    achieved_interval = backcheck$selected$achieved,
    constraint_residuals = list(
      raw = raw,
      scaled = raw / scale,
      scale = scale
    ),
    tolerances = list(constraint = tolerance, pmf_sum = .TOL_PMF_SUM),
    verification = list(
      method = paste0(
        "separate max-shift PMF reconstruction plus direct support, ",
        "group-mass, and mean postcondition audit"
      ),
      performed = TRUE,
      passed = verified,
      numerical_solver_passed = solver_passed,
      positive_requested_groups_represented = structural_mass_passed,
      positive_requested_groups_representable = positive_mass_representable,
      requested_group_masses = requested_masses,
      selected_group_masses = selected_achieved,
      verification_group_masses = verification_achieved,
      pmf_max_abs_delta = max(abs(pmf - verification_pmf)),
      implied = verification_implied,
      mean_residual = mean_verified_residual,
      interval_backcheck = backcheck$verification
    ),
    status = status,
    usable = status %in% c("converged", "boundary"),
    verified = verified,
    message = if (identical(status, "boundary")) {
      "Verified maximum-entropy target on an analytic feasibility boundary."
    } else if (identical(status, "converged")) {
      "Verified maximum-entropy interval target."
    } else if (!solver_passed) {
      "Maximum-entropy candidate was retained, but its root solver did not pass the declared convergence postcondition."
    } else if (!structural_mass_passed) {
      "Maximum-entropy candidate lost a requested positive-mass group at floating-point resolution."
    } else {
      "Maximum-entropy candidate did not pass the declared constraints."
    },
    parameters = list(
      common_tilt = if (is.null(boundary)) tilt else NULL,
      boundary = boundary,
      group_masses = groups$masses
    ),
    attempts = list(attempt),
    .canonical_controls = root_control,
    provenance = list(
      source = "K_interval",
      requested_family = "maxent",
      selected_method = if (is.null(boundary)) {
        "analytic_groups_plus_common_tilt_root"
      } else {
        "analytic_sparse_boundary"
      },
      pmf_identified = TRUE,
      normalization = "analytic_group_masses_and_logsumexp",
      verification_normalization = "separate_max_shift_weight_sum",
      fallback_used = FALSE,
      projection_used = FALSE,
      feasibility = feasibility
    )
  ))
}

.dp_target_K_interval <- function(J, K_interval, mu_K, request, tolerance,
                                  root_control) {
  interval <- .dp_validate_K_interval(K_interval, J)
  interval_mu <- interval$mu_K
  if (!is.null(mu_K)) {
    mu_K <- .dpprior_validate_scalar(
      mu_K, "mu_K", lower = 1, upper = J,
      .subclass = "dpprior_interval_mean_error"
    )
  }
  if (!is.null(mu_K) && !is.null(interval_mu)) {
    scale <- max(1, abs(mu_K), abs(interval_mu))
    if (abs(mu_K - interval_mu) > tolerance * scale) {
      .dp_interval_abort(
        "Top-level mu_K conflicts with K_interval$mu_K.",
        c("dpprior_interval_mean_conflict", "dpprior_conflicting_input"),
        "mu_K", mu_K, sprintf("match K_interval$mu_K=%.12g", interval_mu),
        "mean_conflict",
        interval_mu_K = interval_mu,
        tolerance = tolerance
      )
    }
  }
  resolved_mu <- if (!is.null(mu_K)) mu_K else interval_mu
  interval$mu_K <- resolved_mu

  if (identical(interval$type, "central_mass") && !is.null(resolved_mu) &&
      (resolved_mu < interval$lower || resolved_mu > interval$upper)) {
    .dp_interval_abort(
      paste0(
        "For type='central_mass', mu_K is a location anchor and must lie ",
        "inside the inclusive interval."
      ),
      c("dpprior_interval_mean_error", "dpprior_bounds_error"),
      "mu_K", resolved_mu,
      sprintf("in [%d, %d]", interval$lower, interval$upper),
      "central_mean_outside_interval"
    )
  }

  if (!identical(interval$family, "maxent")) {
    .dp_interval_abort(
      sprintf(
        paste0(
          "K_interval family '%s' is deferred in the initial Phase 7 ",
          "constructor. Use family='maxent' or supply a strict target_pmf; ",
          "no fallback family was selected."
        ),
        interval$family
      ),
      c("dpprior_interval_family_deferred", "dpprior_interval_family_error"),
      "K_interval$family", interval$family,
      "maxent or strict target_pmf", "family_deferred",
      deferred_family = interval$family,
      fallback_selected = FALSE
    )
  }

  .dp_target_from_interval_maxent(
    J, interval, resolved_mu, request, tolerance, root_control
  )
}

.dp_backcheck_pmf <- function(pmf, J, name) {
  pmf <- .dpprior_validate_plain_vector(pmf, name, "dpprior_pmf_error")
  input_length <- length(pmf)
  if (is.numeric(pmf) && input_length == J + 1L) {
    if (!identical(as.numeric(pmf[1L]), 0)) {
      .dp_interval_abort(
        sprintf("%s has a nonzero K=0 entry.", name),
        c("dpprior_pmf_support_error", "dpprior_pmf_error"),
        paste0(name, "[1]"), pmf[1L], "exactly zero", "support"
      )
    }
    pmf <- pmf[-1L]
  }
  pmf <- .dpprior_validate_pmf(
    pmf, name, expected_length = J, require_sum_one = TRUE,
    tolerance = .TOL_PMF_SUM
  )
  unname(as.numeric(pmf))
}

.dp_backcheck_one <- function(pmf, interval, tolerance) {
  achieved <- unlist(.dp_interval_metrics(pmf, interval), use.names = TRUE)
  requested <- .dp_interval_requested_masses(interval)
  residual <- achieved - requested
  constrained <- is.finite(requested)
  component_passed <- rep(NA, length(requested))
  names(component_passed) <- names(requested)
  component_passed[constrained] <-
    abs(residual[constrained]) <= tolerance
  passed <- all(component_passed[constrained])
  scalar_list <- function(x) {
    stats::setNames(lapply(seq_along(x), function(index) unname(x[index])),
                    names(x))
  }
  list(
    requested = scalar_list(requested),
    achieved = scalar_list(achieved),
    residuals = scalar_list(residual),
    constrained = scalar_list(constrained),
    component_passed = scalar_list(component_passed),
    passed = passed
  )
}

# Recompute requested-versus-achieved interval behavior from a selected PMF
# and, when supplied, a separately evaluated verification PMF. The helper does
# not assume that moment matching implies interval matching.
.dp_backcheck_K_interval <- function(pmf, interval, verification_pmf = NULL,
                                     J = NULL, tolerance = 1e-9,
                                     source = "calibrated_prior") {
  if (is.null(J)) {
    if (!is.numeric(pmf)) {
      .dpprior_abort_invalid(
        "J must be supplied when pmf is not a numeric vector.",
        c("dpprior_interval_backcheck_error", "dpprior_type_error"),
        "J", J, "positive integer", "missing"
      )
    }
    # Without an explicit J, a leading structural zero is a valid probability
    # at K=1 and must not be guessed to be a legacy K=0 entry.
    J <- length(pmf)
  }
  assert_valid_J(J)
  J <- as.integer(J)
  tolerance <- .dpprior_validate_control(
    tolerance, "tolerance", "numeric", lower = 0, lower_open = TRUE
  )
  if (!is.character(source) || !is.null(dim(source)) || is.object(source) ||
      length(source) != 1L || is.na(source) || !nzchar(source)) {
    .dpprior_abort_invalid(
      "source must be one non-empty character value.",
      c("dpprior_interval_backcheck_error", "dpprior_type_error"),
      "source", source, "character scalar", "type"
    )
  }

  # Constructor-returned canonical intervals include coverage=1 for hard
  # bounds. Raw user requests are validated through the public input semantics,
  # where coverage is forbidden for hard bounds.
  canonical <- is.list(interval) &&
    all(c("lower", "upper", "type", "coverage", "family",
          "support", "endpoints") %in% names(interval))
  if (canonical) {
    canonical_names_ok <- !is.null(names(interval)) &&
      !anyNA(names(interval)) && all(nzchar(names(interval))) &&
      !anyDuplicated(names(interval))
    support_ok <- .dpprior_is_plain_numeric(interval$support) &&
      length(interval$support) == 2L && !anyNA(interval$support) &&
      all(is.finite(interval$support)) &&
      identical(
        unname(as.numeric(interval$support)), as.numeric(c(1L, J))
      )
    endpoints_ok <- identical(interval$endpoints, "inclusive")
    hard_coverage_ok <- !identical(interval$type, "hard_bounds") || (
      .dpprior_is_plain_numeric(interval$coverage) &&
        length(interval$coverage) == 1L &&
        !is.na(interval$coverage) && is.finite(interval$coverage) &&
        interval$coverage == 1
    )
    if (!canonical_names_ok || !support_ok || !endpoints_ok ||
        !hard_coverage_ok) {
      .dp_interval_abort(
        "Canonical interval is inconsistent with the requested support.",
        c("dpprior_interval_backcheck_error",
          "dpprior_interval_support_error"),
        value = interval,
        expected = sprintf("valid inclusive interval on 1:%d", J),
        code = "canonical_interval_invalid"
      )
    }
    raw_interval <- list(
      lower = interval$lower,
      upper = interval$upper,
      type = interval$type,
      family = interval$family
    )
    if (!identical(interval$type, "hard_bounds")) {
      raw_interval$coverage <- interval$coverage
    }
    if ("mu_K" %in% names(interval) && !is.null(interval$mu_K)) {
      raw_interval$mu_K <- interval$mu_K
    }
    interval <- .dp_validate_K_interval(raw_interval, J)
    if (identical(interval$type, "central_mass") &&
        !is.null(interval$mu_K) &&
        (interval$mu_K < interval$lower ||
         interval$mu_K > interval$upper)) {
      .dp_interval_abort(
        "Canonical central-mass interval has a location anchor outside its inclusive endpoints.",
        c("dpprior_interval_backcheck_error",
          "dpprior_interval_mean_error"),
        "interval$mu_K", interval$mu_K,
        sprintf("in [%d, %d]", interval$lower, interval$upper),
        "canonical_interval_mean"
      )
    }
  } else {
    interval <- .dp_validate_K_interval(interval, J)
  }

  selected_pmf <- .dp_backcheck_pmf(pmf, J, "pmf")
  selected <- .dp_backcheck_one(selected_pmf, interval, tolerance)
  if (is.null(verification_pmf)) {
    stability_delta <- c(
      left_mass = NA_real_, inside_mass = NA_real_, right_mass = NA_real_
    )
    verification <- list(
      performed = FALSE,
      reason = "independent_verification_not_supplied",
      method = NULL,
      requested = selected$requested,
      achieved = NULL,
      residuals = NULL,
      constrained = selected$constrained,
      component_passed = NULL,
      passed = FALSE,
      stability_delta = as.list(stability_delta),
      stability_max_abs = NA_real_
    )
    passed <- FALSE
  } else {
    verification_values <- .dp_backcheck_pmf(
      verification_pmf, J, "verification_pmf"
    )
    verification_core <- .dp_backcheck_one(
      verification_values, interval, tolerance
    )
    stability_delta <- unlist(verification_core$achieved) -
      unlist(selected$achieved)
    verification <- c(list(
      performed = TRUE,
      reason = NULL,
      method = "separately supplied verification PMF"
    ), verification_core, list(
      stability_delta = as.list(stability_delta),
      stability_max_abs = max(abs(stability_delta))
    ))
    passed <- isTRUE(selected$passed) && isTRUE(verification_core$passed) &&
      max(abs(stability_delta)) <= tolerance
  }
  status <- if (passed) "converged" else "approximate"

  structure(list(
    kind = "K_interval_backcheck",
    source = source,
    J = J,
    support = seq_len(J),
    interval = interval,
    selected = selected,
    verification = verification,
    tolerances = list(
      interval_mass_abs = tolerance,
      selected_verification_abs = tolerance,
      pmf_sum_abs = .TOL_PMF_SUM
    ),
    status = status,
    usable = passed,
    verified = passed,
    message = if (passed) {
      "Requested interval constraints were independently reproduced."
    } else if (is.null(verification_pmf)) {
      paste0(
        "Selected-order interval values were computed, but independent ",
        "verification was not supplied."
      )
    } else {
      "Achieved interval behavior misses at least one declared tolerance."
    },
    provenance = list(
      requested_interval_immutable = TRUE,
      selected_and_verification_separate = TRUE,
      moment_matching_assumed_interval_match = FALSE
    )
  ), class = c("dpprior_K_interval_backcheck", "list"))
}
