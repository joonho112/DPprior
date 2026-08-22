# =============================================================================
# Module 20: canonical elicitation sensitivity
# =============================================================================
#
# This internal runner evaluates how prior quantities change across explicit
# elicitation requests.  Every scenario is normalized transactionally before
# calibration starts, and every scientific value is derived from a canonical
# producer or freshly recomputed from retained canonical fit parameters.
# Compatibility views are never a scientific authority here.
# =============================================================================

.dp_sensitivity_closed_statuses <- c(
  "converged", "boundary", "approximate", "infeasible", "failed"
)


# --- Deterministic scientific identity -------------------------------------

.dp_sensitivity_encode_character <- function(value) {
  if (is.na(value)) return("NA")
  encoded <- enc2utf8(as.character(value))
  paste0(nchar(encoded, type = "bytes"), ":", encoded)
}


.dp_sensitivity_encode_atomic <- function(value) {
  type <- typeof(value)
  values <- switch(
    type,
    logical = ifelse(is.na(value), "NA", ifelse(value, "T", "F")),
    integer = ifelse(is.na(value), "NA", as.character(value)),
    double = vapply(value, function(one) {
      if (is.na(one) && !is.nan(one)) return("NA")
      if (is.nan(one)) return("NaN")
      if (is.infinite(one)) return(if (one > 0) "+Inf" else "-Inf")
      sprintf("%a", one)
    }, character(1L)),
    character = vapply(
      value, .dp_sensitivity_encode_character, character(1L)
    ),
    complex = vapply(value, function(one) {
      if (is.na(one)) return("NA")
      paste0(sprintf("%a", Re(one)), "+", sprintf("%a", Im(one)), "i")
    }, character(1L)),
    raw = sprintf("%02x", as.integer(value)),
    stop(sprintf("Unsupported atomic scenario type: %s", type), call. = FALSE)
  )
  value_names <- names(value)
  name_record <- if (is.null(value_names)) {
    "unnamed"
  } else {
    paste(vapply(
      value_names, .dp_sensitivity_encode_character, character(1L)
    ), collapse = ",")
  }
  value_attributes <- attributes(value)
  if (!is.null(value_attributes)) {
    value_attributes <- value_attributes[
      setdiff(names(value_attributes), "names")
    ]
  }
  attribute_record <- if (is.null(value_attributes) ||
                          !length(value_attributes)) {
    "noattrs"
  } else {
    .dp_sensitivity_canonical_string(value_attributes)
  }
  paste0(
    "atomic<", type, ">[", length(value), "]{", name_record, "}<",
    attribute_record, ">(", paste(values, collapse = ","), ")"
  )
}


.dp_sensitivity_canonical_string <- function(value) {
  if (is.null(value)) return("NULL")
  if (is.atomic(value)) return(.dp_sensitivity_encode_atomic(value))
  if (!is.list(value)) {
    stop(
      sprintf("Unsupported scenario content class: %s", class(value)[[1L]]),
      call. = FALSE
    )
  }
  value_names <- names(value)
  if (is.null(value_names) || anyNA(value_names) || any(!nzchar(value_names)) ||
      anyDuplicated(value_names)) {
    stop("Scenario lists must have unique, non-empty names", call. = FALSE)
  }
  ordering <- order(value_names, method = "radix")
  value_names <- value_names[ordering]
  value <- value[ordering]
  value_attributes <- attributes(value)
  if (!is.null(value_attributes)) {
    value_attributes <- value_attributes[
      setdiff(names(value_attributes), "names")
    ]
  }
  attribute_record <- if (is.null(value_attributes) ||
                          !length(value_attributes)) {
    "noattrs"
  } else {
    .dp_sensitivity_canonical_string(value_attributes)
  }
  entries <- vapply(seq_along(value), function(index) {
    paste0(
      .dp_sensitivity_encode_character(value_names[[index]]), "=",
      .dp_sensitivity_canonical_string(value[[index]])
    )
  }, character(1L))
  paste0("list<", attribute_record, ">{", paste(entries, collapse = ";"), "}")
}


.dp_sensitivity_semantic_number <- function(value, integer = FALSE) {
  if (!is.numeric(value) || length(value) != 1L || !is.null(dim(value)) ||
      !all(class(value) %in% c("numeric", "integer"))) return(value)
  value <- unname(value)
  if (integer && is.finite(value) && value == floor(value)) {
    return(as.integer(value))
  }
  as.numeric(value)
}


.dp_sensitivity_identity_request <- function(request) {
  normalized <- request
  for (field in c("mu_K", "var_K", "cv_K")) {
    if (field %in% names(normalized)) {
      normalized[[field]] <- .dp_sensitivity_semantic_number(
        normalized[[field]]
      )
    }
  }
  for (field in c("J", "M")) {
    if (field %in% names(normalized)) {
      normalized[[field]] <- .dp_sensitivity_semantic_number(
        normalized[[field]], integer = TRUE
      )
    }
  }
  for (field in c("method", "confidence")) {
    value <- normalized[[field]]
    if (is.character(value) && length(value) == 1L &&
        all(class(value) == "character")) {
      normalized[[field]] <- unname(value)
    }
  }
  if (is.numeric(normalized[["target_pmf"]]) &&
      !is.object(normalized[["target_pmf"]]) &&
      is.null(dim(normalized[["target_pmf"]]))) {
    normalized[["target_pmf"]] <- unname(as.numeric(
      normalized[["target_pmf"]]
    ))
  }
  interval <- normalized[["K_interval"]]
  if (is.list(interval)) {
    for (field in c("lower", "upper")) {
      if (field %in% names(interval)) {
        interval[[field]] <- .dp_sensitivity_semantic_number(
          interval[[field]], integer = TRUE
        )
      }
    }
    for (field in c("coverage", "mu_K")) {
      if (field %in% names(interval)) {
        interval[[field]] <- .dp_sensitivity_semantic_number(interval[[field]])
      }
    }
    normalized[["K_interval"]] <- interval
  }
  normalized
}


.dp_sensitivity_mod_hash <- function(text, base, modulus) {
  bytes <- as.integer(charToRaw(enc2utf8(text)))
  hash <- 0
  for (byte in bytes) hash <- (hash * base + byte + 1) %% modulus
  as.integer(hash)
}


.dp_sensitivity_content_key <- function(canonical_string) {
  first <- .dp_sensitivity_mod_hash(canonical_string, 131, 2147483629)
  second <- .dp_sensitivity_mod_hash(canonical_string, 137, 2147483587)
  paste0("scn_", sprintf("%08x", first), sprintf("%08x", second))
}


# --- Conditions and exact canonical views ----------------------------------

.dp_sensitivity_new_condition <- function(message, classes, code, ...) {
  .dpprior_new_condition(
    message = message,
    classes = unique(c(
      classes, "dpprior_sensitivity_error", "dpprior_error", "error"
    )),
    code = code, ...
  )
}


.dp_sensitivity_abort <- function(message, classes, code, ...) {
  stop(.dp_sensitivity_new_condition(message, classes, code, ...))
}


.dp_sensitivity_backend_condition <- function(message, result = NULL) {
  .dpprior_new_condition(
    message = message,
    classes = c(
      "dpprior_sensitivity_backend_contract_error",
      "dpprior_backend_contract_error", "dpprior_calibration_error",
      "dpprior_sensitivity_error", "dpprior_error", "error"
    ),
    code = "sensitivity_backend_contract", result = result
  )
}


.dp_sensitivity_diagnostic_condition <- function(message, result = NULL) {
  .dpprior_new_condition(
    message = message,
    classes = c(
      "dpprior_sensitivity_diagnostic_contract_error",
      "dpprior_diagnostics_error", "dpprior_sensitivity_error",
      "dpprior_error", "error"
    ),
    code = "sensitivity_diagnostic_contract", result = result
  )
}


.dp_sensitivity_capture <- function(expression) {
  warnings <- list()
  condition <- NULL
  value <- withCallingHandlers(
    tryCatch(expression, error = function(error) {
      condition <<- error
      NULL
    }),
    warning = function(warning) {
      warnings[[length(warnings) + 1L]] <<- warning
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, condition = condition, warnings = warnings)
}


.dp_sensitivity_condition_message <- function(condition) {
  raw <- tryCatch(unclass(condition), error = function(error) NULL)
  message <- tryCatch(
    raw[["message", exact = TRUE]], error = function(error) NULL
  )
  if (!.dp_sensitivity_plain_scalar_character(message)) {
    return("Backend condition did not retain one ordinary message")
  }
  unname(message)
}


.dp_sensitivity_condition_summary <- function(condition, warning = FALSE) {
  if (is.null(condition)) return(NULL)
  classes <- class(condition)
  raw <- tryCatch(unclass(condition), error = function(error) NULL)
  code <- tryCatch(
    raw[["code", exact = TRUE]], error = function(error) NULL
  )
  if (!is.character(code) || is.object(code) || !is.null(dim(code)) ||
      length(code) != 1L || is.na(code) || !nzchar(code)) {
    code <- if (warning) "sensitivity_backend_warning" else
      "sensitivity_backend_contract"
  }
  if (warning && !("warning" %in% classes)) classes <- c(classes, "warning")
  list(
    class = classes[[1L]], classes = classes, code = unname(code),
    message = .dp_sensitivity_condition_message(condition)
  )
}


.dp_sensitivity_empty_conditions <- function() {
  list(
    calibration = NULL, calibration_warnings = list(), diagnostics = NULL,
    diagnostic_warnings = list(), target = NULL, interval = NULL
  )
}


.dp_sensitivity_gate <- function(value, kind) {
  if (is.null(value)) return(NULL)
  tryCatch(
    .dpprior_require_schema(value, kind = kind, allow_legacy = FALSE),
    error = function(error) NULL
  )
}


.dp_sensitivity_target_projection <- function(target) {
  target <- .dp_sensitivity_gate(target, "target")
  if (is.null(target)) return(NULL)
  raw <- unclass(target)
  list(
    kind = raw[["kind", exact = TRUE]],
    J = raw[["J", exact = TRUE]],
    request = raw[["request", exact = TRUE]],
    used = raw[["used", exact = TRUE]]
  )
}


.dp_sensitivity_target_call <- function(request) {
  request_names <- names(request)
  args <- list(J = request[["J", exact = TRUE]])
  if ("target_pmf" %in% request_names) {
    args[["target_pmf"]] <- request[["target_pmf", exact = TRUE]]
  } else if ("K_interval" %in% request_names) {
    interval <- request[["K_interval", exact = TRUE]]
    if (is.list(interval) && !is.object(interval)) {
      interval <- interval[intersect(
        c("lower", "upper", "type", "coverage", "family", "mu_K"),
        names(interval)
      )]
      if (identical(interval[["type", exact = TRUE]], "hard_bounds") &&
          identical(interval[["coverage", exact = TRUE]], 1)) {
        interval[["coverage"]] <- NULL
      }
    }
    args[["K_interval"]] <- interval
    args[["mu_K"]] <- request[["mu_K", exact = TRUE]]
  } else if ("cv_K" %in% request_names) {
    args[["mu_K"]] <- request[["mu_K", exact = TRUE]]
    args[["cv_K"]] <- request[["cv_K", exact = TRUE]]
  } else if ("confidence" %in% request_names) {
    args[["mu_K"]] <- request[["mu_K", exact = TRUE]]
    args[["confidence"]] <- request[["confidence", exact = TRUE]]
  } else {
    args[["mu_K"]] <- request[["mu_K", exact = TRUE]]
    args[["var_K"]] <- request[["var_K", exact = TRUE]]
  }
  args
}


.dp_sensitivity_canonicalize_scientific_target <- function(
    scientific_request) {
  if (typeof(scientific_request) != "list" ||
      !is.list(scientific_request) || is.object(scientific_request) ||
      is.null(dim(scientific_request)) == FALSE ||
      is.null(names(scientific_request)) || anyNA(names(scientific_request)) ||
      any(!nzchar(names(scientific_request))) ||
      anyDuplicated(names(scientific_request))) {
    return(list(
      target = NULL,
      violation = "scientific request must be one ordinary uniquely named list"
    ))
  }
  attempt <- .dp_sensitivity_capture(do.call(
    .dp_target_K, .dp_sensitivity_target_call(scientific_request)
  ))
  candidate <- attempt$value
  if (is.null(candidate) && !is.null(attempt$condition)) {
    candidate <- .dp_sensitivity_condition_result(attempt$condition)
  }
  candidate <- .dp_sensitivity_gate(candidate, "target")
  if (is.null(candidate)) {
    detail <- if (is.null(attempt$condition)) {
      "constructor returned no canonical target"
    } else {
      .dp_sensitivity_condition_message(attempt$condition)
    }
    return(list(
      target = NULL,
      violation = paste("scientific request could not be canonicalized:", detail)
    ))
  }
  list(target = candidate, violation = NULL)
}


# --- Transactional scenario preflight --------------------------------------

.dp_sensitivity_plain_scalar_character <- function(x, choices = NULL) {
  valid <- is.character(x) && !is.object(x) && is.null(dim(x)) &&
    length(x) == 1L && !is.na(x) && nzchar(x)
  valid && (is.null(choices) || x %in% choices)
}


.dp_sensitivity_plain_scalar_logical <- function(x) {
  is.logical(x) && !is.object(x) && is.null(dim(x)) &&
    length(x) == 1L && !is.na(x)
}


.dp_sensitivity_plain_numeric <- function(x, length = NULL) {
  is.numeric(x) && !is.object(x) && is.null(dim(x)) &&
    (is.null(length) || base::length(x) == length)
}


.dp_sensitivity_preflight_abort <- function(message, code, index = NULL, ...) {
  .dp_sensitivity_abort(
    message,
    c(
      "dpprior_sensitivity_scenarios_error",
      "dpprior_sensitivity_scenario_error", "dpprior_invalid_input"
    ),
    code, scenario_index = index, ...
  )
}


.dp_sensitivity_as_scenario_list <- function(scenarios) {
  if (is.data.frame(scenarios)) {
    if (nrow(scenarios) < 1L) {
      .dp_sensitivity_preflight_abort(
        "scenarios must contain at least one row", "empty_scenario_container"
      )
    }
    return(lapply(seq_len(nrow(scenarios)), function(index) {
      stats::setNames(
        lapply(scenarios, function(column) column[[index]]), names(scenarios)
      )
    }))
  }
  if (typeof(scenarios) != "list" || !is.list(scenarios) ||
      is.object(scenarios) || is.null(dim(scenarios)) == FALSE ||
      length(scenarios) < 1L ||
      !all(vapply(scenarios, function(one) {
        typeof(one) == "list" && is.list(one) && !is.object(one) &&
          is.null(dim(one))
      }, logical(1L)))) {
    .dp_sensitivity_preflight_abort(
      "scenarios must be a non-empty ordinary list of ordinary lists",
      "scenario_container"
    )
  }
  scenarios
}


.dp_sensitivity_normalize_count <- function(value, name, index,
                                             minimum, maximum) {
  if (!is.numeric(value) || is.object(value) || !is.null(dim(value)) ||
      length(value) != 1L || is.na(value) || !is.finite(value) ||
      value != floor(value) || value < minimum || value > maximum) {
    .dp_sensitivity_preflight_abort(
      sprintf("%s must be one integer in [%d, %d]", name, minimum, maximum),
      paste0("invalid_", name), index, field = name, value = value
    )
  }
  as.integer(value)
}


.dp_sensitivity_normalize_numeric <- function(value, name, index,
                                               lower = -Inf,
                                               upper = Inf,
                                               open_lower = FALSE) {
  if (!.dp_sensitivity_plain_numeric(value, 1L) || is.na(value) ||
      !is.finite(value) ||
      (if (open_lower) value <= lower else value < lower) || value > upper) {
    .dp_sensitivity_preflight_abort(
      sprintf("%s must be one finite in-domain number", name),
      paste0("invalid_", name), index, field = name, value = value
    )
  }
  as.numeric(value)
}


.dp_sensitivity_target_route <- function(scenario, index) {
  sources <- c(
    direct_variance = "var_K",
    qualitative_confidence = "confidence",
    coefficient_of_variation = "cv_K",
    interval = "K_interval",
    strict_pmf = "target_pmf"
  )
  present <- names(sources)[vapply(unname(sources), function(field) {
    field %in% names(scenario) && !is.null(scenario[[field, exact = TRUE]])
  }, logical(1L))]
  if (length(present) == 0L) return("qualitative_confidence")
  if (length(present) != 1L) {
    .dp_sensitivity_preflight_abort(
      "Each scenario must select exactly one uncertainty route",
      "multiple_uncertainty_sources", index, routes = present
    )
  }
  present[[1L]]
}


.dp_sensitivity_normalize_interval <- function(value, J, mu_K, index) {
  attempt <- .dp_sensitivity_capture(.dp_target_K(
    J = J, mu_K = mu_K, K_interval = value
  ))
  target <- attempt$value
  if (is.null(target) && !is.null(attempt$condition)) {
    target <- .dp_sensitivity_condition_result(attempt$condition)
  }
  target <- .dp_sensitivity_gate(target, "target")
  if (is.null(target)) {
    .dp_sensitivity_preflight_abort(
      if (is.null(attempt$condition)) {
        "K_interval did not produce a canonical target"
      } else {
        .dp_sensitivity_condition_message(attempt$condition)
      },
      "invalid_interval_request", index, field = "K_interval"
    )
  }
  raw <- unclass(target)
  interval <- raw[["used", exact = TRUE]][["interval", exact = TRUE]]
  if (is.null(interval)) {
    interval <- raw[["interval", exact = TRUE]]
  }
  list(target = target, interval = interval)
}


.dp_sensitivity_normalize_one <- function(
    scenario, index, J, method, M, check_diagnostics) {
  scenario_names <- names(scenario)
  if (is.null(scenario_names) || anyNA(scenario_names) ||
      any(!nzchar(scenario_names)) || anyDuplicated(scenario_names)) {
    .dp_sensitivity_preflight_abort(
      "Each scenario must have unique non-empty names",
      "scenario_names", index
    )
  }
  allowed <- c(
    "J", "mu_K", "var_K", "confidence", "cv_K", "K_interval",
    "target_pmf", "method", "M", "check_diagnostics", "scenario_label",
    "label", "weight_target"
  )
  unknown <- setdiff(scenario_names, allowed)
  if (length(unknown)) {
    .dp_sensitivity_preflight_abort(
      sprintf("Unknown sensitivity scenario field(s): %s",
              paste(unknown, collapse = ", ")),
      "unknown_scenario_fields", index, fields = unknown
    )
  }
  if ("weight_target" %in% scenario_names &&
      !is.null(scenario[["weight_target", exact = TRUE]])) {
    .dp_sensitivity_preflight_abort(
      "weight_target sensitivity requires a future typed authority",
      "weight_target_authority_deferred", index, field = "weight_target"
    )
  }
  scenario_J <- if ("J" %in% scenario_names) {
    .dp_sensitivity_normalize_count(
      scenario[["J", exact = TRUE]], "J", index, 1L, .MAX_J_DEFAULT
    )
  } else {
    J
  }
  if (!identical(scenario_J, J)) {
    .dp_sensitivity_preflight_abort(
      "Every scenario J must exactly equal the runner J",
      "scenario_J_mismatch", index, value = scenario_J, expected = J
    )
  }
  scenario_M <- if ("M" %in% scenario_names) {
    .dp_sensitivity_normalize_count(
      scenario[["M", exact = TRUE]], "M", index, 10L, 256L
    )
  } else {
    M
  }
  diagnostics_requested <- if ("check_diagnostics" %in% scenario_names) {
    value <- scenario[["check_diagnostics", exact = TRUE]]
    if (!.dp_sensitivity_plain_scalar_logical(value)) {
      .dp_sensitivity_preflight_abort(
        "check_diagnostics must be one ordinary non-missing logical",
        "scenario_diagnostics_control", index,
        field = "check_diagnostics", value = value
      )
    }
    isTRUE(value)
  } else {
    check_diagnostics
  }

  route <- .dp_sensitivity_target_route(scenario, index)
  method_explicit <- "method" %in% scenario_names || !is.null(method)
  requested_method <- if ("method" %in% scenario_names) {
    scenario[["method", exact = TRUE]]
  } else if (!is.null(method)) {
    method
  } else if (route %in% c("interval", "strict_pmf")) {
    "A2-KL"
  } else {
    "A2-MN"
  }
  if (!.dp_sensitivity_plain_scalar_character(
    requested_method, c("A2-MN", "A1", "A2-KL")
  )) {
    .dp_sensitivity_preflight_abort(
      "method must be exactly A2-MN, A1, or A2-KL",
      "sensitivity_method", index, field = "method",
      value = requested_method
    )
  }
  if (route %in% c("interval", "strict_pmf") &&
      !identical(requested_method, "A2-KL")) {
    .dp_sensitivity_preflight_abort(
      "Interval and strict-PMF routes require method = 'A2-KL'",
      "target_pmf_wrong_method", index
    )
  }
  if (identical(route, "coefficient_of_variation") &&
      !requested_method %in% c("A1", "A2-MN")) {
    .dp_sensitivity_preflight_abort(
      "The coefficient-of-variation route supports A1 or A2-MN",
      "cv_method_unsupported", index
    )
  }

  mu_present <- "mu_K" %in% scenario_names &&
    !is.null(scenario[["mu_K", exact = TRUE]])
  mu_K <- if (mu_present) {
    .dp_sensitivity_normalize_numeric(
      scenario[["mu_K", exact = TRUE]], "mu_K", index, 1, J
    )
  } else {
    NULL
  }
  confidence_explicit <- identical(route, "qualitative_confidence") &&
    "confidence" %in% scenario_names

  request <- switch(
    route,
    direct_variance = {
      if (is.null(mu_K)) {
        .dp_sensitivity_preflight_abort(
          "Direct variance scenarios require mu_K", "missing_mu_K", index
        )
      }
      variance <- .dp_sensitivity_normalize_numeric(
        scenario[["var_K", exact = TRUE]], "var_K", index, 0, Inf
      )
      list(
        J = J, mu_K = mu_K, var_K = variance,
        method = unname(requested_method), M = scenario_M
      )
    },
    qualitative_confidence = {
      if (is.null(mu_K)) {
        .dp_sensitivity_preflight_abort(
          "Confidence scenarios require mu_K", "missing_mu_K", index
        )
      }
      confidence <- if (confidence_explicit) {
        scenario[["confidence", exact = TRUE]]
      } else {
        "medium"
      }
      if (!.dp_sensitivity_plain_scalar_character(
        confidence, c("low", "medium", "high")
      )) {
        .dp_sensitivity_preflight_abort(
          "confidence must be low, medium, or high",
          "invalid_confidence", index, value = confidence
        )
      }
      list(
        J = J, mu_K = mu_K, confidence = unname(confidence),
        method = unname(requested_method), M = scenario_M
      )
    },
    coefficient_of_variation = {
      if (is.null(mu_K)) {
        .dp_sensitivity_preflight_abort(
          "Coefficient-of-variation scenarios require mu_K",
          "missing_mu_K", index
        )
      }
      cv <- .dp_sensitivity_normalize_numeric(
        scenario[["cv_K", exact = TRUE]], "cv_K", index, 0, Inf,
        open_lower = TRUE
      )
      list(
        J = J, mu_K = mu_K, cv_K = cv,
        method = unname(requested_method), M = scenario_M
      )
    },
    interval = {
      interval <- .dp_sensitivity_normalize_interval(
        scenario[["K_interval", exact = TRUE]], J, mu_K, index
      )
      list(
        J = J, K_interval = interval$interval,
        mu_K = interval$interval[["mu_K", exact = TRUE]],
        method = unname(requested_method), M = scenario_M
      )
    },
    strict_pmf = {
      pmf <- scenario[["target_pmf", exact = TRUE]]
      if (!.dp_sensitivity_plain_numeric(pmf) || is.null(names(pmf)) == FALSE ||
          !length(pmf) %in% c(J, J + 1L) || anyNA(pmf) ||
          any(!is.finite(pmf)) || any(pmf < 0) ||
          (length(pmf) == J + 1L && !identical(as.numeric(pmf[[1L]]), 0))) {
        .dp_sensitivity_preflight_abort(
          "target_pmf must be an unnamed finite nonnegative J or J+1 vector",
          "invalid_target_pmf", index, field = "target_pmf"
        )
      }
      list(
        J = J, target_pmf = unname(as.numeric(pmf)),
        method = unname(requested_method), M = scenario_M
      )
    }
  )

  allowed_route_fields <- switch(
    route,
    direct_variance = c("mu_K", "var_K"),
    qualitative_confidence = c("mu_K", "confidence"),
    coefficient_of_variation = c("mu_K", "cv_K"),
    interval = c("mu_K", "K_interval"),
    strict_pmf = "target_pmf"
  )
  scientific_present <- intersect(
    c("mu_K", "var_K", "confidence", "cv_K", "K_interval", "target_pmf"),
    scenario_names
  )
  extraneous <- setdiff(scientific_present, allowed_route_fields)
  if (length(extraneous)) {
    .dp_sensitivity_preflight_abort(
      sprintf("Route %s does not accept: %s", route,
              paste(extraneous, collapse = ", ")),
      "route_field_conflict", index, fields = extraneous
    )
  }

  canonical <- .dp_sensitivity_canonicalize_scientific_target(request)
  if (!is.null(canonical$violation)) {
    .dp_sensitivity_preflight_abort(
      canonical$violation, "target_preflight_failed", index
    )
  }
  target <- canonical$target
  target_projection <- .dp_sensitivity_target_projection(target)
  target_raw <- unclass(target)
  if (target_raw[["status", exact = TRUE]] %in% c("infeasible", "failed")) {
    # No fit exists to carry diagnostics on a target-only terminal route.
    diagnostics_requested <- FALSE
  }
  if (identical(J, 1L) && identical(route, "strict_pmf")) {
    diagnostics_requested <- FALSE
  }
  content <- .dp_sensitivity_canonical_string(list(
    request = .dp_sensitivity_identity_request(request),
    diagnostics_requested = diagnostics_requested,
    weight_target = NULL
  ))
  list(
    request = request, target = target, target_projection = target_projection,
    diagnostics_requested = diagnostics_requested,
    input_provenance = list(
      method_explicit = method_explicit,
      confidence_explicit = confidence_explicit,
      requested_method = unname(requested_method),
      selected_method = unname(requested_method),
      is_fallback = FALSE, target_route = route
    ),
    canonical_content = content
  )
}


.dp_sensitivity_checked_key <- function(key_fun, canonical_string) {
  attempt <- tryCatch(key_fun(canonical_string), error = function(error) error)
  if (inherits(attempt, "condition")) {
    .dp_sensitivity_abort(
      paste(
        "Scenario key adapter failed:",
        .dp_sensitivity_condition_message(attempt)
      ),
      c("dpprior_sensitivity_key_error", "dpprior_sensitivity_adapter_error"),
      "scenario_key_function"
    )
  }
  if (!.dp_sensitivity_plain_scalar_character(attempt) ||
      !grepl("^scn_[0-9a-f]{16}$", attempt)) {
    .dp_sensitivity_abort(
      "Scenario key adapter must return one canonical scn_<16 hex> key",
      c("dpprior_sensitivity_key_error", "dpprior_sensitivity_adapter_error"),
      "scenario_key_value", value = attempt
    )
  }
  expected <- .dp_sensitivity_content_key(canonical_string)
  if (!identical(unname(attempt), expected)) {
    .dp_sensitivity_abort(
      "Scenario key adapter must preserve the canonical content hash",
      c("dpprior_sensitivity_key_error", "dpprior_sensitivity_adapter_error"),
      "scenario_key_noncanonical", value = unname(attempt),
      expected = expected
    )
  }
  expected
}


.dp_sensitivity_prepare_scenarios <- function(
    J, scenarios, method, M, check_diagnostics, key_fun) {
  raw <- .dp_sensitivity_as_scenario_list(scenarios)
  records <- lapply(seq_along(raw), function(index) {
    .dp_sensitivity_normalize_one(
      raw[[index]], index, J, method, M, check_diagnostics
    )
  })
  canonical <- vapply(
    records, `[[`, character(1L), "canonical_content"
  )
  if (anyDuplicated(canonical)) {
    duplicates <- which(duplicated(canonical) |
                          duplicated(canonical, fromLast = TRUE))
    .dp_sensitivity_abort(
      "Duplicate normalized scientific scenarios are not allowed",
      c(
        "dpprior_sensitivity_duplicate_scenario",
        "dpprior_conflicting_input"
      ),
      "duplicate_scenario_content", duplicate_indices = duplicates
    )
  }
  base_keys <- vapply(
    canonical, function(one) .dp_sensitivity_checked_key(key_fun, one),
    character(1L)
  )
  final_keys <- base_keys
  collisions <- unique(base_keys[
    duplicated(base_keys) | duplicated(base_keys, fromLast = TRUE)
  ])
  for (base_key in collisions) {
    members <- which(base_keys == base_key)
    lexical <- order(canonical[members], method = "radix")
    assigned <- character()
    for (position in seq_along(lexical)) {
      member <- members[[lexical[[position]]]]
      suffix <- position
      repeat {
        candidate <- paste0(base_key, "_c", sprintf("%03d", suffix))
        if (!(candidate %in% base_keys) && !(candidate %in% final_keys[
          setdiff(seq_along(final_keys), members)
        ]) && !(candidate %in% assigned)) break
        suffix <- suffix + 1L
      }
      final_keys[[member]] <- candidate
      assigned <- c(assigned, candidate)
    }
  }
  for (index in seq_along(records)) {
    records[[index]][["base_key"]] <- base_keys[[index]]
    records[[index]][["scenario_key"]] <- final_keys[[index]]
  }
  records[order(final_keys, method = "radix")]
}


.dp_sensitivity_validate_fit_adapter <- function(fit_fun, records) {
  if (!is.function(fit_fun)) {
    .dp_sensitivity_abort(
      ".fit_fun must be a function",
      c("dpprior_sensitivity_adapter_error", "dpprior_type_error"),
      "adapter_type"
    )
  }
  formal_names <- names(formals(fit_fun))
  if (is.null(formal_names) || "..." %in% formal_names) return(invisible(TRUE))
  required <- unique(unlist(lapply(records, function(record) {
    c(names(record$request), "check_diagnostics", "verbose")
  }), use.names = FALSE))
  unavailable <- setdiff(required, formal_names)
  if (length(unavailable)) {
    .dp_sensitivity_abort(
      sprintf("The fit adapter does not accept: %s",
              paste(unavailable, collapse = ", ")),
      c("dpprior_sensitivity_adapter_error", "dpprior_type_error"),
      "fit_adapter_formals", unavailable = unavailable
    )
  }
  invisible(TRUE)
}


# --- Canonical fit retention and fresh fixed-order evidence ----------------

.dp_sensitivity_condition_result <- function(condition) {
  if (is.null(condition)) return(NULL)
  raw <- tryCatch(unclass(condition), error = function(error) NULL)
  tryCatch(
    raw[["result", exact = TRUE]], error = function(error) NULL
  )
}


.dp_sensitivity_condition_code <- function(condition) {
  if (is.null(condition)) return(NULL)
  raw <- tryCatch(unclass(condition), error = function(error) NULL)
  code <- tryCatch(
    raw[["code", exact = TRUE]], error = function(error) NULL
  )
  if (.dp_sensitivity_plain_scalar_character(code)) unname(code) else NULL
}


.dp_sensitivity_expected_fit_mode <- function(method) {
  switch(method, A1 = "a1_proxy", `A2-MN` = "a2_moment", `A2-KL` = "a2_kl")
}


.dp_sensitivity_fit_cross_bound <- function(fit, record) {
  fit <- .dp_sensitivity_gate(fit, "fit")
  if (is.null(fit)) return(FALSE)
  raw <- unclass(fit)
  requested <- record$input_provenance$requested_method
  selected <- raw[["method", exact = TRUE]]
  method_ok <- if (identical(requested, "A2-MN")) {
    selected %in% c("A2-MN", "A2-MN+NM")
  } else {
    identical(selected, requested)
  }
  target <- raw[["target", exact = TRUE]][["K", exact = TRUE]]
  target_projection <- .dp_sensitivity_target_projection(target)
  identical(raw[["object_type", exact = TRUE]], "fit") &&
    identical(raw[["mode", exact = TRUE]],
              .dp_sensitivity_expected_fit_mode(requested)) &&
    method_ok && identical(raw[["J", exact = TRUE]], record$request$J) &&
    identical(target_projection, record$target_projection)
}


.dp_sensitivity_fresh_snapshots <- function(parameters, J, M) {
  M_verification <- as.integer(max(2L * M, M + 40L))
  attempt <- tryCatch({
    logS <- compute_log_stirling(J)
    selected <- pmf_K_marginal(
      J, parameters[["a", exact = TRUE]], parameters[["b", exact = TRUE]],
      logS, M = M, M_verify = M_verification,
      abs_tol = 1e-10, rel_tol = 1e-8, strict = FALSE
    )
    verifier <- attr(
      selected, ".marginal_verification_pmf", exact = TRUE
    )
    selected <- unname(as.numeric(selected[-1L]))
    verifier <- unname(as.numeric(verifier[-1L]))
    selected_moments <- .dpprior_target_pmf_moments(selected)
    verifier_moments <- .dpprior_target_pmf_moments(verifier)
    list(
      evaluator = list(
        method = "gauss_laguerre_marginal_pmf_and_moments",
        M_selected = M, M_verification = M_verification,
        absolute_tolerance = 1e-10, relative_tolerance = 1e-8,
        W_max_point_policy =
          "unavailable_without_retained_typed_backend_evidence"
      ),
      selected_snapshot = list(
        M = M,
        K = list(
          mean = unname(selected_moments[["mean", exact = TRUE]]),
          variance = unname(selected_moments[["variance", exact = TRUE]]),
          pmf = selected
        ),
        finite = TRUE, source = "selected_order"
      ),
      verifier_snapshot = list(
        M = M_verification,
        K = list(
          mean = unname(verifier_moments[["mean", exact = TRUE]]),
          variance = unname(verifier_moments[["variance", exact = TRUE]]),
          pmf = verifier
        ),
        finite = TRUE, source = "independent_verifier"
      )
    )
  }, error = function(error) error)
  attempt
}


.dp_sensitivity_evaluator <- function(M) {
  list(
    method = "gauss_laguerre_marginal_pmf_and_moments",
    M_selected = M,
    M_verification = as.integer(max(2L * M, M + 40L)),
    absolute_tolerance = 1e-10, relative_tolerance = 1e-8,
    W_max_point_policy =
      "unavailable_without_retained_typed_backend_evidence"
  )
}


.dp_sensitivity_native_diagnostics_present <- function(raw) {
  aliases <- raw[["compatibility", exact = TRUE]][[
    "top_level_aliases", exact = TRUE
  ]]
  alias_names <- names(aliases)
  if (is.null(alias_names)) alias_names <- character()
  "diagnostics" %in% names(raw) && !("diagnostics" %in% alias_names)
}


.dp_sensitivity_failed_evidence <- function(record, condition,
                                             conditions = NULL) {
  if (is.null(conditions)) conditions <- .dp_sensitivity_empty_conditions()
  conditions$calibration <- .dp_sensitivity_condition_summary(condition)
  provenance <- record$input_provenance
  list(
    request = record$request, input_provenance = provenance,
    weight_target = NULL, target = record$target_projection,
    method = provenance$selected_method,
    status = "failed", usable = FALSE, verified = FALSE,
    parameters = NULL, evaluator = .dp_sensitivity_evaluator(record$request$M),
    selected_snapshot = NULL, verifier_snapshot = NULL,
    condition_evidence = conditions,
    source = "retained_canonical_fit_evidence"
  )
}


.dp_sensitivity_execute_one <- function(record, fit_fun) {
  call_args <- .dp_sensitivity_target_call(record$request)
  call_args[["method"]] <- record$request[["method", exact = TRUE]]
  call_args[["M"]] <- record$request[["M", exact = TRUE]]
  call_args[["check_diagnostics"]] <- record$diagnostics_requested
  call_args[["verbose"]] <- FALSE
  calibration <- .dp_sensitivity_capture(do.call(fit_fun, call_args))
  condition <- calibration$condition
  condition_result <- .dp_sensitivity_condition_result(condition)

  returned_fit <- .dp_sensitivity_gate(calibration$value, "fit")
  returned_target <- .dp_sensitivity_gate(calibration$value, "target")
  condition_fit <- .dp_sensitivity_gate(condition_result, "fit")
  condition_target <- .dp_sensitivity_gate(condition_result, "target")
  fit <- if (!is.null(returned_fit)) returned_fit else condition_fit
  target <- if (!is.null(returned_target)) returned_target else condition_target
  conditions <- .dp_sensitivity_empty_conditions()
  warning_summaries <- lapply(
    calibration$warnings,
    .dp_sensitivity_condition_summary, warning = TRUE
  )
  conditions$calibration_warnings <- warning_summaries

  if (!is.null(fit) && !.dp_sensitivity_fit_cross_bound(fit, record)) {
    contract <- .dp_sensitivity_backend_condition(
      "Canonical fit did not match the normalized route, method, J, or target",
      result = fit
    )
    evidence <- .dp_sensitivity_failed_evidence(record, contract, conditions)
    return(list(evidence = evidence, row_failure = TRUE))
  }
  if (!is.null(target) && !identical(
    .dp_sensitivity_target_projection(target), record$target_projection
  )) {
    contract <- .dp_sensitivity_backend_condition(
      "Condition target did not match the normalized scenario target",
      result = target
    )
    evidence <- .dp_sensitivity_failed_evidence(record, contract, conditions)
    return(list(evidence = evidence, row_failure = TRUE))
  }

  condition_code <- .dp_sensitivity_condition_code(condition)
  if (is.null(fit)) {
    target_raw <- if (is.null(target)) NULL else unclass(target)
    target_status <- if (is.null(target_raw)) NULL else
      target_raw[["status", exact = TRUE]]
    recognized_target <- !is.null(target) &&
      target_status %in% c("infeasible", "failed")
    recognized_J1 <- !is.null(target) &&
      identical(condition_code, "calibration_nonidentifiable_j1") &&
      identical(record$request$J, 1L) &&
      identical(record$input_provenance$target_route, "strict_pmf")
    if (!recognized_target && !recognized_J1) {
      contract <- .dp_sensitivity_backend_condition(
        if (is.null(condition)) {
          "Fit adapter returned no canonical fit"
        } else {
          paste("Fit adapter failed without usable canonical fit evidence:",
                .dp_sensitivity_condition_message(condition))
        },
        result = target
      )
      evidence <- .dp_sensitivity_failed_evidence(record, contract, conditions)
      return(list(evidence = evidence, row_failure = TRUE))
    }
    conditions$calibration <- .dp_sensitivity_condition_summary(condition)
    if (recognized_target) {
      conditions$target <- conditions$calibration
    }
    evidence_status <- if (recognized_J1) "infeasible" else target_status
    record$input_provenance$selected_method <-
      record$input_provenance$requested_method
    evidence <- list(
      request = record$request,
      input_provenance = record$input_provenance,
      weight_target = NULL, target = record$target_projection,
      method = record$input_provenance$selected_method,
      status = evidence_status, usable = FALSE,
      verified = identical(evidence_status, "infeasible"),
      parameters = NULL,
      evaluator = .dp_sensitivity_evaluator(record$request$M),
      selected_snapshot = NULL, verifier_snapshot = NULL,
      condition_evidence = conditions,
      source = "retained_canonical_fit_evidence"
    )
    return(list(evidence = evidence, row_failure = FALSE))
  }

  raw <- unclass(fit)
  selected_method <- raw[["method", exact = TRUE]]
  record$input_provenance$selected_method <- selected_method
  record$input_provenance$is_fallback <-
    identical(record$input_provenance$requested_method, "A2-MN") &&
    identical(selected_method, "A2-MN+NM")
  parameters <- raw[["parameters", exact = TRUE]]
  snapshots <- .dp_sensitivity_fresh_snapshots(
    parameters, record$request$J, record$request$M
  )
  if (inherits(snapshots, "condition")) {
    contract <- .dp_sensitivity_backend_condition(
      paste("Fresh fixed-order fit evaluation failed:",
            .dp_sensitivity_condition_message(snapshots)),
      result = fit
    )
    evidence <- .dp_sensitivity_failed_evidence(record, contract, conditions)
    return(list(evidence = evidence, row_failure = TRUE))
  }

  diagnostic_failure <- NULL
  if (record$diagnostics_requested && is.null(condition) &&
      !.dp_sensitivity_native_diagnostics_present(raw)) {
    diagnostic_failure <- .dp_sensitivity_diagnostic_condition(
      "Requested diagnostics were not attached to the canonical fit",
      result = fit
    )
  }
  if (!is.null(condition)) {
    if (identical(condition_code, "calibration_unusable")) {
      conditions$calibration <- .dp_sensitivity_condition_summary(condition)
    } else if (identical(condition_code, "fit_diagnostics_approximate") &&
               identical(selected_method, "A1") &&
               identical(raw[["status", exact = TRUE]], "approximate")) {
      conditions$diagnostics <- .dp_sensitivity_condition_summary(condition)
    } else if (identical(condition_code, "fit_diagnostics_approximate")) {
      diagnostic_failure <- .dp_sensitivity_diagnostic_condition(
        .dp_sensitivity_condition_message(condition), result = fit
      )
    } else {
      contract <- .dp_sensitivity_backend_condition(
        paste("Unexpected condition accompanied a canonical fit:",
              .dp_sensitivity_condition_message(condition)),
        result = fit
      )
      evidence <- .dp_sensitivity_failed_evidence(record, contract, conditions)
      return(list(evidence = evidence, row_failure = TRUE))
    }
  }
  if (length(warning_summaries)) {
    diagnostic_failure <- .dp_sensitivity_diagnostic_condition(
      "Warnings accompanied sensitivity calibration; claims were quarantined",
      result = fit
    )
  }
  if (!is.null(diagnostic_failure)) {
    conditions$diagnostics <-
      .dp_sensitivity_condition_summary(diagnostic_failure)
  }

  evidence <- list(
    request = record$request,
    input_provenance = record$input_provenance,
    weight_target = NULL, target = record$target_projection,
    method = selected_method,
    status = raw[["status", exact = TRUE]],
    usable = raw[["usable", exact = TRUE]],
    verified = raw[["verified", exact = TRUE]],
    parameters = parameters,
    evaluator = snapshots$evaluator,
    selected_snapshot = snapshots$selected_snapshot,
    verifier_snapshot = snapshots$verifier_snapshot,
    condition_evidence = conditions,
    source = "retained_canonical_fit_evidence"
  )
  list(evidence = evidence, row_failure = !is.null(diagnostic_failure))
}


# --- Metrics and interval audit --------------------------------------------

.dp_sensitivity_empty_metrics <- function() {
  stats::setNames(rep(NA_real_, length(.DPPRIOR_SENSITIVITY_METRICS)),
                  .DPPRIOR_SENSITIVITY_METRICS)
}


.dp_sensitivity_interval_masses <- function(pmf, interval) {
  support <- seq_along(pmf)
  inside <- support >= interval[["lower", exact = TRUE]] &
    support <= interval[["upper", exact = TRUE]]
  c(
    coverage = sum(pmf[inside]),
    lower_tail = sum(pmf[support < interval[["lower", exact = TRUE]]]),
    upper_tail = sum(pmf[support > interval[["upper", exact = TRUE]]])
  )
}


.dp_sensitivity_interval_constraint_pass <- function(
    masses, interval, tolerance) {
  type <- interval[["type", exact = TRUE]]
  coverage <- interval[["coverage", exact = TRUE]]
  if (identical(type, "hard_bounds")) {
    abs(masses[["lower_tail"]]) <= tolerance &&
      abs(masses[["coverage"]] - 1) <= tolerance &&
      abs(masses[["upper_tail"]]) <= tolerance
  } else if (identical(type, "equal_tail")) {
    tail <- (1 - coverage) / 2
    abs(masses[["lower_tail"]] - tail) <= tolerance &&
      abs(masses[["coverage"]] - coverage) <= tolerance &&
      abs(masses[["upper_tail"]] - tail) <= tolerance
  } else {
    abs(masses[["coverage"]] - coverage) <= tolerance
  }
}


.dp_sensitivity_interval_audit <- function(evidence) {
  route <- evidence$input_provenance$target_route
  if (!identical(route, "interval")) {
    return(list(
      requested = NULL, selected = NULL, verification = NULL,
      status = NULL, source = NULL, usable = FALSE, verified = FALSE,
      reason = "not_interval_scenario"
    ))
  }
  interval <- evidence$target$used[["interval", exact = TRUE]]
  if (is.null(evidence$parameters)) {
    source <- if (identical(evidence$status, "infeasible")) {
      "target_infeasibility_certificate"
    } else {
      "target_construction_failure"
    }
    reason <- evidence$condition_evidence$calibration$message
    return(list(
      requested = interval, selected = NULL, verification = NULL,
      status = evidence$status, source = source,
      usable = FALSE, verified = identical(evidence$status, "infeasible"),
      reason = reason
    ))
  }
  selected <- .dp_sensitivity_interval_masses(
    evidence$selected_snapshot$K$pmf, interval
  )
  verifier <- .dp_sensitivity_interval_masses(
    evidence$verifier_snapshot$K$pmf, interval
  )
  selected <- c(
    as.list(selected),
    list(coverage_residual = unname(
      selected[["coverage"]] - interval[["coverage", exact = TRUE]]
    ))
  )
  verifier <- c(
    as.list(verifier),
    list(coverage_residual = unname(
      verifier[["coverage"]] - interval[["coverage", exact = TRUE]]
    ))
  )
  tolerance <- 1e-8
  stable <- all(abs(c(
    coverage = selected$coverage - verifier$coverage,
    lower_tail = selected$lower_tail - verifier$lower_tail,
    upper_tail = selected$upper_tail - verifier$upper_tail
  )) <= tolerance)
  passed <- .dp_sensitivity_interval_constraint_pass(
    unlist(selected[1:3], use.names = TRUE), interval, tolerance
  ) && .dp_sensitivity_interval_constraint_pass(
    unlist(verifier[1:3], use.names = TRUE), interval, tolerance
  ) && stable
  status <- if (passed) "converged" else "approximate"
  list(
    requested = interval, selected = selected,
    verification = c(verifier, list(
      tolerance = tolerance, passed = passed,
      source = "independent_interval_backcheck"
    )),
    status = status, source = "wrapper_backcheck_selected",
    usable = passed, verified = passed,
    reason = if (passed) "" else
      "Selected and verifier interval constraints did not pass fixed tolerance"
  )
}


.dp_sensitivity_metric_values <- function(
    evidence, interval, diagnostics_requested) {
  metrics <- .dp_sensitivity_empty_metrics()
  parameters <- evidence$parameters
  if (!is.null(parameters)) {
    a <- parameters[["a", exact = TRUE]]
    b <- parameters[["b", exact = TRUE]]
    K <- evidence$selected_snapshot$K
    metrics[c(
      "a", "b", "E_alpha", "CV_alpha", "E_K_J", "Var_K_J", "CV_K_J"
    )] <- c(
      a, b, a / b, 1 / sqrt(a), K$mean, K$variance,
      sqrt(K$variance) / K$mean
    )
    if (isTRUE(diagnostics_requested)) {
      wmax_50 <- unclass(wmax_tail_bounds(0.5, a = a, b = b))
      wmax_90 <- unclass(wmax_tail_bounds(0.9, a = a, b = b))
      metrics[c(
        "E_W_SB", "P_W_SB_gt_50", "P_W_SB_gt_90",
        "P_W_max_gt_50_lower_bound", "P_W_max_gt_50_upper_bound",
        "P_W_max_gt_90_lower_bound", "P_W_max_gt_90_upper_bound", "E_rho"
      )] <- c(
        as.numeric(mean_w1(a, b, evidence$evaluator$M_selected)),
        as.numeric(.diagnostic_wsb_tail(0.5, a, b)),
        as.numeric(.diagnostic_wsb_tail(0.9, a, b)),
        wmax_50[["lower_bound", exact = TRUE]],
        wmax_50[["upper_bound", exact = TRUE]],
        wmax_90[["lower_bound", exact = TRUE]],
        wmax_90[["upper_bound", exact = TRUE]],
        as.numeric(mean_rho(a, b, evidence$evaluator$M_selected))
      )
    }
  }
  if (!is.null(interval$status)) {
    metrics[["interval_requested"]] <- interval$requested$coverage
    if (!is.null(interval$selected)) {
      metrics[c(
        "interval_achieved", "interval_residual",
        "interval_left_tail", "interval_right_tail"
      )] <- c(
        interval$selected$coverage,
        interval$selected$coverage_residual,
        interval$selected$lower_tail,
        interval$selected$upper_tail
      )
    }
  }
  metrics
}


.dp_sensitivity_metric_reasons <- function(
    metrics, evidence, interval, diagnostics_requested) {
  reasons <- stats::setNames(rep(NA_character_, length(metrics)), names(metrics))
  unavailable <- is.na(metrics) & !is.nan(metrics)
  for (metric in names(metrics)[unavailable]) {
    reasons[[metric]] <- if (metric %in% c(
      "P_W_max_gt_50", "P_W_max_gt_90"
    )) {
      evidence$evaluator$W_max_point_policy
    } else if (metric %in% c(
      "interval_requested", "interval_achieved", "interval_residual",
      "interval_left_tail", "interval_right_tail"
    )) {
      if (is.null(interval$status)) "not_interval_scenario" else
        "interval_pmf_unavailable"
    } else if (metric %in% c(
      "E_W_SB", "P_W_SB_gt_50", "P_W_SB_gt_90",
      "P_W_max_gt_50_lower_bound", "P_W_max_gt_50_upper_bound",
      "P_W_max_gt_90_lower_bound", "P_W_max_gt_90_upper_bound", "E_rho"
    ) && !isTRUE(diagnostics_requested)) {
      "diagnostics_not_requested"
    } else {
      code <- evidence$condition_evidence$calibration$code
      if (is.null(code)) "fit_parameters_unavailable" else code
    }
  }
  reasons
}


.dp_sensitivity_worse_status <- function(...) {
  statuses <- unlist(list(...), use.names = FALSE)
  statuses <- statuses[!is.na(statuses) & nzchar(statuses)]
  rank <- c(
    converged = 1L, boundary = 2L, approximate = 3L,
    infeasible = 4L, failed = 5L
  )
  statuses[[which.max(unname(rank[statuses]))]]
}


# --- Canonical sensitivity reconciliation ----------------------------------

.dp_sensitivity_empty_local <- function() {
  data.frame(
    scenario_key = character(), axis = character(), axis_value = numeric(),
    settings_key = character(), lower_scenario_key = character(),
    upper_scenario_key = character(), lower_value = numeric(),
    upper_value = numeric(), metric = character(), component = character(),
    derivative = numeric(), method = character(), reason = character(),
    stringsAsFactors = FALSE
  )
}


.dp_sensitivity_reconciliation_check <- function(value) {
  .dpprior_new_check(
    value = value,
    reference = stats::setNames(rep(TRUE, length(value)), names(value)),
    tolerance = NULL, operator = "identical",
    source = "sensitivity_reconciliation"
  )
}


.dp_sensitivity_truth_check <- function(value) {
  .dpprior_new_check(
    value = value, reference = TRUE, tolerance = NULL,
    operator = "identical", source = "sensitivity_reconciliation"
  )
}


.dp_sensitivity_computation <- function() {
  settings <- list(
    method = "elicitation_sensitivity", controls = list(),
    parameterization = "none"
  )
  .dpprior_new_computation(
    request = settings, used = settings,
    orders = .dpprior_new_orders(
      M_requested = NULL, M_selected = NULL,
      M_verification_required = NULL, M_verification_used = NULL,
      requested_reason = "not_applicable",
      selected_reason = "not_applicable",
      verification_required_reason = "not_applicable",
      verification_used_reason = "not_applicable"
    ),
    scaling = .dpprior_new_scaling(),
    attempts = list(), candidate_evaluations = list(),
    selected_candidate_id = NULL, selected_attempt_id = NULL,
    fallback = .dpprior_new_fallback(),
    termination = .dpprior_new_termination(
      code = "deterministic", source = "sensitivity_reconciliation",
      iterations = NULL
    ),
    resources = list()
  )
}


.dp_sensitivity_provenance <- function(status) {
  approximate <- identical(status, "approximate")
  source_commit <- getOption("DPprior.source_commit", NULL)
  if (!.dp_sensitivity_plain_scalar_character(source_commit)) {
    source_commit <- NULL
  }
  .dpprior_new_provenance(
    requested_method = "elicitation_sensitivity",
    selected_method = "elicitation_sensitivity", is_fallback = FALSE,
    approximation = list(
      active = approximate, opt_in = FALSE,
      kind = if (approximate) "scenario_evidence_approximate" else NULL,
      warning_code = if (approximate) {
        "sensitivity_contains_approximate_scenarios"
      } else {
        NULL
      }
    ),
    projection = list(
      applied = FALSE, opt_in = FALSE, policy = NULL, record = NULL
    ),
    parameterization = "none",
    backend = list(
      package = "DPprior",
      package_version = tryCatch(
        as.character(utils::packageVersion("DPprior")),
        error = function(condition) "development"
      ),
      implementation =
        "R/20_elicitation_sensitivity.R:.dpprior_run_elicitation_sensitivity",
      source_commit = source_commit
    ),
    input_fit = NULL,
    migration = list(
      source_schema = "native", adapter = "none", lossless = TRUE,
      missing_evidence = character(), warnings = character()
    ),
    legacy = list(
      active = FALSE, contract = NULL, deprecation_stage = NULL
    )
  )
}


.dp_sensitivity_build_result <- function(J, records, executions) {
  keys <- vapply(records, `[[`, character(1L), "scenario_key")
  scenario_rows <- vector("list", length(records))
  result_rows <- vector("list", length(records))
  evidence <- conditions <- interval_checks <-
    stats::setNames(vector("list", length(records)), keys)
  metric_rows <- vector("list", length(records))

  for (index in seq_along(records)) {
    record <- records[[index]]
    key <- record$scenario_key
    execution <- executions[[index]]
    one_evidence <- execution$evidence
    interval <- .dp_sensitivity_interval_audit(one_evidence)
    row_status <- .dp_sensitivity_worse_status(
      one_evidence$status,
      if (!is.null(interval$status)) interval$status else character(),
      if (isTRUE(execution$row_failure)) "failed" else character()
    )
    row_usable <- if (row_status %in% c("converged", "boundary")) {
      isTRUE(one_evidence$usable) &&
        (is.null(interval$status) || isTRUE(interval$usable))
    } else if (identical(row_status, "approximate")) {
      isTRUE(one_evidence$usable) &&
        (is.null(interval$status) || isTRUE(interval$usable))
    } else {
      FALSE
    }
    row_verified <- if (row_status %in% c("converged", "boundary")) {
      isTRUE(one_evidence$verified) &&
        (is.null(interval$status) || isTRUE(interval$verified))
    } else {
      identical(row_status, "infeasible") && isTRUE(one_evidence$verified)
    }
    metrics <- .dp_sensitivity_metric_values(
      one_evidence, interval, record$diagnostics_requested
    )
    reasons <- .dp_sensitivity_metric_reasons(
      metrics, one_evidence, interval, record$diagnostics_requested
    )

    scenario_rows[[index]] <- data.frame(
      scenario_key = key, scenario_label = NA_character_,
      canonical_content = record$canonical_content,
      base_key = record$base_key,
      diagnostics_requested = record$diagnostics_requested,
      method_explicit = one_evidence$input_provenance$method_explicit,
      confidence_explicit = one_evidence$input_provenance$confidence_explicit,
      effective_method = one_evidence$input_provenance$requested_method,
      effective_confidence = if (identical(
        one_evidence$input_provenance$target_route,
        "qualitative_confidence"
      )) one_evidence$request$confidence else NA_character_,
      stringsAsFactors = FALSE
    )
    row <- data.frame(
      scenario_key = key, status = row_status,
      usable = row_usable, verified = row_verified,
      stringsAsFactors = FALSE
    )
    for (metric in .DPPRIOR_SENSITIVITY_METRICS) {
      row[[metric]] <- unname(metrics[[metric]])
    }
    result_rows[[index]] <- row
    metric_rows[[index]] <- data.frame(
      scenario_key = rep(key, length(.DPPRIOR_SENSITIVITY_METRICS)),
      metric = .DPPRIOR_SENSITIVITY_METRICS,
      value = unname(metrics), reason = unname(reasons),
      source = unname(.DPPRIOR_SENSITIVITY_METRIC_SOURCE),
      component = unname(.DPPRIOR_SENSITIVITY_METRIC_COMPONENT),
      status = rep(row_status, length(.DPPRIOR_SENSITIVITY_METRICS)),
      usable = rep(row_usable, length(.DPPRIOR_SENSITIVITY_METRICS)),
      verified = rep(row_verified, length(.DPPRIOR_SENSITIVITY_METRICS)),
      stringsAsFactors = FALSE
    )
    evidence[[key]] <- one_evidence
    conditions[[key]] <- one_evidence$condition_evidence
    interval_checks[[key]] <- interval
  }

  scenarios <- do.call(rbind, scenario_rows)
  scenario_results <- do.call(rbind, result_rows)
  metrics_long <- do.call(rbind, metric_rows)
  rownames(scenarios) <- rownames(scenario_results) <-
    rownames(metrics_long) <- NULL
  row_status <- scenario_results$status
  row_usable <- scenario_results$usable
  row_verified <- scenario_results$verified
  precedence <- c(
    converged = 1L, boundary = 2L, approximate = 3L,
    infeasible = 4L, failed = 5L
  )
  status <- row_status[[which.max(unname(precedence[row_status]))]]
  verified <- status %in% c("converged", "boundary", "infeasible") &&
    all(row_verified)
  usable <- status %in% c("converged", "boundary") && all(row_usable)
  global <- list(
    scenario_count = as.integer(length(keys)),
    converged_count = as.integer(sum(row_status %in% c(
      "converged", "boundary"
    ))),
    failed_count = as.integer(sum(row_status %in% c(
      "failed", "infeasible"
    ))),
    metric_count = as.integer(nrow(metrics_long))
  )
  sensitivity <- list(
    scenarios = scenarios, scenario_results = scenario_results,
    fit_evidence = evidence, conditions = conditions,
    interval_checks = interval_checks, metrics_long = metrics_long,
    local = .dp_sensitivity_empty_local(), global = global,
    metadata = list(J = J)
  )
  achieved <- list(scenario_count = as.integer(length(keys)))
  residuals <- list(
    missing_metric_count = as.integer(sum(is.na(metrics_long$value)))
  )
  tolerances <- list(
    expected_metric_count = as.integer(
      length(keys) * length(.DPPRIOR_SENSITIVITY_METRICS)
    )
  )
  selected_snapshot <- .dpprior_new_snapshot(
    parameters = NULL, M = NULL, achieved = achieved,
    residuals = residuals, tolerances = tolerances, finite = TRUE,
    source = "sensitivity_tables"
  )
  verifier_snapshot <- .dpprior_new_snapshot(
    parameters = NULL, M = NULL, achieved = achieved,
    residuals = residuals, tolerances = tolerances, finite = TRUE,
    source = "sensitivity_reconciliation"
  )
  truth <- c(
    scenario_key_identity = identical(scenario_results$scenario_key, keys),
    condition_key_identity = identical(names(conditions), keys),
    interval_key_identity = identical(names(interval_checks), keys),
    metric_grid_identity = identical(
      metrics_long$scenario_key,
      rep(keys, each = length(.DPPRIOR_SENSITIVITY_METRICS))
    ) && identical(
      metrics_long$metric,
      rep(.DPPRIOR_SENSITIVITY_METRICS, times = length(keys))
    ),
    row_status_identity = all(vapply(seq_along(keys), function(index) {
      rows <- metrics_long$scenario_key == keys[[index]]
      all(metrics_long$status[rows] == row_status[[index]]) &&
        all(metrics_long$usable[rows] == row_usable[[index]]) &&
        all(metrics_long$verified[rows] == row_verified[[index]])
    }, logical(1L))),
    top_status_identity = identical(
      status, row_status[[which.max(unname(precedence[row_status]))]]
    ),
    global_summary_identity = identical(
      global$scenario_count, as.integer(length(keys))
    ) && identical(global$metric_count, as.integer(nrow(metrics_long)))
  )
  invariants <- c(
    unique_keys = !anyDuplicated(keys),
    lexicographic_order = identical(keys, sort(keys, method = "radix")),
    finite_or_reason = all(
      is.finite(metrics_long$value) |
        (is.na(metrics_long$value) & !is.nan(metrics_long$value) &
           !is.na(metrics_long$reason) & nzchar(metrics_long$reason))
    ),
    failure_preservation = all(vapply(
      keys[row_status %in% c("failed", "infeasible")], function(key) {
        rows <- metrics_long$scenario_key == key
        key %in% names(conditions) && any(
          !is.na(metrics_long$reason[rows]) &
            nzchar(metrics_long$reason[rows])
        )
      }, logical(1L)
    ))
  )
  verification <- .dpprior_new_verification(
    method = "reconciliation", performed = TRUE, passed = verified,
    reason = if (verified) {
      "all scenario keys, conditions, metrics, and summaries reconciled"
    } else {
      "one or more retained scenario outcomes are not verified"
    },
    settings = list(
      scenario_count = as.integer(length(keys)), scenario_keys = keys,
      metric_names = .DPPRIOR_SENSITIVITY_METRICS
    ),
    selected_snapshot = selected_snapshot,
    verifier_snapshot = verifier_snapshot, stability = NULL,
    components = list(
      reconciliation = .dp_sensitivity_reconciliation_check(truth)
    ),
    invariants = stats::setNames(lapply(
      names(invariants), function(name) {
        .dp_sensitivity_truth_check(unname(invariants[[name]]))
      }
    ), names(invariants))
  )
  output <- .dpprior_new_sensitivity(
    method = "elicitation_sensitivity", J = J, status = status,
    usable = usable, verified = verified,
    message = if (verified) {
      "Canonical elicitation sensitivity completed"
    } else {
      "Canonical elicitation sensitivity retained non-verified outcomes"
    },
    target = list(defaults = list(J = J)),
    achieved = achieved, residuals = residuals, tolerances = tolerances,
    computation = .dp_sensitivity_computation(),
    verification = verification,
    provenance = .dp_sensitivity_provenance(status),
    sensitivity = sensitivity,
    compatibility = .dpprior_new_compatibility()
  )
  .dpprior_validate_result_v1(output)
  output
}


#' Run a canonical elicitation-sensitivity grid
#'
#' This internal API normalizes all scenarios before the first calibration,
#' evaluates actual `DPprior_fit()` routes, and returns one canonical
#' `dpprior.result/1` sensitivity object.  It describes prior-input
#' sensitivity only; posterior robustness requires a simulation study.
#'
#' @param J Common positive integer sample size.
#' @param scenarios Non-empty list of named scenario lists, or a data frame.
#' @param method Optional common A2-MN, A1, or A2-KL method.
#' @param M Common selected quadrature order in 10:256.
#' @param check_diagnostics Whether fresh diagnostic metrics are requested.
#' @param weight_thresholds Fixed c(0.5, 0.9) output thresholds.
#' @param .fit_fun Internal canonical fit adapter.
#' @param .key_fun Internal deterministic scenario-key adapter.
#' @return A canonical dpprior.result/1 sensitivity object.
#' @keywords internal
.dpprior_run_elicitation_sensitivity <- function(
    J, scenarios, method = NULL, M = .QUAD_NODES_DEFAULT,
    check_diagnostics = TRUE, weight_thresholds = c(0.5, 0.9),
    .fit_fun = DPprior_fit, .key_fun = .dp_sensitivity_content_key) {
  J <- .dp_sensitivity_normalize_count(
    J, "J", NULL, 1L, .MAX_J_DEFAULT
  )
  M <- .dp_sensitivity_normalize_count(M, "M", NULL, 10L, 256L)
  if (!is.null(method) && !.dp_sensitivity_plain_scalar_character(
    method, c("A2-MN", "A1", "A2-KL")
  )) {
    .dp_sensitivity_preflight_abort(
      "method must be NULL, A2-MN, A1, or A2-KL",
      "sensitivity_method"
    )
  }
  if (!.dp_sensitivity_plain_scalar_logical(check_diagnostics)) {
    .dp_sensitivity_preflight_abort(
      "check_diagnostics must be one ordinary non-missing logical",
      "sensitivity_diagnostics_control"
    )
  }
  if (!.dp_sensitivity_plain_numeric(weight_thresholds) ||
      length(weight_thresholds) != 2L || anyNA(weight_thresholds) ||
      any(!is.finite(weight_thresholds)) ||
      !setequal(as.numeric(weight_thresholds), c(0.5, 0.9))) {
    .dp_sensitivity_preflight_abort(
      "weight_thresholds must be exactly c(0.5, 0.9)",
      "weight_threshold_schema"
    )
  }
  if (!is.function(.key_fun)) {
    .dp_sensitivity_abort(
      ".key_fun must be a function",
      c("dpprior_sensitivity_adapter_error", "dpprior_type_error"),
      "adapter_type"
    )
  }

  # This completes for every scenario before .fit_fun can be called.
  records <- .dp_sensitivity_prepare_scenarios(
    J, scenarios, method, M, isTRUE(check_diagnostics), .key_fun
  )
  .dp_sensitivity_validate_fit_adapter(.fit_fun, records)
  executions <- lapply(records, .dp_sensitivity_execute_one,
                       fit_fun = .fit_fun)
  .dp_sensitivity_build_result(J, records, executions)
}
