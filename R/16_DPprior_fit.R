# =============================================================================
# Module 16: canonical DPprior_fit() wrapper
# =============================================================================
#
# DPprior_fit() is intentionally a thin producer adapter.  R/10--R/12 own the
# numerical candidate, optimization, fallback, status, and verification
# evidence.  This module owns only public target normalization, dispatch,
# exact cross-binding, optional fit-attached diagnostics, and compatibility
# aliases.  Scientific consumers must read the canonical dpprior.result/1
# spine; compatibility views are never an input to a decision below.


#' Convert a qualitative confidence level to a variance inflation factor
#'
#' @param confidence One of low, medium, or high.
#' @return The corresponding variance inflation factor.
#' @keywords internal
confidence_to_vif_fit <- function(confidence) {
  choices <- c("low", "medium", "high")
  confidence <- .dpprior_fit_match_arg(
    confidence, choices, "confidence", "dpprior_confidence_error"
  )
  switch(confidence, low = 5, medium = 2.5, high = 1.5)
}


#' Convert a variance inflation factor to a K variance
#'
#' @param mu_K Target mean.
#' @param vif Variance inflation factor.
#' @return vif times mu_K minus one.
#' @keywords internal
vif_to_variance_fit <- function(mu_K, vif) {
  if (!.dpprior_is_plain_numeric(mu_K) || length(mu_K) != 1L ||
      !is.finite(mu_K)) {
    .dpprior_abort_legacy_numeric(
      "mu_K must be a finite numeric scalar", mu_K, "mu_K",
      "finite numeric scalar", "dpprior_moment_target_error", scalar = TRUE
    )
  }
  if (!.dpprior_is_plain_numeric(vif) || length(vif) != 1L ||
      !is.finite(vif)) {
    .dpprior_abort_legacy_numeric(
      "vif must be a finite numeric scalar", vif, "vif",
      "finite numeric scalar", "dpprior_moment_target_error", scalar = TRUE
    )
  }
  if (vif <= 1) {
    warning(
      "VIF should be > 1 for overdispersion; using VIF = 1.01",
      call. = FALSE
    )
    vif <- 1.01
  }
  vif * (mu_K - 1)
}


.dpprior_fit_match_arg <- function(value, choices, name,
                                   subclass = "dpprior_choice_error") {
  ordinary <- is.character(value) && !is.object(value) &&
    is.null(dim(value)) && !anyNA(value)
  allowed_length <- length(value) == 1L || identical(value, choices)
  if (!ordinary || !allowed_length) {
    .dpprior_abort_invalid(
      sprintf("%s must be one ordinary character scalar", name),
      c(subclass, "dpprior_choice_error", "dpprior_type_error"),
      name, value, paste(choices, collapse = ", "), "choice"
    )
  }
  tryCatch(
    match.arg(value, choices),
    error = function(error) {
      .dpprior_abort_invalid(
        conditionMessage(error),
        c(subclass, "dpprior_choice_error"),
        name, value, paste(choices, collapse = ", "), "choice"
      )
    }
  )
}


.dpprior_fit_abort <- function(message, code, subclass,
                               result = NULL, action = NULL,
                               field = NULL, actual = NULL,
                               expected = NULL, cause = NULL) {
  stop(.dpprior_new_condition(
    message = message,
    classes = c(
      subclass, "dpprior_fit_error", "dpprior_calibration_error",
      "dpprior_error", "error"
    ),
    code = code,
    action = action,
    field = field,
    actual = actual,
    expected = expected,
    cause = cause,
    result = result
  ))
}


.dpprior_fit_exact_raw <- function(x) {
  if (is.object(x)) unclass(x) else x
}


.dpprior_fit_path_exists <- function(x, path) {
  if (!is.character(path) || is.object(path) || !is.null(dim(path)) ||
      length(path) != 1L || is.na(path) || !nzchar(path)) {
    return(FALSE)
  }
  value <- x
  for (part in strsplit(path, ".", fixed = TRUE)[[1L]]) {
    if (typeof(value) != "list" || !is.list(value)) return(FALSE)
    value <- .dpprior_fit_exact_raw(value)
    value_names <- names(value)
    if (is.null(value_names) || anyDuplicated(value_names) ||
        !(part %in% value_names)) {
      return(FALSE)
    }
    value <- value[[part, exact = TRUE]]
  }
  TRUE
}


.dpprior_fit_validate_warning_policy <- function(policy) {
  if (is.null(policy)) return(NULL)
  canonical <- c(
    "estimand", "direction", "weight_threshold", "action_threshold"
  )
  legacy <- c("estimand", "threshold", "direction", "action_threshold")
  ordinary <- typeof(policy) == "list" && is.list(policy) &&
    !is.object(policy) && !isS4(policy) && is.null(dim(policy)) &&
    identical(names(attributes(policy)), "names") &&
    !anyDuplicated(names(policy))
  canonical_record <- ordinary && identical(names(policy), canonical)
  legacy_record <- ordinary &&
    identical(sort(names(policy)), sort(legacy))
  if (!canonical_record && !legacy_record) {
    .dpprior_abort_invalid(
      paste(
        "warning_policy must be an ordinary canonical list with exact fields",
        "estimand, direction, weight_threshold, and action_threshold, or",
        "the legacy threshold-field record"
      ),
      c("dpprior_warning_policy_error", "dpprior_type_error"),
      "warning_policy", policy,
      paste(canonical, collapse = ", "),
      "invalid_policy_record"
    )
  }
  validate_choice <- function(field, choices) {
    value <- policy[[field, exact = TRUE]]
    valid <- is.character(value) && !is.object(value) && is.null(dim(value)) &&
      length(value) == 1L && !is.na(value) && value %in% choices
    if (!valid) {
      .dpprior_abort_invalid(
        sprintf(
          "warning_policy$%s must be exactly one of: %s",
          field, paste(choices, collapse = ", ")
        ),
        c(
          "dpprior_warning_policy_error", "dpprior_choice_error",
          "dpprior_type_error"
        ),
        paste0("warning_policy$", field), value,
        paste(choices, collapse = ", "), "choice"
      )
    }
    as.character(value)
  }
  estimand <- validate_choice("estimand", c("W_SB", "W_max"))
  direction <- validate_choice("direction", c("above", "below"))
  threshold_field <- if (canonical_record) "weight_threshold" else "threshold"
  weight_threshold <- .dpprior_validate_probability(
    policy[[threshold_field, exact = TRUE]],
    paste0("warning_policy$", threshold_field), scalar = TRUE, open = TRUE
  )
  action_threshold <- .dpprior_validate_probability(
    policy[["action_threshold", exact = TRUE]],
    "warning_policy$action_threshold", scalar = TRUE, open = FALSE
  )
  list(
    estimand = estimand,
    direction = direction,
    weight_threshold = weight_threshold,
    action_threshold = action_threshold
  )
}


.dpprior_fit_target_route <- function(target) {
  raw <- .dpprior_fit_exact_raw(target)
  kind <- raw[["kind", exact = TRUE]]
  request_names <- names(raw[["request", exact = TRUE]])
  if (identical(kind, "pmf")) return("strict_pmf")
  if (identical(kind, "interval")) return("interval")
  if (identical(kind, "cv")) return("coefficient_of_variation")
  if ("confidence" %in% request_names) return("qualitative_confidence")
  "direct_variance"
}


.dpprior_fit_target_implied <- function(target) {
  raw <- .dpprior_fit_exact_raw(target)
  raw[["implied", exact = TRUE]]
}


.dpprior_fit_require_target <- function(target, J) {
  validated <- .dpprior_require_schema(
    target, kind = "target", allow_legacy = FALSE
  )
  raw <- .dpprior_fit_exact_raw(validated)
  if (!identical(raw[["J", exact = TRUE]], as.integer(J))) {
    .dpprior_fit_abort(
      "target_K$J must be exactly identical to the DPprior_fit J argument.",
      "target_J_mismatch", "dpprior_target_integrity_error",
      result = validated, action = "supply_one_J_authority",
      field = "target_K$J", actual = raw[["J", exact = TRUE]],
      expected = as.integer(J)
    )
  }
  if (!isTRUE(raw[["usable", exact = TRUE]])) {
    stop(.dp_target_K_unusable_condition(validated))
  }
  validated
}


.dpprior_fit_target_projection <- function(target) {
  raw <- .dpprior_fit_exact_raw(target)
  raw[["provenance", exact = TRUE]][["projection", exact = TRUE]]
}


.dpprior_fit_target_equal_objective <- function(backend, requested) {
  backend_raw <- .dpprior_fit_exact_raw(backend)
  requested_raw <- .dpprior_fit_exact_raw(requested)
  identical(backend_raw[["J", exact = TRUE]],
            requested_raw[["J", exact = TRUE]]) &&
    identical(backend_raw[["implied", exact = TRUE]],
              requested_raw[["implied", exact = TRUE]]) &&
    identical(backend_raw[["pmf", exact = TRUE]],
              requested_raw[["pmf", exact = TRUE]])
}


.dpprior_fit_bind_A2_KL_objective <- function(target, backend_target) {
  target_raw <- .dpprior_fit_exact_raw(target)
  if (!is.null(target_raw[["pmf", exact = TRUE]])) return(target)
  if (!identical(target_raw[["kind", exact = TRUE]], "moments")) {
    .dpprior_fit_abort(
      paste(
        "A2-KL without an authoritative PMF requires a canonical moments",
        "target; the cv target kind cannot carry the A2-KL objective PMF."
      ),
      "a2_kl_target_kind_unsupported", "dpprior_target_method_conflict",
      result = target, action = "use_A2_MN_or_supply_a_strict_PMF",
      field = "target_K$kind", actual = target_raw[["kind", exact = TRUE]],
      expected = "moments or a target with an authoritative PMF"
    )
  }
  backend_raw <- .dpprior_fit_exact_raw(backend_target)
  objective <- backend_raw[["derivation", exact = TRUE]][[
    "request_to_normalized", exact = TRUE
  ]][["evidence", exact = TRUE]][["A2_KL_objective", exact = TRUE]]
  if (is.null(objective)) {
    .dpprior_fit_abort(
      "A2-KL backend omitted its canonical objective-distribution evidence.",
      "backend_objective_evidence_missing", "dpprior_backend_contract_error",
      result = target, action = "refit_with_a_current_A2_KL_backend"
    )
  }
  derivation <- target_raw[["derivation", exact = TRUE]]
  request_step <- derivation[["request_to_normalized", exact = TRUE]]
  evidence <- request_step[["evidence", exact = TRUE]]
  evidence[["A2_KL_objective"]] <- objective
  request_step[["evidence"]] <- evidence
  derivation[["request_to_normalized"]] <- request_step
  target_raw[["derivation"]] <- derivation
  class(target_raw) <- class(target)
  .dpprior_validate_target_v1(target_raw)
  target_raw
}


.dpprior_fit_wrapper_resource <- function(
    target, target_argument, method_was_missing, requested_method,
    selected_dispatch, M, check_diagnostics, warning_policy,
    a1_projection, a1_projection_explicit) {
  target_raw <- .dpprior_fit_exact_raw(target)
  list(
    api = "DPprior_fit",
    contract = "canonical_dpprior.result/1_passthrough",
    target_argument = target_argument,
    target_route = .dpprior_fit_target_route(target),
    target_schema = "dpprior.target/1",
    target_kind = target_raw[["kind", exact = TRUE]],
    target_authority = "result.target.K",
    requested_method = requested_method,
    selected_dispatch = selected_dispatch,
    method_was_missing = method_was_missing,
    J = target_raw[["J", exact = TRUE]],
    M_requested = as.integer(M),
    check_diagnostics = check_diagnostics,
    warning_policy = warning_policy,
    a1_projection = list(
      explicit = a1_projection_explicit,
      policy = a1_projection
    )
  )
}


.dpprior_fit_rebuild <- function(fit, target_K, wrapper_resource,
                                 diagnostics = NULL) {
  fit <- .dpprior_require_schema(fit, kind = "fit", allow_legacy = FALSE)
  raw <- .dpprior_fit_exact_raw(fit)
  target_K <- .dpprior_fit_require_target(target_K, raw[["J", exact = TRUE]])

  compatibility <- raw[["compatibility", exact = TRUE]]
  aliases <- compatibility[["top_level_aliases", exact = TRUE]]
  alias_names <- names(aliases)
  if (is.null(alias_names)) alias_names <- character()

  target <- raw[["target", exact = TRUE]]
  target[["K"]] <- target_K
  computation <- raw[["computation", exact = TRUE]]
  resources <- computation[["resources", exact = TRUE]]
  if ("wrapper" %in% names(resources) &&
      !identical(resources[["wrapper", exact = TRUE]], wrapper_resource)) {
    .dpprior_fit_abort(
      "Backend computation resources already contain a conflicting wrapper record.",
      "backend_wrapper_resource_conflict", "dpprior_backend_contract_error",
      result = fit, action = "use_an_unwrapped_native_producer_result"
    )
  }
  resources[["wrapper"]] <- wrapper_resource
  computation[["resources"]] <- resources
  provenance <- raw[["provenance", exact = TRUE]]
  provenance[["projection"]] <- .dpprior_fit_target_projection(target_K)

  existing_diagnostics <- "diagnostics" %in% names(raw) &&
    !("diagnostics" %in% alias_names)
  if (is.null(diagnostics) && existing_diagnostics) {
    diagnostics <- raw[["diagnostics", exact = TRUE]]
  }
  extension <- list()
  if ("proxy" %in% names(raw)) {
    extension[["proxy"]] <- raw[["proxy", exact = TRUE]]
  }
  if (!is.null(diagnostics)) extension[["diagnostics"]] <- diagnostics

  prospective <- raw
  if (length(alias_names)) prospective[alias_names] <- NULL
  prospective[["target"]] <- target
  prospective[["computation"]] <- computation
  prospective[["provenance"]] <- provenance
  if (!is.null(diagnostics)) prospective[["diagnostics"]] <- diagnostics
  if (!is.null(diagnostics) && "diagnostics" %in% names(aliases)) {
    aliases <- aliases[names(aliases) != "diagnostics"]
  }
  if (length(aliases)) {
    aliases <- aliases[vapply(
      aliases,
      function(path) .dpprior_fit_path_exists(prospective, path),
      logical(1)
    )]
  }

  rebuilt <- .dpprior_new_fit(
    mode = raw[["mode", exact = TRUE]],
    method = raw[["method", exact = TRUE]],
    J = raw[["J", exact = TRUE]],
    status = raw[["status", exact = TRUE]],
    usable = raw[["usable", exact = TRUE]],
    verified = raw[["verified", exact = TRUE]],
    message = raw[["message", exact = TRUE]],
    parameters = raw[["parameters", exact = TRUE]],
    target = target,
    achieved = raw[["achieved", exact = TRUE]],
    residuals = raw[["residuals", exact = TRUE]],
    tolerances = raw[["tolerances", exact = TRUE]],
    computation = computation,
    verification = raw[["verification", exact = TRUE]],
    provenance = provenance,
    compatibility = .dpprior_new_compatibility(),
    extension = extension
  )
  rebuilt <- .dpprior_append_compatibility_v2(
    rebuilt,
    aliases = aliases,
    views = compatibility[["views", exact = TRUE]],
    deprecations = compatibility[["deprecations", exact = TRUE]]
  )
  .dpprior_validate_result_v1(rebuilt)
  rebuilt
}


.dpprior_fit_bind_backend <- function(fit, requested_target,
                                      requested_dispatch,
                                      wrapper_resource) {
  fit <- .dpprior_require_schema(fit, kind = "fit", allow_legacy = FALSE)
  raw <- .dpprior_fit_exact_raw(fit)
  expected_mode <- switch(
    requested_dispatch,
    A1 = "a1_proxy", `A2-MN` = "a2_moment", `A2-KL` = "a2_kl"
  )
  method_ok <- if (identical(requested_dispatch, "A2-MN")) {
    raw[["method", exact = TRUE]] %in% c("A2-MN", "A2-MN+NM")
  } else {
    identical(raw[["method", exact = TRUE]], requested_dispatch)
  }
  if (!identical(raw[["mode", exact = TRUE]], expected_mode) || !method_ok ||
      !identical(raw[["J", exact = TRUE]],
                 .dpprior_fit_exact_raw(requested_target)[["J", exact = TRUE]])) {
    .dpprior_fit_abort(
      "Calibration backend result does not match the requested mode, method, and J.",
      "backend_dispatch_mismatch", "dpprior_backend_contract_error",
      result = fit, action = "use_the_matching_current_native_producer",
      actual = list(
        mode = raw[["mode", exact = TRUE]],
        method = raw[["method", exact = TRUE]],
        J = raw[["J", exact = TRUE]]
      ),
      expected = list(
        mode = expected_mode, method = requested_dispatch,
        J = .dpprior_fit_exact_raw(requested_target)[["J", exact = TRUE]]
      )
    )
  }

  backend_target <- raw[["target", exact = TRUE]][["K", exact = TRUE]]
  backend_target <- .dpprior_fit_require_target(
    backend_target, raw[["J", exact = TRUE]]
  )
  backend_target_raw <- .dpprior_fit_exact_raw(backend_target)
  projection <- backend_target_raw[["provenance", exact = TRUE]][[
    "projection", exact = TRUE
  ]]

  retained_target <- requested_target
  if (identical(requested_dispatch, "A1") &&
      isTRUE(projection[["applied", exact = TRUE]])) {
    requested_raw <- .dpprior_fit_exact_raw(requested_target)
    projection_record <- projection[["record", exact = TRUE]]
    before <- if (is.null(projection_record)) NULL else {
      projection_record[["before", exact = TRUE]]
    }
    if (!identical(.dpprior_fit_target_route(requested_target),
                   "direct_variance") ||
        !identical(before, requested_raw[["used", exact = TRUE]]) ||
        !identical(projection[["policy", exact = TRUE]], "nearest") ||
        !isTRUE(projection[["opt_in", exact = TRUE]])) {
      .dpprior_fit_abort(
        paste(
          "A1 projection can be retained only when the canonical direct",
          "target is exactly the recorded pre-projection target."
        ),
        "a1_projection_target_mismatch", "dpprior_target_integrity_error",
        result = fit, action = "use_a_direct_variance_target_or_disable_projection"
      )
    }
    retained_target <- backend_target
  } else {
    if (!.dpprior_fit_target_equal_objective(
      backend_target, requested_target
    )) {
      .dpprior_fit_abort(
        "Backend and wrapper targets do not have exactly identical objectives.",
        "backend_target_objective_mismatch", "dpprior_target_integrity_error",
        result = fit, action = "refit_from_one_canonical_target_authority"
      )
    }
    if (identical(requested_dispatch, "A2-KL")) {
      retained_target <- .dpprior_fit_bind_A2_KL_objective(
        requested_target, backend_target
      )
    }
  }
  .dpprior_fit_rebuild(fit, retained_target, wrapper_resource)
}


.dpprior_fit_public_residual <- function(fit) {
  raw <- .dpprior_fit_exact_raw(fit)
  residuals <- raw[["residuals", exact = TRUE]]
  if (typeof(residuals) != "list" || !is.list(residuals)) return(NA_real_)
  residuals <- .dpprior_fit_exact_raw(residuals)

  distribution <- residuals[["distribution", exact = TRUE]]
  if (typeof(distribution) == "list" && is.list(distribution)) {
    distribution <- .dpprior_fit_exact_raw(distribution)
    kl <- distribution[["kl", exact = TRUE]]
    if (.dpprior_is_plain_numeric(kl) && length(kl) == 1L &&
        is.finite(kl)) {
      return(as.numeric(kl))
    }
  }

  K <- residuals[["K", exact = TRUE]]
  if (typeof(K) == "list" && is.list(K)) {
    K <- .dpprior_fit_exact_raw(K)
    components <- list(
      K[["mean", exact = TRUE]], K[["variance", exact = TRUE]]
    )
    finite <- vapply(
      components,
      function(value) {
        .dpprior_is_plain_numeric(value) && length(value) == 1L &&
          is.finite(value)
      },
      logical(1)
    )
    if (any(finite)) {
      values <- vapply(components[finite], as.numeric, numeric(1))
      return(sqrt(sum(values^2)))
    }
  }
  NA_real_
}


.dpprior_fit_stop_unusable_backend <- function(fit, requested_method) {
  fit <- .dpprior_require_schema(fit, kind = "fit", allow_legacy = FALSE)
  raw <- .dpprior_fit_exact_raw(fit)
  status <- raw[["status", exact = TRUE]]
  usable <- raw[["usable", exact = TRUE]]
  verified <- raw[["verified", exact = TRUE]]
  backend_usable <- if (status %in% c("converged", "boundary")) {
    isTRUE(usable) && isTRUE(verified)
  } else if (identical(status, "approximate")) {
    identical(requested_method, "A1") &&
      identical(raw[["mode", exact = TRUE]], "a1_proxy") &&
      identical(raw[["method", exact = TRUE]], "A1") &&
      isTRUE(usable) && !isTRUE(verified)
  } else {
    FALSE
  }
  if (backend_usable) return(fit)

  residual <- .dpprior_fit_public_residual(fit)
  stop(.dpprior_new_condition(
    message = sprintf(
      paste0(
        "Calibration did not produce a usable result for method %s ",
        "(status: %s, residual: %.4g). Inspect condition$result and ",
        "choose a supported target or solver policy."
      ),
      requested_method, status, residual
    ),
    classes = c(
      "dpprior_calibration_unusable", "dpprior_calibration_error",
      "dpprior_error", "error"
    ),
    code = "calibration_unusable",
    method = requested_method,
    status = status,
    residual = residual,
    result = fit,
    action = "inspect_condition_result_and_choose_supported_target_or_solver_policy"
  ))
}


.dpprior_fit_diagnostics_extension <- function(
    fit, M_selected = NULL, warning_policy = NULL,
    allow_approximate = FALSE) {
  raw <- .dpprior_fit_exact_raw(
    .dpprior_require_schema(fit, kind = "fit", allow_legacy = FALSE)
  )
  parameters <- raw[["parameters", exact = TRUE]]
  if (is.null(parameters)) {
    .dpprior_fit_abort(
      "Fit-attached diagnostics require a public canonical parameter candidate.",
      "diagnostic_candidate_unavailable", "dpprior_diagnostics_input_error",
      result = fit, action = "refit_to_obtain_a_public_candidate"
    )
  }
  orders <- raw[["computation", exact = TRUE]][["orders", exact = TRUE]]
  if (is.null(M_selected)) {
    M_selected <- if (identical(raw[["mode", exact = TRUE]], "a1_proxy")) {
      .QUAD_NODES_DEFAULT
    } else {
      orders[["M_selected", exact = TRUE]]
    }
  }
  M_selected <- as.integer(M_selected)
  M_required <- as.integer(max(2L * M_selected, M_selected + 40L))
  M_verification <- if (identical(raw[["mode", exact = TRUE]], "a1_proxy")) {
    M_required
  } else {
    orders[["M_verification_used", exact = TRUE]]
  }
  if (is.na(M_selected) || M_selected < 10L ||
      M_selected > .QUADRATURE_MAX_NODES ||
      is.null(M_verification) || M_required > .QUADRATURE_MAX_NODES ||
      M_verification > .QUADRATURE_MAX_NODES) {
    .dpprior_fit_abort(
      paste(
        "Fit-attached diagnostics require selected and independent",
        "verification orders within the supported quadrature ceiling."
      ),
      "diagnostic_order_unavailable", "dpprior_diagnostics_order_error",
      result = fit, action = "use_M_at_most_256_or_disable_diagnostics",
      actual = list(
        M_selected = M_selected, M_verification_required = M_required,
        M_verification_used = M_verification
      ), expected = sprintf("orders in [10,%d]", .QUADRATURE_MAX_NODES)
    )
  }
  M_verification <- as.integer(M_verification)
  authority <- list(
    method = "fresh_component_specific_diagnostics",
    M_selected = M_selected,
    M_verification_required = M_required,
    M_verification_used = M_verification,
    absolute_tolerance = 1e-10,
    relative_tolerance = 1e-8,
    pmf_mass_tolerance = .TOL_PMF_SUM,
    warning_policy = warning_policy,
    allow_approximate = allow_approximate
  )
  a <- parameters[["a", exact = TRUE]]
  b <- parameters[["b", exact = TRUE]]
  J <- raw[["J", exact = TRUE]]
  K_evidence <- .get_K_pmf_support(
    J, a, b, M = M_selected, M_verify = M_verification,
    abs_tol = authority[["absolute_tolerance", exact = TRUE]],
    rel_tol = authority[["relative_tolerance", exact = TRUE]]
  )
  selected_pmf <- unname(as.numeric(K_evidence[["pmf", exact = TRUE]]))
  verifier_pmf <- unname(as.numeric(
    K_evidence[["verification_pmf", exact = TRUE]]
  ))
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
    coclustering.mean = abs(
      selected_rho[["mean"]] - verifier_rho[["mean"]]
    ),
    coclustering.variance = abs(
      selected_rho[["variance"]] - verifier_rho[["variance"]]
    )
  )
  scalar_tolerance <- function(selected, verifier) {
    authority[["absolute_tolerance", exact = TRUE]] +
      authority[["relative_tolerance", exact = TRUE]] *
      max(abs(selected), abs(verifier), 1)
  }
  tolerance <- c(
    K.mean = scalar_tolerance(
      selected_K[["mean"]], verifier_K[["mean"]]
    ),
    K.variance = scalar_tolerance(
      selected_K[["variance"]], verifier_K[["variance"]]
    ),
    K.pmf_l1 = authority[["absolute_tolerance", exact = TRUE]] +
      authority[["relative_tolerance", exact = TRUE]],
    weights.mean = scalar_tolerance(selected_weight, verifier_weight),
    coclustering.mean = scalar_tolerance(
      selected_rho[["mean"]], verifier_rho[["mean"]]
    ),
    coclustering.variance = scalar_tolerance(
      selected_rho[["variance"]], verifier_rho[["variance"]]
    )
  )
  passed <- c(
    alpha = TRUE,
    K = all(delta[c("K.mean", "K.variance", "K.pmf_l1")] <=
              tolerance[c("K.mean", "K.variance", "K.pmf_l1")]) &&
      abs(sum(selected_pmf) - 1) <=
        authority[["pmf_mass_tolerance", exact = TRUE]] &&
      abs(sum(verifier_pmf) - 1) <=
        authority[["pmf_mass_tolerance", exact = TRUE]],
    weights = delta[["weights.mean"]] <= tolerance[["weights.mean"]],
    coclustering = all(delta[c(
      "coclustering.mean", "coclustering.variance"
    )] <= tolerance[c(
      "coclustering.mean", "coclustering.variance"
    )])
  )
  status <- ifelse(passed, "converged", "approximate")
  usable <- passed | (!passed & allow_approximate)

  policy_results <- if (is.null(warning_policy)) {
    list()
  } else if (identical(
    warning_policy[["estimand", exact = TRUE]], "W_SB"
  )) {
    value <- as.numeric(.diagnostic_wsb_tail(
      warning_policy[["weight_threshold", exact = TRUE]], a, b
    ))
    outcome <- if (identical(
      warning_policy[["direction", exact = TRUE]], "above"
    )) {
      if (value > warning_policy[["action_threshold", exact = TRUE]]) {
        "triggered"
      } else "not_triggered"
    } else if (value <
               warning_policy[["action_threshold", exact = TRUE]]) {
      "triggered"
    } else "not_triggered"
    list(list(
      estimand = "W_SB",
      direction = warning_policy[["direction", exact = TRUE]],
      threshold = warning_policy[["action_threshold", exact = TRUE]],
      value = value,
      lower = NULL,
      upper = NULL,
      outcome = outcome,
      basis = "exact_tail_probability"
    ))
  } else {
    list(list(
      estimand = "W_max",
      direction = warning_policy[["direction", exact = TRUE]],
      threshold = warning_policy[["action_threshold", exact = TRUE]],
      value = NULL,
      lower = NULL,
      upper = NULL,
      outcome = "indeterminate",
      basis = "backend_unavailable"
    ))
  }
  triggered <- length(policy_results) == 1L && identical(
    policy_results[[1L]][["outcome", exact = TRUE]], "triggered"
  )
  warning_message <- if (triggered) {
    sprintf(
      "Explicit W_SB policy triggered: P(W_SB %s %.6g) is %.6g against action threshold %.6g.",
      if (identical(
        warning_policy[["direction", exact = TRUE]], "above"
      )) ">" else "<",
      warning_policy[["weight_threshold", exact = TRUE]],
      policy_results[[1L]][["value", exact = TRUE]],
      warning_policy[["action_threshold", exact = TRUE]]
    )
  } else character()

  list(
    authority = authority,
    policy_results = policy_results,
    warnings = warning_message,
    alpha = list(
      status = unname(status[["alpha"]]),
      usable = unname(usable[["alpha"]]),
      verified = unname(passed[["alpha"]]),
      mean = a / b,
      CV = 1 / sqrt(a)
    ),
    K = list(
      status = unname(status[["K"]]),
      usable = unname(usable[["K"]]),
      verified = unname(passed[["K"]]),
      mean = unname(selected_K[["mean"]]),
      variance = unname(selected_K[["variance"]]),
      pmf = selected_pmf,
      M = M_selected
    ),
    weights = list(
      status = unname(status[["weights"]]),
      usable = unname(usable[["weights"]]),
      verified = unname(passed[["weights"]]),
      mean = selected_weight
    ),
    coclustering = list(
      status = unname(status[["coclustering"]]),
      usable = unname(usable[["coclustering"]]),
      verified = unname(passed[["coclustering"]]),
      mean = unname(selected_rho[["mean"]]),
      variance = unname(selected_rho[["variance"]])
    )
  )
}


.dpprior_fit_attach_diagnostics <- function(fit, M, warning_policy) {
  raw <- .dpprior_fit_exact_raw(
    .dpprior_require_schema(fit, kind = "fit", allow_legacy = FALSE)
  )
  M_selected <- if (identical(raw[["mode", exact = TRUE]], "a1_proxy")) {
    as.integer(M)
  } else {
    raw[["computation", exact = TRUE]][["orders", exact = TRUE]][[
      "M_selected", exact = TRUE
    ]]
  }
  diagnostics <- tryCatch(
    .dpprior_fit_diagnostics_extension(
      fit, M_selected = M_selected, warning_policy = warning_policy,
      allow_approximate = FALSE
    ),
    error = function(error) error
  )
  if (inherits(diagnostics, "condition")) {
    if (inherits(diagnostics, "dpprior_condition")) stop(diagnostics)
    .dpprior_fit_abort(
      paste("Fit-attached diagnostics failed:", conditionMessage(diagnostics)),
      "fit_diagnostics_computation_failed",
      "dpprior_diagnostics_computation_error",
      result = fit, action = "inspect_parameters_and_quadrature_orders",
      cause = list(
        class = class(diagnostics)[[1L]], message = conditionMessage(diagnostics)
      )
    )
  }
  wrapper_resource <- raw[["computation", exact = TRUE]][[
    "resources", exact = TRUE
  ]][["wrapper", exact = TRUE]]
  target_K <- raw[["target", exact = TRUE]][["K", exact = TRUE]]
  with_diagnostics <- .dpprior_fit_rebuild(
    fit, target_K, wrapper_resource, diagnostics = diagnostics
  )
  component_status <- vapply(
    c("alpha", "K", "weights", "coclustering"),
    function(name) diagnostics[[name, exact = TRUE]][["status", exact = TRUE]],
    character(1)
  )
  if (any(component_status == "approximate")) {
    .dpprior_fit_abort(
      paste(
        "Fit-attached diagnostics are approximate; the canonical fit is",
        "retained in condition$result for review."
      ),
      "fit_diagnostics_approximate",
      "dpprior_diagnostics_approximation_error",
      result = with_diagnostics,
      action = "increase_M_or_set_check_diagnostics_FALSE"
    )
  }
  with_diagnostics
}


.dpprior_fit_resignal_backend <- function(condition, requested_target,
                                          requested_dispatch,
                                          wrapper_resource) {
  condition_result <- tryCatch(
    condition[["result", exact = TRUE]], error = function(error) NULL
  )
  canonical_result <- if (is.null(condition_result)) NULL else tryCatch(
    .dpprior_require_schema(
      condition_result, kind = "fit", allow_legacy = FALSE
    ),
    error = function(error) NULL
  )
  if (!is.null(canonical_result)) {
    rebound <- tryCatch(
      .dpprior_fit_bind_backend(
        canonical_result, requested_target, requested_dispatch,
        wrapper_resource
      ),
      error = function(error) error
    )
    if (inherits(rebound, "condition")) {
      .dpprior_fit_abort(
        paste(
          "A backend condition carried a canonical result that failed",
          "wrapper cross-binding."
        ),
        "backend_condition_result_mismatch", "dpprior_backend_contract_error",
        result = requested_target,
        action = "refit_with_the_matching_current_backend",
        cause = list(
          class = class(rebound)[[1L]],
          code = rebound[["code", exact = TRUE]] %||% "unknown",
          message = conditionMessage(rebound)
        )
      )
    }
    condition_raw <- unclass(condition)
    condition_raw[["result"]] <- rebound
    class(condition_raw) <- class(condition)
    stop(condition_raw)
  }

  if (inherits(condition, "dpprior_condition")) {
    condition_raw <- unclass(condition)
    condition_raw[["result"]] <- requested_target
    if (is.null(condition_raw[["action", exact = TRUE]])) {
      condition_raw[["action"]] <- "revise_the_canonical_target_or_method"
    }
    class(condition_raw) <- class(condition)
    stop(condition_raw)
  }

  .dpprior_fit_abort(
    paste("Calibration backend failed:", conditionMessage(condition)),
    "backend_raw_error", "dpprior_backend_contract_error",
    result = requested_target,
    action = "inspect_the_typed_backend_contract",
    cause = list(
      class = class(condition)[[1L]], message = conditionMessage(condition)
    )
  )
}


#' Fit a Gamma hyperprior from an occupied-cluster target
#'
#' \code{DPprior_fit()} constructs one canonical \eqn{K_J} target, or accepts
#' one from \code{\link{DPprior_target_K}}, dispatches to A1, A2-MN, or A2-KL,
#' and returns the numerical producer's native canonical result without
#' flattening or reinterpreting its evidence.
#'
#' @param J Integer design size. When \code{target_K} is supplied, \code{J}
#'   must match the target exactly.
#' @param mu_K,var_K Optional direct target mean and variance.
#' @param confidence Optional \code{"low"}, \code{"medium"}, or
#'   \code{"high"} qualitative uncertainty route.
#' @param method One of \code{"A2-MN"}, \code{"A1"}, or \code{"A2-KL"}.
#' @param target_pmf Optional strict target PMF on \eqn{1,\ldots,J}.
#' @param check_diagnostics Whether to attach freshly verified canonical
#'   diagnostics.
#' @param warn_dominance Deprecated first-size-biased-weight warning adapter.
#' @param M Selected quadrature order.
#' @param verbose Whether the numerical producer prints progress.
#' @param warning_policy Optional ordinary list with exact canonical fields
#'   \code{estimand}, \code{direction}, \code{weight_threshold}, and
#'   \code{action_threshold}. The legacy \code{threshold} field is accepted
#'   here and normalized immediately for compatibility.
#' @param a1_projection A1 projection policy: \code{"error"} or an explicit
#'   opt-in to \code{"nearest"}.
#' @param cv_K Optional coefficient of variation \eqn{SD(K_J)/E(K_J)}.
#' @param K_interval Optional canonical interval-target specification.
#' @param target_K Optional validated \code{dpprior.target/1} K target. It
#'   cannot be combined with scalar target arguments.
#' @param ... Must be empty.
#'
#' @return A canonical \code{dpprior.result/1} \code{DPprior_fit} object.
#'   Use \code{parameters$a} and \code{parameters$b} for the Gamma shape and
#'   rate; \code{target$K}, \code{achieved}, and \code{residuals} for the
#'   scientific claim; and \code{computation}, \code{verification}, and
#'   \code{provenance} for its audit trail. A result is decision-ready only
#'   when its mode-specific contract and \code{status}, \code{usable}, and
#'   \code{verified} fields permit that use. If a backend signals a typed
#'   calibration condition, the canonical target or result is retained in
#'   \code{condition$result}.
#'
#' @details
#' Supply exactly one uncertainty route: direct \code{var_K},
#' \code{confidence}, \code{cv_K}, \code{K_interval}, or
#' \code{target_pmf}. A pre-built \code{target_K} is an alternative to all
#' scalar target arguments. Compatibility aliases are non-authoritative; new
#' code should use the canonical nested fields.
#'
#' @examples
#' target <- DPprior_target_K(J = 50, mu_K = 5, var_K = 8)
#' fit <- DPprior_fit(
#'   J = 50, target_K = target, method = "A2-MN",
#'   check_diagnostics = FALSE
#' )
#' fit$parameters
#' fit[c("status", "usable", "verified")]
#' @family elicitation
#' @export
DPprior_fit <- function(J, mu_K = NULL, var_K = NULL,
                        confidence = c("medium", "low", "high"),
                        method = c("A2-MN", "A1", "A2-KL"),
                        target_pmf = NULL,
                        check_diagnostics = TRUE,
                        warn_dominance = NULL,
                        M = .QUAD_NODES_DEFAULT,
                        verbose = FALSE,
                        warning_policy = NULL,
                        a1_projection = c("error", "nearest"),
                        cv_K = NULL, K_interval = NULL,
                        target_K = NULL, ...) {
  method_was_missing <- missing(method)
  confidence_explicit <- !missing(confidence)
  projection_explicit <- !missing(a1_projection)
  requested_method <- if (method_was_missing) NULL else method

  dots <- match.call(expand.dots = FALSE)[["...", exact = TRUE]]
  if (!is.null(dots) && length(dots)) {
    dot_names <- names(dots)
    if (is.null(dot_names)) dot_names <- rep("<unnamed>", length(dots))
    dot_names[!nzchar(dot_names)] <- "<unnamed>"
    .dpprior_abort_invalid(
      sprintf(
        "Unknown DPprior_fit argument(s) in ...: %s.",
        paste(unique(dot_names), collapse = ", ")
      ),
      c("dpprior_unknown_control_error", "dpprior_fit_control_error"),
      "...", dot_names, "no additional arguments", "unknown_control"
    )
  }

  assert_valid_J(J)
  J <- as.integer(J)
  check_diagnostics <- .dpprior_validate_control(
    check_diagnostics, "check_diagnostics", type = "logical"
  )
  verbose <- .dpprior_validate_control(verbose, "verbose", type = "logical")
  M <- .as_integer_scalar(M, "M", min = 10L)
  method <- .dpprior_fit_match_arg(
    method, c("A2-MN", "A1", "A2-KL"), "method",
    "dpprior_method_error"
  )
  a1_projection <- .dpprior_fit_match_arg(
    a1_projection, c("error", "nearest"), "a1_projection",
    "dpprior_a1_projection_policy_error"
  )

  if (!is.null(warn_dominance)) {
    warn_dominance <- .dpprior_validate_control(
      warn_dominance, "warn_dominance", type = "logical"
    )
    if (isTRUE(warn_dominance)) {
      if (!is.null(warning_policy)) {
        .dpprior_abort_invalid(
          "warning_policy conflicts with deprecated warn_dominance=TRUE",
          c("dpprior_warning_policy_error", "dpprior_conflicting_input"),
          "warning_policy/warn_dominance",
          list(warning_policy = warning_policy,
               warn_dominance = warn_dominance),
          "one warning policy", "conflicting_warning_policy"
        )
      }
      .dpprior_warn(
        paste(
          "warn_dominance is deprecated; translating only to the explicit",
          "legacy W_SB policy P(W_SB > 0.5) > 0.4."
        ),
        "dpprior_deprecated_warning", "warn_dominance", TRUE,
        "warning_policy with estimand='W_SB'", "deprecated_warn_dominance"
      )
      warning_policy <- list(
        estimand = "W_SB", direction = "above",
        weight_threshold = 0.5, action_threshold = 0.4
      )
    }
  }
  warning_policy <- .dpprior_fit_validate_warning_policy(warning_policy)
  if (!is.null(warning_policy) && !check_diagnostics) {
    .dpprior_abort_invalid(
      "warning_policy requires check_diagnostics = TRUE",
      "dpprior_warning_policy_error", "warning_policy", warning_policy,
      "check_diagnostics = TRUE", "diagnostics_disabled"
    )
  }

  target_argument <- if (is.null(target_K)) "scalar_arguments" else "target_K"
  if (!is.null(target_K)) {
    scalar_conflicts <- c(
      if (!is.null(mu_K)) "mu_K",
      if (!is.null(var_K)) "var_K",
      if (confidence_explicit) "confidence",
      if (!is.null(cv_K)) "cv_K",
      if (!is.null(K_interval)) "K_interval",
      if (!is.null(target_pmf)) "target_pmf"
    )
    if (length(scalar_conflicts)) {
      .dpprior_fit_abort(
        sprintf(
          "target_K cannot be combined with scalar target argument(s): %s.",
          paste(scalar_conflicts, collapse = ", ")
        ),
        "target_K_scalar_conflict", "dpprior_target_method_conflict",
        result = tryCatch(
          .dpprior_require_schema(target_K, kind = "target"),
          error = function(error) NULL
        ),
        action = "supply_either_target_K_or_scalar_target_arguments",
        actual = scalar_conflicts, expected = "one target authority"
      )
    }
    target_K <- .dpprior_fit_require_target(target_K, J)
  } else {
    source_names <- c(
      if (!is.null(var_K) && is.null(target_pmf)) "var_K",
      if (confidence_explicit) "confidence",
      if (!is.null(cv_K)) "cv_K",
      if (!is.null(K_interval)) "K_interval",
      if (!is.null(target_pmf)) "target_pmf"
    )
    if (length(source_names) == 0L) {
      confidence <- "medium"
      source_names <- "confidence"
    }
    if (length(source_names) != 1L) {
      .dpprior_abort_invalid(
        sprintf(
          "Supply exactly one uncertainty source; received: %s.",
          paste(source_names, collapse = ", ")
        ),
        c("dpprior_uncertainty_source_conflict", "dpprior_conflicting_input"),
        "uncertainty", source_names, "exactly one source",
        "multiple_uncertainty_sources"
      )
    }
    confidence_used <- if (identical(source_names, "confidence")) {
      .dpprior_fit_match_arg(
        confidence, c("medium", "low", "high"), "confidence",
        "dpprior_confidence_error"
      )
    } else NULL
    target_K <- .dp_target_K(
      J = J, mu_K = mu_K, var_K = var_K,
      confidence = confidence_used, cv_K = cv_K,
      K_interval = K_interval, target_pmf = target_pmf
    )
    target_K <- .dp_target_K_stop_unusable(target_K)
  }

  target_raw <- .dpprior_fit_exact_raw(target_K)
  target_kind <- target_raw[["kind", exact = TRUE]]
  has_target_pmf <- !is.null(target_raw[["pmf", exact = TRUE]])
  if (method_was_missing && has_target_pmf) method <- "A2-KL"
  if (has_target_pmf && !identical(method, "A2-KL")) {
    .dpprior_fit_abort(
      "A canonical PMF or interval target requires method = 'A2-KL'.",
      "target_pmf_wrong_method", "dpprior_target_method_conflict",
      result = target_K, action = "use_A2_KL",
      field = "method", actual = method, expected = "A2-KL"
    )
  }
  if (!has_target_pmf && identical(method, "A2-KL") &&
      !identical(target_kind, "moments")) {
    .dpprior_fit_abort(
      "A2-KL without a PMF requires a canonical moments target.",
      "a2_kl_target_kind_unsupported", "dpprior_target_method_conflict",
      result = target_K, action = "use_A2_MN_or_supply_a_strict_PMF",
      field = "target_K$kind", actual = target_kind, expected = "moments"
    )
  }
  if (projection_explicit && !identical(method, "A1")) {
    .dpprior_abort_invalid(
      "a1_projection can only be supplied with method = 'A1'",
      c("dpprior_a1_projection_policy_error", "dpprior_conflicting_input"),
      "a1_projection", a1_projection, "method = 'A1'",
      "a1_projection_wrong_method"
    )
  }
  implied <- target_raw[["implied", exact = TRUE]]
  implied_mean <- implied[["mean", exact = TRUE]]
  implied_variance <- implied[["variance", exact = TRUE]]
  if (identical(method, "A1") && identical(a1_projection, "nearest") &&
      implied_variance <= implied_mean - 1 &&
      !identical(.dpprior_fit_target_route(target_K), "direct_variance")) {
    .dpprior_fit_abort(
      paste(
        "A1 projection of a derived confidence/CV route cannot preserve both",
        "derivation steps in the canonical target; use a direct variance target."
      ),
      "a1_projection_derived_route_unsupported",
      "dpprior_a1_projection_policy_error",
      result = target_K, action = "construct_a_direct_variance_target"
    )
  }
  if (J == 1L) {
    .dpprior_fit_abort(
      paste(
        "K_1 is identically one for every positive Gamma prior, so a and b",
        "are not identified."
      ),
      "calibration_nonidentifiable_j1",
      "dpprior_calibration_nonidentifiable",
      result = target_K, action = "increase_J_or_specify_alpha_directly"
    )
  }

  wrapper_resource <- .dpprior_fit_wrapper_resource(
    target_K, target_argument, method_was_missing, requested_method,
    method, M, check_diagnostics, warning_policy,
    a1_projection, projection_explicit
  )

  backend <- tryCatch(
    switch(
      method,
      A1 = {
        if (verbose) message("Using A1 closed-form approximation")
        DPprior_a1(
          J, implied_mean, implied_variance, projection = a1_projection
        )
      },
      `A2-MN` = {
        if (verbose) message("Using A2-MN Newton refinement")
        DPprior_a2_newton(
          J, implied_mean, implied_variance, M = M, verbose = verbose
        )
      },
      `A2-KL` = {
        if (verbose) message("Using A2-KL divergence minimization")
        if (has_target_pmf) {
          DPprior_a2_kl(
            J, target = target_raw[["pmf", exact = TRUE]],
            method = "pmf", M = M, verbose = verbose
          )
        } else {
          DPprior_a2_kl(
            J, target = list(
              mu_K = implied_mean, var_K = implied_variance
            ),
            method = "chisq", M = M, verbose = verbose
          )
        }
      }
    ),
    error = function(error) error
  )
  if (inherits(backend, "condition")) {
    .dpprior_fit_resignal_backend(
      backend, target_K, method, wrapper_resource
    )
  }

  fit <- tryCatch(
    .dpprior_fit_bind_backend(backend, target_K, method, wrapper_resource),
    error = function(error) error
  )
  if (inherits(fit, "condition")) {
    if (inherits(fit, "dpprior_fit_error")) stop(fit)
    .dpprior_fit_abort(
      paste("Backend canonical cross-binding failed:", conditionMessage(fit)),
      "backend_schema_binding_failed", "dpprior_backend_contract_error",
      result = target_K, action = "use_a_current_canonical_backend",
      cause = list(
        class = class(fit)[[1L]],
        code = fit[["code", exact = TRUE]] %||% "unknown",
        message = conditionMessage(fit)
      )
    )
  }

  fit <- .dpprior_fit_stop_unusable_backend(fit, method)

  if (check_diagnostics) {
    fit <- .dpprior_fit_attach_diagnostics(fit, M, warning_policy)
  }
  .dpprior_validate_result_v1(fit)

  if (check_diagnostics) {
    diagnostics <- .dpprior_fit_exact_raw(fit)[["diagnostics", exact = TRUE]]
    warnings <- diagnostics[["warnings", exact = TRUE]]
    if (length(warnings)) {
      .dpprior_warn(
        warnings[[1L]],
        c("dpprior_diagnostic_policy_warning", "dpprior_estimand_warning"),
        "warning_policy", warning_policy,
        "explicit estimand-specific user policy", "policy_triggered"
      )
    }
  }
  fit
}


#' Verify the canonical DPprior_fit wrapper
#'
#' @param verbose Whether to print a success message.
#' @return Invisibly TRUE.
#' @keywords internal
verify_DPprior_fit <- function(verbose = TRUE) {
  verbose <- .dpprior_validate_control(verbose, "verbose", type = "logical")
  fit <- DPprior_fit(
    J = 20L, mu_K = 4, var_K = 8, method = "A1",
    check_diagnostics = FALSE, verbose = FALSE
  )
  .dpprior_validate_result_v1(fit)
  if (verbose) message("DPprior_fit canonical wrapper verification passed")
  invisible(TRUE)
}
