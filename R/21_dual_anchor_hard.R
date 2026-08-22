# =============================================================================
# Module 21: Dual-anchor v2 shared primitives and hard inequality calibration
# =============================================================================

.DPPRIOR_V2_STATUS_CODES <- c(
  "converged", "boundary", "approximate", "infeasible", "failed"
)


.dpprior_v2_abort_invalid <- function(message, argument = NULL, value = NULL,
                                      expected = NULL, code = "invalid",
                                      subclass = character(), call = NULL) {
  .dpprior_abort_invalid(
    message = message,
    subclass = c(subclass, "dpprior_dual_input_error"),
    argument = argument,
    value = value,
    expected = expected,
    code = code,
    call = call
  )
}


.dpprior_v2_plain_scalar <- function(x, mode = c("numeric", "character",
                                                 "logical")) {
  mode <- match.arg(mode)
  predicate <- switch(
    mode,
    numeric = is.numeric,
    character = is.character,
    logical = is.logical
  )
  predicate(x) && is.null(dim(x)) && !is.object(x) && length(x) == 1L &&
    !is.na(x) && (mode != "numeric" || is.finite(x))
}


.dpprior_v2_exit_zero <- function(x) {
  .dpprior_v2_plain_scalar(x, "numeric") && x == 0
}


.dpprior_v2_optimizer_evidence_allowed <- function(source) {
  .dpprior_v2_plain_scalar(source, "character") &&
    source %in% c("K_only_L-BFGS-B", "constrained_profile_optimizer")
}


.dpprior_v2_hard_intrinsic_usable <- function(status) {
  .dpprior_v2_plain_scalar(status, "character") &&
    status %in% c("converged", "boundary")
}


.dpprior_v2_hard_return_policy <- function(status, allow_approximate) {
  allow_approximate <- .dpprior_validate_control(
    allow_approximate, "allow_approximate", type = "logical"
  )
  usable <- .dpprior_v2_hard_intrinsic_usable(status)
  list(
    usable = usable,
    return_result = usable ||
      (identical(status, "approximate") && allow_approximate)
  )
}


.dpprior_v2_require_decision_ready_fit <- function(
    fit_info, mode = c("hard", "soft")) {
  mode <- match.arg(mode)
  if (!is.list(fit_info) || !isTRUE(fit_info$usable) ||
      !isTRUE(fit_info$verified)) {
    .dpprior_v2_abort_invalid(
      paste(
        "decision-ready dual calibration requires an input K-only fit with",
        "usable = TRUE and verified = TRUE; unverified A1/legacy proxies",
        "cannot be implicitly superseded"
      ),
      "fit",
      if (is.list(fit_info)) {
        list(status = fit_info$status, usable = fit_info$usable,
             verified = fit_info$verified, mode = mode)
      } else {
        fit_info
      },
      "a usable and independently verified K-only DPprior_fit",
      "input_fit_not_decision_ready",
      c("dpprior_dual_fit_error", "dpprior_status_contract_error")
    )
  }
  invisible(fit_info)
}


.dpprior_v2_validate_named_list <- function(x, name, allowed = NULL,
                                            required = character()) {
  if (!is.list(x) || is.object(x)) {
    .dpprior_v2_abort_invalid(
      sprintf("%s must be an ordinary named list", name),
      name, x, "ordinary named list", "type",
      c("dpprior_dual_type_error", "dpprior_type_error")
    )
  }
  nms <- names(x)
  if (length(x) && (is.null(nms) || anyNA(nms) || any(!nzchar(nms)))) {
    .dpprior_v2_abort_invalid(
      sprintf("every %s component must be named", name),
      name, x, "uniquely named components", "names",
      c("dpprior_dual_name_error", "dpprior_type_error")
    )
  }
  if (anyDuplicated(nms)) {
    .dpprior_v2_abort_invalid(
      sprintf("%s component names must be unique", name),
      name, nms, "unique names", "duplicate_names",
      "dpprior_dual_name_error"
    )
  }
  if (!is.null(allowed)) {
    unknown <- setdiff(nms, allowed)
    if (length(unknown)) {
      .dpprior_v2_abort_invalid(
        sprintf("unknown %s component(s): %s", name,
                paste(unknown, collapse = ", ")),
        name, unknown, paste(allowed, collapse = ", "), "unknown_component",
        c("dpprior_dual_control_error", "dpprior_unknown_control_error")
      )
    }
  }
  missing_required <- setdiff(required, nms)
  if (length(missing_required)) {
    .dpprior_v2_abort_invalid(
      sprintf("%s must contain: %s", name,
              paste(missing_required, collapse = ", ")),
      name, nms, paste(required, collapse = ", "), "missing_component",
      c("dpprior_dual_specification_error", "dpprior_missing_error")
    )
  }
  x
}


.dpprior_v2_exact_field <- function(x, name) {
  if (!is.list(x) || is.object(x)) return(NULL)
  nms <- names(x)
  if (is.null(nms) || anyDuplicated(nms) || !(name %in% nms)) return(NULL)
  x[[name, exact = TRUE]]
}


.dpprior_v2_exact_record <- function(x) {
  if (!is.list(x)) return(NULL)
  record <- if (is.object(x)) unclass(x) else x
  nms <- names(record)
  if (is.null(nms) || anyDuplicated(nms)) return(NULL)
  record
}


.dpprior_v2_request_moment_pair <- function(record) {
  record <- .dpprior_v2_exact_record(record)
  if (is.null(record)) return(NULL)
  mu <- .dpprior_v2_exact_field(record, "mu_K") %||%
    .dpprior_v2_exact_field(record, "mean")
  variance <- .dpprior_v2_exact_field(record, "var_K") %||%
    .dpprior_v2_exact_field(record, "variance")
  cv <- .dpprior_v2_exact_field(record, "cv_K") %||%
    .dpprior_v2_exact_field(record, "cv")
  if (is.null(variance) &&
      .dpprior_v2_plain_scalar(mu, "numeric") &&
      .dpprior_v2_plain_scalar(cv, "numeric")) {
    variance <- (mu * cv)^2
  }
  pair <- c(mean = mu, variance = variance)
  if (length(pair) != 2L || any(vapply(
    pair, is.null, logical(1)
  ))) return(NULL)
  .dpprior_v2_moment_pair(pair)
}


.dpprior_v2_scalar_close <- function(left, right, abs_tol = 1e-12,
                                     rel_tol = 1e-10) {
  .dpprior_v2_plain_scalar(left, "numeric") &&
    .dpprior_v2_plain_scalar(right, "numeric") &&
    abs(left - right) <=
      abs_tol + rel_tol * max(abs(left), abs(right), 1)
}


.dpprior_v2_nullable_scalar_identity <- function(left, right) {
  if (is.null(left) || is.null(right)) return(is.null(left) && is.null(right))
  .dpprior_v2_scalar_close(left, right)
}


.dpprior_v2_wrapper_route_nulls <- function(wrapper, fields) {
  all(vapply(fields, function(field) {
    is.null(.dpprior_v2_exact_field(wrapper, field))
  }, logical(1)))
}


.dpprior_v2_canonical_confidence <- function(raw, used) {
  choices <- c("medium", "low", "high")
  if (!.dpprior_v2_plain_scalar(used, "character") ||
      !(used %in% choices)) return(NULL)
  if (is.null(raw)) {
    return(if (identical(used, "medium")) used else NULL)
  }
  valid_raw <- is.character(raw) && is.null(dim(raw)) && !is.object(raw) &&
    (length(raw) == 1L || identical(raw, choices)) && !anyNA(raw)
  if (!valid_raw) return(NULL)
  normalized <- tryCatch(match.arg(raw, choices), error = function(error) NULL)
  if (identical(normalized, used)) normalized else NULL
}


.dpprior_v2_canonicalize_pmf_request <- function(pmf, J, rule) {
  if (!is.numeric(pmf) || is.null(pmf) || !is.null(dim(pmf)) ||
      is.object(pmf) || anyNA(pmf) || any(!is.finite(pmf))) return(NULL)
  pmf <- unname(as.numeric(pmf))
  valid_rule <- .dpprior_v2_plain_scalar(rule, "character") &&
    rule %in% c("validate_strict_pmf", "drop_structural_k0_zero")
  expected_length <- if (identical(rule, "drop_structural_k0_zero")) {
    J + 1L
  } else {
    J
  }
  if (!valid_rule ||
      !.dpprior_v2_plain_probability_vector(pmf, expected_length)) return(NULL)
  normalized <- if (identical(rule, "drop_structural_k0_zero")) {
    if (!identical(pmf[[1L]], 0)) return(NULL)
    pmf[-1L]
  } else {
    pmf
  }
  if (!.dpprior_v2_plain_probability_vector(normalized, J)) return(NULL)
  list(request = pmf, used = unname(normalized))
}


.dpprior_v2_canonicalize_wrapper_interval <- function(wrapper, J) {
  raw_interval <- .dpprior_v2_exact_field(wrapper, "K_interval")
  interval <- tryCatch(
    .dp_validate_K_interval(raw_interval, J),
    error = function(error) NULL
  )
  if (is.null(interval)) return(NULL)
  wrapper_mean <- .dpprior_v2_exact_field(wrapper, "mu_K")
  interval_mean <- .dpprior_v2_exact_field(interval, "mu_K")
  if (!is.null(wrapper_mean) &&
      (!.dpprior_v2_plain_scalar(wrapper_mean, "numeric") ||
       wrapper_mean < 1 || wrapper_mean > J)) return(NULL)
  if (!is.null(wrapper_mean) && !is.null(interval_mean) &&
      !.dpprior_v2_scalar_close(
        wrapper_mean, interval_mean, abs_tol = 0, rel_tol = 1e-9
      )) return(NULL)
  resolved_mean <- wrapper_mean %||% interval_mean
  if (is.null(resolved_mean) &&
      identical(interval$type, "hard_bounds") &&
      identical(interval$lower, interval$upper)) {
    resolved_mean <- as.numeric(interval$lower)
  }
  interval["mu_K"] <- list(resolved_mean)
  interval
}


.dpprior_v2_interval_identity <- function(left, right, J) {
  left <- .dpprior_v2_exact_record(left)
  right <- .dpprior_v2_exact_record(right)
  fields <- c(
    "lower", "upper", "type", "coverage", "family", "mu_K", "support",
    "endpoints"
  )
  if (is.null(left) || is.null(right) ||
      !identical(names(left), fields) || !identical(names(right), fields)) {
    return(FALSE)
  }
  scalar_fields <- c("lower", "upper", "coverage")
  scalar_identity <- all(vapply(scalar_fields, function(field) {
    .dpprior_v2_scalar_close(left[[field]], right[[field]])
  }, logical(1)))
  mean_identity <- if (is.null(left$mu_K) || is.null(right$mu_K)) {
    is.null(left$mu_K) && is.null(right$mu_K)
  } else {
    .dpprior_v2_scalar_close(left$mu_K, right$mu_K)
  }
  scalar_identity && mean_identity &&
    identical(left$type, right$type) &&
    identical(left$family, right$family) &&
    identical(left$support, c(lower = 1L, upper = as.integer(J))) &&
    identical(right$support, left$support) &&
    identical(left$endpoints, "inclusive") &&
    identical(right$endpoints, left$endpoints)
}


.dpprior_v2_validate_retained_target_identity <- function(
    raw_fit, J, target, method) {
  wrapper <- .dpprior_v2_exact_record(
    .dpprior_v2_exact_field(raw_fit, "wrapper_request")
  )
  elicitation_object <- .dpprior_v2_exact_field(
    raw_fit, "elicitation_target"
  )
  elicitation <- .dpprior_v2_exact_record(
    elicitation_object
  )
  if (is.null(wrapper) && is.null(elicitation)) {
    return(list(
      available = FALSE, passed = TRUE,
      reason = "direct_backend_has_no_wrapper_target_sources",
      authoritative_pmf = NULL
    ))
  }
  family <- if (method %in% c("A2-MN", "A2-MN+NM")) {
    "A2-MN"
  } else {
    method
  }
  request <- .dpprior_v2_exact_record(
    .dpprior_v2_exact_field(elicitation, "request")
  )
  normalized <- .dpprior_v2_exact_record(
    .dpprior_v2_exact_field(elicitation, "normalized")
  )
  used <- .dpprior_v2_exact_record(
    .dpprior_v2_exact_field(elicitation, "used")
  )
  implied <- .dpprior_v2_moment_pair(
    .dpprior_v2_exact_field(elicitation, "implied")
  )
  target_pmf <- .dpprior_v2_exact_field(elicitation, "pmf")
  target_interval <- .dpprior_v2_exact_field(elicitation, "interval")
  target_family <- .dpprior_v2_exact_field(elicitation, "family")
  derivation <- .dpprior_v2_exact_record(
    .dpprior_v2_exact_field(elicitation, "derivation")
  )
  first_derivation <- .dpprior_v2_exact_record(
    .dpprior_v2_exact_field(derivation, "request_to_normalized")
  )
  first_rule <- .dpprior_v2_exact_field(first_derivation, "rule")
  target_kind <- .dpprior_v2_exact_field(elicitation, "kind")
  route <- if (identical(first_rule, "canonicalize_direct_moments") &&
               identical(target_kind, "moments")) {
    "direct"
  } else if (identical(first_rule, "canonicalize_confidence_target") &&
             identical(target_kind, "moments")) {
    "confidence"
  } else if (identical(first_rule, "canonicalize_cv_target") &&
             identical(target_kind, "cv")) {
    "cv"
  } else if (first_rule %in%
             c("validate_strict_pmf", "drop_structural_k0_zero") &&
             identical(target_kind, "pmf")) {
    "pmf"
  } else if (identical(first_rule, "canonicalize_interval_request") &&
             identical(target_kind, "interval")) {
    "interval"
  } else {
    NULL
  }
  canonical_valid <- isTRUE(tryCatch({
    .dpprior_validate_target_v1(elicitation_object)
    TRUE
  }, error = function(error) FALSE))
  target_pair <- c(mean = target$mu_K, variance = target$var_K)
  wrapper_J <- .dpprior_v2_exact_field(wrapper, "J")
  elicitation_J <- .dpprior_v2_exact_field(elicitation, "J")
  selected_method <- .dpprior_v2_exact_field(wrapper, "selected_method")
  basic_identity <-
    !is.null(wrapper) && !is.null(elicitation) && !is.null(request) &&
    !is.null(normalized) && !is.null(used) && !is.null(route) &&
    canonical_valid &&
    .dpprior_v2_plain_positive_integer(wrapper_J) && wrapper_J == J &&
    .dpprior_v2_plain_positive_integer(elicitation_J) &&
    elicitation_J == J &&
    .dpprior_v2_plain_scalar(selected_method, "character") &&
    identical(selected_method, family) &&
    .dpprior_v2_pairs_close(
      implied, target_pair, abs_tol = 1e-12, rel_tol = 1e-10
    )

  request_identity <- FALSE
  authoritative_pmf <- NULL
  if (identical(route, "direct")) {
    wrapper_pair <- .dpprior_v2_request_moment_pair(wrapper)
    request_identity <-
      .dpprior_v2_wrapper_route_nulls(
        wrapper,
        c("cv_K", "K_interval", "confidence", "confidence_used", "target_pmf")
      ) &&
      .dpprior_v2_pairs_close(
        wrapper_pair, .dpprior_v2_request_moment_pair(request),
        abs_tol = 1e-12, rel_tol = 1e-10
      ) &&
      .dpprior_v2_pairs_close(
        wrapper_pair, .dpprior_v2_request_moment_pair(normalized),
        abs_tol = 1e-12, rel_tol = 1e-10
      ) &&
      .dpprior_v2_pairs_close(
        wrapper_pair, .dpprior_v2_request_moment_pair(used),
        abs_tol = 1e-12, rel_tol = 1e-10
      ) &&
      .dpprior_v2_pairs_close(
        wrapper_pair, target_pair, abs_tol = 1e-12, rel_tol = 1e-10
      ) && is.null(target_pmf)
  } else if (identical(route, "confidence")) {
    wrapper_mean <- .dpprior_v2_exact_field(wrapper, "mu_K")
    confidence <- .dpprior_v2_canonical_confidence(
      .dpprior_v2_exact_field(wrapper, "confidence"),
      .dpprior_v2_exact_field(wrapper, "confidence_used")
    )
    vif <- c(low = 5, medium = 2.5, high = 1.5)
    expected_pair <- if (!is.null(confidence) &&
        .dpprior_v2_plain_scalar(wrapper_mean, "numeric")) {
      c(
        mean = wrapper_mean,
        variance = unname(vif[[confidence]]) * (wrapper_mean - 1)
      )
    } else {
      NULL
    }
    request_identity <-
      .dpprior_v2_wrapper_route_nulls(
        wrapper, c("var_K", "cv_K", "K_interval", "target_pmf")
      ) && !is.null(confidence) &&
      identical(.dpprior_v2_exact_field(request, "confidence"), confidence) &&
      identical(.dpprior_v2_exact_field(normalized, "confidence"), confidence) &&
      .dpprior_v2_scalar_close(
        wrapper_mean, .dpprior_v2_exact_field(request, "mean")
      ) &&
      .dpprior_v2_scalar_close(
        wrapper_mean, .dpprior_v2_exact_field(normalized, "mean")
      ) &&
      .dpprior_v2_pairs_close(
        expected_pair, .dpprior_v2_request_moment_pair(used),
        abs_tol = 1e-12, rel_tol = 1e-10
      ) &&
      .dpprior_v2_pairs_close(
        expected_pair, target_pair, abs_tol = 1e-12, rel_tol = 1e-10
      ) && is.null(target_pmf)
  } else if (identical(route, "cv")) {
    wrapper_mean <- .dpprior_v2_exact_field(wrapper, "mu_K")
    wrapper_cv <- .dpprior_v2_exact_field(wrapper, "cv_K")
    expected_pair <- if (
      .dpprior_v2_plain_scalar(wrapper_mean, "numeric") &&
        .dpprior_v2_plain_scalar(wrapper_cv, "numeric") && wrapper_cv > 0
    ) {
      c(mean = wrapper_mean, variance = (wrapper_mean * wrapper_cv)^2)
    } else {
      NULL
    }
    request_identity <-
      .dpprior_v2_wrapper_route_nulls(
        wrapper,
        c("var_K", "K_interval", "confidence", "confidence_used", "target_pmf")
      ) &&
      .dpprior_v2_scalar_close(
        wrapper_mean, .dpprior_v2_exact_field(request, "mean")
      ) &&
      .dpprior_v2_scalar_close(
        wrapper_mean, .dpprior_v2_exact_field(normalized, "mean")
      ) &&
      .dpprior_v2_scalar_close(
        wrapper_cv, .dpprior_v2_exact_field(request, "cv")
      ) &&
      .dpprior_v2_scalar_close(
        wrapper_cv, .dpprior_v2_exact_field(normalized, "cv")
      ) &&
      .dpprior_v2_pairs_close(
        expected_pair, .dpprior_v2_request_moment_pair(used),
        abs_tol = 1e-12, rel_tol = 1e-10
      ) &&
      .dpprior_v2_pairs_close(
        expected_pair, target_pair, abs_tol = 1e-12, rel_tol = 1e-10
      ) && is.null(target_pmf)
  } else if (identical(route, "pmf")) {
    wrapper_request_pmf <- .dpprior_v2_canonicalize_pmf_request(
      .dpprior_v2_exact_field(wrapper, "target_pmf"), J, first_rule
    )
    canonical_request_pmf <- .dpprior_v2_canonicalize_pmf_request(
      .dpprior_v2_exact_field(request, "pmf"), J, first_rule
    )
    normalized_pmf <- .dpprior_v2_exact_field(normalized, "pmf")
    used_pmf <- .dpprior_v2_exact_field(used, "pmf")
    wrapper_mean <- .dpprior_v2_exact_field(wrapper, "mu_K")
    wrapper_variance <- .dpprior_v2_exact_field(wrapper, "var_K")
    assertions_pass <-
      (is.null(wrapper_mean) ||
         .dpprior_v2_scalar_close(
           wrapper_mean, target_pair[["mean"]], abs_tol = 0, rel_tol = 1e-9
         )) &&
      (is.null(wrapper_variance) ||
         .dpprior_v2_scalar_close(
           wrapper_variance, target_pair[["variance"]],
           abs_tol = 0, rel_tol = 1e-9
         ))
    request_identity <-
      .dpprior_v2_wrapper_route_nulls(
        wrapper, c("cv_K", "K_interval", "confidence", "confidence_used")
      ) && assertions_pass && !is.null(wrapper_request_pmf) &&
      !is.null(canonical_request_pmf) &&
      .dpprior_v2_vectors_close(
        wrapper_request_pmf$request, canonical_request_pmf$request,
        abs_tol = 1e-12, rel_tol = 1e-12
      ) &&
      .dpprior_v2_vectors_close(
        wrapper_request_pmf$used, canonical_request_pmf$used,
        abs_tol = 1e-12, rel_tol = 1e-12
      ) &&
      .dpprior_v2_plain_probability_vector(normalized_pmf, J) &&
      .dpprior_v2_plain_probability_vector(used_pmf, J) &&
      .dpprior_v2_plain_probability_vector(target_pmf, J) &&
      .dpprior_v2_vectors_close(
        wrapper_request_pmf$used, normalized_pmf,
        abs_tol = 1e-12, rel_tol = 1e-12
      ) &&
      .dpprior_v2_vectors_close(
        normalized_pmf, used_pmf, abs_tol = 1e-12, rel_tol = 1e-12
      ) &&
      .dpprior_v2_vectors_close(
        used_pmf, target_pmf,
        abs_tol = 1e-12, rel_tol = 1e-12
      )
    if (request_identity) authoritative_pmf <- unname(target_pmf)
  } else if (identical(route, "interval")) {
    wrapper_interval <- .dpprior_v2_canonicalize_wrapper_interval(wrapper, J)
    request_interval <- .dpprior_v2_exact_field(request, "K_interval")
    normalized_interval <- .dpprior_v2_exact_field(normalized, "interval")
    used_interval <- .dpprior_v2_exact_field(used, "interval")
    normalized_family <- .dpprior_v2_exact_field(normalized, "family")
    used_family <- .dpprior_v2_exact_field(used, "family")
    normalized_pmf <- .dpprior_v2_exact_field(normalized, "pmf")
    used_pmf <- .dpprior_v2_exact_field(used, "pmf")
    request_identity <-
      .dpprior_v2_wrapper_route_nulls(
        wrapper,
        c("var_K", "cv_K", "confidence", "confidence_used", "target_pmf")
      ) &&
      .dpprior_v2_interval_identity(wrapper_interval, request_interval, J) &&
      .dpprior_v2_interval_identity(wrapper_interval, normalized_interval, J) &&
      .dpprior_v2_interval_identity(wrapper_interval, used_interval, J) &&
      .dpprior_v2_interval_identity(wrapper_interval, target_interval, J) &&
      .dpprior_v2_nullable_scalar_identity(
        .dpprior_v2_exact_field(wrapper_interval, "mu_K"),
        .dpprior_v2_exact_field(request, "mu_K")
      ) && identical(normalized_family, target_family) &&
      identical(used_family, target_family) &&
      identical(
        .dpprior_v2_exact_field(wrapper_interval, "family"),
        .dpprior_v2_exact_field(target_family, "name")
      ) && is.null(normalized_pmf) &&
      .dpprior_v2_plain_probability_vector(used_pmf, J) &&
      .dpprior_v2_plain_probability_vector(target_pmf, J) &&
      .dpprior_v2_vectors_close(
        used_pmf, target_pmf, abs_tol = 1e-12, rel_tol = 1e-12
      )
    if (request_identity) authoritative_pmf <- unname(target_pmf)
  }
  if (!is.null(authoritative_pmf)) {
    k <- seq_len(J)
    pmf_mean <- sum(k * authoritative_pmf)
    pmf_variance <- sum((k - pmf_mean)^2 * authoritative_pmf)
    request_identity <- request_identity && .dpprior_v2_pairs_close(
      c(mean = pmf_mean, variance = pmf_variance), target_pair,
      abs_tol = 1e-12, rel_tol = 1e-10
    )
  }
  passed <- basic_identity && request_identity
  if (!passed) {
    .dpprior_v2_abort_invalid(
      paste(
        "retained wrapper request, elicitation target, and K target",
        "do not describe one immutable target"
      ),
      "fit retained target sources",
      list(
        method = method, target = target_pair,
        wrapper_request = wrapper, elicitation_target = elicitation,
        canonical_target_valid = canonical_valid, route = route,
        basic_identity = basic_identity,
        request_identity = request_identity
      ),
      "exactly reconciled wrapper, elicitation, and backend target sources",
      "conflicting_retained_target_identity",
      c("dpprior_dual_fit_error", "dpprior_conflicting_input")
    )
  }
  list(
    available = TRUE, passed = TRUE,
    reason = "wrapper_elicitation_and_backend_target_reconciled",
    target_kind = if (is.null(authoritative_pmf)) "moments" else "pmf",
    authoritative_pmf = authoritative_pmf
  )
}


.dpprior_v2_plain_positive_integer <- function(x, maximum = Inf) {
  .dpprior_v2_plain_scalar(x, "numeric") && x == floor(x) && x >= 1 &&
    x <= maximum
}


.dpprior_v2_moment_pair <- function(x) {
  if (is.numeric(x) && is.null(dim(x)) && !is.object(x) &&
      length(x) >= 2L && !is.null(names(x)) &&
      !anyDuplicated(names(x)) && all(c("mean", "variance") %in% names(x))) {
    pair <- x[c("mean", "variance")]
  } else if (is.list(x) && !is.object(x) && !anyDuplicated(names(x))) {
    mean_value <- .dpprior_v2_exact_field(x, "mean") %||%
      .dpprior_v2_exact_field(x, "mu_K")
    variance_value <- .dpprior_v2_exact_field(x, "variance") %||%
      .dpprior_v2_exact_field(x, "var_K")
    if (is.null(mean_value) || is.null(variance_value)) return(NULL)
    pair <- c(mean = mean_value, variance = variance_value)
  } else {
    return(NULL)
  }
  if (!is.numeric(pair) || length(pair) != 2L || is.object(pair) ||
      anyNA(pair) || any(!is.finite(pair))) return(NULL)
  stats::setNames(as.numeric(pair), c("mean", "variance"))
}


.dpprior_v2_pairs_close <- function(left, right, abs_tol = 1e-10,
                                    rel_tol = 1e-8) {
  left <- .dpprior_v2_moment_pair(left)
  right <- .dpprior_v2_moment_pair(right)
  !is.null(left) && !is.null(right) && all(
    abs(left - right) <=
      abs_tol + rel_tol * pmax(abs(left), abs(right), 1)
  )
}


.dpprior_v2_plain_probability_vector <- function(x, length_required) {
  is.numeric(x) && is.null(dim(x)) && !is.object(x) &&
    length(x) == length_required && !anyNA(x) && all(is.finite(x)) &&
    all(x >= 0) && all(x <= 1) &&
    abs(sum(x) - 1) <= 1e-10
}


.dpprior_v2_vectors_close <- function(left, right, abs_tol = 1e-10,
                                      rel_tol = 1e-8) {
  is.numeric(left) && is.null(dim(left)) && !is.object(left) &&
    is.numeric(right) && is.null(dim(right)) && !is.object(right) &&
    length(left) == length(right) && length(left) > 0L &&
    !anyNA(left) && !anyNA(right) &&
    all(is.finite(left)) && all(is.finite(right)) &&
    all(abs(left - right) <=
      abs_tol + rel_tol * pmax(abs(left), abs(right), 1e-8))
}


.dpprior_v2_KL_metrics_close <- function(recorded, fresh) {
  fields <- c(
    "kl", "l1", "mean", "variance", "mean_residual",
    "variance_residual", "target_mean", "target_variance"
  )
  if (!is.list(recorded) || is.object(recorded) ||
      anyDuplicated(names(recorded))) return(FALSE)
  recorded_values <- vapply(fields, function(name) {
    value <- .dpprior_v2_exact_field(recorded, name)
    if (.dpprior_v2_plain_scalar(value, "numeric")) value else NA_real_
  }, numeric(1))
  fresh_values <- vapply(fields, function(name) fresh[[name]], numeric(1))
  all(is.finite(recorded_values)) &&
    .dpprior_v2_vectors_close(
      unname(recorded_values), unname(fresh_values),
      abs_tol = 1e-12, rel_tol = 1e-10
    ) &&
    .dpprior_v2_vectors_close(
      .dpprior_v2_exact_field(recorded, "pmf"), fresh$pmf,
      abs_tol = 1e-12, rel_tol = 1e-10
    )
}


.dpprior_v2_KL_adequacy_close <- function(recorded, fresh) {
  if (!is.list(recorded) || is.object(recorded) ||
      anyDuplicated(names(recorded)) ||
      !isTRUE(.dpprior_v2_exact_field(recorded, "passed"))) return(FALSE)
  numeric_fields <- c(
    "raw", "scales", "scaled", "scaled_tolerances",
    "raw_tolerances", "ratios"
  )
  all(vapply(numeric_fields, function(name) {
    .dpprior_v2_vectors_close(
      .dpprior_v2_exact_field(recorded, name), fresh[[name]],
      abs_tol = 1e-12, rel_tol = 1e-10
    )
  }, logical(1))) && identical(
    .dpprior_v2_exact_field(recorded, "component_passed"),
    fresh$component_passed
  )
}


.dpprior_v2_validate_A2_KL_input_evidence <- function(
    raw_fit, verification, wrapper, a, b, J, target,
    retained_target_identity, M_selected, M_refined) {
  target_record <- .dpprior_v2_exact_field(raw_fit, "target")
  target_type <- .dpprior_v2_exact_field(target_record, "type")
  target_pmf <- .dpprior_v2_exact_field(target_record, "pmf")
  wrapper_target <- .dpprior_v2_exact_field(wrapper, "target_pmf")
  authoritative_pmf <- retained_target_identity$authoritative_pmf
  wrapper_flags <- c(
    numerical_verification_passed = isTRUE(.dpprior_v2_exact_field(
      wrapper, "numerical_verification_passed"
    )),
    tolerance_policy_passed = isTRUE(.dpprior_v2_exact_field(
      wrapper, "tolerance_policy_passed"
    )),
    public_metadata_passed = isTRUE(.dpprior_v2_exact_field(
      wrapper, "public_metadata_passed"
    ))
  )
  target_container_valid <-
    .dpprior_v2_plain_scalar(target_type, "character") &&
    target_type %in% c("pmf", "chisq") &&
    .dpprior_v2_plain_probability_vector(target_pmf, J) &&
    .dpprior_v2_plain_probability_vector(wrapper_target, J) &&
    .dpprior_v2_vectors_close(
      target_pmf, wrapper_target, abs_tol = 1e-12, rel_tol = 1e-12
    )
  if (!is.null(authoritative_pmf)) {
    target_container_valid <- target_container_valid &&
      .dpprior_v2_plain_probability_vector(authoritative_pmf, J) &&
      .dpprior_v2_vectors_close(
        target_pmf, authoritative_pmf,
        abs_tol = 1e-12, rel_tol = 1e-12
      )
  }
  target_identity_pass <- FALSE
  if (target_container_valid && identical(target_type, "pmf")) {
    k <- seq_len(J)
    pmf_mean <- sum(k * target_pmf)
    pmf_variance <- sum((k - pmf_mean)^2 * target_pmf)
    target_identity_pass <- .dpprior_v2_pairs_close(
      c(mean = pmf_mean, variance = pmf_variance),
      c(mean = target$mu_K, variance = target$var_K),
      abs_tol = 1e-12, rel_tol = 1e-10
    )
  } else if (target_container_valid && identical(target_type, "chisq") &&
             target$var_K > 0) {
    expected_df <- 2 * target$mu_K^2 / target$var_K
    expected_scale <- target$var_K / (2 * target$mu_K)
    recorded_df <- .dpprior_v2_exact_field(target_record, "df")
    recorded_scale <- .dpprior_v2_exact_field(target_record, "scale")
    rebuilt_target <- tryCatch(
      discretize_chisq(J, df = expected_df, scale = expected_scale),
      error = identity
    )
    target_identity_pass <-
      .dpprior_v2_plain_scalar(recorded_df, "numeric") &&
      .dpprior_v2_plain_scalar(recorded_scale, "numeric") &&
      abs(recorded_df - expected_df) <=
        1e-12 + 1e-10 * max(abs(recorded_df), abs(expected_df), 1) &&
      abs(recorded_scale - expected_scale) <=
        1e-12 + 1e-10 * max(abs(recorded_scale), abs(expected_scale), 1) &&
      !inherits(rebuilt_target, "condition") &&
      .dpprior_v2_vectors_close(
        target_pmf, rebuilt_target, abs_tol = 1e-12, rel_tol = 1e-10
      )
  }
  target_valid <- target_container_valid && target_identity_pass
  if (is.null(wrapper) || !all(wrapper_flags) || !target_valid) {
    return(list(
      passed = FALSE, reason = "missing_canonical_A2_KL_wrapper_or_target_PMF",
      target_valid = target_valid,
      target_container_valid = target_container_valid,
      target_identity_pass = target_identity_pass,
      wrapper_flags = wrapper_flags
    ))
  }

  bundle <- tryCatch(
    pmf_K_marginal(
      J = J, a = a, b = b, logS = compute_log_stirling(J),
      M = M_selected, M_verify = M_refined, strict = FALSE
    ),
    error = identity
  )
  if (inherits(bundle, "condition")) {
    return(list(
      passed = FALSE, reason = "fresh_A2_KL_PMF_recomputation_failed",
      error = conditionMessage(bundle)
    ))
  }
  refined_bundle <- attr(
    bundle, ".marginal_verification_pmf", exact = TRUE
  )
  metadata <- attr(bundle, "marginal_metadata", exact = TRUE)
  metadata_verification <- .dpprior_v2_exact_field(metadata, "verification")
  selected_pmf <- if (is.numeric(bundle) && length(bundle) == J + 1L) {
    unname(bundle[-1L])
  } else {
    NULL
  }
  refined_pmf <- if (is.numeric(refined_bundle) &&
                     length(refined_bundle) == J + 1L) {
    unname(refined_bundle[-1L])
  } else {
    NULL
  }
  numerical_PMF_pass <-
    .dpprior_v2_plain_probability_vector(selected_pmf, J) &&
    .dpprior_v2_plain_probability_vector(refined_pmf, J) &&
    isTRUE(.dpprior_v2_exact_field(metadata_verification, "performed")) &&
    isTRUE(.dpprior_v2_exact_field(metadata_verification, "passed"))
  if (!numerical_PMF_pass) {
    return(list(
      passed = FALSE, reason = "fresh_A2_KL_PMF_error_control_failed"
    ))
  }

  selected_metrics <- tryCatch(
    .a2_kl_pmf_metrics(
      target_pmf, ifelse(selected_pmf > 0, log(selected_pmf), -Inf)
    ),
    error = identity
  )
  refined_metrics <- tryCatch(
    .a2_kl_pmf_metrics(
      target_pmf, ifelse(refined_pmf > 0, log(refined_pmf), -Inf)
    ),
    error = identity
  )
  if (inherits(selected_metrics, "condition") ||
      inherits(refined_metrics, "condition")) {
    return(list(
      passed = FALSE, reason = "fresh_A2_KL_metric_recomputation_failed"
    ))
  }
  fixed_tolerances <- list(
    kl = 0.015, l1 = 0.11, mean_scaled = 0.01,
    variance_scaled = 0.065
  )
  selected_adequacy <- .a2_kl_assess_adequacy(
    selected_metrics, fixed_tolerances
  )
  refined_adequacy <- .a2_kl_assess_adequacy(
    refined_metrics, fixed_tolerances
  )

  attempts <- .dpprior_v2_exact_field(raw_fit, "attempts")
  selected_optimizer_attempt <- NULL
  if (is.list(attempts) && !is.object(attempts)) {
    selected_indices <- which(vapply(attempts, function(attempt) {
      attempt <- .dpprior_v2_exact_record(attempt)
      !is.null(attempt) &&
        identical(.dpprior_v2_exact_field(attempt, "stage"), "optimizer") &&
        isTRUE(.dpprior_v2_exact_field(attempt, "selected"))
    }, logical(1)))
    if (length(selected_indices) == 1L) {
      selected_optimizer_attempt <- .dpprior_v2_exact_record(
        attempts[[selected_indices]]
      )
    }
  }
  attempt_candidate <- .dpprior_v2_exact_field(
    selected_optimizer_attempt, "candidate"
  )
  attempt_candidate_pass <-
    is.numeric(attempt_candidate) && is.null(dim(attempt_candidate)) &&
    !is.object(attempt_candidate) && length(attempt_candidate) == 2L &&
    !is.null(names(attempt_candidate)) &&
    all(c("a", "b") %in% names(attempt_candidate)) &&
    .dpprior_v2_vectors_close(
      unname(attempt_candidate[c("a", "b")]), c(a, b),
      abs_tol = 1e-12, rel_tol = 1e-10
    )
  attempt_objective <- .dpprior_v2_exact_field(
    selected_optimizer_attempt, "candidate_objective"
  )
  attempt_evidence_pass <-
    .dpprior_v2_exit_zero(.dpprior_v2_exact_field(
      selected_optimizer_attempt, "exit_code"
    )) && attempt_candidate_pass &&
    .dpprior_v2_plain_scalar(attempt_objective, "numeric") &&
    abs(attempt_objective - selected_metrics$kl) <=
      1e-12 + 1e-10 * max(abs(attempt_objective),
                          abs(selected_metrics$kl), 1)

  objective_at <- function(eta) {
    value <- tryCatch({
      local_bundle <- pmf_K_marginal(
        J = J, a = exp(eta[[1L]]), b = exp(eta[[2L]]),
        logS = compute_log_stirling(J), M = M_selected, strict = FALSE
      )
      local_pmf <- unname(local_bundle[-1L])
      .a2_kl_pmf_metrics(
        target_pmf, ifelse(local_pmf > 0, log(local_pmf), -Inf)
      )$kl
    }, error = function(condition) NA_real_)
    if (.dpprior_v2_plain_scalar(value, "numeric")) value else NA_real_
  }
  eta <- log(c(a, b))
  local_step <- 1e-5
  local_gradient <- rep(NA_real_, 2L)
  local_values <- vector("list", 2L)
  for (dimension in seq_len(2L)) {
    lower_eta <- upper_eta <- eta
    lower_eta[[dimension]] <- lower_eta[[dimension]] - local_step
    upper_eta[[dimension]] <- upper_eta[[dimension]] + local_step
    lower_value <- if (lower_eta[[dimension]] >= -15) {
      objective_at(lower_eta)
    } else {
      NA_real_
    }
    upper_value <- if (upper_eta[[dimension]] <= 15) {
      objective_at(upper_eta)
    } else {
      NA_real_
    }
    local_values[[dimension]] <- c(lower = lower_value, upper = upper_value)
    if (is.finite(lower_value) && is.finite(upper_value)) {
      local_gradient[[dimension]] <-
        (upper_value - lower_value) / (2 * local_step)
    } else if (is.finite(upper_value)) {
      local_gradient[[dimension]] <-
        (upper_value - selected_metrics$kl) / local_step
    } else if (is.finite(lower_value)) {
      local_gradient[[dimension]] <-
        (selected_metrics$kl - lower_value) / local_step
    }
  }
  projected_gradient <- local_gradient
  at_lower <- eta <= -15 + local_step
  at_upper <- eta >= 15 - local_step
  projected_gradient[at_lower] <- pmin(projected_gradient[at_lower], 0)
  projected_gradient[at_upper] <- pmax(projected_gradient[at_upper], 0)
  input_status <- .dpprior_v2_exact_field(raw_fit, "status")
  boundary_distance <- min(eta + 15, 15 - eta)
  verified_boundary_claim <- identical(input_status, "boundary") &&
    is.finite(boundary_distance) && boundary_distance <= 1e-6
  bound_coordinate <- verified_boundary_claim & (at_lower | at_upper)
  local_optimality_tolerance <- ifelse(bound_coordinate, 1e-3, 1e-4)
  local_optimality_pass <- all(is.finite(projected_gradient)) &&
    (!identical(input_status, "boundary") || verified_boundary_claim) &&
    all(abs(projected_gradient) <= local_optimality_tolerance)

  recorded_selected_metrics <- .dpprior_v2_exact_field(
    verification, "selected_metrics"
  )
  recorded_refined_metrics <- .dpprior_v2_exact_field(
    verification, "verification_metrics"
  )
  achieved <- .dpprior_v2_exact_field(raw_fit, "achieved")
  fit_summary <- .dpprior_v2_exact_field(raw_fit, "fit")
  achieved_scalar_pass <- all(vapply(c("kl", "l1"), function(name) {
    value <- .dpprior_v2_exact_field(achieved, name)
    .dpprior_v2_plain_scalar(value, "numeric") &&
      abs(value - selected_metrics[[name]]) <=
        1e-12 + 1e-10 * max(abs(value), abs(selected_metrics[[name]]), 1)
  }, logical(1)))
  fit_summary_pass <- all(vapply(
    c("mu_K", "var_K", "kl", "l1"), function(name) {
      value <- .dpprior_v2_exact_field(fit_summary, name)
      fresh_value <- switch(
        name, mu_K = selected_metrics$mean,
        var_K = selected_metrics$variance,
        kl = selected_metrics$kl, l1 = selected_metrics$l1
      )
      .dpprior_v2_plain_scalar(value, "numeric") &&
        abs(value - fresh_value) <=
          1e-12 + 1e-10 * max(abs(value), abs(fresh_value), 1)
    }, logical(1)
  ))
  adequacy_record_pass <-
    .dpprior_v2_KL_adequacy_close(
      .dpprior_v2_exact_field(wrapper, "selected_adequacy"),
      selected_adequacy
    ) &&
    .dpprior_v2_KL_adequacy_close(
      .dpprior_v2_exact_field(wrapper, "verification_adequacy"),
      refined_adequacy
    ) &&
    .dpprior_v2_KL_adequacy_close(
      .dpprior_v2_exact_field(verification, "adequacy"),
      refined_adequacy
    )
  recorded_pass <-
    .dpprior_v2_KL_metrics_close(
      recorded_selected_metrics, selected_metrics
    ) &&
    .dpprior_v2_KL_metrics_close(
      recorded_refined_metrics, refined_metrics
    ) &&
    .dpprior_v2_vectors_close(
      .dpprior_v2_exact_field(wrapper, "selected_pmf"), selected_pmf,
      abs_tol = 1e-12, rel_tol = 1e-10
    ) &&
    .dpprior_v2_vectors_close(
      .dpprior_v2_exact_field(wrapper, "verification_pmf"), refined_pmf,
      abs_tol = 1e-12, rel_tol = 1e-10
    ) &&
    .dpprior_v2_vectors_close(
      .dpprior_v2_exact_field(achieved, "pmf"), selected_pmf,
      abs_tol = 1e-12, rel_tol = 1e-10
    ) && achieved_scalar_pass && fit_summary_pass && adequacy_record_pass
  fixed_adequacy_pass <- isTRUE(selected_adequacy$passed) &&
    isTRUE(refined_adequacy$passed)
  passed <- numerical_PMF_pass && fixed_adequacy_pass && recorded_pass &&
    attempt_evidence_pass && local_optimality_pass
  list(
    passed = passed,
    reason = if (passed) {
      "fresh_selected_refined_PMF_and_fixed_adequacy_passed"
    } else if (!fixed_adequacy_pass) {
      "fixed_A2_KL_adequacy_failed"
    } else if (!attempt_evidence_pass) {
      "A2_KL_selected_attempt_or_objective_identity_failed"
    } else if (!local_optimality_pass) {
      "A2_KL_fresh_local_optimality_failed"
    } else {
      "cached_A2_KL_evidence_disagrees_with_fresh_recomputation"
    },
    target_type = target_type,
    target_pmf = target_pmf,
    selected_pmf = selected_pmf,
    refined_pmf = refined_pmf,
    selected_metrics = selected_metrics,
    refined_metrics = refined_metrics,
    selected_adequacy = selected_adequacy,
    refined_adequacy = refined_adequacy,
    fixed_tolerances = fixed_tolerances,
    numerical_PMF_pass = numerical_PMF_pass,
    fixed_adequacy_pass = fixed_adequacy_pass,
    recorded_pass = recorded_pass,
    attempt_evidence_pass = attempt_evidence_pass,
    local_step = local_step,
    local_values = local_values,
    local_gradient = stats::setNames(local_gradient, c("log_a", "log_b")),
    projected_gradient = stats::setNames(
      projected_gradient, c("log_a", "log_b")
    ),
    boundary_distance = boundary_distance,
    verified_boundary_claim = verified_boundary_claim,
    local_optimality_tolerance = stats::setNames(
      local_optimality_tolerance, c("log_a", "log_b")
    ),
    local_optimality_pass = local_optimality_pass
  )
}


.dpprior_v2_reconcile_order <- function(values, label) {
  values <- Filter(Negate(is.null), values)
  valid <- length(values) > 0L && all(vapply(
    values, .dpprior_v2_plain_positive_integer, logical(1), maximum = 512L
  ))
  if (!valid || length(unique(as.integer(unlist(values)))) != 1L) {
    .dpprior_v2_abort_invalid(
      sprintf("fit %s order metadata is missing, invalid, or contradictory",
              label),
      sprintf("fit verification %s order", label), values,
      "one consistent integer quadrature order in [1, 512]",
      "invalid_input_fit_verification_order",
      c("dpprior_dual_fit_error", "dpprior_status_contract_error")
    )
  }
  as.integer(values[[1L]])
}


# Cached status flags on an S3 object are not verification evidence.  Before a
# usable, verified K fit can define either dual endpoint, recompute its returned
# Gamma parameters at both the selected and refined quadrature orders and bind
# those values to the method-specific evidence retained by the K-only fit.
.dpprior_v2_validate_input_fit_evidence <- function(
    raw_fit, a, b, J, target, retained_target_identity,
    status, usable, verified, method) {
  if (!isTRUE(verified)) {
    return(list(
      required = FALSE, performed = FALSE, passed = FALSE,
      reason = "input_fit_does_not_make_a_verified_claim"
    ))
  }

  family <- if (method %in% c("A2-MN", "A2-MN+NM")) {
    "A2-MN"
  } else if (identical(method, "A2-KL")) {
    "A2-KL"
  } else {
    NA_character_
  }
  verification <- .dpprior_v2_exact_field(raw_fit, "verification")
  if (is.na(family) || !is.list(verification) || is.object(verification) ||
      anyDuplicated(names(verification)) ||
      !isTRUE(.dpprior_v2_exact_field(verification, "performed")) ||
      !isTRUE(.dpprior_v2_exact_field(verification, "passed"))) {
    .dpprior_v2_abort_invalid(
      paste(
        "a verified input fit requires exact method-appropriate",
        "performed-and-passed verification evidence"
      ),
      "fit$verification", verification,
      "substantive A2-MN or A2-KL verification evidence",
      "missing_input_fit_verification",
      c("dpprior_dual_fit_error", "dpprior_status_contract_error")
    )
  }

  wrapper <- .dpprior_v2_exact_field(
    verification, "wrapper_parameter_reverification"
  )
  if (!is.null(wrapper) &&
      (!is.list(wrapper) || is.object(wrapper) ||
       anyDuplicated(names(wrapper)) ||
       !isTRUE(.dpprior_v2_exact_field(wrapper, "performed")) ||
       !isTRUE(.dpprior_v2_exact_field(wrapper, "passed")))) {
    .dpprior_v2_abort_invalid(
      "fit wrapper parameter reverification is malformed or did not pass",
      "fit$verification$wrapper_parameter_reverification", wrapper,
      "performed-and-passed wrapper reverification",
      "invalid_input_fit_wrapper_verification",
      c("dpprior_dual_fit_error", "dpprior_status_contract_error")
    )
  }

  fit_M <- .dpprior_v2_exact_field(raw_fit, "M")
  if (identical(family, "A2-MN")) {
    selected_record <- .dpprior_v2_exact_field(verification, "selected")
    refined_record <- .dpprior_v2_exact_field(verification, "recomputed")
    verification_selected_M <- .dpprior_v2_exact_field(
      verification, "M_selected"
    )
    verification_refined_M <- .dpprior_v2_exact_field(
      verification, "M_verification"
    )
  } else {
    settings <- .dpprior_v2_exact_field(verification, "settings")
    selected_metrics <- .dpprior_v2_exact_field(
      verification, "selected_metrics"
    )
    refined_metrics <- .dpprior_v2_exact_field(
      verification, "verification_metrics"
    )
    selected_record <- .dpprior_v2_moment_pair(selected_metrics)
    refined_record <- .dpprior_v2_moment_pair(refined_metrics)
    verification_selected_M <- .dpprior_v2_exact_field(
      settings, "M_selected"
    )
    verification_refined_M <- .dpprior_v2_exact_field(
      settings, "M_verification"
    )
  }

  wrapper_selected_M <- .dpprior_v2_exact_field(wrapper, "selected_order")
  wrapper_refined_M <- .dpprior_v2_exact_field(wrapper, "verification_order")
  M_selected <- .dpprior_v2_reconcile_order(
    list(fit_M, verification_selected_M, wrapper_selected_M), "selected"
  )
  M_refined <- .dpprior_v2_reconcile_order(
    list(verification_refined_M, wrapper_refined_M), "refined"
  )
  required_refined <- as.integer(max(2L * M_selected, M_selected + 40L))
  if (M_refined < required_refined || M_refined <= M_selected) {
    .dpprior_v2_abort_invalid(
      "fit verification order is not independently refined",
      "fit verification order",
      list(selected = M_selected, refined = M_refined,
           required = required_refined),
      sprintf("refined order at least %d and greater than selected order",
              required_refined),
      "insufficient_input_fit_verification_order",
      c("dpprior_dual_fit_error", "dpprior_status_contract_error")
    )
  }

  fresh <- tryCatch(
    list(
      selected = exact_K_moments(J, a, b, M = M_selected),
      refined = exact_K_moments(J, a, b, M = M_refined)
    ),
    error = function(condition) condition
  )
  if (inherits(fresh, "condition")) {
    .dpprior_v2_abort_invalid(
      "independent input-fit parameter recomputation failed",
      "fit parameters", list(a = a, b = b, J = J),
      "parameters supporting finite selected/refined K moments",
      "input_fit_recomputation_failed",
      c("dpprior_dual_fit_error", "dpprior_status_contract_error")
    )
  }
  fresh_selected <- c(mean = fresh$selected$mean,
                      variance = fresh$selected$var)
  fresh_refined <- c(mean = fresh$refined$mean,
                     variance = fresh$refined$var)
  target_pair <- c(mean = target$mu_K, variance = target$var_K)
  achieved_record <- .dpprior_v2_moment_pair(
    .dpprior_v2_exact_field(raw_fit, "achieved")
  )

  recorded_pass <-
    .dpprior_v2_pairs_close(fresh_selected, achieved_record) &&
    .dpprior_v2_pairs_close(fresh_selected, selected_record) &&
    .dpprior_v2_pairs_close(fresh_refined, refined_record)
  wrapper_pass <- TRUE
  if (!is.null(wrapper) && identical(family, "A2-MN")) {
    wrapper_pass <-
      .dpprior_v2_pairs_close(
        fresh_selected, .dpprior_v2_exact_field(wrapper, "selected")
      ) &&
      .dpprior_v2_pairs_close(
        fresh_refined, .dpprior_v2_exact_field(wrapper, "verification")
      )
  }

  stability_tolerance <- 1e-10 + 1e-8 * pmax(
    abs(fresh_selected), abs(fresh_refined), 1
  )
  stability_pass <- all(
    abs(fresh_selected - fresh_refined) <= stability_tolerance
  )
  target_tolerance <- 1e-8 + 1e-8 * pmax(abs(target_pair), 1)
  target_pass <- identical(family, "A2-KL") || (
    all(abs(fresh_selected - target_pair) <= target_tolerance) &&
      all(abs(fresh_refined - target_pair) <= target_tolerance)
  )

  KL_evidence <- NULL
  method_evidence_pass <- if (identical(family, "A2-MN")) {
    residuals <- .dpprior_v2_exact_field(verification, "residuals")
    selected_residual <- .dpprior_v2_exact_field(residuals, "selected")
    refined_residual <- .dpprior_v2_exact_field(residuals, "verification")
    residual_record_pass <- function(record, fresh_pair) {
      raw <- .dpprior_v2_moment_pair(
        .dpprior_v2_exact_field(record, "raw")
      )
      tolerance <- .dpprior_v2_moment_pair(
        .dpprior_v2_exact_field(record, "tolerance")
      )
      isTRUE(.dpprior_v2_exact_field(record, "passed")) &&
        !is.null(raw) && !is.null(tolerance) && all(tolerance >= 0) &&
        .dpprior_v2_pairs_close(raw, fresh_pair - target_pair,
                                abs_tol = 1e-12, rel_tol = 1e-10) &&
        all(abs(raw) <= tolerance)
    }
    residual_record_pass(selected_residual, fresh_selected) &&
      residual_record_pass(refined_residual, fresh_refined) &&
      isTRUE(.dpprior_v2_exact_field(verification, "available")) &&
      identical(.dpprior_v2_exact_field(verification, "status"),
                "converged")
  } else {
    KL_evidence <- .dpprior_v2_validate_A2_KL_input_evidence(
      raw_fit = raw_fit, verification = verification, wrapper = wrapper,
      a = a, b = b, J = J, target = target,
      retained_target_identity = retained_target_identity,
      M_selected = M_selected,
      M_refined = M_refined
    )
    isTRUE(KL_evidence$passed)
  }

  passed <- recorded_pass && wrapper_pass && stability_pass && target_pass &&
    method_evidence_pass
  if (!passed) {
    .dpprior_v2_abort_invalid(
      paste(
        "input fit verification evidence does not agree with fresh",
        "selected/refined recomputation of the returned parameters"
      ),
      "fit verification evidence",
      list(
        method_family = family, selected_order = M_selected,
        refined_order = M_refined, fresh_selected = fresh_selected,
        fresh_refined = fresh_refined, target = target_pair,
        recorded_pass = recorded_pass, wrapper_pass = wrapper_pass,
        stability_pass = stability_pass, target_pass = target_pass,
        method_evidence_pass = method_evidence_pass,
        method_evidence_reason = KL_evidence$reason %||% NA_character_
      ),
      "method-appropriate evidence bound to returned a, b, target, and orders",
      "input_fit_reverification_failed",
      c("dpprior_dual_fit_error", "dpprior_status_contract_error")
    )
  }

  list(
    required = TRUE, performed = TRUE, passed = TRUE,
    method_family = family, M_selected = M_selected,
    M_refined = M_refined, selected = fresh_selected,
    refined = fresh_refined, target_tolerance = target_tolerance,
    stability_tolerance = stability_tolerance,
    method_evidence_reason = KL_evidence$reason %||%
      "fresh_A2_MN_residual_evidence_passed",
    reason = "fresh_method_appropriate_input_fit_verification_passed"
  )
}


# Canonical consumer helpers deliberately unclass only after the shared schema
# validator has accepted the object.  This prevents S3 `$`/`[[` methods and
# compatibility aliases from becoming scientific input authority.
.dpprior_v2_canonical_record <- function(x, path) {
  if (!is.list(x)) {
    .dpprior_v2_abort_invalid(
      sprintf("%s must be a canonical list record", path), path, x,
      "a canonical ordinary list record", "canonical_record_type",
      c("dpprior_dual_fit_error", "dpprior_type_error")
    )
  }
  raw <- unclass(x)
  nms <- names(raw)
  if (is.null(nms) || anyDuplicated(nms)) {
    .dpprior_v2_abort_invalid(
      sprintf("%s must have unique exact field names", path), path, nms,
      "unique exact canonical field names", "canonical_record_names",
      c("dpprior_dual_fit_error", "dpprior_dual_name_error")
    )
  }
  raw
}


.dpprior_v2_canonical_target_identity <- function(target_K, J, method) {
  target_K <- .dpprior_require_schema(
    target_K, kind = "target", allow_legacy = FALSE
  )
  raw <- .dpprior_v2_canonical_record(target_K, "fit$target$K")
  if (!identical(raw[["J", exact = TRUE]], as.integer(J))) {
    .dpprior_v2_abort_invalid(
      "fit target J is not identical to fit J", "fit$target$K$J",
      raw[["J", exact = TRUE]], as.integer(J),
      "conflicting_retained_target_identity",
      c("dpprior_dual_fit_error", "dpprior_conflicting_input")
    )
  }

  implied <- .dpprior_v2_canonical_record(
    raw[["implied", exact = TRUE]], "fit$target$K$implied"
  )
  request <- .dpprior_v2_canonical_record(
    raw[["request", exact = TRUE]], "fit$target$K$request"
  )
  used <- .dpprior_v2_canonical_record(
    raw[["used", exact = TRUE]], "fit$target$K$used"
  )
  kind <- raw[["kind", exact = TRUE]]
  route <- if (identical(kind, "pmf")) {
    "strict_pmf"
  } else if (identical(kind, "interval")) {
    "interval"
  } else if (identical(kind, "family")) {
    "family"
  } else if (identical(kind, "cv")) {
    "coefficient_of_variation"
  } else if ("confidence" %in% names(request)) {
    "qualitative_confidence"
  } else {
    "direct_variance"
  }
  target_kind <- if (kind %in% c("pmf", "interval", "family")) {
    "pmf"
  } else {
    "moments"
  }
  authoritative_pmf <- raw[["pmf", exact = TRUE]]

  list(
    available = TRUE,
    passed = TRUE,
    reason = "canonical_target_schema_and_derivation_validated",
    route = route,
    target_kind = target_kind,
    kind = kind,
    method = method,
    request = request,
    used = used,
    implied = implied,
    authoritative_pmf = authoritative_pmf
  )
}


.dpprior_v2_canonical_input_fit_reference <- function(raw_fit, target_K) {
  target_K <- .dpprior_require_schema(
    target_K, kind = "target", allow_legacy = FALSE
  )
  target_raw <- .dpprior_v2_canonical_record(target_K, "fit$target$K")
  parameters <- raw_fit[["parameters", exact = TRUE]]
  verification <- .dpprior_v2_canonical_record(
    raw_fit[["verification", exact = TRUE]], "fit$verification"
  )
  selected <- verification[["selected_snapshot", exact = TRUE]]
  verifier <- verification[["verifier_snapshot", exact = TRUE]]
  if (is.null(parameters) || is.null(selected) || is.null(verifier)) {
    return(NULL)
  }
  parameters <- .dpprior_v2_canonical_record(parameters, "fit$parameters")
  target_schema <- .dpprior_v2_canonical_record(
    target_raw[["schema", exact = TRUE]], "fit$target$K$schema"
  )
  target_used <- .dpprior_v2_canonical_record(
    target_raw[["used", exact = TRUE]], "fit$target$K$used"
  )
  target_implied <- .dpprior_v2_canonical_record(
    target_raw[["implied", exact = TRUE]], "fit$target$K$implied"
  )

  snapshot_reference <- function(snapshot, path) {
    snapshot <- .dpprior_v2_canonical_record(snapshot, path)
    achieved <- .dpprior_v2_canonical_record(
      snapshot[["achieved", exact = TRUE]], paste0(path, "$achieved")
    )
    K <- .dpprior_v2_canonical_record(
      achieved[["K", exact = TRUE]], paste0(path, "$achieved$K")
    )
    list(
      parameters = .dpprior_v2_canonical_record(
        snapshot[["parameters", exact = TRUE]], paste0(path, "$parameters")
      ),
      M = snapshot[["M", exact = TRUE]],
      achieved_K = K,
      finite = snapshot[["finite", exact = TRUE]],
      source = snapshot[["source", exact = TRUE]]
    )
  }

  decision_evidence <- if (identical(
    raw_fit[["mode", exact = TRUE]], "a2_kl"
  )) {
    tolerances <- .dpprior_v2_canonical_record(
      raw_fit[["tolerances", exact = TRUE]], "fit$tolerances"
    )
    distribution <- .dpprior_v2_canonical_record(
      tolerances[["distribution", exact = TRUE]],
      "fit$tolerances$distribution"
    )
    list(
      target_K = target_K,
      distribution_tolerances = distribution
    )
  } else {
    NULL
  }

  list(
    schema = "dpprior.result/1",
    mode = raw_fit[["mode", exact = TRUE]],
    method = raw_fit[["method", exact = TRUE]],
    J = raw_fit[["J", exact = TRUE]],
    status = raw_fit[["status", exact = TRUE]],
    usable = raw_fit[["usable", exact = TRUE]],
    verified = raw_fit[["verified", exact = TRUE]],
    parameters = parameters,
    target = list(
      schema = target_schema,
      kind = target_raw[["kind", exact = TRUE]],
      J = target_raw[["J", exact = TRUE]],
      used = target_used,
      implied = target_implied
    ),
    decision_evidence = decision_evidence,
    selected_snapshot = snapshot_reference(
      selected, "fit$verification$selected_snapshot"
    ),
    verifier_snapshot = snapshot_reference(
      verifier, "fit$verification$verifier_snapshot"
    )
  )
}


.dpprior_v2_canonical_input_fit_evidence <- function(
    raw_fit, a, b, J, target, retained_target_identity) {
  verified <- raw_fit[["verified", exact = TRUE]]
  if (!isTRUE(verified)) {
    return(list(
      required = FALSE, performed = FALSE, passed = FALSE,
      reason = "input_fit_does_not_make_a_verified_claim"
    ))
  }

  computation <- .dpprior_v2_canonical_record(
    raw_fit[["computation", exact = TRUE]], "fit$computation"
  )
  orders <- .dpprior_v2_canonical_record(
    computation[["orders", exact = TRUE]], "fit$computation$orders"
  )
  verification <- .dpprior_v2_canonical_record(
    raw_fit[["verification", exact = TRUE]], "fit$verification"
  )
  selected <- .dpprior_v2_canonical_record(
    verification[["selected_snapshot", exact = TRUE]],
    "fit$verification$selected_snapshot"
  )
  verifier <- .dpprior_v2_canonical_record(
    verification[["verifier_snapshot", exact = TRUE]],
    "fit$verification$verifier_snapshot"
  )
  M_selected <- orders[["M_selected", exact = TRUE]]
  M_refined <- orders[["M_verification_used", exact = TRUE]]
  required <- orders[["M_verification_required", exact = TRUE]]

  selected_achieved <- .dpprior_v2_canonical_record(
    .dpprior_v2_canonical_record(
      selected[["achieved", exact = TRUE]],
      "fit$verification$selected_snapshot$achieved"
    )[["K", exact = TRUE]],
    "fit$verification$selected_snapshot$achieved$K"
  )
  refined_achieved <- .dpprior_v2_canonical_record(
    .dpprior_v2_canonical_record(
      verifier[["achieved", exact = TRUE]],
      "fit$verification$verifier_snapshot$achieved"
    )[["K", exact = TRUE]],
    "fit$verification$verifier_snapshot$achieved$K"
  )
  fresh <- tryCatch(
    list(
      selected = exact_K_moments(J, a, b, M = M_selected),
      refined = exact_K_moments(J, a, b, M = M_refined)
    ),
    error = function(condition) condition
  )
  if (inherits(fresh, "condition")) {
    .dpprior_v2_abort_invalid(
      "independent canonical input-fit recomputation failed",
      "fit$parameters", list(a = a, b = b, J = J),
      "finite selected/refined K recomputation",
      "input_fit_recomputation_failed",
      c("dpprior_dual_fit_error", "dpprior_status_contract_error")
    )
  }
  fresh_selected <- c(mean = fresh[["selected"]][["mean"]],
                      variance = fresh[["selected"]][["var"]])
  fresh_refined <- c(mean = fresh[["refined"]][["mean"]],
                     variance = fresh[["refined"]][["var"]])
  recorded_selected <- c(
    mean = selected_achieved[["mean", exact = TRUE]],
    variance = selected_achieved[["variance", exact = TRUE]]
  )
  recorded_refined <- c(
    mean = refined_achieved[["mean", exact = TRUE]],
    variance = refined_achieved[["variance", exact = TRUE]]
  )
  target_pair <- c(mean = target[["mu_K"]], variance = target[["var_K"]])
  target_tolerance <- 1e-8 + 1e-8 * pmax(abs(target_pair), 1)
  stability_tolerance <- 1e-10 + 1e-8 * pmax(
    abs(fresh_selected), abs(fresh_refined), 1
  )
  required_formula <- as.integer(max(
    2L * as.integer(M_selected), as.integer(M_selected) + 40L
  ))
  target_pass <- identical(raw_fit[["mode", exact = TRUE]], "a2_kl") ||
    (all(abs(fresh_selected - target_pair) <= target_tolerance) &&
       all(abs(fresh_refined - target_pair) <= target_tolerance))
  passed <- isTRUE(verification[["performed", exact = TRUE]]) &&
    isTRUE(verification[["passed", exact = TRUE]]) &&
    identical(selected[["M", exact = TRUE]], M_selected) &&
    identical(verifier[["M", exact = TRUE]], M_refined) &&
    identical(required, required_formula) &&
    .dpprior_v2_plain_positive_integer(
      M_refined, maximum = .QUADRATURE_MAX_NODES
    ) && M_refined >= required &&
    .dpprior_v2_pairs_close(fresh_selected, recorded_selected) &&
    .dpprior_v2_pairs_close(fresh_refined, recorded_refined) &&
    target_pass
  if (!passed) {
    .dpprior_v2_abort_invalid(
      paste(
        "canonical input-fit verification does not agree with fresh",
        "selected/refined recomputation"
      ),
      "fit$verification", list(
        selected_order = M_selected, refined_order = M_refined,
        required_order = required, fresh_selected = fresh_selected,
        fresh_refined = fresh_refined, target = target_pair,
        retained_target_identity = retained_target_identity
      ),
      "validated canonical evidence bound to parameters, target, and orders",
      "input_fit_reverification_failed",
      c("dpprior_dual_fit_error", "dpprior_status_contract_error")
    )
  }

  list(
    required = TRUE, performed = TRUE, passed = TRUE,
    method_family = raw_fit[["mode", exact = TRUE]],
    M_selected = M_selected, M_refined = M_refined,
    selected = fresh_selected, refined = fresh_refined,
    target_tolerance = target_tolerance,
    stability_tolerance = stability_tolerance,
    method_evidence_reason = "canonical_schema_and_fresh_snapshots_passed",
    reason = "fresh_canonical_input_fit_verification_passed"
  )
}


# Normalize the K-only fit without mutating or relabeling the input object.
.dpprior_v2_normalize_fit <- function(fit) {
  fit <- .dpprior_require_schema(fit, kind = "fit", allow_legacy = FALSE)
  fit_class <- class(fit)
  raw_fit <- .dpprior_v2_canonical_record(fit, "fit")
  fit_J <- raw_fit[["J", exact = TRUE]]
  fit_status <- raw_fit[["status", exact = TRUE]]
  fit_usable <- raw_fit[["usable", exact = TRUE]]
  fit_verified <- raw_fit[["verified", exact = TRUE]]
  fit_method <- raw_fit[["method", exact = TRUE]]
  fit_mode <- raw_fit[["mode", exact = TRUE]]
  J <- .dpprior_validate_count(
    fit_J, "fit$J", minimum = 2L, maximum = .MAX_J_DEFAULT,
    .subclass = "dpprior_dual_fit_error"
  )

  target_bundle <- .dpprior_v2_canonical_record(
    raw_fit[["target", exact = TRUE]], "fit$target"
  )
  fit_target <- target_bundle[["K", exact = TRUE]]
  target_raw <- .dpprior_v2_canonical_record(fit_target, "fit$target$K")
  implied <- .dpprior_v2_canonical_record(
    target_raw[["implied", exact = TRUE]], "fit$target$K$implied"
  )
  mu_K <- .dpprior_validate_scalar(
    implied[["mean", exact = TRUE]], "fit$target$K$implied$mean",
    lower = 1, upper = J, .subclass = "dpprior_dual_fit_error"
  )
  var_K <- .dpprior_validate_scalar(
    implied[["variance", exact = TRUE]],
    "fit$target$K$implied$variance", lower = 0,
    .subclass = "dpprior_dual_fit_error"
  )
  retained_target_identity <- .dpprior_v2_canonical_target_identity(
    fit_target, J, fit_method
  )

  parameters <- raw_fit[["parameters", exact = TRUE]]
  a <- b <- NULL
  if (!is.null(parameters)) {
    parameters <- .dpprior_v2_canonical_record(parameters, "fit$parameters")
    a <- .dpprior_validate_scalar(
      parameters[["a", exact = TRUE]], "fit$parameters$a",
      lower = 0, lower_open = TRUE,
      .subclass = "dpprior_dual_fit_error"
    )
    b <- .dpprior_validate_scalar(
      parameters[["b", exact = TRUE]], "fit$parameters$b",
      lower = 0, lower_open = TRUE,
      .subclass = "dpprior_dual_fit_error"
    )
  }

  input_fit_evidence <- if (is.null(parameters)) {
    list(
      required = isTRUE(fit_verified), performed = FALSE, passed = FALSE,
      reason = "input_fit_has_no_public_parameters"
    )
  } else {
    .dpprior_v2_canonical_input_fit_evidence(
      raw_fit = raw_fit, a = a, b = b, J = J,
      target = list(mu_K = mu_K, var_K = var_K),
      retained_target_identity = retained_target_identity
    )
  }
  input_fit_reference <- .dpprior_v2_canonical_input_fit_reference(
    raw_fit, fit_target
  )

  list(
    a = a,
    b = b,
    J = as.integer(J),
    target_K = list(
      mu_K = mu_K,
      var_K = var_K,
      units = list(mu_K = "clusters", var_K = "clusters^2"),
      canonical = fit_target,
      source = "fit$target$K"
    ),
    status = fit_status,
    usable = fit_usable,
    verified = fit_verified,
    source = list(
      schema = "dpprior.result/1",
      class = fit_class,
      mode = fit_mode,
      method = fit_method,
      status = fit_status,
      usable = fit_usable,
      verified = fit_verified,
      target_source = "fit$target$K",
      retained_target_identity = retained_target_identity[
        setdiff(names(retained_target_identity), "authoritative_pmf")
      ],
      input_fit_evidence = input_fit_evidence,
      canonical_reference = input_fit_reference
    )
  )
}


.dpprior_v2_normalize_weight_spec <- function(
    spec, mode = c("hard", "soft")) {
  mode <- match.arg(mode)
  required_value <- if (mode == "hard") "bound" else "value"
  allowed <- c(
    "metric", "relation", required_value, "threshold", "probability", "prob"
  )
  spec <- .dpprior_v2_validate_named_list(
    spec, if (mode == "hard") "constraint" else "target",
    allowed = allowed, required = c("metric", "relation", required_value)
  )

  if (!.dpprior_v2_plain_scalar(spec$metric, "character")) {
    .dpprior_v2_abort_invalid(
      "weight metric must be one named character scalar", "metric",
      spec$metric, "named metric", "type",
      c("dpprior_dual_metric_error", "dpprior_type_error")
    )
  }
  metrics <- c(
    "wsb_tail", "wsb_mean", "wsb_quantile", "wmax_tail_upper"
  )
  if (!(spec$metric %in% metrics)) {
    .dpprior_v2_abort_invalid(
      sprintf("metric must be one of: %s", paste(metrics, collapse = ", ")),
      "metric", spec$metric, paste(metrics, collapse = ", "),
      "unsupported_metric", "dpprior_dual_metric_error"
    )
  }
  if (!.dpprior_v2_plain_scalar(spec$relation, "character")) {
    .dpprior_v2_abort_invalid(
      "relation must be one character scalar", "relation", spec$relation,
      "a supported named relation", "type",
      c("dpprior_dual_relation_error", "dpprior_type_error")
    )
  }
  relation_map <- c(
    "<=" = "at_most", "at_most" = "at_most",
    ">=" = "at_least", "at_least" = "at_least",
    "target" = "target"
  )
  relation <- unname(relation_map[spec$relation])
  permitted <- if (mode == "hard") c("at_most", "at_least") else {
    c("target", "at_most", "at_least")
  }
  if (length(relation) != 1L || is.na(relation) || !(relation %in% permitted)) {
    .dpprior_v2_abort_invalid(
      sprintf("relation must be one of: %s", paste(permitted, collapse = ", ")),
      "relation", spec$relation, paste(permitted, collapse = ", "),
      "unsupported_relation", "dpprior_dual_relation_error"
    )
  }

  raw_value <- spec[[required_value]]
  value <- .dpprior_validate_probability(
    raw_value, required_value, scalar = TRUE, open = FALSE,
    .subclass = "dpprior_dual_weight_value_error"
  )
  threshold <- probability <- NULL

  if (spec$metric %in% c("wsb_tail", "wmax_tail_upper")) {
    if (is.null(spec$threshold)) {
      .dpprior_v2_abort_invalid(
        sprintf("%s requires threshold", spec$metric), "threshold", NULL,
        "one probability in [0,1]", "missing_threshold",
        c("dpprior_dual_metric_error", "dpprior_missing_error")
      )
    }
    threshold <- .dpprior_validate_probability(
      spec$threshold, "threshold", scalar = TRUE, open = FALSE,
      .subclass = "dpprior_dual_metric_error"
    )
  } else if (!is.null(spec$threshold)) {
    .dpprior_v2_abort_invalid(
      sprintf("threshold is not defined for metric '%s'", spec$metric),
      "threshold", spec$threshold, "omit threshold", "irrelevant_component",
      "dpprior_dual_metric_error"
    )
  }

  if (identical(spec$metric, "wsb_quantile")) {
    probability_component <- spec[["probability", exact = TRUE]]
    prob_component <- spec[["prob", exact = TRUE]]
    if (!is.null(probability_component) && !is.null(prob_component)) {
      .dpprior_v2_abort_invalid(
        "supply only one of probability and prob", "probability/prob",
        list(probability = probability_component, prob = prob_component),
        "exactly one quantile probability", "conflicting_probability",
        c("dpprior_dual_metric_error", "dpprior_conflicting_input")
      )
    }
    probability_value <- probability_component %||% prob_component
    if (is.null(probability_value)) {
      .dpprior_v2_abort_invalid(
        "wsb_quantile requires probability", "probability", NULL,
        "one probability in [0,1]", "missing_probability",
        c("dpprior_dual_metric_error", "dpprior_missing_error")
      )
    }
    probability <- .dpprior_validate_probability(
      probability_value, "probability", scalar = TRUE, open = FALSE,
      .subclass = "dpprior_dual_metric_error"
    )
  } else if (!is.null(spec[["probability", exact = TRUE]]) ||
             !is.null(spec[["prob", exact = TRUE]])) {
    .dpprior_v2_abort_invalid(
      sprintf("probability is not defined for metric '%s'", spec$metric),
      "probability", spec[["probability", exact = TRUE]] %||%
        spec[["prob", exact = TRUE]],
      "omit probability", "irrelevant_component", "dpprior_dual_metric_error"
    )
  }

  if (mode == "hard" && identical(spec$metric, "wmax_tail_upper") &&
      identical(relation, "at_least")) {
    .dpprior_v2_abort_invalid(
      paste(
        "wmax_tail_upper supports only an at_most hard constraint; an upper",
        "bound cannot certify a lower-tail requirement"
      ),
      "relation", spec$relation, "at_most or <=", "unsafe_relation",
      "dpprior_dual_metric_error"
    )
  }

  estimand <- switch(
    spec$metric,
    wsb_tail = "P(W_SB > threshold)",
    wsb_mean = "E(W_SB)",
    wsb_quantile = "Q_probability(W_SB)",
    wmax_tail_upper = "certified upper bound for P(W_max > threshold)"
  )
  out <- list(
    metric = spec$metric,
    relation = relation,
    operator = switch(relation, at_most = "<=", at_least = ">=",
                      target = "target"),
    value = value,
    threshold = threshold,
    probability = probability,
    estimand = estimand,
    units = "probability",
    raw_input = spec
  )
  if (mode == "hard") out$bound <- value
  out
}


.dpprior_v2_failed_metric <- function(spec, message, method = NA_character_) {
  list(
    metric = spec$metric,
    estimand = spec$estimand,
    value = NA_real_,
    units = spec$units,
    threshold = spec$threshold,
    probability = spec$probability,
    status = "failed",
    usable = FALSE,
    verified = FALSE,
    method = method,
    error_bound = NA_real_,
    certification = list(certified = FALSE, reason = message),
    provenance = list(error = message)
  )
}


# Evaluate only named estimands. In particular, wmax_tail_upper is the
# certified endpoint, never a direct W_max point estimate in disguise.
.dpprior_v2_eval_metric <- function(spec, a, b, J = NULL,
                                    M = .QUAD_NODES_DEFAULT) {
  if (!is.list(spec) || is.null(spec$metric) || is.null(spec$estimand)) {
    .dpprior_v2_abort_invalid(
      "spec must be a normalized weight specification", "spec", spec,
      ".dpprior_v2_normalize_weight_spec() output", "unnormalized_spec",
      "dpprior_dual_metric_error"
    )
  }
  a <- .dpprior_validate_scalar(
    a, "a", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_dual_parameter_error"
  )
  b <- .dpprior_validate_scalar(
    b, "b", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_dual_parameter_error"
  )
  M <- .dpprior_validate_count(
    M, "M", minimum = 1L, maximum = .QUADRATURE_MAX_NODES,
    .subclass = "dpprior_dual_control_error"
  )

  evaluated <- tryCatch({
    switch(
      spec$metric,
      wsb_tail = list(
        value = prob_wsb_exceeds(spec$threshold, a, b),
        method = "closed_form_beta_gamma_survival",
        independently_certified = TRUE,
        details = list(formula = "(b/(b-log(1-threshold)))^a")
      ),
      wsb_mean = list(
        value = mean_w1(a, b, M = M),
        method = "generalized_gauss_laguerre",
        independently_certified = FALSE,
        details = list(M = as.integer(M))
      ),
      wsb_quantile = list(
        value = quantile_w1(spec$probability, a, b),
        method = "closed_form_beta_gamma_quantile",
        independently_certified = TRUE,
        details = list(probability = spec$probability)
      ),
      wmax_tail_upper = {
        bounds <- wmax_tail_bounds(spec$threshold, a = a, b = b)
        list(
          value = bounds$upper_bound,
          method = "certified_size_biased_mass_upper_bound",
          independently_certified = isTRUE(bounds$certified),
          details = list(
            lower_bound = bounds$lower_bound,
            upper_bound = bounds$upper_bound,
            log_upper_bound = bounds$log_upper_bound,
            size_biased_tail = bounds$size_biased_tail,
            bounds_method = bounds$method,
            outward_rounded = bounds$upper_outward_rounded,
            probability_scale_underflow =
              bounds$probability_scale_underflow,
            upper_scale_underflow = bounds$upper_scale_underflow
          )
        )
      },
      stop("unsupported normalized metric")
    )
  }, error = identity)
  if (inherits(evaluated, "condition")) {
    return(.dpprior_v2_failed_metric(
      spec, conditionMessage(evaluated), "metric_evaluation"
    ))
  }

  value <- evaluated$value
  finite_probability <- .dpprior_v2_plain_scalar(value, "numeric") &&
    value >= 0 && value <= 1
  if (!finite_probability) {
    return(.dpprior_v2_failed_metric(
      spec,
      sprintf("metric evaluation returned value outside [0,1]: %s",
              paste(value, collapse = ", ")),
      evaluated$method
    ))
  }
  rounding_bound <- 64 * .Machine$double.eps * max(1, abs(value))
  certified <- isTRUE(evaluated$independently_certified)
  list(
    metric = spec$metric,
    estimand = spec$estimand,
    value = as.numeric(value),
    units = spec$units,
    threshold = spec$threshold,
    probability = spec$probability,
    status = if (certified) "converged" else "approximate",
    usable = certified,
    verified = certified,
    method = evaluated$method,
    error_bound = rounding_bound,
    certification = list(
      certified = certified,
      reason = if (certified) {
        "closed_form_or_certified_bound"
      } else {
        "fixed_order_requires_refined_order_comparison"
      }
    ),
    provenance = c(
      list(
        gamma_parameterization = "shape_rate",
        a = a,
        b = b,
        M = as.integer(M)
      ),
      evaluated$details
    )
  )
}


.dpprior_v2_extract_K <- function(x, name) {
  if (!is.list(x)) {
    .dpprior_v2_abort_invalid(
      sprintf("%s must be a list", name), name, x,
      "list containing K mean and variance", "type", "dpprior_dual_K_error"
    )
  }
  mu <- x$mu_K %||% x$mean
  variance <- x$var_K %||% x$variance %||% x$var
  mu <- .dpprior_validate_scalar(
    mu, paste0(name, "$mu_K"), .subclass = "dpprior_dual_K_error"
  )
  variance <- .dpprior_validate_scalar(
    variance, paste0(name, "$var_K"), lower = 0,
    .subclass = "dpprior_dual_K_error"
  )
  c(mu_K = mu, var_K = variance)
}


.dpprior_v2_k_loss <- function(achieved, target, scales = NULL) {
  achieved_values <- .dpprior_v2_extract_K(achieved, "achieved")
  target_values <- .dpprior_v2_extract_K(target, "target")
  if (is.null(scales)) {
    scales <- c(
      mu_K = max(abs(target_values[["mu_K"]]), 1),
      var_K = max(abs(target_values[["var_K"]]), 1)
    )
  } else {
    if (!is.numeric(scales) || is.object(scales) || !is.null(dim(scales)) ||
        length(scales) != 2L || anyNA(scales) || any(!is.finite(scales)) ||
        any(scales <= 0)) {
      .dpprior_v2_abort_invalid(
        "scales must contain two finite positive numeric values", "scales",
        scales, "positive c(mu_K, var_K)", "invalid_scales",
        "dpprior_dual_K_error"
      )
    }
    if (!is.null(names(scales)) && all(c("mu_K", "var_K") %in% names(scales))) {
      scales <- scales[c("mu_K", "var_K")]
    } else {
      names(scales) <- c("mu_K", "var_K")
    }
  }
  raw <- achieved_values - target_values
  scaled <- raw / scales
  list(
    value = as.numeric(sum(scaled^2)),
    raw = raw,
    scaled = scaled,
    scales = scales,
    scale_formula = c(
      mu_K = "max(abs(target_mu_K), 1)",
      var_K = "max(abs(target_var_K), 1)"
    )
  )
}


.dpprior_v2_constraint_values <- function(spec, achieved, tolerance) {
  value <- if (is.list(achieved)) achieved$value else achieved
  if (!.dpprior_v2_plain_scalar(value, "numeric")) {
    return(list(
      raw = NA_real_, slack = NA_real_, tolerance = tolerance,
      satisfied = FALSE
    ))
  }
  raw <- switch(
    spec$relation,
    at_most = value - spec$value,
    at_least = spec$value - value,
    .dpprior_v2_abort_invalid(
      "hard residual requires at_most or at_least", "relation",
      spec$relation, "at_most or at_least", "unsupported_relation",
      "dpprior_dual_relation_error"
    )
  )
  list(
    raw = as.numeric(raw),
    slack = as.numeric(-raw),
    tolerance = as.numeric(tolerance),
    satisfied = is.finite(raw) && raw <= tolerance
  )
}


.dpprior_v2_make_attempt <- function(
    method, start = NULL, bounds = NULL, control = list(),
    exit_code = NA_integer_, message = NA_character_,
    counts = list(iterations = NA_integer_, evaluations = NA_integer_),
    objective = NA_real_, elapsed = NA_real_, warning = NA_character_,
    error = NA_character_, candidate = NULL) {
  if (!is.list(counts)) counts <- list(evaluations = counts)
  iterations <- counts$iterations %||% NA_integer_
  evaluations <- counts$evaluations %||% counts[["function"]] %||% NA_integer_
  normalized_exit <- if (
    .dpprior_v2_plain_scalar(exit_code, "numeric") &&
    exit_code == as.integer(exit_code)
  ) {
    as.integer(exit_code)
  } else {
    NA_integer_
  }
  list(
    method = as.character(method),
    start = start,
    bounds = bounds,
    control = control,
    exit_code = normalized_exit,
    exit = normalized_exit,
    message = if (length(message) == 1L) as.character(message) else NA_character_,
    iterations = if (length(iterations) == 1L) as.integer(iterations) else NA_integer_,
    evaluations = if (length(evaluations) == 1L) as.integer(evaluations) else NA_integer_,
    counts = list(
      iterations = if (length(iterations) == 1L) as.integer(iterations) else NA_integer_,
      evaluations = if (length(evaluations) == 1L) as.integer(evaluations) else NA_integer_
    ),
    candidate_objective = if (length(objective) == 1L) {
      as.numeric(objective)
    } else {
      NA_real_
    },
    objective = if (length(objective) == 1L) as.numeric(objective) else NA_real_,
    elapsed = if (length(elapsed) == 1L) as.numeric(elapsed) else NA_real_,
    warning = if (length(warning) == 1L) as.character(warning) else NA_character_,
    error = if (length(error) == 1L) as.character(error) else NA_character_,
    candidate = candidate
  )
}


.dpprior_v2_capture <- function(expr) {
  warnings <- character()
  started <- proc.time()[["elapsed"]]
  value <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = identity
  )
  elapsed <- proc.time()[["elapsed"]] - started
  list(
    value = value,
    warning = if (length(warnings)) paste(unique(warnings), collapse = " | ") else NA_character_,
    error = if (inherits(value, "condition")) conditionMessage(value) else NA_character_,
    elapsed = as.numeric(elapsed)
  )
}


.dpprior_v2_normalize_hard_controls <- function(
    constraint_tol, M, M_verify, log_bounds, control, allow_approximate) {
  constraint_tol <- .dpprior_v2_validate_named_list(
    constraint_tol, "constraint_tol", allowed = c("abs", "rel")
  )
  tol_abs <- .dpprior_validate_scalar(
    constraint_tol$abs %||% 1e-6, "constraint_tol$abs", lower = 0,
    .subclass = "dpprior_dual_control_error"
  )
  tol_rel <- .dpprior_validate_scalar(
    constraint_tol$rel %||% 1e-6, "constraint_tol$rel", lower = 0,
    .subclass = "dpprior_dual_control_error"
  )
  M <- .dpprior_validate_count(
    M, "M", minimum = 10L,
    maximum = as.integer(floor(.QUADRATURE_MAX_NODES / 2)),
    .subclass = "dpprior_dual_control_error"
  )
  M_required <- .quadrature_verification_required_order(M)
  if (is.null(M_verify)) M_verify <- M_required
  M_verify <- .dpprior_validate_count(
    M_verify, "M_verify", minimum = M_required,
    maximum = .QUADRATURE_MAX_NODES,
    .subclass = "dpprior_dual_verification_error"
  )
  if (!is.numeric(log_bounds) || is.object(log_bounds) ||
      !is.null(dim(log_bounds)) || length(log_bounds) != 2L ||
      anyNA(log_bounds) || any(!is.finite(log_bounds)) ||
      log_bounds[[1L]] >= log_bounds[[2L]] ||
      log_bounds[[1L]] < -.EXP_MAX || log_bounds[[2L]] > .EXP_MAX) {
    .dpprior_v2_abort_invalid(
      "log_bounds must be two finite increasing values within the safe exponential range",
      "log_bounds", log_bounds,
      sprintf("c(lower, upper), -%g <= lower < upper <= %g", .EXP_MAX, .EXP_MAX),
      "invalid_bounds", c("dpprior_dual_control_error", "dpprior_bounds_error")
    )
  }
  allow_approximate <- .dpprior_validate_control(
    allow_approximate, "allow_approximate", type = "logical"
  )

  defaults <- list(
    maxit = 250L,
    scan_points = 17L,
    scan_keep = 12L,
    profile_starts = 4L,
    root_tol = 1e-10,
    optim_reltol = 1e-10,
    penalty = 1e7,
    boundary_tol = 1e-5,
    verification_abs_tol = 1e-8,
    verification_rel_tol = 1e-6,
    perturbation_step = 1e-6,
    perturbation_abs_tol = 1e-4
  )
  control <- .dpprior_v2_validate_named_list(
    control, "control", allowed = names(defaults)
  )
  supplied_control_names <- names(control)
  control <- utils::modifyList(defaults, control, keep.null = FALSE)
  control$maxit <- .dpprior_validate_count(
    control$maxit, "control$maxit", minimum = 10L, maximum = 100000L,
    .subclass = "dpprior_dual_control_error"
  )
  control$scan_points <- .dpprior_validate_count(
    control$scan_points, "control$scan_points", minimum = 5L, maximum = 101L,
    .subclass = "dpprior_dual_control_error"
  )
  control$scan_keep <- .dpprior_validate_count(
    control$scan_keep, "control$scan_keep", minimum = 1L, maximum = 101L,
    .subclass = "dpprior_dual_control_error"
  )
  control$profile_starts <- .dpprior_validate_count(
    control$profile_starts, "control$profile_starts", minimum = 1L,
    maximum = 20L, .subclass = "dpprior_dual_control_error"
  )
  for (nm in c(
    "root_tol", "optim_reltol", "penalty", "boundary_tol",
    "verification_abs_tol", "verification_rel_tol", "perturbation_step",
    "perturbation_abs_tol"
  )) {
    control[[nm]] <- .dpprior_validate_scalar(
      control[[nm]], paste0("control$", nm), lower = 0,
      lower_open = nm %in% c(
        "root_tol", "optim_reltol", "penalty", "perturbation_step"
      ),
      .subclass = "dpprior_dual_control_error"
    )
  }
  auto_cap <- list(
    scan_keep = control$scan_keep > control$scan_points &&
      !("scan_keep" %in% supplied_control_names),
    profile_starts = control$profile_starts > control$scan_points &&
      !("profile_starts" %in% supplied_control_names)
  )
  if (auto_cap$scan_keep) {
    control$scan_keep <- control$scan_points
  }
  if (auto_cap$profile_starts) {
    control$profile_starts <- control$scan_points
  }
  if (control$scan_keep > control$scan_points ||
      control$profile_starts > control$scan_points) {
    .dpprior_v2_abort_invalid(
      "scan_keep and profile_starts cannot exceed scan_points",
      "control",
      control[c("scan_points", "scan_keep", "profile_starts")],
      "scan_keep <= scan_points and profile_starts <= scan_points",
      "inconsistent_scan_controls",
      c("dpprior_dual_control_error", "dpprior_status_contract_error")
    )
  }
  central_upper <- c(
    root_tol = defaults$root_tol,
    optim_reltol = defaults$optim_reltol,
    verification_abs_tol = defaults$verification_abs_tol,
    verification_rel_tol = defaults$verification_rel_tol,
    perturbation_abs_tol = defaults$perturbation_abs_tol
  )
  inflated <- names(central_upper)[vapply(
    names(central_upper),
    function(nm) control[[nm]] > central_upper[[nm]], logical(1)
  )]
  fixed_control_changed <- c(
    boundary_tol = !identical(
      control$boundary_tol, defaults$boundary_tol
    ),
    perturbation_step = !identical(
      control$perturbation_step, defaults$perturbation_step
    )
  )
  if (length(inflated) || any(fixed_control_changed) ||
      control$penalty > 1e12) {
    fixed_names <- names(fixed_control_changed)[fixed_control_changed]
    .dpprior_v2_abort_invalid(
      paste(
        "internal numerical verification controls may only be made stricter;",
        "relaxed exploratory controls cannot mint a verified hard result"
      ),
      "control",
      control[unique(c(
        inflated, fixed_names,
        if (control$penalty > 1e12) "penalty"
      ))],
      paste(
        "root_tol <= 1e-10, optim_reltol <= 1e-10,",
        "verification_abs_tol <= 1e-8, verification_rel_tol <= 1e-6,",
        "perturbation_abs_tol <= 1e-4, boundary_tol = 1e-5,",
        "perturbation_step = 1e-6, penalty <= 1e12"
      ),
      "noncentral_verification_tolerance",
      c("dpprior_dual_control_error", "dpprior_status_contract_error")
    )
  }
  list(
    tolerance = list(abs = tol_abs, rel = tol_rel),
    M = as.integer(M),
    M_verify = as.integer(M_verify),
    M_verify_required = as.integer(M_required),
    log_bounds = as.numeric(log_bounds),
    control = control,
    auto_cap = auto_cap,
    allow_approximate = allow_approximate
  )
}


.dpprior_v2_effective_constraint_tolerance <- function(spec, tolerance) {
  tolerance$abs + tolerance$rel * max(abs(spec$value), 1e-8)
}


.dpprior_v2_eval_metric_safe <- function(spec, eta, M) {
  if (!is.numeric(eta) || length(eta) != 2L || any(!is.finite(eta)) ||
      any(eta < -.EXP_MAX) || any(eta > .EXP_MAX)) {
    return(.dpprior_v2_failed_metric(spec, "invalid log-parameter candidate"))
  }
  tryCatch(
    .dpprior_v2_eval_metric(spec, exp(eta[[1L]]), exp(eta[[2L]]), M = M),
    error = function(e) .dpprior_v2_failed_metric(spec, conditionMessage(e))
  )
}


.dpprior_v2_K_moments_safe <- function(J, eta, M) {
  if (!is.numeric(eta) || length(eta) != 2L || any(!is.finite(eta))) {
    return(list(ok = FALSE, error = "invalid log-parameter candidate"))
  }
  out <- tryCatch(
    exact_K_moments(J, exp(eta[[1L]]), exp(eta[[2L]]), M = M),
    error = identity
  )
  if (inherits(out, "condition")) {
    return(list(ok = FALSE, error = conditionMessage(out)))
  }
  valid <- .dpprior_v2_plain_scalar(out$mean, "numeric") &&
    .dpprior_v2_plain_scalar(out$var, "numeric") && out$mean >= 1 &&
    out$mean <= J && out$var >= 0
  if (!valid) {
    return(list(ok = FALSE, error = "K moments violated finite support invariants"))
  }
  list(ok = TRUE, mean = out$mean, var = out$var, raw = out)
}


.dpprior_v2_hard_candidate <- function(
    eta, source, fit_info, spec, M, effective_tolerance,
    optimizer_exit_code = NA_integer_, attempt_index = NA_integer_) {
  metric <- .dpprior_v2_eval_metric_safe(spec, eta, M)
  moments <- .dpprior_v2_K_moments_safe(fit_info$J, eta, M)
  if (!isTRUE(moments$ok) || !is.finite(metric$value)) {
    return(list(
      valid = FALSE,
      eta = eta,
      a = if (length(eta) >= 1L && is.finite(eta[[1L]])) exp(eta[[1L]]) else NA_real_,
      b = if (length(eta) >= 2L && is.finite(eta[[2L]])) exp(eta[[2L]]) else NA_real_,
      source = source,
      optimizer_exit_code = optimizer_exit_code,
      attempt_index = attempt_index,
      error = paste(
        c(if (!isTRUE(moments$ok)) moments$error,
          if (!is.finite(metric$value)) metric$certification$reason),
        collapse = " | "
      )
    ))
  }
  K <- .dpprior_v2_k_loss(
    list(mu_K = moments$mean, var_K = moments$var), fit_info$target_K
  )
  constraint <- .dpprior_v2_constraint_values(
    spec, metric, effective_tolerance
  )
  list(
    valid = TRUE,
    eta = as.numeric(eta),
    a = exp(eta[[1L]]),
    b = exp(eta[[2L]]),
    source = source,
    optimizer_exit_code = optimizer_exit_code,
    attempt_index = attempt_index,
    moments = moments,
    metric = metric,
    K = K,
    constraint = constraint,
    error = NA_character_
  )
}


.dpprior_v2_metric_extrema <- function(spec, log_bounds, M, M_verify,
                                       effective_tolerance) {
  # Every supported metric is non-increasing in Gamma shape a and
  # non-decreasing in Gamma rate b. Therefore these two corners are global
  # extrema on the declared rectangle, not merely sampled grid points.
  eta_minimum <- c(log_bounds[[2L]], log_bounds[[1L]])
  eta_maximum <- c(log_bounds[[1L]], log_bounds[[2L]])
  evaluate_corner <- function(eta) {
    selected <- .dpprior_v2_eval_metric_safe(spec, eta, M)
    refined <- .dpprior_v2_eval_metric_safe(spec, eta, M_verify)
    finite <- is.finite(selected$value) && is.finite(refined$value)
    uncertainty <- if (finite) {
      max(
        abs(selected$value - refined$value),
        selected$error_bound %||% 0,
        refined$error_bound %||% 0,
        64 * .Machine$double.eps
      )
    } else {
      Inf
    }
    list(
      eta = eta,
      parameters = list(a = exp(eta[[1L]]), b = exp(eta[[2L]])),
      selected = selected,
      refined = refined,
      finite = finite,
      uncertainty = uncertainty,
      lower = if (finite) max(0, refined$value - uncertainty) else NA_real_,
      upper = if (finite) min(1, refined$value + uncertainty) else NA_real_
    )
  }
  minimum <- evaluate_corner(eta_minimum)
  maximum <- evaluate_corner(eta_maximum)
  stable <- minimum$finite && maximum$finite &&
    minimum$uncertainty <= 1e-6 && maximum$uncertainty <= 1e-6
  # Agreement of two quadrature orders is a stability diagnostic, not a
  # certified enclosure. Only closed-form metrics and the outward-rounded
  # W_max upper endpoint can support a global infeasibility certificate.
  certified_extrema <- stable &&
    isTRUE(minimum$refined$certification$certified) &&
    isTRUE(maximum$refined$certification$certified)
  certified_infeasible <- FALSE
  feasible_corner <- NULL
  if (certified_extrema && identical(spec$relation, "at_most")) {
    certified_infeasible <- minimum$lower > spec$value + effective_tolerance
    if (!certified_infeasible &&
        minimum$refined$value - spec$value <= effective_tolerance) {
      feasible_corner <- minimum
    }
  } else if (certified_extrema && identical(spec$relation, "at_least")) {
    certified_infeasible <- maximum$upper < spec$value - effective_tolerance
    if (!certified_infeasible &&
        spec$value - maximum$refined$value <= effective_tolerance) {
      feasible_corner <- maximum
    }
  } else if (stable && identical(spec$relation, "at_most") &&
             minimum$refined$value - spec$value <= effective_tolerance) {
    # This is a search start only. It is not a feasibility or infeasibility
    # certificate until the independent candidate verifier passes.
    feasible_corner <- minimum
  } else if (stable && identical(spec$relation, "at_least") &&
             spec$value - maximum$refined$value <= effective_tolerance) {
    feasible_corner <- maximum
  }
  classification <- if (certified_infeasible) {
    "certified_infeasible"
  } else if (!certified_extrema) {
    "unknown"
  } else if (!is.null(feasible_corner)) {
    "feasible_corner_certified"
  } else {
    "unknown"
  }
  list(
    classification = classification,
    certified_infeasible = certified_infeasible,
    feasibility_unknown = identical(classification, "unknown"),
    feasible_corner = feasible_corner,
    minimum = minimum,
    maximum = maximum,
    certificate = list(
      type = if (certified_extrema) {
        "analytic_global_monotonicity_with_refined_order_corner_enclosures"
      } else if (stable) {
        "analytic_monotonicity_with_uncertified_two_order_stability"
      } else {
        "corner_evaluation_unstable_or_failed"
      },
      domain = list(
        log_a = log_bounds, log_b = log_bounds,
        a = exp(log_bounds), b = exp(log_bounds)
      ),
      monotonicity = list(a = "non_increasing", b = "non_decreasing"),
      M = as.integer(M),
      M_verify = as.integer(M_verify),
      effective_tolerance = effective_tolerance,
      certified = certified_extrema
    )
  )
}


.dpprior_v2_feasible_b_interval <- function(
    log_a, spec, log_bounds, M, root_tol) {
  lo <- log_bounds[[1L]]
  hi <- log_bounds[[2L]]
  value_at <- function(log_b) {
    .dpprior_v2_eval_metric_safe(spec, c(log_a, log_b), M)$value
  }
  v_lo <- value_at(lo)
  v_hi <- value_at(hi)
  if (!is.finite(v_lo) || !is.finite(v_hi)) {
    return(list(feasible = FALSE, unknown = TRUE, lower = NA_real_,
                upper = NA_real_, error = "nonfinite metric endpoint"))
  }

  if (identical(spec$relation, "at_most")) {
    if (v_lo > spec$value) {
      return(list(feasible = FALSE, unknown = FALSE, lower = NA_real_,
                  upper = NA_real_, endpoint_values = c(v_lo, v_hi)))
    }
    if (v_hi <= spec$value) {
      return(list(feasible = TRUE, unknown = FALSE, lower = lo, upper = hi,
                  endpoint_values = c(v_lo, v_hi), active = FALSE))
    }
    root <- tryCatch(
      stats::uniroot(
        function(log_b) value_at(log_b) - spec$value,
        interval = c(lo, hi), tol = root_tol
      )$root,
      error = function(e) NA_real_
    )
    return(list(
      feasible = is.finite(root), unknown = !is.finite(root), lower = lo,
      upper = root, endpoint_values = c(v_lo, v_hi), active = TRUE
    ))
  }

  if (v_hi < spec$value) {
    return(list(feasible = FALSE, unknown = FALSE, lower = NA_real_,
                upper = NA_real_, endpoint_values = c(v_lo, v_hi)))
  }
  if (v_lo >= spec$value) {
    return(list(feasible = TRUE, unknown = FALSE, lower = lo, upper = hi,
                endpoint_values = c(v_lo, v_hi), active = FALSE))
  }
  root <- tryCatch(
    stats::uniroot(
      function(log_b) value_at(log_b) - spec$value,
      interval = c(lo, hi), tol = root_tol
    )$root,
    error = function(e) NA_real_
  )
  list(
    feasible = is.finite(root), unknown = !is.finite(root), lower = root,
    upper = hi, endpoint_values = c(v_lo, v_hi), active = TRUE
  )
}


.dpprior_v2_profile_candidate <- function(
    log_a, fit_info, spec, settings, effective_tolerance,
    source = "candidate_scan", optimizer_exit_code = NA_integer_,
    attempt_index = NA_integer_) {
  interval <- .dpprior_v2_feasible_b_interval(
    log_a, spec, settings$log_bounds, settings$M,
    settings$control$root_tol
  )
  if (!isTRUE(interval$feasible)) {
    return(list(
      valid = FALSE, eta = c(log_a, NA_real_), source = source,
      optimizer_exit_code = optimizer_exit_code,
      attempt_index = attempt_index,
      error = if (isTRUE(interval$unknown)) {
        "feasible log-b interval could not be evaluated"
      } else {
        "no exact feasible log-b value at this log-a"
      }
    ))
  }

  objective_b <- function(log_b) {
    moments <- .dpprior_v2_K_moments_safe(
      fit_info$J, c(log_a, log_b), settings$M
    )
    if (!isTRUE(moments$ok)) return(.PENALTY_INF)
    .dpprior_v2_k_loss(
      list(mu_K = moments$mean, var_K = moments$var), fit_info$target_K
    )$value
  }
  width <- interval$upper - interval$lower
  b_values <- c(interval$lower, interval$upper)
  if (fit_info$b > 0) {
    input_log_b <- log(fit_info$b)
    if (input_log_b >= interval$lower && input_log_b <= interval$upper) {
      b_values <- c(b_values, input_log_b)
    }
  }
  inner_exit <- NA_integer_
  if (is.finite(width) && width > settings$control$root_tol) {
    inner <- tryCatch(
      stats::optimize(
        objective_b, interval = c(interval$lower, interval$upper),
        tol = settings$control$optim_reltol
      ),
      error = identity
    )
    if (!inherits(inner, "condition") && is.finite(inner$objective)) {
      b_values <- c(b_values, inner$minimum)
      inner_exit <- 0L
    } else {
      inner_exit <- 1L
    }
  } else if (is.finite(width) && width >= 0) {
    inner_exit <- 0L
  }
  objectives <- vapply(b_values, objective_b, numeric(1))
  if (!any(is.finite(objectives) & objectives < .PENALTY_INF)) {
    return(list(
      valid = FALSE, eta = c(log_a, NA_real_), source = source,
      optimizer_exit_code = optimizer_exit_code,
      attempt_index = attempt_index,
      error = "all feasible log-b profile evaluations failed"
    ))
  }
  best <- which.min(objectives)
  candidate <- .dpprior_v2_hard_candidate(
    c(log_a, b_values[[best]]), source, fit_info, spec, settings$M,
    effective_tolerance,
    # An inner log-b minimum at one fixed log-a is scan evidence, not a
    # successful two-parameter constrained solve. Only callers performing the
    # outer profile optimization may supply optimizer-exit evidence.
    optimizer_exit_code = optimizer_exit_code,
    attempt_index = attempt_index
  )
  candidate$profile_interval <- interval
  candidate$inner_optimizer_exit_code <- inner_exit
  candidate
}


.dpprior_v2_verify_candidate <- function(
    candidate, fit_info, spec, settings, effective_tolerance) {
  if (!is.list(candidate) || !isTRUE(candidate$valid)) {
    return(list(
      performed = FALSE, passed = FALSE, constraint_satisfied = FALSE,
      reason = "invalid_fit_order_candidate", candidate = candidate
    ))
  }
  eta <- candidate$eta
  inside_domain <- length(eta) == 2L && all(is.finite(eta)) &&
    all(eta >= settings$log_bounds[[1L]]) &&
    all(eta <= settings$log_bounds[[2L]])
  if (!inside_domain) {
    return(list(
      performed = TRUE, passed = FALSE, constraint_satisfied = FALSE,
      reason = "candidate_outside_declared_log_domain", eta = eta
    ))
  }

  moments_verify <- .dpprior_v2_K_moments_safe(
    fit_info$J, eta, settings$M_verify
  )
  metric_verify <- .dpprior_v2_eval_metric_safe(
    spec, eta, settings$M_verify
  )
  if (!isTRUE(moments_verify$ok) || !is.finite(metric_verify$value)) {
    return(list(
      performed = TRUE, passed = FALSE, constraint_satisfied = FALSE,
      reason = "independent_recomputation_failed",
      M_fit = settings$M, M_verify = settings$M_verify,
      moments_error = moments_verify$error %||% NA_character_,
      metric_error = metric_verify$certification$reason %||% NA_character_
    ))
  }

  K_verify <- .dpprior_v2_k_loss(
    list(mu_K = moments_verify$mean, var_K = moments_verify$var),
    fit_info$target_K,
    scales = candidate$K$scales
  )
  constraint_refined <- .dpprior_v2_constraint_values(
    spec, metric_verify, effective_tolerance
  )
  constraint_selected <- candidate$constraint
  K_differences <- c(
    mu_K = abs(candidate$moments$mean - moments_verify$mean),
    var_K = abs(candidate$moments$var - moments_verify$var)
  )
  K_tolerances <- settings$control$verification_abs_tol +
    settings$control$verification_rel_tol * pmax(
      abs(c(candidate$moments$mean, candidate$moments$var)),
      abs(c(moments_verify$mean, moments_verify$var)), 1e-8
    )
  names(K_tolerances) <- c("mu_K", "var_K")
  K_stable <- all(is.finite(K_differences)) &&
    all(K_differences <= K_tolerances)

  metric_difference <- abs(candidate$metric$value - metric_verify$value)
  metric_tolerance <- settings$control$verification_abs_tol +
    settings$control$verification_rel_tol * max(
      abs(candidate$metric$value), abs(metric_verify$value), 1e-8
    )
  metric_order_stable <- is.finite(metric_difference) &&
    metric_difference <= metric_tolerance
  metric_verification_kind <- if (identical(spec$metric, "wsb_mean")) {
    "selected_refined_order_agreement_not_global_enclosure"
  } else {
    "closed_form_or_certified_bound"
  }
  metric_verification_passed <- if (identical(spec$metric, "wsb_mean")) {
    metric_order_stable
  } else {
    isTRUE(metric_verify$certification$certified)
  }
  metric_selected_published <- candidate$metric
  if (identical(spec$metric, "wsb_mean") && metric_verification_passed) {
    backend_status <- metric_selected_published$status
    metric_selected_published$status <- "converged"
    metric_selected_published$usable <- TRUE
    metric_selected_published$verified <- TRUE
    metric_selected_published$certification <- list(
      certified = TRUE,
      scope = "selected_point_two_order_error_control",
      global_enclosure_certified = FALSE,
      supports_global_infeasibility = FALSE,
      reason = "selected_and_refined_order_agreement"
    )
    metric_selected_published$provenance$backend_status <- backend_status
    metric_selected_published$provenance$verification_promotion <-
      "selected_point_only_not_global_domain_certificate"
  }

  perturbations <- list()
  perturbation_values <- numeric()
  h <- settings$control$perturbation_step
  k <- 0L
  for (dimension in seq_len(2L)) {
    for (direction in c(-1, 1)) {
      perturbed <- eta
      perturbed[[dimension]] <- perturbed[[dimension]] + direction * h
      if (perturbed[[dimension]] < settings$log_bounds[[1L]] ||
          perturbed[[dimension]] > settings$log_bounds[[2L]]) next
      k <- k + 1L
      evaluated <- .dpprior_v2_eval_metric_safe(
        spec, perturbed, settings$M_verify
      )
      perturbations[[k]] <- list(
        dimension = c("log_a", "log_b")[[dimension]],
        direction = direction,
        step = h,
        value = evaluated$value,
        finite = is.finite(evaluated$value)
      )
      perturbation_values[[k]] <- evaluated$value
    }
  }
  perturbation_delta <- if (length(perturbation_values) &&
                            all(is.finite(perturbation_values))) {
    max(abs(perturbation_values - metric_verify$value))
  } else {
    Inf
  }
  perturbation_tolerance <- settings$control$perturbation_abs_tol +
    settings$control$verification_rel_tol * max(abs(metric_verify$value), 1e-8)
  perturbation_stable <- length(perturbation_values) >= 2L &&
    is.finite(perturbation_delta) &&
    perturbation_delta <= perturbation_tolerance

  probability_invariant <- metric_verify$value >= 0 &&
    metric_verify$value <= 1
  moments_invariant <- moments_verify$mean >= 1 &&
    moments_verify$mean <= fit_info$J && moments_verify$var >= 0
  numerical_passed <- K_stable && metric_order_stable &&
    metric_verification_passed && perturbation_stable &&
    probability_invariant && moments_invariant
  selected_constraint_passed <- isTRUE(constraint_selected$satisfied)
  refined_constraint_passed <- isTRUE(constraint_refined$satisfied)
  both_constraint_orders_passed <- selected_constraint_passed &&
    refined_constraint_passed && metric_verification_passed &&
    metric_order_stable
  passed <- both_constraint_orders_passed && numerical_passed
  failed_checks <- c(
    if (!selected_constraint_passed) "constraint_residual_selected_order",
    if (!refined_constraint_passed) "constraint_residual_refined_order",
    if (!K_stable) "K_order_stability",
    if (!metric_order_stable) "metric_order_stability",
    if (!metric_verification_passed) "metric_error_control",
    if (!perturbation_stable) "parameter_perturbation_stability",
    if (!probability_invariant) "probability_invariant",
    if (!moments_invariant) "K_support_invariant"
  )
  list(
    performed = TRUE,
    passed = passed,
    numerical_passed = numerical_passed,
    constraint_satisfied = both_constraint_orders_passed,
    reason = if (passed) "independent_verification_passed" else {
      paste("failed", paste(failed_checks, collapse = ", "))
    },
    method = "fresh_higher_order_recomputation_and_local_perturbation",
    settings = list(
      M_fit = settings$M,
      M_verify = settings$M_verify,
      M_verify_required = settings$M_verify_required,
      log_bounds = settings$log_bounds,
      verification_abs_tol = settings$control$verification_abs_tol,
      verification_rel_tol = settings$control$verification_rel_tol,
      perturbation_step = h,
      perturbation_abs_tol = settings$control$perturbation_abs_tol
    ),
    achieved_K_selected = list(
      mu_K = candidate$moments$mean,
      var_K = candidate$moments$var,
      units = list(mu_K = "clusters", var_K = "clusters^2"),
      M = settings$M
    ),
    achieved_K_refined = list(
      mu_K = moments_verify$mean,
      var_K = moments_verify$var,
      units = list(mu_K = "clusters", var_K = "clusters^2"),
      M = settings$M_verify
    ),
    achieved_weight_selected = metric_selected_published,
    achieved_weight_refined = metric_verify,
    K_loss = candidate$K,
    K_loss_refined = K_verify,
    constraint = constraint_selected,
    constraint_refined = constraint_refined,
    constraint_satisfaction_basis = list(
      required = c(
        "selected_order_residual", "refined_order_residual",
        "selected_refined_stability", "metric_certification"
      ),
      selected_order_passed = selected_constraint_passed,
      refined_order_passed = refined_constraint_passed,
      stability_passed = metric_order_stable,
      metric_verification_passed = metric_verification_passed,
      metric_verification_kind = metric_verification_kind,
      global_enclosure_certified = isTRUE(
        metric_verify$certification$certified
      ),
      all_passed = both_constraint_orders_passed
    ),
    order_stability = list(
      K_difference = K_differences,
      K_tolerance = K_tolerances,
      K_passed = K_stable,
      metric_difference = metric_difference,
      metric_tolerance = metric_tolerance,
      metric_passed = metric_order_stable
    ),
    perturbation = list(
      evaluations = perturbations,
      maximum_metric_delta = perturbation_delta,
      tolerance = perturbation_tolerance,
      passed = perturbation_stable
    ),
    invariants = list(
      probability = probability_invariant,
      K_support = moments_invariant,
      finite_parameters_inside_domain = inside_domain
    )
  )
}


.dpprior_v2_classify_hard_candidate <- function(
    verification, optimizer_exit_code, eta, log_bounds, boundary_tol,
    source = NULL) {
  constraint_ok <- is.list(verification) &&
    isTRUE(verification$constraint_satisfied)
  numerical_ok <- is.list(verification) &&
    (isTRUE(verification$numerical_passed) || isTRUE(verification$passed))
  verified <- is.list(verification) && isTRUE(verification$passed) &&
    constraint_ok
  optimizer_ok <- .dpprior_v2_exit_zero(optimizer_exit_code) &&
    .dpprior_v2_optimizer_evidence_allowed(source)
  if (!verified) {
    if (numerical_ok && !constraint_ok) return("approximate")
    return("failed")
  }
  if (!optimizer_ok) return("approximate")
  domain_distance <- min(
    eta - log_bounds[[1L]], log_bounds[[2L]] - eta
  )
  active_constraint <- abs(verification$constraint$raw) <=
    verification$constraint$tolerance
  if (domain_distance <= boundary_tol || active_constraint) {
    "boundary"
  } else {
    "converged"
  }
}


.dpprior_v2_hard_target <- function(spec) {
  request <- list(
    metric = spec$metric,
    relation = spec$raw_input$relation,
    bound = spec$value,
    threshold = spec$threshold,
    probability = spec$probability
  )
  normalized <- list(
    metric = spec$metric,
    relation = spec$relation,
    value = spec$value,
    threshold = spec$threshold,
    probability = spec$probability
  )
  .dpprior_new_weight_target(
    request = request,
    normalized = normalized,
    used = normalized,
    metric = spec$metric,
    relation = spec$relation,
    operator = spec$operator,
    value = spec$value,
    threshold = spec$threshold,
    probability = spec$probability,
    estimand = spec$estimand,
    units = spec$units,
    certification = if (identical(spec$metric, "wmax_tail_upper")) {
      list(
        kind = "upper_bound",
        certified = TRUE,
        passed = TRUE,
        method = "certified_size_biased_mass_upper_bound",
        source = "wmax_tail_bounds"
      )
    } else {
      list()
    },
    provenance = list(
      source = "DPprior_dual_hard request canonicalization",
      transformation = list(
        rule = "canonicalize_hard_weight_target",
        opt_in = FALSE,
        before = request,
        after = normalized,
        evidence = list(
          mode = "hard",
          value_field = "bound",
          relation_from = request$relation,
          relation_to = normalized$relation,
          probability_field = if (is.null(request$probability)) {
            "none"
          } else {
            "probability"
          }
        )
      ),
      selection = NULL
    )
  )
}


.dpprior_v2_hard_controls <- function(settings) {
  list(
    maxit = settings$control$maxit,
    scan_points = settings$control$scan_points,
    scan_keep = settings$control$scan_keep,
    profile_starts = settings$control$profile_starts,
    log_bounds = settings$log_bounds,
    root_tol = settings$control$root_tol,
    optim_reltol = settings$control$optim_reltol,
    penalty = settings$control$penalty,
    boundary_tol = settings$control$boundary_tol,
    constraint_abs_tol = settings$tolerance$abs,
    constraint_rel_tol = settings$tolerance$rel,
    verification_abs_tol = settings$control$verification_abs_tol,
    verification_rel_tol = settings$control$verification_rel_tol,
    perturbation_step = settings$control$perturbation_step,
    perturbation_abs_tol = settings$control$perturbation_abs_tol,
    auto_cap = settings$auto_cap
  )
}


.dpprior_v2_hard_tolerances <- function(spec, settings) {
  constraint <- list(
    absolute = settings$tolerance$abs,
    relative = settings$tolerance$rel,
    effective = .dpprior_v2_effective_constraint_tolerance(
      spec, settings$tolerance
    )
  )
  K <- list(
    absolute = settings$control$verification_abs_tol,
    relative = settings$control$verification_rel_tol,
    scale_floor = 1e-8
  )
  list(
    constraint = constraint,
    K = K,
    weight = K,
    perturbation = list(
      absolute = settings$control$perturbation_abs_tol,
      relative = settings$control$verification_rel_tol,
      scale_floor = 1e-8,
      step = settings$control$perturbation_step,
      min_evaluations = 2L
    ),
    boundary = settings$control$boundary_tol,
    certificate = list(
      corner_uncertainty = 1e-6,
      rounding_floor = 64 * .Machine$double.eps
    )
  )
}


.dpprior_v2_hard_parameters <- function(candidate) {
  if (!is.list(candidate) || !isTRUE(candidate$valid) ||
      !.dpprior_v2_plain_scalar(candidate$a, "numeric") ||
      !.dpprior_v2_plain_scalar(candidate$b, "numeric") ||
      candidate$a <= 0 || candidate$b <= 0) {
    return(NULL)
  }
  .dpprior_new_parameters(candidate$a, candidate$b, "log_ab")
}


.dpprior_v2_hard_achieved_weight <- function(spec, metric, source) {
  list(
    metric = spec$metric,
    value = as.numeric(metric$value),
    source = source
  )
}


.dpprior_v2_hard_snapshot <- function(
    candidate, verification, fit_info, spec, settings, tolerances,
    refined = FALSE) {
  parameters <- .dpprior_v2_hard_parameters(candidate)
  if (is.null(parameters)) return(NULL)

  if (refined) {
    if (!is.list(verification) || !isTRUE(verification$performed) ||
        is.null(verification$achieved_K_refined) ||
        is.null(verification$achieved_weight_refined)) {
      return(NULL)
    }
    old_K <- verification$achieved_K_refined
    old_weight <- verification$achieved_weight_refined
    M <- settings$M_verify
    snapshot_source <- "independent_verifier"
    achieved_source <- "independent_higher_order_recomputation"
  } else {
    old_K <- list(
      mu_K = candidate$moments$mean,
      var_K = candidate$moments$var
    )
    old_weight <- if (is.list(verification) &&
                      !is.null(verification$achieved_weight_selected)) {
      verification$achieved_weight_selected
    } else {
      candidate$metric
    }
    M <- settings$M
    snapshot_source <- "selected_order"
    achieved_source <- "selected_order_recomputation"
  }
  if (!.dpprior_v2_plain_scalar(old_K$mu_K, "numeric") ||
      !.dpprior_v2_plain_scalar(old_K$var_K, "numeric") ||
      !.dpprior_v2_plain_scalar(old_weight$value, "numeric")) {
    return(NULL)
  }

  achieved <- list(
    K = list(
      mean = as.numeric(old_K$mu_K),
      variance = as.numeric(old_K$var_K),
      estimand = "K_J",
      source = achieved_source,
      M = as.integer(M)
    ),
    weight = .dpprior_v2_hard_achieved_weight(
      spec, old_weight, achieved_source
    )
  )
  target_moments <- unclass(fit_info$target_K$canonical)[[
    "implied", exact = TRUE
  ]]
  raw_weight <- achieved$weight$value - spec$value
  directed_weight <- if (identical(spec$relation, "at_most")) {
    raw_weight
  } else {
    -raw_weight
  }
  residuals <- list(
    K = list(
      mean = achieved$K$mean - target_moments[["mean", exact = TRUE]],
      variance = achieved$K$variance -
        target_moments[["variance", exact = TRUE]]
    ),
    weight = list(raw = raw_weight, directed = directed_weight)
  )
  finite <- all(is.finite(c(
    parameters$a, parameters$b, achieved$K$mean, achieved$K$variance,
    achieved$weight$value, unlist(residuals, use.names = FALSE)
  )))
  .dpprior_new_snapshot(
    parameters = parameters,
    M = as.integer(M),
    achieved = achieved,
    residuals = residuals,
    tolerances = tolerances,
    finite = finite,
    source = snapshot_source
  )
}


.dpprior_v2_hard_candidate_checks <- function(
    id, selected_snapshot, verifier_snapshot, verification, fit_info, spec,
    tolerances) {
  source <- paste0("candidate:", id)
  selected_K <- selected_snapshot$achieved$K
  selected_weight <- selected_snapshot$achieved$weight$value
  if (is.null(verifier_snapshot)) {
    return(list(
      candidate_domain = .dpprior_new_check(
        value = c(
          K_support = selected_K$mean >= 1 && selected_K$mean <= fit_info$J,
          weight_support = selected_weight >= 0 && selected_weight <= 1
        ),
        reference = c(K_support = TRUE, weight_support = TRUE),
        tolerance = NULL,
        operator = "identical",
        source = source
      )
    ))
  }

  residual_for <- function(value) {
    if (identical(spec$relation, "at_most")) value - spec$value else
      spec$value - value
  }
  selected_values <- c(
    K.mean = selected_snapshot$achieved$K$mean,
    K.variance = selected_snapshot$achieved$K$variance,
    weight.value = selected_weight
  )
  verifier_values <- c(
    K.mean = verifier_snapshot$achieved$K$mean,
    K.variance = verifier_snapshot$achieved$K$variance,
    weight.value = verifier_snapshot$achieved$weight$value
  )
  delta <- abs(selected_values - verifier_values)
  stability_tolerance <- c(
    K.mean = tolerances$K$absolute + tolerances$K$relative * max(
      abs(selected_values[["K.mean"]]),
      abs(verifier_values[["K.mean"]]), tolerances$K$scale_floor
    ),
    K.variance = tolerances$K$absolute + tolerances$K$relative * max(
      abs(selected_values[["K.variance"]]),
      abs(verifier_values[["K.variance"]]), tolerances$K$scale_floor
    ),
    weight.value = tolerances$weight$absolute +
      tolerances$weight$relative * max(
        abs(selected_values[["weight.value"]]),
        abs(verifier_values[["weight.value"]]),
        tolerances$weight$scale_floor
      )
  )
  metric_identity <- c(
    selected_metric = identical(
      selected_snapshot$achieved$weight$metric, spec$metric
    ),
    refined_metric = identical(
      verifier_snapshot$achieved$weight$metric, spec$metric
    ),
    selected_finite = is.finite(selected_weight),
    refined_finite = is.finite(verifier_values[["weight.value"]])
  )
  invariants <- c(
    selected_finite = selected_snapshot$finite,
    refined_finite = verifier_snapshot$finite,
    parameter_identity = identical(
      selected_snapshot$parameters, verifier_snapshot$parameters
    ),
    support_valid = selected_values[["K.mean"]] >= 1 &&
      selected_values[["K.mean"]] <= fit_info$J &&
      verifier_values[["K.mean"]] >= 1 &&
      verifier_values[["K.mean"]] <= fit_info$J
  )
  maximum_metric_delta <- if (is.list(verification$perturbation) &&
      .dpprior_v2_plain_scalar(
        verification$perturbation$maximum_metric_delta, "numeric"
      )) {
    verification$perturbation$maximum_metric_delta
  } else {
    Inf
  }
  perturbation_tolerance <- tolerances$perturbation$absolute +
    tolerances$perturbation$relative * max(
      abs(selected_weight), tolerances$perturbation$scale_floor
    )
  list(
    constraint_selected = .dpprior_new_check(
      value = residual_for(selected_weight), reference = 0,
      tolerance = tolerances$constraint$effective,
      operator = "lte", source = source
    ),
    constraint_refined = .dpprior_new_check(
      value = residual_for(verifier_values[["weight.value"]]), reference = 0,
      tolerance = tolerances$constraint$effective,
      operator = "lte", source = source
    ),
    order_stability = .dpprior_new_check(
      value = delta,
      reference = setNames(rep(0, length(delta)), names(delta)),
      tolerance = stability_tolerance,
      operator = "lte", source = source
    ),
    metric_certification = .dpprior_new_check(
      value = metric_identity,
      reference = setNames(rep(TRUE, length(metric_identity)),
                           names(metric_identity)),
      tolerance = NULL,
      operator = "identical", source = source
    ),
    perturbation = .dpprior_new_check(
      value = c(maximum_metric_delta = maximum_metric_delta),
      reference = 0,
      tolerance = c(maximum_metric_delta = perturbation_tolerance),
      operator = "lte", source = source
    ),
    invariants = .dpprior_new_check(
      value = invariants,
      reference = setNames(rep(TRUE, length(invariants)), names(invariants)),
      tolerance = NULL,
      operator = "identical", source = source
    )
  )
}


.dpprior_v2_hard_candidate_parent <- function(candidate, attempts) {
  source <- candidate$source
  if (source %in% c(
    "deterministic_profile_scan", "input_fit", "feasibility_extreme"
  )) {
    return(NA_integer_)
  }
  if (.dpprior_v2_plain_scalar(candidate$attempt_index, "numeric") &&
      candidate$attempt_index == as.integer(candidate$attempt_index) &&
      candidate$attempt_index >= 1L &&
      candidate$attempt_index <= length(attempts)) {
    return(as.integer(candidate$attempt_index))
  }
  NA_integer_
}


.dpprior_v2_hard_candidate_route <- function(candidate, parent_index) {
  source <- candidate$source
  if (identical(source, "deterministic_profile_scan")) {
    return(list(
      method = "deterministic_profile_scan",
      generator = "deterministic_profile_scan"
    ))
  }
  if (identical(source, "input_fit")) {
    return(list(method = "input_fit", generator = "input_fit"))
  }
  if (identical(source, "feasibility_extreme")) {
    return(list(
      method = "analytic_feasibility_extreme",
      generator = "feasibility_extreme"
    ))
  }
  if (identical(source, "constrained_profile_optimizer")) {
    return(list(
      method = "constrained_profile_optimize", generator = "direct_attempt"
    ))
  }
  if (identical(source, "penalty_L-BFGS-B_diagnostic")) {
    return(list(
      method = "penalty_L-BFGS-B_diagnostic",
      generator = "derived_diagnostic_attempt"
    ))
  }
  list(method = "K_only_L-BFGS-B", generator = "direct_attempt")
}


.dpprior_v2_hard_candidate_descriptors <- function(internal, tolerances) {
  candidates <- internal$candidates
  attempts <- internal$attempts
  verified <- internal$verified_candidates
  verified_map <- setNames(
    verified,
    vapply(verified, function(entry) as.character(entry$index), character(1))
  )
  selected_index <- if (is.list(internal$selected)) {
    internal$selected$index
  } else {
    NA_integer_
  }
  out <- list()
  for (index in seq_along(candidates)) {
    candidate <- candidates[[index]]
    if (!isTRUE(candidate$valid) || !is.finite(candidate$K$value) ||
        is.null(.dpprior_v2_hard_parameters(candidate))) next
    id <- sprintf("candidate-%03d", index)
    verification_entry <- verified_map[[as.character(index)]]
    verification <- if (is.null(verification_entry)) NULL else
      verification_entry$verification
    selected_snapshot <- .dpprior_v2_hard_snapshot(
      candidate, verification, internal$fit_info, internal$spec,
      internal$settings, tolerances, refined = FALSE
    )
    verifier_snapshot <- .dpprior_v2_hard_snapshot(
      candidate, verification, internal$fit_info, internal$spec,
      internal$settings, tolerances, refined = TRUE
    )
    if (is.null(selected_snapshot)) next
    parent_index <- .dpprior_v2_hard_candidate_parent(candidate, attempts)
    route <- .dpprior_v2_hard_candidate_route(candidate, parent_index)
    checks <- .dpprior_v2_hard_candidate_checks(
      id, selected_snapshot, verifier_snapshot, verification,
      internal$fit_info, internal$spec, tolerances
    )
    check_pass <- vapply(
      checks, function(check) isTRUE(check$passed), logical(1)
    )
    diagnostic_eligible <- !is.null(verifier_snapshot) &&
      all(check_pass[c(
        "order_stability", "metric_certification", "perturbation",
        "invariants"
      )]) &&
      !all(check_pass[c("constraint_selected", "constraint_refined")])
    out[[length(out) + 1L]] <- list(
      id = id,
      raw_index = as.integer(index),
      parent_index = parent_index,
      method = route$method,
      generator = route$generator,
      parameters = selected_snapshot$parameters,
      selected_snapshot = selected_snapshot,
      verifier_snapshot = verifier_snapshot,
      checks = checks,
      fresh_objective = as.numeric(candidate$K$value),
      diagnostic_eligible = diagnostic_eligible,
      selected = identical(as.integer(index), as.integer(selected_index))
    )
  }
  out
}


.dpprior_v2_hard_attempt_value <- function(value, kind = c("record", "count")) {
  kind <- match.arg(kind)
  if (kind == "count") {
    if (.dpprior_v2_plain_scalar(value, "numeric") && value >= 0 &&
        value == as.integer(value)) return(as.integer(value))
    return(NULL)
  }
  if (is.null(value)) return(NULL)
  if (is.numeric(value) && !is.object(value) && is.null(dim(value)) &&
      !anyNA(value) && all(is.finite(value))) return(value)
  if (is.list(value) && !is.object(value) && !anyDuplicated(names(value))) {
    valid <- vapply(value, function(component) {
      is.numeric(component) && !is.object(component) &&
        is.null(dim(component)) && !anyNA(component) &&
        all(is.finite(component))
    }, logical(1))
    if (all(valid)) return(value)
  }
  NULL
}


.dpprior_v2_hard_attempts <- function(
    old_attempts, descriptors, status, certificate = NULL) {
  if (identical(status, "infeasible")) {
    active <- if (identical(certificate$relation, "at_most")) {
      certificate$minimum
    } else {
      certificate$maximum
    }
    return(list(.dpprior_new_attempt(
      id = "attempt-feasibility-probe",
      stage = "feasibility",
      method = "analytic_monotonicity_feasibility_probe",
      start = NULL,
      bounds = list(
        log_a = certificate$domain$log_a,
        log_b = certificate$domain$log_b
      ),
      control = list(
        M = certificate$M_selected,
        M_verify = certificate$M_verification
      ),
      exit_code = 0L,
      message = "certified_infeasible",
      iterations = 0L,
      evaluations = list(function_count = 2L),
      candidate_parameters = NULL,
      candidate_objective = active$refined$value,
      elapsed_seconds = 0,
      warnings = character(),
      error = NULL,
      selected = FALSE,
      reason_code = "globally_infeasible_by_certificate",
      unavailable = c(
        start = "analytic global corner proof has no optimizer start",
        candidate_parameters = "certificate route has no public candidate"
      )
    )))
  }

  direct_descriptor <- function(index) {
    matches <- which(vapply(descriptors, function(descriptor) {
      identical(descriptor$parent_index, as.integer(index)) &&
        descriptor$generator %in% c(
          "direct_attempt", "derived_diagnostic_attempt"
        )
    }, logical(1)))
    if (length(matches)) descriptors[[matches[[1L]]]] else NULL
  }
  out <- vector("list", length(old_attempts))
  for (index in seq_along(old_attempts)) {
    old <- old_attempts[[index]]
    method <- old$method
    stage <- switch(
      method,
      analytic_monotonicity_feasibility_probe = "feasibility",
      `K_only_L-BFGS-B` = "optimizer",
      deterministic_feasible_profile_scan = "profile_scan",
      constrained_profile_optimize = "optimizer",
      `penalty_L-BFGS-B_diagnostic` = "diagnostic",
      "optimizer"
    )
    aggregate_scan <- identical(method, "deterministic_feasible_profile_scan")
    probe <- identical(method, "analytic_monotonicity_feasibility_probe")
    descriptor <- direct_descriptor(index)
    candidate_parameters <- if (!aggregate_scan && !probe &&
        !is.null(descriptor)) descriptor$parameters else NULL
    candidate_objective <- if (!aggregate_scan && !probe &&
        !is.null(descriptor) &&
        .dpprior_v2_plain_scalar(old$candidate_objective, "numeric")) {
      as.numeric(old$candidate_objective)
    } else {
      NULL
    }
    start <- .dpprior_v2_hard_attempt_value(old$start)
    bounds <- .dpprior_v2_hard_attempt_value(old$bounds)
    control <- if (is.list(old$control) && !is.object(old$control) &&
        !anyDuplicated(names(old$control))) old$control else NULL
    exit_code <- .dpprior_v2_hard_attempt_value(old$exit_code, "count")
    iterations <- .dpprior_v2_hard_attempt_value(old$iterations, "count")
    evaluation_count <- .dpprior_v2_hard_attempt_value(
      old$evaluations, "count"
    )
    evaluations <- if (is.null(evaluation_count)) NULL else
      list(function_count = evaluation_count)
    elapsed <- if (.dpprior_v2_plain_scalar(old$elapsed, "numeric") &&
        old$elapsed >= 0) as.numeric(old$elapsed) else NULL
    warning <- old$warning
    warnings <- if (is.character(warning) && length(warning) == 1L &&
        !is.na(warning) && nzchar(warning)) strsplit(warning, " \\| ")[[1L]] else
      character()
    error_message <- old$error
    error <- if (is.character(error_message) && length(error_message) == 1L &&
        !is.na(error_message) && nzchar(error_message)) {
      list(
        class = "simpleError", code = "optimizer_error",
        message = error_message
      )
    } else {
      NULL
    }
    selected <- !is.null(descriptor) && isTRUE(descriptor$selected) &&
      identical(descriptor$generator, "direct_attempt")
    reason_code <- if (selected) {
      "selected"
    } else if (aggregate_scan) {
      "diagnostic_only"
    } else if (!is.null(descriptor) && all(vapply(
      descriptor$checks, function(check) isTRUE(check$passed), logical(1)
    ))) {
      "eligible_not_selected"
    } else if (!is.null(error)) {
      "optimizer_error"
    } else if (!is.null(exit_code) && exit_code != 0L) {
      "optimizer_exit_nonzero"
    } else if (is.null(candidate_parameters)) {
      "nonfinite_candidate"
    } else if (isTRUE(descriptor$diagnostic_eligible)) {
      "constraint_verification_failed"
    } else {
      "candidate_eligibility_failed"
    }
    nullable <- list(
      start = start, bounds = bounds, control = control,
      exit_code = exit_code, iterations = iterations,
      evaluations = evaluations, candidate_parameters = candidate_parameters,
      candidate_objective = candidate_objective,
      elapsed_seconds = elapsed
    )
    missing <- names(nullable)[vapply(nullable, is.null, logical(1))]
    unavailable <- setNames(
      paste("source attempt did not retain", gsub("_", " ", missing)),
      missing
    )
    out[[index]] <- .dpprior_new_attempt(
      id = sprintf("attempt-%03d", index),
      stage = stage,
      method = method,
      start = start,
      bounds = bounds,
      control = control,
      exit_code = exit_code,
      message = if (is.character(old$message) && length(old$message) == 1L &&
          !is.na(old$message)) old$message else "",
      iterations = iterations,
      evaluations = evaluations,
      candidate_parameters = candidate_parameters,
      candidate_objective = candidate_objective,
      elapsed_seconds = elapsed,
      warnings = warnings,
      error = error,
      selected = selected,
      reason_code = reason_code,
      unavailable = unavailable
    )
  }
  out
}


.dpprior_v2_hard_candidate_evaluations <- function(
    descriptors, attempts, tolerances) {
  lapply(descriptors, function(descriptor) {
    attempt_id <- if (is.na(descriptor$parent_index)) NULL else
      sprintf("attempt-%03d", descriptor$parent_index)
    attempt <- if (is.null(attempt_id)) NULL else
      attempts[[descriptor$parent_index]]
    recorded_objective <- if (is.null(attempt)) NULL else
      attempt$candidate_objective
    recorded_kind <- if (is.null(recorded_objective)) {
      NULL
    } else if (identical(descriptor$generator,
                         "derived_diagnostic_attempt")) {
      "penalized_diagnostic"
    } else {
      "K_loss"
    }
    execution_success <- if (is.null(attempt)) TRUE else
      identical(attempt$exit_code, 0L) && is.null(attempt$error)
    optimizer_supported <- !is.null(attempt) &&
      identical(descriptor$generator, "direct_attempt") &&
      execution_success && descriptor$method %in% c(
        "dual_anchor_hard_inequality", "L-BFGS-B",
        "K_only_L-BFGS-B", "constrained_profile_optimize"
      ) && attempt$stage %in% c("primary", "optimizer", "profile")
    .dpprior_new_candidate_evaluation(
      id = descriptor$id,
      attempt_id = attempt_id,
      method = descriptor$method,
      generator = descriptor$generator,
      parameters = descriptor$parameters,
      objective_kind = "K_loss",
      recorded_objective_kind = recorded_kind,
      selection_objective_kind = if (
        is.null(descriptor$verifier_snapshot)
      ) NULL else "K_loss",
      recorded_objective = recorded_objective,
      recorded_objective_reason = if (is.null(recorded_objective)) {
        "source route retained no comparable per-candidate objective"
      } else {
        NULL
      },
      fresh_objective = descriptor$fresh_objective,
      selection_objective = if (is.null(descriptor$verifier_snapshot)) {
        NULL
      } else {
        descriptor$fresh_objective
      },
      objective_tolerance = tolerances$K$absolute,
      selected_snapshot = descriptor$selected_snapshot,
      verifier_snapshot = descriptor$verifier_snapshot,
      checks = descriptor$checks,
      execution_success = execution_success,
      optimizer_supported = optimizer_supported,
      diagnostic_eligible = descriptor$diagnostic_eligible &&
        execution_success,
      selected = descriptor$selected,
      source = "R/21_dual_anchor_hard.R:candidate_ledger"
    )
  })
}


.dpprior_v2_hard_stability <- function(selected, verifier, tolerances) {
  selected_values <- c(
    K.mean = selected$achieved$K$mean,
    K.variance = selected$achieved$K$variance,
    weight.value = selected$achieved$weight$value
  )
  verifier_values <- c(
    K.mean = verifier$achieved$K$mean,
    K.variance = verifier$achieved$K$variance,
    weight.value = verifier$achieved$weight$value
  )
  delta <- abs(selected_values - verifier_values)
  tolerance <- c(
    K.mean = tolerances$K$absolute + tolerances$K$relative * max(
      abs(selected_values[["K.mean"]]),
      abs(verifier_values[["K.mean"]]), tolerances$K$scale_floor
    ),
    K.variance = tolerances$K$absolute + tolerances$K$relative * max(
      abs(selected_values[["K.variance"]]),
      abs(verifier_values[["K.variance"]]), tolerances$K$scale_floor
    ),
    weight.value = tolerances$weight$absolute +
      tolerances$weight$relative * max(
        abs(selected_values[["weight.value"]]),
        abs(verifier_values[["weight.value"]]),
        tolerances$weight$scale_floor
      )
  )
  .dpprior_new_stability(
    delta = delta,
    tolerance = tolerance,
    formula = setNames(
      rep("absolute_plus_relative_max", length(delta)), names(delta)
    ),
    scale_floor = c(
      K.mean = tolerances$K$scale_floor,
      K.variance = tolerances$K$scale_floor,
      weight.value = tolerances$weight$scale_floor
    ),
    source = "independent_verifier"
  )
}


.dpprior_v2_hard_unavailable_optimality <- function(reason) {
  list(
    performed = FALSE,
    passed = FALSE,
    selection_rule = NULL,
    K_scales = NULL,
    selected_K_loss = NULL,
    minimum_K_loss = NULL,
    tie_tolerance = NULL,
    perturbation_passed = NULL,
    source = "independent_candidate_ledger",
    unavailable_reason = reason
  )
}


.dpprior_v2_hard_certificate_evaluation <- function(metric, M) {
  list(
    metric = metric$metric,
    value = as.numeric(metric$value),
    method = metric$method,
    error_bound = as.numeric(metric$error_bound),
    certified = isTRUE(metric$certification$certified),
    M = as.integer(M),
    source = "independent_metric_evaluator"
  )
}


.dpprior_v2_hard_certificate_corner <- function(corner, M, M_verify) {
  eta <- c(log_a = corner$eta[[1L]], log_b = corner$eta[[2L]])
  list(
    eta = eta,
    parameters = list(a = corner$parameters$a, b = corner$parameters$b),
    selected = .dpprior_v2_hard_certificate_evaluation(
      corner$selected, M
    ),
    refined = .dpprior_v2_hard_certificate_evaluation(
      corner$refined, M_verify
    ),
    finite = TRUE,
    uncertainty = as.numeric(corner$uncertainty),
    lower = as.numeric(corner$lower),
    upper = as.numeric(corner$upper)
  )
}


.dpprior_v2_hard_certificate <- function(
    feasibility, fit_info, spec, settings, effective_tolerance) {
  minimum <- .dpprior_v2_hard_certificate_corner(
    feasibility$minimum, settings$M, settings$M_verify
  )
  maximum <- .dpprior_v2_hard_certificate_corner(
    feasibility$maximum, settings$M, settings$M_verify
  )
  list(
    method = "analytic_global_monotonicity_with_refined_order_corner_enclosures",
    version = "1",
    J = fit_info$J,
    metric = spec$metric,
    relation = spec$relation,
    target_value = spec$value,
    threshold = spec$threshold,
    probability = spec$probability,
    support = c(0, 1),
    domain = list(
      log_a = settings$log_bounds,
      log_b = settings$log_bounds,
      a = exp(settings$log_bounds),
      b = exp(settings$log_bounds)
    ),
    monotonicity = list(a = "non_increasing", b = "non_decreasing"),
    M_selected = settings$M,
    M_verification = settings$M_verify,
    effective_tolerance = effective_tolerance,
    minimum = minimum,
    maximum = maximum,
    lower_bound = minimum$lower,
    upper_bound = maximum$upper,
    tolerance = effective_tolerance,
    certified = TRUE,
    source = "phase8_metric_extrema_corner_enclosures"
  )
}


.dpprior_v2_hard_provenance <- function(
    fit_info, status, settings, diagnostic_selection = FALSE) {
  target_raw <- unclass(fit_info$target_K$canonical)
  target_provenance <- target_raw[["provenance", exact = TRUE]]
  approximation <- identical(status, "approximate")
  .dpprior_new_provenance(
    requested_method = "dual_anchor_hard_inequality",
    selected_method = "dual_anchor_hard_inequality",
    is_fallback = FALSE,
    approximation = list(
      active = approximation,
      opt_in = FALSE,
      kind = if (!approximation) NULL else if (diagnostic_selection) {
        "signed_unsatisfied_hard_diagnostic"
      } else {
        "optimizer_unsupported_hard_candidate"
      },
      warning_code = if (approximation) "dual_hard_approximate" else NULL
    ),
    projection = target_provenance[["projection", exact = TRUE]],
    parameterization = "log_ab",
    backend = list(
      package = "DPprior",
      package_version = tryCatch(
        as.character(utils::packageVersion("DPprior")),
        error = function(condition) "development"
      ),
      implementation = "R/21_dual_anchor_hard.R:DPprior_dual_hard",
      source_commit = NULL
    ),
    input_fit = fit_info$source$canonical_reference,
    migration = list(
      source_schema = "native",
      adapter = "none",
      lossless = TRUE,
      missing_evidence = character(),
      warnings = character()
    ),
    legacy = list(
      active = FALSE, contract = NULL, deprecation_stage = NULL
    )
  )
}


.dpprior_v2_hard_canonical_result <- function(old_result) {
  old <- unclass(old_result)
  internal <- old[[".internal", exact = TRUE]]
  fit_info <- internal$fit_info
  spec <- internal$spec
  settings <- internal$settings
  status <- old[["status", exact = TRUE]]
  message <- old[["message", exact = TRUE]]
  tolerances <- .dpprior_v2_hard_tolerances(spec, settings)
  weight_target <- .dpprior_v2_hard_target(spec)
  target <- list(K = fit_info$target_K$canonical, weight = weight_target)
  scales <- list(K = list(
    mean = max(abs(fit_info$target_K$mu_K), 1),
    variance = max(abs(fit_info$target_K$var_K), 1)
  ))
  controls <- .dpprior_v2_hard_controls(settings)
  computation_settings <- list(
    method = "dual_anchor_hard_inequality",
    controls = controls,
    parameterization = "log_ab"
  )
  scaling <- .dpprior_new_scaling(
    requested = scales,
    used = scales,
    formula = "fixed_from_input_target_max_abs_one",
    values = scales,
    fixed_from_input = TRUE
  )

  if (identical(status, "infeasible")) {
    effective <- tolerances$constraint$effective
    certificate <- .dpprior_v2_hard_certificate(
      old$feasibility, fit_info, spec, settings, effective
    )
    attempts <- .dpprior_v2_hard_attempts(
      internal$attempts, list(), status, certificate
    )
    computation <- .dpprior_new_computation(
      request = computation_settings,
      used = computation_settings,
      orders = .dpprior_new_orders(
        M_requested = NULL, M_selected = NULL,
        M_verification_required = NULL, M_verification_used = NULL,
        requested_reason = "not_applicable",
        selected_reason = "no_candidate",
        verification_required_reason = "certificate_only",
        verification_used_reason = "certificate_only"
      ),
      scaling = scaling,
      attempts = attempts,
      candidate_evaluations = list(),
      selected_candidate_id = NULL,
      selected_attempt_id = NULL,
      fallback = .dpprior_new_fallback(),
      termination = .dpprior_new_termination(
        code = "certified_infeasible",
        message = message,
        source = "analytic_certificate",
        iterations = NULL,
        boundary_reason = NULL
      ),
      resources = list(
        requested_M = settings$M,
        requested_M_verify = settings$M_verify,
        log_bounds = settings$log_bounds
      )
    )
    verification <- .dpprior_new_verification(
      method = "analytic_global_monotonicity_corner_certificate",
      performed = TRUE,
      passed = TRUE,
      reason = "global domain infeasibility certified",
      settings = list(),
      selected_snapshot = NULL,
      verifier_snapshot = NULL,
      stability = NULL,
      components = list(
        infeasibility_certificate = .dpprior_new_check(
          value = if (identical(spec$relation, "at_most")) {
            certificate$lower_bound
          } else {
            certificate$upper_bound
          },
          reference = certificate$target_value,
          tolerance = certificate$tolerance,
          operator = if (identical(spec$relation, "at_most")) "gt" else "lt",
          source = "analytic_certificate"
        )
      ),
      invariants = list(
        domain_monotonicity = .dpprior_new_check(
          value = TRUE, reference = TRUE, tolerance = NULL,
          operator = "identical", source = "analytic_certificate"
        )
      )
    )
    constraint <- .dpprior_new_constraint(
      relation = spec$relation,
      operator = spec$operator,
      residual = NULL,
      slack = NULL,
      tolerance = NULL,
      satisfied = NULL,
      active = NULL,
      feasibility = list(
        classification = "certified_infeasible",
        certified_infeasible = TRUE,
        feasibility_unknown = FALSE,
        certificate = certificate,
        candidate_count = 0L,
        verified_candidate_count = 0L,
        feasible_candidate_count = 0L
      ),
      optimality = .dpprior_v2_hard_unavailable_optimality(
        "certificate-only route has no finite candidate"
      )
    )
    return(.dpprior_new_fit(
      mode = "dual_hard", method = "dual_anchor_hard_inequality",
      J = fit_info$J, status = "infeasible", usable = FALSE,
      verified = TRUE, message = message, parameters = NULL,
      target = target, achieved = list(), residuals = list(),
      tolerances = tolerances, computation = computation,
      verification = verification,
      provenance = .dpprior_v2_hard_provenance(
        fit_info, "infeasible", settings
      ),
      extension = list(constraint = constraint)
    ))
  }

  internal$attempts <- old[["attempts", exact = TRUE]]
  descriptors <- .dpprior_v2_hard_candidate_descriptors(
    internal, tolerances
  )
  has_selected <- any(vapply(
    descriptors, function(descriptor) isTRUE(descriptor$selected), logical(1)
  ))
  if (!has_selected) {
    descriptors <- lapply(descriptors, function(descriptor) {
      descriptor$selected <- FALSE
      descriptor
    })
  }
  attempts <- .dpprior_v2_hard_attempts(
    internal$attempts, descriptors, status
  )
  evaluations <- .dpprior_v2_hard_candidate_evaluations(
    descriptors, attempts, tolerances
  )
  selected_index <- which(vapply(
    evaluations, function(evaluation) isTRUE(evaluation$selected), logical(1)
  ))
  selected_evaluation <- if (length(selected_index) == 1L) {
    evaluations[[selected_index]]
  } else {
    NULL
  }
  finite_result <- !is.null(selected_evaluation)
  verified_count <- sum(vapply(
    evaluations,
    function(evaluation) !is.null(evaluation$verifier_snapshot), logical(1)
  ))
  feasible_count <- sum(vapply(
    evaluations,
    function(evaluation) isTRUE(evaluation$selection_eligible), logical(1)
  ))

  if (!finite_result) {
    has_finite_candidates <- length(evaluations) > 0L
    orders <- if (has_finite_candidates) {
      .dpprior_new_orders(
        M_requested = settings$M,
        M_selected = settings$M,
        M_verification_required = settings$M_verify_required,
        M_verification_used = settings$M_verify,
        requested_reason = "public_argument",
        selected_reason = "candidate_evaluation_order",
        verification_required_reason = "independent_order_policy",
        verification_used_reason = "public_argument_or_required_default"
      )
    } else {
      .dpprior_new_orders(
        M_requested = NULL, M_selected = NULL,
        M_verification_required = NULL, M_verification_used = NULL,
        requested_reason = "not_applicable",
        selected_reason = "no_candidate",
        verification_required_reason = "not_applicable",
        verification_used_reason = "not_applicable"
      )
    }
    computation <- .dpprior_new_computation(
      request = computation_settings,
      used = computation_settings,
      orders = orders,
      scaling = scaling,
      attempts = attempts,
      candidate_evaluations = evaluations,
      selected_candidate_id = NULL,
      selected_attempt_id = NULL,
      fallback = .dpprior_new_fallback(),
      termination = .dpprior_new_termination(
        code = "no_candidate", message = message,
        source = "no_candidate", iterations = NULL,
        boundary_reason = NULL
      ),
      resources = list(
        generated_candidate_count = as.integer(length(evaluations)),
        verifier_invocation_count = as.integer(verified_count),
        requested_M = settings$M,
        requested_M_verify = settings$M_verify,
        log_bounds = settings$log_bounds
      )
    )
    verification <- .dpprior_new_verification(
      method = "no_candidate", performed = FALSE, passed = FALSE,
      reason = "no candidate passed the public selection contract",
      settings = list(), selected_snapshot = NULL,
      verifier_snapshot = NULL, stability = NULL,
      components = list(),
      invariants = list(
        no_public_candidate = .dpprior_new_check(
          value = TRUE, reference = TRUE, tolerance = NULL,
          operator = "identical", source = "independent_verifier"
        )
      )
    )
    constraint <- .dpprior_new_constraint(
      relation = spec$relation, operator = spec$operator,
      residual = NULL, slack = NULL, tolerance = NULL,
      satisfied = NULL, active = NULL,
      feasibility = list(
        classification = "unknown", certified_infeasible = FALSE,
        feasibility_unknown = TRUE, certificate = NULL,
        candidate_count = as.integer(length(evaluations)),
        verified_candidate_count = as.integer(verified_count),
        feasible_candidate_count = as.integer(feasible_count)
      ),
      optimality = .dpprior_v2_hard_unavailable_optimality(
        "no finite public candidate passed the selection contract"
      )
    )
    return(.dpprior_new_fit(
      mode = "dual_hard", method = "dual_anchor_hard_inequality",
      J = fit_info$J, status = "failed", usable = FALSE,
      verified = FALSE, message = message, parameters = NULL,
      target = target, achieved = list(), residuals = list(),
      tolerances = tolerances, computation = computation,
      verification = verification,
      provenance = .dpprior_v2_hard_provenance(
        fit_info, "failed", settings
      ),
      extension = list(constraint = constraint)
    ))
  }

  selected_descriptor <- descriptors[[selected_index]]
  selected_snapshot <- selected_evaluation$selected_snapshot
  verifier_snapshot <- selected_evaluation$verifier_snapshot
  selected_attempt_id <- selected_evaluation$attempt_id
  owns_attempt <- !is.null(selected_attempt_id) &&
    identical(selected_evaluation$generator, "direct_attempt")
  if (!owns_attempt) selected_attempt_id <- NULL
  selected_attempt <- if (is.null(selected_attempt_id)) NULL else
    attempts[[match(
      selected_attempt_id,
      vapply(attempts, function(attempt) attempt$id, character(1))
    )]]
  pool <- which(vapply(
    evaluations,
    function(evaluation) isTRUE(evaluation$selection_eligible), logical(1)
  ))
  if (!length(pool)) {
    pool <- which(vapply(
      evaluations,
      function(evaluation) isTRUE(evaluation$diagnostic_eligible), logical(1)
    ))
  }
  pool_objectives <- vapply(
    evaluations[pool],
    function(evaluation) evaluation$selection_objective,
    numeric(1)
  )
  minimum_loss <- min(pool_objectives)
  tie_tolerance <- max(
    tolerances$K$absolute,
    64 * .Machine$double.eps * max(1, abs(minimum_loss))
  )
  candidate_checks <- selected_evaluation$checks
  optimality <- list(
    performed = TRUE,
    passed = selected_evaluation$selection_objective <=
      minimum_loss + tie_tolerance &&
      isTRUE(candidate_checks$perturbation$passed),
    selection_rule = "minimum_K_loss",
    K_scales = scales$K,
    selected_K_loss = selected_evaluation$selection_objective,
    minimum_K_loss = minimum_loss,
    tie_tolerance = tie_tolerance,
    perturbation_passed = isTRUE(candidate_checks$perturbation$passed),
    source = "independent_candidate_ledger",
    unavailable_reason = NULL
  )
  selected_weight <- selected_snapshot$achieved$weight$value
  refined_weight <- verifier_snapshot$achieved$weight$value
  residual_for <- function(value) {
    if (identical(spec$relation, "at_most")) value - spec$value else
      spec$value - value
  }
  selected_residual <- residual_for(selected_weight)
  refined_residual <- residual_for(refined_weight)
  satisfied <- selected_residual <= tolerances$constraint$effective &&
    refined_residual <= tolerances$constraint$effective
  diagnostic_selection <- isTRUE(selected_evaluation$diagnostic_eligible)
  top_verified <- status %in% c("converged", "boundary")
  top_usable <- top_verified
  stability <- .dpprior_v2_hard_stability(
    selected_snapshot, verifier_snapshot, tolerances
  )
  components <- list(
    constraint_selected = .dpprior_new_check(
      value = selected_residual, reference = 0,
      tolerance = tolerances$constraint$effective,
      operator = "lte", source = "independent_verifier"
    ),
    constraint_refined = .dpprior_new_check(
      value = refined_residual, reference = 0,
      tolerance = tolerances$constraint$effective,
      operator = "lte", source = "independent_verifier"
    ),
    order_stability = .dpprior_new_check(
      value = stability$delta[["weight.value"]], reference = 0,
      tolerance = stability$tolerance[["weight.value"]],
      operator = "lte", source = "independent_verifier"
    ),
    metric_certification = .dpprior_new_check(
      value = candidate_checks$metric_certification$value,
      reference = candidate_checks$metric_certification$reference,
      tolerance = NULL, operator = "identical",
      source = "independent_verifier"
    ),
    candidate_selection = .dpprior_new_check(
      value = selected_evaluation$selection_objective - minimum_loss,
      reference = 0, tolerance = tie_tolerance,
      operator = "lte", source = "independent_verifier"
    ),
    perturbation = .dpprior_new_check(
      value = isTRUE(candidate_checks$perturbation$passed),
      reference = TRUE, tolerance = NULL,
      operator = "identical", source = "independent_verifier"
    )
  )
  invariants <- list(
    probability = .dpprior_new_check(
      value = selected_weight >= 0 && selected_weight <= 1,
      reference = TRUE, tolerance = NULL,
      operator = "identical", source = "independent_verifier"
    ),
    K_support = .dpprior_new_check(
      value = selected_snapshot$achieved$K$mean >= 1 &&
        selected_snapshot$achieved$K$mean <= fit_info$J &&
        verifier_snapshot$achieved$K$mean >= 1 &&
        verifier_snapshot$achieved$K$mean <= fit_info$J,
      reference = TRUE, tolerance = NULL,
      operator = "identical", source = "independent_verifier"
    ),
    finite_parameters_inside_domain = .dpprior_new_check(
      value = all(
        log(unlist(selected_snapshot$parameters[c("a", "b")])) >=
          settings$log_bounds[[1L]] &
          log(unlist(selected_snapshot$parameters[c("a", "b")])) <=
            settings$log_bounds[[2L]]
      ),
      reference = TRUE, tolerance = NULL,
      operator = "identical", source = "independent_verifier"
    )
  )
  verification <- .dpprior_new_verification(
    method = "fresh_higher_order_recomputation_and_local_perturbation",
    performed = TRUE,
    passed = top_verified,
    reason = if (top_verified) {
      "selected candidate passed all decision checks"
    } else if (diagnostic_selection) {
      "selected candidate is a signed unsatisfied diagnostic"
    } else {
      "selected candidate lacks ordinary optimizer support"
    },
    settings = list(
      M_fit = settings$M,
      M_verify = settings$M_verify,
      M_verify_required = settings$M_verify_required,
      log_bounds = settings$log_bounds,
      verification_abs_tol = settings$control$verification_abs_tol,
      verification_rel_tol = settings$control$verification_rel_tol,
      perturbation_step = settings$control$perturbation_step,
      perturbation_abs_tol = settings$control$perturbation_abs_tol
    ),
    selected_snapshot = selected_snapshot,
    verifier_snapshot = verifier_snapshot,
    stability = stability,
    components = components,
    invariants = invariants
  )
  feasibility <- list(
    classification = if (diagnostic_selection) "unknown" else
      "feasible_candidate",
    certified_infeasible = FALSE,
    feasibility_unknown = diagnostic_selection,
    certificate = NULL,
    candidate_count = as.integer(length(evaluations)),
    verified_candidate_count = as.integer(verified_count),
    feasible_candidate_count = as.integer(feasible_count)
  )
  constraint <- .dpprior_new_constraint(
    relation = spec$relation, operator = spec$operator,
    residual = selected_residual, slack = -selected_residual,
    tolerance = tolerances$constraint,
    satisfied = satisfied,
    active = abs(selected_residual) <= tolerances$constraint$effective,
    feasibility = feasibility,
    optimality = optimality
  )
  termination_iterations <- if (is.null(selected_attempt)) NULL else
    selected_attempt$iterations
  computation <- .dpprior_new_computation(
    request = computation_settings,
    used = computation_settings,
    orders = .dpprior_new_orders(
      M_requested = settings$M,
      M_selected = settings$M,
      M_verification_required = settings$M_verify_required,
      M_verification_used = settings$M_verify,
      requested_reason = "public_argument",
      selected_reason = "selected_order_recomputation",
      verification_required_reason = "independent_order_policy",
      verification_used_reason = "public_argument_or_required_default"
    ),
    scaling = scaling,
    attempts = attempts,
    candidate_evaluations = evaluations,
    selected_candidate_id = selected_evaluation$id,
    selected_attempt_id = selected_attempt_id,
    fallback = .dpprior_new_fallback(),
    termination = .dpprior_new_termination(
      code = status,
      message = message,
      source = if (identical(status, "approximate")) {
        "candidate_evaluation"
      } else {
        "optimizer"
      },
      iterations = termination_iterations,
      boundary_reason = if (identical(status, "boundary")) {
        if (min(
          log(unlist(selected_snapshot$parameters[c("a", "b")])) -
            settings$log_bounds[[1L]],
          settings$log_bounds[[2L]] -
            log(unlist(selected_snapshot$parameters[c("a", "b")]))
        ) <= settings$control$boundary_tol) {
          "declared_log_parameter_boundary"
        } else {
          "active_hard_constraint"
        }
      } else {
        NULL
      }
    ),
    resources = list(
      generated_candidate_count = as.integer(length(evaluations)),
      verifier_invocation_count = as.integer(verified_count),
      feasible_candidate_count = as.integer(feasible_count),
      requested_M = settings$M,
      requested_M_verify = settings$M_verify,
      log_bounds = settings$log_bounds,
      probability_alias_used = "prob" %in% names(spec$raw_input)
    )
  )
  .dpprior_new_fit(
    mode = "dual_hard", method = "dual_anchor_hard_inequality",
    J = fit_info$J, status = status, usable = top_usable,
    verified = top_verified, message = message,
    parameters = selected_snapshot$parameters,
    target = target,
    achieved = selected_snapshot$achieved,
    residuals = selected_snapshot$residuals,
    tolerances = tolerances,
    computation = computation,
    verification = verification,
    provenance = .dpprior_v2_hard_provenance(
      fit_info, status, settings, diagnostic_selection
    ),
    extension = list(constraint = constraint)
  )
}


.dpprior_v2_condition_from_result <- function(result) {
  result <- .dpprior_require_schema(
    result, kind = "fit", allow_legacy = FALSE
  )
  raw <- unclass(result)
  status <- raw[["status", exact = TRUE]]
  message <- raw[["message", exact = TRUE]]
  subclass <- switch(
    status,
    infeasible = "dpprior_dual_infeasible",
    approximate = "dpprior_dual_approximation_error",
    failed = "dpprior_dual_hard_error",
    "dpprior_dual_hard_error"
  )
  .dpprior_new_condition(
    message = message,
    classes = c(subclass, "dpprior_dual_numerical_error",
                "dpprior_numerical_error", "dpprior_error", "error"),
    code = paste0("dual_hard_", status),
    status = status,
    result = result
  )
}


.dpprior_v2_empty_hard_result <- function(
    fit_info, spec, settings, feasibility, attempts, status, message,
    candidates = list(), verified_candidates = list()) {
  verified <- identical(status, "infeasible") &&
    isTRUE(feasibility$certified_infeasible) &&
    isTRUE(feasibility$certificate$certified)
  structure(
    list(
      mode = "hard_inequality",
      method = "dual_anchor_hard_inequality",
      a = NA_real_,
      b = NA_real_,
      J = fit_info$J,
      parameters = NULL,
      target_K = fit_info$target_K,
      target_weight = spec,
      achieved_K = NULL,
      achieved_weight = NULL,
      K_loss = NA_real_,
      constraint_residual = NA_real_,
      constraint_slack = NA_real_,
      constraint_tolerance = list(
        abs = settings$tolerance$abs,
        rel = settings$tolerance$rel,
        effective = .dpprior_v2_effective_constraint_tolerance(
          spec, settings$tolerance
        ),
        formula = "abs + rel * max(abs(bound), 1e-8)"
      ),
      constraint_satisfied = FALSE,
      status = status,
      usable = FALSE,
      verified = verified,
      converged = FALSE,
      message = message,
      target = list(K = fit_info$target_K, weight = spec),
      achieved = NULL,
      residuals = list(
        constraint = list(
          raw = NA_real_, scaled = NA_real_, slack = NA_real_,
          units = spec$units,
          reason = "no verified candidate"
        )
      ),
      tolerances = list(
        constraint = list(
          absolute = settings$tolerance$abs,
          relative = settings$tolerance$rel,
          effective = .dpprior_v2_effective_constraint_tolerance(
            spec, settings$tolerance
          ),
          formula = "abs + rel * max(abs(bound), 1e-8)"
        ),
        verification = list(
          absolute = settings$control$verification_abs_tol,
          relative = settings$control$verification_rel_tol
        )
      ),
      attempts = attempts,
      verification = list(
        performed = verified,
        passed = verified,
        reason = if (verified) {
          "global_domain_infeasibility_certificate"
        } else {
          "candidate_verification_unavailable"
        },
        feasibility = feasibility
      ),
      feasibility = feasibility,
      feasibility_unknown = !isTRUE(feasibility$certified_infeasible),
      provenance = list(
        requested_method = "dual_anchor_hard_inequality",
        selected_method = NA_character_,
        is_fallback = FALSE,
        legacy = FALSE,
        lambda_used = FALSE,
        parameterization = "log(shape), log(rate)",
        log_bounds = settings$log_bounds,
        quadrature = list(
          M_fit = settings$M,
          M_verify = settings$M_verify,
          M_verify_required = settings$M_verify_required
        ),
        K_scaling = list(
          mu_K = max(abs(fit_info$target_K$mu_K), 1),
          var_K = max(abs(fit_info$target_K$var_K), 1),
          fixed_from_input = TRUE
        ),
        input_fit = fit_info$source
      ),
      .internal = list(
        fit_info = fit_info,
        spec = spec,
        settings = settings,
        candidates = candidates,
        verified_candidates = verified_candidates,
        selected = NULL
      )
    ),
    class = c("DPprior_dual_hard", "DPprior_fit")
  )
}


.dpprior_v2_hard_backend <- function(fit_info, spec, settings) {
  effective_tolerance <- .dpprior_v2_effective_constraint_tolerance(
    spec, settings$tolerance
  )
  attempts <- list()
  candidates <- list()
  add_candidate <- function(candidate) {
    if (is.list(candidate)) candidates[[length(candidates) + 1L]] <<- candidate
    invisible(candidate)
  }

  probe_capture <- .dpprior_v2_capture(
    .dpprior_v2_metric_extrema(
      spec, settings$log_bounds, settings$M, settings$M_verify,
      effective_tolerance
    )
  )
  if (inherits(probe_capture$value, "condition")) {
    feasibility <- list(
      classification = "unknown",
      certified_infeasible = FALSE,
      feasibility_unknown = TRUE,
      certificate = list(
        certified = FALSE,
        type = "feasibility_probe_error",
        error = probe_capture$error
      )
    )
    attempts[[1L]] <- .dpprior_v2_make_attempt(
      "analytic_monotonicity_feasibility_probe",
      bounds = list(log_a = settings$log_bounds,
                    log_b = settings$log_bounds),
      control = list(M = settings$M, M_verify = settings$M_verify),
      exit_code = 1L, message = "feasibility probe failed",
      elapsed = probe_capture$elapsed, warning = probe_capture$warning,
      error = probe_capture$error
    )
  } else {
    feasibility <- probe_capture$value
    attempts[[1L]] <- .dpprior_v2_make_attempt(
      "analytic_monotonicity_feasibility_probe",
      bounds = list(log_a = settings$log_bounds,
                    log_b = settings$log_bounds),
      control = list(M = settings$M, M_verify = settings$M_verify),
      exit_code = if (isTRUE(feasibility$certificate$certified)) 0L else 1L,
      message = feasibility$classification,
      objective = if (identical(spec$relation, "at_most")) {
        feasibility$minimum$refined$value %||% NA_real_
      } else {
        feasibility$maximum$refined$value %||% NA_real_
      },
      elapsed = probe_capture$elapsed,
      warning = probe_capture$warning,
      candidate = feasibility$feasible_corner$parameters %||% NULL
    )
  }

  if (isTRUE(feasibility$certified_infeasible)) {
    result <- .dpprior_v2_empty_hard_result(
      fit_info, spec, settings, feasibility, attempts, "infeasible",
      paste(
        "The requested hard inequality is certified infeasible on the",
        "declared log-parameter domain."
      ),
      candidates = candidates,
      verified_candidates = list()
    )
    return(result)
  }
  if (!is.null(feasibility$feasible_corner)) {
    add_candidate(.dpprior_v2_hard_candidate(
      feasibility$feasible_corner$eta, "feasibility_extreme", fit_info,
      spec, settings$M, effective_tolerance
    ))
  }

  start_eta <- log(c(fit_info$a, fit_info$b))
  K_objective <- function(eta) {
    moments <- .dpprior_v2_K_moments_safe(fit_info$J, eta, settings$M)
    if (!isTRUE(moments$ok)) return(.PENALTY_INF)
    .dpprior_v2_k_loss(
      list(mu_K = moments$mean, var_K = moments$var), fit_info$target_K
    )$value
  }

  unconstrained_capture <- .dpprior_v2_capture(stats::optim(
    par = start_eta,
    fn = K_objective,
    method = "L-BFGS-B",
    lower = rep(settings$log_bounds[[1L]], 2L),
    upper = rep(settings$log_bounds[[2L]], 2L),
    control = list(
      maxit = settings$control$maxit,
      factr = max(1, settings$control$optim_reltol / .Machine$double.eps)
    )
  ))
  if (inherits(unconstrained_capture$value, "condition")) {
    attempts[[length(attempts) + 1L]] <- .dpprior_v2_make_attempt(
      "K_only_L-BFGS-B",
      start = start_eta,
      bounds = list(lower = rep(settings$log_bounds[[1L]], 2L),
                    upper = rep(settings$log_bounds[[2L]], 2L)),
      control = list(maxit = settings$control$maxit,
                     reltol = settings$control$optim_reltol),
      exit_code = 1L, message = "optimizer raised an error",
      elapsed = unconstrained_capture$elapsed,
      warning = unconstrained_capture$warning,
      error = unconstrained_capture$error
    )
  } else {
    opt <- unconstrained_capture$value
    attempt_index <- length(attempts) + 1L
    unconstrained_candidate <- .dpprior_v2_hard_candidate(
      opt$par, "K_only_L-BFGS-B", fit_info, spec, settings$M,
      effective_tolerance, opt$convergence, attempt_index
    )
    attempts[[attempt_index]] <- .dpprior_v2_make_attempt(
      "K_only_L-BFGS-B",
      start = start_eta,
      bounds = list(lower = rep(settings$log_bounds[[1L]], 2L),
                    upper = rep(settings$log_bounds[[2L]], 2L)),
      control = list(maxit = settings$control$maxit,
                     reltol = settings$control$optim_reltol),
      exit_code = opt$convergence,
      message = opt$message %||% if (opt$convergence == 0L) "converged" else "not converged",
      counts = list(evaluations = unname(opt$counts[["function"]] %||% NA_integer_)),
      objective = opt$value,
      elapsed = unconstrained_capture$elapsed,
      warning = unconstrained_capture$warning,
      candidate = list(eta = opt$par, a = exp(opt$par[[1L]]),
                       b = exp(opt$par[[2L]]))
    )
    add_candidate(unconstrained_candidate)
  }

  # Always retain the untouched input candidate. It can demonstrate immediate
  # feasibility, but it has no optimizer-exit evidence of its own.
  add_candidate(.dpprior_v2_hard_candidate(
    start_eta, "input_fit", fit_info, spec, settings$M,
    effective_tolerance
  ))

  scan_started <- proc.time()[["elapsed"]]
  log_a_grid <- sort(unique(c(
    seq(settings$log_bounds[[1L]], settings$log_bounds[[2L]],
        length.out = settings$control$scan_points),
    start_eta[[1L]]
  )))
  scan_candidates <- lapply(log_a_grid, function(log_a) {
    .dpprior_v2_profile_candidate(
      log_a, fit_info, spec, settings, effective_tolerance,
      source = "deterministic_profile_scan"
    )
  })
  for (candidate in scan_candidates) add_candidate(candidate)
  finite_scan <- vapply(
    scan_candidates,
    function(x) isTRUE(x$valid) && is.finite(x$K$value), logical(1)
  )
  best_scan <- if (any(finite_scan)) {
    finite_scan_indices <- which(finite_scan)
    best_finite_position <- which.min(vapply(
      scan_candidates[finite_scan_indices], function(x) x$K$value, numeric(1)
    ))
    scan_candidates[[finite_scan_indices[[best_finite_position]]]]
  } else {
    NULL
  }
  attempts[[length(attempts) + 1L]] <- .dpprior_v2_make_attempt(
    "deterministic_feasible_profile_scan",
    start = start_eta,
    bounds = list(log_a = settings$log_bounds,
                  log_b = settings$log_bounds),
    control = list(
      scan_points = settings$control$scan_points,
      root_tol = settings$control$root_tol
    ),
    exit_code = if (any(finite_scan)) 0L else 1L,
    message = if (any(finite_scan)) "at least one feasible profile candidate" else "no feasible profile candidate",
    counts = list(evaluations = length(log_a_grid)),
    objective = best_scan$K$value %||% NA_real_,
    elapsed = proc.time()[["elapsed"]] - scan_started,
    candidate = if (is.null(best_scan)) NULL else {
      list(eta = best_scan$eta, a = best_scan$a, b = best_scan$b)
    }
  )

  # Refine several deterministic scan basins. The inner log-b optimization is
  # constrained to the exact feasible interval at every log-a evaluation.
  if (any(finite_scan)) {
    finite_indices <- which(finite_scan)
    ranked <- finite_indices[order(vapply(
      scan_candidates[finite_indices], function(x) x$K$value, numeric(1)
    ))]
    ranked <- utils::head(ranked, settings$control$profile_starts)
    for (idx in ranked) {
      grid_idx <- match(log_a_grid[[idx]], log_a_grid)
      bracket <- c(
        log_a_grid[[max(1L, grid_idx - 1L)]],
        log_a_grid[[min(length(log_a_grid), grid_idx + 1L)]]
      )
      if (diff(bracket) <= settings$control$root_tol) next
      evaluations <- 0L
      objective_a <- function(log_a) {
        evaluations <<- evaluations + 1L
        candidate <- .dpprior_v2_profile_candidate(
          log_a, fit_info, spec, settings, effective_tolerance,
          source = "constrained_profile_optimizer"
        )
        if (isTRUE(candidate$valid)) candidate$K$value else .PENALTY_INF
      }
      captured <- .dpprior_v2_capture(stats::optimize(
        objective_a, interval = bracket,
        tol = settings$control$optim_reltol
      ))
      attempt_index <- length(attempts) + 1L
      if (inherits(captured$value, "condition")) {
        attempts[[attempt_index]] <- .dpprior_v2_make_attempt(
          "constrained_profile_optimize",
          start = log_a_grid[[idx]],
          bounds = list(log_a = bracket, log_b = settings$log_bounds),
          control = list(
            root_tol = settings$control$root_tol,
            reltol = settings$control$optim_reltol
          ),
          exit_code = 1L, message = "profile optimizer raised an error",
          counts = list(evaluations = evaluations),
          elapsed = captured$elapsed, warning = captured$warning,
          error = captured$error
        )
      } else {
        prof <- captured$value
        exit_code <- if (is.finite(prof$objective) &&
                         prof$objective < .PENALTY_INF) 0L else 1L
        candidate <- .dpprior_v2_profile_candidate(
          prof$minimum, fit_info, spec, settings, effective_tolerance,
          source = "constrained_profile_optimizer",
          optimizer_exit_code = exit_code,
          attempt_index = attempt_index
        )
        attempts[[attempt_index]] <- .dpprior_v2_make_attempt(
          "constrained_profile_optimize",
          start = log_a_grid[[idx]],
          bounds = list(log_a = bracket, log_b = settings$log_bounds),
          control = list(
            root_tol = settings$control$root_tol,
            reltol = settings$control$optim_reltol
          ),
          exit_code = exit_code,
          message = if (exit_code == 0L) "finite constrained profile minimum" else "nonfinite constrained profile minimum",
          counts = list(evaluations = evaluations),
          objective = prof$objective,
          elapsed = captured$elapsed,
          warning = captured$warning,
          candidate = if (isTRUE(candidate$valid)) {
            list(eta = candidate$eta, a = candidate$a, b = candidate$b)
          } else NULL
        )
        add_candidate(candidate)
      }
    }
  }

  # A separately named penalty attempt is useful diagnostic evidence, but its
  # optimizer exit code never substitutes for verified inequality residuals.
  penalty_start <- if (!is.null(best_scan) && isTRUE(best_scan$valid)) {
    best_scan$eta
  } else {
    start_eta
  }
  penalty_objective <- function(eta) {
    K_value <- K_objective(eta)
    if (!is.finite(K_value) || K_value >= .PENALTY_INF) return(.PENALTY_INF)
    metric <- .dpprior_v2_eval_metric_safe(spec, eta, settings$M)
    residual <- .dpprior_v2_constraint_values(spec, metric, 0)$raw
    if (!is.finite(residual)) return(.PENALTY_INF)
    K_value + settings$control$penalty * max(residual, 0)^2
  }
  penalty_capture <- .dpprior_v2_capture(stats::optim(
    par = penalty_start,
    fn = penalty_objective,
    method = "L-BFGS-B",
    lower = rep(settings$log_bounds[[1L]], 2L),
    upper = rep(settings$log_bounds[[2L]], 2L),
    control = list(maxit = settings$control$maxit)
  ))
  attempt_index <- length(attempts) + 1L
  if (inherits(penalty_capture$value, "condition")) {
    attempts[[attempt_index]] <- .dpprior_v2_make_attempt(
      "penalty_L-BFGS-B_diagnostic",
      start = penalty_start,
      bounds = list(lower = rep(settings$log_bounds[[1L]], 2L),
                    upper = rep(settings$log_bounds[[2L]], 2L)),
      control = list(maxit = settings$control$maxit,
                     penalty = settings$control$penalty),
      exit_code = 1L, message = "penalty optimizer raised an error",
      elapsed = penalty_capture$elapsed, warning = penalty_capture$warning,
      error = penalty_capture$error
    )
  } else {
    opt <- penalty_capture$value
    penalty_candidate <- .dpprior_v2_hard_candidate(
      opt$par, "penalty_L-BFGS-B_diagnostic", fit_info, spec, settings$M,
      effective_tolerance, NA_integer_, attempt_index
    )
    penalty_candidate$diagnostic_optimizer_exit_code <- opt$convergence
    attempts[[attempt_index]] <- .dpprior_v2_make_attempt(
      "penalty_L-BFGS-B_diagnostic",
      start = penalty_start,
      bounds = list(lower = rep(settings$log_bounds[[1L]], 2L),
                    upper = rep(settings$log_bounds[[2L]], 2L)),
      control = list(maxit = settings$control$maxit,
                     penalty = settings$control$penalty),
      exit_code = opt$convergence,
      message = opt$message %||% if (opt$convergence == 0L) "converged" else "not converged",
      counts = list(evaluations = unname(opt$counts[["function"]] %||% NA_integer_)),
      objective = opt$value,
      elapsed = penalty_capture$elapsed,
      warning = penalty_capture$warning,
      candidate = list(eta = opt$par, a = exp(opt$par[[1L]]),
                       b = exp(opt$par[[2L]]))
    )
    add_candidate(penalty_candidate)
  }

  valid <- vapply(
    candidates,
    function(x) isTRUE(x$valid) && is.finite(x$K$value), logical(1)
  )
  if (!any(valid)) {
    feasibility$classification <- "unknown"
    feasibility$feasibility_unknown <- TRUE
    result <- .dpprior_v2_empty_hard_result(
      fit_info, spec, settings, feasibility, attempts, "failed",
      paste(
        "No finite hard-inequality candidate was produced; feasibility",
        "remains unknown."
      ),
      candidates = candidates,
      verified_candidates = list()
    )
    return(result)
  }

  valid_indices <- which(valid)
  ranked <- valid_indices[order(vapply(
    candidates[valid_indices], function(x) x$K$value, numeric(1)
  ))]
  optimizer_indices <- valid_indices[vapply(
    candidates[valid_indices],
    function(x) !is.na(x$optimizer_exit_code), logical(1)
  )]
  verify_indices <- unique(c(
    utils::head(ranked, settings$control$scan_keep), optimizer_indices
  ))
  verified_candidates <- lapply(verify_indices, function(index) {
    candidate <- candidates[[index]]
    verification <- .dpprior_v2_verify_candidate(
      candidate, fit_info, spec, settings, effective_tolerance
    )
    list(index = index, candidate = candidate, verification = verification)
  })
  candidate_execution_success <- function(candidate) {
    if (candidate$source %in% c(
      "deterministic_profile_scan", "input_fit", "feasibility_extreme"
    )) {
      return(TRUE)
    }
    attempt_index <- candidate$attempt_index
    if (!.dpprior_v2_plain_scalar(attempt_index, "numeric") ||
        attempt_index != as.integer(attempt_index) || attempt_index < 1L ||
        attempt_index > length(attempts)) {
      return(FALSE)
    }
    attempt <- attempts[[as.integer(attempt_index)]]
    error_message <- attempt$error
    has_error <- is.character(error_message) && length(error_message) == 1L &&
      !is.na(error_message) && nzchar(error_message)
    .dpprior_v2_exit_zero(attempt$exit_code) && !has_error
  }
  feasible_verified <- vapply(
    verified_candidates,
    function(x) isTRUE(x$verification$passed) &&
      isTRUE(x$verification$constraint_satisfied), logical(1)
  )
  diagnostic_verified <- vapply(
    verified_candidates,
    function(x) {
      verification <- x$verification
      candidate_execution_success(x$candidate) &&
        isTRUE(verification$performed) &&
        isTRUE(verification$numerical_passed) &&
        !isTRUE(verification$constraint_satisfied)
    },
    logical(1)
  )
  if (!any(feasible_verified) && !any(diagnostic_verified)) {
    feasibility$classification <- "unknown"
    feasibility$feasibility_unknown <- TRUE
    feasibility$candidate_verifications <- verified_candidates
    result <- .dpprior_v2_empty_hard_result(
      fit_info, spec, settings, feasibility, attempts, "failed",
      paste(
        "No candidate passed independent higher-order constraint",
        "verification; feasibility remains unknown."
      ),
      candidates = candidates,
      verified_candidates = verified_candidates
    )
    return(result)
  }

  diagnostic_selection <- !any(feasible_verified)
  eligible <- verified_candidates[if (diagnostic_selection) {
    diagnostic_verified
  } else {
    feasible_verified
  }]
  eligible_losses <- vapply(
    eligible, function(x) x$candidate$K$value, numeric(1)
  )
  minimum_candidate_loss <- min(eligible_losses)
  selection_tie_tolerance <- max(
    settings$control$verification_abs_tol,
    64 * .Machine$double.eps * max(1, abs(minimum_candidate_loss))
  )
  numerically_tied <- eligible_losses <=
    minimum_candidate_loss + selection_tie_tolerance
  selection_pool <- eligible[numerically_tied]
  selection_order <- order(
    -as.integer(vapply(
      selection_pool,
      function(x) {
        .dpprior_v2_exit_zero(x$candidate$optimizer_exit_code) &&
          .dpprior_v2_optimizer_evidence_allowed(x$candidate$source)
      },
      logical(1)
    )),
    vapply(selection_pool, function(x) x$index, integer(1))
  )
  selected <- selection_pool[[selection_order[[1L]]]]
  candidate <- selected$candidate
  verification <- selected$verification
  status <- if (diagnostic_selection) {
    "approximate"
  } else {
    .dpprior_v2_classify_hard_candidate(
      verification, candidate$optimizer_exit_code, candidate$eta,
      settings$log_bounds, settings$control$boundary_tol,
      source = candidate$source
    )
  }
  usable <- .dpprior_v2_hard_return_policy(
    status, settings$allow_approximate
  )$usable
  selected_attempt <- if (!is.na(candidate$attempt_index) &&
                          candidate$attempt_index <= length(attempts)) {
    attempts[[candidate$attempt_index]]$method
  } else {
    candidate$source
  }
  # Phase 6 selected-order quarantine: the public identity is tied to the
  # fit-order computation at returned (a,b). Refined values remain only in the
  # verification record and gate satisfaction/status jointly with this one.
  constraint <- candidate$constraint
  K <- candidate$K
  achieved_K_selected <- list(
    mu_K = candidate$moments$mean,
    var_K = candidate$moments$var,
    units = list(mu_K = "clusters", var_K = "clusters^2"),
    M = settings$M
  )
  achieved_weight_selected <- verification$achieved_weight_selected
  active_constraint <- abs(constraint$raw) <= constraint$tolerance
  domain_distance <- min(
    candidate$eta - settings$log_bounds[[1L]],
    settings$log_bounds[[2L]] - candidate$eta
  )
  message <- switch(
    status,
    converged = paste(
      "Hard inequality solution passed independent verification and",
      "the recorded optimality check."
    ),
    boundary = if (domain_distance <= settings$control$boundary_tol) {
      "Hard inequality solution passed verification at a declared log-parameter boundary."
    } else {
      "Hard inequality solution passed verification with the constraint active."
    },
    approximate = if (diagnostic_selection) {
      paste(
        "No independently verified feasible candidate was found; the",
        "returned finite candidate is a signed unsatisfied diagnostic only."
      )
    } else {
      paste(
        "A verified feasible candidate was found, but no successful constrained",
        "optimizer exit supports ordinary convergence."
      )
    },
    failed = "Candidate verification failed."
  )
  feasibility$classification <- if (diagnostic_selection) {
    "unknown"
  } else {
    "feasible_candidate"
  }
  feasibility$feasibility_unknown <- diagnostic_selection
  feasibility$candidate_count <- length(candidates)
  feasibility$verified_candidate_count <- length(verified_candidates)

  structure(
    list(
      mode = "hard_inequality",
      method = "dual_anchor_hard_inequality",
      a = candidate$a,
      b = candidate$b,
      J = fit_info$J,
      parameters = list(a = candidate$a, b = candidate$b),
      target_K = fit_info$target_K,
      target_weight = spec,
      achieved_K = achieved_K_selected,
      achieved_weight = achieved_weight_selected,
      K_loss = K$value,
      constraint_residual = constraint$raw,
      constraint_slack = constraint$slack,
      constraint_tolerance = list(
        abs = settings$tolerance$abs,
        rel = settings$tolerance$rel,
        effective = constraint$tolerance,
        formula = "abs + rel * max(abs(bound), 1e-8)"
      ),
      constraint_satisfied = isTRUE(verification$constraint_satisfied),
      status = status,
      usable = usable,
      verified = isTRUE(verification$passed),
      converged = identical(status, "converged"),
      message = message,
      target = list(K = fit_info$target_K, weight = spec),
      achieved = list(
        K = achieved_K_selected,
        weight = achieved_weight_selected
      ),
      residuals = list(
        K = list(
          raw = K$raw,
          scaled = K$scaled,
          scales = K$scales,
          scale_formula = K$scale_formula,
          units = c(mu_K = "clusters", var_K = "clusters^2")
        ),
        constraint = list(
          raw = constraint$raw,
          scaled = constraint$raw,
          slack = constraint$slack,
          units = spec$units,
          sign_convention = paste0(
            if (identical(spec$relation, "at_most")) {
              "achieved - bound"
            } else {
              "bound - achieved"
            }, "; feasible when raw <= tolerance"
          )
        )
      ),
      tolerances = list(
        constraint = list(
          absolute = settings$tolerance$abs,
          relative = settings$tolerance$rel,
          effective = constraint$tolerance,
          formula = "abs + rel * max(abs(bound), 1e-8)"
        ),
        verification = list(
          absolute = settings$control$verification_abs_tol,
          relative = settings$control$verification_rel_tol,
          perturbation_absolute = settings$control$perturbation_abs_tol
        )
      ),
      attempts = attempts,
      verification = c(
        verification,
        list(
          candidate_index = selected$index,
          candidates_checked = length(verified_candidates),
          verified_feasible_candidates = sum(feasible_verified),
          selection_rule = paste(
            "minimum selected-order fixed-input-scale K_loss among candidates",
            "whose selected and refined constraints both passed verification;",
            "successful optimizer evidence breaks numerical ties"
          ),
          minimum_candidate_K_loss = minimum_candidate_loss,
          selection_tie_tolerance = selection_tie_tolerance,
          selected_K_loss_delta_from_minimum =
            candidate$K$value - minimum_candidate_loss,
          all_candidate_results = verified_candidates
        )
      ),
      feasibility = feasibility,
      feasibility_unknown = FALSE,
      provenance = list(
        requested_method = "dual_anchor_hard_inequality",
        selected_method = selected_attempt,
        is_fallback = !(selected_attempt %in% c(
          "K_only_L-BFGS-B", "constrained_profile_optimize"
        )),
        fallback_reason = if (!(selected_attempt %in% c(
          "K_only_L-BFGS-B", "constrained_profile_optimize"
        ))) {
          "minimum verified feasible candidate came from another declared attempt"
        } else {
          NA_character_
        },
        legacy = FALSE,
        lambda_used = FALSE,
        parameterization = "log(shape), log(rate)",
        log_bounds = settings$log_bounds,
        boundary_distance_log_scale = domain_distance,
        constraint_active = active_constraint,
        optimality_evidence = if (
          .dpprior_v2_exit_zero(candidate$optimizer_exit_code) &&
          .dpprior_v2_optimizer_evidence_allowed(candidate$source)
        ) {
          "successful_selected_optimizer_exit_plus_independent_verification"
        } else {
          "unrefined_verified_feasible_candidate"
        },
        quadrature = list(
          M_fit = settings$M,
          M_verify = settings$M_verify,
          M_verify_required = settings$M_verify_required
        ),
        K_scaling = list(
          mu_K = K$scales[["mu_K"]],
          var_K = K$scales[["var_K"]],
          fixed_from_input = TRUE
        ),
        input_fit = fit_info$source
      ),
      .internal = list(
        fit_info = fit_info,
        spec = spec,
        settings = settings,
        candidates = candidates,
        verified_candidates = verified_candidates,
        selected = selected,
        diagnostic_selection = diagnostic_selection,
        minimum_candidate_loss = minimum_candidate_loss,
        selection_tie_tolerance = selection_tie_tolerance
      )
    ),
    class = c("DPprior_dual_hard", "DPprior_fit")
  )
}


#' Dual-anchor calibration with a verified hard inequality
#'
#' Minimizes the fixed-input-scale cluster-count loss subject to one named
#' weight inequality. This hard method has no lambda: lambda belongs only to
#' soft trade-off calibration.
#'
#' @param fit A `DPprior_fit` object retaining positive `a`, positive `b`, `J`,
#'   and the requested K mean and variance.
#' @param constraint A named list with `metric`, `relation`, and `bound`.
#'   Supported metrics are `"wsb_tail"` (requires `threshold`),
#'   `"wsb_mean"`, `"wsb_quantile"` (requires `probability`), and
#'   `"wmax_tail_upper"` (requires `threshold`). Relations may be `"<="` or
#'   `"at_most"`, and `">="` or `"at_least"` except that a certified W-max
#'   upper bound supports only an at-most safety constraint.
#' @param constraint_tol Named list containing non-negative `abs` and `rel`
#'   tolerances. Independent satisfaction uses
#'   `abs + rel * max(abs(bound), 1e-8)`.
#' @param M Fit-order Gauss--Laguerre quadrature order, no greater than 256 so
#'   that a required higher-order verification remains available.
#' @param M_verify Independent verification order. It must be at least
#'   `max(2*M, M+40)` and no greater than 512. The default is that required
#'   minimum.
#' @param log_bounds Two finite increasing bounds shared by `log(a)` and
#'   `log(b)`.
#' @param control Named solver-control list. Supported components are
#'   `maxit`, `scan_points`, `scan_keep`, `profile_starts`, `root_tol`,
#'   `optim_reltol`, `penalty`, `boundary_tol`,
#'   `verification_abs_tol`, `verification_rel_tol`, `perturbation_step`, and
#'   `perturbation_abs_tol`.
#' @param allow_approximate Whether to return an independently verified
#'   feasible scan candidate without successful optimizer evidence. The
#'   result remains visibly `status = "approximate"` and intrinsically
#'   `usable = FALSE`; this flag changes return-versus-condition behavior only.
#' @param ... Must be empty. In particular, supplying `lambda` raises a typed
#'   invalid-input condition before fitting.
#'
#' @return A canonical \code{dpprior.result/1} \code{DPprior_dual_hard}
#'   object for an approved solution. The \code{constraint} extension contains
#'   the named inequality, achieved metric, signed residual, tolerance, and
#'   independent satisfaction certificate. Certified infeasibility, unknown
#'   feasibility, and an unapproved approximate result raise a typed condition
#'   whose \code{condition$result} contains the complete canonical object.
#'
#' @details
#' The signed residual is `achieved - bound` for an at-most constraint and
#' `bound - achieved` for an at-least constraint. Thus feasibility always has
#' the auditable form `constraint_residual <= constraint_tolerance$effective`.
#' `constraint_satisfied` is set only from a fresh higher-order calculation;
#' optimizer exit codes are never treated as constraint evidence.
#' Decision-ready use requires the returned object's mode-specific contract
#' together with \code{usable = TRUE} and \code{verified = TRUE}.
#'
#' @family elicitation
#' @export
DPprior_dual_hard <- function(
    fit,
    constraint = list(
      metric = "wsb_tail", threshold = 0.5,
      relation = "<=", bound = 0.30
    ),
    constraint_tol = list(abs = 1e-6, rel = 1e-6),
    M = .QUAD_NODES_DEFAULT,
    M_verify = NULL,
    log_bounds = .LOG_BOUNDS_DEFAULT,
    control = list(),
    allow_approximate = FALSE,
    ...) {
  dots <- match.call(expand.dots = FALSE)$...
  if (!is.null(dots) && length(dots)) {
    dot_names <- names(dots)
    if (is.null(dot_names)) dot_names <- rep("<unnamed>", length(dots))
    dot_names[!nzchar(dot_names)] <- "<unnamed>"
    if ("lambda" %in% dot_names) {
      .dpprior_v2_abort_invalid(
        "lambda is not defined for hard inequality calibration; use DPprior_dual_soft() for a trade-off",
        "lambda", NULL, "omit lambda in hard mode", "lambda_in_hard_mode",
        c("dpprior_dual_lambda_error", "dpprior_conflicting_input")
      )
    }
    .dpprior_v2_abort_invalid(
      sprintf("unknown DPprior_dual_hard argument(s): %s",
              paste(unique(dot_names), collapse = ", ")),
      "...", dot_names, "no additional arguments", "unknown_control",
      c("dpprior_dual_control_error", "dpprior_unknown_control_error")
    )
  }

  fit_info <- .dpprior_v2_normalize_fit(fit)
  .dpprior_v2_require_decision_ready_fit(fit_info, mode = "hard")
  spec <- .dpprior_v2_normalize_weight_spec(constraint, mode = "hard")
  settings <- .dpprior_v2_normalize_hard_controls(
    constraint_tol, M, M_verify, log_bounds, control, allow_approximate
  )
  start_eta <- log(c(fit_info$a, fit_info$b))
  if (any(start_eta < settings$log_bounds[[1L]]) ||
      any(start_eta > settings$log_bounds[[2L]])) {
    .dpprior_v2_abort_invalid(
      paste(
        "fit parameters lie outside the declared log_bounds; hard mode does",
        "not silently project or nudge the input fit"
      ),
      "fit", list(a = fit_info$a, b = fit_info$b),
      sprintf("log(a), log(b) in [%g,%g]", settings$log_bounds[[1L]],
              settings$log_bounds[[2L]]),
      "input_outside_solver_domain",
      c("dpprior_dual_fit_error", "dpprior_bounds_error")
    )
  }

  backend_result <- .dpprior_v2_hard_backend(fit_info, spec, settings)
  result <- .dpprior_v2_hard_canonical_result(backend_result)
  raw_result <- unclass(result)
  return_policy <- .dpprior_v2_hard_return_policy(
    raw_result[["status", exact = TRUE]], settings$allow_approximate
  )
  if (isTRUE(return_policy$return_result)) {
    return(result)
  }
  stop(.dpprior_v2_condition_from_result(result))
}
