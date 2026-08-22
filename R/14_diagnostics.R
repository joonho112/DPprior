# =============================================================================
# Module 14: Estimand-labelled diagnostics and provenance
# =============================================================================
#
# A failed K_J PMF is never replaced by invented data. Weight diagnostics
# distinguish the first size-biased DP weight (W_SB) from the largest
# population weight (W_max), and v2 objects contain no categorical
# "dominance risk" label.
# =============================================================================

.DIAGNOSTIC_STATUS_CODES <- c(
  "converged", "boundary", "approximate", "infeasible", "failed"
)

.diagnostic_abort_numerical <- function(message, result = NULL,
                                        cause = NULL,
                                        code = "diagnostic_failure",
                                        subclass = character()) {
  stop(.dpprior_new_condition(
    message = message,
    classes = c(
      subclass,
      "dpprior_diagnostics_error", "dpprior_numerical_error",
      "dpprior_error", "error"
    ),
    result = result, cause = cause, code = code
  ))
}

.diagnostic_validate_status <- function(status, name = "status") {
  if (!is.character(status) || length(status) != 1L || is.na(status) ||
      !status %in% .DIAGNOSTIC_STATUS_CODES) {
    .dpprior_abort_invalid(
      sprintf("%s must be one of: %s", name,
              paste(.DIAGNOSTIC_STATUS_CODES, collapse = ", ")),
      c("dpprior_diagnostic_status_error", "dpprior_bounds_error"),
      name, status, paste(.DIAGNOSTIC_STATUS_CODES, collapse = ", "),
      "unknown_status"
    )
  }
  status
}

.diagnostic_combine_status <- function(statuses) {
  statuses <- vapply(
    statuses, .diagnostic_validate_status, character(1),
    name = "component status"
  )
  if (any(statuses == "failed")) return("failed")
  if (any(statuses == "infeasible")) return("infeasible")
  if (any(statuses == "approximate")) return("approximate")
  if (any(statuses == "boundary")) return("boundary")
  "converged"
}

.diagnostic_controls <- function(M, M_verify = NULL,
                                 abs_tol = 1e-10, rel_tol = 1e-8) {
  controls <- .marginal_verification_controls(
    M, M_verify = M_verify, abs_tol = abs_tol, rel_tol = rel_tol,
    strict = FALSE
  )
  if (is.null(M_verify) && isTRUE(controls$verification_available)) {
    controls <- .marginal_verification_controls(
      M, M_verify = controls$M_verification_required,
      abs_tol = abs_tol, rel_tol = rel_tol, strict = FALSE
    )
  }
  controls
}

.diagnostic_scalar_audit <- function(selected, verification,
                                     abs_tol, rel_tol) {
  difference <- if (is.null(verification)) NA_real_ else {
    abs(selected - verification)
  }
  tolerance <- if (is.null(verification)) NA_real_ else {
    abs_tol + rel_tol * max(abs(selected), abs(verification))
  }
  list(
    selected = as.numeric(selected),
    verification = if (is.null(verification)) NA_real_ else {
      as.numeric(verification)
    },
    difference = difference,
    tolerance = tolerance,
    passed = if (is.null(verification)) NA else {
      is.finite(difference) && is.finite(tolerance) &&
        difference <= tolerance
    }
  )
}

.diagnostic_attempt <- function(name, method, fn) {
  started <- proc.time()[["elapsed"]]
  caught <- NULL
  value <- tryCatch(fn(), error = function(e) {
    caught <<- e
    NULL
  })
  elapsed <- proc.time()[["elapsed"]] - started
  if (!is.null(caught)) {
    return(list(
      ok = FALSE, value = NULL, error = caught,
      record = list(
        component = name, method = method, status = "failed",
        elapsed_seconds = as.numeric(elapsed),
        message = conditionMessage(caught),
        condition_classes = class(caught)
      )
    ))
  }
  status <- if (is.list(value) && !is.null(value$status)) {
    .diagnostic_validate_status(value$status, paste0(name, "$status"))
  } else {
    "converged"
  }
  list(
    ok = TRUE, value = value, error = NULL,
    record = list(
      component = name, method = method, status = status,
      elapsed_seconds = as.numeric(elapsed),
      message = if (is.list(value) && !is.null(value$message)) {
        value$message
      } else {
        "completed"
      },
      condition_classes = character(0)
    )
  )
}

.diagnostic_tail_names <- function(thresholds, prefix = "threshold_") {
  paste0(prefix, format(
    thresholds, scientific = FALSE, trim = TRUE, digits = 15L
  ))
}

.diagnostic_wsb_tail <- function(threshold, a, b) {
  # Adapter seam for the canonical alias introduced by the weight backend.
  if (exists("prob_wsb_exceeds", mode = "function", inherits = TRUE)) {
    return(prob_wsb_exceeds(threshold, a, b))
  }
  prob_w1_exceeds(threshold, a, b)
}

.diagnostic_wmax_backend <- function(thresholds, a, b) {
  # G8 deliberately excludes a direct W_max point-probability claim from the
  # canonical diagnostics API. Certified upper bounds remain available from
  # the dedicated wmax_tail_bounds() API; this compatibility component only
  # records that no W_SB value was substituted for W_max.
  list(
    requested = FALSE, estimand = "W_max",
    label = "Largest population DP weight",
    status = "failed", usable = FALSE, verified = FALSE,
    message = paste(
      "W_max point-probability calibration is not part of canonical prior",
      "diagnostics; use wmax_tail_bounds() for certified bounds."
    ),
    values = NULL,
    provenance = list(
      requested_method = "canonical_diagnostics",
      selected_method = NULL,
      is_fallback = FALSE,
      policy = "backend_unavailable",
      certified_upper_bound_api = "wmax_tail_bounds"
    )
  )
}

# =============================================================================
# Alpha diagnostics
# =============================================================================

#' Alpha Distribution Diagnostics
#'
#' @param a Numeric Gamma shape.
#' @param b Numeric Gamma rate.
#' @return A provenance-labelled diagnostic component.
#' @keywords internal
compute_alpha_diagnostics <- function(a, b) {
  assert_positive(a, "a")
  assert_positive(b, "b")
  probs <- c(0.05, 0.25, 0.50, 0.75, 0.95)
  quantiles <- stats::qgamma(probs, shape = a, rate = b)
  names(quantiles) <- paste0(
    "q", formatC(100 * probs, format = "f", digits = 0)
  )
  list(
    status = "converged", usable = TRUE, verified = TRUE,
    message = "Gamma summaries were evaluated from closed-form identities.",
    estimand = "alpha", label = "DP concentration parameter",
    mean = a / b, sd = sqrt(a) / b, cv = 1 / sqrt(a),
    median = unname(quantiles[["q50"]]), quantiles = quantiles,
    method = "closed_form_and_stats_qgamma",
    provenance = list(
      requested_method = "closed_form",
      selected_method = "closed_form_and_stats_qgamma",
      is_fallback = FALSE, approximation = FALSE
    )
  )
}

# =============================================================================
# K_J diagnostics: no invented PMF fallback
# =============================================================================

#' Get K_J PMF and status metadata
#'
#' @param J Integer sample size.
#' @param a,b Numeric Gamma hyperparameters.
#' @param M Selected quadrature order.
#' @param M_verify Optional verification order.
#' @param abs_tol,rel_tol Comparison tolerances.
#' @return PMF, support, optional verification PMF, and provenance.
#' @keywords internal
.get_K_pmf_support <- function(
    J, a, b, M = .QUAD_NODES_DEFAULT, M_verify = NULL,
    abs_tol = 1e-10, rel_tol = 1e-8) {
  assert_valid_J(J)
  assert_positive(a, "a")
  assert_positive(b, "b")
  J <- as.integer(J)
  controls <- .diagnostic_controls(M, M_verify, abs_tol, rel_tol)
  if (!exists("pmf_K_marginal", mode = "function", inherits = TRUE)) {
    .diagnostic_abort_numerical(
      "pmf_K_marginal() is unavailable; K_J PMF diagnostics cannot be computed",
      code = "pmf_backend_unavailable"
    )
  }
  if (!exists("compute_log_stirling", mode = "function", inherits = TRUE)) {
    .diagnostic_abort_numerical(
      "compute_log_stirling() is unavailable; K_J PMF diagnostics cannot be computed",
      code = "stirling_backend_unavailable"
    )
  }
  logS <- compute_log_stirling(J)
  pmf0 <- pmf_K_marginal(
    J, a, b, logS, M = controls$M,
    M_verify = controls$M_verify,
    abs_tol = controls$abs_tol, rel_tol = controls$rel_tol,
    strict = FALSE
  )
  metadata <- attr(pmf0, "marginal_metadata", exact = TRUE)
  verification0 <- attr(
    pmf0, ".marginal_verification_pmf", exact = TRUE
  )
  if (length(pmf0) != J + 1L) {
    .diagnostic_abort_numerical(
      sprintf("pmf_K_marginal() returned length %d; expected %d",
              length(pmf0), J + 1L),
      code = "pmf_length"
    )
  }
  .dpprior_validate_pmf(pmf0, expected_length = J + 1L)
  if (!identical(as.numeric(pmf0[[1L]]), 0)) {
    .diagnostic_abort_numerical(
      "K_J PMF assigned nonzero mass to impossible support value K_J=0",
      code = "pmf_support"
    )
  }
  if (is.null(metadata) || is.null(metadata$status)) {
    .diagnostic_abort_numerical(
      "K_J PMF did not carry required marginal status metadata",
      code = "pmf_metadata"
    )
  }
  .diagnostic_validate_status(metadata$status, "PMF metadata status")
  if (!is.null(verification0)) {
    .dpprior_validate_pmf(verification0, expected_length = J + 1L)
    if (!identical(as.numeric(verification0[[1L]]), 0)) {
      .diagnostic_abort_numerical(
        "verification PMF assigned nonzero mass to K_J=0",
        code = "verification_pmf_support"
      )
    }
  }
  list(
    pmf = unname(as.numeric(pmf0[-1L])),
    verification_pmf = if (is.null(verification0)) NULL else {
      unname(as.numeric(verification0[-1L]))
    },
    support = seq_len(J), status = metadata$status,
    metadata = metadata, controls = controls
  )
}

.diagnostic_discrete_K <- function(pmf, support, probs) {
  cdf <- cumsum(pmf)
  quantiles <- vapply(probs, function(p) {
    index <- which(cdf >= p)[1L]
    if (is.na(index)) {
      .diagnostic_abort_numerical(
        sprintf("K_J PMF CDF did not reach probability %.17g", p),
        code = "pmf_quantile"
      )
    }
    as.integer(support[[index]])
  }, integer(1))
  names(quantiles) <- paste0(
    "q", formatC(100 * probs, format = "f", digits = 0)
  )
  list(
    quantiles = quantiles,
    median = as.integer(quantiles[["q50"]]),
    mode = as.integer(support[[which.max(pmf)]])
  )
}

#' K Distribution Diagnostics
#'
#' @param J Integer sample size.
#' @param a,b Numeric Gamma hyperparameters.
#' @param M Selected quadrature order.
#' @param M_verify Optional verification order.
#' @param abs_tol,rel_tol Comparison tolerances.
#' @return A canonical K_J diagnostic component.
#' @keywords internal
compute_K_diagnostics <- function(
    J, a, b, M = .QUAD_NODES_DEFAULT, M_verify = NULL,
    abs_tol = 1e-10, rel_tol = 1e-8) {
  assert_valid_J(J)
  assert_positive(a, "a")
  assert_positive(b, "b")
  J <- as.integer(J)
  controls <- .diagnostic_controls(M, M_verify, abs_tol, rel_tol)
  moments <- exact_K_moments(
    J, a, b, M = controls$M, M_verify = controls$M_verify,
    abs_tol = controls$abs_tol, rel_tol = controls$rel_tol,
    strict = FALSE
  )
  pmf_obj <- .get_K_pmf_support(
    J, a, b, M = controls$M, M_verify = controls$M_verify,
    abs_tol = controls$abs_tol, rel_tol = controls$rel_tol
  )
  probs <- c(0.05, 0.25, 0.50, 0.75, 0.95)
  selected_discrete <- .diagnostic_discrete_K(
    pmf_obj$pmf, pmf_obj$support, probs
  )
  verification_discrete <- if (is.null(pmf_obj$verification_pmf)) NULL else {
    .diagnostic_discrete_K(
      pmf_obj$verification_pmf, pmf_obj$support, probs
    )
  }
  discrete_passed <- if (is.null(verification_discrete)) NA else {
    identical(
      unname(selected_discrete$quantiles),
      unname(verification_discrete$quantiles)
    ) && identical(selected_discrete$mode, verification_discrete$mode)
  }
  moment_verified <- isTRUE(moments$quadrature$verification_passed)
  pmf_verified <- isTRUE(pmf_obj$metadata$verification$passed)
  discrete_verified <- isTRUE(discrete_passed)
  verification_complete <-
    moment_verified && pmf_verified && discrete_verified
  status <- .diagnostic_combine_status(list(moments$status, pmf_obj$status))
  if (status %in% c("converged", "boundary") &&
      !verification_complete) {
    status <- "approximate"
  }
  verified <- status %in% c("converged", "boundary") &&
    verification_complete
  list(
    status = status,
    usable = status %in% c("converged", "boundary"),
    verified = verified,
    message = if (verified) {
      "K_J moments, PMF, and discrete summaries passed higher-order verification."
    } else {
      paste(
        "K_J diagnostics are approximate; inspect moment, PMF, and discrete",
        "verification records."
      )
    },
    estimand = "K_J",
    label = sprintf("Occupied cluster count among J=%d exchangeable units", J),
    mean = as.numeric(moments$mean), var = as.numeric(moments$var),
    sd = as.numeric(moments$sd), cv = as.numeric(moments$cv),
    mode = selected_discrete$mode, median = selected_discrete$median,
    quantiles = selected_discrete$quantiles, support = pmf_obj$support,
    pmf = pmf_obj$pmf,
    method = "gauss-laguerre-marginal-pmf-and-moments",
    verification = list(
      moments = moments$quadrature,
      pmf = pmf_obj$metadata$verification,
      discrete = list(
        performed = !is.null(verification_discrete),
        passed = discrete_passed, probabilities = probs,
        selected = selected_discrete,
        verification = verification_discrete
      )
    ),
    residuals = list(
      pmf_mass_error = abs(sum(pmf_obj$pmf) - 1),
      mean_refinement_difference = moments$quadrature$mean_difference,
      variance_refinement_difference = moments$quadrature$variance_difference
    ),
    tolerances = list(
      pmf_mass = .TOL_PMF_SUM,
      absolute = controls$abs_tol, relative = controls$rel_tol
    ),
    provenance = list(
      requested_method = "marginal_K_diagnostic",
      selected_method = "gauss-laguerre-marginal-pmf-and-moments",
      is_fallback = FALSE,
      approximations = if (status == "approximate") {
        c(
          moments$quadrature$reason, pmf_obj$metadata$reason,
          if (!moment_verified) "moment_verification_incomplete" else NULL,
          if (!pmf_verified) "pmf_verification_incomplete" else NULL,
          if (!discrete_verified) "discrete_verification_incomplete" else NULL
        )
      } else character(0),
      selected_order = controls$M,
      verification_order = controls$M_verify,
      verification_order_required = controls$M_verification_required,
      verification_available = controls$verification_available
    )
  )
}

# =============================================================================
# Weight diagnostics: W_SB and W_max are different estimands
# =============================================================================

#' First size-biased and largest-weight diagnostics
#'
#' @param a,b Numeric Gamma hyperparameters.
#' @param thresholds Numeric vector in (0,1).
#' @param M Selected quadrature order for E(W_SB).
#' @param M_verify Optional verification order.
#' @param abs_tol,rel_tol Comparison tolerances.
#' @return Separate \code{size_biased} and \code{maximum} components. Numeric
#'   compatibility aliases refer only to W_SB; no risk category is returned.
#' @family diagnostics
#' @export
compute_weight_diagnostics <- function(
    a, b, thresholds = c(0.3, 0.5, 0.7, 0.9),
    M = .QUAD_NODES_DEFAULT, M_verify = NULL,
    abs_tol = 1e-10, rel_tol = 1e-8) {
  assert_positive(a, "a")
  assert_positive(b, "b")
  thresholds <- .dpprior_validate_probability(
    thresholds, "thresholds", scalar = FALSE, open = TRUE
  )
  controls <- .diagnostic_controls(M, M_verify, abs_tol, rel_tol)
  mean_selected <- mean_w1(a, b, controls$M)
  mean_verification <- if (is.null(controls$M_verify)) NULL else {
    mean_w1(a, b, controls$M_verify)
  }
  mean_audit <- .diagnostic_scalar_audit(
    mean_selected, mean_verification,
    controls$abs_tol, controls$rel_tol
  )
  size_biased_status <- if (isTRUE(mean_audit$passed)) {
    "converged"
  } else {
    "approximate"
  }
  probs <- c(0.05, 0.25, 0.50, 0.75, 0.95)
  quantiles <- quantile_w1(probs, a, b)
  names(quantiles) <- paste0(
    "q", formatC(100 * probs, format = "f", digits = 0)
  )
  tail_probability <- vapply(
    thresholds, .diagnostic_wsb_tail, numeric(1), a = a, b = b
  )
  names(tail_probability) <- .diagnostic_tail_names(thresholds)
  size_biased <- list(
    status = size_biased_status,
    usable = size_biased_status == "converged",
    verified = isTRUE(mean_audit$passed),
    message = if (isTRUE(mean_audit$passed)) {
      "W_SB closed-form distribution and refined quadrature mean are available."
    } else {
      "W_SB mean is unverified or disagrees at the required refinement order."
    },
    estimand = "W_SB", label = "First size-biased DP weight",
    interpretation = paste(
      "Mass of the population cluster containing a randomly selected unit;",
      "not the largest population weight"
    ),
    mean = as.numeric(mean_selected),
    median = unname(quantiles[["q50"]]), quantiles = quantiles,
    thresholds = thresholds, tail_probability = tail_probability,
    units = "probability",
    method = "closed-form-distribution-plus-refined-quadrature-mean",
    verification = list(mean = mean_audit),
    provenance = list(
      requested_method = "W_SB_marginal",
      selected_method = "closed_form_and_gauss_laguerre",
      is_fallback = FALSE,
      selected_order = controls$M,
      verification_order = controls$M_verify,
      verification_order_required = controls$M_verification_required,
      verification_available = controls$verification_available
    )
  )
  maximum <- .diagnostic_wmax_backend(thresholds, a, b)
  statuses <- list(size_biased$status)
  if (isTRUE(maximum$requested)) statuses <- c(statuses, list(maximum$status))
  status <- .diagnostic_combine_status(statuses)
  legacy_prob_names <- paste0(
    "prob_gt_",
    format(thresholds, scientific = FALSE, trim = TRUE, digits = 15L)
  )
  compatibility_prob <- unname(tail_probability)
  names(compatibility_prob) <- legacy_prob_names
  list(
    status = status,
    usable = status %in% c("converged", "boundary"),
    verified = isTRUE(size_biased$verified) &&
      (!isTRUE(maximum$requested) || isTRUE(maximum$verified)),
    message = if (status == "converged") {
      "Requested weight diagnostics passed their method contracts."
    } else {
      "One or more weight diagnostics are approximate or unavailable."
    },
    size_biased = size_biased, maximum = maximum,
    mean = size_biased$mean, median = size_biased$median,
    quantiles = size_biased$quantiles, prob_exceeds = compatibility_prob,
    provenance = list(
      schema_version = 2L, estimand_separation = TRUE,
      is_fallback = FALSE,
      removed_ambiguous_category = "dominance_risk",
      compatibility_alias_owner = "W_SB"
    )
  )
}

# =============================================================================
# Co-clustering diagnostics
# =============================================================================

#' Co-Clustering Diagnostics (rho)
#'
#' @param a,b Numeric Gamma hyperparameters.
#' @param M Selected quadrature order.
#' @param M_verify Optional verification order.
#' @param abs_tol,rel_tol Comparison tolerances.
#' @return Status-aware rho component without a qualitative category.
#' @keywords internal
compute_coclustering_diagnostics <- function(
    a, b, M = .QUAD_NODES_DEFAULT, M_verify = NULL,
    abs_tol = 1e-10, rel_tol = 1e-8) {
  assert_positive(a, "a")
  assert_positive(b, "b")
  controls <- .diagnostic_controls(M, M_verify, abs_tol, rel_tol)
  mean_selected <- mean_rho(a, b, controls$M)
  var_selected <- var_rho(a, b, controls$M)
  if (!is.finite(mean_selected) || mean_selected < 0 || mean_selected > 1 ||
      !is.finite(var_selected) || var_selected < 0) {
    .diagnostic_abort_numerical(
      "rho diagnostics violated probability or non-negative variance bounds",
      code = "rho_numerical_contract"
    )
  }
  mean_verification <- var_verification <- NULL
  if (!is.null(controls$M_verify)) {
    mean_verification <- mean_rho(a, b, controls$M_verify)
    var_verification <- var_rho(a, b, controls$M_verify)
  }
  mean_audit <- .diagnostic_scalar_audit(
    mean_selected, mean_verification, controls$abs_tol, controls$rel_tol
  )
  var_audit <- .diagnostic_scalar_audit(
    var_selected, var_verification, controls$abs_tol, controls$rel_tol
  )
  verified <- isTRUE(mean_audit$passed) && isTRUE(var_audit$passed)
  status <- if (verified) "converged" else "approximate"
  list(
    status = status, usable = status == "converged", verified = verified,
    message = if (verified) {
      "rho mean and variance passed higher-order quadrature verification."
    } else {
      "rho moments are unverified or disagree at the required refinement order."
    },
    estimand = "rho",
    label = "Conditional pairwise co-clustering probability",
    mean = as.numeric(mean_selected), var = as.numeric(var_selected),
    sd = sqrt(as.numeric(var_selected)), units = "probability",
    method = "gauss-laguerre-marginal-moments",
    verification = list(mean = mean_audit, variance = var_audit),
    provenance = list(
      requested_method = "rho_marginal_moments",
      selected_method = "gauss_laguerre", is_fallback = FALSE,
      selected_order = controls$M,
      verification_order = controls$M_verify,
      verification_order_required = controls$M_verification_required,
      verification_available = controls$verification_available
    )
  )
}

# =============================================================================
# Explicit warning policy and legacy adapter
# =============================================================================

.diagnostic_validate_warning_policy <- function(policy) {
  if (is.null(policy)) return(NULL)
  required <- c(
    "estimand", "direction", "weight_threshold", "action_threshold"
  )
  ordinary <- typeof(policy) == "list" && is.list(policy) &&
    !is.object(policy) && !isS4(policy) && is.null(dim(policy)) &&
    identical(names(attributes(policy)), "names") &&
    identical(names(policy), required)
  if (!ordinary) {
    .dpprior_abort_invalid(
      paste(
        "warning_policy must be NULL or an ordinary list with exact fields:",
        paste(required, collapse = ", ")
      ),
      c(
        "dpprior_warning_policy_error", "dpprior_diagnostics_input_error",
        "dpprior_type_error"
      ),
      "warning_policy", policy,
      paste(required, collapse = ", "), "warning_policy_shape"
    )
  }
  estimand <- policy[["estimand", exact = TRUE]]
  direction <- policy[["direction", exact = TRUE]]
  if (!is.character(estimand) || is.object(estimand) ||
      !is.null(dim(estimand)) || length(estimand) != 1L ||
      is.na(estimand) || !estimand %in% c("W_SB", "W_max")) {
    .dpprior_abort_invalid(
      "warning_policy$estimand must be 'W_SB' or 'W_max'",
      c(
        "dpprior_warning_policy_error", "dpprior_diagnostics_input_error",
        "dpprior_bounds_error"
      ),
      "warning_policy$estimand", estimand, "W_SB or W_max",
      "unknown_estimand"
    )
  }
  if (!is.character(direction) || is.object(direction) ||
      !is.null(dim(direction)) || length(direction) != 1L ||
      is.na(direction) || !direction %in% c("above", "below")) {
    .dpprior_abort_invalid(
      "warning_policy$direction must be 'above' or 'below'",
      c(
        "dpprior_warning_policy_error", "dpprior_diagnostics_input_error",
        "dpprior_bounds_error"
      ),
      "warning_policy$direction", direction, "above or below",
      "unknown_direction"
    )
  }
  validate_probability <- function(field, open) {
    value <- policy[[field, exact = TRUE]]
    valid <- is.numeric(value) && !is.object(value) && is.null(dim(value)) &&
      length(value) == 1L && !is.na(value) && is.finite(value) &&
      if (open) value > 0 && value < 1 else value >= 0 && value <= 1
    if (!valid) {
      .dpprior_abort_invalid(
        sprintf(
          "warning_policy$%s must be one ordinary probability %s.",
          field, if (open) "strictly between 0 and 1" else "in [0,1]"
        ),
        c(
          "dpprior_warning_policy_error", "dpprior_diagnostics_input_error",
          "dpprior_probability_error", "dpprior_bounds_error"
        ),
        paste0("warning_policy$", field), value,
        if (open) "ordinary scalar in (0,1)" else "ordinary scalar in [0,1]",
        "warning_policy_probability"
      )
    }
    as.numeric(value)
  }
  weight_threshold <- validate_probability("weight_threshold", TRUE)
  action_threshold <- validate_probability("action_threshold", FALSE)
  list(
    estimand = estimand,
    direction = direction,
    weight_threshold = weight_threshold,
    action_threshold = action_threshold
  )
}

.diagnostic_policy_evidence <- function(policy, a, b) {
  if (is.null(policy)) {
    return(list(records = list(), messages = character(0)))
  }
  estimand <- policy[["estimand", exact = TRUE]]
  direction <- policy[["direction", exact = TRUE]]
  weight_threshold <- policy[["weight_threshold", exact = TRUE]]
  action_threshold <- policy[["action_threshold", exact = TRUE]]
  if (identical(estimand, "W_SB")) {
    value <- as.numeric(.diagnostic_wsb_tail(weight_threshold, a, b))
    outcome <- if (identical(direction, "above")) {
      if (value > action_threshold) "triggered" else "not_triggered"
    } else if (value < action_threshold) {
      "triggered"
    } else {
      "not_triggered"
    }
    record <- list(
      estimand = "W_SB", direction = direction,
      threshold = action_threshold, value = value,
      lower = NULL, upper = NULL, outcome = outcome,
      basis = "exact_tail_probability"
    )
  } else {
    record <- list(
      estimand = "W_max", direction = direction,
      threshold = action_threshold, value = NULL,
      lower = NULL, upper = NULL, outcome = "indeterminate",
      basis = "backend_unavailable"
    )
  }
  messages <- if (identical(record[["outcome", exact = TRUE]],
                            "triggered")) {
    messages <- sprintf(
      paste(
        "Explicit W_SB policy triggered: P(W_SB > %.6g) is %.6g and",
        "the requested comparison is '%s' action threshold %.6g."
      ),
      weight_threshold, record[["value", exact = TRUE]], direction,
      action_threshold
    )
  } else character()
  list(records = list(record), messages = messages)
}

#' Deprecated Ambiguous Weight-Risk Check
#'
#' This compatibility wrapper evaluates the explicit policy
#' \code{P(W_SB > threshold) > risk_level}; it does not evaluate W_max or
#' assign a category.
#'
#' @param a,b Numeric Gamma hyperparameters.
#' @param threshold W_SB mass threshold in (0,1).
#' @param risk_level Action threshold on the probability scale.
#' @return Logical policy result, with a typed deprecation warning.
#' @export
check_dominance_risk <- function(a, b, threshold = 0.5, risk_level = 0.3) {
  assert_positive(a, "a")
  assert_positive(b, "b")
  threshold <- .dpprior_validate_probability(
    threshold, "threshold", scalar = TRUE, open = TRUE
  )
  risk_level <- .dpprior_validate_probability(
    risk_level, "risk_level", scalar = TRUE, open = FALSE
  )
  .dpprior_warn(
    paste(
      "check_dominance_risk() is deprecated because its name does not identify",
      "an estimand. This call evaluates P(W_SB > threshold) > risk_level;",
      "it does not evaluate W_max."
    ),
    "dpprior_deprecated_warning",
    "check_dominance_risk", NULL,
    "an explicit DPprior_diagnostics(warning_policy=...) policy",
    "deprecated_ambiguous_name"
  )
  .diagnostic_wsb_tail(threshold, a, b) > risk_level
}

# =============================================================================
# Canonical diagnostic bundle
# =============================================================================

.diagnostic_abort_order <- function(message, fit, code, actual, expected) {
  stop(.dpprior_new_condition(
    message = message,
    classes = c(
      "dpprior_diagnostics_order_error", "dpprior_diagnostics_error",
      "dpprior_invalid_input", "dpprior_error", "error"
    ),
    code = code, argument = "M_verify", value = actual,
    expected = expected, result = fit
  ))
}

.diagnostic_validate_reserved_thresholds <- function(thresholds) {
  validated <- .dpprior_validate_probability(
    thresholds, "thresholds", scalar = FALSE, open = TRUE
  )
  if (!identical(validated, c(0.5, 0.9))) {
    .dpprior_abort_invalid(
      paste(
        "thresholds is reserved for the fixed summary view and must remain",
        "c(0.5, 0.9); use prob_wsb_exceeds() for W_SB tails or",
        "wmax_tail_bounds() for certified W_max bounds."
      ),
      c("dpprior_diagnostics_input_error", "dpprior_control_error"),
      "thresholds", thresholds, "c(0.5, 0.9)", "reserved_thresholds"
    )
  }
  invisible(validated)
}

.diagnostic_canonical_orders <- function(fit, raw, M_verify) {
  fit_orders <- raw[["computation", exact = TRUE]][["orders", exact = TRUE]]
  a1 <- identical(raw[["mode", exact = TRUE]], "a1_proxy")
  M_selected <- if (a1) {
    as.integer(.QUAD_NODES_DEFAULT)
  } else {
    fit_orders[["M_selected", exact = TRUE]]
  }
  selected_ok <- is.integer(M_selected) && !is.object(M_selected) &&
    is.null(dim(M_selected)) && length(M_selected) == 1L &&
    !is.na(M_selected) && M_selected >= 1L &&
    M_selected <= .QUADRATURE_MAX_NODES
  if (!selected_ok) {
    .diagnostic_abort_order(
      "The canonical input fit has no usable selected diagnostic order.",
      fit, "diagnostic_selected_order_unavailable", M_selected,
      sprintf("integer in [1,%d]", .QUADRATURE_MAX_NODES)
    )
  }
  M_required <- as.integer(.quadrature_verification_required_order(M_selected))
  if (M_required > .QUADRATURE_MAX_NODES) {
    .diagnostic_abort_order(
      sprintf(
        paste(
          "Independent diagnostics verification is unavailable for selected",
          "order %d because required order %d exceeds %d."
        ),
        M_selected, M_required, .QUADRATURE_MAX_NODES
      ),
      fit, "diagnostic_verification_order_unavailable", M_required,
      sprintf("required order <= %d", .QUADRATURE_MAX_NODES)
    )
  }

  explicit <- !is.null(M_verify)
  if (explicit) {
    valid_M_verify <- is.numeric(M_verify) && !is.object(M_verify) &&
      is.null(attributes(M_verify)) && length(M_verify) == 1L &&
      !is.na(M_verify) && is.finite(M_verify) &&
      M_verify == floor(M_verify) && M_verify >= 1 &&
      M_verify <= .QUADRATURE_MAX_NODES
    if (!valid_M_verify) {
      .diagnostic_abort_order(
        sprintf(
          "M_verify must be one ordinary integer in [1,%d].",
          .QUADRATURE_MAX_NODES
        ),
        fit, "diagnostic_verification_order_bounds", M_verify,
        sprintf("ordinary integer in [1,%d]", .QUADRATURE_MAX_NODES)
      )
    }
    M_used <- as.integer(M_verify)
    if (M_used < M_required || M_used <= M_selected) {
      .diagnostic_abort_order(
        sprintf(
          "M_verify must be an integer in [%d,%d] and exceed M_selected=%d.",
          M_required, .QUADRATURE_MAX_NODES, M_selected
        ),
        fit, "insufficient_diagnostic_verification_order", M_used,
        sprintf("integer in [%d,%d]", M_required, .QUADRATURE_MAX_NODES)
      )
    }
    used_reason <- "explicit_diagnostics_verifier_order"
  } else {
    retained <- if (a1) NULL else {
      fit_orders[["M_verification_used", exact = TRUE]]
    }
    retained_ok <- is.integer(retained) && !is.object(retained) &&
      is.null(dim(retained)) && length(retained) == 1L &&
      !is.na(retained) && retained >= M_required &&
      retained <= .QUADRATURE_MAX_NODES
    M_used <- if (retained_ok) retained else M_required
    used_reason <- if (retained_ok) {
      "input_fit_verifier_order"
    } else {
      "canonical_required_verifier_order"
    }
  }
  list(
    M_selected = as.integer(M_selected),
    M_required = M_required,
    M_used = as.integer(M_used),
    requested_reason = if (a1) {
      "a1_diagnostics_default_order"
    } else {
      "input_fit_selected_order"
    },
    selected_reason = if (a1) {
      "fresh_diagnostics_default_selected_order"
    } else {
      "fresh_diagnostics_input_fit_selected_order"
    },
    verification_used_reason = used_reason
  )
}

.diagnostic_canonical_truth <- function(
    parameters, J, M_selected, M_verification, absolute, relative) {
  a <- parameters[["a", exact = TRUE]]
  b <- parameters[["b", exact = TRUE]]
  K <- .get_K_pmf_support(
    J, a, b, M = M_selected, M_verify = M_verification,
    abs_tol = absolute, rel_tol = relative
  )
  selected_pmf <- unname(as.numeric(K[["pmf", exact = TRUE]]))
  verifier_pmf <- unname(as.numeric(
    K[["verification_pmf", exact = TRUE]]
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
      abs(sum(selected_pmf) - 1) <= .TOL_PMF_SUM &&
      abs(sum(verifier_pmf) - 1) <= .TOL_PMF_SUM,
    weights = delta[["weights.mean"]] <= refinement[["weights.mean"]],
    coclustering = all(delta[c(
      "coclustering.mean", "coclustering.variance"
    )] <= refinement[c(
      "coclustering.mean", "coclustering.variance"
    )])
  )
  selected <- list(
    alpha = list(mean = a / b, CV = 1 / sqrt(a)),
    K = list(
      mean = unname(selected_K[["mean"]]),
      variance = unname(selected_K[["variance"]]),
      estimand = "K_J", source = "fresh_diagnostics_selected_order",
      M = M_selected, pmf = selected_pmf
    ),
    weights = list(mean = selected_weight),
    coclustering = as.list(selected_rho)
  )
  verifier <- list(
    alpha = list(mean = a / b, CV = 1 / sqrt(a)),
    K = list(
      mean = unname(verifier_K[["mean"]]),
      variance = unname(verifier_K[["variance"]]),
      estimand = "K_J", source = "fresh_diagnostics_verifier_evidence",
      M = M_verification, pmf = verifier_pmf
    ),
    weights = list(mean = verifier_weight),
    coclustering = as.list(verifier_rho)
  )
  list(
    controls = list(
      absolute_tolerance = absolute,
      relative_tolerance = relative,
      pmf_mass_tolerance = .TOL_PMF_SUM
    ),
    tolerances = list(diagnostics = list(
      absolute = absolute, relative = relative,
      pmf_mass = .TOL_PMF_SUM, refinement = refinement
    )),
    residuals = list(diagnostics = delta),
    selected = selected, verifier = verifier,
    component_pass = component_pass
  )
}

.diagnostic_canonical_attempts <- function(component_pass) {
  lapply(seq_along(.DPPRIOR_DIAGNOSTIC_COMPONENTS), function(index) {
    component <- .DPPRIOR_DIAGNOSTIC_COMPONENTS[[index]]
    status <- if (component_pass[[component]]) "converged" else "approximate"
    .dpprior_new_attempt(
      id = paste0("diagnostic-", component),
      stage = "diagnostic_component",
      method = unname(.DPPRIOR_DIAGNOSTIC_ATTEMPT_METHODS[[component]]),
      start = NULL, bounds = NULL,
      control = list(component = component),
      exit_code = 0L,
      message = "Component diagnostic was freshly recomputed.",
      iterations = 0L,
      evaluations = list(function_count = 1L),
      candidate_parameters = NULL, candidate_objective = NULL,
      elapsed_seconds = 0, warnings = character(), error = NULL,
      selected = FALSE, reason_code = paste0("component_", status),
      unavailable = c(
        start = "component diagnostic has no optimizer start",
        bounds = "component diagnostic has no optimizer bounds",
        candidate_parameters = "diagnostics do not select fit parameters",
        candidate_objective = "diagnostics do not optimize an objective"
      )
    )
  })
}

.diagnostic_canonical_computation <- function(truth, orders) {
  setting <- list(
    method = "canonical_prior_diagnostics",
    controls = truth[["controls", exact = TRUE]],
    parameterization = "Gamma(shape=a, rate=b)"
  )
  .dpprior_new_computation(
    request = setting, used = setting,
    orders = .dpprior_new_orders(
      M_requested = orders[["M_selected", exact = TRUE]],
      M_selected = orders[["M_selected", exact = TRUE]],
      M_verification_required = orders[["M_required", exact = TRUE]],
      M_verification_used = orders[["M_used", exact = TRUE]],
      requested_reason = orders[["requested_reason", exact = TRUE]],
      selected_reason = orders[["selected_reason", exact = TRUE]],
      verification_required_reason =
        "canonical_independent_quadrature_contract",
      verification_used_reason = orders[[
        "verification_used_reason", exact = TRUE
      ]]
    ),
    scaling = .dpprior_new_scaling(),
    attempts = .diagnostic_canonical_attempts(
      truth[["component_pass", exact = TRUE]]
    ),
    candidate_evaluations = list(),
    selected_candidate_id = NULL, selected_attempt_id = NULL,
    fallback = .dpprior_new_fallback(),
    termination = .dpprior_new_termination(
      code = "diagnostics_recomputed",
      message = "All diagnostic components were freshly recomputed.",
      source = "component_aggregation", iterations = NULL
    ),
    trace = NULL,
    resources = list(component_elapsed_seconds = stats::setNames(
      rep(0, length(.DPPRIOR_DIAGNOSTIC_COMPONENTS)),
      .DPPRIOR_DIAGNOSTIC_COMPONENTS
    ))
  )
}

.diagnostic_canonical_verification <- function(truth, parameters, orders) {
  selected <- .dpprior_new_snapshot(
    parameters = parameters,
    M = orders[["M_selected", exact = TRUE]],
    achieved = truth[["selected", exact = TRUE]],
    residuals = truth[["residuals", exact = TRUE]],
    tolerances = truth[["tolerances", exact = TRUE]],
    finite = TRUE, source = "fresh_diagnostics_selected_order"
  )
  verifier <- .dpprior_new_snapshot(
    parameters = parameters,
    M = orders[["M_used", exact = TRUE]],
    achieved = truth[["verifier", exact = TRUE]],
    residuals = truth[["residuals", exact = TRUE]],
    tolerances = truth[["tolerances", exact = TRUE]],
    finite = TRUE, source = "fresh_diagnostics_verifier_evidence"
  )
  component_pass <- truth[["component_pass", exact = TRUE]]
  component_reference <- stats::setNames(
    rep(TRUE, length(component_pass)), names(component_pass)
  )
  check <- function(value, reference = TRUE) {
    .dpprior_new_check(
      value = value, reference = reference, tolerance = NULL,
      operator = "identical", source = "fresh_component_specific_checks"
    )
  }
  .dpprior_new_verification(
    method = "fresh_component_specific_diagnostics",
    performed = TRUE, passed = all(component_pass),
    reason = "Fresh component-specific verification completed.",
    settings = list(
      M_selected = orders[["M_selected", exact = TRUE]],
      M_verification = orders[["M_used", exact = TRUE]]
    ),
    selected_snapshot = selected, verifier_snapshot = verifier,
    stability = NULL,
    components = list(component_aggregation = check(
      component_pass, component_reference
    )),
    invariants = list(
      fixed_parameters = check(TRUE),
      dominance_category_removed = check(TRUE)
    )
  )
}

.diagnostic_canonical_provenance <- function(status, allow_approximate) {
  approximate <- identical(status, "approximate")
  source_commit <- getOption("DPprior.source_commit", NULL)
  if (!is.character(source_commit) || is.object(source_commit) ||
      !is.null(dim(source_commit)) || length(source_commit) != 1L ||
      is.na(source_commit) || !nzchar(source_commit)) {
    source_commit <- NULL
  }
  .dpprior_new_provenance(
    requested_method = "canonical_prior_diagnostics",
    selected_method = "canonical_prior_diagnostics",
    is_fallback = FALSE,
    approximation = list(
      active = approximate,
      opt_in = approximate && allow_approximate,
      kind = if (approximate) {
        "diagnostic_refinement_not_verified"
      } else NULL,
      warning_code = if (approximate) "diagnostics_approximate" else NULL
    ),
    projection = list(
      applied = FALSE, opt_in = FALSE, policy = NULL, record = NULL
    ),
    parameterization = "Gamma(shape=a, rate=b)",
    backend = list(
      package = "DPprior",
      package_version = tryCatch(
        as.character(utils::packageVersion("DPprior")),
        error = function(condition) "development"
      ),
      implementation = "R/14_diagnostics.R:DPprior_diagnostics",
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

.diagnostic_new_bundle <- function(
    raw, orders, truth, warning_policy, warning_result,
    allow_approximate) {
  component_pass <- truth[["component_pass", exact = TRUE]]
  component_status <- ifelse(component_pass, "converged", "approximate")
  status <- if (all(component_pass)) "converged" else "approximate"
  usable <- identical(status, "converged")
  verified <- usable && all(component_pass)
  diagnostics <- list(
    policy_results = warning_result[["records", exact = TRUE]],
    warnings = warning_result[["messages", exact = TRUE]],
    alpha = c(
      list(
        status = unname(component_status[["alpha"]]),
        usable = unname(component_pass[["alpha"]]),
        verified = unname(component_pass[["alpha"]])
      ),
      truth[["selected", exact = TRUE]][["alpha", exact = TRUE]]
    ),
    K = c(
      list(
        status = unname(component_status[["K"]]),
        usable = unname(component_pass[["K"]]),
        verified = unname(component_pass[["K"]])
      ),
      truth[["selected", exact = TRUE]][["K", exact = TRUE]][
        c("mean", "variance", "pmf", "M")
      ]
    ),
    weights = c(
      list(
        status = unname(component_status[["weights"]]),
        usable = unname(component_pass[["weights"]]),
        verified = unname(component_pass[["weights"]])
      ),
      truth[["selected", exact = TRUE]][["weights", exact = TRUE]]
    ),
    coclustering = c(
      list(
        status = unname(component_status[["coclustering"]]),
        usable = unname(component_pass[["coclustering"]]),
        verified = unname(component_pass[["coclustering"]])
      ),
      truth[["selected", exact = TRUE]][[
        "coclustering", exact = TRUE
      ]]
    )
  )
  parameters <- .dpprior_new_parameters(
    raw[["parameters", exact = TRUE]][["a", exact = TRUE]],
    raw[["parameters", exact = TRUE]][["b", exact = TRUE]],
    "Gamma(shape=a, rate=b)"
  )
  result <- .dpprior_new_diagnostics(
    method = "canonical_prior_diagnostics",
    J = raw[["J", exact = TRUE]],
    status = status, usable = usable, verified = verified,
    message = if (verified) {
      paste(
        "All canonical diagnostic components passed independent",
        "selected-versus-verifier checks."
      )
    } else {
      paste(
        "At least one canonical diagnostic component did not pass the",
        "fixed refinement contract."
      )
    },
    parameters = parameters,
    target = list(
      requested_components = .DPPRIOR_DIAGNOSTIC_COMPONENTS,
      warning_policy = warning_policy
    ),
    achieved = truth[["selected", exact = TRUE]],
    residuals = truth[["residuals", exact = TRUE]],
    tolerances = truth[["tolerances", exact = TRUE]],
    computation = .diagnostic_canonical_computation(truth, orders),
    verification = .diagnostic_canonical_verification(
      truth, parameters, orders
    ),
    provenance = .diagnostic_canonical_provenance(
      status, allow_approximate
    ),
    diagnostics = diagnostics,
    compatibility = .dpprior_new_compatibility()
  )
  .dpprior_validate_object(result)
  result
}

#' Comprehensive Prior Diagnostics
#'
#' @param fit A canonical \code{dpprior.result/1} fit. Legacy and flat lists
#'   are rejected; call \code{upgrade_DPprior_object()} explicitly first.
#' @param thresholds Reserved fixed summary-view thresholds. The only accepted
#'   value is \code{c(0.5, 0.9)}. Use \code{prob_wsb_exceeds()} for arbitrary
#'   W_SB tails or \code{wmax_tail_bounds()} for certified W_max bounds.
#' @param warning_policy Optional ordinary list with exact fields
#'   \code{estimand}, \code{direction}, \code{weight_threshold}, and
#'   \code{action_threshold}. W_max policies are retained as indeterminate
#'   because canonical diagnostics make no W_max point-probability claim.
#' @param M_verify Optional independent quadrature order.
#' @param abs_tol,rel_tol Comparison tolerances.
#' @param allow_approximate Return an approximate bundle when TRUE; otherwise
#'   signal a typed condition carrying the complete bundle.
#' @return A \code{DPprior_diagnostics} object using the shared result schema.
#' @family diagnostics
#' @export
DPprior_diagnostics <- function(
    fit, thresholds = c(0.5, 0.9), warning_policy = NULL,
    M_verify = NULL, abs_tol = 1e-10, rel_tol = 1e-8,
    allow_approximate = FALSE) {
  fit <- .dpprior_require_schema(fit, kind = "fit", allow_legacy = FALSE)
  raw <- unclass(fit)
  .diagnostic_validate_reserved_thresholds(thresholds)
  allow_approximate <- .dpprior_validate_control(
    allow_approximate, "allow_approximate", type = "logical"
  )
  abs_tol <- .dpprior_validate_scalar(
    abs_tol, "abs_tol", lower = 0,
    .subclass = "dpprior_diagnostics_input_error"
  )
  rel_tol <- .dpprior_validate_scalar(
    rel_tol, "rel_tol", lower = 0,
    .subclass = "dpprior_diagnostics_input_error"
  )
  if (abs_tol > 1e-10 || rel_tol > 1e-8) {
    .dpprior_abort_invalid(
      "Diagnostic tolerances cannot exceed absolute=1e-10 or relative=1e-8.",
      c("dpprior_diagnostics_input_error", "dpprior_control_error"),
      "abs_tol/rel_tol", c(abs_tol = abs_tol, rel_tol = rel_tol),
      "abs_tol <= 1e-10 and rel_tol <= 1e-8",
      "diagnostic_tolerance_bounds"
    )
  }
  warning_policy <- .diagnostic_validate_warning_policy(warning_policy)
  parameters <- raw[["parameters", exact = TRUE]]
  if (is.null(parameters)) {
    .diagnostic_abort_numerical(
      paste(
        "Canonical diagnostics require a public finite parameter candidate;",
        "the validated input fit is retained in condition$result."
      ),
      result = fit, code = "diagnostic_candidate_unavailable",
      subclass = c(
        "dpprior_diagnostics_input_error",
        "dpprior_diagnostics_computation_error"
      )
    )
  }
  orders <- .diagnostic_canonical_orders(fit, raw, M_verify)
  caught_warning <- NULL
  truth <- tryCatch(
    withCallingHandlers(
      .diagnostic_canonical_truth(
        .dpprior_new_parameters(
          parameters[["a", exact = TRUE]],
          parameters[["b", exact = TRUE]],
          "Gamma(shape=a, rate=b)"
        ),
        raw[["J", exact = TRUE]],
        orders[["M_selected", exact = TRUE]],
        orders[["M_used", exact = TRUE]], abs_tol, rel_tol
      ),
      warning = function(warning) {
        caught_warning <<- warning
        invokeRestart("muffleWarning")
      }
    ),
    error = function(error) error
  )
  if (inherits(truth, "condition") || !is.null(caught_warning)) {
    cause <- if (inherits(truth, "condition")) truth else caught_warning
    .diagnostic_abort_numerical(
      paste(
        "Fresh canonical diagnostics recomputation failed; the validated",
        "input fit is retained in condition$result."
      ),
      result = fit, cause = cause, code = "diagnostic_recomputation_failed",
      subclass = "dpprior_diagnostics_computation_error"
    )
  }
  warning_result <- .diagnostic_policy_evidence(
    warning_policy,
    parameters[["a", exact = TRUE]], parameters[["b", exact = TRUE]]
  )
  result <- .diagnostic_new_bundle(
    raw, orders, truth, warning_policy, warning_result,
    allow_approximate
  )
  result_raw <- unclass(result)
  if (identical(result_raw[["status", exact = TRUE]], "approximate") &&
      !allow_approximate) {
    .diagnostic_abort_numerical(
      paste(
        result_raw[["message", exact = TRUE]],
        "Set allow_approximate=TRUE only to return the unusable review object."
      ),
      result = result, code = "approximation_not_accepted",
      subclass = "dpprior_diagnostics_approximation_error"
    )
  }
  messages <- warning_result[["messages", exact = TRUE]]
  if (length(messages)) {
    warning(.dpprior_new_condition(
      message = messages[[1L]],
      classes = c(
        "dpprior_diagnostic_policy_warning", "dpprior_estimand_warning",
        "dpprior_warning", "warning"
      ),
      code = "policy_triggered", argument = "warning_policy",
      value = warning_policy,
      expected = "explicit estimand-specific user policy",
      result = result
    ))
  }
  result
}

# =============================================================================
# S3 methods
# =============================================================================

.dpprior_s3_diagnostics_gate <- function(x) {
  validated <- .dpprior_require_schema(
    x, kind = "diagnostics", allow_legacy = FALSE
  )
  list(
    canonical = TRUE,
    detected = .DPPRIOR_RESULT_SCHEMA_V1,
    raw = unclass(validated)
  )
}


.dpprior_s3_diagnostics_scalar <- function(component, field,
                                            default = NA_real_) {
  if (is.null(component)) {
    return(default)
  }
  value <- component[[field, exact = TRUE]]
  if (is.null(value) || length(value) != 1L) default else value
}


.dpprior_s3_canonical_diagnostics_view <- function(raw) {
  parameters <- raw[["parameters", exact = TRUE]]
  diagnostics <- raw[["diagnostics", exact = TRUE]]
  alpha <- diagnostics[["alpha", exact = TRUE]]
  K <- diagnostics[["K", exact = TRUE]]
  weights <- diagnostics[["weights", exact = TRUE]]
  coclustering <- diagnostics[["coclustering", exact = TRUE]]
  provenance <- raw[["provenance", exact = TRUE]]
  migration <- provenance[["migration", exact = TRUE]]
  migrated <- !is.null(migration) && identical(
    migration[["adapter", exact = TRUE]], "upgrade_DPprior_object"
  )
  a <- parameters[["a", exact = TRUE]]
  b <- parameters[["b", exact = TRUE]]
  pmf <- if (is.null(K)) NULL else K[["pmf", exact = TRUE]]
  K_mode <- if (is.null(pmf)) NA_integer_ else as.integer(which.max(pmf))
  K_median <- if (is.null(pmf)) NA_integer_ else {
    as.integer(which(cumsum(pmf) >= 0.5)[[1L]])
  }

  list(
    schema = .DPPRIOR_RESULT_SCHEMA_V1,
    migrated = migrated,
    status = raw[["status", exact = TRUE]],
    usable = raw[["usable", exact = TRUE]],
    verified = raw[["verified", exact = TRUE]],
    method = raw[["method", exact = TRUE]],
    message = raw[["message", exact = TRUE]],
    J = raw[["J", exact = TRUE]],
    a = a,
    b = b,
    alpha = alpha,
    K = K,
    weights = weights,
    coclustering = coclustering,
    alpha_sd = if (is.null(alpha)) NA_real_ else sqrt(a) / b,
    alpha_median = if (is.null(alpha)) NA_real_ else {
      stats::qgamma(0.5, shape = a, rate = b)
    },
    K_sd = if (is.null(K)) NA_real_ else {
      sqrt(K[["variance", exact = TRUE]])
    },
    K_mode = K_mode,
    K_median = K_median,
    wsb_median = if (is.null(weights)) NA_real_ else quantile_w1(0.5, a, b),
    wsb_tail_50 = if (is.null(weights)) NA_real_ else {
      .diagnostic_wsb_tail(0.5, a, b)
    },
    wsb_tail_90 = if (is.null(weights)) NA_real_ else {
      .diagnostic_wsb_tail(0.9, a, b)
    },
    rho_sd = if (is.null(coclustering)) NA_real_ else {
      sqrt(coclustering[["variance", exact = TRUE]])
    },
    policy_results = diagnostics[["policy_results", exact = TRUE]],
    warnings = diagnostics[["warnings", exact = TRUE]]
  )
}


#' Print Method for DPprior_diagnostics Objects
#'
#' @param x A \code{DPprior_diagnostics} object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.DPprior_diagnostics <- function(x, ...) {
  original <- x
  gate <- .dpprior_s3_diagnostics_gate(x)
  if (isTRUE(gate[["canonical", exact = TRUE]])) {
    view <- .dpprior_s3_canonical_diagnostics_view(
      gate[["raw", exact = TRUE]]
    )
    cat("DPprior Prior Diagnostics\n")
    cat(strrep("=", 60), "\n")
    cat(sprintf("Schema: %s\n", view[["schema", exact = TRUE]]))
    cat(sprintf("Status: %s\n", toupper(view[["status", exact = TRUE]])))
    cat(sprintf("Usable: %s; verified: %s\n",
                if (isTRUE(view[["usable", exact = TRUE]])) "yes" else "no",
                if (isTRUE(view[["verified", exact = TRUE]])) "yes" else "no"))
    cat(sprintf("Method used: %s\n", view[["method", exact = TRUE]]))
    cat(sprintf("Message: %s\n\n", view[["message", exact = TRUE]]))
    cat(sprintf(
      "Prior: alpha ~ Gamma(shape=%.6g, rate=%.6g); J=%d\n\n",
      view[["a", exact = TRUE]], view[["b", exact = TRUE]],
      view[["J", exact = TRUE]]
    ))

    alpha <- view[["alpha", exact = TRUE]]
    cat("alpha (DP concentration parameter)\n")
    if (is.null(alpha)) {
      cat("  Unavailable: canonical alpha component was not produced.\n\n")
    } else {
      cat(sprintf(
        "  Status %s; mean %.4g; SD %.4g; CV %.4g; median %.4g\n\n",
        alpha[["status", exact = TRUE]], alpha[["mean", exact = TRUE]],
        view[["alpha_sd", exact = TRUE]], alpha[["CV", exact = TRUE]],
        view[["alpha_median", exact = TRUE]]
      ))
    }

    K <- view[["K", exact = TRUE]]
    cat(sprintf("K_J (occupied clusters among J=%d units)\n",
                view[["J", exact = TRUE]]))
    if (is.null(K)) {
      cat("  Unavailable: canonical K_J component was not produced.\n\n")
    } else {
      cat(sprintf(
        "  Status %s; E[K_J] %.4g; SD %.4g; median %d; mode %d; M=%d\n\n",
        K[["status", exact = TRUE]], K[["mean", exact = TRUE]],
        view[["K_sd", exact = TRUE]], view[["K_median", exact = TRUE]],
        view[["K_mode", exact = TRUE]], K[["M", exact = TRUE]]
      ))
    }

    weights <- view[["weights", exact = TRUE]]
    if (is.null(weights)) {
      cat("W_SB: unavailable; canonical weight component was not produced.\n\n")
    } else {
      cat("W_SB (first size-biased DP weight)\n")
      cat(sprintf("  Status %s; mean %.4g; median %.4g\n",
                  weights[["status", exact = TRUE]],
                  weights[["mean", exact = TRUE]],
                  view[["wsb_median", exact = TRUE]]))
      cat(sprintf("  P(W_SB > 0.5) = %.6g\n",
                  view[["wsb_tail_50", exact = TRUE]]))
      cat(sprintf("  P(W_SB > 0.9) = %.6g\n",
                  view[["wsb_tail_90", exact = TRUE]]))
      cat("  W_SB is not W_max.\n\n")
    }
    if (isTRUE(view[["migrated", exact = TRUE]])) {
      cat(paste0(
        "W_max: unavailable (migration freshly reconstructed W_SB only; no ",
        "W_max estimate was substituted).\n\n"
      ))
    } else {
      cat(paste0(
        "W_max: unavailable (this canonical diagnostics bundle contains no ",
        "W_max estimate; no W_SB value was substituted).\n\n"
      ))
    }

    coclustering <- view[["coclustering", exact = TRUE]]
    cat("rho (conditional pairwise co-clustering probability)\n")
    if (is.null(coclustering)) {
      cat("  Unavailable: canonical co-clustering component was not produced.\n\n")
    } else {
      cat(sprintf("  Status %s; mean %.4g; SD %.4g\n\n",
                  coclustering[["status", exact = TRUE]],
                  coclustering[["mean", exact = TRUE]],
                  view[["rho_sd", exact = TRUE]]))
    }
    policy_results <- view[["policy_results", exact = TRUE]]
    if (length(policy_results)) {
      cat(sprintf("Explicit warning policy: %d recorded result(s).\n",
                  length(policy_results)))
    } else {
      cat("No warning policy requested; no categorical warning was computed.\n")
    }
    return(invisible(original))
  }
  stop("unreachable canonical diagnostics print state", call. = FALSE)
}

#' Summary Method for DPprior_diagnostics Objects
#'
#' @param object A \code{DPprior_diagnostics} object.
#' @param ... Additional arguments (ignored).
#' @return One-row data frame of labelled numerical metrics.
#' @export
summary.DPprior_diagnostics <- function(object, ...) {
  gate <- .dpprior_s3_diagnostics_gate(object)
  if (isTRUE(gate[["canonical", exact = TRUE]])) {
    view <- .dpprior_s3_canonical_diagnostics_view(
      gate[["raw", exact = TRUE]]
    )
    alpha <- view[["alpha", exact = TRUE]]
    K <- view[["K", exact = TRUE]]
    weights <- view[["weights", exact = TRUE]]
    coclustering <- view[["coclustering", exact = TRUE]]
    return(data.frame(
      schema = view[["schema", exact = TRUE]],
      migrated = view[["migrated", exact = TRUE]],
      status = view[["status", exact = TRUE]],
      usable = view[["usable", exact = TRUE]],
      verified = view[["verified", exact = TRUE]],
      message = view[["message", exact = TRUE]],
      method = view[["method", exact = TRUE]],
      J = view[["J", exact = TRUE]],
      a = view[["a", exact = TRUE]],
      b = view[["b", exact = TRUE]],
      alpha_available = !is.null(alpha),
      E_alpha = .dpprior_s3_diagnostics_scalar(alpha, "mean"),
      CV_alpha = .dpprior_s3_diagnostics_scalar(alpha, "CV"),
      K_available = !is.null(K),
      E_K_J = .dpprior_s3_diagnostics_scalar(K, "mean"),
      SD_K_J = view[["K_sd", exact = TRUE]],
      Mode_K_J = view[["K_mode", exact = TRUE]],
      weights_available = !is.null(weights),
      E_W_SB = .dpprior_s3_diagnostics_scalar(weights, "mean"),
      P_W_SB_gt_50 = view[["wsb_tail_50", exact = TRUE]],
      P_W_SB_gt_90 = view[["wsb_tail_90", exact = TRUE]],
      W_max_available = FALSE,
      P_W_max_gt_50 = NA_real_,
      P_W_max_gt_90 = NA_real_,
      W_max_guidance = if (isTRUE(view[["migrated", exact = TRUE]])) {
        paste(
          "W_max was not reconstructed during migration; no W_SB value was",
          "substituted. Use wmax_tail_bounds() for certified bounds."
        )
      } else {
        paste(
          "This canonical bundle contains no W_max estimate; no W_SB value",
          "was substituted. Use wmax_tail_bounds() for certified bounds."
        )
      },
      coclustering_available = !is.null(coclustering),
      E_rho = .dpprior_s3_diagnostics_scalar(coclustering, "mean"),
      n_policy_warnings = length(view[["warnings", exact = TRUE]]),
      row.names = NULL,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ))
  }
  stop("unreachable canonical diagnostics summary state", call. = FALSE)
}

#' Data-Frame Method for DPprior_diagnostics Objects
#'
#' @param x A canonical \code{DPprior_diagnostics} object.
#' @param row.names Optional row names; only \code{NULL} is supported.
#' @param optional Ignored compatibility argument from
#'   \code{as.data.frame()}.
#' @param ... Additional arguments passed to the canonical summary method.
#' @return The one-row canonical diagnostics summary data frame.
#' @export
as.data.frame.DPprior_diagnostics <- function(
    x, row.names = NULL, optional = FALSE, ...) {
  if (!is.null(row.names)) {
    .dpprior_abort_invalid(
      "row.names must be NULL for a one-row canonical diagnostics summary.",
      c("dpprior_diagnostics_input_error", "dpprior_control_error"),
      "row.names", row.names, "NULL", "unsupported_row_names"
    )
  }
  .dpprior_validate_control(
    optional, "optional", type = "logical"
  )
  summary.DPprior_diagnostics(x, ...)
}

#' Compare Diagnostics Across Multiple Fits
#'
#' @param ... Named DPprior fit objects.
#' @param M Reserved historical argument. It must remain at the canonical
#'   default because each fit's retained selected order is authoritative.
#' @return Data frame of estimand-labelled diagnostics.
#' @keywords internal
compare_diagnostics <- function(..., M = .QUAD_NODES_DEFAULT) {
  M <- .dpprior_validate_count(
    M, "M", minimum = 1L, maximum = .QUADRATURE_MAX_NODES,
    .subclass = "dpprior_diagnostics_input_error"
  )
  if (!identical(M, as.integer(.QUAD_NODES_DEFAULT))) {
    .dpprior_abort_invalid(
      paste(
        "compare_diagnostics() uses each canonical fit's retained order;",
        "the historical M argument must remain at its default."
      ),
      c("dpprior_diagnostics_input_error", "dpprior_control_error"),
      "M", M, as.integer(.QUAD_NODES_DEFAULT), "reserved_order"
    )
  }
  fits <- list(...)
  if (!length(fits)) {
    .dpprior_abort_invalid(
      "At least one fit object is required",
      c("dpprior_diagnostics_input_error", "dpprior_length_error"),
      "...", fits, "one or more fits", "length"
    )
  }
  fit_names <- names(fits)
  if (is.null(fit_names) || any(!nzchar(fit_names))) {
    fit_names <- paste0("Fit_", seq_along(fits))
  }
  results <- lapply(seq_along(fits), function(i) {
    diagnostic <- DPprior_diagnostics(
      fits[[i]], thresholds = c(0.5, 0.9)
    )
    output <- summary(diagnostic)
    output$fit <- fit_names[[i]]
    output[, c("fit", setdiff(names(output), "fit")), drop = FALSE]
  })
  do.call(rbind, results)
}

#' Verify Diagnostics Module
#'
#' @param verbose Logical; print progress when TRUE.
#' @return Invisibly returns TRUE.
#' @keywords internal
verify_diagnostics <- function(verbose = TRUE) {
  verbose <- .dpprior_validate_control(verbose, "verbose", type = "logical")
  cases <- list(
    c(a = 0.5, b = 0.5), c(a = 1, b = 1), c(a = 2, b = 1),
    c(a = 1.6, b = 1.22), c(a = 5, b = 2)
  )
  for (parameters in cases) {
    stopifnot(abs(
      mean_w1(parameters[["a"]], parameters[["b"]]) -
        mean_rho(parameters[["a"]], parameters[["b"]])
    ) < 1e-6)
  }
  weight <- compute_weight_diagnostics(
    1.6, 1.22, thresholds = c(0.5, 0.9)
  )
  stopifnot(is.null(weight$dominance_risk))
  stopifnot(weight$size_biased$estimand == "W_SB")
  stopifnot(weight$size_biased$tail_probability[[1L]] > 0.4)
  stopifnot(!isTRUE(weight$maximum$requested))
  stopifnot(is.null(weight$maximum$values))
  stopifnot(abs(.diagnostic_wsb_tail(0.5, 1.6, 1.22) - 0.4868311) < 1e-5)
  fit <- DPprior_fit(
    20L, 4, 8, method = "A1", check_diagnostics = FALSE
  )
  diagnostic <- DPprior_diagnostics(fit)
  stopifnot(inherits(diagnostic, "dpprior_result"))
  .dpprior_validate_result_v1(diagnostic)
  if (isTRUE(verbose)) cat("Diagnostics verification passed.\n")
  invisible(TRUE)
}
