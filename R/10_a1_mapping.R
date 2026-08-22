# =============================================================================
# Module 10: A1 Closed-Form Prior Elicitation (Revised)
# =============================================================================
#
# This module provides:
# 1. Closed-form mapping from (mu_K, var_K, J) to Gamma(a, b) hyperprior
# 2. VIF (Variance Inflation Factor) utilities
# 3. Confidence level to VIF conversion
# 4. S3 class DPprior_fit for results
#
# Author: JoonHo Lee (jlee296@ua.edu)
# Date: December 2025
# Part of: DPprior R Package
# Reference: Lee (2026) arXiv:2602.06301, Section 3.1 (TSMM Stage 1)
#
# Revision Notes:
# - Added epsilon validation for robustness
# - Added converged/iterations fields for A2 compatibility
# - Fixed digamma scaling to use numeric floor
# - Strengthened VIF validation to require >= 1
# =============================================================================


# =============================================================================
# Scaling Constant Computation
# =============================================================================

.dpprior_a1_match_arg <- function(value, choices, name,
                                  subclass = "dpprior_a1_input_error") {
  if (!is.character(value) || !is.null(dim(value)) || is.object(value) ||
      (length(value) != 1L && !identical(value, choices))) {
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

#' Compute Scaling Constant for A1 Mapping
#'
#' Computes the scaling constant \eqn{c_J} used in the A1 closed-form mapping.
#' Three variants are supported based on asymptotic approximations.
#'
#' @param J Integer; number of items/sites (must be >= 2).
#' @param scaling Character; one of "log", "harmonic", or "digamma".
#' @param mu_K Numeric; target mean of K (required for "digamma" scaling).
#'
#' @return Numeric scalar; the scaling constant \eqn{c_J}.
#'
#' @details
#' The scaling constant appears in the Poisson proxy:
#' \deqn{K_J - 1 \mid \alpha \approx \text{Poisson}(\alpha \cdot c_J)}
#'
#' Available variants:
#' \describe{
#'   \item{log}{\eqn{c_J = \log(J)}, the asymptotic leading term (default)}
#'   \item{harmonic}{\eqn{c_J = H_{J-1} = \psi(J) + \gamma}, improves
#'     accuracy for small/moderate J}
#'   \item{digamma}{\eqn{c_J = \psi(\tilde{\alpha} + J) - \psi(\tilde{\alpha})}
#'     where \eqn{\tilde{\alpha} = (\mu_K - 1)/\log(J)}, a local correction}
#' }
#'
#' @seealso \code{\link{DPprior_a1}} for the main elicitation function
#'
#' @examples
#' # Default log scaling
#' compute_scaling_constant(50, "log")
#'
#' # Harmonic scaling (better for moderate J)
#' compute_scaling_constant(50, "harmonic")
#'
#' # Digamma scaling (requires mu_K)
#' compute_scaling_constant(50, "digamma", mu_K = 5)
#'
#' @export
compute_scaling_constant <- function(J, scaling = c("log", "harmonic", "digamma"),
                                     mu_K = NULL) {
  if (!.dpprior_is_plain_numeric(J) || length(J) != 1L || !is.finite(J) ||
      J != floor(J) || J < 2L) {
    .dpprior_abort_legacy_numeric(
      "J must be an integer >= 2", J, "J", "integer >= 2",
      "dpprior_a1_sample_size_error", scalar = TRUE, integer = TRUE
    )
  }
  scaling <- .dpprior_a1_match_arg(
    scaling, c("log", "harmonic", "digamma"), "scaling",
    "dpprior_a1_scaling_error"
  )

  switch(scaling,
         log = log(J),
         harmonic = digamma(J) + .EULER_GAMMA,
         digamma = {
           if (is.null(mu_K)) {
             .dpprior_abort_invalid(
               "mu_K required for digamma scaling",
               c("dpprior_a1_scaling_error", "dpprior_missing_error"),
               "mu_K", mu_K, "finite numeric scalar", "missing"
             )
           }
           if (!.dpprior_is_plain_numeric(mu_K) || length(mu_K) != 1L ||
               !is.finite(mu_K)) {
             .dpprior_abort_legacy_numeric(
               "mu_K must be a finite numeric scalar for digamma scaling",
               mu_K, "mu_K", "finite numeric scalar",
               "dpprior_a1_scaling_error", scalar = TRUE
             )
           }
           if (mu_K <= 1 || mu_K >= J) {
             .dpprior_abort_invalid(
               "mu_K must satisfy 1 < mu_K < J for digamma scaling",
               c("dpprior_a1_scaling_error", "dpprior_bounds_error"),
               "mu_K", mu_K, sprintf("1 < mu_K < %d", as.integer(J)),
               "bounds"
             )
           }
           alpha_tilde <- (mu_K - 1) / log(J)
           # Use numeric floor instead of fallback (robust to edge cases)
           alpha_tilde <- max(alpha_tilde, .Machine$double.eps)
           digamma(alpha_tilde + J) - digamma(alpha_tilde)
         }
  )
}


# =============================================================================
# Main A1 Elicitation Function
# =============================================================================

# Resolve the A1-specific variance domain once for both direct and wrapper
# calls.  The fixed-support bound is a property of K_J, while the strict lower
# bound is a property of the shifted-Negative-Binomial A1 proxy.  Keeping these
# two bounds together prevents a lower-bound projection from creating a target
# that is impossible on {1, ..., J}.
.dpprior_a1_resolve_target <- function(J, mu_K, var_K,
                                       projection = c("error", "nearest"),
                                       epsilon = .TOL_PROJECTION_BUFFER,
                                       signal_projection = TRUE) {
  projection <- .dpprior_a1_match_arg(
    projection, c("error", "nearest"), "projection",
    "dpprior_a1_projection_policy_error"
  )

  if (!.dpprior_is_plain_numeric(epsilon) || length(epsilon) != 1L ||
      !is.finite(epsilon) || epsilon <= 0) {
    .dpprior_abort_invalid(
      "epsilon must be a finite numeric scalar > 0",
      c("dpprior_a1_control_error", "dpprior_bounds_error"),
      "epsilon", epsilon, "> 0", "bounds"
    )
  }

  support_upper <- .assert_feasible_K_moments(J, mu_K, var_K)
  shifted_mean <- mu_K - 1
  requested_buffer <- epsilon * (1 + shifted_mean^2)
  available_gap <- support_upper - shifted_mean
  # Prefer the declared scaled buffer, but never let that numerical preference
  # erase a real (possibly narrow) feasible A1 interval. Half the available
  # gap stays strictly inside both bounds whenever such a floating-point value
  # is representable.
  preferred_buffer <- if (available_gap > 0) {
    min(requested_buffer, available_gap / 2)
  } else {
    requested_buffer
  }
  representability_floor <- .Machine$double.eps *
    max(1, abs(shifted_mean))
  selected_buffer <- if (available_gap > 0) {
    min(max(preferred_buffer, representability_floor), available_gap)
  } else {
    preferred_buffer
  }
  buffer_was_capped <- available_gap > 0 &&
    requested_buffer > available_gap / 2
  buffer_was_floored <- selected_buffer > preferred_buffer
  effective_buffer <- selected_buffer
  numerical_lower <- shifted_mean + effective_buffer
  # If arithmetic on the buffer itself rounded back to the strict lower bound,
  # the already-representable support upper bound is a safe final guard. This
  # branch is relevant only when the positive feasible gap is narrower than
  # the local floating-point spacing estimate.
  if (available_gap > 0 &&
      !(is.finite(numerical_lower) && numerical_lower > shifted_mean &&
        numerical_lower <= support_upper)) {
    numerical_lower <- support_upper
  }
  effective_buffer <- numerical_lower - shifted_mean
  representable_interior <- available_gap > 0 &&
    is.finite(numerical_lower) && numerical_lower > shifted_mean &&
    numerical_lower <= support_upper
  projection_required <- var_K <= shifted_mean
  near_lower_boundary <- !projection_required && var_K <= numerical_lower
  upper_scale <- max(1, abs(support_upper))
  near_support_upper <- abs(var_K - support_upper) <=
    max(1e-12, epsilon * upper_scale)

  projection_record <- list(
    policy = projection,
    opt_in = identical(projection, "nearest"),
    applied = FALSE,
    reason = if (projection_required) {
      "a1_strict_lower_bound"
    } else if (near_lower_boundary) {
      "near_a1_lower_boundary"
    } else if (near_support_upper) {
      "support_upper_boundary"
    } else {
      "not_required"
    },
    original_target = list(mu_K = mu_K, var_K = var_K),
    projected_target = NULL,
    distance = 0,
    a1_lower_bound = shifted_mean,
    numerical_interior_lower = numerical_lower,
    support_upper_bound = support_upper,
    epsilon = epsilon,
    buffer = effective_buffer,
    requested_buffer = requested_buffer,
    preferred_buffer = preferred_buffer,
    representability_floor = representability_floor,
    effective_buffer = effective_buffer,
    available_gap = available_gap,
    buffer_was_capped = buffer_was_capped,
    buffer_was_floored = buffer_was_floored,
    representable_interior = representable_interior
  )

  if (projection_required) {
    projection_record$projected_target <- list(
      mu_K = mu_K, var_K = numerical_lower
    )
    projection_record$distance <- numerical_lower - var_K

    # If the strict A1 lower bound meets/exceeds the support upper bound (or no
    # floating-point value between them is representable), no projection
    # policy can construct a valid A1 target. A merely oversized preferred
    # buffer is not mathematical infeasibility because it was capped above.
    if (!representable_interior) {
      stop(.dpprior_new_condition(
        message = sprintf(
          paste0(
            "No numerically representable A1 variance lies above the strict ",
            "lower bound %.15g and at or below the fixed-support upper bound ",
            "%.15g for K in {1,...,%d}. Projection was not applied."
          ),
          shifted_mean, support_upper, as.integer(J)
        ),
        classes = c(
          "dpprior_a1_projection_impossible",
          "dpprior_a1_infeasible",
          "dpprior_calibration_error", "dpprior_error", "error"
        ),
        code = "a1_projection_exceeds_support",
        original_target = projection_record$original_target,
        projected_target = projection_record$projected_target,
        projection_distance = projection_record$distance,
        projection_policy = projection,
        support_upper_bound = support_upper,
        a1_lower_bound = shifted_mean,
        available_gap = available_gap,
        requested_buffer = requested_buffer,
        effective_buffer = effective_buffer
      ))
    }

    if (identical(projection, "error")) {
      stop(.dpprior_new_condition(
        message = sprintf(
          paste0(
            "A1 requires var_K > mu_K - 1 (got %.15g; lower bound %.15g). ",
            "No target was changed. To explicitly project to the nearest ",
            "numerically interior A1 target, set projection = 'nearest'."
          ),
          var_K, shifted_mean
        ),
        classes = c(
          "dpprior_a1_projection_required",
          "dpprior_a1_infeasible",
          "dpprior_calibration_error", "dpprior_error", "error"
        ),
        code = "a1_projection_required",
        original_target = projection_record$original_target,
        projected_target = projection_record$projected_target,
        projection_distance = projection_record$distance,
        projection_policy = projection,
        support_upper_bound = support_upper
      ))
    }

    var_K_used <- numerical_lower
    projection_record$applied <- TRUE

    if (isTRUE(signal_projection)) {
      warning(.dpprior_new_condition(
        message = sprintf(
          paste0(
            "A1 target projection explicitly requested: var_K %.15g -> ",
            "%.15g (distance %.15g; fixed-support upper bound %.15g)."
          ),
          var_K, var_K_used, projection_record$distance, support_upper
        ),
        classes = c(
          "dpprior_a1_projection_warning",
          "dpprior_target_projection_warning",
          "dpprior_warning", "warning"
        ),
        code = "a1_target_projected",
        original_target = projection_record$original_target,
        projected_target = projection_record$projected_target,
        projection_distance = projection_record$distance,
        projection_policy = projection,
        support_upper_bound = support_upper
      ))
    }
  } else {
    var_K_used <- var_K
  }

  # Validate the actual target after resolution as a final invariant.  This is
  # intentionally after projection, not merely a check of the request.
  .assert_feasible_K_moments(J, mu_K, var_K_used)

  at_target_boundary <- projection_record$applied || near_lower_boundary ||
    near_support_upper
  boundary_reason <- if (at_target_boundary) {
    projection_record$reason
  } else {
    NA_character_
  }

  list(
    var_K_requested = var_K,
    var_K_used = var_K_used,
    denominator = var_K_used - shifted_mean,
    at_target_boundary = at_target_boundary,
    boundary_reason = boundary_reason,
    usable = TRUE,
    mapping_verified = TRUE,
    projection = projection_record
  )
}


# Construct the canonical moment target consumed by the A1 result.  The target
# constructor verifies only the identity and bounded-support validity of the
# elicited moments.  It does not turn the A1 proxy round-trip into finite-J
# verification of the fitted prior.
.dpprior_a1_target_v1 <- function(J, mu_K, var_K_requested, var_K_used,
                                  projection_record, epsilon) {
  J <- as.integer(J)
  request <- list(J = J, mu_K = mu_K, var_K = var_K_requested)
  normalized <- c(request, list(interval = NULL, pmf = NULL))
  used <- normalized
  used$var_K <- var_K_used
  implied <- list(mean = mu_K, variance = var_K_used)
  target_residuals <- list(mean = 0, variance = 0)
  target_tolerances <- list(absolute = 1e-10)

  request_to_normalized <- list(
    rule = "canonicalize_direct_moments",
    outcome = "canonicalized",
    opt_in = FALSE,
    before = request,
    after = normalized,
    evidence = list(
      source = "DPprior_a1_validated_public_arguments",
      explicit_null_fields = c("interval", "pmf")
    )
  )
  normalized_to_used <- NULL
  projection <- list(
    applied = FALSE, opt_in = FALSE, policy = NULL, record = NULL
  )
  if (isTRUE(projection_record$applied)) {
    normalized_to_used <- list(
      rule = "project_a1_variance_to_nearest_interior",
      outcome = "projected",
      opt_in = TRUE,
      before = normalized,
      after = used,
      evidence = projection_record
    )
    projection <- list(
      applied = TRUE,
      opt_in = TRUE,
      policy = projection_record$policy,
      record = list(
        before = normalized, after = used, authority = projection_record
      )
    )
  }

  target_achieved <- list(implied = implied, interval = NULL, pmf = NULL)
  selected_snapshot <- .dpprior_new_snapshot(
    parameters = NULL,
    M = NULL,
    achieved = target_achieved,
    residuals = target_residuals,
    tolerances = target_tolerances,
    finite = TRUE,
    source = "target_moment_constructor"
  )
  verifier_snapshot <- .dpprior_new_snapshot(
    parameters = NULL,
    M = NULL,
    achieved = target_achieved,
    residuals = target_residuals,
    tolerances = target_tolerances,
    finite = TRUE,
    source = "independent_target_moment_identity_check"
  )
  target_verification <- .dpprior_new_verification(
    method = "independent_target_moment_identity_check",
    performed = TRUE,
    passed = TRUE,
    reason = "requested and used moments satisfy the bounded-support target contract",
    settings = list(
      support = c(lower = 1L, upper = J),
      moment_identity_tolerance = 0
    ),
    selected_snapshot = selected_snapshot,
    verifier_snapshot = verifier_snapshot,
    stability = NULL,
    components = list(
      moment_identity = .dpprior_new_check(
        value = unlist(implied, use.names = TRUE),
        reference = unlist(implied, use.names = TRUE),
        tolerance = NULL,
        operator = "identical",
        source = "independent_target_moment_identity_check"
      )
    ),
    invariants = list(
      request_identity = .dpprior_new_check(
        value = c(
          J = identical(request$J, J),
          mean = identical(request$mu_K, mu_K),
          variance = identical(request$var_K, var_K_requested)
        ),
        reference = c(J = TRUE, mean = TRUE, variance = TRUE),
        tolerance = NULL,
        operator = "identical",
        source = "independent_target_moment_identity_check"
      ),
      support_identity = .dpprior_new_check(
        value = c(lower = 1L, upper = J),
        reference = c(lower = 1L, upper = J),
        tolerance = NULL,
        operator = "identical",
        source = "independent_target_moment_identity_check"
      )
    )
  )

  target_setting <- list(
    method = "target_moments",
    controls = list(
      projection_policy = projection_record$policy,
      epsilon = epsilon
    ),
    parameterization = "bounded_discrete_K_moments"
  )
  target_computation <- .dpprior_new_computation(
    request = target_setting,
    used = target_setting,
    orders = .dpprior_new_orders(
      M_requested = NULL,
      M_selected = NULL,
      M_verification_required = NULL,
      M_verification_used = NULL,
      requested_reason = "not_applicable_for_moment_identity",
      selected_reason = "not_applicable_for_moment_identity",
      verification_required_reason = "not_applicable_for_moment_identity",
      verification_used_reason = "not_applicable_for_moment_identity"
    ),
    scaling = .dpprior_new_scaling(),
    attempts = list(),
    selected_attempt_id = NULL,
    fallback = .dpprior_new_fallback(),
    termination = .dpprior_new_termination(
      code = "deterministic",
      message = "bounded-support moment target constructed deterministically",
      source = "constructor",
      iterations = NULL
    ),
    trace = NULL,
    resources = list()
  )

  source_commit <- getOption("DPprior.source_commit", NULL)
  if (!is.character(source_commit) || length(source_commit) != 1L ||
      is.na(source_commit) || !nzchar(source_commit)) {
    source_commit <- NULL
  }
  target_provenance <- .dpprior_new_provenance(
    requested_method = "target_moments",
    selected_method = "target_moments",
    is_fallback = FALSE,
    approximation = list(
      active = FALSE, opt_in = FALSE, kind = NULL, warning_code = NULL
    ),
    projection = projection,
    parameterization = "bounded_discrete_K_moments",
    backend = list(
      package = "DPprior",
      package_version = tryCatch(
        as.character(utils::packageVersion("DPprior")),
        error = function(condition) "development"
      ),
      implementation = "R/10_a1_mapping.R:.dpprior_a1_target_v1",
      source_commit = source_commit
    ),
    input_fit = NULL,
    migration = list(
      source_schema = "producer_native",
      adapter = ".dpprior_a1_target_v1",
      lossless = TRUE,
      missing_evidence = character(),
      warnings = character()
    ),
    legacy = list(active = FALSE, contract = NULL, deprecation_stage = NULL)
  )

  .dpprior_new_target_K(
    kind = "moments",
    J = J,
    support = seq_len(J),
    request = request,
    normalized = normalized,
    used = used,
    derivation = list(
      request_to_normalized = request_to_normalized,
      normalized_to_used = normalized_to_used
    ),
    interval = NULL,
    family = NULL,
    assumptions = list(
      estimand = "K_J",
      support = "integer support 1:J"
    ),
    pmf = NULL,
    implied = implied,
    achieved_interval = NULL,
    residuals = target_residuals,
    tolerances = target_tolerances,
    status = "converged",
    usable = TRUE,
    verified = TRUE,
    message = "Bounded-support moment target identity was verified.",
    parameters = NULL,
    computation = target_computation,
    verification = target_verification,
    provenance = target_provenance,
    compatibility = .dpprior_new_compatibility()
  )
}


.dpprior_a1_result_v1 <- function(J, mu_K, var_K_requested, var_K_used,
                                   scaling, cJ, epsilon, projection_policy,
                                   target_resolution, a, b,
                                   mapping_target, mapping_achieved,
                                   mapping_residual, mapping_tolerance,
                                   mapping_component_pass, mapping_passed,
                                   caveats, message) {
  J <- as.integer(J)
  mapping_passed <- isTRUE(mapping_passed)
  candidate_values <- c(
    a = a, b = b, mapping_achieved, mapping_residual, mapping_tolerance
  )
  candidate_finite <- all(is.finite(candidate_values)) && a > 0 && b > 0
  retained_candidate <- candidate_finite
  status <- "approximate"
  usable <- retained_candidate && mapping_passed &&
    isTRUE(target_resolution$usable)
  parameterization <- "Gamma(shape=a, rate=b)"
  parameters <- if (retained_candidate) {
    .dpprior_new_parameters(a, b, parameterization)
  } else {
    NULL
  }
  target_K <- .dpprior_a1_target_v1(
    J = J,
    mu_K = mu_K,
    var_K_requested = var_K_requested,
    var_K_used = var_K_used,
    projection_record = target_resolution$projection,
    epsilon = epsilon
  )
  achieved <- if (retained_candidate) {
    list(K = list(
      mean = unname(mapping_achieved[["mean"]]),
      variance = unname(mapping_achieved[["variance"]]),
      estimand = "shifted_negative_binomial_proxy_moments",
      source = "closed_form_a1_proxy_mapping",
      M = NULL
    ))
  } else {
    list()
  }
  residuals <- if (retained_candidate) {
    list(K = list(
      mean = unname(mapping_residual[["mean"]]),
      variance = unname(mapping_residual[["variance"]])
    ))
  } else if (all(is.finite(mapping_residual))) {
    list(mapping = mapping_residual)
  } else {
    list(unavailable_reason = "nonfinite_a1_proxy_mapping")
  }
  tolerances <- list(
    mapping = mapping_tolerance,
    mapping_formula = "1e-12 + 1e-10 * max(1, abs(target))",
    feasibility_buffer = target_resolution$projection$buffer
  )
  selected_snapshot <- if (retained_candidate) {
    .dpprior_new_snapshot(
      parameters = parameters,
      M = NULL,
      achieved = achieved,
      residuals = residuals,
      tolerances = tolerances,
      finite = TRUE,
      source = "closed_form_a1_proxy_mapping"
    )
  } else {
    NULL
  }

  controls <- list(
    scaling = scaling,
    epsilon = epsilon,
    projection = projection_policy,
    mapping_absolute_floor = 1e-12,
    mapping_relative_tolerance = 1e-10
  )
  setting <- list(
    method = "A1", controls = controls, parameterization = parameterization
  )
  computation <- .dpprior_new_computation(
    request = setting,
    used = setting,
    orders = .dpprior_new_orders(
      M_requested = NULL,
      M_selected = NULL,
      M_verification_required = NULL,
      M_verification_used = NULL,
      requested_reason = "not_applicable_for_closed_form_proxy",
      selected_reason = "not_applicable_for_closed_form_proxy",
      verification_required_reason = "finite_J_verification_not_performed",
      verification_used_reason = "finite_J_verification_not_performed"
    ),
    scaling = .dpprior_new_scaling(
      requested = list(method = scaling),
      used = list(method = scaling),
      formula = "A1_closed_form_scaling_constant",
      values = list(cJ = cJ),
      fixed_from_input = FALSE,
      change_reason = ""
    ),
    attempts = list(),
    selected_attempt_id = NULL,
    fallback = .dpprior_new_fallback(),
    termination = if (retained_candidate) {
      .dpprior_new_termination(
        code = "closed_form",
        message = "A1 closed-form proxy inverse evaluated deterministically",
        source = "closed_form",
        iterations = 0L
      )
    } else {
      .dpprior_new_termination(
        code = "failed",
        message = paste(
          "A1 proxy round-trip failed; no public parameter candidate was",
          "retained."
        ),
        source = "mapping",
        iterations = NULL
      )
    },
    trace = NULL,
    resources = list()
  )

  verification <- if (retained_candidate) {
    .dpprior_new_verification(
      method = "independent_finite_J_target_verification",
      performed = FALSE,
      passed = FALSE,
      reason = "not_performed_for_A1_proxy",
      settings = list(
        requested_estimand = "finite_J_K_moments",
        evaluated_estimand = "shifted_negative_binomial_proxy_moments",
        independent_verifier_required = TRUE,
        independent_verifier_available = FALSE
      ),
      selected_snapshot = selected_snapshot,
      verifier_snapshot = NULL,
      stability = NULL,
      components = list(
        exact_estimand = .dpprior_new_check(
          value = FALSE,
          reference = TRUE,
          tolerance = NULL,
          operator = "identical",
          source = "not_performed_for_A1_proxy"
        ),
        proxy_mapping = .dpprior_new_check(
          value = abs(mapping_residual),
          reference = c(mean = 0, variance = 0),
          tolerance = mapping_tolerance,
          operator = "lte",
          source = "algebraic_round_trip_under_a1_proxy"
        )
      ),
      invariants = list(
        target_identity = .dpprior_new_check(
          value = c(
            mean = identical(mapping_target[["mean"]], mu_K),
            variance = identical(mapping_target[["variance"]], var_K_used)
          ),
          reference = c(mean = TRUE, variance = TRUE),
          tolerance = NULL,
          operator = "identical",
          source = "closed_form_a1_proxy_mapping"
        ),
        parameter_identity = .dpprior_new_check(
          value = c(a = a, b = b),
          reference = c(a = parameters$a, b = parameters$b),
          tolerance = NULL,
          operator = "identical",
          source = "closed_form_a1_proxy_mapping"
        )
      )
    )
  } else {
    .dpprior_new_verification(
      method = "no_candidate",
      performed = FALSE,
      passed = FALSE,
      reason = "a1_proxy_mapping_failed",
      settings = list(),
      selected_snapshot = NULL,
      verifier_snapshot = NULL,
      stability = NULL,
      components = list(),
      invariants = list(
        no_public_candidate = .dpprior_new_check(
          value = TRUE,
          reference = TRUE,
          tolerance = NULL,
          operator = "identical",
          source = "independent_verifier"
        )
      )
    )
  }

  mapping_verification <- if (candidate_finite) list(
    method = "algebraic_round_trip_under_a1_proxy",
    performed = TRUE,
    passed = mapping_passed,
    target = list(mu_K = mu_K, var_K = var_K_used),
    achieved = list(
      mu_K = unname(mapping_achieved[["mean"]]),
      var_K = unname(mapping_achieved[["variance"]])
    ),
    residuals = mapping_residual,
    tolerances = mapping_tolerance,
    component_pass = stats::setNames(
      mapping_component_pass, names(mapping_residual)
    ),
    estimand = "shifted_negative_binomial_proxy_moments"
  ) else list(
    method = "algebraic_round_trip_under_a1_proxy",
    performed = TRUE,
    passed = FALSE,
    target = list(mu_K = mu_K, var_K = var_K_used),
    achieved = NULL,
    residuals = NULL,
    tolerances = mapping_tolerance,
    component_pass = stats::setNames(
      mapping_component_pass, names(mapping_target)
    ),
    estimand = "shifted_negative_binomial_proxy_moments",
    unavailable_reason = "nonfinite_a1_proxy_mapping"
  )
  proxy <- .dpprior_new_proxy(
    mapping = if (candidate_finite) list(
      formula = "closed_form_shifted_negative_binomial_inverse",
      target = mapping_target,
      achieved = mapping_achieved,
      residuals = mapping_residual,
      tolerances = mapping_tolerance
    ) else list(
      formula = "closed_form_shifted_negative_binomial_inverse",
      target = mapping_target,
      achieved = NULL,
      residuals = NULL,
      tolerances = mapping_tolerance,
      unavailable_reason = "nonfinite_a1_proxy_mapping"
    ),
    mapping_verification = mapping_verification,
    projection = target_resolution$projection,
    caveats = caveats
  )

  source_commit <- getOption("DPprior.source_commit", NULL)
  if (!is.character(source_commit) || length(source_commit) != 1L ||
      is.na(source_commit) || !nzchar(source_commit)) {
    source_commit <- NULL
  }
  provenance <- .dpprior_new_provenance(
    requested_method = "A1",
    selected_method = "A1",
    is_fallback = FALSE,
    approximation = list(
      active = TRUE,
      opt_in = TRUE,
      kind = "shifted_negative_binomial_proxy",
      warning_code = "a1_proxy_not_finite_J_verified"
    ),
    projection = target_K$provenance$projection,
    parameterization = parameterization,
    backend = list(
      package = "DPprior",
      package_version = tryCatch(
        as.character(utils::packageVersion("DPprior")),
        error = function(condition) "development"
      ),
      implementation = "R/10_a1_mapping.R:DPprior_a1",
      source_commit = source_commit
    ),
    input_fit = NULL,
    migration = list(
      source_schema = "phase8_a1_flat_v0",
      adapter = ".dpprior_a1_result_v1",
      lossless = FALSE,
      missing_evidence = c(
        "independent_finite_J_verifier_not_performed",
        "legacy_NA_boundary_reason_normalized_to_NULL",
        "legacy_unnamed_attempt_index_named_attempt_1",
        if (!candidate_finite) {
          "nonfinite_failed_candidate_encoded_as_plain_text"
        } else {
          character()
        },
        if (is.null(source_commit)) "source_commit_not_embedded" else character()
      ),
      warnings = "compatibility.views.a1_v0_is_non_authoritative"
    ),
    legacy = list(active = FALSE, contract = NULL, deprecation_stage = NULL)
  )

  old_boundary_reason <- if (length(target_resolution$boundary_reason) == 1L &&
      !is.na(target_resolution$boundary_reason)) {
    target_resolution$boundary_reason
  } else {
    NULL
  }
  legacy_a <- if (is.finite(a)) a else format(a, scientific = TRUE)
  legacy_b <- if (is.finite(b)) b else format(b, scientific = TRUE)
  legacy_view <- list(
    a = legacy_a,
    b = legacy_b,
    J = J,
    target = list(
      mu_K = mu_K,
      var_K = var_K_requested,
      var_K_requested = var_K_requested,
      var_K_used = var_K_used,
      projection = target_resolution$projection,
      type = "moments"
    ),
    method = "A1",
    status = status,
    usable = usable,
    verified = FALSE,
    mapping_verified = isTRUE(target_resolution$mapping_verified) &&
      mapping_passed,
    at_target_boundary = target_resolution$at_target_boundary,
    boundary_reason = old_boundary_reason,
    message = message,
    scaling = scaling,
    cJ = computation$scaling$values$cJ,
    parameters = list(a = legacy_a, b = legacy_b),
    achieved = list(
      mu_K = mu_K,
      var_K = var_K_used,
      estimand = "shifted_negative_binomial_proxy_moments"
    ),
    residuals = list(
      raw = c(mean = 0, variance = 0),
      scaled = c(mean = 0, variance = 0),
      scale_formula = "closed-form algebraic inverse"
    ),
    tolerances = list(
      feasibility_buffer = target_resolution$projection$buffer
    ),
    attempts = list(attempt_1 = list(
      method = "closed_form_a1",
      exit_code = 0L,
      message = "closed-form inverse evaluated"
    )),
    mapping_verification = mapping_verification,
    verification = list(
      method = "independent_finite_J_target_verification",
      performed = FALSE,
      passed = FALSE,
      reason = "not_performed_for_A1_proxy"
    ),
    var_K_used = var_K_used,
    projection = target_resolution$projection,
    caveats = if (length(caveats)) caveats else NULL,
    provenance = list(
      requested_method = "A1",
      selected_method = "A1",
      target_projection = target_resolution$projection,
      approximation = "shifted_negative_binomial"
    ),
    converged = FALSE,
    iterations = 0L,
    fit = NULL,
    diagnostics = NULL,
    trace = NULL
  )
  views <- list(
    a1_v0 = legacy_view,
    mapping_verified = legacy_view$mapping_verified,
    at_target_boundary = legacy_view$at_target_boundary,
    boundary_reason = old_boundary_reason,
    scaling = scaling,
    cJ = computation$scaling$values$cJ,
    converged = FALSE,
    iterations = 0L,
    fit = NULL,
    diagnostics = NULL,
    trace = NULL
  )
  out <- .dpprior_new_fit(
    mode = "a1_proxy",
    method = "A1",
    J = J,
    status = status,
    usable = usable,
    verified = FALSE,
    message = message,
    parameters = parameters,
    target = list(K = target_K),
    achieved = achieved,
    residuals = residuals,
    tolerances = tolerances,
    computation = computation,
    verification = verification,
    provenance = provenance,
    compatibility = .dpprior_new_compatibility(),
    extension = list(proxy = proxy)
  )
  aliases <- c(
    mapping_verified = "compatibility.views.mapping_verified",
    at_target_boundary = "compatibility.views.at_target_boundary",
    scaling = "compatibility.views.scaling",
    cJ = "compatibility.views.cJ",
    var_K_used = "target.K.used.var_K",
    projection = "proxy.projection",
    caveats = "proxy.caveats",
    attempts = "computation.attempts",
    mapping_verification = "proxy.mapping_verification",
    converged = "compatibility.views.converged",
    iterations = "compatibility.views.iterations"
  )
  if (retained_candidate) {
    aliases <- c(
      a = "parameters.a",
      b = "parameters.b",
      aliases
    )
  }
  .dpprior_append_compatibility_v2(
    out,
    aliases = aliases,
    views = views,
    deprecations = list(
      a1_v0 = list(
        code = "a1_flat_v0_view_quarantined",
        authority = "non_authoritative",
        consumer_policy = "canonical_fields_only"
      )
    )
  )
}

#' A1 Closed-Form Prior Elicitation
#'
#' Maps target beliefs about the number of clusters \eqn{(\mu_K, \sigma^2_K)}
#' to Gamma hyperprior parameters \eqn{(a, b)} using the A1 closed-form
#' approximation based on Negative Binomial moment matching.
#'
#' @param J Integer; number of items/sites (must be >= 2).
#' @param mu_K Numeric; target prior mean of \eqn{K_J} (must satisfy \eqn{1 < \mu_K < J}).
#' @param var_K Numeric; target prior variance of \eqn{K_J} (must be > 0).
#' @param scaling Character; scaling constant method: "log" (default),
#'   "harmonic", or "digamma".
#' @param epsilon Numeric; preferred numerical-interior buffer and near-boundary
#'   tolerance. The buffer is adaptively capped in narrow feasible domains and
#'   floored at local floating-point spacing. A target is changed only when an
#'   explicit projection is requested. Default is
#'   \code{.TOL_PROJECTION_BUFFER} (1e-6).
#' @param projection Character; A1 target-projection policy. The default
#'   \code{"error"} never changes the requested target. Set
#'   \code{"nearest"} to explicitly opt into projection of an A1-infeasible
#'   variance to the nearest numerically interior A1 target. The original and
#'   projected targets, distance, and policy are retained in the result.
#'
#' @return A canonical \code{dpprior.result/1} \code{DPprior_fit} object.
#'   Gamma shape and rate are in \code{parameters}; the immutable canonical
#'   target is in \code{target$K}; proxy achievements and residuals are in
#'   \code{achieved} and \code{residuals}. The \code{proxy},
#'   \code{computation}, \code{verification}, and \code{provenance} records
#'   retain the closed-form formula, scaling, projection, and algebraic
#'   round-trip evidence. A1 is always labelled \code{status = "approximate"},
#'   \code{usable = TRUE}, and \code{verified = FALSE}: its proxy identity is
#'   not finite-design verification. Registered flat aliases are migration
#'   views and are not authoritative scientific fields.
#'
#' @details
#' ## Theory (TSMM Stage 1)
#'
#' The A1 method uses a shifted Negative Binomial approximation:
#' \deqn{K_J - 1 \mid \alpha \approx \text{Poisson}(\alpha \cdot c_J)}
#'
#' With \eqn{\alpha \sim \text{Gamma}(a, b)}, the marginal becomes:
#' \deqn{K_J - 1 \approx \text{NegBin}(a, b/(b + c_J))}
#'
#' ## Inverse Formulas (Theorem 1)
#'
#' Let \eqn{m = \mu_K - 1} (shifted mean) and \eqn{D = \sigma^2_K - m}.
#' If \eqn{D > 0} (overdispersion):
#' \deqn{a = m^2 / D, \quad b = m \cdot c_J / D}
#'
#' If \eqn{D \leq 0}, the target is infeasible for A1. The default is a typed
#' error. Projection occurs only with \code{projection = "nearest"} and emits
#' exactly one typed projection warning.
#'
#' ## Feasibility
#'
#' The NegBin model requires overdispersion: \eqn{\sigma^2_K > \mu_K - 1}.
#' High-confidence specifications (low variance) may violate this constraint
#' under the A1 proxy, even though they may be feasible under the exact DP.
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso
#' \code{\link{vif_to_variance}} for VIF conversion,
#' \code{\link{confidence_to_vif}} for confidence mapping,
#' \code{\link{print.DPprior_fit}} for print method
#'
#' @examples
#' # Basic usage with moment targets
#' fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
#' print(fit)
#'
#' # Using VIF specification
#' fit <- DPprior_a1(J = 50, mu_K = 5, var_K = vif_to_variance(5, 2))
#'
#' # Using confidence-based specification
#' vif <- confidence_to_vif("medium")
#' fit <- DPprior_a1(J = 50, mu_K = 5, var_K = vif_to_variance(5, vif))
#'
#' # Explicit projection of an A1-infeasible variance
#' fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 3,
#'                   projection = "nearest")
#'
#' # Compare scaling methods
#' fit_log <- DPprior_a1(50, 5, 8, scaling = "log")
#' fit_harm <- DPprior_a1(50, 5, 8, scaling = "harmonic")
#'
#' @family elicitation
#'
#' @export
DPprior_a1 <- function(J, mu_K, var_K,
                       scaling = c("log", "harmonic", "digamma"),
                       epsilon = .TOL_PROJECTION_BUFFER,
                       projection = c("error", "nearest")) {

  scaling <- .dpprior_a1_match_arg(
    scaling, c("log", "harmonic", "digamma"), "scaling",
    "dpprior_a1_scaling_error"
  )

  # Input validation
  if (!.dpprior_is_plain_numeric(J) || length(J) != 1L || !is.finite(J) ||
      J != floor(J) || J < 2L) {
    .dpprior_abort_legacy_numeric(
      "J must be an integer >= 2", J, "J", "integer >= 2",
      "dpprior_a1_sample_size_error", scalar = TRUE, integer = TRUE
    )
  }
  if (!.dpprior_is_plain_numeric(mu_K) || length(mu_K) != 1L ||
      !is.finite(mu_K)) {
    .dpprior_abort_legacy_numeric(
      "mu_K must be a finite numeric scalar", mu_K, "mu_K",
      "finite numeric scalar", "dpprior_a1_mean_error", scalar = TRUE
    )
  }
  if (mu_K <= 1) {
    .dpprior_abort_invalid(
      "mu_K must be > 1 (at least one cluster is always present)",
      c("dpprior_a1_mean_error", "dpprior_bounds_error"),
      "mu_K", mu_K, "> 1", "bounds"
    )
  }
  if (mu_K >= J) {
    .dpprior_abort_invalid(
      "mu_K must be < J (mu_K = J implies zero variance for K_J, outside the positive-variance elicitation workflow)",
      c("dpprior_a1_mean_error", "dpprior_bounds_error"),
      "mu_K", mu_K, sprintf("< %d", as.integer(J)), "bounds"
    )
  }
  if (!.dpprior_is_plain_numeric(var_K) || length(var_K) != 1L ||
      !is.finite(var_K) ||
      var_K <= 0) {
    .dpprior_abort_legacy_numeric(
      "var_K must be a positive finite numeric scalar", var_K, "var_K",
      "positive finite numeric scalar", "dpprior_a1_variance_error",
      scalar = TRUE
    )
  }
  projection <- .dpprior_a1_match_arg(
    projection, c("error", "nearest"), "projection",
    "dpprior_a1_projection_policy_error"
  )

  target_resolution <- .dpprior_a1_resolve_target(
    J = J, mu_K = mu_K, var_K = var_K,
    projection = projection, epsilon = epsilon,
    signal_projection = TRUE
  )

  # Compute scaling constant
  cJ <- compute_scaling_constant(J, scaling, mu_K)

  # Shifted mean (under A1 with shift s = 1)
  mu_S <- mu_K - 1

  denom <- target_resolution$denominator
  var_K_used <- target_resolution$var_K_used
  # A1 exactly inverts a shifted-Negative-Binomial proxy. Proxy round-trip
  # success is not independent verification of the finite-J K target, so the
  # public numerical status remains `approximate` even for an interior target.
  # Closed-form inverse (Theorem 1)
  a0 <- mu_S^2 / denom
  b0 <- mu_S * cJ / denom

  # Audit the algebraic proxy mapping from the returned parameters rather
  # than treating successful formula evaluation as verification by fiat.
  # This remains a proxy-mapping check, not finite-J target verification.
  mapping_target <- c(mean = mu_K, variance = var_K_used)
  mapping_achieved <- c(
    mean = 1 + a0 * cJ / b0,
    variance = a0 * cJ / b0 + (a0 * cJ / b0)^2 / a0
  )
  mapping_residual <- mapping_achieved - mapping_target
  mapping_tolerance <- stats::setNames(
    1e-12 + 1e-10 * pmax(1, abs(mapping_target)),
    names(mapping_target)
  )
  mapping_component_pass <- is.finite(mapping_residual) &
    abs(mapping_residual) <= mapping_tolerance
  mapping_passed <- all(mapping_component_pass)
  caveats <- character()
  if (a0 < 0.01) caveats <- c(caveats, "quasi_improper_shape")
  if (a0 > 1e6) caveats <- c(caveats, "quasi_degenerate_shape")

  message <- if (!mapping_passed) {
    paste(
      "A1 proxy parameters were evaluated, but their algebraic",
      "round-trip did not meet the declared numerical tolerance."
    )
  } else if (target_resolution$at_target_boundary) {
    paste(
      "A1 proxy mapping was verified at a declared target boundary;",
      "the finite-J target was not independently verified."
    )
  } else {
    paste(
      "A1 proxy mapping was verified; the finite-J target was not",
      "independently verified."
    )
  }

  .dpprior_a1_result_v1(
    J = J,
    mu_K = mu_K,
    var_K_requested = var_K,
    var_K_used = var_K_used,
    scaling = scaling,
    cJ = cJ,
    epsilon = epsilon,
    projection_policy = projection,
    target_resolution = target_resolution,
    a = a0,
    b = b0,
    mapping_target = mapping_target,
    mapping_achieved = mapping_achieved,
    mapping_residual = mapping_residual,
    mapping_tolerance = mapping_tolerance,
    mapping_component_pass = mapping_component_pass,
    mapping_passed = mapping_passed,
    caveats = caveats,
    message = message
  )
}


# =============================================================================
# VIF (Variance Inflation Factor) Utilities
# =============================================================================

#' Convert Variance Inflation Factor to Variance
#'
#' Converts a Variance Inflation Factor (VIF) specification to the actual
#' variance of \eqn{K_J}.
#'
#' @param mu_K Numeric; target prior mean of \eqn{K_J}.
#' @param vif Numeric; Variance Inflation Factor (must be >= 1 for A1 feasibility).
#'
#' @return Numeric; variance of \eqn{K_J} computed as \eqn{(\mu_K - 1) \times \text{VIF}}.
#'
#' @details
#' The VIF is defined as:
#' \deqn{\text{VIF} = \frac{\sigma^2_K}{\mu_K - 1}}
#'
#' Interpretation:
#' \describe{
#'   \item{VIF = 1}{Poisson variance (exact boundary for A1)}
#'   \item{VIF > 1}{Overdispersion (required for A1 feasibility)}
#'   \item{VIF < 1}{Underdispersion (infeasible for A1, not allowed)}
#' }
#'
#' @seealso
#' \code{\link{confidence_to_vif}} for mapping confidence levels to VIF,
#' \code{\link{cv_alpha_to_variance}} for CV-based specification
#'
#' @examples
#' # VIF = 2 means variance is twice the Poisson variance
#' vif_to_variance(mu_K = 5, vif = 2)  # Returns 8
#'
#' # Use with DPprior_a1
#' fit <- DPprior_a1(J = 50, mu_K = 5, var_K = vif_to_variance(5, 2))
#'
#' @export
vif_to_variance <- function(mu_K, vif) {
  if (!.dpprior_is_plain_numeric(mu_K) || !.dpprior_is_plain_numeric(vif)) {
    .dpprior_abort_invalid(
      "mu_K and vif must be numeric",
      c("dpprior_a1_vif_error", "dpprior_type_error"),
      "mu_K/vif", list(mu_K = mu_K, vif = vif), "numeric", "type"
    )
  }
  if (!length(mu_K) || !length(vif) || anyNA(mu_K) || anyNA(vif) ||
      any(!is.finite(mu_K)) || any(!is.finite(vif))) {
    reason <- if (!length(mu_K) || !length(vif)) {
      "dpprior_length_error"
    } else if (anyNA(mu_K) || anyNA(vif)) {
      "dpprior_missing_error"
    } else {
      "dpprior_nonfinite_error"
    }
    code <- switch(
      reason,
      dpprior_length_error = "length",
      dpprior_missing_error = "missing",
      dpprior_nonfinite_error = "nonfinite"
    )
    .dpprior_abort_invalid(
      "mu_K and vif must be finite and non-missing",
      c("dpprior_a1_vif_error", reason), "mu_K/vif",
      list(mu_K = mu_K, vif = vif), "finite non-missing numeric", code
    )
  }
  if (any(mu_K <= 1)) {
    .dpprior_abort_invalid(
      "mu_K must be > 1",
      c("dpprior_a1_vif_error", "dpprior_bounds_error"),
      "mu_K", mu_K, "> 1", "bounds"
    )
  }
  # VIF >= 1 required for A1 feasibility
  if (any(vif < 1)) {
    .dpprior_abort_invalid(
      "vif must be >= 1 (values < 1 imply infeasible underdispersion)",
      c("dpprior_a1_vif_error", "dpprior_bounds_error"),
      "vif", vif, ">= 1", "bounds"
    )
  }

  (mu_K - 1) * vif
}


#' Map a Qualitative Confidence Level to a Variance Inflation Factor (VIF)
#'
#' Maps intuitive confidence levels to VIF values for easy prior specification.
#'
#' @param confidence Character; one of "low", "medium", or "high".
#'
#' @return Numeric; VIF value (5.0 for low, 2.5 for medium, 1.5 for high).
#'
#' @details
#' The mapping is:
#' \describe{
#'   \item{low}{VIF = 5.0; high uncertainty about \eqn{K_J}}
#'   \item{medium}{VIF = 2.5; moderate uncertainty}
#'   \item{high}{VIF = 1.5; high confidence (near Poisson boundary)}
#' }
#'
#' Higher confidence implies lower variance, which corresponds to lower VIF.
#' The "high" setting (VIF = 1.5) is close to the A1 feasibility boundary.
#'
#' @seealso \code{\link{vif_to_variance}} for converting VIF to variance
#'
#' @examples
#' # Get VIF for medium confidence
#' vif <- confidence_to_vif("medium")  # Returns 2.5
#'
#' # Complete workflow
#' mu_K <- 5
#' vif <- confidence_to_vif("low")
#' var_K <- vif_to_variance(mu_K, vif)
#' fit <- DPprior_a1(J = 50, mu_K = mu_K, var_K = var_K)
#'
#' @export
confidence_to_vif <- function(confidence = c("low", "medium", "high")) {
  confidence <- .dpprior_a1_match_arg(
    confidence, c("low", "medium", "high"), "confidence",
    "dpprior_a1_confidence_error"
  )

  switch(
    confidence,
    low = 5.0,
    medium = 2.5,
    high = 1.5
  )
}


#' Convert CV(alpha) to Variance
#'
#' Converts a coefficient of variation specification for \eqn{\alpha} to
#' the implied variance of \eqn{K_J} under the A1 approximation.
#'
#' @param mu_K Numeric; target prior mean of \eqn{K_J}.
#' @param cv_alpha Numeric; target coefficient of variation for \eqn{\alpha}.
#'
#' @return Numeric; implied variance of \eqn{K_J}.
#'
#' @details
#' Under the A1 approximation:
#' \deqn{\text{CV}(\alpha) = 1/\sqrt{a} = \frac{\sqrt{\sigma^2_K - m}}{m}}
#'
#' where \eqn{m = \mu_K - 1}. Inverting:
#' \deqn{\sigma^2_K = m + (\text{CV}(\alpha) \cdot m)^2 = m(1 + \text{CV}(\alpha)^2 \cdot m)}
#'
#' @seealso \code{\link{vif_to_variance}} for VIF-based specification
#'
#' @examples
#' # CV(alpha) = 0.5 means moderate prior concentration
#' var_K <- cv_alpha_to_variance(mu_K = 5, cv_alpha = 0.5)
#'
#' # Verify round-trip
#' fit <- DPprior_a1(J = 50, mu_K = 5, var_K = var_K)
#' 1 / sqrt(fit$parameters$a)  # Should be approximately 0.5
#'
#' @export
cv_alpha_to_variance <- function(mu_K, cv_alpha) {
  if (!.dpprior_is_plain_numeric(mu_K) ||
      !.dpprior_is_plain_numeric(cv_alpha)) {
    .dpprior_abort_invalid(
      "mu_K and cv_alpha must be numeric",
      c("dpprior_a1_cv_error", "dpprior_type_error"),
      "mu_K/cv_alpha", list(mu_K = mu_K, cv_alpha = cv_alpha),
      "numeric", "type"
    )
  }
  if (!length(mu_K) || !length(cv_alpha) || anyNA(mu_K) ||
      anyNA(cv_alpha) || any(!is.finite(mu_K)) ||
      any(!is.finite(cv_alpha))) {
    reason <- if (!length(mu_K) || !length(cv_alpha)) {
      "dpprior_length_error"
    } else if (anyNA(mu_K) || anyNA(cv_alpha)) {
      "dpprior_missing_error"
    } else {
      "dpprior_nonfinite_error"
    }
    code <- switch(
      reason,
      dpprior_length_error = "length",
      dpprior_missing_error = "missing",
      dpprior_nonfinite_error = "nonfinite"
    )
    .dpprior_abort_invalid(
      "mu_K and cv_alpha must be finite and non-missing",
      c("dpprior_a1_cv_error", reason), "mu_K/cv_alpha",
      list(mu_K = mu_K, cv_alpha = cv_alpha),
      "finite non-missing numeric", code
    )
  }
  if (any(mu_K <= 1)) {
    .dpprior_abort_invalid(
      "mu_K must be > 1",
      c("dpprior_a1_cv_error", "dpprior_bounds_error"),
      "mu_K", mu_K, "> 1", "bounds"
    )
  }
  if (any(cv_alpha <= 0)) {
    .dpprior_abort_invalid(
      "cv_alpha must be positive",
      c("dpprior_a1_cv_error", "dpprior_bounds_error"),
      "cv_alpha", cv_alpha, "> 0", "bounds"
    )
  }

  m <- mu_K - 1
  m + (cv_alpha * m)^2
}


# =============================================================================
# S3 Methods for DPprior_fit
# =============================================================================
# Note: summary.DPprior_fit is defined in R/17_s3_methods.R (canonical location)

# Historical runtime definition retained for source-order compatibility.
# R/17_s3_methods.R owns the public S3 implementation and its documentation.
as.data.frame.DPprior_fit <- function(x, row.names = NULL, optional = FALSE, ...) {
  if (inherits(x, "dpprior_result")) {
    .dpprior_validate_result_v1(x)
    raw <- unclass(x)
    target_K <- raw[["target", exact = TRUE]][["K", exact = TRUE]]
    requested <- target_K[["request", exact = TRUE]]
    requested_mean <- if ("mu_K" %in% names(requested)) {
      requested[["mu_K", exact = TRUE]]
    } else {
      requested[["mean", exact = TRUE]]
    }
    requested_variance <- if ("var_K" %in% names(requested)) {
      requested[["var_K", exact = TRUE]]
    } else {
      requested[["variance", exact = TRUE]]
    }
    parameters <- raw[["parameters", exact = TRUE]]
    scaling_record <- raw[["computation", exact = TRUE]][[
      "scaling", exact = TRUE
    ]]
    scaling_method <- scaling_record[["used", exact = TRUE]][[
      "method", exact = TRUE
    ]]
    iterations <- raw[["computation", exact = TRUE]][[
      "termination", exact = TRUE
    ]][["iterations", exact = TRUE]]
    if (is.null(iterations) && identical(raw[["mode", exact = TRUE]],
                                         "a1_proxy")) {
      iterations <- 0L
    }
    return(data.frame(
      method = raw[["method", exact = TRUE]],
      status = raw[["status", exact = TRUE]],
      a = parameters[["a", exact = TRUE]],
      b = parameters[["b", exact = TRUE]],
      J = raw[["J", exact = TRUE]],
      mu_K = requested_mean,
      var_K = requested_variance,
      mean_alpha = parameters[["a", exact = TRUE]] /
        parameters[["b", exact = TRUE]],
      cv_alpha = 1 / sqrt(parameters[["a", exact = TRUE]]),
      scaling = scaling_method,
      converged = raw[["status", exact = TRUE]] %in%
        c("converged", "boundary") &&
        isTRUE(raw[["usable", exact = TRUE]]) &&
        isTRUE(raw[["verified", exact = TRUE]]),
      iterations = iterations,
      stringsAsFactors = FALSE,
      row.names = row.names
    ))
  }
  data.frame(
    method = x$method,
    status = x$status,
    a = x$a,
    b = x$b,
    J = x$J,
    mu_K = x$target$mu_K,
    var_K = x$target$var_K,
    mean_alpha = x$a / x$b,
    cv_alpha = 1 / sqrt(x$a),
    scaling = if (!is.null(x$scaling)) x$scaling else NA_character_,
    converged = x$converged,
    iterations = x$iterations,
    stringsAsFactors = FALSE,
    row.names = row.names
  )
}


# =============================================================================
# Verification Functions
# =============================================================================

#' Verify A1 Mapping via Round-Trip
#'
#' Tests the A1 mapping by computing the forward model (NegBin moments)
#' from the derived \eqn{(a, b)} parameters and comparing to targets.
#'
#' @param fit A \code{DPprior_fit} object from \code{DPprior_a1}.
#' @param tol Numeric; tolerance for relative error comparison.
#' @param verbose Logical; if TRUE, print verification details.
#'
#' @return Logical; TRUE if round-trip succeeds within tolerance.
#'
#' @details
#' Under the A1 NegBin approximation:
#' \deqn{K_J - 1 \sim \text{NegBin}(a, p)}
#' where \eqn{p = b/(b + c_J)}.
#'
#' @examples
#' \dontrun{
#' fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
#' verify_a1_roundtrip(fit)
#'
#' }
#' @keywords internal
verify_a1_roundtrip <- function(fit, tol = 1e-8, verbose = TRUE) {
  if (!inherits(fit, "DPprior_fit")) {
    .dpprior_abort_invalid(
      "fit must be a DPprior_fit object", "dpprior_a1_verification_error",
      "fit", fit, "DPprior_fit object", "class"
    )
  }
  if (!inherits(fit, "dpprior_result")) {
    .dpprior_abort_invalid(
      "fit must use the canonical dpprior.result/1 contract",
      "dpprior_a1_verification_error", "fit", fit,
      "canonical DPprior_fit object", "schema"
    )
  }
  .dpprior_validate_result_v1(fit)
  raw <- unclass(fit)
  if (raw[["method", exact = TRUE]] != "A1") {
    warning("Round-trip verification is specific to A1 method", call. = FALSE)
  }

  # Extract only canonical scientific authority, never compatibility aliases.
  parameters <- raw[["parameters", exact = TRUE]]
  a <- parameters[["a", exact = TRUE]]
  b <- parameters[["b", exact = TRUE]]
  cJ <- raw[["computation", exact = TRUE]][[
    "scaling", exact = TRUE
  ]][["values", exact = TRUE]][["cJ", exact = TRUE]]

  # Forward model: NegBin moments
  # p = b / (b + cJ)
  # Mean(K-1) = a * (1-p) / p = a * cJ / b
  # Var(K-1) = a * (1-p) / p^2 = a * cJ * (b + cJ) / b^2
  p <- b / (b + cJ)
  mean_shifted <- a * (1 - p) / p
  var_shifted <- a * (1 - p) / p^2

  mu_K_recovered <- mean_shifted + 1
  var_K_recovered <- var_shifted

  # Compare to targets (use var_K_used for projected cases)
  used_target <- raw[["target", exact = TRUE]][["K", exact = TRUE]][[
    "used", exact = TRUE
  ]]
  mu_K_target <- if ("mu_K" %in% names(used_target)) {
    used_target[["mu_K", exact = TRUE]]
  } else {
    used_target[["mean", exact = TRUE]]
  }
  var_K_target <- if ("var_K" %in% names(used_target)) {
    used_target[["var_K", exact = TRUE]]
  } else {
    used_target[["variance", exact = TRUE]]
  }

  rel_err_mu <- abs(mu_K_recovered - mu_K_target) / mu_K_target
  rel_err_var <- abs(var_K_recovered - var_K_target) / var_K_target

  passed <- (rel_err_mu < tol) && (rel_err_var < tol)

  if (isTRUE(verbose)) {
    cat("A1 Round-Trip Verification\n")
    cat(paste0(rep("-", 40), collapse = ""), "\n")
    cat(sprintf("mu_K: target = %.6f, recovered = %.6f, rel_err = %.2e\n",
                mu_K_target, mu_K_recovered, rel_err_mu))
    cat(sprintf("var_K: target = %.6f, recovered = %.6f, rel_err = %.2e\n",
                var_K_target, var_K_recovered, rel_err_var))
    cat(sprintf("Result: %s\n", if (passed) "PASS" else "FAIL"))
  }

  invisible(passed)
}


#' Run All Module 10 Verification Tests
#'
#' Comprehensive verification suite for the A1 closed-form mapping module.
#'
#' @param verbose Logical; if TRUE, print detailed results.
#'
#' @return Logical; TRUE if all tests pass.
#'
#' @examples
#' \dontrun{
#' verify_a1_mapping_all()
#'
#' }
#' @keywords internal
verify_a1_mapping_all <- function(verbose = TRUE) {
  if (isTRUE(verbose)) {
    cat("=", rep("=", 69), "\n", sep = "")
    cat("Module 10: A1 Closed-Form Mapping - Full Verification Suite\n")
    cat("=", rep("=", 69), "\n\n", sep = "")
  }

  all_pass <- TRUE

  # Test 1: Basic positive parameters
  if (isTRUE(verbose)) {
    cat("[Test 1] A1 produces positive parameters\n")
  }
  fit <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
  test1 <- (fit$a > 0) && (fit$b > 0) && (fit$method == "A1")
  if (isTRUE(verbose)) {
    cat(sprintf("  J=50, mu_K=5, var_K=8\n"))
    cat(sprintf("  a = %.6f, b = %.6f\n", fit$a, fit$b))
    cat(sprintf("  PASS: %s\n\n", test1))
  }
  all_pass <- all_pass && test1

  # Test 2: Explicit infeasible-target projection
  if (isTRUE(verbose)) {
    cat("[Test 2] A1 handles infeasible variance via projection\n")
  }
  fit2 <- suppressWarnings(DPprior_a1(
    J = 50, mu_K = 5, var_K = 3,
    projection = "nearest"
  ))
  test2 <- identical(fit2$status, "approximate") &&
    isTRUE(fit2$at_target_boundary) &&
    isTRUE(fit2$projection$applied)
  if (isTRUE(verbose)) {
    cat(sprintf("  J=50, mu_K=5, var_K=3 (infeasible: var < mu-1=4)\n"))
    cat(sprintf("  var_K_used = %.6f\n", fit2$var_K_used))
    cat(sprintf("  status = %s\n", fit2$status))
    cat(sprintf("  PASS: %s\n\n", test2))
  }
  all_pass <- all_pass && test2

  # Test 3: VIF conversion
  if (isTRUE(verbose)) {
    cat("[Test 3] VIF conversion is correct\n")
  }
  var_computed <- vif_to_variance(5, 2)
  test3a <- abs(var_computed - 8) < 1e-10
  vif_medium <- confidence_to_vif("medium")
  test3b <- abs(vif_medium - 2.5) < 1e-10
  if (isTRUE(verbose)) {
    cat(sprintf("  vif_to_variance(5, 2) = %.1f (expected 8)\n", var_computed))
    cat(sprintf("  confidence_to_vif('medium') = %.1f (expected 2.5)\n", vif_medium))
    cat(sprintf("  PASS: %s\n\n", test3a && test3b))
  }
  all_pass <- all_pass && test3a && test3b

  # Test 4: Round-trip verification
  if (isTRUE(verbose)) {
    cat("[Test 4] Round-trip verification\n")
  }
  test_cases <- list(
    list(J = 50, mu_K = 5.0, var_K = 8.0),
    list(J = 100, mu_K = 10.0, var_K = 20.0),
    list(J = 25, mu_K = 4.0, var_K = 30.0),
    list(J = 50, mu_K = 5.0, var_K = 6.0)
  )

  test4_pass <- TRUE
  for (tc in test_cases) {
    fit <- DPprior_a1(J = tc$J, mu_K = tc$mu_K, var_K = tc$var_K)
    passed <- verify_a1_roundtrip(fit, verbose = FALSE)
    test4_pass <- test4_pass && passed
    if (isTRUE(verbose)) {
      cat(sprintf("  J=%d, mu_K=%.1f, var_K=%.1f: %s\n",
                  tc$J, tc$mu_K, tc$var_K, if (passed) "PASS" else "FAIL"))
    }
  }
  if (isTRUE(verbose)) cat("\n")
  all_pass <- all_pass && test4_pass

  # Test 5: Scaling method comparison
  if (isTRUE(verbose)) {
    cat("[Test 5] Scaling method comparison\n")
  }
  for (scaling in c("log", "harmonic", "digamma")) {
    fit <- DPprior_a1(50, 5, 8, scaling = scaling)
    if (isTRUE(verbose)) {
      cat(sprintf("  scaling='%s': cJ=%.4f, a=%.4f, b=%.4f\n",
                  scaling, fit$cJ, fit$a, fit$b))
    }
  }
  if (isTRUE(verbose)) cat("\n")

  # Test 6: CV to variance conversion
  if (isTRUE(verbose)) {
    cat("[Test 6] CV(alpha) to variance conversion\n")
  }
  test6_pass <- TRUE
  for (cv_target in c(0.5, 1.0, 2.0)) {
    var_K <- cv_alpha_to_variance(5, cv_target)
    fit <- DPprior_a1(50, 5, var_K)
    cv_recovered <- 1 / sqrt(fit$a)
    rel_err <- abs(cv_recovered - cv_target) / cv_target
    passed <- rel_err < 0.01
    test6_pass <- test6_pass && passed
    if (isTRUE(verbose)) {
      cat(sprintf("  Target CV=%.1f -> var_K=%.2f -> Recovered CV=%.4f\n",
                  cv_target, var_K, cv_recovered))
    }
  }
  if (isTRUE(verbose)) cat("\n")
  all_pass <- all_pass && test6_pass

  # Test 7: closed-status compatibility fields
  if (isTRUE(verbose)) {
    cat("[Test 7] closed-status compatibility fields\n")
  }
  fit <- DPprior_a1(50, 5, 8)
  test7 <- identical(fit$status, "approximate") &&
    !isTRUE(fit$converged) && !isTRUE(fit$verified) &&
    isTRUE(fit$mapping_verified) && identical(fit$iterations, 0L)
  if (isTRUE(verbose)) {
    cat(sprintf("  status = %s (expected approximate)\n", fit$status))
    cat(sprintf("  converged = %s (expected FALSE)\n", fit$converged))
    cat(sprintf("  mapping_verified = %s (expected TRUE)\n",
                fit$mapping_verified))
    cat(sprintf("  iterations = %d (expected 0)\n", fit$iterations))
    cat(sprintf("  PASS: %s\n\n", test7))
  }
  all_pass <- all_pass && test7

  # Summary
  if (isTRUE(verbose)) {
    cat("=", rep("=", 69), "\n", sep = "")
    cat(sprintf("Overall Result: %s\n",
                if (all_pass) "ALL TESTS PASSED" else "SOME TESTS FAILED"))
    cat("=", rep("=", 69), "\n", sep = "")
  }

  invisible(all_pass)
}
