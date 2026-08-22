# =============================================================================
# Canonical schema v1 kernel: constructors, trust boundary, and mutations
# =============================================================================

.schema23_truth_controls <- function(method) {
  switch(
    method,
    `A2-MN` = list(
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
    `A2-KL` = list(
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
    dual_anchor_hard_inequality = list(
      maxit = 250L, scan_points = 17L, scan_keep = 12L,
      profile_starts = 4L, log_bounds = c(-15, 15),
      root_tol = 1e-10, optim_reltol = 1e-10,
      penalty = 1e7, boundary_tol = 1e-5,
      constraint_abs_tol = 1e-6, constraint_rel_tol = 1e-6,
      verification_abs_tol = 1e-8, verification_rel_tol = 1e-6,
      perturbation_step = 1e-6, perturbation_abs_tol = 1e-4,
      auto_cap = list(scan_keep = FALSE, profile_starts = FALSE)
    ),
    `dual-soft` = list(
      max_iter = 100L, log_bounds = c(-15, 15),
      primary = list(
        maxit = 100L, fnscale = 1,
        parscale = c(log_shape = 1, log_rate = 1),
        ndeps = c(log_shape = 1e-3, log_rate = 1e-3)
      ),
      fallback = list(
        maxit = 200L, reltol = sqrt(.Machine$double.eps), fnscale = 1,
        parscale = c(log_shape = 1, log_rate = 1),
        ndeps = c(log_shape = 1e-3, log_rate = 1e-3)
      ),
      boundary_tol = 1e-6, verification_abs_tol = 1e-8,
      verification_rel_tol = 1e-6, objective_abs_tol = 1e-10,
      objective_rel_tol = 1e-8, stationarity_step = 1e-5,
      stationarity_tol = 1e-4, local_neighbor_step = 1e-3,
      selection_tolerance = 0, optimizer_adapter = "stats::optim"
    ),
    list(tolerance = 1e-8, selection_tolerance = 0)
  )
}


.schema23_settings <- function(method = "A2-MN",
                               parameterization = "log_ab") {
  list(
    method = method, controls = .schema23_truth_controls(method),
    parameterization = parameterization
  )
}


.schema23_check <- function(value = 0, reference = 0, tolerance = 1e-6,
                            operator = "abs_lte",
                            source = "independent_fixture_check") {
  .dpprior_new_check(
    value = value, reference = reference, tolerance = tolerance,
    operator = operator, source = source
  )
}


.schema23_hard_optimality <- function(available = TRUE,
                                      scales = list(mean = 5, variance = 2),
                                      tie_tolerance = 1e-8) {
  if (available) {
    return(list(
      performed = TRUE, passed = TRUE,
      selection_rule = "minimum_K_loss",
      K_scales = scales,
      selected_K_loss = 0, minimum_K_loss = 0,
      tie_tolerance = tie_tolerance, perturbation_passed = TRUE,
      source = "independent_fixture_check", unavailable_reason = NULL
    ))
  }
  list(
    performed = FALSE, passed = FALSE, selection_rule = NULL,
    K_scales = NULL, selected_K_loss = NULL, minimum_K_loss = NULL,
    tie_tolerance = NULL, perturbation_passed = NULL,
    source = "independent_fixture_check",
    unavailable_reason = "no finite feasible candidate"
  )
}


.schema23_soft_optimality <- function(available = TRUE) {
  if (available) {
    return(list(
      performed = TRUE, passed = TRUE,
      recorded_objective = 0, recomputed_objective = 0,
      objective_difference = 0, objective_tolerance = 1e-10 + 1e-8,
      objective_passed = TRUE,
      gradient = c(log_shape = 0, log_rate = 0),
      gradient_method = c(
        log_shape = "central finite-difference interior gradient",
        log_rate = "central finite-difference interior gradient"
      ),
      bound_state = c(log_shape = "interior", log_rate = "interior"),
      boundary_tolerance = 1e-6,
      stationarity_operator = c(
        log_shape = "abs_lte", log_rate = "abs_lte"
      ),
      stationarity_tolerance = 1e-4,
      component_pass = c(log_shape = TRUE, log_rate = TRUE),
      stationarity_passed = TRUE,
      local_base_objective = 0,
      neighbor_objectives = c(
        log_shape_plus = 0, log_shape_minus = 0,
        log_rate_plus = 0, log_rate_minus = 0,
        diagonal_pp = 0, diagonal_pm = 0,
        diagonal_mp = 0, diagonal_mm = 0
      ),
      neighbor_tolerances = c(
        log_shape_plus = 1e-10 + 1e-8,
        log_shape_minus = 1e-10 + 1e-8,
        log_rate_plus = 1e-10 + 1e-8,
        log_rate_minus = 1e-10 + 1e-8,
        diagonal_pp = 1e-10 + 1e-8,
        diagonal_pm = 1e-10 + 1e-8,
        diagonal_mp = 1e-10 + 1e-8,
        diagonal_mm = 1e-10 + 1e-8
      ),
      local_minimum_passed = TRUE, start_objective = 0,
      candidate_objective = 0, start_tolerance = 1e-10 + 1e-8,
      no_worse_start = TRUE, selection_tolerance = 0,
      source = "independent_refined_objective_verification",
      unavailable_reason = NULL
    ))
  }
  list(
    performed = FALSE, passed = FALSE, recorded_objective = NULL,
    recomputed_objective = NULL, objective_difference = NULL,
    objective_tolerance = NULL, objective_passed = NULL, gradient = NULL,
    gradient_method = NULL, bound_state = NULL, boundary_tolerance = NULL,
    stationarity_operator = NULL,
    stationarity_tolerance = NULL, component_pass = NULL,
    stationarity_passed = NULL,
    local_base_objective = NULL, neighbor_objectives = NULL,
    neighbor_tolerances = NULL, local_minimum_passed = NULL,
    start_objective = NULL, candidate_objective = NULL,
    start_tolerance = NULL, no_worse_start = NULL,
    selection_tolerance = NULL,
    source = "independent_fixture_check",
    unavailable_reason = "no finite optimizer candidate"
  )
}


.schema23_provenance <- function(method = "A2-MN",
                                 parameterization = "log_ab",
                                 approximation = FALSE,
                                 opt_in = FALSE,
                                 legacy = FALSE,
                                 input_fit = NULL) {
  .dpprior_new_provenance(
    requested_method = method,
    selected_method = method,
    is_fallback = FALSE,
    approximation = list(
      active = approximation,
      opt_in = opt_in,
      kind = if (approximation) "test_approximation" else NULL,
      warning_code = if (approximation) "test_approximation" else NULL
    ),
    projection = list(
      applied = FALSE, opt_in = FALSE, policy = NULL, record = NULL
    ),
    parameterization = parameterization,
    backend = list(
      package = "DPprior", package_version = "test",
      implementation = "schema_fixture", source_commit = NULL
    ),
    input_fit = input_fit,
    migration = list(
      source_schema = "native", adapter = "none", lossless = TRUE,
      missing_evidence = character(), warnings = character()
    ),
    legacy = list(
      active = legacy,
      contract = if (legacy) "DPprior_dual_v1" else NULL,
      deprecation_stage = if (legacy) "v2_retained" else NULL
    )
  )
}


.schema23_computation <- function(method = "A2-MN",
                                  parameters = NULL,
                                  M_selected = NULL,
                                  M_verification = NULL,
                                  scaling = .dpprior_new_scaling()) {
  has_candidate <- !is.null(parameters)
  attempts <- if (has_candidate) {
    list(.dpprior_new_attempt(
      id = "attempt-1", stage = "primary", method = method,
      start = list(a = 1, b = 1),
      bounds = list(lower = c(-10, -10), upper = c(10, 10)),
      control = list(maxit = 100L), exit_code = 0L,
      message = "optimizer completed", iterations = 5L,
      evaluations = list(function_count = 12L, gradient_count = 6L),
      candidate_parameters = parameters, candidate_objective = 0,
      elapsed_seconds = 0.01, warnings = character(), error = NULL,
      selected = TRUE, reason_code = "selected", unavailable = character()
    ))
  } else {
    list()
  }
  orders <- .dpprior_new_orders(
    M_requested = M_selected,
    M_selected = M_selected,
    M_verification_required = M_verification,
    M_verification_used = M_verification,
    requested_reason = if (is.null(M_selected)) "not_applicable" else "input",
    selected_reason = if (is.null(M_selected)) "not_applicable" else "calibration",
    verification_required_reason = if (is.null(M_verification)) {
      "not_applicable"
    } else {
      "policy"
    },
    verification_used_reason = if (is.null(M_verification)) {
      "not_applicable"
    } else {
      "independent_verifier"
    }
  )
  .dpprior_new_computation(
    request = .schema23_settings(method),
    used = .schema23_settings(method),
    orders = orders,
    scaling = scaling,
    attempts = attempts,
    selected_attempt_id = if (has_candidate) "attempt-1" else NULL,
    fallback = .dpprior_new_fallback(),
    termination = .dpprior_new_termination(
      code = if (has_candidate) "selected" else "deterministic",
      source = if (has_candidate) "optimizer" else "constructor",
      iterations = if (has_candidate) 5L else NULL
    ),
    trace = NULL,
    resources = list(elapsed_seconds = if (has_candidate) 0.01 else 0)
  )
}


.schema23_target <- function(kind = "moments", mean = 5, variance = 2,
                             pmf_override = NULL) {
  stopifnot(kind %in% c("moments", "pmf"))
  pmf <- if (identical(kind, "pmf")) {
    if (is.null(pmf_override)) {
      out <- numeric(20L)
      out[c(3L, 5L, 7L)] <- c(0.25, 0.5, 0.25)
      out
    } else {
      unname(as.numeric(pmf_override))
    }
  } else {
    NULL
  }
  if (!is.null(pmf_override)) {
    moments <- .dpprior_target_pmf_moments(pmf)
    mean <- unname(moments[["mean"]])
    variance <- unname(moments[["variance"]])
  }
  snapshot <- .dpprior_new_snapshot(
    parameters = NULL, M = NULL,
    achieved = list(
      implied = list(mean = mean, variance = variance), interval = NULL,
      pmf = pmf
    ),
    residuals = list(mean = 0, variance = 0),
    tolerances = list(absolute = 1e-10),
    finite = TRUE, source = "target_constructor"
  )
  verification <- .dpprior_new_verification(
    method = "structural_target_check", performed = TRUE, passed = TRUE,
    reason = "verified", settings = list(tolerance = 1e-10),
    selected_snapshot = snapshot, verifier_snapshot = snapshot,
    stability = NULL,
    components = list(support = TRUE, moments = TRUE),
    invariants = list(identity = TRUE)
  )
  request <- if (identical(kind, "pmf")) {
    list(J = 20L, pmf = pmf)
  } else {
    list(J = 20L, mu_K = mean, var_K = variance)
  }
  canonical_request <- if (identical(kind, "pmf")) {
    list(J = 20L, pmf = pmf, interval = NULL)
  } else {
    c(request, list(interval = NULL, pmf = NULL))
  }
  target_method <- if (identical(kind, "pmf")) "target_pmf" else
    "target_moments"
  .dpprior_new_target_K(
    kind = kind, J = 20L, request = request,
    normalized = canonical_request, used = canonical_request,
    derivation = list(
      request_to_normalized = list(
        rule = if (identical(kind, "pmf")) {
          "validate_strict_pmf"
        } else {
          "canonicalize_direct_moments"
        },
        outcome = "canonicalized", opt_in = FALSE,
        before = request, after = canonical_request,
        evidence = list(
          source = "schema_fixture",
          explicit_null_fields = "interval"
        )
      ),
      normalized_to_used = NULL
    ),
    interval = NULL, family = NULL,
    assumptions = list(estimand = "K_J"), pmf = pmf,
    implied = list(mean = mean, variance = variance), achieved_interval = NULL,
    residuals = list(mean = 0, variance = 0),
    tolerances = list(absolute = 1e-10),
    status = "converged", usable = TRUE, verified = TRUE,
    parameters = NULL,
    computation = .schema23_computation(
      method = target_method, parameters = NULL
    ),
    verification = verification,
    provenance = .schema23_provenance(
      method = target_method, parameterization = "none"
    )
  )
}


.schema23_target_truth <- function(kind, constraint = 1e-9) {
  list(
    tolerances = list(
      constraint = constraint, pmf_l1 = .TOL_PMF_SUM,
      moment_relative = 1e-8, moment_scale_floor = 1
    ),
    controls = list(
      constraint_tolerance = constraint,
      root_tolerance = if (identical(kind, "interval")) {
        max(.Machine$double.eps, min(1e-12, constraint / 10))
      } else {
        NULL
      },
      max_iterations = if (identical(kind, "interval")) 1000L else NULL,
      pmf_l1_tolerance = .TOL_PMF_SUM,
      moment_relative_tolerance = 1e-8,
      moment_scale_floor = 1
    )
  )
}


.schema23_constructed_target_computation <- function(method, truth, status) {
  computation <- .schema23_computation(method = method, parameters = NULL)
  computation$request$controls <- truth$controls
  computation$used$controls <- truth$controls
  if (identical(status, "boundary")) {
    computation$termination <- .dpprior_new_termination(
      code = "boundary", source = "constructor", iterations = NULL,
      boundary_reason = "analytic_feasibility_boundary"
    )
  }
  computation
}


.schema23_constructed_target_verification <- function(
    kind, interval, tolerances, controls, pmf, verifier_pmf = pmf,
    implied = .dpprior_target_pmf_moments(pmf),
    verifier_implied = .dpprior_target_pmf_moments(verifier_pmf),
    residuals = list(), source = "target_constructor") {
  achieved_interval <- if (identical(kind, "interval")) {
    as.list(.dpprior_target_interval_masses(pmf, interval))
  } else {
    NULL
  }
  verifier_interval <- if (identical(kind, "interval")) {
    as.list(.dpprior_target_interval_masses(verifier_pmf, interval))
  } else {
    NULL
  }
  selected <- .dpprior_new_snapshot(
    parameters = NULL, M = NULL,
    achieved = list(
      implied = as.list(implied), interval = achieved_interval, pmf = pmf
    ),
    residuals = residuals, tolerances = tolerances, finite = TRUE,
    source = source
  )
  verifier <- .dpprior_new_snapshot(
    parameters = NULL, M = NULL,
    achieved = list(
      implied = as.list(verifier_implied), interval = verifier_interval,
      pmf = verifier_pmf
    ),
    residuals = residuals, tolerances = tolerances, finite = TRUE,
    source = "independent_target_reconstruction"
  )
  provisional <- list(
    kind = kind, interval = interval, tolerances = tolerances,
    verification = list(
      selected_snapshot = selected, verifier_snapshot = verifier
    )
  )
  stability <- .dpprior_expected_target_stability(provisional)
  list(
    achieved_interval = achieved_interval,
    verification = .dpprior_new_verification(
      method = "independent_target_reconstruction", performed = TRUE,
      passed = TRUE, reason = "verified",
      settings = list(
        pmf_l1_tolerance = tolerances$pmf_l1,
        moment_relative_tolerance = tolerances$moment_relative,
        moment_scale_floor = tolerances$moment_scale_floor,
        constraint_tolerance = tolerances$constraint,
        root_tolerance = controls$root_tolerance,
        max_iterations = controls$max_iterations
      ),
      selected_snapshot = selected, verifier_snapshot = verifier,
      stability = stability,
      components = list(
        target_reconstruction = .schema23_check(
          value = sum(abs(verifier_pmf - pmf)), reference = 0,
          tolerance = tolerances$pmf_l1,
          operator = "lte",
          source = "independent_target_reconstruction"
        ),
        order_stability = .schema23_check(
          value = stability$delta,
          reference = setNames(
            rep(0, length(stability$delta)), names(stability$delta)
          ),
          tolerance = stability$tolerance, operator = "lte",
          source = "independent_target_reconstruction"
        )
      ),
      invariants = list(
        support_identity = .schema23_check(
          value = c(
            target_support = TRUE, selected_length = TRUE,
            verifier_length = TRUE
          ),
          reference = c(
            target_support = TRUE, selected_length = TRUE,
            verifier_length = TRUE
          ), tolerance = NULL, operator = "identical",
          source = "independent_target_reconstruction"
        ),
        authority_identity = .schema23_check(
          value = c(
            selected_pmf = TRUE, used_pmf = TRUE,
            normalized_pmf_unidentified = TRUE, request_J = TRUE
          ),
          reference = c(
            selected_pmf = TRUE, used_pmf = TRUE,
            normalized_pmf_unidentified = TRUE, request_J = TRUE
          ), tolerance = NULL, operator = "identical",
          source = "independent_target_reconstruction"
        )
      )
    )
  )
}


.schema23_maxent_pmf <- function(J, groups, masses, tilt = NULL,
                                 boundary = NULL) {
  pmf <- numeric(J)
  for (group in names(groups)) {
    support <- groups[[group]]
    mass <- masses[[group]]
    if (mass == 0) next
    if (!is.null(boundary)) {
      point <- if (identical(boundary, "minimum")) min(support) else max(support)
      pmf[[point]] <- pmf[[point]] + mass
    } else {
      logits <- tilt * support
      weight <- exp(logits - max(logits))
      pmf[support] <- mass * weight / sum(weight)
    }
  }
  pmf
}


.schema23_interval_target <- function(
    type = c("equal_tail", "central_mass", "hard_bounds")) {
  type <- match.arg(type)
  if (identical(type, "hard_bounds")) {
    J <- 5L
    lower <- 2L
    upper <- 4L
    coverage <- 1
    mean_constraint <- 2
  } else {
    J <- 20L
    lower <- 3L
    upper <- 10L
    coverage <- 0.8
    mean_constraint <- 6.5
  }
  support <- seq_len(J)
  groups <- switch(
    type,
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
  masses <- switch(
    type,
    hard_bounds = c(inside = 1),
    equal_tail = c(
      left = (1 - coverage) / 2, inside = coverage,
      right = (1 - coverage) / 2
    ),
    central_mass = c(inside = coverage, outside = 1 - coverage)
  )
  lower_mean <- sum(vapply(
    names(groups), function(group) masses[[group]] * min(groups[[group]]),
    numeric(1)
  ))
  upper_mean <- sum(vapply(
    names(groups), function(group) masses[[group]] * max(groups[[group]]),
    numeric(1)
  ))
  feasibility_tolerance <- 64 * .Machine$double.eps * max(
    1, abs(lower_mean), abs(upper_mean), abs(mean_constraint)
  )
  boundary <- if (abs(mean_constraint - lower_mean) <= feasibility_tolerance) {
    "minimum"
  } else if (abs(mean_constraint - upper_mean) <= feasibility_tolerance) {
    "maximum"
  } else {
    NULL
  }
  tilt <- if (is.null(boundary)) {
    stats::uniroot(
      function(value) {
        pmf <- .schema23_maxent_pmf(J, groups, masses, tilt = value)
        sum(support * pmf) - mean_constraint
      },
      interval = c(-10, 10), tol = 1e-12, maxiter = 1000L
    )$root
  } else {
    NULL
  }
  pmf <- .schema23_maxent_pmf(
    J, groups, masses, tilt = tilt, boundary = boundary
  )
  interval <- list(
    lower = lower, upper = upper, type = type, coverage = coverage,
    family = "maxent", mu_K = mean_constraint,
    support = c(lower = 1L, upper = J), endpoints = "inclusive"
  )
  family <- list(
    name = "maxent",
    parameterization = "common exponential tilt within fixed-mass groups",
    explicit = TRUE, reference_measure = "counting measure on 1:J"
  )
  request <- list(J = J, K_interval = interval, mu_K = mean_constraint)
  normalized <- list(
    J = J, interval = interval, family = family, pmf = NULL
  )
  used <- list(J = J, interval = interval, family = family, pmf = pmf)
  truth <- .schema23_target_truth("interval")
  residuals <- list(
    mean = sum(support * pmf) - mean_constraint,
    interval = as.list(
      .dpprior_target_interval_masses(pmf, interval) -
        c(
          left_mass = sum(pmf[support < lower]),
          inside_mass = coverage,
          right_mass = sum(pmf[support > upper])
        )
    )
  )
  verification_parts <- .schema23_constructed_target_verification(
    "interval", interval, truth$tolerances, truth$controls, pmf,
    residuals = residuals
  )
  method <- "target_interval_maxent"
  status <- if (is.null(boundary)) "converged" else "boundary"
  .dpprior_new_target_K(
    kind = "interval", J = J, request = request,
    normalized = normalized, used = used,
    derivation = list(
      request_to_normalized = list(
        rule = "canonicalize_interval_request", outcome = "canonicalized",
        opt_in = FALSE, before = request, after = normalized,
        evidence = list(source = "schema_fixture")
      ),
      normalized_to_used = list(
        rule = switch(
          type,
          hard_bounds = "construct_maxent_hard_bounds_pmf",
          equal_tail = "construct_maxent_equal_tail_pmf",
          central_mass = "construct_maxent_central_mass_pmf"
        ),
        outcome = "derived", opt_in = FALSE,
        before = normalized, after = used,
        evidence = list(
          constructor_method = if (is.null(boundary)) "uniroot" else "boundary",
          group_masses = as.list(masses), common_tilt = tilt,
          boundary = list(active = !is.null(boundary), side = boundary),
          mean_constraint = mean_constraint,
          feasible_mean_lower = lower_mean,
          feasible_mean_upper = upper_mean,
          feasibility_tolerance = feasibility_tolerance,
          root_tolerance = truth$controls$root_tolerance
        )
      )
    ),
    interval = interval, family = family,
    assumptions = list(estimand = "K_J", maximum_entropy = TRUE),
    pmf = pmf, implied = as.list(.dpprior_target_pmf_moments(pmf)),
    achieved_interval = verification_parts$achieved_interval,
    residuals = residuals, tolerances = truth$tolerances,
    status = status, usable = TRUE, verified = TRUE,
    parameters = NULL,
    computation = .schema23_constructed_target_computation(
      method, truth, status
    ),
    verification = verification_parts$verification,
    provenance = .schema23_provenance(method, parameterization = "none")
  )
}


.schema23_family_target <- function() {
  J <- 20L
  requested_mean <- 5
  requested_variance <- 2
  df <- 2 * requested_mean^2 / requested_variance
  scale <- requested_variance / (2 * requested_mean)
  lower <- (seq_len(J) - 0.5) / scale
  upper <- (seq_len(J) + 0.5) / scale
  raw_mass <- pmax(
    stats::pchisq(upper, df = df) - stats::pchisq(lower, df = df), 0
  )
  retained_mass <- sum(raw_mass)
  pmf <- raw_mass / retained_mass
  family <- list(
    name = "scaled_chisq",
    parameterization = "df_and_scale_from_requested_moments",
    explicit = TRUE,
    reference_measure = paste(
      "continuity-corrected chi-square bins conditioned on 1:J"
    )
  )
  request <- list(
    J = J, mean = requested_mean, variance = requested_variance,
    family = family
  )
  normalized <- list(
    J = J, mean = requested_mean, variance = requested_variance,
    interval = NULL, family = family, pmf = NULL
  )
  used <- list(J = J, interval = NULL, family = family, pmf = pmf)
  truth <- .schema23_target_truth("family")
  residuals <- list(mean = 0, variance = 0)
  verification_parts <- .schema23_constructed_target_verification(
    "family", NULL, truth$tolerances, truth$controls, pmf,
    residuals = residuals
  )
  method <- "target_scaled_chisq"
  .dpprior_new_target_K(
    kind = "family", J = J, request = request,
    normalized = normalized, used = used,
    derivation = list(
      request_to_normalized = list(
        rule = "canonicalize_scaled_chisq_family_request",
        outcome = "canonicalized", opt_in = FALSE,
        before = request, after = normalized,
        evidence = list(source = "schema_fixture")
      ),
      normalized_to_used = list(
        rule = "construct_scaled_chisq_conditioned_pmf",
        outcome = "derived", opt_in = FALSE,
        before = normalized, after = used,
        evidence = list(
          df = df, scale = scale,
          binning = "continuity_corrected_half_integer_bins",
          retained_mass_before_normalization = retained_mass,
          omitted_mass = 1 - retained_mass,
          normalization = "explicit_support_conditioning",
          support = c(lower = 1L, upper = J)
        )
      )
    ),
    interval = NULL, family = family,
    assumptions = list(estimand = "K_J", support_conditioned = TRUE),
    pmf = pmf, implied = as.list(.dpprior_target_pmf_moments(pmf)),
    achieved_interval = NULL, residuals = residuals,
    tolerances = truth$tolerances,
    status = "converged", usable = TRUE, verified = TRUE,
    parameters = NULL,
    computation = .schema23_constructed_target_computation(
      method, truth, "converged"
    ),
    verification = verification_parts$verification,
    provenance = .schema23_provenance(method, parameterization = "none")
  )
}


.schema23_infeasible_interval_target <- function() {
  J <- 5L
  interval <- list(
    lower = 2L, upper = 4L, type = "hard_bounds", coverage = 1,
    family = "maxent", mu_K = 5,
    support = c(lower = 1L, upper = J), endpoints = "inclusive"
  )
  family <- list(
    name = "maxent",
    parameterization = "common exponential tilt within fixed-mass groups",
    explicit = TRUE, reference_measure = "counting measure on 1:J"
  )
  request <- list(J = J, K_interval = interval, mu_K = 5)
  normalized <- list(
    J = J, interval = interval, family = family, pmf = NULL
  )
  assumptions <- list(estimand = "K_J", construction = "maxent")
  lower_bound <- 2
  upper_bound <- 4
  feasibility_tolerance <- 64 * .Machine$double.eps * 5
  outside_distance <- 1
  certificate <- list(
    version = "1",
    method = "analytic_interval_group_support_feasibility",
    kind = "mean_outside_group_mass_hull",
    assumptions = assumptions, request = request, J = J,
    support = seq_len(J), interval = interval, family = family,
    group_masses = list(inside = 1),
    group_support_counts = list(inside = 3L), empty_groups = NULL,
    requested_mean = 5, feasible_mean_lower = lower_bound,
    feasible_mean_upper = upper_bound,
    feasibility_tolerance = feasibility_tolerance,
    side = "above", outside_distance = outside_distance,
    certified = TRUE, source = "analytic_maxent_interval_feasibility"
  )
  verification <- .dpprior_new_verification(
    method = "analytic_interval_infeasibility_certificate",
    performed = TRUE, passed = TRUE, reason = "certified_infeasible",
    settings = list(certificate = certificate),
    selected_snapshot = NULL, verifier_snapshot = NULL, stability = NULL,
    components = list(infeasibility_certificate = .schema23_check(
      value = outside_distance, reference = 0,
      tolerance = feasibility_tolerance, operator = "gt",
      source = "analytic_maxent_interval_feasibility"
    )),
    invariants = list(
      request_identity = .schema23_check(
        value = c(J = TRUE, interval = TRUE, mean = TRUE),
        reference = c(J = TRUE, interval = TRUE, mean = TRUE),
        tolerance = NULL, operator = "identical",
        source = "analytic_maxent_interval_feasibility"
      ),
      support_identity = .schema23_check(
        value = c(top = TRUE, interval = TRUE),
        reference = c(top = TRUE, interval = TRUE),
        tolerance = NULL, operator = "identical",
        source = "analytic_maxent_interval_feasibility"
      )
    )
  )
  truth <- .schema23_target_truth("interval")
  attempt <- .dpprior_new_attempt(
    id = "attempt-interval-feasibility", stage = "feasibility",
    method = "analytic_interval_group_support_feasibility",
    start = list(
      J = J, requested_mean = 5, interval_lower = 2L,
      interval_upper = 4L, coverage = 1
    ),
    bounds = list(
      feasible_mean_lower = lower_bound,
      feasible_mean_upper = upper_bound
    ),
    control = list(feasibility_tolerance = feasibility_tolerance),
    exit_code = 0L, message = "analytic interval infeasibility certified",
    iterations = 0L, evaluations = list(function_count = 1L),
    candidate_parameters = NULL, candidate_objective = NULL,
    elapsed_seconds = 0, warnings = character(), error = NULL,
    selected = FALSE, reason_code = "globally_infeasible_by_certificate",
    unavailable = c(
      candidate_parameters = "certificate route has no candidate parameters",
      candidate_objective = "certificate route has no candidate objective"
    )
  )
  computation <- .schema23_constructed_target_computation(
    "target_interval_maxent", truth, "converged"
  )
  computation$attempts <- list(attempt)
  computation$termination <- .dpprior_new_termination(
    code = "certified_infeasible", source = "analytic_certificate",
    iterations = NULL
  )
  .dpprior_new_target_K(
    kind = "interval", J = J, request = request,
    normalized = normalized, used = normalized,
    derivation = list(
      request_to_normalized = list(
        rule = "canonicalize_interval_request", outcome = "canonicalized",
        opt_in = FALSE, before = request, after = normalized,
        evidence = list(source = "schema_fixture")
      ),
      normalized_to_used = NULL
    ),
    interval = interval, family = family, assumptions = assumptions,
    pmf = NULL, implied = NULL, achieved_interval = NULL,
    residuals = list(
      unavailable_reason = "analytic interval construction infeasible"
    ),
    tolerances = truth$tolerances,
    status = "infeasible", usable = FALSE, verified = TRUE,
    parameters = NULL, computation = computation,
    verification = verification,
    provenance = .schema23_provenance(
      "target_interval_maxent", parameterization = "none"
    )
  )
}


.schema23_a1_projected_target <- function() {
  target <- .schema23_target()
  resolved <- suppressWarnings(.dpprior_a1_resolve_target(
    20L, 5, 2, projection = "nearest", signal_projection = FALSE
  ))
  target$used$var_K <- resolved$var_K_used
  target$implied$variance <- resolved$var_K_used
  for (snapshot_name in c("selected_snapshot", "verifier_snapshot")) {
    target$verification[[snapshot_name]]$achieved$implied$variance <-
      resolved$var_K_used
  }
  target$derivation$normalized_to_used <- list(
    rule = "project_a1_variance_to_nearest_interior",
    outcome = "projected", opt_in = TRUE,
    before = target$normalized, after = target$used,
    evidence = resolved$projection
  )
  target$provenance$projection <- list(
    applied = TRUE, opt_in = TRUE, policy = "nearest",
    record = list(
      before = target$normalized, after = target$used,
      authority = resolved$projection
    )
  )
  .dpprior_validate_target_v1(target)
  target
}


.schema23_weight_target <- function(relation = "target",
                                    metric = "wsb_mean", value = 0.4) {
  operator <- switch(relation, at_most = "<=", at_least = ">=", target = "target")
  threshold <- if (metric %in% c("wsb_tail", "wmax_tail_upper")) 0.5 else NULL
  probability <- if (identical(metric, "wsb_quantile")) 0.9 else NULL
  estimand <- switch(
    metric,
    wsb_tail = "P(W_SB > threshold)",
    wsb_mean = "E(W_SB)",
    wsb_quantile = "Q_probability(W_SB)",
    wmax_tail_upper = "certified upper bound for P(W_max > threshold)"
  )
  spec <- list(
    metric = metric, relation = relation, value = value,
    threshold = threshold, probability = probability
  )
  .dpprior_new_weight_target(
    request = spec, normalized = spec, used = spec,
    metric = metric, relation = relation, operator = operator, value = value,
    threshold = threshold, probability = probability, estimand = estimand,
    units = "probability",
    certification = if (identical(metric, "wmax_tail_upper")) {
      list(
        kind = "upper_bound", certified = TRUE, passed = TRUE,
        method = "certified_size_biased_mass_upper_bound",
        source = "wmax_tail_bounds"
      )
    } else {
      list()
    },
    provenance = list(source = "schema_fixture")
  )
}


.schema23_transformed_weight_target <- function(route = c("hard", "soft")) {
  route <- match.arg(route)
  if (identical(route, "hard")) {
    request <- list(
      metric = "wsb_tail", relation = "<=", bound = 0.4,
      threshold = 0.5, probability = NULL
    )
    normalized <- list(
      metric = "wsb_tail", relation = "at_most", value = 0.4,
      threshold = 0.5, probability = NULL
    )
    rule <- "canonicalize_hard_weight_target"
    evidence <- list(
      mode = "hard", value_field = "bound", relation_from = "<=",
      relation_to = "at_most", probability_field = "none"
    )
    metric <- "wsb_tail"
    relation <- "at_most"
    operator <- "<="
    threshold <- 0.5
    estimand <- "P(W_SB > threshold)"
  } else {
    request <- list(
      metric = "wsb_mean", relation = "target", value = 0.4,
      threshold = NULL, probability = NULL
    )
    normalized <- list(
      metric = "wsb_mean", relation = "target", value = 0.4,
      threshold = NULL, probability = NULL
    )
    rule <- "canonicalize_soft_weight_target"
    evidence <- list(
      mode = "soft", value_field = "value", relation_from = "target",
      relation_to = "target", probability_field = "none"
    )
    metric <- "wsb_mean"
    relation <- operator <- "target"
    threshold <- NULL
    estimand <- "E(W_SB)"
  }
  transformation <- if (identical(route, "hard")) {
    list(
      rule = rule, opt_in = FALSE, before = request, after = normalized,
      evidence = evidence
    )
  } else {
    NULL
  }
  .dpprior_new_weight_target(
    request = request, normalized = normalized, used = normalized,
    metric = metric, relation = relation, operator = operator, value = 0.4,
    threshold = threshold, probability = NULL, estimand = estimand,
    units = "probability", certification = list(),
    provenance = list(
      source = "schema_fixture",
      transformation = transformation,
      selection = NULL
    )
  )
}


.schema23_fit_parts <- function(method = "A2-MN",
                                M_selected = 80L,
                                M_verification = 160L) {
  parameters <- .dpprior_new_parameters(2, 3, "log_ab")
  achieved <- list(K = list(
    mean = 5, variance = 2, estimand = "K_J", source = "selected", M = M_selected
  ))
  residuals <- list(K = list(mean = 0, variance = 0))
  tolerances <- switch(
    method,
    `A2-MN` = list(
      K_adequacy = list(
        absolute = 1e-8, relative = 1e-8,
        scale_formula = "max(abs(target),1)"
      ),
      K_stability = list(absolute = 1e-10, relative = 1e-8, scale_floor = 1),
      step = 1e-10, boundary = 1e-6
    ),
    `A2-KL` = list(
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
    ),
    dual_anchor_hard_inequality = list(
      constraint = list(
        absolute = 1e-6, relative = 1e-6,
        effective = 1e-6 + 1e-6 * max(0.4, 1e-8)
      ),
      K = list(absolute = 1e-8, relative = 1e-6, scale_floor = 1e-8),
      weight = list(
        absolute = 1e-8, relative = 1e-6, scale_floor = 1e-8
      ),
      perturbation = list(
        absolute = 1e-4, relative = 1e-6, scale_floor = 1e-8,
        step = 1e-6, min_evaluations = 2L
      ),
      boundary = 1e-5,
      certificate = list(
        corner_uncertainty = 1e-6,
        rounding_floor = 64 * .Machine$double.eps
      )
    ),
    `dual-soft` = list(
      K = list(absolute = 1e-8, relative = 1e-6, scale_floor = 1),
      weight = list(absolute = 1e-8, relative = 1e-6, scale_floor = 1),
      objective = list(absolute = 1e-10, relative = 1e-8, scale_floor = 1),
      boundary = 1e-6,
      stationarity = list(tolerance = 1e-4, step = 1e-5),
      neighborhood = list(step = 1e-3), selection = 0
    ),
    list(K = list(absolute = 1e-6, relative = 1e-6))
  )
  selected <- .dpprior_new_snapshot(
    parameters, M_selected, achieved, residuals, tolerances, TRUE,
    "selected_order"
  )
  verifier_achieved <- achieved
  verifier_achieved$K$mean <- 5 + 1e-10
  verifier_achieved$K$M <- M_verification
  verifier_residuals <- residuals
  verifier_residuals$K$mean <- (5 + 1e-10) - 5
  verifier <- .dpprior_new_snapshot(
    parameters, M_verification, verifier_achieved, verifier_residuals,
    tolerances,
    TRUE, "independent_verifier"
  )
  stability_delta <- c(
    K.mean = abs(5 - (5 + 1e-10)),
    K.variance = 0
  )
  stability_tolerance <- if (identical(method, "A2-MN")) {
    c(
      K.mean = 1e-10 + 1e-8 * max(abs(5), abs(5 + 1e-10), 1),
      K.variance = 1e-10 + 1e-8 * max(abs(2), abs(2), 1)
    )
  } else if (identical(method, "dual_anchor_hard_inequality")) {
    c(
      K.mean = 1e-8 + 1e-6 * max(abs(5), abs(5 + 1e-10), 1e-8),
      K.variance = 1e-8 + 1e-6 * max(abs(2), abs(2), 1e-8)
    )
  } else {
    c(
      K.mean = 1e-8 + 1e-6 * max(abs(5), abs(5 + 1e-10), 1),
      K.variance = 1e-8 + 1e-6 * max(abs(2), abs(2), 1)
    )
  }
  stability_floor <- if (method %in% c("A2-MN", "A2-KL", "dual-soft")) {
    1
  } else {
    1e-8
  }
  list(
    parameters = parameters,
    achieved = achieved,
    residuals = residuals,
    tolerances = tolerances,
    computation = .schema23_computation(
      method, parameters, M_selected, M_verification
    ),
    verification = .dpprior_new_verification(
      method = "higher_order", performed = TRUE, passed = TRUE,
      reason = "all checks passed",
      settings = list(
        M_selected = M_selected, M_verification = M_verification
      ),
      selected_snapshot = selected, verifier_snapshot = verifier,
      stability = .dpprior_new_stability(
        delta = stability_delta,
        tolerance = stability_tolerance,
        formula = c(
          K.mean = "absolute_plus_relative_max",
          K.variance = "absolute_plus_relative_max"
        ),
        scale_floor = c(
          K.mean = stability_floor, K.variance = stability_floor
        ),
        source = "independent_verifier"
      ),
      components = list(adequacy = TRUE, order_stability = TRUE),
      invariants = list(identity = TRUE, finite = TRUE)
    )
  )
}


.schema23_legacy_fit <- function(lambda = 0.5,
                                 loss_type = "relative",
                                 max_iter = 50L,
                                 J = 20L,
                                 mu_K = 5,
                                 var_K = 8,
                                 weight_target = list(prob = list(
                                   threshold = 0.5, value = 0.3
                                 ))) {
  K_only <- DPprior_fit(
    J = J, mu_K = mu_K, var_K = var_K,
    method = "A2-MN", M = 80L, check_diagnostics = FALSE
  )
  suppressWarnings(DPprior_dual(
    K_only, weight_target, lambda = lambda, max_iter = max_iter,
    M = 80L, loss_type = loss_type
  ))
}


.schema23_migrated_dual <- function(allow_legacy = TRUE, verify = FALSE) {
  source <- structure(list(
    a = 1.27223464009856,
    b = 0.58075377814199,
    J = 20,
    target = list(
      mu_K = 5, var_K = 8, var_K_used = 8,
      confidence = NULL, type = "moments"
    ),
    method = "dual-anchor",
    status = "success",
    converged = TRUE,
    iterations = 7L,
    termination = "residual",
    fit = list(
      mu_K = 5.00000000001361,
      var_K = 7.99999999987637,
      residual = 1.24381729905704e-10
    ),
    solver_diagnostics = list(
      a0 = 4, b0 = 2.99573227355399, tol_F = 1e-8,
      tol_step = 1e-10, M = 80L, fallback_used = FALSE
    ),
    trace = data.frame(
      iteration = 1L, residual = 1.24381729905704e-10
    ),
    dual_anchor = list(
      w1_target = list(prob = list(threshold = 0.5, value = 0.3)),
      lambda = 1,
      loss_type = "relative",
      w1_achieved = list(
        mean = 0.42, prob_gt_50 = 0.368, prob_gt_90 = 0.13
      ),
      K_loss = 0,
      init = list(a = 1.27223464009856, b = 0.58075377814199),
      note = "lambda = 1 returns the K-only solution"
    )
  ), class = "DPprior_fit")
  suppressWarnings(upgrade_DPprior_object(
    source, verify = verify, allow_legacy = allow_legacy
  ))
}


.schema23_expect_legacy_reject <- function(x, mutate) {
  forged <- unserialize(serialize(x, NULL, xdr = TRUE))
  forged <- mutate(forged)
  expect_error(
    .dpprior_validate_result_v1(forged), class = "dpprior_schema_error"
  )
}


.schema23_set_migration_destination_version <- function(x, version) {
  x$provenance$backend$package_version <- version
  x$target$K$provenance$backend$package_version <- version
  x
}


.schema23_swap_migration_source_and_destination <- function(x) {
  x <- .schema23_set_migration_destination_version(x, "1.1.0")
  x$provenance$migration$source_schema <- "DPprior/2.0/fit"
  x$target$K$provenance$migration$source_schema <- "DPprior/2.0/fit"
  x$target$K$derivation$request_to_normalized$evidence$source_schema <-
    "DPprior/2.0/fit"
  x$target$weight$provenance$source_schema <- "DPprior/2.0/fit"
  x$compatibility$views$source$source_schema <- "DPprior/2.0/fit"
  x$compatibility$views$source$source_package_version <- "2.0.0"
  x
}


.schema23_fit <- function(mode = "a2_moment") {
  if (identical(mode, "dual_legacy")) {
    return(.schema23_legacy_fit())
  }
  method <- switch(
    mode,
    a1_proxy = "A1",
    a2_moment = "A2-MN",
    a2_kl = "A2-KL",
    dual_hard = "dual_anchor_hard_inequality",
    dual_soft = "dual-soft",
    dual_legacy = "dual-anchor"
  )
  parts <- .schema23_fit_parts(method)
  target <- list(K = .schema23_target())
  status <- "converged"
  usable <- verified <- TRUE
  approximation <- opt_in <- legacy <- FALSE
  input_fit <- NULL
  extension <- list()

  if (identical(mode, "a2_moment")) {
    adequacy_tolerance <- c(
      selected.mean = 1e-8 + 1e-8 * 5,
      selected.variance = 1e-8 + 1e-8 * 2,
      refined.mean = 1e-8 + 1e-8 * 5,
      refined.variance = 1e-8 + 1e-8 * 2,
      selected.standardized_norm = 1,
      refined.standardized_norm = 1
    )
    parts$verification$components <- list(
      residual_adequacy = .schema23_check(
        value = c(
          selected.mean = 0, selected.variance = 0,
          refined.mean = abs((5 + 1e-10) - 5), refined.variance = 0,
          selected.standardized_norm = 0,
          refined.standardized_norm = sqrt(mean(c(
            ((5 + 1e-10) - 5) /
              (1e-8 + 1e-8 * max(abs(5), 1)),
            0
          )^2))
        ),
        reference = c(
          selected.mean = 0, selected.variance = 0,
          refined.mean = 0, refined.variance = 0,
          selected.standardized_norm = 0,
          refined.standardized_norm = 0
        ),
        tolerance = adequacy_tolerance, operator = "lte"
      ),
      order_stability = .schema23_check(
        value = parts$verification$stability$delta,
        reference = c(K.mean = 0, K.variance = 0),
        tolerance = parts$verification$stability$tolerance,
        operator = "lte"
      ),
      parameter_identity = .schema23_check(
        value = c(a = 2, b = 3), reference = c(a = 2, b = 3),
        tolerance = NULL, operator = "identical"
      )
    )
  }

  if (identical(mode, "a2_kl")) {
    target$K <- .schema23_target("pmf")
    target_pmf <- target$K$pmf
    distribution_residuals <- list(
      kl = 0, l1 = 0, mean = 0, variance = 0
    )
    distribution_tolerances <- parts$tolerances$distribution
    parts$achieved$K$pmf <- target_pmf
    parts$residuals$distribution <- distribution_residuals
    parts$tolerances$distribution <- distribution_tolerances
    parts$verification$selected_snapshot$achieved <- parts$achieved
    parts$verification$selected_snapshot$residuals <- parts$residuals
    parts$verification$selected_snapshot$tolerances <- parts$tolerances
    parts$verification$verifier_snapshot$achieved <- parts$achieved
    parts$verification$verifier_snapshot$residuals <- parts$residuals
    parts$verification$verifier_snapshot$achieved$K$M <-
      parts$verification$verifier_snapshot$M
    parts$verification$verifier_snapshot$tolerances <- parts$tolerances
    parts$verification$stability <- .dpprior_new_stability(
      delta = c(pmf.l1 = 0),
      tolerance = c(pmf.l1 = 1e-10 + 1e-8),
      formula = c(pmf.l1 = "direct_pmf_l1_tolerance"),
      scale_floor = c(pmf.l1 = 0),
      source = "independent_verifier"
    )
    adequacy_value <- c(
      selected.kl = 0, selected.l1 = 0,
      selected.mean_scaled = 0, selected.variance_scaled = 0,
      refined.kl = 0, refined.l1 = 0,
      refined.mean_scaled = 0, refined.variance_scaled = 0
    )
    adequacy_tolerance <- c(
      selected.kl = 0.015, selected.l1 = 0.11,
      selected.mean_scaled = 0.01, selected.variance_scaled = 0.065,
      refined.kl = 0.015, refined.l1 = 0.11,
      refined.mean_scaled = 0.01, refined.variance_scaled = 0.065
    )
    parts$verification$components <- list(
      target_identity = .schema23_check(
        value = target_pmf, reference = target_pmf,
        tolerance = NULL, operator = "identical"
      ),
      pmf_adequacy = .schema23_check(
        value = adequacy_value,
        reference = setNames(rep(0, length(adequacy_value)),
                             names(adequacy_value)),
        tolerance = adequacy_tolerance, operator = "lte"
      ),
      order_stability = .schema23_check(
        value = c(pmf.l1 = 0), reference = c(pmf.l1 = 0),
        tolerance = c(pmf.l1 = 1e-10 + 1e-8),
        operator = "lte"
      ),
      candidate_selection = .schema23_check(
        tolerance = 0, operator = "lte"
      )
    )
  }

  if (identical(mode, "a1_proxy")) {
    status <- "approximate"
    verified <- FALSE
    approximation <- opt_in <- TRUE
    parts$verification <- .dpprior_new_verification(
      method = "exact_estimand_not_performed", performed = FALSE,
      passed = FALSE, reason = "A1 is a proxy",
      selected_snapshot = parts$verification$selected_snapshot,
      verifier_snapshot = NULL, stability = NULL,
      components = list(exact_estimand = FALSE),
      invariants = list(identity = TRUE)
    )
    extension <- list(proxy = .dpprior_new_proxy(
      mapping = list(formula = "closed_form"),
      mapping_verification = list(passed = TRUE),
      projection = list(applied = FALSE),
      caveats = "proxy mapping is not exact estimand verification"
    ))
  }

  if (mode %in% c("dual_hard", "dual_soft", "dual_legacy")) {
    relation <- if (identical(mode, "dual_soft")) "target" else "at_most"
    soft_selected <- soft_refined <- NULL
    if (mode %in% c("dual_hard", "dual_soft")) {
      provisional_fit <- list(J = 20L)
      metric_spec <- list(
        metric = "wsb_mean", estimand = "E(W_SB)", units = "probability",
        threshold = NULL, probability = NULL
      )
      soft_selected <- list(
        K = exact_K_moments(20L, 2, 3, 80L),
        weight = if (identical(mode, "dual_soft")) {
          .dpprior_v2_eval_metric(metric_spec, 2, 3, J = 20L, M = 80L)
        } else NULL
      )
      soft_refined <- list(
        K = exact_K_moments(20L, 2, 3, 160L),
        weight = if (identical(mode, "dual_soft")) {
          .dpprior_v2_eval_metric(metric_spec, 2, 3, J = 20L, M = 160L)
        } else NULL
      )
      target$K <- .schema23_target(
        mean = soft_selected$K$mean, variance = soft_selected$K$var
      )
      parts$achieved$K$mean <- soft_selected$K$mean
      parts$achieved$K$variance <- soft_selected$K$var
      parts$residuals$K <- list(mean = 0, variance = 0)
      parts$verification$selected_snapshot$achieved <- parts$achieved
      parts$verification$selected_snapshot$residuals <- parts$residuals
      parts$verification$verifier_snapshot$achieved$K$mean <-
        soft_refined$K$mean
      parts$verification$verifier_snapshot$achieved$K$variance <-
        soft_refined$K$var
      parts$verification$verifier_snapshot$residuals$K <- list(
        mean = soft_refined$K$mean - soft_selected$K$mean,
        variance = soft_refined$K$var - soft_selected$K$var
      )
      parts$verification$verifier_snapshot$residuals$weight <-
        list(raw = 0, directed = 0)
      parts$verification$verifier_snapshot$residuals$weight <-
        list(raw = 0, directed = 0)
    }
    target$weight <- .schema23_weight_target(
      relation = relation,
      value = if (identical(mode, "dual_soft")) {
        soft_selected$weight$value
      } else {
        0.4
      }
    )
    parts$achieved$weight <- list(
      metric = "wsb_mean",
      value = if (identical(mode, "dual_soft")) {
        soft_selected$weight$value
      } else {
        0.4
      },
      source = "selected"
    )
    parts$residuals$weight <- list(raw = 0, directed = 0)
    parts$verification$selected_snapshot$achieved <- parts$achieved
    parts$verification$selected_snapshot$residuals <- parts$residuals
    parts$verification$selected_snapshot$tolerances <- parts$tolerances
    parts$verification$verifier_snapshot$achieved <- parts$achieved
    if (identical(mode, "dual_soft")) {
      parts$verification$verifier_snapshot$achieved$K$mean <-
        soft_refined$K$mean
      parts$verification$verifier_snapshot$achieved$K$variance <-
        soft_refined$K$var
      parts$verification$verifier_snapshot$achieved$weight$value <-
        soft_refined$weight$value
      parts$verification$verifier_snapshot$residuals$K <- list(
        mean = soft_refined$K$mean - soft_selected$K$mean,
        variance = soft_refined$K$var - soft_selected$K$var
      )
      parts$verification$verifier_snapshot$residuals$weight <- list(
        raw = soft_refined$weight$value - soft_selected$weight$value,
        directed = soft_refined$weight$value - soft_selected$weight$value
      )
    } else if (identical(mode, "dual_hard")) {
      parts$verification$verifier_snapshot$achieved$K$mean <-
        soft_refined$K$mean
      parts$verification$verifier_snapshot$achieved$K$variance <-
        soft_refined$K$var
      parts$verification$verifier_snapshot$achieved$K$M <-
        parts$verification$verifier_snapshot$M
      parts$verification$verifier_snapshot$residuals$K <- list(
        mean = soft_refined$K$mean - soft_selected$K$mean,
        variance = soft_refined$K$var - soft_selected$K$var
      )
      parts$verification$verifier_snapshot$residuals$weight <-
        list(raw = 0, directed = 0)
    } else {
      parts$verification$verifier_snapshot$achieved$K$mean <- 5 + 1e-10
      parts$verification$verifier_snapshot$achieved$K$M <-
        parts$verification$verifier_snapshot$M
    }
    if (identical(mode, "dual_legacy")) {
      parts$verification$verifier_snapshot$residuals <- parts$residuals
    }
    parts$verification$verifier_snapshot$tolerances <- parts$tolerances
    parts$verification$stability <- .dpprior_new_stability(
      delta = c(parts$verification$stability$delta, weight.value = 0),
      tolerance = c(
        parts$verification$stability$tolerance,
        weight.value = if (identical(mode, "dual_hard")) {
          1e-8 + 1e-6 * max(0.4, 0.4, 1e-8)
        } else {
          1e-8 + 1e-6 * max(
            parts$verification$selected_snapshot$achieved$weight$value,
            parts$verification$verifier_snapshot$achieved$weight$value, 1
          )
        }
      ),
      formula = c(
        K.mean = "absolute_plus_relative_max",
        K.variance = "absolute_plus_relative_max",
        weight.value = "absolute_plus_relative_max"
      ),
      scale_floor = setNames(
        rep(if (identical(mode, "dual_soft")) 1 else 1e-8, 3L),
        c("K.mean", "K.variance", "weight.value")
      ),
      source = "independent_verifier"
    )
    if (identical(mode, "dual_soft")) {
      parts$verification$stability <- .dpprior_new_stability(
        delta = c(
          K.mean = abs(soft_selected$K$mean - soft_refined$K$mean),
          K.variance = abs(soft_selected$K$var - soft_refined$K$var),
          weight.value = abs(
            soft_selected$weight$value - soft_refined$weight$value
          )
        ),
        tolerance = c(
          K.mean = 1e-8 + 1e-6 * max(
            abs(soft_selected$K$mean), abs(soft_refined$K$mean), 1
          ),
          K.variance = 1e-8 + 1e-6 * max(
            abs(soft_selected$K$var), abs(soft_refined$K$var), 1
          ),
          weight.value = 1e-8 + 1e-6 * max(
            abs(soft_selected$weight$value),
            abs(soft_refined$weight$value), 1
          )
        ),
        formula = setNames(
          rep("absolute_plus_relative_max", 3L),
          c("K.mean", "K.variance", "weight.value")
        ),
        scale_floor = setNames(
          rep(1, 3L), c("K.mean", "K.variance", "weight.value")
        ),
        source = "independent_verifier"
      )
    } else if (identical(mode, "dual_hard")) {
      hard_delta <- c(
        K.mean = abs(soft_selected$K$mean - soft_refined$K$mean),
        K.variance = abs(soft_selected$K$var - soft_refined$K$var),
        weight.value = 0
      )
      parts$verification$stability <- .dpprior_new_stability(
        delta = hard_delta,
        tolerance = c(
          K.mean = 1e-8 + 1e-6 * max(
            abs(soft_selected$K$mean), abs(soft_refined$K$mean), 1e-8
          ),
          K.variance = 1e-8 + 1e-6 * max(
            abs(soft_selected$K$var), abs(soft_refined$K$var), 1e-8
          ),
          weight.value = 1e-8 + 1e-6 * max(0.4, 1e-8)
        ),
        formula = setNames(
          rep("absolute_plus_relative_max", 3L), names(hard_delta)
        ),
        scale_floor = setNames(rep(1e-8, 3L), names(hard_delta)),
        source = "independent_verifier"
      )
    }
    parts$achieved$K$M <- parts$verification$selected_snapshot$M
    parts$verification$selected_snapshot$achieved$K$M <-
      parts$verification$selected_snapshot$M
    parts$verification$verifier_snapshot$achieved$K$M <-
      parts$verification$verifier_snapshot$M
    input_snapshot_reference <- function(snapshot) {
      list(
        parameters = snapshot$parameters, M = snapshot$M,
        achieved_K = snapshot$achieved$K, finite = snapshot$finite,
        source = snapshot$source
      )
    }
    input_fit <- list(
      schema = "dpprior.result/1", mode = "a2_moment", method = "A2-MN",
      J = 20L, status = "converged", usable = TRUE, verified = TRUE,
      parameters = parts$parameters,
      target = list(
        schema = target$K$schema, kind = target$K$kind, J = target$K$J,
        used = target$K$used, implied = target$K$implied
      ),
      decision_evidence = NULL,
      selected_snapshot = input_snapshot_reference(
        parts$verification$selected_snapshot
      ),
      verifier_snapshot = input_snapshot_reference(
        parts$verification$verifier_snapshot
      )
    )
    if (identical(mode, "dual_legacy")) input_fit <- NULL
  }

  if (identical(mode, "dual_hard")) {
    parts$verification$components <- list(
      constraint_selected = .schema23_check(
        tolerance = parts$tolerances$constraint$effective, operator = "lte"
      ),
      constraint_refined = .schema23_check(
        tolerance = parts$tolerances$constraint$effective, operator = "lte"
      ),
      order_stability = .schema23_check(
        value = 0, reference = 0,
        tolerance = parts$verification$stability$tolerance[["weight.value"]],
        operator = "lte"
      ),
      metric_certification = .schema23_check(),
      candidate_selection = .schema23_check(
        tolerance = 1e-8, operator = "lte"
      ),
      perturbation = .schema23_check(
        value = TRUE, reference = TRUE, tolerance = NULL,
        operator = "identical"
      )
    )
    extension <- list(constraint = .dpprior_new_constraint(
      relation = "at_most", operator = "<=", residual = 0, slack = 0,
      tolerance = parts$tolerances$constraint,
      satisfied = TRUE, active = TRUE,
      feasibility = list(
        classification = "feasible_candidate", certified_infeasible = FALSE,
        feasibility_unknown = FALSE, certificate = NULL,
        candidate_count = 1L, verified_candidate_count = 1L,
        feasible_candidate_count = 1L
      ),
        optimality = .schema23_hard_optimality(scales = list(
          mean = max(abs(target$K$implied$mean), 1),
          variance = max(abs(target$K$implied$variance), 1)
        ))
    ))
    hard_scales <- list(K = list(
      mean = max(abs(target$K$implied$mean), 1),
      variance = max(abs(target$K$implied$variance), 1)
    ))
    parts$computation$scaling <- .dpprior_new_scaling(
      requested = hard_scales, used = hard_scales,
      formula = "fixed_from_input_target_max_abs_one",
      values = hard_scales, fixed_from_input = TRUE
    )
  }

  if (identical(mode, "dual_soft")) {
    scales <- list(
      K = list(
        mean = max(abs(target$K$implied$mean), 1),
        variance = max(abs(target$K$implied$variance), 1)
      ),
      weight = 1
    )
    parts$computation$scaling <- .dpprior_new_scaling(
      requested = scales, used = scales,
      formula = "fixed_scaled_squared_loss", values = scales,
      fixed_from_input = TRUE
    )
    provisional <- list(
      J = 20L, parameters = parts$parameters, target = target,
      tolerances = parts$tolerances, computation = parts$computation,
      provenance = list(input_fit = input_fit),
      tradeoff = list(lambda = 0.5, scales = scales)
    )
    bounds <- parts$computation$attempts[[1L]]$bounds
    soft_truth <- .dpprior_expected_soft_optimality(
      provisional, bounds$lower, bounds$upper
    )
    objective_policy <- parts$tolerances$objective
    objective_tolerance_for <- function(left, right) {
      objective_policy$absolute + objective_policy$relative * max(
        abs(left), abs(right), objective_policy$scale_floor
      )
    }
    objective_tolerance <- objective_tolerance_for(
      soft_truth$selected_objective, soft_truth$selected_objective
    )
    neighbor_tolerances <- vapply(
      soft_truth$neighbor_objectives,
      function(value) objective_tolerance_for(
        soft_truth$refined_objective, value
      ), numeric(1)
    )
    start_tolerance <- objective_tolerance_for(
      soft_truth$refined_objective, soft_truth$start_objective
    )
    stationarity_components <- vapply(
      seq_along(soft_truth$gradient), function(index) switch(
        soft_truth$stationarity_operator[[index]],
        abs_lte = abs(soft_truth$gradient[[index]]) <=
          parts$tolerances$stationarity$tolerance,
        gte = soft_truth$gradient[[index]] >=
          -parts$tolerances$stationarity$tolerance,
        lte = soft_truth$gradient[[index]] <=
          parts$tolerances$stationarity$tolerance
      ), logical(1)
    )
    names(stationarity_components) <- names(soft_truth$gradient)
    local_pass <- all(
      soft_truth$refined_objective <=
        soft_truth$neighbor_objectives + neighbor_tolerances
    )
    start_pass <- soft_truth$refined_objective <=
      soft_truth$start_objective + start_tolerance
    optimality <- list(
      performed = TRUE,
      passed = all(stationarity_components) && local_pass && start_pass,
      recorded_objective = soft_truth$selected_objective,
      recomputed_objective = soft_truth$selected_objective,
      objective_difference = 0, objective_tolerance = objective_tolerance,
      objective_passed = TRUE, gradient = soft_truth$gradient,
      gradient_method = soft_truth$gradient_method,
      bound_state = soft_truth$bound_state,
      boundary_tolerance = parts$tolerances$boundary,
      stationarity_operator = soft_truth$stationarity_operator,
      stationarity_tolerance = parts$tolerances$stationarity$tolerance,
      component_pass = stationarity_components,
      stationarity_passed = all(stationarity_components),
      local_base_objective = soft_truth$refined_objective,
      neighbor_objectives = soft_truth$neighbor_objectives,
      neighbor_tolerances = neighbor_tolerances,
      local_minimum_passed = local_pass,
      start_objective = soft_truth$start_objective,
      candidate_objective = soft_truth$selected_objective,
      start_tolerance = start_tolerance, no_worse_start = start_pass,
      selection_tolerance = 0,
      source = "independent_refined_objective_verification",
      unavailable_reason = NULL
    )
    parts$computation$attempts[[1L]]$candidate_objective <-
      soft_truth$selected_objective
    parts$verification$components <- list(
      objective_recomputation = .schema23_check(
        value = 0, reference = 0, tolerance = objective_tolerance,
        operator = "abs_lte"
      ),
      order_stability = .schema23_check(),
      local_optimality = .schema23_check(
        value = stationarity_components,
        reference = setNames(
          rep(TRUE, length(stationarity_components)),
          names(stationarity_components)
        ), tolerance = NULL, operator = "identical"
      ),
      candidate_selection = .schema23_check(
        tolerance = 0, operator = "lte"
      )
    )
    extension <- list(tradeoff = .dpprior_new_tradeoff(
      lambda = 0.5, K_loss = 0, weight_loss = 0, total_loss = 0,
      target_residual = 0, directed_residual = 0, scales = scales,
      endpoint = FALSE, optimality = optimality
    ))
  }

  if (identical(mode, "dual_legacy")) {
    status <- "approximate"
    verified <- FALSE
    approximation <- opt_in <- legacy <- TRUE
    parts$verification <- .dpprior_new_verification(
      method = "legacy_unverified", performed = FALSE, passed = FALSE,
      reason = "legacy approximation",
      selected_snapshot = parts$verification$selected_snapshot,
      components = list(exact_estimand = FALSE),
      invariants = list(identity = TRUE)
    )
    extension <- list(legacy = .dpprior_new_legacy_details(
      contract = "DPprior_dual_v1", lambda = 0.5,
      losses = list(K = 0.2, weight = 0.1, total = 0.15),
      approximation_opt_in = TRUE,
      warning_code = "legacy_dual_approximation"
    ))
  }

  if (mode %in% c("a2_moment", "a2_kl", "dual_hard", "dual_soft")) {
    orders <- parts$computation$orders
    controls <- parts$computation$used$controls
    contract <- switch(
      mode,
      a2_moment = list(
        method = "independent_higher_order_moment_recomputation",
        settings = list(
          M_selected = orders$M_selected,
          M_verification_required = orders$M_verification_required,
          M_verification = orders$M_verification_used
        ),
        invariants = c("finite_parameters", "K_support")
      ),
      a2_kl = list(
        method = "fresh higher-order marginal PMF and direct-moment audit",
        settings = list(
          M_selected = orders$M_selected,
          M_verification = orders$M_verification_used,
          M_verification_required = orders$M_verification_required,
          pmf_abs_tol = parts$tolerances$distribution$order$pmf_absolute,
          pmf_rel_tol = parts$tolerances$distribution$order$pmf_relative
        ),
        invariants = c(
          "finite_parameters", "parameter_identity", "pmf_probability",
          "support_identity"
        )
      ),
      dual_hard = list(
        method = "fresh_higher_order_recomputation_and_local_perturbation",
        settings = list(
          M_fit = orders$M_selected,
          M_verify = orders$M_verification_used,
          M_verify_required = orders$M_verification_required,
          log_bounds = c(-15, 15),
          verification_abs_tol = controls$verification_abs_tol,
          verification_rel_tol = controls$verification_rel_tol,
          perturbation_step = controls$perturbation_step,
          perturbation_abs_tol = controls$perturbation_abs_tol
        ),
        invariants = c(
          "probability", "K_support", "finite_parameters_inside_domain"
        )
      ),
      dual_soft = list(
        method = "fresh higher-order K moments and named weight metric",
        settings = list(
          M_fit = orders$M_selected,
          M_verify = orders$M_verification_used,
          log_bounds = controls$log_bounds,
          verification_abs_tol = controls$verification_abs_tol,
          verification_rel_tol = controls$verification_rel_tol,
          objective_abs_tol = controls$objective_abs_tol,
          objective_rel_tol = controls$objective_rel_tol,
          boundary_tol = controls$boundary_tol,
          stationarity_step = controls$stationarity_step,
          stationarity_tol = controls$stationarity_tol,
          neighborhood_step = controls$local_neighbor_step
        ),
        invariants = c(
          "probability", "K_support", "finite_parameters_inside_domain",
          "fixed_input_scales"
        )
      )
    )
    parts$verification$method <- contract$method
    parts$verification$settings <- contract$settings
    parts$verification$components <- lapply(
      parts$verification$components,
      function(check) {
        check$source <- "independent_verifier"
        check
      }
    )
    parts$verification$invariants <- setNames(
      lapply(contract$invariants, function(name) {
        .schema23_check(
          value = TRUE, reference = TRUE, tolerance = NULL,
          operator = "identical", source = "independent_verifier"
        )
      }),
      contract$invariants
    )
  }

  if (mode %in% c("a2_moment", "a2_kl", "dual_hard", "dual_soft")) {
    candidate_source <- "candidate:candidate-1"
    candidate_checks <- switch(
      mode,
      a2_moment = list(
        candidate_finite = .schema23_check(
          value = c(parameters = TRUE, snapshot = TRUE, objective = TRUE),
          reference = c(parameters = TRUE, snapshot = TRUE, objective = TRUE),
          tolerance = NULL, operator = "identical", source = candidate_source
        )
      ),
      a2_kl = list(
        candidate_distribution = .schema23_check(
          value = c(pmf_mass_error = 0, pmf_minimum_violation = 0),
          reference = c(pmf_mass_error = 0, pmf_minimum_violation = 0),
          tolerance = c(
            pmf_mass_error = .TOL_PMF_SUM, pmf_minimum_violation = 0
          ),
          operator = "lte", source = candidate_source
        )
      ),
      dual_hard = list(
        constraint_selected = .schema23_check(
          tolerance = parts$tolerances$constraint$effective,
          operator = "lte", source = candidate_source
        ),
        constraint_refined = .schema23_check(
          tolerance = parts$tolerances$constraint$effective,
          operator = "lte", source = candidate_source
        ),
        order_stability = .schema23_check(
          value = parts$verification$stability$delta,
          reference = setNames(
            rep(0, length(parts$verification$stability$delta)),
            names(parts$verification$stability$delta)
          ),
          tolerance = parts$verification$stability$tolerance,
          operator = "lte", source = candidate_source
        ),
        metric_certification = .schema23_check(
          value = c(
            selected_metric = TRUE, refined_metric = TRUE,
            selected_finite = TRUE, refined_finite = TRUE
          ),
          reference = c(
            selected_metric = TRUE, refined_metric = TRUE,
            selected_finite = TRUE, refined_finite = TRUE
          ), tolerance = NULL,
          operator = "identical", source = candidate_source
        ),
        perturbation = .schema23_check(
          value = c(maximum_metric_delta = 0), reference = 0,
          tolerance = c(
            maximum_metric_delta = 1e-4 + 1e-6 * max(0.4, 1e-8)
          ), operator = "lte",
          source = candidate_source
        ),
        invariants = .schema23_check(
          value = c(
            selected_finite = TRUE, refined_finite = TRUE,
            parameter_identity = TRUE, support_valid = TRUE
          ),
          reference = c(
            selected_finite = TRUE, refined_finite = TRUE,
            parameter_identity = TRUE, support_valid = TRUE
          ), tolerance = NULL,
          operator = "identical", source = candidate_source
        )
      ),
      dual_soft = list(
        candidate_domain = .schema23_check(
          value = c(
            K_support = TRUE, K_variance = TRUE, weight_support = TRUE
          ),
          reference = c(
            K_support = TRUE, K_variance = TRUE, weight_support = TRUE
          ), tolerance = NULL,
          operator = "identical", source = candidate_source
        )
      )
    )
    candidate_recorded_objective <- if (identical(mode, "dual_soft")) {
      parts$computation$attempts[[1L]]$candidate_objective
    } else {
      0
    }
    candidate_fresh_objective <- if (identical(mode, "dual_soft")) {
      soft_truth$selected_objective
    } else {
      0
    }
    candidate <- .dpprior_new_candidate_evaluation(
      id = "candidate-1", attempt_id = "attempt-1", method = method,
      generator = "direct_attempt", parameters = parts$parameters,
      objective_kind = switch(
        mode, a2_moment = "standardized_residual", a2_kl = "kl",
        dual_hard = "K_loss", dual_soft = "soft_tradeoff"
      ),
      recorded_objective = candidate_recorded_objective,
      fresh_objective = candidate_fresh_objective,
      selection_objective = candidate_fresh_objective,
      objective_tolerance = switch(
        mode, a2_moment = 0, a2_kl = 0,
        dual_hard = 1e-8, dual_soft = 1e-10 + 1e-8
      ),
      selected_snapshot = parts$verification$selected_snapshot,
      verifier_snapshot = if (identical(mode, "dual_hard")) {
        parts$verification$verifier_snapshot
      } else {
        NULL
      },
      checks = candidate_checks, execution_success = TRUE,
      optimizer_supported = TRUE, selected = TRUE,
      source = "schema_fixture_candidate_ledger"
    )
    parts$computation$candidate_evaluations <- list(candidate)
    parts$computation$selected_candidate_id <- "candidate-1"
    .dpprior_validate_computation(parts$computation)
  }

  .dpprior_new_fit(
    mode = mode, method = method, J = 20L, status = status,
    usable = usable, verified = verified,
    message = if (verified) "" else "explicit approximation",
    parameters = parts$parameters, target = target,
    achieved = parts$achieved, residuals = parts$residuals,
    tolerances = parts$tolerances, computation = parts$computation,
    verification = parts$verification,
    provenance = .schema23_provenance(
      method = method, approximation = approximation, opt_in = opt_in,
      legacy = legacy, input_fit = input_fit
    ),
    extension = extension
  )
}


.schema23_input_K_record <- function(J, parameters, M, pmf = NULL) {
  moments <- if (is.null(pmf)) {
    fresh <- exact_K_moments(J, parameters$a, parameters$b, M = M)
    c(mean = fresh$mean, variance = fresh$var)
  } else {
    .dpprior_target_pmf_moments(pmf)
  }
  out <- list(
    mean = unname(moments[["mean"]]),
    variance = unname(moments[["variance"]]),
    estimand = "K_J", source = "selected", M = as.integer(M)
  )
  if (!is.null(pmf)) out$pmf <- unname(as.numeric(pmf))
  out
}


.schema23_with_A2_MN_input_orders <- function(
    fit, M_selected = 80L, M_verification = 200L) {
  input_fit <- fit$provenance$input_fit
  input_fit$selected_snapshot$M <- as.integer(M_selected)
  input_fit$selected_snapshot$achieved_K <- .schema23_input_K_record(
    fit$J, input_fit$parameters, M_selected
  )
  input_fit$verifier_snapshot$M <- as.integer(M_verification)
  input_fit$verifier_snapshot$achieved_K <- .schema23_input_K_record(
    fit$J, input_fit$parameters, M_verification
  )
  fit$provenance$input_fit <- input_fit
  .dpprior_validate_result_v1(fit)
  fit
}


.schema23_with_A2_KL_input_fit <- function(
    fit, M_selected = 80L, M_verification = 200L, validate = TRUE) {
  parameters <- fit$provenance$input_fit$parameters
  distribution <- .schema23_fit_parts("A2-KL")$tolerances$distribution
  fresh <- .get_K_pmf_support(
    fit$J, parameters$a, parameters$b,
    M = M_selected, M_verify = M_verification,
    abs_tol = distribution$order$pmf_absolute,
    rel_tol = distribution$order$pmf_relative
  )
  selected_pmf <- unname(as.numeric(fresh$pmf))
  verifier_pmf <- unname(as.numeric(fresh$verification_pmf))
  original_implied <- fit$target$K$implied
  target_K <- .schema23_target(
    "pmf", pmf_override = selected_pmf
  )
  # The PMF-derived and exact-moment paths agree within the fixed numerical
  # identity tolerance.  Preserve the dual target's exact public moment record
  # while retaining the objective PMF as the A2-KL decision authority.
  target_K$implied <- original_implied
  target_K$verification$selected_snapshot$achieved$implied <- original_implied
  target_K$verification$verifier_snapshot$achieved$implied <- original_implied
  fit$target$K <- target_K
  compact_target <- list(
    schema = target_K$schema, kind = target_K$kind, J = target_K$J,
    used = target_K$used, implied = target_K$implied
  )
  snapshot <- function(M, pmf, source) {
    list(
      parameters = parameters, M = as.integer(M),
      achieved_K = .schema23_input_K_record(fit$J, parameters, M, pmf),
      finite = TRUE, source = source
    )
  }
  fit$provenance$input_fit <- list(
    schema = "dpprior.result/1", mode = "a2_kl", method = "A2-KL",
    J = fit$J, status = "converged", usable = TRUE, verified = TRUE,
    parameters = parameters, target = compact_target,
    decision_evidence = list(
      target_K = target_K, distribution_tolerances = distribution
    ),
    selected_snapshot = snapshot(M_selected, selected_pmf, "selected_order"),
    verifier_snapshot = snapshot(
      M_verification, verifier_pmf, "independent_verifier"
    )
  )
  if (isTRUE(validate)) .dpprior_validate_result_v1(fit)
  fit
}


.schema23_no_candidate <- function(mode, status = "failed") {
  stopifnot(mode %in% c("dual_hard", "dual_soft"))
  fit <- .schema23_fit(mode)
  fit$status <- status
  fit$usable <- FALSE
  fit$verified <- identical(status, "infeasible")
  fit$message <- if (identical(status, "infeasible")) {
    "certified infeasible"
  } else {
    "no finite candidate"
  }
  fit["parameters"] <- list(NULL)
  fit$achieved <- list()
  fit$residuals <- list()
  fit$computation$attempts <- list()
  fit$computation$candidate_evaluations <- list()
  fit$computation["selected_candidate_id"] <- list(NULL)
  fit$computation$orders <- .dpprior_new_orders(
    M_requested = NULL, M_selected = NULL,
    M_verification_required = NULL, M_verification_used = NULL,
    requested_reason = "not_applicable",
    selected_reason = "no_candidate",
    verification_required_reason = "certificate_only",
    verification_used_reason = "certificate_only"
  )
  fit$computation["selected_attempt_id"] <- list(NULL)
  fit$computation$termination$code <- if (identical(status, "infeasible")) {
    "certified_infeasible"
  } else {
    "no_candidate"
  }
  fit$computation$termination$source <- if (identical(status, "infeasible")) {
    "analytic_certificate"
  } else {
    "no_candidate"
  }
  fit$computation$termination["iterations"] <- list(NULL)

  if (identical(status, "infeasible")) {
    fit$target$weight <- .schema23_weight_target(
      relation = "at_most", metric = "wsb_tail"
    )
    log_bounds <- c(-2, -1)
    fit$computation$request$controls$log_bounds <- log_bounds
    fit$computation$used$controls$log_bounds <- log_bounds
    metric_spec <- list(
      metric = fit$target$weight$metric,
      estimand = fit$target$weight$estimand,
      units = fit$target$weight$units,
      threshold = fit$target$weight$threshold,
      probability = fit$target$weight$probability
    )
    certificate_evaluation <- function(parameters, M) {
      evaluated <- .dpprior_v2_eval_metric(
        metric_spec, parameters$a, parameters$b, J = fit$J, M = M
      )
      list(
        metric = fit$target$weight$metric, value = evaluated$value,
        method = evaluated$method, error_bound = evaluated$error_bound,
        certified = evaluated$certification$certified, M = as.integer(M),
        source = "independent_metric_evaluator"
      )
    }
    certificate_corner <- function(eta) {
      parameters <- list(a = exp(eta[[1L]]), b = exp(eta[[2L]]))
      selected <- certificate_evaluation(parameters, 80L)
      refined <- certificate_evaluation(parameters, 160L)
      uncertainty <- max(
        abs(selected$value - refined$value), selected$error_bound,
        refined$error_bound, 64 * .Machine$double.eps
      )
      list(
        eta = c(log_a = eta[[1L]], log_b = eta[[2L]]),
        parameters = parameters, selected = selected, refined = refined,
        finite = TRUE, uncertainty = uncertainty,
        lower = max(0, refined$value - uncertainty),
        upper = min(1, refined$value + uncertainty)
      )
    }
    minimum_corner <- certificate_corner(c(log_bounds[[2L]], log_bounds[[1L]]))
    maximum_corner <- certificate_corner(c(log_bounds[[1L]], log_bounds[[2L]]))
    certificate <- list(
      method = paste0(
        "analytic_global_monotonicity_with_refined_order_",
        "corner_enclosures"
      ),
      version = "1", J = fit$J, metric = "wsb_tail", relation = "at_most",
      target_value = 0.4, threshold = fit$target$weight$threshold,
      probability = fit$target$weight$probability, support = c(0, 1),
      domain = list(
        log_a = log_bounds, log_b = log_bounds,
        a = exp(log_bounds), b = exp(log_bounds)
      ),
      monotonicity = list(a = "non_increasing", b = "non_decreasing"),
      M_selected = 80L, M_verification = 160L,
      effective_tolerance = fit$tolerances$constraint$effective,
      minimum = minimum_corner, maximum = maximum_corner,
      lower_bound = minimum_corner$lower,
      upper_bound = maximum_corner$upper,
      tolerance = fit$tolerances$constraint$effective,
      certified = TRUE, source = "phase8_metric_extrema_corner_enclosures"
    )
    fit$computation$attempts <- list(.dpprior_new_attempt(
      id = "attempt-feasibility-probe", stage = "feasibility",
      method = "analytic_monotonicity_feasibility_probe", start = NULL,
      bounds = list(log_a = log_bounds, log_b = log_bounds),
      control = list(M = 80L, M_verify = 160L), exit_code = 0L,
      message = "certified_infeasible", iterations = 0L,
      evaluations = list(function_count = 2L),
      candidate_parameters = NULL,
      candidate_objective = minimum_corner$refined$value,
      elapsed_seconds = 0.001, warnings = character(), error = NULL,
      selected = FALSE, reason_code = "globally_infeasible_by_certificate",
      unavailable = c(
        start = "analytic global corner proof has no optimizer start",
        candidate_parameters = "certificate route has no public candidate"
      )
    ))
    fit$verification <- .dpprior_new_verification(
      method = "analytic_global_monotonicity_corner_certificate",
      performed = TRUE,
      passed = TRUE, reason = "certified",
      components = list(infeasibility_certificate = .schema23_check(
        value = certificate$lower_bound,
        reference = certificate$target_value,
        tolerance = certificate$tolerance, operator = "gt",
        source = "analytic_certificate"
      )),
      invariants = list(domain_monotonicity = .schema23_check(
        value = TRUE, reference = TRUE, tolerance = NULL,
        operator = "identical", source = "analytic_certificate"
      ))
    )
  } else {
    fit$verification <- .dpprior_new_verification(
      method = "no_candidate", performed = FALSE, passed = FALSE,
      reason = "all optimizer attempts failed",
      components = list(),
      invariants = list(no_public_candidate = .schema23_check(
        value = TRUE, reference = TRUE, tolerance = NULL,
        operator = "identical", source = "independent_verifier"
      ))
    )
  }

  if (identical(mode, "dual_hard")) {
    fit$constraint[c(
      "residual", "slack", "tolerance", "satisfied", "active"
    )] <- rep(list(NULL), 5L)
    fit$constraint$feasibility <- if (identical(status, "infeasible")) {
      list(
        classification = "certified_infeasible",
        certified_infeasible = TRUE,
        feasibility_unknown = FALSE,
        certificate = certificate,
        candidate_count = 0L,
        verified_candidate_count = 0L,
        feasible_candidate_count = 0L
      )
    } else {
      list(
        classification = "unknown",
        certified_infeasible = FALSE,
        feasibility_unknown = TRUE,
        certificate = NULL,
        candidate_count = 0L,
        verified_candidate_count = 0L,
        feasible_candidate_count = 0L
      )
    }
    fit$constraint$optimality <- .schema23_hard_optimality(FALSE)
  } else {
    fit$tradeoff[c(
      "K_loss", "weight_loss", "total_loss", "target_residual",
      "directed_residual"
    )] <- rep(list(NULL), 5L)
    fit$tradeoff$optimality <- .schema23_soft_optimality(FALSE)
  }
  .dpprior_validate_result_v1(fit)
  fit
}


.schema23_soft_failed_endpoint <- function(
    message = "soft path returned no canonical candidate") {
  fit <- .schema23_no_candidate("dual_soft", "failed")
  fit$message <- message
  fit$tradeoff$lambda <- 1
  fit$tradeoff$endpoint <- TRUE
  .dpprior_validate_result_v1(fit)
  fit
}


.schema23_hard_failed_rejected_candidates <- function() {
  template <- .schema23_fit("dual_hard")
  fit <- .schema23_no_candidate("dual_hard", "failed")
  fit$computation$scaling <- .dpprior_new_scaling(
    requested = list(K = list(mean = 5, variance = 2)),
    used = list(K = list(mean = 5, variance = 2)),
    formula = "fixed_from_input_target_max_abs_one",
    values = list(K = list(mean = 5, variance = 2)),
    fixed_from_input = TRUE
  )

  make_candidate <- function(index, parameters, method, stage) {
    attempt <- template$computation$attempts[[1L]]
    attempt$id <- paste0("attempt-", index)
    attempt$method <- method
    attempt$stage <- stage
    attempt$candidate_parameters <- parameters
    attempt$candidate_objective <- 0
    attempt$selected <- FALSE
    attempt$reason_code <- "order_stability_failed"

    selected <- template$verification$selected_snapshot
    verifier <- template$verification$verifier_snapshot
    selected$parameters <- parameters
    verifier$parameters <- parameters
    selected$tolerances <- fit$tolerances
    verifier$tolerances <- fit$tolerances
    selected$achieved$weight$value <- 0.4000009
    verifier$achieved$weight$value <- 0.4000011
    verifier$achieved$K$mean <- 5.1
    verifier$residuals$K$mean <- 0.1
    selected_residual <- selected$achieved$weight$value - 0.4
    refined_residual <- verifier$achieved$weight$value - 0.4
    selected$residuals$weight <- list(
      raw = selected_residual, directed = selected_residual
    )
    verifier$residuals$weight <- list(
      raw = refined_residual, directed = refined_residual
    )

    source <- paste0("candidate:candidate-", index)
    stability_delta <- c(
      K.mean = abs(
        selected$achieved$K$mean - verifier$achieved$K$mean
      ),
      K.variance = abs(
        selected$achieved$K$variance - verifier$achieved$K$variance
      ),
      weight.value = abs(
        selected$achieved$weight$value - verifier$achieved$weight$value
      )
    )
    stability_tolerance <- c(
      K.mean = 1e-8 + 1e-6 * max(
        abs(selected$achieved$K$mean), abs(verifier$achieved$K$mean), 1e-8
      ),
      K.variance = 1e-8 + 1e-6 * max(
        abs(selected$achieved$K$variance),
        abs(verifier$achieved$K$variance), 1e-8
      ),
      weight.value = fit$tolerances$weight$absolute +
        fit$tolerances$weight$relative * max(
          abs(selected$achieved$weight$value),
          abs(verifier$achieved$weight$value),
          fit$tolerances$weight$scale_floor
        )
    )
    checks <- list(
      constraint_selected = .schema23_check(
        value = selected_residual, reference = 0,
        tolerance = fit$tolerances$constraint$effective,
        operator = "lte", source = source
      ),
      constraint_refined = .schema23_check(
        value = refined_residual, reference = 0,
        tolerance = fit$tolerances$constraint$effective,
        operator = "lte", source = source
      ),
      order_stability = .schema23_check(
        value = stability_delta,
        reference = setNames(rep(0, 3L), names(stability_delta)),
        tolerance = stability_tolerance, operator = "lte", source = source
      ),
      metric_certification = .schema23_check(
        value = c(
          selected_metric = TRUE, refined_metric = TRUE,
          selected_finite = TRUE, refined_finite = TRUE
        ),
        reference = c(
          selected_metric = TRUE, refined_metric = TRUE,
          selected_finite = TRUE, refined_finite = TRUE
        ),
        tolerance = NULL, operator = "identical", source = source
      ),
      perturbation = .schema23_check(
        value = c(maximum_metric_delta = 0), reference = 0,
        tolerance = c(
          maximum_metric_delta = fit$tolerances$perturbation$absolute +
            fit$tolerances$perturbation$relative * max(
              abs(selected$achieved$weight$value),
              fit$tolerances$perturbation$scale_floor
            )
        ), operator = "lte",
        source = source
      ),
      invariants = .schema23_check(
        value = c(
          selected_finite = TRUE, refined_finite = TRUE,
          parameter_identity = TRUE, support_valid = TRUE
        ),
        reference = c(
          selected_finite = TRUE, refined_finite = TRUE,
          parameter_identity = TRUE, support_valid = TRUE
        ),
        tolerance = NULL, operator = "identical", source = source
      )
    )
    candidate <- .dpprior_new_candidate_evaluation(
      id = paste0("candidate-", index), attempt_id = attempt$id,
      method = method, generator = "direct_attempt",
      parameters = parameters, objective_kind = "K_loss",
      recorded_objective = 0, fresh_objective = 0,
      selection_objective = 0, objective_tolerance = 1e-8,
      selected_snapshot = selected, verifier_snapshot = verifier,
      checks = checks, execution_success = TRUE,
      optimizer_supported = TRUE, selected = FALSE,
      source = "retained_rejected_candidate"
    )
    list(attempt = attempt, candidate = candidate)
  }

  first <- make_candidate(
    1L, .dpprior_new_parameters(2, 3, "log_ab"),
    "L-BFGS-B", "primary"
  )
  second <- make_candidate(
    2L, .dpprior_new_parameters(2.1, 3, "log_ab"),
    "K_only_L-BFGS-B", "optimizer"
  )
  fit$computation$attempts <- list(first$attempt, second$attempt)
  fit$computation$candidate_evaluations <- list(
    first$candidate, second$candidate
  )
  fit$computation$orders <- template$computation$orders
  fit$constraint$feasibility$candidate_count <- 2L
  fit$constraint$feasibility$verified_candidate_count <- 2L
  fit$constraint$feasibility$feasible_candidate_count <- 0L
  .dpprior_validate_result_v1(fit)
  fit
}


.schema23_hard_approximate_scan <- function() {
  fit <- .schema23_fit("dual_hard")
  fit$status <- "approximate"
  fit$usable <- FALSE
  fit$verified <- FALSE
  fit$message <- "scientific candidate lacks ordinary optimizer convergence"
  fit$verification$passed <- FALSE
  fit$computation$attempts <- list()
  fit$computation["selected_attempt_id"] <- list(NULL)
  old <- fit$computation$candidate_evaluations[[1L]]
  fit$computation$candidate_evaluations <- list(
    .dpprior_new_candidate_evaluation(
      id = "candidate-1", attempt_id = NULL,
      method = "deterministic_profile_scan",
      generator = "deterministic_profile_scan",
      parameters = fit$parameters, objective_kind = "K_loss",
      recorded_objective = NULL,
      recorded_objective_reason = "generated candidate has no raw objective",
      fresh_objective = 0, selection_objective = 0,
      objective_tolerance = 1e-12,
      selected_snapshot = fit$verification$selected_snapshot,
      verifier_snapshot = fit$verification$verifier_snapshot,
      checks = old$checks, execution_success = TRUE,
      optimizer_supported = FALSE, selected = TRUE,
      source = "deterministic_profile_scan_candidate"
    )
  )
  fit$computation$termination$code <- "approximate"
  fit$computation$termination$source <- "candidate_evaluation"
  fit$computation$termination["iterations"] <- list(NULL)
  .dpprior_validate_result_v1(fit)
  fit
}


.schema23_hard_unsatisfied_diagnostic <- function() {
  fit <- .schema23_fit("dual_hard")
  fit$status <- "approximate"
  fit$usable <- FALSE
  fit$verified <- FALSE
  fit$message <- "finite hard candidate retained as an unsatisfied diagnostic"
  fit$verification$passed <- FALSE
  fit$verification$reason <- "two-order hard constraint is not satisfied"

  selected_weight <- 0.5
  refined_weight <- 0.5000001
  selected_residual <- selected_weight - fit$target$weight$value
  refined_residual <- refined_weight - fit$target$weight$value
  fit$achieved$weight$value <- selected_weight
  fit$residuals$weight <- list(
    raw = selected_residual, directed = selected_residual
  )
  fit$verification$selected_snapshot$achieved <- fit$achieved
  fit$verification$selected_snapshot$residuals <- fit$residuals
  fit$verification$verifier_snapshot$achieved$weight$value <- refined_weight
  fit$verification$verifier_snapshot$residuals$weight <- list(
    raw = refined_residual, directed = refined_residual
  )
  expected_stability <- .dpprior_expected_result_stability(fit)
  fit$verification$stability <- .dpprior_new_stability(
    delta = expected_stability$delta,
    tolerance = expected_stability$tolerance,
    formula = expected_stability$formula,
    scale_floor = expected_stability$scale_floor,
    source = "independent_verifier"
  )
  fit$constraint$residual <- selected_residual
  fit$constraint$slack <- -selected_residual
  fit$constraint$satisfied <- FALSE
  fit$constraint$active <- FALSE
  fit$constraint$feasibility <- list(
    classification = "unknown", certified_infeasible = FALSE,
    feasibility_unknown = TRUE, certificate = NULL,
    candidate_count = 1L, verified_candidate_count = 1L,
    feasible_candidate_count = 0L
  )
  fit$constraint$optimality$passed <- TRUE
  fit$constraint$optimality$selected_K_loss <- 0
  fit$constraint$optimality$minimum_K_loss <- 0

  fit$verification$components$constraint_selected <- .schema23_check(
    value = selected_residual, reference = 0,
    tolerance = fit$tolerances$constraint$effective,
    operator = "lte", source = "independent_verifier"
  )
  fit$verification$components$constraint_refined <- .schema23_check(
    value = refined_residual, reference = 0,
    tolerance = fit$tolerances$constraint$effective,
    operator = "lte", source = "independent_verifier"
  )
  fit$verification$components$order_stability <- .schema23_check(
    value = expected_stability$delta[["weight.value"]], reference = 0,
    tolerance = expected_stability$tolerance[["weight.value"]],
    operator = "lte", source = "independent_verifier"
  )

  source <- "candidate:candidate-1"
  candidate_checks <- fit$computation$candidate_evaluations[[1L]]$checks
  candidate_checks$constraint_selected <- .schema23_check(
    value = selected_residual, reference = 0,
    tolerance = fit$tolerances$constraint$effective,
    operator = "lte", source = source
  )
  candidate_checks$constraint_refined <- .schema23_check(
    value = refined_residual, reference = 0,
    tolerance = fit$tolerances$constraint$effective,
    operator = "lte", source = source
  )
  candidate_checks$order_stability <- .schema23_check(
    value = expected_stability$delta,
    reference = setNames(rep(0, length(expected_stability$delta)),
                         names(expected_stability$delta)),
    tolerance = expected_stability$tolerance,
    operator = "lte", source = source
  )
  candidate_checks$perturbation <- .schema23_check(
    value = c(maximum_metric_delta = 0), reference = 0,
    tolerance = c(
      maximum_metric_delta = fit$tolerances$perturbation$absolute +
        fit$tolerances$perturbation$relative * max(
          selected_weight, fit$tolerances$perturbation$scale_floor
        )
    ), operator = "lte", source = source
  )
  fit$computation$candidate_evaluations[[1L]] <-
    .dpprior_new_candidate_evaluation(
      id = "candidate-1", attempt_id = "attempt-1",
      method = "dual_anchor_hard_inequality", generator = "direct_attempt",
      parameters = fit$parameters, objective_kind = "K_loss",
      recorded_objective = 0, fresh_objective = 0,
      selection_objective = 0, objective_tolerance = 1e-8,
      selected_snapshot = fit$verification$selected_snapshot,
      verifier_snapshot = fit$verification$verifier_snapshot,
      checks = candidate_checks, execution_success = TRUE,
      optimizer_supported = TRUE, diagnostic_eligible = TRUE,
      selected = TRUE, source = "signed_hard_diagnostic_fixture"
    )
  fit$computation$termination <- .dpprior_new_termination(
    code = "approximate", message = fit$message,
    source = "candidate_evaluation", iterations = 5L
  )
  .dpprior_validate_result_v1(fit)
  fit
}


.schema23_hard_feasible_beats_diagnostic <- function() {
  fit <- .schema23_fit("dual_hard")
  target_mean <- fit$target$K$implied$mean
  mean_scale <- fit$constraint$optimality$K_scales$mean
  selected_mean <- target_mean + mean_scale
  fit$achieved$K$mean <- selected_mean
  fit$residuals$K$mean <- selected_mean - target_mean
  fit$verification$selected_snapshot$achieved <- fit$achieved
  fit$verification$selected_snapshot$residuals <- fit$residuals
  fit$verification$verifier_snapshot$achieved$K$mean <- selected_mean + 1e-10
  fit$verification$verifier_snapshot$residuals$K$mean <-
    selected_mean + 1e-10 - target_mean
  expected_stability <- .dpprior_expected_result_stability(fit)
  fit$verification$stability <- .dpprior_new_stability(
    delta = expected_stability$delta,
    tolerance = expected_stability$tolerance,
    formula = expected_stability$formula,
    scale_floor = expected_stability$scale_floor,
    source = "independent_verifier"
  )
  fit$verification$components$order_stability <- .schema23_check(
    value = expected_stability$delta[["weight.value"]], reference = 0,
    tolerance = expected_stability$tolerance[["weight.value"]],
    operator = "lte", source = "independent_verifier"
  )
  fit$constraint$optimality$selected_K_loss <- 1
  fit$constraint$optimality$minimum_K_loss <- 1
  fit$computation$attempts[[1L]]$candidate_objective <- 1

  selected_source <- "candidate:candidate-1"
  selected_checks <- fit$computation$candidate_evaluations[[1L]]$checks
  selected_checks$order_stability <- .schema23_check(
    value = expected_stability$delta,
    reference = setNames(rep(0, length(expected_stability$delta)),
                         names(expected_stability$delta)),
    tolerance = expected_stability$tolerance,
    operator = "lte", source = selected_source
  )
  fit$computation$candidate_evaluations[[1L]] <-
    .dpprior_new_candidate_evaluation(
      id = "candidate-1", attempt_id = "attempt-1",
      method = "dual_anchor_hard_inequality", generator = "direct_attempt",
      parameters = fit$parameters, objective_kind = "K_loss",
      recorded_objective = 1, fresh_objective = 1,
      selection_objective = 1, objective_tolerance = 1e-8,
      selected_snapshot = fit$verification$selected_snapshot,
      verifier_snapshot = fit$verification$verifier_snapshot,
      checks = selected_checks, execution_success = TRUE,
      optimizer_supported = TRUE, selected = TRUE,
      source = "feasible_priority_fixture"
    )

  diagnostic <- .schema23_hard_unsatisfied_diagnostic()
  diagnostic_attempt <- diagnostic$computation$attempts[[1L]]
  diagnostic_attempt$id <- "attempt-diagnostic"
  diagnostic_attempt$candidate_parameters <- .dpprior_new_parameters(
    2.1, 3, "log_ab"
  )
  diagnostic_attempt$selected <- FALSE
  diagnostic_attempt$reason_code <- "constraint_verification_failed"
  diagnostic_candidate <-
    diagnostic$computation$candidate_evaluations[[1L]]
  diagnostic_candidate$id <- "candidate-diagnostic"
  diagnostic_candidate$attempt_id <- diagnostic_attempt$id
  diagnostic_candidate$parameters <- diagnostic_attempt$candidate_parameters
  diagnostic_candidate$selected_snapshot$parameters <-
    diagnostic_attempt$candidate_parameters
  diagnostic_candidate$verifier_snapshot$parameters <-
    diagnostic_attempt$candidate_parameters
  diagnostic_candidate$selected <- FALSE
  diagnostic_candidate$outcome <- "rejected"
  names(diagnostic_candidate$checks) <- names(
    diagnostic$computation$candidate_evaluations[[1L]]$checks
  )
  for (check_name in names(diagnostic_candidate$checks)) {
    diagnostic_candidate$checks[[check_name]]$source <-
      "candidate:candidate-diagnostic"
  }
  .dpprior_validate_candidate_evaluation(diagnostic_candidate)
  fit$computation$attempts <- c(
    fit$computation$attempts, list(diagnostic_attempt)
  )
  fit$computation$candidate_evaluations <- c(
    fit$computation$candidate_evaluations, list(diagnostic_candidate)
  )
  fit$constraint$feasibility$candidate_count <- 2L
  fit$constraint$feasibility$verified_candidate_count <- 2L
  fit$constraint$feasibility$feasible_candidate_count <- 1L
  .dpprior_validate_result_v1(fit)
  fit
}


.schema23_a2_kl_approximate_initializer <- function() {
  fit <- .schema23_fit("a2_kl")
  fit$status <- "approximate"
  fit$usable <- FALSE
  fit$verified <- FALSE
  fit$message <- "initializer retained without optimizer convergence"
  fit$verification$passed <- FALSE
  attempt <- fit$computation$attempts[[1L]]
  attempt$stage <- "initialization"
  attempt$method <- "A2-MN"
  attempt["exit_code"] <- list(NULL)
  attempt["candidate_objective"] <- list(NULL)
  attempt$unavailable <- c(
    exit_code = "initializer has no optimizer exit",
    candidate_objective = "initializer has no recorded KL objective"
  )
  fit$computation$attempts[[1L]] <- attempt
  old <- fit$computation$candidate_evaluations[[1L]]
  fit$computation$candidate_evaluations <- list(
    .dpprior_new_candidate_evaluation(
      id = "candidate-1", attempt_id = "attempt-1", method = "A2-MN",
      generator = "initialization", parameters = fit$parameters,
      objective_kind = "kl", recorded_objective = NULL,
      recorded_objective_reason = "initializer did not record KL",
      fresh_objective = 0, selection_objective = 0,
      objective_tolerance = 1e-12,
      selected_snapshot = fit$verification$selected_snapshot,
      checks = old$checks, execution_success = TRUE,
      optimizer_supported = FALSE, selected = TRUE,
      source = "A2_MN_initialization_candidate"
    )
  )
  fit$computation$termination$code <- "approximate"
  fit$computation$termination$source <- "candidate_evaluation"
  .dpprior_validate_result_v1(fit)
  fit
}


.schema23_soft_approximate_diagnostic <- function() {
  fit <- .schema23_fit("dual_soft")
  fit$status <- "approximate"
  fit$usable <- FALSE
  fit$verified <- FALSE
  fit$message <- "selected diagnostic candidate failed local optimality"
  fit$verification$passed <- FALSE
  fit$computation$request$controls$stationarity_tol <- 1e-12
  fit$computation$used$controls$stationarity_tol <- 1e-12
  fit$tolerances$stationarity$tolerance <- 1e-12
  fit$verification$selected_snapshot$tolerances <- fit$tolerances
  fit$verification$verifier_snapshot$tolerances <- fit$tolerances
  fit$computation$candidate_evaluations[[1L]]$selected_snapshot$tolerances <-
    fit$tolerances
  fit$verification$settings$stationarity_tol <- 1e-12
  fit$tradeoff$optimality$stationarity_tolerance <- 1e-12
  fit$tradeoff$optimality$component_pass <-
    abs(fit$tradeoff$optimality$gradient) <= 1e-12
  fit$tradeoff$optimality$stationarity_passed <- FALSE
  fit$tradeoff$optimality$passed <- FALSE
  fit$verification$components$local_optimality <- .schema23_check(
    value = fit$tradeoff$optimality$component_pass,
    reference = c(log_shape = TRUE, log_rate = TRUE),
    tolerance = NULL, operator = "identical",
    source = "independent_verifier"
  )
  fit$computation$termination$code <- "approximate"
  fit$computation$termination$source <- "candidate_evaluation"
  .dpprior_validate_result_v1(fit)
  fit
}


.schema23_soft_endpoint <- function() {
  fit <- .schema23_fit("dual_soft")
  fit$tradeoff$lambda <- 1
  fit$tradeoff$endpoint <- TRUE
  fit$tradeoff$optimality <- .schema23_soft_optimality(FALSE)
  fit$computation$attempts <- list()
  fit$computation$candidate_evaluations <- list()
  fit$computation["selected_attempt_id"] <- list(NULL)
  fit$computation["selected_candidate_id"] <- list(NULL)
  fit$computation$termination <- .dpprior_new_termination(
    code = "endpoint", source = "endpoint", iterations = NULL
  )
  input_fit <- fit$provenance$input_fit
  target_K <- fit$target$K
  expected_input_target <- list(
    schema = target_K$schema, kind = target_K$kind, J = target_K$J,
    used = target_K$used, implied = target_K$implied
  )
  snapshot_reference <- function(snapshot) {
    list(
      parameters = snapshot$parameters, M = snapshot$M,
      achieved_K = snapshot$achieved$K, finite = snapshot$finite,
      source = snapshot$source
    )
  }
  endpoint_identity <- c(
    schema = identical(input_fit$schema, "dpprior.result/1"),
    mode = input_fit$mode %in% c("a2_moment", "a2_kl", "dual_hard", "dual_soft"),
    method = input_fit$method %in% .DPPRIOR_MODE_METHODS[[input_fit$mode]],
    J = identical(input_fit$J, fit$J),
    status = input_fit$status %in% c("converged", "boundary"),
    usable = isTRUE(input_fit$usable), verified = isTRUE(input_fit$verified),
    parameters = identical(input_fit$parameters, fit$parameters),
    target = identical(input_fit$target, expected_input_target),
    decision_evidence = is.null(input_fit$decision_evidence),
    selected_snapshot = identical(
      input_fit$selected_snapshot,
      snapshot_reference(fit$verification$selected_snapshot)
    ),
    verifier_snapshot = identical(
      input_fit$verifier_snapshot,
      snapshot_reference(fit$verification$verifier_snapshot)
    )
  )
  fit$verification$components <- list(
    endpoint_input_identity = .schema23_check(
      value = endpoint_identity,
      reference = setNames(rep(TRUE, length(endpoint_identity)),
                           names(endpoint_identity)),
      tolerance = NULL, operator = "identical",
      source = "independent_verifier"
    ),
    order_stability = .schema23_check(
      value = fit$verification$stability$delta,
      reference = setNames(
        rep(0, length(fit$verification$stability$delta)),
        names(fit$verification$stability$delta)
      ),
      tolerance = fit$verification$stability$tolerance, operator = "lte",
      source = "independent_verifier"
    )
  )
  .dpprior_validate_result_v1(fit)
  fit
}


.schema23_diagnostics <- function() {
  parameters <- .dpprior_new_parameters(2, 3, "log_ab")
  absolute_tolerance <- 1e-10
  relative_tolerance <- 1e-8
  pmf_mass_tolerance <- .TOL_PMF_SUM
  K_evidence <- .get_K_pmf_support(
    20L, parameters$a, parameters$b, M = 80L, M_verify = 160L,
    abs_tol = absolute_tolerance, rel_tol = relative_tolerance
  )
  selected_pmf <- unname(as.numeric(K_evidence$pmf))
  verifier_pmf <- unname(as.numeric(K_evidence$verification_pmf))
  selected_K <- .dpprior_target_pmf_moments(selected_pmf)
  verifier_K <- .dpprior_target_pmf_moments(verifier_pmf)
  selected_weight <- as.numeric(mean_w1(parameters$a, parameters$b, 80L))
  verifier_weight <- as.numeric(mean_w1(parameters$a, parameters$b, 160L))
  selected_rho <- c(
    mean = as.numeric(mean_rho(parameters$a, parameters$b, 80L)),
    variance = as.numeric(var_rho(parameters$a, parameters$b, 80L))
  )
  verifier_rho <- c(
    mean = as.numeric(mean_rho(parameters$a, parameters$b, 160L)),
    variance = as.numeric(var_rho(parameters$a, parameters$b, 160L))
  )
  refinement_delta <- c(
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
  refinement_tolerance <- c(
    K.mean = absolute_tolerance + relative_tolerance * max(
      abs(selected_K[["mean"]]), abs(verifier_K[["mean"]]), 1
    ),
    K.variance = absolute_tolerance + relative_tolerance * max(
      abs(selected_K[["variance"]]), abs(verifier_K[["variance"]]), 1
    ),
    K.pmf_l1 = absolute_tolerance + relative_tolerance,
    weights.mean = absolute_tolerance + relative_tolerance * max(
      abs(selected_weight), abs(verifier_weight), 1
    ),
    coclustering.mean = absolute_tolerance + relative_tolerance * max(
      abs(selected_rho[["mean"]]), abs(verifier_rho[["mean"]]), 1
    ),
    coclustering.variance = absolute_tolerance + relative_tolerance * max(
      abs(selected_rho[["variance"]]), abs(verifier_rho[["variance"]]), 1
    )
  )
  diagnostics <- list(
    policy_results = list(), warnings = character(),
    alpha = list(
      status = "converged", usable = TRUE, verified = TRUE,
      mean = 2 / 3, CV = 1 / sqrt(2)
    ),
    K = list(
      status = "converged", usable = TRUE, verified = TRUE,
      mean = unname(selected_K[["mean"]]),
      variance = unname(selected_K[["variance"]]),
      pmf = selected_pmf, M = 80L
    ),
    weights = list(
      status = "converged", usable = TRUE, verified = TRUE,
      mean = selected_weight
    ),
    coclustering = list(
      status = "converged", usable = TRUE, verified = TRUE,
      mean = selected_rho[["mean"]], variance = selected_rho[["variance"]]
    )
  )
  achieved <- list(
    alpha = diagnostics$alpha[c("mean", "CV")],
    K = list(
      mean = unname(selected_K[["mean"]]),
      variance = unname(selected_K[["variance"]]), estimand = "K_J",
      source = "fresh_diagnostics_selected_order", M = 80L,
      pmf = selected_pmf
    ),
    weights = diagnostics$weights["mean"],
    coclustering = diagnostics$coclustering[c("mean", "variance")]
  )
  residuals <- list(diagnostics = refinement_delta)
  tolerances <- list(diagnostics = list(
    absolute = absolute_tolerance, relative = relative_tolerance,
    pmf_mass = pmf_mass_tolerance, refinement = refinement_tolerance
  ))
  verifier_achieved <- list(
    alpha = achieved$alpha,
    K = list(
      mean = unname(verifier_K[["mean"]]),
      variance = unname(verifier_K[["variance"]]), estimand = "K_J",
      source = "fresh_diagnostics_verifier_evidence", M = 160L,
      pmf = verifier_pmf
    ),
    weights = list(mean = verifier_weight),
    coclustering = as.list(verifier_rho)
  )
  component_attempt <- function(name, method) {
    .dpprior_new_attempt(
      id = paste0("diagnostic-", name), stage = "diagnostic_component",
      method = method, start = NULL, bounds = NULL,
      control = list(component = name), exit_code = 0L,
      message = "fresh diagnostic component recomputed", iterations = 0L,
      evaluations = list(function_count = 1L), candidate_parameters = NULL,
      candidate_objective = NULL, elapsed_seconds = 0.001,
      warnings = character(), error = NULL, selected = FALSE,
      reason_code = "component_converged",
      unavailable = c(
        start = "component diagnostic has no optimizer start",
        bounds = "component diagnostic has no optimizer bounds",
        candidate_parameters = "diagnostics do not select fit parameters",
        candidate_objective = "diagnostics do not optimize an objective"
      )
    )
  }
  computation <- .schema23_computation(
    "canonical_prior_diagnostics", parameters = NULL,
    M_selected = 80L, M_verification = 160L
  )
  diagnostic_controls <- list(
    absolute_tolerance = absolute_tolerance,
    relative_tolerance = relative_tolerance,
    pmf_mass_tolerance = pmf_mass_tolerance
  )
  computation$request$controls <- diagnostic_controls
  computation$used$controls <- diagnostic_controls
  computation$attempts <- Map(
    component_attempt, .DPPRIOR_DIAGNOSTIC_COMPONENTS,
    unname(.DPPRIOR_DIAGNOSTIC_ATTEMPT_METHODS)
  )
  computation$termination <- .dpprior_new_termination(
    code = "diagnostics_recomputed", source = "component_aggregation",
    iterations = NULL
  )
  selected_snapshot <- .dpprior_new_snapshot(
    parameters, 80L, achieved, residuals, tolerances, TRUE,
    "fresh_diagnostics_selected_order"
  )
  verifier_snapshot <- .dpprior_new_snapshot(
    parameters, 160L, verifier_achieved, residuals, tolerances, TRUE,
    "fresh_diagnostics_verifier_evidence"
  )
  component_truth <- setNames(rep(TRUE, 4L), .DPPRIOR_DIAGNOSTIC_COMPONENTS)
  verification <- .dpprior_new_verification(
    method = "fresh_component_specific_diagnostics", performed = TRUE,
    passed = TRUE, reason = "all four diagnostic components verified",
    settings = list(M_selected = 80L, M_verification = 160L),
    selected_snapshot = selected_snapshot,
    verifier_snapshot = verifier_snapshot, stability = NULL,
    components = list(component_aggregation = .schema23_check(
      value = component_truth, reference = component_truth,
      tolerance = NULL, operator = "identical",
      source = "fresh_component_specific_checks"
    )),
    invariants = list(
      fixed_parameters = .schema23_check(
        value = TRUE, reference = TRUE, tolerance = NULL,
        operator = "identical", source = "fresh_component_specific_checks"
      ),
      dominance_category_removed = .schema23_check(
        value = TRUE, reference = TRUE, tolerance = NULL,
        operator = "identical", source = "fresh_component_specific_checks"
      )
    )
  )
  .dpprior_new_diagnostics(
    method = "canonical_prior_diagnostics", J = 20L,
    status = "converged",
    usable = TRUE, verified = TRUE,
    parameters = parameters,
    target = list(
      requested_components = .DPPRIOR_DIAGNOSTIC_COMPONENTS,
      warning_policy = NULL
    ),
    achieved = achieved, residuals = residuals, tolerances = tolerances,
    computation = computation, verification = verification,
    provenance = .schema23_provenance(
      method = "canonical_prior_diagnostics"
    ),
    diagnostics = diagnostics
  )
}


.schema23_fit_diagnostics_extension <- function(
    fit, M_selected = NULL, warning_policy = NULL,
    allow_approximate = FALSE) {
  raw <- unclass(fit)
  if (is.null(M_selected)) {
    M_selected <- if (identical(raw$mode, "a1_proxy")) {
      80L
    } else raw$computation$orders$M_selected
  }
  M_selected <- as.integer(M_selected)
  M_required <- as.integer(max(2L * M_selected, M_selected + 40L))
  M_verification <- if (identical(raw$mode, "a1_proxy")) {
    M_required
  } else raw$computation$orders$M_verification_used
  authority <- list(
    method = "fresh_component_specific_diagnostics",
    M_selected = M_selected,
    M_verification_required = M_required,
    M_verification_used = M_verification,
    absolute_tolerance = 1e-10, relative_tolerance = 1e-8,
    pmf_mass_tolerance = .TOL_PMF_SUM,
    warning_policy = warning_policy,
    allow_approximate = allow_approximate
  )
  a <- raw$parameters$a
  b <- raw$parameters$b
  K_evidence <- .get_K_pmf_support(
    raw$J, a, b, M = M_selected, M_verify = M_verification,
    abs_tol = authority$absolute_tolerance,
    rel_tol = authority$relative_tolerance
  )
  selected_pmf <- unname(as.numeric(K_evidence$pmf))
  verifier_pmf <- unname(as.numeric(K_evidence$verification_pmf))
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
    authority$absolute_tolerance + authority$relative_tolerance *
      max(abs(selected), abs(verifier), 1)
  }
  tolerance <- c(
    K.mean = scalar_tolerance(selected_K[["mean"]], verifier_K[["mean"]]),
    K.variance = scalar_tolerance(
      selected_K[["variance"]], verifier_K[["variance"]]
    ),
    K.pmf_l1 = authority$absolute_tolerance + authority$relative_tolerance,
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
      abs(sum(selected_pmf) - 1) <= authority$pmf_mass_tolerance &&
      abs(sum(verifier_pmf) - 1) <= authority$pmf_mass_tolerance,
    weights = delta[["weights.mean"]] <= tolerance[["weights.mean"]],
    coclustering = all(delta[c(
      "coclustering.mean", "coclustering.variance"
    )] <= tolerance[c("coclustering.mean", "coclustering.variance")])
  )
  status <- ifelse(passed, "converged", "approximate")
  usable <- passed | (!passed & allow_approximate)
  policy_results <- if (is.null(warning_policy)) {
    list()
  } else if (identical(warning_policy$estimand, "W_SB")) {
    value <- as.numeric(.diagnostic_wsb_tail(
      warning_policy$weight_threshold, a, b
    ))
    outcome <- if (identical(warning_policy$direction, "above")) {
      if (value > warning_policy$action_threshold) {
        "triggered"
      } else "not_triggered"
    } else if (value < warning_policy$action_threshold) {
      "triggered"
    } else "not_triggered"
    list(list(
      estimand = "W_SB", direction = warning_policy$direction,
      threshold = warning_policy$action_threshold,
      value = value, lower = NULL, upper = NULL,
      outcome = outcome, basis = "exact_tail_probability"
    ))
  } else {
    list(list(
      estimand = "W_max", direction = warning_policy$direction,
      threshold = warning_policy$action_threshold,
      value = NULL, lower = NULL, upper = NULL,
      outcome = "indeterminate", basis = "backend_unavailable"
    ))
  }
  list(
    authority = authority,
    policy_results = policy_results,
    warnings = if (length(policy_results) == 1L &&
                   identical(policy_results[[1L]]$outcome, "triggered")) {
      "canonical fit-attached diagnostic warning"
    } else character(),
    alpha = list(
      status = unname(status[["alpha"]]),
      usable = unname(usable[["alpha"]]),
      verified = unname(passed[["alpha"]]),
      mean = a / b, CV = 1 / sqrt(a)
    ),
    K = list(
      status = unname(status[["K"]]),
      usable = unname(usable[["K"]]),
      verified = unname(passed[["K"]]),
      mean = unname(selected_K[["mean"]]),
      variance = unname(selected_K[["variance"]]),
      pmf = selected_pmf, M = M_selected
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


.schema23_attach_fit_diagnostics <- function(
    fit, M_selected = NULL, warning_policy = NULL,
    allow_approximate = FALSE) {
  diagnostics <- .schema23_fit_diagnostics_extension(
    fit, M_selected = M_selected, warning_policy = warning_policy,
    allow_approximate = allow_approximate
  )
  raw <- unclass(fit)
  aliases <- raw$compatibility$top_level_aliases
  alias_names <- names(aliases)
  alias_values <- raw[alias_names]
  if (length(alias_names) > 0L) raw[alias_names] <- NULL
  if ("diagnostics" %in% alias_names) {
    aliases <- aliases[names(aliases) != "diagnostics"]
    alias_names <- names(aliases)
    alias_values <- alias_values[alias_names]
  }
  raw$compatibility$top_level_aliases <- aliases
  raw$diagnostics <- diagnostics
  for (alias in alias_names) raw[[alias]] <- alias_values[[alias]]
  class(raw) <- class(fit)
  .dpprior_validate_result_v1(raw)
  raw
}


.schema23_empty_condition_evidence <- function() {
  list(
    calibration = NULL, calibration_warnings = list(), diagnostics = NULL,
    diagnostic_warnings = list(), target = NULL, interval = NULL
  )
}


.schema23_sensitivity_condition <- function(code, message = code) {
  contract <- switch(
    code,
    calibration_unusable = list(
      class = "dpprior_calibration_unusable",
      classes = c(
        "dpprior_calibration_unusable", "dpprior_calibration_error",
        "dpprior_error", "error", "dpprior_condition", "condition"
      )
    ),
    fit_diagnostics_approximate = list(
      class = "dpprior_diagnostics_approximation_error",
      classes = c(
        "dpprior_diagnostics_approximation_error", "dpprior_fit_error",
        "dpprior_calibration_error", "dpprior_error", "error",
        "dpprior_condition", "condition"
      )
    ),
    calibration_nonidentifiable_j1 = list(
      class = "dpprior_calibration_nonidentifiable",
      classes = c(
        "dpprior_calibration_nonidentifiable", "dpprior_fit_error",
        "dpprior_calibration_error", "dpprior_error", "error",
        "dpprior_condition", "condition"
      )
    ),
    sensitivity_diagnostic_contract = list(
      class = "dpprior_sensitivity_diagnostic_contract_error",
      classes = c(
        "dpprior_sensitivity_diagnostic_contract_error",
        "dpprior_diagnostics_error", "dpprior_sensitivity_error",
        "dpprior_error", "error", "dpprior_condition", "condition"
      )
    ),
    stop("unknown sensitivity fixture condition")
  )
  c(contract, list(code = code, message = message))
}


.schema23_sensitivity_fit_evidence <- function(
    request = NULL, status = "converged", usable = TRUE, verified = TRUE,
    parameters = .dpprior_new_parameters(2, 3, "log_ab"),
    condition_evidence = .schema23_empty_condition_evidence(),
    input_provenance = NULL) {
  J <- if (is.null(request)) 20L else request[["J", exact = TRUE]]
  M_selected <- if (is.null(request)) 80L else request[["M", exact = TRUE]]
  M_verification <- max(2L * M_selected, M_selected + 40L)
  evaluator <- list(
    method = "gauss_laguerre_marginal_pmf_and_moments",
    M_selected = M_selected, M_verification = M_verification,
    absolute_tolerance = 1e-10, relative_tolerance = 1e-8,
    W_max_point_policy =
      "unavailable_without_retained_typed_backend_evidence"
  )
  selected_snapshot <- verifier_snapshot <- NULL
  if (!is.null(parameters)) {
    logS <- compute_log_stirling(J)
    selected_pmf <- pmf_K_marginal(
      J, parameters$a, parameters$b, logS, M = M_selected,
      M_verify = M_verification, abs_tol = 1e-10, rel_tol = 1e-8,
      strict = FALSE
    )
    verifier_pmf <- attr(
      selected_pmf, ".marginal_verification_pmf", exact = TRUE
    )
    selected_pmf <- unname(as.numeric(selected_pmf[-1L]))
    verifier_pmf <- unname(as.numeric(verifier_pmf[-1L]))
    selected_moments <- .dpprior_target_pmf_moments(selected_pmf)
    verifier_moments <- .dpprior_target_pmf_moments(verifier_pmf)
    if (is.null(request)) {
      request <- list(
        J = J, mu_K = unname(selected_moments[["mean"]]),
        var_K = unname(selected_moments[["variance"]]),
        method = "A2-MN", M = M_selected
      )
    }
    selected_snapshot <- list(
      M = M_selected,
      K = list(
        mean = unname(selected_moments[["mean"]]),
        variance = unname(selected_moments[["variance"]]), pmf = selected_pmf
      ),
      finite = TRUE, source = "selected_order"
    )
    verifier_snapshot <- list(
      M = M_verification,
      K = list(
        mean = unname(verifier_moments[["mean"]]),
        variance = unname(verifier_moments[["variance"]]), pmf = verifier_pmf
      ),
      finite = TRUE, source = "independent_verifier"
    )
  }
  stopifnot(!is.null(request))
  target_route <- if ("K_interval" %in% names(request)) {
    "interval"
  } else if ("target_pmf" %in% names(request)) {
    "strict_pmf"
  } else if ("cv_K" %in% names(request)) {
    "coefficient_of_variation"
  } else if ("confidence" %in% names(request)) {
    "qualitative_confidence"
  } else {
    "direct_variance"
  }
  if (is.null(input_provenance)) {
    input_provenance <- list(
      method_explicit = TRUE,
      confidence_explicit = identical(target_route, "qualitative_confidence"),
      requested_method = request$method, selected_method = request$method,
      is_fallback = FALSE, target_route = target_route
    )
  }
  canonical <- .dpprior_sensitivity_canonical_target(
    request, target_route, "fixture.request"
  )
  target_raw <- unclass(canonical)
  target <- list(
    kind = target_raw$kind, J = target_raw$J,
    request = target_raw$request, used = target_raw$used
  )
  list(
    request = request, input_provenance = input_provenance,
    weight_target = NULL, target = target,
    method = input_provenance$selected_method,
    status = status, usable = usable,
    verified = verified, parameters = parameters, evaluator = evaluator,
    selected_snapshot = selected_snapshot, verifier_snapshot = verifier_snapshot,
    condition_evidence = condition_evidence,
    source = "retained_canonical_fit_evidence"
  )
}


.schema23_sensitivity <- function() {
  evidence <- .schema23_sensitivity_fit_evidence()
  canonical_content <- .dp_sensitivity_canonical_string(list(
    request = .dp_sensitivity_identity_request(evidence$request),
    diagnostics_requested = TRUE, weight_target = evidence$weight_target
  ))
  keys <- .dp_sensitivity_content_key(canonical_content)
  scenarios <- data.frame(
    scenario_key = keys, scenario_label = NA_character_,
    canonical_content = canonical_content,
    base_key = keys, diagnostics_requested = TRUE,
    method_explicit = evidence$input_provenance$method_explicit,
    confidence_explicit = evidence$input_provenance$confidence_explicit,
    effective_method = evidence$input_provenance$requested_method,
    effective_confidence = if (identical(
      evidence$input_provenance$target_route, "qualitative_confidence"
    )) evidence$request$confidence else NA_character_,
    stringsAsFactors = FALSE
  )
  wmax_50 <- unclass(wmax_tail_bounds(0.5, a = 2, b = 3))
  wmax_90 <- unclass(wmax_tail_bounds(0.9, a = 2, b = 3))
  selected_K <- evidence$selected_snapshot$K
  scenario_results <- data.frame(
    scenario_key = keys, status = "converged", usable = TRUE,
    verified = TRUE, stringsAsFactors = FALSE
  )
  metric_values <- c(
    a = 2, b = 3, E_alpha = 2 / 3, CV_alpha = 1 / sqrt(2),
    E_K_J = selected_K$mean, Var_K_J = selected_K$variance,
    CV_K_J = sqrt(selected_K$variance) / selected_K$mean,
    interval_requested = NA_real_, interval_achieved = NA_real_,
    interval_residual = NA_real_, interval_left_tail = NA_real_,
    interval_right_tail = NA_real_,
    E_W_SB = as.numeric(mean_w1(2, 3, 80L)),
    P_W_SB_gt_50 = as.numeric(.diagnostic_wsb_tail(0.5, 2, 3)),
    P_W_SB_gt_90 = as.numeric(.diagnostic_wsb_tail(0.9, 2, 3)),
    P_W_max_gt_50 = NA_real_, P_W_max_gt_90 = NA_real_,
    P_W_max_gt_50_lower_bound = wmax_50$lower_bound,
    P_W_max_gt_50_upper_bound = wmax_50$upper_bound,
    P_W_max_gt_90_lower_bound = wmax_90$lower_bound,
    P_W_max_gt_90_upper_bound = wmax_90$upper_bound,
    E_rho = as.numeric(mean_rho(2, 3, 80L))
  )
  stopifnot(identical(names(metric_values), .DPPRIOR_SENSITIVITY_METRICS))
  for (metric in .DPPRIOR_SENSITIVITY_METRICS) {
    scenario_results[[metric]] <- unname(metric_values[[metric]])
  }
  conditions <- setNames(list(evidence$condition_evidence), keys)
  interval_checks <- setNames(list(list(
    requested = NULL, selected = NULL, verification = NULL,
    status = NULL, source = NULL, usable = FALSE, verified = FALSE,
    reason = "not_interval_scenario"
  )), keys)
  metric_count <- length(.DPPRIOR_SENSITIVITY_METRICS)
  metrics_long <- data.frame(
    scenario_key = rep(keys, each = metric_count),
    metric = .DPPRIOR_SENSITIVITY_METRICS,
    value = unname(metric_values),
    reason = ifelse(
      is.na(metric_values), "not_available_for_base_scenario", NA_character_
    ),
    source = unname(.DPPRIOR_SENSITIVITY_METRIC_SOURCE),
    component = unname(.DPPRIOR_SENSITIVITY_METRIC_COMPONENT),
    status = rep("converged", metric_count),
    usable = rep(TRUE, metric_count), verified = rep(TRUE, metric_count),
    stringsAsFactors = FALSE
  )
  local <- data.frame(
    scenario_key = character(), axis = character(), axis_value = numeric(),
    settings_key = character(), lower_scenario_key = character(),
    upper_scenario_key = character(), lower_value = numeric(),
    upper_value = numeric(), metric = character(), component = character(),
    derivative = numeric(), method = character(), reason = character(),
    stringsAsFactors = FALSE
  )
  sensitivity <- list(
    scenarios = scenarios, scenario_results = scenario_results,
    fit_evidence = setNames(list(evidence), keys),
    conditions = conditions, interval_checks = interval_checks,
    metrics_long = metrics_long, local = local,
    global = list(
      scenario_count = 1L, converged_count = 1L, failed_count = 0L,
      metric_count = as.integer(metric_count)
    ),
    metadata = list(J = 20L)
  )
  achieved <- list(scenario_count = 1L)
  residuals <- list(
    missing_metric_count = as.integer(sum(is.na(metric_values)))
  )
  tolerances <- list(expected_metric_count = as.integer(metric_count))
  selected_snapshot <- .dpprior_new_snapshot(
    NULL, NULL, achieved, residuals, tolerances, TRUE, "sensitivity_tables"
  )
  verifier_snapshot <- .dpprior_new_snapshot(
    NULL, NULL, achieved, residuals, tolerances, TRUE,
    "sensitivity_reconciliation"
  )
  truth <- c(
    scenario_key_identity = TRUE, condition_key_identity = TRUE,
    interval_key_identity = TRUE, metric_grid_identity = TRUE,
    row_status_identity = TRUE, top_status_identity = TRUE,
    global_summary_identity = TRUE
  )
  verification <- .dpprior_new_verification(
    method = "reconciliation", performed = TRUE, passed = TRUE,
    reason = "all keys reconciled",
    settings = list(
      scenario_count = 1L, scenario_keys = keys,
      metric_names = .DPPRIOR_SENSITIVITY_METRICS
    ),
    selected_snapshot = selected_snapshot, verifier_snapshot = verifier_snapshot,
    stability = NULL,
    components = list(reconciliation = .schema23_check(
      value = truth, reference = setNames(rep(TRUE, length(truth)), names(truth)),
      tolerance = NULL, operator = "identical",
      source = "sensitivity_reconciliation"
    )),
    invariants = setNames(lapply(
      c(
        "unique_keys", "lexicographic_order", "finite_or_reason",
        "failure_preservation"
      ),
      function(name) .schema23_check(
        value = TRUE, reference = TRUE, tolerance = NULL,
        operator = "identical", source = "sensitivity_reconciliation"
      )
    ), c(
      "unique_keys", "lexicographic_order", "finite_or_reason",
      "failure_preservation"
    ))
  )
  computation <- .schema23_computation("elicitation_sensitivity")
  computation$termination <- .dpprior_new_termination(
    code = "deterministic", source = "sensitivity_reconciliation",
    iterations = NULL
  )
  .dpprior_new_sensitivity(
    method = "elicitation_sensitivity", J = 20L,
    status = "converged",
    usable = TRUE, verified = TRUE,
    target = list(defaults = list(J = 20L)),
    achieved = achieved, residuals = residuals, tolerances = tolerances,
    computation = computation,
    verification = verification,
    provenance = .schema23_provenance(
      "elicitation_sensitivity", parameterization = "none"
    ),
    sensitivity = sensitivity
  )
}


.schema23_sensitivity_route_evidence <- function(
    route,
    method_explicit = TRUE,
    confidence_explicit = FALSE,
    confidence = "medium",
    fallback = FALSE,
    structural_K0 = FALSE,
    A1_diagnostic_condition = FALSE) {
  base <- .schema23_sensitivity_fit_evidence()
  selected <- base$selected_snapshot$K
  request <- switch(
    route,
    direct_variance = list(
      J = 20L, mu_K = selected$mean, var_K = selected$variance,
      method = "A2-MN", M = 80L
    ),
    qualitative_confidence = list(
      J = 20L, mu_K = selected$mean, confidence = confidence,
      method = "A1", M = 80L
    ),
    coefficient_of_variation = list(
      J = 20L, mu_K = selected$mean,
      cv_K = sqrt(selected$variance) / selected$mean,
      method = "A2-MN", M = 80L
    ),
    strict_pmf = list(
      J = 20L,
      target_pmf = if (structural_K0) c(0, selected$pmf) else selected$pmf,
      method = "A2-KL", M = 80L
    ),
    stop("unknown sensitivity fixture route")
  )
  selected_method <- if (fallback) "A2-MN+NM" else request$method
  input_provenance <- list(
    method_explicit = method_explicit,
    confidence_explicit = confidence_explicit,
    requested_method = request$method,
    selected_method = selected_method,
    is_fallback = fallback,
    target_route = route
  )
  A1 <- identical(selected_method, "A1")
  condition_evidence <- .schema23_empty_condition_evidence()
  if (A1_diagnostic_condition) {
    stopifnot(A1)
    condition_evidence$diagnostics <- .schema23_sensitivity_condition(
      "fit_diagnostics_approximate",
      "A1 fit retained with approximate diagnostics evidence"
    )
  }
  .schema23_sensitivity_fit_evidence(
    request = request,
    status = if (A1) "approximate" else "converged",
    usable = TRUE,
    verified = !A1,
    condition_evidence = condition_evidence,
    input_provenance = input_provenance
  )
}


.schema23_sensitivity_interval_evidence <- function(
    type = c("equal_tail", "central_mass", "hard_bounds")) {
  type <- match.arg(type)
  raw_interval <- switch(
    type,
    equal_tail = list(
      lower = 3L, upper = 10L, type = type, coverage = 0.8,
      family = "maxent"
    ),
    central_mass = list(
      lower = 3L, upper = 10L, type = type, coverage = 0.8,
      family = "maxent"
    ),
    hard_bounds = list(
      lower = 5L, upper = 5L, type = type, family = "maxent"
    )
  )
  target <- .dp_target_K(
    J = 20L, mu_K = if (identical(type, "central_mass")) 6.5 else NULL,
    K_interval = raw_interval
  )
  interval <- unclass(target)$used$interval
  request <- list(
    J = 20L, K_interval = interval, mu_K = interval$mu_K,
    method = "A2-KL", M = 80L
  )
  conditions <- .schema23_empty_condition_evidence()
  conditions$calibration <- .schema23_sensitivity_condition(
    "calibration_unusable",
    "A2-KL interval fit retained as an unusable approximation"
  )
  .schema23_sensitivity_fit_evidence(
    request = request, status = "approximate", usable = FALSE,
    verified = FALSE, condition_evidence = conditions
  )
}


.schema23_sensitivity_from_evidence <- function(evidence, validate = TRUE) {
  stopifnot(!identical(evidence$input_provenance$target_route, "interval"))
  out <- .schema23_sensitivity()
  old_key <- out$sensitivity$scenarios$scenario_key[[1L]]
  canonical_content <- .dp_sensitivity_canonical_string(list(
    request = .dp_sensitivity_identity_request(evidence$request),
    diagnostics_requested = TRUE, weight_target = evidence$weight_target
  ))
  key <- .dp_sensitivity_content_key(canonical_content)
  row_status <- evidence$status
  row_usable <- evidence$usable
  row_verified <- evidence$verified
  top_usable <- row_status %in% c("converged", "boundary") && row_usable
  top_verified <- row_status %in%
    c("converged", "boundary", "infeasible") && row_verified

  out$sensitivity$scenarios$scenario_key <- key
  out$sensitivity$scenarios$base_key <- key
  out$sensitivity$scenarios$canonical_content <- canonical_content
  out$sensitivity$scenarios$method_explicit <-
    evidence$input_provenance$method_explicit
  out$sensitivity$scenarios$confidence_explicit <-
    evidence$input_provenance$confidence_explicit
  out$sensitivity$scenarios$effective_method <-
    evidence$input_provenance$requested_method
  out$sensitivity$scenarios$effective_confidence <- if (identical(
    evidence$input_provenance$target_route, "qualitative_confidence"
  )) evidence$request$confidence else NA_character_
  out$sensitivity$scenario_results$scenario_key <- key
  out$sensitivity$scenario_results$status <- row_status
  out$sensitivity$scenario_results$usable <- row_usable
  out$sensitivity$scenario_results$verified <- row_verified
  out$sensitivity$fit_evidence <- setNames(list(evidence), key)
  out$sensitivity$conditions <- setNames(
    list(evidence$condition_evidence), key
  )
  names(out$sensitivity$interval_checks) <- key
  out$sensitivity$metrics_long$scenario_key <- key
  out$sensitivity$metrics_long$status <- row_status
  out$sensitivity$metrics_long$usable <- row_usable
  out$sensitivity$metrics_long$verified <- row_verified
  out$sensitivity$global$converged_count <- as.integer(
    row_status %in% c("converged", "boundary")
  )
  out$sensitivity$global$failed_count <- as.integer(
    row_status %in% c("failed", "infeasible")
  )
  out$status <- row_status
  out$usable <- top_usable
  out$verified <- top_verified
  out$verification$settings$scenario_keys <- key
  out$verification$passed <- top_verified
  out$verification$reason <- if (top_verified) {
    "all keys reconciled"
  } else {
    "retained route is not decision-ready"
  }
  stopifnot(!identical(old_key, key) || identical(
    evidence$request, .schema23_sensitivity_fit_evidence()$request
  ))
  if (validate) .dpprior_validate_result_v1(out)
  out
}


.schema23_diagnostics_approximate <- function() {
  out <- .schema23_diagnostics()
  out$tolerances$diagnostics$absolute <- 0
  out$tolerances$diagnostics$relative <- 0
  out$tolerances$diagnostics$refinement[] <- 0
  out$computation$request$controls$absolute_tolerance <- 0
  out$computation$request$controls$relative_tolerance <- 0
  out$computation$used$controls$absolute_tolerance <- 0
  out$computation$used$controls$relative_tolerance <- 0
  expected <- .dpprior_expected_diagnostics_evidence(out)
  stopifnot(any(!expected$component_pass))
  component_status <- ifelse(
    expected$component_pass, "converged", "approximate"
  )
  for (name in .DPPRIOR_DIAGNOSTIC_COMPONENTS) {
    out$diagnostics[[name]]$status <- unname(component_status[[name]])
    out$diagnostics[[name]]$usable <-
      unname(expected$component_pass[[name]])
    out$diagnostics[[name]]$verified <-
      unname(expected$component_pass[[name]])
  }
  out$status <- "approximate"
  out$usable <- FALSE
  out$verified <- FALSE
  out$achieved <- expected$selected
  out$residuals <- list(diagnostics = expected$delta)
  out$verification$passed <- FALSE
  out$verification$reason <- "one or more fixed-tolerance checks failed"
  out$verification$selected_snapshot$achieved <- expected$selected
  out$verification$selected_snapshot$residuals <- out$residuals
  out$verification$selected_snapshot$tolerances <- out$tolerances
  out$verification$verifier_snapshot$achieved <- expected$verifier
  out$verification$verifier_snapshot$residuals <- out$residuals
  out$verification$verifier_snapshot$tolerances <- out$tolerances
  component_truth <- expected$component_pass
  out$verification$components$component_aggregation <- .schema23_check(
    value = component_truth,
    reference = setNames(rep(TRUE, length(component_truth)),
                         names(component_truth)),
    tolerance = NULL, operator = "identical",
    source = "fresh_component_specific_checks"
  )
  for (index in seq_along(out$computation$attempts)) {
    out$computation$attempts[[index]]$reason_code <- paste0(
      "component_", component_status[[index]]
    )
  }
  .dpprior_validate_result_v1(out)
  out
}


.schema23_sensitivity_failed <- function() {
  out <- .schema23_sensitivity()
  key <- out$sensitivity$scenarios$scenario_key[[1L]]
  failure_condition <- .schema23_sensitivity_condition(
    "calibration_unusable", "calibration backend returned a failed fit"
  )
  condition_evidence <- .schema23_empty_condition_evidence()
  condition_evidence$calibration <- failure_condition
  evidence <- out$sensitivity$fit_evidence[[key]]
  evidence$status <- "failed"
  evidence$usable <- FALSE
  evidence$verified <- FALSE
  evidence[c("parameters", "selected_snapshot", "verifier_snapshot")] <-
    list(NULL, NULL, NULL)
  evidence$condition_evidence <- condition_evidence
  out$sensitivity$fit_evidence[[key]] <- evidence
  out$sensitivity$conditions[[key]] <- condition_evidence
  out$status <- "failed"
  out$usable <- FALSE
  out$verified <- FALSE
  out$sensitivity$scenario_results$status <- "failed"
  out$sensitivity$scenario_results$usable <- FALSE
  out$sensitivity$scenario_results$verified <- FALSE
  for (metric in .DPPRIOR_SENSITIVITY_METRICS) {
    out$sensitivity$scenario_results[[metric]] <- NA_real_
  }
  out$sensitivity$metrics_long$value <- NA_real_
  out$sensitivity$metrics_long$reason <- "calibration_unusable"
  out$sensitivity$metrics_long$status <- "failed"
  out$sensitivity$metrics_long$usable <- FALSE
  out$sensitivity$metrics_long$verified <- FALSE
  out$sensitivity$global$converged_count <- 0L
  out$sensitivity$global$failed_count <- 1L
  out$residuals$missing_metric_count <- as.integer(
    length(.DPPRIOR_SENSITIVITY_METRICS)
  )
  out$verification$selected_snapshot$residuals <- out$residuals
  out$verification$verifier_snapshot$residuals <- out$residuals
  out$verification$passed <- FALSE
  out$verification$reason <- "failed scenario evidence preserved"
  .dpprior_validate_result_v1(out)
  out
}


.schema23_sensitivity_infeasible <- function() {
  out <- .schema23_sensitivity()
  request <- list(J = 1L, target_pmf = 1, method = "A2-KL", M = 20L)
  certificate <- .schema23_sensitivity_condition(
    "calibration_nonidentifiable_j1",
    "K_1 is identically one for every positive Gamma prior, so a and b are not identified."
  )
  condition_evidence <- .schema23_empty_condition_evidence()
  condition_evidence$calibration <- certificate
  evidence <- .schema23_sensitivity_fit_evidence(
    request = request, status = "infeasible", usable = FALSE,
    verified = TRUE, parameters = NULL,
    condition_evidence = condition_evidence
  )
  canonical_content <- .dp_sensitivity_canonical_string(list(
    request = .dp_sensitivity_identity_request(request),
    diagnostics_requested = FALSE, weight_target = NULL
  ))
  key <- .dp_sensitivity_content_key(canonical_content)
  out$J <- 1L
  out$target$defaults$J <- 1L
  out$sensitivity$metadata$J <- 1L
  out$sensitivity$scenarios$canonical_content <- canonical_content
  out$sensitivity$scenarios$base_key <- key
  out$sensitivity$scenarios$scenario_key <- key
  out$sensitivity$scenarios$diagnostics_requested <- FALSE
  out$sensitivity$scenarios$effective_method <- "A2-KL"
  out$sensitivity$scenario_results$scenario_key <- key
  out$sensitivity$scenario_results$status <- "infeasible"
  out$sensitivity$scenario_results$usable <- FALSE
  out$sensitivity$scenario_results$verified <- TRUE
  for (metric in .DPPRIOR_SENSITIVITY_METRICS) {
    out$sensitivity$scenario_results[[metric]] <- NA_real_
  }
  out$sensitivity$fit_evidence <- setNames(list(evidence), key)
  out$sensitivity$conditions <- setNames(list(condition_evidence), key)
  names(out$sensitivity$interval_checks) <- key
  out$sensitivity$metrics_long$scenario_key <- key
  out$sensitivity$metrics_long$value <- NA_real_
  out$sensitivity$metrics_long$reason <- "calibration_nonidentifiable_j1"
  out$sensitivity$metrics_long$status <- "infeasible"
  out$sensitivity$metrics_long$usable <- FALSE
  out$sensitivity$metrics_long$verified <- TRUE
  out$sensitivity$global$converged_count <- 0L
  out$sensitivity$global$failed_count <- 1L
  out$status <- "infeasible"
  out$usable <- FALSE
  out$verified <- TRUE
  out$residuals$missing_metric_count <- as.integer(
    length(.DPPRIOR_SENSITIVITY_METRICS)
  )
  out$verification$settings$scenario_keys <- key
  out$verification$selected_snapshot$residuals <- out$residuals
  out$verification$verifier_snapshot$residuals <- out$residuals
  out$verification$reason <- "certified J=1 nonidentifiability reconciled"
  .dpprior_validate_result_v1(out)
  out
}


.schema23_sensitivity_interval <- function() {
  out <- .schema23_sensitivity()
  old_key <- out$sensitivity$scenarios$scenario_key[[1L]]
  raw_requested <- list(
    lower = 3L, upper = 10L, type = "equal_tail", coverage = 0.8,
    family = "maxent"
  )
  target_attempt <- .dp_target_K(J = 20L, K_interval = raw_requested)
  requested <- unclass(target_attempt)$interval
  request <- list(
    J = 20L, K_interval = requested, mu_K = requested$mu_K,
    method = "A2-KL", M = 80L
  )
  condition_evidence <- .schema23_empty_condition_evidence()
  condition_evidence$calibration <- .schema23_sensitivity_condition(
    "calibration_unusable",
    "A2-KL interval fit retained as an unusable approximation"
  )
  evidence <- .schema23_sensitivity_fit_evidence(
    request = request, status = "approximate", usable = FALSE,
    verified = FALSE, condition_evidence = condition_evidence
  )
  canonical_content <- .dp_sensitivity_canonical_string(list(
    request = .dp_sensitivity_identity_request(request),
    diagnostics_requested = TRUE, weight_target = NULL
  ))
  key <- .dp_sensitivity_content_key(canonical_content)
  out$sensitivity$scenarios$canonical_content <- canonical_content
  out$sensitivity$scenarios$base_key <- key
  out$sensitivity$scenarios$scenario_key <- key
  out$sensitivity$scenarios$effective_method <- "A2-KL"
  out$sensitivity$scenario_results$scenario_key <- key
  out$sensitivity$scenario_results$status <- "approximate"
  out$sensitivity$scenario_results$usable <- FALSE
  out$sensitivity$scenario_results$verified <- FALSE
  out$sensitivity$metrics_long$scenario_key <- key
  out$sensitivity$metrics_long$status <- "approximate"
  out$sensitivity$metrics_long$usable <- FALSE
  out$sensitivity$metrics_long$verified <- FALSE
  out$sensitivity$fit_evidence <- setNames(list(evidence), key)
  out$sensitivity$conditions <- setNames(list(condition_evidence), key)
  out$verification$settings$scenario_keys <- key
  support <- seq_len(20L)
  interval_masses <- function(pmf) {
    coverage <- sum(pmf[support >= 3L & support <= 10L])
    list(
      coverage = coverage, lower_tail = sum(pmf[support < 3L]),
      upper_tail = sum(pmf[support > 10L]),
      coverage_residual = coverage - 0.8
    )
  }
  selected <- interval_masses(evidence$selected_snapshot$K$pmf)
  verifier <- interval_masses(evidence$verifier_snapshot$K$pmf)
  out$sensitivity$interval_checks <- setNames(list(list(
    requested = requested, selected = selected,
    verification = c(verifier, list(
      tolerance = 1e-8, passed = FALSE,
      source = "independent_interval_backcheck"
    )),
    status = "approximate", source = "wrapper_backcheck_selected",
    usable = FALSE, verified = FALSE,
    reason = paste(
      "Selected and verifier interval constraints did not pass",
      "fixed tolerance"
    )
  )), key)
  interval_metrics <- c(
    interval_requested = 0.8, interval_achieved = selected$coverage,
    interval_residual = selected$coverage_residual,
    interval_left_tail = selected$lower_tail,
    interval_right_tail = selected$upper_tail
  )
  for (metric in names(interval_metrics)) {
    out$sensitivity$scenario_results[[metric]] <-
      unname(interval_metrics[[metric]])
    row <- match(metric, out$sensitivity$metrics_long$metric)
    out$sensitivity$metrics_long$value[[row]] <-
      unname(interval_metrics[[metric]])
    out$sensitivity$metrics_long$reason[[row]] <- NA_character_
  }
  out$residuals$missing_metric_count <- as.integer(sum(
    is.na(out$sensitivity$metrics_long$value)
  ))
  out$verification$selected_snapshot$residuals <- out$residuals
  out$verification$verifier_snapshot$residuals <- out$residuals
  out$sensitivity$global$converged_count <- 0L
  out$sensitivity$global$failed_count <- 0L
  out$status <- "approximate"
  out$usable <- FALSE
  out$verified <- FALSE
  out$verification$passed <- FALSE
  out$verification$reason <- "interval scenario retained as diagnostic only"
  stopifnot(!identical(old_key, key))
  .dpprior_validate_result_v1(out)
  out
}


.schema23_sensitivity_interval_infeasible <- function(
    certificate = c("empty_group", "mean_hull")) {
  certificate <- match.arg(certificate)
  out <- .schema23_sensitivity_infeasible()
  J <- if (identical(certificate, "empty_group")) 10L else 20L
  mu_K <- if (identical(certificate, "empty_group")) NULL else 2
  raw_interval <- if (identical(certificate, "empty_group")) {
    list(
      lower = 1L, upper = 4L, type = "equal_tail", coverage = 0.8,
      family = "maxent"
    )
  } else {
    list(
      lower = 3L, upper = 10L, type = "equal_tail", coverage = 0.8,
      family = "maxent"
    )
  }
  target_attempt <- tryCatch(
    .dp_target_K(J = J, mu_K = mu_K, K_interval = raw_interval),
    error = function(condition) condition[["result", exact = TRUE]]
  )
  stopifnot(inherits(target_attempt, "dpprior_K_target"))
  target_raw <- unclass(target_attempt)
  normalized_interval <- target_raw$used$interval
  request <- list(
    J = J, K_interval = normalized_interval,
    mu_K = normalized_interval$mu_K, method = "A2-KL", M = 80L
  )
  canonical_target <- .dpprior_sensitivity_canonical_target(
    request, "interval", "fixture.interval_infeasible.request"
  )
  target_condition <- tryCatch(
    .dp_target_K_stop_unusable(canonical_target),
    error = function(condition) condition
  )
  stopifnot(inherits(target_condition, "condition"))
  condition_summary <- list(
    class = class(target_condition)[[1L]], classes = class(target_condition),
    code = target_condition[["code", exact = TRUE]],
    message = conditionMessage(target_condition)
  )
  condition_evidence <- .schema23_empty_condition_evidence()
  condition_evidence$calibration <- condition_summary
  condition_evidence$target <- condition_summary
  input_provenance <- list(
    method_explicit = TRUE, confidence_explicit = FALSE,
    requested_method = "A2-KL", selected_method = "A2-KL",
    is_fallback = FALSE, target_route = "interval"
  )
  evidence <- .schema23_sensitivity_fit_evidence(
    request = request, status = "infeasible", usable = FALSE,
    verified = TRUE, parameters = NULL,
    condition_evidence = condition_evidence,
    input_provenance = input_provenance
  )
  canonical_content <- .dp_sensitivity_canonical_string(list(
    request = .dp_sensitivity_identity_request(request),
    diagnostics_requested = FALSE, weight_target = NULL
  ))
  key <- .dp_sensitivity_content_key(canonical_content)
  old_key <- out$sensitivity$scenarios$scenario_key[[1L]]
  out$J <- J
  out$target$defaults$J <- J
  out$sensitivity$metadata$J <- J
  out$sensitivity$scenarios$scenario_key <- key
  out$sensitivity$scenarios$base_key <- key
  out$sensitivity$scenarios$canonical_content <- canonical_content
  out$sensitivity$scenarios$diagnostics_requested <- FALSE
  out$sensitivity$scenarios$method_explicit <- TRUE
  out$sensitivity$scenarios$confidence_explicit <- FALSE
  out$sensitivity$scenarios$effective_method <- "A2-KL"
  out$sensitivity$scenarios$effective_confidence <- NA_character_
  out$sensitivity$scenario_results$scenario_key <- key
  out$sensitivity$fit_evidence <- setNames(list(evidence), key)
  out$sensitivity$conditions <- setNames(list(condition_evidence), key)
  out$sensitivity$interval_checks <- setNames(list(list(
    requested = normalized_interval, selected = NULL, verification = NULL,
    status = "infeasible", source = "target_infeasibility_certificate",
    usable = FALSE, verified = TRUE, reason = condition_summary$message
  )), key)
  out$sensitivity$metrics_long$scenario_key <- key
  out$verification$settings$scenario_keys <- key
  out$sensitivity$scenario_results$interval_requested <-
    normalized_interval$coverage
  interval_requested_row <- match(
    "interval_requested", out$sensitivity$metrics_long$metric
  )
  out$sensitivity$metrics_long$value[[interval_requested_row]] <-
    normalized_interval$coverage
  out$sensitivity$metrics_long$reason[[interval_requested_row]] <- NA_character_
  out$residuals$missing_metric_count <- as.integer(
    length(.DPPRIOR_SENSITIVITY_METRICS) - 1L
  )
  out$verification$selected_snapshot$residuals <- out$residuals
  out$verification$verifier_snapshot$residuals <- out$residuals
  stopifnot(!identical(old_key, key))
  .dpprior_validate_result_v1(out)
  out
}


.schema23_legacy_A2_source <- function(mode, wrapper = FALSE,
                                       include_diagnostics = FALSE) {
  if (identical(mode, "a2_moment")) {
    raw <- list(
      a = 1.27223464009856, b = 0.58075377814199, J = 20,
      target = list(mu_K = 5, var_K = 8, type = "moments"),
      method = "A2-MN", status = "success", converged = TRUE,
      iterations = 7L, termination = "residual",
      fit = list(
        mu_K = 5.00000000001361, var_K = 7.99999999987637,
        residual = 1.24381729905704e-10
      ),
      diagnostics = list(
        a0 = 4, b0 = 2.99573227355399, tol_F = 1e-8,
        tol_step = 1e-10, M = 80L, fallback_used = FALSE
      ),
      trace = data.frame(
        iteration = 1L, residual = 1.24381729905704e-10
      )
    )
    if (wrapper) {
      raw$target <- list(
        mu_K = 5, var_K = 8, var_K_used = 8,
        confidence = NULL, type = "moments"
      )
      raw$solver_diagnostics <- raw$diagnostics
      raw$diagnostics <- NULL
      raw <- raw[c(
        "a", "b", "J", "target", "method", "status", "converged",
        "iterations", "termination", "fit", "solver_diagnostics", "trace"
      )]
      if (include_diagnostics) {
        raw$diagnostics <- list(retained_nested_bundle = TRUE)
      }
    }
    return(structure(raw, class = "DPprior_fit"))
  }
  pmf <- rep(1 / 20, 20L)
  support <- seq_along(pmf)
  mean <- sum(support * pmf)
  variance <- sum((support - mean)^2 * pmf)
  structure(list(
    a = 1.5, b = 0.7, J = 20L,
    target = list(
      type = "chisq", pmf = pmf, mu_K = mean, var_K = variance,
      df = 6.25, scale = 0.8,
      mu_K_discrete = mean, var_K_discrete = variance
    ),
    method = "A2-KL", status = "success", converged = TRUE,
    iterations = 8L, termination = "optim_converged",
    fit = list(mu_K = 5, var_K = 8, kl = 0.02, residual = 0.02),
    diagnostics = list(M = 80L, fallback_used = FALSE),
    trace = data.frame(evaluation = 1L, kl = 0.02)
  ), class = "DPprior_fit")
}


.schema23_migrated_A2 <- function(mode, verify = FALSE, wrapper = FALSE,
                                  include_diagnostics = FALSE,
                                  M_verify = NULL) {
  suppressWarnings(upgrade_DPprior_object(
    .schema23_legacy_A2_source(
      mode, wrapper = wrapper, include_diagnostics = include_diagnostics
    ),
    verify = verify, M_verify = M_verify, allow_legacy = TRUE
  ))
}


.schema23_mutate_plain_leaf <- function(x) {
  if (typeof(x) == "list" && is.list(x)) {
    stopifnot(length(x) > 0L)
    first <- names(x)[[1L]]
    x[[first]] <- .schema23_mutate_plain_leaf(x[[first, exact = TRUE]])
    return(x)
  }
  if (is.logical(x)) {
    x[[1L]] <- !x[[1L]]
    return(x)
  }
  if (is.numeric(x)) {
    x[[1L]] <- x[[1L]] + if (is.integer(x)) 1L else 0.25
    return(x)
  }
  if (is.character(x)) {
    x[[1L]] <- paste0(x[[1L]], "_forged")
    return(x)
  }
  stop("fixture mutation helper received an unsupported leaf")
}


.schema23_migration_decision_and_consumer_state <- function(x) {
  target_K <- unclass(x$target$K)
  target_K$compatibility <- NULL
  canonical_calculation <- if (is.null(x$parameters)) {
    NULL
  } else {
    c(
      E_alpha = x$parameters$a / x$parameters$b,
      CV_alpha = 1 / sqrt(x$parameters$a)
    )
  }
  list(
    canonical = list(
      mode = x$mode, method = x$method, J = x$J, status = x$status,
      usable = x$usable, verified = x$verified, message = x$message,
      parameters = x$parameters, target_K = target_K,
      achieved = x$achieved, residuals = x$residuals,
      tolerances = x$tolerances, computation = x$computation,
      verification = x$verification, provenance = x$provenance
    ),
    fixed_candidate_audit = x$compatibility$views[
      "fixed_candidate_recomputation"
    ],
    canonical_calculation = canonical_calculation,
    summary = unclass(summary(x, print_output = FALSE)),
    data_frame = as.data.frame(x),
    printed = capture.output(print(x))
  )
}


test_that("schema/status constructors reject coercive scalar inputs", {
  expect_identical(
    .dpprior_schema("result"),
    list(name = "dpprior.result", version = 1L)
  )
  expect_error(
    .dpprior_schema(matrix("result", 1, 1)),
    class = "dpprior_schema_error"
  )
  expect_error(
    .dpprior_schema(structure("result", class = "forged")),
    class = "dpprior_schema_error"
  )
  expect_error(
    .dpprior_new_status(matrix("converged", 1, 1), TRUE, TRUE),
    class = "dpprior_schema_error"
  )
  expect_error(
    .dpprior_new_status("converged", matrix(TRUE, 1, 1), TRUE),
    class = "dpprior_schema_error"
  )
  expect_error(
    .dpprior_new_status("infeasible", FALSE, FALSE),
    class = "dpprior_schema_error"
  )
  expect_error(
    .dpprior_new_target_K(
      kind = "moments", J = 2.5, request = list(), normalized = list(),
      used = list(), status = "failed", usable = FALSE, verified = FALSE,
      computation = .schema23_computation("target_moments"),
      verification = .dpprior_new_verification(
        "none", FALSE, FALSE, reason = "not run"
      ),
      provenance = .schema23_provenance(
        "target_moments", parameterization = "none"
      )
    ),
    class = "dpprior_schema_error"
  )
})


test_that("schema errors expose stable machine-readable fields", {
  condition <- tryCatch(
    .dpprior_new_status("converged", TRUE, FALSE, "bad quartet"),
    dpprior_schema_error = identity
  )
  expect_s3_class(condition, "dpprior_schema_error")
  expect_s3_class(condition, "dpprior_error")
  expect_s3_class(condition, "error")
  expect_identical(condition$code, "status_quartet")
  expect_identical(condition$path, "status_record")
  expect_true(!is.null(condition$expected))
  expect_null(condition$call)
})


test_that("canonical leaf constructors use exact stable shapes", {
  parameters <- .dpprior_new_parameters(2, 3, "log_ab")
  attempt <- .schema23_computation("A2-MN", parameters, 80L, 160L)$attempts[[1L]]
  expect_identical(names(attempt), .DPPRIOR_ATTEMPT_FIELDS)
  expect_identical(names(.schema23_computation()), .DPPRIOR_COMPUTATION_FIELDS)
  expect_identical(
    names(.schema23_fit_parts()$verification), .DPPRIOR_VERIFICATION_FIELDS
  )
  expect_identical(
    names(.schema23_provenance()), .DPPRIOR_PROVENANCE_FIELDS
  )
})


test_that("K and weight target constructors validate exact v1 contracts", {
  target <- .schema23_target()
  expect_s3_class(target, "dpprior_K_target")
  expect_identical(names(unclass(target)), .DPPRIOR_TARGET_FIELDS)
  expect_invisible(.dpprior_validate_target_v1(target))
  expect_null(target$pmf)

  for (metric in c("wsb_mean", "wsb_tail", "wsb_quantile",
                   "wmax_tail_upper")) {
    expect_invisible(.dpprior_validate_weight_target_v1(
      .schema23_weight_target(
        relation = if (metric == "wmax_tail_upper") "at_most" else "target",
        metric = metric
      )
    ))
  }
})


test_that("all result-mode constructors share one canonical spine", {
  modes <- c(
    "a1_proxy", "a2_moment", "a2_kl", "dual_hard", "dual_soft",
    "dual_legacy"
  )
  fits <- lapply(modes, .schema23_fit)
  for (fit in fits) {
    expect_identical(
      names(unclass(fit))[seq_along(.DPPRIOR_RESULT_COMMON_FIELDS)],
      .DPPRIOR_RESULT_COMMON_FIELDS
    )
    expect_invisible(.dpprior_validate_result_v1(fit))
    expect_invisible(.dpprior_validate_object(fit))
  }
  expect_invisible(.dpprior_validate_object(.schema23_diagnostics()))
  expect_invisible(.dpprior_validate_object(.schema23_sensitivity()))
})


test_that("required fields, exact names, and nested ordinary lists fail closed", {
  fit <- .schema23_fit()
  missing <- fit
  missing$target <- NULL
  expect_error(
    .dpprior_validate_result_v1(missing), class = "dpprior_schema_error"
  )

  unknown <- unclass(fit)
  unknown$surprise <- 1
  class(unknown) <- class(fit)
  expect_error(
    .dpprior_validate_result_v1(unknown), class = "dpprior_schema_error"
  )

  duplicated <- unclass(fit)
  duplicated <- c(duplicated, list(status = "converged"))
  class(duplicated) <- class(fit)
  expect_error(
    .dpprior_validate_result_v1(duplicated), class = "dpprior_schema_error"
  )

  forged_attempt <- fit
  class(forged_attempt$computation$attempts[[1L]]) <- c("forged", "list")
  expect_error(
    .dpprior_validate_result_v1(forged_attempt),
    class = "dpprior_schema_error"
  )
})


test_that("top-level class spoofing fails before forged accessors dispatch", {
  fit <- .schema23_fit()
  forged <- unclass(fit)
  class(forged) <- c("schema23_forged", "dpprior_result", "list")
  old <- get0("$.schema23_forged", envir = .GlobalEnv, inherits = FALSE)
  assign("$.schema23_forged", function(x, name) {
    stop("forged accessor dispatched")
  }, envir = .GlobalEnv)
  on.exit({
    if (is.null(old)) {
      rm("$.schema23_forged", envir = .GlobalEnv)
    } else {
      assign("$.schema23_forged", old, envir = .GlobalEnv)
    }
  }, add = TRUE)
  expect_error(
    .dpprior_validate_result_v1(forged), class = "dpprior_schema_error"
  )
})


test_that("validated nested targets cannot dispatch forged accessors", {
  fits <- lapply(
    c("a1_proxy", "a2_moment", "a2_kl", "dual_hard", "dual_soft",
      "dual_legacy"),
    .schema23_fit
  )
  method_names <- c(
    "[[.dpprior_K_target", "$.dpprior_K_target",
    "[[.dpprior_weight_target", "$.dpprior_weight_target"
  )
  old_methods <- lapply(
    method_names,
    get0, envir = .GlobalEnv, inherits = FALSE
  )
  names(old_methods) <- method_names
  for (method_name in method_names) {
    assign(method_name, function(...) {
      stop("forged nested target accessor dispatched")
    }, envir = .GlobalEnv)
  }
  on.exit({
    for (method_name in method_names) {
      old <- old_methods[[method_name]]
      if (is.null(old)) {
        rm(list = method_name, envir = .GlobalEnv)
      } else {
        assign(method_name, old, envir = .GlobalEnv)
      }
    }
  }, add = TRUE)

  for (fit in fits) {
    expect_invisible(.dpprior_validate_result_v1(fit))
  }
})


test_that("partial-name and malicious nested target forgeries are rejected", {
  target <- .schema23_target()
  for (value in list(environment(), quote(x + 1), c(20L, 20L),
                     matrix(20L, 1, 1), structure(20L, class = "forged"))) {
    bad <- target
    bad$request$J <- value
    expect_error(
      .dpprior_validate_target_v1(bad), class = "dpprior_schema_error"
    )
  }

  weight <- .schema23_weight_target("at_most", "wmax_tail_upper")
  weight$certification <- list(kind_fake = "upper_bound")
  expect_error(
    .dpprior_validate_weight_target_v1(weight),
    class = "dpprior_schema_error"
  )

  fit <- .schema23_fit()
  bad_provenance <- fit
  bad_provenance$provenance$approximation_fake <-
    bad_provenance$provenance$approximation
  bad_provenance$provenance$approximation <- NULL
  expect_error(
    .dpprior_validate_result_v1(bad_provenance),
    class = "dpprior_schema_error"
  )

  bad_verification <- fit
  bad_verification$verification$passed_fake <- TRUE
  bad_verification$verification$passed <- NULL
  expect_error(
    .dpprior_validate_result_v1(bad_verification),
    class = "dpprior_schema_error"
  )
})


test_that("target authority, kind, class, and verification are immutable", {
  target <- .schema23_target()

  verification_mismatch <- target
  verification_mismatch$status <- "approximate"
  verification_mismatch$usable <- FALSE
  verification_mismatch$verified <- FALSE
  expect_error(.dpprior_validate_target_v1(verification_mismatch),
               class = "dpprior_schema_error")

  snapshot_mismatch <- target
  snapshot_mismatch$verification$selected_snapshot$residuals$mean <- 0.1
  expect_error(.dpprior_validate_target_v1(snapshot_mismatch),
               class = "dpprior_schema_error")

  missing_interval <- target
  missing_interval$kind <- "interval"
  expect_error(.dpprior_validate_target_v1(missing_interval),
               class = "dpprior_schema_error")

  forged_class <- target
  class(forged_class) <- c("forged_target", class(target))
  expect_error(.dpprior_validate_target_v1(forged_class),
               class = "dpprior_schema_error")

  opaque_request <- target
  opaque_request$request$payload <- new.env(parent = emptyenv())
  expect_error(.dpprior_validate_target_v1(opaque_request),
               class = "dpprior_schema_error")

  explicit_null <- target
  explicit_null$request["var_K"] <- list(NULL)
  expect_error(.dpprior_validate_target_v1(explicit_null),
               class = "dpprior_schema_error")
  expect_identical(
    unserialize(serialize(target, NULL, version = 3L)), target
  )
})


test_that("weight targets bind canonical authority and certification", {
  weight <- .schema23_weight_target()
  changed_used <- weight
  changed_used$used$value <- 0.3
  expect_error(.dpprior_validate_weight_target_v1(changed_used),
               class = "dpprior_schema_error")

  changed_normalized <- weight
  changed_normalized$normalized$relation <- "at_most"
  expect_error(.dpprior_validate_weight_target_v1(changed_normalized),
               class = "dpprior_schema_error")

  certified <- .schema23_weight_target("at_most", "wmax_tail_upper")
  false_certificate <- certified
  false_certificate$certification$passed <- FALSE
  expect_error(.dpprior_validate_weight_target_v1(false_certificate),
               class = "dpprior_schema_error")

  forged_class <- weight
  class(forged_class) <- c("forged_weight", class(weight))
  expect_error(.dpprior_validate_weight_target_v1(forged_class),
               class = "dpprior_schema_error")

  opaque_used <- weight
  opaque_used$used$payload <- quote(system("unexpected"))
  expect_error(.dpprior_validate_weight_target_v1(opaque_used),
               class = "dpprior_schema_error")
})


test_that("hard and soft weight request canonicalization is exact and typed", {
  hard <- .schema23_transformed_weight_target("hard")
  soft <- .schema23_transformed_weight_target("soft")
  expect_invisible(.dpprior_validate_weight_target_v1(hard))
  expect_invisible(.dpprior_validate_weight_target_v1(soft))

  wrong_bound <- hard
  wrong_bound$request$bound <- 0.7
  wrong_bound$provenance$transformation$before <- wrong_bound$request
  expect_error(.dpprior_validate_weight_target_v1(wrong_bound),
               class = "dpprior_schema_error")

  wrong_relation <- hard
  wrong_relation$provenance$transformation$evidence$relation_to <- "at_least"
  expect_error(.dpprior_validate_weight_target_v1(wrong_relation),
               class = "dpprior_schema_error")

  wrong_mode <- soft
  wrong_mode$provenance$transformation <- list(
    rule = "canonicalize_soft_weight_target", opt_in = FALSE,
    before = soft$request, after = soft$normalized,
    evidence = list(
      mode = "hard", value_field = "value", relation_from = "target",
      relation_to = "target", probability_field = "none"
    )
  )
  expect_error(.dpprior_validate_weight_target_v1(wrong_mode),
               class = "dpprior_schema_error")

  pairlist_request <- hard
  pairlist_request$request <- as.pairlist(pairlist_request$request)
  expect_error(.dpprior_validate_weight_target_v1(pairlist_request),
               class = "dpprior_schema_error")
})


test_that("target transformations and PMFs cannot be silently repaired", {
  target <- .schema23_target()
  transformed <- target
  transformed$used$mu_K <- 5.1
  expect_error(
    .dpprior_validate_target_v1(transformed),
    class = "dpprior_schema_error"
  )

  target_snapshot <- target$verification$selected_snapshot
  expect_error(
    .dpprior_new_target_K(
      kind = "pmf", J = 2L,
      request = list(J = 2L, target_pmf = c(0.2, 0.7)),
      normalized = list(J = 2L, target_pmf = c(0.2, 0.7)),
      used = list(J = 2L, target_pmf = c(0.2, 0.7)),
      pmf = c(0.2, 0.7), status = "converged", usable = TRUE,
      verified = TRUE, computation = .schema23_computation("target_pmf"),
      verification = .dpprior_new_verification(
        "strict_pmf", TRUE, TRUE, selected_snapshot = target_snapshot,
        verifier_snapshot = target_snapshot,
        components = list(mass = TRUE), invariants = list(support = TRUE)
      ),
      provenance = .schema23_provenance(
        "target_pmf", parameterization = "none"
      )
    ),
    class = "dpprior_schema_error"
  )
})


test_that("attempt selection, fallback, and M evidence cross-check", {
  fit <- .schema23_fit()
  bad_id <- fit
  bad_id$computation$selected_attempt_id <- "missing"
  expect_error(.dpprior_validate_result_v1(bad_id),
               class = "dpprior_schema_error")

  missing_reason <- fit
  missing_reason$computation$attempts[[1L]]$elapsed_seconds <- NULL
  expect_error(.dpprior_validate_result_v1(missing_reason),
               class = "dpprior_schema_error")

  low_order <- fit
  low_order$computation$orders$M_verification_used <- 100L
  expect_error(.dpprior_validate_result_v1(low_order),
               class = "dpprior_schema_error")

  changed_method <- fit
  changed_method$computation$used$method <- "fallback"
  expect_error(.dpprior_validate_result_v1(changed_method),
               class = "dpprior_schema_error")

  no_selection <- fit
  no_selection$computation$attempts[[1L]]$selected <- FALSE
  no_selection$computation$selected_attempt_id <- NULL
  expect_error(.dpprior_validate_result_v1(no_selection),
               class = "dpprior_schema_error")

  nonzero_exit <- fit
  nonzero_exit$computation$attempts[[1L]]$exit_code <- 2L
  expect_error(.dpprior_validate_result_v1(nonzero_exit),
               class = "dpprior_schema_error")

  optimizer_without_attempt <- fit
  optimizer_without_attempt$computation$attempts <- list()
  optimizer_without_attempt$computation$selected_attempt_id <- NULL
  expect_error(.dpprior_validate_result_v1(optimizer_without_attempt),
               class = "dpprior_schema_error")

  boundary_without_reason <- fit
  boundary_without_reason$status <- "boundary"
  expect_error(.dpprior_validate_result_v1(boundary_without_reason),
               class = "dpprior_schema_error")

  contradictory_code <- fit
  contradictory_code$computation$termination$code <- "failed"
  expect_error(.dpprior_validate_result_v1(contradictory_code),
               class = "dpprior_schema_error")

  contradictory_iterations <- fit
  contradictory_iterations$computation$termination$iterations <- 6L
  expect_error(.dpprior_validate_result_v1(contradictory_iterations),
               class = "dpprior_schema_error")
})


test_that("attempt records reject malformed and dispatch-capable evidence", {
  fit <- .schema23_fit()
  mutations <- list(
    function(x) {
      x$computation$attempts[[1L]]$start <- matrix(c(1, 1), 1, 2)
      x
    },
    function(x) {
      x$computation$attempts[[1L]]$bounds <- structure(
        list(c(-1, -1), c(1, 1)), names = c("lower", "lower")
      )
      x
    },
    function(x) {
      class(x$computation$attempts[[1L]]$control) <- "forged_control"
      x
    },
    function(x) {
      x$computation$attempts[[1L]]$exit_code <- 0.5
      x
    },
    function(x) {
      x$computation$attempts[[1L]]$evaluations$function_count <- 1.5
      x
    },
    function(x) {
      class(x$computation$attempts) <- "forged_attempts"
      x
    },
    function(x) {
      x$computation$attempts <- pairlist(x$computation$attempts[[1L]])
      x
    },
    function(x) {
      x$computation$candidate_evaluations <- pairlist(
        x$computation$candidate_evaluations[[1L]]
      )
      x
    },
    function(x) {
      x$computation$attempts[[1L]]$error <- pairlist(
        class = "optimizer_error", code = "failed", message = "failure"
      )
      x
    },
    function(x) {
      x$computation$attempts[[1L]]$control <- pairlist(maxit = 100L)
      x
    },
    function(x) {
      x$computation$attempts[[1L]]$unavailable <- c(
        control = "declared unavailable despite retained value"
      )
      x
    },
    function(x) {
      x$computation$attempts[[1L]]$control <- NULL
      x
    }
  )
  for (mutate in mutations) {
    expect_error(
      .dpprior_validate_result_v1(mutate(fit)),
      class = "dpprior_schema_error"
    )
  }

  diagnostics <- .schema23_diagnostics()
  diagnostics$computation$attempts[[1L]]$selected <- FALSE
  diagnostics$computation["selected_attempt_id"] <- list(NULL)
  expect_invisible(.dpprior_validate_result_v1(diagnostics))
})


test_that("finite-support achieved K records are internally authoritative", {
  fit <- .schema23_fit()
  impossible_mean <- fit
  impossible_mean$achieved$K$mean <- 999
  expect_error(.dpprior_validate_result_v1(impossible_mean),
               class = "dpprior_schema_error")

  impossible_variance <- fit
  impossible_variance$achieved$K$variance <- 999
  expect_error(.dpprior_validate_result_v1(impossible_variance),
               class = "dpprior_schema_error")

  wrong_length <- fit
  wrong_length$achieved$K$pmf <- c(0.5, 0.5)
  expect_error(.dpprior_validate_result_v1(wrong_length),
               class = "dpprior_schema_error")

  wrong_moments <- fit
  wrong_moments$achieved$K$pmf <- c(rep(0, 4), 1, rep(0, 15))
  expect_error(.dpprior_validate_result_v1(wrong_moments),
               class = "dpprior_schema_error")
})


test_that("selected-order values and verifier candidate remain quarantined", {
  fit <- .schema23_fit()
  public_from_verifier <- fit
  public_from_verifier$achieved <- fit$verification$verifier_snapshot$achieved
  expect_error(
    .dpprior_validate_result_v1(public_from_verifier),
    class = "dpprior_schema_error"
  )

  changed_candidate <- fit
  changed_candidate$verification$verifier_snapshot$parameters$a <- 99
  expect_error(
    .dpprior_validate_result_v1(changed_candidate),
    class = "dpprior_schema_error"
  )

  changed_M <- fit
  changed_M$achieved$K$M <- 160L
  changed_M$verification$selected_snapshot$achieved$K$M <- 160L
  expect_error(
    .dpprior_validate_result_v1(changed_M),
    class = "dpprior_schema_error"
  )

  changed_stability <- fit
  changed_stability$verification$stability$delta[["K.mean"]] <- 0
  expect_error(
    .dpprior_validate_result_v1(changed_stability),
    class = "dpprior_schema_error"
  )

  malformed_verifier <- fit
  malformed_verifier$verification$verifier_snapshot$achieved$K$mean <-
    matrix(5, 1, 1)
  expect_error(
    .dpprior_validate_result_v1(malformed_verifier),
    class = "dpprior_schema_error"
  )
})


test_that("hard, soft, and legacy semantics cannot leak across modes", {
  hard <- .schema23_fit("dual_hard")
  hard$constraint$optimality$lambda <- 0.5
  expect_error(.dpprior_validate_result_v1(hard),
               class = "dpprior_schema_error")

  hard_check <- .schema23_fit("dual_hard")
  hard_check$verification$components$constraint_refined <- FALSE
  expect_error(.dpprior_validate_result_v1(hard_check),
               class = "dpprior_schema_error")

  soft <- .schema23_fit("dual_soft")
  soft$tradeoff$optimality$constraint <- list(satisfied = TRUE)
  expect_error(.dpprior_validate_result_v1(soft),
               class = "dpprior_schema_error")

  soft_loss <- .schema23_fit("dual_soft")
  soft_loss$tradeoff$total_loss <- 0.16
  expect_error(.dpprior_validate_result_v1(soft_loss),
               class = "dpprior_schema_error")

  soft_scale <- .schema23_fit("dual_soft")
  soft_scale$tradeoff$scales$K <- 2
  expect_error(.dpprior_validate_result_v1(soft_scale),
               class = "dpprior_schema_error")

  hard_residual <- .schema23_fit("dual_hard")
  hard_residual$constraint$residual <- 5e-7
  hard_residual$constraint$slack <- -5e-7
  expect_error(.dpprior_validate_result_v1(hard_residual),
               class = "dpprior_schema_error")

  soft_recomputed <- .schema23_fit("dual_soft")
  soft_recomputed$tradeoff$K_loss <- 0.2
  soft_recomputed$tradeoff$total_loss <- 0.1
  expect_error(.dpprior_validate_result_v1(soft_recomputed),
               class = "dpprior_schema_error")

  legacy <- .schema23_fit("dual_legacy")
  legacy$legacy$losses$constraint <- list(satisfied = TRUE)
  expect_error(.dpprior_validate_result_v1(legacy),
               class = "dpprior_schema_error")
})


test_that("legacy lambda retains zero while soft lambda remains open", {
  legacy_zero <- .schema23_legacy_fit(lambda = 0)
  expect_invisible(.dpprior_validate_result_v1(legacy_zero))

  expect_silent(.dpprior_new_legacy_details(
    contract = "DPprior_dual_v1", lambda = 0,
    losses = list(K = 0.2, weight = 0.1, total = 0.1),
    approximation_opt_in = TRUE,
    warning_code = "legacy_dual_approximation"
  ))

  for (lambda in c(-.Machine$double.eps, 1 + .Machine$double.eps)) {
    outside <- .schema23_fit("dual_legacy")
    outside$legacy$lambda <- lambda
    expect_error(
      .dpprior_validate_result_v1(outside), class = "dpprior_schema_error"
    )
  }

  soft_zero <- .schema23_fit("dual_soft")
  soft_zero$tradeoff$lambda <- 0
  expect_error(
    .dpprior_validate_result_v1(soft_zero), class = "dpprior_schema_error"
  )
})


test_that("dual_legacy accepts every frozen producer loss and lambda route", {
  for (loss_type in c("relative", "adaptive", "absolute")) {
    for (lambda in c(0, 0.5, 1)) {
      fit <- .schema23_legacy_fit(lambda = lambda, loss_type = loss_type)
      expect_invisible(.dpprior_validate_result_v1(fit))
      expect_identical(fit$legacy$lambda, lambda)
      expect_identical(fit$legacy$losses$loss_type, loss_type)
      expected_attempts <- if (identical(lambda, 1)) {
        0L
      } else if (identical(loss_type, "adaptive")) {
        2L
      } else {
        1L
      }
      expect_identical(length(fit$computation$attempts), expected_attempts)
    }
  }

  targets <- list(
    list(prob = list(threshold = 0.37, value = 0.3)),
    list(mean = 0.3),
    list(quantile = list(prob = 0.73, value = 0.4))
  )
  for (target in targets) {
    fit <- .schema23_legacy_fit(lambda = 1, weight_target = target)
    expect_invisible(.dpprior_validate_result_v1(fit))
  }
})


test_that("dual_legacy binds natural and forced fallback execution routes", {
  natural <- .schema23_legacy_fit(
    lambda = 0.5, loss_type = "adaptive", J = 50L,
    mu_K = 3, var_K = 10,
    weight_target = list(prob = list(threshold = 0.5, value = 0.25))
  )
  expect_true(natural$computation$fallback$used)
  expect_identical(
    vapply(natural$computation$attempts, `[[`, character(1), "stage"),
    c("scaling", "primary", "fallback")
  )
  expect_invisible(.dpprior_validate_result_v1(natural))
  for (attempt_index in seq_along(natural$computation$attempts)) {
    .schema23_expect_legacy_reject(natural, function(x) {
      x$computation$attempts[[attempt_index]]$candidate_objective <-
        x$computation$attempts[[attempt_index]]$candidate_objective + 1
      x
    })
  }

  mock_optim <- function(par, fn, method, control, ...) {
    if (identical(method, "BFGS")) {
      return(list(
        par = par, value = fn(par),
        counts = c("function" = 1, gradient = 1),
        convergence = 1L, message = "forced primary exit"
      ))
    }
    stats::optim(par, fn, method = method, control = control, ...)
  }
  testthat::local_mocked_bindings(
    .dpprior_legacy_dual_optim = mock_optim, .package = "DPprior"
  )
  forced <- .schema23_legacy_fit(lambda = 0.5, loss_type = "relative")
  expect_true(forced$computation$fallback$used)
  expect_identical(
    forced$computation$selected_attempt_id, "attempt-fallback-001"
  )
  expect_invisible(.dpprior_validate_result_v1(forced))
})


test_that("dual_legacy accepts a freshly bound A2-KL K-only baseline", {
  K_only <- DPprior_fit(
    J = 50L, mu_K = 5, var_K = 8, method = "A2-KL", M = 80L,
    check_diagnostics = TRUE
  )
  fit <- suppressWarnings(DPprior_dual(
    K_only, list(mean = 0.3), lambda = 1, loss_type = "adaptive", M = 80L
  ))
  expect_identical(
    fit$computation$resources$K_only_baseline$mode, "a2_kl"
  )
  expect_true(!is.null(
    fit$computation$resources$K_only_baseline$selected_snapshot$achieved_K$pmf
  ))
  expect_invisible(.dpprior_validate_result_v1(fit))
})


test_that("dual_legacy closes literals, provenance, settings, and resources", {
  base <- .schema23_legacy_fit(lambda = 0.5, loss_type = "relative")
  mutations <- list(
    function(x) { x$legacy$contract <- "forged"; x },
    function(x) { x$provenance$legacy$contract <- "forged"; x },
    function(x) { x$provenance$legacy$deprecation_stage <- "removed"; x },
    function(x) { x$legacy$approximation_opt_in <- FALSE; x },
    function(x) { x$legacy$warning_code <- "forged"; x },
    function(x) { x$provenance$approximation$active <- FALSE; x },
    function(x) { x$provenance$approximation$opt_in <- FALSE; x },
    function(x) { x$provenance$approximation$kind <- "forged"; x },
    function(x) { x$provenance$approximation$warning_code <- "forged"; x },
    function(x) { x$provenance$backend$implementation <- "forged"; x },
    function(x) { x$provenance$migration$adapter <- "forged"; x },
    function(x) { x$computation$request$method <- "forged"; x },
    function(x) { x$computation$used$parameterization <- "forged"; x },
    function(x) {
      forged_parameterization <- "log_ab"
      x$parameters$parameterization <- forged_parameterization
      x$provenance$parameterization <- forged_parameterization
      x$computation$request$parameterization <- forged_parameterization
      x$computation$used$parameterization <- forged_parameterization
      x$verification$selected_snapshot$parameters$parameterization <-
        forged_parameterization
      x$computation$resources$K_only_baseline$parameters$parameterization <-
        forged_parameterization
      x$computation$resources$K_only_baseline$selected_snapshot$
        parameters$parameterization <- forged_parameterization
      x$computation$attempts <- lapply(
        x$computation$attempts,
        function(attempt) {
          attempt$candidate_parameters$parameterization <-
            forged_parameterization
          attempt
        }
      )
      x
    },
    function(x) { x$computation$request$controls$lambda <- 0.7; x },
    function(x) { x$computation$used$controls$loss_type <- "absolute"; x },
    function(x) { x$computation$orders$M_requested <- 81L; x },
    function(x) { x$computation$orders$M_verification_used <- 160L; x },
    function(x) {
      x$computation$resources$optimizer_controls$primary$maxit <- 999L
      x
    },
    function(x) { x$computation$resources$objective$lambda <- 0.7; x },
    function(x) { x$legacy$lambda <- 0.7; x },
    function(x) {
      x$computation$resources$legacy_weight_request$prob$value <- 0.4
      x
    },
    function(x) {
      x$target$weight$provenance$legacy_request$prob$value <- 0.4
      x
    },
    function(x) { x$target$weight$request$value <- 0.4; x },
    function(x) {
      x$computation$resources$K_only_baseline$target$implied$mean <- 9
      x
    },
    function(x) {
      x$computation$resources$K_only_baseline$selected_snapshot$M <- 81L
      x
    },
    function(x) {
      x$computation$resources$K_only_baseline$selected_snapshot$achieved_K$mean <-
        9
      x
    },
    function(x) { x$verification$method <- "rubber_stamp"; x },
    function(x) { x$verification$components$forged <- TRUE; x },
    function(x) { x$message <- "forged"; x }
  )
  for (mutate in mutations) {
    .schema23_expect_legacy_reject(base, mutate)
  }

  overbound <- unserialize(serialize(base, NULL, xdr = TRUE))
  too_large <- as.numeric(floor(.Machine$integer.max / 2) + 1)
  overbound$computation$request$controls$max_iter <- too_large
  overbound$computation$used$controls$max_iter <- too_large
  overbound$computation$resources$optimizer_controls$primary$maxit <- too_large
  overbound$computation$resources$optimizer_controls$fallback$maxit <-
    2 * too_large
  expect_error(
    .dpprior_validate_result_v1(overbound), class = "dpprior_schema_error"
  )
})


test_that("dual_legacy freshly binds scaling, losses, and every attempt", {
  adaptive <- .schema23_legacy_fit(lambda = 0.5, loss_type = "adaptive")
  mutations <- list(
    function(x) { x$computation$scaling$formula <- "forged"; x },
    function(x) { x$computation$scaling$fixed_from_input <- TRUE; x },
    function(x) { x$computation$scaling$values$L_K_scale <- 999; x },
    function(x) { x$computation$scaling$values$L_w_scale <- 999; x },
    function(x) { x$legacy$losses$scaling$L_K_scale <- 999; x },
    function(x) { x$legacy$losses$K_loss <- 999; x },
    function(x) { x$legacy$losses$weight_loss <- 999; x },
    function(x) { x$legacy$losses$total_loss <- 999; x },
    function(x) { x$computation$resources$objective$K_loss <- 999; x },
    function(x) { x$computation$resources$objective$weight_loss <- 999; x },
    function(x) { x$computation$resources$objective$total_loss <- 999; x },
    function(x) { x$computation$attempts[[1L]]$stage <- "primary"; x },
    function(x) { x$computation$attempts[[1L]]$method <- "BFGS"; x },
    function(x) { x$computation$attempts[[1L]]$start[[1L]] <- 999; x },
    function(x) { x$computation$attempts[[1L]]$control$maxit <- 999L; x },
    function(x) { x$computation$attempts[[1L]]$candidate_objective <- 999; x },
    function(x) { x$computation$attempts[[2L]]$candidate_objective <- 999; x },
    function(x) { x$computation$attempts[[1L]]$iterations <- 999L; x },
    function(x) { x$computation$termination$source <- "legacy_adapter"; x },
    function(x) { x$provenance$is_fallback <- TRUE; x }
  )
  for (mutate in mutations) {
    .schema23_expect_legacy_reject(adaptive, mutate)
  }

  endpoint <- .schema23_legacy_fit(lambda = 1, loss_type = "adaptive")
  .schema23_expect_legacy_reject(endpoint, function(x) {
    x$computation$scaling$values$L_K_scale <- 1
    x$legacy$losses$scaling$L_K_scale <- 1
    x
  })
  .schema23_expect_legacy_reject(endpoint, function(x) {
    x$computation$termination$code <- "selected"
    x
  })

  coordinated <- unserialize(serialize(adaptive, NULL, xdr = TRUE))
  new_mean <- coordinated$achieved$K$mean + 0.1
  target_mean <- coordinated$target$K$implied$mean
  target_variance <- coordinated$target$K$implied$variance
  new_residual <- new_mean - target_mean
  variance_residual <- coordinated$residuals$K$variance
  raw_K <- (new_residual / target_mean)^2 +
    (variance_residual / target_variance)^2
  new_K_loss <- raw_K / coordinated$computation$scaling$values$L_K_scale
  new_total <- coordinated$legacy$lambda * new_K_loss +
    (1 - coordinated$legacy$lambda) *
      coordinated$legacy$losses$weight_loss
  coordinated$achieved$K$mean <- new_mean
  coordinated$residuals$K$mean <- new_residual
  coordinated$verification$selected_snapshot$achieved$K$mean <- new_mean
  coordinated$verification$selected_snapshot$residuals$K$mean <- new_residual
  coordinated$legacy$losses$K_loss <- new_K_loss
  coordinated$legacy$losses$total_loss <- new_total
  coordinated$computation$resources$objective$K_loss <- new_K_loss
  coordinated$computation$resources$objective$total_loss <- new_total
  selected_index <- match(
    coordinated$computation$selected_attempt_id,
    vapply(coordinated$computation$attempts, `[[`, character(1), "id")
  )
  coordinated$computation$attempts[[selected_index]]$candidate_objective <-
    new_total
  expect_error(
    .dpprior_validate_result_v1(coordinated), class = "dpprior_schema_error"
  )
})


test_that("dual_legacy compatibility is quarantined and serialization-safe", {
  fit <- .schema23_legacy_fit()
  restored <- unserialize(serialize(fit, NULL, xdr = TRUE))
  expect_identical(restored, fit)
  expect_invisible(.dpprior_validate_result_v1(restored))

  .schema23_expect_legacy_reject(fit, function(x) {
    x$compatibility$views$legacy_dual_v2$authority <- "authoritative"
    x
  })
  .schema23_expect_legacy_reject(fit, function(x) {
    x$compatibility$deprecations$legacy_dual_v2$code <- "forged"
    x
  })
  .schema23_expect_legacy_reject(fit, function(x) {
    x$compatibility$views$dual_anchor$consumer_policy <-
      "scientific_consumer_allowed"
    x
  })
  .schema23_expect_legacy_reject(fit, function(x) {
    x$compatibility$views$scientific_decision <- list(value = 1)
    x
  })
})


test_that("dual_legacy preserves only the exact finite v1.1 migration route", {
  migrated <- .schema23_migrated_dual(allow_legacy = TRUE, verify = FALSE)
  migrated_audit <- .schema23_migrated_dual(
    allow_legacy = TRUE, verify = TRUE
  )
  migrated_no_opt_in <- .schema23_migrated_dual(
    allow_legacy = FALSE, verify = FALSE
  )
  expect_invisible(.dpprior_validate_result_v1(migrated))
  expect_invisible(.dpprior_validate_result_v1(migrated_audit))
  expect_invisible(.dpprior_validate_result_v1(migrated_no_opt_in))
  expect_true(migrated$usable)
  expect_false(migrated_no_opt_in$usable)
  expect_identical(migrated$provenance$backend$package_version, "2.0.0")
  expect_identical(
    migrated$target$K$provenance$backend$package_version, "2.0.0"
  )
  expect_identical(
    migrated$compatibility$views$source$source_package_version, "1.1.0"
  )

  mutations <- list(
    function(x) { x$provenance$migration$source_schema <- "native"; x },
    function(x) { x$provenance$migration$adapter <- "none"; x },
    function(x) { x$provenance$migration$lossless <- TRUE; x },
    function(x) { x$provenance$migration$missing_evidence <- character(); x },
    function(x) { x$provenance$backend$package <- "forged"; x },
    function(x) { x$provenance$backend$package_version <- "1.1.0"; x },
    function(x) { x$provenance$backend$implementation <- "forged"; x },
    function(x) { x$provenance$backend$source_commit <- "forged"; x },
    function(x) {
      x$target$K$provenance$migration$source_schema <- "dpprior.result/1"
      x
    },
    function(x) {
      x$target$K$derivation$request_to_normalized$evidence$source_schema <-
        "dpprior.result/1"
      x
    },
    function(x) {
      x$target$weight$provenance$source_schema <- "dpprior.result/1"
      x
    },
    function(x) { x$provenance$legacy$active <- FALSE; x },
    function(x) { x$provenance$legacy$contract <- "forged"; x },
    function(x) { x$provenance$legacy$deprecation_stage <- "forged"; x },
    function(x) { x$provenance$approximation$kind <- "forged"; x },
    function(x) { x$provenance$approximation$warning_code <- "forged"; x },
    function(x) { x$legacy$contract <- "forged"; x },
    function(x) { x$legacy$warning_code <- "forged"; x },
    function(x) { x$verification$method <- "rubber_stamp"; x },
    function(x) { x$computation$termination$code <- "selected"; x },
    function(x) { x$computation$termination$source <- "optimizer"; x },
    function(x) { x$computation$attempts <- list(list()); x },
    function(x) { x$computation$candidate_evaluations <- list(list()); x },
    function(x) { x$computation$selected_attempt_id <- "attempt-1"; x },
    function(x) { x$computation$selected_candidate_id <- "candidate-1"; x },
    function(x) { x$provenance$input_fit <- list(); x },
    function(x) {
      x$compatibility$deprecations$legacy_schema$code <- "forged"
      x
    },
    function(x) {
      x$compatibility$deprecations$legacy_schema$removal_floor <- "removed"
      x
    },
    function(x) { x$compatibility$views$scientific_decision <- list(); x },
    function(x) { x$usable <- FALSE; x }
  )
  for (mutate in mutations) {
    .schema23_expect_legacy_reject(migrated, mutate)
  }

  native_relabel <- .schema23_legacy_fit()
  .schema23_expect_legacy_reject(native_relabel, function(x) {
    x$provenance$backend <- list(
      package = "DPprior", package_version = "1.1.0",
      implementation = "schema_upgrade_v1", source_commit = NULL
    )
    x$provenance$migration <- list(
      source_schema = "DPprior/1.1/fit",
      adapter = "upgrade_DPprior_object", lossless = FALSE,
      missing_evidence = c(
        "requested_controls", "M_requested", "M_selected",
        "M_verification_required", "M_verification_used", "scaling",
        "canonical_attempts", "candidate_selection", "fallback_lineage",
        "independent_verifier_snapshot", "source_commit",
        "legacy_objective_scaling", "legacy_optimizer_lineage"
      ),
      warnings = "legacy_object_upgraded"
    )
    x$provenance$legacy <- list(
      active = TRUE, contract = "v1.1_path_scaled_soft_equality_loss",
      deprecation_stage = "v2_migration"
    )
    x$provenance$approximation <- list(
      active = TRUE, opt_in = TRUE, kind = "legacy_schema_migration",
      warning_code = "legacy_object_upgraded"
    )
    x$legacy$contract <- "v1.1_path_scaled_soft_equality_loss"
    x$legacy$warning_code <- "legacy_object_upgraded"
    x$verification$method <- "legacy_evidence_quarantine"
    x$computation$attempts <- list()
    x$computation$candidate_evaluations <- list()
    x$computation[c("selected_candidate_id", "selected_attempt_id")] <-
      list(NULL, NULL)
    x$computation$termination$code <- "legacy_migration_no_selection"
    x$computation$termination$source <- "constructor"
    x$compatibility$deprecations <- list(legacy_schema = list(
      code = "legacy_object_upgraded",
      first_deprecated_version = "2.0.0",
      removal_floor = "not_scheduled"
    ))
    x
  })
})


test_that("migration keeps v1.1 source truth and a retained v2.x destination", {
  migrated <- list(
    dual_legacy = .schema23_migrated_dual(),
    a2_moment = .schema23_migrated_A2("a2_moment"),
    a2_kl = .schema23_migrated_A2("a2_kl")
  )
  invalid_versions <- list(
    "1.1.0", "0.9.0", "3.0.0", "2", "2.x.0", "2..0",
    "development", "v2.0.0", " 2.0.0", "2.0.0 ", "2.0.0\n", "",
    NA_character_,
    structure("2.1.0", class = "retained_package_version")
  )

  for (fit in migrated) {
    expect_identical(fit$provenance$backend$package_version, "2.0.0")
    expect_identical(
      fit$target$K$provenance$backend$package_version, "2.0.0"
    )
    expect_identical(
      fit$provenance$migration$source_schema, "DPprior/1.1/fit"
    )
    expect_identical(
      fit$compatibility$views$source$source_package_version, "1.1.0"
    )

    retained_later_v2 <- .schema23_set_migration_destination_version(
      fit, "2.12.3.9000"
    )
    expect_invisible(.dpprior_validate_result_v1(retained_later_v2))
    expect_invisible(.dpprior_validate_result_v1(unserialize(serialize(
      retained_later_v2, NULL, xdr = TRUE
    ))))

    for (version in invalid_versions) {
      forged <- .schema23_set_migration_destination_version(fit, version)
      expect_error(
        .dpprior_validate_result_v1(forged), class = "dpprior_schema_error"
      )
    }

    swapped <- .schema23_swap_migration_source_and_destination(fit)
    expect_error(
      .dpprior_validate_result_v1(swapped), class = "dpprior_schema_error"
    )
  }
})


test_that("dual_legacy truth binding never dispatches nested target accessors", {
  fits <- list(
    .schema23_legacy_fit(
      lambda = 0.5,
      weight_target = list(prob = list(threshold = 0.37, value = 0.3))
    ),
    .schema23_legacy_fit(lambda = 1, weight_target = list(mean = 0.3)),
    .schema23_legacy_fit(
      lambda = 0,
      weight_target = list(quantile = list(prob = 0.73, value = 0.4))
    )
  )
  method_names <- c(
    "[[.dpprior_K_target", "$.dpprior_K_target",
    "[[.dpprior_weight_target", "$.dpprior_weight_target"
  )
  old_methods <- lapply(
    method_names, get0, envir = .GlobalEnv, inherits = FALSE
  )
  names(old_methods) <- method_names
  dispatch_count <- 0L
  hostile <- function(...) {
    dispatch_count <<- dispatch_count + 1L
    stop("hostile target accessor dispatched")
  }
  for (method_name in method_names) {
    assign(method_name, hostile, envir = .GlobalEnv)
  }
  on.exit({
    for (method_name in method_names) {
      old <- old_methods[[method_name]]
      if (is.null(old)) {
        rm(list = method_name, envir = .GlobalEnv)
      } else {
        assign(method_name, old, envir = .GlobalEnv)
      }
    }
  }, add = TRUE)
  for (fit in fits) {
    expect_invisible(.dpprior_validate_result_v1(fit))
  }
  expect_identical(dispatch_count, 0L)
})


test_that("mode, verification evidence, and provenance claims are closed", {
  mismatched <- .schema23_fit("a2_kl")
  mismatched$method <- "A1"
  mismatched$computation$request$method <- "A1"
  mismatched$computation$used$method <- "A1"
  mismatched$provenance$requested_method <- "A1"
  mismatched$provenance$selected_method <- "A1"
  expect_error(.dpprior_validate_result_v1(mismatched),
               class = "dpprior_schema_error")

  for (value in list(matrix("A2-MN", 1, 1),
                     structure("A2-MN", class = "forged_method"))) {
    bad <- .schema23_fit()
    bad$method <- value
    expect_error(.dpprior_validate_result_v1(bad),
                 class = "dpprior_schema_error")
  }

  empty_evidence <- .schema23_fit()
  empty_evidence$verification$components <- list()
  expect_error(.dpprior_validate_result_v1(empty_evidence),
               class = "dpprior_schema_error")

  missing_component <- .schema23_fit("a2_kl")
  missing_component$verification$components$pmf_adequacy <- NULL
  expect_error(.dpprior_validate_result_v1(missing_component),
               class = "dpprior_schema_error")

  bare_decision <- .schema23_fit()
  bare_decision$verification$components$residual_adequacy <- TRUE
  expect_error(.dpprior_validate_result_v1(bare_decision),
               class = "dpprior_schema_error")

  false_decision <- .schema23_fit()
  false_decision$verification$components$residual_adequacy$value <- 999
  false_decision$verification$components$residual_adequacy$tolerance <- 0
  expect_error(.dpprior_validate_result_v1(false_decision),
               class = "dpprior_schema_error")

  soft_missing <- .schema23_fit("dual_soft")
  soft_missing$verification$components$local_optimality <- NULL
  expect_error(.dpprior_validate_result_v1(soft_missing),
               class = "dpprior_schema_error")

  hard_without_input <- .schema23_fit("dual_hard")
  hard_without_input$provenance["input_fit"] <- list(NULL)
  expect_error(.dpprior_validate_result_v1(hard_without_input),
               class = "dpprior_schema_error")

  soft_bad_input <- .schema23_fit("dual_soft")
  soft_bad_input$provenance$input_fit$status <- "failed"
  expect_error(.dpprior_validate_result_v1(soft_bad_input),
               class = "dpprior_schema_error")

  false_lossless <- .schema23_fit()
  false_lossless$provenance$migration$warnings <- "evidence missing"
  expect_error(.dpprior_validate_result_v1(false_lossless),
               class = "dpprior_schema_error")

  opaque_control <- .schema23_fit()
  opaque_control$computation$used$controls$payload <-
    new.env(parent = emptyenv())
  expect_error(.dpprior_validate_result_v1(opaque_control),
               class = "dpprior_schema_error")

  null_control <- .schema23_fit()
  null_control$computation$request$controls["optional"] <- list(NULL)
  expect_invisible(.dpprior_schema_validate_plain_record_value(
    null_control$computation$request$controls,
    "result.computation.request.controls"
  ))
  expect_error(.dpprior_validate_result_v1(null_control),
               class = "dpprior_schema_error")

  used_unknown <- null_control
  used_unknown$computation$used$controls["optional"] <- list(NULL)
  expect_error(.dpprior_validate_result_v1(used_unknown),
               class = "dpprior_schema_error")
})


test_that("verification check records reject shape and S3 forgeries", {
  fit <- .schema23_fit()
  forged_checks <- list(
    matrix(TRUE, 1, 1),
    structure(list(passed = TRUE), class = c("forged_check", "list")),
    structure(list(TRUE, FALSE), names = c("passed", "passed")),
    structure(list(passed = TRUE), dim = c(1L))
  )
  old <- get0("$.forged_check", envir = .GlobalEnv, inherits = FALSE)
  assign("$.forged_check", function(x, name) {
    stop("forged accessor dispatched")
  }, envir = .GlobalEnv)
  on.exit({
    if (is.null(old)) {
      rm("$.forged_check", envir = .GlobalEnv)
    } else {
      assign("$.forged_check", old, envir = .GlobalEnv)
    }
  }, add = TRUE)
  for (check in forged_checks) {
    bad <- fit
    bad$verification$components$residual_adequacy <- check
    expect_error(.dpprior_validate_result_v1(bad),
                 class = "dpprior_schema_error")
  }
})


test_that("hard status and constraint evidence cannot contradict", {
  unsatisfied_success <- .schema23_fit("dual_hard")
  unsatisfied_success$constraint$residual <- 0.1
  unsatisfied_success$constraint$slack <- -0.1
  unsatisfied_success$constraint$satisfied <- FALSE
  expect_error(.dpprior_validate_result_v1(unsatisfied_success),
               class = "dpprior_schema_error")

  failed_satisfied <- .schema23_fit("dual_hard")
  failed_satisfied$status <- "failed"
  failed_satisfied$usable <- FALSE
  failed_satisfied$verified <- FALSE
  failed_satisfied["parameters"] <- list(NULL)
  failed_satisfied$verification$passed <- FALSE
  failed_satisfied$verification$reason <- "failed"
  failed_satisfied$computation$attempts[[1L]]$selected <- FALSE
  failed_satisfied$computation["selected_attempt_id"] <- list(NULL)
  failed_satisfied$computation$termination$code <- "failed"
  failed_satisfied$computation$attempts[[1L]]$exit_code <- 1L
  expect_error(.dpprior_validate_result_v1(failed_satisfied),
               class = "dpprior_schema_error")
})


test_that("hard and soft no-candidate states retain typed unavailability", {
  infeasible <- .schema23_no_candidate("dual_hard", "infeasible")
  expect_null(infeasible$parameters)
  expect_length(infeasible$achieved, 0L)
  expect_null(infeasible$constraint$residual)
  expect_identical(
    infeasible$constraint$feasibility$classification,
    "certified_infeasible"
  )
  expect_invisible(.dpprior_validate_result_v1(infeasible))

  unknown <- .schema23_no_candidate("dual_hard", "failed")
  expect_true(unknown$constraint$feasibility$feasibility_unknown)
  expect_null(unknown$constraint$satisfied)
  expect_invisible(.dpprior_validate_result_v1(unknown))

  soft_failed <- .schema23_no_candidate("dual_soft", "failed")
  expect_null(soft_failed$tradeoff$K_loss)
  expect_identical(
    soft_failed$tradeoff$optimality$unavailable_reason,
    "no finite optimizer candidate"
  )
  expect_invisible(.dpprior_validate_result_v1(soft_failed))

  fake_hard <- unknown
  fake_hard$constraint$residual <- 0
  expect_error(.dpprior_validate_result_v1(fake_hard),
               class = "dpprior_schema_error")

  fake_soft <- soft_failed
  fake_soft$tradeoff$K_loss <- 0
  expect_error(.dpprior_validate_result_v1(fake_soft),
               class = "dpprior_schema_error")
})


test_that("compatibility aliases are derived and identity checked", {
  fit <- .schema23_fit()
  aliased <- .dpprior_append_compatibility_v2(
    fit,
    aliases = c(
      a = "parameters.a",
      b = "parameters.b",
      achieved_K = "achieved.K",
      converged = "compatibility.views.converged"
    ),
    views = list(converged = TRUE),
    deprecations = list(a = list(removal_floor = "3.0.0"))
  )
  expect_identical(aliased$a, fit$parameters$a)
  expect_identical(aliased$achieved_K, fit$achieved$K)
  expect_invisible(.dpprior_validate_result_v1(aliased))

  diverged <- aliased
  diverged$a <- 999
  expect_error(
    .dpprior_validate_result_v1(diverged), class = "dpprior_schema_error"
  )
  expect_error(
    .dpprior_append_compatibility_v2(
      fit, aliases = c(method = "parameters.a")
    ),
    class = "dpprior_schema_error"
  )

  boundary <- .schema23_fit()
  boundary$status <- "boundary"
  boundary$computation$termination$boundary_reason <- "fixture_boundary"
  boundary_alias <- .dpprior_append_compatibility_v2(
    boundary,
    aliases = c(converged = "compatibility.views.converged"),
    views = list(converged = TRUE)
  )
  expect_true(boundary_alias$converged)

  infeasible <- .schema23_no_candidate("dual_hard", "infeasible")
  infeasible_alias <- .dpprior_append_compatibility_v2(
    infeasible,
    aliases = c(converged = "compatibility.views.converged"),
    views = list(converged = FALSE)
  )
  expect_false(infeasible_alias$converged)
})


test_that("collect mode reports violations without signaling", {
  fit <- .schema23_fit()
  valid <- .dpprior_validate_object(fit, collect = TRUE)
  expect_true(valid$valid)
  expect_length(valid$errors, 0L)

  fit$method <- "wrong"
  invalid <- .dpprior_validate_object(fit, collect = TRUE)
  expect_false(invalid$valid)
  expect_length(invalid$errors, 1L)
  expect_s3_class(invalid$errors[[1L]], "dpprior_schema_error")
})


test_that("hard verification binds refined residuals to the central constraint", {
  hard <- .schema23_fit("dual_hard")
  hard$verification$verifier_snapshot$achieved$weight$value <- 0.9
  hard$tolerances$weight$absolute <- 1
  hard$verification$selected_snapshot$tolerances$weight$absolute <- 1
  hard$verification$verifier_snapshot$tolerances$weight$absolute <- 1
  hard$constraint$tolerance$absolute <- 1
  hard$constraint$tolerance$effective <- 1
  hard$verification$stability <- .dpprior_new_stability(
    delta = c(K.mean = 1e-10, K.variance = 0, weight.value = 0.5),
    tolerance = c(K.mean = 6e-6, K.variance = 3e-6, weight.value = 1),
    formula = c(
      K.mean = "absolute_plus_relative_max",
      K.variance = "absolute_plus_relative_max",
      weight.value = "absolute_plus_relative_max"
    ),
    scale_floor = c(K.mean = 1e-8, K.variance = 1e-8,
                    weight.value = 1e-8),
    source = "forged_loose_stability"
  )
  hard$verification$components$order_stability <- .schema23_check(
    value = 0.5, tolerance = 1, operator = "lte"
  )
  hard$verification$components$constraint_selected <- .schema23_check(
    value = 0, tolerance = 1, operator = "lte"
  )
  hard$verification$components$constraint_refined <- .schema23_check(
    value = 0, tolerance = 1, operator = "lte"
  )
  expect_error(.dpprior_validate_result_v1(hard),
               class = "dpprior_schema_error")
})


test_that("soft optimality cannot contradict a verified scientific claim", {
  soft <- .schema23_fit("dual_soft")
  soft$tradeoff$optimality$gradient[] <- c(999, -999)
  expect_error(.dpprior_validate_result_v1(soft),
               class = "dpprior_schema_error")

  soft <- .schema23_fit("dual_soft")
  soft$tradeoff$optimality$local_minimum_passed <- FALSE
  expect_error(.dpprior_validate_result_v1(soft),
               class = "dpprior_schema_error")
})


test_that("soft lambda-one endpoint is exact input-fit identity", {
  endpoint <- .schema23_soft_endpoint()
  expect_invisible(.dpprior_validate_result_v1(endpoint))
  expect_true(endpoint$tradeoff$endpoint)
  expect_false(endpoint$tradeoff$optimality$performed)
  expect_false(endpoint$tradeoff$optimality$passed)
  expect_length(endpoint$computation$attempts, 0L)
  expect_length(endpoint$computation$candidate_evaluations, 0L)

  wrong_input <- endpoint
  wrong_input$provenance$input_fit$parameters <-
    .dpprior_new_parameters(2.1, 3, "log_ab")
  expect_error(.dpprior_validate_result_v1(wrong_input),
               class = "dpprior_schema_error")

  arbitrary_identity <- endpoint
  arbitrary_identity$provenance$input_fit$identity <- "trust-me"
  expect_error(.dpprior_validate_result_v1(arbitrary_identity),
               class = "dpprior_schema_error")

  wrong_method <- endpoint
  wrong_method$provenance$input_fit$method <- "A2-KL"
  expect_error(.dpprior_validate_result_v1(wrong_method),
               class = "dpprior_schema_error")

  wrong_J <- endpoint
  wrong_J$provenance$input_fit$J <- 19L
  wrong_J$provenance$input_fit$target$J <- 19L
  expect_error(.dpprior_validate_result_v1(wrong_J),
               class = "dpprior_schema_error")

  wrong_target <- endpoint
  wrong_target$provenance$input_fit$target$implied$mean <- 6
  expect_error(.dpprior_validate_result_v1(wrong_target),
               class = "dpprior_schema_error")

  wrong_snapshot <- endpoint
  wrong_snapshot$provenance$input_fit$selected_snapshot$achieved_K$mean <- 6
  expect_error(.dpprior_validate_result_v1(wrong_snapshot),
               class = "dpprior_schema_error")

  wrong_code <- endpoint
  wrong_code$computation$termination$code <- "closed_form"
  expect_error(.dpprior_validate_result_v1(wrong_code),
               class = "dpprior_schema_error")

  wrong_source <- endpoint
  wrong_source$computation$termination$source <- "constructor"
  expect_error(.dpprior_validate_result_v1(wrong_source),
               class = "dpprior_schema_error")

  active_optimality <- endpoint
  active_optimality$tradeoff$optimality$performed <- TRUE
  expect_error(.dpprior_validate_result_v1(active_optimality),
               class = "dpprior_schema_error")
})


test_that("failed soft lambda-one endpoint is an exact no-candidate route", {
  no_candidate <- .schema23_soft_failed_endpoint()
  backend_error <- .schema23_soft_failed_endpoint(
    "soft path backend failed: boom"
  )
  for (failed in list(no_candidate, backend_error)) {
    expect_invisible(.dpprior_validate_result_v1(failed))
    expect_identical(failed$status, "failed")
    expect_false(failed$usable)
    expect_false(failed$verified)
    expect_null(failed$parameters)
    expect_true(failed$tradeoff$endpoint)
    expect_identical(failed$tradeoff$lambda, 1)
    expect_length(failed$computation$attempts, 0L)
    expect_length(failed$computation$candidate_evaluations, 0L)
    expect_null(failed$computation$selected_attempt_id)
    expect_null(failed$computation$selected_candidate_id)
    expect_identical(failed$computation$termination$code, "no_candidate")
    expect_identical(failed$computation$termination$source, "no_candidate")
    expect_identical(failed$verification$method, "no_candidate")
    expect_false(failed$verification$performed)
    expect_false(failed$verification$passed)
    expect_length(failed$verification$components, 0L)
    expect_null(failed$verification$selected_snapshot)
    expect_null(failed$verification$verifier_snapshot)
    expect_null(failed$verification$stability)
  }

  finite_endpoint <- .schema23_soft_endpoint()
  expect_invisible(.dpprior_validate_result_v1(finite_endpoint))

  attacks <- list(
    function(x) {
      x$verification$components$endpoint_input_identity <- .schema23_check(
        value = TRUE, reference = TRUE, tolerance = NULL,
        operator = "identical", source = "independent_verifier"
      )
      x
    },
    function(x) {
      x$verification$stability <- finite_endpoint$verification$stability
      x
    },
    function(x) {
      x$verification$selected_snapshot <-
        finite_endpoint$verification$selected_snapshot
      x
    },
    function(x) {
      x$verification <- finite_endpoint$verification
      x
    },
    function(x) {
      x$parameters <- finite_endpoint$parameters
      x
    },
    function(x) {
      x$computation$attempts <-
        .schema23_fit("dual_soft")$computation$attempts
      x
    },
    function(x) {
      x$computation$candidate_evaluations <-
        .schema23_fit("dual_soft")$computation$candidate_evaluations
      x
    },
    function(x) {
      x$computation$selected_attempt_id <- "attempt-1"
      x
    },
    function(x) {
      x$computation$selected_candidate_id <- "candidate-1"
      x
    },
    function(x) {
      x$computation$termination$code <- "endpoint"
      x$computation$termination$source <- "endpoint"
      x
    },
    function(x) {
      x$computation$termination$code <- "endpoint"
      x
    },
    function(x) {
      x$computation$termination$source <- "endpoint"
      x
    },
    function(x) {
      x$tradeoff$optimality$performed <- TRUE
      x
    },
    function(x) {
      x$tradeoff$endpoint <- FALSE
      x
    },
    function(x) {
      x$status <- "approximate"
      x
    }
  )
  for (mutate in attacks) {
    expect_error(
      .dpprior_validate_result_v1(mutate(no_candidate)),
      class = "dpprior_schema_error"
    )
  }
})


test_that("A2 moment and KL claims are recomputed from target authority", {
  moment <- .schema23_fit("a2_moment")
  moment$achieved$K$mean <- 10
  moment$residuals$K$mean <- 5
  moment$verification$selected_snapshot$achieved$K$mean <- 10
  moment$verification$selected_snapshot$residuals$K$mean <- 5
  moment$verification$verifier_snapshot$achieved$K$mean <- 10 + 1e-10
  moment$verification$verifier_snapshot$residuals$K$mean <- 5 + 1e-10
  moment$verification$stability <- .dpprior_new_stability(
    delta = c(K.mean = 1e-10, K.variance = 0),
    tolerance = c(K.mean = 11e-6 + 1e-16, K.variance = 3e-6),
    formula = c(
      K.mean = "absolute_plus_relative_max",
      K.variance = "absolute_plus_relative_max"
    ),
    scale_floor = c(K.mean = 1, K.variance = 1),
    source = "forged_self_consistent"
  )
  moment$verification$components$order_stability <- .schema23_check(
    value = moment$verification$stability$delta,
    reference = c(K.mean = 0, K.variance = 0),
    tolerance = moment$verification$stability$tolerance,
    operator = "lte"
  )
  expect_error(.dpprior_validate_result_v1(moment),
               class = "dpprior_schema_error")

  missing_objective <- .schema23_fit("a2_kl")
  missing_objective$target$K <- .schema23_target()
  expect_error(.dpprior_validate_result_v1(missing_objective),
               class = "dpprior_schema_error")

  chisq <- .schema23_fit("a2_kl")
  objective_df <- 25
  objective_scale <- 0.2
  lower_edges <- (1:20 - 0.5) / objective_scale
  upper_edges <- (1:20 + 0.5) / objective_scale
  raw_mass <- stats::pchisq(upper_edges, df = objective_df) -
    stats::pchisq(lower_edges, df = objective_df)
  retained_mass <- sum(raw_mass)
  objective_pmf <- raw_mass / retained_mass
  objective_mean <- sum((1:20) * objective_pmf)
  objective_variance <- sum(((1:20) - objective_mean)^2 * objective_pmf)
  moment_target <- .schema23_target()
  moment_target$derivation$request_to_normalized$evidence$A2_KL_objective <-
    list(
      method = "chisq", df = objective_df, scale = objective_scale,
      binning = "continuity_corrected_half_integer_bins",
      normalization = "explicit_support_conditioning",
      support = c(lower = 1L, upper = 20L),
      retained_mass_before_normalization = retained_mass,
      omitted_mass = max(0, 1 - retained_mass), pmf = objective_pmf,
      mu_K_discrete = objective_mean,
      var_K_discrete = objective_variance,
      source = "A2_KL_backend_target"
    )
  expect_invisible(.dpprior_validate_target_v1(moment_target))
  chisq$target$K <- moment_target
  chisq$achieved$K$pmf <- objective_pmf
  chisq$achieved$K$mean <- objective_mean
  chisq$achieved$K$variance <- objective_variance
  chisq_mean_residual <- objective_mean - 5
  chisq_variance_residual <- objective_variance - 2
  chisq$residuals$K <- list(
    mean = chisq_mean_residual, variance = chisq_variance_residual
  )
  chisq$residuals$distribution <- list(
    kl = 0, l1 = 0, mean = 0, variance = 0
  )
  for (snapshot_name in c("selected_snapshot", "verifier_snapshot")) {
    chisq$verification[[snapshot_name]]$achieved$K <- chisq$achieved$K
    chisq$verification[[snapshot_name]]$achieved$K$M <-
      chisq$verification[[snapshot_name]]$M
    chisq$verification[[snapshot_name]]$residuals <- chisq$residuals
  }
  chisq$verification$components$target_identity <- .schema23_check(
    value = objective_pmf, reference = objective_pmf,
    tolerance = NULL, operator = "identical",
    source = "independent_verifier"
  )
  expected_chisq_stability <- .dpprior_expected_result_stability(chisq)
  chisq$verification$stability <- .dpprior_new_stability(
    delta = expected_chisq_stability$delta,
    tolerance = expected_chisq_stability$tolerance,
    formula = expected_chisq_stability$formula,
    scale_floor = expected_chisq_stability$scale_floor,
    source = "independent_verifier"
  )
  chisq$verification$components$order_stability <- .schema23_check(
    value = expected_chisq_stability$delta,
    reference = setNames(
      rep(0, length(expected_chisq_stability$delta)),
      names(expected_chisq_stability$delta)
    ),
    tolerance = expected_chisq_stability$tolerance,
    operator = "lte", source = "independent_verifier"
  )
  chisq_candidate <- chisq$computation$candidate_evaluations[[1L]]
  chisq_candidate$selected_snapshot <- chisq$verification$selected_snapshot
  chisq_candidate$checks$candidate_distribution <- .schema23_check(
    value = c(
      pmf_mass_error = abs(sum(objective_pmf) - 1),
      pmf_minimum_violation = max(0, -min(objective_pmf))
    ),
    reference = c(pmf_mass_error = 0, pmf_minimum_violation = 0),
    tolerance = c(
      pmf_mass_error = .TOL_PMF_SUM, pmf_minimum_violation = 0
    ),
    operator = "lte", source = "candidate:candidate-1"
  )
  chisq$computation$candidate_evaluations[[1L]] <- chisq_candidate
  expect_invisible(.dpprior_validate_result_v1(chisq))

  same_moments <- chisq
  null_direction <- c(-1, 3, -3, 1)
  epsilon <- min(
    objective_pmf[[4L]], objective_pmf[[6L]] / 3
  ) / 2
  same_moments$target$K$derivation$request_to_normalized$evidence$
    A2_KL_objective$pmf[4:7] <-
    objective_pmf[4:7] + epsilon * null_direction
  expect_error(.dpprior_validate_result_v1(same_moments),
               class = "dpprior_schema_error")

  chisq$target$K$derivation$request_to_normalized$evidence$
    A2_KL_objective$df <- 999
  expect_error(.dpprior_validate_result_v1(chisq),
               class = "dpprior_schema_error")
})


test_that("target derivations and verifier moments cannot silently change", {
  target <- .schema23_target()
  target$normalized$mu_K <- 10
  target$normalized$var_K <- 3
  target$used <- target$normalized
  target$implied <- list(mean = 10, variance = 3)
  target$verification$selected_snapshot$achieved$implied <- target$implied
  target$verification$verifier_snapshot$achieved$implied <- target$implied
  target$derivation$request_to_normalized$after <- target$normalized
  expect_error(.dpprior_validate_target_v1(target),
               class = "dpprior_schema_error")

  forged_verifier <- .schema23_target()
  forged_verifier$verification$verifier_snapshot$achieved$implied$mean <- 999
  expect_error(.dpprior_validate_target_v1(forged_verifier),
               class = "dpprior_schema_error")
})


test_that("each scientific derivation rule recomputes its retained evidence", {
  make_confidence <- function(level, vif, variance = vif * 4) {
    object <- unclass(.schema23_target())
    object$request <- list(J = 20L, mean = 5, confidence = level)
    object$normalized <- c(
      object$request, list(interval = NULL, pmf = NULL)
    )
    object$used <- list(
      J = 20L, mean = 5, variance = variance,
      interval = NULL, pmf = NULL
    )
    object$derivation <- list(
      request_to_normalized = list(
        rule = "canonicalize_confidence_target", outcome = "canonicalized",
        opt_in = FALSE, before = object$request,
        after = object$normalized, evidence = list(source = "fixture")
      ),
      normalized_to_used = list(
        rule = "derive_variance_from_confidence_vif", outcome = "derived",
        opt_in = FALSE, before = object$normalized, after = object$used,
        evidence = list(
          confidence = level, vif = vif,
          formula = "variance = vif * (mean - 1)"
        )
      )
    )
    object
  }
  for (level in names(c(low = 5, medium = 2.5, high = 1.5))) {
    vif <- c(low = 5, medium = 2.5, high = 1.5)[[level]]
    expect_invisible(.dpprior_validate_target_derivation(
      make_confidence(level, vif)
    ))
  }
  for (bad_level in list(0.9, "unknown", structure("medium", class = "evil"))) {
    expect_error(
      .dpprior_validate_target_derivation(
        make_confidence(bad_level, 2.5)
      ),
      class = "dpprior_schema_error"
    )
  }

  confidence <- unclass(.schema23_target())
  confidence$request <- list(J = 20L, mean = 5, confidence = "medium")
  confidence$normalized <- c(
    confidence$request, list(interval = NULL, pmf = NULL)
  )
  confidence$used <- list(
    J = 20L, mean = 5, variance = 999, interval = NULL, pmf = NULL
  )
  confidence$derivation <- list(
    request_to_normalized = list(
      rule = "canonicalize_confidence_target", outcome = "canonicalized",
      opt_in = FALSE, before = confidence$request,
      after = confidence$normalized, evidence = list(source = "fixture")
    ),
    normalized_to_used = list(
      rule = "derive_variance_from_confidence_vif", outcome = "derived",
      opt_in = FALSE, before = confidence$normalized,
      after = confidence$used,
      evidence = list(
        confidence = "medium", vif = 2.5,
        formula = "variance = vif * (mean - 1)"
      )
    )
  )
  expect_error(.dpprior_validate_target_derivation(confidence),
               class = "dpprior_schema_error")

  cv <- unclass(.schema23_target())
  cv$kind <- "cv"
  cv$request <- list(J = 20L, mean = 5, cv = 0.2)
  cv$normalized <- c(cv$request, list(interval = NULL, pmf = NULL))
  cv$used <- list(
    J = 20L, mean = 5, variance = 9, interval = NULL, pmf = NULL
  )
  cv$derivation <- list(
    request_to_normalized = list(
      rule = "canonicalize_cv_target", outcome = "canonicalized",
      opt_in = FALSE, before = cv$request, after = cv$normalized,
      evidence = list(source = "fixture")
    ),
    normalized_to_used = list(
      rule = "derive_variance_from_cv", outcome = "derived",
      opt_in = FALSE, before = cv$normalized, after = cv$used,
      evidence = list(
        cv = 0.2, definition = "SD(K_J) / E(K_J)",
        formula = "variance = (cv * mean)^2"
      )
    )
  )
  expect_error(.dpprior_validate_target_derivation(cv),
               class = "dpprior_schema_error")

  interval <- unclass(.schema23_target("pmf"))
  interval$kind <- "interval"
  interval$interval <- list(
    lower = 4L, upper = 6L, type = "equal_tail", coverage = 0.5
  )
  interval$family <- list(name = "maxent")
  interval$request <- list(J = 20L, interval = interval$interval)
  interval$normalized <- list(
    J = 20L, interval = interval$interval, family = interval$family,
    pmf = NULL
  )
  interval$used <- list(
    J = 20L, interval = interval$interval, family = interval$family,
    pmf = interval$pmf
  )
  interval$derivation <- list(
    request_to_normalized = list(
      rule = "canonicalize_interval_request", outcome = "canonicalized",
      opt_in = FALSE, before = interval$request,
      after = interval$normalized, evidence = list(source = "fixture")
    ),
    normalized_to_used = list(
      rule = "construct_maxent_equal_tail_pmf", outcome = "derived",
      opt_in = FALSE, before = interval$normalized, after = interval$used,
      evidence = list(
        constructor_method = "analytic", group_masses = list(
          left = 0, inside = 1, right = 0
        ), common_tilt = 0, boundary = list(active = FALSE),
        mean_constraint = 5, feasible_mean_lower = 4,
        feasible_mean_upper = 6, feasibility_tolerance = 1e-8,
        root_tolerance = 1e-10
      )
    )
  )
  expect_error(.dpprior_validate_target_derivation(interval),
               class = "dpprior_schema_error")

  family <- unclass(.schema23_target("pmf"))
  family$kind <- "family"
  family$family <- list(name = "scaled_chisq")
  family$request <- list(J = 20L, mean = 5, variance = 2)
  family$normalized <- list(
    J = 20L, mean = 5, variance = 2, interval = NULL,
    family = family$family, pmf = NULL
  )
  family$used <- list(
    J = 20L, interval = NULL, family = family$family, pmf = family$pmf
  )
  family$derivation <- list(
    request_to_normalized = list(
      rule = "canonicalize_scaled_chisq_family_request",
      outcome = "canonicalized", opt_in = FALSE,
      before = family$request, after = family$normalized,
      evidence = list(source = "fixture")
    ),
    normalized_to_used = list(
      rule = "construct_scaled_chisq_conditioned_pmf",
      outcome = "derived", opt_in = FALSE,
      before = family$normalized, after = family$used,
      evidence = list(
        df = 999, scale = 0.2, binning = "continuity_corrected",
        retained_mass_before_normalization = 0.9, omitted_mass = 0.1,
        normalization = "condition_on_1_to_J", support = 1:20
      )
    )
  )
  expect_error(.dpprior_validate_target_derivation(family),
               class = "dpprior_schema_error")

  projection <- unclass(.schema23_target())
  projection$used$var_K <- 3
  projection$derivation$normalized_to_used <- list(
    rule = "project_a1_variance_to_nearest_interior",
    outcome = "projected", opt_in = TRUE,
    before = projection$normalized, after = projection$used,
    evidence = list(
      policy = "nearest_interior", reason = "fixture", distance = 1,
      a1_lower_bound = 1, numerical_interior_lower = 1.1,
      support_upper_bound = 75, epsilon = 1e-8,
      requested_buffer = 0.1, preferred_buffer = 0.1,
      representability_floor = 1e-8, effective_buffer = 0.1,
      available_gap = 74, buffer_was_capped = FALSE,
      buffer_was_floored = FALSE, representable_interior = TRUE
    )
  )
  expect_error(.dpprior_validate_target_derivation(projection),
               class = "dpprior_schema_error")
})


test_that("A1 projection retains and recomputes the exact R10 authority", {
  expect_error(
    .dpprior_a1_resolve_target(20L, 5, 2, "nearest", FALSE),
    class = "dpprior_error"
  )
  target <- .schema23_a1_projected_target()
  expect_invisible(.dpprior_validate_target_v1(target))
  evidence <- target$derivation$normalized_to_used$evidence
  expect_identical(
    evidence$numerical_interior_lower,
    evidence$a1_lower_bound + evidence$effective_buffer
  )
  expect_identical(
    evidence$available_gap,
    evidence$support_upper_bound - evidence$a1_lower_bound
  )

  double_buffer <- target
  forged_variance <- evidence$numerical_interior_lower +
    evidence$effective_buffer
  double_buffer$used$var_K <- forged_variance
  double_buffer$implied$variance <- forged_variance
  for (snapshot_name in c("selected_snapshot", "verifier_snapshot")) {
    double_buffer$verification[[snapshot_name]]$achieved$implied$variance <-
      forged_variance
  }
  double_buffer$derivation$normalized_to_used$after <- double_buffer$used
  double_buffer$derivation$normalized_to_used$evidence$
    numerical_interior_lower <- forged_variance
  double_buffer$derivation$normalized_to_used$evidence$
    projected_target$var_K <- forged_variance
  double_buffer$derivation$normalized_to_used$evidence$distance <-
    forged_variance - 2
  double_buffer$provenance$projection$record <- list(
    before = double_buffer$normalized, after = double_buffer$used,
    authority = double_buffer$derivation$normalized_to_used$evidence
  )
  expect_error(.dpprior_validate_target_v1(double_buffer),
               class = "dpprior_schema_error")

  atomic_mutations <- list(
    epsilon = evidence$epsilon * 2,
    buffer = evidence$buffer * 2,
    requested_buffer = evidence$requested_buffer * 2,
    preferred_buffer = evidence$preferred_buffer * 2,
    representability_floor = evidence$representability_floor * 2,
    effective_buffer = evidence$effective_buffer * 2,
    available_gap = evidence$available_gap + 1,
    a1_lower_bound = evidence$a1_lower_bound + 1,
    support_upper_bound = evidence$support_upper_bound + 1,
    distance = evidence$distance + 1,
    buffer_was_capped = !evidence$buffer_was_capped,
    buffer_was_floored = !evidence$buffer_was_floored,
    representable_interior = !evidence$representable_interior,
    reason = "near_a1_lower_boundary",
    policy = "error",
    opt_in = FALSE,
    applied = FALSE
  )
  for (field in names(atomic_mutations)) {
    forged <- target
    forged$derivation$normalized_to_used$evidence[[field]] <-
      atomic_mutations[[field]]
    forged$provenance$projection$record$authority[[field]] <-
      atomic_mutations[[field]]
    if (identical(field, "policy")) {
      forged$provenance$projection$policy <- atomic_mutations[[field]]
    }
    expect_error(.dpprior_validate_target_v1(forged),
                 class = "dpprior_schema_error")
  }

  false_claim <- .schema23_target()
  false_claim$provenance$projection <- list(
    applied = TRUE, opt_in = TRUE, policy = "nearest",
    record = list(fake = TRUE)
  )
  expect_error(.dpprior_validate_target_v1(false_claim),
               class = "dpprior_schema_error")

  partial_opt_in <- .schema23_target()
  partial_opt_in$provenance$projection <- list(
    applied = FALSE, opt_in = TRUE, policy = "nearest", record = NULL
  )
  expect_error(.dpprior_validate_target_v1(partial_opt_in),
               class = "dpprior_schema_error")

  result_claim <- .schema23_fit()
  result_claim$provenance$projection <- list(
    applied = TRUE, opt_in = TRUE, policy = "nearest",
    record = list(fake = TRUE)
  )
  expect_error(.dpprior_validate_result_v1(result_claim),
               class = "dpprior_schema_error")
})


test_that("constructed target requests are exact scientific authority", {
  interval <- .schema23_interval_target("equal_tail")
  family <- .schema23_family_target()
  expect_invisible(.dpprior_validate_target_v1(interval))
  expect_invisible(.dpprior_validate_target_v1(family))

  interval_mutations <- list(
    missing_mean = function(x) {
      x$request$mu_K <- NULL
      x
    },
    garbage = function(x) {
      x$request$garbage <- "ignored"
      x
    },
    changed_mean = function(x) {
      x$request$mu_K <- 7.5
      x
    },
    changed_coverage = function(x) {
      x$request$K_interval$coverage <- 0.7
      x
    }
  )
  for (mutate in interval_mutations) {
    forged <- mutate(interval)
    forged$derivation$request_to_normalized$before <- forged$request
    expect_error(.dpprior_validate_target_v1(forged),
                 class = "dpprior_schema_error")
  }

  family_mutations <- list(
    missing_family = function(x) {
      x$request$family <- NULL
      x
    },
    changed_mean = function(x) {
      x$request$mean <- 9
      x
    },
    changed_variance = function(x) {
      x$request$variance <- 9
      x
    }
  )
  for (mutate in family_mutations) {
    forged <- mutate(family)
    forged$derivation$request_to_normalized$before <- forged$request
    expect_error(.dpprior_validate_target_v1(forged),
                 class = "dpprior_schema_error")
  }

  relabeled <- family
  for (record in c("request", "normalized", "used")) {
    relabeled[[record]]$family$name <- "poisson"
  }
  relabeled$family$name <- "poisson"
  relabeled$derivation$request_to_normalized$before <- relabeled$request
  relabeled$derivation$request_to_normalized$after <- relabeled$normalized
  relabeled$derivation$normalized_to_used$before <- relabeled$normalized
  relabeled$derivation$normalized_to_used$after <- relabeled$used
  expect_error(.dpprior_validate_target_v1(relabeled),
               class = "dpprior_schema_error")

  confidence <- unclass(.schema23_target())
  confidence$request <- list(J = 20L, mean = 5, confidence = "medium")
  confidence$normalized <- c(
    confidence$request, list(interval = NULL, pmf = NULL)
  )
  confidence$normalized$variance <- 999
  confidence$derivation$request_to_normalized <- list(
    rule = "canonicalize_confidence_target", outcome = "canonicalized",
    opt_in = FALSE, before = confidence$request,
    after = confidence$normalized, evidence = list(source = "fixture")
  )
  expect_error(.dpprior_validate_target_derivation(confidence),
               class = "dpprior_schema_error")

  cv <- confidence
  cv$kind <- "cv"
  cv$request <- list(J = 20L, mean = 5, cv = 0.2)
  cv$normalized <- c(cv$request, list(interval = NULL, pmf = NULL))
  cv$normalized$variance <- 999
  cv$derivation$request_to_normalized <- list(
    rule = "canonicalize_cv_target", outcome = "canonicalized",
    opt_in = FALSE, before = cv$request, after = cv$normalized,
    evidence = list(source = "fixture")
  )
  expect_error(.dpprior_validate_target_derivation(cv),
               class = "dpprior_schema_error")
})


test_that("interval and family targets retain independent reconstruction truth", {
  interval <- .schema23_interval_target("equal_tail")
  central <- .schema23_interval_target("central_mass")
  boundary <- .schema23_interval_target("hard_bounds")
  family <- .schema23_family_target()
  for (target in list(interval, central, boundary, family)) {
    expect_invisible(.dpprior_validate_target_v1(target))
    expect_true(target$verification$stability$passed)
  }
  expect_identical(
    names(central$derivation$normalized_to_used$evidence$group_masses),
    c("inside", "outside")
  )
  expect_identical(
    boundary$derivation$normalized_to_used$evidence$boundary,
    list(active = TRUE, side = "minimum")
  )
  expect_null(boundary$derivation$normalized_to_used$evidence$common_tilt)

  forged_family <- family
  verifier_pmf <- rev(forged_family$verification$verifier_snapshot$achieved$pmf)
  forged_family$verification$verifier_snapshot$achieved$pmf <- verifier_pmf
  forged_family$verification$verifier_snapshot$achieved$implied <-
    as.list(.dpprior_target_pmf_moments(verifier_pmf))
  expect_error(.dpprior_validate_target_v1(forged_family),
               class = "dpprior_schema_error")

  forged_implied <- family
  forged_implied$verification$verifier_snapshot$achieved$implied <- list(
    mean = 999, variance = 999
  )
  expect_error(.dpprior_validate_target_v1(forged_implied),
               class = "dpprior_schema_error")

  forged_interval <- interval
  verifier_pmf <- rev(forged_interval$verification$verifier_snapshot$achieved$pmf)
  forged_interval$verification$verifier_snapshot$achieved$pmf <- verifier_pmf
  forged_interval$verification$verifier_snapshot$achieved$implied <-
    as.list(.dpprior_target_pmf_moments(verifier_pmf))
  forged_interval$verification$verifier_snapshot$achieved$interval <-
    as.list(.dpprior_target_interval_masses(verifier_pmf, interval$interval))
  expect_error(.dpprior_validate_target_v1(forged_interval),
               class = "dpprior_schema_error")

  inflated_stability <- interval
  inflated_stability$verification$stability$tolerance[["pmf.l1"]] <- 1
  expect_error(.dpprior_validate_target_v1(inflated_stability),
               class = "dpprior_schema_error")
})


test_that("MaxEnt evidence is recomputed from interval type and fixed controls", {
  target <- .schema23_interval_target("equal_tail")
  evidence <- target$derivation$normalized_to_used$evidence
  mutations <- list(
    group_masses = within(evidence$group_masses, left <- left + 0.01),
    feasible_mean_lower = evidence$feasible_mean_lower - 1,
    feasible_mean_upper = evidence$feasible_mean_upper + 1,
    feasibility_tolerance = 999,
    root_tolerance = 999,
    boundary = list(active = TRUE, side = "minimum"),
    mean_constraint = 7,
    common_tilt = evidence$common_tilt + 0.1,
    constructor_method = "analytic"
  )
  for (field in names(mutations)) {
    forged <- target
    forged$derivation$normalized_to_used$evidence[[field]] <- mutations[[field]]
    expect_error(.dpprior_validate_target_v1(forged),
                 class = "dpprior_schema_error")
  }

  wrong_constraint_authority <- target
  wrong_constraint_authority$tolerances$constraint <- 1
  expect_error(.dpprior_validate_target_v1(wrong_constraint_authority),
               class = "dpprior_schema_error")

  wrong_requested_control <- target
  wrong_requested_control$computation$request$controls$root_tolerance <- 1e-6
  expect_error(.dpprior_validate_target_v1(wrong_requested_control),
               class = "dpprior_schema_error")

  all_mirrors <- target
  all_mirrors$tolerances$constraint <- 1
  for (record in c("request", "used")) {
    all_mirrors$computation[[record]]$controls$constraint_tolerance <- 1
  }
  for (snapshot in c("selected_snapshot", "verifier_snapshot")) {
    all_mirrors$verification[[snapshot]]$tolerances$constraint <- 1
  }
  all_mirrors$verification$settings$constraint_tolerance <- 1
  all_mirrors$verification$stability <-
    .dpprior_expected_target_stability(all_mirrors)
  all_mirrors$verification$components$order_stability <- .schema23_check(
    value = all_mirrors$verification$stability$delta,
    reference = setNames(
      rep(0, length(all_mirrors$verification$stability$delta)),
      names(all_mirrors$verification$stability$delta)
    ), tolerance = all_mirrors$verification$stability$tolerance,
    operator = "lte", source = "independent_target_reconstruction"
  )
  expect_error(.dpprior_validate_target_v1(all_mirrors),
               class = "dpprior_schema_error")

  loose_mean <- .schema23_interval_target("central_mass")
  support <- seq_len(loose_mean$J)
  interval_support <- support[
    support >= loose_mean$interval$lower &
      support <= loose_mean$interval$upper
  ]
  groups <- list(
    inside = interval_support,
    outside = setdiff(support, interval_support)
  )
  masses <- c(inside = 0.8, outside = 0.2)
  loose_tilt <- stats::uniroot(
    function(value) {
      candidate <- .schema23_maxent_pmf(
        loose_mean$J, groups, masses, tilt = value
      )
      sum(support * candidate) - 3.53
    }, interval = c(-10, 10), tol = 1e-12
  )$root
  loose_pmf <- .schema23_maxent_pmf(
    loose_mean$J, groups, masses, tilt = loose_tilt
  )
  loose_mean$pmf <- loose_pmf
  loose_mean$used$pmf <- loose_pmf
  loose_mean$derivation$normalized_to_used$after <- loose_mean$used
  loose_mean$derivation$normalized_to_used$evidence$common_tilt <- loose_tilt
  loose_mean$implied <- as.list(.dpprior_target_pmf_moments(loose_pmf))
  loose_mean$achieved_interval <- as.list(
    .dpprior_target_interval_masses(loose_pmf, loose_mean$interval)
  )
  loose_mean$residuals$mean <- loose_mean$implied$mean - 6.5
  loose_mean$tolerances$constraint <- 1
  for (record in c("request", "used")) {
    loose_mean$computation[[record]]$controls$constraint_tolerance <- 1
  }
  for (snapshot in c("selected_snapshot", "verifier_snapshot")) {
    loose_mean$verification[[snapshot]]$achieved <- list(
      implied = loose_mean$implied,
      interval = loose_mean$achieved_interval,
      pmf = loose_pmf
    )
    loose_mean$verification[[snapshot]]$residuals <- loose_mean$residuals
    loose_mean$verification[[snapshot]]$tolerances <- loose_mean$tolerances
  }
  loose_mean$verification$settings$constraint_tolerance <- 1
  loose_mean$verification$stability <-
    .dpprior_expected_target_stability(loose_mean)
  loose_mean$verification$components$target_reconstruction <- .schema23_check(
    value = 0, reference = 0, tolerance = .TOL_PMF_SUM, operator = "lte",
    source = "independent_target_reconstruction"
  )
  loose_mean$verification$components$order_stability <- .schema23_check(
    value = loose_mean$verification$stability$delta,
    reference = setNames(
      rep(0, length(loose_mean$verification$stability$delta)),
      names(loose_mean$verification$stability$delta)
    ), tolerance = loose_mean$verification$stability$tolerance,
    operator = "lte", source = "independent_target_reconstruction"
  )
  expect_equal(loose_mean$implied$mean, 3.53, tolerance = 1e-8)
  expect_error(.dpprior_validate_target_v1(loose_mean),
               class = "dpprior_schema_error")

  for (inflated_root in c(1, 1e6)) {
    root_inflation <- target
    for (record in c("request", "used")) {
      root_inflation$computation[[record]]$controls$root_tolerance <-
        inflated_root
    }
    root_inflation$derivation$normalized_to_used$evidence$root_tolerance <-
      inflated_root
    root_inflation$verification$settings$root_tolerance <- inflated_root
    expect_error(.dpprior_validate_target_v1(root_inflation),
                 class = "dpprior_schema_error")
  }

  iteration_inflation <- target
  for (record in c("request", "used")) {
    iteration_inflation$computation[[record]]$controls$max_iterations <-
      1000000L
  }
  iteration_inflation$verification$settings$max_iterations <- 1000000L
  expect_error(.dpprior_validate_target_v1(iteration_inflation),
               class = "dpprior_schema_error")

  for (field in c("method", "settings", "component", "invariant", "source")) {
    forged <- target
    if (identical(field, "method")) {
      forged$verification$method <- "rubber_stamp"
    } else if (identical(field, "settings")) {
      forged$verification$settings$policy <- "rubber_stamp"
    } else if (identical(field, "component")) {
      names(forged$verification$components)[[1L]] <- "rubber_stamp"
    } else if (identical(field, "invariant")) {
      names(forged$verification$invariants)[[1L]] <- "rubber_stamp"
    } else {
      forged$verification$components$target_reconstruction$source <-
        "rubber_stamp"
    }
    expect_error(.dpprior_validate_target_v1(forged),
                 class = "dpprior_schema_error")
  }

  relabeled <- target
  relabeled$interval$type <- "central_mass"
  relabeled$normalized$interval <- relabeled$interval
  relabeled$used$interval <- relabeled$interval
  relabeled$derivation$request_to_normalized$after <- relabeled$normalized
  relabeled$derivation$normalized_to_used$before <- relabeled$normalized
  relabeled$derivation$normalized_to_used$after <- relabeled$used
  relabeled$request$K_interval <- relabeled$interval
  relabeled$derivation$request_to_normalized$before <- relabeled$request
  expect_error(.dpprior_validate_target_v1(relabeled),
               class = "dpprior_schema_error")
})


test_that("target infeasibility requires a request-bound analytic certificate", {
  target <- .schema23_infeasible_interval_target()
  expect_invisible(.dpprior_validate_target_v1(target))
  expect_null(target$pmf)
  expect_null(target$verification$selected_snapshot)

  mutations <- list(
    generic_check = function(x) {
      x$verification$components$infeasibility_certificate <- .schema23_check(
        value = 1, reference = 0, tolerance = 0, operator = "gt",
        source = "analytic_maxent_interval_feasibility"
      )
      x
    },
    wrong_bound = function(x) {
      x$verification$settings$certificate$feasible_mean_upper <- 999
      x
    },
    wrong_source = function(x) {
      x$verification$settings$certificate$source <- "fabricated"
      x
    },
    wrong_termination = function(x) {
      x$computation$termination$source <- "constructor"
      x
    },
    wrong_attempt = function(x) {
      x$computation$attempts[[1L]]$method <- "rubber_stamp"
      x
    },
    invented_pmf = function(x) {
      x$pmf <- rep(0.2, 5L)
      x
    }
  )
  for (mutate in mutations) {
    expect_error(.dpprior_validate_target_v1(mutate(target)),
                 class = "dpprior_schema_error")
  }
})


test_that("infeasibility certificates are substantive and target bound", {
  infeasible <- .schema23_no_candidate("dual_hard", "infeasible")
  expect_length(infeasible$computation$attempts, 1L)
  expect_identical(
    infeasible$computation$attempts[[1L]]$method,
    "analytic_monotonicity_feasibility_probe"
  )
  expect_identical(infeasible$computation$attempts[[1L]]$exit_code, 0L)
  expect_length(infeasible$computation$candidate_evaluations, 0L)
  expect_identical(
    infeasible$constraint$feasibility[c(
      "candidate_count", "verified_candidate_count", "feasible_candidate_count"
    )],
    list(
      candidate_count = 0L, verified_candidate_count = 0L,
      feasible_candidate_count = 0L
    )
  )
  empty <- infeasible
  empty$constraint$feasibility$certificate <- list()
  expect_error(.dpprior_validate_result_v1(empty),
               class = "dpprior_schema_error")

  fabricated <- infeasible
  fabricated$constraint$feasibility$certificate$method <- "fabricated"
  expect_error(.dpprior_validate_result_v1(fabricated),
               class = "dpprior_schema_error")

  wrong_target <- infeasible
  wrong_target$constraint$feasibility$certificate$target_value <- 0.2
  expect_error(.dpprior_validate_result_v1(wrong_target),
               class = "dpprior_schema_error")

  wrong_bound <- infeasible
  wrong_bound$constraint$feasibility$certificate$lower_bound <- 0.7
  expect_error(.dpprior_validate_result_v1(wrong_bound),
               class = "dpprior_schema_error")

  wrong_corner <- infeasible
  wrong_corner$constraint$feasibility$certificate$minimum$refined$value <-
    0.7
  expect_error(.dpprior_validate_result_v1(wrong_corner),
               class = "dpprior_schema_error")

  wrong_source <- infeasible
  wrong_source$constraint$feasibility$certificate$source <-
    "unverified_claim"
  expect_error(.dpprior_validate_result_v1(wrong_source),
               class = "dpprior_schema_error")

  false_corner_claim <- infeasible
  false_certificate <- false_corner_claim$constraint$feasibility$certificate
  false_certificate$domain <- list(
    log_a = c(-2, 2), log_b = c(-2, 2),
    a = exp(c(-2, 2)), b = exp(c(-2, 2))
  )
  reset_corner <- function(corner, eta, value) {
    corner$eta <- c(log_a = eta[[1L]], log_b = eta[[2L]])
    corner$parameters <- list(a = exp(eta[[1L]]), b = exp(eta[[2L]]))
    corner$selected$value <- value
    corner$refined$value <- value
    corner$uncertainty <- 64 * .Machine$double.eps
    corner$lower <- value - corner$uncertainty
    corner$upper <- value + corner$uncertainty
    corner
  }
  false_certificate$minimum <- reset_corner(
    false_certificate$minimum, c(2, -2), 0.95
  )
  false_certificate$maximum <- reset_corner(
    false_certificate$maximum, c(-2, 2), 0.99
  )
  false_certificate$lower_bound <- false_certificate$minimum$lower
  false_certificate$upper_bound <- false_certificate$maximum$upper
  false_corner_claim$constraint$feasibility$certificate <- false_certificate
  false_corner_claim$verification$components$infeasibility_certificate <-
    .schema23_check(
      value = false_certificate$lower_bound,
      reference = false_certificate$target_value,
      tolerance = false_certificate$tolerance, operator = "gt"
    )
  false_corner_claim$computation$attempts[[1L]]$bounds <- list(
    log_a = c(-2, 2), log_b = c(-2, 2)
  )
  false_corner_claim$computation$attempts[[1L]]$candidate_objective <- 0.95
  expect_error(.dpprior_validate_result_v1(false_corner_claim),
               class = "dpprior_schema_error")
})


test_that("soft optimality applies componentwise one-sided boundary KKT", {
  lower <- .schema23_soft_optimality()
  lower$gradient[["log_shape"]] <- 100
  lower$gradient_method[["log_shape"]] <-
    "forward feasible-direction KKT gradient"
  lower$bound_state[["log_shape"]] <- "lower"
  lower$stationarity_operator[["log_shape"]] <- "gte"
  expect_invisible(.dpprior_validate_soft_optimality(lower, available = TRUE))

  upper <- .schema23_soft_optimality()
  upper$gradient[["log_rate"]] <- -100
  upper$gradient_method[["log_rate"]] <-
    "backward feasible-direction KKT gradient"
  upper$bound_state[["log_rate"]] <- "upper"
  upper$stationarity_operator[["log_rate"]] <- "lte"
  expect_invisible(.dpprior_validate_soft_optimality(upper, available = TRUE))

  forged_sign <- lower
  forged_sign$gradient[["log_shape"]] <- -100
  expect_error(
    .dpprior_validate_soft_optimality(forged_sign, available = TRUE),
    class = "dpprior_schema_error"
  )

  forged_operator <- lower
  forged_operator$stationarity_operator[["log_shape"]] <- "abs_lte"
  expect_error(
    .dpprior_validate_soft_optimality(forged_operator, available = TRUE),
    class = "dpprior_schema_error"
  )

  result_state <- .schema23_fit("dual_soft")
  result_state$tradeoff$optimality$bound_state[["log_shape"]] <- "lower"
  result_state$tradeoff$optimality$gradient_method[["log_shape"]] <-
    "forward feasible-direction KKT gradient"
  result_state$tradeoff$optimality$stationarity_operator[["log_shape"]] <-
    "gte"
  expect_error(.dpprior_validate_result_v1(result_state),
               class = "dpprior_schema_error")
})


test_that("candidate selection and termination reject forged optimizer evidence", {
  for (mode in c("a2_kl", "dual_hard", "dual_soft")) {
    fit <- .schema23_fit(mode)
    better <- fit$computation$attempts[[1L]]
    better$id <- "attempt-better"
    better$selected <- FALSE
    better$reason_code <- "eligible_not_selected"
    better$candidate_objective <- -1
    fit$computation$attempts <- c(list(better), fit$computation$attempts)
    expect_error(.dpprior_validate_result_v1(fit),
                 class = "dpprior_schema_error")
  }

  method <- .schema23_fit("a2_kl")
  method$computation$attempts[[1L]]$method <- "different_algorithm"
  expect_error(.dpprior_validate_result_v1(method),
               class = "dpprior_schema_error")

  unofficial_hard_alias <- .schema23_fit("dual_hard")
  unofficial_hard_alias$computation$attempts[[1L]]$method <-
    "K_only_L_BFGS_B"
  unofficial_hard_alias$computation$candidate_evaluations[[1L]]$method <-
    "K_only_L_BFGS_B"
  expect_error(.dpprior_validate_result_v1(unofficial_hard_alias),
               class = "dpprior_schema_error")

  moment_method <- .schema23_fit("a2_moment")
  moment_method$computation$attempts[[1L]]$method <- "different_algorithm"
  expect_error(.dpprior_validate_result_v1(moment_method),
               class = "dpprior_schema_error")

  reason <- .schema23_fit("a2_kl")
  reason$computation$attempts[[1L]]$reason_code <- "failed_rejected"
  expect_error(.dpprior_validate_result_v1(reason),
               class = "dpprior_schema_error")

  moment_reason <- .schema23_fit("a2_moment")
  moment_reason$computation$attempts[[1L]]$reason_code <- "failed_rejected"
  expect_error(.dpprior_validate_result_v1(moment_reason),
               class = "dpprior_schema_error")

  termination <- .schema23_fit("a2_kl")
  termination$computation$termination$code <- "closed_form"
  termination$computation$termination$source <- "fatal_error"
  expect_error(.dpprior_validate_result_v1(termination),
               class = "dpprior_schema_error")

  optimizer_closed_form <- .schema23_fit("a2_kl")
  optimizer_closed_form$computation$termination$code <- "closed_form"
  expect_error(.dpprior_validate_result_v1(optimizer_closed_form),
               class = "dpprior_schema_error")

  empty_optimizer <- .schema23_no_candidate("dual_hard", "failed")
  empty_optimizer$computation$termination$code <- "failed"
  empty_optimizer$computation$termination$source <- "optimizer"
  expect_error(.dpprior_validate_result_v1(empty_optimizer),
               class = "dpprior_schema_error")

  rejected_optimizer <- .schema23_hard_failed_rejected_candidates()
  rejected_optimizer$computation$termination$code <- "failed"
  rejected_optimizer$computation$termination$source <- "optimizer"
  expect_error(.dpprior_validate_result_v1(rejected_optimizer),
               class = "dpprior_schema_error")

  certificate_code <- .schema23_no_candidate("dual_hard", "infeasible")
  certificate_code$computation$termination$code <- "infeasible"
  expect_error(.dpprior_validate_result_v1(certificate_code),
               class = "dpprior_schema_error")
})


test_that("parameterless A2 migration quarantine has one closed route", {
  codes <- c(
    "selected", "approximate", "deterministic", "closed_form",
    "legacy_migration_no_selection"
  )
  sources <- c(
    "constructor", "optimizer", "legacy_schema_upgrade", "fabricated_source"
  )
  for (mode in c("a2_moment", "a2_kl")) {
    migrated <- .schema23_migrated_A2(mode)
    expect_invisible(.dpprior_validate_result_v1(migrated))
    expect_identical(
      migrated$computation$termination[c("code", "source")],
      list(code = "legacy_migration_no_selection", source = "constructor")
    )

    for (code in codes) {
      for (source in sources) {
        if (identical(code, "legacy_migration_no_selection") &&
            identical(source, "constructor")) {
          next
        }
        forged <- migrated
        forged$computation$termination$code <- code
        forged$computation$termination$source <- source
        expect_error(.dpprior_validate_result_v1(forged),
                     class = "dpprior_schema_error")
      }
    }

    for (mutation in c(
      "parameters", "achieved", "residuals", "orders", "scaling",
      "resources", "verification_reason", "verification_settings",
      "approximation", "backend", "migration", "compatibility"
    )) {
      forged <- migrated
      if (identical(mutation, "parameters")) {
        forged$parameters <- .dpprior_new_parameters(
          1, 1, "Gamma(shape=a, rate=b)"
        )
      } else if (identical(mutation, "achieved")) {
        forged$achieved$migration_note <- "fabricated public evidence"
      } else if (identical(mutation, "residuals")) {
        forged$residuals$migration_note <- "fabricated public evidence"
      } else if (identical(mutation, "orders")) {
        forged$computation$orders$M_requested <- 80L
        forged$computation$orders$requested_reason <- "forged order"
      } else if (identical(mutation, "scaling")) {
        forged$computation$scaling$formula <- "forged scaling"
      } else if (identical(mutation, "resources")) {
        forged$computation$resources <- list(worker = "forged")
      } else if (identical(mutation, "verification_reason")) {
        forged$verification$reason <- "fabricated quarantine reason"
      } else if (identical(mutation, "verification_settings")) {
        forged$verification$settings <- list(rubber_stamp = TRUE)
      } else if (identical(mutation, "approximation")) {
        forged$provenance$approximation$kind <- "fabricated migration"
      } else if (identical(mutation, "backend")) {
        forged$provenance$backend$implementation <- "fabricated_adapter"
      } else if (identical(mutation, "migration")) {
        forged$provenance$migration$missing_evidence <-
          forged$provenance$migration$missing_evidence[-1L]
      } else {
        forged$compatibility$views$source$required_action <- "none"
      }
      expect_error(.dpprior_validate_result_v1(forged),
                   class = "dpprior_schema_error")
    }

    failed_rewrite <- migrated
    failed_rewrite$status <- "failed"
    failed_rewrite$message <- "fabricated no-candidate failure"
    failed_rewrite$provenance$approximation <- list(
      active = FALSE, opt_in = FALSE, kind = NULL, warning_code = NULL
    )
    failed_rewrite$verification <- .dpprior_new_verification(
      method = "no_candidate", performed = FALSE, passed = FALSE,
      reason = "fabricated no candidate", settings = list(),
      selected_snapshot = NULL, verifier_snapshot = NULL, stability = NULL,
      components = list(), invariants = list(
        no_public_candidate = .schema23_check(
          value = TRUE, reference = TRUE, tolerance = NULL,
          operator = "identical", source = "independent_verifier"
        )
      )
    )
    failed_rewrite$computation$termination$code <- "no_candidate"
    failed_rewrite$computation$termination$source <- "no_candidate"
    expect_error(.dpprior_validate_result_v1(failed_rewrite),
                 class = "dpprior_schema_error")

    erased_markers <- failed_rewrite
    erased_markers$provenance$backend$implementation <- "native_current_fit"
    erased_markers$provenance$migration <- list(
      source_schema = "dpprior.result/1", adapter = "direct_constructor",
      lossless = FALSE, missing_evidence = character(), warnings = character()
    )
    erased_markers$compatibility$views$source$source_schema <-
      "dpprior.result/1"
    erased_markers$compatibility$views$source$
      public_candidate_quarantined <- FALSE
    erased_markers$compatibility$views$source$required_action <- "none"
    erased_markers$target$K$provenance$approximation <- list(
      active = FALSE, opt_in = FALSE, kind = NULL, warning_code = NULL
    )
    erased_markers$target$K$provenance$backend$implementation <-
      "native_current_target"
    erased_markers$target$K$provenance$migration <- list(
      source_schema = "dpprior.target/1", adapter = "direct_constructor",
      lossless = FALSE, missing_evidence = character(), warnings = character()
    )
    erased_markers$target$K$verification$method <- "unverified_target"
    erased_markers$target$K$computation$termination$code <- "approximate"
    expect_error(.dpprior_validate_result_v1(erased_markers),
                 class = "dpprior_schema_error")

    retained_failure <- erased_markers
    retained_failure$message <- "retained optimizer error; no public candidate"
    retained_failure$computation$termination$code <- "failed"
    retained_failure$computation$termination$source <- "optimizer"
    retained_failure$computation$attempts <- list(.dpprior_new_attempt(
      id = "attempt-1", stage = "primary", method = if (
        identical(mode, "a2_moment")
      ) "A2-MN" else "A2-KL",
      start = NULL, bounds = NULL, control = NULL, exit_code = 2L,
      message = "optimizer failed", iterations = 1L,
      evaluations = list(function_count = 1L),
      candidate_parameters = NULL, candidate_objective = NULL,
      elapsed_seconds = NULL, warnings = character(),
      error = list(
        class = "simpleError", code = "optimizer_error",
        message = "optimizer failed"
      ),
      selected = FALSE, reason_code = "optimizer_error",
      unavailable = c(
        start = "not retained", bounds = "not retained",
        control = "not retained", candidate_parameters = "optimizer failed",
        candidate_objective = "optimizer failed",
        elapsed_seconds = "not retained"
      )
    ))
    expect_invisible(.dpprior_validate_result_v1(retained_failure))
  }
})


test_that("A2 migration target and optional fixed audit are fully bound", {
  quarantine_boundary <- list(
    authority = "non_authoritative", lossy = TRUE,
    consumer_policy = "ignored_by_scientific_and_decision_consumers"
  )
  positive_cases <- list(
    direct_moment = list(mode = "a2_moment", wrapper = FALSE,
                         include_diagnostics = FALSE),
    direct_kl = list(mode = "a2_kl", wrapper = FALSE,
                     include_diagnostics = FALSE),
    wrapper_moment = list(mode = "a2_moment", wrapper = TRUE,
                          include_diagnostics = FALSE),
    wrapper_diagnostics = list(mode = "a2_moment", wrapper = TRUE,
                               include_diagnostics = TRUE)
  )
  for (case in positive_cases) {
    without_audit <- .schema23_migrated_A2(
      case$mode, verify = FALSE, wrapper = case$wrapper,
      include_diagnostics = case$include_diagnostics
    )
    with_audit <- .schema23_migrated_A2(
      case$mode, verify = TRUE, wrapper = case$wrapper,
      include_diagnostics = case$include_diagnostics
    )
    expect_invisible(.dpprior_validate_result_v1(without_audit))
    expect_invisible(.dpprior_validate_result_v1(with_audit))
    expect_identical(names(without_audit$compatibility$views), "source")
    expect_identical(
      names(with_audit$compatibility$views),
      c("source", "fixed_candidate_recomputation")
    )
    expect_identical(
      without_audit$compatibility$views$source[names(quarantine_boundary)],
      quarantine_boundary
    )
    expect_identical(
      without_audit$compatibility$deprecations$legacy_schema[
        names(quarantine_boundary)
      ],
      quarantine_boundary
    )
    expect_identical(
      without_audit$target$K$compatibility$deprecations$
        legacy_target_view[names(quarantine_boundary)],
      quarantine_boundary
    )
  }

  for (mode in c("a2_moment", "a2_kl")) {
    custom_order <- .schema23_migrated_A2(
      mode, verify = TRUE, M_verify = 200L
    )
    expect_invisible(.dpprior_validate_result_v1(custom_order))
    expect_identical(
      custom_order$compatibility$views$fixed_candidate_recomputation$
        M_verification,
      200L
    )
  }

  codes <- c(
    "selected", "approximate", "deterministic", "closed_form",
    "legacy_migration_no_selection"
  )
  sources <- c(
    "constructor", "optimizer", "legacy_schema_upgrade", "fabricated_source"
  )
  for (mode in c("a2_moment", "a2_kl")) {
    migrated <- .schema23_migrated_A2(mode, verify = FALSE)
    for (code in codes) {
      for (source in sources) {
        if (identical(code, "legacy_migration_no_selection") &&
            identical(source, "constructor")) {
          next
        }
        forged <- migrated
        forged$target$K$computation$termination$code <- code
        forged$target$K$computation$termination$source <- source
        expect_error(.dpprior_validate_result_v1(forged),
                     class = "dpprior_schema_error")
      }
    }

    target_mutations <- c(
      "message", "computation_method", "computation_controls", "orders",
      "scaling", "resources", "attempts", "verification_reason",
      "verification_settings", "verification_components",
      "approximation_opt_in", "approximation_kind", "warning_code",
      "parameterization", "backend_package", "backend_version",
      "backend_implementation", "migration_missing", "migration_warning",
      "legacy", "compatibility_extra", "legacy_target_fields"
    )
    for (mutation in target_mutations) {
      forged <- migrated
      K <- forged$target$K
      if (identical(mutation, "message")) {
        K$message <- "fabricated target quarantine"
      } else if (identical(mutation, "computation_method")) {
        K$computation$request$method <- "legacy_target:fabricated"
      } else if (identical(mutation, "computation_controls")) {
        K$computation$request$controls <- list(rubber_stamp = TRUE)
        K$computation$used$controls <- list(rubber_stamp = TRUE)
      } else if (identical(mutation, "orders")) {
        K$computation$orders$M_selected <- 80L
        K$computation$orders$selected_reason <- "fabricated order"
      } else if (identical(mutation, "scaling")) {
        K$computation$scaling$formula <- "fabricated scaling"
      } else if (identical(mutation, "resources")) {
        K$computation$resources <- list(worker = "fabricated")
      } else if (identical(mutation, "attempts")) {
        K$computation$attempts <- list(fabricated = list(value = TRUE))
      } else if (identical(mutation, "verification_reason")) {
        K$verification$reason <- "fabricated target verification reason"
      } else if (identical(mutation, "verification_settings")) {
        K$verification$settings <- list(rubber_stamp = TRUE)
      } else if (identical(mutation, "verification_components")) {
        K$verification$components <- list(rubber_stamp = TRUE)
      } else if (identical(mutation, "approximation_opt_in")) {
        K$provenance$approximation$opt_in <- TRUE
      } else if (identical(mutation, "approximation_kind")) {
        K$provenance$approximation$kind <- "fabricated_migration"
      } else if (identical(mutation, "warning_code")) {
        K$provenance$approximation$warning_code <- "fabricated_warning"
      } else if (identical(mutation, "parameterization")) {
        K$provenance$parameterization <- "fabricated_parameterization"
      } else if (identical(mutation, "backend_package")) {
        K$provenance$backend$package <- "fabricated_package"
      } else if (identical(mutation, "backend_version")) {
        K$provenance$backend$package_version <- "999.0.0"
      } else if (identical(mutation, "backend_implementation")) {
        K$provenance$backend$implementation <- "fabricated_adapter"
      } else if (identical(mutation, "migration_missing")) {
        K$provenance$migration$missing_evidence <-
          K$provenance$migration$missing_evidence[-1L]
      } else if (identical(mutation, "migration_warning")) {
        K$provenance$migration$warnings <- "fabricated_warning"
      } else if (identical(mutation, "legacy")) {
        K$provenance$legacy$active <- TRUE
      } else if (identical(mutation, "compatibility_extra")) {
        K$compatibility$views$fabricated <- list(value = TRUE)
      } else {
        K$compatibility$views$legacy_target_fields <- "fabricated"
      }
      forged$target$K <- K
      expect_error(.dpprior_validate_result_v1(forged),
                   class = "dpprior_schema_error")
    }

    top_mutations <- c(
      "message", "tolerances", "controls", "parameterization",
      "backend_version", "source_authority", "source_lossy",
      "source_consumer_policy", "deprecation_authority",
      "target_deprecation_authority", "extra_view", "extra_alias"
    )
    for (mutation in top_mutations) {
      forged <- migrated
      if (identical(mutation, "message")) {
        forged$message <- "fabricated migration message"
      } else if (identical(mutation, "tolerances")) {
        forged$tolerances$boundary <- 1
      } else if (identical(mutation, "controls")) {
        forged$computation$request$controls$selection_tolerance <- 1
        forged$computation$used$controls$selection_tolerance <- 1
      } else if (identical(mutation, "parameterization")) {
        forged$provenance$parameterization <- "fabricated_parameterization"
      } else if (identical(mutation, "backend_version")) {
        forged$provenance$backend$package_version <- "999.0.0"
      } else if (identical(mutation, "source_authority")) {
        forged$compatibility$views$source$authority <- "authoritative"
      } else if (identical(mutation, "source_lossy")) {
        forged$compatibility$views$source$lossy <- FALSE
      } else if (identical(mutation, "source_consumer_policy")) {
        forged$compatibility$views$source$consumer_policy <-
          "may_drive_scientific_decisions"
      } else if (identical(mutation, "deprecation_authority")) {
        forged$compatibility$deprecations$legacy_schema$authority <-
          "authoritative"
      } else if (identical(mutation, "target_deprecation_authority")) {
        forged$target$K$compatibility$deprecations$
          legacy_target_view$authority <- "authoritative"
      } else if (identical(mutation, "extra_view")) {
        forged$compatibility$views$fabricated <- list(value = TRUE)
      } else {
        forged$compatibility$top_level_aliases <- c(fake = "fabricated")
        forged$fake <- "fabricated"
      }
      expect_error(.dpprior_validate_result_v1(forged),
                   class = "dpprior_schema_error")
    }
  }

  for (mode in c("a2_moment", "a2_kl")) {
    migrated <- .schema23_migrated_A2(mode, verify = TRUE)
    audit <- migrated$compatibility$views$fixed_candidate_recomputation
    for (field in names(audit)) {
      forged <- migrated
      forged$compatibility$views$fixed_candidate_recomputation[[field]] <-
        .schema23_mutate_plain_leaf(audit[[field, exact = TRUE]])
      expect_error(.dpprior_validate_result_v1(forged),
                   class = "dpprior_schema_error")

      removed <- migrated
      removed$compatibility$views$fixed_candidate_recomputation[[field]] <- NULL
      expect_error(.dpprior_validate_result_v1(removed),
                   class = "dpprior_schema_error")
    }
    extra_field <- migrated
    extra_field$compatibility$views$fixed_candidate_recomputation$
      fabricated <- TRUE
    expect_error(.dpprior_validate_result_v1(extra_field),
                 class = "dpprior_schema_error")

    atomic_audit <- migrated
    atomic_audit$compatibility$views$fixed_candidate_recomputation <- 1
    expect_error(.dpprior_validate_result_v1(atomic_audit),
                 class = "dpprior_schema_error")

    wrong_source_candidate <- migrated
    wrong_source_candidate$compatibility$views$source$source_parameters$a <-
      wrong_source_candidate$compatibility$views$source$source_parameters$a +
      0.25
    wrong_source_candidate$compatibility$views$
      fixed_candidate_recomputation$candidate$a <-
      wrong_source_candidate$compatibility$views$source$source_parameters$a
    expect_error(.dpprior_validate_result_v1(wrong_source_candidate),
                 class = "dpprior_schema_error")

    # No independent verify-request bit exists in R24. Removing the complete
    # optional audit is therefore indistinguishable from, and safely demotes to,
    # the canonical verify=FALSE quarantine. It cannot raise any status flag.
    conservative_downgrade <- migrated
    conservative_downgrade$compatibility$views$
      fixed_candidate_recomputation <- NULL
    expect_invisible(.dpprior_validate_result_v1(conservative_downgrade))
    expect_false(conservative_downgrade$usable)
    expect_false(conservative_downgrade$verified)
  }
})


test_that("opaque migration metadata cannot influence canonical consumers", {
  base <- .schema23_migrated_A2("a2_kl", verify = FALSE)
  base_state <- .schema23_migration_decision_and_consumer_state(base)
  mutations <- list(
    legacy_mu = function(x) {
      x$target$K$compatibility$views$legacy_target$mu_K <-
        x$target$K$compatibility$views$legacy_target$mu_K + 0.25
      x
    },
    legacy_variance = function(x) {
      x$target$K$compatibility$views$legacy_target$var_K <-
        x$target$K$compatibility$views$legacy_target$var_K + 0.25
      x
    },
    legacy_df = function(x) {
      x$target$K$compatibility$views$legacy_target$df <-
        x$target$K$compatibility$views$legacy_target$df + 0.25
      x
    },
    legacy_scale = function(x) {
      x$target$K$compatibility$views$legacy_target$scale <-
        x$target$K$compatibility$views$legacy_target$scale + 0.25
      x
    },
    source_status = function(x) {
      x$compatibility$views$source$source_status <- "historical_unknown"
      x
    },
    source_converged = function(x) {
      x$compatibility$views$source$source_converged <- FALSE
      x
    },
    source_digest = function(x) {
      x$compatibility$views$source$source_digest <- paste(rep("a", 32L),
                                                          collapse = "")
      x
    },
    quarantined_a = function(x) {
      x$compatibility$views$source$source_parameters$a <-
        x$compatibility$views$source$source_parameters$a + 0.25
      x
    },
    quarantined_b = function(x) {
      x$compatibility$views$source$source_parameters$b <-
        x$compatibility$views$source$source_parameters$b + 0.25
      x
    },
    discarded_metadata = function(x) {
      x$compatibility$views$source[
        "discarded_legacy_fields"
      ] <- list(NULL)
      x
    }
  )
  for (mutate in mutations) {
    opaque <- mutate(base)
    expect_invisible(.dpprior_validate_result_v1(opaque))
    expect_identical(
      .schema23_migration_decision_and_consumer_state(opaque),
      base_state
    )
  }

  wrapper <- .schema23_migrated_A2(
    "a2_moment", verify = FALSE, wrapper = TRUE
  )
  wrapper_state <- .schema23_migration_decision_and_consumer_state(wrapper)
  wrapper$target$K$compatibility$views$legacy_target$confidence <- "low"
  expect_invisible(.dpprior_validate_result_v1(wrapper))
  expect_identical(
    .schema23_migration_decision_and_consumer_state(wrapper),
    wrapper_state
  )

  audited <- .schema23_migrated_A2("a2_kl", verify = TRUE)
  audited_state <- .schema23_migration_decision_and_consumer_state(audited)
  audited$target$K$compatibility$views$legacy_target$df <-
    audited$target$K$compatibility$views$legacy_target$df + 0.25
  audited$compatibility$views$source$source_status <- "historical_unknown"
  audited$compatibility$views$source$source_digest <- paste(
    rep("b", 32L), collapse = ""
  )
  expect_invisible(.dpprior_validate_result_v1(audited))
  expect_identical(
    .schema23_migration_decision_and_consumer_state(audited),
    audited_state
  )
})


test_that("native producer compatibility views stay non-authoritative", {
  fits <- list(
    a1 = DPprior_a1(50L, 5, 8),
    a2_moment = DPprior_a2_newton(50L, 5, 8),
    a2_kl_explicit = DPprior_a2_kl(
      12L, rep(1 / 12, 12L), method = "pmf", M = 40L
    ),
    a2_kl_chisq = DPprior_a2_kl(
      20L, list(mu_K = 4, var_K = 5), method = "chisq", M = 40L
    )
  )
  markers <- c(
    a1 = "a1_v0", a2_moment = "legacy_v2",
    a2_kl_explicit = "a2_kl_v0", a2_kl_chisq = "a2_kl_v0"
  )
  expected_boundary <- list(
    authority = "non_authoritative", lossy = TRUE,
    consumer_policy = "ignored_by_scientific_and_decision_consumers"
  )
  expected_codes <- c(
    a1 = "a1_flat_v0_view_quarantined",
    a2_moment = "a2_moment_legacy_view",
    a2_kl_explicit = "a2_kl_flat_v0_view_quarantined",
    a2_kl_chisq = "a2_kl_flat_v0_view_quarantined"
  )

  for (name in names(fits)) {
    fit <- fits[[name]]
    marker <- unname(markers[[name]])
    expect_invisible(.dpprior_validate_result_v1(fit))
    expect_identical(
      fit$compatibility$views[[marker]][names(expected_boundary)],
      expected_boundary
    )
    expect_identical(
      fit$compatibility$deprecations[[marker]][names(expected_boundary)],
      expected_boundary
    )
    expect_identical(
      fit$compatibility$deprecations[[marker]]$code,
      unname(expected_codes[[name]])
    )

    baseline <- .schema23_migration_decision_and_consumer_state(fit)
    opaque <- fit
    if (identical(marker, "a1_v0")) {
      opaque$compatibility$views[[marker]]$a <- 999
      opaque$compatibility$views[[marker]]$achieved$mu_K <- 999
    } else if (identical(marker, "legacy_v2")) {
      opaque$compatibility$views[[marker]]$numerical_candidate$a <- 999
      opaque$compatibility$views[[marker]]$source_status <- "failed"
    } else {
      opaque$compatibility$views[[marker]]$a <- 999
      opaque$compatibility$views[[marker]]$achieved$mu_K <- 999
    }
    expect_invisible(.dpprior_validate_result_v1(opaque))
    expect_identical(
      .schema23_migration_decision_and_consumer_state(opaque), baseline
    )

    for (field in names(expected_boundary)) {
      forged_view <- fit
      forged_view$compatibility$views[[marker]][[field]] <-
        if (identical(field, "lossy")) FALSE else "forged"
      expect_error(.dpprior_validate_result_v1(forged_view),
                   class = "dpprior_schema_error")

      forged_deprecation <- fit
      forged_deprecation$compatibility$deprecations[[marker]][[field]] <-
        if (identical(field, "lossy")) FALSE else "forged"
      expect_error(.dpprior_validate_result_v1(forged_deprecation),
                   class = "dpprior_schema_error")
    }

    wrong_code <- fit
    wrong_code$compatibility$deprecations[[marker]]$code <- "forged"
    expect_error(.dpprior_validate_result_v1(wrong_code),
                 class = "dpprior_schema_error")

    extra_view <- fit
    extra_view$compatibility$views$decision_ready <- TRUE
    expect_error(.dpprior_validate_result_v1(extra_view),
                 class = "dpprior_schema_error")
    extra_science <- fit
    extra_science$compatibility$views$scientific_result <- list(
      status = "converged", usable = TRUE, verified = TRUE
    )
    expect_error(.dpprior_validate_result_v1(extra_science),
                 class = "dpprior_schema_error")
    extra_alias <- fit
    extra_alias$compatibility$top_level_aliases <- c(
      extra_alias$compatibility$top_level_aliases,
      decision_ready = "compatibility.views.converged"
    )
    extra_alias$decision_ready <-
      extra_alias$compatibility$views$converged
    expect_error(.dpprior_validate_result_v1(extra_alias),
                 class = "dpprior_schema_error")
    extra_deprecation <- fit
    extra_deprecation$compatibility$deprecations$decision_ready <- list(
      code = "forged", authority = "authoritative"
    )
    expect_error(.dpprior_validate_result_v1(extra_deprecation),
                 class = "dpprior_schema_error")

    reordered <- fit
    reordered$compatibility$views <-
      reordered$compatibility$views[rev(names(reordered$compatibility$views))]
    expect_error(.dpprior_validate_result_v1(reordered),
                 class = "dpprior_schema_error")

    escalated <- fit
    escalated$compatibility$views$converged <-
      !escalated$compatibility$views$converged
    escalated$converged <- escalated$compatibility$views$converged
    expect_error(.dpprior_validate_result_v1(escalated),
                 class = "dpprior_schema_error")

    wrong_backend <- fit
    wrong_backend$provenance$backend$implementation <-
      "R/99_forged.R:scientific_consumer"
    expect_error(.dpprior_validate_result_v1(wrong_backend),
                 class = "dpprior_schema_error")

    stripped <- unclass(fit)
    alias_names <- names(stripped$compatibility$top_level_aliases)
    stripped[alias_names] <- NULL
    stripped$compatibility <- .dpprior_new_compatibility()
    class(stripped) <- class(fit)
    expect_error(.dpprior_validate_result_v1(stripped),
                 class = "dpprior_schema_error")
  }

  cross_route <- fits$a1
  cross_route$compatibility$views$a2_kl_v0 <-
    fits$a2_kl_explicit$compatibility$views$a2_kl_v0
  cross_route$compatibility$deprecations$a2_kl_v0 <-
    fits$a2_kl_explicit$compatibility$deprecations$a2_kl_v0
  expect_error(.dpprior_validate_result_v1(cross_route),
               class = "dpprior_schema_error")
})


test_that("candidate ledger binds execution, objectives, counts, and fallback", {
  empty_ledger <- .schema23_fit("a2_kl")
  empty_ledger$computation$candidate_evaluations <- list()
  empty_ledger$computation$selected_candidate_id <- NULL
  expect_error(.dpprior_validate_result_v1(empty_ledger),
               class = "dpprior_schema_error")

  execution_flip <- .schema23_fit("a2_kl")
  candidate <- execution_flip$computation$candidate_evaluations[[1L]]
  candidate$execution_success <- FALSE
  candidate$optimizer_supported <- FALSE
  candidate$decision_eligible <- FALSE
  candidate$outcome <- "selected_diagnostic"
  candidate$rejection_codes <- c("execution_failed", "optimizer_unsupported")
  execution_flip$computation$candidate_evaluations[[1L]] <- candidate
  expect_error(.dpprior_validate_result_v1(execution_flip),
               class = "dpprior_schema_error")

  optimizer_flip <- .schema23_fit("a2_kl")
  candidate <- optimizer_flip$computation$candidate_evaluations[[1L]]
  candidate$optimizer_supported <- FALSE
  candidate$decision_eligible <- FALSE
  candidate$outcome <- "selected_diagnostic"
  candidate$rejection_codes <- "optimizer_unsupported"
  optimizer_flip$computation$candidate_evaluations[[1L]] <- candidate
  expect_error(.dpprior_validate_result_v1(optimizer_flip),
               class = "dpprior_schema_error")

  forged_objective <- .schema23_fit("a2_kl")
  forged_objective$computation$attempts[[1L]]$candidate_objective <- 1
  candidate <- forged_objective$computation$candidate_evaluations[[1L]]
  candidate$recorded_objective <- 1
  candidate$fresh_objective <- 1
  candidate$selection_objective <- 1
  forged_objective$computation$candidate_evaluations[[1L]] <- candidate
  expect_error(.dpprior_validate_result_v1(forged_objective),
               class = "dpprior_schema_error")

  negative_moment_objective <- .schema23_fit("a2_moment")
  negative_moment_objective$computation$attempts[[1L]]$
    candidate_objective <- -999
  negative_moment_objective$computation$candidate_evaluations[[1L]]$
    recorded_objective <- -999
  expect_error(.dpprior_validate_result_v1(negative_moment_objective),
               class = "dpprior_schema_error")

  wrong_count <- .schema23_fit("dual_hard")
  wrong_count$constraint$feasibility$candidate_count <- 999L
  expect_error(.dpprior_validate_result_v1(wrong_count),
               class = "dpprior_schema_error")

  wrong_feasible_count <- .schema23_fit("dual_hard")
  wrong_feasible_count$constraint$feasibility$feasible_candidate_count <- 0L
  expect_error(.dpprior_validate_result_v1(wrong_feasible_count),
               class = "dpprior_schema_error")

  no_execution <- .schema23_no_candidate("dual_hard", "failed")
  no_execution$computation$termination$iterations <- 77L
  expect_error(.dpprior_validate_result_v1(no_execution),
               class = "dpprior_schema_error")

  fallback_mismatch <- .schema23_fit("dual_soft")
  fallback_mismatch$computation$attempts[[1L]]$method <- "Nelder-Mead"
  fallback_mismatch$computation$attempts[[1L]]$stage <- "fallback"
  fallback_mismatch$computation$candidate_evaluations[[1L]]$method <-
    "Nelder-Mead"
  expect_error(.dpprior_validate_result_v1(fallback_mismatch),
               class = "dpprior_schema_error")
})


test_that("unused fallback evidence identifies its trigger and later attempt", {
  fit <- .schema23_fit("a2_kl")
  fallback_attempt <- fit$computation$attempts[[1L]]
  fallback_attempt$id <- "attempt-2"
  fallback_attempt$stage <- "fallback"
  fallback_attempt$method <- "nlminb"
  fallback_attempt$exit_code <- 2L
  fallback_attempt["candidate_parameters"] <- list(NULL)
  fallback_attempt["candidate_objective"] <- list(NULL)
  fallback_attempt$selected <- FALSE
  fallback_attempt$reason_code <- "optimizer_exit_nonzero"
  fallback_attempt$unavailable <- c(
    candidate_parameters = "fallback returned no finite parameters",
    candidate_objective = "fallback returned no finite objective"
  )
  fit$computation$attempts <- c(
    fit$computation$attempts, list(fallback_attempt)
  )
  fit$computation$fallback <- .dpprior_new_fallback(
    attempted = TRUE, used = FALSE, trigger_attempt_id = "attempt-1",
    selected_attempt_id = NULL, reason_code = "primary_selected",
    outcome = "attempted_not_selected"
  )
  expect_invisible(.dpprior_validate_result_v1(fit))

  selected_id <- fit
  selected_id$computation$fallback$selected_attempt_id <- "attempt-2"
  expect_error(.dpprior_validate_result_v1(selected_id),
               class = "dpprior_schema_error")

  no_reason <- fit
  no_reason$computation$fallback["reason_code"] <- list(NULL)
  expect_error(.dpprior_validate_result_v1(no_reason),
               class = "dpprior_schema_error")

  open_outcome <- fit
  open_outcome$computation$fallback$outcome <- "maybe"
  expect_error(.dpprior_validate_result_v1(open_outcome),
               class = "dpprior_schema_error")
})


test_that("objective determinism is scoped to recorded and selection kinds", {
  fit <- .schema23_fit("a2_kl")
  second_attempt <- fit$computation$attempts[[1L]]
  second_attempt$id <- "attempt-2"
  second_attempt$stage <- "optimizer"
  second_attempt$method <- "nlminb"
  second_attempt$selected <- FALSE
  second_attempt$reason_code <- "eligible_not_selected"
  second_attempt$candidate_objective <- 0

  second_candidate <- .dpprior_new_candidate_evaluation(
    id = "candidate-2", attempt_id = "attempt-2", method = "nlminb",
    generator = "direct_attempt", parameters = fit$parameters,
    objective_kind = "kl", recorded_objective_kind = "kl",
    recorded_objective = 0,
    fresh_objective = 0, selection_objective = 0,
    objective_tolerance = 1e-12,
    selected_snapshot = fit$verification$selected_snapshot,
    checks = list(candidate_distribution = .schema23_check(
      value = c(pmf_mass_error = 0, pmf_minimum_violation = 0),
      reference = c(pmf_mass_error = 0, pmf_minimum_violation = 0),
      tolerance = c(
        pmf_mass_error = .TOL_PMF_SUM, pmf_minimum_violation = 0
      ),
      operator = "lte", source = "candidate:candidate-2"
    )),
    execution_success = TRUE, optimizer_supported = TRUE,
    selected = FALSE, source = "second_optimizer_fixture"
  )
  fit$computation$attempts <- c(
    fit$computation$attempts, list(second_attempt)
  )
  fit$computation$candidate_evaluations <- c(
    fit$computation$candidate_evaluations, list(second_candidate)
  )
  expect_invisible(.dpprior_validate_result_v1(fit))

  inconsistent_recorded <- fit
  inconsistent_recorded$computation$attempts[[2L]]$candidate_objective <- 1
  candidate <- inconsistent_recorded$computation$candidate_evaluations[[2L]]
  candidate$recorded_objective <- 1
  candidate$objective_passed <- FALSE
  candidate$rejection_codes <- "recorded_objective_mismatch"
  inconsistent_recorded$computation$candidate_evaluations[[2L]] <- candidate
  expect_error(.dpprior_validate_result_v1(inconsistent_recorded),
               class = "dpprior_schema_error")

  inconsistent_fresh <- fit
  candidate <- inconsistent_fresh$computation$candidate_evaluations[[2L]]
  candidate$fresh_objective <- 1
  candidate$selection_objective <- 1
  inconsistent_fresh$computation$candidate_evaluations[[2L]] <- candidate
  expect_error(.dpprior_validate_result_v1(inconsistent_fresh),
               class = "dpprior_schema_error")
})


test_that("producer-negative candidates remain complete without decision claims", {
  rejected <- .schema23_hard_failed_rejected_candidates()
  expect_invisible(.dpprior_validate_result_v1(rejected))
  expect_length(rejected$computation$attempts, 2L)
  expect_length(rejected$computation$candidate_evaluations, 2L)
  for (candidate in rejected$computation$candidate_evaluations) {
    expect_true(candidate$execution_success)
    expect_false(candidate$selection_eligible)
    expect_identical(
      candidate$rejection_codes, "check_failed:order_stability"
    )
  }
  expect_identical(
    rejected$constraint$feasibility[c(
      "candidate_count", "verified_candidate_count", "feasible_candidate_count"
    )],
    list(
      candidate_count = 2L, verified_candidate_count = 2L,
      feasible_candidate_count = 0L
    )
  )

  heterogeneous <- .schema23_fit("dual_hard")
  target_mean <- heterogeneous$target$K$implied$mean
  mean_scale <- heterogeneous$constraint$optimality$K_scales$mean
  selected_mean <- target_mean + 0.2 * mean_scale
  heterogeneous$achieved$K$mean <- selected_mean
  heterogeneous$residuals$K$mean <- selected_mean - target_mean
  heterogeneous$verification$selected_snapshot$achieved$K$mean <- selected_mean
  heterogeneous$verification$selected_snapshot$residuals$K$mean <-
    selected_mean - target_mean
  heterogeneous$verification$verifier_snapshot$achieved$K$mean <-
    selected_mean + 1e-10
  heterogeneous$verification$verifier_snapshot$residuals$K$mean <-
    selected_mean + 1e-10 - target_mean
  stability <- .dpprior_expected_result_stability(heterogeneous)
  heterogeneous$verification$stability <- .dpprior_new_stability(
    delta = stability$delta, tolerance = stability$tolerance,
    formula = stability$formula, scale_floor = stability$scale_floor,
    source = "independent_verifier"
  )
  heterogeneous$constraint$optimality$selected_K_loss <- 0.04
  heterogeneous$constraint$optimality$minimum_K_loss <- 0.04
  heterogeneous$computation$attempts[[1L]]$candidate_objective <- 0.04
  selected_candidate <-
    heterogeneous$computation$candidate_evaluations[[1L]]
  selected_candidate$recorded_objective <- 0.04
  selected_candidate$fresh_objective <- 0.04
  selected_candidate$selection_objective <- 0.04
  selected_candidate$selected_snapshot <-
    heterogeneous$verification$selected_snapshot
  selected_candidate$verifier_snapshot <-
    heterogeneous$verification$verifier_snapshot
  candidate_stability <- .dpprior_candidate_expected_stability(
    heterogeneous, selected_candidate
  )
  selected_candidate$checks$order_stability <- .schema23_check(
    value = candidate_stability$delta,
    reference = setNames(
      rep(0, length(candidate_stability$delta)),
      names(candidate_stability$delta)
    ),
    tolerance = candidate_stability$tolerance, operator = "lte",
    source = "candidate:candidate-1"
  )
  heterogeneous$computation$candidate_evaluations[[1L]] <-
    selected_candidate

  penalty_parameters <- .dpprior_new_parameters(2.1, 3, "log_ab")
  penalty_attempt <- heterogeneous$computation$attempts[[1L]]
  penalty_attempt$id <- "attempt-penalty"
  penalty_attempt$stage <- "diagnostic"
  penalty_attempt$method <- "penalty_L-BFGS-B_diagnostic"
  penalty_attempt$control$penalty <- 100
  penalty_attempt$candidate_parameters <- penalty_parameters
  penalty_attempt$candidate_objective <- 0
  penalty_attempt$selected <- FALSE
  penalty_attempt$reason_code <- "diagnostic_only"
  penalty_snapshot <- heterogeneous$verification$selected_snapshot
  penalty_snapshot$parameters <- penalty_parameters
  penalty_snapshot$achieved$K$mean <- target_mean
  penalty_snapshot$residuals$K$mean <- 0
  penalty_candidate <- .dpprior_new_candidate_evaluation(
    id = "candidate-penalty", attempt_id = "attempt-penalty",
    method = "penalty_L-BFGS-B_diagnostic",
    generator = "derived_diagnostic_attempt",
    parameters = penalty_parameters, objective_kind = "K_loss",
    recorded_objective_kind = "penalized_diagnostic",
    recorded_objective = 0, fresh_objective = 0,
    selection_objective = NULL, objective_tolerance = 1e-8,
    selected_snapshot = penalty_snapshot,
    checks = list(candidate_domain = .schema23_check(
      value = c(K_support = TRUE, weight_support = TRUE),
      reference = c(K_support = TRUE, weight_support = TRUE),
      tolerance = NULL, operator = "identical",
      source = "candidate:candidate-penalty"
    )),
    execution_success = TRUE, optimizer_supported = FALSE,
    selected = FALSE, source = "penalty_diagnostic_fixture"
  )
  heterogeneous$computation$attempts <- c(
    heterogeneous$computation$attempts, list(penalty_attempt)
  )
  heterogeneous$computation$candidate_evaluations <- c(
    heterogeneous$computation$candidate_evaluations,
    list(penalty_candidate)
  )
  heterogeneous$constraint$feasibility$candidate_count <- 2L
  heterogeneous$constraint$feasibility$verified_candidate_count <- 1L
  heterogeneous$constraint$feasibility$feasible_candidate_count <- 1L
  expect_lt(
    penalty_attempt$candidate_objective,
    heterogeneous$computation$attempts[[1L]]$candidate_objective
  )
  expect_null(penalty_candidate$selection_objective)
  expect_invisible(.dpprior_validate_result_v1(heterogeneous))

  same_parameters_different_recorded_kinds <- heterogeneous
  shared_parameters <- heterogeneous$parameters
  same_parameters_different_recorded_kinds$computation$attempts[[2L]]$
    candidate_parameters <- shared_parameters
  candidate <- same_parameters_different_recorded_kinds$computation$
    candidate_evaluations[[2L]]
  candidate$parameters <- shared_parameters
  candidate$selected_snapshot$parameters <- shared_parameters
  same_parameters_different_recorded_kinds$computation$
    candidate_evaluations[[2L]] <- candidate
  expect_identical(candidate$objective_kind, "K_loss")
  expect_identical(candidate$recorded_objective_kind, "penalized_diagnostic")
  expect_invisible(.dpprior_validate_result_v1(
    same_parameters_different_recorded_kinds
  ))

  scan <- .schema23_hard_approximate_scan()
  expect_invisible(.dpprior_validate_result_v1(scan))
  expect_length(scan$computation$attempts, 0L)
  expect_null(scan$computation$selected_attempt_id)
  expect_true(scan$constraint$satisfied)
  expect_false(scan$usable)
  expect_false(scan$verified)

  false_optimizer <- scan
  candidate <- false_optimizer$computation$candidate_evaluations[[1L]]
  candidate$optimizer_supported <- TRUE
  candidate$decision_eligible <- TRUE
  candidate$outcome <- "selected"
  candidate$rejection_codes <- character()
  false_optimizer$computation$candidate_evaluations[[1L]] <- candidate
  expect_error(.dpprior_validate_result_v1(false_optimizer),
               class = "dpprior_schema_error")

  false_scan_source <- scan
  false_scan_source$computation$termination$source <- "analytic_certificate"
  expect_error(.dpprior_validate_result_v1(false_scan_source),
               class = "dpprior_schema_error")

  parented_scan <- scan
  scan_parent <- .dpprior_new_attempt(
    id = "attempt-scan", stage = "profile_scan",
    method = "deterministic_feasible_profile_scan",
    start = list(log_a = -5), bounds = list(lower = -15, upper = 15),
    control = list(scan_points = 17L), exit_code = 0L,
    message = "deterministic profile scan completed", iterations = 17L,
    evaluations = list(function_count = 17L), candidate_parameters = NULL,
    candidate_objective = NULL, elapsed_seconds = 0.01,
    warnings = character(), error = NULL, selected = FALSE,
    reason_code = "diagnostic_only",
    unavailable = c(
      candidate_parameters = "aggregate scan has per-row candidates",
      candidate_objective = "per-row fresh K-loss is in the candidate ledger"
    )
  )
  parented_scan$computation$attempts <- list(scan_parent)
  parented_scan$computation$candidate_evaluations[[1L]]$attempt_id <-
    "attempt-scan"
  expect_invisible(.dpprior_validate_result_v1(parented_scan))

  aggregate_label <- parented_scan
  aggregate_label$computation$candidate_evaluations[[1L]]$generator <-
    "aggregate_profile_scan"
  expect_error(.dpprior_validate_result_v1(aggregate_label),
               class = "dpprior_schema_error")

  feasibility_parent <- scan
  probe <- .dpprior_new_attempt(
    id = "attempt-probe", stage = "feasibility",
    method = "analytic_monotonicity_feasibility_probe",
    start = list(log_a = -15, log_b = -15),
    bounds = list(lower = c(-15, -15), upper = c(15, 15)),
    control = list(M = 80L), exit_code = 0L,
    message = "analytic feasibility extreme evaluated", iterations = 0L,
    evaluations = list(function_count = 1L),
    candidate_parameters = scan$parameters, candidate_objective = 0.4,
    elapsed_seconds = 0.001, warnings = character(), error = NULL,
    selected = FALSE, reason_code = "diagnostic_only", unavailable = character()
  )
  feasibility_parent$computation$attempts <- list(probe)
  candidate <- feasibility_parent$computation$candidate_evaluations[[1L]]
  candidate$attempt_id <- "attempt-probe"
  candidate$method <- "analytic_feasibility_extreme"
  candidate$generator <- "feasibility_extreme"
  candidate$recorded_objective_kind <- "weight_metric_probe"
  candidate$recorded_objective <- 0.4
  candidate$recorded_objective_available <- TRUE
  candidate$recorded_objective_reason <-
    "source_objective_kind_mismatch:weight_metric_probe->K_loss"
  candidate["objective_passed"] <- list(NULL)
  candidate$rejection_codes <- c(
    "recorded_objective_noncomparable", "optimizer_unsupported"
  )
  feasibility_parent$computation$candidate_evaluations[[1L]] <- candidate
  expect_invisible(.dpprior_validate_result_v1(feasibility_parent))

  wrong_probe_objective <- feasibility_parent
  wrong_probe_objective$computation$attempts[[1L]]$candidate_objective <- 0.5
  expect_error(.dpprior_validate_result_v1(wrong_probe_objective),
               class = "dpprior_schema_error")

  wrong_probe_parent <- feasibility_parent
  wrong_probe_parent$computation$attempts[[1L]]$method <-
    "deterministic_feasible_profile_scan"
  expect_error(.dpprior_validate_result_v1(wrong_probe_parent),
               class = "dpprior_schema_error")

  initializer <- .schema23_a2_kl_approximate_initializer()
  expect_invisible(.dpprior_validate_result_v1(initializer))
  expect_null(
    initializer$computation$attempts[[1L]]$candidate_objective
  )
  expect_true(
    initializer$computation$candidate_evaluations[[1L]]$selection_eligible
  )
  expect_false(
    initializer$computation$candidate_evaluations[[1L]]$optimizer_supported
  )
  wrong_initializer_generator <- initializer
  wrong_initializer_generator$computation$candidate_evaluations[[1L]]$
    generator <- "direct_attempt"
  expect_error(.dpprior_validate_result_v1(wrong_initializer_generator),
               class = "dpprior_schema_error")

  soft_diagnostic <- .schema23_soft_approximate_diagnostic()
  expect_invisible(.dpprior_validate_result_v1(soft_diagnostic))
  expect_false(soft_diagnostic$tradeoff$optimality$passed)
  expect_false(soft_diagnostic$usable)
  expect_false(soft_diagnostic$verified)

  false_soft_source <- soft_diagnostic
  false_soft_source$computation$termination$source <- "analytic_certificate"
  expect_error(.dpprior_validate_result_v1(false_soft_source),
               class = "dpprior_schema_error")

  expect_error(
    .schema23_no_candidate("dual_soft", "infeasible"),
    class = "dpprior_schema_error"
  )
})


test_that("signed hard unsatisfied candidates remain diagnostic, never usable", {
  diagnostic <- .schema23_hard_unsatisfied_diagnostic()
  expect_invisible(.dpprior_validate_result_v1(diagnostic))
  expect_identical(diagnostic$status, "approximate")
  expect_false(diagnostic$usable)
  expect_false(diagnostic$verified)
  expect_false(diagnostic$constraint$satisfied)
  expect_true(diagnostic$constraint$feasibility$feasibility_unknown)
  expect_true(
    diagnostic$computation$candidate_evaluations[[1L]]$diagnostic_eligible
  )
  expect_identical(
    diagnostic$computation$candidate_evaluations[[1L]]$outcome,
    "selected_diagnostic"
  )

  forged <- diagnostic
  forged$computation$candidate_evaluations[[1L]]$diagnostic_eligible <- FALSE
  expect_error(.dpprior_validate_result_v1(forged),
               class = "dpprior_schema_error")

  promoted <- diagnostic
  promoted$status <- "converged"
  promoted$usable <- TRUE
  promoted$verified <- TRUE
  promoted$verification$passed <- TRUE
  promoted$computation$termination$code <- "converged"
  promoted$computation$termination$source <- "optimizer"
  expect_error(.dpprior_validate_result_v1(promoted),
               class = "dpprior_schema_error")

  priority <- .schema23_hard_feasible_beats_diagnostic()
  expect_invisible(.dpprior_validate_result_v1(priority))
  expect_identical(priority$constraint$optimality$selected_K_loss, 1)
  expect_identical(
    priority$computation$candidate_evaluations[[2L]]$selection_objective, 0
  )
  expect_true(
    priority$computation$candidate_evaluations[[2L]]$diagnostic_eligible
  )
  expect_identical(priority$computation$selected_candidate_id, "candidate-1")
  expect_identical(
    priority$constraint$feasibility[c(
      "candidate_count", "verified_candidate_count", "feasible_candidate_count"
    )],
    list(
      candidate_count = 2L, verified_candidate_count = 2L,
      feasible_candidate_count = 1L
    )
  )
})


test_that("truth-making tolerances and scales cannot be inflated in place", {
  hard_constraint <- .schema23_hard_unsatisfied_diagnostic()
  hard_constraint$computation$request$controls$constraint_abs_tol <- 1
  hard_constraint$computation$request$controls$constraint_rel_tol <- 0
  hard_constraint$computation$used$controls$constraint_abs_tol <- 1
  hard_constraint$computation$used$controls$constraint_rel_tol <- 0
  hard_constraint$tolerances$constraint <- list(
    absolute = 1, relative = 0, effective = 1
  )
  hard_constraint$constraint$tolerance <-
    hard_constraint$tolerances$constraint
  for (snapshot_name in c("selected_snapshot", "verifier_snapshot")) {
    hard_constraint$verification[[snapshot_name]]$tolerances <-
      hard_constraint$tolerances
    hard_constraint$computation$candidate_evaluations[[1L]][[
      snapshot_name
    ]]$tolerances <- hard_constraint$tolerances
  }
  expect_error(.dpprior_validate_result_v1(hard_constraint),
               class = "dpprior_schema_error")

  hard_perturbation <- .schema23_fit("dual_hard")
  hard_perturbation$computation$request$controls$perturbation_abs_tol <- 999
  hard_perturbation$computation$used$controls$perturbation_abs_tol <- 999
  hard_perturbation$tolerances$perturbation$absolute <- 999
  for (snapshot_name in c("selected_snapshot", "verifier_snapshot")) {
    hard_perturbation$verification[[snapshot_name]]$tolerances <-
      hard_perturbation$tolerances
    hard_perturbation$computation$candidate_evaluations[[1L]][[
      snapshot_name
    ]]$tolerances <- hard_perturbation$tolerances
  }
  hard_perturbation$computation$candidate_evaluations[[1L]]$checks$
    perturbation <- .schema23_check(
      value = c(maximum_metric_delta = 999), reference = 0,
      tolerance = c(maximum_metric_delta = 999), operator = "lte",
      source = "candidate:candidate-1"
    )
  expect_error(.dpprior_validate_result_v1(hard_perturbation),
               class = "dpprior_schema_error")

  hard_scales <- .schema23_fit("dual_hard")
  reversed <- list(K = list(mean = 2, variance = 5))
  hard_scales$computation$scaling$requested <- reversed
  hard_scales$computation$scaling$used <- reversed
  hard_scales$computation$scaling$values <- reversed
  hard_scales$constraint$optimality$K_scales <- reversed$K
  expect_error(.dpprior_validate_result_v1(hard_scales),
               class = "dpprior_schema_error")

  hard_tie <- .schema23_fit("dual_hard")
  hard_tie$constraint$optimality$tie_tolerance <- 1
  hard_tie$verification$components$candidate_selection$tolerance <- 1
  expect_error(.dpprior_validate_result_v1(hard_tie),
               class = "dpprior_schema_error")

  certificate_gap <- .schema23_no_candidate("dual_hard", "infeasible")
  certificate <- certificate_gap$constraint$feasibility$certificate
  inside_value <- certificate$target_value +
    certificate$tolerance / 2 + 64 * .Machine$double.eps
  certificate$minimum$selected$value <- inside_value
  certificate$minimum$refined$value <- inside_value
  certificate$minimum$uncertainty <- 64 * .Machine$double.eps
  certificate$minimum$lower <- inside_value - 64 * .Machine$double.eps
  certificate$minimum$upper <- inside_value + 64 * .Machine$double.eps
  certificate$lower_bound <- certificate$minimum$lower
  certificate_gap$constraint$feasibility$certificate <- certificate
  expect_error(.dpprior_validate_result_v1(certificate_gap),
               class = "dpprior_schema_error")

  soft_gradient <- .schema23_fit("dual_soft")
  soft_gradient$computation$request$controls$stationarity_tol <- 1000
  soft_gradient$computation$used$controls$stationarity_tol <- 1000
  soft_gradient$tolerances$stationarity$tolerance <- 1000
  soft_gradient$tradeoff$optimality$stationarity_tolerance <- 1000
  soft_gradient$tradeoff$optimality$gradient[] <- c(999, -999)
  soft_gradient$tradeoff$optimality$component_pass[] <- TRUE
  soft_gradient$tradeoff$optimality$stationarity_passed <- TRUE
  expect_error(.dpprior_validate_result_v1(soft_gradient),
               class = "dpprior_schema_error")

  soft_objective <- .schema23_fit("dual_soft")
  soft_objective$computation$request$controls$objective_abs_tol <- 100
  soft_objective$computation$request$controls$objective_rel_tol <- 0
  soft_objective$computation$used$controls$objective_abs_tol <- 100
  soft_objective$computation$used$controls$objective_rel_tol <- 0
  soft_objective$tolerances$objective <- list(
    absolute = 100, relative = 0, scale_floor = 1
  )
  soft_objective$tradeoff$optimality$recorded_objective <- 100
  soft_objective$tradeoff$optimality$objective_difference <- 100
  soft_objective$tradeoff$optimality$objective_tolerance <- 100
  soft_objective$tradeoff$optimality$objective_passed <- TRUE
  expect_error(.dpprior_validate_result_v1(soft_objective),
               class = "dpprior_schema_error")

  soft_negative_start <- .schema23_fit("dual_soft")
  soft_negative_start$tradeoff$optimality$start_objective <- -100
  expect_error(.dpprior_validate_result_v1(soft_negative_start),
               class = "dpprior_schema_error")

  soft_negative_neighbor <- .schema23_fit("dual_soft")
  soft_negative_neighbor$tradeoff$optimality$neighbor_objectives[[1L]] <- -100
  expect_error(.dpprior_validate_result_v1(soft_negative_neighbor),
               class = "dpprior_schema_error")

  soft_selection <- .schema23_fit("dual_soft")
  soft_selection$computation$request$controls$selection_tolerance <- 1
  soft_selection$computation$used$controls$selection_tolerance <- 1
  soft_selection$tolerances$selection <- 1
  soft_selection$tradeoff$optimality$selection_tolerance <- 1
  soft_selection$verification$components$candidate_selection$tolerance <- 1
  expect_error(.dpprior_validate_result_v1(soft_selection),
               class = "dpprior_schema_error")

  moment <- .schema23_fit("a2_moment")
  moment$tolerances$K_adequacy$absolute <- 100
  moment$verification$selected_snapshot$tolerances <- moment$tolerances
  moment$verification$verifier_snapshot$tolerances <- moment$tolerances
  moment$computation$candidate_evaluations[[1L]]$selected_snapshot$tolerances <-
    moment$tolerances
  expect_error(.dpprior_validate_result_v1(moment),
               class = "dpprior_schema_error")

  kl <- .schema23_fit("a2_kl")
  kl$tolerances$distribution$adequacy$kl <- 100
  kl$verification$selected_snapshot$tolerances <- kl$tolerances
  kl$verification$verifier_snapshot$tolerances <- kl$tolerances
  kl$computation$candidate_evaluations[[1L]]$selected_snapshot$tolerances <-
    kl$tolerances
  expect_error(.dpprior_validate_result_v1(kl),
               class = "dpprior_schema_error")
})


test_that("extensions, resources, and traces are serialization-safe evidence", {
  proxy <- .schema23_fit("a1_proxy")
  proxy$proxy$mapping$evil <- new.env(parent = emptyenv())
  expect_error(.dpprior_validate_result_v1(proxy),
               class = "dpprior_schema_error")

  diagnostics <- .schema23_diagnostics()
  diagnostics$diagnostics$K$evil <- new.env(parent = emptyenv())
  expect_error(.dpprior_validate_result_v1(diagnostics),
               class = "dpprior_schema_error")

  resources <- .schema23_fit()
  resources$computation$resources$evil <- new.env(parent = emptyenv())
  expect_error(.dpprior_validate_result_v1(resources),
               class = "dpprior_schema_error")

  trace <- .schema23_fit()
  trace$computation$trace <- data.frame(x = 1)
  trace$computation$trace$evil <- I(list(new.env(parent = emptyenv())))
  expect_error(.dpprior_validate_result_v1(trace),
               class = "dpprior_schema_error")

  duplicate_trace <- .schema23_fit()
  duplicate_trace$computation$trace <- data.frame(x = 1, y = 2)
  names(duplicate_trace$computation$trace) <- c("x", "x")
  expect_error(.dpprior_validate_result_v1(duplicate_trace),
               class = "dpprior_schema_error")

  top_attribute <- .schema23_fit()
  attr(top_attribute, "payload") <- new.env(parent = emptyenv())
  expect_error(.dpprior_validate_result_v1(top_attribute),
               class = "dpprior_schema_error")

  attempts_attribute <- .schema23_fit()
  attr(attempts_attribute$computation$attempts, "payload") <-
    new.env(parent = emptyenv())
  expect_error(.dpprior_validate_result_v1(attempts_attribute),
               class = "dpprior_schema_error")

  pmf_attribute <- .schema23_fit("a2_kl")
  attr(pmf_attribute$target$K$pmf, "payload") <-
    new.env(parent = emptyenv())
  expect_error(.dpprior_validate_result_v1(pmf_attribute),
               class = "dpprior_schema_error")

  duplicate_pmf_names <- .schema23_fit("a2_kl")
  names(duplicate_pmf_names$target$K$pmf) <- rep("mass", 20L)
  expect_error(.dpprior_validate_result_v1(duplicate_pmf_names),
               class = "dpprior_schema_error")

  duplicate_support_names <- .schema23_no_candidate(
    "dual_hard", "infeasible"
  )
  names(duplicate_support_names$constraint$feasibility$certificate$support) <-
    c("bound", "bound")
  expect_error(.dpprior_validate_result_v1(duplicate_support_names),
               class = "dpprior_schema_error")
})


test_that("weight certification uses the frozen metric-specific contract", {
  certified <- .schema23_weight_target("at_most", "wmax_tail_upper")
  certified$certification$method <- "fabricated_bound"
  expect_error(.dpprior_validate_weight_target_v1(certified),
               class = "dpprior_schema_error")

  ordinary <- .schema23_weight_target()
  ordinary$certification <- list(
    kind = "upper_bound", certified = TRUE, passed = TRUE,
    method = "certified_size_biased_mass_upper_bound",
    source = "wmax_tail_bounds"
  )
  expect_error(.dpprior_validate_weight_target_v1(ordinary),
               class = "dpprior_schema_error")
})


test_that("all-passing candidate checks do not mint phantom rejection codes", {
  candidate <- .schema23_fit("a2_moment")$computation$
    candidate_evaluations[[1L]]
  expect_identical(candidate$rejection_codes, character())

  forged <- candidate
  forged$rejection_codes <- "check_failed:"
  expect_error(
    .dpprior_validate_candidate_evaluation(forged),
    class = "dpprior_schema_error"
  )
})


test_that("diagnostic truth and aggregation cannot be self-certified", {
  diagnostics <- .schema23_diagnostics()
  expect_invisible(.dpprior_validate_result_v1(diagnostics))

  coordinated_K_forgery <- diagnostics
  coordinated_K_forgery$diagnostics$K$mean <- 999
  coordinated_K_forgery$achieved$K$mean <- 999
  coordinated_K_forgery$verification$selected_snapshot$achieved$K$mean <- 999
  expect_error(.dpprior_validate_result_v1(coordinated_K_forgery),
               class = "dpprior_schema_error")

  coordinated_rho_forgery <- diagnostics
  coordinated_rho_forgery$diagnostics$coclustering$mean <- 0.01
  coordinated_rho_forgery$achieved$coclustering$mean <- 0.01
  coordinated_rho_forgery$verification$selected_snapshot$achieved$
    coclustering$mean <- 0.01
  expect_error(.dpprior_validate_result_v1(coordinated_rho_forgery),
               class = "dpprior_schema_error")

  inflated_refinement <- diagnostics
  inflated_refinement$tolerances$diagnostics$refinement[] <- 1
  inflated_refinement$verification$selected_snapshot$tolerances <-
    inflated_refinement$tolerances
  inflated_refinement$verification$verifier_snapshot$tolerances <-
    inflated_refinement$tolerances
  expect_error(.dpprior_validate_result_v1(inflated_refinement),
               class = "dpprior_schema_error")

  reordered_attempts <- diagnostics
  reordered_attempts$computation$attempts <-
    rev(reordered_attempts$computation$attempts)
  expect_error(.dpprior_validate_result_v1(reordered_attempts),
               class = "dpprior_schema_error")

  false_execution <- diagnostics
  false_execution$computation$attempts[[1L]]$exit_code <- 1L
  expect_error(.dpprior_validate_result_v1(false_execution),
               class = "dpprior_schema_error")

  wrong_status <- diagnostics
  wrong_status$diagnostics$alpha$status <- "boundary"
  expect_error(.dpprior_validate_result_v1(wrong_status),
               class = "dpprior_schema_error")

  self_report <- diagnostics
  self_report$verification$verifier_snapshot$source <- "optimizer_self_report"
  expect_error(.dpprior_validate_result_v1(self_report),
               class = "dpprior_schema_error")

  placeholder_invariant <- diagnostics
  placeholder_invariant$verification$invariants$fixed_parameters$source <-
    "placeholder"
  expect_error(.dpprior_validate_result_v1(placeholder_invariant),
               class = "dpprior_schema_error")

  false_policy <- diagnostics
  false_policy$diagnostics$policy_results <- list(list(
    estimand = "W_SB", direction = "above", threshold = 0.5,
    value = 0.9, lower = NULL, upper = NULL,
    outcome = "not_triggered", basis = "exact_tail_probability"
  ))
  expect_error(.dpprior_validate_result_v1(false_policy),
               class = "dpprior_schema_error")
})


test_that("sensitivity reconciliation binds every key and metric claim", {
  sensitivity <- .schema23_sensitivity()
  expect_invisible(.dpprior_validate_result_v1(sensitivity))

  duplicate_key <- sensitivity
  duplicate_key$sensitivity$scenarios <- rbind(
    duplicate_key$sensitivity$scenarios,
    duplicate_key$sensitivity$scenarios
  )
  expect_error(.dpprior_validate_result_v1(duplicate_key),
               class = "dpprior_schema_error")

  forged_content_key <- sensitivity
  forged_content_key$sensitivity$scenarios$canonical_content <-
    "J=999|method=A1|diagnostics=FALSE"
  expect_error(.dpprior_validate_result_v1(forged_content_key),
               class = "dpprior_schema_error")

  missing_metric <- sensitivity
  missing_metric$sensitivity$metrics_long <-
    missing_metric$sensitivity$metrics_long[-1L, , drop = FALSE]
  expect_error(.dpprior_validate_result_v1(missing_metric),
               class = "dpprior_schema_error")

  infinite_metric <- sensitivity
  infinite_metric$sensitivity$metrics_long$value[[1L]] <- Inf
  infinite_metric$sensitivity$scenario_results$a[[1L]] <- Inf
  expect_error(.dpprior_validate_result_v1(infinite_metric),
               class = "dpprior_schema_error")

  unavailable_without_reason <- sensitivity
  unavailable_without_reason$sensitivity$metrics_long$value[[1L]] <- NA_real_
  unavailable_without_reason$sensitivity$scenario_results$a[[1L]] <- NA_real_
  expect_error(.dpprior_validate_result_v1(unavailable_without_reason),
               class = "dpprior_schema_error")

  wide_long_mismatch <- sensitivity
  wide_long_mismatch$sensitivity$metrics_long$value[[1L]] <- 999
  expect_error(.dpprior_validate_result_v1(wide_long_mismatch),
               class = "dpprior_schema_error")

  coordinated_alpha_forgery <- sensitivity
  alpha_row <- match(
    "E_alpha", coordinated_alpha_forgery$sensitivity$metrics_long$metric
  )
  coordinated_alpha_forgery$sensitivity$scenario_results$E_alpha <- 999
  coordinated_alpha_forgery$sensitivity$metrics_long$value[[alpha_row]] <- 999
  expect_error(.dpprior_validate_result_v1(coordinated_alpha_forgery),
               class = "dpprior_schema_error")

  wrong_metric_source <- sensitivity
  wrong_metric_source$sensitivity$metrics_long$source[[1L]] <- "rubber_stamp"
  expect_error(.dpprior_validate_result_v1(wrong_metric_source),
               class = "dpprior_schema_error")

  malformed_condition <- sensitivity
  key <- malformed_condition$sensitivity$scenarios$scenario_key[[1L]]
  malformed_condition$sensitivity$conditions[[key]]$calibration <-
    list(class = "condition")
  expect_error(.dpprior_validate_result_v1(malformed_condition),
               class = "dpprior_schema_error")

  malformed_interval <- sensitivity
  key <- malformed_interval$sensitivity$scenarios$scenario_key[[1L]]
  malformed_interval$sensitivity$interval_checks[[key]]$reason <- "fabricated"
  expect_error(.dpprior_validate_result_v1(malformed_interval),
               class = "dpprior_schema_error")

  fabricated_interval <- sensitivity
  key <- fabricated_interval$sensitivity$scenarios$scenario_key[[1L]]
  fabricated_interval$sensitivity$interval_checks[[key]] <- list(
    requested = list(
      lower = 1, upper = 5, type = "hard_bounds", coverage = 1,
      family = "maxent"
    ),
    selected = list(
      coverage = 1, lower_tail = 0, upper_tail = 0,
      coverage_residual = 0
    ),
    verification = NULL, status = "approximate", source = "rubber_stamp",
    usable = FALSE, verified = FALSE, reason = "fabricated"
  )
  expect_error(.dpprior_validate_result_v1(fabricated_interval),
               class = "dpprior_schema_error")

  wrong_local_shape <- sensitivity
  wrong_local_shape$sensitivity$local <- list()
  expect_error(.dpprior_validate_result_v1(wrong_local_shape),
               class = "dpprior_schema_error")

  fabricated_local <- sensitivity
  key <- fabricated_local$sensitivity$scenarios$scenario_key[[1L]]
  fabricated_local$sensitivity$local <- data.frame(
    scenario_key = key, axis = "mu_K", axis_value = 5,
    settings_key = "settings", lower_scenario_key = key,
    upper_scenario_key = key, lower_value = 4, upper_value = 6,
    metric = "E_alpha", component = "alpha", derivative = 999,
    method = "bracketed_secant_across_nearest_same-setting_neighbors",
    reason = NA_character_, stringsAsFactors = FALSE
  )
  expect_error(.dpprior_validate_result_v1(fabricated_local),
               class = "dpprior_schema_error")

  wrong_global <- sensitivity
  wrong_global$sensitivity$global$scenario_count <- 999L
  expect_error(.dpprior_validate_result_v1(wrong_global),
               class = "dpprior_schema_error")

  wrong_row_precedence <- sensitivity
  wrong_row_precedence$sensitivity$scenario_results$status <- "boundary"
  wrong_row_precedence$sensitivity$metrics_long$status <- "boundary"
  expect_error(.dpprior_validate_result_v1(wrong_row_precedence),
               class = "dpprior_schema_error")

  self_report <- sensitivity
  self_report$verification$verifier_snapshot$source <- "optimizer_self_report"
  expect_error(.dpprior_validate_result_v1(self_report),
               class = "dpprior_schema_error")

  placeholder <- sensitivity
  placeholder$verification$components$reconciliation$source <- "placeholder"
  expect_error(.dpprior_validate_result_v1(placeholder),
               class = "dpprior_schema_error")
})


test_that("fit verifier sources and candidate orders are authoritative", {
  for (mode in c("a2_moment", "a2_kl", "dual_hard", "dual_soft")) {
    self_report <- .schema23_fit(mode)
    self_report$verification$verifier_snapshot$source <-
      "optimizer_self_report"
    expect_error(.dpprior_validate_result_v1(self_report),
                 class = "dpprior_schema_error")

    placeholder_method <- .schema23_fit(mode)
    placeholder_method$verification$method <- "placeholder_verification"
    expect_error(.dpprior_validate_result_v1(placeholder_method),
                 class = "dpprior_schema_error")

    placeholder_invariant <- .schema23_fit(mode)
    placeholder_invariant$verification$invariants[[1L]]$source <- "placeholder"
    expect_error(.dpprior_validate_result_v1(placeholder_invariant),
                 class = "dpprior_schema_error")

    wrong_selected_M <- .schema23_fit(mode)
    wrong_selected_M$computation$candidate_evaluations[[1L]]$
      selected_snapshot$M <- 1L
    expect_error(.dpprior_validate_result_v1(wrong_selected_M),
                 class = "dpprior_schema_error")

    wrong_selected_K_M <- .schema23_fit(mode)
    wrong_selected_K_M$computation$candidate_evaluations[[1L]]$
      selected_snapshot$achieved$K$M <- 1L
    expect_error(.dpprior_validate_result_v1(wrong_selected_K_M),
                 class = "dpprior_schema_error")

    wrong_verifier_M <- .schema23_fit(mode)
    wrong_verifier_M$computation$candidate_evaluations[[1L]]$
      verifier_snapshot$M <- 999L
    expect_error(.dpprior_validate_result_v1(wrong_verifier_M),
                 class = "dpprior_schema_error")

    wrong_public_M <- .schema23_fit(mode)
    wrong_public_M$verification$selected_snapshot$M <- 1L
    expect_error(.dpprior_validate_result_v1(wrong_public_M),
                 class = "dpprior_schema_error")

    for (forged_M in list(NULL, 1L, 999L)) {
      wrong_verifier_K_M <- .schema23_fit(mode)
      wrong_verifier_K_M$verification$verifier_snapshot$achieved$K[
        "M"
      ] <- list(forged_M)
      expect_error(.dpprior_validate_result_v1(wrong_verifier_K_M),
                   class = "dpprior_schema_error")
    }
  }
})


test_that("dual input identity and retained no-candidate facts stay coherent", {
  for (mode in c("dual_hard", "dual_soft")) {
    wrong_J <- .schema23_fit(mode)
    wrong_J$provenance$input_fit$J <- 11L
    wrong_J$provenance$input_fit$target$J <- 11L
    expect_error(.dpprior_validate_result_v1(wrong_J),
                 class = "dpprior_schema_error")

    wrong_parameters <- .schema23_fit(mode)
    replacement <- .dpprior_new_parameters(2.1, 3, "log_ab")
    wrong_parameters$provenance$input_fit$parameters <- replacement
    wrong_parameters$provenance$input_fit$selected_snapshot$parameters <-
      replacement
    expect_error(.dpprior_validate_result_v1(wrong_parameters),
                 class = "dpprior_schema_error")

    coordinated_parameters <- .schema23_fit(mode)
    coordinated_parameters$provenance$input_fit$parameters <- replacement
    coordinated_parameters$provenance$input_fit$selected_snapshot$parameters <-
      replacement
    coordinated_parameters$provenance$input_fit$verifier_snapshot$parameters <-
      replacement
    expect_error(.dpprior_validate_result_v1(coordinated_parameters),
                 class = "dpprior_schema_error")

    swapped_sources <- .schema23_fit(mode)
    swapped_sources$provenance$input_fit$selected_snapshot$source <-
      "independent_verifier"
    swapped_sources$provenance$input_fit$verifier_snapshot$source <-
      "selected_order"
    expect_error(.dpprior_validate_result_v1(swapped_sources),
                 class = "dpprior_schema_error")

    forged_orders <- .schema23_fit(mode)
    forged_orders$provenance$input_fit$selected_snapshot$M <- 999L
    forged_orders$provenance$input_fit$selected_snapshot$achieved_K$M <- 999L
    forged_orders$provenance$input_fit$verifier_snapshot$M <- 1998L
    forged_orders$provenance$input_fit$verifier_snapshot$achieved_K$M <- 1998L
    expect_error(.dpprior_validate_result_v1(forged_orders),
                 class = "dpprior_schema_error")
  }

  contradictory_attempt <- .schema23_hard_failed_rejected_candidates()
  contradictory_attempt$computation$attempts[[1L]][
    "candidate_parameters"
  ] <- list(NULL)
  contradictory_attempt$computation$attempts[[1L]]$reason_code <-
    "diagnostic_only"
  contradictory_attempt$computation$attempts[[1L]]$unavailable <- c(
    candidate_parameters = "fabricated missing candidate"
  )
  expect_error(.dpprior_validate_result_v1(contradictory_attempt),
               class = "dpprior_schema_error")

  satisfied_without_invariant <- .schema23_hard_approximate_scan()
  source <- satisfied_without_invariant$verification$invariants[[1L]]$source
  satisfied_without_invariant$verification$invariants[[1L]] <-
    .schema23_check(
      value = FALSE, reference = TRUE, tolerance = NULL,
      operator = "identical", source = source
    )
  expect_error(.dpprior_validate_result_v1(satisfied_without_invariant),
               class = "dpprior_schema_error")
})


test_that("hard log bounds and input-fit decision authority are exact", {
  default_hard <- .schema23_fit("dual_hard")
  expect_identical(
    default_hard$computation$used$controls$log_bounds, c(-15, 15)
  )
  expect_invisible(.dpprior_validate_result_v1(default_hard))

  custom_hard <- default_hard
  custom_bounds <- c(-2, 2)
  custom_hard$computation$request$controls$log_bounds <- custom_bounds
  custom_hard$computation$used$controls$log_bounds <- custom_bounds
  custom_hard$verification$settings$log_bounds <- custom_bounds
  expect_invisible(.dpprior_validate_result_v1(custom_hard))

  request_mismatch <- custom_hard
  request_mismatch$computation$request$controls$log_bounds <- c(-3, 2)
  expect_error(
    .dpprior_validate_result_v1(request_mismatch),
    class = "dpprior_schema_error"
  )
  settings_mismatch <- custom_hard
  settings_mismatch$verification$settings$log_bounds <- c(-3, 2)
  expect_error(
    .dpprior_validate_result_v1(settings_mismatch),
    class = "dpprior_schema_error"
  )
  invalid_bounds <- custom_hard
  invalid_bounds$computation$request$controls$log_bounds <-
    invalid_bounds$computation$used$controls$log_bounds <- c(-Inf, 2)
  expect_error(
    .dpprior_validate_result_v1(invalid_bounds),
    class = "dpprior_schema_error"
  )
  forged_domain_invariant <- custom_hard
  forged_domain_invariant$verification$invariants$
    finite_parameters_inside_domain$value <- FALSE
  forged_domain_invariant$verification$invariants$
    finite_parameters_inside_domain$passed <- FALSE
  expect_error(
    .dpprior_validate_result_v1(forged_domain_invariant),
    class = "dpprior_schema_error"
  )

  certificate <- .schema23_no_candidate("dual_hard", "infeasible")
  expect_identical(
    certificate$constraint$feasibility$certificate$domain$log_a,
    certificate$computation$used$controls$log_bounds
  )
  expect_invisible(.dpprior_validate_result_v1(certificate))
  certificate_domain_mismatch <- certificate
  certificate_domain_mismatch$computation$request$controls$log_bounds <-
    certificate_domain_mismatch$computation$used$controls$log_bounds <-
      c(-2, 0)
  expect_error(
    .dpprior_validate_result_v1(certificate_domain_mismatch),
    class = "dpprior_schema_error"
  )

  MN_80_200 <- .schema23_with_A2_MN_input_orders(
    .schema23_fit("dual_hard"), 80L, 200L
  )
  MN_40_160 <- .schema23_with_A2_MN_input_orders(
    .schema23_fit("dual_hard"), 40L, 160L
  )
  expect_invisible(.dpprior_validate_result_v1(MN_80_200))
  expect_invisible(.dpprior_validate_result_v1(MN_40_160))

  below_required <- MN_80_200
  below_required$provenance$input_fit$verifier_snapshot$M <- 159L
  below_required$provenance$input_fit$verifier_snapshot$achieved_K <-
    .schema23_input_K_record(
      below_required$J, below_required$provenance$input_fit$parameters, 159L
    )
  expect_error(
    .dpprior_validate_result_v1(below_required),
    class = "dpprior_schema_error"
  )
  above_ceiling <- MN_80_200
  above_ceiling$provenance$input_fit$verifier_snapshot$M <- 513L
  above_ceiling$provenance$input_fit$verifier_snapshot$achieved_K$M <- 513L
  expect_error(
    .dpprior_validate_result_v1(above_ceiling),
    class = "dpprior_schema_error"
  )
  forged_actual_order <- MN_80_200
  forged_actual_order$provenance$input_fit$verifier_snapshot$achieved_K$mean <-
    forged_actual_order$provenance$input_fit$verifier_snapshot$achieved_K$mean +
      0.1
  expect_error(
    .dpprior_validate_result_v1(forged_actual_order),
    class = "dpprior_schema_error"
  )

  KL_input <- .schema23_with_A2_KL_input_fit(
    .schema23_fit("dual_hard"), 80L, 200L
  )
  expect_invisible(.dpprior_validate_result_v1(KL_input))
  expect_identical(
    KL_input$provenance$input_fit$decision_evidence$target_K,
    KL_input$target$K
  )
  KL_roundtrip <- unserialize(serialize(KL_input, NULL, version = 3L))
  expect_identical(KL_roundtrip, KL_input)
  expect_invisible(.dpprior_validate_result_v1(KL_roundtrip))

  missing_decision <- KL_input
  missing_decision$provenance$input_fit["decision_evidence"] <- list(NULL)
  expect_error(
    .dpprior_validate_result_v1(missing_decision),
    class = "dpprior_schema_error"
  )
  non_KL_decision <- MN_80_200
  non_KL_decision$provenance$input_fit$decision_evidence <- list(
    target_K = non_KL_decision$target$K,
    distribution_tolerances = .schema23_fit_parts("A2-KL")$
      tolerances$distribution
  )
  expect_error(
    .dpprior_validate_result_v1(non_KL_decision),
    class = "dpprior_schema_error"
  )
  inflated_adequacy <- KL_input
  inflated_adequacy$provenance$input_fit$decision_evidence$
    distribution_tolerances$adequacy$kl <- 1
  expect_error(
    .dpprior_validate_result_v1(inflated_adequacy),
    class = "dpprior_schema_error"
  )
  inflated_order <- KL_input
  inflated_order$provenance$input_fit$decision_evidence$
    distribution_tolerances$order$pmf_relative <- 1e-6
  inflated_order$provenance$input_fit$decision_evidence$
    distribution_tolerances$order$pmf_l1 <- 1e-10 + 1e-6
  inflated_order$provenance$input_fit$decision_evidence$
    distribution_tolerances$order$direct_moment_relative <- 1e-6
  expect_error(
    .dpprior_validate_result_v1(inflated_order),
    class = "dpprior_schema_error"
  )
  forged_PMF <- KL_input
  changed <- forged_PMF$provenance$input_fit$selected_snapshot$achieved_K$pmf
  changed[[1L]] <- changed[[1L]] + 1e-4
  changed[[2L]] <- changed[[2L]] - 1e-4
  changed_moments <- .dpprior_target_pmf_moments(changed)
  forged_PMF$provenance$input_fit$selected_snapshot$achieved_K$pmf <- changed
  forged_PMF$provenance$input_fit$selected_snapshot$achieved_K$mean <-
    unname(changed_moments[["mean"]])
  forged_PMF$provenance$input_fit$selected_snapshot$achieved_K$variance <-
    unname(changed_moments[["variance"]])
  expect_error(
    .dpprior_validate_result_v1(forged_PMF),
    class = "dpprior_schema_error"
  )
  unbound_full_target <- KL_input
  unbound_full_target$provenance$input_fit$decision_evidence$target_K$message <-
    "scientifically unbound copy"
  expect_error(
    .dpprior_validate_result_v1(unbound_full_target),
    class = "dpprior_schema_error"
  )
  unclassed_decision_target <- KL_input
  unclassed_decision_target$provenance$input_fit$decision_evidence$target_K <-
    unclass(unclassed_decision_target$provenance$input_fit$
      decision_evidence$target_K)
  expect_error(
    .dpprior_validate_result_v1(unclassed_decision_target),
    class = "dpprior_schema_error"
  )

  endpoint <- .schema23_soft_endpoint()
  expect_null(endpoint$provenance$input_fit$decision_evidence)
  forged_endpoint <- endpoint
  forged_endpoint$provenance$input_fit$decision_evidence <- list(
    target_K = endpoint$target$K,
    distribution_tolerances = .schema23_fit_parts("A2-KL")$
      tolerances$distribution
  )
  expect_error(
    .dpprior_validate_result_v1(forged_endpoint),
    class = "dpprior_schema_error"
  )
})


test_that("direct weight, objective-kind, and soft-neighbor vocabularies close", {
  for (alias in c("bound", "mode", "prob")) {
    target <- .schema23_weight_target()
    target$request[[alias]] <- if (identical(alias, "mode")) "soft" else 0.4
    target$normalized <- target$request
    target$used <- target$request
    expect_error(.dpprior_validate_weight_target_v1(target),
                 class = "dpprior_schema_error")
  }

  wrong_kind <- .schema23_fit("dual_hard")
  wrong_kind$computation$candidate_evaluations[[1L]]$
    recorded_objective_kind <- "penalized_diagnostic"
  expect_error(.dpprior_validate_result_v1(wrong_kind),
               class = "dpprior_schema_error")

  missing_neighbor <- .schema23_fit("dual_soft")
  missing_neighbor$tradeoff$optimality$neighbor_objectives <-
    missing_neighbor$tradeoff$optimality$neighbor_objectives[-1L]
  expect_error(.dpprior_validate_result_v1(missing_neighbor),
               class = "dpprior_schema_error")

  renamed_neighbor <- .schema23_fit("dual_soft")
  names(renamed_neighbor$tradeoff$optimality$neighbor_objectives)[[1L]] <-
    "fabricated_direction"
  expect_error(.dpprior_validate_result_v1(renamed_neighbor),
               class = "dpprior_schema_error")

  arbitrary_start <- .schema23_fit("dual_soft")
  arbitrary_start$tradeoff$optimality$start_objective <-
    arbitrary_start$tradeoff$optimality$start_objective + 1
  expect_error(.dpprior_validate_result_v1(arbitrary_start),
               class = "dpprior_schema_error")
})


test_that("diagnostic authority fails closed and preserves approximate truth", {
  expect_invisible(
    .dpprior_validate_result_v1(.schema23_diagnostics_approximate())
  )

  missing_parameters <- .schema23_diagnostics()
  missing_parameters$parameters <- NULL
  expect_error(.dpprior_validate_result_v1(missing_parameters),
               class = "dpprior_schema_error")

  missing_selected_M <- .schema23_diagnostics()
  missing_selected_M$computation$orders$M_selected <- NULL
  expect_error(.dpprior_validate_result_v1(missing_selected_M),
               class = "dpprior_schema_error")

  missing_verifier_M <- .schema23_diagnostics()
  missing_verifier_M$computation$orders$M_verification_used <- NULL
  expect_error(.dpprior_validate_result_v1(missing_verifier_M),
               class = "dpprior_schema_error")

  forbidden_target <- .schema23_diagnostics()
  forbidden_target$target$dominance_risk <- "high"
  expect_error(.dpprior_validate_result_v1(forbidden_target),
               class = "dpprior_schema_error")

  policy <- .schema23_diagnostics()
  weight_threshold <- 0.5
  action_threshold <- 0.1
  value <- as.numeric(.diagnostic_wsb_tail(weight_threshold, 2, 3))
  outcome <- if (value > action_threshold) "triggered" else "not_triggered"
  policy$target$warning_policy <- list(
    estimand = "W_SB", direction = "above",
    weight_threshold = weight_threshold,
    action_threshold = action_threshold
  )
  policy$diagnostics$policy_results <- list(list(
    estimand = "W_SB", direction = "above", threshold = action_threshold,
    value = value, lower = NULL, upper = NULL, outcome = outcome,
    basis = "exact_tail_probability"
  ))
  policy$diagnostics$warnings <- if (identical(outcome, "triggered")) {
    "canonical explicit warning policy triggered"
  } else character()
  expect_invisible(.dpprior_validate_result_v1(policy))

  forged_policy <- policy
  forged_policy$diagnostics$policy_results[[1L]]$value <-
    min(1, value + 0.25)
  forged_policy$diagnostics$policy_results[[1L]]$outcome <-
    if (forged_policy$diagnostics$policy_results[[1L]]$value >
        action_threshold) "triggered" else "not_triggered"
  expect_error(.dpprior_validate_result_v1(forged_policy),
               class = "dpprior_schema_error")

  for (mutation in c("control", "evaluations", "iterations", "warnings",
                     "unavailable")) {
    forged_attempt <- .schema23_diagnostics()
    if (identical(mutation, "control")) {
      forged_attempt$computation$attempts[[1L]]$control$component <- "K"
    } else if (identical(mutation, "evaluations")) {
      forged_attempt$computation$attempts[[1L]]$evaluations$gradient_count <-
        999L
    } else if (identical(mutation, "iterations")) {
      forged_attempt$computation$attempts[[1L]]$iterations <- 999L
    } else if (identical(mutation, "warnings")) {
      forged_attempt$computation$attempts[[1L]]$warnings <- "fabricated"
    } else {
      forged_attempt$computation$attempts[[1L]]$unavailable[["start"]] <-
        "fabricated reason"
    }
    expect_error(.dpprior_validate_result_v1(forged_attempt),
                 class = "dpprior_schema_error")
  }
})


test_that("sensitivity route provenance and target authority are exact", {
  route_evidence <- list(
    direct = .schema23_sensitivity_route_evidence("direct_variance"),
    direct_method_implicit = .schema23_sensitivity_route_evidence(
      "direct_variance", method_explicit = FALSE
    ),
    confidence_low = .schema23_sensitivity_route_evidence(
      "qualitative_confidence", confidence_explicit = TRUE,
      confidence = "low"
    ),
    confidence_medium_explicit = .schema23_sensitivity_route_evidence(
      "qualitative_confidence", confidence_explicit = TRUE,
      confidence = "medium"
    ),
    confidence_medium_default = .schema23_sensitivity_route_evidence(
      "qualitative_confidence", confidence_explicit = FALSE,
      confidence = "medium"
    ),
    confidence_high = .schema23_sensitivity_route_evidence(
      "qualitative_confidence", confidence_explicit = TRUE,
      confidence = "high"
    ),
    cv = .schema23_sensitivity_route_evidence(
      "coefficient_of_variation"
    ),
    pmf_J = .schema23_sensitivity_route_evidence("strict_pmf"),
    pmf_J_plus_one = .schema23_sensitivity_route_evidence(
      "strict_pmf", structural_K0 = TRUE
    ),
    fallback = .schema23_sensitivity_route_evidence(
      "direct_variance", fallback = TRUE
    ),
    A1_diagnostics = .schema23_sensitivity_route_evidence(
      "qualitative_confidence", confidence_explicit = TRUE,
      A1_diagnostic_condition = TRUE
    )
  )
  route_objects <- lapply(route_evidence, .schema23_sensitivity_from_evidence)
  for (name in names(route_evidence)) {
    expect_type(
      .dpprior_sensitivity_validate_fit_evidence(
        route_evidence[[name]], paste0("scn_", name),
        paste0("fixture.", name)
      ),
      "list"
    )
    expect_invisible(.dpprior_validate_result_v1(route_objects[[name]]))
  }
  for (scenario in list(
    list(mu_K = 5, var_K = 8, method = "A2-KL"),
    list(mu_K = 5, confidence = "medium", method = "A2-KL")
  )) {
    producer <- .dpprior_run_elicitation_sensitivity(
      50L, list(scenario), M = 80L, check_diagnostics = FALSE
    )
    producer_raw <- unclass(producer)
    expect_identical(producer_raw$status, "converged")
    expect_identical(
      producer_raw$sensitivity$fit_evidence[[1L]]$method, "A2-KL"
    )
    expect_invisible(.dpprior_validate_result_v1(producer))
  }
  A1_diagnostic_failure <- route_objects$confidence_medium_default
  key <- A1_diagnostic_failure$sensitivity$scenarios$scenario_key[[1L]]
  diagnostic_condition <- .schema23_sensitivity_condition(
    "sensitivity_diagnostic_contract",
    "diagnostic adapter returned noncanonical evidence"
  )
  A1_diagnostic_failure$sensitivity$conditions[[key]]$diagnostics <-
    diagnostic_condition
  A1_diagnostic_failure$sensitivity$fit_evidence[[key]]$condition_evidence$
    diagnostics <- diagnostic_condition
  A1_diagnostic_failure$sensitivity$scenario_results$status <- "failed"
  A1_diagnostic_failure$sensitivity$scenario_results$usable <- FALSE
  A1_diagnostic_failure$sensitivity$scenario_results$verified <- FALSE
  A1_diagnostic_failure$sensitivity$metrics_long$status <- "failed"
  A1_diagnostic_failure$sensitivity$metrics_long$usable <- FALSE
  A1_diagnostic_failure$sensitivity$metrics_long$verified <- FALSE
  A1_diagnostic_failure$sensitivity$global$failed_count <- 1L
  A1_diagnostic_failure$status <- "failed"
  A1_diagnostic_failure$usable <- FALSE
  A1_diagnostic_failure$verified <- FALSE
  A1_diagnostic_failure$verification$passed <- FALSE
  A1_diagnostic_failure$verification$reason <-
    "diagnostic contract failure preserved"
  expect_invisible(.dpprior_validate_result_v1(A1_diagnostic_failure))
  interval_evidence <- lapply(
    c("equal_tail", "central_mass", "hard_bounds"),
    .schema23_sensitivity_interval_evidence
  )
  names(interval_evidence) <- c("equal_tail", "central_mass", "hard_bounds")
  for (name in names(interval_evidence)) {
    expect_type(
      .dpprior_sensitivity_validate_fit_evidence(
        interval_evidence[[name]], paste0("scn_interval_", name),
        paste0("fixture.interval.", name)
      ),
      "list"
    )
  }
  expect_identical(
    interval_evidence$hard_bounds$request$K_interval$coverage, 1
  )
  expect_identical(
    interval_evidence$hard_bounds$request$K_interval$mu_K, 5
  )
  expect_identical(
    interval_evidence$central_mass$request$K_interval$mu_K, 6.5
  )
  forged_verified_interval <- interval_evidence$equal_tail
  forged_verified_interval$status <- "converged"
  forged_verified_interval$usable <- TRUE
  forged_verified_interval$verified <- TRUE
  forged_verified_interval$condition_evidence <-
    .schema23_empty_condition_evidence()
  expect_error(
    .dpprior_sensitivity_validate_fit_evidence(
      forged_verified_interval, "scn_forged_verified_interval",
      "fixture.forged_verified_interval"
    ),
    class = "dpprior_schema_error"
  )

  explicit_medium <- route_objects$confidence_medium_explicit
  default_medium <- route_objects$confidence_medium_default
  expect_identical(
    explicit_medium$sensitivity$scenarios$canonical_content,
    default_medium$sensitivity$scenarios$canonical_content
  )
  expect_identical(
    explicit_medium$sensitivity$scenarios$scenario_key,
    default_medium$sensitivity$scenarios$scenario_key
  )
  expect_false(default_medium$sensitivity$scenarios$confidence_explicit)
  expect_true(explicit_medium$sensitivity$scenarios$confidence_explicit)
  expect_identical(
    route_objects$direct$sensitivity$scenarios$canonical_content,
    route_objects$direct_method_implicit$sensitivity$scenarios$
      canonical_content
  )
  expect_false(identical(
    route_objects$pmf_J$sensitivity$scenarios$scenario_key,
    route_objects$pmf_J_plus_one$sensitivity$scenarios$scenario_key
  ))

  confidence_vif <- c(low = 5, medium = 2.5, high = 1.5)
  for (level in names(confidence_vif)) {
    evidence <- route_evidence[[if (identical(level, "medium")) {
      "confidence_medium_explicit"
    } else {
      paste0("confidence_", level)
    }]]
    expect_equal(
      evidence$target$used$variance,
      confidence_vif[[level]] * (evidence$request$mu_K - 1),
      tolerance = 0
    )
  }
  expect_identical(route_evidence$cv$target$kind, "cv")
  expect_equal(
    route_evidence$cv$target$used$variance,
    (route_evidence$cv$request$mu_K * route_evidence$cv$request$cv_K)^2,
    tolerance = 0
  )
  expect_length(route_evidence$pmf_J$target$used$pmf, 20L)
  expect_length(route_evidence$pmf_J_plus_one$target$used$pmf, 20L)
  expect_identical(
    route_evidence$pmf_J$target$used$pmf,
    route_evidence$pmf_J_plus_one$target$used$pmf
  )
  expect_identical(route_evidence$fallback$method, "A2-MN+NM")

  direct <- route_objects$direct
  direct_key <- direct$sensitivity$scenarios$scenario_key[[1L]]
  provenance_mutations <- list(
    method_explicit = function(x) {
      x$sensitivity$fit_evidence[[direct_key]]$input_provenance$
        method_explicit <- FALSE
      x
    },
    confidence_explicit = function(x) {
      x$sensitivity$fit_evidence[[direct_key]]$input_provenance$
        confidence_explicit <- TRUE
      x
    },
    requested_method = function(x) {
      x$sensitivity$fit_evidence[[direct_key]]$input_provenance$
        requested_method <- "A1"
      x
    },
    selected_method = function(x) {
      x$sensitivity$fit_evidence[[direct_key]]$input_provenance$
        selected_method <- "A1"
      x
    },
    is_fallback = function(x) {
      x$sensitivity$fit_evidence[[direct_key]]$input_provenance$is_fallback <-
        TRUE
      x
    },
    target_route = function(x) {
      x$sensitivity$fit_evidence[[direct_key]]$input_provenance$target_route <-
        "coefficient_of_variation"
      x
    }
  )
  for (mutate in provenance_mutations) {
    expect_error(
      .dpprior_validate_result_v1(mutate(direct)),
      class = "dpprior_schema_error"
    )
  }
  reordered_provenance <- direct
  reordered_provenance$sensitivity$fit_evidence[[direct_key]]$
    input_provenance <- reordered_provenance$sensitivity$fit_evidence[[
      direct_key
    ]]$input_provenance[c(
      "confidence_explicit", "method_explicit", "requested_method",
      "selected_method", "is_fallback", "target_route"
    )]
  expect_error(.dpprior_validate_result_v1(reordered_provenance),
               class = "dpprior_schema_error")

  wrong_scenario_method <- direct
  wrong_scenario_method$sensitivity$scenarios$effective_method <- "A1"
  expect_error(.dpprior_validate_result_v1(wrong_scenario_method),
               class = "dpprior_schema_error")
  wrong_scenario_explicit <- direct
  wrong_scenario_explicit$sensitivity$scenarios$method_explicit <- FALSE
  expect_error(.dpprior_validate_result_v1(wrong_scenario_explicit),
               class = "dpprior_schema_error")
  coordinated_implicit_A1 <- route_objects$confidence_medium_default
  key <- coordinated_implicit_A1$sensitivity$scenarios$scenario_key[[1L]]
  coordinated_implicit_A1$sensitivity$fit_evidence[[key]]$input_provenance$
    method_explicit <- FALSE
  coordinated_implicit_A1$sensitivity$scenarios$method_explicit <- FALSE
  expect_error(.dpprior_validate_result_v1(coordinated_implicit_A1),
               class = "dpprior_schema_error")

  forged_confidence <- route_objects$confidence_low
  key <- forged_confidence$sensitivity$scenarios$scenario_key[[1L]]
  forged_confidence$sensitivity$fit_evidence[[key]]$target$used$variance <-
    forged_confidence$sensitivity$fit_evidence[[key]]$target$used$variance + 1
  expect_error(.dpprior_validate_result_v1(forged_confidence),
               class = "dpprior_schema_error")
  forged_cv <- route_objects$cv
  key <- forged_cv$sensitivity$scenarios$scenario_key[[1L]]
  forged_cv$sensitivity$fit_evidence[[key]]$target$used$variance <-
    forged_cv$sensitivity$fit_evidence[[key]]$target$used$variance + 1
  expect_error(.dpprior_validate_result_v1(forged_cv),
               class = "dpprior_schema_error")

  forged_K0 <- route_objects$pmf_J_plus_one
  key <- forged_K0$sensitivity$scenarios$scenario_key[[1L]]
  forged_K0$sensitivity$fit_evidence[[key]]$request$target_pmf[[1L]] <- 1e-12
  expect_error(.dpprior_validate_result_v1(forged_K0),
               class = "dpprior_schema_error")
  retained_K0 <- route_objects$pmf_J_plus_one
  key <- retained_K0$sensitivity$scenarios$scenario_key[[1L]]
  retained_K0$sensitivity$fit_evidence[[key]]$target$used$pmf <-
    retained_K0$sensitivity$fit_evidence[[key]]$request$target_pmf
  expect_error(.dpprior_validate_result_v1(retained_K0),
               class = "dpprior_schema_error")
  forged_pmf_targets <- list(
    uniform = rep(1 / 20, 20L),
    point_mass = c(1, rep(0, 19L))
  )
  for (forged_pmf in forged_pmf_targets) {
    forged_evidence <- route_evidence$pmf_J
    forged_evidence$request$target_pmf <- forged_pmf
    forged_target <- .dpprior_sensitivity_canonical_target(
      forged_evidence$request, "strict_pmf", "fixture.forged_pmf"
    )
    forged_target_raw <- unclass(forged_target)
    forged_evidence$target <- list(
      kind = forged_target_raw$kind, J = forged_target_raw$J,
      request = forged_target_raw$request, used = forged_target_raw$used
    )
    expect_error(
      .dpprior_sensitivity_validate_fit_evidence(
        forged_evidence, "scn_forged_pmf", "fixture.forged_pmf"
      ),
      class = "dpprior_schema_error"
    )
    expect_error(
      .dpprior_validate_result_v1(.schema23_sensitivity_from_evidence(
        forged_evidence, validate = FALSE
      )),
      class = "dpprior_schema_error"
    )
  }

  false_A1_status <- route_objects$confidence_medium_default
  key <- false_A1_status$sensitivity$scenarios$scenario_key[[1L]]
  false_A1_status$sensitivity$fit_evidence[[key]]$status <- "converged"
  false_A1_status$sensitivity$fit_evidence[[key]]$verified <- TRUE
  expect_error(.dpprior_validate_result_v1(false_A1_status),
               class = "dpprior_schema_error")
  false_A2_usable <- .schema23_sensitivity_interval()
  key <- false_A2_usable$sensitivity$scenarios$scenario_key[[1L]]
  false_A2_usable$sensitivity$fit_evidence[[key]]$usable <- TRUE
  expect_error(.dpprior_validate_result_v1(false_A2_usable),
               class = "dpprior_schema_error")
  missing_A2_condition <- .schema23_sensitivity_interval()
  key <- missing_A2_condition$sensitivity$scenarios$scenario_key[[1L]]
  missing_A2_condition$sensitivity$conditions[[key]]["calibration"] <-
    list(NULL)
  missing_A2_condition$sensitivity$fit_evidence[[key]]$condition_evidence[
    "calibration"
  ] <- list(NULL)
  expect_error(.dpprior_validate_result_v1(missing_A2_condition),
               class = "dpprior_schema_error")
  misplaced_A1_condition <- route_objects$A1_diagnostics
  key <- misplaced_A1_condition$sensitivity$scenarios$scenario_key[[1L]]
  condition <- misplaced_A1_condition$sensitivity$conditions[[key]]$diagnostics
  misplaced_A1_condition$sensitivity$conditions[[key]][
    c("calibration", "diagnostics")
  ] <- list(condition, NULL)
  misplaced_A1_condition$sensitivity$fit_evidence[[key]]$condition_evidence[
    c("calibration", "diagnostics")
  ] <- list(condition, NULL)
  expect_error(.dpprior_validate_result_v1(misplaced_A1_condition),
               class = "dpprior_schema_error")

  dispatch_count <- 0L
  hostile_accessor <- function(...) {
    dispatch_count <<- dispatch_count + 1L
    stop("hostile sensitivity request accessor dispatched")
  }
  assign("[[.schema23_hostile_sensitivity_request", hostile_accessor,
         envir = .GlobalEnv)
  assign("$.schema23_hostile_sensitivity_request", hostile_accessor,
         envir = .GlobalEnv)
  assign("[[.schema23_hostile_sensitivity_evidence", hostile_accessor,
         envir = .GlobalEnv)
  assign("$.schema23_hostile_sensitivity_evidence", hostile_accessor,
         envir = .GlobalEnv)
  assign("[<-.schema23_hostile_sensitivity_evidence", hostile_accessor,
         envir = .GlobalEnv)
  assign("[[.schema23_hostile_sensitivity_condition", hostile_accessor,
         envir = .GlobalEnv)
  assign("$.schema23_hostile_sensitivity_condition", hostile_accessor,
         envir = .GlobalEnv)
  assign("[<-.schema23_hostile_sensitivity_condition", hostile_accessor,
         envir = .GlobalEnv)
  on.exit({
    rm("[[.schema23_hostile_sensitivity_request", envir = .GlobalEnv)
    rm("$.schema23_hostile_sensitivity_request", envir = .GlobalEnv)
    rm("[[.schema23_hostile_sensitivity_evidence", envir = .GlobalEnv)
    rm("$.schema23_hostile_sensitivity_evidence", envir = .GlobalEnv)
    rm("[<-.schema23_hostile_sensitivity_evidence", envir = .GlobalEnv)
    rm("[[.schema23_hostile_sensitivity_condition", envir = .GlobalEnv)
    rm("$.schema23_hostile_sensitivity_condition", envir = .GlobalEnv)
    rm("[<-.schema23_hostile_sensitivity_condition", envir = .GlobalEnv)
  }, add = TRUE)
  hostile_request <- direct
  class(hostile_request$sensitivity$fit_evidence[[direct_key]]$request) <-
    "schema23_hostile_sensitivity_request"
  expect_error(.dpprior_validate_result_v1(hostile_request),
               class = "dpprior_schema_error")
  expect_identical(dispatch_count, 0L)
  hostile_evidence <- direct
  class(hostile_evidence$sensitivity$fit_evidence[[direct_key]]) <-
    "schema23_hostile_sensitivity_evidence"
  expect_error(.dpprior_validate_result_v1(hostile_evidence),
               class = "dpprior_schema_error")
  expect_identical(dispatch_count, 0L)
  hostile_condition <- direct
  class(hostile_condition$sensitivity$conditions[[direct_key]]) <-
    "schema23_hostile_sensitivity_condition"
  expect_error(.dpprior_validate_result_v1(hostile_condition),
               class = "dpprior_schema_error")
  expect_identical(dispatch_count, 0L)

  for (object in route_objects) {
    restored <- unserialize(serialize(object, NULL, version = 3L))
    expect_identical(restored, object)
    expect_invisible(.dpprior_validate_result_v1(restored))
  }
})


test_that("fit-attached diagnostics are bound to fresh canonical authority", {
  source_fits <- list(
    a1 = DPprior_a1(20L, 5, 8),
    a2_moment = DPprior_a2_newton(20L, 5, 8, M = 80L),
    a2_kl = DPprior_a2_kl(
      20L, list(mu_K = 4, var_K = 5), method = "chisq", M = 40L
    )
  )
  attached <- lapply(source_fits, .schema23_attach_fit_diagnostics)
  for (fit in attached) {
    expect_identical(
      names(fit$diagnostics),
      c(
        "authority", "policy_results", "warnings", "alpha", "K",
        "weights", "coclustering"
      )
    )
    expect_invisible(.dpprior_validate_result_v1(fit))
    restored <- unserialize(serialize(fit, NULL, version = 3L))
    expect_identical(restored, fit)
    expect_invisible(.dpprior_validate_result_v1(restored))
  }
  expect_false(
    "diagnostics" %in%
      names(attached$a2_moment$compatibility$top_level_aliases)
  )
  expect_true(
    "legacy_v2" %in% names(attached$a2_moment$compatibility$views)
  )

  WSB_policy <- list(
    estimand = "W_SB", direction = "above", weight_threshold = 0.5,
    action_threshold = 0.2
  )
  WSB <- .schema23_attach_fit_diagnostics(
    source_fits$a2_moment, warning_policy = WSB_policy
  )
  expect_identical(WSB$diagnostics$policy_results[[1L]]$basis,
                   "exact_tail_probability")
  expect_identical(length(WSB$diagnostics$warnings), 1L)
  expect_invisible(.dpprior_validate_result_v1(WSB))

  Wmax <- .schema23_attach_fit_diagnostics(
    source_fits$a2_kl,
    warning_policy = list(
      estimand = "W_max", direction = "above", weight_threshold = 0.5,
      action_threshold = 0.2
    )
  )
  expect_identical(Wmax$diagnostics$policy_results[[1L]]$basis,
                   "backend_unavailable")
  expect_identical(Wmax$diagnostics$policy_results[[1L]]$outcome,
                   "indeterminate")
  expect_length(Wmax$diagnostics$warnings, 0L)
  expect_invisible(.dpprior_validate_result_v1(Wmax))

  approximate_closed <- .schema23_attach_fit_diagnostics(
    source_fits$a1, M_selected = 10L, allow_approximate = FALSE
  )
  expect_false(all(vapply(
    approximate_closed$diagnostics[.DPPRIOR_DIAGNOSTIC_COMPONENTS],
    `[[`, logical(1), "usable"
  )))
  expect_invisible(.dpprior_validate_result_v1(approximate_closed))

  approximate <- .schema23_attach_fit_diagnostics(
    source_fits$a1, M_selected = 10L, allow_approximate = TRUE
  )
  expect_identical(
    unname(vapply(
      approximate$diagnostics[.DPPRIOR_DIAGNOSTIC_COMPONENTS],
      `[[`, character(1), "status"
    )),
    c("converged", "approximate", "approximate", "approximate")
  )
  expect_true(all(vapply(
    approximate$diagnostics[.DPPRIOR_DIAGNOSTIC_COMPONENTS],
    `[[`, logical(1), "usable"
  )))
  expect_identical(approximate$status, source_fits$a1$status)
  expect_identical(approximate$usable, source_fits$a1$usable)
  expect_identical(approximate$verified, source_fits$a1$verified)
  expect_invisible(.dpprior_validate_result_v1(approximate))

  standalone <- .schema23_diagnostics()
  expect_invisible(.dpprior_validate_result_v1(standalone))
  standalone$diagnostics <- c(
    list(authority = attached$a2_moment$diagnostics$authority),
    standalone$diagnostics
  )
  expect_error(.dpprior_validate_result_v1(standalone),
               class = "dpprior_schema_error")

  missing_authority <- attached$a2_moment
  missing_authority$diagnostics["authority"] <- list(NULL)
  expect_error(.dpprior_validate_result_v1(missing_authority),
               class = "dpprior_schema_error")

  numeric_forgers <- list(
    alpha = function(x) {
      x$diagnostics$alpha$mean <- 9
      x
    },
    K = function(x) {
      x$diagnostics$K$mean <- 9
      x
    },
    weights = function(x) {
      x$diagnostics$weights$mean <- 0.9
      x
    },
    coclustering = function(x) {
      x$diagnostics$coclustering$mean <- 0.9
      x
    }
  )
  for (forge in numeric_forgers) {
    expect_error(.dpprior_validate_result_v1(forge(attached$a2_moment)),
                 class = "dpprior_schema_error")
  }
  coordinated_pmf <- attached$a2_moment
  coordinated_pmf$diagnostics$K$pmf <-
    rev(coordinated_pmf$diagnostics$K$pmf)
  coordinated_moments <- .dpprior_target_pmf_moments(
    coordinated_pmf$diagnostics$K$pmf
  )
  coordinated_pmf$diagnostics$K$mean <-
    unname(coordinated_moments[["mean"]])
  coordinated_pmf$diagnostics$K$variance <-
    unname(coordinated_moments[["variance"]])
  expect_error(.dpprior_validate_result_v1(coordinated_pmf),
               class = "dpprior_schema_error")

  authority_forgers <- list(
    method = function(x) {
      x$diagnostics$authority$method <- "optimizer_self_report"
      x
    },
    M_selected = function(x) {
      x$diagnostics$authority$M_selected <- 81L
      x
    },
    M_required = function(x) {
      x$diagnostics$authority$M_verification_required <- 161L
      x
    },
    M_used = function(x) {
      x$diagnostics$authority$M_verification_used <- 161L
      x
    },
    absolute = function(x) {
      x$diagnostics$authority$absolute_tolerance <- 1e-3
      x
    },
    relative = function(x) {
      x$diagnostics$authority$relative_tolerance <- 1e-3
      x
    },
    pmf_mass = function(x) {
      x$diagnostics$authority$pmf_mass_tolerance <- 1e-3
      x
    },
    K_M = function(x) {
      x$diagnostics$K$M <- 160L
      x
    },
    K_status = function(x) {
      x$diagnostics$K$status <- "approximate"
      x$diagnostics$K$usable <- FALSE
      x$diagnostics$K$verified <- FALSE
      x
    }
  )
  for (forge in authority_forgers) {
    expect_error(.dpprior_validate_result_v1(forge(attached$a2_moment)),
                 class = "dpprior_schema_error")
  }

  wrong_A1_order <- attached$a1
  wrong_A1_order$diagnostics$authority$M_verification_required <- 161L
  wrong_A1_order$diagnostics$authority$M_verification_used <- 161L
  expect_error(.dpprior_validate_result_v1(wrong_A1_order),
               class = "dpprior_schema_error")

  approximate_without_opt_in <- approximate
  approximate_without_opt_in$diagnostics$authority$allow_approximate <- FALSE
  expect_error(.dpprior_validate_result_v1(approximate_without_opt_in),
               class = "dpprior_schema_error")

  forged_WSB <- WSB
  forged_WSB$diagnostics$policy_results[[1L]]$value <- 0.99
  forged_WSB$diagnostics$policy_results[[1L]]$outcome <- "triggered"
  expect_error(.dpprior_validate_result_v1(forged_WSB),
               class = "dpprior_schema_error")
  wrong_WSB_warning_count <- WSB
  wrong_WSB_warning_count$diagnostics$warnings <- character()
  expect_error(.dpprior_validate_result_v1(wrong_WSB_warning_count),
               class = "dpprior_schema_error")
  wrong_policy_authority <- WSB
  wrong_policy_authority$diagnostics$authority$warning_policy$weight_threshold <-
    0.6
  expect_error(.dpprior_validate_result_v1(wrong_policy_authority),
               class = "dpprior_schema_error")
  invalid_policy_vocabulary <- WSB
  invalid_policy_vocabulary$diagnostics$authority$warning_policy$estimand <-
    "largest_weight"
  expect_error(.dpprior_validate_result_v1(invalid_policy_vocabulary),
               class = "dpprior_schema_error")

  Wmax_escalation <- Wmax
  Wmax_escalation$diagnostics$policy_results[[1L]]$basis <-
    "certified_bounds"
  Wmax_escalation$diagnostics$policy_results[[1L]]$lower <- 0.1
  Wmax_escalation$diagnostics$policy_results[[1L]]$upper <- 0.9
  expect_error(.dpprior_validate_result_v1(Wmax_escalation),
               class = "dpprior_schema_error")

  parameterless_raw <- unclass(attached$a2_moment)
  parameterless_raw["parameters"] <- list(NULL)
  expect_error(
    .dpprior_validate_diagnostics_extension(
      parameterless_raw$diagnostics, fit_raw = parameterless_raw
    ),
    class = "dpprior_schema_error"
  )
  parameterless_fit <- unclass(.schema23_migrated_A2(
    "a2_moment", verify = FALSE
  ))
  parameterless_fit$diagnostics <- attached$a2_moment$diagnostics
  class(parameterless_fit) <- class(attached$a2_moment)
  expect_error(.dpprior_validate_result_v1(parameterless_fit),
               class = "dpprior_schema_error")

  attributed_authority <- attached$a2_moment
  attr(attributed_authority$diagnostics$authority, "payload") <- new.env()
  expect_error(.dpprior_validate_result_v1(attributed_authority),
               class = "dpprior_schema_error")
  attributed_policy <- WSB
  attr(
    attributed_policy$diagnostics$authority$warning_policy, "payload"
  ) <- new.env()
  expect_error(.dpprior_validate_result_v1(attributed_policy),
               class = "dpprior_schema_error")
})


test_that("sensitivity content, intervals, metrics, and failure rows reconcile", {
  expect_invisible(.dpprior_validate_result_v1(.schema23_sensitivity_failed()))
  expect_invisible(
    .dpprior_validate_result_v1(.schema23_sensitivity_infeasible())
  )
  expect_invisible(.dpprior_validate_result_v1(.schema23_sensitivity_interval()))
  expect_invisible(.dpprior_validate_result_v1(
    .schema23_sensitivity_interval_infeasible()
  ))
  expect_invisible(.dpprior_validate_result_v1(
    .schema23_sensitivity_interval_infeasible("mean_hull")
  ))

  forged_infeasibility <- .schema23_sensitivity_infeasible()
  key <- forged_infeasibility$sensitivity$scenarios$scenario_key[[1L]]
  forged_infeasibility$sensitivity$conditions[[key]]$calibration$code <-
    "fabricated_certificate"
  expect_error(.dpprior_validate_result_v1(forged_infeasibility),
               class = "dpprior_schema_error")

  forged_interval_certificate <- .schema23_sensitivity_interval_infeasible()
  key <- forged_interval_certificate$sensitivity$scenarios$scenario_key[[1L]]
  forged_interval_certificate$sensitivity$conditions[[key]]$target$message <-
    "fabricated target certificate"
  expect_error(.dpprior_validate_result_v1(forged_interval_certificate),
               class = "dpprior_schema_error")
  missing_interval_chain <- .schema23_sensitivity_interval_infeasible()
  key <- missing_interval_chain$sensitivity$scenarios$scenario_key[[1L]]
  missing_interval_chain$sensitivity$conditions[[key]]["target"] <- list(NULL)
  missing_interval_chain$sensitivity$fit_evidence[[key]]$condition_evidence[
    "target"
  ] <- list(NULL)
  expect_error(.dpprior_validate_result_v1(missing_interval_chain),
               class = "dpprior_schema_error")
  duplicated_interval_chain <- .schema23_sensitivity_interval_infeasible()
  key <- duplicated_interval_chain$sensitivity$scenarios$scenario_key[[1L]]
  duplicated_interval_chain$sensitivity$conditions[[key]]$interval <-
    duplicated_interval_chain$sensitivity$conditions[[key]]$target
  duplicated_interval_chain$sensitivity$fit_evidence[[key]]$
    condition_evidence$interval <-
    duplicated_interval_chain$sensitivity$conditions[[key]]$target
  expect_error(.dpprior_validate_result_v1(duplicated_interval_chain),
               class = "dpprior_schema_error")
  terminal_target_with_parameters <-
    .schema23_sensitivity_interval_infeasible()
  key <- terminal_target_with_parameters$sensitivity$scenarios$
    scenario_key[[1L]]
  numerical <- .schema23_sensitivity_fit_evidence()
  terminal_target_with_parameters$sensitivity$fit_evidence[[key]][c(
    "parameters", "selected_snapshot", "verifier_snapshot"
  )] <- numerical[c("parameters", "selected_snapshot", "verifier_snapshot")]
  terminal_target_with_parameters$sensitivity$fit_evidence[[key]]$status <-
    "converged"
  terminal_target_with_parameters$sensitivity$fit_evidence[[key]]$usable <- TRUE
  terminal_target_with_parameters$sensitivity$fit_evidence[[key]]$verified <-
    TRUE
  expect_error(.dpprior_validate_result_v1(terminal_target_with_parameters),
               class = "dpprior_schema_error")

  J1_available <- .schema23_sensitivity_infeasible()
  key <- J1_available$sensitivity$scenarios$scenario_key[[1L]]
  J1_evidence <- J1_available$sensitivity$fit_evidence[[key]]
  numerical <- .schema23_sensitivity_fit_evidence(
    request = J1_evidence$request, status = "converged", usable = TRUE,
    verified = TRUE,
    input_provenance = J1_evidence$input_provenance
  )
  expect_error(
    .dpprior_sensitivity_validate_fit_evidence(
      numerical, key, "fixture.J1_available"
    ),
    class = "dpprior_schema_error"
  )
  J1_interval_target <- .dp_target_K(
    J = 1L,
    K_interval = list(
      lower = 1L, upper = 1L, type = "hard_bounds", family = "maxent"
    )
  )
  J1_interval <- unclass(J1_interval_target)$used$interval
  J1_interval_request <- list(
    J = 1L, K_interval = J1_interval, mu_K = J1_interval$mu_K,
    method = "A2-KL", M = 80L
  )
  J1_interval_evidence <- .schema23_sensitivity_fit_evidence(
    request = J1_interval_request, status = "converged", usable = TRUE,
    verified = TRUE,
    input_provenance = list(
      method_explicit = TRUE, confidence_explicit = FALSE,
      requested_method = "A2-KL", selected_method = "A2-KL",
      is_fallback = FALSE, target_route = "interval"
    )
  )
  expect_error(
    .dpprior_sensitivity_validate_fit_evidence(
      J1_interval_evidence, "scn_J1_interval", "fixture.J1_interval"
    ),
    class = "dpprior_schema_error"
  )

  wrong_default_J <- .schema23_sensitivity()
  wrong_default_J$target$defaults$J <- 999L
  expect_error(.dpprior_validate_result_v1(wrong_default_J),
               class = "dpprior_schema_error")

  wrong_diagnostic_flag <- .schema23_sensitivity()
  wrong_diagnostic_flag$sensitivity$scenarios$diagnostics_requested <- FALSE
  expect_error(.dpprior_validate_result_v1(wrong_diagnostic_flag),
               class = "dpprior_schema_error")

  wrong_effective_method <- .schema23_sensitivity()
  wrong_effective_method$sensitivity$scenarios$effective_method <- "A1"
  expect_error(.dpprior_validate_result_v1(wrong_effective_method),
               class = "dpprior_schema_error")

  fabricated_label <- .schema23_sensitivity()
  fabricated_label$sensitivity$scenarios$scenario_label <- "fabricated label"
  expect_error(.dpprior_validate_result_v1(fabricated_label),
               class = "dpprior_schema_error")

  condition_without_evidence <- .schema23_sensitivity_failed()
  key <- condition_without_evidence$sensitivity$scenarios$scenario_key[[1L]]
  condition_without_evidence$sensitivity$conditions[[key]]$calibration$code <-
    "fabricated_condition"
  expect_error(.dpprior_validate_result_v1(condition_without_evidence),
               class = "dpprior_schema_error")

  interval_condition_on_direct <- .schema23_sensitivity_failed()
  key <- interval_condition_on_direct$sensitivity$scenarios$scenario_key[[1L]]
  interval_solver_condition <- list(
    class = "dpprior_interval_solver_failed",
    classes = c(
      "dpprior_interval_solver_failed", "dpprior_numerical_error",
      "dpprior_calibration_error", "dpprior_error", "error",
      "dpprior_condition", "condition"
    ),
    code = "target_failed", message = "fabricated interval solver failure"
  )
  interval_condition_on_direct$sensitivity$conditions[[key]]$calibration <-
    interval_solver_condition
  interval_condition_on_direct$sensitivity$fit_evidence[[key]]$
    condition_evidence$calibration <- interval_solver_condition
  expect_error(.dpprior_validate_result_v1(interval_condition_on_direct),
               class = "dpprior_schema_error")

  retained_warning <- .schema23_sensitivity_failed()
  key <- retained_warning$sensitivity$scenarios$scenario_key[[1L]]
  warning_record <- list(
    class = "simpleWarning", classes = c("simpleWarning", "warning", "condition"),
    code = "retained_backend_warning", message = "retained backend warning"
  )
  retained_warning$sensitivity$conditions[[key]]$calibration_warnings <-
    list(warning_record)
  retained_warning$sensitivity$conditions[[key]]$diagnostics <-
    .schema23_sensitivity_condition(
      "sensitivity_diagnostic_contract",
      "warnings accompanied sensitivity calibration; claims quarantined"
    )
  retained_warning$sensitivity$fit_evidence[[key]]$condition_evidence$
    calibration_warnings <- list(warning_record)
  retained_warning$sensitivity$fit_evidence[[key]]$condition_evidence$
    diagnostics <- retained_warning$sensitivity$conditions[[key]]$diagnostics
  expect_invisible(.dpprior_validate_result_v1(retained_warning))

  misplaced_warning <- retained_warning
  misplaced_warning$sensitivity$conditions[[key]]$calibration_warnings <-
    list()
  misplaced_warning$sensitivity$conditions[[key]]$diagnostic_warnings <-
    list(warning_record)
  misplaced_warning$sensitivity$fit_evidence[[key]]$condition_evidence$
    calibration_warnings <- list()
  misplaced_warning$sensitivity$fit_evidence[[key]]$condition_evidence$
    diagnostic_warnings <- list(warning_record)
  expect_error(.dpprior_validate_result_v1(misplaced_warning),
               class = "dpprior_schema_error")

  warning_without_quarantine <- retained_warning
  warning_without_quarantine$sensitivity$conditions[[key]]$diagnostics <- NULL
  warning_without_quarantine$sensitivity$fit_evidence[[key]]$
    condition_evidence$diagnostics <- NULL
  expect_error(.dpprior_validate_result_v1(warning_without_quarantine),
               class = "dpprior_schema_error")

  A2_condition <- tryCatch(
    DPprior_fit(
      20L, 4, 8, method = "A2-KL", M = 80L,
      check_diagnostics = FALSE, verbose = FALSE
    ),
    error = identity
  )
  expect_s3_class(A2_condition, "dpprior_calibration_unusable")
  A2_warning <- .dpprior_run_elicitation_sensitivity(
    20L,
    list(list(mu_K = 4, var_K = 8, method = "A2-KL", M = 80L)),
    check_diagnostics = FALSE,
    .fit_fun = function(...) {
      warning("adapter warning")
      stop(A2_condition)
    }
  )
  A2_warning_raw <- unclass(A2_warning)
  A2_warning_evidence <- A2_warning_raw$sensitivity$fit_evidence[[1L]]
  expect_identical(A2_warning_raw$status, "failed")
  expect_identical(A2_warning_evidence$status, "approximate")
  expect_identical(
    A2_warning_evidence$condition_evidence$calibration$code,
    "calibration_unusable"
  )
  expect_identical(
    A2_warning_evidence$condition_evidence$diagnostics$code,
    "sensitivity_diagnostic_contract"
  )
  expect_length(
    A2_warning_evidence$condition_evidence$calibration_warnings, 1L
  )
  expect_invisible(.dpprior_validate_result_v1(A2_warning))

  fabricated_warning <- .schema23_sensitivity_failed()
  key <- fabricated_warning$sensitivity$scenarios$scenario_key[[1L]]
  warning_record <- list(
    class = "dpprior_calibration_error",
    classes = c("dpprior_calibration_error", "error", "condition"),
    code = "fabricated_warning", message = "fabricated warning"
  )
  fabricated_warning$sensitivity$conditions[[key]]$calibration_warnings <-
    list(warning_record)
  fabricated_warning$sensitivity$fit_evidence[[key]]$condition_evidence$
    calibration_warnings <- list(warning_record)
  expect_error(.dpprior_validate_result_v1(fabricated_warning),
               class = "dpprior_schema_error")

  fake_snapshot_parameters <- .schema23_sensitivity()
  parameters <- .dpprior_new_parameters(2, 3, "log_ab")
  fake_snapshot_parameters$verification$selected_snapshot$parameters <-
    parameters
  fake_snapshot_parameters$verification$verifier_snapshot$parameters <-
    parameters
  expect_error(.dpprior_validate_result_v1(fake_snapshot_parameters),
               class = "dpprior_schema_error")

  fake_snapshot_M <- .schema23_sensitivity()
  fake_snapshot_M$verification$selected_snapshot$M <- 1L
  expect_error(.dpprior_validate_result_v1(fake_snapshot_M),
               class = "dpprior_schema_error")

  fake_public_parameters <- .schema23_sensitivity()
  fake_public_parameters$parameters <- parameters
  fake_public_parameters$provenance$parameterization <- "log_ab"
  expect_error(.dpprior_validate_result_v1(fake_public_parameters),
               class = "dpprior_schema_error")

  for (order_name in c(
    "M_requested", "M_selected", "M_verification_required",
    "M_verification_used"
  )) {
    fake_order <- .schema23_sensitivity()
    fake_order$computation$orders[[order_name]] <- 1L
    expect_error(.dpprior_validate_result_v1(fake_order),
                 class = "dpprior_schema_error")
  }

  impossible_probability <- .schema23_sensitivity()
  probability_row <- match(
    "P_W_SB_gt_50", impossible_probability$sensitivity$metrics_long$metric
  )
  impossible_probability$sensitivity$scenario_results$P_W_SB_gt_50 <- 2
  impossible_probability$sensitivity$metrics_long$value[[probability_row]] <- 2
  expect_error(.dpprior_validate_result_v1(impossible_probability),
               class = "dpprior_schema_error")

  wrong_probability_formula <- .schema23_sensitivity()
  probability_row <- match(
    "P_W_SB_gt_50", wrong_probability_formula$sensitivity$metrics_long$metric
  )
  wrong_probability_formula$sensitivity$scenario_results$P_W_SB_gt_50 <- 0.2
  wrong_probability_formula$sensitivity$metrics_long$value[[probability_row]] <-
    0.2
  expect_error(.dpprior_validate_result_v1(wrong_probability_formula),
               class = "dpprior_schema_error")

  for (metric in c(
    "E_K_J", "Var_K_J", "E_W_SB", "E_rho",
    "P_W_max_gt_50_lower_bound"
  )) {
    coordinated_metric <- .schema23_sensitivity()
    current <- coordinated_metric$sensitivity$scenario_results[[metric]][[1L]]
    forged <- if (metric %in% c("E_K_J", "Var_K_J")) {
      current + 0.25
    } else {
      current / 2
    }
    coordinated_metric$sensitivity$scenario_results[[metric]][[1L]] <- forged
    row <- match(metric, coordinated_metric$sensitivity$metrics_long$metric)
    coordinated_metric$sensitivity$metrics_long$value[[row]] <- forged
    expect_error(.dpprior_validate_result_v1(coordinated_metric),
                 class = "dpprior_schema_error")
  }

  forged_fit_snapshot <- .schema23_sensitivity()
  key <- forged_fit_snapshot$sensitivity$scenarios$scenario_key[[1L]]
  forged_fit_snapshot$sensitivity$fit_evidence[[key]]$selected_snapshot$K$
    mean <- forged_fit_snapshot$sensitivity$fit_evidence[[key]]$
      selected_snapshot$K$mean + 0.25
  expect_error(.dpprior_validate_result_v1(forged_fit_snapshot),
               class = "dpprior_schema_error")

  impossible_partial_K <- .schema23_sensitivity_failed()
  impossible_partial_K$sensitivity$scenario_results$E_K_J <- 999
  impossible_partial_K$sensitivity$scenario_results$Var_K_J <- NA_real_
  impossible_partial_K$sensitivity$scenario_results$CV_K_J <- NA_real_
  for (metric in c("E_K_J", "Var_K_J", "CV_K_J")) {
    row <- match(metric, impossible_partial_K$sensitivity$metrics_long$metric)
    impossible_partial_K$sensitivity$metrics_long$value[[row]] <-
      impossible_partial_K$sensitivity$scenario_results[[metric]]
    impossible_partial_K$sensitivity$metrics_long$reason[[row]] <-
      if (identical(metric, "E_K_J")) NA_character_ else "unavailable"
  }
  expect_error(.dpprior_validate_result_v1(impossible_partial_K),
               class = "dpprior_schema_error")

  impossible_partial_variance <- .schema23_sensitivity_failed()
  impossible_partial_variance$sensitivity$scenario_results$E_K_J <- NA_real_
  impossible_partial_variance$sensitivity$scenario_results$Var_K_J <- 999
  impossible_partial_variance$sensitivity$scenario_results$CV_K_J <- NA_real_
  for (metric in c("E_K_J", "Var_K_J", "CV_K_J")) {
    row <- match(
      metric, impossible_partial_variance$sensitivity$metrics_long$metric
    )
    impossible_partial_variance$sensitivity$metrics_long$value[[row]] <-
      impossible_partial_variance$sensitivity$scenario_results[[metric]]
    impossible_partial_variance$sensitivity$metrics_long$reason[[row]] <-
      if (identical(metric, "Var_K_J")) NA_character_ else "unavailable"
  }
  expect_error(.dpprior_validate_result_v1(impossible_partial_variance),
               class = "dpprior_schema_error")

  unpaired_Wmax_bound <- .schema23_sensitivity_failed()
  unpaired_Wmax_bound$sensitivity$scenario_results$
    P_W_max_gt_50_lower_bound <- 0.2
  row <- match(
    "P_W_max_gt_50_lower_bound",
    unpaired_Wmax_bound$sensitivity$metrics_long$metric
  )
  unpaired_Wmax_bound$sensitivity$metrics_long$value[[row]] <- 0.2
  unpaired_Wmax_bound$sensitivity$metrics_long$reason[[row]] <- NA_character_
  expect_error(.dpprior_validate_result_v1(unpaired_Wmax_bound),
               class = "dpprior_schema_error")

  disconnected_interval <- .schema23_sensitivity()
  key <- disconnected_interval$sensitivity$scenarios$scenario_key[[1L]]
  disconnected_interval$sensitivity$interval_checks[[key]] <- list(
    requested = list(
      lower = 3L, upper = 10L, type = "equal_tail", coverage = 0.8,
      family = "maxent"
    ),
    selected = list(
      coverage = 0.79, lower_tail = 0.1, upper_tail = 0.11,
      coverage_residual = -0.01
    ),
    verification = NULL, status = "approximate",
    source = "selected_diagnostic_pmf_recalculation",
    usable = FALSE, verified = FALSE,
    reason = "diagnostic recalculation is not decision-ready"
  )
  expect_error(.dpprior_validate_result_v1(disconnected_interval),
               class = "dpprior_schema_error")

  opaque_interval_verifier <- .schema23_sensitivity_interval()
  key <- opaque_interval_verifier$sensitivity$scenarios$scenario_key[[1L]]
  opaque_interval_verifier$sensitivity$interval_checks[[key]]$verification <-
    list(forged = "opaque")
  expect_error(.dpprior_validate_result_v1(opaque_interval_verifier),
               class = "dpprior_schema_error")

  interval_authority_mutations <- list(
    top_mu = function(x, key) {
      x$sensitivity$fit_evidence[[key]]$request$mu_K <- 7
      x
    },
    nested_mu = function(x, key) {
      x$sensitivity$fit_evidence[[key]]$request$K_interval$mu_K <- 7
      x
    },
    support = function(x, key) {
      x$sensitivity$fit_evidence[[key]]$request$K_interval$support <-
        c(lower = 0L, upper = 20L)
      x
    },
    endpoints = function(x, key) {
      x$sensitivity$fit_evidence[[key]]$request$K_interval$endpoints <-
        "exclusive"
      x
    },
    family = function(x, key) {
      x$sensitivity$fit_evidence[[key]]$request$K_interval$family <-
        "fabricated"
      x
    }
  )
  for (mutate in interval_authority_mutations) {
    forged <- .schema23_sensitivity_interval()
    key <- forged$sensitivity$scenarios$scenario_key[[1L]]
    expect_error(.dpprior_validate_result_v1(mutate(forged, key)),
                 class = "dpprior_schema_error")
  }

  failed_backcheck <- .schema23_sensitivity_interval()
  key <- failed_backcheck$sensitivity$scenarios$scenario_key[[1L]]
  failed_backcheck$sensitivity$interval_checks[[key]]$status <- "failed"
  expect_error(.dpprior_validate_result_v1(failed_backcheck),
               class = "dpprior_schema_error")
  for (forged_tolerance in c(0, 1e-12, 1e-9)) {
    forged_backcheck_tolerance <- .schema23_sensitivity_interval()
    key <- forged_backcheck_tolerance$sensitivity$scenarios$
      scenario_key[[1L]]
    forged_backcheck_tolerance$sensitivity$interval_checks[[key]]$
      verification$tolerance <- forged_tolerance
    expect_error(
      .dpprior_validate_result_v1(forged_backcheck_tolerance),
      class = "dpprior_schema_error"
    )
  }
  passing_interval <- .dpprior_run_elicitation_sensitivity(
    20L,
    list(list(
      K_interval = list(
        lower = 1L, upper = 20L, type = "hard_bounds",
        family = "maxent"
      ),
      method = "A2-KL", M = 80L
    )),
    check_diagnostics = FALSE
  )
  passing_interval_raw <- unclass(passing_interval)
  expect_identical(
    passing_interval_raw$sensitivity$interval_checks[[1L]]$reason, ""
  )
  expect_invisible(.dpprior_validate_result_v1(passing_interval))
  passing_interval_raw$sensitivity$interval_checks[[1L]]$reason <-
    "fabricated successful audit message"
  class(passing_interval_raw) <- class(passing_interval)
  expect_error(.dpprior_validate_result_v1(passing_interval_raw),
               class = "dpprior_schema_error")

  divergent_interval_orders <- .schema23_sensitivity_interval()
  key <- divergent_interval_orders$sensitivity$scenarios$scenario_key[[1L]]
  divergent_interval_orders$sensitivity$interval_checks[[key]]$selected <- list(
    coverage = 0.2, lower_tail = 0.4, upper_tail = 0.4,
    coverage_residual = -0.6
  )
  changed_interval_metrics <- c(
    interval_achieved = 0.2, interval_residual = -0.6,
    interval_left_tail = 0.4, interval_right_tail = 0.4
  )
  for (metric in names(changed_interval_metrics)) {
    divergent_interval_orders$sensitivity$scenario_results[[metric]] <-
      unname(changed_interval_metrics[[metric]])
    row <- match(
      metric, divergent_interval_orders$sensitivity$metrics_long$metric
    )
    divergent_interval_orders$sensitivity$metrics_long$value[[row]] <-
      unname(changed_interval_metrics[[metric]])
  }
  expect_error(.dpprior_validate_result_v1(divergent_interval_orders),
               class = "dpprior_schema_error")

  coordinated_interval_forgery <- .schema23_sensitivity_interval()
  key <- coordinated_interval_forgery$sensitivity$scenarios$
    scenario_key[[1L]]
  forged_selected <- list(
    coverage = 0.8, lower_tail = 0.05, upper_tail = 0.15,
    coverage_residual = 0
  )
  forged_verifier <- c(
    forged_selected,
    list(
      tolerance = 1e-9, passed = FALSE,
      source = "independent_interval_backcheck"
    )
  )
  coordinated_interval_forgery$sensitivity$interval_checks[[key]]$selected <-
    forged_selected
  coordinated_interval_forgery$sensitivity$interval_checks[[key]]$
    verification <- forged_verifier
  changed_interval_metrics <- c(
    interval_achieved = 0.8, interval_residual = 0,
    interval_left_tail = 0.05, interval_right_tail = 0.15
  )
  for (metric in names(changed_interval_metrics)) {
    value <- unname(changed_interval_metrics[[metric]])
    coordinated_interval_forgery$sensitivity$scenario_results[[metric]] <-
      value
    row <- match(
      metric, coordinated_interval_forgery$sensitivity$metrics_long$metric
    )
    coordinated_interval_forgery$sensitivity$metrics_long$value[[row]] <- value
  }
  expect_error(.dpprior_validate_result_v1(coordinated_interval_forgery),
               class = "dpprior_schema_error")

  suppressed_interval <- .schema23_sensitivity_interval()
  key <- suppressed_interval$sensitivity$scenarios$scenario_key[[1L]]
  suppressed_interval$sensitivity$interval_checks[[key]] <- list(
    requested = NULL, selected = NULL, verification = NULL, status = NULL,
    source = NULL, usable = FALSE, verified = FALSE,
    reason = "not_interval_scenario"
  )
  for (metric in c(
    "interval_requested", "interval_achieved", "interval_residual",
    "interval_left_tail", "interval_right_tail"
  )) {
    suppressed_interval$sensitivity$scenario_results[[metric]] <- NA_real_
    row <- match(metric, suppressed_interval$sensitivity$metrics_long$metric)
    suppressed_interval$sensitivity$metrics_long$value[[row]] <- NA_real_
    suppressed_interval$sensitivity$metrics_long$reason[[row]] <-
      "not_interval_scenario"
  }
  expect_error(.dpprior_validate_result_v1(suppressed_interval),
               class = "dpprior_schema_error")

  unavailable_local <- .schema23_sensitivity()
  key <- unavailable_local$sensitivity$scenarios$scenario_key[[1L]]
  unavailable_local$sensitivity$local <- data.frame(
    scenario_key = key, axis = "mu_K", axis_value = 5,
    settings_key = .dp_sensitivity_content_key("settings-without-mu-K"),
    lower_scenario_key = NA_character_,
    upper_scenario_key = NA_character_, lower_value = NA_real_,
    upper_value = NA_real_, metric = "E_alpha", component = "alpha",
    derivative = NA_real_,
    method = "bracketed_secant_across_nearest_same-setting_neighbors",
    reason = "insufficient_neighbors", stringsAsFactors = FALSE
  )
  expect_error(.dpprior_validate_result_v1(unavailable_local),
               class = "dpprior_schema_error")

  forged_local <- unavailable_local
  forged_local$sensitivity$local$lower_scenario_key <- "ghost-lower"
  forged_local$sensitivity$local$upper_scenario_key <- "ghost-upper"
  forged_local$sensitivity$local$lower_value <- -Inf
  forged_local$sensitivity$local$upper_value <- Inf
  forged_local$sensitivity$local$component <- "weights"
  expect_error(.dpprior_validate_result_v1(forged_local),
               class = "dpprior_schema_error")
})


test_that("canonical objects survive serialization without shape drift", {
  objects <- c(
    list(
      .schema23_target(), .schema23_interval_target(),
      .schema23_family_target(), .schema23_infeasible_interval_target()
    ),
    lapply(
      c("a1_proxy", "a2_moment", "a2_kl", "dual_hard", "dual_soft",
        "dual_legacy"),
      .schema23_fit
    ),
    list(
      .schema23_hard_failed_rejected_candidates(),
      .schema23_hard_approximate_scan(),
      .schema23_a2_kl_approximate_initializer(),
      .schema23_soft_approximate_diagnostic(), .schema23_soft_endpoint()
    ),
    list(
      .schema23_diagnostics(), .schema23_diagnostics_approximate(),
      .schema23_sensitivity(), .schema23_sensitivity_failed(),
      .schema23_sensitivity_infeasible(), .schema23_sensitivity_interval(),
      .schema23_sensitivity_interval_infeasible(),
      .schema23_sensitivity_from_evidence(
        .schema23_sensitivity_route_evidence(
          "qualitative_confidence", confidence_explicit = FALSE
        )
      ),
      .schema23_sensitivity_from_evidence(
        .schema23_sensitivity_route_evidence(
          "coefficient_of_variation"
        )
      ),
      .schema23_sensitivity_from_evidence(
        .schema23_sensitivity_route_evidence(
          "strict_pmf", structural_K0 = TRUE
        )
      )
    )
  )
  for (object in objects) {
    restored <- unserialize(serialize(object, NULL, version = 3L))
    expect_identical(restored, object)
    expect_invisible(.dpprior_validate_object(restored))
  }
})
