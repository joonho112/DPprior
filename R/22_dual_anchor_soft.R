# =============================================================================
# Module 22: verified soft Dual-Anchor trade-off calibration
# =============================================================================

.DPPRIOR_SOFT_STATUS <- c(
  "converged", "boundary", "approximate", "infeasible", "failed"
)

.dpprior_soft_abort_invalid <- function(message, argument = NULL,
                                        value = NULL, expected = NULL,
                                        code = "dual_soft_invalid") {
  .dpprior_abort_invalid(
    message,
    c("dpprior_dual_soft_invalid_input", "dpprior_dual_anchor_error"),
    argument, value, expected, code
  )
}

.dpprior_soft_plain_character <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x) &&
    is.null(dim(x)) && !is.object(x)
}

.dpprior_soft_plain_logical <- function(x) {
  is.logical(x) && length(x) == 1L && !is.na(x) &&
    is.null(dim(x)) && !is.object(x)
}

.dpprior_soft_safe_condition_message <- function(
    condition, fallback = "captured condition without a safe message") {
  if (is.null(condition)) return(fallback)
  record <- if (is.list(condition)) {
    tryCatch(unclass(condition), error = function(e) NULL)
  } else {
    NULL
  }
  message <- if (is.list(record)) {
    record[["message", exact = TRUE]]
  } else {
    NULL
  }
  if (.dpprior_soft_plain_character(message)) message else fallback
}

.dpprior_soft_scalar <- function(x, name, lower = NULL, upper = NULL,
                                 lower_open = FALSE, upper_open = FALSE) {
  .dpprior_validate_scalar(
    x, name, lower = lower, upper = upper,
    lower_open = lower_open, upper_open = upper_open,
    .subclass = "dpprior_dual_soft_scalar_error"
  )
}

.dpprior_soft_integer <- function(x, name, lower = 1L) {
  value <- .dpprior_soft_scalar(x, name, lower = lower)
  if (value != floor(value)) {
    .dpprior_soft_abort_invalid(
      sprintf("%s must be an integer", name), name, x,
      sprintf("integer >= %d", lower), "dual_soft_integer"
    )
  }
  as.integer(value)
}

.dpprior_soft_control <- function(control, max_iter) {
  allowed <- c(
    "primary", "fallback", "boundary_tol", "verification_abs_tol",
    "verification_rel_tol", "objective_abs_tol", "objective_rel_tol",
    "stationarity_step", "stationarity_tol", ".optim_fun"
  )
  control <- .dpprior_v2_validate_named_list(
    control, "control", allowed = allowed
  )
  get_control <- function(name) control[[name, exact = TRUE]]
  primary <- get_control("primary") %||% list()
  fallback <- get_control("fallback") %||% list()
  primary <- .dpprior_v2_validate_named_list(
    primary, "control$primary",
    allowed = c(
      "maxit", "factr", "pgtol", "trace", "REPORT", "fnscale",
      "parscale", "ndeps"
    )
  )
  fallback <- .dpprior_v2_validate_named_list(
    fallback, "control$fallback",
    allowed = c(
      "maxit", "reltol", "alpha", "beta", "gamma", "trace", "REPORT",
      "fnscale", "parscale", "ndeps"
    )
  )
  validate_optimizer_geometry <- function(value, label) {
    fnscale <- value[["fnscale", exact = TRUE]]
    if (!is.null(fnscale) &&
        (!.dpprior_is_plain_numeric(fnscale) || length(fnscale) != 1L ||
         !is.finite(fnscale) || !identical(as.numeric(fnscale), 1))) {
      .dpprior_soft_abort_invalid(
        paste(label, "fnscale must be exactly 1; objective inversion is forbidden"),
        paste0(label, "$fnscale"), fnscale, "1",
        "dual_soft_truth_altering_optimizer_control"
      )
    }
    for (field in c("parscale", "ndeps")) {
      component <- value[[field, exact = TRUE]]
      if (!is.null(component) &&
          (!.dpprior_is_plain_numeric(component) || length(component) != 2L ||
           any(!is.finite(component)) || any(component <= 0))) {
        .dpprior_soft_abort_invalid(
          paste(label, field, "must be two positive plain numeric values"),
          paste0(label, "$", field), component,
          "positive finite numeric vector of length 2",
          "dual_soft_optimizer_geometry_control"
        )
      }
    }
  }
  validate_optimizer_geometry(primary, "control$primary")
  validate_optimizer_geometry(fallback, "control$fallback")
  normalize_scalar_control <- function(container, field, label, lower = NULL,
                                       upper = NULL, lower_open = FALSE,
                                       upper_open = FALSE,
                                       integer = FALSE) {
    value <- container[[field, exact = TRUE]]
    if (is.null(value)) return(container)
    normalized <- if (isTRUE(integer)) {
      .dpprior_soft_integer(value, paste0(label, "$", field),
                            lower = as.integer(lower %||% 0))
    } else {
      .dpprior_soft_scalar(
        value, paste0(label, "$", field), lower = lower, upper = upper,
        lower_open = lower_open, upper_open = upper_open
      )
    }
    container[[field]] <- normalized
    container
  }
  primary <- normalize_scalar_control(
    primary, "factr", "control$primary", lower = 0, lower_open = TRUE
  )
  primary <- normalize_scalar_control(
    primary, "pgtol", "control$primary", lower = 0
  )
  primary <- normalize_scalar_control(
    primary, "trace", "control$primary", lower = 0, integer = TRUE
  )
  primary <- normalize_scalar_control(
    primary, "REPORT", "control$primary", lower = 1, integer = TRUE
  )
  fallback <- normalize_scalar_control(
    fallback, "reltol", "control$fallback", lower = 0, lower_open = TRUE
  )
  fallback <- normalize_scalar_control(
    fallback, "alpha", "control$fallback", lower = 0, lower_open = TRUE
  )
  fallback <- normalize_scalar_control(
    fallback, "beta", "control$fallback", lower = 0, upper = 1,
    lower_open = TRUE, upper_open = TRUE
  )
  fallback <- normalize_scalar_control(
    fallback, "gamma", "control$fallback", lower = 1, lower_open = TRUE
  )
  fallback <- normalize_scalar_control(
    fallback, "trace", "control$fallback", lower = 0, integer = TRUE
  )
  fallback <- normalize_scalar_control(
    fallback, "REPORT", "control$fallback", lower = 1, integer = TRUE
  )
  optim_fun <- get_control(".optim_fun") %||% stats::optim
  if (!is.function(optim_fun)) {
    .dpprior_soft_abort_invalid(
      "control$.optim_fun must be a function when supplied",
      "control$.optim_fun", optim_fun, "function",
      "dual_soft_optimizer_adapter"
    )
  }
  primary[["maxit"]] <- primary[["maxit", exact = TRUE]] %||% max_iter
  fallback[["maxit"]] <- fallback[["maxit", exact = TRUE]] %||%
    (2L * max_iter)
  primary[["maxit"]] <- .dpprior_soft_integer(
    primary[["maxit", exact = TRUE]], "control$primary$maxit"
  )
  fallback[["maxit"]] <- .dpprior_soft_integer(
    fallback[["maxit", exact = TRUE]], "control$fallback$maxit"
  )
  boundary_tol <- .dpprior_soft_scalar(
    get_control("boundary_tol") %||% 1e-6,
    "control$boundary_tol", lower = 0
  )
  verification_abs_tol <- .dpprior_soft_scalar(
    get_control("verification_abs_tol") %||% 1e-8,
    "control$verification_abs_tol", lower = 0
  )
  verification_rel_tol <- .dpprior_soft_scalar(
    get_control("verification_rel_tol") %||% 1e-6,
    "control$verification_rel_tol", lower = 0
  )
  objective_abs_tol <- .dpprior_soft_scalar(
    get_control("objective_abs_tol") %||% 1e-10,
    "control$objective_abs_tol", lower = 0
  )
  objective_rel_tol <- .dpprior_soft_scalar(
    get_control("objective_rel_tol") %||% 1e-8,
    "control$objective_rel_tol", lower = 0
  )
  stationarity_step <- .dpprior_soft_scalar(
    get_control("stationarity_step") %||% 1e-5,
    "control$stationarity_step", lower = 0, lower_open = TRUE
  )
  stationarity_tol <- .dpprior_soft_scalar(
    get_control("stationarity_tol") %||% 1e-4,
    "control$stationarity_tol", lower = 0
  )
  policy_maxima <- c(
    boundary_tol = 1e-4,
    verification_abs_tol = 1e-8,
    verification_rel_tol = 1e-6,
    objective_abs_tol = 1e-10,
    objective_rel_tol = 1e-8,
    stationarity_step = 1e-3,
    stationarity_tol = 1e-4
  )
  supplied <- c(
    boundary_tol = boundary_tol,
    verification_abs_tol = verification_abs_tol,
    verification_rel_tol = verification_rel_tol,
    objective_abs_tol = objective_abs_tol,
    objective_rel_tol = objective_rel_tol,
    stationarity_step = stationarity_step,
    stationarity_tol = stationarity_tol
  )
  relaxed <- names(supplied)[supplied > policy_maxima[names(supplied)]]
  if (!identical(stationarity_step, 1e-5)) {
    relaxed <- unique(c(relaxed, "stationarity_step"))
  }
  if (length(relaxed)) {
    .dpprior_soft_abort_invalid(
      paste(
        "soft verification controls must remain within the fixed truth-making",
        "policy (stationarity_step is exactly 1e-5):",
        paste(relaxed, collapse = ", ")
      ),
      "control", control, policy_maxima,
      "dual_soft_noncanonical_tolerance"
    )
  }
  list(
    max_iter = as.integer(max_iter),
    primary = primary,
    fallback = fallback,
    boundary_tol = boundary_tol,
    verification_abs_tol = verification_abs_tol,
    verification_rel_tol = verification_rel_tol,
    objective_abs_tol = objective_abs_tol,
    objective_rel_tol = objective_rel_tol,
    stationarity_step = stationarity_step,
    stationarity_tol = stationarity_tol,
    policy_maxima = policy_maxima,
    optim_fun = optim_fun,
    optim_fun_label = if (identical(optim_fun, stats::optim)) {
      "stats::optim"
    } else {
      "injected_optimizer_adapter"
    }
  )
}

.dpprior_soft_normalize_fit <- function(fit) {
  if (!exists(".dpprior_v2_normalize_fit", mode = "function", inherits = TRUE)) {
    stop("Phase 8 shared fit normalizer is unavailable", call. = FALSE)
  }
  normalized <- .dpprior_v2_normalize_fit(fit)
  required <- c("a", "b", "J", "target_K", "status", "usable", "verified")
  if (!is.list(normalized) || !all(required %in% names(normalized))) {
    .dpprior_soft_abort_invalid(
      "fit did not satisfy the Phase 8 normalized-fit contract", "fit", fit,
      paste(required, collapse = ", "), "dual_soft_fit_contract"
    )
  }
  normalized
}

.dpprior_soft_normalize_target <- function(target) {
  if (!exists(".dpprior_v2_normalize_weight_spec", mode = "function",
              inherits = TRUE)) {
    stop("Phase 8 shared weight-target normalizer is unavailable", call. = FALSE)
  }
  normalized <- .dpprior_v2_normalize_weight_spec(target, mode = "soft")
  if (identical(normalized$metric, "wmax_tail")) {
    .dpprior_soft_abort_invalid(
      paste(
        "metric='wmax_tail' is not available for soft calibration until a",
        "direct estimator passes its independent error-control contract; use",
        "metric='wmax_tail_upper' for the certified conservative bound"
      ),
      "target$metric", normalized$metric,
      "wsb_tail, wsb_mean, wsb_quantile, or wmax_tail_upper",
      "dual_soft_wmax_direct_deferred"
    )
  }
  normalized
}

.dpprior_soft_eval <- function(a, b, fit_info, target, M) {
  moments <- tryCatch(
    exact_K_moments(fit_info$J, a, b, M),
    error = identity
  )
  if (inherits(moments, "condition") || !is.list(moments) ||
      !all(c("mean", "var") %in% names(moments)) ||
      any(!is.finite(c(moments$mean, moments$var)))) {
    return(list(
      ok = FALSE,
      condition = if (inherits(moments, "condition")) moments else NULL,
      reason = "K moment evaluation failed"
    ))
  }
  metric <- tryCatch(
    .dpprior_v2_eval_metric(
      target, a = a, b = b, J = fit_info$J, M = M
    ),
    error = identity
  )
  if (inherits(metric, "condition") || !is.list(metric) ||
      !is.numeric(metric$value) || length(metric$value) != 1L ||
      !is.finite(metric$value)) {
    return(list(
      ok = FALSE,
      condition = if (inherits(metric, "condition")) metric else NULL,
      reason = "weight metric evaluation failed"
    ))
  }
  list(
    ok = TRUE,
    K = list(mean = as.numeric(moments$mean),
             variance = as.numeric(moments$var), M = as.integer(M)),
    weight = metric
  )
}

.dpprior_soft_weight_loss <- function(achieved, target) {
  raw <- achieved - target$value
  directed <- switch(
    target$relation,
    target = raw,
    at_most = max(0, raw),
    at_least = max(0, -raw),
    stop("unreachable soft relation", call. = FALSE)
  )
  list(
    value = directed^2,
    residual = raw,
    directed_residual = directed,
    scale = 1,
    scaled = directed,
    scale_formula = "fixed probability/weight scale s_D = 1",
    loss_formula = switch(
      target$relation,
      target = "((achieved - value) / s_D)^2",
      at_most = "(max(0, achieved - value) / s_D)^2",
      at_least = "(max(0, value - achieved) / s_D)^2"
    )
  )
}

.dpprior_soft_losses <- function(evaluation, fit_info, target, lambda,
                                 scales) {
  K <- .dpprior_v2_k_loss(
    list(mean = evaluation$K$mean, variance = evaluation$K$variance),
    list(mean = fit_info$target_K$mu_K,
         variance = fit_info$target_K$var_K),
    scales = scales
  )
  weight <- .dpprior_soft_weight_loss(evaluation$weight$value, target)
  list(
    K = K,
    weight = weight,
    total = lambda * K$value + (1 - lambda) * weight$value
  )
}

.dpprior_soft_objective <- function(fit_info, target, lambda, M, scales,
                                    log_bounds) {
  force(fit_info)
  force(target)
  force(lambda)
  force(M)
  force(scales)
  force(log_bounds)
  function(eta) {
    if (!.dpprior_is_plain_numeric(eta) || length(eta) != 2L ||
        any(!is.finite(eta)) || any(eta < log_bounds[1L]) ||
        any(eta > log_bounds[2L])) {
      return(.PENALTY_INF)
    }
    evaluation <- .dpprior_soft_eval(
      exp(eta[1L]), exp(eta[2L]), fit_info, target, M
    )
    if (!isTRUE(evaluation$ok)) {
      return(.PENALTY_INF)
    }
    losses <- .dpprior_soft_losses(
      evaluation, fit_info, target, lambda, scales
    )
    if (is.finite(losses$total)) losses$total else .PENALTY_INF
  }
}

.dpprior_soft_call_optim <- function(optim_fun, method, start, objective,
                                     log_bounds, control) {
  warnings <- character()
  started <- proc.time()[["elapsed"]]
  value <- withCallingHandlers(
    tryCatch({
      if (identical(method, "L-BFGS-B")) {
        optim_fun(
          par = start, fn = objective, method = method,
          lower = rep(log_bounds[1L], 2L),
          upper = rep(log_bounds[2L], 2L), control = control
        )
      } else {
        optim_fun(par = start, fn = objective, method = method,
                  control = control)
      }
    }, error = identity),
    warning = function(w) {
      warnings <<- c(warnings, .dpprior_soft_safe_condition_message(w))
      invokeRestart("muffleWarning")
    }
  )
  elapsed <- proc.time()[["elapsed"]] - started
  if (inherits(value, "condition")) {
    return(list(
      result = NULL, condition = value, warnings = warnings,
      elapsed = elapsed
    ))
  }
  list(result = value, condition = NULL, warnings = warnings, elapsed = elapsed)
}

.dpprior_soft_attempt <- function(call, method, start, bounds, control,
                                  objective) {
  opt <- call$result
  candidate <- NULL
  ordinary_opt <- if (is.list(opt) && !is.object(opt)) opt else NULL
  opt_names <- if (is.list(ordinary_opt)) names(ordinary_opt) else NULL
  list_schema <- is.list(ordinary_opt) &&
    !is.null(opt_names) && !anyNA(opt_names) && all(nzchar(opt_names)) &&
    !anyDuplicated(opt_names) &&
    all(c("par", "value", "convergence") %in% opt_names) &&
    !length(setdiff(opt_names, c(
      "par", "value", "counts", "convergence", "message", "hessian"
    )))
  par <- if (is.list(ordinary_opt)) ordinary_opt[["par", exact = TRUE]] else NULL
  reported_value <- if (is.list(ordinary_opt)) {
    ordinary_opt[["value", exact = TRUE]]
  } else NULL
  convergence <- if (is.list(ordinary_opt)) {
    ordinary_opt[["convergence", exact = TRUE]]
  } else NULL
  message <- if (is.list(ordinary_opt)) {
    ordinary_opt[["message", exact = TRUE]]
  } else NULL
  counts_record <- if (is.list(ordinary_opt)) {
    ordinary_opt[["counts", exact = TRUE]]
  } else NULL
  counts_names <- if (!is.null(counts_record)) names(counts_record) else NULL
  counts_schema <- is.null(counts_record) ||
    (.dpprior_is_plain_numeric(counts_record) && length(counts_record) <= 2L &&
       all(is.na(counts_record) |
           (is.finite(counts_record) & counts_record >= 0 &
              counts_record == floor(counts_record))) &&
       !is.null(counts_names) && !anyNA(counts_names) &&
       all(nzchar(counts_names)) && !anyDuplicated(counts_names) &&
       !length(setdiff(counts_names, c("function", "gradient"))))
  message_schema <- is.null(message) || .dpprior_soft_plain_character(message)
  valid_schema <- list_schema && counts_schema && message_schema &&
    .dpprior_is_plain_numeric(par) && length(par) == 2L &&
    all(is.finite(par)) && .dpprior_is_plain_numeric(reported_value) &&
    length(reported_value) == 1L && is.finite(reported_value) &&
    .dpprior_is_plain_numeric(convergence) && length(convergence) == 1L &&
    is.finite(convergence) &&
    convergence == floor(convergence)
  if (valid_schema) {
    candidate <- c(
      log_shape = unname(par[1L]), log_rate = unname(par[2L]),
      shape = exp(unname(par[1L])), rate = exp(unname(par[2L]))
    )
  }
  objective_value <- if (valid_schema) {
    as.numeric(reported_value)
  } else if (!is.null(candidate)) {
    objective(candidate[c("log_shape", "log_rate")])
  } else {
    NA_real_
  }
  counts <- if (is.list(counts_record) && !is.object(counts_record)) {
    counts_record
  } else if (.dpprior_is_plain_numeric(counts_record)) {
    as.list(counts_record)
  } else {
    list()
  }
  valid_count <- function(value) {
    .dpprior_is_plain_numeric(value) && length(value) == 1L &&
      is.finite(value) && value >= 0 && value == floor(value) &&
      value <= .Machine$integer.max
  }
  iteration_count <- counts[["function", exact = TRUE]]
  evaluation_count <- counts[["gradient", exact = TRUE]]
  .dpprior_v2_make_attempt(
    method = method,
    start = stats::setNames(as.numeric(start), c("log_shape", "log_rate")),
    bounds = list(
      lower = stats::setNames(rep(bounds[1L], 2L),
                              c("log_shape", "log_rate")),
      upper = stats::setNames(rep(bounds[2L], 2L),
                              c("log_shape", "log_rate")),
      parameterization = "log(shape), log(rate)"
    ),
    control = control,
    exit_code = if (valid_schema) as.integer(convergence) else NA_integer_,
    message = if (!is.null(call$condition)) {
      .dpprior_soft_safe_condition_message(call$condition)
    } else if (valid_schema && .dpprior_soft_plain_character(message)) {
      message
    } else if (!is.null(opt) && !valid_schema) {
      "optimizer returned an invalid result schema"
    } else {
      NA_character_
    },
    counts = list(
      iterations = if (valid_count(iteration_count)) {
        as.integer(iteration_count)
      } else {
        NA_integer_
      },
      evaluations = if (valid_count(evaluation_count)) {
        as.integer(evaluation_count)
      } else {
        NA_integer_
      }
    ),
    objective = objective_value,
    elapsed = as.numeric(call$elapsed),
    warning = if (length(call$warnings)) paste(call$warnings, collapse = "; ")
      else NA_character_,
    error = if (!is.null(call$condition)) {
      .dpprior_soft_safe_condition_message(call$condition)
    } else if (!is.null(opt) && !valid_schema) {
      "invalid_optimizer_result"
    } else NA_character_,
    candidate = candidate
  )
}

.dpprior_soft_attempt_candidate <- function(attempt) {
  candidate <- attempt[["candidate", exact = TRUE]]
  if (!.dpprior_is_plain_numeric(candidate) ||
      !all(c("log_shape", "log_rate", "shape", "rate") %in%
           names(candidate)) || any(!is.finite(candidate))) {
    return(NULL)
  }
  list(
    eta = as.numeric(candidate[c("log_shape", "log_rate")]),
    a = as.numeric(candidate[["shape"]]),
    b = as.numeric(candidate[["rate"]]),
    objective = as.numeric(
      attempt[["candidate_objective", exact = TRUE]] %||%
        attempt[["objective", exact = TRUE]] %||% NA_real_
    ),
    exit_code = as.integer(
      attempt[["exit_code", exact = TRUE]] %||%
        attempt[["exit", exact = TRUE]] %||% NA_integer_
    ),
    method = attempt[["method", exact = TRUE]]
  )
}

.dpprior_soft_stability <- function(selected, verification, abs_tol, rel_tol) {
  selected_values <- c(
    K_mean = selected$K$mean,
    K_variance = selected$K$variance,
    weight_metric = selected$weight$value
  )
  verification_values <- c(
    K_mean = verification$K$mean,
    K_variance = verification$K$variance,
    weight_metric = verification$weight$value
  )
  difference <- abs(verification_values - selected_values)
  tolerance <- abs_tol + rel_tol * pmax(
    abs(selected_values), abs(verification_values), 1
  )
  component_pass <- is.finite(difference) & difference <= tolerance
  list(
    selected = selected_values,
    verification = verification_values,
    absolute_difference = difference,
    tolerance = tolerance,
    tolerance_formula =
      "abs_tol + rel_tol * max(abs(selected), abs(verification), 1)",
    component_pass = component_pass,
    passed = length(component_pass) == 3L && !anyNA(component_pass) &&
      all(component_pass)
  )
}

.dpprior_soft_weight_verified <- function(metric) {
  finite <- is.list(metric) &&
    .dpprior_is_plain_numeric(metric$value) && length(metric$value) == 1L &&
    is.finite(metric$value) && metric$value >= 0 && metric$value <= 1
  if (!finite || !.dpprior_soft_plain_character(metric$status)) {
    return(FALSE)
  }
  if (identical(metric$metric, "wsb_mean")) {
    # One quadrature order is deliberately labelled approximate by the shared
    # evaluator. The soft verifier certifies the mean only by comparing the
    # selected and independently refined orders below.
    return(metric$status %in% c("approximate", "converged", "boundary"))
  }
  metric$status %in% c("converged", "boundary") &&
    .dpprior_soft_plain_logical(metric$usable) && isTRUE(metric$usable) &&
    .dpprior_soft_plain_logical(metric$verified) && isTRUE(metric$verified)
}

.dpprior_soft_optimality <- function(candidate, fit_info, target, lambda,
                                     M_fit, M_verify, scales, controls,
                                     log_bounds, endpoint) {
  if (isTRUE(endpoint)) {
    return(list(
      performed = FALSE,
      passed = FALSE,
      recorded_objective = NULL,
      recomputed_objective = NULL,
      objective_difference = NULL,
      objective_tolerance = NULL,
      objective_passed = NULL,
      gradient = NULL,
      gradient_method = NULL,
      bound_state = NULL,
      boundary_tolerance = NULL,
      stationarity_operator = NULL,
      stationarity_tolerance = NULL,
      component_pass = NULL,
      stationarity_passed = NULL,
      local_base_objective = NULL,
      neighbor_objectives = NULL,
      neighbor_tolerances = NULL,
      local_minimum_passed = NULL,
      start_objective = NULL,
      candidate_objective = NULL,
      start_tolerance = NULL,
      no_worse_start = NULL,
      selection_tolerance = NULL,
      source = "independent_refined_objective_verification",
      unavailable_reason = "lambda=1 endpoint performs no optimization"
    ))
  }
  objective_selected <- .dpprior_soft_objective(
    fit_info, target, lambda, M_fit, scales, log_bounds
  )
  objective <- .dpprior_soft_objective(
    fit_info, target, lambda, M_verify, scales, log_bounds
  )
  eta <- log(c(candidate$a, candidate$b))
  names(eta) <- c("log_shape", "log_rate")
  base <- objective(eta)
  selected_base <- objective_selected(eta)
  recorded <- as.numeric(candidate$objective)
  objective_tolerance <- controls$objective_abs_tol +
    controls$objective_rel_tol * max(abs(recorded), abs(selected_base), 1)
  objective_bound <- is.finite(recorded) && is.finite(selected_base) &&
    abs(recorded - selected_base) <= objective_tolerance
  gradient <- component_pass <- rep(NA, 2L)
  names(gradient) <- names(component_pass) <- c("log_shape", "log_rate")
  method <- bound_state <- stationarity_operator <- character(2L)
  names(method) <- names(bound_state) <- names(stationarity_operator) <-
    names(gradient)
  for (index in 1:2) {
    lower_distance <- eta[index] - log_bounds[1L]
    upper_distance <- log_bounds[2L] - eta[index]
    step <- min(
      controls$stationarity_step,
      max(lower_distance, 0) / 2,
      max(upper_distance, 0) / 2
    )
    on_lower <- lower_distance <= controls$boundary_tol
    on_upper <- upper_distance <= controls$boundary_tol
    if (!on_lower && !on_upper && is.finite(step) && step > 0) {
      plus <- minus <- eta
      plus[index] <- plus[index] + step
      minus[index] <- minus[index] - step
      gradient[index] <- (objective(plus) - objective(minus)) / (2 * step)
      component_pass[index] <- is.finite(gradient[index]) &&
        abs(gradient[index]) <= controls$stationarity_tol
      method[index] <- "central finite-difference interior gradient"
      bound_state[index] <- "interior"
      stationarity_operator[index] <- "abs_lte"
    } else if (on_lower) {
      step <- min(controls$stationarity_step,
                  max(upper_distance, 0) / 2)
      plus <- eta
      plus[index] <- plus[index] + step
      gradient[index] <- (objective(plus) - base) / step
      component_pass[index] <- is.finite(gradient[index]) &&
        gradient[index] >= -controls$stationarity_tol
      method[index] <- "forward feasible-direction KKT gradient"
      bound_state[index] <- "lower"
      stationarity_operator[index] <- "gte"
    } else {
      step <- min(controls$stationarity_step,
                  max(lower_distance, 0) / 2)
      minus <- eta
      minus[index] <- minus[index] - step
      gradient[index] <- (base - objective(minus)) / step
      component_pass[index] <- is.finite(gradient[index]) &&
        gradient[index] <= controls$stationarity_tol
      method[index] <- "backward feasible-direction KKT gradient"
      bound_state[index] <- "upper"
      stationarity_operator[index] <- "lte"
    }
  }
  stationarity_passed <- !anyNA(component_pass) && all(component_pass)
  local_step <- 1e-3
  directions <- rbind(
    log_shape_plus = c(1, 0), log_shape_minus = c(-1, 0),
    log_rate_plus = c(0, 1), log_rate_minus = c(0, -1),
    diagonal_pp = c(1, 1) / sqrt(2),
    diagonal_pm = c(1, -1) / sqrt(2),
    diagonal_mp = c(-1, 1) / sqrt(2),
    diagonal_mm = c(-1, -1) / sqrt(2)
  )
  local_values <- rep(NA_real_, nrow(directions))
  names(local_values) <- rownames(directions)
  for (index in seq_len(nrow(directions))) {
    proposal <- eta + local_step * directions[index, ]
    if (all(proposal >= log_bounds[1L]) &&
        all(proposal <= log_bounds[2L])) {
      local_values[index] <- objective(proposal)
    }
  }
  local_values <- local_values[is.finite(local_values)]
  local_tolerance <- controls$objective_abs_tol +
    controls$objective_rel_tol * pmax(abs(base), abs(local_values), 1)
  names(local_tolerance) <- names(local_values)
  local_component_pass <- base <= local_values + local_tolerance
  local_minimum_passed <- length(local_values) > 0L &&
    all(local_component_pass)
  reference_eta <- log(c(fit_info$a, fit_info$b))
  reference_objective <- objective(reference_eta)
  reference_tolerance <- controls$objective_abs_tol +
    controls$objective_rel_tol * max(abs(base), abs(reference_objective), 1)
  reference_passed <- is.finite(reference_objective) && is.finite(base) &&
    base <= reference_objective + reference_tolerance
  list(
    performed = TRUE,
    passed = objective_bound && stationarity_passed && local_minimum_passed &&
      reference_passed,
    recorded_objective = recorded,
    recomputed_objective = selected_base,
    objective_difference = abs(recorded - selected_base),
    objective_tolerance = objective_tolerance,
    objective_passed = objective_bound,
    gradient = gradient,
    gradient_method = method,
    bound_state = bound_state,
    boundary_tolerance = controls$boundary_tol,
    stationarity_operator = stationarity_operator,
    stationarity_tolerance = controls$stationarity_tol,
    component_pass = component_pass,
    stationarity_passed = stationarity_passed,
    local_base_objective = base,
    neighbor_objectives = local_values,
    neighbor_tolerances = local_tolerance,
    local_minimum_passed = local_minimum_passed,
    start_objective = reference_objective,
    candidate_objective = selected_base,
    start_tolerance = reference_tolerance,
    no_worse_start = reference_passed,
    selection_tolerance = 0,
    source = "independent_refined_objective_verification",
    unavailable_reason = NULL
  )
}

.dpprior_soft_verify <- function(candidate, fit_info, target, lambda,
                                 M_fit, M_verify, scales, controls,
                                 log_bounds, endpoint = FALSE) {
  selected <- .dpprior_soft_eval(
    candidate$a, candidate$b, fit_info, target, M_fit
  )
  verification <- .dpprior_soft_eval(
    candidate$a, candidate$b, fit_info, target, M_verify
  )
  performed <- isTRUE(selected$ok) && isTRUE(verification$ok)
  if (!performed) {
    return(list(
      performed = TRUE, passed = FALSE,
      selected_order = M_fit, verification_order = M_verify,
      selected = selected, verification = verification,
      stability = NULL,
      reason = "selected or independent verification evaluation failed"
    ))
  }
  selected_losses <- .dpprior_soft_losses(
    selected, fit_info, target, lambda, scales
  )
  verification_losses <- .dpprior_soft_losses(
    verification, fit_info, target, lambda, scales
  )
  stability <- .dpprior_soft_stability(
    selected, verification,
    controls$verification_abs_tol, controls$verification_rel_tol
  )
  optimality <- .dpprior_soft_optimality(
    candidate, fit_info, target, lambda, M_fit, M_verify, scales, controls,
    log_bounds, endpoint
  )
  endpoint_target <- if (isTRUE(endpoint)) {
    input_reference <- fit_info[["source", exact = TRUE]][[
      "canonical_reference", exact = TRUE
    ]]
    selected_reference <- input_reference[[
      "selected_snapshot", exact = TRUE
    ]][["achieved_K", exact = TRUE]]
    verifier_reference <- input_reference[[
      "verifier_snapshot", exact = TRUE
    ]][["achieved_K", exact = TRUE]]
    selected_values <- c(
      selected.mean = selected$K$mean,
      selected.variance = selected$K$variance,
      verifier.mean = verification$K$mean,
      verifier.variance = verification$K$variance
    )
    reference_values <- c(
      selected.mean = selected_reference[["mean", exact = TRUE]],
      selected.variance = selected_reference[["variance", exact = TRUE]],
      verifier.mean = verifier_reference[["mean", exact = TRUE]],
      verifier.variance = verifier_reference[["variance", exact = TRUE]]
    )
    difference <- abs(selected_values - reference_values)
    tolerance <- controls$verification_abs_tol +
      controls$verification_rel_tol * pmax(abs(reference_values), 1)
    component_pass <- is.finite(difference) & difference <= tolerance
    reference_pass <- identical(
      input_reference[["status", exact = TRUE]] %in%
        c("converged", "boundary"), TRUE
    ) && isTRUE(input_reference[["usable", exact = TRUE]]) &&
      isTRUE(input_reference[["verified", exact = TRUE]])
    list(
      performed = TRUE, selected = selected_values,
      target = reference_values,
      absolute_difference = difference, tolerance = tolerance,
      component_pass = component_pass,
      passed = reference_pass && !anyNA(component_pass) && all(component_pass),
      reason = paste(
        "lambda=1 requires fresh selected/refined K moments to match the",
        "validated canonical input-fit snapshots"
      )
    )
  } else {
    list(
      performed = FALSE, passed = TRUE,
      reason = "not applicable away from the exact K-only endpoint"
    )
  }
  passed <- isTRUE(stability$passed) &&
    .dpprior_soft_weight_verified(selected$weight) &&
    .dpprior_soft_weight_verified(verification$weight) &&
    is.finite(selected_losses$total) && is.finite(verification_losses$total) &&
    (isTRUE(endpoint) || isTRUE(optimality$passed)) &&
    isTRUE(endpoint_target$passed)
  list(
    method = "fresh higher-order K moments and named weight metric",
    performed = TRUE,
    passed = passed,
    selected_order = M_fit,
    verification_order = M_verify,
    selected = selected,
    verification = verification,
    selected_losses = selected_losses,
    verification_losses = verification_losses,
    stability = stability,
    optimality = optimality,
    endpoint_target = endpoint_target,
    reason = if (passed) "independent soft-tradeoff verification passed"
      else "independent soft-tradeoff verification failed"
  )
}

.dpprior_soft_legacy_achieved <- function(a, b, M) {
  list(
    mean = mean_w1(a, b, M),
    prob_gt_50 = prob_wsb_exceeds(0.5, a, b),
    prob_gt_90 = prob_wsb_exceeds(0.9, a, b)
  )
}

.dpprior_soft_canonical_controls <- function(controls, max_iter, log_bounds) {
  normalize_pair <- function(value, default) {
    value <- value %||% default
    stats::setNames(as.numeric(value), c("log_shape", "log_rate"))
  }
  primary <- controls[["primary", exact = TRUE]]
  fallback <- controls[["fallback", exact = TRUE]]
  list(
    max_iter = as.integer(max_iter),
    log_bounds = unname(as.numeric(log_bounds)),
    primary = list(
      maxit = as.integer(primary[["maxit", exact = TRUE]]),
      fnscale = as.numeric(primary[["fnscale", exact = TRUE]] %||% 1),
      parscale = normalize_pair(
        primary[["parscale", exact = TRUE]], c(1, 1)
      ),
      ndeps = normalize_pair(
        primary[["ndeps", exact = TRUE]], c(1e-3, 1e-3)
      )
    ),
    fallback = list(
      maxit = as.integer(fallback[["maxit", exact = TRUE]]),
      reltol = as.numeric(
        fallback[["reltol", exact = TRUE]] %||%
          sqrt(.Machine$double.eps)
      ),
      fnscale = as.numeric(fallback[["fnscale", exact = TRUE]] %||% 1),
      parscale = normalize_pair(
        fallback[["parscale", exact = TRUE]], c(1, 1)
      ),
      ndeps = normalize_pair(
        fallback[["ndeps", exact = TRUE]], c(1e-3, 1e-3)
      )
    ),
    boundary_tol = controls[["boundary_tol", exact = TRUE]],
    verification_abs_tol = controls[["verification_abs_tol", exact = TRUE]],
    verification_rel_tol = controls[["verification_rel_tol", exact = TRUE]],
    objective_abs_tol = controls[["objective_abs_tol", exact = TRUE]],
    objective_rel_tol = controls[["objective_rel_tol", exact = TRUE]],
    stationarity_step = controls[["stationarity_step", exact = TRUE]],
    stationarity_tol = controls[["stationarity_tol", exact = TRUE]],
    local_neighbor_step = 1e-3,
    selection_tolerance = 0,
    optimizer_adapter = controls[["optim_fun_label", exact = TRUE]]
  )
}

.dpprior_soft_canonical_tolerances <- function(controls) {
  stability <- list(
    absolute = controls[["verification_abs_tol", exact = TRUE]],
    relative = controls[["verification_rel_tol", exact = TRUE]],
    scale_floor = 1
  )
  list(
    K = stability,
    weight = stability,
    objective = list(
      absolute = controls[["objective_abs_tol", exact = TRUE]],
      relative = controls[["objective_rel_tol", exact = TRUE]],
      scale_floor = 1
    ),
    boundary = controls[["boundary_tol", exact = TRUE]],
    stationarity = list(
      tolerance = controls[["stationarity_tol", exact = TRUE]],
      step = controls[["stationarity_step", exact = TRUE]]
    ),
    neighborhood = list(step = 1e-3),
    selection = 0
  )
}

.dpprior_soft_canonical_weight_target <- function(target) {
  raw_input <- target[["raw_input", exact = TRUE]]
  relation_from <- raw_input[["relation", exact = TRUE]]
  request <- list(
    metric = target[["metric", exact = TRUE]],
    relation = relation_from,
    value = target[["value", exact = TRUE]],
    threshold = target[["threshold", exact = TRUE]],
    probability = target[["probability", exact = TRUE]]
  )
  normalized <- list(
    metric = target[["metric", exact = TRUE]],
    relation = target[["relation", exact = TRUE]],
    value = target[["value", exact = TRUE]],
    threshold = target[["threshold", exact = TRUE]],
    probability = target[["probability", exact = TRUE]]
  )
  transformed <- !identical(request, normalized)
  .dpprior_new_weight_target(
    request = request,
    normalized = normalized,
    used = normalized,
    metric = target[["metric", exact = TRUE]],
    relation = target[["relation", exact = TRUE]],
    operator = target[["operator", exact = TRUE]],
    value = target[["value", exact = TRUE]],
    threshold = target[["threshold", exact = TRUE]],
    probability = target[["probability", exact = TRUE]],
    estimand = target[["estimand", exact = TRUE]],
    units = target[["units", exact = TRUE]],
    certification = if (identical(
      target[["metric", exact = TRUE]], "wmax_tail_upper"
    )) {
      list(
        kind = "upper_bound", certified = TRUE, passed = TRUE,
        method = "certified_size_biased_mass_upper_bound",
        source = "wmax_tail_bounds"
      )
    } else {
      list()
    },
    provenance = list(
      source = "DPprior_dual_soft request canonicalization",
      transformation = if (transformed) {
        list(
          rule = "canonicalize_soft_weight_target",
          opt_in = FALSE,
          before = request,
          after = normalized,
          evidence = list(
            mode = "soft", value_field = "value",
            relation_from = relation_from,
            relation_to = target[["relation", exact = TRUE]],
            probability_field = if (is.null(
              target[["probability", exact = TRUE]]
            )) "none" else "probability"
          )
        )
      } else {
        NULL
      },
      selection = NULL,
      probability_alias_used = "prob" %in% names(raw_input)
    )
  )
}

.dpprior_soft_achieved <- function(evaluation, target, M,
                                   K_override = NULL) {
  K <- if (is.null(K_override)) {
    list(
      mean = evaluation[["K", exact = TRUE]][["mean", exact = TRUE]],
      variance = evaluation[["K", exact = TRUE]][["variance", exact = TRUE]],
      estimand = "K_J", source = "selected", M = as.integer(M)
    )
  } else {
    K_override
  }
  list(
    K = K,
    weight = list(
      metric = target[["metric", exact = TRUE]],
      value = evaluation[["weight", exact = TRUE]][["value", exact = TRUE]],
      source = "selected"
    )
  )
}

.dpprior_soft_residuals <- function(achieved, fit_info, target) {
  raw_weight <- achieved[["weight", exact = TRUE]][["value", exact = TRUE]] -
    target[["value", exact = TRUE]]
  directed_weight <- switch(
    target[["relation", exact = TRUE]],
    target = raw_weight,
    at_most = max(0, raw_weight),
    at_least = max(0, -raw_weight)
  )
  list(
    K = list(
      mean = achieved[["K", exact = TRUE]][["mean", exact = TRUE]] -
        fit_info[["target_K", exact = TRUE]][["mu_K", exact = TRUE]],
      variance = achieved[["K", exact = TRUE]][["variance", exact = TRUE]] -
        fit_info[["target_K", exact = TRUE]][["var_K", exact = TRUE]]
    ),
    weight = list(raw = raw_weight, directed = directed_weight)
  )
}

.dpprior_soft_snapshot <- function(evaluation, parameters, fit_info, target,
                                   M, tolerances, source,
                                   K_override = NULL) {
  achieved <- .dpprior_soft_achieved(
    evaluation, target, M, K_override = K_override
  )
  .dpprior_new_snapshot(
    parameters = parameters,
    M = as.integer(M),
    achieved = achieved,
    residuals = .dpprior_soft_residuals(achieved, fit_info, target),
    tolerances = tolerances,
    finite = TRUE,
    source = source
  )
}

.dpprior_soft_canonical_stability <- function(selected, verifier,
                                              tolerances) {
  selected_values <- c(
    K.mean = selected[["achieved", exact = TRUE]][["K", exact = TRUE]][[
      "mean", exact = TRUE
    ]],
    K.variance = selected[["achieved", exact = TRUE]][["K", exact = TRUE]][[
      "variance", exact = TRUE
    ]],
    weight.value = selected[["achieved", exact = TRUE]][[
      "weight", exact = TRUE
    ]][["value", exact = TRUE]]
  )
  verifier_values <- c(
    K.mean = verifier[["achieved", exact = TRUE]][["K", exact = TRUE]][[
      "mean", exact = TRUE
    ]],
    K.variance = verifier[["achieved", exact = TRUE]][["K", exact = TRUE]][[
      "variance", exact = TRUE
    ]],
    weight.value = verifier[["achieved", exact = TRUE]][[
      "weight", exact = TRUE
    ]][["value", exact = TRUE]]
  )
  delta <- abs(selected_values - verifier_values)
  tolerance <- c(
    K.mean = tolerances$K$absolute + tolerances$K$relative * max(
      abs(selected_values[["K.mean"]]), abs(verifier_values[["K.mean"]]),
      tolerances$K$scale_floor
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
    scale_floor = c(K.mean = 1, K.variance = 1, weight.value = 1),
    source = "independent_verifier"
  )
}

.dpprior_soft_unavailable_optimality <- function(reason) {
  out <- setNames(
    rep(list(NULL), length(.DPPRIOR_SOFT_OPTIMALITY_FIELDS)),
    .DPPRIOR_SOFT_OPTIMALITY_FIELDS
  )
  out[["performed"]] <- FALSE
  out[["passed"]] <- FALSE
  out[["source"]] <- "independent_refined_objective_verification"
  out[["unavailable_reason"]] <- reason
  out
}

.dpprior_soft_canonical_attempts <- function(attempts, selected_index = NULL,
                                             parameterization,
                                             retain_candidates = TRUE) {
  lapply(seq_along(attempts), function(index) {
    old <- attempts[[index]]
    raw_candidate <- old[["candidate", exact = TRUE]]
    candidate_parameters <- if (isTRUE(retain_candidates) &&
      .dpprior_is_plain_numeric(raw_candidate) &&
        all(c("shape", "rate") %in% names(raw_candidate)) &&
        all(is.finite(raw_candidate[c("shape", "rate")])) &&
        all(raw_candidate[c("shape", "rate")] > 0)
    ) {
      .dpprior_new_parameters(
        as.numeric(raw_candidate[["shape"]]),
        as.numeric(raw_candidate[["rate"]]), parameterization
      )
    } else {
      NULL
    }
    scalar_or_null <- function(value, count = FALSE) {
      valid <- .dpprior_is_plain_numeric(value) && length(value) == 1L &&
        !is.na(value) && is.finite(value) && (!count || value >= 0)
      if (!valid) NULL else if (count) as.integer(value) else as.numeric(value)
    }
    exit_code <- scalar_or_null(old[["exit_code", exact = TRUE]], TRUE)
    iterations <- scalar_or_null(old[["iterations", exact = TRUE]], TRUE)
    evaluations_value <- scalar_or_null(
      old[["evaluations", exact = TRUE]], TRUE
    )
    evaluations <- if (is.null(evaluations_value)) NULL else
      list(gradient_count = evaluations_value)
    candidate_objective <- if (isTRUE(retain_candidates)) {
      scalar_or_null(old[["candidate_objective", exact = TRUE]])
    } else {
      NULL
    }
    if (!is.null(candidate_objective) && candidate_objective < 0) {
      candidate_objective <- NULL
    }
    elapsed <- scalar_or_null(old[["elapsed", exact = TRUE]])
    old_error <- old[["error", exact = TRUE]]
    error <- if (.dpprior_soft_plain_character(old_error)) {
      list(
        class = "optimizer_error",
        code = if (identical(old_error, "invalid_optimizer_result")) {
          "invalid_optimizer_result"
        } else {
          "optimizer_error"
        },
        message = old_error
      )
    } else {
      NULL
    }
    selected <- !is.null(selected_index) && identical(index, selected_index)
    reason <- if (selected) {
      "selected"
    } else if (!is.null(error)) {
      "optimizer_error"
    } else if (!is.null(exit_code) && exit_code != 0L) {
      "optimizer_exit_nonzero"
    } else if (!isTRUE(retain_candidates)) {
      "candidate_evaluation_failed"
    } else if (is.null(candidate_parameters)) {
      "nonfinite_candidate"
    } else if (is.null(candidate_objective)) {
      "nonfinite_objective"
    } else {
      "eligible_not_selected"
    }
    start <- old[["start", exact = TRUE]]
    if (!.dpprior_is_plain_numeric(start) || length(start) != 2L ||
        any(!is.finite(start))) start <- NULL
    bounds_old <- old[["bounds", exact = TRUE]]
    bounds <- if (is.list(bounds_old) && !is.object(bounds_old) &&
                  .dpprior_is_plain_numeric(bounds_old[["lower", exact = TRUE]]) &&
                  .dpprior_is_plain_numeric(bounds_old[["upper", exact = TRUE]])) {
      list(
        lower = unname(as.numeric(bounds_old[["lower", exact = TRUE]])),
        upper = unname(as.numeric(bounds_old[["upper", exact = TRUE]]))
      )
    } else {
      NULL
    }
    control <- old[["control", exact = TRUE]]
    if (!is.list(control) || is.object(control)) control <- NULL
    warnings <- old[["warning", exact = TRUE]]
    warnings <- if (.dpprior_soft_plain_character(warnings)) warnings else
      character()
    message <- old[["message", exact = TRUE]]
    if (!.dpprior_soft_plain_character(message)) message <- ""
    evidence <- list(
      start = start, bounds = bounds, control = control,
      exit_code = exit_code, iterations = iterations,
      evaluations = evaluations, candidate_parameters = candidate_parameters,
      candidate_objective = candidate_objective, elapsed_seconds = elapsed
    )
    missing <- names(evidence)[vapply(evidence, is.null, logical(1))]
    unavailable <- if (length(missing)) {
      setNames(
        paste("source attempt did not retain", gsub("_", " ", missing)),
        missing
      )
    } else {
      character()
    }
    .dpprior_new_attempt(
      id = sprintf("attempt-%03d", index),
      stage = if (identical(old[["method", exact = TRUE]], "Nelder-Mead")) {
        "fallback"
      } else {
        "primary"
      },
      method = old[["method", exact = TRUE]],
      start = start, bounds = bounds, control = control,
      exit_code = exit_code, message = message,
      iterations = iterations, evaluations = evaluations,
      candidate_parameters = candidate_parameters,
      candidate_objective = candidate_objective,
      elapsed_seconds = elapsed, warnings = warnings, error = error,
      selected = selected, reason_code = reason, unavailable = unavailable
    )
  })
}

.dpprior_soft_provenance <- function(fit_info, status, fallback_used,
                                     parameterization) {
  target_raw <- unclass(fit_info[["target_K", exact = TRUE]][[
    "canonical", exact = TRUE
  ]])
  approximation <- identical(status, "approximate")
  .dpprior_new_provenance(
    requested_method = "dual-soft",
    selected_method = "dual-soft",
    is_fallback = fallback_used,
    approximation = list(
      active = approximation,
      opt_in = FALSE,
      kind = if (approximation) "soft_candidate_not_decision_ready" else NULL,
      warning_code = if (approximation) "dual_soft_approximate" else NULL
    ),
    projection = target_raw[["provenance", exact = TRUE]][[
      "projection", exact = TRUE
    ]],
    parameterization = parameterization,
    backend = list(
      package = "DPprior",
      package_version = tryCatch(
        as.character(utils::packageVersion("DPprior")),
        error = function(condition) "development"
      ),
      implementation = "R/22_dual_anchor_soft.R:DPprior_dual_soft",
      source_commit = NULL
    ),
    input_fit = fit_info[["source", exact = TRUE]][[
      "canonical_reference", exact = TRUE
    ]],
    migration = list(
      source_schema = "native", adapter = "none", lossless = TRUE,
      missing_evidence = character(), warnings = character()
    ),
    legacy = list(
      active = FALSE, contract = NULL, deprecation_stage = NULL
    )
  )
}

.dpprior_soft_build_result <- function(fit, fit_info, target, lambda,
                                       candidate, verification, attempts,
                                       scales, log_bounds, controls,
                                       endpoint = FALSE,
                                       allow_approximate = FALSE,
                                       warm_start = NULL) {
  tolerances <- .dpprior_soft_canonical_tolerances(controls)
  input_reference <- fit_info[["source", exact = TRUE]][[
    "canonical_reference", exact = TRUE
  ]]
  parameterization <- input_reference[["parameters", exact = TRUE]][[
    "parameterization", exact = TRUE
  ]]
  parameters <- .dpprior_new_parameters(
    candidate$a, candidate$b, parameterization
  )
  selected <- verification$selected
  refined <- verification$verification
  selected_K_override <- verifier_K_override <- NULL
  if (isTRUE(endpoint)) {
    selected_K_override <- input_reference[["selected_snapshot", exact = TRUE]][[
      "achieved_K", exact = TRUE
    ]]
    verifier_K_override <- input_reference[["verifier_snapshot", exact = TRUE]][[
      "achieved_K", exact = TRUE
    ]]
  }
  selected_snapshot <- .dpprior_soft_snapshot(
    selected, parameters, fit_info, target, verification$selected_order,
    tolerances, "selected_order", K_override = selected_K_override
  )
  verifier_snapshot <- .dpprior_soft_snapshot(
    refined, parameters, fit_info, target, verification$verification_order,
    tolerances, "independent_verifier", K_override = verifier_K_override
  )
  stability <- .dpprior_soft_canonical_stability(
    selected_snapshot, verifier_snapshot, tolerances
  )
  losses <- .dpprior_soft_losses(
    list(
      K = list(
        mean = selected_snapshot$achieved$K$mean,
        variance = selected_snapshot$achieved$K$variance
      ),
      weight = list(value = selected_snapshot$achieved$weight$value)
    ), fit_info, target, lambda, scales
  )
  eta <- log(c(candidate$a, candidate$b))
  boundary_distance <- min(
    eta - log_bounds[1L], log_bounds[2L] - eta
  )
  selected_attempt <- if (length(attempts)) {
    which(vapply(
      attempts,
      function(x) identical(x$method, candidate$method) &&
        isTRUE(all.equal(
          as.numeric(x$candidate[c("shape", "rate")]),
          c(candidate$a, candidate$b), tolerance = 0
        )),
      logical(1)
    ))[1L]
  } else {
    NA_integer_
  }
  optimizer_passed <- isTRUE(endpoint) ||
    (length(selected_attempt) == 1L && !is.na(selected_attempt) &&
       identical(as.integer(attempts[[selected_attempt]]$exit_code), 0L))
  status <- if (!isTRUE(verification$performed) ||
                !isTRUE(verification$passed)) {
    "approximate"
  } else if (!isTRUE(optimizer_passed)) {
    "approximate"
  } else if (!isTRUE(endpoint) && boundary_distance <= controls$boundary_tol) {
    "boundary"
  } else {
    "converged"
  }
  stopifnot(status %in% .DPPRIOR_SOFT_STATUS)
  verified <- status %in% c("converged", "boundary")
  usable <- verified
  message <- if (isTRUE(endpoint) && isTRUE(verified)) {
      "lambda=1 returned the independently verified K-only fit without optimization"
    } else if (isTRUE(endpoint)) {
      paste(
        "lambda=1 retained the exact K-only input without optimization, but",
        "the independent weight/order stability contract did not pass"
      )
    } else if (identical(status, "converged")) {
      "soft trade-off optimizer and independent verification passed"
    } else if (identical(status, "boundary")) {
      "verified soft trade-off candidate is on the declared solver boundary"
    } else {
      "finite soft trade-off candidate did not pass the complete convergence contract"
    }
  canonical_controls <- .dpprior_soft_canonical_controls(
    controls, controls$max_iter, log_bounds
  )
  settings <- list(
    method = "dual-soft", controls = canonical_controls,
    parameterization = parameterization
  )
  scaling_values <- list(
    K = list(mean = unname(scales[["mean"]]),
             variance = unname(scales[["variance"]])),
    weight = 1
  )
  scaling <- .dpprior_new_scaling(
    requested = scaling_values, used = scaling_values,
    formula = "fixed_scaled_squared_loss", values = scaling_values,
    fixed_from_input = TRUE
  )
  canonical_attempts <- if (isTRUE(endpoint)) list() else
    .dpprior_soft_canonical_attempts(
      attempts, selected_attempt, parameterization
    )
  selected_attempt_id <- if (isTRUE(endpoint)) NULL else
    sprintf("attempt-%03d", selected_attempt)
  selected_candidate_id <- if (isTRUE(endpoint)) NULL else
    sprintf("candidate-%03d", selected_attempt)
  candidate_evaluations <- list()
  if (!isTRUE(endpoint)) {
    candidate_evaluations <- lapply(seq_along(canonical_attempts), function(i) {
      attempt <- canonical_attempts[[i]]
      candidate_parameters <- attempt[["candidate_parameters", exact = TRUE]]
      if (is.null(candidate_parameters)) return(NULL)
      evaluation <- .dpprior_soft_eval(
        candidate_parameters$a, candidate_parameters$b, fit_info, target,
        verification$selected_order
      )
      if (!isTRUE(evaluation$ok)) return(NULL)
      snapshot <- .dpprior_soft_snapshot(
        evaluation, candidate_parameters, fit_info, target,
        verification$selected_order, tolerances, "selected_order"
      )
      candidate_losses <- .dpprior_soft_losses(
        evaluation, fit_info, target, lambda, scales
      )
      fresh <- candidate_losses$total
      recorded <- attempt[["candidate_objective", exact = TRUE]]
      objective_tolerance <- controls$objective_abs_tol +
        controls$objective_rel_tol * max(abs(recorded), abs(fresh), 1)
      id <- sprintf("candidate-%03d", i)
      .dpprior_new_candidate_evaluation(
        id = id, attempt_id = attempt$id, method = attempt$method,
        generator = "direct_attempt", parameters = candidate_parameters,
        objective_kind = "soft_tradeoff",
        recorded_objective_kind = if (is.null(recorded)) NULL else
          "soft_tradeoff",
        selection_objective_kind = "soft_tradeoff",
        recorded_objective = recorded,
        recorded_objective_reason = if (is.null(recorded)) {
          "source optimizer did not retain a finite objective"
        } else {
          NULL
        },
        fresh_objective = fresh, selection_objective = fresh,
        objective_tolerance = objective_tolerance,
        selected_snapshot = snapshot, verifier_snapshot = NULL,
        checks = list(candidate_domain = .dpprior_new_check(
          value = c(
            K_support = evaluation$K$mean >= 1 &&
              evaluation$K$mean <= fit_info$J,
            K_variance = evaluation$K$variance >= 0,
            weight_support = evaluation$weight$value >= 0 &&
              evaluation$weight$value <= 1
          ),
          reference = c(
            K_support = TRUE, K_variance = TRUE, weight_support = TRUE
          ), tolerance = NULL, operator = "identical",
          source = paste0("candidate:", id)
        )),
        execution_success = identical(attempt$exit_code, 0L) &&
          is.null(attempt$error),
        optimizer_supported = identical(attempt$exit_code, 0L) &&
          is.null(attempt$error),
        diagnostic_eligible = FALSE,
        selected = identical(i, selected_attempt),
        source = "R/22_dual_anchor_soft.R:candidate_ledger"
      )
    })
    candidate_evaluations <- Filter(Negate(is.null), candidate_evaluations)
    for (evaluation in candidate_evaluations) {
      evaluation_raw <- unclass(evaluation)
      attempt_id <- evaluation_raw[["attempt_id", exact = TRUE]]
      attempt_index <- match(
        attempt_id,
        vapply(
          canonical_attempts,
          function(x) unclass(x)[["id", exact = TRUE]],
          character(1)
        )
      )
      if (is.na(attempt_index) ||
          isTRUE(evaluation_raw[["selected", exact = TRUE]])) {
        next
      }
      attempt <- unclass(canonical_attempts[[attempt_index]])
      attempt[["reason_code"]] <- if (isTRUE(
        evaluation_raw[["selection_eligible", exact = TRUE]]
      )) {
        "eligible_not_selected"
      } else if (identical(
        evaluation_raw[["objective_passed", exact = TRUE]], FALSE
      )) {
        "objective_recomputation_failed"
      } else {
        "candidate_eligibility_failed"
      }
      canonical_attempts[[attempt_index]] <- attempt
    }
  }
  fallback_attempted <- !isTRUE(endpoint) && length(canonical_attempts) > 1L
  fallback_used <- fallback_attempted && selected_attempt > 1L
  fallback <- .dpprior_new_fallback(
    attempted = fallback_attempted, used = fallback_used,
    trigger_attempt_id = if (fallback_attempted) "attempt-001" else NULL,
    selected_attempt_id = if (fallback_used) selected_attempt_id else NULL,
    reason_code = if (fallback_attempted) {
      "primary_attempt_not_accepted"
    } else {
      NULL
    },
    message = if (fallback_attempted) {
      "fallback attempted after the primary convergence contract did not pass"
    } else {
      ""
    },
    outcome = if (fallback_used) "selected" else if (fallback_attempted) {
      "attempted_not_selected"
    } else {
      "not_attempted"
    }
  )
  optimality <- if (isTRUE(endpoint)) {
    .dpprior_soft_unavailable_optimality(
      "lambda=1 endpoint performs no optimization"
    )
  } else {
    verification$optimality
  }
  components <- if (isTRUE(endpoint)) {
    endpoint_reference <- function(snapshot) {
      list(
        parameters = snapshot$parameters, M = snapshot$M,
        achieved_K = snapshot$achieved$K, finite = snapshot$finite,
        source = snapshot$source
      )
    }
    target_raw <- unclass(fit_info$target_K$canonical)
    endpoint_identity <- c(
      schema = identical(input_reference$schema, "dpprior.result/1"),
      mode = input_reference$mode %in%
        c("a2_moment", "a2_kl", "dual_hard", "dual_soft"),
      method = input_reference$method %in%
        .DPPRIOR_MODE_METHODS[[input_reference$mode]],
      J = identical(input_reference$J, fit_info$J),
      status = input_reference$status %in% c("converged", "boundary"),
      usable = isTRUE(input_reference$usable),
      verified = isTRUE(input_reference$verified),
      parameters = identical(input_reference$parameters, parameters),
      target = identical(input_reference$target, list(
        schema = target_raw$schema, kind = target_raw$kind, J = target_raw$J,
        used = target_raw$used, implied = target_raw$implied
      )),
      decision_evidence = if (identical(input_reference$mode, "a2_kl")) {
        !is.null(input_reference$decision_evidence) && identical(
          unclass(input_reference$decision_evidence$target_K), target_raw
        )
      } else {
        is.null(input_reference$decision_evidence)
      },
      selected_snapshot = identical(
        input_reference$selected_snapshot,
        endpoint_reference(selected_snapshot)
      ),
      verifier_snapshot = identical(
        input_reference$verifier_snapshot,
        endpoint_reference(verifier_snapshot)
      )
    )
    list(
      endpoint_input_identity = .dpprior_new_check(
        value = endpoint_identity,
        reference = setNames(rep(TRUE, length(endpoint_identity)),
                             names(endpoint_identity)),
        tolerance = NULL, operator = "identical",
        source = "independent_verifier"
      ),
      order_stability = .dpprior_new_check(
        value = stability$delta,
        reference = setNames(rep(0, length(stability$delta)),
                             names(stability$delta)),
        tolerance = stability$tolerance, operator = "lte",
        source = "independent_verifier"
      )
    )
  } else {
    selected_objectives <- vapply(
      candidate_evaluations,
      function(x) x[["selection_objective", exact = TRUE]], numeric(1)
    )
    selected_objective <- candidate_evaluations[[
      match(selected_candidate_id, vapply(
        candidate_evaluations, `[[`, character(1), "id"
      ))
    ]][["selection_objective", exact = TRUE]]
    list(
      objective_recomputation = .dpprior_new_check(
        value = optimality$objective_difference, reference = 0,
        tolerance = optimality$objective_tolerance, operator = "abs_lte",
        source = "independent_verifier"
      ),
      order_stability = .dpprior_new_check(
        value = stability$delta,
        reference = setNames(rep(0, length(stability$delta)),
                             names(stability$delta)),
        tolerance = stability$tolerance, operator = "lte",
        source = "independent_verifier"
      ),
      local_optimality = .dpprior_new_check(
        value = optimality$component_pass,
        reference = setNames(
          rep(TRUE, length(optimality$component_pass)),
          names(optimality$component_pass)
        ), tolerance = NULL, operator = "identical",
        source = "independent_verifier"
      ),
      candidate_selection = .dpprior_new_check(
        value = selected_objective - min(selected_objectives), reference = 0,
        tolerance = 0, operator = "lte", source = "independent_verifier"
      )
    )
  }
  invariants <- list(
    probability = .dpprior_new_check(
      value = selected_snapshot$achieved$weight$value >= 0 &&
        selected_snapshot$achieved$weight$value <= 1,
      reference = TRUE, tolerance = NULL, operator = "identical",
      source = "independent_verifier"
    ),
    K_support = .dpprior_new_check(
      value = selected_snapshot$achieved$K$mean >= 1 &&
        selected_snapshot$achieved$K$mean <= fit_info$J &&
        selected_snapshot$achieved$K$variance >= 0,
      reference = TRUE, tolerance = NULL, operator = "identical",
      source = "independent_verifier"
    ),
    finite_parameters_inside_domain = .dpprior_new_check(
      value = all(is.finite(eta)) && all(eta >= log_bounds[1L]) &&
        all(eta <= log_bounds[2L]),
      reference = TRUE, tolerance = NULL, operator = "identical",
      source = "independent_verifier"
    ),
    fixed_input_scales = .dpprior_new_check(
      value = TRUE, reference = TRUE, tolerance = NULL,
      operator = "identical", source = "independent_verifier"
    )
  )
  verification_passed <- verified &&
    all(vapply(components, function(x) isTRUE(x$passed), logical(1))) &&
    all(vapply(invariants, function(x) isTRUE(x$passed), logical(1))) &&
    isTRUE(stability$passed)
  if (isTRUE(endpoint) && !isTRUE(
    components[["endpoint_input_identity", exact = TRUE]][[
      "passed", exact = TRUE
    ]]
  )) {
    .dpprior_soft_abort_invalid(
      paste(
        "lambda=1 requires M_fit/M_verify and fresh snapshots identical to",
        "the canonical input-fit evidence"
      ),
      "M_fit/M_verify", c(M_fit = verification$selected_order,
                            M_verify = verification$verification_order),
      c(
        M_fit = input_reference$selected_snapshot$M,
        M_verify = input_reference$verifier_snapshot$M
      ), "dual_soft_endpoint_input_identity"
    )
  }
  canonical_verification <- .dpprior_new_verification(
    method = "fresh higher-order K moments and named weight metric",
    performed = TRUE, passed = verification_passed,
    reason = if (verification_passed) {
      "independent soft-tradeoff verification passed"
    } else {
      "independent soft-tradeoff verification did not establish decision readiness"
    },
    settings = list(
      M_fit = as.integer(verification$selected_order),
      M_verify = as.integer(verification$verification_order),
      log_bounds = unname(as.numeric(log_bounds)),
      verification_abs_tol = controls$verification_abs_tol,
      verification_rel_tol = controls$verification_rel_tol,
      objective_abs_tol = controls$objective_abs_tol,
      objective_rel_tol = controls$objective_rel_tol,
      boundary_tol = controls$boundary_tol,
      stationarity_step = controls$stationarity_step,
      stationarity_tol = controls$stationarity_tol,
      neighborhood_step = 1e-3
    ),
    selected_snapshot = selected_snapshot,
    verifier_snapshot = verifier_snapshot,
    stability = stability,
    components = components,
    invariants = invariants
  )
  if (!identical(verification_passed, verified)) {
    status <- "approximate"
    usable <- verified <- FALSE
    message <- paste(
      "finite soft trade-off candidate did not pass the complete",
      "independent convergence contract"
    )
  }
  computation <- .dpprior_new_computation(
    request = settings, used = settings,
    orders = .dpprior_new_orders(
      M_requested = as.integer(verification$selected_order),
      M_selected = as.integer(verification$selected_order),
      M_verification_required = .quadrature_verification_required_order(
        verification$selected_order
      ),
      M_verification_used = as.integer(verification$verification_order),
      requested_reason = "public_argument",
      selected_reason = "selected_order_recomputation",
      verification_required_reason = "independent_order_policy",
      verification_used_reason = "public_argument_or_required_default"
    ),
    scaling = scaling,
    attempts = canonical_attempts,
    candidate_evaluations = candidate_evaluations,
    selected_candidate_id = selected_candidate_id,
    selected_attempt_id = selected_attempt_id,
    fallback = fallback,
    termination = .dpprior_new_termination(
      code = if (isTRUE(endpoint)) "endpoint" else status,
      message = message,
      source = if (isTRUE(endpoint)) "endpoint" else if (
        identical(status, "approximate")
      ) {
        "candidate_evaluation"
      } else if (fallback_used) {
        "fallback_optimizer"
      } else {
        "optimizer"
      },
      iterations = if (isTRUE(endpoint)) NULL else
        canonical_attempts[[selected_attempt]]$iterations,
      boundary_reason = if (identical(status, "boundary")) {
        "declared_log_parameter_boundary"
      } else {
        NULL
      }
    ),
    resources = list(
      warm_start = warm_start,
      allow_approximate_return = isTRUE(allow_approximate),
      lambda_policy = "0 < lambda <= 1; 1-lambda must be representable",
      raw_optimizer_control = list(
        primary = controls$primary, fallback = controls$fallback
      )
    )
  )
  tradeoff <- .dpprior_new_tradeoff(
    lambda = lambda,
    K_loss = losses$K$value,
    weight_loss = losses$weight$value,
    total_loss = losses$total,
    target_residual = losses$weight$residual,
    directed_residual = losses$weight$directed_residual,
    scales = scaling_values,
    endpoint = isTRUE(endpoint), optimality = optimality
  )
  .dpprior_new_fit(
    mode = "dual_soft", method = "dual-soft", J = fit_info$J,
    status = status, usable = usable, verified = verified,
    message = message, parameters = parameters,
    target = list(
      K = fit_info$target_K$canonical,
      weight = .dpprior_soft_canonical_weight_target(target)
    ),
    achieved = selected_snapshot$achieved,
    residuals = selected_snapshot$residuals,
    tolerances = tolerances,
    computation = computation,
    verification = canonical_verification,
    provenance = .dpprior_soft_provenance(
      fit_info, status, fallback_used, parameterization
    ),
    extension = list(tradeoff = tradeoff)
  )
}

.dpprior_soft_failed_result <- function(fit_info, target, lambda, attempts,
                                        scales, log_bounds, controls,
                                        M_fit, M_verify, message) {
  input_reference <- fit_info[["source", exact = TRUE]][[
    "canonical_reference", exact = TRUE
  ]]
  parameterization <- input_reference[["parameters", exact = TRUE]][[
    "parameterization", exact = TRUE
  ]]
  canonical_controls <- .dpprior_soft_canonical_controls(
    controls, controls[["max_iter", exact = TRUE]], log_bounds
  )
  settings <- list(
    method = "dual-soft", controls = canonical_controls,
    parameterization = parameterization
  )
  tolerances <- .dpprior_soft_canonical_tolerances(controls)
  scaling_values <- list(
    K = list(
      mean = unname(scales[["mean"]]),
      variance = unname(scales[["variance"]])
    ),
    weight = 1
  )
  canonical_attempts <- .dpprior_soft_canonical_attempts(
    attempts, selected_index = NULL, parameterization = parameterization,
    retain_candidates = FALSE
  )
  fallback_attempted <- length(canonical_attempts) > 1L
  fallback <- .dpprior_new_fallback(
    attempted = fallback_attempted, used = FALSE,
    trigger_attempt_id = if (fallback_attempted) "attempt-001" else NULL,
    selected_attempt_id = NULL,
    reason_code = if (fallback_attempted) {
      "primary_attempt_not_accepted"
    } else {
      NULL
    },
    message = if (fallback_attempted) {
      "fallback was attempted but produced no public candidate"
    } else {
      ""
    },
    outcome = if (fallback_attempted) "failed" else "not_attempted"
  )
  computation <- .dpprior_new_computation(
    request = settings, used = settings,
    orders = .dpprior_new_orders(
      M_requested = NULL, M_selected = NULL,
      M_verification_required = NULL, M_verification_used = NULL,
      requested_reason = "not_applicable",
      selected_reason = "no_candidate",
      verification_required_reason = "not_applicable",
      verification_used_reason = "not_applicable"
    ),
    scaling = .dpprior_new_scaling(
      requested = scaling_values, used = scaling_values,
      formula = "fixed_scaled_squared_loss", values = scaling_values,
      fixed_from_input = TRUE
    ),
    attempts = canonical_attempts,
    candidate_evaluations = list(),
    selected_candidate_id = NULL, selected_attempt_id = NULL,
    fallback = fallback,
    termination = .dpprior_new_termination(
      code = "no_candidate", message = message,
      source = "no_candidate", iterations = NULL
    ),
    resources = list(
      requested_M = as.integer(M_fit),
      requested_M_verify = as.integer(M_verify),
      raw_candidate_count = as.integer(sum(vapply(
        attempts,
        function(x) !is.null(x[["candidate", exact = TRUE]]), logical(1)
      ))),
      allow_approximate_return = FALSE
    )
  )
  verification <- .dpprior_new_verification(
    method = "no_candidate", performed = FALSE, passed = FALSE,
    reason = "no candidate passed the public selection contract",
    settings = list(), selected_snapshot = NULL, verifier_snapshot = NULL,
    stability = NULL, components = list(),
    invariants = list(no_public_candidate = .dpprior_new_check(
      value = TRUE, reference = TRUE, tolerance = NULL,
      operator = "identical", source = "independent_verifier"
    ))
  )
  tradeoff <- .dpprior_new_tradeoff(
    lambda = lambda, K_loss = NULL, weight_loss = NULL, total_loss = NULL,
    target_residual = NULL, directed_residual = NULL,
    scales = scaling_values, endpoint = identical(lambda, 1),
    optimality = .dpprior_soft_unavailable_optimality(
      "no finite optimizer candidate passed the selection contract"
    )
  )
  .dpprior_new_fit(
    mode = "dual_soft", method = "dual-soft", J = fit_info$J,
    status = "failed", usable = FALSE, verified = FALSE,
    message = message, parameters = NULL,
    target = list(
      K = fit_info$target_K$canonical,
      weight = .dpprior_soft_canonical_weight_target(target)
    ),
    achieved = list(), residuals = list(), tolerances = tolerances,
    computation = computation, verification = verification,
    provenance = .dpprior_soft_provenance(
      fit_info, "failed", FALSE, parameterization
    ),
    extension = list(tradeoff = tradeoff)
  )
}

.dpprior_soft_signal_unusable <- function(result) {
  result <- .dpprior_require_schema(
    result, kind = "fit", allow_legacy = FALSE
  )
  raw <- unclass(result)
  status <- raw[["status", exact = TRUE]]
  message <- raw[["message", exact = TRUE]]
  subclass <- switch(
    status,
    approximate = "dpprior_dual_soft_approximation_error",
    infeasible = "dpprior_dual_soft_infeasible_error",
    failed = "dpprior_dual_soft_computation_error",
    "dpprior_dual_soft_computation_error"
  )
  action <- if (identical(status, "approximate")) {
    "review_then_explicitly_allow_approximate"
  } else {
    "refit_or_change_controls"
  }
  guidance <- if (identical(status, "approximate")) {
    paste(
      "Inspect condition$result; set allow_approximate=TRUE only after",
      "reviewing its residuals, attempts, and verification."
    )
  } else {
    paste(
      "Inspect condition$result and refit or revise the declared controls;",
      "allow_approximate does not return failed/no-candidate results."
    )
  }
  stop(.dpprior_new_condition(
    message = paste(message, guidance),
    classes = c(
      subclass, "dpprior_dual_soft_error", "dpprior_calibration_error",
      "dpprior_error", "error"
    ),
    result = result,
    status = status,
    action = action,
    code = paste0("dual_soft_", status)
  ))
}

#' Calibrate a fixed-scale soft Dual-Anchor trade-off
#'
#' @param fit A K-only \code{DPprior_fit} object.
#' @param target Named weight target with \code{metric}, \code{relation}, and
#'   \code{value}; threshold or probability is required when the metric needs
#'   it. Relations are \code{"target"}, \code{"at_most"}, or
#'   \code{"at_least"}.
#' @param lambda Mandatory finite scalar satisfying \code{0 < lambda <= 1}.
#' @param max_iter Maximum iterations for each optimizer attempt.
#' @param M_fit Selected quadrature order.
#' @param M_verify Independent verification order; by default at least twice
#'   \code{M_fit} and at least \code{M_fit + 40}.
#' @param log_bounds Two finite log-parameter bounds.
#' @param control Optimizer and verification controls.
#' @param allow_approximate Return a finite but unusable approximate candidate
#'   when \code{TRUE}; otherwise signal a typed condition retaining the result.
#' @param start Optional positive \code{c(a,b)} warm start.
#' @param ... Must be empty. Unknown or legacy hard-mode arguments are rejected
#'   as typed invalid input without evaluating their expressions.
#'
#' @return A canonical \code{dpprior.result/1} \code{DPprior_dual_soft}
#'   object. The \code{tradeoff} extension records the fixed input-derived K
#'   scales, lambda, weight target, component losses, and total objective.
#'   Soft objects never contain \code{constraint_satisfied}. Unapproved
#'   approximate, infeasible, or failed outcomes are signalled as typed
#'   conditions whose \code{condition$result} retains the complete object;
#'   \code{allow_approximate} changes return policy, not scientific status.
#'
#' @details A soft result is decision-ready only when its mode-specific
#'   contract and \code{status}, \code{usable}, and \code{verified} fields
#'   permit that use. Lambda is a trade-off weight, not a constraint
#'   probability or a hard-satisfaction certificate.
#'
#' @family elicitation
#' @export
DPprior_dual_soft <- function(
    fit, target, lambda, max_iter = 100L,
    M_fit = .QUAD_NODES_DEFAULT, M_verify = NULL,
    log_bounds = .LOG_BOUNDS_DEFAULT, control = list(),
    allow_approximate = FALSE, start = NULL, ...) {
  dots <- match.call(expand.dots = FALSE)$...
  if (!is.null(dots) && length(dots)) {
    dot_names <- names(dots)
    if (is.null(dot_names)) dot_names <- rep("<unnamed>", length(dots))
    dot_names[is.na(dot_names) | !nzchar(dot_names)] <- "<unnamed>"
    .dpprior_soft_abort_invalid(
      sprintf(
        "unknown DPprior_dual_soft argument(s): %s",
        paste(unique(dot_names), collapse = ", ")
      ),
      "...", dot_names, "no additional arguments",
      "dual_soft_unknown_argument"
    )
  }
  if (missing(target)) {
    .dpprior_soft_abort_invalid(
      "target is mandatory in soft mode", "target", NULL,
      "named weight target", "dual_soft_missing_target"
    )
  }
  if (missing(lambda)) {
    .dpprior_soft_abort_invalid(
      "lambda is mandatory in soft mode", "lambda", NULL,
      "plain finite scalar in (0,1]", "dual_soft_missing_lambda"
    )
  }
  fit_info <- .dpprior_soft_normalize_fit(fit)
  if (!.dpprior_is_plain_numeric(fit_info$J) || length(fit_info$J) != 1L ||
      !is.finite(fit_info$J) || fit_info$J < 2 ||
      fit_info$J != floor(fit_info$J)) {
    .dpprior_soft_abort_invalid(
      "soft dual-anchor calibration requires an identified design with J >= 2",
      "fit$J", fit_info$J, "integer J >= 2", "dual_soft_nonidentifiable_J"
    )
  }
  if (!isTRUE(fit_info$usable) || !isTRUE(fit_info$verified)) {
    .dpprior_soft_abort_invalid(
      paste(
        "decision-ready v2 soft calibration requires an input K-only fit",
        "with usable=TRUE and verified=TRUE before any optimization"
      ),
      "fit", fit, "verified usable K-only fit",
      "dual_soft_unverified_input_fit"
    )
  }
  target <- .dpprior_soft_normalize_target(target)
  lambda <- .dpprior_soft_scalar(
    lambda, "lambda", lower = 0, upper = 1, lower_open = TRUE
  )
  if (lambda < 1 && identical(1 - lambda, 1)) {
    .dpprior_soft_abort_invalid(
      paste(
        "lambda is positive mathematically but is not distinguishable from",
        "zero in the finite-double soft objective because 1-lambda equals 1"
      ),
      "lambda", lambda,
      "value in (0,1] with 1-lambda numerically distinct from 1",
      "dual_soft_lambda_not_representable"
    )
  }
  max_iter <- .dpprior_soft_integer(max_iter, "max_iter")
  if (max_iter > 100000L) {
    .dpprior_soft_abort_invalid(
      "max_iter must be no greater than 100000", "max_iter", max_iter,
      "integer in [1,100000]", "dual_soft_max_iter"
    )
  }
  M_fit <- .dpprior_validate_count(
    M_fit, "M_fit", minimum = 10L,
    maximum = as.integer(floor(.QUADRATURE_MAX_NODES / 2)),
    .subclass = "dpprior_dual_soft_control_error"
  )
  required <- .quadrature_verification_required_order(M_fit)
  if (is.null(M_verify)) {
    M_verify <- required
  } else {
    M_verify <- .dpprior_validate_count(
      M_verify, "M_verify", minimum = required,
      maximum = .QUADRATURE_MAX_NODES,
      .subclass = "dpprior_dual_soft_verification_error"
    )
  }
  if (!.dpprior_is_plain_numeric(log_bounds) || length(log_bounds) != 2L ||
      any(!is.finite(log_bounds)) || log_bounds[1L] >= log_bounds[2L] ||
      log_bounds[1L] < -.EXP_MAX || log_bounds[2L] > .EXP_MAX) {
    .dpprior_soft_abort_invalid(
      "log_bounds must be two increasing finite ordinary numeric values",
      "log_bounds", log_bounds, "finite c(lower, upper) with lower < upper",
      "dual_soft_log_bounds"
    )
  }
  allow_approximate <- .dpprior_validate_control(
    allow_approximate, "allow_approximate", type = "logical"
  )
  controls <- .dpprior_soft_control(control, max_iter)
  scales <- c(
    mean = max(abs(fit_info$target_K$mu_K), 1),
    variance = max(abs(fit_info$target_K$var_K), 1)
  )
  if (is.null(start)) {
    start_ab <- c(a = fit_info$a, b = fit_info$b)
    warm_start <- NULL
  } else {
    if (!.dpprior_is_plain_numeric(start) || length(start) != 2L ||
        any(!is.finite(start)) || any(start <= 0)) {
      .dpprior_soft_abort_invalid(
        "start must be two positive finite ordinary numeric values c(a,b)",
        "start", start, "positive finite c(a,b)", "dual_soft_start"
      )
    }
    start_ab <- stats::setNames(as.numeric(start), c("a", "b"))
    warm_start <- as.list(start_ab)
  }
  start_eta <- log(start_ab)
  if (any(start_eta < log_bounds[1L]) || any(start_eta > log_bounds[2L])) {
    .dpprior_soft_abort_invalid(
      "start parameters must lie inside log_bounds", "start", start,
      "positive c(a,b) inside exp(log_bounds)", "dual_soft_start_bounds"
    )
  }

  if (identical(lambda, 1)) {
    if (!is.null(start)) {
      .dpprior_soft_abort_invalid(
        "start must be omitted when lambda=1 because no optimizer is run",
        "start", start, "NULL for the exact K-only endpoint",
        "dual_soft_endpoint_irrelevant_start"
      )
    }
    endpoint_eta <- log(c(fit_info$a, fit_info$b))
    if (any(endpoint_eta < log_bounds[1L]) ||
        any(endpoint_eta > log_bounds[2L])) {
      .dpprior_soft_abort_invalid(
        "the lambda=1 input fit parameters lie outside log_bounds",
        "fit", c(a = fit_info$a, b = fit_info$b),
        "K-only fit parameters inside exp(log_bounds)",
        "dual_soft_endpoint_fit_bounds"
      )
    }
    candidate <- list(
      eta = log(c(fit_info$a, fit_info$b)), a = fit_info$a, b = fit_info$b,
      objective = NA_real_, exit_code = 0L,
      method = "K-only exact short-circuit"
    )
    verification <- .dpprior_soft_verify(
      candidate, fit_info, target, lambda, M_fit, M_verify, scales, controls,
      log_bounds, endpoint = TRUE
    )
    result <- .dpprior_soft_build_result(
      fit, fit_info, target, lambda, candidate, verification, list(), scales,
      log_bounds, controls, endpoint = TRUE,
      allow_approximate = allow_approximate, warm_start = NULL
    )
    result_raw <- unclass(result)
    if (!isTRUE(result_raw[["usable", exact = TRUE]]) &&
        !isTRUE(allow_approximate)) {
      .dpprior_soft_signal_unusable(result)
    }
    return(result)
  }

  objective <- .dpprior_soft_objective(
    fit_info, target, lambda, M_fit, scales, log_bounds
  )
  primary_call <- .dpprior_soft_call_optim(
    controls$optim_fun, "L-BFGS-B", start_eta, objective, log_bounds,
    controls$primary
  )
  attempts <- list(.dpprior_soft_attempt(
    primary_call, "L-BFGS-B", start_eta, log_bounds,
    controls$primary, objective
  ))
  primary_candidate <- .dpprior_soft_attempt_candidate(attempts[[1L]])
  primary_fresh_objective <- if (!is.null(primary_candidate)) {
    objective(primary_candidate$eta)
  } else NA_real_
  primary_objective_tolerance <- controls$objective_abs_tol +
    controls$objective_rel_tol * max(
      abs(primary_candidate$objective %||% NA_real_),
      abs(primary_fresh_objective), 1, na.rm = TRUE
    )
  primary_objective_bound <- !is.null(primary_candidate) &&
    is.finite(primary_candidate$objective) &&
    is.finite(primary_fresh_objective) &&
    abs(primary_candidate$objective - primary_fresh_objective) <=
      primary_objective_tolerance
  primary_passed <- !is.null(primary_candidate) &&
    identical(primary_candidate$exit_code, 0L) &&
    isTRUE(primary_objective_bound) &&
    primary_fresh_objective < .PENALTY_INF

  if (!primary_passed) {
    fallback_start <- if (!is.null(primary_candidate)) {
      pmin(pmax(primary_candidate$eta, log_bounds[1L]), log_bounds[2L])
    } else {
      start_eta
    }
    bounded_objective <- function(eta) {
      if (any(eta < log_bounds[1L]) || any(eta > log_bounds[2L])) {
        return(.PENALTY_INF + sum(pmax(log_bounds[1L] - eta, 0)^2) +
                 sum(pmax(eta - log_bounds[2L], 0)^2))
      }
      objective(eta)
    }
    fallback_call <- .dpprior_soft_call_optim(
      controls$optim_fun, "Nelder-Mead", fallback_start,
      bounded_objective, log_bounds, controls$fallback
    )
    attempts[[2L]] <- .dpprior_soft_attempt(
      fallback_call, "Nelder-Mead", fallback_start, log_bounds,
      controls$fallback, objective
    )
  }
  candidates <- Filter(Negate(is.null), lapply(
    attempts, .dpprior_soft_attempt_candidate
  ))
  candidates <- Filter(function(x) {
    all(x$eta >= log_bounds[1L]) && all(x$eta <= log_bounds[2L])
  }, candidates)
  candidates <- lapply(candidates, function(x) {
    x$fresh_objective <- objective(x$eta)
    x
  })
  candidates <- Filter(function(x) {
    if (!is.finite(x$objective) || !is.finite(x$fresh_objective) ||
        x$fresh_objective >= .PENALTY_INF) {
      return(FALSE)
    }
    objective_tolerance <- controls$objective_abs_tol +
      controls$objective_rel_tol * max(
        abs(x$objective), abs(x$fresh_objective), 1
      )
    abs(x$objective - x$fresh_objective) <= objective_tolerance
  }, candidates)
  if (!length(candidates)) {
    failed <- .dpprior_soft_failed_result(
      fit_info = fit_info, target = target, lambda = lambda,
      attempts = attempts, scales = scales, log_bounds = log_bounds,
      controls = controls, M_fit = M_fit, M_verify = M_verify,
      message = "all declared soft trade-off optimizer attempts failed"
    )
    .dpprior_soft_signal_unusable(failed)
  }
  candidate <- candidates[[which.min(vapply(
    candidates, `[[`, numeric(1), "fresh_objective"
  ))]]
  verification <- .dpprior_soft_verify(
    candidate, fit_info, target, lambda, M_fit, M_verify, scales, controls,
    log_bounds, endpoint = FALSE
  )
  result <- .dpprior_soft_build_result(
    fit, fit_info, target, lambda, candidate, verification, attempts, scales,
    log_bounds, controls, endpoint = FALSE,
    allow_approximate = allow_approximate, warm_start = warm_start
  )
  result_raw <- unclass(result)
  if (!isTRUE(result_raw[["usable", exact = TRUE]]) &&
      !isTRUE(allow_approximate)) {
    .dpprior_soft_signal_unusable(result)
  }
  result
}

.dpprior_soft_legacy_target <- function(w1_target, relation = "target") {
  if (!is.list(w1_target) || is.object(w1_target)) {
    .dpprior_soft_abort_invalid(
      "w1_target must be an ordinary list", "w1_target", w1_target,
      "one legacy prob, mean, or quantile target", "dual_soft_legacy_target"
    )
  }
  allowed <- c("prob", "mean", "quantile")
  if (is.null(names(w1_target)) || anyNA(names(w1_target)) ||
      any(!nzchar(names(w1_target))) || anyDuplicated(names(w1_target)) ||
      length(setdiff(names(w1_target), allowed))) {
    .dpprior_soft_abort_invalid(
      "w1_target must have unique names drawn from prob, mean, quantile",
      "w1_target", w1_target, paste(allowed, collapse = ", "),
      "dual_soft_legacy_target_names"
    )
  }
  prob_target <- w1_target[["prob", exact = TRUE]]
  mean_target <- w1_target[["mean", exact = TRUE]]
  quantile_target <- w1_target[["quantile", exact = TRUE]]
  present <- c(
    prob = !is.null(prob_target),
    mean = !is.null(mean_target),
    quantile = !is.null(quantile_target)
  )
  if (sum(present) != 1L) {
    .dpprior_soft_abort_invalid(
      "w1_target must specify exactly one of prob, mean, or quantile",
      "w1_target", w1_target, "exactly one legacy target",
      "dual_soft_legacy_target"
    )
  }
  if (present[["prob"]]) {
    prob_target <- .dpprior_v2_validate_named_list(
      prob_target, "w1_target$prob", allowed = c("threshold", "value"),
      required = c("threshold", "value")
    )
    return(list(
      metric = "wsb_tail", relation = relation,
      threshold = prob_target[["threshold", exact = TRUE]],
      value = prob_target[["value", exact = TRUE]]
    ))
  }
  if (present[["mean"]]) {
    return(list(
      metric = "wsb_mean", relation = relation,
      value = mean_target
    ))
  }
  quantile_target <- .dpprior_v2_validate_named_list(
    quantile_target, "w1_target$quantile", allowed = c("prob", "value"),
    required = c("prob", "value")
  )
  list(
    metric = "wsb_quantile", relation = relation,
    probability = quantile_target[["prob", exact = TRUE]],
    value = quantile_target[["value", exact = TRUE]]
  )
}

.dpprior_soft_format_number <- function(x) {
  bytes <- writeBin(as.double(x), raw(), size = 8L, endian = "big")
  paste(sprintf("%02x", as.integer(bytes)), collapse = "")
}

.dpprior_soft_point_id <- function(J, target_K, target, lambda) {
  fields <- c(
    "soft", paste0("J=", J),
    paste0("mu=", .dpprior_soft_format_number(target_K$mu_K)),
    paste0("var=", .dpprior_soft_format_number(target_K$var_K)),
    paste0("metric=", target$metric), paste0("relation=", target$relation),
    paste0("value=", .dpprior_soft_format_number(target$value)),
    paste0("threshold=", if (is.null(target$threshold)) "NA" else
      .dpprior_soft_format_number(target$threshold)),
    paste0("probability=", if (is.null(target$probability)) "NA" else
      .dpprior_soft_format_number(target$probability)),
    paste0("lambda=", .dpprior_soft_format_number(lambda))
  )
  paste(fields, collapse = "|")
}

.dpprior_soft_condition_fields <- function(condition) {
  if (is.null(condition)) {
    return(list(class = NA_character_, code = NA_character_, message = NA_character_))
  }
  condition_record <- if (is.list(condition)) unclass(condition) else NULL
  code_value <- if (is.list(condition_record)) {
    condition_record[["code", exact = TRUE]]
  } else {
    NULL
  }
  condition_classes <- class(condition)
  condition_class <- if (length(condition_classes) &&
                         .dpprior_soft_plain_character(
                           unname(condition_classes[[1L]])
                         )) {
    unname(condition_classes[[1L]])
  } else {
    "condition"
  }
  list(
    class = condition_class,
    code = if (.dpprior_soft_plain_character(code_value)) code_value
      else NA_character_,
    message = .dpprior_soft_safe_condition_message(condition)
  )
}

.dpprior_soft_same_scalar <- function(left, right, tolerance = 1e-12) {
  .dpprior_is_plain_numeric(left) && length(left) == 1L && is.finite(left) &&
    .dpprior_is_plain_numeric(right) && length(right) == 1L &&
    is.finite(right) && abs(left - right) <= tolerance
}

.dpprior_soft_curve_contract_violations <- function(
    point, J, target, lambda, base_info, M_fit, M_verify,
    log_bounds, controls) {
  point <- .dpprior_require_schema(
    point, kind = "fit", allow_legacy = FALSE
  )
  violations <- character()
  if (!identical(
    class(point),
    c("DPprior_dual_soft", "DPprior_fit", "dpprior_result", "list")
  )) {
    violations <- c(
      violations,
      "returned point must have the exact canonical dual-soft class vector"
    )
  }
  raw <- unclass(point)
  get <- function(name) raw[[name, exact = TRUE]]
  computation <- unclass(get("computation"))
  orders <- unclass(computation[["orders", exact = TRUE]])
  used <- unclass(computation[["used", exact = TRUE]])
  used_controls <- unclass(used[["controls", exact = TRUE]])
  resources <- unclass(computation[["resources", exact = TRUE]])
  tradeoff <- unclass(get("tradeoff"))
  provenance <- unclass(get("provenance"))
  target_bundle <- unclass(get("target"))
  status <- get("status")
  finite_status <- status %in% c("converged", "boundary", "approximate")

  if (!identical(get("mode"), "dual_soft") ||
      !identical(get("method"), "dual-soft")) {
    violations <- c(violations, "mode/method must be dual_soft/dual-soft")
  }
  if (!identical(get("J"), as.integer(J))) {
    violations <- c(
      violations, "J disagrees with the authoritative path request"
    )
  }
  if (!identical(tradeoff[["lambda", exact = TRUE]], lambda)) {
    violations <- c(
      violations, "lambda disagrees with the authoritative path request"
    )
  }
  if (!identical(
    target_bundle[["K", exact = TRUE]],
    base_info[["target_K", exact = TRUE]][["canonical", exact = TRUE]]
  )) {
    violations <- c(
      violations, "K target differs from the canonical input-fit target"
    )
  }
  if (!identical(
    target_bundle[["weight", exact = TRUE]],
    .dpprior_soft_canonical_weight_target(target)
  )) {
    violations <- c(
      violations, "weight target differs from the canonical path request"
    )
  }
  if (!identical(
    provenance[["input_fit", exact = TRUE]],
    base_info[["source", exact = TRUE]][["canonical_reference", exact = TRUE]]
  )) {
    violations <- c(
      violations, "provenance.input_fit differs from the canonical input fit"
    )
  }
  if (!identical(
    used_controls[["log_bounds", exact = TRUE]],
    unname(as.numeric(log_bounds))
  ) || !identical(
    used_controls[["max_iter", exact = TRUE]],
    as.integer(controls[["max_iter", exact = TRUE]])
  )) {
    violations <- c(
      violations, "used log_bounds/max_iter disagree with path controls"
    )
  }
  if (finite_status) {
    if (!identical(
      orders[["M_selected", exact = TRUE]], as.integer(M_fit)
    ) || !identical(
      orders[["M_verification_used", exact = TRUE]], as.integer(M_verify)
    ) || !identical(
      orders[["M_verification_required", exact = TRUE]],
      .quadrature_verification_required_order(M_fit)
    )) {
      violations <- c(
        violations, "quadrature orders disagree with path controls"
      )
    }
  } else if (identical(status, "failed")) {
    if (!identical(
      resources[["requested_M", exact = TRUE]], as.integer(M_fit)
    ) || !identical(
      resources[["requested_M_verify", exact = TRUE]], as.integer(M_verify)
    )) {
      violations <- c(
        violations, "failed point does not retain the requested orders"
      )
    }
  }
  unique(violations)
}
.dpprior_soft_curve_failed_result <- function(
    fit_info, target, lambda, message, M_fit, M_verify,
    log_bounds, controls) {
  scales <- c(
    mean = max(abs(fit_info[["target_K", exact = TRUE]][[
      "mu_K", exact = TRUE
    ]]), 1),
    variance = max(abs(fit_info[["target_K", exact = TRUE]][[
      "var_K", exact = TRUE
    ]]), 1)
  )
  .dpprior_soft_failed_result(
    fit_info = fit_info, target = target, lambda = lambda,
    attempts = list(), scales = scales, log_bounds = log_bounds,
    controls = controls, M_fit = M_fit, M_verify = M_verify,
    message = message
  )
}
.dpprior_soft_curve_row <- function(fit, condition, point_id, lambda,
                                    target, warm_start_from) {
  fit <- .dpprior_require_schema(
    fit, kind = "fit", allow_legacy = FALSE
  )
  raw <- unclass(fit)
  if (!identical(
    class(fit),
    c("DPprior_dual_soft", "DPprior_fit", "dpprior_result", "list")
  ) || !identical(raw[["mode", exact = TRUE]], "dual_soft")) {
    .dpprior_soft_abort_invalid(
      "curve rows require a canonical dual-soft result",
      "fit", class(fit), "canonical DPprior_dual_soft result",
      "dual_soft_curve_row_schema"
    )
  }
  fields <- .dpprior_soft_condition_fields(condition)
  parameters <- raw[["parameters", exact = TRUE]]
  parameters <- if (is.null(parameters)) NULL else unclass(parameters)
  achieved <- unclass(raw[["achieved", exact = TRUE]])
  achieved_K <- achieved[["K", exact = TRUE]]
  achieved_K <- if (is.null(achieved_K)) NULL else unclass(achieved_K)
  achieved_weight_record <- achieved[["weight", exact = TRUE]]
  achieved_weight_record <- if (is.null(achieved_weight_record)) {
    NULL
  } else {
    unclass(achieved_weight_record)
  }
  tradeoff <- unclass(raw[["tradeoff", exact = TRUE]])
  computation <- unclass(raw[["computation", exact = TRUE]])
  orders <- unclass(computation[["orders", exact = TRUE]])
  attempts <- computation[["attempts", exact = TRUE]]
  provenance <- unclass(raw[["provenance", exact = TRUE]])
  status <- raw[["status", exact = TRUE]]
  a <- if (is.null(parameters)) NA_real_ else
    parameters[["a", exact = TRUE]]
  b <- if (is.null(parameters)) NA_real_ else
    parameters[["b", exact = TRUE]]
  finite_parameters <- .dpprior_is_plain_numeric(a) && length(a) == 1L &&
    is.finite(a) && a > 0 &&
    .dpprior_is_plain_numeric(b) && length(b) == 1L &&
    is.finite(b) && b > 0
  selected_order <- orders[["M_selected", exact = TRUE]]
  wsb50 <- if (finite_parameters) {
    prob_wsb_exceeds(0.5, a, b)
  } else {
    NA_real_
  }
  wsb_mean <- if (finite_parameters) {
    mean_w1(a, b, selected_order %||% .QUAD_NODES_DEFAULT)
  } else {
    NA_real_
  }
  scalar_or_na <- function(value) {
    if (.dpprior_is_plain_numeric(value) && length(value) == 1L &&
        !is.na(value) && is.finite(value)) value else NA_real_
  }
  data.frame(
    point_id = point_id,
    lambda = lambda,
    mode = raw[["mode", exact = TRUE]],
    metric = target[["metric", exact = TRUE]],
    relation = target[["relation", exact = TRUE]],
    target_value = target[["value", exact = TRUE]],
    status = status,
    usable = isTRUE(raw[["usable", exact = TRUE]]),
    verified = isTRUE(raw[["verified", exact = TRUE]]),
    converged = identical(status, "converged"),
    outcome = if (is.null(condition)) "fit_returned" else
      "condition_retained",
    condition_class = fields[["class", exact = TRUE]],
    condition_code = fields[["code", exact = TRUE]],
    condition_message = fields[["message", exact = TRUE]],
    a = a,
    b = b,
    mu_K = scalar_or_na(if (is.null(achieved_K)) NULL else
      achieved_K[["mean", exact = TRUE]]),
    var_K = scalar_or_na(if (is.null(achieved_K)) NULL else
      achieved_K[["variance", exact = TRUE]]),
    achieved_weight = scalar_or_na(
      if (is.null(achieved_weight_record)) NULL else
        achieved_weight_record[["value", exact = TRUE]]
    ),
    target_residual = scalar_or_na(
      tradeoff[["target_residual", exact = TRUE]]
    ),
    K_loss = scalar_or_na(tradeoff[["K_loss", exact = TRUE]]),
    weight_loss = scalar_or_na(tradeoff[["weight_loss", exact = TRUE]]),
    total_loss = scalar_or_na(tradeoff[["total_loss", exact = TRUE]]),
    attempt_count = as.integer(length(attempts)),
    selected_method = provenance[["selected_method", exact = TRUE]],
    warm_start_from = warm_start_from %||% NA_character_,
    w_loss = scalar_or_na(tradeoff[["weight_loss", exact = TRUE]]),
    w1_prob_gt_50 = wsb50,
    E_w1 = wsb_mean,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}
#' Compute a status-aware soft trade-off path
#'
#' @param J Design size.
#' @param K_target List with \code{mu_K} and \code{var_K}.
#' @param w1_target Optional legacy W_SB target, translated to an explicit
#'   equality target.
#' @param lambda_seq Unique values satisfying \code{0 < lambda <= 1}.
#' @param max_iter,M,M_verify,log_bounds,control Soft solver controls.
#' @param verbose Show path progress.
#' @param loss_type Deprecated curve argument. Only \code{NULL} is accepted;
#'   the v2 curve always uses fixed input-derived scales.
#' @param target Preferred explicit named weight target.
#' @param allow_approximate Retain approximate points without making them
#'   usable warm starts.
#'
#' @return A deterministic \code{dpprior_tradeoff_curve} data frame with
#'   retained fit and condition lists as attributes. No lambda is selected as
#'   best.
#'
#' @export
compute_tradeoff_curve <- function(
    J, K_target, w1_target = NULL,
    lambda_seq = seq(0.1, 1, by = 0.1), max_iter = 100L,
    M = .QUAD_NODES_DEFAULT, verbose = FALSE, loss_type = NULL,
    target = NULL, M_verify = NULL, log_bounds = .LOG_BOUNDS_DEFAULT,
    control = list(), allow_approximate = TRUE) {
  if (!is.null(loss_type)) {
    .dpprior_soft_abort_invalid(
      paste(
        "loss_type is not a v2 trade-off-curve argument; the soft path uses",
        "fixed input-derived scales. Use DPprior_dual() only to reproduce",
        "legacy relative/adaptive/absolute equality-loss behavior."
      ),
      "loss_type", loss_type, "NULL", "dual_soft_legacy_loss_type"
    )
  }
  if (!is.null(target) && !is.null(w1_target)) {
    .dpprior_soft_abort_invalid(
      "supply exactly one of target or legacy w1_target",
      "target/w1_target", list(target = target, w1_target = w1_target),
      "one weight target", "dual_soft_conflicting_target"
    )
  }
  if (is.null(target)) {
    if (is.null(w1_target)) {
      .dpprior_soft_abort_invalid(
        "target is required", "target", NULL,
        "named soft weight target", "dual_soft_missing_target"
      )
    }
    target <- .dpprior_soft_legacy_target(w1_target, relation = "target")
  }
  target <- .dpprior_soft_normalize_target(target)
  target_for_fit <- list(
    metric = target[["metric", exact = TRUE]],
    relation = target[["relation", exact = TRUE]],
    value = target[["value", exact = TRUE]]
  )
  if (!is.null(target[["threshold", exact = TRUE]])) {
    target_for_fit[["threshold"]] <- target[["threshold", exact = TRUE]]
  }
  if (!is.null(target[["probability", exact = TRUE]])) {
    target_for_fit[["probability"]] <- target[["probability", exact = TRUE]]
  }
  assert_valid_J(J)
  if (J < 2) {
    .dpprior_soft_abort_invalid(
      "soft trade-off curves require an identified design with J >= 2",
      "J", J, "integer J >= 2", "dual_soft_nonidentifiable_J"
    )
  }
  K_target <- .dpprior_v2_validate_named_list(
    K_target, "K_target", allowed = c("mu_K", "var_K"),
    required = c("mu_K", "var_K")
  )
  mu_K <- .dpprior_soft_scalar(
    K_target[["mu_K", exact = TRUE]], "K_target$mu_K",
                               lower = 1, upper = J)
  var_K <- .dpprior_soft_scalar(
    K_target[["var_K", exact = TRUE]], "K_target$var_K", lower = 0
  )
  if (!.dpprior_is_plain_numeric(lambda_seq) || !length(lambda_seq) ||
      any(!is.finite(lambda_seq)) || any(lambda_seq <= 0) ||
      any(lambda_seq > 1)) {
    .dpprior_soft_abort_invalid(
      "lambda_seq must be a nonempty ordinary numeric vector in (0,1]",
      "lambda_seq", lambda_seq, "unique finite values in (0,1]",
      "dual_soft_lambda_path"
    )
  }
  lambda_seq <- as.numeric(lambda_seq)
  indistinguishable_zero <- vapply(
    lambda_seq,
    function(value) value < 1 && identical(1 - value, 1),
    logical(1)
  )
  if (any(indistinguishable_zero)) {
    .dpprior_soft_abort_invalid(
      "lambda_seq contains a positive value not representable in the soft objective",
      "lambda_seq", lambda_seq[indistinguishable_zero],
      "values whose 1-lambda is numerically distinct from 1",
      "dual_soft_lambda_not_representable"
    )
  }
  if (anyDuplicated(lambda_seq)) {
    .dpprior_soft_abort_invalid(
      "lambda_seq must not contain duplicate scientific points", "lambda_seq",
      lambda_seq, "unique values in (0,1]", "dual_soft_duplicate_lambda"
    )
  }
  verbose <- .dpprior_validate_control(verbose, "verbose", type = "logical")
  allow_approximate <- .dpprior_validate_control(
    allow_approximate, "allow_approximate", type = "logical"
  )
  max_iter <- .dpprior_soft_integer(max_iter, "max_iter")
  if (max_iter > 100000L) {
    .dpprior_soft_abort_invalid(
      "max_iter must be no greater than 100000", "max_iter", max_iter,
      "integer in [1,100000]", "dual_soft_max_iter"
    )
  }
  M <- .dpprior_validate_count(
    M, "M", minimum = 10L,
    maximum = as.integer(floor(.QUADRATURE_MAX_NODES / 2)),
    .subclass = "dpprior_dual_soft_control_error"
  )
  required_M_verify <- .quadrature_verification_required_order(M)
  if (is.null(M_verify)) {
    M_verify <- required_M_verify
  } else {
    M_verify <- .dpprior_validate_count(
      M_verify, "M_verify", minimum = required_M_verify,
      maximum = .QUADRATURE_MAX_NODES,
      .subclass = "dpprior_dual_soft_verification_error"
    )
  }
  if (!.dpprior_is_plain_numeric(log_bounds) || length(log_bounds) != 2L ||
      any(!is.finite(log_bounds)) || log_bounds[1L] >= log_bounds[2L] ||
      log_bounds[1L] < -.EXP_MAX || log_bounds[2L] > .EXP_MAX) {
    .dpprior_soft_abort_invalid(
      "log_bounds must be two increasing finite ordinary numeric values",
      "log_bounds", log_bounds,
      "finite c(lower, upper) with lower < upper", "dual_soft_log_bounds"
    )
  }
  control <- .dpprior_v2_validate_named_list(
    control, "control",
    allowed = c(
      "primary", "fallback", "boundary_tol", "verification_abs_tol",
      "verification_rel_tol", "objective_abs_tol", "objective_rel_tol",
      "stationarity_step", "stationarity_tol", ".optim_fun", ".fit_fun"
    )
  )
  fit_fun <- control[[".fit_fun", exact = TRUE]] %||% DPprior_dual_soft
  if (!is.function(fit_fun)) {
    .dpprior_soft_abort_invalid(
      "control$.fit_fun must be a function when supplied", "control$.fit_fun",
      fit_fun, "function", "dual_soft_curve_adapter"
    )
  }
  solver_control <- control
  solver_control[[".fit_fun"]] <- NULL
  curve_controls <- .dpprior_soft_control(solver_control, max_iter)
  base_fit <- tryCatch(
    DPprior_a2_newton(
      J = J, mu_K = mu_K, var_K = var_K,
      M = M, M_verify = M_verify, verbose = FALSE
    ),
    error = function(e) {
      condition_data <- if (is.list(e)) unclass(e) else NULL
      retained <- if (is.list(condition_data)) {
        condition_data[["result", exact = TRUE]]
      } else {
        NULL
      }
      if (!is.null(retained)) {
        return(.dpprior_require_schema(
          retained, kind = "fit", allow_legacy = FALSE
        ))
      }
      stop(e)
    }
  )
  base_info <- .dpprior_soft_normalize_fit(base_fit)
  if (!isTRUE(base_info$usable) || !isTRUE(base_info$verified)) {
    .dpprior_soft_abort_invalid(
      "K_target did not produce a verified usable K-only starting fit",
      "K_target", K_target, "verified usable K-only fit",
      "dual_soft_K_fit_unusable"
    )
  }

  evaluation_order <- sort(lambda_seq, decreasing = TRUE)
  fits <- list()
  conditions <- list()
  rows <- list()
  warm <- c(a = base_info$a, b = base_info$b)
  warm_id <- "K-only input fit"
  for (index in seq_along(evaluation_order)) {
    lambda <- evaluation_order[[index]]
    point_id <- .dpprior_soft_point_id(
      J, list(mu_K = mu_K, var_K = var_K), target, lambda
    )
    if (isTRUE(verbose)) {
      message(sprintf("soft trade-off lambda=%g", lambda))
    }
    condition <- NULL
    fit_arguments <- list(
        fit = base_fit, target = target_for_fit, lambda = lambda,
        max_iter = max_iter, M_fit = M, M_verify = M_verify,
        log_bounds = log_bounds, control = solver_control,
        allow_approximate = allow_approximate
      )
    if (!identical(lambda, 1)) fit_arguments$start <- warm
    point <- tryCatch(
      do.call(fit_fun, fit_arguments),
      error = function(e) {
        condition <<- e
        condition_data <- if (is.list(e)) unclass(e) else NULL
        if (is.list(condition_data)) {
          condition_data[["result", exact = TRUE]] %||% NULL
        } else {
          NULL
        }
      }
    )
    if (is.null(point)) {
      failure_message <- if (is.null(condition)) {
        "soft path returned NULL"
      } else {
        .dpprior_soft_safe_condition_message(condition)
      }
      upstream_condition <- condition
      point <- .dpprior_soft_curve_failed_result(
        base_info, target, lambda, failure_message,
        M, M_verify, log_bounds, curve_controls
      )
      condition <- .dpprior_new_condition(
        paste("soft path backend failed:", failure_message),
        c(
          "dpprior_dual_soft_backend_contract_error",
          "dpprior_dual_soft_error", "dpprior_calibration_error",
          "dpprior_error", "error"
        ),
        result = point,
        upstream_condition = upstream_condition,
        violations = "backend returned no canonical result",
        code = "dual_soft_backend_contract"
      )
    } else {
      validator_condition <- NULL
      violations <- tryCatch(
        .dpprior_soft_curve_contract_violations(
          point, J, target, lambda, base_info, M, M_verify,
          log_bounds, curve_controls
        ),
        error = function(e) {
          validator_condition <<- e
          paste(
            "result validator failed closed:",
            .dpprior_soft_safe_condition_message(e)
          )
        }
      )
      if (length(violations)) {
        upstream_condition <- condition %||% validator_condition
        failure_message <- paste(
          "soft path adapter violated the result contract:",
          paste(violations, collapse = "; ")
        )
        point <- .dpprior_soft_curve_failed_result(
          base_info, target, lambda, failure_message,
          M, M_verify, log_bounds, curve_controls
        )
        condition <- .dpprior_new_condition(
          failure_message,
          c(
            "dpprior_dual_soft_backend_contract_error",
            "dpprior_dual_soft_error", "dpprior_calibration_error",
            "dpprior_error", "error"
          ),
          result = point,
          upstream_condition = upstream_condition,
          violations = violations,
          code = "dual_soft_backend_contract"
        )
      } else if (!is.null(condition)) {
        condition_data <- if (is.list(condition)) unclass(condition) else NULL
        retained <- if (is.list(condition_data)) {
          condition_data[["result", exact = TRUE]]
        } else {
          NULL
        }
        if (!identical(retained, point)) {
          upstream_condition <- condition
          condition <- .dpprior_new_condition(
            "soft path condition did not retain its canonical point result",
            c(
              "dpprior_dual_soft_backend_contract_error",
              "dpprior_dual_soft_error", "dpprior_calibration_error",
              "dpprior_error", "error"
            ),
            result = point,
            upstream_condition = upstream_condition,
            violations = "condition$result identity mismatch",
            code = "dual_soft_backend_contract"
          )
        }
      }
    }
    fits[[point_id]] <- point
    conditions[point_id] <- list(condition)
    rows[[point_id]] <- .dpprior_soft_curve_row(
      point, condition, point_id, lambda, target, warm_id
    )
    point <- .dpprior_require_schema(
      point, kind = "fit", allow_legacy = FALSE
    )
    point_record <- unclass(point)
    point_parameters <- point_record[["parameters", exact = TRUE]]
    point_parameters <- if (is.null(point_parameters)) NULL else
      unclass(point_parameters)
    point_usable <- point_record[["usable", exact = TRUE]]
    point_verified <- point_record[["verified", exact = TRUE]]
    point_a <- if (is.null(point_parameters)) NULL else
      point_parameters[["a", exact = TRUE]]
    point_b <- if (is.null(point_parameters)) NULL else
      point_parameters[["b", exact = TRUE]]
    if (isTRUE(point_usable) && isTRUE(point_verified) &&
        .dpprior_soft_same_scalar(point_a, point_a) && point_a > 0 &&
        .dpprior_soft_same_scalar(point_b, point_b) && point_b > 0) {
      warm <- c(a = point_a, b = point_b)
      warm_id <- point_id
    }
  }
  curve <- do.call(rbind, unname(rows))
  curve <- curve[order(curve$lambda, curve$point_id), , drop = FALSE]
  rownames(curve) <- NULL
  fits <- fits[curve$point_id]
  conditions <- conditions[curve$point_id]
  attr(curve, "fits") <- fits
  attr(curve, "conditions") <- conditions
  attr(curve, "metadata") <- list(
    mode = "soft_tradeoff",
    evaluation_order = evaluation_order,
    output_order = "ascending lambda then point_id",
    warm_start_policy = "descending lambda; verified usable predecessor only",
    failure_policy = "all requested points retained",
    selection_policy = "descriptive path; no best lambda is selected",
    target_K = list(mu_K = mu_K, var_K = var_K),
    target_weight = target,
    fixed_scales = list(
      K_mean = max(abs(mu_K), 1),
      K_variance = max(abs(var_K), 1), weight = 1
    )
  )
  class(curve) <- c("dpprior_tradeoff_curve", "data.frame")
  curve
}
