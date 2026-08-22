# =============================================================================
# Module 15: legacy Dual-Anchor compatibility adapter
# =============================================================================
#
# The historical function in this module is retained throughout v2.x, with no
# removal before v3.0 and a migration review. It reproduces the old
# equality-loss calculation, but its result is explicitly legacy,
# approximate, and unverified. New work uses DPprior_dual_soft() or
# DPprior_dual_hard().

#' Legacy Dual-Anchor loss
#'
#' @param a,b Positive Gamma shape and rate.
#' @param J Design size.
#' @param K_target Legacy list with \code{mu_K} and \code{var_K}.
#' @param w1_target Legacy W_SB target.
#' @param lambda Legacy mixing weight in \code{[0,1]}.
#' @param M Quadrature order.
#' @param loss_type Historical scaling rule.
#' @param L_K_scale,L_w_scale Historical adaptive scale factors.
#'
#' @return Historical scalar equality loss. This helper is not the v2 soft
#'   contract because adaptive scales are path-derived.
#'
#' @keywords internal
dual_anchor_loss <- function(
    a, b, J, K_target, w1_target, lambda,
    M = .QUAD_NODES_DEFAULT,
    loss_type = c("relative", "adaptive", "absolute"),
    L_K_scale = NULL, L_w_scale = NULL) {
  loss_type <- match.arg(loss_type)
  if (!.dpprior_is_plain_numeric(c(a, b)) || length(c(a, b)) != 2L ||
      any(!is.finite(c(a, b))) || any(c(a, b) <= 0)) {
    return(.PENALTY_INF)
  }
  moments <- tryCatch(exact_K_moments(J, a, b, M), error = function(e) NULL)
  if (is.null(moments)) {
    return(.PENALTY_INF)
  }
  if (identical(loss_type, "absolute")) {
    L_K <- (moments$mean - K_target$mu_K)^2 +
      (moments$var - K_target$var_K)^2
  } else {
    L_K <- ((moments$mean - K_target$mu_K) / K_target$mu_K)^2 +
      ((moments$var - K_target$var_K) / K_target$var_K)^2
  }
  L_w <- tryCatch({
    if (!is.null(w1_target$quantile)) {
      (quantile_w1(w1_target$quantile$prob, a, b) -
         w1_target$quantile$value)^2
    } else if (!is.null(w1_target$prob)) {
      (prob_wsb_exceeds(w1_target$prob$threshold, a, b) -
         w1_target$prob$value)^2
    } else if (!is.null(w1_target$mean)) {
      (mean_w1(a, b, M) - w1_target$mean)^2
    } else {
      stop("missing legacy W_SB target", call. = FALSE)
    }
  }, error = function(e) .PENALTY_INF)
  if (!is.finite(L_w)) {
    return(.PENALTY_INF)
  }
  if (identical(loss_type, "adaptive")) {
    if (!is.null(L_K_scale) && is.finite(L_K_scale) && L_K_scale > 0) {
      L_K <- L_K / L_K_scale
    }
    if (!is.null(L_w_scale) && is.finite(L_w_scale) && L_w_scale > 0) {
      L_w <- L_w / L_w_scale
    }
  }
  lambda * L_K + (1 - lambda) * L_w
}

.dpprior_legacy_dual_abort <- function(message, argument, value, expected,
                                       code) {
  .dpprior_abort_invalid(
    message,
    c("dpprior_dual_legacy_invalid_input", "dpprior_dual_anchor_error"),
    argument, value, expected, code
  )
}

.dpprior_legacy_dual_input <- function(fit, argument = "fit") {
  if (!inherits(fit, "DPprior_fit")) {
    .dpprior_legacy_dual_abort(
      sprintf("%s must be a canonical DPprior_fit object", argument),
      argument, fit, "canonical K-only DPprior_fit", "legacy_dual_fit"
    )
  }
  validated <- .dpprior_require_schema(
    fit, kind = "fit", allow_legacy = FALSE
  )
  raw <- unclass(validated)
  mode <- raw[["mode", exact = TRUE]]
  if (!mode %in% c("a2_moment", "a2_kl") ||
      !raw[["status", exact = TRUE]] %in% c("converged", "boundary") ||
      !isTRUE(raw[["usable", exact = TRUE]]) ||
      !isTRUE(raw[["verified", exact = TRUE]])) {
    .dpprior_legacy_dual_abort(
      sprintf("%s must be a decision-ready K-only fit", argument),
      argument, fit,
      "verified usable a2_moment or a2_kl fit",
      "legacy_dual_K_only_contract"
    )
  }
  J <- raw[["J", exact = TRUE]]
  if (!is.integer(J) || length(J) != 1L || is.na(J) || J < 2L) {
    .dpprior_legacy_dual_abort(
      sprintf("%s$J must be an integer at least 2", argument),
      argument, J, "integer >= 2", "legacy_dual_J"
    )
  }
  parameters <- raw[["parameters", exact = TRUE]]
  target_bundle <- raw[["target", exact = TRUE]]
  target_K <- target_bundle[["K", exact = TRUE]]
  target_raw <- unclass(target_K)
  implied <- target_raw[["implied", exact = TRUE]]
  target_loss <- list(
    mu_K = implied[["mean", exact = TRUE]],
    var_K = implied[["variance", exact = TRUE]]
  )
  target_reference <- list(
    schema = target_raw[["schema", exact = TRUE]],
    kind = target_raw[["kind", exact = TRUE]],
    J = target_raw[["J", exact = TRUE]],
    used = target_raw[["used", exact = TRUE]],
    implied = implied
  )
  input_snapshot <- raw[["verification", exact = TRUE]][[
    "selected_snapshot", exact = TRUE
  ]]
  selected_snapshot_reference <- list(
    parameters = input_snapshot[["parameters", exact = TRUE]],
    M = input_snapshot[["M", exact = TRUE]],
    achieved_K = input_snapshot[["achieved", exact = TRUE]][[
      "K", exact = TRUE
    ]],
    finite = input_snapshot[["finite", exact = TRUE]],
    source = input_snapshot[["source", exact = TRUE]]
  )
  baseline <- list(
    schema = "dpprior.result/1",
    mode = mode,
    method = raw[["method", exact = TRUE]],
    J = J,
    status = raw[["status", exact = TRUE]],
    usable = raw[["usable", exact = TRUE]],
    verified = raw[["verified", exact = TRUE]],
    parameters = parameters,
    target = target_reference,
    selected_snapshot = selected_snapshot_reference
  )
  list(
    a = parameters[["a", exact = TRUE]],
    b = parameters[["b", exact = TRUE]],
    J = J,
    parameterization = parameters[["parameterization", exact = TRUE]],
    target_K = target_K,
    target_loss = target_loss,
    target_reference = target_reference,
    projection = target_raw[["provenance", exact = TRUE]][[
      "projection", exact = TRUE
    ]],
    baseline = baseline
  )
}

.dpprior_legacy_dual_validate <- function(fit, w1_target, lambda,
                                          max_iter, M, verbose, loss_type) {
  allowed_loss_types <- c("relative", "adaptive", "absolute")
  if (!.dpprior_schema_is_scalar_character(loss_type) ||
      !loss_type %in% allowed_loss_types) {
    .dpprior_abort_invalid(
      "loss_type must be one plain scalar legacy loss type",
      c(
        "dpprior_dual_legacy_loss_type_error",
        "dpprior_dual_legacy_invalid_input",
        "dpprior_dual_anchor_error"
      ),
      "loss_type", loss_type,
      "one of relative, adaptive, or absolute",
      "legacy_dual_loss_type"
    )
  }
  lambda <- as.numeric(.dpprior_validate_scalar(
    lambda, "lambda", lower = 0, upper = 1,
    .subclass = "dpprior_dual_legacy_lambda_error"
  ))
  max_iter <- .dpprior_validate_count(
    max_iter, "max_iter", minimum = 1L,
    maximum = floor(.Machine$integer.max / 2),
    .subclass = "dpprior_dual_legacy_control_error"
  )
  M <- .dpprior_validate_count(
    M, "M", minimum = 1L, maximum = .QUADRATURE_MAX_NODES,
    .subclass = "dpprior_dual_legacy_control_error"
  )
  verbose <- .dpprior_validate_control(verbose, "verbose", type = "logical")
  canonical <- .dpprior_soft_legacy_target(w1_target, relation = "target")
  normalized <- .dpprior_soft_normalize_target(canonical)
  safe_legacy_target <- switch(
    normalized[["metric", exact = TRUE]],
    wsb_tail = list(prob = list(
      threshold = normalized[["threshold", exact = TRUE]],
      value = normalized[["value", exact = TRUE]]
    )),
    wsb_mean = list(mean = normalized[["value", exact = TRUE]]),
    wsb_quantile = list(quantile = list(
      prob = normalized[["probability", exact = TRUE]],
      value = normalized[["value", exact = TRUE]]
    ))
  )
  fit_info <- .dpprior_legacy_dual_input(fit)
  list(
    lambda = lambda,
    loss_type = loss_type,
    max_iter = max_iter,
    M = M,
    verbose = verbose,
    canonical_target = canonical,
    normalized_target = normalized,
    legacy_target = safe_legacy_target,
    fit_info = fit_info
  )
}

.dpprior_legacy_dual_warning <- function() {
  .dpprior_warn(
    paste(
      "DPprior_dual() is deprecated and reproduces the legacy soft equality",
      "loss with historical scaling. It is not a hard constraint and does",
      "not certify target satisfaction. Use DPprior_dual_soft() with a named",
      "metric/relation or DPprior_dual_hard()."
    ),
    c("dpprior_dual_legacy_warning", "dpprior_deprecated_warning"),
    "DPprior_dual", "legacy equality-loss adapter",
    "DPprior_dual_soft() or DPprior_dual_hard()",
    "deprecated_legacy_dual_anchor"
  )
}

.dpprior_legacy_dual_wsb <- function(a, b, M) {
  list(
    mean = mean_w1(a, b, M),
    prob_gt_50 = prob_wsb_exceeds(0.5, a, b),
    prob_gt_90 = prob_wsb_exceeds(0.9, a, b)
  )
}

.dpprior_legacy_dual_weight_value <- function(a, b, target, M) {
  switch(
    target[["metric", exact = TRUE]],
    wsb_tail = prob_wsb_exceeds(
      target[["threshold", exact = TRUE]], a, b
    ),
    wsb_mean = mean_w1(a, b, M),
    wsb_quantile = quantile_w1(
      target[["probability", exact = TRUE]], a, b
    ),
    stop("unsupported legacy weight metric", call. = FALSE)
  )
}

.dpprior_legacy_dual_weight_target <- function(target, legacy_target) {
  authority <- list(
    metric = target[["metric", exact = TRUE]],
    relation = target[["relation", exact = TRUE]],
    value = target[["value", exact = TRUE]],
    threshold = target[["threshold", exact = TRUE]],
    probability = target[["probability", exact = TRUE]]
  )
  .dpprior_new_weight_target(
    request = authority,
    normalized = authority,
    used = authority,
    metric = target[["metric", exact = TRUE]],
    relation = target[["relation", exact = TRUE]],
    operator = target[["operator", exact = TRUE]],
    value = target[["value", exact = TRUE]],
    threshold = target[["threshold", exact = TRUE]],
    probability = target[["probability", exact = TRUE]],
    estimand = target[["estimand", exact = TRUE]],
    units = target[["units", exact = TRUE]],
    certification = list(),
    provenance = list(
      source = "R/15_dual_anchor.R:legacy_weight_target",
      transformation = NULL,
      selection = NULL,
      legacy_request = legacy_target
    )
  )
}

# A narrow indirection retained so fallback behavior can be tested without
# changing the public signature or the historical optimizer calls.
.dpprior_legacy_dual_optim <- function(...) stats::optim(...)

.dpprior_legacy_dual_optim_call <- function(...) {
  warnings <- character()
  value <- withCallingHandlers(
    tryCatch(.dpprior_legacy_dual_optim(...), error = identity),
    warning = function(condition) {
      warnings <<- c(
        warnings,
        .dpprior_soft_safe_condition_message(condition, "optimizer warning")
      )
      invokeRestart("muffleWarning")
    }
  )
  list(
    result = if (inherits(value, "condition")) NULL else value,
    condition = if (inherits(value, "condition")) value else NULL,
    warnings = warnings
  )
}

.dpprior_legacy_dual_optim_info <- function(call) {
  result <- call[["result", exact = TRUE]]
  ordinary <- if (is.list(result) && !is.object(result)) result else NULL
  result_names <- if (is.null(ordinary)) NULL else names(ordinary)
  schema_ok <- !is.null(ordinary) && !is.null(result_names) &&
    !anyNA(result_names) && all(nzchar(result_names)) &&
    !anyDuplicated(result_names) &&
    all(c("par", "value", "convergence") %in% result_names)
  par <- if (schema_ok) ordinary[["par", exact = TRUE]] else NULL
  value <- if (schema_ok) ordinary[["value", exact = TRUE]] else NULL
  convergence <- if (schema_ok) {
    ordinary[["convergence", exact = TRUE]]
  } else {
    NULL
  }
  valid <- schema_ok && .dpprior_is_plain_numeric(par) &&
    length(par) == 2L && all(is.finite(par)) &&
    .dpprior_is_plain_numeric(value) && length(value) == 1L &&
    is.finite(value) && .dpprior_is_plain_numeric(convergence) &&
    length(convergence) == 1L && is.finite(convergence) &&
    convergence == floor(convergence)
  if (valid) {
    natural_parameters <- exp(unname(par))
    valid <- all(is.finite(natural_parameters)) &&
      all(natural_parameters > 0)
  }
  counts <- if (schema_ok) ordinary[["counts", exact = TRUE]] else NULL
  counts <- if (.dpprior_is_plain_numeric(counts) &&
                !is.null(names(counts)) && !anyDuplicated(names(counts))) {
    counts
  } else {
    numeric()
  }
  optimizer_message <- if (schema_ok) {
    ordinary[["message", exact = TRUE]]
  } else {
    NULL
  }
  if (!is.character(optimizer_message) || length(optimizer_message) != 1L ||
      is.na(optimizer_message) || is.object(optimizer_message)) {
    optimizer_message <- ""
  }
  condition <- call[["condition", exact = TRUE]]
  error <- if (!is.null(condition)) {
    list(
      class = class(condition)[[1L]],
      code = "optimizer_error",
      message = .dpprior_soft_safe_condition_message(condition)
    )
  } else if (!valid) {
    list(
      class = "optimizer_result_error",
      code = "invalid_optimizer_result",
      message = "optimizer returned an invalid result schema"
    )
  } else {
    NULL
  }
  list(
    valid = valid,
    eta = if (valid) unname(as.numeric(par)) else NULL,
    a = if (valid) exp(unname(par[[1L]])) else NULL,
    b = if (valid) exp(unname(par[[2L]])) else NULL,
    value = if (valid) as.numeric(value) else NULL,
    exit_code = if (valid) as.integer(convergence) else NULL,
    message = optimizer_message,
    counts = counts,
    warnings = call[["warnings", exact = TRUE]],
    error = error
  )
}

.dpprior_legacy_dual_count <- function(x) {
  if (.dpprior_is_plain_numeric(x) && length(x) == 1L &&
      is.finite(x) && x >= 0 && x == floor(x) &&
      x <= .Machine$integer.max) {
    as.integer(x)
  } else {
    NULL
  }
}

.dpprior_legacy_dual_attempt <- function(id, stage, method, start, control,
                                         info, parameterization, selected,
                                         reason_code) {
  counts <- info[["counts", exact = TRUE]]
  iterations <- .dpprior_legacy_dual_count(
    if ("function" %in% names(counts)) counts[["function"]] else NULL
  )
  function_count <- iterations
  gradient_count <- .dpprior_legacy_dual_count(
    if ("gradient" %in% names(counts)) counts[["gradient"]] else NULL
  )
  evaluations <- list()
  if (!is.null(function_count)) evaluations$function_count <- function_count
  if (!is.null(gradient_count)) evaluations$gradient_count <- gradient_count
  if (!length(evaluations)) evaluations <- NULL
  candidate_parameters <- if (isTRUE(info[["valid", exact = TRUE]])) {
    .dpprior_new_parameters(
      info[["a", exact = TRUE]], info[["b", exact = TRUE]],
      parameterization
    )
  } else {
    NULL
  }
  evidence <- list(
    start = start,
    bounds = NULL,
    control = control,
    exit_code = info[["exit_code", exact = TRUE]],
    iterations = iterations,
    evaluations = evaluations,
    candidate_parameters = candidate_parameters,
    candidate_objective = info[["value", exact = TRUE]],
    elapsed_seconds = NULL
  )
  missing <- names(evidence)[vapply(evidence, is.null, logical(1))]
  unavailable <- if (length(missing)) {
    stats::setNames(
      paste("legacy optimizer did not retain", gsub("_", " ", missing)),
      missing
    )
  } else {
    character()
  }
  .dpprior_new_attempt(
    id = id,
    stage = stage,
    method = method,
    start = start,
    bounds = NULL,
    control = control,
    exit_code = info[["exit_code", exact = TRUE]],
    message = info[["message", exact = TRUE]],
    iterations = iterations,
    evaluations = evaluations,
    candidate_parameters = candidate_parameters,
    candidate_objective = info[["value", exact = TRUE]],
    elapsed_seconds = NULL,
    warnings = info[["warnings", exact = TRUE]],
    error = info[["error", exact = TRUE]],
    selected = selected,
    reason_code = reason_code,
    unavailable = unavailable
  )
}

.dpprior_legacy_dual_old_attempt <- function(method, info) {
  counts <- info[["counts", exact = TRUE]]
  count_names <- names(counts)
  keep <- !is.null(count_names) && !anyNA(count_names) &&
    all(nzchar(count_names))
  counts <- if (keep) {
    counts[is.finite(counts) & counts >= 0 & counts == floor(counts)]
  } else {
    numeric()
  }
  counts <- if (length(counts)) as.list(counts) else list()
  list(
    method = method,
    exit_code = info[["exit_code", exact = TRUE]],
    message = info[["message", exact = TRUE]],
    candidate_objective = info[["value", exact = TRUE]],
    candidate = if (isTRUE(info[["valid", exact = TRUE]])) {
      c(a = info[["a", exact = TRUE]], b = info[["b", exact = TRUE]])
    } else {
      NULL
    },
    counts = counts
  )
}

.dpprior_legacy_dual_result <- function(
    validated, loss_type, a, b, moments, K_loss, w_loss, total_loss,
    L_K_scale = NULL, L_w_scale = NULL, canonical_attempts = list(),
    old_attempts = list(), selected_attempt_id = NULL,
    selected_info = NULL, fallback_attempted = FALSE,
    fallback_used = FALSE, endpoint = FALSE) {
  fit_info <- validated[["fit_info", exact = TRUE]]
  target <- validated[["normalized_target", exact = TRUE]]
  lambda <- validated[["lambda", exact = TRUE]]
  M <- validated[["M", exact = TRUE]]
  parameterization <- fit_info[["parameterization", exact = TRUE]]
  parameters <- .dpprior_new_parameters(a, b, parameterization)
  weight_value <- .dpprior_legacy_dual_weight_value(a, b, target, M)
  achieved <- list(
    K = list(
      mean = moments[["mean", exact = TRUE]],
      variance = moments[["var", exact = TRUE]],
      estimand = "K_J",
      source = "legacy_selected_order",
      M = M
    ),
    weight = list(
      metric = target[["metric", exact = TRUE]],
      value = weight_value,
      source = "legacy_selected_order"
    )
  )
  residuals <- list(
    K = list(
      mean = moments[["mean", exact = TRUE]] -
        fit_info[["target_loss", exact = TRUE]][["mu_K", exact = TRUE]],
      variance = moments[["var", exact = TRUE]] -
        fit_info[["target_loss", exact = TRUE]][["var_K", exact = TRUE]]
    ),
    weight = list(
      raw = weight_value - target[["value", exact = TRUE]],
      directed = weight_value - target[["value", exact = TRUE]]
    )
  )
  tolerances <- list()
  selected_snapshot <- .dpprior_new_snapshot(
    parameters = parameters,
    M = M,
    achieved = achieved,
    residuals = residuals,
    tolerances = tolerances,
    finite = TRUE,
    source = "legacy_selected_order"
  )
  settings <- list(
    method = "dual-anchor",
    controls = list(
      lambda = lambda,
      loss_type = loss_type,
      max_iter = validated[["max_iter", exact = TRUE]],
      M = M
    ),
    parameterization = parameterization
  )
  scale_values <- if (identical(loss_type, "adaptive")) {
    list(L_K_scale = L_K_scale, L_w_scale = L_w_scale)
  } else {
    list()
  }
  scaling <- .dpprior_new_scaling(
    requested = list(loss_type = loss_type),
    used = list(loss_type = loss_type),
    formula = switch(
      loss_type,
      relative = "legacy_relative_squared_loss",
      adaptive = "legacy_path_derived_adaptive_loss",
      absolute = "legacy_absolute_squared_loss"
    ),
    values = scale_values,
    fixed_from_input = FALSE,
    change_reason = ""
  )
  fallback <- if (isTRUE(fallback_attempted)) {
    .dpprior_new_fallback(
      attempted = TRUE,
      used = fallback_used,
      trigger_attempt_id = "attempt-primary-001",
      selected_attempt_id = if (fallback_used) selected_attempt_id else NULL,
      reason_code = "primary_optimizer_exit_nonzero",
      message = "fallback was attempted after the primary optimizer",
      outcome = if (fallback_used) "selected" else
        "attempted_not_selected"
    )
  } else {
    .dpprior_new_fallback()
  }
  iterations <- if (isTRUE(endpoint)) {
    0L
  } else {
    selected_index <- match(
      selected_attempt_id,
      vapply(canonical_attempts, `[[`, character(1), "id")
    )
    canonical_attempts[[selected_index]][["iterations", exact = TRUE]]
  }
  termination <- if (isTRUE(endpoint)) {
    .dpprior_new_termination(
      code = "deterministic",
      message = "lambda=1 retained the K-only parameters",
      source = "legacy_adapter",
      iterations = 0L
    )
  } else {
    .dpprior_new_termination(
      code = if (identical(
        selected_info[["exit_code", exact = TRUE]], 0L
      )) "selected" else "approximate",
      message = paste(
        "legacy equality-loss candidate retained; optimizer exit is",
        "descriptive evidence only"
      ),
      source = if (fallback_used) "fallback_optimizer" else "optimizer",
      iterations = iterations
    )
  }
  computation <- .dpprior_new_computation(
    request = settings,
    used = settings,
    orders = .dpprior_new_orders(
      M_requested = M,
      M_selected = M,
      M_verification_required = NULL,
      M_verification_used = NULL,
      requested_reason = "public_argument",
      selected_reason = "legacy_selected_order",
      verification_required_reason =
        "legacy_adapter_did_not_request_independent_order",
      verification_used_reason =
        "legacy_adapter_did_not_perform_independent_order"
    ),
    scaling = scaling,
    attempts = canonical_attempts,
    candidate_evaluations = list(),
    selected_candidate_id = NULL,
    selected_attempt_id = selected_attempt_id,
    fallback = fallback,
    termination = termination,
    trace = NULL,
    resources = list(
      K_only_baseline = fit_info[["baseline", exact = TRUE]],
      optimizer_controls = list(
        adaptive = list(maxit = 50L),
        primary = list(maxit = validated[["max_iter", exact = TRUE]]),
        fallback = list(
          maxit = 2L * validated[["max_iter", exact = TRUE]]
        )
      ),
      legacy_weight_request = validated[["legacy_target", exact = TRUE]],
      objective = list(
        lambda = lambda,
        loss_type = loss_type,
        K_loss = K_loss,
        weight_loss = w_loss,
        total_loss = total_loss
      )
    )
  )
  verification <- .dpprior_new_verification(
    method = "legacy_unverified",
    performed = FALSE,
    passed = FALSE,
    reason = "legacy adapter retains no independent verification",
    settings = list(),
    selected_snapshot = selected_snapshot,
    verifier_snapshot = NULL,
    stability = NULL,
    components = list(),
    invariants = list()
  )
  provenance <- .dpprior_new_provenance(
    requested_method = "dual-anchor",
    selected_method = "dual-anchor",
    is_fallback = fallback_used,
    approximation = list(
      active = TRUE,
      opt_in = TRUE,
      kind = "legacy_soft_equality_loss",
      warning_code = "deprecated_legacy_dual_anchor"
    ),
    projection = fit_info[["projection", exact = TRUE]],
    parameterization = parameterization,
    backend = list(
      package = "DPprior",
      package_version = tryCatch(
        as.character(utils::packageVersion("DPprior")),
        error = function(condition) "development"
      ),
      implementation = "R/15_dual_anchor.R:DPprior_dual",
      source_commit = NULL
    ),
    input_fit = NULL,
    migration = list(
      source_schema = "native",
      adapter = "none",
      lossless = TRUE,
      missing_evidence = character(),
      warnings = character()
    ),
    legacy = list(
      active = TRUE,
      contract = "DPprior_dual_v1_path_scaled_soft_equality_loss",
      deprecation_stage = "compatibility_window"
    )
  )
  legacy <- .dpprior_new_legacy_details(
    contract = "DPprior_dual_v1_path_scaled_soft_equality_loss",
    lambda = lambda,
    losses = list(
      loss_type = loss_type,
      K_loss = K_loss,
      weight_loss = w_loss,
      total_loss = total_loss,
      scaling = scale_values
    ),
    approximation_opt_in = TRUE,
    warning_code = "deprecated_legacy_dual_anchor"
  )
  message <- if (isTRUE(endpoint)) {
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
  result <- .dpprior_new_fit(
    mode = "dual_legacy",
    method = "dual-anchor",
    J = fit_info[["J", exact = TRUE]],
    status = "approximate",
    usable = TRUE,
    verified = FALSE,
    message = message,
    parameters = parameters,
    target = list(
      K = fit_info[["target_K", exact = TRUE]],
      weight = .dpprior_legacy_dual_weight_target(
        target, validated[["legacy_target", exact = TRUE]]
      )
    ),
    achieved = achieved,
    residuals = residuals,
    tolerances = tolerances,
    computation = computation,
    verification = verification,
    provenance = provenance,
    extension = list(legacy = legacy)
  )
  achieved_fit <- list(
    mu_K = moments[["mean", exact = TRUE]],
    var_K = moments[["var", exact = TRUE]],
    residual = sqrt(
      residuals[["K", exact = TRUE]][["mean", exact = TRUE]]^2 +
        residuals[["K", exact = TRUE]][["variance", exact = TRUE]]^2
    )
  )
  selected_counts <- if (is.null(selected_info)) {
    list()
  } else {
    retained_counts <- selected_info[["counts", exact = TRUE]]
    count_names <- names(retained_counts)
    retained_counts <- if (!is.null(count_names) && !anyNA(count_names) &&
      all(nzchar(count_names))) {
      retained_counts[
        is.finite(retained_counts) & retained_counts >= 0 &
          retained_counts == floor(retained_counts)
      ]
    } else {
      numeric()
    }
    if (length(retained_counts)) as.list(retained_counts) else list()
  }
  dual_anchor <- list(
    w1_target = validated[["legacy_target", exact = TRUE]],
    canonical_target = validated[["canonical_target", exact = TRUE]],
    lambda = lambda,
    loss_type = loss_type,
    w1_achieved = .dpprior_legacy_dual_wsb(a, b, M),
    K_loss = K_loss,
    w_loss = w_loss,
    total_loss = total_loss,
    init = list(
      a = fit_info[["a", exact = TRUE]],
      b = fit_info[["b", exact = TRUE]]
    ),
    scaling = if (identical(loss_type, "adaptive")) scale_values else NULL,
    optim = list(
      convergence = if (isTRUE(endpoint)) 0L else
        selected_info[["exit_code", exact = TRUE]],
      value = if (isTRUE(endpoint)) NULL else
        selected_info[["value", exact = TRUE]],
      counts = selected_counts
    ),
    legacy_optimizer_converged = isTRUE(endpoint) || identical(
      selected_info[["exit_code", exact = TRUE]], 0L
    ),
    note = if (isTRUE(endpoint)) {
      "legacy lambda=1 K-only short-circuit"
    } else {
      "legacy equality-loss optimization"
    }
  )
  boundary <- .dpprior_compatibility_quarantine_boundary()
  dual_view <- c(dual_anchor, boundary)
  deprecation <- c(list(
    code = "legacy_dual_flat_view_quarantined",
    first_deprecated_version = "2.0.0",
    removal_floor = "not_scheduled"
  ), boundary)
  result <- .dpprior_append_compatibility_v2(
    result,
    aliases = c(
      a = "parameters.a",
      b = "parameters.b",
      converged = "compatibility.views.converged",
      iterations = "compatibility.views.iterations",
      fit = "compatibility.views.fit",
      attempts = "compatibility.views.attempts",
      dual_anchor = "compatibility.views.dual_anchor"
    ),
    views = list(
      legacy_dual_v2 = dual_view,
      converged = FALSE,
      iterations = iterations %||% 0L,
      fit = achieved_fit,
      attempts = old_attempts,
      dual_anchor = dual_view
    ),
    deprecations = list(legacy_dual_v2 = deprecation)
  )
  .dpprior_validate_result_v1(result)
  result
}

#' Legacy Dual-Anchor equality-loss adapter
#'
#' @description
#' \code{DPprior_dual()} is retained throughout v2.x, with no removal before
#' v3.0 and a migration review. It emits one typed lifecycle warning and
#' reproduces the historical equality-loss objective, including the old
#' path-derived adaptive scaling. Its result is always labelled legacy,
#' approximate, and unverified. It never represents a hard inequality or
#' target-satisfaction certificate.
#'
#' @param fit A K-only \code{DPprior_fit}.
#' @param w1_target Historical W_SB target using \code{prob}, \code{mean}, or
#'   \code{quantile} nesting.
#' @param lambda Historical mixing weight in \code{[0,1]}. The new soft API
#'   rejects zero; this compatibility adapter retains it only to reproduce old
#'   analyses.
#' @param max_iter Maximum iterations per historical optimizer.
#' @param verbose Print historical optimizer progress.
#' @param M Quadrature order.
#' @param loss_type Historical \code{"relative"}, \code{"adaptive"}, or
#'   \code{"absolute"} scaling.
#'
#' @return A canonical \code{dpprior.result/1} legacy \code{DPprior_fit} with
#'   \code{mode = "dual_legacy"}, \code{status = "approximate"},
#'   \code{usable = TRUE}, and \code{verified = FALSE}. The authoritative
#'   equality-loss, optimizer, scaling, and lifecycle records are
#'   \code{legacy}, \code{provenance}, and \code{computation}. The
#'   \code{dual_anchor} alias is a quarantined, non-authoritative compatibility
#'   view only. Use \code{\link{DPprior_dual_hard}} for a verified inequality
#'   or \code{\link{DPprior_dual_soft}} for a current fixed-scale trade-off.
#'
#' @family elicitation
#' @export
DPprior_dual <- function(
    fit, w1_target, lambda = 0.5, max_iter = 100L,
    verbose = FALSE, M = .QUAD_NODES_DEFAULT,
    loss_type = c("relative", "adaptive", "absolute")) {
  if (missing(loss_type)) loss_type <- "relative"
  validated <- .dpprior_legacy_dual_validate(
    fit, w1_target, lambda, max_iter, M, verbose, loss_type
  )
  loss_type <- validated[["loss_type", exact = TRUE]]
  if (!identical(loss_type, "absolute") &&
      validated[["fit_info", exact = TRUE]][[
        "target_loss", exact = TRUE
      ]][["var_K", exact = TRUE]] <= 0) {
    .dpprior_legacy_dual_abort(
      "relative legacy K loss requires a positive K variance target",
      "fit$target$K$implied$variance",
      validated[["fit_info", exact = TRUE]][[
        "target_loss", exact = TRUE
      ]][["var_K", exact = TRUE]],
      "positive variance for relative or adaptive loss",
      "legacy_dual_relative_variance"
    )
  }
  .dpprior_legacy_dual_warning()
  lambda <- validated[["lambda", exact = TRUE]]
  max_iter <- validated[["max_iter", exact = TRUE]]
  M <- validated[["M", exact = TRUE]]
  verbose <- validated[["verbose", exact = TRUE]]
  w1_target <- validated[["legacy_target", exact = TRUE]]
  fit_info <- validated[["fit_info", exact = TRUE]]
  J <- fit_info[["J", exact = TRUE]]
  K_target <- fit_info[["target_loss", exact = TRUE]]
  a_init <- fit_info[["a", exact = TRUE]]
  b_init <- fit_info[["b", exact = TRUE]]
  L_K_scale <- NULL
  L_w_scale <- NULL

  if (identical(lambda, 1)) {
    moments <- exact_K_moments(J, a_init, b_init, M)
    K_loss <- dual_anchor_loss(
      a_init, b_init, J, K_target, w1_target, 1, M, loss_type
    )
    w_loss <- dual_anchor_loss(
      a_init, b_init, J, K_target, w1_target, 0, M, loss_type
    )
    return(.dpprior_legacy_dual_result(
      validated = validated,
      loss_type = loss_type,
      a = a_init,
      b = b_init,
      moments = moments,
      K_loss = K_loss,
      w_loss = w_loss,
      total_loss = lambda * K_loss + (1 - lambda) * w_loss,
      canonical_attempts = list(),
      old_attempts = list(),
      selected_attempt_id = NULL,
      selected_info = NULL,
      fallback_attempted = FALSE,
      fallback_used = FALSE,
      endpoint = TRUE
    ))
  }

  scaling_attempt <- NULL
  scaling_info <- NULL
  if (identical(loss_type, "adaptive")) {
    L_w_init <- dual_anchor_loss(
      a_init, b_init, J, K_target, w1_target, 0, M, "relative"
    )
    w_only <- function(eta) {
      dual_anchor_loss(
        exp(eta[1L]), exp(eta[2L]), J, K_target, w1_target,
        0, M, "relative"
      )
    }
    scaling_call <- .dpprior_legacy_dual_optim_call(
      c(log(a_init), log(b_init)), w_only, method = "BFGS",
      control = list(maxit = 50L)
    )
    scaling_info <- .dpprior_legacy_dual_optim_info(scaling_call)
    a_w <- if (isTRUE(scaling_info[["valid", exact = TRUE]])) {
      scaling_info[["a", exact = TRUE]]
    } else {
      a_init
    }
    b_w <- if (isTRUE(scaling_info[["valid", exact = TRUE]])) {
      scaling_info[["b", exact = TRUE]]
    } else {
      b_init
    }
    L_K_at_w <- dual_anchor_loss(
      a_w, b_w, J, K_target, w1_target, 1, M, "relative"
    )
    L_K_scale <- max(L_K_at_w, 0.01)
    L_w_scale <- max(L_w_init, 0.001)
    if (isTRUE(verbose)) {
      cat(sprintf("Legacy adaptive L_K_scale = %.4e\n", L_K_scale))
      cat(sprintf("Legacy adaptive L_w_scale = %.4e\n", L_w_scale))
    }
    scaling_attempt <- .dpprior_legacy_dual_attempt(
      id = "attempt-scaling-001",
      stage = "scaling",
      method = "BFGS-weight-scaling",
      start = c(log_shape = log(a_init), log_rate = log(b_init)),
      control = list(maxit = 50L),
      info = scaling_info,
      parameterization = fit_info[["parameterization", exact = TRUE]],
      selected = FALSE,
      reason_code = if (isTRUE(scaling_info[["valid", exact = TRUE]])) {
        "scaling_evidence"
      } else {
        "optimizer_error"
      }
    )
  }

  objective <- function(eta) {
    eta <- pmin(pmax(eta, -.EXP_MAX), .EXP_MAX)
    dual_anchor_loss(
      exp(eta[1L]), exp(eta[2L]), J, K_target, w1_target,
      lambda, M, loss_type, L_K_scale, L_w_scale
    )
  }
  start <- c(log_shape = log(a_init), log_rate = log(b_init))
  primary_call <- .dpprior_legacy_dual_optim_call(
    start, objective, method = "BFGS", control = list(maxit = max_iter)
  )
  primary_info <- .dpprior_legacy_dual_optim_info(primary_call)
  fallback_info <- NULL
  fallback_attempted <- !isTRUE(primary_info[["valid", exact = TRUE]]) ||
    !identical(primary_info[["exit_code", exact = TRUE]], 0L)
  if (fallback_attempted) {
    fallback_call <- .dpprior_legacy_dual_optim_call(
      start, objective, method = "Nelder-Mead",
      control = list(maxit = 2L * max_iter)
    )
    fallback_info <- .dpprior_legacy_dual_optim_info(fallback_call)
  }
  fallback_used <- fallback_attempted &&
    isTRUE(fallback_info[["valid", exact = TRUE]]) &&
    (!isTRUE(primary_info[["valid", exact = TRUE]]) ||
       fallback_info[["value", exact = TRUE]] <
         primary_info[["value", exact = TRUE]])
  selected_info <- if (fallback_used) fallback_info else primary_info
  if (!isTRUE(selected_info[["valid", exact = TRUE]])) {
    .dpprior_legacy_dual_abort(
      "legacy optimizers did not retain a finite candidate",
      "optimizer", list(primary = primary_info, fallback = fallback_info),
      "finite optimizer candidate", "legacy_dual_optimizer_no_candidate"
    )
  }
  selected_attempt_id <- if (fallback_used) {
    "attempt-fallback-001"
  } else {
    "attempt-primary-001"
  }
  primary_attempt <- .dpprior_legacy_dual_attempt(
    id = "attempt-primary-001",
    stage = "primary",
    method = "BFGS",
    start = start,
    control = list(maxit = max_iter),
    info = primary_info,
    parameterization = fit_info[["parameterization", exact = TRUE]],
    selected = !fallback_used,
    reason_code = if (!fallback_used) {
      "selected"
    } else if (!is.null(primary_info[["error", exact = TRUE]])) {
      "optimizer_error"
    } else {
      "optimizer_exit_nonzero"
    }
  )
  canonical_attempts <- Filter(
    Negate(is.null), list(scaling_attempt, primary_attempt)
  )
  old_attempts <- list(
    primary = .dpprior_legacy_dual_old_attempt("BFGS", primary_info)
  )
  if (fallback_attempted) {
    fallback_attempt <- .dpprior_legacy_dual_attempt(
      id = "attempt-fallback-001",
      stage = "fallback",
      method = "Nelder-Mead",
      start = start,
      control = list(maxit = 2L * max_iter),
      info = fallback_info,
      parameterization = fit_info[["parameterization", exact = TRUE]],
      selected = fallback_used,
      reason_code = if (fallback_used) {
        "selected"
      } else if (!is.null(fallback_info[["error", exact = TRUE]])) {
        "optimizer_error"
      } else if (!identical(
        fallback_info[["exit_code", exact = TRUE]], 0L
      )) {
        "optimizer_exit_nonzero"
      } else {
        "eligible_not_selected"
      }
    )
    canonical_attempts[[length(canonical_attempts) + 1L]] <- fallback_attempt
    old_attempts[["fallback"]] <- .dpprior_legacy_dual_old_attempt(
      "Nelder-Mead", fallback_info
    )
  }
  a_opt <- selected_info[["a", exact = TRUE]]
  b_opt <- selected_info[["b", exact = TRUE]]
  moments <- exact_K_moments(J, a_opt, b_opt, M)
  K_loss <- dual_anchor_loss(
    a_opt, b_opt, J, K_target, w1_target, 1, M,
    loss_type, L_K_scale, L_w_scale
  )
  w_loss <- dual_anchor_loss(
    a_opt, b_opt, J, K_target, w1_target, 0, M,
    loss_type, L_K_scale, L_w_scale
  )
  .dpprior_legacy_dual_result(
    validated = validated,
    loss_type = loss_type,
    a = a_opt,
    b = b_opt,
    moments = moments,
    K_loss = K_loss,
    w_loss = w_loss,
    total_loss = selected_info[["value", exact = TRUE]],
    L_K_scale = L_K_scale,
    L_w_scale = L_w_scale,
    canonical_attempts = canonical_attempts,
    old_attempts = old_attempts,
    selected_attempt_id = selected_attempt_id,
    selected_info = selected_info,
    fallback_attempted = fallback_attempted,
    fallback_used = fallback_used,
    endpoint = FALSE
  )
}

#' Compare a Dual-Anchor fit with its K-only start
#'
#' @param fit_dual A legacy or v2 soft Dual-Anchor fit.
#' @param fit_K_only Optional K-only fit.
#' @param M Quadrature order.
#'
#' @return A comparison data frame with explicit W_SB labels.
#'
#' @export
dual_anchor_diagnostics <- function(
    fit_dual, fit_K_only = NULL, M = .QUAD_NODES_DEFAULT) {
  diagnostic_abort <- function(message, argument, value, expected, code) {
    .dpprior_abort_invalid(
      message,
      c("dpprior_dual_diagnostics_invalid_input", "dpprior_dual_anchor_error"),
      argument, value, expected, code
    )
  }
  require_fit <- function(value, argument) {
    if (!inherits(value, "DPprior_fit")) {
      diagnostic_abort(
        sprintf("%s must be a canonical DPprior_fit object", argument),
        argument, value, "canonical DPprior_fit",
        "dual_diagnostics_fit_class"
      )
    }
    tryCatch(
      .dpprior_require_schema(value, kind = "fit", allow_legacy = FALSE),
      error = function(condition) {
        diagnostic_abort(
          sprintf("%s failed canonical schema validation", argument),
          argument, value, "valid dpprior.result/1 fit",
          "dual_diagnostics_schema"
        )
      }
    )
  }

  fit_dual <- require_fit(fit_dual, "fit_dual")
  dual_raw <- unclass(fit_dual)
  mode <- dual_raw[["mode", exact = TRUE]]
  if (!mode %in% c("dual_legacy", "dual_soft")) {
    diagnostic_abort(
      "fit_dual must be a canonical legacy or soft Dual-Anchor fit",
      "fit_dual", fit_dual, "mode dual_legacy or dual_soft",
      "dual_diagnostics_mode_contract"
    )
  }
  parameters <- dual_raw[["parameters", exact = TRUE]]
  if (is.null(parameters)) {
    diagnostic_abort(
      "fit_dual must retain finite public parameters",
      "fit_dual", fit_dual, "finite canonical parameters",
      "dual_diagnostics_parameters"
    )
  }
  target_K <- dual_raw[["target", exact = TRUE]][["K", exact = TRUE]]
  target_raw <- unclass(target_K)
  implied <- target_raw[["implied", exact = TRUE]]
  expected_target_reference <- list(
    schema = target_raw[["schema", exact = TRUE]],
    kind = target_raw[["kind", exact = TRUE]],
    J = target_raw[["J", exact = TRUE]],
    used = target_raw[["used", exact = TRUE]],
    implied = implied
  )

  baseline <- if (identical(mode, "dual_legacy")) {
    dual_raw[["computation", exact = TRUE]][["resources", exact = TRUE]][[
      "K_only_baseline", exact = TRUE
    ]]
  } else {
    dual_raw[["provenance", exact = TRUE]][["input_fit", exact = TRUE]]
  }
  baseline_names <- if (is.list(baseline)) names(baseline) else NULL
  baseline_parameters <- if (is.list(baseline)) {
    baseline[["parameters", exact = TRUE]]
  } else {
    NULL
  }
  baseline_target <- if (is.list(baseline)) {
    baseline[["target", exact = TRUE]]
  } else {
    NULL
  }
  baseline_snapshot <- if (is.list(baseline)) {
    baseline[["selected_snapshot", exact = TRUE]]
  } else {
    NULL
  }
  snapshot_parameters <- if (is.list(baseline_snapshot)) {
    baseline_snapshot[["parameters", exact = TRUE]]
  } else {
    NULL
  }
  snapshot_M <- if (is.list(baseline_snapshot)) {
    baseline_snapshot[["M", exact = TRUE]]
  } else {
    NULL
  }
  snapshot_K <- if (is.list(baseline_snapshot)) {
    baseline_snapshot[["achieved_K", exact = TRUE]]
  } else {
    NULL
  }
  baseline_ok <- is.list(baseline) && !is.object(baseline) &&
    !is.null(baseline_names) && !anyNA(baseline_names) &&
    all(nzchar(baseline_names)) && !anyDuplicated(baseline_names) &&
    identical(baseline[["schema", exact = TRUE]], "dpprior.result/1") &&
    identical(baseline[["J", exact = TRUE]], dual_raw[["J", exact = TRUE]]) &&
    baseline[["status", exact = TRUE]] %in% c("converged", "boundary") &&
    isTRUE(baseline[["usable", exact = TRUE]]) &&
    isTRUE(baseline[["verified", exact = TRUE]]) &&
    is.list(baseline_parameters) && !is.object(baseline_parameters) &&
    identical(names(baseline_parameters), c("a", "b", "parameterization")) &&
    identical(baseline_target, expected_target_reference) &&
    is.list(baseline_snapshot) && !is.object(baseline_snapshot) &&
    identical(snapshot_parameters, baseline_parameters) &&
    .dpprior_is_plain_numeric(snapshot_M) && length(snapshot_M) == 1L &&
    is.finite(snapshot_M) && snapshot_M >= 1 && snapshot_M == floor(snapshot_M) &&
    is.list(snapshot_K) && !is.object(snapshot_K) &&
    all(c("mean", "variance") %in% names(snapshot_K))
  if (!isTRUE(baseline_ok)) {
    diagnostic_abort(
      "fit_dual lacks its exact decision-ready K-only baseline reference",
      "fit_dual", fit_dual, "ordinary canonical K-only baseline reference",
      "dual_diagnostics_baseline_reference"
    )
  }
  fresh_baseline <- tryCatch(
    exact_K_moments(
      dual_raw[["J", exact = TRUE]],
      baseline_parameters[["a", exact = TRUE]],
      baseline_parameters[["b", exact = TRUE]],
      as.integer(snapshot_M)
    ),
    error = function(condition) NULL
  )
  close <- function(recorded, recomputed) {
    .dpprior_is_plain_numeric(recorded) && length(recorded) == 1L &&
      is.finite(recorded) && is.finite(recomputed) &&
      abs(recorded - recomputed) <=
        1e-8 + 1e-8 * max(abs(recorded), abs(recomputed), 1)
  }
  if (is.null(fresh_baseline) ||
      !close(snapshot_K[["mean", exact = TRUE]], fresh_baseline$mean) ||
      !close(snapshot_K[["variance", exact = TRUE]], fresh_baseline$var)) {
    diagnostic_abort(
      "fit_dual K-only baseline failed fresh fixed-order recomputation",
      "fit_dual", fit_dual, "fresh identity with retained baseline snapshot",
      "dual_diagnostics_baseline_recomputation"
    )
  }

  if (!is.null(fit_K_only)) {
    explicit <- tryCatch(
      .dpprior_legacy_dual_input(fit_K_only, "fit_K_only"),
      error = function(condition) {
        diagnostic_abort(
          "fit_K_only failed the decision-ready K-only contract",
          "fit_K_only", fit_K_only, "decision-ready canonical K-only fit",
          "dual_diagnostics_K_contract"
        )
      }
    )
    if (!identical(
          explicit[["J", exact = TRUE]], dual_raw[["J", exact = TRUE]]
        ) || !identical(explicit[["target_K", exact = TRUE]], target_K)) {
      diagnostic_abort(
        "fit_K_only must retain the exact Dual-Anchor K target",
        "fit_K_only", fit_K_only, "exact J and canonical K-target identity",
        "dual_diagnostics_target_mismatch"
      )
    }
    if (!identical(
      explicit[["baseline", exact = TRUE]][["parameters", exact = TRUE]],
      baseline_parameters
    )) {
      diagnostic_abort(
        "fit_K_only parameters must equal the retained baseline parameters",
        "fit_K_only", fit_K_only, "exact retained baseline parameters",
        "dual_diagnostics_baseline_mismatch"
      )
    }
    baseline_parameters <- explicit[["baseline", exact = TRUE]][[
      "parameters", exact = TRUE
    ]]
  }

  M <- .dpprior_validate_count(
    M, "M", minimum = 1L, maximum = .QUADRATURE_MAX_NODES,
    .subclass = "dpprior_dual_diagnostics_invalid_input"
  )
  J <- dual_raw[["J", exact = TRUE]]
  a_K <- baseline_parameters[["a", exact = TRUE]]
  b_K <- baseline_parameters[["b", exact = TRUE]]
  a_dual <- parameters[["a", exact = TRUE]]
  b_dual <- parameters[["b", exact = TRUE]]
  moments_K <- exact_K_moments(J, a_K, b_K, M)
  moments_dual <- exact_K_moments(J, a_dual, b_dual, M)
  data.frame(
    Metric = c(
      "a", "b", "E[K_J]", "Var[K_J]", "P(W_SB > 0.5)",
      "P(W_SB > 0.9)", "E[W_SB]"
    ),
    K_only = c(
      a_K, b_K, moments_K$mean, moments_K$var,
      prob_wsb_exceeds(0.5, a_K, b_K),
      prob_wsb_exceeds(0.9, a_K, b_K),
      mean_w1(a_K, b_K, M)
    ),
    Dual_anchor = c(
      a_dual, b_dual, moments_dual$mean, moments_dual$var,
      prob_wsb_exceeds(0.5, a_dual, b_dual),
      prob_wsb_exceeds(0.9, a_dual, b_dual),
      mean_w1(a_dual, b_dual, M)
    ),
    stringsAsFactors = FALSE
  )
}

#' Verify the retained legacy adapter
#'
#' @param verbose Print a short result.
#'
#' @return Invisibly \code{TRUE}.
#'
#' @keywords internal
verify_dual_anchor <- function(verbose = TRUE) {
  fit <- DPprior_fit(
    J = 50, mu_K = 5, var_K = 8,
    method = "A2-MN", M = 80, check_diagnostics = FALSE
  )
  result <- suppressWarnings(DPprior_dual(
    fit,
    list(prob = list(threshold = 0.5, value = 0.3)),
    lambda = 0.5, M = 80, loss_type = "relative"
  ))
  .dpprior_validate_result_v1(result)
  raw <- unclass(result)
  input_raw <- unclass(fit)
  stopifnot(
    identical(raw[["mode", exact = TRUE]], "dual_legacy"),
    identical(raw[["method", exact = TRUE]], "dual-anchor"),
    identical(raw[["status", exact = TRUE]], "approximate"),
    isTRUE(raw[["usable", exact = TRUE]]),
    identical(raw[["verified", exact = TRUE]], FALSE),
    isTRUE(raw[["provenance", exact = TRUE]][[
      "legacy", exact = TRUE
    ]][["active", exact = TRUE]]),
    isTRUE(raw[["provenance", exact = TRUE]][[
      "approximation", exact = TRUE
    ]][["opt_in", exact = TRUE]]),
    is.null(raw[["provenance", exact = TRUE]][[
      "input_fit", exact = TRUE
    ]]),
    identical(
      raw[["target", exact = TRUE]][["K", exact = TRUE]],
      input_raw[["target", exact = TRUE]][["K", exact = TRUE]]
    ),
    identical(
      raw[["compatibility", exact = TRUE]][[
        "views", exact = TRUE
      ]][["legacy_dual_v2", exact = TRUE]][
        names(.dpprior_compatibility_quarantine_boundary())
      ],
      .dpprior_compatibility_quarantine_boundary()
    )
  )
  if (isTRUE(verbose)) {
    cat("Legacy Dual-Anchor adapter verification passed\n")
  }
  invisible(TRUE)
}
