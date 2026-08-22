# =============================================================================
# Module 11: A2-MN Newton Solver for Exact Moment Matching
# =============================================================================
#
# This module implements the A2-MN (exact-moment Newton) algorithm for
# calibrating Gamma hyperpriors on the Dirichlet process concentration
# parameter alpha to match target moments of K_J.
#
# Theory Background (Lee, 2026, Section 3.2):
# ------------------------------------
# Given target moments (mu_K, var_K), the A2-MN algorithm finds (a*, b*)
# such that the induced marginal distribution of K_J under alpha ~ Gamma(a*, b*)
# has:
#   E[K_J] = mu_K
#   Var(K_J) = var_K
#
# Algorithm:
# 1. Initialize from A1 closed-form: (a0, b0) <- DPprior_a1(J, mu_K, var_K)
# 2. Log-parameterize for positivity: eta = (log a, log b)
# 3. Newton iteration with backtracking line search
# 4. Return (a, b) = exp(eta)
#
# Key Features:
# - Uses score-based Jacobian (Module 07) for exact derivatives
# - Log-parameterization ensures positivity of (a, b)
# - Damped Newton with backtracking for global convergence
# - Optional Nelder-Mead fallback for difficult cases
#
# Author: JoonHo Lee (jlee296@ua.edu)
# Date: December 2025
# Part of: DPprior R Package
# Reference: Lee (2026), Section 3.2
# Dependencies: Modules 00, 02, 03, 05, 07, 10
# =============================================================================


# =============================================================================
# Main A2-MN Newton Solver
# =============================================================================

# Construct the componentwise residual contract used by both Newton and the
# fallback.  Natural-unit residuals remain available for interpretation, while
# scaling prevents the variance component from dominating merely because it is
# measured on a larger numerical scale.  Passing the aggregate norm never
# overrides a failed component.
.a2_residual_contract <- function(observed, target, abs_tol, rel_tol) {
  observed <- stats::setNames(as.numeric(observed), c("mean", "variance"))
  target <- stats::setNames(as.numeric(target), c("mean", "variance"))
  scale <- pmax(abs(target), 1)
  tolerance <- abs_tol + rel_tol * scale
  raw <- observed - target
  scaled <- raw / scale
  scaled_tolerance <- tolerance / scale
  standardized <- raw / tolerance
  component_pass <- is.finite(standardized) & abs(standardized) <= 1
  standardized_norm <- sqrt(mean(standardized^2))

  list(
    raw = raw,
    scale = stats::setNames(scale, names(raw)),
    scale_formula = "max(abs(target), 1)",
    scaled = scaled,
    tolerance = stats::setNames(tolerance, names(raw)),
    tolerance_formula = "absolute + relative * max(abs(target), 1)",
    scaled_tolerance = stats::setNames(scaled_tolerance, names(raw)),
    standardized = standardized,
    standardized_norm = standardized_norm,
    norm = "root-mean-square of residual/tolerance",
    component_pass = stats::setNames(component_pass, names(raw)),
    norm_pass = is.finite(standardized_norm) && standardized_norm <= 1,
    passed = all(component_pass) &&
      is.finite(standardized_norm) && standardized_norm <= 1
  )
}


# Fixed-order moment callback used by the optimizers.  The independent final
# verifier below deliberately does not reuse an optimizer-returned objective or
# cached moment vector.
.a2_moments_at_order <- function(J, a, b, M) {
  result <- exact_K_moments(J, a, b, M = M)
  c(mean = result$mean, variance = result$var)
}


.a2_required_verification_order <- function(M) {
  as.integer(max(2L * M, M + 40L))
}


# Independently recompute the selected candidate at the fit order and a
# contractually distinct higher order.  The fit-order candidate is retained by
# the caller; the higher-order values are evidence and never replace it.
.a2_verify_candidate <- function(
    J, a, b, M, M_verify, target, abs_tol, rel_tol,
    stability_abs_tol, stability_rel_tol) {
  required <- .a2_required_verification_order(M)
  available <- required <= .QUADRATURE_MAX_NODES

  selected <- tryCatch(
    .a2_moments_at_order(J, a, b, M),
    error = function(error) error
  )
  if (inherits(selected, "error")) {
    return(list(
      status = "failed", passed = FALSE, available = available,
      performed = FALSE, reason = "selected_order_recomputation_failed",
      message = conditionMessage(selected), M_selected = as.integer(M),
      M_required = required, M_verification = if (is.null(M_verify)) {
        NA_integer_
      } else {
        as.integer(M_verify)
      },
      selected = NULL, recomputed = NULL, residuals = NULL,
      stability = NULL
    ))
  }

  selected_contract <- .a2_residual_contract(
    selected, target, abs_tol, rel_tol
  )
  if (!available || is.null(M_verify) || is.na(M_verify)) {
    return(list(
      status = "approximate", passed = FALSE, available = FALSE,
      performed = FALSE, reason = "verification_order_exceeds_ceiling",
      message = sprintf(
        paste(
          "independent verification requires M_verify=%d, which exceeds",
          "the supported ceiling %d"
        ),
        required, .QUADRATURE_MAX_NODES
      ),
      M_selected = as.integer(M), M_required = required,
      M_verification = NA_integer_, selected = selected,
      recomputed = NULL,
      residuals = list(selected = selected_contract, verification = NULL),
      stability = NULL
    ))
  }

  recomputed <- tryCatch(
    .a2_moments_at_order(J, a, b, M_verify),
    error = function(error) error
  )
  if (inherits(recomputed, "error")) {
    return(list(
      status = "failed", passed = FALSE, available = TRUE,
      performed = TRUE, reason = "higher_order_recomputation_failed",
      message = conditionMessage(recomputed), M_selected = as.integer(M),
      M_required = required, M_verification = as.integer(M_verify),
      selected = selected, recomputed = NULL,
      residuals = list(selected = selected_contract, verification = NULL),
      stability = NULL
    ))
  }

  verification_contract <- .a2_residual_contract(
    recomputed, target, abs_tol, rel_tol
  )
  difference <- abs(selected - recomputed)
  stability_scale <- pmax(abs(selected), abs(recomputed), 1)
  stability_tolerance <- stability_abs_tol +
    stability_rel_tol * stability_scale
  stability_pass <- is.finite(difference) &
    difference <= stability_tolerance
  passed <- isTRUE(selected_contract$passed) &&
    isTRUE(verification_contract$passed) && all(stability_pass)
  reason <- if (!isTRUE(selected_contract$passed)) {
    "selected_residual_tolerance_not_met"
  } else if (!isTRUE(verification_contract$passed)) {
    "verification_residual_tolerance_not_met"
  } else if (!all(stability_pass)) {
    "quadrature_order_disagreement"
  } else {
    "independent_verification_passed"
  }

  list(
    status = if (passed) "converged" else "approximate",
    passed = passed, available = TRUE, performed = TRUE,
    reason = reason,
    message = if (passed) {
      "independent higher-order moment verification passed"
    } else {
      "independent higher-order moment verification did not meet every tolerance"
    },
    M_selected = as.integer(M), M_required = required,
    M_verification = as.integer(M_verify), selected = selected,
    recomputed = recomputed,
    residuals = list(
      selected = selected_contract,
      verification = verification_contract
    ),
    stability = list(
      absolute_difference = stats::setNames(difference, names(selected)),
      scale = stats::setNames(stability_scale, names(selected)),
      tolerance = stats::setNames(stability_tolerance, names(selected)),
      tolerance_formula = paste(
        "stability_absolute + stability_relative *",
        "max(abs(selected), abs(verification), 1)"
      ),
      component_pass = stats::setNames(stability_pass, names(selected)),
      passed = all(stability_pass)
    )
  )
}


# Kept as a small helper so failure-path tests can inject a malformed or failed
# fallback without replacing stats namespace bindings.  Parameters are always
# log(shape), log(rate); raw positive parameters are never optimized.
.a2_nelder_mead <- function(start, objective, maxit = 1000L) {
  stats::optim(
    par = start, fn = objective, method = "Nelder-Mead",
    control = list(maxit = maxit, reltol = 1e-12)
  )
}


# Fixed log-parameter starts are independent of the requested target.  They
# span diffuse/concentrated and low/high mean concentration regimes while
# remaining strictly inside the declared solver box.  Target residuals are
# used only to choose among these predeclared starts.
.a2_initialization_grid <- function() {
  parameters <- rbind(
    c(shape = 1, rate = 1),
    c(shape = 0.25, rate = 1),
    c(shape = 4, rate = 1),
    c(shape = 1, rate = 0.25),
    c(shape = 4, rate = 0.1),
    c(shape = 10, rate = 0.01),
    c(shape = 1, rate = 0.01),
    c(shape = 0.1, rate = 10),
    c(shape = 10, rate = 0.1),
    c(shape = 100, rate = 0.1)
  )
  log(parameters)
}


.a2_select_initialization_grid <- function(
    J, target, M, abs_tol, rel_tol) {
  grid <- .a2_initialization_grid()
  records <- vector("list", nrow(grid))
  objectives <- rep(Inf, nrow(grid))
  for (index in seq_len(nrow(grid))) {
    eta <- grid[index, ]
    evaluated <- tryCatch(
      .a2_moments_at_order(J, exp(eta[1L]), exp(eta[2L]), M),
      error = function(error) error
    )
    if (inherits(evaluated, "error")) {
      records[[index]] <- list(
        method = "fixed_log_parameter_grid", index = index,
        log_parameters = stats::setNames(
          as.numeric(eta), c("log_shape", "log_rate")
        ),
        parameters = stats::setNames(
          exp(as.numeric(eta)), c("shape", "rate")
        ),
        status = "failed", objective = NA_real_,
        error = conditionMessage(evaluated)
      )
      next
    }
    residual <- .a2_residual_contract(
      evaluated, target, abs_tol, rel_tol
    )
    objective <- sum(residual$standardized^2)
    if (is.finite(objective)) objectives[index] <- objective
    records[[index]] <- list(
      method = "fixed_log_parameter_grid", index = index,
      log_parameters = stats::setNames(
        as.numeric(eta), c("log_shape", "log_rate")
      ),
      parameters = stats::setNames(
        exp(as.numeric(eta)), c("shape", "rate")
      ),
      status = if (is.finite(objective)) "evaluated" else "failed",
      objective = if (is.finite(objective)) objective else NA_real_,
      error = if (is.finite(objective)) NA_character_ else "non-finite objective"
    )
  }
  if (!any(is.finite(objectives))) {
    .dpprior_abort_invalid(
      "every fixed A2-MN initialization candidate failed numerical evaluation",
      c("dpprior_a2_initialization_error", "dpprior_numerical_error"),
      "initialization_grid", records, "at least one finite candidate",
      "initialization_grid_failed"
    )
  }
  selected <- which.min(objectives)
  list(
    eta = stats::setNames(
      as.numeric(grid[selected, ]), c("log_shape", "log_rate")
    ),
    selected_index = selected,
    selected_objective = objectives[selected],
    selection_rule = paste(
      "minimum standardized residual objective among a fixed",
      "target-independent log-parameter grid"
    ),
    candidates = records
  )
}


.a2_trace_row <- function(
    iter, eta, observed = NULL, residual = NULL,
    step = NA_real_, step_norm = NA_real_, line_iterations = 0L,
    accepted = FALSE, determinant = NA_real_, reciprocal_condition = NA_real_,
    jacobian_status = NA_character_, derivative_status = NA_character_,
    reason = NA_character_) {
  if (is.null(observed)) observed <- c(mean = NA_real_, variance = NA_real_)
  if (is.null(residual)) {
    raw <- scaled <- c(mean = NA_real_, variance = NA_real_)
    residual_norm <- standardized_norm <- max_ratio <- NA_real_
  } else {
    raw <- residual$raw
    scaled <- residual$scaled
    residual_norm <- sqrt(sum(raw^2))
    standardized_norm <- residual$standardized_norm
    max_ratio <- max(abs(residual$standardized))
  }
  data.frame(
    iter = as.integer(iter), a = unname(exp(eta[1L])),
    b = unname(exp(eta[2L])),
    M1 = observed[["mean"]], V = observed[["variance"]],
    residual = residual_norm, residual_mean = raw[["mean"]],
    residual_variance = raw[["variance"]],
    scaled_mean = scaled[["mean"]],
    scaled_variance = scaled[["variance"]],
    standardized_norm = standardized_norm,
    max_budget_ratio = max_ratio, step = step, step_norm = step_norm,
    line_search_iterations = as.integer(line_iterations),
    accepted = isTRUE(accepted), det_Jlog = determinant,
    reciprocal_condition = reciprocal_condition,
    jacobian_status = jacobian_status,
    derivative_status = derivative_status,
    reason_code = reason, stringsAsFactors = FALSE
  )
}


.a2_scaled_jacobian_condition <- function(jacobian) {
  values <- tryCatch(
    svd(jacobian, nu = 0L, nv = 0L)$d,
    error = function(error) c(NA_real_, NA_real_)
  )
  reciprocal <- if (all(is.finite(values)) && max(values) > 0) {
    min(values) / max(values)
  } else {
    0
  }
  status <- if (!is.finite(reciprocal) ||
                reciprocal <= .JACOBIAN_RCOND_SINGULAR) {
    "singular"
  } else if (reciprocal <= .JACOBIAN_RCOND_ILL) {
    "ill_conditioned"
  } else {
    "well_conditioned"
  }
  list(
    parameterization = "standardized residual by log(shape),log(rate)",
    singular_values = values, reciprocal_condition = reciprocal,
    status = status,
    thresholds = c(
      singular = .JACOBIAN_RCOND_SINGULAR,
      ill_conditioned = .JACOBIAN_RCOND_ILL
    )
  )
}


.a2_run_newton <- function(
    J, target, eta_start, M, abs_tol, rel_tol, tol_step,
    max_iter, damping, verbose) {
  eta <- eta_start
  rows <- list()
  best <- list(
    eta = eta, moments = NULL, residual = NULL, objective = Inf
  )
  solved <- FALSE
  reason <- "max_iterations"
  error_message <- NA_character_
  last_conditioning <- NULL
  last_derivative_status <- NA_character_
  started <- proc.time()[["elapsed"]]

  consider <- function(eta_value, moments) {
    residual <- .a2_residual_contract(moments, target, abs_tol, rel_tol)
    objective <- sum(residual$standardized^2)
    list(residual = residual, objective = objective)
  }

  for (iter in seq_len(max_iter)) {
    a <- exp(eta[1L])
    b <- exp(eta[2L])
    evaluation <- tryCatch(
      moments_with_jacobian(J, a, b, M),
      error = function(error) error
    )
    if (inherits(evaluation, "error")) {
      reason <- "moment_or_jacobian_evaluation_failed"
      error_message <- conditionMessage(evaluation)
      rows[[length(rows) + 1L]] <- .a2_trace_row(
        iter, eta, reason = reason
      )
      break
    }

    moments <- c(mean = evaluation$mean, variance = evaluation$var)
    assessment <- consider(eta, moments)
    residual <- assessment$residual
    objective <- assessment$objective
    if (is.finite(objective) && objective < best$objective) {
      best <- list(
        eta = eta, moments = moments, residual = residual,
        objective = objective
      )
    }
    jacobian_log <- evaluation$jacobian %*% diag(c(a, b))
    jacobian_standardized <- sweep(
      jacobian_log, 1L, residual$tolerance, "/"
    )
    conditioning <- .a2_scaled_jacobian_condition(jacobian_standardized)
    last_conditioning <- conditioning
    derivative_status <- evaluation$derivative_diagnostics$status %||%
      NA_character_
    last_derivative_status <- derivative_status
    determinant <- suppressWarnings(as.numeric(det(jacobian_log)))

    if (isTRUE(residual$passed)) {
      solved <- TRUE
      reason <- "residual_tolerance_met"
      rows[[length(rows) + 1L]] <- .a2_trace_row(
        iter, eta, moments, residual, determinant = determinant,
        reciprocal_condition = conditioning$reciprocal_condition,
        jacobian_status = conditioning$status,
        derivative_status = derivative_status, reason = reason
      )
      break
    }
    if (!identical(conditioning$status, "well_conditioned")) {
      reason <- paste0(conditioning$status, "_scaled_jacobian")
      rows[[length(rows) + 1L]] <- .a2_trace_row(
        iter, eta, moments, residual, determinant = determinant,
        reciprocal_condition = conditioning$reciprocal_condition,
        jacobian_status = conditioning$status,
        derivative_status = derivative_status, reason = reason
      )
      break
    }

    delta <- tryCatch(
      -solve(jacobian_standardized, residual$standardized),
      error = function(error) error
    )
    if (inherits(delta, "error") || length(delta) != 2L ||
        any(!is.finite(delta))) {
      reason <- "newton_linear_solve_failed"
      error_message <- if (inherits(delta, "error")) {
        conditionMessage(delta)
      } else {
        "Newton linear solve returned a non-finite step"
      }
      rows[[length(rows) + 1L]] <- .a2_trace_row(
        iter, eta, moments, residual, determinant = determinant,
        reciprocal_condition = conditioning$reciprocal_condition,
        jacobian_status = conditioning$status,
        derivative_status = derivative_status, reason = reason
      )
      break
    }
    full_step_norm <- sqrt(sum(delta^2))
    if (full_step_norm < tol_step) {
      reason <- "stagnation_small_step"
      rows[[length(rows) + 1L]] <- .a2_trace_row(
        iter, eta, moments, residual, step = 0,
        step_norm = full_step_norm, determinant = determinant,
        reciprocal_condition = conditioning$reciprocal_condition,
        jacobian_status = conditioning$status,
        derivative_status = derivative_status, reason = reason
      )
      break
    }

    max_line <- if (isTRUE(damping)) 20L else 1L
    factor <- 1
    accepted <- FALSE
    trial_eta <- eta
    line_iterations <- 0L
    for (line_iter in seq_len(max_line)) {
      line_iterations <- line_iter
      candidate_eta <- eta + factor * delta
      inside <- all(is.finite(candidate_eta)) &&
        all(candidate_eta >= .LOG_BOUNDS_DEFAULT[1L]) &&
        all(candidate_eta <= .LOG_BOUNDS_DEFAULT[2L])
      candidate <- if (inside) {
        tryCatch(
          .a2_moments_at_order(
            J, exp(candidate_eta[1L]), exp(candidate_eta[2L]), M
          ),
          error = function(error) NULL
        )
      } else {
        NULL
      }
      if (!is.null(candidate)) {
        candidate_objective <- consider(candidate_eta, candidate)$objective
        if (is.finite(candidate_objective) && candidate_objective < objective) {
          accepted <- TRUE
          trial_eta <- candidate_eta
          break
        }
      }
      factor <- factor / 2
    }
    reason_row <- if (accepted) {
      "step_accepted"
    } else if (isTRUE(damping)) {
      "line_search_failed"
    } else {
      "undamped_step_rejected"
    }
    rows[[length(rows) + 1L]] <- .a2_trace_row(
      iter, eta, moments, residual, step = factor,
      step_norm = factor * full_step_norm,
      line_iterations = line_iterations, accepted = accepted,
      determinant = determinant,
      reciprocal_condition = conditioning$reciprocal_condition,
      jacobian_status = conditioning$status,
      derivative_status = derivative_status, reason = reason_row
    )
    if (!accepted) {
      reason <- reason_row
      break
    }
    eta <- trial_eta
    if (isTRUE(verbose)) {
      cat(sprintf(
        "iter %d: standardized RMS %.3g; step %.3g; rcond %.3g\n",
        iter, residual$standardized_norm, factor,
        conditioning$reciprocal_condition
      ))
    }
  }

  # A final accepted step can occur on the last permitted iteration.
  post <- tryCatch(
    .a2_moments_at_order(J, exp(eta[1L]), exp(eta[2L]), M),
    error = function(error) NULL
  )
  if (!is.null(post)) {
    assessment <- consider(eta, post)
    if (is.finite(assessment$objective) &&
        assessment$objective < best$objective) {
      best <- list(
        eta = eta, moments = post, residual = assessment$residual,
        objective = assessment$objective
      )
    }
    if (!solved && isTRUE(assessment$residual$passed)) {
      solved <- TRUE
      reason <- "post_loop_residual_tolerance_met"
    }
  }

  history <- if (length(rows)) do.call(rbind, rows) else data.frame()
  attempt <- list(
    method = "scaled_log_newton", start = eta_start,
    bounds = list(
      parameterization = "log(shape), log(rate)",
      lower = rep(.LOG_BOUNDS_DEFAULT[1L], 2L),
      upper = rep(.LOG_BOUNDS_DEFAULT[2L], 2L)
    ),
    control = list(
      max_iter = max_iter, damping = damping, tol_step = tol_step,
      residual_absolute_tolerance = abs_tol,
      residual_relative_tolerance = rel_tol
    ),
    exit_code = if (solved) 0L else 1L,
    message = reason, iterations = nrow(history), evaluations = NA_integer_,
    candidate_objective = best$objective,
    elapsed_seconds = unname(proc.time()[["elapsed"]] - started),
    warning = NA_character_, error = error_message,
    reason_code = reason,
    candidate = c(
      log_shape = unname(best$eta[1L]),
      log_rate = unname(best$eta[2L]),
      shape = unname(exp(best$eta[1L])),
      rate = unname(exp(best$eta[2L]))
    )
  )
  list(
    solved = solved, reason = reason, best = best, trace = history,
    attempt = attempt, conditioning = last_conditioning,
    derivative_status = last_derivative_status
  )
}


.a2_run_fallback <- function(
    J, target, start, M, abs_tol, rel_tol, incumbent_objective) {
  warnings <- character()
  objective <- function(log_parameters) {
    valid <- is.numeric(log_parameters) && length(log_parameters) == 2L &&
      all(is.finite(log_parameters)) &&
      all(log_parameters >= .LOG_BOUNDS_DEFAULT[1L]) &&
      all(log_parameters <= .LOG_BOUNDS_DEFAULT[2L])
    if (!valid) return(.Machine$double.xmax / 1024)
    moments <- tryCatch(
      .a2_moments_at_order(
        J, exp(log_parameters[1L]), exp(log_parameters[2L]), M
      ),
      error = function(error) NULL
    )
    if (is.null(moments)) return(.Machine$double.xmax / 1024)
    residual <- .a2_residual_contract(moments, target, abs_tol, rel_tol)
    value <- sum(residual$standardized^2)
    if (is.finite(value)) value else .Machine$double.xmax / 1024
  }
  started <- proc.time()[["elapsed"]]
  result <- tryCatch(
    withCallingHandlers(
      .a2_nelder_mead(start, objective),
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(error) error
  )
  elapsed <- proc.time()[["elapsed"]] - started
  reason <- "optimizer_error"
  error_message <- NA_character_
  exit_code <- iterations <- evaluations <- NA_integer_
  candidate_value <- NA_real_
  optimizer_value <- NA_real_
  candidate <- NULL
  candidate_record <- NULL
  selected <- FALSE
  success <- FALSE

  if (inherits(result, "error")) {
    error_message <- conditionMessage(result)
  } else if (!is.list(result) || !is.numeric(result$par) ||
             length(result$par) != 2L || any(!is.finite(result$par))) {
    reason <- "invalid_optimizer_result"
    error_message <- "Nelder-Mead returned a malformed log-parameter candidate"
  } else {
    eta <- as.numeric(result$par)
    exit_code <- as.integer(result$convergence %||% NA_integer_)
    iterations <- as.integer(result$counts[["function"]] %||% NA_integer_)
    evaluations <- iterations
    optimizer_value <- as.numeric(result$value %||% NA_real_)
    inside <- all(eta >= .LOG_BOUNDS_DEFAULT[1L]) &&
      all(eta <= .LOG_BOUNDS_DEFAULT[2L])
    moments <- if (inside) {
      tryCatch(
        .a2_moments_at_order(J, exp(eta[1L]), exp(eta[2L]), M),
        error = function(error) NULL
      )
    } else {
      NULL
    }
    if (is.null(moments)) {
      reason <- if (inside) {
        "fallback_candidate_evaluation_failed"
      } else {
        "fallback_candidate_outside_bounds"
      }
    } else {
      residual <- .a2_residual_contract(moments, target, abs_tol, rel_tol)
      assessed_objective <- sum(residual$standardized^2)
      candidate_value <- assessed_objective
      candidate_record <- list(
        eta = eta, moments = moments, residual = residual,
        objective = assessed_objective
      )
      candidate <- c(
        log_shape = unname(eta[1L]), log_rate = unname(eta[2L]),
        shape = unname(exp(eta[1L])), rate = unname(exp(eta[2L]))
      )
      selected <- is.finite(assessed_objective) &&
        assessed_objective < incumbent_objective
      exit_ok <- length(exit_code) == 1L && !is.na(exit_code) &&
        exit_code == 0L
      success <- exit_ok && isTRUE(residual$passed) && selected
      reason <- if (!exit_ok) {
        "optimizer_exit_nonzero"
      } else if (!isTRUE(residual$passed)) {
        "fallback_residual_tolerance_not_met"
      } else if (!selected) {
        "fallback_candidate_not_better"
      } else {
        "fallback_residual_tolerance_met"
      }
    }
  }

  attempt <- list(
    method = "nelder_mead_log", start = start,
    bounds = list(
      parameterization = "log(shape), log(rate)",
      lower = rep(.LOG_BOUNDS_DEFAULT[1L], 2L),
      upper = rep(.LOG_BOUNDS_DEFAULT[2L], 2L),
      enforcement = "finite objective penalty"
    ),
    control = list(maxit = 1000L, reltol = 1e-12),
    exit_code = exit_code, message = reason, iterations = iterations,
    evaluations = evaluations, candidate_objective = candidate_value,
    optimizer_reported_objective = optimizer_value,
    elapsed_seconds = unname(elapsed),
    warning = if (length(warnings)) paste(warnings, collapse = "; ") else NA_character_,
    error = error_message, reason_code = reason, candidate = candidate
  )
  list(
    success = success, selected = selected, candidate = candidate_record,
    reason = reason, attempt = attempt
  )
}


# Canonical dpprior.result/1 adapter -----------------------------------------

# A2-MN keeps its historical solver controls, but decision evidence is capped
# at the frozen central truth policy.  A relaxed stopping tolerance may select
# a numerical candidate; it cannot relax the public adequacy claim.
.a2_schema_controls <- function(max_iter, damping, use_fallback, tol_step) {
  list(
    max_iter = as.integer(max_iter),
    damping = damping,
    use_fallback = use_fallback,
    tol_step = tol_step,
    log_bounds = unname(as.numeric(.LOG_BOUNDS_DEFAULT)),
    boundary_tol = 1e-6,
    line_search_max = 20L,
    jacobian_rcond_singular = .JACOBIAN_RCOND_SINGULAR,
    jacobian_rcond_ill = .JACOBIAN_RCOND_ILL,
    fallback = list(
      maxit = 1000L, reltol = 1e-12,
      finite_penalty = .Machine$double.xmax / 1024
    ),
    selection_tolerance = 0
  )
}


.a2_schema_tolerances <- function(tol_F, tol_rel, tol_step,
                                  verification_abs_tol,
                                  verification_rel_tol) {
  list(
    K_adequacy = list(
      absolute = min(tol_F, 1e-8),
      relative = min(tol_rel, 1e-8),
      scale_formula = "max(abs(target),1)"
    ),
    K_stability = list(
      absolute = min(verification_abs_tol, 1e-10),
      relative = min(verification_rel_tol, 1e-8),
      scale_floor = 1
    ),
    step = min(tol_step, 1e-10),
    boundary = 1e-6
  )
}


.a2_schema_parameters <- function(a, b) {
  .dpprior_new_parameters(a, b, "log(shape), log(rate)")
}


.a2_schema_achieved <- function(moments, M, source) {
  list(K = list(
    mean = unname(as.numeric(moments[["mean"]])),
    variance = unname(as.numeric(moments[["variance"]])),
    estimand = "K_J", source = source, M = as.integer(M)
  ))
}


.a2_schema_residuals <- function(moments, target) {
  list(K = list(
    mean = unname(as.numeric(moments[["mean"]] - target[["mean"]])),
    variance = unname(as.numeric(
      moments[["variance"]] - target[["variance"]]
    ))
  ))
}


.a2_schema_fresh_evidence <- function(J, a, b, M, M_verify, target,
                                      tolerances) {
  selected <- tryCatch(
    exact_K_moments(J, a, b, M = M),
    error = function(condition) condition
  )
  if (inherits(selected, "condition")) {
    return(list(
      ok = FALSE, code = "selected_order_recomputation_failed",
      message = conditionMessage(selected), cause = selected
    ))
  }
  selected_values <- c(mean = selected$mean, variance = selected$var)
  adequacy <- tolerances$K_adequacy
  selected_contract <- .a2_residual_contract(
    selected_values, target, adequacy$absolute, adequacy$relative
  )
  if (is.null(M_verify) || is.na(M_verify)) {
    return(list(
      ok = FALSE, code = "verification_order_exceeds_ceiling",
      message = paste(
        "The required independent A2-MN verification order is unavailable;",
        "refit with M <= 256 so max(2*M, M+40) is supported."
      ),
      selected = selected_values, selected_contract = selected_contract
    ))
  }
  verifier <- tryCatch(
    exact_K_moments(J, a, b, M = M_verify),
    error = function(condition) condition
  )
  if (inherits(verifier, "condition")) {
    return(list(
      ok = FALSE, code = "higher_order_recomputation_failed",
      message = conditionMessage(verifier), cause = verifier
    ))
  }
  verifier_values <- c(mean = verifier$mean, variance = verifier$var)
  verifier_contract <- .a2_residual_contract(
    verifier_values, target, adequacy$absolute, adequacy$relative
  )
  stability_policy <- tolerances$K_stability
  stability_delta <- c(
    K.mean = abs(selected_values[["mean"]] - verifier_values[["mean"]]),
    K.variance = abs(
      selected_values[["variance"]] - verifier_values[["variance"]]
    )
  )
  stability_tolerance <- stability_policy$absolute +
    stability_policy$relative * pmax(
      abs(c(
        K.mean = selected_values[["mean"]],
        K.variance = selected_values[["variance"]]
      )),
      abs(c(
        K.mean = verifier_values[["mean"]],
        K.variance = verifier_values[["variance"]]
      )),
      stability_policy$scale_floor
    )
  stability <- .dpprior_new_stability(
    delta = stability_delta,
    tolerance = stability_tolerance,
    formula = c(
      K.mean = "absolute_plus_relative_max",
      K.variance = "absolute_plus_relative_max"
    ),
    scale_floor = c(
      K.mean = stability_policy$scale_floor,
      K.variance = stability_policy$scale_floor
    ),
    source = "independent_verifier"
  )
  list(
    ok = TRUE,
    selected = selected_values,
    verifier = verifier_values,
    selected_contract = selected_contract,
    verifier_contract = verifier_contract,
    stability = stability,
    passed = isTRUE(selected_contract$passed) &&
      isTRUE(verifier_contract$passed) && isTRUE(stability$passed),
    reason = if (!isTRUE(selected_contract$passed)) {
      "selected_residual_tolerance_not_met"
    } else if (!isTRUE(verifier_contract$passed)) {
      "verification_residual_tolerance_not_met"
    } else if (!isTRUE(stability$passed)) {
      "quadrature_order_disagreement"
    } else {
      "independent_verification_passed"
    }
  )
}


.a2_schema_snapshot <- function(parameters, moments, M, target, tolerances,
                                source) {
  achieved <- .a2_schema_achieved(moments, M, source)
  residuals <- .a2_schema_residuals(moments, target)
  .dpprior_new_snapshot(
    parameters = parameters, M = as.integer(M), achieved = achieved,
    residuals = residuals, tolerances = tolerances, finite = TRUE,
    source = source
  )
}


.a2_schema_attempt_error <- function(attempt) {
  value <- attempt[["error"]]
  if (is.null(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
    return(NULL)
  }
  list(
    class = "simpleError", code = attempt[["reason_code"]] %||%
      "optimizer_error", message = as.character(value)
  )
}


.a2_schema_attempt_parameters <- function(attempt) {
  candidate <- attempt[["candidate"]]
  if (!is.numeric(candidate) || is.object(candidate) || is.null(names(candidate)) ||
      !all(c("shape", "rate") %in% names(candidate)) ||
      any(!is.finite(candidate[c("shape", "rate")])) ||
      any(candidate[c("shape", "rate")] <= 0)) {
    return(NULL)
  }
  .a2_schema_parameters(
    unname(candidate[["shape"]]), unname(candidate[["rate"]])
  )
}


.a2_schema_attempt <- function(attempt, id, stage, selected,
                               candidate_parameters) {
  scalar_or_null <- function(value, integer = FALSE) {
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value)) return(NULL)
    if (integer) as.integer(value) else unname(as.numeric(value))
  }
  exit_code <- scalar_or_null(attempt[["exit_code"]], integer = TRUE)
  iterations <- scalar_or_null(attempt[["iterations"]], integer = TRUE)
  evaluations_value <- scalar_or_null(
    attempt[["evaluations"]], integer = TRUE
  )
  evaluations <- if (is.null(evaluations_value)) NULL else
    list(function_count = evaluations_value)
  elapsed <- scalar_or_null(attempt[["elapsed_seconds"]])
  warnings <- attempt[["warning"]]
  warnings <- if (is.character(warnings) && length(warnings) == 1L &&
      !is.na(warnings) && nzchar(warnings)) warnings else character()
  error <- .a2_schema_attempt_error(attempt)
  start <- attempt[["start"]]
  if (!is.numeric(start) || any(!is.finite(start))) start <- NULL
  bounds <- attempt[["bounds"]]
  bounds <- if (is.list(bounds) && all(c("lower", "upper") %in% names(bounds))) {
    list(
      lower = unname(as.numeric(bounds[["lower"]])),
      upper = unname(as.numeric(bounds[["upper"]]))
    )
  } else {
    NULL
  }
  control <- attempt[["control"]]
  if (!is.list(control)) control <- NULL
  nullable <- list(
    start = start, bounds = bounds, control = control,
    exit_code = exit_code, iterations = iterations,
    evaluations = evaluations, candidate_parameters = candidate_parameters,
    candidate_objective = NULL, elapsed_seconds = elapsed
  )
  unavailable <- vapply(
    names(nullable)[vapply(nullable, is.null, logical(1))],
    function(field) switch(
      field,
      candidate_objective = paste(
        "source solver recorded a mean-square residual objective; canonical",
        "selection uses a fresh sum-of-squares objective"
      ),
      evaluations = "source attempt did not retain named evaluation counts",
      "source attempt did not retain this evidence"
    ),
    character(1)
  )
  reason <- if (selected) {
    "selected"
  } else if (!is.null(candidate_parameters)) {
    "eligible_not_selected"
  } else if (!is.null(error)) {
    "optimizer_error"
  } else {
    "no_candidate"
  }
  .dpprior_new_attempt(
    id = id, stage = stage, method = attempt[["method"]],
    start = start, bounds = bounds, control = control,
    exit_code = exit_code, message = attempt[["message"]] %||% reason,
    iterations = iterations, evaluations = evaluations,
    candidate_parameters = candidate_parameters,
    candidate_objective = NULL, elapsed_seconds = elapsed,
    warnings = warnings, error = error, selected = selected,
    reason_code = reason, unavailable = unavailable
  )
}


.a2_schema_candidate <- function(id, attempt_id, attempt, parameters,
                                 snapshot, target, tolerances, selected) {
  target_tolerance <- c(
    mean = tolerances$K_adequacy$absolute +
      tolerances$K_adequacy$relative * max(abs(target[["mean"]]), 1),
    variance = tolerances$K_adequacy$absolute +
      tolerances$K_adequacy$relative * max(abs(target[["variance"]]), 1)
  )
  residual <- unlist(snapshot$residuals$K, use.names = TRUE)
  objective <- sum((residual / target_tolerance)^2)
  source <- paste0("candidate:", id)
  exit_code <- attempt[["exit_code"]]
  execution_success <- is.numeric(exit_code) && length(exit_code) == 1L &&
    !is.na(exit_code) && identical(as.integer(exit_code), 0L) &&
    is.null(.a2_schema_attempt_error(attempt))
  .dpprior_new_candidate_evaluation(
    id = id, attempt_id = attempt_id, method = attempt[["method"]],
    generator = "direct_attempt", parameters = parameters,
    objective_kind = "standardized_residual",
    recorded_objective = NULL, fresh_objective = objective,
    selection_objective = objective, objective_tolerance = 0,
    selected_snapshot = snapshot, verifier_snapshot = NULL,
    checks = list(candidate_finite = .dpprior_new_check(
      value = c(parameters = TRUE, snapshot = TRUE, objective = TRUE),
      reference = c(parameters = TRUE, snapshot = TRUE, objective = TRUE),
      tolerance = NULL, operator = "identical", source = source
    )),
    execution_success = execution_success,
    optimizer_supported = execution_success,
    selected = selected, source = "A2-MN fresh selected-order candidate audit"
  )
}


.a2_schema_verification <- function(parameters, evidence, M, M_verify,
                                    target, tolerances, passed, reason) {
  selected_snapshot <- .a2_schema_snapshot(
    parameters, evidence$selected, M, target, tolerances, "selected_order"
  )
  verifier_snapshot <- .a2_schema_snapshot(
    parameters, evidence$verifier, M_verify, target, tolerances,
    "independent_verifier"
  )
  selected_residual <- unlist(
    selected_snapshot$residuals$K, use.names = TRUE
  )
  verifier_residual <- unlist(
    verifier_snapshot$residuals$K, use.names = TRUE
  )
  adequacy_tolerance <- c(
    selected.mean = evidence$selected_contract$tolerance[["mean"]],
    selected.variance = evidence$selected_contract$tolerance[["variance"]],
    refined.mean = evidence$verifier_contract$tolerance[["mean"]],
    refined.variance = evidence$verifier_contract$tolerance[["variance"]],
    selected.standardized_norm = 1,
    refined.standardized_norm = 1
  )
  adequacy_value <- c(
    selected.mean = abs(selected_residual[["mean"]]),
    selected.variance = abs(selected_residual[["variance"]]),
    refined.mean = abs(verifier_residual[["mean"]]),
    refined.variance = abs(verifier_residual[["variance"]]),
    selected.standardized_norm = evidence$selected_contract$standardized_norm,
    refined.standardized_norm = evidence$verifier_contract$standardized_norm
  )
  components <- list(
    residual_adequacy = .dpprior_new_check(
      value = adequacy_value,
      reference = stats::setNames(rep(0, length(adequacy_value)),
                                  names(adequacy_value)),
      tolerance = adequacy_tolerance, operator = "lte",
      source = "independent_verifier"
    ),
    order_stability = .dpprior_new_check(
      value = evidence$stability$delta,
      reference = stats::setNames(
        rep(0, length(evidence$stability$delta)),
        names(evidence$stability$delta)
      ),
      tolerance = evidence$stability$tolerance, operator = "lte",
      source = "independent_verifier"
    ),
    parameter_identity = .dpprior_new_check(
      value = unlist(parameters[c("a", "b")]),
      reference = unlist(parameters[c("a", "b")]),
      tolerance = NULL, operator = "identical",
      source = "independent_verifier"
    )
  )
  invariants <- list(
    finite_parameters = .dpprior_new_check(
      value = all(is.finite(unlist(parameters[c("a", "b")]))),
      reference = TRUE, tolerance = NULL, operator = "identical",
      source = "independent_verifier"
    ),
    K_support = .dpprior_new_check(
      value = TRUE, reference = TRUE, tolerance = NULL,
      operator = "identical", source = "independent_verifier"
    )
  )
  .dpprior_new_verification(
    method = "independent_higher_order_moment_recomputation",
    performed = TRUE, passed = passed, reason = reason,
    settings = list(
      M_selected = as.integer(M),
      M_verification_required = .a2_required_verification_order(M),
      M_verification = as.integer(M_verify)
    ),
    selected_snapshot = selected_snapshot,
    verifier_snapshot = verifier_snapshot,
    stability = evidence$stability,
    components = components, invariants = invariants
  )
}


.a2_schema_compatibility_view <- function(
    a, b, J, method, status, usable, verified, primary, attempts,
    termination, diagnostics, evidence, initialization,
    required_action = "none") {
  safe_diagnostics <- list(
    a0 = diagnostics$a0, b0 = diagnostics$b0,
    tol_F = diagnostics$tol_F, tol_rel = diagnostics$tol_rel,
    tol_step = diagnostics$tol_step, M = as.integer(diagnostics$M),
    M_verify = if (is.na(diagnostics$M_verify)) NULL else
      as.integer(diagnostics$M_verify),
    M_verify_required = as.integer(diagnostics$M_verify_required),
    verification_available = diagnostics$verification_available,
    fallback_used = diagnostics$fallback_used,
    fallback_attempted = diagnostics$fallback_attempted,
    derivative_status = if (
      is.character(diagnostics$derivative_status) &&
        length(diagnostics$derivative_status) == 1L &&
        !is.na(diagnostics$derivative_status) &&
        nzchar(diagnostics$derivative_status)
    ) diagnostics$derivative_status else "unavailable",
    quasi_improper = diagnostics$quasi_improper,
    boundary_distance_log_scale = diagnostics$boundary_distance_log_scale
  )
  selected_fit <- if (
    is.numeric(evidence$selected) && length(evidence$selected) == 2L &&
      all(is.finite(evidence$selected)) &&
      is.list(evidence$selected_contract)
  ) {
    list(
      mu_K = evidence$selected[["mean"]],
      var_K = evidence$selected[["variance"]],
      residual = sqrt(sum(evidence$selected_contract$raw^2)),
      scaled_residual = max(abs(evidence$selected_contract$standardized)),
      residual_components = evidence$selected_contract$raw
    )
  } else {
    NULL
  }
  initialization_attempts <- initialization$attempts
  initialization_summary <- list(
    method = initialization$method,
    projection_policy = initialization$projection_policy,
    target_was_not_projected_for_A2 =
      initialization$target_was_not_projected_for_A2,
    requested_A2_target = initialization$requested_A2_target,
    target_projection = initialization$target_projection,
    attempt_methods = vapply(
      initialization_attempts, function(attempt) attempt$method, character(1)
    ),
    attempt_statuses = vapply(
      initialization_attempts, function(attempt) attempt$status, character(1)
    ),
    attempt_reason_codes = vapply(
      initialization_attempts,
      function(attempt) attempt$reason_code,
      character(1)
    )
  )
  list(
    contract = "R11_pre_dpprior.result/1",
    authority = "non_authoritative", lossy = TRUE,
    consumer_policy = "ignored_by_scientific_and_decision_consumers",
    source_status = status, source_usable = usable,
    source_verified = verified, source_method = method,
    numerical_candidate = list(a = a, b = b, J = as.integer(J)),
    selected_fit = selected_fit,
    solver_diagnostics = safe_diagnostics,
    iterations = as.integer(nrow(primary$trace)),
    termination = termination,
    source_attempt_methods = vapply(
      attempts, function(attempt) attempt$method, character(1)
    ),
    source_attempt_reason_codes = vapply(
      attempts, function(attempt) attempt$reason_code, character(1)
    ),
    initialization = initialization_summary,
    required_action = required_action
  )
}


.a2_schema_provenance <- function(method, fallback_used, status, target_K) {
  source_commit <- getOption("DPprior.source_commit", NULL)
  if (!is.character(source_commit) || is.object(source_commit) ||
      is.null(dim(source_commit)) == FALSE || length(source_commit) != 1L ||
      is.na(source_commit) || !nzchar(source_commit)) {
    source_commit <- NULL
  }
  approximate <- identical(status, "approximate")
  .dpprior_new_provenance(
    requested_method = "A2-MN", selected_method = method,
    is_fallback = fallback_used,
    approximation = list(
      active = approximate, opt_in = FALSE,
      kind = if (approximate) "numerical_candidate_not_decision_ready" else NULL,
      warning_code = if (approximate) "a2_moment_approximate" else NULL
    ),
    projection = target_K$provenance$projection,
    parameterization = "log(shape), log(rate)",
    backend = list(
      package = "DPprior",
      package_version = tryCatch(
        as.character(utils::packageVersion("DPprior")),
        error = function(condition) "development"
      ),
      implementation = "R/11_a2_newton.R:DPprior_a2_newton",
      source_commit = source_commit
    ),
    input_fit = NULL,
    migration = list(
      source_schema = "native", adapter = "none", lossless = TRUE,
      missing_evidence = character(), warnings = character()
    ),
    legacy = list(active = FALSE, contract = NULL, deprecation_stage = NULL)
  )
}


.a2_schema_scaling <- function(target) {
  values <- list(
    mean = max(abs(target[["mean"]]), 1),
    variance = max(abs(target[["variance"]]), 1)
  )
  .dpprior_new_scaling(
    requested = list(K = values), used = list(K = values),
    formula = "max(abs(target),1)", values = values,
    fixed_from_input = TRUE, change_reason = ""
  )
}


.a2_schema_settings <- function(method, controls) {
  list(
    method = method, controls = controls,
    parameterization = "log(shape), log(rate)"
  )
}


.a2_schema_compatibility <- function(view, converged) {
  .dpprior_new_compatibility(
    views = list(converged = converged, legacy_v2 = view),
    deprecations = list(legacy_v2 = list(
      code = "a2_moment_legacy_view",
      first_deprecated_version = "2.0.0",
      removal_floor = "not_scheduled",
      authority = "non_authoritative", lossy = TRUE,
      consumer_policy = "ignored_by_scientific_and_decision_consumers"
    ))
  )
}


.a2_schema_append_finite_aliases <- function(result) {
  .dpprior_append_compatibility_v2(
    result,
    aliases = c(
      a = "parameters.a", b = "parameters.b",
      converged = "compatibility.views.converged",
      iterations = "compatibility.views.legacy_v2.iterations",
      termination = "compatibility.views.legacy_v2.termination",
      fit = "compatibility.views.legacy_v2.selected_fit",
      diagnostics = "compatibility.views.legacy_v2.solver_diagnostics",
      trace = "computation.trace", attempts = "computation.attempts"
    ),
    views = result$compatibility$views,
    deprecations = result$compatibility$deprecations
  )
}


.a2_schema_finite_result <- function(
    J, target, a, b, method, status, message, evidence, M, M_verify,
    M_required, controls, tolerances, source_attempts, primary,
    fallback_attempted, fallback_selected, fallback_reason,
    initialization, legacy_status, legacy_usable, legacy_verified,
    legacy_termination, legacy_diagnostics) {
  parameters <- .a2_schema_parameters(a, b)
  verification <- .a2_schema_verification(
    parameters, evidence, M, M_verify, target, tolerances,
    passed = status %in% c("converged", "boundary"), reason = evidence$reason
  )
  selected_snapshot <- verification$selected_snapshot
  selected_index <- if (fallback_selected) length(source_attempts) else 1L
  canonical_attempts <- candidate_evaluations <- list()
  for (index in seq_along(source_attempts)) {
    source_attempt <- source_attempts[[index]]
    candidate_parameters <- .a2_schema_attempt_parameters(source_attempt)
    candidate_snapshot <- NULL
    if (!is.null(candidate_parameters)) {
      candidate_moments <- tryCatch(
        exact_K_moments(
          J, candidate_parameters$a, candidate_parameters$b, M = M
        ),
        error = function(condition) NULL
      )
      if (is.list(candidate_moments)) {
        candidate_values <- c(
          mean = candidate_moments$mean, variance = candidate_moments$var
        )
        candidate_snapshot <- .a2_schema_snapshot(
          candidate_parameters, candidate_values, M, target, tolerances,
          "selected_order"
        )
      } else {
        candidate_parameters <- NULL
      }
    }
    selected_attempt <- identical(index, selected_index)
    if (selected_attempt) {
      candidate_parameters <- parameters
      candidate_snapshot <- selected_snapshot
    }
    stage <- if (identical(index, 1L)) "primary" else "fallback"
    canonical_attempts[[index]] <- .a2_schema_attempt(
      source_attempt, paste0("attempt-", index), stage,
      selected_attempt, candidate_parameters
    )
    if (!is.null(candidate_parameters)) {
      candidate_evaluations[[length(candidate_evaluations) + 1L]] <-
        .a2_schema_candidate(
          paste0("candidate-", index), paste0("attempt-", index),
          source_attempt, candidate_parameters, candidate_snapshot,
          target, tolerances, selected_attempt
        )
    }
  }
  fallback <- if (!fallback_attempted) {
    .dpprior_new_fallback()
  } else if (fallback_selected) {
    .dpprior_new_fallback(
      attempted = TRUE, used = TRUE,
      trigger_attempt_id = "attempt-1",
      selected_attempt_id = paste0("attempt-", selected_index),
      reason_code = fallback_reason,
      message = source_attempts[[selected_index]]$message,
      outcome = "selected"
    )
  } else {
    last <- source_attempts[[length(source_attempts)]]
    .dpprior_new_fallback(
      attempted = TRUE, used = FALSE,
      trigger_attempt_id = "attempt-1", selected_attempt_id = NULL,
      reason_code = fallback_reason,
      message = last$message,
      outcome = if (is.null(.a2_schema_attempt_parameters(last))) {
        "failed"
      } else {
        "rejected"
      }
    )
  }
  orders <- .dpprior_new_orders(
    M_requested = as.integer(M), M_selected = as.integer(M),
    M_verification_required = as.integer(M_required),
    M_verification_used = as.integer(M_verify),
    requested_reason = "public M argument",
    selected_reason = "selected candidate fit order",
    verification_required_reason = "max(2*M_selected,M_selected+40)",
    verification_used_reason = "independent higher-order verifier"
  )
  termination_source <- if (identical(status, "approximate")) {
    "candidate_evaluation"
  } else if (fallback_selected) {
    "fallback_optimizer"
  } else {
    "optimizer"
  }
  termination_code <- switch(
    status, converged = "converged", boundary = "boundary",
    approximate = "approximate", "failed"
  )
  termination <- .dpprior_new_termination(
    code = termination_code, message = message,
    source = termination_source,
    iterations = canonical_attempts[[selected_index]]$iterations,
    boundary_reason = if (identical(status, "boundary")) {
      "selected log-parameter is within the fixed boundary tolerance"
    } else {
      NULL
    }
  )
  computation <- .dpprior_new_computation(
    request = .a2_schema_settings("A2-MN", controls),
    used = .a2_schema_settings(method, controls),
    orders = orders, scaling = .a2_schema_scaling(target),
    attempts = canonical_attempts,
    candidate_evaluations = candidate_evaluations,
    selected_candidate_id = paste0("candidate-", selected_index),
    selected_attempt_id = paste0("attempt-", selected_index),
    fallback = fallback, termination = termination,
    trace = primary$trace,
    resources = list(
      implementation = "scaled log-Newton with optional log-NM fallback",
      initialization_method = initialization$method
    )
  )
  target_K <- DPprior_target_K(
    J = J, mu_K = target[["mean"]], var_K = target[["variance"]]
  )
  legacy_view <- .a2_schema_compatibility_view(
    a, b, J, method, legacy_status, legacy_usable, legacy_verified,
    primary, source_attempts, legacy_termination, legacy_diagnostics,
    evidence, initialization
  )
  result <- .dpprior_new_fit(
    mode = "a2_moment", method = method, J = J, status = status,
    usable = status %in% c("converged", "boundary"),
    verified = status %in% c("converged", "boundary"),
    message = message, parameters = parameters,
    target = list(K = target_K),
    achieved = selected_snapshot$achieved,
    residuals = selected_snapshot$residuals,
    tolerances = tolerances, computation = computation,
    verification = verification,
    provenance = .a2_schema_provenance(
      method, fallback_selected, status, target_K
    ),
    compatibility = .a2_schema_compatibility(
      legacy_view,
      status %in% c("converged", "boundary")
    )
  )
  .a2_schema_append_finite_aliases(result)
}


.a2_schema_failed_result <- function(
    J, target, a, b, source_method, source_status, source_usable,
    source_verified, message, code, action, evidence, M, M_required,
    controls, tolerances, source_attempts, primary, initialization,
    legacy_termination, legacy_diagnostics) {
  target_K <- DPprior_target_K(
    J = J, mu_K = target[["mean"]], var_K = target[["variance"]]
  )
  orders <- .dpprior_new_orders(
    M_requested = as.integer(M), M_selected = as.integer(M),
    M_verification_required = as.integer(M_required),
    M_verification_used = NULL,
    requested_reason = "public M argument",
    selected_reason = "non-authoritative source candidate fit order",
    verification_required_reason = "max(2*M_selected,M_selected+40)",
    verification_used_reason = "unavailable above quadrature ceiling"
  )
  failure_attempt <- .dpprior_new_attempt(
    id = "verification-gate-1", stage = "primary", method = "A2-MN",
    start = list(
      M_selected = as.integer(M),
      M_verification_required = as.integer(M_required)
    ),
    bounds = NULL, control = controls, exit_code = 1L,
    message = message, iterations = 0L,
    evaluations = list(function_count = 1L),
    candidate_parameters = NULL, candidate_objective = NULL,
    elapsed_seconds = 0, warnings = character(),
    error = list(
      class = "dpprior_a2_verification_error", code = code,
      message = message
    ),
    selected = FALSE, reason_code = "no_candidate",
    unavailable = c(
      bounds = "verification gate has no optimizer bounds",
      candidate_parameters = paste(
        "source parameters are quarantined in compatibility and are not",
        "a public canonical candidate"
      ),
      candidate_objective = paste(
        "no canonical objective is available without the required",
        "independent verification order"
      )
    )
  )
  computation <- .dpprior_new_computation(
    request = .a2_schema_settings("A2-MN", controls),
    used = .a2_schema_settings("A2-MN", controls), orders = orders,
    scaling = .a2_schema_scaling(target), attempts = list(failure_attempt),
    candidate_evaluations = list(), selected_candidate_id = NULL,
    selected_attempt_id = NULL, fallback = .dpprior_new_fallback(),
    termination = .dpprior_new_termination(
      code = "no_candidate", message = message, source = "no_candidate"
    ),
    trace = NULL,
    resources = list(
      implementation = "A2-MN verification fail-closed quarantine",
      required_action = action
    )
  )
  verification <- .dpprior_new_verification(
    method = "no_candidate", performed = FALSE, passed = FALSE,
    reason = message, settings = list(), selected_snapshot = NULL,
    verifier_snapshot = NULL, stability = NULL, components = list(),
    invariants = list(no_public_candidate = .dpprior_new_check(
      value = TRUE, reference = TRUE, tolerance = NULL,
      operator = "identical", source = "independent_verifier"
    ))
  )
  legacy_view <- .a2_schema_compatibility_view(
    a, b, J, source_method, source_status, source_usable, source_verified,
    primary, source_attempts, legacy_termination, legacy_diagnostics,
    evidence, initialization, required_action = action
  )
  result <- .dpprior_new_fit(
    mode = "a2_moment", method = "A2-MN", J = J,
    status = "failed", usable = FALSE, verified = FALSE,
    message = message, parameters = NULL, target = list(K = target_K),
    achieved = list(), residuals = list(), tolerances = tolerances,
    computation = computation, verification = verification,
    provenance = .a2_schema_provenance("A2-MN", FALSE, "failed", target_K),
    compatibility = .a2_schema_compatibility(legacy_view, FALSE)
  )
  result <- .dpprior_append_compatibility_v2(
    result,
    aliases = c(converged = "compatibility.views.converged"),
    views = result$compatibility$views,
    deprecations = result$compatibility$deprecations
  )
  condition <- .dpprior_new_condition(
    message,
    classes = c(
      "dpprior_a2_no_candidate", "dpprior_a2_verification_error",
      "dpprior_calibration_error", "dpprior_error", "error"
    ),
    code = code, action = action, result = result
  )
  stop(condition)
}


#' A2-MN Calibrated-Moment Newton Solver
#'
#' Seeks Gamma(a, b) hyperprior parameters that match target moments for the
#' number of clusters K_J under a Dirichlet process prior, and reports
#' convergence only after an independent higher-order check.
#'
#' @param J Integer; sample size (number of observations/sites). Must be >= 2.
#' @param mu_K Numeric; target prior mean \eqn{E[K_J]}. Must satisfy \eqn{1 < \mu_K < J}.
#' @param var_K Numeric; target prior variance \eqn{\mathrm{Var}(K_J)}. Must be positive.
#' @param a0 Numeric or NULL; initial shape parameter. If NULL, computed via
#'   \code{\link{DPprior_a1}}.
#' @param b0 Numeric or NULL; initial rate parameter. If NULL, computed via
#'   \code{\link{DPprior_a1}}.
#' @param tol_F Positive numeric; absolute part of the componentwise moment
#'   tolerance. Retained under its historical name for compatibility.
#' @param tol_rel Non-negative numeric; relative part of the componentwise
#'   tolerance. Each component uses
#'   \code{tol_F + tol_rel * max(abs(target), 1)}.
#' @param tol_step Numeric; stopping tolerance for Newton step size.
#'   Default: 1e-10.
#' @param max_iter Integer; maximum Newton iterations. Default: 20.
#' @param damping Logical; if TRUE, use backtracking line search for
#'   damped Newton updates. Default: TRUE.
#' @param use_fallback Logical; if TRUE, use Nelder-Mead fallback when Newton
#'   fails to converge. Default: TRUE.
#' @param M Integer; number of quadrature nodes for moment computation.
#'   Default: 80.
#' @param M_verify Optional integer; independent verification order. It must be
#'   at least \code{max(2*M, M+40)} and no greater than 512. The required order
#'   is selected automatically when available. If the required order exceeds
#'   512, the function fails closed with a typed no-candidate condition. The
#'   numerical source candidate is then available only in the condition's
#'   non-authoritative compatibility view, with guidance to refit using
#'   \code{M <= 256}.
#' @param verification_abs_tol,verification_rel_tol Non-negative absolute and
#'   relative tolerances for selected-versus-verification order stability. At
#'   least one must be positive.
#' @param verbose Logical; if TRUE, print iteration progress. Default: FALSE.
#'
#' @return A canonical \code{dpprior.result/1} \code{DPprior_fit} object with
#'   authoritative nested components:
#'   \describe{
#'     \item{\code{parameters}}{Gamma shape, rate, and parameterization.}
#'     \item{\code{J}}{Integer; sample size}
#'     \item{\code{target}}{Canonical immutable \code{K} target.}
#'     \item{\code{achieved},\code{residuals}}{Selected-order public
#'       moment evidence.}
#'     \item{\code{method}}{Character; "A2-MN" or "A2-MN+NM" if fallback was used}
#'     \item{\code{status}}{One of \code{"converged"}, \code{"boundary"},
#'       \code{"approximate"}, or \code{"failed"}.}
#'     \item{\code{usable},\code{verified}}{Logical status fields derived from
#'       the numerical and independent-verification gates.}
#'     \item{\code{computation}}{Orders, fixed scaling, typed attempts,
#'       candidate ledger, fallback, termination, and trace.}
#'     \item{\code{verification}}{Independent higher-order recomputation,
#'       componentwise adequacy, and order-stability evidence.}
#'     \item{\code{provenance},\code{compatibility}}{Canonical provenance
#'       plus a quarantined non-authoritative legacy view.}
#'   }
#'   Finite results also retain exact compatibility aliases \code{a}, \code{b},
#'   \code{converged}, \code{iterations}, \code{termination}, \code{fit},
#'   \code{diagnostics}, \code{trace}, and \code{attempts}. If no independent
#'   verifier order is available, the function signals a
#'   \code{dpprior_a2_no_candidate} condition containing a valid parameterless
#'   failed result in its \code{result} field.
#'
#' @details
#' This implements TSMM Stage 2 (A2-MN) from Lee (2026).
#' The A2-MN algorithm uses Newton's method in log-scale to ensure positivity
#' of the Gamma parameters. The Jacobian uses score-function identities and
#' their order-refinement diagnostics. Residual components are scaled by
#' \code{max(abs(target), 1)} and checked individually under an explicit
#' absolute-plus-relative budget; a small aggregate norm cannot hide a failed
#' component.
#'
#' \strong{Algorithm Steps:}
#' \enumerate{
#'   \item Initialize: \eqn{(a_0, b_0)} from A1 closed-form or user-provided
#'   \item Log-parameterize: \eqn{\eta = (\log a, \log b)}
#'   \item For each iteration:
#'     \itemize{
#'       \item Compute moments \eqn{(M_1, V)} and Jacobian \eqn{J_F}
#'       \item Compute natural-unit and scaled component residuals
#'       \item Transform Jacobian to log-scale: \eqn{J_{\log} = J_F \cdot \text{diag}(a, b)}
#'       \item Newton step: \eqn{\Delta = -J_{\log}^{-1} F}
#'       \item Backtracking line search (if damping enabled)
#'       \item Update: \eqn{\eta \leftarrow \eta + \lambda \Delta}
#'     }
#'   \item If needed, run the predeclared Nelder--Mead fallback on
#'     \eqn{(\log a,\log b)}, never on raw positive parameters
#'   \item Independently recompute the selected candidate at a distinct
#'     higher quadrature order; retain the fit-order candidate unchanged
#' }
#'
#' \strong{Termination and status:}
#' \itemize{
#'   \item Optimizer exit, a small step, or line-search exhaustion is never by
#'     itself convergence.
#'   \item \code{converged} requires the selected- and higher-order residual
#'     gates plus the quadrature-stability gate.
#'   \item \code{approximate} retains a finite but unusable candidate when a
#'     solver or available independent-verification gate is unmet.
#'   \item An unavailable required verifier order produces a typed
#'     \code{failed}/no-candidate result; it never exposes the numerical source
#'     candidate as public scientific output.
#'   \item Singular/ill-conditioned Jacobians, line-search failure,
#'     stagnation, fallback outcome, and verification failure have stable
#'     reason codes in \code{computation} and the compatibility audit view.
#' }
#'
#' A1 projection may be requested explicitly to obtain an initialization only.
#' Its projection record is retained in the non-authoritative initialization
#' audit view, while the canonical requested A2 target is immutable and is
#' never silently projected.
#'
#' @seealso
#' \code{\link{DPprior_a1}} for closed-form initialization,
#' \code{\link{moments_with_jacobian}} for Jacobian computation,
#' \code{\link{exact_K_moments}} for moment verification
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @examples
#' # Basic usage
#' fit <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8)
#' print(fit)
#'
#' # Inspect the selected-order match and independent verification
#' achieved <- exact_K_moments(
#'   50, fit$parameters$a, fit$parameters$b
#' )
#' cat(sprintf("Target E[K]=5, Achieved E[K]=%.10f\n", achieved$mean))
#' cat(sprintf("Target Var=8, Achieved Var=%.10f\n", achieved$var))
#'
#' # Compare A1 vs A2 accuracy
#' a1 <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
#' a1_mom <- exact_K_moments(
#'   50, a1$parameters$a, a1$parameters$b
#' )
#' a2_mom <- exact_K_moments(
#'   50, fit$parameters$a, fit$parameters$b
#' )
#' cat(sprintf("A1 mean error: %.6f\n", abs(a1_mom$mean - 5)))
#' cat(sprintf("A2 mean error: %.2e\n", abs(a2_mom$mean - 5)))
#'
#' # View iteration trace (includes step size and Jacobian determinant)
#' head(fit$computation$trace)
#'
#' @family elicitation
#'
#' @export
DPprior_a2_newton <- function(J, mu_K, var_K,
                              a0 = NULL, b0 = NULL,
                              tol_F = .TOL_NEWTON,
                              tol_rel = 1e-8,
                              tol_step = 1e-10,
                              max_iter = 20L,
                              damping = TRUE,
                              use_fallback = TRUE,
                              M = .QUAD_NODES_DEFAULT,
                              M_verify = NULL,
                              verification_abs_tol = 1e-10,
                              verification_rel_tol = 1e-8,
                              verbose = FALSE) {
  assert_valid_J(J)
  mu_K <- .dpprior_validate_scalar(
    mu_K, "mu_K", lower = 1, upper = J,
    lower_open = TRUE, upper_open = TRUE,
    .subclass = "dpprior_moment_target_error"
  )
  var_K <- .dpprior_validate_scalar(
    var_K, "var_K", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_moment_target_error"
  )
  .assert_feasible_K_moments(J, mu_K, var_K)
  tol_F <- .dpprior_validate_scalar(
    tol_F, "tol_F", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_control_error"
  )
  tol_rel <- .dpprior_validate_scalar(
    tol_rel, "tol_rel", lower = 0,
    .subclass = "dpprior_control_error"
  )
  tol_step <- .dpprior_validate_scalar(
    tol_step, "tol_step", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_control_error"
  )
  verification_abs_tol <- .dpprior_validate_scalar(
    verification_abs_tol, "verification_abs_tol", lower = 0,
    .subclass = "dpprior_control_error"
  )
  verification_rel_tol <- .dpprior_validate_scalar(
    verification_rel_tol, "verification_rel_tol", lower = 0,
    .subclass = "dpprior_control_error"
  )
  if (verification_abs_tol == 0 && verification_rel_tol == 0) {
    .dpprior_abort_invalid(
      "at least one verification tolerance must be positive",
      "dpprior_control_error", "verification_abs_tol",
      verification_abs_tol, "positive absolute or relative tolerance",
      "zero_tolerances"
    )
  }
  max_iter <- .dpprior_validate_count(
    max_iter, "max_iter", minimum = 1L,
    .subclass = "dpprior_control_error"
  )
  M <- .dpprior_validate_count(
    M, "M", minimum = 10L, maximum = .QUADRATURE_MAX_NODES,
    .subclass = "dpprior_control_error"
  )
  damping <- .dpprior_validate_control(damping, "damping", type = "logical")
  use_fallback <- .dpprior_validate_control(
    use_fallback, "use_fallback", type = "logical"
  )
  verbose <- .dpprior_validate_control(verbose, "verbose", type = "logical")

  M_required <- .a2_required_verification_order(M)
  verification_available <- M_required <= .QUADRATURE_MAX_NODES
  if (is.null(M_verify)) {
    M_verify <- if (verification_available) M_required else NA_integer_
  } else {
    if (!verification_available) {
      .dpprior_abort_invalid(
        sprintf(
          paste(
            "A2-MN verification is unavailable for M=%d:",
            "required M_verify=%d exceeds the supported ceiling %d"
          ),
          M, M_required, .QUADRATURE_MAX_NODES
        ),
        c("dpprior_a2_verification_error", "dpprior_bounds_error"),
        "M_verify", M_verify,
        sprintf("required order <= %d", .QUADRATURE_MAX_NODES),
        "verification_unavailable"
      )
    }
    M_verify <- .dpprior_validate_count(
      M_verify, "M_verify", minimum = M_required,
      maximum = .QUADRATURE_MAX_NODES,
      .subclass = "dpprior_a2_verification_error"
    )
  }

  if (xor(is.null(a0), is.null(b0))) {
    .dpprior_abort_invalid(
      "a0 and b0 must either both be supplied or both be NULL",
      "dpprior_initialization_error", "a0", a0,
      "both a0 and b0, or neither", "incomplete_initialization"
    )
  }
  target <- c(mean = mu_K, variance = var_K)
  initialization_warning <- NA_character_
  if (is.null(a0)) {
    # A1 projection is allowed solely to obtain a finite starting point.  The
    # requested A2 target remains unchanged and is what every residual uses.
    # A1 infeasibility is only a failed initialization attempt; it is not
    # evidence that the exact A2 target is infeasible.
    a1_error <- NULL
    init <- tryCatch(
      withCallingHandlers(
        DPprior_a1(J, mu_K, var_K, projection = "nearest"),
        warning = function(warning) {
          if (inherits(warning, "dpprior_a1_projection_warning")) {
            initialization_warning <<- conditionMessage(warning)
            invokeRestart("muffleWarning")
          }
        }
      ),
      error = function(error) {
        a1_error <<- error
        NULL
      }
    )
    if (!is.null(init)) {
      a0 <- init$a
      b0 <- init$b
      initialization <- list(
        method = "A1", projection_policy = "nearest_for_start_only",
        target_projection = init$projection %||%
          init$target$projection %||% NULL,
        requested_A2_target = list(mu_K = mu_K, var_K = var_K),
        warning = initialization_warning,
        target_was_not_projected_for_A2 = TRUE,
        attempts = list(list(
          method = "A1_nearest_start", status = "selected",
          reason_code = "a1_initialization_available",
          projection_policy = "nearest_for_start_only",
          candidate = c(a = unname(a0), b = unname(b0)),
          warning = initialization_warning, error = NA_character_
        ))
      )
    } else {
      grid <- .a2_select_initialization_grid(
        J, target, M, tol_F, tol_rel
      )
      a0 <- exp(grid$eta[["log_shape"]])
      b0 <- exp(grid$eta[["log_rate"]])
      initialization <- list(
        method = "fixed_log_parameter_grid",
        projection_policy = "A1_failed_no_A2_projection",
        target_projection = NULL,
        requested_A2_target = list(mu_K = mu_K, var_K = var_K),
        warning = initialization_warning,
        target_was_not_projected_for_A2 = TRUE,
        attempts = list(
          list(
            method = "A1_nearest_start", status = "failed",
            reason_code = a1_error$code %||%
              "a1_initialization_failed",
            projection_policy = "nearest_for_start_only",
            candidate = NULL, warning = initialization_warning,
            error = conditionMessage(a1_error),
            error_classes = class(a1_error)
          ),
          list(
            method = "fixed_log_parameter_grid", status = "selected",
            reason_code = "fixed_grid_candidate_selected",
            selection_rule = grid$selection_rule,
            selected_index = grid$selected_index,
            candidate = c(a = unname(a0), b = unname(b0)),
            candidate_objective = grid$selected_objective,
            candidates = grid$candidates,
            warning = NA_character_, error = NA_character_
          )
        )
      )
    }
  } else {
    initialization <- list(
      method = "user", projection_policy = "not_applicable",
      target_projection = NULL,
      requested_A2_target = list(mu_K = mu_K, var_K = var_K),
      warning = NA_character_, target_was_not_projected_for_A2 = TRUE,
      attempts = list(list(
        method = "user_supplied", status = "selected",
        reason_code = "user_initialization",
        candidate = c(a = unname(a0), b = unname(b0)),
        warning = NA_character_, error = NA_character_
      ))
    )
  }
  a0 <- .dpprior_validate_scalar(
    a0, "a0", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_initialization_error"
  )
  b0 <- .dpprior_validate_scalar(
    b0, "b0", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_initialization_error"
  )
  eta_start <- c(log_shape = log(a0), log_rate = log(b0))
  if (any(eta_start < .LOG_BOUNDS_DEFAULT[1L]) ||
      any(eta_start > .LOG_BOUNDS_DEFAULT[2L])) {
    .dpprior_abort_invalid(
      "A2-MN initialization lies outside the declared log-parameter bounds",
      "dpprior_initialization_error", "a0/b0", c(a0 = a0, b0 = b0),
      sprintf(
        "log parameters in [%g, %g]",
        .LOG_BOUNDS_DEFAULT[1L], .LOG_BOUNDS_DEFAULT[2L]
      ),
      "initialization_outside_bounds"
    )
  }
  if (isTRUE(verbose)) {
    cat("A2-MN calibrated-moment solver\n")
    cat(sprintf("Target: E[K]=%.6g, Var(K)=%.6g\n", mu_K, var_K))
    cat(sprintf(
      "Residual tolerance: absolute %.3g + relative %.3g * scale\n",
      tol_F, tol_rel
    ))
    cat(sprintf(
      "Fit M=%d; independent verification M=%s\n",
      M, if (is.na(M_verify)) "unavailable" else as.character(M_verify)
    ))
  }
  primary <- .a2_run_newton(
    J, target, eta_start, M, tol_F, tol_rel, tol_step,
    max_iter, damping, verbose
  )
  selected <- primary$best
  attempts <- list(primary$attempt)
  selected_method <- "scaled_log_newton"
  fallback_attempted <- FALSE
  fallback_selected <- FALSE
  solver_acceptable <- primary$solved
  if (!primary$solved && isTRUE(use_fallback)) {
    fallback_attempted <- TRUE
    fallback <- .a2_run_fallback(
      J, target, selected$eta, M, tol_F, tol_rel,
      selected$objective
    )
    attempts[[2L]] <- fallback$attempt
    if (isTRUE(fallback$selected)) {
      selected <- fallback$candidate
      selected_method <- "nelder_mead_log"
      fallback_selected <- TRUE
    }
    solver_acceptable <- isTRUE(fallback$success) && fallback_selected
  }

  a_star <- unname(exp(selected$eta[1L]))
  b_star <- unname(exp(selected$eta[2L]))
  if (is.null(selected$moments) || is.null(selected$residual)) {
    selected$moments <- tryCatch(
      .a2_moments_at_order(J, a_star, b_star, M),
      error = function(error) c(mean = NA_real_, variance = NA_real_)
    )
    selected$residual <- .a2_residual_contract(
      selected$moments, target, tol_F, tol_rel
    )
  }
  # The historical verification verdict is retained only in the explicitly
  # non-authoritative compatibility view.  Canonical decision evidence is
  # freshly recomputed below with the frozen truth tolerances.
  legacy_verification <- .a2_verify_candidate(
    J, a_star, b_star, M,
    if (is.na(M_verify)) NULL else M_verify,
    target, tol_F, tol_rel,
    verification_abs_tol, verification_rel_tol
  )
  boundary_distance <- min(
    selected$eta - .LOG_BOUNDS_DEFAULT[1L],
    .LOG_BOUNDS_DEFAULT[2L] - selected$eta
  )
  at_boundary <- is.finite(boundary_distance) && boundary_distance <= 1e-6
  legacy_status <- if (identical(legacy_verification$status, "failed")) {
    "failed"
  } else if (solver_acceptable && isTRUE(legacy_verification$passed) &&
             at_boundary) {
    "boundary"
  } else if (solver_acceptable && isTRUE(legacy_verification$passed)) {
    "converged"
  } else {
    "approximate"
  }
  legacy_usable <- legacy_status %in% c("converged", "boundary")
  legacy_verified <- isTRUE(legacy_verification$passed) && legacy_usable
  legacy_termination <- if (fallback_selected && solver_acceptable) {
    "nelder_mead"
  } else if (primary$solved) {
    "residual"
  } else if (identical(primary$reason, "stagnation_small_step")) {
    "step"
  } else {
    "max_iter"
  }
  method <- if (fallback_selected) "A2-MN+NM" else "A2-MN"
  controls <- .a2_schema_controls(
    max_iter, damping, use_fallback, tol_step
  )
  tolerances <- .a2_schema_tolerances(
    tol_F, tol_rel, tol_step,
    verification_abs_tol, verification_rel_tol
  )
  fresh_evidence <- .a2_schema_fresh_evidence(
    J, a_star, b_star, M,
    if (is.na(M_verify)) NULL else M_verify,
    target, tolerances
  )
  legacy_diagnostics <- list(
    a0 = a0, b0 = b0, tol_F = tol_F, tol_rel = tol_rel,
    tol_step = tol_step, M = M,
    M_verify = if (is.na(M_verify)) NA_integer_ else as.integer(M_verify),
    M_verify_required = M_required,
    verification_available = verification_available,
    fallback_used = fallback_selected,
    fallback_attempted = fallback_attempted,
    conditioning = primary$conditioning,
    derivative_status = primary$derivative_status,
    quasi_improper = a_star < 0.1,
    boundary_distance_log_scale = boundary_distance
  )

  if (!isTRUE(fresh_evidence$ok)) {
    action <- if (identical(
      fresh_evidence$code, "verification_order_exceeds_ceiling"
    )) {
      "refit_with_smaller_M_at_or_below_256"
    } else {
      "refit_after_resolving_independent_verification_failure"
    }
    failure_message <- paste(
      fresh_evidence$message,
      "No public finite candidate is exposed; the source candidate is",
      "retained only in a non-authoritative compatibility view.",
      "Required action:", action
    )
    .a2_schema_failed_result(
      J = J, target = target, a = a_star, b = b_star,
      source_method = method, source_status = legacy_status,
      source_usable = legacy_usable, source_verified = legacy_verified,
      message = failure_message,
      code = paste0("a2_", fresh_evidence$code), action = action,
      evidence = fresh_evidence, M = M, M_required = M_required,
      controls = controls, tolerances = tolerances,
      source_attempts = attempts, primary = primary,
      initialization = initialization,
      legacy_termination = legacy_termination,
      legacy_diagnostics = legacy_diagnostics
    )
  }

  status <- if (solver_acceptable && isTRUE(fresh_evidence$passed) &&
                at_boundary) {
    "boundary"
  } else if (solver_acceptable && isTRUE(fresh_evidence$passed)) {
    "converged"
  } else {
    "approximate"
  }
  message <- switch(
    status,
    converged = paste(
      "scaled component residuals and independent higher-order",
      "verification passed"
    ),
    boundary = "verified candidate is on the declared log-parameter boundary",
    approximate = sprintf(
      "finite candidate retained, but convergence was withheld (%s; %s)",
      if (solver_acceptable) {
        "solver tolerance passed"
      } else {
        "solver tolerance or exit failed"
      },
      fresh_evidence$reason
    )
  )

  if (isTRUE(verbose)) {
    cat(sprintf(
      "Final status: %s; method: %s; verification: %s\n",
      status, selected_method, fresh_evidence$reason
    ))
  }
  .a2_schema_finite_result(
    J = J, target = target, a = a_star, b = b_star,
    method = method, status = status, message = message,
    evidence = fresh_evidence, M = M, M_verify = M_verify,
    M_required = M_required, controls = controls, tolerances = tolerances,
    source_attempts = attempts, primary = primary,
    fallback_attempted = fallback_attempted,
    fallback_selected = fallback_selected,
    fallback_reason = if (fallback_attempted) primary$reason else "not_attempted",
    initialization = initialization,
    legacy_status = legacy_status, legacy_usable = legacy_usable,
    legacy_verified = legacy_verified,
    legacy_termination = legacy_termination,
    legacy_diagnostics = legacy_diagnostics
  )
}


# =============================================================================
# Verification Functions
# =============================================================================

# Every R11 consumer crosses the same fail-closed boundary.  Returning a small
# classless view keeps downstream helper logic from accidentally consulting
# deprecated aliases or the explicitly non-authoritative compatibility layer.
.a2_consumer_canonical_raw <- function(
    fit, expected_mode, J, mu_K, var_K) {
  validated <- .dpprior_require_schema(
    fit, kind = "fit", allow_legacy = FALSE
  )
  raw <- unclass(validated)
  .dpprior_schema_require(
    identical(raw[["mode", exact = TRUE]], expected_mode),
    "consumer_mode", "result.mode", expected_mode,
    raw[["mode", exact = TRUE]]
  )
  .dpprior_schema_require(
    identical(raw[["J", exact = TRUE]], as.integer(J)),
    "consumer_J", "result.J", "identity with the helper request",
    raw[["J", exact = TRUE]]
  )
  target <- raw[["target", exact = TRUE]][["K", exact = TRUE]]
  target_raw <- unclass(target)
  used <- target_raw[["used", exact = TRUE]]
  recorded_target <- c(
    mean = used[["mu_K", exact = TRUE]],
    variance = used[["var_K", exact = TRUE]]
  )
  expected_target <- c(
    mean = as.numeric(mu_K), variance = as.numeric(var_K)
  )
  .dpprior_schema_require(
    identical(recorded_target, expected_target),
    "consumer_target_identity", "result.target.K.used",
    "exact identity with the helper moment request", recorded_target
  )
  raw
}


.a2_consumer_view <- function(fit, J, mu_K, var_K) {
  raw <- .a2_consumer_canonical_raw(
    fit, "a2_moment", J, mu_K, var_K
  )
  parameters <- raw[["parameters", exact = TRUE]]
  .dpprior_schema_require(
    !is.null(parameters), "consumer_parameters", "result.parameters",
    "a finite canonical A2-MN candidate", NULL
  )
  achieved <- raw[["achieved", exact = TRUE]][["K", exact = TRUE]]
  residuals <- raw[["residuals", exact = TRUE]][["K", exact = TRUE]]
  computation <- raw[["computation", exact = TRUE]]
  termination <- computation[["termination", exact = TRUE]]
  residual_vector <- c(
    mean = residuals[["mean", exact = TRUE]],
    variance = residuals[["variance", exact = TRUE]]
  )
  status <- raw[["status", exact = TRUE]]
  usable <- raw[["usable", exact = TRUE]]
  verified <- raw[["verified", exact = TRUE]]
  list(
    status = status, usable = usable, verified = verified,
    decision_ready = status %in% c("converged", "boundary") &&
      usable && verified,
    a = parameters[["a", exact = TRUE]],
    b = parameters[["b", exact = TRUE]],
    mean = achieved[["mean", exact = TRUE]],
    variance = achieved[["variance", exact = TRUE]],
    residual = sqrt(sum(residual_vector^2)),
    residual_components = residual_vector,
    iterations = termination[["iterations", exact = TRUE]],
    termination_code = termination[["code", exact = TRUE]],
    termination_source = termination[["source", exact = TRUE]],
    trace = computation[["trace", exact = TRUE]]
  )
}


.a2_consumer_a1_parameters <- function(fit, J, mu_K, var_K) {
  raw <- .a2_consumer_canonical_raw(
    fit, "a1_proxy", J, mu_K, var_K
  )
  parameters <- raw[["parameters", exact = TRUE]]
  .dpprior_schema_require(
    !is.null(parameters), "consumer_parameters", "result.parameters",
    "finite canonical A1 proxy parameters", NULL
  )
  c(
    a = parameters[["a", exact = TRUE]],
    b = parameters[["b", exact = TRUE]]
  )
}

#' Verify A2-MN Moment Matching
#'
#' Tests that the A2-MN solver achieves exact moment matching.
#'
#' @param J Integer; sample size.
#' @param mu_K Numeric; target mean.
#' @param var_K Numeric; target variance.
#' @param tol Numeric; tolerance for verification.
#' @param verbose Logical; if TRUE, print results.
#'
#' @return Logical; TRUE if verification passes.
#'
#' @examples
#' \dontrun{
#' verify_a2_moment_matching(J = 50, mu_K = 5, var_K = 8)
#'
#' }
#' @keywords internal
verify_a2_moment_matching <- function(J, mu_K, var_K, tol = 1e-6, verbose = TRUE) {
  fit <- DPprior_a2_newton(J, mu_K, var_K, verbose = FALSE)
  view <- .a2_consumer_view(fit, J, mu_K, var_K)
  achieved_mean <- view[["mean", exact = TRUE]]
  achieved_variance <- view[["variance", exact = TRUE]]
  pass <- view[["decision_ready", exact = TRUE]] &&
    abs(achieved_mean - mu_K) < tol &&
    abs(achieved_variance - var_K) < tol

  if (isTRUE(verbose)) {
    cat(sprintf("A2-MN Moment Matching Verification\n"))
    cat(strrep("-", 50), "\n")
    cat(sprintf("Target: E[K]=%.4f, Var(K)=%.4f\n", mu_K, var_K))
    cat(sprintf("Achieved: E[K]=%.10f, Var(K)=%.10f\n",
                achieved_mean, achieved_variance))
    cat(sprintf("Mean error: %.2e\n", abs(achieved_mean - mu_K)))
    cat(sprintf("Var error: %.2e\n", abs(achieved_variance - var_K)))
    cat(sprintf(
      "Termination: %s (%s)\n",
      view[["termination_code", exact = TRUE]],
      view[["termination_source", exact = TRUE]]
    ))
    cat(sprintf("Status: %s\n", if (pass) "PASS" else "FAIL"))
  }

  invisible(pass)
}


#' Compare A1 vs A2 Accuracy
#'
#' Compares the accuracy of A1 closed-form and A2 Newton methods.
#'
#' @param J Integer; sample size.
#' @param mu_K Numeric; target mean.
#' @param var_K Numeric; target variance.
#' @param verbose Logical; if TRUE, print comparison.
#'
#' @return A list with A1 and A2 results and error comparison.
#'
#' @examples
#' compare_a1_a2(J = 50, mu_K = 5, var_K = 8)
#'
#' @export
compare_a1_a2 <- function(J, mu_K, var_K, verbose = TRUE) {
  # A1 solution
  a1 <- DPprior_a1(J, mu_K, var_K)
  a1_parameters <- .a2_consumer_a1_parameters(
    a1, J, mu_K, var_K
  )
  a1_a <- a1_parameters[["a", exact = TRUE]]
  a1_b <- a1_parameters[["b", exact = TRUE]]
  a1_mom <- exact_K_moments(J, a1_a, a1_b)
  a1_mean <- a1_mom[["mean", exact = TRUE]]
  a1_variance <- a1_mom[["var", exact = TRUE]]
  a1_residual <- sqrt(
    (a1_mean - mu_K)^2 + (a1_variance - var_K)^2
  )

  # A2 solution
  a2 <- DPprior_a2_newton(J, mu_K, var_K, verbose = FALSE)
  a2_view <- .a2_consumer_view(a2, J, mu_K, var_K)
  a2_residual <- a2_view[["residual", exact = TRUE]]

  improvement_ratio <- a1_residual / max(a2_residual, 1e-15)

  if (isTRUE(verbose)) {
    cat(sprintf("A1 vs A2 Comparison (J=%d, mu_K=%.2f, var_K=%.2f)\n", J, mu_K, var_K))
    cat(strrep("-", 60), "\n")
    cat(sprintf("%-20s %12s %12s\n", "", "A1", "A2"))
    cat(strrep("-", 60), "\n")
    cat(sprintf(
      "%-20s %12.6f %12.6f\n", "Shape (a)", a1_a,
      a2_view[["a", exact = TRUE]]
    ))
    cat(sprintf(
      "%-20s %12.6f %12.6f\n", "Rate (b)", a1_b,
      a2_view[["b", exact = TRUE]]
    ))
    cat(sprintf(
      "%-20s %12.6f %12.10f\n", "E[K] achieved", a1_mean,
      a2_view[["mean", exact = TRUE]]
    ))
    cat(sprintf(
      "%-20s %12.6f %12.10f\n", "Var achieved", a1_variance,
      a2_view[["variance", exact = TRUE]]
    ))
    cat(sprintf(
      "%-20s %12.6f %12.2e\n", "Residual", a1_residual, a2_residual
    ))
    cat(strrep("-", 60), "\n")
    cat(sprintf("Improvement ratio: %.0fx\n", improvement_ratio))
  }

  invisible(list(
    a1 = list(a = a1_a, b = a1_b, mean = a1_mean, var = a1_variance,
              residual = a1_residual),
    a2 = list(
      a = a2_view[["a", exact = TRUE]],
      b = a2_view[["b", exact = TRUE]],
      mean = a2_view[["mean", exact = TRUE]],
      var = a2_view[["variance", exact = TRUE]],
      residual = a2_residual
    ),
    improvement_ratio = improvement_ratio
  ))
}


#' Run All A2-MN Verification Tests
#'
#' Comprehensive verification suite for the A2-MN Newton solver.
#'
#' @param verbose Logical; if TRUE, print detailed results.
#'
#' @return Logical; TRUE if all tests pass.
#'
#' @examples
#' \dontrun{
#' verify_a2_all()
#'
#' }
#' @keywords internal
verify_a2_all <- function(verbose = TRUE) {
  if (isTRUE(verbose)) {
    cat(strrep("=", 70), "\n")
    cat("Module 11: A2-MN Newton Solver - Full Verification Suite\n")
    cat(strrep("=", 70), "\n\n")
  }

  all_pass <- TRUE

  # Test 1: Basic convergence
  if (isTRUE(verbose)) {
    cat("[Test 1] Basic convergence (J=50, mu_K=5, var_K=8)\n")
    cat(strrep("-", 50), "\n")
  }

  fit <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8, verbose = verbose)
  fit_view <- .a2_consumer_view(fit, 50, 5, 8)
  test1_pass <- fit_view[["decision_ready", exact = TRUE]] &&
    abs(fit_view[["mean", exact = TRUE]] - 5) < 1e-6 &&
    abs(fit_view[["variance", exact = TRUE]] - 8) < 1e-6

  if (isTRUE(verbose)) {
    cat(sprintf("\nResult: %s\n\n", if (test1_pass) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test1_pass

  # Test 2: Fast convergence
  if (isTRUE(verbose)) {
    cat("[Test 2] Fast convergence (< 10 iterations)\n")
    cat(strrep("-", 50), "\n")
  }

  fit <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8, verbose = FALSE)
  fit_view <- .a2_consumer_view(fit, 50, 5, 8)
  test2_pass <- fit_view[["iterations", exact = TRUE]] < 10

  if (isTRUE(verbose)) {
    cat(sprintf("Iterations: %d\n", fit_view[["iterations", exact = TRUE]]))
    cat(sprintf(
      "Termination: %s (%s)\n",
      fit_view[["termination_code", exact = TRUE]],
      fit_view[["termination_source", exact = TRUE]]
    ))
    cat(sprintf("Result: %s\n\n", if (test2_pass) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test2_pass

  # Test 3: A2 corrects A1 error
  if (isTRUE(verbose)) {
    cat("[Test 3] A2 improves over A1\n")
    cat(strrep("-", 50), "\n")
  }

  comparison <- compare_a1_a2(J = 50, mu_K = 5, var_K = 8, verbose = verbose)
  test3_pass <- comparison[["a2", exact = TRUE]][["residual", exact = TRUE]] <
    comparison[["a1", exact = TRUE]][["residual", exact = TRUE]]

  if (isTRUE(verbose)) {
    cat(sprintf("\nResult: %s\n\n", if (test3_pass) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test3_pass

  # Test 4: Various target scenarios
  if (isTRUE(verbose)) {
    cat("[Test 4] Various target scenarios\n")
    cat(strrep("-", 50), "\n")
  }

  test_cases <- list(
    list(J = 30, mu_K = 3, var_K = 5),
    list(J = 50, mu_K = 10, var_K = 15),
    list(J = 100, mu_K = 5, var_K = 8),
    list(J = 50, mu_K = 25, var_K = 50)
  )

  for (tc in test_cases) {
    fit <- DPprior_a2_newton(
      tc[["J", exact = TRUE]], tc[["mu_K", exact = TRUE]],
      tc[["var_K", exact = TRUE]], verbose = FALSE
    )
    fit_view <- .a2_consumer_view(
      fit, tc[["J", exact = TRUE]], tc[["mu_K", exact = TRUE]],
      tc[["var_K", exact = TRUE]]
    )
    case_pass <- fit_view[["decision_ready", exact = TRUE]] &&
      fit_view[["residual", exact = TRUE]] < 1e-6
    status <- if (case_pass) "PASS" else "FAIL"
    if (isTRUE(verbose)) {
      cat(sprintf("  J=%3d, mu_K=%2d, var_K=%2d: %s (iter=%d, term=%s, res=%.2e)\n",
                  tc[["J", exact = TRUE]], tc[["mu_K", exact = TRUE]],
                  tc[["var_K", exact = TRUE]], status,
                  fit_view[["iterations", exact = TRUE]],
                  fit_view[["termination_code", exact = TRUE]],
                  fit_view[["residual", exact = TRUE]]))
    }
    all_pass <- all_pass && case_pass
  }

  if (isTRUE(verbose)) {
    cat("\n")
  }

  # Test 5: Edge case with high VIF (quasi-improper prior)
  if (isTRUE(verbose)) {
    cat("[Test 5] Edge case with high VIF (quasi-improper prior)\n")
    cat(strrep("-", 50), "\n")
  }

  # This case requires very small a
  fit <- DPprior_a2_newton(J = 50, mu_K = 3, var_K = 10, verbose = FALSE)
  fit_view <- .a2_consumer_view(fit, 50, 3, 10)
  test5_pass <- fit_view[["decision_ready", exact = TRUE]] &&
    fit_view[["residual", exact = TRUE]] < 1e-5

  if (isTRUE(verbose)) {
    cat(sprintf("Target: mu_K=3, var_K=10 (VIF=%.1f)\n", 10 / 2))
    cat(sprintf(
      "Decision ready: %s, Iterations: %d\n",
      fit_view[["decision_ready", exact = TRUE]],
      fit_view[["iterations", exact = TRUE]]
    ))
    cat(sprintf(
      "Termination: %s (%s)\n",
      fit_view[["termination_code", exact = TRUE]],
      fit_view[["termination_source", exact = TRUE]]
    ))
    cat(sprintf(
      "Achieved: a=%.6f (quasi-improper: %s)\n",
      fit_view[["a", exact = TRUE]], fit_view[["a", exact = TRUE]] < 0.1
    ))
    cat(sprintf("Residual: %.2e\n", fit_view[["residual", exact = TRUE]]))
    cat(sprintf("Result: %s\n\n", if (test5_pass) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test5_pass

  # Test 6: Termination field consistency
  if (isTRUE(verbose)) {
    cat("[Test 6] Termination field consistency\n")
    cat(strrep("-", 50), "\n")
  }

  fit <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8, verbose = FALSE)
  fit_view <- .a2_consumer_view(fit, 50, 5, 8)
  test6_pass <- fit_view[["termination_code", exact = TRUE]] %in%
    c("converged", "boundary", "approximate")

  if (isTRUE(verbose)) {
    cat(sprintf(
      "Termination field: '%s'\n",
      fit_view[["termination_code", exact = TRUE]]
    ))
    cat(sprintf("Valid termination: %s\n", test6_pass))
    cat(sprintf("Result: %s\n\n", if (test6_pass) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test6_pass

  # Test 7: Trace contains required diagnostics
  if (isTRUE(verbose)) {
    cat("[Test 7] Trace contains required diagnostics\n")
    cat(strrep("-", 50), "\n")
  }

  required_cols <- c("iter", "a", "b", "M1", "V", "residual", "step", "det_Jlog")
  trace <- fit_view[["trace", exact = TRUE]]
  test7_pass <- all(required_cols %in% names(trace))

  if (isTRUE(verbose)) {
    cat(sprintf("Required columns: %s\n", paste(required_cols, collapse = ", ")))
    cat(sprintf("Present columns:  %s\n", paste(names(trace), collapse = ", ")))
    cat(sprintf("Result: %s\n\n", if (test7_pass) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test7_pass

  # Summary
  if (isTRUE(verbose)) {
    cat(strrep("=", 70), "\n")
    cat(sprintf("Overall Result: %s\n",
                if (all_pass) "ALL TESTS PASSED" else "SOME TESTS FAILED"))
    cat(strrep("=", 70), "\n")
  }

  invisible(all_pass)
}
