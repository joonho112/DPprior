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
# Author: JoonHo Lee
# Date: December 2025
# Part of: DPprior R Package
# Reference: Lee (2026), Section 3.2
# Dependencies: Modules 00, 02, 03, 05, 07, 10
# =============================================================================


# =============================================================================
# Main A2-MN Newton Solver
# =============================================================================

#' A2-MN Exact-Moment Newton Solver
#'
#' Finds Gamma(a, b) hyperprior parameters that exactly match target moments
#' for the number of clusters K_J under a Dirichlet process prior.
#'
#' @param J Integer; sample size (number of observations/sites). Must be >= 2.
#' @param mu_K Numeric; target prior mean \eqn{E[K_J]}. Must satisfy \eqn{1 < \mu_K < J}.
#' @param var_K Numeric; target prior variance \eqn{\mathrm{Var}(K_J)}. Must be positive.
#' @param a0 Numeric or NULL; initial shape parameter. If NULL, computed via
#'   \code{\link{DPprior_a1}}.
#' @param b0 Numeric or NULL; initial rate parameter. If NULL, computed via
#'   \code{\link{DPprior_a1}}.
#' @param tol_F Numeric; stopping tolerance for the residual norm
#'   ||F|| = sqrt((M1 - mu_K)^2 + (V - var_K)^2). Default: 1e-8.
#' @param tol_step Numeric; stopping tolerance for Newton step size.
#'   Default: 1e-10.
#' @param max_iter Integer; maximum Newton iterations. Default: 20.
#' @param damping Logical; if TRUE, use backtracking line search for
#'   damped Newton updates. Default: TRUE.
#' @param use_fallback Logical; if TRUE, use Nelder-Mead fallback when Newton
#'   fails to converge. Default: TRUE.
#' @param M Integer; number of quadrature nodes for moment computation.
#'   Default: 80.
#' @param verbose Logical; if TRUE, print iteration progress. Default: FALSE.
#'
#' @return A \code{DPprior_fit} object (S3 class) with components:
#'   \describe{
#'     \item{\code{a}}{Numeric; optimal shape parameter}
#'     \item{\code{b}}{Numeric; optimal rate parameter}
#'     \item{\code{J}}{Integer; sample size}
#'     \item{\code{target}}{List with \code{mu_K}, \code{var_K}, and \code{type = "moments"}}
#'     \item{\code{method}}{Character; "A2-MN" or "A2-MN+NM" if fallback was used}
#'     \item{\code{status}}{Character; convergence status ("success", "stagnated", "max_iter")}
#'     \item{\code{converged}}{Logical; whether algorithm converged to target tolerance}
#'     \item{\code{iterations}}{Integer; number of iterations}
#'     \item{\code{termination}}{Character; reason for termination ("residual", "step", "max_iter", "nelder_mead")}
#'     \item{\code{fit}}{List with achieved \code{mu_K}, \code{var_K}, and \code{residual}}
#'     \item{\code{diagnostics}}{List with diagnostic information}
#'     \item{\code{trace}}{Data frame with iteration history}
#'   }
#'
#' @details
#' This implements TSMM Stage 2 (A2-MN) from Lee (2026).
#' The A2-MN algorithm uses Newton's method in log-scale to ensure positivity
#' of the Gamma parameters. The Jacobian is computed exactly using score
#' function identities (Lee, 2026, Section 3.2, Corollary 1), avoiding finite difference
#' approximations.
#'
#' \strong{Algorithm Steps:}
#' \enumerate{
#'   \item Initialize: \eqn{(a_0, b_0)} from A1 closed-form or user-provided
#'   \item Log-parameterize: \eqn{\eta = (\log a, \log b)}
#'   \item For each iteration:
#'     \itemize{
#'       \item Compute moments \eqn{(M_1, V)} and Jacobian \eqn{J_F}
#'       \item Compute residual \eqn{F = (M_1 - \mu_K, V - \sigma^2_K)}
#'       \item Transform Jacobian to log-scale: \eqn{J_{\log} = J_F \cdot \text{diag}(a, b)}
#'       \item Newton step: \eqn{\Delta = -J_{\log}^{-1} F}
#'       \item Backtracking line search (if damping enabled)
#'       \item Update: \eqn{\eta \leftarrow \eta + \lambda \Delta}
#'     }
#'   \item Return \eqn{(a, b) = \exp(\eta)}
#' }
#'
#' \strong{Convergence Behavior:}
#' \itemize{
#'   \item For typical targets, converges in 3-8 iterations
#'   \item For targets requiring very small \eqn{a} (quasi-improper priors),
#'         convergence may be slower; the Nelder-Mead fallback handles these cases
#'   \item Machine-precision accuracy (residual < 1e-10) is typically achieved
#' }
#'
#' \strong{Termination Conditions:}
#' \itemize{
#'   \item \code{residual}: ||F|| < tol_F (success)
#'   \item \code{step}: ||delta|| < tol_step (stagnation - may not have converged)
#'   \item \code{max_iter}: maximum iterations reached
#'   \item \code{nelder_mead}: Nelder-Mead fallback succeeded
#' }
#'
#' \strong{Error Handling:}
#' \itemize{
#'   \item If Jacobian becomes singular, falls back to gradient descent
#'   \item If Newton fails after \code{max_iter}, uses Nelder-Mead if enabled
#'   \item Warns for quasi-improper priors (\eqn{a < 0.1})
#' }
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
#' # Verify exact moment matching
#' achieved <- exact_K_moments(50, fit$a, fit$b)
#' cat(sprintf("Target E[K]=5, Achieved E[K]=%.10f\n", achieved$mean))
#' cat(sprintf("Target Var=8, Achieved Var=%.10f\n", achieved$var))
#'
#' # Compare A1 vs A2 accuracy
#' a1 <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
#' a1_mom <- exact_K_moments(50, a1$a, a1$b)
#' a2_mom <- exact_K_moments(50, fit$a, fit$b)
#' cat(sprintf("A1 mean error: %.6f\n", abs(a1_mom$mean - 5)))
#' cat(sprintf("A2 mean error: %.2e\n", abs(a2_mom$mean - 5)))
#'
#' # View iteration trace (includes step size and Jacobian determinant)
#' head(fit$trace)
#'
#' @family elicitation
#'
#' @export
DPprior_a2_newton <- function(J, mu_K, var_K,
                              a0 = NULL, b0 = NULL,
                              tol_F = .TOL_NEWTON,
                              tol_step = 1e-10,
                              max_iter = 20L,
                              damping = TRUE,
                              use_fallback = TRUE,
                              M = .QUAD_NODES_DEFAULT,
                              verbose = FALSE) {

  # -------------------------------------------------------------------------
  # Input Validation
  # -------------------------------------------------------------------------
  assert_valid_J(J)

  if (!is.numeric(mu_K) || length(mu_K) != 1L || !is.finite(mu_K)) {
    stop("mu_K must be a finite numeric scalar", call. = FALSE)
  }
  if (mu_K <= 1) {
    stop("mu_K must be > 1 (at least one cluster is always present)", call. = FALSE)
  }
  if (mu_K > J) {
    stop("mu_K must be <= J", call. = FALSE)
  }

  if (!is.numeric(var_K) || length(var_K) != 1L || !is.finite(var_K) ||
      var_K <= 0) {
    stop("var_K must be a positive finite numeric scalar", call. = FALSE)
  }

  if (!is.numeric(tol_F) || length(tol_F) != 1L || tol_F <= 0) {
    stop("tol_F must be a positive numeric scalar", call. = FALSE)
  }

  if (!is.numeric(tol_step) || length(tol_step) != 1L || tol_step <= 0) {
    stop("tol_step must be a positive numeric scalar", call. = FALSE)
  }

  max_iter <- as.integer(max_iter)
  if (max_iter < 1L) {
    stop("max_iter must be a positive integer", call. = FALSE)
  }

  M <- as.integer(M)
  if (M < 10L) {
    stop("M must be at least 10", call. = FALSE)
  }

  # -------------------------------------------------------------------------
  # Initialization
  # -------------------------------------------------------------------------
  if (is.null(a0) || is.null(b0)) {
    init <- DPprior_a1(J, mu_K, var_K)
    a0 <- init$a
    b0 <- init$b
    if (isTRUE(verbose)) {
      cat(sprintf("A2-MN Newton Solver\n"))
      cat(sprintf("Target: E[K]=%.4f, Var(K)=%.4f\n", mu_K, var_K))
      cat(sprintf("A1 initialization: a0=%.6f, b0=%.6f\n", a0, b0))
      cat(strrep("-", 80), "\n")
      cat(sprintf("%4s | %10s | %10s | %10s | %10s | %10s | %8s | %10s\n",
                  "Iter", "a", "b", "E[K]", "Var(K)", "||F||", "step", "det(J)"))
      cat(strrep("-", 80), "\n")
    }
  } else {
    if (!is.numeric(a0) || a0 <= 0) {
      stop("a0 must be positive", call. = FALSE)
    }
    if (!is.numeric(b0) || b0 <= 0) {
      stop("b0 must be positive", call. = FALSE)
    }
  }

  # Log-parameterization for positivity constraint
  eta <- c(log(a0), log(b0))

  # History storage (pre-allocate data.frame structure)
  history <- data.frame(
    iter = integer(0),
    a = numeric(0),
    b = numeric(0),
    M1 = numeric(0),
    V = numeric(0),
    residual = numeric(0),
    step = numeric(0),
    det_Jlog = numeric(0),
    stringsAsFactors = FALSE
  )

  converged <- FALSE
  termination <- "max_iter"
  status <- "max_iter"

  # -------------------------------------------------------------------------
  # Newton Iteration
  # -------------------------------------------------------------------------
  for (iter in seq_len(max_iter)) {
    # Current parameters
    a <- exp(eta[1L])
    b <- exp(eta[2L])

    # Compute moments and Jacobian
    mj <- moments_with_jacobian(J, a, b, M)

    # Residual vector F = (M1 - mu_K, V - var_K)
    F_val <- c(mj$mean - mu_K, mj$var - var_K)
    residual <- sqrt(sum(F_val^2))

    # Transform Jacobian to log-scale
    # Chain rule: d/d(log a) = d/da * a, d/d(log b) = d/db * b
    # J_log = J_F %*% diag(a, b)
    J_log <- mj$jacobian %*% diag(c(a, b))
    det_Jlog <- as.numeric(det(J_log))

    # Record history (step will be updated after backtracking)
    history <- rbind(
      history,
      data.frame(
        iter = iter, a = a, b = b,
        M1 = mj$mean, V = mj$var,
        residual = residual,
        step = NA_real_, det_Jlog = det_Jlog,
        stringsAsFactors = FALSE
      )
    )

    # Check convergence on residual
    if (residual < tol_F) {
      converged <- TRUE
      termination <- "residual"
      status <- "success"
      history$step[nrow(history)] <- NA_real_  # No step taken
      if (isTRUE(verbose)) {
        cat(sprintf("%4d | %10.6f | %10.6f | %10.6f | %10.6f | %10.2e | %8s | %10.2e\n",
                    iter, a, b, mj$mean, mj$var, residual, "---", det_Jlog))
        cat(sprintf("\nConverged: ||F|| = %.2e < %.2e\n", residual, tol_F))
      }
      break
    }

    # Newton step
    delta <- NULL
    if (!is.finite(det_Jlog) || abs(det_Jlog) < .TOL_SINGULAR_JACOBIAN) {
      # Singular Jacobian: use gradient descent fallback
      delta <- -0.1 * crossprod(J_log, F_val)
      if (isTRUE(verbose)) {
        cat(sprintf("  Warning: Near-singular Jacobian (det=%.2e), using gradient descent\n",
                    det_Jlog))
      }
    } else {
      delta <- tryCatch(
        -solve(J_log, F_val),
        error = function(e) NULL
      )
      if (is.null(delta)) {
        delta <- -0.1 * F_val
      }
    }

    # Damping via backtracking line search
    step_size <- 1.0
    if (isTRUE(damping)) {
      for (k in seq_len(10L)) {
        eta_new <- eta + step_size * delta
        a_new <- exp(eta_new[1L])
        b_new <- exp(eta_new[2L])

        if (!is.finite(a_new) || !is.finite(b_new) || a_new <= 0 || b_new <= 0) {
          step_size <- step_size * 0.5
          next
        }

        mj_new <- exact_K_moments(J, a_new, b_new, M)
        F_new <- c(mj_new$mean - mu_K, mj_new$var - var_K)
        residual_new <- sqrt(sum(F_new^2))

        if (residual_new < residual) {
          break
        }
        step_size <- step_size * 0.5
      }

      # Guard: if no improvement was found after all halvings, prevent bad step
      old_norm <- residual
      eta_trial <- eta + step_size * delta
      a_trial <- exp(eta_trial[1L])
      b_trial <- exp(eta_trial[2L])
      if (is.finite(a_trial) && is.finite(b_trial) && a_trial > 0 && b_trial > 0) {
        mj_trial <- exact_K_moments(J, a_trial, b_trial, M)
        F_trial <- c(mj_trial$mean - mu_K, mj_trial$var - var_K)
        new_norm <- sqrt(sum(F_trial^2))
      } else {
        new_norm <- old_norm
      }
      if (new_norm >= old_norm) {
        delta <- rep(0, length(delta))
      } else {
        delta <- step_size * delta
      }
    }

    # Store step size in trace
    history$step[nrow(history)] <- step_size

    if (isTRUE(verbose)) {
      cat(sprintf("%4d | %10.6f | %10.6f | %10.6f | %10.6f | %10.2e | %8.4f | %10.2e\n",
                  iter, a, b, mj$mean, mj$var, residual, step_size, det_Jlog))
    }

    # Check for stagnation (step size too small)
    # NOTE: This is NOT convergence - residual may still be above tolerance
    step_norm <- sqrt(sum(delta^2))
    if (step_norm < tol_step) {
      termination <- "step"
      status <- "stagnated"
      if (isTRUE(verbose)) {
        cat(sprintf("\nStagnated: ||delta|| = %.2e < %.2e (residual still %.2e)\n",
                    step_norm, tol_step, residual))
      }
      break
    }

    # Update eta
    eta <- eta + delta
  }

  # -------------------------------------------------------------------------
  # Post-loop Convergence Check
  # -------------------------------------------------------------------------
  # Check if final residual meets tolerance even if loop ended for other reasons
  a_current <- exp(eta[1L])
  b_current <- exp(eta[2L])
  final_check <- exact_K_moments(J, a_current, b_current, M)
  final_residual <- sqrt((final_check$mean - mu_K)^2 + (final_check$var - var_K)^2)

  if (!converged && final_residual < tol_F) {
    converged <- TRUE
    status <- "success"
    termination <- "residual"
  }

  # -------------------------------------------------------------------------
  # Fallback to Nelder-Mead if Newton didn't converge
  # -------------------------------------------------------------------------
  method <- "A2-MN"
  fallback_used <- FALSE

  if (!converged && isTRUE(use_fallback)) {
    if (isTRUE(verbose)) {
      cat("\nNewton did not converge, using Nelder-Mead fallback...\n")
    }

    # Objective function for optimization
    objective <- function(params) {
      a_opt <- params[1L]
      b_opt <- params[2L]
      if (a_opt <= 0 || b_opt <= 0) {
        return(.PENALTY_INF)
      }
      mom <- exact_K_moments(J, a_opt, b_opt, M)
      (mom$mean - mu_K)^2 + (mom$var - var_K)^2
    }

    # Start from Newton's last point
    a_last <- exp(eta[1L])
    b_last <- exp(eta[2L])

    opt_result <- optim(
      par = c(a_last, b_last),
      fn = objective,
      method = "Nelder-Mead",
      control = list(maxit = 1000, reltol = 1e-12)
    )

    if (opt_result$value < tol_F^2) {
      converged <- TRUE
      termination <- "nelder_mead"
      status <- "success"
      method <- "A2-MN+NM"
      fallback_used <- TRUE
      eta <- log(opt_result$par)

      if (isTRUE(verbose)) {
        cat(sprintf("Nelder-Mead converged: objective = %.2e\n", opt_result$value))
      }
    }
  }

  # -------------------------------------------------------------------------
  # Final Result
  # -------------------------------------------------------------------------
  a_star <- exp(eta[1L])
  b_star <- exp(eta[2L])
  final <- exact_K_moments(J, a_star, b_star, M)
  final_residual <- sqrt((final$mean - mu_K)^2 + (final$var - var_K)^2)

  # Build status message with warnings
  if (converged) {
    if (a_star < 0.1) {
      status <- paste(status, "[warning: quasi-improper prior, a < 0.1]")
    }
    if (fallback_used) {
      status <- paste(status, "[Nelder-Mead fallback used]")
    }
  } else {
    if (final_residual < 0.1) {
      status <- paste(status, "[acceptable accuracy]")
    }
  }

  if (isTRUE(verbose) && !converged) {
    cat(sprintf("\nDid not converge after %d iterations (final residual: %.2e)\n",
                nrow(history), final_residual))
  }

  # Construct DPprior_fit object
  structure(
    list(
      a = a_star,
      b = b_star,
      J = J,
      target = list(
        mu_K = mu_K,
        var_K = var_K,
        type = "moments"
      ),
      method = method,
      status = status,
      converged = converged,
      iterations = nrow(history),
      termination = termination,
      fit = list(
        mu_K = final$mean,
        var_K = final$var,
        residual = final_residual
      ),
      diagnostics = list(
        a0 = a0,
        b0 = b0,
        tol_F = tol_F,
        tol_step = tol_step,
        M = M,
        fallback_used = fallback_used
      ),
      trace = history
    ),
    class = "DPprior_fit"
  )
}


# =============================================================================
# Verification Functions
# =============================================================================

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

  pass <- fit$converged &&
    abs(fit$fit$mu_K - mu_K) < tol &&
    abs(fit$fit$var_K - var_K) < tol

  if (isTRUE(verbose)) {
    cat(sprintf("A2-MN Moment Matching Verification\n"))
    cat(strrep("-", 50), "\n")
    cat(sprintf("Target: E[K]=%.4f, Var(K)=%.4f\n", mu_K, var_K))
    cat(sprintf("Achieved: E[K]=%.10f, Var(K)=%.10f\n",
                fit$fit$mu_K, fit$fit$var_K))
    cat(sprintf("Mean error: %.2e\n", abs(fit$fit$mu_K - mu_K)))
    cat(sprintf("Var error: %.2e\n", abs(fit$fit$var_K - var_K)))
    cat(sprintf("Termination: %s\n", fit$termination))
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
  a1_mom <- exact_K_moments(J, a1$a, a1$b)
  a1_residual <- sqrt((a1_mom$mean - mu_K)^2 + (a1_mom$var - var_K)^2)

  # A2 solution
  a2 <- DPprior_a2_newton(J, mu_K, var_K, verbose = FALSE)

  improvement_ratio <- a1_residual / max(a2$fit$residual, 1e-15)

  if (isTRUE(verbose)) {
    cat(sprintf("A1 vs A2 Comparison (J=%d, mu_K=%.2f, var_K=%.2f)\n", J, mu_K, var_K))
    cat(strrep("-", 60), "\n")
    cat(sprintf("%-20s %12s %12s\n", "", "A1", "A2"))
    cat(strrep("-", 60), "\n")
    cat(sprintf("%-20s %12.6f %12.6f\n", "Shape (a)", a1$a, a2$a))
    cat(sprintf("%-20s %12.6f %12.6f\n", "Rate (b)", a1$b, a2$b))
    cat(sprintf("%-20s %12.6f %12.10f\n", "E[K] achieved", a1_mom$mean, a2$fit$mu_K))
    cat(sprintf("%-20s %12.6f %12.10f\n", "Var achieved", a1_mom$var, a2$fit$var_K))
    cat(sprintf("%-20s %12.6f %12.2e\n", "Residual", a1_residual, a2$fit$residual))
    cat(strrep("-", 60), "\n")
    cat(sprintf("Improvement ratio: %.0fx\n", improvement_ratio))
  }

  invisible(list(
    a1 = list(a = a1$a, b = a1$b, mean = a1_mom$mean, var = a1_mom$var,
              residual = a1_residual),
    a2 = list(a = a2$a, b = a2$b, mean = a2$fit$mu_K, var = a2$fit$var_K,
              residual = a2$fit$residual),
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
  test1_pass <- fit$converged &&
    abs(fit$fit$mu_K - 5) < 1e-6 &&
    abs(fit$fit$var_K - 8) < 1e-6

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
  test2_pass <- fit$iterations < 10

  if (isTRUE(verbose)) {
    cat(sprintf("Iterations: %d\n", fit$iterations))
    cat(sprintf("Termination: %s\n", fit$termination))
    cat(sprintf("Result: %s\n\n", if (test2_pass) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test2_pass

  # Test 3: A2 corrects A1 error
  if (isTRUE(verbose)) {
    cat("[Test 3] A2 improves over A1\n")
    cat(strrep("-", 50), "\n")
  }

  comparison <- compare_a1_a2(J = 50, mu_K = 5, var_K = 8, verbose = verbose)
  test3_pass <- comparison$a2$residual < comparison$a1$residual

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
    fit <- DPprior_a2_newton(tc$J, tc$mu_K, tc$var_K, verbose = FALSE)
    case_pass <- fit$converged && fit$fit$residual < 1e-6
    status <- if (case_pass) "PASS" else "FAIL"
    if (isTRUE(verbose)) {
      cat(sprintf("  J=%3d, mu_K=%2d, var_K=%2d: %s (iter=%d, term=%s, res=%.2e)\n",
                  tc$J, tc$mu_K, tc$var_K, status, fit$iterations,
                  fit$termination, fit$fit$residual))
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
  test5_pass <- fit$converged && fit$fit$residual < 1e-5

  if (isTRUE(verbose)) {
    cat(sprintf("Target: mu_K=3, var_K=10 (VIF=%.1f)\n", 10 / 2))
    cat(sprintf("Converged: %s, Iterations: %d\n", fit$converged, fit$iterations))
    cat(sprintf("Termination: %s\n", fit$termination))
    cat(sprintf("Achieved: a=%.6f (quasi-improper: %s)\n", fit$a, fit$a < 0.1))
    cat(sprintf("Residual: %.2e\n", fit$fit$residual))
    cat(sprintf("Result: %s\n\n", if (test5_pass) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test5_pass

  # Test 6: Termination field consistency
  if (isTRUE(verbose)) {
    cat("[Test 6] Termination field consistency\n")
    cat(strrep("-", 50), "\n")
  }

  fit <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8, verbose = FALSE)
  test6_pass <- fit$termination %in% c("residual", "step", "max_iter", "nelder_mead")

  if (isTRUE(verbose)) {
    cat(sprintf("Termination field: '%s'\n", fit$termination))
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
  test7_pass <- all(required_cols %in% names(fit$trace))

  if (isTRUE(verbose)) {
    cat(sprintf("Required columns: %s\n", paste(required_cols, collapse = ", ")))
    cat(sprintf("Present columns:  %s\n", paste(names(fit$trace), collapse = ", ")))
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
