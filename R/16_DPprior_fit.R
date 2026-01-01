# =============================================================================
# Module 16: DPprior_fit() Main Wrapper
# =============================================================================
#
# This module provides the primary user-facing function for DPprior prior
# elicitation. It serves as a one-click entry point that:
#   1. Handles user-friendly parameter specification (confidence levels)
#   2. Validates inputs and handles edge cases
#   3. Dispatches to appropriate calibration methods
#   4. Computes diagnostics and issues warnings
#
# Author: JoonHo Lee
# Date: December 2025
# Part of: DPprior R Package
#
# Key improvements integrated:
#   - var_K upper bound check: var_K <= (J-1)^2/4
#   - A1-only projection (not at wrapper level for A2-MN)
#   - Updated VIF values: low=5.0, medium=2.5, high=1.5
#   - mu_K > 1 constraint (trivial case excluded)
#   - solver_diagnostics preservation
#   - A2-KL dual-mode support (pmf, chisq)
#
# References:
#   - RN-03: Closed-Form Mapping for the Gamma Hyperprior
#   - RN-04: Small-J Correction via Newton Refinement
#   - RN-07: Unintended Prior Diagnostic
# =============================================================================


# =============================================================================
# Confidence Level Utilities
# =============================================================================

#' Convert Confidence Level to Variance Inflation Factor
#'
#' Maps qualitative confidence levels to quantitative VIF values.
#' The VIF determines prior variance as: Var(K) = VIF * (mu_K - 1).
#'
#' @param confidence Character; one of "low", "medium", or "high".
#'
#' @return Numeric; the corresponding VIF value.
#'
#' @details
#' The mapping is:
#' \itemize{
#'   \item "low": VIF = 5.0 (very high uncertainty about cluster count)
#'   \item "medium": VIF = 2.5 (moderate uncertainty)
#'   \item "high": VIF = 1.5 (strong prior belief)
#' }
#'
#' These values are chosen to provide a wide range of uncertainty levels.
#' VIF = 5 corresponds to "I have little idea how many clusters there are."
#'
#' @seealso \code{\link{vif_to_variance_fit}}, \code{\link{DPprior_fit}}
#'
#' @examples
#' \dontrun{
#' confidence_to_vif_fit("low")     # 5.0
#' confidence_to_vif_fit("medium")  # 2.5
#' confidence_to_vif_fit("high")    # 1.5
#' }
#'
#' @keywords internal
confidence_to_vif_fit <- function(confidence) {
  valid <- c("low", "medium", "high")
  if (!is.character(confidence) || length(confidence) != 1L ||
      !(confidence %in% valid)) {
    stop(sprintf("confidence must be one of: %s",
                 paste(valid, collapse = ", ")), call. = FALSE)
  }

  switch(confidence,
         "low" = 5.0,
         "medium" = 2.5,
         "high" = 1.5)
}


#' Convert VIF to Target Variance
#'
#' Computes the target variance from a Variance Inflation Factor.
#'
#' @param mu_K Numeric; target expected number of clusters.
#' @param vif Numeric; Variance Inflation Factor (must be > 1).
#'
#' @return Numeric; target variance: \code{vif * (mu_K - 1)}.
#'
#' @details
#' The VIF relates variance to the shifted mean (mu_K - 1):
#' \deqn{Var(K_J) = VIF \cdot (\mu_K - 1)}
#'
#' This parameterization is natural for the shifted NegBin approximation.
#'
#' @seealso \code{\link{confidence_to_vif_fit}}, \code{\link{DPprior_fit}}
#'
#' @examples
#' \dontrun{
#' vif_to_variance_fit(5, 2.5)  # 10 = 2.5 * (5 - 1)
#' }
#'
#' @keywords internal
vif_to_variance_fit <- function(mu_K, vif) {
  if (vif <= 1) {
    warning("VIF should be > 1 for overdispersion; using VIF = 1.01",
            call. = FALSE)
    vif <- 1.01
  }
  vif * (mu_K - 1)
}


# =============================================================================
# Null-Coalescing Operator
# =============================================================================

#' Null-Coalescing Operator
#'
#' Returns the left operand if not NULL, otherwise the right operand.
#'
#' @param x Left operand.
#' @param y Right operand (default value).
#'
#' @return \code{x} if not NULL, else \code{y}.
#'
#' @noRd
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


# =============================================================================
# Main Entry Point: DPprior_fit()
# =============================================================================

#' Fit a Gamma Hyperprior for DP Concentration Parameter
#'
#' Main entry point for DPprior prior elicitation. Takes user-specified
#' cluster count expectations and returns calibrated Gamma(a, b) parameters.
#'
#' This function provides a convenient interface for specifying beliefs about
#' the number of clusters and uncertainty, then calibrates a Gamma hyperprior
#' for the Dirichlet Process concentration parameter alpha.
#'
#' @param J Integer; sample size (number of sites/units). Must be positive
#'   and not exceed the maximum supported value (default: 500).
#' @param mu_K Numeric; target expected number of clusters \code{E[K_J]}. Must be
#'   in the range (1, J]. Note: mu_K = 1 is trivial (single cluster).
#' @param var_K Numeric; target variance of \code{K_J}. If NULL, computed from
#'   the \code{confidence} argument. Must satisfy var_K <= (J-1)^2/4
#'   (maximum possible variance for K in \{1,...,J\}).
#' @param confidence Character; alternative to var_K for specifying uncertainty.
#'   One of:
#'   \itemize{
#'     \item "low": VIF = 5.0 (very high uncertainty, wide prior)
#'     \item "medium": VIF = 2.5 (moderate uncertainty)
#'     \item "high": VIF = 1.5 (low uncertainty, concentrated prior)
#'   }
#'   Only used if \code{var_K} is NULL.
#' @param method Character; calibration method. One of:
#'   \itemize{
#'     \item "A2-MN" (default): Exact moment matching via Newton's method.
#'           Recommended for accurate calibration.
#'     \item "A1": Fast closed-form approximation based on shifted NegBin.
#'           Good for initial exploration or large J. Note: may project
#'           var_K to feasible region if var_K < mu_K - 1.
#'     \item "A2-KL": KL divergence minimization. Supports both moment
#'           target (default) and custom target PMF via \code{target_pmf}.
#'   }
#' @param target_pmf Numeric vector; optional custom target PMF for A2-KL.
#'   If provided, A2-KL uses PMF matching mode. Length must equal J.
#' @param check_diagnostics Logical; if TRUE (default), compute comprehensive
#'   diagnostics including weight distribution analysis.
#' @param warn_dominance Logical; if TRUE (default), issue a warning if the
#'   prior implies high probability of a dominant cluster (P(w1 > 0.5) > 40%).
#' @param M Integer; number of quadrature nodes for numerical integration.
#'   Default is 80, which provides good accuracy for most cases.
#' @param verbose Logical; if TRUE, print progress messages during calibration.
#' @param ... Additional arguments passed to method-specific functions.
#'   Note: only matching formal arguments are forwarded to backends.
#'
#' @return An S3 object of class "DPprior_fit" containing:
#'   \describe{
#'     \item{a}{Numeric; shape parameter of the calibrated Gamma hyperprior.}
#'     \item{b}{Numeric; rate parameter of the calibrated Gamma hyperprior.}
#'     \item{J}{Integer; sample size used for calibration.}
#'     \item{target}{List; target specification including mu_K, var_K,
#'       confidence level (if used), and target type.}
#'     \item{method}{Character; calibration method used.}
#'     \item{converged}{Logical; whether the calibration converged.}
#'     \item{iterations}{Integer; number of iterations (for iterative methods).}
#'     \item{fit}{List; achieved moments (mu_K, var_K) and residual.}
#'     \item{diagnostics}{List; comprehensive RN-07 diagnostics (if requested).}
#'     \item{solver_diagnostics}{List; backend solver details (if applicable).}
#'     \item{trace}{Data frame; iteration history (for iterative methods).}
#'   }
#'
#' @details
#' \subsection{Variance Constraints}{
#' The target variance must satisfy two constraints:
#' \enumerate{
#'   \item Upper bound: \code{var_K <= (J-1)^2/4} (maximum possible variance
#'         for a distribution on \{1,...,J\})
#'   \item Lower bound (A1 only): \code{var_K >= mu_K - 1} (NegBin feasibility).
#'         For A1, infeasible variance is projected to the boundary with a warning.
#'         A2-MN does not have this constraint.
#' }
#' }
#'
#' \subsection{Diagnostics}{
#' When \code{check_diagnostics = TRUE}, the function computes:
#' \itemize{
#'   \item Alpha distribution: \code{E[alpha]}, \code{Var(alpha)}, \code{CV(alpha)}
#'   \item K distribution: achieved moments, mode, quantiles
#'   \item Weight distribution: \code{E[w_1]}, \code{P(w_1 > 0.5)}, dominance risk
#'   \item Co-clustering: \code{E[rho]}, interpretation
#' }
#' Backend solver details are preserved in \code{solver_diagnostics}.
#' }
#'
#' @seealso
#' \code{\link{DPprior_a1}} for A1 closed-form approximation,
#' \code{\link{DPprior_a2_newton}} for A2 Newton refinement,
#' \code{\link{DPprior_a2_kl}} for A2-KL divergence minimization,
#' \code{\link{DPprior_diagnostics}} for detailed diagnostics,
#' \code{\link{DPprior_dual}} for dual-anchor calibration with weight constraints.
#'
#' @examples
#' # Basic usage with target moments
#' fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8)
#' print(fit)
#'
#' # Using confidence level instead of explicit variance
#' fit_medium <- DPprior_fit(J = 50, mu_K = 5, confidence = "medium")
#' fit_low <- DPprior_fit(J = 50, mu_K = 5, confidence = "low")
#'
#' # Quick approximation for exploration
#' fit_quick <- DPprior_fit(J = 50, mu_K = 5, var_K = 8, method = "A1")
#'
#' # Verbose output for debugging
#' fit_verbose <- DPprior_fit(J = 50, mu_K = 5, var_K = 8, verbose = TRUE)
#'
#' # Access results
#' fit$a  # Gamma shape
#' fit$b  # Gamma rate
#' fit$fit$mu_K  # Achieved mean
#' fit$diagnostics$weights$dominance_risk  # Dominance risk level
#'
#' @references
#' Lee, J. (2025). Prior Elicitation for the Dirichlet Process Concentration
#' Parameter in Low-Information Settings. \emph{Working Paper}.
#'
#' RN-03: Closed-Form Mapping for the Gamma Hyperprior.
#' RN-04: Small-J Correction via Newton Refinement.
#' RN-07: Unintended Prior Diagnostic.
#'
#' @export
DPprior_fit <- function(J, mu_K, var_K = NULL,
                        confidence = c("medium", "low", "high"),
                        method = c("A2-MN", "A1", "A2-KL"),
                        target_pmf = NULL,
                        check_diagnostics = TRUE,
                        warn_dominance = TRUE,
                        M = .QUAD_NODES_DEFAULT,
                        verbose = FALSE, ...) {

  # ---------------------------------------------------------------------------
  # Input Validation
  # ---------------------------------------------------------------------------

  # Validate J
  assert_valid_J(J)

  # Validate M
  if (!is.numeric(M) || length(M) != 1L || M < 10L || M != floor(M)) {
    stop("M must be an integer >= 10", call. = FALSE)
  }
  M <- as.integer(M)

  # Validate mu_K: must be finite, positive, and in (1, J]
  if (!is.numeric(mu_K) || length(mu_K) != 1L || !is.finite(mu_K)) {
    stop("mu_K must be a finite numeric scalar", call. = FALSE)
  }

  # GPT improvement: mu_K = 1 is trivial (single cluster)
  if (mu_K <= 1) {
    stop("mu_K must be > 1 (mu_K = 1 implies trivial single-cluster case)",
         call. = FALSE)
  }

  if (mu_K > J) {
    stop(sprintf("mu_K must be <= J (got mu_K = %.2f, J = %d)", mu_K, J),
         call. = FALSE)
  }

  # Match method argument
  method <- match.arg(method)

  # ---------------------------------------------------------------------------
  # Handle var_K vs Confidence
  # ---------------------------------------------------------------------------

  confidence_used <- NULL

  if (is.null(var_K)) {
    # Use confidence level
    confidence <- match.arg(confidence)
    vif <- confidence_to_vif_fit(confidence)
    var_K <- vif_to_variance_fit(mu_K, vif)
    confidence_used <- confidence

    if (isTRUE(verbose)) {
      message(sprintf("Using confidence='%s' -> VIF=%.1f -> var_K=%.2f",
                      confidence, vif, var_K))
    }
  } else {
    # Explicit var_K provided
    if (!is.numeric(var_K) || length(var_K) != 1L || !is.finite(var_K)) {
      stop("var_K must be a finite numeric scalar", call. = FALSE)
    }
    if (var_K <= 0) {
      stop("var_K must be positive", call. = FALSE)
    }
  }

  # ---------------------------------------------------------------------------
  # GPT improvement: var_K Upper Bound Check (Universal)
  # ---------------------------------------------------------------------------

  # Hard upper bound: K_J in {1, ..., J} implies Var(K) <= (J-1)^2 / 4
  # This is the variance of a distribution at the endpoints (Bernoulli-like)
  var_upper <- (J - 1)^2 / 4

  if (var_K > var_upper + 1e-12) {
    stop(sprintf(
      "var_K = %.2f exceeds maximum possible variance (%.2f) for K in {1,...,%d}.\n  Reduce var_K or use a lower confidence level.",
      var_K, var_upper, J), call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Validate target_pmf (if provided)
  # ---------------------------------------------------------------------------

  if (!is.null(target_pmf)) {
    if (!is.numeric(target_pmf) || length(target_pmf) != J) {
      stop(sprintf("target_pmf must be a numeric vector of length J = %d", J),
           call. = FALSE)
    }
    if (any(target_pmf < 0) || abs(sum(target_pmf) - 1) > 1e-6) {
      stop("target_pmf must be non-negative and sum to 1", call. = FALSE)
    }
  }

  # ---------------------------------------------------------------------------
  # Method Dispatch
  # ---------------------------------------------------------------------------

  # Store original var_K for reporting
  var_K_original <- var_K
  var_K_used <- var_K

  fit <- switch(method,
                "A1" = {
                  if (isTRUE(verbose)) {
                    message("Using A1 closed-form approximation")
                  }

                  # GPT improvement: A1-only feasibility projection
                  # NegBin requires var_K >= mu_K - 1
                  # A2-MN does NOT have this constraint
                  min_var <- mu_K - 1 + .TOL_PROJECTION_BUFFER
                  if (var_K < min_var) {
                    warning(sprintf(
                      "A1 method: var_K=%.4f < mu_K-1=%.4f (infeasible for NegBin proxy).\n  Projecting to feasible boundary: %.6f",
                      var_K, mu_K - 1, min_var), call. = FALSE)
                    var_K_used <- min_var
                  }

                  DPprior_a1(J, mu_K, var_K_used)
                },

                "A2-MN" = {
                  if (isTRUE(verbose)) {
                    message("Using A2-MN Newton refinement")
                  }

                  # GPT insight: A2-MN does NOT require the NegBin feasibility constraint
                  # It matches exact DP moments, which can handle any valid (mu_K, var_K)

                  DPprior_a2_newton(J, mu_K, var_K_used, M = M, verbose = verbose)
                },

                "A2-KL" = {
                  if (isTRUE(verbose)) {
                    message("Using A2-KL divergence minimization")
                  }

                  # GPT improvement: dual-mode dispatch for A2-KL
                  # API: DPprior_a2_kl(J, target, method = c("pmf", "chisq"), ...)
                  if (!is.null(target_pmf)) {
                    # PMF matching mode: target is numeric vector
                    DPprior_a2_kl(J, target = target_pmf, method = "pmf",
                                  M = M, verbose = verbose)
                  } else {
                    # Moment-based chi-squared mode: target is list(mu_K, var_K)
                    DPprior_a2_kl(J, target = list(mu_K = mu_K, var_K = var_K_used),
                                  method = "chisq", M = M, verbose = verbose)
                  }
                }
  )

  # ---------------------------------------------------------------------------
  # Build Consistent Output Structure
  # ---------------------------------------------------------------------------

  # Compute achieved moments if not provided by backend
  achieved <- fit$fit %||% {
    mom <- exact_K_moments(J, fit$a, fit$b, M)
    list(
      mu_K = mom$mean,
      var_K = mom$var,
      residual = sqrt((mom$mean - mu_K)^2 + (mom$var - var_K_used)^2)
    )
  }

  # GPT improvement: preserve solver-specific diagnostics before overwriting
  solver_diag <- NULL
  if (!is.null(fit$diagnostics) &&
      !("weights" %in% names(fit$diagnostics))) {
    solver_diag <- fit$diagnostics
  }

  result <- list(
    a = fit$a,
    b = fit$b,
    J = J,
    target = list(
      mu_K = mu_K,
      var_K = var_K_original,   # Store original (before any projection)
      var_K_used = var_K_used,  # Actual value used (may be projected for A1)
      confidence = confidence_used,
      type = if (!is.null(target_pmf)) "pmf" else "moments"
    ),
    method = fit$method %||% method,
    status = fit$status %||% "success",
    converged = fit$converged %||% TRUE,
    iterations = fit$iterations %||% NA_integer_,
    termination = fit$termination %||% NA_character_,
    fit = achieved,
    solver_diagnostics = solver_diag,
    trace = fit$trace %||% NULL
  )

  # ---------------------------------------------------------------------------
  # Optional Diagnostics (RN-07)
  # ---------------------------------------------------------------------------

  if (isTRUE(check_diagnostics)) {
    # Check if DPprior_diagnostics is available
    if (exists("DPprior_diagnostics", mode = "function")) {
      result$diagnostics <- DPprior_diagnostics(result)

      # Check for high dominance risk and warn
      if (isTRUE(warn_dominance)) {
        if (!is.null(result$diagnostics$weights$dominance_risk) &&
            result$diagnostics$weights$dominance_risk == "high") {

          # Get the actual probability for the warning message
          p_gt_50 <- prob_w1_exceeds(0.5, result$a, result$b)

          warning(sprintf(
            "HIGH DOMINANCE RISK: P(w1 > 0.5) = %.1f%% exceeds 40%%.\n",
            100 * p_gt_50),
            "  This may indicate unintended prior behavior (RN-07).\n",
            "  Consider using DPprior_dual() for weight-constrained elicitation.\n",
            "  See ?DPprior_diagnostics for interpretation.",
            call. = FALSE)

          result$dominance_warning <- TRUE
        }
      }
    } else {
      warning("DPprior_diagnostics() not available. Skipping diagnostics.",
              call. = FALSE)
      result$diagnostics <- NULL
    }
  }

  # ---------------------------------------------------------------------------
  # Set Class and Return
  # ---------------------------------------------------------------------------

  class(result) <- "DPprior_fit"
  result
}


# =============================================================================
# S3 Methods for DPprior_fit
# =============================================================================

#' @export
print.DPprior_fit <- function(x, ...) {
  cat("DPprior Elicitation Result\n")
  cat(strrep("=", 55), "\n")

  cat(sprintf("Method: %s\n", x$method))
  cat(sprintf("Sample size J: %d\n", x$J))

  cat("\nTarget:\n")
  cat(sprintf("  E[K_J]   = %.4f\n", x$target$mu_K))
  if (!is.null(x$target$var_K_used) &&
      abs(x$target$var_K - x$target$var_K_used) > 1e-10) {
    cat(sprintf("  Var(K_J) = %.4f (requested), %.4f (used after projection)\n",
                x$target$var_K, x$target$var_K_used))
  } else {
    cat(sprintf("  Var(K_J) = %.4f\n", x$target$var_K))
  }
  if (!is.null(x$target$confidence)) {
    cat(sprintf("  (from confidence = '%s')\n", x$target$confidence))
  }

  cat("\nOptimal Gamma(a, b) hyperprior:\n")
  cat(sprintf("  a (shape) = %.8f\n", x$a))
  cat(sprintf("  b (rate)  = %.8f\n", x$b))
  cat(sprintf("  E[alpha]  = %.4f\n", x$a / x$b))
  cat(sprintf("  CV(alpha) = %.4f\n", 1 / sqrt(x$a)))

  cat("\nConvergence:\n")
  if (!is.null(x$status)) {
    cat(sprintf("  Status:      %s\n", x$status))
  }
  cat(sprintf("  Converged:   %s\n", x$converged))
  if (!is.na(x$iterations)) {
    cat(sprintf("  Iterations:  %d\n", x$iterations))
  }
  if (!is.null(x$termination) && !is.na(x$termination)) {
    cat(sprintf("  Termination: %s\n", x$termination))
  }

  if (!is.null(x$fit)) {
    cat("\nAchieved fit:\n")
    cat(sprintf("  E[K_J]   = %.10f\n", x$fit$mu_K))
    cat(sprintf("  Var(K_J) = %.10f\n", x$fit$var_K))
    if (!is.null(x$fit$residual) && !is.na(x$fit$residual)) {
      cat(sprintf("  Residual = %.2e\n", x$fit$residual))
    }
  }

  # Print diagnostic summary if available
  if (!is.null(x$diagnostics)) {
    cat("\nDiagnostics Summary:\n")
    if (!is.null(x$diagnostics$weights)) {
      risk <- x$diagnostics$weights$dominance_risk
      risk_symbol <- switch(risk,
                            "low" = "\u2713",      # checkmark
                            "moderate" = "\u26A0", # warning
                            "high" = "\u2718",     # X mark
                            "?")
      cat(sprintf("  Dominance risk: %s (%s)\n", toupper(risk), risk_symbol))
    }
    cat("  (Run DPprior_diagnostics(fit) for full report)\n")
  }

  invisible(x)
}


#' @export
summary.DPprior_fit <- function(object, ...) {
  # Helper to safely extract a scalar value
  safe_scalar <- function(x, default = NA) {
    if (is.null(x) || length(x) == 0L) return(default)
    if (length(x) == 1L) return(x)
    return(x[1L])
  }

  # Return list structure matching test-10_a1_mapping.R expectations
  # Expected: c("method", "status", "gamma_prior", "alpha_summary",
  #             "target", "scaling", "converged", "iterations")
  result <- list(
    method = safe_scalar(object$method, NA_character_),
    status = safe_scalar(object$status, "success"),
    gamma_prior = list(
      a = safe_scalar(object$a, NA_real_),
      b = safe_scalar(object$b, NA_real_)
    ),
    alpha_summary = list(
      E_alpha = safe_scalar(object$a / object$b, NA_real_),
      Var_alpha = safe_scalar(object$a / object$b^2, NA_real_),
      CV_alpha = safe_scalar(1 / sqrt(object$a), NA_real_)
    ),
    target = list(
      mu_K = safe_scalar(object$target$mu_K, NA_real_),
      var_K = safe_scalar(object$target$var_K, NA_real_),
      var_K_used = safe_scalar(object$target$var_K_used %||% object$target$var_K, NA_real_),
      confidence = safe_scalar(object$target$confidence, NA_character_)
    ),
    scaling = list(
      J = safe_scalar(object$J, NA_integer_)
    ),
    converged = safe_scalar(object$converged, NA),
    iterations = safe_scalar(object$iterations, NA_integer_)
  )

  class(result) <- "summary.DPprior_fit"
  result
}


# =============================================================================
# Verification Function
# =============================================================================

#' Verify DPprior_fit Module
#'
#' Runs comprehensive verification tests for the DPprior_fit module.
#'
#' @param verbose Logical; if TRUE, print detailed test output.
#'
#' @return Invisibly returns TRUE if all tests pass.
#'
#' @examples
#' verify_DPprior_fit()
#'
#' @export
verify_DPprior_fit <- function(verbose = TRUE) {

  if (isTRUE(verbose)) {
    cat(strrep("=", 60), "\n")
    cat("Module 16: DPprior_fit() - Verification Suite\n")
    cat(strrep("=", 60), "\n\n")
  }

  all_pass <- TRUE

  # -------------------------------------------------------------------------
  # Test 1: Confidence to Variance Conversion (Updated VIF values)
  # -------------------------------------------------------------------------
  if (isTRUE(verbose)) {
    cat("[Test 1] Confidence to Variance Conversion (GPT VIF values)\n")
    cat(strrep("-", 50), "\n")
  }

  mu_K <- 5.0
  expected_vifs <- c("low" = 5.0, "medium" = 2.5, "high" = 1.5)

  for (conf in names(expected_vifs)) {
    vif <- confidence_to_vif_fit(conf)
    var_K <- vif_to_variance_fit(mu_K, vif)
    expected <- expected_vifs[conf] * (mu_K - 1)

    test_pass <- abs(vif - expected_vifs[conf]) < 1e-10 &&
      abs(var_K - expected) < 1e-10

    if (isTRUE(verbose)) {
      cat(sprintf("  confidence='%s': VIF=%.1f, var_K=%.2f [%s]\n",
                  conf, vif, var_K, if (test_pass) "PASS" else "FAIL"))
    }
    all_pass <- all_pass && test_pass
  }
  if (isTRUE(verbose)) cat("\n")

  # -------------------------------------------------------------------------
  # Test 2: var_K Upper Bound Check (GPT improvement)
  # -------------------------------------------------------------------------
  if (isTRUE(verbose)) {
    cat("[Test 2] var_K Upper Bound Check (GPT improvement)\n")
    cat(strrep("-", 50), "\n")
  }

  J <- 50
  var_upper <- (J - 1)^2 / 4  # 600.25

  test_error <- tryCatch({
    DPprior_fit(J = J, mu_K = 25, var_K = 700,  # Exceeds (49)^2/4
                check_diagnostics = FALSE)
    FALSE  # Should not reach here
  }, error = function(e) {
    grepl("exceeds maximum", e$message)
  })

  if (isTRUE(verbose)) {
    cat(sprintf("  Max variance for J=%d: %.2f\n", J, var_upper))
    cat(sprintf("  Error on var_K=700: %s [%s]\n",
                test_error, if (test_error) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test_error
  if (isTRUE(verbose)) cat("\n")

  # -------------------------------------------------------------------------
  # Test 3: A1-only Feasibility Projection (GPT insight)
  # -------------------------------------------------------------------------
  if (isTRUE(verbose)) {
    cat("[Test 3] A1-only Feasibility Projection (GPT insight)\n")
    cat(strrep("-", 50), "\n")
  }

  # A1 should project var_K < mu_K - 1
  fit_a1 <- suppressWarnings(
    DPprior_fit(J = 50, mu_K = 5, var_K = 3,  # 3 < 4 = mu_K - 1
                method = "A1", check_diagnostics = FALSE)
  )

  a1_projected <- fit_a1$target$var_K_used >= (5 - 1)

  # A2-MN with valid var_K (>= mu_K - 1) should not need projection
  # Use var_K = 5 which satisfies var_K >= mu_K - 1 = 4
  fit_a2 <- suppressWarnings(
    DPprior_fit(J = 50, mu_K = 5, var_K = 5,
                method = "A2-MN", check_diagnostics = FALSE,
                warn_dominance = FALSE)
  )

  # A2 uses the original var_K when feasible
  a2_no_project <- abs(fit_a2$target$var_K_used - 5) < 1e-6

  if (isTRUE(verbose)) {
    cat(sprintf("  A1 with var_K=3: projected to %.6f [%s]\n",
                fit_a1$target$var_K_used,
                if (a1_projected) "PASS" else "FAIL"))
    cat(sprintf("  A2-MN with var_K=5: var_K_used=%.6f [%s]\n",
                fit_a2$target$var_K_used,
                if (a2_no_project) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && a1_projected && a2_no_project
  if (isTRUE(verbose)) cat("\n")

  # -------------------------------------------------------------------------
  # Test 4: mu_K Boundary Check (GPT improvement)
  # -------------------------------------------------------------------------
  if (isTRUE(verbose)) {
    cat("[Test 4] mu_K Boundary Check (GPT improvement)\n")
    cat(strrep("-", 50), "\n")
  }

  # mu_K = 1 should error
  test_mu1_error <- tryCatch({
    DPprior_fit(J = 50, mu_K = 1, var_K = 8, check_diagnostics = FALSE)
    FALSE
  }, error = function(e) {
    grepl("mu_K must be > 1", e$message)
  })

  # mu_K > J should error
  test_mu_gt_J <- tryCatch({
    DPprior_fit(J = 50, mu_K = 51, var_K = 8, check_diagnostics = FALSE)
    FALSE
  }, error = function(e) {
    grepl("mu_K must be <= J", e$message)
  })

  if (isTRUE(verbose)) {
    cat(sprintf("  mu_K=1 rejected: %s [%s]\n",
                test_mu1_error, if (test_mu1_error) "PASS" else "FAIL"))
    cat(sprintf("  mu_K>J rejected: %s [%s]\n",
                test_mu_gt_J, if (test_mu_gt_J) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test_mu1_error && test_mu_gt_J
  if (isTRUE(verbose)) cat("\n")

  # -------------------------------------------------------------------------
  # Test 5: Method Dispatch
  # -------------------------------------------------------------------------
  if (isTRUE(verbose)) {
    cat("[Test 5] Method Dispatch\n")
    cat(strrep("-", 50), "\n")
  }

  for (meth in c("A1", "A2-MN")) {
    fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8, method = meth,
                       check_diagnostics = FALSE, warn_dominance = FALSE)

    test_pass <- fit$a > 0 && fit$b > 0 &&
      grepl(substr(meth, 1, 2), fit$method)

    if (isTRUE(verbose)) {
      cat(sprintf("  method='%s': a=%.6f, b=%.6f [%s]\n",
                  meth, fit$a, fit$b, if (test_pass) "PASS" else "FAIL"))
    }
    all_pass <- all_pass && test_pass
  }
  if (isTRUE(verbose)) cat("\n")

  # -------------------------------------------------------------------------
  # Test 6: A2 Improves A1
  # -------------------------------------------------------------------------
  if (isTRUE(verbose)) {
    cat("[Test 6] A2 Improves A1 Accuracy\n")
    cat(strrep("-", 50), "\n")
  }

  mu_K <- 5; var_K <- 8

  fit_a1 <- DPprior_fit(J = 50, mu_K = mu_K, var_K = var_K, method = "A1",
                        check_diagnostics = FALSE)
  fit_a2 <- DPprior_fit(J = 50, mu_K = mu_K, var_K = var_K, method = "A2-MN",
                        check_diagnostics = FALSE)

  # Compute A1 residual
  a1_mom <- exact_K_moments(50, fit_a1$a, fit_a1$b)
  a1_residual <- sqrt((a1_mom$mean - mu_K)^2 + (a1_mom$var - var_K)^2)

  # A2 residual from fit
  a2_residual <- fit_a2$fit$residual %||% {
    a2_mom <- exact_K_moments(50, fit_a2$a, fit_a2$b)
    sqrt((a2_mom$mean - mu_K)^2 + (a2_mom$var - var_K)^2)
  }

  test_pass <- a2_residual < a1_residual

  if (isTRUE(verbose)) {
    cat(sprintf("  A1 residual: %.6e\n", a1_residual))
    cat(sprintf("  A2 residual: %.6e\n", a2_residual))
    cat(sprintf("  Improvement: %.0fx [%s]\n",
                a1_residual / max(a2_residual, 1e-15),
                if (test_pass) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test_pass
  if (isTRUE(verbose)) cat("\n")

  # -------------------------------------------------------------------------
  # Test 7: Diagnostics and solver_diagnostics (GPT improvement)
  # -------------------------------------------------------------------------
  if (isTRUE(verbose)) {
    cat("[Test 7] Diagnostics Structure (GPT: solver_diagnostics)\n")
    cat(strrep("-", 50), "\n")
  }

  fit <- suppressWarnings(
    DPprior_fit(J = 50, mu_K = 5, var_K = 8, method = "A2-MN",
                check_diagnostics = TRUE, warn_dominance = FALSE)
  )

  has_diagnostics <- !is.null(fit$diagnostics)
  has_weights <- !is.null(fit$diagnostics$weights)

  if (isTRUE(verbose)) {
    cat(sprintf("  diagnostics present: %s\n", has_diagnostics))
    cat(sprintf("  weights diagnostics: %s\n", has_weights))
    if (has_weights) {
      cat(sprintf("  dominance_risk: %s\n",
                  fit$diagnostics$weights$dominance_risk))
    }
    cat(sprintf("  [%s]\n", if (has_diagnostics && has_weights) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && has_diagnostics && has_weights
  if (isTRUE(verbose)) cat("\n")

  # -------------------------------------------------------------------------
  # Test 8: Output Structure
  # -------------------------------------------------------------------------
  if (isTRUE(verbose)) {
    cat("[Test 8] Output Structure\n")
    cat(strrep("-", 50), "\n")
  }

  required_fields <- c("a", "b", "J", "target", "method",
                       "converged", "fit")

  fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8,
                     check_diagnostics = FALSE, warn_dominance = FALSE)

  test8_pass <- TRUE
  for (field in required_fields) {
    present <- field %in% names(fit)
    if (isTRUE(verbose)) {
      cat(sprintf("  $%s: %s\n", field, if (present) "present" else "MISSING"))
    }
    test8_pass <- test8_pass && present
  }

  # Check nested structure (including GPT's var_K_used)
  target_ok <- all(c("mu_K", "var_K", "var_K_used") %in% names(fit$target))
  fit_ok <- all(c("mu_K", "var_K") %in% names(fit$fit))

  if (isTRUE(verbose)) {
    cat(sprintf("  $target structure (incl. var_K_used): %s\n",
                if (target_ok) "OK" else "FAIL"))
    cat(sprintf("  $fit structure: %s\n", if (fit_ok) "OK" else "FAIL"))
  }
  all_pass <- all_pass && test8_pass && target_ok && fit_ok
  if (isTRUE(verbose)) cat("\n")

  # -------------------------------------------------------------------------
  # Summary
  # -------------------------------------------------------------------------
  if (isTRUE(verbose)) {
    cat(strrep("=", 60), "\n")
    cat(sprintf("Overall Result: %s\n",
                if (all_pass) "ALL TESTS PASSED" else "SOME TESTS FAILED"))
    cat(strrep("=", 60), "\n")
  }

  invisible(all_pass)
}
