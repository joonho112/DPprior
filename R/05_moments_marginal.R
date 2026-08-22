# =============================================================================
# Module 05: Marginal Moments of K_J under alpha ~ Gamma(a, b)
# =============================================================================
#
# This module provides computation of marginal moments E[K_J] and Var(K_J)
# when the concentration parameter alpha has a Gamma(a, b) prior distribution.
#
# Theory Background (Lee, 2026, Section 3.2):
# ------------------------------------
# Using the Law of Total Expectation and Variance:
#
#   M_1(a,b) = E[K_J | a,b] = E_{alpha ~ Gamma(a,b)}[mu_J(alpha)]
#
#   V(a,b) = Var(K_J | a,b) = E[v_J(alpha)] + Var(mu_J(alpha))
#                           = E[v_J(alpha)] + E[mu_J(alpha)^2] - M_1^2
#
# where:
#   mu_J(alpha) = E[K_J | alpha] = alpha * (psi(alpha+J) - psi(alpha))
#   v_J(alpha)  = Var(K_J | alpha) = mu_J(alpha) - alpha^2 * (psi1(alpha) - psi1(alpha+J))
#
# The expectations are computed via Gauss-Laguerre quadrature (Module 02).
#
# Key Properties:
# ---------------
# 1. Within-alpha and between-alpha variance components are non-negative
# 2. Mean bounds: 1 <= E[K_J] <= J
# 3. Total variance: Var(K_J) = E[Var(K_J | alpha)] +
#    Var(E[K_J | alpha])
#
# Author: JoonHo Lee (jlee296@ua.edu)
# Date: December 2025
# Part of: DPprior R Package
# Reference: Lee (2026), Sections 2--3 and Section 3.2
# Dependencies: Module 02 (quadrature), Module 03 (conditional moments)
# =============================================================================


# =============================================================================
# Core Exported Functions
# =============================================================================

# Validate fixed-order and optional higher-order verification controls shared by
# the marginal moment and PMF modules. A missing M_verify is permitted only in
# non-strict mode and is recorded as an unverified approximation by callers.
.marginal_verification_controls <- function(
    M, M_verify = NULL, abs_tol = 1e-10, rel_tol = 1e-8,
    strict = FALSE) {
  M <- .dpprior_validate_count(
    M, "M", minimum = 1L, maximum = .QUADRATURE_MAX_NODES,
    .subclass = "dpprior_control_error"
  )
  abs_tol <- .dpprior_validate_scalar(
    abs_tol, "abs_tol", lower = 0,
    .subclass = "dpprior_control_error"
  )
  rel_tol <- .dpprior_validate_scalar(
    rel_tol, "rel_tol", lower = 0,
    .subclass = "dpprior_control_error"
  )
  strict <- .dpprior_validate_control(strict, "strict", type = "logical")
  M_verification_required <-
    .quadrature_verification_required_order(M)
  verification_available <-
    M_verification_required <= .QUADRATURE_MAX_NODES

  if (!verification_available && (!is.null(M_verify) || strict)) {
    .dpprior_abort_invalid(
      sprintf(
        paste(
          "marginal verification is unavailable for M=%d:",
          "required M_verify=%d exceeds the supported ceiling (%d)"
        ),
        M, M_verification_required, .QUADRATURE_MAX_NODES
      ),
      c(
        "dpprior_marginal_verification_error",
        "dpprior_marginal_convergence_error",
        "dpprior_bounds_error"
      ),
      "M_verify", M_verification_required,
      sprintf("required order <= %d", .QUADRATURE_MAX_NODES),
      "verification_unavailable"
    )
  }

  if (verification_available && !is.null(M_verify)) {
    M_verify <- .dpprior_validate_count(
      M_verify, "M_verify", minimum = 1L,
      maximum = .QUADRATURE_MAX_NODES,
      .subclass = "dpprior_marginal_verification_error"
    )
    if (M_verify < M_verification_required) {
      .dpprior_abort_invalid(
        sprintf(
          "M_verify must be at least %d for selected order M=%d",
          M_verification_required, M
        ),
        c("dpprior_marginal_verification_error", "dpprior_bounds_error"),
        "M_verify", M_verify,
        sprintf(
          "integer in [%d, %d]",
          M_verification_required, .QUADRATURE_MAX_NODES
        ),
        "insufficient_verification_order"
      )
    }
  } else if (verification_available && strict) {
    .dpprior_abort_invalid(
      "strict marginal computation requires an explicit M_verify",
      "dpprior_marginal_convergence_error", "M_verify", NULL,
      sprintf(
        "integer in [%d, %d]",
        M_verification_required, .QUADRATURE_MAX_NODES
      ),
      "verification_required"
    )
  }

  list(
    M = M,
    M_verify = M_verify,
    abs_tol = abs_tol,
    rel_tol = rel_tol,
    strict = strict,
    M_verification_required = M_verification_required,
    verification_available = verification_available
  )
}


.marginal_moments_fixed <- function(J, a, b, M) {
  quad <- build_gamma_quadrature(a, b, M)
  alphas <- quad$alpha_nodes
  weights <- quad$weights_normalized
  mu <- mean_K_given_alpha(J, alphas)
  variance <- var_K_given_alpha(J, alphas)

  mean_K <- .quadrature_weighted_sum(weights, mu)
  within_alpha <- .quadrature_weighted_sum(weights, variance)
  between_alpha <- .quadrature_weighted_sum(weights, (mu - mean_K)^2)
  var_K <- within_alpha + between_alpha

  support_tolerance <- 64 * .Machine$double.eps * max(1, J)
  if (!is.finite(mean_K) || mean_K < 1 - support_tolerance ||
      mean_K > J + support_tolerance ||
      !is.finite(within_alpha) || within_alpha < 0 ||
      !is.finite(between_alpha) || between_alpha < 0 ||
      !is.finite(var_K) || var_K < 0) {
    .dpprior_abort_invalid(
      paste(
        "marginal moment quadrature violated finite support or",
        "non-negative variance-decomposition requirements"
      ),
      c("dpprior_marginal_moment_error", "dpprior_numerical_error"),
      "moments",
      c(
        mean = mean_K, var = var_K,
        within_alpha = within_alpha, between_alpha = between_alpha
      ),
      sprintf(
        "finite mean in [1,%d] and non-negative finite variance components",
        J
      ),
      "numerical_contract"
    )
  }

  list(
    mean = as.numeric(mean_K),
    var = as.numeric(var_K),
    within_alpha = as.numeric(within_alpha),
    between_alpha = as.numeric(between_alpha),
    node_rule = .quadrature_metadata(quad)
  )
}

#' Marginal Moments of K_J under Gamma Prior
#'
#' Computes the marginal mean \eqn{E[K_J]} and variance \eqn{Var(K_J)}
#' when the DP concentration parameter follows a Gamma(a, b) prior.
#'
#' @param J Integer; sample size (positive integer >= 1).
#' @param a Numeric; shape parameter of Gamma prior (> 0).
#' @param b Numeric; rate parameter of Gamma prior (> 0).
#' @param M Integer; number of quadrature nodes (default: 80).
#' @param M_verify Optional integer quadrature order at least
#'   \code{max(2*M, M+40)} and no greater than 512. When supplied, the
#'   selected-order moments are compared with a fresh higher-order computation.
#'   When omitted, the result is explicitly marked as an unverified
#'   approximation. For \code{M > 256}, no admissible verification order exists;
#'   only omitted \code{M_verify} with \code{strict = FALSE} is allowed.
#' @param abs_tol,rel_tol Non-negative absolute and relative tolerances used
#'   for selected-versus-verification comparisons.
#' @param strict Logical; if \code{TRUE}, require higher-order verification and
#'   raise a typed convergence error when it is absent or fails.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{\code{mean}}{Marginal mean \eqn{E[K_J | a, b]}}
#'     \item{\code{var}}{Marginal variance \eqn{Var(K_J | a, b)}}
#'     \item{\code{sd}}{Marginal standard deviation}
#'     \item{\code{cv}}{Coefficient of variation (sd/mean)}
#'     \item{\code{status}}{Either \code{"converged"} after successful
#'       higher-order verification or \code{"approximate"}.}
#'     \item{\code{decomposition}}{Within-alpha and between-alpha components
#'       of the law of total variance.}
#'     \item{\code{quadrature}}{Selected/verification order, tolerances,
#'       discrepancies, and node-rule provenance.}
#'   }
#'
#' @details
#' Uses Gauss-Laguerre quadrature to numerically evaluate:
#' \deqn{M_1(a,b) = E_{\alpha \sim \Gamma(a,b)}[\mu_J(\alpha)]}
#' \deqn{V(a,b) = E[v_J(\alpha)] + E[\mu_J(\alpha)^2] - M_1^2}
#'
#' where \eqn{\mu_J(\alpha)} and \eqn{v_J(\alpha)} are the conditional
#' mean and variance from Module 03.
#'
#' The Law of Total Variance decomposes the marginal variance into:
#' \itemize{
#'   \item Within-alpha variance: \eqn{E[v_J(\alpha)]}
#'   \item Between-alpha variance: \eqn{Var(\mu_J(\alpha))}
#' }
#'
#' The function name is retained for compatibility. The returned values use a
#' fixed Gauss--Laguerre rule and are therefore numerical approximations.
#' Supply \code{M_verify} to obtain an explicit higher-order agreement check;
#' the selected-order result is always returned and is never silently replaced.
#' The quadrature metadata reports the required verification order and whether
#' that order is available under the 512-node implementation ceiling.
#'
#' \strong{Key properties:}
#' \itemize{
#'   \item The mean is bounded: \eqn{1 \leq E[K_J] \leq J}
#'   \item Gamma mixing adds the non-negative between-\eqn{\alpha} component
#'     to the average conditional variance; this identity alone does not imply
#'     a universal marginal variance-to-mean ratio
#'   \item The marginal variance equals the average conditional variance plus
#'     the non-negative between-\eqn{\alpha} variance of the conditional mean
#' }
#'
#' @examples
#' # Example: J=50, Gamma(1.5, 0.5) prior
#' result <- exact_K_moments(50, 1.5, 0.5)
#' print(result)
#'
#' # Verify the law-of-total-variance decomposition
#' sum(unlist(result$decomposition)) - result$var
#'
#' # Verify mean bounds
#' 1 <= result$mean && result$mean <= 50  # TRUE
#'
#' @seealso \code{\link{K_moments}} for convenience wrapper,
#'   \code{\link{mean_K_given_alpha}}, \code{\link{var_K_given_alpha}}
#'
#' @references
#' Antoniak, C. E. (1974). Mixtures of Dirichlet Processes.
#' \emph{The Annals of Statistics}, 2(6), 1152-1174.
#'
#' @family marginal_K
#'
#' @export
exact_K_moments <- function(
    J, a, b, M = .QUAD_NODES_DEFAULT,
    M_verify = NULL, abs_tol = 1e-10, rel_tol = 1e-8,
    strict = FALSE) {
  # Input validation
  assert_valid_J(J)
  a <- .dpprior_validate_scalar(
    a, "a", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  b <- .dpprior_validate_scalar(
    b, "b", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  controls <- .marginal_verification_controls(
    M, M_verify, abs_tol, rel_tol, strict
  )

  selected <- .marginal_moments_fixed(J, a, b, controls$M)
  verification <- NULL
  verification_passed <- NA
  reason <- "fixed_order_unverified"
  mean_difference <- NA_real_
  variance_difference <- NA_real_
  mean_tolerance <- NA_real_
  variance_tolerance <- NA_real_

  if (!is.null(controls$M_verify)) {
    verification <- .marginal_moments_fixed(
      J, a, b, controls$M_verify
    )
    mean_difference <- abs(selected$mean - verification$mean)
    variance_difference <- abs(selected$var - verification$var)
    mean_tolerance <- controls$abs_tol + controls$rel_tol *
      max(abs(selected$mean), abs(verification$mean))
    variance_tolerance <- controls$abs_tol + controls$rel_tol *
      max(abs(selected$var), abs(verification$var))
    verification_passed <- is.finite(mean_difference) &&
      is.finite(variance_difference) &&
      mean_difference <= mean_tolerance &&
      variance_difference <= variance_tolerance
    reason <- if (verification_passed) {
      "higher_order_agreement"
    } else {
      "higher_order_disagreement"
    }
  }

  status <- if (isTRUE(verification_passed)) "converged" else "approximate"
  if (controls$strict && !isTRUE(verification_passed)) {
    .dpprior_abort_invalid(
      sprintf(
        "marginal moments did not meet higher-order tolerance (mean difference %.3g; variance difference %.3g)",
        mean_difference, variance_difference
      ),
      "dpprior_marginal_convergence_error", "M_verify",
      controls$M_verify, "selected and higher-order agreement", reason
    )
  }

  sd_K <- sqrt(selected$var)
  cv_K <- sd_K / selected$mean

  list(
    mean = selected$mean,
    var = selected$var,
    sd = sd_K,
    cv = cv_K,
    status = status,
    decomposition = list(
      within_alpha = selected$within_alpha,
      between_alpha = selected$between_alpha
    ),
    quadrature = list(
      schema_version = 1L,
      engine = "gauss-laguerre",
      M_selected = controls$M,
      M_verification = if (is.null(controls$M_verify)) {
        NA_integer_
      } else {
        controls$M_verify
      },
      M_verification_required = controls$M_verification_required,
      verification_available = controls$verification_available,
      status = status,
      reason = reason,
      verification_performed = !is.null(controls$M_verify),
      verification_passed = verification_passed,
      absolute_tolerance = controls$abs_tol,
      relative_tolerance = controls$rel_tol,
      mean_difference = mean_difference,
      variance_difference = variance_difference,
      mean_tolerance = mean_tolerance,
      variance_tolerance = variance_tolerance,
      selected = selected$node_rule,
      verification = if (is.null(verification)) NULL else verification$node_rule
    )
  )
}


#' Convenience Wrapper for Marginal Moments
#'
#' Returns marginal mean and variance as a named numeric vector.
#'
#' @param J Integer; sample size (positive integer >= 1).
#' @param a Numeric; shape parameter of Gamma prior (> 0).
#' @param b Numeric; rate parameter of Gamma prior (> 0).
#' @param M Integer; number of quadrature nodes (default: 80).
#' @param M_verify Optional quadrature order satisfying the same verification
#'   contract as \code{exact_K_moments()}: at least \code{max(2*M, M+40)} and
#'   no greater than 512.
#' @param abs_tol,rel_tol Non-negative selected-versus-verification tolerances.
#' @param strict Logical; require successful higher-order verification.
#'
#' @return Named numeric vector \code{c(mean = ..., var = ...)} carrying a
#'   \code{"marginal_metadata"} attribute with status, variance decomposition,
#'   and quadrature verification details.
#'
#' @examples
#' K_moments(50, 2.0, 1.0)
#'
#' @seealso \code{\link{exact_K_moments}} for full output
#'
#' @family marginal_K
#'
#' @export
K_moments <- function(
    J, a, b, M = .QUAD_NODES_DEFAULT,
    M_verify = NULL, abs_tol = 1e-10, rel_tol = 1e-8,
    strict = FALSE) {
  result <- exact_K_moments(
    J, a, b, M, M_verify, abs_tol, rel_tol, strict
  )
  moments <- c(mean = result$mean, var = result$var)
  attr(moments, "marginal_metadata") <- list(
    schema_version = 1L,
    status = result$status,
    decomposition = result$decomposition,
    quadrature = result$quadrature
  )
  moments
}


# =============================================================================
# Diagnostic Functions
# =============================================================================

#' Variance Inflation Ratio
#'
#' Computes the ratio \eqn{Var(K_J) / E[K_J]} as an overdispersion measure.
#'
#' @param J Integer; sample size.
#' @param a Numeric; shape parameter of Gamma prior.
#' @param b Numeric; rate parameter of Gamma prior.
#' @param M Integer; number of quadrature nodes (default: 80).
#'
#' @return Numeric; the variance inflation ratio (VIR).
#'
#' @details
#' For a Poisson distribution, this ratio equals 1.
#' For the marginal distribution of \eqn{K_J} under a Gamma prior on
#' \eqn{\alpha}, this ratio is typically > 1, indicating overdispersion.
#'
#' This ratio is useful for:
#' \itemize{
#'   \item Diagnosing the appropriateness of Poisson approximations
#'   \item Comparing different prior specifications
#'   \item Understanding the "spread" induced by uncertainty in \eqn{\alpha}
#' }
#'
#' @examples
#' \dontrun{
#' # Typical overdispersion
#' vir <- variance_inflation_ratio(50, 2.0, 1.0)
#' vir > 1  # TRUE: overdispersed
#'
#' # Compare across prior specifications
#' variance_inflation_ratio(50, 2.0, 0.5)  # Higher uncertainty in alpha
#' variance_inflation_ratio(50, 8.0, 4.0)  # Same mean, lower variance in alpha
#'
#' }
#' @keywords internal
variance_inflation_ratio <- function(J, a, b, M = .QUAD_NODES_DEFAULT) {
  result <- exact_K_moments(J, a, b, M)
  if (result$mean > 0) {
    return(result$var / result$mean)
  } else {
    return(Inf)
  }
}


#' Compare Exact Moments to NegBin Approximation
#'
#' Compares exact marginal moments (via quadrature) to the Negative Binomial
#' approximation from the A1 method for error analysis.
#'
#' @param J Integer; sample size.
#' @param a Numeric; shape parameter of Gamma prior.
#' @param b Numeric; rate parameter of Gamma prior.
#' @param M Integer; number of quadrature nodes (default: 80).
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{exact}}{List with exact mean and var}
#'     \item{\code{negbin}}{List with NegBin approximation mean and var}
#'     \item{\code{abs_error}}{Absolute errors (negbin - exact)}
#'     \item{\code{rel_error}}{Relative errors}
#'   }
#'
#' @details
#' The NegBin approximation (A1 method from Lee, 2026, Section 3.1) assumes:
#' \deqn{K_J - 1 | \alpha \approx \text{Poisson}(\alpha \cdot c_J)}
#'
#' where \eqn{c_J = \log(J)}. With \eqn{\alpha \sim \text{Gamma}(a, b)}:
#' \deqn{E[K_J] \approx 1 + (a/b) \cdot c_J}
#' \deqn{Var(K_J) \approx m \cdot (1 + m/a), \quad m = (a/b) \cdot c_J}
#'
#' This comparison helps diagnose when the A1 approximation is insufficient
#' and exact A2 moment matching is needed.
#'
#' @examples
#' \dontrun{
#' # Large approximation error for small J
#' compare_to_negbin(50, 1.5, 0.5)
#'
#' # Error decreases with J
#' compare_to_negbin(300, 1.5, 0.5)
#'
#' }
#' @keywords internal
compare_to_negbin <- function(J, a, b, M = .QUAD_NODES_DEFAULT) {
  # Exact moments via quadrature
  exact <- exact_K_moments(J, a, b, M)

  # NegBin approximation (A1)
  c_J <- log(J)
  m <- (a / b) * c_J  # Mean of shifted distribution

  negbin_mean <- 1 + m
  negbin_var <- m * (1 + m / a)

  # Errors
  abs_error_mean <- negbin_mean - exact$mean
  abs_error_var <- negbin_var - exact$var

  rel_error_mean <- abs_error_mean / exact$mean
  rel_error_var <- if (exact$var > 0) abs_error_var / exact$var else Inf

  list(
    exact = list(mean = exact$mean, var = exact$var),
    negbin = list(mean = negbin_mean, var = negbin_var),
    abs_error = list(mean = abs_error_mean, var = abs_error_var),
    rel_error = list(mean = rel_error_mean, var = rel_error_var)
  )
}


# =============================================================================
# Advanced Functions for Jacobian/Newton Support
# =============================================================================

#' Marginal Moments with Jacobian
#'
#' Computes marginal moments and their Jacobian matrix with respect to
#' the Gamma hyperparameters (a, b). Used for Newton-type optimization.
#'
#' @param J Integer; sample size.
#' @param a Numeric; shape parameter of Gamma prior.
#' @param b Numeric; rate parameter of Gamma prior.
#' @param M Integer; number of quadrature nodes (default: 80).
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{mean}}{Marginal mean}
#'     \item{\code{var}}{Marginal variance}
#'     \item{\code{jacobian}}{2x2 matrix of partial derivatives}
#'     \item{\code{derivative_diagnostics}}{Diagnostics from the canonical
#'       centered-score implementation.}
#'     \item{\code{conditioning}}{Raw- and log-scale conditioning metadata.}
#'   }
#'
#' @details
#' This compatibility adapter delegates all numerical work to
#' \code{\link{moments_with_jacobian}}. It only restores the historical
#' row and column labels used by this internal helper; no second derivative
#' implementation is maintained here.
#'
#' The Jacobian matrix is:
#' \deqn{J = \begin{pmatrix}
#'   \partial M_1 / \partial a & \partial M_1 / \partial b \\
#'   \partial V / \partial a & \partial V / \partial b
#' \end{pmatrix}}
#'
#' Derivatives are computed using the score function identity:
#' \deqn{\frac{\partial}{\partial \theta} E[f(\alpha)] = E[f(\alpha) \cdot s_\theta(\alpha)]}
#'
#' where \eqn{s_\theta(\alpha) = \partial \log p(\alpha | a, b) / \partial \theta}.
#'
#' For Gamma(a, b):
#' \itemize{
#'   \item \eqn{s_a(\alpha) = \log(b) - \psi(a) + \log(\alpha)}
#'   \item \eqn{s_b(\alpha) = a/b - \alpha}
#' }
#'
#' @examples
#' \dontrun{
#' result <- marginal_moments_with_jacobian(50, 2.0, 1.0)
#' result$jacobian
#'
#' }
#' @seealso \code{\link{exact_K_moments}}, Module 07 (Jacobian)
#'
#' @keywords internal
marginal_moments_with_jacobian <- function(J, a, b, M = .QUAD_NODES_DEFAULT) {
  result <- moments_with_jacobian(J, a, b, M)
  dimnames(result$jacobian) <- list(c("mean", "var"), c("a", "b"))
  result
}


# =============================================================================
# Verification Functions
# =============================================================================

#' Verify Marginal Moments Properties
#'
#' Runs verification tests on the marginal moment computations.
#'
#' @param J Integer; sample size to test.
#' @param a Numeric; shape parameter to test.
#' @param b Numeric; rate parameter to test.
#' @param verbose Logical; if \code{TRUE}, print detailed results.
#'
#' @return Logical; \code{TRUE} if all verifications pass.
#'
#' @examples
#' \dontrun{
#' verify_marginal_moments(50, 2.0, 1.0)
#'
#' }
#' @keywords internal
verify_marginal_moments <- function(J, a, b, verbose = TRUE) {
  all_pass <- TRUE

  # Get marginal moments
  result <- exact_K_moments(J, a, b, M = 100)

  if (isTRUE(verbose)) {
    cat(sprintf("Marginal Moments Verification (J=%d, a=%.2f, b=%.2f):\n",
                J, a, b))
    cat(sprintf("  E[K_J] = %.6f\n", result$mean))
    cat(sprintf("  Var(K_J) = %.6f\n", result$var))
    cat(sprintf("  SD(K_J) = %.6f\n", result$sd))
    cat(sprintf("  CV(K_J) = %.6f\n\n", result$cv))
  }

  # Test 1: Mean bounds (1 <= E[K] <= J)
  test1 <- (result$mean >= 1 - 1e-10) && (result$mean <= J + 1e-10)
  if (isTRUE(verbose)) {
    cat(sprintf("  Test 1 (Mean bounds): E[K] in [1, %d]? %s\n",
                J, if (test1) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test1

  # Test 2: Variance non-negative
  test2 <- result$var >= 0
  if (isTRUE(verbose)) {
    cat(sprintf("  Test 2 (Var non-negative): Var(K) >= 0? %s\n",
                if (test2) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test2

  # Test 3: law-of-total-variance decomposition
  within <- result$decomposition$within_alpha
  between <- result$decomposition$between_alpha
  decomposition_error <- abs(result$var - within - between)
  test3 <- within >= 0 && between >= 0 && decomposition_error <= 1e-10
  if (isTRUE(verbose)) {
    cat(sprintf(
      paste0(
        "  Test 3 (Total variance): within=%.4f + between=%.4f, ",
        "error=%.2e? %s\n"
      ),
      within, between, decomposition_error, if (test3) "PASS" else "FAIL"
    ))
  }
  all_pass <- all_pass && test3

  # Test 4: CV is finite and non-negative (zero at J=1)
  test4 <- is.finite(result$cv) && result$cv >= 0
  if (isTRUE(verbose)) {
    cat(sprintf("  Test 4 (CV valid): finite CV >= 0? %s\n",
                if (test4) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test4

  if (isTRUE(verbose)) {
    cat(sprintf("\n  Overall: %s\n", if (all_pass) "PASS" else "FAIL"))
  }

  invisible(all_pass)
}


#' Report Quadrature Discrepancies for Marginal Moments
#'
#' Reports successive-order marginal-moment discrepancies and classifies each
#' selected/refined pair against an explicit mixed absolute/relative budget.
#' Successive discrepancies are not assumed to decrease monotonically because
#' quadrature error depends on the integrand and parameter regime.
#'
#' @param J Integer; sample size.
#' @param a Numeric; shape parameter.
#' @param b Numeric; rate parameter.
#' @param M_values Integer vector; numbers of quadrature nodes to test.
#' @param verbose Logical; if \code{TRUE}, print detailed results.
#' @param abs_tol,rel_tol Non-negative mixed-error tolerances.
#'
#' @return Data frame of successive-order discrepancies, mixed-error budgets,
#'   and pairwise budget classifications.
#'
#' @examples
#' \dontrun{
#' verify_quadrature_convergence(50, 1.5, 0.5)
#'
#' }
#' @keywords internal
verify_quadrature_convergence <- function(J, a, b,
                                          M_values = c(20, 40, 60, 80, 100, 120),
                                          verbose = TRUE,
                                          abs_tol = 1e-10,
                                          rel_tol = 1e-8) {
  abs_tol <- .dpprior_validate_scalar(
    abs_tol, "abs_tol", lower = 0,
    .subclass = "dpprior_control_error"
  )
  rel_tol <- .dpprior_validate_scalar(
    rel_tol, "rel_tol", lower = 0,
    .subclass = "dpprior_control_error"
  )
  results <- data.frame(
    M = integer(0),
    mean = numeric(0),
    var = numeric(0),
    mean_change = numeric(0),
    var_change = numeric(0),
    mean_tolerance = numeric(0),
    var_tolerance = numeric(0),
    within_budget = logical(0)
  )

  prev_mean <- NA
  prev_var <- NA

  for (M in M_values) {
    moments <- exact_K_moments(J, a, b, M)

    mean_change <- if (is.na(prev_mean)) NA else abs(moments$mean - prev_mean)
    var_change <- if (is.na(prev_var)) NA else abs(moments$var - prev_var)
    mean_tolerance <- if (is.na(prev_mean)) {
      NA_real_
    } else {
      abs_tol + rel_tol * max(abs(moments$mean), abs(prev_mean))
    }
    var_tolerance <- if (is.na(prev_var)) {
      NA_real_
    } else {
      abs_tol + rel_tol * max(abs(moments$var), abs(prev_var))
    }
    within_budget <- if (is.na(mean_change)) {
      NA
    } else {
      mean_change <= mean_tolerance && var_change <= var_tolerance
    }

    results <- rbind(results, data.frame(
      M = M,
      mean = moments$mean,
      var = moments$var,
      mean_change = mean_change,
      var_change = var_change,
      mean_tolerance = mean_tolerance,
      var_tolerance = var_tolerance,
      within_budget = within_budget
    ))

    prev_mean <- moments$mean
    prev_var <- moments$var
  }

  if (isTRUE(verbose)) {
    cat(sprintf("Quadrature Order Discrepancies (J=%d, a=%.2f, b=%.2f):\n",
                J, a, b))
    print(results, row.names = FALSE)
  }

  invisible(results)
}


#' Run All Module 05 Verification Tests
#'
#' Comprehensive verification suite for the marginal moments module.
#'
#' @param verbose Logical; if \code{TRUE}, print detailed results.
#'
#' @return Logical; \code{TRUE} if all tests pass.
#'
#' @examples
#' \dontrun{
#' verify_moments_marginal_all()
#'
#' }
#' @keywords internal
verify_moments_marginal_all <- function(verbose = TRUE) {
  if (isTRUE(verbose)) {
    cat("=" , rep("=", 69), "\n", sep = "")
    cat("Module 05: Marginal Moments - Full Verification Suite\n")
    cat("=" , rep("=", 69), "\n\n", sep = "")
  }

  all_pass <- TRUE

  # Test cases
  test_cases <- list(
    list(J = 50, a = 1.5, b = 0.5),
    list(J = 50, a = 2.0, b = 1.0),
    list(J = 100, a = 2.0, b = 1.0),
    list(J = 10, a = 0.5, b = 0.5),
    list(J = 10, a = 5.0, b = 2.0)
  )

  for (tc in test_cases) {
    result <- verify_marginal_moments(tc$J, tc$a, tc$b, verbose = verbose)
    all_pass <- all_pass && result
    if (isTRUE(verbose)) cat("\n")
  }

  # NegBin comparison
  if (isTRUE(verbose)) {
    cat("-", rep("-", 69), "\n", sep = "")
    cat("NegBin Approximation Error Analysis:\n")
    cat("-", rep("-", 69), "\n", sep = "")

    for (J in c(50, 100, 300)) {
      comp <- compare_to_negbin(J, 1.5, 0.5)
      cat(sprintf("  J=%3d: E[K] exact=%.2f, NegBin=%.2f, rel_err=%.1f%%\n",
                  J, comp$exact$mean, comp$negbin$mean,
                  comp$rel_error$mean * 100))
    }
    cat("\n")
  }

  # Convergence test
  if (isTRUE(verbose)) {
    cat("-", rep("-", 69), "\n", sep = "")
    cat("Quadrature Convergence Test:\n")
    cat("-", rep("-", 69), "\n", sep = "")
    verify_quadrature_convergence(50, 1.5, 0.5, verbose = TRUE)
    cat("\n")
  }

  if (isTRUE(verbose)) {
    cat("=" , rep("=", 69), "\n", sep = "")
    cat(sprintf("Overall Result: %s\n", if (all_pass) "ALL TESTS PASSED" else "SOME TESTS FAILED"))
    cat("=" , rep("=", 69), "\n", sep = "")
  }

  invisible(all_pass)
}
