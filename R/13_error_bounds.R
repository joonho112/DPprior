# =============================================================================
# Module 13: Error Bounds Implementation
# =============================================================================
#
# This module provides functions for quantifying approximation errors in the
# A1 large-J approximation, implementing the error quantification framework
# from Lee (2026, Section 3.3).
#
# The A1 approximation (shifted NegBin) differs from the exact K_J distribution
# due to three error sources:
# 1. Poissonization error: Bernoulli sum -> Poisson approximation
# 2. Mean-linearization error: Exact mean lambda_J(alpha) -> alpha*c_J
# 3. Mixing error: Poisson-Gamma identity -> shifted NegBin
#
# Author: JoonHo Lee
# Date: December 2025
# Part of: DPprior R Package
# Reference: Lee (2026), Section 3.3
# =============================================================================


# =============================================================================
# Poissonization Error Bound
# =============================================================================

#' Poissonization Error Bound (Raw Sum of Squared Probabilities)
#'
#' Computes the raw sum of squared Bernoulli probabilities:
#' \deqn{\sum_{i=2}^{J} p_i^2 = \alpha^2 [\psi_1(\alpha+1) - \psi_1(\alpha+J)]}
#' where \eqn{p_i = \alpha / (\alpha + i - 1)}.
#'
#' This quantity represents the "underdispersion gap" between the conditional
#' variance of \eqn{K_J | \alpha} and a Poisson with the same mean.
#'
#' @param J Integer; sample size (number of observations).
#' @param alpha Numeric; concentration parameter (can be vectorized).
#'
#' @return Numeric vector; sum of squared probabilities for each alpha value.
#'
#' @details
#' From the Poisson-binomial representation, \eqn{S_J = K_J - 1 = \sum_{i=2}^{J} I_i}
#' where \eqn{I_i \sim \text{Bernoulli}(p_i)}.
#'
#' This sum equals:
#' \deqn{\sum_{i=2}^{J} p_i^2 = \alpha^2 [\psi_1(\alpha+1) - \psi_1(\alpha+J)]}
#' using the identity for sums of squared reciprocals.
#'
#' @seealso \code{\link{compute_poissonization_bound}} for the full Chen-Stein bound
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @examples
#' \dontrun{
#' # Compute raw sum for J=50, alpha=1
#' compute_sum_p_squared(J = 50, alpha = 1)
#'
#' # Vectorized over alpha
#' compute_sum_p_squared(J = 50, alpha = c(0.5, 1, 2))
#'
#' }
#' @keywords internal
compute_sum_p_squared <- function(J, alpha) {
  assert_valid_J(J)
  assert_positive(alpha, "alpha")

  alpha^2 * (trigamma(alpha + 1) - trigamma(alpha + J))
}


#' Poissonization Error Bound (Chen-Stein / Le Cam Bound)
#'
#' Computes an upper bound on the conditional total variation distance
#' between the shifted cluster count \eqn{S_J = K_J - 1} and a Poisson law
#' with the same mean (Poissonization error).
#'
#' @param J Integer; sample size (number of observations).
#' @param alpha Numeric; concentration parameter (can be vectorized).
#' @param raw Logical; if TRUE, return just sum(p_i^2) without the prefactor.
#'   Default is FALSE.
#'
#' @return Numeric vector; upper bound on
#'   \eqn{d_{TV}(S_J | \alpha, \text{Poisson}(\lambda_J(\alpha)))}.
#'
#' @details
#' Under the CRP representation,
#' \eqn{S_J = \sum_{i=2}^J I_i} where \eqn{I_i \sim \text{Bernoulli}(p_i)} and
#' \eqn{p_i = \alpha / (\alpha + i - 1)}.
#'
#' A standard Chen-Stein/Le Cam bound gives:
#' \deqn{d_{TV}(S_J, \text{Poisson}(\lambda)) \le
#'       \frac{1 - e^{-\lambda}}{\lambda} \sum_{i=2}^J p_i^2}
#' where \eqn{\lambda = \sum_{i=2}^J p_i = E[S_J | \alpha]}.
#'
#' The prefactor \eqn{(1 - e^{-\lambda})/\lambda} is always in (0, 1] and
#' approaches 1 as \eqn{\lambda \to 0}. This provides a tighter bound than
#' simply using \eqn{\sum p_i^2} alone.
#'
#' The returned value is capped at 1 (since total variation is always between 0 and 1).
#'
#' @seealso \code{\link{compute_sum_p_squared}}, \code{\link{compute_linearization_bound}}
#'
#' @references
#' Le Cam, L. (1960). An approximation theorem for the Poisson binomial
#' distribution. \emph{Pacific Journal of Mathematics}, 10(4), 1181-1197.
#'
#' Chen, L. H. Y. (1975). Poisson approximation for dependent trials.
#' \emph{The Annals of Probability}, 3(3), 534-545.
#'
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @examples
#' \dontrun{
#' # Full Chen-Stein bound
#' compute_poissonization_bound(J = 50, alpha = 1)
#'
#' # Raw bound (sum of p_i^2)
#' compute_poissonization_bound(J = 50, alpha = 1, raw = TRUE)
#'
#' # Vectorized
#' compute_poissonization_bound(J = 50, alpha = c(0.5, 1, 2, 5))
#'
#' }
#' @keywords internal
compute_poissonization_bound <- function(J, alpha, raw = FALSE) {
  assert_valid_J(J)
  assert_positive(alpha, "alpha")

  # Sum_{i=2}^J p_i^2 in closed form via trigamma
  sum_p_sq <- compute_sum_p_squared(J, alpha)

  if (raw) {
    return(sum_p_sq)
  }


  # Shifted mean: lambda_J(alpha) = E[S_J | alpha] = E[K_J | alpha] - 1
  lambda <- pmax(0, mean_K_given_alpha(J, alpha) - 1)

  # Chen-Stein prefactor: (1 - exp(-lambda))/lambda
  # This equals 1 when lambda = 0 (by L'Hopital or Taylor expansion)
  # The prefactor is always in (0, 1], providing a tighter bound
  prefactor <- ifelse(lambda > 0, (1 - exp(-lambda)) / lambda, 1)

  # Cap at 1 (TV is bounded by 1)
  pmin(1, prefactor * sum_p_sq)
}


# =============================================================================
# Mean-Linearization Error Bound
# =============================================================================

#' Poisson-Poisson KL Divergence
#'
#' Computes the Kullback-Leibler divergence between two Poisson distributions:
#' \deqn{KL(\text{Poisson}(\lambda) || \text{Poisson}(\lambda')) =
#'       \lambda \log(\lambda/\lambda') + \lambda' - \lambda}
#'
#' @param lambda Numeric; mean of first Poisson distribution.
#' @param lambda_prime Numeric; mean of second Poisson distribution.
#'
#' @return Numeric; KL divergence (non-negative, possibly Inf).
#'
#' @details
#' Special cases:
#' \itemize{
#'   \item If both \eqn{\lambda = 0} and \eqn{\lambda' = 0}: KL = 0
#'   \item If \eqn{\lambda = 0} and \eqn{\lambda' > 0}: KL = \eqn{\lambda'}
#'   \item If \eqn{\lambda > 0} and \eqn{\lambda' = 0}: KL = Inf
#' }
#'
#' @keywords internal
poisson_kl_divergence <- function(lambda, lambda_prime) {
  n <- length(lambda)
  kl <- rep(0, n)

  # Both zero: KL = 0 (already initialized)

  # lambda = 0, lambda' > 0: KL = lambda'
  idx0 <- (lambda == 0) & (lambda_prime > 0)
  kl[idx0] <- lambda_prime[idx0]

  # lambda > 0, lambda' > 0: standard formula
  idxp <- (lambda > 0) & (lambda_prime > 0)
  kl[idxp] <- lambda[idxp] * log(lambda[idxp] / lambda_prime[idxp]) +
    (lambda_prime[idxp] - lambda[idxp])

  # lambda > 0, lambda' = 0: KL = Inf
  idx_inf <- (lambda > 0) & (lambda_prime == 0)
  kl[idx_inf] <- Inf

  kl
}


#' Mean-Linearization Error Bound
#'
#' Computes an upper bound on the TV distance between two Poisson distributions,
#' \eqn{\text{Poisson}(\lambda_J(\alpha))} and \eqn{\text{Poisson}(\alpha c_J)},
#' using the Poisson-Poisson KL divergence together with Pinsker's inequality.
#'
#' @param J Integer; sample size.
#' @param alpha Numeric; concentration parameter (vectorized).
#' @param cJ Numeric; scaling constant (default: log(J)).
#'
#' @return Numeric; upper bound via Pinsker's inequality.
#'
#' @details
#' Let \eqn{\lambda = \lambda_J(\alpha)} (exact shifted mean) and
#' \eqn{\lambda' = \alpha c_J} (A1 approximate mean).
#'
#' The KL divergence is:
#' \deqn{KL(\text{Poisson}(\lambda) || \text{Poisson}(\lambda')) =
#'       \lambda \log(\lambda/\lambda') + \lambda' - \lambda}
#'
#' By Pinsker's inequality:
#' \deqn{d_{TV}(\text{Poisson}(\lambda), \text{Poisson}(\lambda')) \le \sqrt{KL/2}}
#'
#' Numerical safeguards handle edge cases where \eqn{\lambda} or \eqn{c_J} is zero.
#'
#' @seealso \code{\link{compute_poissonization_bound}}, \code{\link{compute_total_tv_bound}}
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @examples
#' \dontrun{
#' # Linearization bound for J=50, alpha=1
#' compute_linearization_bound(J = 50, alpha = 1)
#'
#' # Effect of J on linearization bound (should decrease)
#' sapply(c(25, 50, 100, 200), function(J)
#'   compute_linearization_bound(J, alpha = 2))
#'
#' }
#' @keywords internal
compute_linearization_bound <- function(J, alpha, cJ = log(J)) {
  assert_valid_J(J)
  assert_positive(alpha, "alpha")

  # Validate cJ
  if (!is.numeric(cJ) || length(cJ) != 1L || !is.finite(cJ) || cJ < 0) {
    stop("cJ must be a finite non-negative scalar", call. = FALSE)
  }

  # Exact conditional mean for shifted count S_J = K_J - 1
  lambda_exact <- pmax(0, mean_K_given_alpha(J, alpha) - 1)
  lambda_approx <- alpha * cJ

  # Handle all edge cases properly
  both_zero <- (lambda_exact == 0) & (lambda_approx == 0)

  # Compute KL divergence with proper edge case handling
  kl_div <- poisson_kl_divergence(lambda_exact, lambda_approx)

  # Pinsker's inequality: d_TV <= sqrt(KL/2)
  out <- sqrt(0.5 * pmax(0, kl_div))

  # Both zero means identical distributions: TV = 0
  out[both_zero] <- 0

  # Cap at 1
  pmin(1, out)
}


# =============================================================================
# Total TV Bound (Conditional)
# =============================================================================

#' Total TV Error Bound (Conditional)
#'
#' Computes the combined conditional TV bound using the triangle inequality:
#' \deqn{d_{TV}(K_J | \alpha, 1 + \text{NegBin}) \le B_{\text{Pois}} + B_{\text{lin}}}
#'
#' The result is capped at 1 since TV distance is bounded by 1.
#'
#' @param J Integer; sample size.
#' @param alpha Numeric; concentration parameter (vectorized).
#' @param cJ Numeric; scaling constant (default: log(J)).
#'
#' @return Numeric; upper bound on total TV error (capped at 1).
#'
#' @details
#' From Lee (2026, Section 3.3, Theorem 1), the total conditional TV error decomposes as:
#' \enumerate{
#'   \item Poissonization error: \eqn{S_J | \alpha} vs \eqn{\text{Poisson}(\lambda_J(\alpha))}
#'   \item Linearization error: \eqn{\text{Poisson}(\lambda_J(\alpha))} vs \eqn{\text{Poisson}(\alpha c_J)}
#' }
#'
#' @seealso \code{\link{compute_poissonization_bound}}, \code{\link{compute_linearization_bound}},
#'   \code{\link{expected_tv_bound}}
#'
#' @examples
#' \dontrun{
#' # Total bound at alpha = E[alpha] under Gamma(2, 1)
#' compute_total_tv_bound(J = 50, alpha = 2)
#'
#' # Vectorized
#' compute_total_tv_bound(J = 50, alpha = c(0.5, 1, 2, 5))
#'
#' }
#' @keywords internal
compute_total_tv_bound <- function(J, alpha, cJ = log(J)) {
  assert_valid_J(J)
  assert_positive(alpha, "alpha")

  B_pois <- compute_poissonization_bound(J, alpha, raw = FALSE)
  B_lin <- compute_linearization_bound(J, alpha, cJ)

  # TV distance is bounded by 1
  pmin(1, B_pois + B_lin)
}


# =============================================================================
# A1 Moment Error
# =============================================================================

#' A1 Approximation Moment Errors
#'
#' Computes the discrepancy between the A1 (shifted NegBin) approximation
#' and exact marginal moments of \eqn{K_J}.
#'
#' @param J Integer; sample size.
#' @param a,b Numeric; Gamma hyperparameters (shape, rate).
#' @param cJ Numeric; scaling constant (default: log(J)).
#' @param M Integer; number of quadrature nodes (default: 80).
#'
#' @return A list with components:
#'   \describe{
#'     \item{exact_mean, exact_var}{Exact moments via Gauss-Laguerre quadrature}
#'     \item{a1_mean, a1_var}{A1 (shifted NegBin) approximation moments}
#'     \item{error_mean_abs, error_var_abs}{Absolute errors}
#'     \item{error_mean_rel, error_var_rel}{Relative errors (percentage)}
#'   }
#'
#' @details
#' The A1 approximation models \eqn{K_J \approx 1 + \text{NegBin}(a, p_J)} where
#' \eqn{p_J = b / (b + c_J)}.
#'
#' The NegBin(a, p) moments are:
#' \itemize{
#'   \item Mean: \eqn{a(1-p)/p}
#'   \item Variance: \eqn{a(1-p)/p^2}
#' }
#'
#' @seealso \code{\link{exact_K_moments}}, \code{\link{DPprior_error_bounds}}
#'
#' @examples
#' # Moment errors for J=50, Gamma(2, 1) prior
#' errors <- a1_moment_error(J = 50, a = 2, b = 1)
#' print(errors)
#'
#' # Compare A1 accuracy at different J values
#' sapply(c(25, 50, 100, 200), function(J) {
#'   err <- a1_moment_error(J, a = 2, b = 1)
#'   c(mean_err = err$error_mean_rel, var_err = err$error_var_rel)
#' })
#'
#' # Compare different scaling constants
#' a1_moment_error(J = 50, a = 2, b = 1, cJ = log(50))
#' a1_moment_error(J = 50, a = 2, b = 1, cJ = digamma(50) + 0.5772)
#'
#' @export
a1_moment_error <- function(J, a, b, cJ = log(J), M = .QUAD_NODES_DEFAULT) {
  assert_valid_J(J)
  assert_positive(a, "a")
  assert_positive(b, "b")

  # Validate cJ
  if (!is.numeric(cJ) || length(cJ) != 1L || !is.finite(cJ) || cJ < 0) {
    stop("cJ must be a finite non-negative scalar", call. = FALSE)
  }

  # Exact moments via quadrature
  exact <- exact_K_moments(J, a, b, M)

  # A1 approximation moments (shifted NegBin from Poisson-Gamma identity)
  # S := K_J - 1 ~ NegBin(a, pJ) with pJ = b / (b + cJ)
  pJ <- b / (b + cJ)

  # NegBin(a, pJ) moments: mean = a*(1-pJ)/pJ, var = a*(1-pJ)/pJ^2
  mu_S <- a * (1 - pJ) / pJ
  var_S <- a * (1 - pJ) / pJ^2

  a1_mean <- 1 + mu_S
  a1_var <- var_S

  # Compute errors
  error_mean_abs <- abs(exact$mean - a1_mean)
  error_var_abs <- abs(exact$var - a1_var)
  error_mean_rel <- if (exact$mean > .Machine$double.eps) {
    100 * error_mean_abs / exact$mean
  } else NA_real_
  error_var_rel <- if (exact$var > .Machine$double.eps) {
    100 * error_var_abs / exact$var
  } else NA_real_

  list(
    exact_mean = exact$mean,
    exact_var = exact$var,
    a1_mean = a1_mean,
    a1_var = a1_var,
    error_mean_abs = error_mean_abs,
    error_var_abs = error_var_abs,
    error_mean_rel = error_mean_rel,
    error_var_rel = error_var_rel
  )
}


# =============================================================================
# Expected (Marginal) TV Bound
# =============================================================================

#' Expected TV Bound Under Gamma Prior
#'
#' Integrates the conditional TV bound over \eqn{\alpha \sim \text{Gamma}(a, b)}
#' to obtain the marginal error bound.
#'
#' @param J Integer; sample size.
#' @param a,b Numeric; Gamma hyperparameters.
#' @param cJ Numeric; scaling constant (default: log(J)).
#' @param M Integer; number of quadrature nodes.
#'
#' @return Numeric; \eqn{E[d_{TV} \text{ bound} | a, b]}.
#'
#' @details
#' From Lee (2026, Section 3.3, Corollary 1), the TV error between the exact prior predictive
#' \eqn{p(S_J | a, b)} and the A1 shifted NegBin proxy is bounded by:
#' \deqn{d_{TV}(P^{\text{exact}}, Q^{A1}) \le E_{\alpha \sim \Gamma(a,b)}[B_{\text{Pois}} + B_{\text{lin}}]}
#'
#' This follows from the mixture contraction property of TV distance.
#'
#' @seealso \code{\link{compute_total_tv_bound}}, \code{\link{integrate_gamma}}
#'
#' @examples
#' \dontrun{
#' # Marginal TV bound for J=100, Gamma(1, 1)
#' expected_tv_bound(J = 100, a = 1, b = 1)
#'
#' }
#' @keywords internal
expected_tv_bound <- function(J, a, b, cJ = log(J), M = .QUAD_NODES_DEFAULT) {
  assert_valid_J(J)
  assert_positive(a, "a")
  assert_positive(b, "b")

  integrate_gamma(
    function(alpha) compute_total_tv_bound(J, alpha, cJ),
    a, b, M
  )
}


# =============================================================================
# Threshold J Finder
# =============================================================================

#' Find Threshold J for A1 Adequacy
#'
#' Determines the minimum sample size J for which the A1 approximation
#' achieves target accuracy in moment matching.
#'
#' @param a,b Numeric; Gamma hyperparameters.
#' @param target_error Numeric; target relative error for mean (default: 5%).
#' @param target_var_error Numeric; target relative error for variance
#'   (default: 2 * target_error).
#' @param J_min,J_max Integer; search range for J.
#' @param step Integer; step size for search.
#'
#' @return Integer; minimum J achieving target accuracy, or NA if not found.
#'
#' @details
#' Searches over J values to find the smallest J where:
#' \itemize{
#'   \item Mean relative error < target_error
#'   \item Variance relative error < target_var_error
#' }
#'
#' Note: For many parameter combinations, especially with high \eqn{E[\alpha]},
#' the A1 approximation may never achieve low errors within practical J ranges.
#' In such cases, A2 refinement is recommended.
#'
#' @seealso \code{\link{a1_moment_error}}, \code{\link{DPprior_error_bounds}}
#'
#' @examples
#' \dontrun{
#' # Find threshold for 5% mean error
#' find_a1_threshold_J(a = 1, b = 2)
#'
#' # For higher E[alpha], threshold may not exist
#' find_a1_threshold_J(a = 2, b = 1)  # Likely returns NA
#' }
#'
#' @keywords internal
find_a1_threshold_J <- function(a, b, target_error = 0.05,
                                target_var_error = NULL,
                                J_min = 10, J_max = 500, step = 10) {
  assert_positive(a, "a")
  assert_positive(b, "b")
  assert_positive(target_error, "target_error")

  if (is.null(target_var_error)) {
    target_var_error <- 2 * target_error
  }

  for (J in seq(J_min, J_max, by = step)) {
    errors <- a1_moment_error(J, a, b)
    if (!is.na(errors$error_mean_rel) && !is.na(errors$error_var_rel) &&
        errors$error_mean_rel / 100 < target_error &&
        errors$error_var_rel / 100 < target_var_error) {
      return(as.integer(J))
    }
  }

  NA_integer_
}


# =============================================================================
# Main User-Facing Function
# =============================================================================

#' Compute A1 Approximation Error Bounds
#'
#' Comprehensive error analysis for the A1 large-J approximation.
#' Implements the error quantification framework from Lee (2026, Section 3.3).
#'
#' @param J Integer; sample size.
#' @param a,b Numeric; Gamma hyperparameters (shape, rate).
#' @param cJ Numeric; scaling constant (default: log(J)).
#' @param M Integer; number of quadrature nodes (default: 80).
#'
#' @return An S3 object of class "DPprior_error_bounds" with components:
#'   \describe{
#'     \item{J, a, b, cJ}{Input parameters}
#'     \item{moment_errors}{List of moment error metrics from \code{a1_moment_error}}
#'     \item{tv_bounds}{List with conditional and marginal TV bounds}
#'     \item{recommendation}{"A1_sufficient" or "A2_recommended"}
#'     \item{threshold_J}{Estimated J threshold for A1 adequacy (or NA)}
#'   }
#'
#' @details
#' The recommendation is based on:
#' \itemize{
#'   \item A1 sufficient if: mean relative error < 5\% AND variance relative error < 10\%
#'   \item A2 recommended otherwise
#' }
#'
#' This function provides:
#' \enumerate{
#'   \item Moment errors: Exact vs A1 approximation for mean and variance
#'   \item TV bounds: Conditional bounds at multiple alpha values, plus marginal bound
#'   \item Recommendation: Whether to use A1 or refine with A2
#'   \item Threshold: Minimum J for A1 adequacy with current prior
#' }
#'
#' @seealso \code{\link{a1_moment_error}}, \code{\link{compute_total_tv_bound}},
#'   \code{\link{DPprior_a1}}, \code{\link{DPprior_a2_newton}}
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @family diagnostics
#'
#' @examples
#' # Check A1 adequacy for J=50, typical prior
#' bounds <- DPprior_error_bounds(J = 50, a = 1.6, b = 1.2)
#' print(bounds)
#'
#' # For larger J, A1 becomes more adequate
#' bounds_200 <- DPprior_error_bounds(J = 200, a = 1.6, b = 1.2)
#' print(bounds_200)
#'
#' @export
DPprior_error_bounds <- function(J, a, b, cJ = log(J), M = .QUAD_NODES_DEFAULT) {
  assert_valid_J(J)
  assert_positive(a, "a")
  assert_positive(b, "b")

  # Validate cJ
  if (!is.numeric(cJ) || length(cJ) != 1L || !is.finite(cJ) || cJ < 0) {
    stop("cJ must be a finite non-negative scalar", call. = FALSE)
  }

  # Moment errors
  moment_errors <- a1_moment_error(J, a, b, cJ = cJ, M = M)

  # TV bounds at selected alpha values (around E[alpha] = a/b)
  E_alpha <- a / b
  alpha_multipliers <- c(0.1, 0.25, 0.5, 1, 2, 3, 5)
  alpha_vals <- alpha_multipliers * E_alpha
  # Filter to reasonable range
  alpha_vals <- alpha_vals[alpha_vals > 0.01 & alpha_vals < 100]

  cond_bounds <- data.frame(
    alpha = alpha_vals,
    poissonization_raw = sapply(alpha_vals, function(al)
      compute_poissonization_bound(J, al, raw = TRUE)),
    poissonization = sapply(alpha_vals, function(al)
      compute_poissonization_bound(J, al, raw = FALSE)),
    linearization = sapply(alpha_vals, function(al)
      compute_linearization_bound(J, al, cJ)),
    total = sapply(alpha_vals, function(al)
      compute_total_tv_bound(J, al, cJ))
  )

  # Marginal bound
  marginal_bound <- expected_tv_bound(J, a, b, cJ, M)

  # Recommendation based on moment errors
  a1_sufficient <- !is.na(moment_errors$error_mean_rel) &&
    !is.na(moment_errors$error_var_rel) &&
    moment_errors$error_mean_rel < 5 &&
    moment_errors$error_var_rel < 10
  recommendation <- if (a1_sufficient) "A1_sufficient" else "A2_recommended"

  # Threshold J estimation
  threshold_J <- find_a1_threshold_J(a, b, target_error = 0.05)

  result <- list(
    J = J,
    a = a,
    b = b,
    cJ = cJ,
    moment_errors = moment_errors,
    tv_bounds = list(
      conditional = cond_bounds,
      marginal = marginal_bound
    ),
    recommendation = recommendation,
    threshold_J = threshold_J
  )

  class(result) <- "DPprior_error_bounds"
  result
}


# =============================================================================
# S3 Print Method
# =============================================================================

#' Print Method for DPprior_error_bounds Objects
#'
#' Displays a formatted summary of the A1 approximation error analysis.
#'
#' @param x An object of class "DPprior_error_bounds".
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.DPprior_error_bounds <- function(x, ...) {
  cat("DPprior A1 Approximation Error Analysis\n")
  cat(strrep("=", 50), "\n\n")

  cat(sprintf("Sample size J = %d, c_J = %.4f\n", x$J, x$cJ))
  cat(sprintf("Gamma prior: alpha ~ Gamma(%.4f, %.4f) [shape-rate]\n", x$a, x$b))
  cat(sprintf("E[alpha] = %.4f, CV(alpha) = %.4f\n\n", x$a / x$b, 1 / sqrt(x$a)))

  cat("Moment Errors (A1 vs Exact):\n")
  cat(strrep("-", 45), "\n")
  cat(sprintf("  E[K_J]:   exact = %8.4f, A1 = %8.4f, error = %6.2f%%\n",
              x$moment_errors$exact_mean, x$moment_errors$a1_mean,
              x$moment_errors$error_mean_rel))
  cat(sprintf("  Var(K_J): exact = %8.4f, A1 = %8.4f, error = %6.2f%%\n",
              x$moment_errors$exact_var, x$moment_errors$a1_var,
              x$moment_errors$error_var_rel))

  cat("\nTV Bounds:\n")
  cat(strrep("-", 45), "\n")
  cat(sprintf("  Marginal E[d_TV] <= %.4f\n", x$tv_bounds$marginal))
  cat(sprintf("  At E[alpha] = %.2f:\n", x$a / x$b))

  # Find row closest to E[alpha]
  E_alpha <- x$a / x$b
  idx <- which.min(abs(x$tv_bounds$conditional$alpha - E_alpha))
  if (length(idx) > 0) {
    row <- x$tv_bounds$conditional[idx, ]
    cat(sprintf("    Poissonization: %.4f (raw: %.4f)\n",
                row$poissonization, row$poissonization_raw))
    cat(sprintf("    Linearization:  %.4f\n", row$linearization))
    cat(sprintf("    Total:          %.4f\n", row$total))
  }

  cat(sprintf("\nRecommendation: %s\n", x$recommendation))
  if (!is.na(x$threshold_J)) {
    cat(sprintf("  (A1 adequate for J >= %d with this prior)\n", x$threshold_J))
  } else {
    cat("  (A1 may not achieve < 5%% mean error for J <= 500 with this prior)\n")
  }

  invisible(x)
}


# =============================================================================
# Summary Method
# =============================================================================

#' Summary Method for DPprior_error_bounds Objects
#'
#' Provides a detailed summary including conditional bounds at multiple alpha values.
#'
#' @param object An object of class "DPprior_error_bounds".
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
summary.DPprior_error_bounds <- function(object, ...) {
  print(object)

  cat("\nConditional TV Bounds at Various alpha:\n")
  cat(strrep("-", 70), "\n")
  cat(sprintf("%10s %12s %12s %12s %12s\n",
              "alpha", "Pois(raw)", "Pois(C-S)", "Linear", "Total"))
  cat(strrep("-", 70), "\n")

  for (i in seq_len(nrow(object$tv_bounds$conditional))) {
    row <- object$tv_bounds$conditional[i, ]
    cat(sprintf("%10.3f %12.4f %12.4f %12.4f %12.4f\n",
                row$alpha, row$poissonization_raw, row$poissonization,
                row$linearization, row$total))
  }

  invisible(object)
}


# =============================================================================
# Utility Function: Error Landscape
# =============================================================================

#' Compute A1 Error Landscape
#'
#' Computes error metrics over a grid of (J, alpha) values for visualization.
#'
#' @param J_seq Numeric vector; sequence of J values.
#' @param alpha_seq Numeric vector; sequence of alpha values.
#' @param cJ_fun Function; scaling constant function (default: log).
#'
#' @return A data frame with columns: J, alpha, lambda_exact, lambda_approx,
#'   pois_raw, pois_bound, lin_bound, total_tv.
#'
#' @details
#' This function is useful for creating error landscape visualizations
#' as shown in Lee (2026, Section 3.3).
#'
#' @examples
#' # Create error landscape
#' landscape <- compute_error_landscape(
#'   J_seq = c(25, 50, 100),
#'   alpha_seq = c(0.5, 1, 2, 5)
#' )
#' print(landscape)
#'
#' @family diagnostics
#'
#' @export
compute_error_landscape <- function(J_seq, alpha_seq, cJ_fun = log) {
  results <- expand.grid(J = J_seq, alpha = alpha_seq)

  results$lambda_exact <- mapply(function(J, alpha) {
    mean_K_given_alpha(J, alpha) - 1
  }, results$J, results$alpha)

  results$lambda_approx <- mapply(function(J, alpha) {
    alpha * cJ_fun(J)
  }, results$J, results$alpha)

  results$pois_raw <- mapply(function(J, alpha) {
    compute_poissonization_bound(J, alpha, raw = TRUE)
  }, results$J, results$alpha)

  results$pois_bound <- mapply(function(J, alpha) {
    compute_poissonization_bound(J, alpha, raw = FALSE)
  }, results$J, results$alpha)

  results$lin_bound <- mapply(function(J, alpha) {
    compute_linearization_bound(J, alpha, cJ_fun(J))
  }, results$J, results$alpha)

  results$total_tv <- mapply(function(J, alpha) {
    compute_total_tv_bound(J, alpha, cJ_fun(J))
  }, results$J, results$alpha)

  results
}
