# =============================================================================
# Module 13: Error Bounds Implementation
# =============================================================================
#
# This module provides functions for quantifying approximation errors in the
# A1 large-J approximation, implementing the error quantification framework
# from Lee (2026, Section 3.3).
#
# The A1 approximation (shifted NegBin after Gamma mixing) differs from the
# exact K_J distribution through two conditional approximation steps:
# 1. Poissonization error: Bernoulli sum -> Poisson approximation
# 2. Mean-linearization error: Exact mean lambda_J(alpha) -> alpha*c_J
# Gamma mixing is exact for the proxy and propagates, rather than creates, the
# two conditional errors.
#
# Author: JoonHo Lee (jlee296@ua.edu)
# Date: December 2025
# Part of: DPprior R Package
# Reference: Lee (2026), Section 3.3
# =============================================================================


# =============================================================================
# Poissonization Error Bound
# =============================================================================

# Evaluate lambda_J(alpha) = sum_{r=1}^{J-1} alpha/(alpha+r) directly.
# Subtracting one from E[K_J|alpha] loses the shifted mean when alpha is tiny;
# the finite sum is also the independent Bernoulli representation used by the
# manuscript theorem.
.a1_shifted_mean <- function(J, alpha) {
  if (J == 1L) return(rep(0, length(alpha)))
  offsets <- seq_len(J - 1L)
  vapply(alpha, function(value) {
    sum(value / (value + offsets))
  }, numeric(1))
}


# TV bounds are mathematical probabilities.  Permit only a small floating
# rounding excursion before projection; clipping an invalid computation must
# never manufacture a seemingly valid bound.
.a1_clamp_tv_bound <- function(x, quantity = "TV bound") {
  rounding_tol <- 64 * .Machine$double.eps
  if (any(!is.finite(x)) || any(x < -rounding_tol) ||
      any(x > 1 + rounding_tol)) {
    stop(.dpprior_new_condition(
      sprintf("%s must be finite and in [0, 1]", quantity),
      classes = c(
        "dpprior_tv_bound_error", "dpprior_numerical_error",
        "dpprior_error", "error"
      ),
      quantity = quantity,
      value = x,
      expected = "finite value in [0, 1]",
      code = "invalid_tv_bound"
    ))
  }
  pmin(1, pmax(0, x))
}


.a1_validate_cJ <- function(cJ) {
  .dpprior_validate_scalar(
    cJ, "cJ", lower = 0,
    .subclass = "dpprior_a1_scaling_error"
  )
}


.a1_logspace_add <- function(log_x, log_y) {
  if (log_x >= log_y) {
    log_x + log1p(exp(log_y - log_x))
  } else {
    log_y + log1p(exp(log_x - log_y))
  }
}


.a1_nonnegative_from_log <- function(log_value, quantity) {
  if (identical(log_value, -Inf)) return(0)
  if (!is.finite(log_value) || log_value > log(.Machine$double.xmax)) {
    stop(.dpprior_new_condition(
      sprintf("%s is outside the representable finite range", quantity),
      classes = c(
        "dpprior_a1_moment_error", "dpprior_numerical_error",
        "dpprior_error", "error"
      ),
      quantity = quantity,
      log_value = log_value,
      code = "unrepresentable_a1_moment"
    ))
  }
  value <- exp(log_value)
  if (!is.finite(value) || value < 0) {
    stop(.dpprior_new_condition(
      sprintf("%s could not be evaluated as a finite non-negative value", quantity),
      classes = c(
        "dpprior_a1_moment_error", "dpprior_numerical_error",
        "dpprior_error", "error"
      ),
      quantity = quantity,
      log_value = log_value,
      value = value,
      code = "invalid_a1_moment"
    ))
  }
  value
}


.a1_proxy_moments <- function(a, b, cJ) {
  if (cJ == 0) return(list(mean_shifted = 0, variance = 0))
  log_a <- log(a)
  log_b <- log(b)
  log_cJ <- log(cJ)
  log_b_plus_cJ <- .a1_logspace_add(log_b, log_cJ)
  list(
    mean_shifted = .a1_nonnegative_from_log(
      log_a + log_cJ - log_b, "A1 shifted mean"
    ),
    variance = .a1_nonnegative_from_log(
      log_a + log_cJ + log_b_plus_cJ - 2 * log_b,
      "A1 variance"
    )
  )
}

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

  if (J == 1L) return(rep(0, length(alpha)))

  # The trigamma identity in the documentation is exact, but subtracting two
  # nearly equal trigamma values can return zero or NaN for large alpha.  The
  # direct finite Bernoulli sum is O(J), with J <= 500 in the supported domain,
  # and keeps every summand in [0, 1].
  offsets <- seq_len(J - 1L)
  vapply(alpha, function(value) {
    probabilities <- value / (value + offsets)
    sum(probabilities * probabilities)
  }, numeric(1))
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
#' It is also no larger than the submitted manuscript Appendix D bound
#' \eqn{\min(1,1/\lambda)\sum p_i^2} (Equation D10), so this implementation is
#' a strengthened Chen--Stein crosswalk rather than a different estimand.
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
  raw <- .dpprior_validate_control(raw, "raw", type = "logical")

  # Sum_{i=2}^J p_i^2 in closed form via trigamma
  sum_p_sq <- compute_sum_p_squared(J, alpha)

  if (raw) {
    return(sum_p_sq)
  }


  # Shifted mean: lambda_J(alpha) = E[S_J | alpha] = E[K_J | alpha] - 1
  lambda <- .a1_shifted_mean(J, alpha)

  # Chen-Stein prefactor: (1 - exp(-lambda))/lambda
  # This equals 1 when lambda = 0 (by L'Hopital or Taylor expansion)
  # The prefactor is always in (0, 1], providing a tighter bound
  prefactor <- rep(1, length(lambda))
  positive <- lambda > 0
  prefactor[positive] <- -expm1(-lambda[positive]) / lambda[positive]

  # Cap at 1 (TV is bounded by 1)
  .a1_clamp_tv_bound(
    pmin(1, prefactor * sum_p_sq), "Poissonization TV bound"
  )
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
  if (!is.numeric(lambda) || !is.numeric(lambda_prime) ||
      length(lambda) == 0L || length(lambda_prime) == 0L ||
      anyNA(lambda) || anyNA(lambda_prime) ||
      any(lambda < 0) || any(lambda_prime < 0) ||
      any(is.nan(lambda)) || any(is.nan(lambda_prime)) ||
      any(lambda == Inf) || any(lambda_prime == -Inf)) {
    .dpprior_abort_invalid(
      "lambda and lambda_prime must be non-empty non-negative numeric values; only lambda_prime may be +Inf",
      "dpprior_poisson_kl_error", "lambda/lambda_prime",
      list(lambda = lambda, lambda_prime = lambda_prime),
      "compatible non-negative Poisson means", "invalid_poisson_mean"
    )
  }
  n <- max(length(lambda), length(lambda_prime))
  if (!(length(lambda) %in% c(1L, n)) ||
      !(length(lambda_prime) %in% c(1L, n))) {
    .dpprior_abort_invalid(
      "lambda and lambda_prime must have equal lengths or one must be scalar",
      c("dpprior_poisson_kl_error", "dpprior_length_error"),
      "lambda/lambda_prime", c(length(lambda), length(lambda_prime)),
      "equal lengths or scalar recycling", "incompatible_lengths"
    )
  }
  lambda <- rep(lambda, length.out = n)
  lambda_prime <- rep(lambda_prime, length.out = n)
  kl <- rep(0, n)

  # Both zero: KL = 0 (already initialized)

  # lambda = 0, lambda' > 0: KL = lambda'
  idx0 <- (lambda == 0) & (lambda_prime > 0)
  kl[idx0] <- lambda_prime[idx0]

  # lambda > 0, finite lambda' > 0: use a centered expression.  The direct
  # formula subtracts nearly equal O(lambda) terms when the two means agree.
  idxp <- (lambda > 0) & is.finite(lambda_prime) & (lambda_prime > 0)
  lambda_positive <- lambda[idxp]
  lambda_prime_positive <- lambda_prime[idxp]
  difference <- lambda_prime_positive - lambda_positive
  near <- abs(difference) < 1e-4 * lambda_positive
  centered <- numeric(length(lambda_positive))
  if (any(near)) {
    u <- difference[near] / lambda_positive[near]
    series <- numeric(length(u))
    power <- u * u
    for (order in 2:12) {
      series <- series + if (order %% 2L == 0L) {
        power / order
      } else {
        -power / order
      }
      power <- power * u
    }
    centered[near] <- series
  }
  if (any(!near)) {
    # Log differences remain finite even when the ratio itself would overflow
    # or when (lambda_prime-lambda)/lambda rounds to exactly -1.
    centered[!near] <-
      log(lambda_positive[!near]) - log(lambda_prime_positive[!near]) +
      lambda_prime_positive[!near] / lambda_positive[!near] - 1
  }
  kl[idxp] <- lambda_positive * centered

  # A finite first mean and infinite second mean have infinite divergence.
  idx_prime_inf <- is.finite(lambda) & is.infinite(lambda_prime)
  kl[idx_prime_inf] <- Inf

  # lambda > 0, lambda' = 0: KL = Inf
  idx_inf <- (lambda > 0) & (lambda_prime == 0)
  kl[idx_inf] <- Inf

  # Protect only against a few ulps below zero from the centered calculation.
  tiny_negative <- kl < 0 & kl >= -64 * .Machine$double.eps *
    pmax(1, lambda, ifelse(is.finite(lambda_prime), lambda_prime, 1))
  kl[tiny_negative] <- 0
  if (any(kl < 0 | is.nan(kl))) {
    stop(.dpprior_new_condition(
      "Poisson KL calculation produced an invalid negative or NaN value",
      classes = c(
        "dpprior_poisson_kl_error", "dpprior_numerical_error",
        "dpprior_error", "error"
      ),
      value = kl,
      code = "invalid_kl"
    ))
  }

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
  cJ <- .a1_validate_cJ(cJ)

  # Exact conditional mean for shifted count S_J = K_J - 1
  lambda_exact <- .a1_shifted_mean(J, alpha)
  lambda_approx <- if (cJ == 0) {
    rep(0, length(alpha))
  } else {
    # Explicitly represent overflow as an infinite proxy mean; the KL and TV
    # limits are then Inf and 1 rather than NaN.
    ifelse(alpha > .Machine$double.xmax / cJ, Inf, alpha * cJ)
  }

  # Handle all edge cases properly
  both_zero <- (lambda_exact == 0) & (lambda_approx == 0)

  # Compute KL divergence with proper edge case handling
  kl_div <- poisson_kl_divergence(lambda_exact, lambda_approx)

  # Pinsker's inequality: d_TV <= sqrt(KL/2)
  out <- sqrt(0.5 * kl_div)

  # Both zero means identical distributions: TV = 0
  out[both_zero] <- 0

  # Infinite KL is a finite, valid but uninformative TV bound of one.
  out[is.infinite(out)] <- 1
  .a1_clamp_tv_bound(pmin(1, out), "mean-linearization TV bound")
}


# =============================================================================
# Total TV Bound (Conditional)
# =============================================================================

#' Total TV Error Bound (Conditional)
#'
#' Computes the combined conditional TV bound using the triangle inequality:
#' \deqn{d_{TV}(\mathcal{L}(S_J | \alpha),
#'              \text{Poisson}(\alpha c_J))
#'       \le B_{\text{Pois}} + B_{\text{lin}}}
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
#' In the submitted manuscript this is Appendix D, Corollary
#' \code{cor:combined-conditional} (Equations D14--D16). The total conditional
#' TV error decomposes as:
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
  cJ <- .a1_validate_cJ(cJ)

  B_pois <- compute_poissonization_bound(J, alpha, raw = FALSE)
  B_lin <- compute_linearization_bound(J, alpha, cJ)

  # TV distance is bounded by 1
  total <- .a1_clamp_tv_bound(
    pmin(1, B_pois + B_lin), "combined conditional TV bound"
  )
  ordering_tol <- 64 * .Machine$double.eps
  if (any(total + ordering_tol < B_pois) ||
      any(total + ordering_tol < B_lin)) {
    stop(.dpprior_new_condition(
      "combined TV bound violated component ordering",
      classes = c(
        "dpprior_tv_bound_error", "dpprior_numerical_error",
        "dpprior_error", "error"
      ),
      poissonization = B_pois,
      linearization = B_lin,
      total = total,
      code = "component_ordering"
    ))
  }
  total
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
#' @param M_verify Optional independent quadrature order for the exact moments.
#'   It must satisfy the package-wide verification-order contract.
#' @param abs_tol,rel_tol Non-negative selected-versus-verification tolerances.
#' @param strict Logical; require verified exact-moment quadrature.
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
a1_moment_error <- function(
    J, a, b, cJ = log(J), M = .QUAD_NODES_DEFAULT,
    M_verify = NULL, abs_tol = 1e-10, rel_tol = 1e-8,
    strict = FALSE) {
  assert_valid_J(J)
  a <- .dpprior_validate_scalar(
    a, "a", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  b <- .dpprior_validate_scalar(
    b, "b", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  cJ <- .a1_validate_cJ(cJ)

  # Exact moments via quadrature
  exact <- exact_K_moments(
    J, a, b, M = M, M_verify = M_verify,
    abs_tol = abs_tol, rel_tol = rel_tol, strict = strict
  )

  # A1 approximation moments (shifted NegBin from Poisson-Gamma identity)
  # S := K_J - 1 ~ NegBin(a, pJ) with pJ = b / (b + cJ)
  # Evaluate the algebraic forms in log space: direct products can overflow or
  # underflow even when the resulting moment is finite and representable.
  # NegBin(a, pJ), pJ=b/(b+cJ):
  # mean = a*cJ/b; variance = a*cJ*(b+cJ)/b^2.
  proxy_moments <- .a1_proxy_moments(a, b, cJ)
  mu_S <- proxy_moments$mean_shifted
  var_S <- proxy_moments$variance

  a1_mean <- 1 + mu_S
  a1_var <- var_S
  computed <- c(
    exact_mean = exact$mean, exact_var = exact$var,
    a1_mean = a1_mean, a1_var = a1_var
  )
  if (any(!is.finite(computed)) || any(computed < 0)) {
    stop(.dpprior_new_condition(
      "A1 moment calculation produced non-finite or negative moments",
      classes = c(
        "dpprior_a1_moment_error", "dpprior_numerical_error",
        "dpprior_error", "error"
      ),
      value = computed,
      code = "invalid_a1_moments"
    ))
  }

  # Compute errors
  error_mean_abs <- abs(exact$mean - a1_mean)
  error_var_abs <- abs(exact$var - a1_var)
  error_mean_rel <- if (exact$mean > .Machine$double.eps) {
    100 * error_mean_abs / exact$mean
  } else NA_real_
  error_var_rel <- if (exact$var > .Machine$double.eps) {
    100 * error_var_abs / exact$var
  } else if (error_var_abs <= .Machine$double.eps) {
    0
  } else {
    Inf
  }

  result <- list(
    exact_mean = exact$mean,
    exact_var = exact$var,
    a1_mean = a1_mean,
    a1_var = a1_var,
    error_mean_abs = error_mean_abs,
    error_var_abs = error_var_abs,
    error_mean_rel = error_mean_rel,
    error_var_rel = error_var_rel
  )
  attr(result, "a1_error_metadata") <- list(
    schema_version = 1L,
    status = exact$status,
    verified = identical(exact$status, "converged"),
    exact_moment_quadrature = exact$quadrature,
    proxy = list(
      conditional = "Poisson(alpha * cJ)",
      marginal = "1 + NegativeBinomial(a, b/(b+cJ))",
      cJ = cJ
    )
  )
  result
}


# =============================================================================
# Expected (Marginal) TV Bound
# =============================================================================

.a1_expected_tv_adaptive <- function(J, a, b, cJ, tail_probability,
                                     rel_tol, abs_tol) {
  lower_alpha <- stats::qgamma(tail_probability, shape = a, rate = b)
  upper_alpha <- stats::qgamma(
    tail_probability, shape = a, rate = b, lower.tail = FALSE
  )
  if (!is.finite(lower_alpha) || lower_alpha < 0 ||
      !is.finite(upper_alpha) || upper_alpha <= 0) {
    stop(.dpprior_new_condition(
      "adaptive marginal TV integration could not construct finite Gamma quantiles",
      classes = c(
        "dpprior_tv_bound_integration_error", "dpprior_numerical_error",
        "dpprior_error", "error"
      ),
      lower_alpha = lower_alpha,
      upper_alpha = upper_alpha,
      tail_probability = tail_probability,
      code = "invalid_integration_quantiles"
    ))
  }
  if (lower_alpha == 0) lower_alpha <- .Machine$double.xmin
  if (lower_alpha >= upper_alpha) {
    stop(.dpprior_new_condition(
      "adaptive marginal TV integration has an empty finite interval",
      classes = c(
        "dpprior_tv_bound_integration_error", "dpprior_numerical_error",
        "dpprior_error", "error"
      ),
      lower_alpha = lower_alpha,
      upper_alpha = upper_alpha,
      code = "empty_integration_interval"
    ))
  }

  log_integrand <- function(log_alpha) {
    alpha <- exp(log_alpha)
    log_measure <- stats::dgamma(
      alpha, shape = a, rate = b, log = TRUE
    ) + log_alpha
    values <- compute_total_tv_bound(J, alpha, cJ) * exp(log_measure)
    values[!is.finite(values) & log_measure == -Inf] <- 0
    values
  }
  integral <- tryCatch(
    stats::integrate(
      log_integrand,
      lower = log(lower_alpha), upper = log(upper_alpha),
      subdivisions = 1000L,
      rel.tol = rel_tol, abs.tol = abs_tol,
      stop.on.error = FALSE
    ),
    error = function(error) error
  )
  if (inherits(integral, "error") ||
      !is.finite(integral$value) || !is.finite(integral$abs.error) ||
      integral$abs.error < 0 ||
      !identical(integral$message, "OK")) {
    stop(.dpprior_new_condition(
      sprintf(
        "adaptive marginal TV integration failed: %s",
        if (inherits(integral, "error")) {
          conditionMessage(integral)
        } else {
          integral$message
        }
      ),
      classes = c(
        "dpprior_tv_bound_integration_error", "dpprior_numerical_error",
        "dpprior_error", "error"
      ),
      result = integral,
      code = "integration_failure"
    ))
  }

  omitted_probability <- stats::pgamma(
    lower_alpha, shape = a, rate = b
  ) + stats::pgamma(
    upper_alpha, shape = a, rate = b, lower.tail = FALSE
  )
  numerical_error_bound <- integral$abs.error + omitted_probability
  upper_bound <- .a1_clamp_tv_bound(
    min(1, integral$value + numerical_error_bound),
    "adaptive marginal expected TV bound"
  )

  list(
    estimate = integral$value,
    upper_bound = upper_bound,
    numerical_error_bound = numerical_error_bound,
    integration_error = integral$abs.error,
    omitted_probability = omitted_probability,
    lower_alpha = lower_alpha,
    upper_alpha = upper_alpha,
    message = integral$message,
    subdivisions = integral$subdivisions
  )
}


.a1_expected_tv_gl_audit <- function(J, a, b, cJ, controls) {
  evaluate <- function(order, label) {
    tryCatch(
      {
        value <- integrate_gamma(
          function(alpha) compute_total_tv_bound(J, alpha, cJ),
          a, b, order
        )
        .a1_clamp_tv_bound(
          value, sprintf("%s Gauss-Laguerre marginal TV audit", label)
        )
      },
      error = function(error) error
    )
  }

  selected <- evaluate(controls$M, "selected-order")
  verification <- if (is.null(controls$M_verify)) {
    NULL
  } else {
    evaluate(controls$M_verify, "verification-order")
  }
  errors <- Filter(
    function(value) inherits(value, "error"),
    list(selected = selected, verification = verification)
  )
  if (length(errors) > 0L) {
    return(list(
      status = "failed",
      verified = FALSE,
      reason = "quadrature_audit_failure",
      method = "Gauss-Laguerre secondary audit",
      M_selected = controls$M,
      M_verification = if (is.null(controls$M_verify)) {
        NA_integer_
      } else {
        controls$M_verify
      },
      selected = if (inherits(selected, "error")) NA_real_ else selected,
      verification = if (inherits(verification, "error")) {
        NA_real_
      } else {
        verification
      },
      absolute_difference = NA_real_,
      tolerance = NA_real_,
      message = paste(vapply(errors, conditionMessage, character(1)),
                      collapse = "; ")
    ))
  }

  if (is.null(verification)) {
    return(list(
      status = "approximate",
      verified = FALSE,
      reason = "fixed_order_unverified",
      method = "Gauss-Laguerre secondary audit",
      M_selected = controls$M,
      M_verification = NA_integer_,
      selected = selected,
      verification = NULL,
      absolute_difference = NA_real_,
      tolerance = NA_real_,
      message = "No independent Gauss-Laguerre audit order was requested"
    ))
  }

  difference <- abs(selected - verification)
  tolerance <- controls$abs_tol + controls$rel_tol *
    max(abs(selected), abs(verification))
  passed <- is.finite(difference) && difference <= tolerance
  list(
    status = if (passed) "converged" else "approximate",
    verified = passed,
    reason = if (passed) {
      "higher_order_agreement"
    } else {
      "higher_order_disagreement"
    },
    method = "Gauss-Laguerre secondary audit",
    M_selected = controls$M,
    M_verification = controls$M_verify,
    selected = selected,
    verification = verification,
    absolute_difference = difference,
    tolerance = tolerance,
    message = if (passed) {
      "Selected and verification Gauss-Laguerre audit orders agree"
    } else {
      "Gauss-Laguerre audit orders disagree; the adaptive result is retained"
    }
  )
}

#' Expected TV Bound Under Gamma Prior
#'
#' Integrates the conditional TV bound over \eqn{\alpha \sim \text{Gamma}(a, b)}
#' to obtain the marginal error bound.
#'
#' @param J Integer; sample size.
#' @param a,b Numeric; Gamma hyperparameters.
#' @param cJ Numeric; scaling constant (default: log(J)).
#' @param M Integer; number of nodes for a secondary Gauss-Laguerre audit.
#' @param M_verify Optional independent order at least
#'   \code{max(2*M, M+40)} and no greater than 512 for the secondary audit.
#' @param abs_tol,rel_tol Non-negative adaptive selected-versus-verification
#'   tolerances. They are also used and reported for the Gauss-Laguerre audit.
#' @param strict Logical; if \code{TRUE}, require successful agreement between
#'   two independently controlled adaptive integrations and otherwise raise a
#'   typed
#'   \code{dpprior_tv_bound_convergence_error}.
#'
#' @return Numeric; \eqn{E[d_{TV} \text{ bound} | a, b]}. The scalar carries
#'   a \code{"tv_bound_metadata"} attribute with status, numerical error and
#'   omitted-tail allowances, discrepancies, a Gauss-Laguerre audit, and
#'   theorem provenance. The returned conservative adaptive upper bound is
#'   never silently replaced by either verification value.
#'
#' @details
#' From the submitted manuscript Appendix D, Theorem
#' \code{thm:marginal-tv} (Equation D17), the TV error between the exact prior
#' predictive \eqn{p(S_J | a, b)} and the A1 Negative-Binomial proxy is bounded by:
#' \deqn{d_{TV}(P^{\text{exact}}, Q^{A1}) \le E_{\alpha \sim \Gamma(a,b)}[B_{\text{Pois}} + B_{\text{lin}}]}
#'
#' This follows from the mixture contraction property of TV distance. The
#' expectation is evaluated on log-alpha by adaptive quadrature over central
#' Gamma quantiles. The returned value adds both the integrator's reported
#' absolute error and the full omitted Gamma-tail probability to the numerical
#' estimate. A tighter independently controlled integration determines the
#' convergence status. Fixed-order Gauss-Laguerre values are retained in the
#' metadata as a reproducibility audit because their convergence can be
#' non-monotone for this capped, non-smooth integrand.
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
expected_tv_bound <- function(
    J, a, b, cJ = log(J), M = .QUAD_NODES_DEFAULT,
    M_verify = NULL, abs_tol = 1e-10, rel_tol = 1e-8,
    strict = FALSE) {
  assert_valid_J(J)
  a <- .dpprior_validate_scalar(
    a, "a", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  b <- .dpprior_validate_scalar(
    b, "b", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  cJ <- .a1_validate_cJ(cJ)
  controls <- .marginal_verification_controls(
    M, M_verify, abs_tol, rel_tol, strict
  )

  selected_adaptive <- .a1_expected_tv_adaptive(
    J, a, b, cJ,
    tail_probability = 1e-11,
    rel_tol = 1e-9,
    abs_tol = 1e-11
  )
  verification_adaptive <- .a1_expected_tv_adaptive(
    J, a, b, cJ,
    tail_probability = 1e-13,
    rel_tol = 1e-11,
    abs_tol = 1e-13
  )
  selected <- selected_adaptive$upper_bound
  verification <- verification_adaptive$upper_bound
  difference <- abs(selected - verification)
  tolerance <- controls$abs_tol + controls$rel_tol *
    max(abs(selected), abs(verification))
  passed <- is.finite(difference) && difference <= tolerance
  reason <- if (passed) {
    "adaptive_integration_agreement"
  } else {
    "adaptive_integration_disagreement"
  }
  quadrature_audit <- .a1_expected_tv_gl_audit(
    J, a, b, cJ, controls
  )

  status <- if (isTRUE(passed)) "converged" else "approximate"
  metadata <- list(
    schema_version = 1L,
    status = status,
    usable = identical(status, "converged"),
    verified = isTRUE(passed),
    reason = reason,
    method = paste(
      "adaptive log-alpha integration of capped conditional TV bound",
      "with numerical-error and omitted-tail allowance"
    ),
    theorem = "submitted manuscript Appendix D, Theorem marginal-tv (D17)",
    M_selected = controls$M,
    M_verification = if (is.null(controls$M_verify)) {
      NA_integer_
    } else {
      controls$M_verify
    },
    M_verification_required = controls$M_verification_required,
    verification_available = controls$verification_available,
    selected = selected,
    verification = verification,
    selected_estimate = selected_adaptive$estimate,
    verification_estimate = verification_adaptive$estimate,
    selected_numerical_error_bound =
      selected_adaptive$numerical_error_bound,
    verification_numerical_error_bound =
      verification_adaptive$numerical_error_bound,
    adaptive_selected = selected_adaptive,
    adaptive_verification = verification_adaptive,
    absolute_difference = difference,
    tolerance = tolerance,
    absolute_tolerance = controls$abs_tol,
    relative_tolerance = controls$rel_tol,
    quadrature_audit = quadrature_audit
  )

  if (controls$strict && !isTRUE(passed)) {
    stop(.dpprior_new_condition(
      sprintf(
        "marginal expected TV bound did not pass adaptive verification (%s)",
        reason
      ),
      classes = c(
        "dpprior_tv_bound_convergence_error", "dpprior_numerical_error",
        "dpprior_error", "error"
      ),
      result = metadata,
      code = reason
    ))
  }

  attr(selected, "tv_bound_metadata") <- metadata
  selected
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
#' @param M Selected quadrature order for every exact-moment evaluation.
#' @param M_verify Optional independent verification order.
#' @param abs_tol,rel_tol Non-negative quadrature comparison tolerances.
#' @param strict Logical; require every scanned exact-moment result to be
#'   verified.
#'
#' @return Integer; minimum scanned J achieving target accuracy, or NA if not
#'   found. A \code{"threshold_metadata"} attribute records every scan point
#'   and whether the conclusion is verified or provisional.
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
                                J_min = 10, J_max = 500, step = 10,
                                M = .QUAD_NODES_DEFAULT, M_verify = NULL,
                                abs_tol = 1e-10, rel_tol = 1e-8,
                                strict = FALSE) {
  a <- .dpprior_validate_scalar(
    a, "a", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  b <- .dpprior_validate_scalar(
    b, "b", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  target_error <- .dpprior_validate_scalar(
    target_error, "target_error", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_control_error"
  )

  if (is.null(target_var_error)) {
    target_var_error <- 2 * target_error
  }
  target_var_error <- .dpprior_validate_scalar(
    target_var_error, "target_var_error", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_control_error"
  )
  J_min <- .dpprior_validate_count(
    J_min, "J_min", minimum = 1L, maximum = .MAX_J_DEFAULT,
    .subclass = "dpprior_control_error"
  )
  J_max <- .dpprior_validate_count(
    J_max, "J_max", minimum = J_min, maximum = .MAX_J_DEFAULT,
    .subclass = "dpprior_control_error"
  )
  step <- .dpprior_validate_count(
    step, "step", minimum = 1L, .subclass = "dpprior_control_error"
  )

  scan <- vector("list", length(seq(J_min, J_max, by = step)))
  candidate <- NA_integer_
  evaluated <- 0L

  for (J in seq(J_min, J_max, by = step)) {
    evaluated <- evaluated + 1L
    errors <- a1_moment_error(
      J, a, b, M = M, M_verify = M_verify,
      abs_tol = abs_tol, rel_tol = rel_tol, strict = strict
    )
    error_metadata <- attr(errors, "a1_error_metadata", exact = TRUE)
    meets <- !is.na(errors$error_mean_rel) &&
      !is.na(errors$error_var_rel) &&
      errors$error_mean_rel / 100 < target_error &&
      errors$error_var_rel / 100 < target_var_error
    scan[[evaluated]] <- data.frame(
      J = as.integer(J),
      error_mean_rel = errors$error_mean_rel,
      error_var_rel = errors$error_var_rel,
      meets_criteria = meets,
      exact_moment_status = error_metadata$status,
      stringsAsFactors = FALSE
    )
    if (meets) {
      candidate <- as.integer(J)
      break
    }
  }

  scan <- do.call(rbind, scan[seq_len(evaluated)])
  verified <- nrow(scan) > 0L && all(scan$exact_moment_status == "converged")
  status <- if (verified) "converged" else "approximate"
  result <- candidate
  attr(result, "threshold_metadata") <- list(
    schema_version = 1L,
    status = status,
    usable = verified,
    verified = verified,
    message = if (verified) {
      "Every evaluated threshold-scan point passed exact-moment verification"
    } else {
      paste(
        "The threshold candidate is provisional because at least one",
        "evaluated exact-moment result is unverified"
      )
    },
    candidate = candidate,
    candidate_is_provisional = !verified,
    target_mean_relative_error = target_error,
    target_variance_relative_error = target_var_error,
    scan = scan,
    quadrature = list(
      M_selected = as.integer(M),
      M_verification = if (is.null(M_verify)) NA_integer_ else {
        as.integer(M_verify)
      },
      absolute_tolerance = abs_tol,
      relative_tolerance = rel_tol
    )
  )
  result
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
#' @param M_verify Optional independent quadrature order. When omitted and an
#'   admissible order exists, the package-wide minimum
#'   \code{max(2*M, M+40)} is used automatically.
#' @param abs_tol,rel_tol Non-negative selected-versus-verification tolerances.
#' @param strict Logical; require successful numerical verification of the
#'   current exact moments, every adequacy-threshold scan point, and the
#'   marginal TV bound.
#'
#' @return An S3 object of class "DPprior_error_bounds" with components:
#'   \describe{
#'     \item{J, a, b, cJ}{Input parameters}
#'     \item{moment_errors}{List of moment error metrics from \code{a1_moment_error}}
#'     \item{tv_bounds}{List with conditional and marginal TV bounds}
#'     \item{recommendation}{"A1_sufficient" or "A2_recommended"}
#'     \item{threshold_J}{First scan-grid J meeting the A1 adequacy criteria
#'       (or NA), with verified/provisional status in \code{verification}}
#'     \item{status, verified, usable}{Numerical verification state for this
#'       error-analysis computation. These fields do not assert that A1 itself
#'       is accurate; that is described by the reported errors and bounds.}
#'     \item{verification, provenance}{Numerical checks, including threshold
#'       scan provenance, and the theorem crosswalk.}
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
DPprior_error_bounds <- function(
    J, a, b, cJ = log(J), M = .QUAD_NODES_DEFAULT,
    M_verify = NULL, abs_tol = 1e-10, rel_tol = 1e-8,
    strict = FALSE) {
  assert_valid_J(J)
  a <- .dpprior_validate_scalar(
    a, "a", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  b <- .dpprior_validate_scalar(
    b, "b", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_gamma_parameter_error"
  )
  cJ <- .a1_validate_cJ(cJ)
  M <- .dpprior_validate_count(
    M, "M", minimum = 1L, maximum = .QUADRATURE_MAX_NODES,
    .subclass = "dpprior_control_error"
  )
  strict <- .dpprior_validate_control(strict, "strict", type = "logical")
  required_order <- .quadrature_verification_required_order(M)
  if (is.null(M_verify) && required_order <= .QUADRATURE_MAX_NODES) {
    M_verify <- required_order
  }

  # Moment errors
  moment_errors <- a1_moment_error(
    J, a, b, cJ = cJ, M = M, M_verify = M_verify,
    abs_tol = abs_tol, rel_tol = rel_tol, strict = strict
  )
  moment_metadata <- attr(moment_errors, "a1_error_metadata", exact = TRUE)

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
  marginal_bound <- expected_tv_bound(
    J, a, b, cJ = cJ, M = M, M_verify = M_verify,
    abs_tol = abs_tol, rel_tol = rel_tol, strict = strict
  )
  marginal_metadata <- attr(
    marginal_bound, "tv_bound_metadata", exact = TRUE
  )

  # Recommendation based on moment errors
  a1_sufficient <- !is.na(moment_errors$error_mean_rel) &&
    !is.na(moment_errors$error_var_rel) &&
    moment_errors$error_mean_rel < 5 &&
    moment_errors$error_var_rel < 10
  recommendation <- if (a1_sufficient) "A1_sufficient" else "A2_recommended"

  # Threshold J estimation
  threshold_J <- find_a1_threshold_J(
    a, b, target_error = 0.05,
    M = M, M_verify = M_verify,
    abs_tol = abs_tol, rel_tol = rel_tol, strict = strict
  )
  threshold_metadata <- attr(
    threshold_J, "threshold_metadata", exact = TRUE
  )
  attr(threshold_J, "threshold_metadata") <- NULL

  component_status <- c(
    exact_moments = moment_metadata$status,
    marginal_tv = marginal_metadata$status,
    adequacy_threshold = threshold_metadata$status
  )
  status <- if (all(component_status == "converged")) {
    "converged"
  } else {
    "approximate"
  }
  message <- if (identical(status, "converged")) {
    "All primary numerical error-analysis components passed verification"
  } else {
    paste(
      "Finite selected-order results are retained, but at least one",
      "primary numerical verification tolerance was not met"
    )
  }

  result <- list(
    status = status,
    usable = identical(status, "converged"),
    verified = identical(status, "converged"),
    message = message,
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
    threshold_J = threshold_J,
    verification = list(
      component_status = component_status,
      exact_moments = moment_metadata$exact_moment_quadrature,
      marginal_tv = marginal_metadata,
      adequacy_threshold = threshold_metadata
    ),
    provenance = list(
      schema_version = 1L,
      selected_method = "A1_error_decomposition",
      conditional_proxy = "Poisson(alpha * cJ)",
      marginal_proxy = "1 + NegativeBinomial(a, b/(b+cJ))",
      cJ = cJ,
      theorem_crosswalk = list(
        poissonization = paste(
          "Appendix D D10-D11; implementation uses the stronger",
          "(1-exp(-lambda))/lambda Chen-Stein prefactor"
        ),
        linearization = "Appendix D D14-D16 (Poisson KL plus Pinsker)",
        mixing = "Appendix D D17 (TV contraction under common Gamma mixing)",
        mixing_is_not_an_additional_approximation = TRUE
      )
    )
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

  cat(sprintf("Status:      %s\n", toupper(x$status)))
  cat(sprintf("Verified:    %s\n", if (isTRUE(x$verified)) "yes" else "no"))
  cat("Method used: A1 error decomposition; adaptive marginal-TV integration\n")
  marginal_metadata <- x$verification$marginal_tv
  gl_audit <- marginal_metadata$quadrature_audit
  if (!is.null(gl_audit)) {
    cat(sprintf("GL audit:    %s", toupper(gl_audit$status)))
    if (is.finite(gl_audit$absolute_difference)) {
      cat(sprintf(" (absolute difference %.3g)",
                  gl_audit$absolute_difference))
    }
    cat("\n")
  }
  if (!isTRUE(x$verified)) {
    cat(sprintf("Caveat:      %s\n", x$message))
  } else if (!is.null(gl_audit) &&
             !identical(gl_audit$status, "converged")) {
    cat(sprintf("Audit note:  %s\n", gl_audit$message))
  }
  cat("\n")

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
  threshold_metadata <- x$verification$adequacy_threshold
  if (isTRUE(threshold_metadata$verified) && !is.na(x$threshold_J)) {
    cat(sprintf("  First verified scan-grid J meeting criteria: %d\n",
                x$threshold_J))
  } else if (isTRUE(threshold_metadata$verified)) {
    cat("  No scanned J met the A1 adequacy criteria.\n")
  } else if (!is.na(x$threshold_J)) {
    cat(sprintf("  Provisional unverified scan-grid candidate: J = %d\n",
                x$threshold_J))
  } else {
    cat("  Threshold scan is unverified; no certified conclusion is available.\n")
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
    assert_valid_J(J)
    assert_positive(alpha, "alpha")
    .a1_shifted_mean(J, alpha)
  }, results$J, results$alpha)

  results$lambda_approx <- mapply(function(J, alpha) {
    cJ <- .a1_validate_cJ(cJ_fun(J))
    if (cJ > 0 && alpha > .Machine$double.xmax / cJ) Inf else alpha * cJ
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
