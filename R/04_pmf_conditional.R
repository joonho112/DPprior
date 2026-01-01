# =============================================================================
# Module 04: Conditional PMF of K_J | alpha (Antoniak Distribution)
# =============================================================================
#
# This module provides numerically stable computation of P(K_J = k | alpha)
# using the Antoniak distribution, which is derived from unsigned Stirling
# numbers of the first kind.
#
# Theory (RN-01 Theorem 2):
# ------------------------
# The exact PMF of K_J given alpha is:
#   P(K_J = k | alpha) = |s(J,k)| * alpha^k / (alpha)_J
#
# where:
#   - |s(J,k)| = unsigned Stirling number of first kind (count of permutations
#                of J elements with exactly k cycles)
#   - (alpha)_J = alpha * (alpha+1) * ... * (alpha+J-1) = Gamma(alpha+J)/Gamma(alpha)
#                 is the rising factorial (Pochhammer symbol)
#
# Log-Space Computation:
# ---------------------
# To avoid numerical overflow, we compute in log-space:
#   log P(K_J = k | alpha) = log|s(J,k)| + k*log(alpha) - log(alpha)_J
#                          = L_{J,k} + k*log(alpha) - [lgamma(alpha+J) - lgamma(alpha)]
#
# Dependencies:
# -------------
# - Module 00: logsumexp, softmax functions
# - Module 01: compute_log_stirling for log|s(J,k)|
# - Module 03: mean_K_given_alpha, var_K_given_alpha for verification
#
# Author: JoonHo Lee
# Date: December 2025
# Part of: DPprior R Package
# Reference: RN-01 Theorem 2, Antoniak (1974)
# =============================================================================


# =============================================================================
# Internal Validation Helpers
# =============================================================================

#' Validate Precomputed Log-Stirling Matrix Size
#'
#' Checks that the log-Stirling matrix is properly formatted and large enough
#' for the given sample size J.
#'
#' @param J Sample size (integer >= 1).
#' @param logS Pre-computed log-Stirling matrix from \code{compute_log_stirling()}.
#'
#' @return Invisible \code{TRUE} if validation passes.
#'
#' @keywords internal
.check_logS_size <- function(J, logS) {

  if (!is.matrix(logS) || nrow(logS) != ncol(logS)) {
    stop("logS must be a square matrix returned by compute_log_stirling()",
         call. = FALSE)
  }
  if (J + 1L > nrow(logS)) {
    stop(sprintf("Stirling matrix too small: J=%d but matrix has J_max=%d",
                 J, nrow(logS) - 1L), call. = FALSE)
  }
  invisible(TRUE)
}


# =============================================================================
# Rising Factorial (Pochhammer Symbol)
# =============================================================================

#' Log Rising Factorial (Pochhammer Symbol)
#'
#' Computes the logarithm of the rising factorial (Pochhammer symbol)
#' \eqn{(\alpha)_J = \alpha(\alpha+1)\cdots(\alpha+J-1)}.
#'
#' @param alpha Numeric; concentration parameter (must be positive scalar).
#' @param J Integer; number of terms (must be >= 1).
#'
#' @return Numeric; \eqn{\log(\alpha)_J}.
#'
#' @details
#' Uses the identity:
#' \deqn{(\alpha)_J = \frac{\Gamma(\alpha+J)}{\Gamma(\alpha)}}
#'
#' Therefore:
#' \deqn{\log(\alpha)_J = \log\Gamma(\alpha+J) - \log\Gamma(\alpha)}
#'
#' This is numerically stable for all \eqn{\alpha > 0} and \eqn{J \geq 1}.
#'
#' @examples
#' # (2)_3 = 2 * 3 * 4 = 24
#' exp(log_rising_factorial(2, 3))
#'
#' # (1)_5 = 5!
#' exp(log_rising_factorial(1, 5))
#' factorial(5)
#'
#' @seealso \code{\link{lgamma}} for the log-gamma function
#'
#' @keywords internal
#' @export
log_rising_factorial <- function(alpha, J) {
  assert_valid_J(J)
  assert_positive(alpha, "alpha")
  if (length(alpha) != 1L) {
    stop("alpha must be a scalar", call. = FALSE)
  }

  lgamma(alpha + J) - lgamma(alpha)
}


# =============================================================================
# Log-PMF Computation
# =============================================================================

#' Log-PMF of K Given Alpha (Antoniak Distribution)
#'
#' Computes \eqn{\log P(K_J = k \mid \alpha)} for \eqn{k = 0, 1, \ldots, J}.
#'
#' @param J Integer; sample size (number of observations, must be >= 1).
#' @param alpha Numeric; DP concentration parameter (must be positive scalar).
#' @param logS Matrix; pre-computed log-Stirling matrix from
#'   \code{\link{compute_log_stirling}}.
#'
#' @return Numeric vector of length \eqn{J+1} containing
#'   \eqn{\log P(K_J = k \mid \alpha)} for \eqn{k = 0, 1, \ldots, J}.
#'   Note that entry \code{[1]} corresponds to \eqn{k=0} and always equals
#'   \code{-Inf} (since \eqn{P(K_J = 0) = 0}).
#'
#' @details
#' Uses the Antoniak distribution formula in log-space:
#' \deqn{\log P(K_J = k \mid \alpha) = \log|s(J,k)| + k\log\alpha - \log(\alpha)_J}
#'
#' where \eqn{|s(J,k)|} is the unsigned Stirling number of the first kind and
#' \eqn{(\alpha)_J} is the rising factorial.
#'
#' This log-space computation is numerically stable for large \eqn{J} where
#' direct computation would overflow.
#'
#' @examples
#' logS <- compute_log_stirling(50)
#' log_pmf <- log_pmf_K_given_alpha(50, 2.0, logS)
#'
#' # Convert to probabilities
#' pmf <- softmax(log_pmf)
#' sum(pmf)  # Should be 1
#'
#' @seealso \code{\link{pmf_K_given_alpha}} for normalized PMF,
#'   \code{\link{compute_log_stirling}} for Stirling computation
#'
#' @keywords internal
#' @export
log_pmf_K_given_alpha <- function(J, alpha, logS) {
  # Input validation
  assert_valid_J(J)
  assert_positive(alpha, "alpha")
  if (length(alpha) != 1L) {
    stop("alpha must be a scalar", call. = FALSE)
  }
  .check_logS_size(J, logS)

  J <- as.integer(J)

  # Compute log of rising factorial: log(alpha)_J
  log_pochhammer <- log_rising_factorial(alpha, J)
  log_alpha <- log(alpha)

  # Initialize log-PMF with -Inf (P(K=k) = 0)
  log_pmf <- rep(-Inf, J + 1L)

  # Compute log P(K = k) for k = 1, ..., J
  for (k in seq_len(J)) {
    # Get log|s(J, k)| from pre-computed matrix
    # Note: R is 1-indexed, so logS[J+1, k+1] = log|s(J, k)|
    log_stirling_k <- logS[J + 1L, k + 1L]

    if (is.finite(log_stirling_k)) {
      log_pmf[k + 1L] <- log_stirling_k + k * log_alpha - log_pochhammer
    }
    # If log_stirling_k is -Inf, log_pmf[k+1] stays -Inf
  }

  log_pmf
}


# =============================================================================
# Normalized PMF
# =============================================================================

#' PMF of K Given Alpha (Antoniak Distribution)
#'
#' Computes \eqn{P(K_J = k \mid \alpha)} for \eqn{k = 0, 1, \ldots, J} using
#' the Antoniak distribution derived from the Dirichlet process.
#'
#' @param J Integer; sample size (number of observations, must be >= 1).
#' @param alpha Numeric; DP concentration parameter (must be positive scalar).
#' @param logS Matrix; pre-computed log-Stirling matrix from
#'   \code{\link{compute_log_stirling}}.
#' @param normalize Logical; if \code{TRUE} (default), use softmax normalization
#'   for numerical stability. If \code{FALSE}, return raw exponentiated values.
#'
#' @return Numeric vector of length \eqn{J+1} containing
#'   \eqn{P(K_J = k \mid \alpha)} for \eqn{k = 0, 1, \ldots, J}.
#'   Entry \code{[1]} corresponds to \eqn{k=0} and always equals 0.
#'   The vector sums to 1 (when \code{normalize = TRUE}).
#'
#' @details
#' The Antoniak distribution gives the exact PMF of the number of occupied
#' clusters \eqn{K_J} in a Dirichlet process with concentration parameter
#' \eqn{\alpha}:
#' \deqn{P(K_J = k \mid \alpha) = |s(J,k)| \frac{\alpha^k}{(\alpha)_J}}
#'
#' where \eqn{|s(J,k)|} is the unsigned Stirling number of the first kind
#' and \eqn{(\alpha)_J} is the rising factorial.
#'
#' Key properties:
#' \itemize{
#'   \item \eqn{P(K_J = 0) = 0} always (at least one cluster exists)
#'   \item \eqn{P(K_J = J) > 0} for all \eqn{\alpha > 0}
#'   \item As \eqn{\alpha \to 0^+}, mass concentrates on \eqn{K_J = 1}
#'   \item As \eqn{\alpha \to \infty}, mass concentrates on \eqn{K_J = J}
#' }
#'
#' Even though the Antoniak formula is theoretically normalized, converting
#' log-probabilities to the probability scale via \code{exp()} can underflow
#' for large \code{J} or extreme \code{alpha}. Setting \code{normalize=TRUE}
#' mitigates this by applying \code{softmax()} to the log-PMF.
#'
#' @examples
#' # Compute PMF for J=50, alpha=2
#' logS <- compute_log_stirling(50)
#' pmf <- pmf_K_given_alpha(50, 2.0, logS)
#'
#' # Verify normalization
#' sum(pmf)  # Should be 1
#'
#' # Verify P(K=0) = 0
#' pmf[1]    # Should be 0
#'
#' # Most likely number of clusters (mode)
#' which.max(pmf) - 1  # Subtract 1 for 0-indexing
#'
#' # Compare with moments
#' k_vals <- 0:50
#' mean_K <- sum(k_vals * pmf)
#' var_K <- sum(k_vals^2 * pmf) - mean_K^2
#'
#' # These should match digamma formulas
#' mean_K_given_alpha(50, 2.0)
#' var_K_given_alpha(50, 2.0)
#'
#' @seealso \code{\link{log_pmf_K_given_alpha}} for log-scale computation,
#'   \code{\link{mode_K_given_alpha}}, \code{\link{cdf_K_given_alpha}},
#'   \code{\link{quantile_K_given_alpha}}
#'
#' @references
#' Antoniak, C. E. (1974). Mixtures of Dirichlet Processes with Applications
#' to Bayesian Nonparametric Problems. \emph{The Annals of Statistics},
#' 2(6), 1152-1174.
#'
#' @export
pmf_K_given_alpha <- function(J, alpha, logS, normalize = TRUE) {
  log_pmf <- log_pmf_K_given_alpha(J, alpha, logS)

  if (isTRUE(normalize)) {
    softmax(log_pmf)
  } else {
    exp(log_pmf)
  }
}


# =============================================================================
# Distribution Summary Functions
# =============================================================================

#' Mode of K Given Alpha
#'
#' Computes the mode (most likely value) of \eqn{K_J \mid \alpha}.
#'
#' @param J Integer; sample size.
#' @param alpha Numeric; DP concentration parameter.
#' @param logS Matrix; pre-computed log-Stirling matrix.
#'
#' @return Integer; the value \eqn{k} that maximizes \eqn{P(K_J = k \mid \alpha)}.
#'
#' @examples
#' logS <- compute_log_stirling(50)
#' mode_K_given_alpha(50, 2.0, logS)
#'
#' @keywords internal
#' @export
mode_K_given_alpha <- function(J, alpha, logS) {
  pmf <- pmf_K_given_alpha(J, alpha, logS, normalize = TRUE)
  as.integer(which.max(pmf) - 1L)
}


#' CDF of K Given Alpha
#'
#' Computes the cumulative distribution function \eqn{P(K_J \leq k \mid \alpha)}
#' for \eqn{k = 0, 1, \ldots, J}.
#'
#' @param J Integer; sample size.
#' @param alpha Numeric; DP concentration parameter.
#' @param logS Matrix; pre-computed log-Stirling matrix.
#'
#' @return Numeric vector of length \eqn{J+1} containing
#'   \eqn{P(K_J \leq k \mid \alpha)} for \eqn{k = 0, 1, \ldots, J}.
#'
#' @details
#' The CDF satisfies:
#' \itemize{
#'   \item \eqn{F(0) = 0} (since \eqn{P(K_J = 0) = 0})
#'   \item \eqn{F(J) = 1}
#'   \item \eqn{F(k)} is non-decreasing in \eqn{k}
#' }
#'
#' @examples
#' logS <- compute_log_stirling(50)
#' cdf <- cdf_K_given_alpha(50, 2.0, logS)
#'
#' # Verify CDF ends at 1
#' cdf[51]  # Should be 1
#'
#' # P(K <= 5)
#' cdf[6]
#'
#' @keywords internal
#' @export
cdf_K_given_alpha <- function(J, alpha, logS) {
  pmf <- pmf_K_given_alpha(J, alpha, logS, normalize = TRUE)
  cumsum(pmf)
}


#' Quantile of K Given Alpha
#'
#' Computes the \eqn{p}-th quantile of \eqn{K_J \mid \alpha}, defined as
#' the smallest \eqn{k} such that \eqn{P(K_J \leq k \mid \alpha) \geq p}.
#'
#' @param p Numeric; probability level(s) in \eqn{[0, 1]}. Can be scalar or vector.
#' @param J Integer; sample size.
#' @param alpha Numeric; DP concentration parameter.
#' @param logS Matrix; pre-computed log-Stirling matrix.
#'
#' @return Integer (or integer vector); the \eqn{p}-th quantile(s) of
#'   \eqn{K_J \mid \alpha}.
#'
#' @details
#' For \eqn{p = 0.5}, this gives the median. Note that for discrete
#' distributions, the quantile is defined as the smallest value where
#' the CDF meets or exceeds \eqn{p}.
#'
#' Edge cases:
#' \itemize{
#'   \item \eqn{p = 0}: returns 0 (the first k where CDF >= 0)
#'   \item \eqn{p = 1}: returns J (the maximum possible K)
#' }
#'
#' @examples
#' logS <- compute_log_stirling(50)
#'
#' # Single quantile (median)
#' quantile_K_given_alpha(0.5, 50, 2.0, logS)
#'
#' # Multiple quantiles
#' quantile_K_given_alpha(c(0.25, 0.5, 0.75), 50, 2.0, logS)
#'
#' # Edge cases
#' quantile_K_given_alpha(0, 50, 2.0, logS)  # Returns 0
#' quantile_K_given_alpha(1, 50, 2.0, logS)  # Returns 50
#'
#' @keywords internal
#' @export
quantile_K_given_alpha <- function(p, J, alpha, logS) {

  assert_probability(p, "p")

  J <- as.integer(J)
  cdf <- cdf_K_given_alpha(J, alpha, logS)

  # Internal function to find quantile for a single p
  # Uses tolerance for floating-point comparison at boundaries
  qfun <- function(pp) {
    # Special case: p = 1 (or very close) always returns J
    # This avoids floating-point issues where CDF[k] ≈ 0.9999999999 >= 1.0
    if (pp >= 1 - 1e-12) {
      return(J)
    }

    # Special case: p = 0 (or very close) returns 0
    if (pp <= 1e-12) {
      return(0L)
    }

    # Find smallest k where CDF >= p
    idx <- which(cdf >= pp)
    if (length(idx) == 0L) {
      # Fallback: return J (should not happen with valid inputs)
      return(J)
    }
    as.integer(min(idx) - 1L)
  }

  # Handle both scalar and vector inputs
  if (length(p) == 1L) {
    qfun(p)
  } else {
    vapply(p, qfun, integer(1))
  }
}


# =============================================================================
# Moments from PMF (for verification)
# =============================================================================

#' Mean of K from PMF
#'
#' Computes \eqn{E[K_J \mid \alpha]} by summing over the PMF.
#' This is primarily for verification against the closed-form digamma formula.
#'
#' @param J Integer; sample size.
#' @param alpha Numeric; DP concentration parameter.
#' @param logS Matrix; pre-computed log-Stirling matrix.
#'
#' @return Numeric; the conditional mean \eqn{E[K_J \mid \alpha]}.
#'
#' @examples
#' logS <- compute_log_stirling(50)
#'
#' # Should match mean_K_given_alpha(50, 2.0)
#' mean_K_from_pmf(50, 2.0, logS)
#' mean_K_given_alpha(50, 2.0)
#'
#' @keywords internal
#' @export
mean_K_from_pmf <- function(J, alpha, logS) {
  pmf <- pmf_K_given_alpha(J, alpha, logS, normalize = TRUE)
  k_vals <- 0:J
  sum(k_vals * pmf)
}


#' Variance of K from PMF
#'
#' Computes \eqn{Var(K_J \mid \alpha)} by summing over the PMF.
#' This is primarily for verification against the closed-form trigamma formula.
#'
#' @param J Integer; sample size.
#' @param alpha Numeric; DP concentration parameter.
#' @param logS Matrix; pre-computed log-Stirling matrix.
#'
#' @return Numeric; the conditional variance \eqn{Var(K_J \mid \alpha)}.
#'
#' @examples
#' logS <- compute_log_stirling(50)
#'
#' # Should match var_K_given_alpha(50, 2.0)
#' var_K_from_pmf(50, 2.0, logS)
#' var_K_given_alpha(50, 2.0)
#'
#' @keywords internal
#' @export
var_K_from_pmf <- function(J, alpha, logS) {
  pmf <- pmf_K_given_alpha(J, alpha, logS, normalize = TRUE)
  k_vals <- 0:J
  mean_k <- sum(k_vals * pmf)
  sum(k_vals^2 * pmf) - mean_k^2
}


# =============================================================================
# Comprehensive Summary
# =============================================================================

#' Summary of Conditional PMF
#'
#' Returns a comprehensive summary of the conditional distribution
#' \eqn{K_J \mid \alpha}.
#'
#' @param J Integer; sample size.
#' @param alpha Numeric; DP concentration parameter.
#' @param logS Matrix; pre-computed log-Stirling matrix.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{J}}{Sample size}
#'     \item{\code{alpha}}{Concentration parameter}
#'     \item{\code{mean}}{Conditional mean \eqn{E[K_J \mid \alpha]}}
#'     \item{\code{var}}{Conditional variance \eqn{Var(K_J \mid \alpha)}}
#'     \item{\code{sd}}{Conditional standard deviation}
#'     \item{\code{mode}}{Most likely value of K}
#'     \item{\code{median}}{Median of K}
#'     \item{\code{quantiles}}{25th, 50th, and 75th percentiles}
#'     \item{\code{pmf}}{Full PMF vector}
#'     \item{\code{cdf}}{Full CDF vector}
#'   }
#'
#' @examples
#' logS <- compute_log_stirling(50)
#' summary <- summary_pmf_K_given_alpha(50, 2.0, logS)
#' print(summary)
#'
#' @export
summary_pmf_K_given_alpha <- function(J, alpha, logS) {
  pmf <- pmf_K_given_alpha(J, alpha, logS, normalize = TRUE)
  cdf <- cumsum(pmf)
  k_vals <- 0:J

  mean_k <- sum(k_vals * pmf)
  var_k <- sum(k_vals^2 * pmf) - mean_k^2

  list(
    J = J,
    alpha = alpha,
    mean = mean_k,
    var = var_k,
    sd = sqrt(var_k),
    mode = as.integer(which.max(pmf) - 1L),
    median = quantile_K_given_alpha(0.5, J, alpha, logS),
    quantiles = c(
      q25 = quantile_K_given_alpha(0.25, J, alpha, logS),
      q50 = quantile_K_given_alpha(0.50, J, alpha, logS),
      q75 = quantile_K_given_alpha(0.75, J, alpha, logS)
    ),
    pmf = pmf,
    cdf = cdf
  )
}


# =============================================================================
# Verification Functions
# =============================================================================

#' Verify PMF Consistency with Moments
#'
#' Verifies that moments computed from the PMF match the closed-form
#' digamma/trigamma formulas.
#'
#' @param J Integer; sample size.
#' @param alpha Numeric; DP concentration parameter.
#' @param logS Matrix; pre-computed log-Stirling matrix.
#' @param tol Numeric; tolerance for comparison (default: 1e-8).
#' @param verbose Logical; if \code{TRUE}, print results.
#'
#' @return Logical; \code{TRUE} if verification passes.
#'
#' @examples
#' logS <- compute_log_stirling(50)
#' verify_pmf_moments(50, 2.0, logS)
#'
#' @export
verify_pmf_moments <- function(J, alpha, logS, tol = 1e-8, verbose = TRUE) {
  # Moments from PMF
  mean_pmf <- mean_K_from_pmf(J, alpha, logS)
  var_pmf <- var_K_from_pmf(J, alpha, logS)

  # Moments from closed-form (Module 03)
  mean_cf <- mean_K_given_alpha(J, alpha)
  var_cf <- var_K_given_alpha(J, alpha)

  # Errors
  err_mean <- abs(mean_pmf - mean_cf)
  err_var <- abs(var_pmf - var_cf)

  pass_mean <- err_mean < tol
  pass_var <- err_var < tol
  all_pass <- pass_mean && pass_var

  if (isTRUE(verbose)) {
    cat(sprintf("PMF-Moments Verification (J=%d, alpha=%.2f):\n", J, alpha))
    cat(sprintf("  Mean: PMF=%.8f, digamma=%.8f, error=%.2e [%s]\n",
                mean_pmf, mean_cf, err_mean,
                if (pass_mean) "PASS" else "FAIL"))
    cat(sprintf("  Var:  PMF=%.8f, trigamma=%.8f, error=%.2e [%s]\n",
                var_pmf, var_cf, err_var,
                if (pass_var) "PASS" else "FAIL"))
    cat(sprintf("  Overall: %s\n", if (all_pass) "PASS" else "FAIL"))
  }

  invisible(all_pass)
}


#' Verify PMF Normalization
#'
#' Verifies that the PMF sums to 1.
#'
#' @param J Integer; sample size.
#' @param alpha Numeric; DP concentration parameter.
#' @param logS Matrix; pre-computed log-Stirling matrix.
#' @param tol Numeric; tolerance (default: 1e-10).
#' @param verbose Logical; if \code{TRUE}, print results.
#'
#' @return Logical; \code{TRUE} if verification passes.
#'
#' @export
verify_pmf_normalization <- function(J, alpha, logS, tol = 1e-10, verbose = TRUE) {
  pmf <- pmf_K_given_alpha(J, alpha, logS, normalize = TRUE)
  total <- sum(pmf)

  passed <- abs(total - 1.0) < tol

  if (isTRUE(verbose)) {
    cat(sprintf("PMF Normalization (J=%d, alpha=%.2f):\n", J, alpha))
    cat(sprintf("  Sum = %.12f (expected 1.0) [%s]\n",
                total, if (passed) "PASS" else "FAIL"))
  }

  invisible(passed)
}


#' Verify Zero Probability at K=0
#'
#' Verifies that P(K_J = 0 | alpha) = 0.
#'
#' @param J Integer; sample size.
#' @param alpha Numeric; DP concentration parameter.
#' @param logS Matrix; pre-computed log-Stirling matrix.
#' @param tol Numeric; tolerance (default: 1e-15).
#' @param verbose Logical; if \code{TRUE}, print results.
#'
#' @return Logical; \code{TRUE} if verification passes.
#'
#' @export
verify_zero_probability <- function(J, alpha, logS, tol = 1e-15, verbose = TRUE) {
  pmf <- pmf_K_given_alpha(J, alpha, logS, normalize = TRUE)

  passed <- pmf[1] < tol

  if (isTRUE(verbose)) {
    cat(sprintf("Zero Probability (J=%d, alpha=%.2f):\n", J, alpha))
    cat(sprintf("  P(K=0) = %.2e (expected 0) [%s]\n",
                pmf[1], if (passed) "PASS" else "FAIL"))
  }

  invisible(passed)
}


#' Verify CDF Properties
#'
#' Verifies that the CDF is non-decreasing and ends at 1.
#'
#' @param J Integer; sample size.
#' @param alpha Numeric; DP concentration parameter.
#' @param logS Matrix; pre-computed log-Stirling matrix.
#' @param tol Numeric; tolerance (default: 1e-10).
#' @param verbose Logical; if \code{TRUE}, print results.
#'
#' @return Logical; \code{TRUE} if verification passes.
#'
#' @export
verify_cdf_properties <- function(J, alpha, logS, tol = 1e-10, verbose = TRUE) {
  cdf <- cdf_K_given_alpha(J, alpha, logS)

  # Check monotonicity (allow small negative due to numerical errors)
  diffs <- diff(cdf)
  monotone <- all(diffs >= -tol)

  # Check final value
  final_ok <- abs(cdf[length(cdf)] - 1.0) < tol

  all_pass <- monotone && final_ok

  if (isTRUE(verbose)) {
    cat(sprintf("CDF Properties (J=%d, alpha=%.2f):\n", J, alpha))
    cat(sprintf("  Monotonicity: min(diff) = %.2e [%s]\n",
                min(diffs), if (monotone) "PASS" else "FAIL"))
    cat(sprintf("  CDF[J] = %.12f (expected 1.0) [%s]\n",
                cdf[length(cdf)], if (final_ok) "PASS" else "FAIL"))
    cat(sprintf("  Overall: %s\n", if (all_pass) "PASS" else "FAIL"))
  }

  invisible(all_pass)
}


#' Run All PMF Verifications
#'
#' Runs comprehensive verification tests for the conditional PMF module.
#'
#' @param J_values Integer vector; sample sizes to test.
#' @param alpha_values Numeric vector; concentration parameters to test.
#' @param verbose Logical; if \code{TRUE}, print detailed results.
#'
#' @return Logical; \code{TRUE} if all verifications pass.
#'
#' @examples
#' verify_pmf_all()
#'
#' @export
verify_pmf_all <- function(J_values = c(10, 50, 100),
                           alpha_values = c(0.5, 1.0, 2.0, 5.0),
                           verbose = TRUE) {
  # Pre-compute Stirling numbers
  J_max <- max(J_values)
  logS <- compute_log_stirling(J_max)

  all_pass <- TRUE

  for (J in J_values) {
    for (alpha in alpha_values) {
      if (isTRUE(verbose)) {
        cat(sprintf("\n%s\n", strrep("=", 60)))
        cat(sprintf("Test case: J=%d, alpha=%.2f\n", J, alpha))
        cat(sprintf("%s\n", strrep("=", 60)))
      }

      pass <- TRUE
      pass <- pass && verify_pmf_normalization(J, alpha, logS, verbose = verbose)
      pass <- pass && verify_zero_probability(J, alpha, logS, verbose = verbose)
      pass <- pass && verify_pmf_moments(J, alpha, logS, verbose = verbose)
      pass <- pass && verify_cdf_properties(J, alpha, logS, verbose = verbose)

      all_pass <- all_pass && pass
    }
  }

  if (isTRUE(verbose)) {
    cat(sprintf("\n%s\n", strrep("=", 60)))
    cat(sprintf("Overall: %s\n", if (all_pass) "ALL TESTS PASSED" else "SOME TESTS FAILED"))
    cat(sprintf("%s\n", strrep("=", 60)))
  }

  invisible(all_pass)
}
