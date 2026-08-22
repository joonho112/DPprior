# =============================================================================
# Module 03: Conditional Moments of K_J | alpha
# =============================================================================
#
# This module provides numerically stable computation of the conditional mean
# and variance of K_J given alpha, where K_J is the number of occupied clusters
# induced by a Dirichlet process (DP) prior with concentration parameter alpha.
#
# Theory (Lee, 2026, Sections 2--3):
#   mu_J(alpha) = E[K_J | alpha] = alpha * (psi(alpha + J) - psi(alpha))
#   v_J(alpha)  = Var(K_J | alpha) =
#                   mu_J(alpha) - alpha^2 * (psi1(alpha) - psi1(alpha + J))
#
# where psi(·) is the digamma function and psi1(·) is the trigamma function.
#
# Numerical notes:
# - Differences of digamma/trigamma values lose accuracy at large alpha.
# - Hard-coded small-alpha limits create discontinuities and give incorrect
#   derivatives. All four quantities are therefore evaluated as finite sums of
#   stable positive/ratio terms over r = 1, ..., J - 1.
# - The supported package domain J <= 500 makes these O(J) sums inexpensive.
#
# Author: JoonHo Lee (jlee296@ua.edu)
# Date: December 2025
# Part of: DPprior R Package
# Reference: Lee (2026), Sections 2--3
# =============================================================================


# =============================================================================
# Stable finite-sum primitives
# =============================================================================

# Evaluate one conditional-moment quantity without subtracting nearly equal
# special-function values. The ratio forms keep intermediate values in [0, 1]
# over the tested alpha domain and avoid squaring/cubing a large alpha.
.conditional_K_scalar <- function(J, alpha,
                                  quantity = c("mean", "var", "dmean", "dvar")) {
  quantity <- match.arg(quantity)

  if (J == 1L) {
    return(switch(quantity, mean = 1.0, var = 0.0,
                  dmean = 0.0, dvar = 0.0))
  }

  r <- seq_len(J - 1L)

  if (quantity == "mean") {
    # Sum the small side of each complement. The two branches are algebraically
    # identical: 1 + sum(alpha/(alpha+r)) = J - sum(r/(alpha+r)).
    if (alpha <= J) {
      return(1.0 + sum(alpha / (alpha + r)))
    }
    return(J - sum(r / (alpha + r)))
  }

  alpha_is_smaller <- alpha <= r
  ratio <- ifelse(alpha_is_smaller, alpha / r, r / alpha)

  if (quantity == "var") {
    # alpha*r/(alpha+r)^2 = x/(1+x)^2 for x=min/max.
    return(sum(ratio / (1.0 + ratio)^2))
  }

  if (quantity == "dmean") {
    # r/(alpha+r)^2, expressed without a large squared denominator.
    terms <- ifelse(
      alpha_is_smaller,
      1.0 / (r * (1.0 + ratio)^2),
      ratio / (alpha * (1.0 + ratio)^2)
    )
    return(sum(terms))
  }

  # r*(r-alpha)/(alpha+r)^3. This tends to 1/r as alpha -> 0+
  # and may be positive or negative away from the limit.
  terms <- ifelse(
    alpha_is_smaller,
    (1.0 - ratio) / (r * (1.0 + ratio)^3),
    ratio * (ratio - 1.0) / (alpha * (1.0 + ratio)^3)
  )
  sum(terms)
}

.conditional_K_vector <- function(J, alpha, quantity) {
  vapply(alpha, function(a) .conditional_K_scalar(J, a, quantity), numeric(1))
}


# =============================================================================
# Core Exported Functions
# =============================================================================

#' Conditional Mean of K_J Given Alpha
#'
#' Computes \eqn{\mathbb{E}[K_J \mid \alpha]} under a Dirichlet process prior.
#'
#' @param J Sample size (integer >= 1).
#' @param alpha Concentration parameter (positive numeric, vectorized).
#'
#' @return Numeric vector of conditional means (same length as \code{alpha}).
#'
#' @details
#' The classical digamma closed form is:
#' \deqn{\mu_J(\alpha) = \alpha\{\psi(\alpha+J)-\psi(\alpha)\}}
#' where \eqn{\psi(\cdot)} is the digamma function.
#'
#' This is equivalent to the direct summation:
#' \deqn{\mu_J(\alpha) = \sum_{i=1}^{J} \frac{\alpha}{\alpha + i - 1}}
#'
#' \strong{Limiting behavior:}
#' \itemize{
#'   \item \eqn{\alpha \to 0^+}: \eqn{\mu_J(\alpha) \to 1}
#'   \item \eqn{\alpha \to \infty}: \eqn{\mu_J(\alpha) \to J}
#' }
#'
#' The implementation evaluates the equivalent finite sum. For small
#' \code{alpha} it adds the small increments to 1; for large \code{alpha} it
#' subtracts the small deficits from \code{J}. This avoids special-function
#' cancellation and a discontinuous small-alpha guard.
#'
#' @examples
#' mean_K_given_alpha(50, 2.0)
#' mean_K_given_alpha(50, c(0.5, 1, 2, 5))
#'
#' # Limiting behavior
#' mean_K_given_alpha(50, 1e-10)  # Returns a value just above 1
#' mean_K_given_alpha(50, 1e6)    # Returns ~50
#'
#' @seealso \code{\link{var_K_given_alpha}}, \code{\link{moments_K_given_alpha}}
#'
#' @family conditional_K
#'
#' @export
mean_K_given_alpha <- function(J, alpha) {
  assert_valid_J(J)
  assert_positive(alpha, "alpha")

  alpha <- as.numeric(alpha)
  out <- .conditional_K_vector(J, alpha, "mean")

  out
}


#' Conditional Variance of K_J Given Alpha
#'
#' Computes \eqn{\mathrm{Var}(K_J \mid \alpha)} under a Dirichlet process prior.
#'
#' @param J Sample size (integer >= 1).
#' @param alpha Concentration parameter (positive numeric, vectorized).
#'
#' @return Numeric vector of conditional variances (same length as \code{alpha}).
#'
#' @details
#' The classical trigamma closed form is:
#' \deqn{v_J(\alpha) = \mu_J(\alpha) - \alpha^2\{\psi_1(\alpha) - \psi_1(\alpha+J)\}}
#' where \eqn{\psi_1(\cdot)} is the trigamma function.
#'
#' This is equivalent to the direct summation:
#' \deqn{v_J(\alpha) = \sum_{i=1}^{J} \frac{\alpha(i-1)}{(\alpha + i - 1)^2}}
#'
#' \strong{Key property:} \eqn{0 \leq v_J(\alpha) < \mu_J(\alpha)} for all
#' \eqn{\alpha > 0}. The lower inequality is strict when \code{J > 1}; at
#' \code{J = 1}, \eqn{K_J=1} deterministically and the variance is zero.
#'
#' \strong{Limiting behavior:}
#' \itemize{
#'   \item \eqn{\alpha \to 0^+}: \eqn{v_J(\alpha) \to 0}
#'   \item \eqn{\alpha \to \infty}: \eqn{v_J(\alpha) \to 0}
#' }
#'
#' The implementation evaluates each Bernoulli variance contribution as a
#' non-negative ratio term. It remains stable at both small and large
#' \code{alpha} without clipping cancellation artifacts to zero.
#'
#' @examples
#' var_K_given_alpha(50, 2.0)
#' var_K_given_alpha(50, c(0.5, 1, 2, 5))
#'
#' # Verify underdispersion
#' J <- 50; alpha <- 2.0
#' mean_K_given_alpha(J, alpha) > var_K_given_alpha(J, alpha)  # TRUE
#'
#' @seealso \code{\link{mean_K_given_alpha}}, \code{\link{moments_K_given_alpha}}
#'
#' @family conditional_K
#'
#' @export
var_K_given_alpha <- function(J, alpha) {
  assert_valid_J(J)
  assert_positive(alpha, "alpha")

  alpha <- as.numeric(alpha)
  out <- .conditional_K_vector(J, alpha, "var")

  out
}


# =============================================================================
# Internal Functions (for Jacobian support)
# =============================================================================

#' Derivative of Conditional Mean w.r.t. Alpha
#'
#' Computes \eqn{\frac{d}{d\alpha} \mathbb{E}[K_J \mid \alpha]}.
#' This is used internally when building Jacobians for Newton-type solvers.
#'
#' @param J Sample size (integer >= 1).
#' @param alpha Concentration parameter (positive numeric, vectorized).
#'
#' @return Numeric vector of derivatives (same length as \code{alpha}).
#'
#' @details
#' A special-function closed form is:
#' \deqn{\frac{d}{d\alpha}\mu_J(\alpha) =
#'       \{\psi(\alpha+J)-\psi(\alpha)\} + \alpha\{\psi_1(\alpha+J)-\psi_1(\alpha)\}}
#'
#' This derivative is positive for \eqn{\alpha > 0} when \code{J > 1},
#' confirming that \eqn{E[K_J | \alpha]} is then strictly increasing in
#' \eqn{\alpha}. At \code{J = 1}, the mean is identically one and the
#' derivative is zero.
#' The implementation uses the stable finite sum
#' \deqn{\sum_{r=1}^{J-1} r/(\alpha+r)^2}. In particular, both the mean and
#' variance derivatives approach \eqn{H_{J-1}} as \eqn{\alpha\to0^+}; they do
#' not approach zero.
#'
#' @examples
#' \dontrun{
#' dmean_dalpha(50, 2.0)
#'
#' # Verify with finite difference
#' J <- 50; alpha <- 2.0; eps <- 1e-6
#' fd <- (mean_K_given_alpha(J, alpha + eps) -
#'        mean_K_given_alpha(J, alpha - eps)) / (2 * eps)
#' abs(fd - dmean_dalpha(J, alpha)) < 1e-5  # TRUE
#'
#' }
#' @keywords internal
dmean_dalpha <- function(J, alpha) {
  assert_valid_J(J)
  assert_positive(alpha, "alpha")

  alpha <- as.numeric(alpha)
  out <- .conditional_K_vector(J, alpha, "dmean")

  out
}


#' Derivative of Conditional Variance w.r.t. Alpha
#'
#' Computes \eqn{\frac{d}{d\alpha} \mathrm{Var}(K_J \mid \alpha)}.
#'
#' @param J Sample size (integer >= 1).
#' @param alpha Concentration parameter (positive numeric, vectorized).
#'
#' @return Numeric vector of derivatives (same length as \code{alpha}).
#'
#' @details
#' Uses the summation form:
#' \deqn{\frac{d}{d\alpha} v_J(\alpha) = \sum_{r=1}^{J-1} \frac{r(r - \alpha)}{(\alpha + r)^3}}
#'
#' This derivative can be positive or negative depending on \eqn{\alpha}.
#'
#' @keywords internal
dvar_dalpha <- function(J, alpha) {
  assert_valid_J(J)
  assert_positive(alpha, "alpha")

  alpha <- as.numeric(alpha)
  out <- .conditional_K_vector(J, alpha, "dvar")

  out
}


# =============================================================================
# Convenience Functions
# =============================================================================

#' Conditional Mean and Variance of K_J Given Alpha
#'
#' Convenience wrapper returning both conditional moments in one call.
#'
#' @param J Sample size (integer >= 1).
#' @param alpha Concentration parameter (positive numeric, vectorized).
#'
#' @return If \code{alpha} is scalar, a named numeric vector
#'   \code{c(mean = ..., var = ...)}. If \code{alpha} is a vector, a numeric
#'   matrix with two columns \code{mean} and \code{var} (one row per element of
#'   \code{alpha}).
#'
#' @examples
#' moments_K_given_alpha(50, 2.0)
#' moments_K_given_alpha(50, c(0.5, 1, 2, 5))
#'
#' @seealso \code{\link{mean_K_given_alpha}}, \code{\link{var_K_given_alpha}}
#'
#' @family conditional_K
#'
#' @export
moments_K_given_alpha <- function(J, alpha) {
  mu <- mean_K_given_alpha(J, alpha)
  v <- var_K_given_alpha(J, alpha)

  if (length(mu) == 1L) {
    return(c(mean = as.numeric(mu), var = as.numeric(v)))
  }

  out <- cbind(mean = as.numeric(mu), var = as.numeric(v))
  rownames(out) <- names(mu)
  out
}


#' Coefficient of Variation for K Given Alpha
#'
#' Computes the coefficient of variation (CV) of \eqn{K_J | \alpha},
#' defined as \eqn{CV = SD(K) / E[K]}.
#'
#' @param J Sample size (integer >= 1).
#' @param alpha Concentration parameter (positive numeric, vectorized).
#'
#' @return Numeric vector of coefficients of variation.
#'
#' @family conditional_K
#'
#' @export
cv_K_given_alpha <- function(J, alpha) {
  mu <- mean_K_given_alpha(J, alpha)
  v <- var_K_given_alpha(J, alpha)
  sqrt(v) / mu
}


#' Dispersion Index for K Given Alpha
#'
#' Computes the dispersion index (variance-to-mean ratio) of \eqn{K_J | \alpha}.
#' Always < 1 for this distribution (underdispersion).
#'
#' @param J Sample size (integer >= 1).
#' @param alpha Concentration parameter (positive numeric, vectorized).
#'
#' @return Numeric vector of dispersion indices.
#'
#' @family conditional_K
#'
#' @export
dispersion_K_given_alpha <- function(J, alpha) {
  var_K_given_alpha(J, alpha) / mean_K_given_alpha(J, alpha)
}


#' Summary of Conditional Moments
#'
#' Returns a comprehensive summary of the conditional distribution of K given alpha.
#'
#' @param J Sample size (integer >= 1).
#' @param alpha Concentration parameter (scalar).
#'
#' @return A list with components: J, alpha, mean, var, sd, cv, dispersion, dmean_dalpha.
#'
#' @examples
#' \dontrun{
#' summary_K_given_alpha(50, 2.0)
#'
#' }
#' @keywords internal
summary_K_given_alpha <- function(J, alpha) {
  if (length(alpha) != 1L) {
    stop("alpha must be a scalar for summary_K_given_alpha()", call. = FALSE)
  }

  mu <- mean_K_given_alpha(J, alpha)
  v <- var_K_given_alpha(J, alpha)

  list(
    J = J,
    alpha = alpha,
    mean = as.numeric(mu),
    var = as.numeric(v),
    sd = sqrt(as.numeric(v)),
    cv = sqrt(as.numeric(v)) / as.numeric(mu),
    dispersion = as.numeric(v) / as.numeric(mu),
    dmean_dalpha = as.numeric(dmean_dalpha(J, alpha))
  )
}


# =============================================================================
# Verification Functions (for development and testing)
# =============================================================================

#' Compute Conditional Mean via Direct Summation
#'
#' For verification purposes only.
#'
#' @param J Sample size.
#' @param alpha Concentration parameter (scalar).
#'
#' @return Numeric; conditional mean computed via summation.
#'
#' @keywords internal
mean_K_given_alpha_sum <- function(J, alpha) {
  assert_valid_J(J)
  assert_positive(alpha, "alpha")
  sum(alpha / (alpha + seq_len(J) - 1))
}


#' Compute Conditional Variance via Direct Summation
#'
#' For verification purposes only.
#'
#' @param J Sample size.
#' @param alpha Concentration parameter (scalar).
#'
#' @return Numeric; conditional variance computed via summation.
#'
#' @keywords internal
var_K_given_alpha_sum <- function(J, alpha) {
  assert_valid_J(J)
  assert_positive(alpha, "alpha")
  i_vals <- seq_len(J)
  sum(alpha * (i_vals - 1) / (alpha + i_vals - 1)^2)
}


#' Validate Conditional Moments Computation
#'
#' Validates polygamma formulas against direct summation.
#'
#' @param J Sample size.
#' @param alpha Concentration parameter.
#' @param tol Tolerance (default: 1e-10).
#' @param verbose Print results if TRUE.
#'
#' @return Logical; TRUE if validation passes.
#'
#' @keywords internal
validate_moments_conditional <- function(J, alpha, tol = 1e-10, verbose = TRUE) {
  mean_poly <- mean_K_given_alpha(J, alpha)
  mean_sum <- mean_K_given_alpha_sum(J, alpha)
  var_poly <- var_K_given_alpha(J, alpha)
  var_sum <- var_K_given_alpha_sum(J, alpha)

  err_mean <- abs(mean_poly - mean_sum)
  err_var <- abs(var_poly - var_sum)

  pass_mean <- err_mean < tol
  pass_var <- err_var < tol
  all_pass <- pass_mean && pass_var

  if (verbose) {
    cat(sprintf("Validation (J=%d, alpha=%.2f):\n", J, alpha))
    cat(sprintf("  Mean error: %.2e [%s]\n", err_mean, if (pass_mean) "PASS" else "FAIL"))
    cat(sprintf("  Var error:  %.2e [%s]\n", err_var, if (pass_var) "PASS" else "FAIL"))
  }

  invisible(all_pass)
}


#' Verify Underdispersion Inequality
#'
#' Verifies that \eqn{0 \leq Var(K_J | \alpha) < E[K_J | \alpha]} for all
#' specified values of \eqn{\alpha}. The variance is zero at \code{J = 1}.
#'
#' @param J Sample size (integer >= 1).
#' @param alpha_values Numeric vector of alpha values to test
#'   (default: c(0.1, 0.5, 1, 2, 5, 10)).
#' @param verbose Logical; if TRUE, print results.
#'
#' @return Logical; TRUE if inequality holds for all alpha values.
#'
#' @examples
#' verify_underdispersion(50)
#' verify_underdispersion(50, c(0.1, 1, 10), verbose = TRUE)
#'
#' @export
verify_underdispersion <- function(J, alpha_values = c(0.1, 0.5, 1, 2, 5, 10),
                                   verbose = TRUE) {
  all_pass <- TRUE

  if (verbose) {
    cat(sprintf("Underdispersion verification (J=%d):\n", J))
  }

  for (alpha in alpha_values) {
    mu <- mean_K_given_alpha(J, alpha)
    v <- var_K_given_alpha(J, alpha)
    ok <- (v >= 0) && (v < mu)
    all_pass <- all_pass && ok

    if (verbose) {
      status <- if (ok) "PASS" else "FAIL"
      cat(sprintf("  alpha=%5.2f: E[K]=%8.4f, Var(K)=%8.4f, D=%.4f [%s]\n",
                  alpha, mu, v, v / mu, status))
    }
  }

  invisible(all_pass)
}


#' Verify Derivative via Finite Difference
#'
#' Verifies the analytic derivative of the conditional mean against
#' finite difference approximation.
#'
#' @param J Sample size (integer >= 1).
#' @param alpha Concentration parameter (positive scalar).
#' @param eps Finite difference step size (default: 1e-6).
#' @param tol Tolerance for comparison (default: 1e-5).
#' @param verbose Logical; if TRUE, print results.
#'
#' @return Logical; TRUE if derivative matches finite difference.
#'
#' @examples
#' \dontrun{
#' verify_derivative(50, 2.0)
#'
#' }
#' @keywords internal
verify_derivative <- function(J, alpha, eps = 1e-6, tol = 1e-5, verbose = TRUE) {
  deriv_analytic <- dmean_dalpha(J, alpha)
  deriv_fd <- (mean_K_given_alpha(J, alpha + eps) -
                 mean_K_given_alpha(J, alpha - eps)) / (2 * eps)

  err <- abs(deriv_analytic - deriv_fd)
  passed <- err < tol

  if (verbose) {
    cat(sprintf("Derivative verification (J=%d, alpha=%.2f):\n", J, alpha))
    cat(sprintf("  Analytic:    %.10f\n", deriv_analytic))
    cat(sprintf("  Finite diff: %.10f\n", deriv_fd))
    cat(sprintf("  Error:       %.2e [%s]\n", err, if (passed) "PASS" else "FAIL"))
  }

  invisible(passed)
}
