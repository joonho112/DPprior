# =============================================================================
# Module 05: Marginal Moments of K_J under alpha ~ Gamma(a, b)
# =============================================================================
#
# This module provides computation of marginal moments E[K_J] and Var(K_J)
# when the concentration parameter alpha has a Gamma(a, b) prior distribution.
#
# Theory Background (RN-04 Section 3):
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
# 1. Marginal overdispersion: Var(K_J) > E[K_J] (despite conditional underdispersion)
# 2. Mean bounds: 1 <= E[K_J] <= J
# 3. Variance inflation: Var(K_J) > Var(K_J | alpha = E[alpha])
#
# Author: JoonHo Lee
# Date: December 2025
# Part of: DPprior R Package
# Reference: RN-01 Corollary 1, RN-04 Section 3
# Dependencies: Module 02 (quadrature), Module 03 (conditional moments)
# =============================================================================


# =============================================================================
# Core Exported Functions
# =============================================================================

#' Exact Marginal Moments of K_J under Gamma Prior
#'
#' Computes the exact marginal mean \eqn{E[K_J]} and variance \eqn{Var(K_J)}
#' when the DP concentration parameter follows a Gamma(a, b) prior.
#'
#' @param J Integer; sample size (positive integer >= 1).
#' @param a Numeric; shape parameter of Gamma prior (> 0).
#' @param b Numeric; rate parameter of Gamma prior (> 0).
#' @param M Integer; number of quadrature nodes (default: 80).
#'
#' @return A named list with components:
#'   \describe{
#'     \item{\code{mean}}{Marginal mean \eqn{E[K_J | a, b]}}
#'     \item{\code{var}}{Marginal variance \eqn{Var(K_J | a, b)}}
#'     \item{\code{sd}}{Marginal standard deviation}
#'     \item{\code{cv}}{Coefficient of variation (sd/mean)}
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
#' \strong{Key properties:}
#' \itemize{
#'   \item The mean is bounded: \eqn{1 \leq E[K_J] \leq J}
#'   \item Despite conditional underdispersion (\eqn{v_J(\alpha) < \mu_J(\alpha)}),
#'         the marginal distribution is typically overdispersed
#'   \item Marginal variance exceeds conditional variance at \eqn{\alpha = E[\alpha]}
#' }
#'
#' @examples
#' # Example from RN-01: J=50, Gamma(1.5, 0.5) prior
#' result <- exact_K_moments(50, 1.5, 0.5)
#' print(result)
#'
#' # Compare with conditional variance at E[alpha] = 3
#' cond_var <- var_K_given_alpha(50, 3.0)
#' result$var > cond_var  # TRUE: marginal > conditional
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
#' @export
exact_K_moments <- function(J, a, b, M = .QUAD_NODES_DEFAULT) {
  # Input validation
  assert_valid_J(J)
  assert_positive(a, "a")
  assert_positive(b, "b")

  if (!is.numeric(M) || length(M) != 1L || !is.finite(M) ||
      M != floor(M) || M < 1L) {
    stop("M must be a positive integer", call. = FALSE)
  }

  M <- as.integer(M)

  # Build quadrature for Gamma(a, b)
  quad <- build_gamma_quadrature(a, b, M)
  alphas <- quad$alpha_nodes
  w <- quad$weights_normalized

  # Compute conditional moments at each quadrature node
  # Using vectorized implementations from Module 03
  mu_vec <- mean_K_given_alpha(J, alphas)
  v_vec <- var_K_given_alpha(J, alphas)

  # Marginal mean: E[K_J] = E[mu_J(alpha)]
  mean_K <- sum(w * mu_vec)

  # Marginal variance via Law of Total Variance:
  # Var(K_J) = E[Var(K_J | alpha)] + Var(E[K_J | alpha])
  #          = E[v_J(alpha)] + E[mu_J(alpha)^2] - E[mu_J(alpha)]^2
  E_mu_sq <- sum(w * mu_vec^2)
  E_v <- sum(w * v_vec)
  var_K <- E_v + E_mu_sq - mean_K^2

  # Numerical safety: enforce non-negativity
  var_K <- max(0, var_K)
  sd_K <- sqrt(var_K)
  cv_K <- if (mean_K > 0) sd_K / mean_K else Inf

  list(
    mean = mean_K,
    var = var_K,
    sd = sd_K,
    cv = cv_K
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
#'
#' @return Named numeric vector \code{c(mean = ..., var = ...)}.
#'
#' @examples
#' K_moments(50, 2.0, 1.0)
#'
#' @seealso \code{\link{exact_K_moments}} for full output
#'
#' @export
K_moments <- function(J, a, b, M = .QUAD_NODES_DEFAULT) {
  result <- exact_K_moments(J, a, b, M)
  c(mean = result$mean, var = result$var)
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
#' # Typical overdispersion
#' vir <- variance_inflation_ratio(50, 2.0, 1.0)
#' vir > 1  # TRUE: overdispersed
#'
#' # Compare across prior specifications
#' variance_inflation_ratio(50, 2.0, 0.5)  # Higher uncertainty in alpha
#' variance_inflation_ratio(50, 8.0, 4.0)  # Same mean, lower variance in alpha
#'
#' @keywords internal
#' @export
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
#' The NegBin approximation (A1 method from RN-03) assumes:
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
#' # Large approximation error for small J
#' compare_to_negbin(50, 1.5, 0.5)
#'
#' # Error decreases with J
#' compare_to_negbin(300, 1.5, 0.5)
#'
#' @keywords internal
#' @export
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
#'   }
#'
#' @details
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
#' result <- marginal_moments_with_jacobian(50, 2.0, 1.0)
#' result$jacobian
#'
#' @seealso \code{\link{exact_K_moments}}, Module 07 (Jacobian)
#'
#' @keywords internal
#' @export
marginal_moments_with_jacobian <- function(J, a, b, M = .QUAD_NODES_DEFAULT) {
  # Input validation
  assert_valid_J(J)
  assert_positive(a, "a")
  assert_positive(b, "b")

  # Build quadrature
  quad <- build_gamma_quadrature(a, b, M)
  alphas <- quad$alpha_nodes
  w <- quad$weights_normalized

  # Conditional moments
  mu_vec <- mean_K_given_alpha(J, alphas)
  v_vec <- var_K_given_alpha(J, alphas)
  m2_vec <- v_vec + mu_vec^2  # E[K^2 | alpha]

  # Marginal moments
  mean_K <- sum(w * mu_vec)
  E_v <- sum(w * v_vec)
  E_mu_sq <- sum(w * mu_vec^2)
  var_K <- max(0, E_v + E_mu_sq - mean_K^2)

  # Score functions for Gamma(a, b)
  # s_a(alpha) = log(b) - psi(a) + log(alpha)
  # s_b(alpha) = a/b - alpha
  score_a <- log(b) - digamma(a) + log(alphas)
  score_b <- a / b - alphas

  # Derivatives of moments via score identity
  # d/d_theta E[f] = E[f * s_theta]

  # Derivatives of mean
  dmean_da <- sum(w * mu_vec * score_a)
  dmean_db <- sum(w * mu_vec * score_b)

  # Derivatives of E[mu^2] and E[v]
  dmu2_da <- sum(w * mu_vec^2 * score_a)
  dmu2_db <- sum(w * mu_vec^2 * score_b)
  dv_da <- sum(w * v_vec * score_a)
  dv_db <- sum(w * v_vec * score_b)

  # Derivatives of variance via chain rule
  # V = E[v] + E[mu^2] - E[mu]^2
  # dV/d_theta = dE[v]/d_theta + dE[mu^2]/d_theta - 2*E[mu]*dE[mu]/d_theta
  dvar_da <- dv_da + dmu2_da - 2 * mean_K * dmean_da
  dvar_db <- dv_db + dmu2_db - 2 * mean_K * dmean_db

  # Construct Jacobian matrix
  jacobian <- matrix(
    c(dmean_da, dmean_db, dvar_da, dvar_db),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("mean", "var"), c("a", "b"))
  )

  list(
    mean = mean_K,
    var = var_K,
    jacobian = jacobian
  )
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
#' verify_marginal_moments(50, 2.0, 1.0)
#'
#' @export
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

  # Test 3: Marginal variance > conditional variance at E[alpha]
  alpha_mean <- a / b
  cond_var <- var_K_given_alpha(J, alpha_mean)
  test3 <- result$var > cond_var - 1e-10  # Small tolerance
  if (isTRUE(verbose)) {
    cat(sprintf("  Test 3 (Var inflation): Var(K)=%.4f > Var(K|E[alpha])=%.4f? %s\n",
                result$var, cond_var, if (test3) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test3

  # Test 4: CV reasonable (between 0 and 1 for typical cases)
  test4 <- result$cv > 0
  if (isTRUE(verbose)) {
    cat(sprintf("  Test 4 (CV positive): CV > 0? %s\n",
                if (test4) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test4

  if (isTRUE(verbose)) {
    cat(sprintf("\n  Overall: %s\n", if (all_pass) "PASS" else "FAIL"))
  }

  invisible(all_pass)
}


#' Verify Quadrature Convergence for Marginal Moments
#'
#' Tests that marginal moments converge as the number of quadrature nodes increases.
#'
#' @param J Integer; sample size.
#' @param a Numeric; shape parameter.
#' @param b Numeric; rate parameter.
#' @param M_values Integer vector; numbers of quadrature nodes to test.
#' @param verbose Logical; if \code{TRUE}, print detailed results.
#'
#' @return Data frame with convergence results.
#'
#' @examples
#' verify_quadrature_convergence(50, 1.5, 0.5)
#'
#' @keywords internal
#' @export
verify_quadrature_convergence <- function(J, a, b,
                                          M_values = c(20, 40, 60, 80, 100, 120),
                                          verbose = TRUE) {
  results <- data.frame(
    M = integer(0),
    mean = numeric(0),
    var = numeric(0),
    mean_change = numeric(0),
    var_change = numeric(0)
  )

  prev_mean <- NA
  prev_var <- NA

  for (M in M_values) {
    moments <- exact_K_moments(J, a, b, M)

    mean_change <- if (is.na(prev_mean)) NA else abs(moments$mean - prev_mean)
    var_change <- if (is.na(prev_var)) NA else abs(moments$var - prev_var)

    results <- rbind(results, data.frame(
      M = M,
      mean = moments$mean,
      var = moments$var,
      mean_change = mean_change,
      var_change = var_change
    ))

    prev_mean <- moments$mean
    prev_var <- moments$var
  }

  if (isTRUE(verbose)) {
    cat(sprintf("Quadrature Convergence (J=%d, a=%.2f, b=%.2f):\n", J, a, b))
    print(results, row.names = FALSE)

    # Check if converged
    final_change <- results$mean_change[nrow(results)]
    if (!is.na(final_change) && final_change < 1e-8) {
      cat("\nMean: CONVERGED to machine precision\n")
    }
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
#' verify_moments_marginal_all()
#'
#' @export
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
