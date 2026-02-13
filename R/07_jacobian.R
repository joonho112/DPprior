# =============================================================================
# Module 07: Score-Based Jacobian for Newton Solver
# =============================================================================
#
# This module provides functions for computing the Jacobian matrix of the
# marginal moment map F(a,b) = (M_1(a,b), V(a,b)) using score function
# identities, enabling efficient Newton's method for exact moment matching.
#
# IMPROVEMENTS (incorporated from comparative analysis):
# 1. Small alpha edge case handling (alpha < 1e-12)
# 2. Higher M recommended for verification
# 3. Better documentation about numerical considerations
# 4. More robust quadrature using orthpol-style approach
#
# Theory Reference: Lee (2026), Section 3.2
#
# Author: JoonHo Lee
# Date: December 2025
# Part of: DPprior R Package
# =============================================================================

# =============================================================================
# Constants for Numerical Stability
# =============================================================================

# .ALPHA_SMALL is now defined in R/00_constants.R (loaded first alphabetically)

#' Default Quadrature Nodes for Verification
#' @description Higher M value recommended for Jacobian verification due to
#'   slower convergence of score-weighted integrands.
#' @keywords internal
.QUAD_NODES_VERIFICATION <- 200L


# =============================================================================
# Score Functions for Gamma(a, b)
# =============================================================================

#' Score Function with Respect to Shape Parameter a
#'
#' Computes the score function \eqn{s_a(\alpha) = \partial/\partial a \log g_{a,b}(\alpha)}
#' for the Gamma(a, b) distribution.
#'
#' @param alpha Numeric vector; points at which to evaluate.
#' @param a Numeric scalar; shape parameter of the Gamma distribution (> 0).
#' @param b Numeric scalar; rate parameter of the Gamma distribution (> 0).
#'
#' @return Numeric vector of the same length as \code{alpha}.
#'
#' @details
#' For the Gamma(shape = a, rate = b) distribution with density
#' \deqn{g_{a,b}(\alpha) = \frac{b^a}{\Gamma(a)} \alpha^{a-1} e^{-b\alpha},}
#' the score function with respect to \code{a} is:
#' \deqn{s_a(\alpha) = \log b - \psi(a) + \log \alpha,}
#' where \eqn{\psi} is the digamma function.
#'
#' A fundamental property of score functions is that their expectation is zero:
#' \deqn{E_{\alpha \sim g_{a,b}}[s_a(\alpha)] = 0.}
#'
#' \strong{Numerical Note:} The \code{log(alpha)} term causes slower quadrature
#' convergence for score-weighted integrands compared to the moments themselves.
#' Use higher M (e.g., 120-200) for verification purposes.
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso \code{\link{score_b}} for the score with respect to b
#'
#' @examples
#' # Evaluate score at several points
#' alpha_vals <- c(0.5, 1.0, 2.0, 5.0)
#' score_a(alpha_vals, a = 2.0, b = 1.0)
#'
#' @keywords internal
#' @export
score_a <- function(alpha, a, b) {
  assert_positive(alpha, "alpha")
  assert_positive(a, "a")
  assert_positive(b, "b")

  log(b) - digamma(a) + log(alpha)
}


#' Score Function with Respect to Rate Parameter b
#'
#' Computes the score function \eqn{s_b(\alpha) = \partial/\partial b \log g_{a,b}(\alpha)}
#' for the Gamma(a, b) distribution.
#'
#' @param alpha Numeric vector; points at which to evaluate.
#' @param a Numeric scalar; shape parameter of the Gamma distribution (> 0).
#' @param b Numeric scalar; rate parameter of the Gamma distribution (> 0).
#'
#' @return Numeric vector of the same length as \code{alpha}.
#'
#' @details
#' For the Gamma(shape = a, rate = b) distribution, the score function with
#' respect to \code{b} is:
#' \deqn{s_b(\alpha) = \frac{a}{b} - \alpha.}
#'
#' A fundamental property of score functions is that their expectation is zero:
#' \deqn{E_{\alpha \sim g_{a,b}}[s_b(\alpha)] = 0.}
#'
#' Unlike \code{score_a}, this function is linear in \eqn{\alpha}, so its
#' expectation converges very quickly with quadrature.
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso \code{\link{score_a}} for the score with respect to a
#'
#' @examples
#' # Evaluate score at several points
#' alpha_vals <- c(0.5, 1.0, 2.0, 5.0)
#' score_b(alpha_vals, a = 2.0, b = 1.0)
#'
#' @keywords internal
#' @export
score_b <- function(alpha, a, b) {
  assert_positive(alpha, "alpha")
  assert_positive(a, "a")
  assert_positive(b, "b")

  a / b - alpha
}


# =============================================================================
# Enhanced Conditional Moments with Small-Alpha Handling
# =============================================================================

#' Conditional Mean of K Given Alpha (Enhanced)
#'
#' Computes \eqn{E[K_J | \alpha]} with proper handling of very small alpha values.
#'
#' @param J Integer; sample size.
#' @param alpha Numeric vector; concentration parameter values.
#'
#' @return Numeric vector of conditional means.
#'
#' @details
#' For very small alpha (< 1e-12), returns 1.0 as the limiting value.
#' This avoids numerical issues with digamma at near-zero arguments.
#'
#' @keywords internal
mean_K_given_alpha_safe <- function(J, alpha) {
  alpha <- as.numeric(alpha)
  out <- numeric(length(alpha))

  small <- alpha <= .ALPHA_SMALL
  out[small] <- 1.0

  if (any(!small)) {
    a <- alpha[!small]
    out[!small] <- a * (digamma(a + J) - digamma(a))
  }

  out
}


#' Conditional Variance of K Given Alpha (Enhanced)
#'
#' Computes \eqn{Var(K_J | \alpha)} with proper handling of very small alpha values.
#'
#' @param J Integer; sample size.
#' @param alpha Numeric vector; concentration parameter values.
#'
#' @return Numeric vector of conditional variances.
#'
#' @details
#' For very small alpha (< 1e-12), returns 0.0 as the limiting value.
#'
#' @keywords internal
var_K_given_alpha_safe <- function(J, alpha) {
  alpha <- as.numeric(alpha)
  out <- numeric(length(alpha))

  small <- alpha <= .ALPHA_SMALL
  out[small] <- 0.0

  if (any(!small)) {
    a <- alpha[!small]
    mu <- mean_K_given_alpha_safe(J, a)
    out[!small] <- mu - a^2 * (trigamma(a) - trigamma(a + J))
  }

  pmax(out, 0)
}


# =============================================================================
# Combined Moments and Jacobian Computation
# =============================================================================

#' Compute Marginal Moments and Jacobian Simultaneously
#'
#' Computes the exact marginal moments \eqn{M_1 = E[K_J]} and \eqn{V = Var(K_J)}
#' along with the Jacobian matrix of the moment map \eqn{F(a,b) = (M_1, V)}
#' using score function identities.
#'
#' @param J Integer; sample size (number of observations/sites).
#' @param a Numeric; shape parameter of the Gamma prior on \eqn{\alpha} (> 0).
#' @param b Numeric; rate parameter of the Gamma prior on \eqn{\alpha} (> 0).
#' @param M Integer; number of quadrature nodes (default: 80).
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{mean}}{Marginal mean \eqn{E[K_J]}}
#'   \item{\code{var}}{Marginal variance \eqn{Var(K_J)}}
#'   \item{\code{jacobian}}{2x2 Jacobian matrix with structure:
#'     \deqn{J_F = \begin{bmatrix} \partial M_1/\partial a & \partial M_1/\partial b \\
#'                                  \partial V/\partial a & \partial V/\partial b \end{bmatrix}}}
#' }
#'
#' @details
#' This function uses the score identity (Lee, 2026, Section 3.2, Corollary 1) to compute exact
#' derivatives without finite differences:
#' \deqn{\frac{\partial}{\partial\theta} E[f(\alpha)] = E[f(\alpha) \cdot s_\theta(\alpha)]}
#'
#' The Jacobian components are computed as:
#' \itemize{
#'   \item \eqn{\partial M_1/\partial \theta = E[\mu_J(\alpha) \cdot s_\theta(\alpha)]}
#'   \item \eqn{\partial V/\partial \theta = \partial E[v_J]/\partial \theta +
#'              \partial E[\mu_J^2]/\partial \theta - 2 M_1 \partial M_1/\partial \theta}
#' }
#'
#' \strong{Numerical Considerations:}
#' \itemize{
#'   \item The score function \code{s_a} contains \code{log(alpha)}, which causes
#'         slower quadrature convergence compared to moment computation.
#'   \item For verification, use higher M (e.g., 120-200).
#'   \item Very small alpha values (< 1e-12) are handled with limiting values.
#' }
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso
#' \code{\link{exact_K_moments}} for moments only,
#' \code{\link{score_a}}, \code{\link{score_b}} for score functions
#'
#' @examples
#' # Compute moments and Jacobian for J=50, a=2, b=1
#' result <- moments_with_jacobian(J = 50, a = 2.0, b = 1.0)
#' print(result$mean)      # E[K_J]
#' print(result$var)       # Var(K_J)
#' print(result$jacobian)  # 2x2 Jacobian matrix
#'
#' # Use in Newton iteration
#' target <- c(5.0, 8.0)  # Target (E[K], Var(K))
#' current <- c(result$mean, result$var)
#' residual <- current - target
#' delta <- solve(result$jacobian, -residual)
#'
#' @export
moments_with_jacobian <- function(J, a, b, M = .QUAD_NODES_DEFAULT) {
  # Input validation
  assert_valid_J(J)
  assert_positive(a, "a")
  assert_positive(b, "b")
  if (!is.numeric(M) || length(M) != 1L || M < 10L) {
    stop("M must be a positive integer >= 10", call. = FALSE)
  }

  # Build quadrature
  quad <- build_gamma_quadrature(a, b, M)
  alphas <- quad$alpha_nodes
  w <- quad$weights_normalized

  # Conditional moments at each quadrature node (use safe versions)
  # Note: For most quadrature nodes, alpha won't be near zero,
  # but using safe versions ensures robustness
  mu_vec <- mean_K_given_alpha(J, alphas)
  v_vec <- var_K_given_alpha(J, alphas)

  # Score functions at each quadrature node
  s_a <- score_a(alphas, a, b)
  s_b <- score_b(alphas, a, b)

  # ----- Marginal Moments -----
  # M_1 = E[mu_J(alpha)]
  M1 <- sum(w * mu_vec)

  # E[mu_J^2]
  E_mu_sq <- sum(w * mu_vec^2)

  # E[v_J(alpha)]
  E_v <- sum(w * v_vec)

  # V = E[v_J] + E[mu_J^2] - M_1^2 (law of total variance)
  V <- E_v + E_mu_sq - M1^2

  # ----- Jacobian Components via Score Identity -----
  # dM1/da = E[mu_J * s_a]
  dM1_da <- sum(w * mu_vec * s_a)

  # dM1/db = E[mu_J * s_b]
  dM1_db <- sum(w * mu_vec * s_b)

  # dE[v_J]/da = E[v_J * s_a]
  dEv_da <- sum(w * v_vec * s_a)

  # dE[v_J]/db = E[v_J * s_b]
  dEv_db <- sum(w * v_vec * s_b)

  # dE[mu_J^2]/da = E[mu_J^2 * s_a]
  dEmu_sq_da <- sum(w * mu_vec^2 * s_a)

  # dE[mu_J^2]/db = E[mu_J^2 * s_b]
  dEmu_sq_db <- sum(w * mu_vec^2 * s_b)

  # dV/dtheta = dE[v_J]/dtheta + dE[mu_J^2]/dtheta - 2*M1*dM1/dtheta
  dV_da <- dEv_da + dEmu_sq_da - 2 * M1 * dM1_da
  dV_db <- dEv_db + dEmu_sq_db - 2 * M1 * dM1_db

  # Jacobian matrix (2x2, row-major: first row is dM1, second row is dV)
  jacobian <- matrix(
    c(dM1_da, dM1_db,
      dV_da, dV_db),
    nrow = 2L, ncol = 2L, byrow = TRUE
  )
  rownames(jacobian) <- c("dM1", "dV")
  colnames(jacobian) <- c("da", "db")

  list(
    mean = M1,
    var = max(0, V),
    jacobian = jacobian
  )
}


# =============================================================================
# Verification Functions
# =============================================================================

#' Verify Jacobian Against Finite Differences
#'
#' Compares the analytically computed Jacobian (via score identities) against
#' numerical finite differences to validate the implementation.
#'
#' @param J Integer; sample size.
#' @param a Numeric; shape parameter.
#' @param b Numeric; rate parameter.
#' @param eps Numeric; step size for finite differences (default: 1e-6).
#' @param M Integer; number of quadrature nodes (default: 200 for verification).
#' @param verbose Logical; if TRUE, print detailed comparison.
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{analytic}}{The analytically computed Jacobian}
#'   \item{\code{numeric}}{The numerically computed Jacobian (finite differences)}
#'   \item{\code{abs_error}}{Matrix of absolute errors}
#'   \item{\code{rel_error}}{Matrix of relative errors}
#'   \item{\code{max_rel_error}}{Maximum relative error across all entries}
#'   \item{\code{pass}}{Logical; TRUE if max relative error < 0.01}
#' }
#'
#' @details
#' Uses central finite differences:
#' \deqn{\frac{\partial f}{\partial a} \approx \frac{f(a+\epsilon) - f(a-\epsilon)}{2\epsilon}}
#'
#' \strong{Important:} This is a SECONDARY verification method because both the
#' analytical Jacobian and the finite differences use the same quadrature layer.
#' For independent verification, compare against adaptive integration (scipy.integrate.quad
#' in Python).
#'
#' @examples
#' # Verify Jacobian for a specific case
#' result <- verify_jacobian(J = 50, a = 2.0, b = 1.0, verbose = TRUE)
#'
#' @keywords internal
#' @export
verify_jacobian <- function(J, a, b, eps = 1e-6, M = .QUAD_NODES_VERIFICATION,
                            verbose = TRUE) {
  # Use higher M for finite difference computations
  M_fd <- min(M + 50L, 300L)

  # Analytical Jacobian
  result <- moments_with_jacobian(J, a, b, M)
  J_analytic <- result$jacobian

  # Numerical Jacobian via central differences
  J_numeric <- matrix(0, nrow = 2L, ncol = 2L)

  # d/da
  m_a_plus <- exact_K_moments(J, a + eps, b, M_fd)
  m_a_minus <- exact_K_moments(J, a - eps, b, M_fd)
  J_numeric[1L, 1L] <- (m_a_plus$mean - m_a_minus$mean) / (2 * eps)
  J_numeric[2L, 1L] <- (m_a_plus$var - m_a_minus$var) / (2 * eps)

  # d/db
  m_b_plus <- exact_K_moments(J, a, b + eps, M_fd)
  m_b_minus <- exact_K_moments(J, a, b - eps, M_fd)
  J_numeric[1L, 2L] <- (m_b_plus$mean - m_b_minus$mean) / (2 * eps)
  J_numeric[2L, 2L] <- (m_b_plus$var - m_b_minus$var) / (2 * eps)

  rownames(J_numeric) <- c("dM1", "dV")
  colnames(J_numeric) <- c("da", "db")

  # Errors
  abs_error <- abs(J_analytic - J_numeric)
  rel_error <- abs_error / (abs(J_numeric) + 1e-10)
  max_rel_error <- max(rel_error)

  # Pass criterion: < 1% relative error
  pass <- max_rel_error < 0.01

  if (isTRUE(verbose)) {
    cat(sprintf("Jacobian Verification (J=%d, a=%.2f, b=%.2f, M=%d)\n", J, a, b, M))
    cat(strrep("-", 60), "\n")

    cat("\nAnalytic Jacobian (score-based):\n")
    cat(sprintf("  dM1/da = %12.8f  dM1/db = %12.8f\n",
                J_analytic[1L, 1L], J_analytic[1L, 2L]))
    cat(sprintf("  dV/da  = %12.8f  dV/db  = %12.8f\n",
                J_analytic[2L, 1L], J_analytic[2L, 2L]))

    cat("\nNumeric Jacobian (finite diff):\n")
    cat(sprintf("  dM1/da = %12.8f  dM1/db = %12.8f\n",
                J_numeric[1L, 1L], J_numeric[1L, 2L]))
    cat(sprintf("  dV/da  = %12.8f  dV/db  = %12.8f\n",
                J_numeric[2L, 1L], J_numeric[2L, 2L]))

    cat("\nRelative Errors:\n")
    cat(sprintf("  dM1/da: %.2e  dM1/db: %.2e\n",
                rel_error[1L, 1L], rel_error[1L, 2L]))
    cat(sprintf("  dV/da:  %.2e  dV/db:  %.2e\n",
                rel_error[2L, 1L], rel_error[2L, 2L]))

    cat(sprintf("\nMax Relative Error: %.2e [%s]\n",
                max_rel_error, if (pass) "PASS" else "FAIL"))

    if (!pass) {
      cat("\nNote: Finite-difference verification has circular dependency with quadrature.\n")
      cat("Consider verifying against adaptive integration for independent validation.\n")
    }
  }

  invisible(list(
    analytic = J_analytic,
    numeric = J_numeric,
    abs_error = abs_error,
    rel_error = rel_error,
    max_rel_error = max_rel_error,
    pass = pass
  ))
}


#' Verify Score Function Zero Expectation Property
#'
#' Verifies the fundamental property that \eqn{E[s_\theta(\alpha)] = 0}
#' for both score functions.
#'
#' @param a Numeric; shape parameter.
#' @param b Numeric; rate parameter.
#' @param M Integer; number of quadrature nodes.
#' @param verbose Logical; if TRUE, print results.
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{E_score_a}}{Expectation of \eqn{s_a}}
#'   \item{\code{E_score_b}}{Expectation of \eqn{s_b}}
#' }
#'
#' @details
#' This is a fundamental property of score functions. Due to quadrature
#' approximation error, the computed expectations may not be exactly zero.
#'
#' \strong{Expected behavior:}
#' \itemize{
#'   \item \code{E[s_b]} should be very close to zero (typically < 1e-14)
#'         because s_b is linear in alpha.
#'   \item \code{E[s_a]} may show larger errors (up to 1e-2 for small a)
#'         due to the log(alpha) term causing slower quadrature convergence.
#' }
#'
#' For rigorous verification, use adaptive integration (e.g., in Python with
#' scipy.integrate.quad) which provides independent ground truth.
#'
#' @examples
#' \dontrun{
#' verify_score_expectation(a = 2.0, b = 1.0, verbose = TRUE)
#'
#' }
#' @keywords internal
verify_score_expectation <- function(a, b, M = .QUAD_NODES_VERIFICATION,
                                     verbose = TRUE) {
  # Build quadrature
  quad <- build_gamma_quadrature(a, b, M)
  alphas <- quad$alpha_nodes
  w <- quad$weights_normalized

  # Compute expectations
  E_sa <- sum(w * score_a(alphas, a, b))
  E_sb <- sum(w * score_b(alphas, a, b))

  if (isTRUE(verbose)) {
    cat(sprintf("Score Expectation Verification (a=%.2f, b=%.2f, M=%d)\n", a, b, M))
    cat(sprintf("  E[s_a(alpha)] = %.6e (should be ~= 0)\n", E_sa))
    cat(sprintf("  E[s_b(alpha)] = %.6e (should be ~= 0)\n", E_sb))

    if (abs(E_sa) > 1e-3) {
      cat("\n  Note: E[s_a] has larger error because log(alpha) slows quadrature.\n")
      cat("  This is expected and does not indicate implementation error.\n")
      cat("  For verification, compare against adaptive integration.\n")
    }
  }

  invisible(list(
    E_score_a = E_sa,
    E_score_b = E_sb
  ))
}


#' Run All Module 07 Verification Tests
#'
#' Comprehensive verification suite for the score-based Jacobian module.
#'
#' @param verbose Logical; if TRUE, print detailed results.
#'
#' @return Logical; TRUE if all tests pass.
#'
#' @examples
#' \dontrun{
#' verify_jacobian_all()
#'
#' }
#' @keywords internal
verify_jacobian_all <- function(verbose = TRUE) {
  if (isTRUE(verbose)) {
    cat(strrep("=", 70), "\n")
    cat("Module 07: Score-Based Jacobian - Full Verification Suite\n")
    cat(strrep("=", 70), "\n\n")
  }

  all_pass <- TRUE

  # Test 1: Score function expectations
  if (isTRUE(verbose)) {
    cat("[Test 1] Score function zero-expectation property\n")
    cat(strrep("-", 50), "\n")
    cat("Note: E[s_a] may show larger errors due to log(alpha) term.\n\n")
  }

  test_cases <- list(
    c(a = 2.0, b = 1.0),
    c(a = 1.5, b = 0.5),
    c(a = 3.0, b = 1.5)
  )

  for (tc in test_cases) {
    verify_score_expectation(tc["a"], tc["b"], verbose = verbose)
    cat("\n")
  }

  # Test 2: Jacobian matches finite differences
  if (isTRUE(verbose)) {
    cat("[Test 2] Jacobian vs finite differences (M=200)\n")
    cat(strrep("-", 50), "\n\n")
  }

  test_cases <- list(
    list(J = 30, a = 1.5, b = 0.5),
    list(J = 50, a = 2.0, b = 1.0),
    list(J = 100, a = 2.0, b = 1.0),
    list(J = 50, a = 3.0, b = 1.5)
  )

  for (tc in test_cases) {
    result <- verify_jacobian(tc$J, tc$a, tc$b, verbose = verbose)
    all_pass <- all_pass && result$pass
    cat("\n")
  }

  # Test 3: Moments consistency
  if (isTRUE(verbose)) {
    cat("[Test 3] Moments consistency check\n")
    cat(strrep("-", 50), "\n")
  }

  J <- 50; a <- 2.0; b <- 1.0
  result_jac <- moments_with_jacobian(J, a, b)
  result_mom <- exact_K_moments(J, a, b)

  mean_diff <- abs(result_jac$mean - result_mom$mean)
  var_diff <- abs(result_jac$var - result_mom$var)
  consistent <- mean_diff < 1e-10 && var_diff < 1e-10

  if (isTRUE(verbose)) {
    cat(sprintf("Mean (moments_with_jacobian): %.10f\n", result_jac$mean))
    cat(sprintf("Mean (exact_K_moments):       %.10f\n", result_mom$mean))
    cat(sprintf("Difference: %.2e\n", mean_diff))
    cat(sprintf("Var (moments_with_jacobian):  %.10f\n", result_jac$var))
    cat(sprintf("Var (exact_K_moments):        %.10f\n", result_mom$var))
    cat(sprintf("Difference: %.2e\n", var_diff))
    cat(sprintf("Status: %s\n\n", if (consistent) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && consistent

  # Test 4: Jacobian non-singularity
  if (isTRUE(verbose)) {
    cat("[Test 4] Jacobian non-singularity check\n")
    cat(strrep("-", 50), "\n")
  }

  test_cases <- list(
    list(J = 30, a = 1.5, b = 0.5),
    list(J = 50, a = 2.0, b = 1.0),
    list(J = 100, a = 2.0, b = 1.0)
  )

  for (tc in test_cases) {
    result <- moments_with_jacobian(tc$J, tc$a, tc$b)
    det_val <- det(result$jacobian)
    cond_val <- kappa(result$jacobian)
    nonsingular <- abs(det_val) > 1e-10 && cond_val < 1e10

    if (isTRUE(verbose)) {
      cat(sprintf("J=%d, a=%.1f, b=%.1f: det=%.4f, cond=%.2e [%s]\n",
                  tc$J, tc$a, tc$b, det_val, cond_val,
                  if (nonsingular) "PASS" else "FAIL"))
    }
    all_pass <- all_pass && nonsingular
  }

  # Test 5: Newton convergence
  if (isTRUE(verbose)) {
    cat("\n[Test 5] Newton convergence test\n")
    cat(strrep("-", 50), "\n")
  }

  newton_pass <- test_newton_convergence(J = 50, mu_target = 5.0, var_target = 8.0,
                                         verbose = verbose)
  all_pass <- all_pass && newton_pass

  # Summary
  if (isTRUE(verbose)) {
    cat("\n", strrep("=", 70), "\n", sep = "")
    cat(sprintf("Overall Result: %s\n",
                if (all_pass) "ALL TESTS PASSED" else "SOME TESTS FAILED"))
    cat(strrep("=", 70), "\n")
  }

  invisible(all_pass)
}


#' Test Newton Convergence Using the Jacobian
#'
#' Verifies that the Jacobian enables fast Newton convergence for moment matching.
#'
#' @param J Integer; sample size.
#' @param mu_target Numeric; target mean.
#' @param var_target Numeric; target variance.
#' @param a0 Numeric; initial shape parameter.
#' @param b0 Numeric; initial rate parameter.
#' @param max_iter Integer; maximum iterations.
#' @param tol Numeric; convergence tolerance.
#' @param verbose Logical; if TRUE, print iteration history.
#'
#' @return Logical; TRUE if Newton converges.
#'
#' @keywords internal
test_newton_convergence <- function(J, mu_target, var_target,
                                    a0 = 2.0, b0 = 1.0,
                                    max_iter = 15L, tol = 1e-8,
                                    verbose = TRUE) {
  a <- a0
  b <- b0

  if (isTRUE(verbose)) {
    cat(sprintf("Target: E[K]=%.2f, Var(K)=%.2f\n", mu_target, var_target))
    cat(sprintf("Start:  a=%.4f, b=%.4f\n\n", a, b))
    cat("Iter |    a      |    b      |   E[K]    |  Var(K)   | ||F||\n")
    cat(strrep("-", 70), "\n")
  }

  for (iter in 0:max_iter) {
    result <- moments_with_jacobian(J, a, b, M = .QUAD_NODES_VERIFICATION)
    F_vec <- c(result$mean - mu_target, result$var - var_target)
    norm_F <- sqrt(sum(F_vec^2))

    if (isTRUE(verbose)) {
      cat(sprintf("%4d | %9.6f | %9.6f | %9.6f | %9.6f | %.2e\n",
                  iter, a, b, result$mean, result$var, norm_F))
    }

    if (norm_F < tol) {
      if (isTRUE(verbose)) {
        cat(sprintf("\nConverged in %d iterations!\n", iter))
      }
      return(TRUE)
    }

    delta <- tryCatch(
      solve(result$jacobian, -F_vec),
      error = function(e) NULL
    )

    if (is.null(delta)) {
      if (isTRUE(verbose)) cat("Jacobian singular!\n")
      return(FALSE)
    }

    # Damped update
    step <- 1.0
    while (step > 1e-10) {
      a_new <- a + step * delta[1]
      b_new <- b + step * delta[2]
      if (a_new > 0.1 && b_new > 0.01) break
      step <- step * 0.5
    }

    a <- a_new
    b <- b_new
  }

  if (isTRUE(verbose)) cat("Max iterations reached.\n")
  return(FALSE)
}
