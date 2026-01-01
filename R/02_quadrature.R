# =============================================================================
# Module 02: Gauss-Laguerre Quadrature for Gamma Distributions
# =============================================================================
#
# This module provides numerical integration against Gamma distributions
# using generalized Gauss-Laguerre quadrature.
#
# Theory Background:
# -----------------
# Standard Laguerre quadrature integrates:
#   int_0^inf f(x) x^beta e^{-x} dx = sum w_m f(x_m)
#
# For Gamma(a, b) prior expectations E[g(alpha)]:
#   E[g(alpha)] = int_0^inf g(alpha) p(alpha|a,b) d_alpha
#
# Using change of variables x = b*alpha:
#   E[g(alpha)] = (1/Gamma(a)) int_0^inf g(x/b) x^{a-1} e^{-x} dx
#
# Generalized Gauss-Laguerre quadrature with parameter alpha_param = a - 1
# provides nodes and weights for exact integration of polynomial * weight.
#
# Author: JoonHo Lee
# Date: December 2025
# Part of: DPprior R Package
# Reference: RN-01 Section 6.2-6.3
# =============================================================================


#' Gauss-Laguerre Quadrature Nodes and Weights
#'
#' Computes generalized Gauss-Laguerre quadrature nodes and weights via
#' eigendecomposition of the Jacobi matrix. These are used to approximate
#' integrals of the form \eqn{\int_0^\infty f(x) x^\beta e^{-x} dx}.
#'
#' @param M Integer; number of quadrature nodes (at least 1, typically 40-120).
#' @param alpha_param Numeric; Laguerre parameter (must be > -1).
#'   For standard Laguerre polynomials, use 0. For integrating against
#'   Gamma(a, b), use \code{alpha_param = a - 1}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{nodes}}{Numeric vector of quadrature nodes \eqn{x_m}.}
#'     \item{\code{weights}}{Numeric vector of quadrature weights \eqn{w_m}.}
#'     \item{\code{weights_log}}{Numeric vector of \eqn{\log(w_m)} for
#'       numerical stability.}
#'   }
#'
#' @details
#' The algorithm constructs a tridiagonal Jacobi matrix \eqn{J} of size
#' \eqn{M \times M}:
#' \itemize{
#'   \item Diagonal: \eqn{a_k = 2k + 1 + \alpha} for \eqn{k = 0, \ldots, M-1}
#'   \item Off-diagonal: \eqn{b_k = \sqrt{k(k + \alpha)}} for \eqn{k = 1, \ldots, M-1}
#' }
#' Eigendecomposition \eqn{J = V D V^T} yields:
#' \itemize{
#'   \item Nodes = eigenvalues (diagonal of \eqn{D})
#'   \item Weights = \eqn{\Gamma(\alpha + 1) \cdot V[1,:]^2}
#' }
#'
#' The generalized Laguerre polynomials \eqn{L_n^{(\alpha)}(x)} are orthogonal
#' with respect to the weight function \eqn{w(x) = x^\alpha e^{-x}} on
#' \eqn{[0, \infty)}.
#'
#' @examples
#' # Standard Laguerre (alpha_param = 0)
#' quad <- gauss_laguerre_nodes(40)
#' sum(quad$weights)  # Should equal Gamma(1) = 1
#'
#' # For Gamma(2.5, b) integration, use alpha_param = 1.5
#' quad <- gauss_laguerre_nodes(80, alpha_param = 1.5)
#' sum(quad$weights)  # Should equal Gamma(2.5)
#'
#' @references
#' Golub, G. H., & Welsch, J. H. (1969). Calculation of Gauss Quadrature Rules.
#' \emph{Mathematics of Computation}, 23(106), 221-230.
#'
#' @seealso \code{\link{build_gamma_quadrature}} for Gamma distribution integration,
#'   \code{\link{integrate_gamma}} for high-level expectation computation
#'
#' @export
gauss_laguerre_nodes <- function(M, alpha_param = 0) {
  # Input validation
  if (!is.numeric(M) || length(M) != 1L || !is.finite(M) ||
      M != floor(M) || M < 1L) {
    stop("M must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(alpha_param) || length(alpha_param) != 1L ||
      !is.finite(alpha_param) || alpha_param <= -1) {
    stop("alpha_param must be a finite number > -1", call. = FALSE)
  }

  M <- as.integer(M)

  # Handle M = 1 special case
  if (M == 1L) {
    node <- alpha_param + 1
    weights_norm <- 1.0
    weight_log <- lgamma(alpha_param + 1)
    weight <- exp(weight_log)
    return(list(
      nodes = node,
      weights = weight,
      weights_log = weight_log,
      weights_norm = weights_norm
    ))
  }

  # Construct tridiagonal Jacobi matrix
  # Diagonal elements: a_k = 2k + 1 + alpha for k = 0, ..., M-1
  n <- 0:(M - 1L)
  diag_elem <- 2 * n + alpha_param + 1

  # Off-diagonal elements: b_k = sqrt(k * (k + alpha)) for k = 1, ..., M-1
  k <- 1:(M - 1L)
  off_diag <- sqrt(k * (k + alpha_param))

  # Build symmetric tridiagonal matrix
  J <- diag(diag_elem)
  for (i in seq_len(M - 1L)) {
    J[i, i + 1L] <- off_diag[i]
    J[i + 1L, i] <- off_diag[i]
  }

  # Eigendecomposition (symmetric matrix)
  eig <- eigen(J, symmetric = TRUE)

  # Sort by eigenvalues (nodes)
  ord <- order(eig$values)
  nodes <- eig$values[ord]
  V <- eig$vectors[, ord, drop = FALSE]

  # Normalized weights are V[1,]^2 (sum to 1 up to rounding)
  weights_norm <- V[1L, ]^2

  # Log weights (allow -Inf if some weights_norm are numerically zero)
  weights_log <- lgamma(alpha_param + 1) + suppressWarnings(log(weights_norm))

  # Raw weights can overflow for large alpha_param; kept for completeness
  weights <- suppressWarnings(exp(weights_log))

  list(
    nodes = nodes,
    weights = weights,
    weights_log = weights_log,
    weights_norm = weights_norm
  )
}


#' Build Quadrature for Gamma(a, b) Integration
#'
#' Transforms standard Gauss-Laguerre quadrature for integration against
#' a Gamma(a, b) distribution with shape \code{a} and rate \code{b}.
#'
#' @param a Numeric; shape parameter of Gamma distribution (must be > 0).
#' @param b Numeric; rate parameter of Gamma distribution (must be > 0).
#' @param M Integer; number of quadrature nodes (default: 80).
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{a}}{Shape parameter.}
#'     \item{\code{b}}{Rate parameter.}
#'     \item{\code{alpha_nodes}}{Numeric vector of transformed nodes
#'       \eqn{\alpha_m = x_m / b} on the \eqn{\alpha} scale.}
#'     \item{\code{weights_normalized}}{Numeric vector of normalized weights
#'       that sum to 1.}
#'   }
#'
#' @details
#' For \eqn{\alpha \sim \text{Gamma}(a, b)}:
#' \deqn{E[g(\alpha)] = \frac{1}{\Gamma(a)} \int_0^\infty g(x/b) x^{a-1} e^{-x} dx}
#'
#' Using generalized Laguerre quadrature with parameter \eqn{\alpha_{\text{param}} = a - 1}:
#' \deqn{E[g(\alpha)] \approx \sum_{m=1}^M \tilde{w}_m g(\alpha_m)}
#'
#' where \eqn{\alpha_m = x_m / b} and \eqn{\tilde{w}_m} are normalized weights
#' summing to 1.
#'
#' The weights are normalized in log-space for numerical stability, which
#' ensures monotone convergence as M increases.
#'
#' @examples
#' # Build quadrature for Gamma(2.5, 1.5)
#' quad <- build_gamma_quadrature(2.5, 1.5)
#'
#' # Check weights sum to 1
#' sum(quad$weights_normalized)
#'
#' # Check mean: E[alpha] should be a/b = 2.5/1.5
#' sum(quad$weights_normalized * quad$alpha_nodes)
#'
#' @seealso \code{\link{gauss_laguerre_nodes}} for raw quadrature computation,
#'   \code{\link{integrate_gamma}} for high-level expectation computation
#'
#' @export
build_gamma_quadrature <- function(a, b, M = .QUAD_NODES_DEFAULT) {
  # Input validation
  assert_positive(a, "a")
  assert_positive(b, "b")

  if (!is.numeric(M) || length(M) != 1L || !is.finite(M) ||
      M != floor(M) || M < 1L) {
    stop("M must be a positive integer", call. = FALSE)
  }

  M <- as.integer(M)

  # Get standard Laguerre nodes and weights with parameter a - 1
  quad <- gauss_laguerre_nodes(M, alpha_param = a - 1)

  # Transform nodes to alpha scale
  alpha_nodes <- quad$nodes / b

  # Normalize weights in log-space for numerical stability
  # This ensures monotone convergence as M increases
  logw <- quad$weights_log
  logw_norm <- logw - logsumexp_vec(logw)
  weights_normalized <- exp(logw_norm)

  list(
    a = a,
    b = b,
    alpha_nodes = alpha_nodes,
    weights_normalized = weights_normalized
  )
}


#' Integrate Function Against Gamma Distribution
#'
#' High-level interface for computing \eqn{E[f(\alpha)]} where
#' \eqn{\alpha \sim \text{Gamma}(a, b)}.
#'
#' @param f Function to integrate; must accept a single numeric argument
#'   and return a single numeric value.
#' @param a Numeric; shape parameter of Gamma distribution (must be > 0).
#' @param b Numeric; rate parameter of Gamma distribution (must be > 0).
#' @param M Integer; number of quadrature nodes (default: 80).
#'
#' @return Numeric; approximation of \eqn{E[f(\alpha)]}.
#'
#' @details
#' Uses Gauss-Laguerre quadrature to approximate:
#' \deqn{E[f(\alpha)] = \int_0^\infty f(\alpha) \frac{b^a}{\Gamma(a)} \alpha^{a-1} e^{-b\alpha} d\alpha}
#'
#' The approximation is:
#' \deqn{E[f(\alpha)] \approx \sum_{m=1}^M w_m f(\alpha_m)}
#'
#' where \eqn{\alpha_m} are quadrature nodes and \eqn{w_m} are normalized weights.
#'
#' For polynomial integrands of degree up to \eqn{2M - 1}, the quadrature is exact.
#' For other smooth functions, accuracy improves rapidly with \eqn{M}.
#'
#' @examples
#' # E[alpha] for Gamma(2.5, 1.5) should be 2.5/1.5
#' integrate_gamma(identity, 2.5, 1.5)
#'
#' # E[alpha^2] for Gamma(2.5, 1.5) should be 2.5*3.5/1.5^2
#' integrate_gamma(function(x) x^2, 2.5, 1.5)
#'
#' # More complex function
#' integrate_gamma(function(x) log(x + 1), 2.5, 1.5)
#'
#' @seealso \code{\link{build_gamma_quadrature}} for the underlying quadrature
#'
#' @export
integrate_gamma <- function(f, a, b, M = .QUAD_NODES_DEFAULT) {
  # Build quadrature
  quad <- build_gamma_quadrature(a, b, M)

  # Evaluate function at quadrature nodes
  f_vals <- sapply(quad$alpha_nodes, f)

  # Weighted sum
  sum(quad$weights_normalized * f_vals)
}


#' Verify Quadrature Accuracy Against Known Gamma Moments
#'
#' Validates quadrature implementation by comparing computed expectations
#' against known closed-form Gamma distribution moments.
#'
#' @param a Numeric; shape parameter of Gamma distribution.
#' @param b Numeric; rate parameter of Gamma distribution.
#' @param M Integer; number of quadrature nodes.
#' @param tol Numeric; tolerance for verification (default: 1e-10).
#' @param verbose Logical; if \code{TRUE}, print detailed results.
#'
#' @return Logical; \code{TRUE} if all moments match within tolerance.
#'
#' @details
#' For \eqn{\alpha \sim \text{Gamma}(a, b)}:
#' \itemize{
#'   \item \eqn{E[\alpha] = a/b}
#'   \item \eqn{E[\alpha^2] = a(a+1)/b^2}
#'   \item \eqn{Var(\alpha) = a/b^2}
#' }
#'
#' @examples
#' # Should return TRUE
#' verify_quadrature(2.5, 1.5, M = 80)
#'
#' # More challenging case
#' verify_quadrature(0.5, 2.0, M = 100, verbose = TRUE)
#'
#' @export
verify_quadrature <- function(a, b, M, tol = 1e-10, verbose = TRUE) {
  assert_positive(a, "a")
  assert_positive(b, "b")

  # Theoretical moments
  E_alpha_true <- a / b
  E_alpha_sq_true <- a * (a + 1) / b^2
  Var_alpha_true <- a / b^2

  # Computed via quadrature
  E_alpha_quad <- integrate_gamma(identity, a, b, M)
  E_alpha_sq_quad <- integrate_gamma(function(x) x^2, a, b, M)
  Var_alpha_quad <- E_alpha_sq_quad - E_alpha_quad^2

  # Errors
  err_mean <- abs(E_alpha_quad - E_alpha_true)
  err_sq <- abs(E_alpha_sq_quad - E_alpha_sq_true)
  err_var <- abs(Var_alpha_quad - Var_alpha_true)

  pass_mean <- err_mean < tol
  pass_sq <- err_sq < tol
  pass_var <- err_var < tol
  all_pass <- pass_mean && pass_sq && pass_var

  if (verbose) {
    cat(sprintf("Quadrature verification (a=%.2f, b=%.2f, M=%d):\n", a, b, M))
    cat(sprintf("  E[alpha]:   true=%.6f, quad=%.10f, error=%.2e [%s]\n",
                E_alpha_true, E_alpha_quad, err_mean,
                if (pass_mean) "PASS" else "FAIL"))
    cat(sprintf("  E[alpha^2]:  true=%.6f, quad=%.10f, error=%.2e [%s]\n",
                E_alpha_sq_true, E_alpha_sq_quad, err_sq,
                if (pass_sq) "PASS" else "FAIL"))
    cat(sprintf("  Var(alpha): true=%.6f, quad=%.10f, error=%.2e [%s]\n",
                Var_alpha_true, Var_alpha_quad, err_var,
                if (pass_var) "PASS" else "FAIL"))
    cat(sprintf("  Overall: %s\n", if (all_pass) "PASS" else "FAIL"))
  }

  invisible(all_pass)
}


#' Get Quadrature Information Summary
#'
#' Returns summary information about the quadrature nodes and weights
#' for a given Gamma(a, b) distribution.
#'
#' @param a Numeric; shape parameter of Gamma distribution.
#' @param b Numeric; rate parameter of Gamma distribution.
#' @param M Integer; number of quadrature nodes.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{n_nodes}}{Number of quadrature nodes.}
#'     \item{\code{alpha_range}}{Range of alpha nodes (min, max).}
#'     \item{\code{weight_range}}{Range of normalized weights (min, max).}
#'     \item{\code{gamma_mean}}{Theoretical mean of Gamma(a, b).}
#'     \item{\code{gamma_sd}}{Theoretical SD of Gamma(a, b).}
#'     \item{\code{coverage}}{Approximate coverage in terms of SD from mean.}
#'   }
#'
#' @examples
#' summary_quadrature(2.5, 1.5, M = 80)
#'
#' @keywords internal
#' @export
summary_quadrature <- function(a, b, M = .QUAD_NODES_DEFAULT) {
  quad <- build_gamma_quadrature(a, b, M)

  gamma_mean <- a / b
  gamma_sd <- sqrt(a) / b

  list(
    n_nodes = M,
    alpha_range = c(min = min(quad$alpha_nodes), max = max(quad$alpha_nodes)),
    weight_range = c(min = min(quad$weights_normalized),
                     max = max(quad$weights_normalized)),
    gamma_mean = gamma_mean,
    gamma_sd = gamma_sd,
    coverage = c(
      lower_sd = (gamma_mean - min(quad$alpha_nodes)) / gamma_sd,
      upper_sd = (max(quad$alpha_nodes) - gamma_mean) / gamma_sd
    )
  )
}


#' Convergence Diagnostic for Quadrature
#'
#' Examines how quadrature accuracy improves with increasing number of nodes.
#'
#' @param f Function to integrate.
#' @param a Numeric; shape parameter of Gamma distribution.
#' @param b Numeric; rate parameter of Gamma distribution.
#' @param M_values Integer vector; number of nodes to test.
#' @param true_value Numeric; known true value (optional).
#'
#' @return A data frame with columns: M, estimate, change, relative_change.
#'
#' @examples
#' # Check convergence for E[alpha]
#' convergence_quadrature(identity, 2.5, 1.5,
#'                        M_values = c(10, 20, 50, 80, 100),
#'                        true_value = 2.5/1.5)
#'
#' @keywords internal
#' @export
convergence_quadrature <- function(f, a, b, M_values = c(10, 20, 50, 80, 100),
                                   true_value = NULL) {
  estimates <- sapply(M_values, function(M) integrate_gamma(f, a, b, M))

  result <- data.frame(
    M = M_values,
    estimate = estimates,
    change = c(NA, diff(estimates)),
    relative_change = c(NA, abs(diff(estimates)) / abs(estimates[-length(estimates)]))
  )

  if (!is.null(true_value)) {
    result$error <- abs(estimates - true_value)
    result$relative_error <- result$error / abs(true_value)
  }

  result
}
