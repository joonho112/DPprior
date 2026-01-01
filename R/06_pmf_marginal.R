# =============================================================================
# Module 06: Marginal PMF of K_J under alpha ~ Gamma(a, b)
# =============================================================================
#
# This module computes the marginal distribution of the number of occupied
# clusters K_J induced by a Dirichlet process (DP) when the concentration
# parameter alpha has a Gamma(a, b) hyperprior (shape-rate).
#
# The key identity is a mixture:
#   P(K_J = k | a, b) = integral P(K_J = k | alpha) g_{a,b}(alpha) d alpha
#                     approx sum_{m=1}^M w_m P(K_J = k | alpha_m),
#
# where (alpha_m, w_m) are Gauss-Laguerre quadrature nodes/weights adapted to
# the Gamma(a, b) distribution (Module 02).
#
# Key Implementation Notes:
# -------------------------
# 1. All mixing is performed in LOG-SPACE for numerical stability. This is
#    critical for large J or extreme alpha values where probabilities can
#    underflow in linear space.
#
# 2. The conditional PMF at each quadrature node is normalized before mixing
#    to ensure proper probability distributions.
#
# 3. A defensive final normalization is applied after mixing.
#
# Important Warning:
# ------------------
# Some older drafts of RN-01 report E[K_J] ~ 10.23 for J=50 and
# alpha ~ Gamma(1.5, 0.5) under the shape-rate parameterization.
# Direct numerical integration of the conditional mean returns
# E[K_50] ~ 8.3555. If you see a mismatch with any hard-coded RN-01
# numbers, treat those as potentially stale and rely on the internal
# consistency checks (PMF <-> moments) implemented in Modules 05-06.
#
# Author: JoonHo Lee
# Date: December 2025
# Part of: DPprior R Package
# Reference: RN-01 Section 6.5
# Dependencies: Module 00 (constants, logsumexp),
#               Module 02 (quadrature),
#               Module 04 (conditional PMF)
# =============================================================================


# =============================================================================
# Internal: Log-Space Marginal PMF (Core Algorithm)
# =============================================================================

#' Log Marginal PMF of K_J under Gamma Hyperprior
#'
#' Computes \eqn{\log P(K_J = k \mid a, b)} for \eqn{k = 0, 1, \ldots, J}
#' using log-space mixing for numerical stability.
#'
#' @param J Integer; sample size (must be >= 1).
#' @param a Numeric; shape parameter of Gamma prior (must be > 0).
#' @param b Numeric; rate parameter of Gamma prior (must be > 0).
#' @param logS Matrix; pre-computed log-Stirling matrix from
#'   \code{\link{compute_log_stirling}}.
#' @param M Integer; number of quadrature nodes (default: \code{.QUAD_NODES_DEFAULT}).
#'
#' @return Numeric vector of length \eqn{J+1} containing log-probabilities
#'   for \eqn{k = 0, 1, \ldots, J}. Entry \code{[1]} corresponds to \eqn{k=0}
#'   and is always \code{-Inf}.
#'
#' @details
#' This routine normalizes \eqn{P(K_J = \cdot \mid \alpha_m)} at each quadrature
#' node before mixing, and then mixes in log-space via \code{logsumexp_vec}:
#' \deqn{\log p_k \approx \log\sum_m \exp\{\log w_m + \log p_{k\mid m}\}.}
#'
#' The log-space computation is essential for numerical stability when:
#' \itemize{
#'   \item J is large (tail probabilities become very small)
#'   \item Alpha values span a wide range (extreme quadrature nodes)
#'   \item Parameters lead to concentrated distributions
#' }
#'
#' @keywords internal
log_pmf_K_marginal <- function(J, a, b, logS, M = .QUAD_NODES_DEFAULT) {
  # Input validation
  assert_valid_J(J)
  assert_positive(a, "a")
  assert_positive(b, "b")
  .check_logS_size(J, logS)

  if (!is.numeric(M) || length(M) != 1L || !is.finite(M) ||
      M != floor(M) || M < 1L) {
    stop("M must be a positive integer", call. = FALSE)
  }

  J <- as.integer(J)
  M <- as.integer(M)

  # Build quadrature for Gamma(a, b)
  quad <- build_gamma_quadrature(a, b, M)
  alpha_nodes <- quad$alpha_nodes
  w <- quad$weights_normalized
  logw <- log(w)

  # Matrix to store normalized log-PMF at each quadrature node
  # logpmf_mat[m, k+1] = log P(K_J = k | alpha_m), normalized over k=1..J
  logpmf_mat <- matrix(-Inf, nrow = M, ncol = J + 1L)

  for (m in seq_len(M)) {
    # Get raw log-PMF from Antoniak distribution
    lp_raw <- log_pmf_K_given_alpha(J, alpha_nodes[m], logS)

    # Normalize across k=1..J (exclude k=0 which is impossible)
    # This ensures each conditional distribution is proper before mixing
    logZ <- logsumexp_vec(lp_raw[-1L])
    logpmf_mat[m, ] <- lp_raw - logZ
  }

  # Mix across quadrature nodes in log-space
  # For each k: log p_k = logsumexp_m(log w_m + log p_{k|m})
  logp <- rep(-Inf, J + 1L)
  for (k_idx in seq_len(J + 1L)) {
    logp[k_idx] <- logsumexp_vec(logw + logpmf_mat[, k_idx])
  }

  # Enforce support constraint: P(K=0) = 0 exactly
  logp[1L] <- -Inf

  # Defensive normalization across k=1..J
  logZ2 <- logsumexp_vec(logp[-1L])
  logp[-1L] <- logp[-1L] - logZ2


  logp
}


# =============================================================================
# Core Exported Functions
# =============================================================================

#' Marginal PMF of K_J under Gamma Hyperprior
#'
#' Computes \eqn{P(K_J = k \mid a, b)} for \eqn{k = 0, 1, \ldots, J} when
#' \eqn{\alpha \sim \mathrm{Gamma}(a, b)} (shape-rate parameterization).
#'
#' @param J Integer; sample size (positive integer >= 1).
#' @param a Numeric; shape parameter of Gamma prior (> 0).
#' @param b Numeric; rate parameter of Gamma prior (> 0).
#' @param logS Matrix; pre-computed log-Stirling matrix from
#'   \code{\link{compute_log_stirling}}.
#' @param M Integer; number of quadrature nodes (default: 80).
#'
#' @return Numeric vector of length \eqn{J+1} containing
#'   \eqn{P(K_J = k \mid a, b)} for \eqn{k = 0, 1, \ldots, J}.
#'   Entry \code{[1]} corresponds to \eqn{k=0} and always equals 0.
#'   The vector sums to 1.
#'
#' @details
#' Uses Gauss-Laguerre quadrature to numerically evaluate:
#' \deqn{P(K_J = k \mid a, b) = \int_0^\infty P(K_J = k \mid \alpha) \cdot g_{a,b}(\alpha) d\alpha}
#' \deqn{\approx \sum_{m=1}^M \tilde{w}_m \cdot P(K_J = k \mid \alpha_m)}
#'
#' where \eqn{P(K_J = k \mid \alpha)} is the Antoniak distribution from Module 04
#' and \eqn{(\alpha_m, \tilde{w}_m)} are the transformed quadrature nodes and
#' normalized weights from Module 02.
#'
#' \strong{Implementation:} All mixing is performed in log-space for numerical
#' stability. This is critical for large J or extreme parameter values.
#'
#' \strong{Key properties:}
#' \itemize{
#'   \item \eqn{P(K_J = 0) = 0} always (at least one cluster exists)
#'   \item The PMF sums to 1
#'   \item Moments from the PMF match \code{exact_K_moments()} within numerical tolerance
#'   \item Mode is typically near \eqn{E[K_J]} but may differ
#' }
#'
#' @examples
#' # Pre-compute Stirling numbers
#' logS <- compute_log_stirling(50)
#'
#' # Compute marginal PMF for J=50, Gamma(1.5, 0.5) prior
#' pmf <- pmf_K_marginal(50, 1.5, 0.5, logS)
#'
#' # Verify normalization
#' sum(pmf)
#'
#' # Most likely number of clusters
#' which.max(pmf) - 1
#'
#' # Compare mean with exact_K_moments
#' k_vals <- 0:50
#' mean_pmf <- sum(k_vals * pmf)
#' exact <- exact_K_moments(50, 1.5, 0.5)
#' abs(mean_pmf - exact$mean)
#'
#' @seealso \code{\link{log_pmf_K_marginal}} for log-scale computation,
#'   \code{\link{pmf_K_given_alpha}} for conditional PMF,
#'   \code{\link{exact_K_moments}} for marginal moments,
#'   \code{\link{cdf_K_marginal}}, \code{\link{quantile_K_marginal}},
#'   \code{\link{mode_K_marginal}}, \code{\link{summary_K_marginal}}
#'
#' @references
#' Antoniak, C. E. (1974). Mixtures of Dirichlet Processes with Applications
#' to Bayesian Nonparametric Problems. \emph{The Annals of Statistics},
#' 2(6), 1152-1174.
#'
#' @export
pmf_K_marginal <- function(J, a, b, logS, M = .QUAD_NODES_DEFAULT) {
  # Compute in log-space for numerical stability
  logp <- log_pmf_K_marginal(J, a, b, logS, M)

  # Convert to linear scale
  pmf <- exp(logp)

  # Enforce P(K=0) = 0 exactly and renormalize
  pmf[1L] <- 0.0
  s <- sum(pmf)
  if (is.finite(s) && s > 0) {
    pmf <- pmf / s
  }

  pmf
}


#' CDF of Marginal K Distribution
#'
#' Computes the cumulative distribution function \eqn{P(K_J \leq k \mid a, b)}
#' for \eqn{k = 0, 1, \ldots, J}.
#'
#' @param J Integer; sample size.
#' @param a Numeric; shape parameter of Gamma prior (> 0).
#' @param b Numeric; rate parameter of Gamma prior (> 0).
#' @param logS Matrix; pre-computed log-Stirling matrix.
#' @param M Integer; number of quadrature nodes (default: 80).
#'
#' @return Numeric vector of length \eqn{J+1} containing
#'   \eqn{P(K_J \leq k \mid a, b)} for \eqn{k = 0, 1, \ldots, J}.
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
#' cdf <- cdf_K_marginal(50, 1.5, 0.5, logS)
#'
#' # Verify CDF ends at 1
#' cdf[51]
#'
#' # P(K <= 10)
#' cdf[11]
#'
#' @seealso \code{\link{pmf_K_marginal}}, \code{\link{quantile_K_marginal}}
#'
#' @export
cdf_K_marginal <- function(J, a, b, logS, M = .QUAD_NODES_DEFAULT) {
  pmf <- pmf_K_marginal(J, a, b, logS, M)
  cumsum(pmf)
}


#' Quantile of Marginal K Distribution
#'
#' Computes the \eqn{p}-th quantile of the marginal distribution of \eqn{K_J}.
#'
#' @param p Numeric; probability level(s) in \eqn{[0, 1]}. Can be scalar or vector.
#' @param J Integer; sample size.
#' @param a Numeric; shape parameter of Gamma prior (> 0).
#' @param b Numeric; rate parameter of Gamma prior (> 0).
#' @param logS Matrix; pre-computed log-Stirling matrix.
#' @param M Integer; number of quadrature nodes (default: 80).
#'
#' @return Integer vector of quantiles (same length as \code{p}).
#'   Each element is the smallest \eqn{k} such that \eqn{P(K_J \leq k) \geq p}.
#'
#' @details
#' This is the standard quantile definition for discrete distributions:
#' \eqn{Q(p) = \min\{k : F(k) \geq p\}}.
#'
#' The function is vectorized over \code{p}, allowing efficient computation
#' of multiple quantiles in a single call.
#'
#' @examples
#' logS <- compute_log_stirling(50)
#'
#' # Single quantile (median)
#' quantile_K_marginal(0.5, 50, 1.5, 0.5, logS)
#'
#' # Multiple quantiles at once
#' quantile_K_marginal(c(0.1, 0.25, 0.5, 0.75, 0.9), 50, 1.5, 0.5, logS)
#'
#' # Interquartile range
#' qs <- quantile_K_marginal(c(0.25, 0.75), 50, 1.5, 0.5, logS)
#' diff(qs)
#'
#' @seealso \code{\link{cdf_K_marginal}}, \code{\link{pmf_K_marginal}}
#'
#' @export
quantile_K_marginal <- function(p, J, a, b, logS, M = .QUAD_NODES_DEFAULT) {
  # Input validation
  assert_probability(p, "p")

  J <- as.integer(J)
  cdf <- cdf_K_marginal(J, a, b, logS, M)

  # Vectorized quantile computation
  as.integer(sapply(p, function(pi) {
    if (pi <= 0) return(0L)
    if (pi >= 1) return(J)
    # Find smallest k such that CDF[k] >= pi
    # Note: cdf[k+1] = P(K <= k) due to 0-indexing
    idx <- which(cdf >= pi)
    if (length(idx) == 0L) return(J)
    min(idx) - 1L
  }))
}


#' Mode of Marginal K Distribution
#'
#' Computes the mode (most likely value) of the marginal distribution of \eqn{K_J}.
#'
#' @param J Integer; sample size.
#' @param a Numeric; shape parameter of Gamma prior (> 0).
#' @param b Numeric; rate parameter of Gamma prior (> 0).
#' @param logS Matrix; pre-computed log-Stirling matrix.
#' @param M Integer; number of quadrature nodes (default: 80).
#'
#' @return Integer; the value \eqn{k} that maximizes \eqn{P(K_J = k \mid a, b)}.
#'
#' @details
#' The mode is always >= 1 since \eqn{P(K_J = 0) = 0}.
#'
#' @examples
#' logS <- compute_log_stirling(50)
#' mode_K_marginal(50, 1.5, 0.5, logS)
#'
#' @seealso \code{\link{pmf_K_marginal}}, \code{\link{summary_K_marginal}}
#'
#' @export
mode_K_marginal <- function(J, a, b, logS, M = .QUAD_NODES_DEFAULT) {
  pmf <- pmf_K_marginal(J, a, b, logS, M)
  as.integer(which.max(pmf) - 1L)
}


#' Summary Statistics for Marginal K Distribution
#'
#' Computes comprehensive summary statistics for the marginal distribution
#' of \eqn{K_J} under a Gamma prior on \eqn{\alpha}.
#'
#' @param J Integer; sample size (positive integer >= 1).
#' @param a Numeric; shape parameter of Gamma prior (> 0).
#' @param b Numeric; rate parameter of Gamma prior (> 0).
#' @param logS Matrix; pre-computed log-Stirling matrix.
#' @param M Integer; number of quadrature nodes (default: 80).
#' @param probs Numeric vector; probability levels for quantiles (default:
#'   \code{c(0.05, 0.25, 0.5, 0.75, 0.95)}).
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{J}}{Sample size}
#'     \item{\code{a}}{Gamma shape parameter}
#'     \item{\code{b}}{Gamma rate parameter}
#'     \item{\code{mean}}{Mean \eqn{E[K_J \mid a, b]}}
#'     \item{\code{var}}{Variance \eqn{Var(K_J \mid a, b)}}
#'     \item{\code{sd}}{Standard deviation}
#'     \item{\code{cv}}{Coefficient of variation (sd/mean)}
#'     \item{\code{mode}}{Mode (most likely value)}
#'     \item{\code{median}}{Median (50th percentile)}
#'     \item{\code{quantiles}}{Named integer vector of quantiles at \code{probs}}
#'     \item{\code{pmf}}{Full PMF vector}
#'     \item{\code{cdf}}{Full CDF vector}
#'   }
#'
#' @details
#' This function provides a complete summary of the marginal distribution,
#' combining PMF-based and CDF-based statistics. The mean and variance
#' computed from the PMF should match \code{exact_K_moments()} within
#' numerical tolerance.
#'
#' The \code{probs} argument allows customization of which quantiles to
#' report, making it flexible for different reporting needs.
#'
#' @examples
#' logS <- compute_log_stirling(50)
#' summary <- summary_K_marginal(50, 1.5, 0.5, logS)
#'
#' # View main statistics
#' summary$mean
#' summary$var
#' summary$mode
#' summary$quantiles
#'
#' # Custom quantiles
#' summary2 <- summary_K_marginal(50, 1.5, 0.5, logS,
#'                                probs = c(0.025, 0.5, 0.975))
#' summary2$quantiles
#'
#' # Compare with exact moments
#' exact <- exact_K_moments(50, 1.5, 0.5)
#' c(summary$mean - exact$mean, summary$var - exact$var)
#'
#' @seealso \code{\link{pmf_K_marginal}}, \code{\link{exact_K_moments}}
#'
#' @export
summary_K_marginal <- function(J, a, b, logS, M = .QUAD_NODES_DEFAULT,
                               probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) {
  # Input validation
  assert_probability(probs, "probs")

  # Compute PMF and CDF
  pmf <- pmf_K_marginal(J, a, b, logS, M)
  cdf <- cumsum(pmf)

  k_vals <- 0:as.integer(J)

  # Mean and variance from PMF
  mean_K <- sum(k_vals * pmf)
  var_K <- sum(k_vals^2 * pmf) - mean_K^2
  var_K <- max(0, var_K)  # Numerical safety
  sd_K <- sqrt(var_K)
  cv_K <- if (mean_K > 0) sd_K / mean_K else Inf

  # Mode
  mode_K <- as.integer(which.max(pmf) - 1L)

  # Quantiles (vectorized)
  qs <- quantile_K_marginal(probs, J, a, b, logS, M)
  names(qs) <- paste0("q", formatC(100 * probs, format = "f", digits = 0))

  # Median (always include even if not in probs)
  median_K <- as.integer(min(which(cdf >= 0.5)) - 1L)

  list(
    J = J,
    a = a,
    b = b,
    mean = as.numeric(mean_K),
    var = as.numeric(var_K),
    sd = sd_K,
    cv = cv_K,
    mode = mode_K,
    median = median_K,
    quantiles = qs,
    pmf = pmf,
    cdf = cdf
  )
}


# =============================================================================
# Convenience Functions
# =============================================================================

#' Mean of Marginal K from PMF
#'
#' Computes \eqn{E[K_J \mid a, b]} from the marginal PMF.
#'
#' @param J Integer; sample size.
#' @param a Numeric; shape parameter of Gamma prior.
#' @param b Numeric; rate parameter of Gamma prior.
#' @param logS Matrix; pre-computed log-Stirling matrix.
#' @param M Integer; number of quadrature nodes.
#'
#' @return Numeric; marginal mean.
#'
#' @examples
#' logS <- compute_log_stirling(50)
#' mean_K_from_marginal_pmf(50, 1.5, 0.5, logS)
#'
#' @keywords internal
#' @export
mean_K_from_marginal_pmf <- function(J, a, b, logS, M = .QUAD_NODES_DEFAULT) {
  pmf <- pmf_K_marginal(J, a, b, logS, M)
  sum((0:J) * pmf)
}


#' Variance of Marginal K from PMF
#'
#' Computes \eqn{Var(K_J \mid a, b)} from the marginal PMF.
#'
#' @param J Integer; sample size.
#' @param a Numeric; shape parameter of Gamma prior.
#' @param b Numeric; rate parameter of Gamma prior.
#' @param logS Matrix; pre-computed log-Stirling matrix.
#' @param M Integer; number of quadrature nodes.
#'
#' @return Numeric; marginal variance.
#'
#' @examples
#' logS <- compute_log_stirling(50)
#' var_K_from_marginal_pmf(50, 1.5, 0.5, logS)
#'
#' @keywords internal
#' @export
var_K_from_marginal_pmf <- function(J, a, b, logS, M = .QUAD_NODES_DEFAULT) {
  pmf <- pmf_K_marginal(J, a, b, logS, M)
  k_vals <- 0:J
  mean_K <- sum(k_vals * pmf)
  max(0, sum(k_vals^2 * pmf) - mean_K^2)
}


# =============================================================================
# Verification Functions
# =============================================================================

#' Verify Marginal PMF Properties
#'
#' Verifies that the marginal PMF satisfies basic probability properties.
#'
#' @param J Integer; sample size.
#' @param a Numeric; shape parameter of Gamma prior.
#' @param b Numeric; rate parameter of Gamma prior.
#' @param logS Matrix; pre-computed log-Stirling matrix.
#' @param M Integer; number of quadrature nodes.
#' @param tol Numeric; tolerance for comparisons.
#' @param verbose Logical; if TRUE, print results.
#'
#' @return Logical; TRUE if all verifications pass.
#'
#' @examples
#' logS <- compute_log_stirling(50)
#' verify_pmf_marginal_properties(50, 1.5, 0.5, logS)
#'
#' @export
verify_pmf_marginal_properties <- function(J, a, b, logS, M = .QUAD_NODES_DEFAULT,
                                           tol = 1e-10, verbose = TRUE) {
  pmf <- pmf_K_marginal(J, a, b, logS, M)

  # Test 1: Sum to 1
  sum_test <- abs(sum(pmf) - 1.0) < tol

  # Test 2: P(K=0) = 0
  zero_test <- pmf[1L] < tol

  # Test 3: All non-negative
  nonneg_test <- all(pmf >= -tol)

  # Test 4: CDF is monotonic
  cdf <- cumsum(pmf)
  monotone_test <- all(diff(cdf) >= -tol)

  # Test 5: CDF ends at 1
  cdf_end_test <- abs(cdf[J + 1L] - 1.0) < tol

  all_pass <- sum_test && zero_test && nonneg_test && monotone_test && cdf_end_test

  if (isTRUE(verbose)) {
    cat(sprintf("Marginal PMF Properties (J=%d, a=%.2f, b=%.2f):\n", J, a, b))
    cat(sprintf("  Sum = 1:          %s (sum = %.12f)\n",
                if (sum_test) "PASS" else "FAIL", sum(pmf)))
    cat(sprintf("  P(K=0) = 0:       %s (P(K=0) = %.2e)\n",
                if (zero_test) "PASS" else "FAIL", pmf[1L]))
    cat(sprintf("  Non-negative:     %s (min = %.2e)\n",
                if (nonneg_test) "PASS" else "FAIL", min(pmf)))
    cat(sprintf("  CDF monotonic:    %s\n",
                if (monotone_test) "PASS" else "FAIL"))
    cat(sprintf("  CDF[J] = 1:       %s (CDF[J] = %.12f)\n",
                if (cdf_end_test) "PASS" else "FAIL", cdf[J + 1L]))
    cat(sprintf("  Overall: %s\n", if (all_pass) "PASS" else "FAIL"))
  }

  invisible(all_pass)
}


#' Verify Moments Consistency
#'
#' Verifies that moments computed from the marginal PMF match those from
#' \code{exact_K_moments()}.
#'
#' @param J Integer; sample size.
#' @param a Numeric; shape parameter of Gamma prior.
#' @param b Numeric; rate parameter of Gamma prior.
#' @param logS Matrix; pre-computed log-Stirling matrix.
#' @param M Integer; number of quadrature nodes.
#' @param tol Numeric; tolerance for comparisons.
#' @param verbose Logical; if TRUE, print results.
#'
#' @return Logical; TRUE if moments match within tolerance.
#'
#' @examples
#' logS <- compute_log_stirling(50)
#' verify_pmf_marginal_moments(50, 1.5, 0.5, logS)
#'
#' @export
verify_pmf_marginal_moments <- function(J, a, b, logS, M = .QUAD_NODES_DEFAULT,
                                        tol = 1e-6, verbose = TRUE) {
  # Moments from PMF
  summary <- summary_K_marginal(J, a, b, logS, M)

  # Moments from quadrature (exact method)
  exact <- exact_K_moments(J, a, b, M)

  # Errors
  mean_err <- abs(summary$mean - exact$mean)
  var_err <- abs(summary$var - exact$var)

  pass_mean <- mean_err < tol
  pass_var <- var_err < tol
  all_pass <- pass_mean && pass_var

  if (isTRUE(verbose)) {
    cat(sprintf("Moments Consistency (J=%d, a=%.2f, b=%.2f):\n", J, a, b))
    cat(sprintf("  Mean (PMF):       %.8f\n", summary$mean))
    cat(sprintf("  Mean (quadrature):%.8f\n", exact$mean))
    cat(sprintf("  Mean error:       %.2e [%s]\n",
                mean_err, if (pass_mean) "PASS" else "FAIL"))
    cat(sprintf("  Var (PMF):        %.8f\n", summary$var))
    cat(sprintf("  Var (quadrature): %.8f\n", exact$var))
    cat(sprintf("  Var error:        %.2e [%s]\n",
                var_err, if (pass_var) "PASS" else "FAIL"))
    cat(sprintf("  Overall: %s\n", if (all_pass) "PASS" else "FAIL"))
  }

  invisible(all_pass)
}


#' Verify Quadrature Convergence
#'
#' Checks that the marginal PMF converges as the number of quadrature
#' nodes increases.
#'
#' @param J Integer; sample size.
#' @param a Numeric; shape parameter of Gamma prior.
#' @param b Numeric; rate parameter of Gamma prior.
#' @param logS Matrix; pre-computed log-Stirling matrix.
#' @param M_values Integer vector; quadrature node counts to test.
#' @param verbose Logical; if TRUE, print results.
#'
#' @return Data frame with convergence results.
#'
#' @examples
#' logS <- compute_log_stirling(50)
#' verify_pmf_marginal_convergence(50, 1.5, 0.5, logS)
#'
#' @export
verify_pmf_marginal_convergence <- function(J, a, b, logS,
                                            M_values = c(20L, 40L, 80L, 120L),
                                            verbose = TRUE) {
  results <- data.frame(
    M = integer(length(M_values)),
    mean = numeric(length(M_values)),
    var = numeric(length(M_values)),
    mean_change = numeric(length(M_values)),
    var_change = numeric(length(M_values)),
    L1_change = numeric(length(M_values))
  )

  prev_pmf <- NULL

  for (i in seq_along(M_values)) {
    M <- M_values[i]
    pmf <- pmf_K_marginal(J, a, b, logS, M)
    summary <- summary_K_marginal(J, a, b, logS, M)

    results$M[i] <- M
    results$mean[i] <- summary$mean
    results$var[i] <- summary$var

    if (i > 1L) {
      results$mean_change[i] <- abs(results$mean[i] - results$mean[i - 1L])
      results$var_change[i] <- abs(results$var[i] - results$var[i - 1L])
      results$L1_change[i] <- sum(abs(pmf - prev_pmf))
    }

    prev_pmf <- pmf
  }

  if (isTRUE(verbose)) {
    cat(sprintf("Quadrature Convergence (J=%d, a=%.2f, b=%.2f):\n", J, a, b))
    cat(sprintf("  %4s %12s %12s %12s %12s %12s\n",
                "M", "mean", "var", "mean_chg", "var_chg", "L1_chg"))
    for (i in seq_len(nrow(results))) {
      cat(sprintf("  %4d %12.8f %12.8f %12.2e %12.2e %12.2e\n",
                  results$M[i], results$mean[i], results$var[i],
                  results$mean_change[i], results$var_change[i],
                  results$L1_change[i]))
    }
  }

  invisible(results)
}


#' Run All Module 06 Verification Tests
#'
#' Comprehensive verification suite for the marginal PMF module.
#'
#' @param verbose Logical; if TRUE, print detailed results.
#'
#' @return Logical; TRUE if all tests pass.
#'
#' @examples
#' verify_pmf_marginal_all()
#'
#' @export
verify_pmf_marginal_all <- function(verbose = TRUE) {
  if (isTRUE(verbose)) {
    cat("=", rep("=", 69), "\n", sep = "")
    cat("Module 06: Marginal PMF - Full Verification Suite\n")
    cat("=", rep("=", 69), "\n\n", sep = "")
  }

  # Pre-compute Stirling numbers
  logS <- compute_log_stirling(100)

  all_pass <- TRUE

  # Test cases
  test_cases <- list(
    list(J = 50, a = 1.5, b = 0.5),
    list(J = 50, a = 2.0, b = 1.0),
    list(J = 100, a = 1.5, b = 0.5),
    list(J = 30, a = 1.0, b = 0.5),
    list(J = 50, a = 3.0, b = 1.5)
  )

  for (tc in test_cases) {
    if (isTRUE(verbose)) {
      cat(sprintf("\n[Test case: J=%d, a=%.2f, b=%.2f]\n",
                  tc$J, tc$a, tc$b))
      cat(strrep("-", 50), "\n")
    }

    pass1 <- verify_pmf_marginal_properties(tc$J, tc$a, tc$b, logS,
                                            verbose = verbose)
    pass2 <- verify_pmf_marginal_moments(tc$J, tc$a, tc$b, logS,
                                         verbose = verbose)

    all_pass <- all_pass && pass1 && pass2
  }

  # Convergence test
  if (isTRUE(verbose)) {
    cat("\n[Convergence test]\n")
    cat(strrep("-", 50), "\n")
  }
  verify_pmf_marginal_convergence(50, 1.5, 0.5, logS, verbose = verbose)

  # Summary
  if (isTRUE(verbose)) {
    cat("\n", strrep("=", 70), "\n", sep = "")
    cat(sprintf("Overall Result: %s\n",
                if (all_pass) "ALL TESTS PASSED" else "SOME TESTS FAILED"))
    cat(strrep("=", 70), "\n")
  }

  invisible(all_pass)
}
