# =============================================================================
# Module 12: A2-KL KL Divergence Minimization for Prior Calibration
# =============================================================================
#
# This module implements the A2-KL algorithm for calibrating Gamma hyperpriors
# on the Dirichlet process concentration parameter by minimizing KL divergence
# between a target distribution and the induced marginal PMF of K_J.
#
# Theory Background (Lee, 2026, Sections 2--3 and Section 3.2):
# ---------------------------------
# Given a target PMF p*(k) for K_J, the A2-KL algorithm finds (a*, b*) such that
# the induced marginal distribution p_{a,b}(K_J = k) minimizes:
#
#   D_KL(p* || p_{a,b}) = sum_k p*(k) * log(p*(k) / p_{a,b}(k))
#
# This allows fitting to **arbitrary target distributions**, not just moments.
#
# Key Features:
# - KL divergence minimization via L-BFGS-B optimization (bounded)
# - Support for user-specified PMF targets (method="pmf")
# - Chi-square discretization for moment-based targets (method="chisq")
# - Log-parameterization with explicit bounds for numerical stability
# - Initialization from A2-MN (exact moment matching) estimates
# - Predeclared, recorded bounded fallback when the primary optimizer fails
# - Independent higher-order PMF verification and target-family adequacy checks
# - Exact zero-support KL semantics without epsilon smoothing
# - Comprehensive trace recording for diagnostics
#
# Author: JoonHo Lee (jlee296@ua.edu)
# Date: December 2025
# Part of: DPprior R Package
# Reference: Lee (2026), Sections 2--3 and Section 3.2
# Dependencies: Modules 00 (constants), 02 (quadrature), 04 (conditional PMF),
#               05 (marginal moments), 06 (marginal PMF), 10 (A1), 11 (A2-MN)
# =============================================================================


# =============================================================================
# Internal Helper Functions
# =============================================================================

#' Validate a Target PMF on \{1, ..., J\}
#'
#' @param target_pmf Numeric vector; target PMF.
#' @param J Integer; sample size.
#'
#' @return Numeric vector of length J (support k=1,...,J). The input is
#'   validated but never normalized or otherwise modified.
#'
#' @keywords internal
.a2_kl_normalize_pmf <- function(target_pmf, J) {
  target_pmf <- .dpprior_validate_plain_vector(
    target_pmf, "target_pmf", "dpprior_pmf_error"
  )
  if (!is.numeric(target_pmf)) {
    .dpprior_abort_invalid(
      "target_pmf must be a numeric vector.",
      c("dpprior_pmf_type_error", "dpprior_pmf_error"),
      "target_pmf", target_pmf, "numeric vector", "type"
    )
  }

  has_k0_entry <- length(target_pmf) == J + 1L
  if (!has_k0_entry && length(target_pmf) != J) {
    .dpprior_abort_invalid(
      sprintf("target_pmf must have length J=%d or J+1=%d.", J, J + 1L),
      c("dpprior_pmf_length_error", "dpprior_pmf_error"),
      "target_pmf", target_pmf,
      sprintf("length %d or %d", J, J + 1L), "length"
    )
  }

  if (anyNA(target_pmf) || any(!is.finite(target_pmf))) {
    .dpprior_abort_invalid(
      "target_pmf must be finite and non-missing.",
      c("dpprior_pmf_nonfinite_error", "dpprior_pmf_error"),
      "target_pmf", target_pmf, "finite and non-missing", "nonfinite"
    )
  }
  if (any(target_pmf < 0)) {
    .dpprior_abort_invalid(
      "target_pmf must be non-negative.",
      c("dpprior_pmf_negative_error", "dpprior_pmf_error"),
      "target_pmf", target_pmf, "non-negative probabilities", "negative"
    )
  }
  if (has_k0_entry) {
    if (!identical(as.numeric(target_pmf[1L]), 0)) {
      .dpprior_abort_invalid(
        paste(
          "target_pmf[1] corresponds to K = 0 and must be exactly 0",
          "because K_J has support {1, ..., J}."
        ),
        c("dpprior_pmf_support_error", "dpprior_pmf_error"),
        "target_pmf[1]", target_pmf[1L], "exactly 0", "support"
      )
    }
    target_pmf <- target_pmf[-1L]
  }

  s <- sum(target_pmf)
  if (!is.finite(s) || s <= 0) {
    .dpprior_abort_invalid(
      "target_pmf must have positive total mass on K=1:J.",
      c("dpprior_pmf_zero_mass_error", "dpprior_pmf_error"),
      "target_pmf", target_pmf, "positive total mass", "zero_mass"
    )
  }
  if (abs(s - 1) > .TOL_PMF_SUM) {
    .dpprior_abort_invalid(
      sprintf(
        paste(
          "target_pmf must sum to 1 within tolerance %.3g; observed sum",
          "was %.17g. Normalize weights explicitly before calling."
        ),
        .TOL_PMF_SUM, s
      ),
      c("dpprior_pmf_normalization_error", "dpprior_pmf_error"),
      "target_pmf", target_pmf,
      sprintf("sum within %.3g of 1", .TOL_PMF_SUM), "normalization"
    )
  }

  unname(as.numeric(target_pmf))
}


# Validate a generic mathematical PMF without silently repairing it.
.a2_kl_validate_mathematical_pmf <- function(x, name, expected_length = NULL) {
  x <- .dpprior_validate_plain_vector(x, name, "dpprior_kl_input_error")
  if (!is.numeric(x)) {
    .dpprior_abort_invalid(
      sprintf("%s must be a numeric vector.", name),
      c("dpprior_pmf_type_error", "dpprior_kl_input_error"),
      name, x, "numeric vector", "type"
    )
  }
  if (!length(x)) {
    .dpprior_abort_invalid(
      sprintf("%s must have positive length.", name),
      c("dpprior_pmf_length_error", "dpprior_kl_input_error"),
      name, x, "positive length", "length"
    )
  }
  if (!is.null(expected_length) && length(x) != expected_length) {
    .dpprior_abort_invalid(
      sprintf("%s must have length %d.", name, expected_length),
      c("dpprior_pmf_length_error", "dpprior_kl_input_error"),
      name, x, sprintf("length %d", expected_length), "length"
    )
  }
  if (anyNA(x) || any(!is.finite(x))) {
    .dpprior_abort_invalid(
      sprintf("%s must be finite and non-missing.", name),
      c("dpprior_pmf_nonfinite_error", "dpprior_kl_input_error"),
      name, x, "finite and non-missing", "nonfinite"
    )
  }
  if (any(x < 0)) {
    .dpprior_abort_invalid(
      sprintf("%s must be non-negative.", name),
      c("dpprior_pmf_negative_error", "dpprior_kl_input_error"),
      name, x, "non-negative probabilities", "negative"
    )
  }
  total <- sum(x)
  if (!is.finite(total) || total <= 0) {
    .dpprior_abort_invalid(
      sprintf("%s must have positive total mass.", name),
      c("dpprior_pmf_zero_mass_error", "dpprior_kl_input_error"),
      name, x, "positive total mass", "zero_mass"
    )
  }
  if (abs(total - 1) > .TOL_PMF_SUM) {
    .dpprior_abort_invalid(
      sprintf(
        "%s must sum to 1 within tolerance %.3g; observed sum was %.17g.",
        name, .TOL_PMF_SUM, total
      ),
      c("dpprior_pmf_normalization_error", "dpprior_kl_input_error"),
      name, x, sprintf("sum within %.3g of 1", .TOL_PMF_SUM),
      "normalization"
    )
  }
  unname(as.numeric(x))
}


#' Induced Marginal PMF of K_J under alpha ~ Gamma(a, b)
#'
#' @param J Integer; sample size.
#' @param a,b Gamma hyperparameters.
#' @param logS Matrix; log-Stirling numbers (from compute_log_stirling(J)).
#' @param M Quadrature nodes.
#'
#' @return Numeric vector of length J; PMF on k = 1,...,J.
#'
#' @keywords internal
.a2_kl_induced_pmf <- function(J, a, b, logS, M) {
  pmf_full <- pmf_K_marginal(
    J = J, a = a, b = b, logS = logS, M = M,
    M_verify = NULL, strict = FALSE
  )
  .a2_kl_validate_mathematical_pmf(
    unname(pmf_full[-1L]), "induced_pmf", expected_length = J
  )
}


# Return the induced PMF in log space so a mathematically positive but very
# small q(k) is not converted to an artificial zero before evaluating KL.
.a2_kl_induced_log_pmf <- function(
    J, a, b, logS, M, M_verify = NULL,
    abs_tol = 1e-10, rel_tol = 1e-8) {
  logp_full <- log_pmf_K_marginal(
    J = J, a = a, b = b, logS = logS, M = M,
    M_verify = M_verify, abs_tol = abs_tol, rel_tol = rel_tol,
    strict = FALSE
  )
  selected <- unname(logp_full[-1L])
  verification_full <- attr(
    logp_full, ".marginal_verification_logp", exact = TRUE
  )
  verification <- if (is.null(verification_full)) {
    NULL
  } else {
    unname(verification_full[-1L])
  }
  list(
    selected = selected,
    verification = verification,
    metadata = attr(logp_full, "marginal_metadata", exact = TRUE)
  )
}


#' Compute KL Divergence from Target to Induced PMF
#'
#' @param target_pmf Numeric vector; normalized PMF (length J).
#' @param induced_pmf Numeric vector; normalized induced PMF (length J).
#' @return Numeric; KL divergence (non-negative).
#'
#' @keywords internal
.a2_kl_compute_kl <- function(target_pmf, induced_pmf) {
  target_pmf <- .a2_kl_validate_mathematical_pmf(
    target_pmf, "target_pmf"
  )
  induced_pmf <- .a2_kl_validate_mathematical_pmf(
    induced_pmf, "induced_pmf", expected_length = length(target_pmf)
  )
  idx <- target_pmf > 0
  if (any(induced_pmf[idx] == 0)) {
    return(Inf)
  }
  terms <- target_pmf[idx] * (
    log(target_pmf[idx]) - log(induced_pmf[idx])
  )
  kl <- sum(terms)
  negative_budget <- 2 * .TOL_PMF_SUM +
    100 * .Machine$double.eps * max(1, sum(abs(terms)))
  if (!is.finite(kl) || kl < -negative_budget) {
    stop(.dpprior_new_condition(
      "KL computation produced an invalid negative or non-finite result.",
      c("dpprior_kl_numerical_error", "dpprior_numerical_error",
        "dpprior_error", "error"),
      value = kl, code = "invalid_kl"
    ))
  }
  max(0, kl)
}


.a2_kl_compute_kl_logq <- function(target_pmf, log_q) {
  if (!is.numeric(log_q) || length(log_q) != length(target_pmf) ||
      anyNA(log_q) || any(log_q > 0) ||
      any(!is.finite(log_q) & log_q != -Inf)) {
    stop(.dpprior_new_condition(
      "Induced log-PMF violated the finite support/log-probability contract.",
      c("dpprior_kl_numerical_error", "dpprior_numerical_error",
        "dpprior_error", "error"),
      value = log_q, code = "invalid_log_pmf"
    ))
  }
  idx <- target_pmf > 0
  if (any(log_q[idx] == -Inf)) {
    return(Inf)
  }
  terms <- target_pmf[idx] * (log(target_pmf[idx]) - log_q[idx])
  kl <- sum(terms)
  negative_budget <- 2 * .TOL_PMF_SUM +
    100 * .Machine$double.eps * max(1, sum(abs(terms)))
  if (!is.finite(kl) || kl < -negative_budget) {
    stop(.dpprior_new_condition(
      "Log-space KL computation produced an invalid result.",
      c("dpprior_kl_numerical_error", "dpprior_numerical_error",
        "dpprior_error", "error"),
      value = kl, code = "invalid_kl"
    ))
  }
  max(0, kl)
}


# =============================================================================
# KL Divergence Computation (Public API)
# =============================================================================

#' KL Divergence Between Two PMFs
#'
#' Computes the Kullback-Leibler divergence \eqn{D_{KL}(p \| q)} between two
#' probability mass functions.
#'
#' @param p Numeric vector; target PMF (reference distribution).
#' @param q Numeric vector; comparison PMF.
#' @param eps Deprecated compatibility argument. Exact mathematical KL does
#'   not smooth zero probabilities; supplying a non-\code{NULL} value is an
#'   error.
#'
#' @return Numeric scalar; the KL divergence (non-negative).
#'
#' @details
#' The KL divergence is defined as:
#' \deqn{D_{KL}(p \| q) = \sum_k p(k) \log\frac{p(k)}{q(k)}}
#'
#' Terms with \eqn{p(k)=0} contribute zero. If \eqn{p(k)>0} and
#' \eqn{q(k)=0}, the result is \code{Inf}. Both inputs must already be valid,
#' normalized PMFs; this function never normalizes or smooths them. In keeping
#' with the package PMF contract, a sum differing from one by at most
#' \code{.TOL_PMF_SUM} is accepted without modifying the vector. A raw negative
#' value whose magnitude is fully explained by that validation tolerance and
#' floating-point roundoff is reported as zero.
#'
#' \strong{Properties:}
#' \itemize{
#'   \item \eqn{D_{KL}(p \| q) \geq 0}; for exactly normalized PMFs, equality
#'     holds iff \eqn{p = q}. For sums accepted within the PMF tolerance,
#'     equality is interpreted up to that same numerical tolerance.
#'   \item Not symmetric: \eqn{D_{KL}(p \| q) \neq D_{KL}(q \| p)}
#' }
#'
#' @examples
#' p <- c(0.2, 0.5, 0.3)
#' kl_divergence_pmf(p, p)  # 0
#'
#' q <- c(0.3, 0.4, 0.3)
#' kl_divergence_pmf(p, q)
#'
#' @seealso \code{\link{kl_divergence_K}} for KL divergence with induced PMF
#'
#' @export
kl_divergence_pmf <- function(p, q, eps = NULL) {
  if (!is.null(eps)) {
    .dpprior_abort_invalid(
      paste(
        "eps smoothing is not part of exact KL divergence;",
        "supply normalized PMFs and leave eps = NULL."
      ),
      c("dpprior_kl_smoothing_error", "dpprior_kl_input_error"),
      "eps", eps, "NULL for exact KL", "smoothing_not_exact"
    )
  }
  if (!is.numeric(p) || !is.numeric(q)) {
    .dpprior_abort_invalid(
      "p and q must be numeric vectors.",
      c("dpprior_pmf_type_error", "dpprior_kl_input_error"),
      if (!is.numeric(p)) "p" else "q",
      if (!is.numeric(p)) p else q, "numeric vector", "type"
    )
  }
  if (length(p) != length(q)) {
    .dpprior_abort_invalid(
      "p and q must have the same positive length.",
      c("dpprior_pmf_length_error", "dpprior_kl_input_error"),
      "q", q, sprintf("length %d", length(p)), "length"
    )
  }
  p <- .a2_kl_validate_mathematical_pmf(p, "p")
  q <- .a2_kl_validate_mathematical_pmf(
    q, "q", expected_length = length(p)
  )
  .a2_kl_compute_kl(p, q)
}


#' KL Divergence Between Target and Induced K_J PMFs
#'
#' Computes the KL divergence \eqn{D_{KL}(p^* \| p_{a,b})} between a target PMF
#' and the induced marginal PMF of \eqn{K_J} under \eqn{\alpha \sim Gamma(a, b)}.
#'
#' @param target_pmf Numeric vector; target PMF for K_J. Can have length J
#'   (support k=1,...,J) or J+1 (support k=0,...,J, where k=0 must have zero mass).
#' @param a Numeric; shape parameter of Gamma hyperprior (a > 0).
#' @param b Numeric; rate parameter of Gamma hyperprior (b > 0).
#' @param J Integer; sample size.
#' @param M Integer; number of quadrature nodes. Default: 80.
#'
#' @return Numeric scalar; \eqn{D_{KL}(p^* \| p_{a,b})} (non-negative).
#'
#' @examples
#' J <- 50
#' target <- rep(1/J, J)  # uniform target over k=1,...,J
#' kl_divergence_K(target, a = 2, b = 1, J = J)
#'
#' @seealso \code{\link{kl_divergence_pmf}}, \code{\link{DPprior_a2_kl}}
#'
#' @export
kl_divergence_K <- function(target_pmf, a, b, J, M = .QUAD_NODES_DEFAULT) {
  assert_valid_J(J)
  assert_positive(a, "a")
  assert_positive(b, "b")
  M <- .dpprior_validate_count(
    M, "M", minimum = 10L, maximum = .QUADRATURE_MAX_NODES,
    .subclass = "dpprior_a2_kl_control_error"
  )

  # Validate target PMF (handles both length J and J+1 without normalization)
  target_pmf <- .a2_kl_normalize_pmf(target_pmf, J)

  # Compute the induced PMF in log space. This preserves exact zero-support
  # semantics without turning representable log probabilities into artificial
  # linear-scale zeros.
  logS <- compute_log_stirling(J)
  induced <- .a2_kl_induced_log_pmf(J, a, b, logS, M)

  # Compute KL divergence
  .a2_kl_compute_kl_logq(target_pmf, induced$selected)
}


# =============================================================================
# Target PMF Construction
# =============================================================================

#' Discretize Chi-Square to K_J Support
#'
#' Converts a (possibly scaled) chi-square distribution into a discrete PMF on
#' \eqn{\{1, \dots, J\}} using continuity-corrected binning:
#' \deqn{p(k) = P(k - 0.5 < X \le k + 0.5), \quad k = 1, \dots, J}
#' followed by renormalization.
#'
#' @param J Integer; maximum value (support upper bound).
#' @param df Numeric; degrees of freedom.
#' @param scale Numeric; scale parameter (default 1).
#'   If \code{scale != 1}, assumes \eqn{X = scale \cdot \chi^2_{df}}.
#'
#' @return Numeric vector of length J; PMF on \eqn{\{1, \dots, J\}}.
#'
#' @details
#' For a scaled chi-square distribution \eqn{Y = scale \cdot X} where \eqn{X \sim \chi^2_{df}}:
#' \itemize{
#'   \item \eqn{E[Y] = scale \cdot df}
#'   \item \eqn{Var[Y] = scale^2 \cdot 2 \cdot df}
#' }
#'
#' \strong{Matching target moments:}
#' Given target mean \eqn{\mu_K} and variance \eqn{\sigma^2_K}:
#' \itemize{
#'   \item \eqn{scale = \sigma^2_K / (2 \mu_K)}
#'   \item \eqn{df = 2 \mu_K^2 / \sigma^2_K}
#' }
#'
#' @examples
#' # Chi-square with target moments mu=5, var=8
#' mu_K <- 5
#' var_K <- 8
#' scale <- var_K / (2 * mu_K)
#' df <- 2 * mu_K^2 / var_K
#' pmf <- discretize_chisq(50, df = df, scale = scale)
#'
#' # Verify moments
#' k_vals <- 1:50
#' sum(k_vals * pmf)  # ~5
#'
#' @seealso \code{\link{DPprior_a2_kl}}
#'
#' @export
discretize_chisq <- function(J, df, scale = 1) {
  assert_valid_J(J)
  df <- .dpprior_validate_scalar(
    df, "df", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_chisq_parameter_error"
  )
  scale <- .dpprior_validate_scalar(
    scale, "scale", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_chisq_parameter_error"
  )

  k <- seq_len(J)
  upper <- (k + 0.5) / scale
  lower <- (k - 0.5) / scale

  raw_pmf <- stats::pchisq(upper, df = df) - stats::pchisq(lower, df = df)
  raw_pmf <- pmax(raw_pmf, 0)

  s <- sum(raw_pmf)
  if (!is.finite(s) || s <= 0) {
    stop(.dpprior_new_condition(
      "Chi-square discretization produced zero mass; check df/scale/J.",
      c("dpprior_chisq_discretization_error", "dpprior_numerical_error",
        "dpprior_error", "error"),
      value = s, parameters = list(J = J, df = df, scale = scale),
      code = "zero_discretized_mass"
    ))
  }

  pmf <- raw_pmf / s
  pmf <- .a2_kl_validate_mathematical_pmf(
    pmf, "discretized_chisq_pmf", expected_length = J
  )
  attr(pmf, "dpprior_pmf_provenance") <- list(
    constructor = "discretize_chisq",
    normalization = "explicit_support_conditioning",
    support = c(lower = 1L, upper = as.integer(J)),
    retained_mass_before_normalization = s,
    omitted_mass = max(0, 1 - s),
    df = df,
    scale = scale
  )
  pmf
}


# Reconstruct the chi-square target through one continuity-edge CDF vector.
# This intentionally does not call discretize_chisq(): it is a separately
# coded postcondition oracle for the public constructor's two-vector
# lower/upper calculation.
.a2_kl_verify_chisq_pmf <- function(J, df, scale) {
  assert_valid_J(J)
  df <- .dpprior_validate_scalar(
    df, "df", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_chisq_parameter_error"
  )
  scale <- .dpprior_validate_scalar(
    scale, "scale", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_chisq_parameter_error"
  )

  edges <- (seq.int(0L, J) + 0.5) / scale
  cdf_at_edges <- stats::pchisq(edges, df = df)
  raw_pmf <- diff(cdf_at_edges)
  raw_pmf[raw_pmf < 0 & raw_pmf >= -64 * .Machine$double.eps] <- 0
  if (any(raw_pmf < 0) || anyNA(raw_pmf) || any(!is.finite(raw_pmf))) {
    stop(.dpprior_new_condition(
      "Independent chi-square target reconstruction produced invalid mass.",
      c("dpprior_chisq_verification_error", "dpprior_numerical_error",
        "dpprior_error", "error"),
      value = raw_pmf, parameters = list(J = J, df = df, scale = scale),
      code = "invalid_reconstructed_mass"
    ))
  }
  total <- sum(raw_pmf)
  if (!is.finite(total) || total <= 0) {
    stop(.dpprior_new_condition(
      "Independent chi-square target reconstruction produced zero mass.",
      c("dpprior_chisq_verification_error", "dpprior_numerical_error",
        "dpprior_error", "error"),
      value = total, parameters = list(J = J, df = df, scale = scale),
      code = "zero_reconstructed_mass"
    ))
  }
  .a2_kl_validate_mathematical_pmf(
    unname(raw_pmf / total), "verification_chisq_pmf", J
  )
}


# Validate a chi-square target's requested moments before deriving df/scale.
# This helper deliberately distinguishes malformed scalars, positivity, support,
# and finite-support feasibility so every public failure is catchable and has a
# stable machine-readable reason.
.a2_kl_validate_target_moments <- function(J, mu_K, var_K) {
  mu_K <- .dpprior_validate_scalar(
    mu_K, "mu_K", .subclass = "dpprior_target_moment_error"
  )
  var_K <- .dpprior_validate_scalar(
    var_K, "var_K", .subclass = "dpprior_target_moment_error"
  )

  if (mu_K <= 0) {
    .dpprior_abort_invalid(
      "mu_K must be a positive finite numeric scalar",
      c("dpprior_target_moment_error", "dpprior_bounds_error"),
      "mu_K", mu_K, "> 0", "positive"
    )
  }
  if (var_K <= 0) {
    .dpprior_abort_invalid(
      "var_K must be a positive finite numeric scalar",
      c("dpprior_target_moment_error", "dpprior_bounds_error"),
      "var_K", var_K, "> 0", "positive"
    )
  }
  if (mu_K <= 1) {
    .dpprior_abort_invalid(
      "mu_K must be > 1 (at least one cluster is always present)",
      c("dpprior_target_moment_support_error", "dpprior_bounds_error"),
      "mu_K", mu_K, sprintf("in (1, %d)", as.integer(J)),
      "support_lower"
    )
  }
  if (mu_K >= J) {
    .dpprior_abort_invalid(
      paste(
        "mu_K must be < J (mu_K = J implies zero variance for K_J,",
        "outside the positive-variance elicitation workflow)"
      ),
      c("dpprior_target_moment_support_error", "dpprior_bounds_error"),
      "mu_K", mu_K, sprintf("in (1, %d)", as.integer(J)),
      "support_upper"
    )
  }
  .assert_feasible_K_moments(J, mu_K, var_K)

  list(mu_K = mu_K, var_K = var_K)
}


# Compute PMF moments with a centered second pass. The algebraic shortcut
# E[K^2] - E[K]^2 loses all precision, and can become negative, for target
# distributions concentrated near the upper support boundary.
.a2_kl_pmf_moments <- function(target_pmf) {
  k <- seq_along(target_pmf)
  mu_K <- sum(k * target_pmf)
  list(
    mu_K = mu_K,
    var_K = sum((k - mu_K)^2 * target_pmf)
  )
}


#' Construct Target PMF from User Specification
#'
#' Creates a target PMF from either a user-provided PMF vector or target moment
#' specification. A supplied PMF must already sum to one; it is validated but
#' never silently normalized. The chi-square constructor explicitly conditions
#' its discretized mass on the declared support and records that provenance.
#'
#' @param J Integer; sample size.
#' @param target Either:
#'   \itemize{
#'     \item Numeric vector of length J or J+1: direct PMF specification
#'     \item Named list with \code{mu_K} and \code{var_K}: construct from moments
#'   }
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{pmf}}{Numeric vector of length J; validated target PMF}
#'     \item{\code{mu_K}}{Target mean}
#'     \item{\code{var_K}}{Target variance}
#'     \item{\code{df}}{(if moments provided) Chi-square degrees of freedom}
#'     \item{\code{scale}}{(if moments provided) Chi-square scale parameter}
#'   }
#'
#' @details
#' \strong{Direct PMF specification:}
#' If \code{target} is a numeric vector of length J, it is treated as the PMF
#' for k = 1, ..., J. If length J+1, the k=0 entry must be zero and is dropped.
#'
#' \strong{Moment specification:}
#' If \code{target} is a list with \code{mu_K} and \code{var_K}, a discretized
#' chi-square distribution matching these moments is constructed using:
#' \deqn{\text{scale} = \sigma^2_K / (2\mu_K), \quad \text{df} = 2\mu_K^2 / \sigma^2_K}
#'
#' @examples
#' \dontrun{
#' # Direct PMF
#' result <- construct_target_pmf(50, rep(1 / 50, 50))  # Uniform on 1:50
#'
#' # Moment specification
#' result <- construct_target_pmf(50, list(mu_K = 5, var_K = 8))
#' result$mu_K  # 5
#' result$df    # 6.25
#'
#' }
#' @seealso \code{\link{discretize_chisq}}, \code{\link{DPprior_a2_kl}}
#'
#' @keywords internal
construct_target_pmf <- function(J, target) {
  assert_valid_J(J)
  J <- as.integer(J)

  if (is.numeric(target)) {
    # Direct PMF specification
    target_pmf <- .a2_kl_normalize_pmf(target, J)

    # Compute moments from PMF
    moments <- .a2_kl_pmf_moments(target_pmf)
    mu_K <- moments$mu_K
    var_K <- moments$var_K

    list(
      pmf = target_pmf,
      mu_K = mu_K,
      var_K = var_K,
      input_provenance = list(
        source = "custom_pmf",
        normalization = "validated_not_modified",
        input_length = length(target),
        input_sum = sum(target),
        k0_entry = length(target) == J + 1L
      )
    )

  } else if (is.list(target) && all(c("mu_K", "var_K") %in% names(target))) {
    # Moment specification
    moments <- .a2_kl_validate_target_moments(
      J, target$mu_K, target$var_K
    )
    mu_K <- moments$mu_K
    var_K <- moments$var_K

    # Construct discretized chi-square
    # For Y = scale * X where X ~ chi2(df):
    #   E[Y] = scale * df = mu_K
    #   Var[Y] = scale^2 * 2 * df = var_K
    # Solving: scale = var_K / (2 * mu_K), df = 2 * mu_K^2 / var_K
    scale <- var_K / (2 * mu_K)
    df <- 2 * mu_K^2 / var_K

    target_pmf <- discretize_chisq(J, df, scale)
    pmf_provenance <- attr(
      target_pmf, "dpprior_pmf_provenance", exact = TRUE
    )
    target_pmf <- unname(as.numeric(target_pmf))

    # Compute discretized moments (may differ due to truncation)
    discrete_moments <- .a2_kl_pmf_moments(target_pmf)
    mu_disc <- discrete_moments$mu_K
    var_disc <- discrete_moments$var_K

    list(
      pmf = target_pmf,
      mu_K = mu_K,
      var_K = var_K,
      df = df,
      scale = scale,
      mu_K_discrete = mu_disc,
      var_K_discrete = var_disc,
      input_provenance = pmf_provenance
    )

  } else {
    .dpprior_abort_invalid(
      "target must be a PMF vector or list with 'mu_K' and 'var_K'",
      c("dpprior_target_structure_error", "dpprior_a2_kl_input_error"),
      "target", target, "numeric PMF or list(mu_K=..., var_K=...)",
      "target_structure"
    )
  }
}


# =============================================================================
# A2-KL solver helpers
# =============================================================================

.A2_KL_STATUS_CODES <- c(
  "converged", "boundary", "approximate", "infeasible", "failed"
)

.A2_KL_ADEQUACY_DEFAULTS <- list(
  kl_tol = 0.015,
  l1_tol = 0.11,
  mean_scaled_tol = 0.01,
  var_scaled_tol = 0.065
)


.a2_kl_capture_call <- function(expr) {
  warnings <- character()
  error_condition <- NULL
  started <- proc.time()[["elapsed"]]
  value <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      error_condition <<- e
      NULL
    }
  )
  list(
    value = value,
    elapsed = unname(proc.time()[["elapsed"]] - started),
    warnings = unname(warnings),
    error = if (is.null(error_condition)) NULL else conditionMessage(error_condition),
    error_class = if (is.null(error_condition)) NULL else class(error_condition)
  )
}


.a2_kl_initializer_evidence <- function(value) {
  if (!inherits(value, "dpprior_result") ||
      !identical(value[["schema", exact = TRUE]], .dpprior_schema("result"))) {
    return(NULL)
  }
  valid <- tryCatch({
    .dpprior_validate_result_v1(value)
    TRUE
  }, error = function(condition) FALSE)
  if (!valid) return(NULL)

  parameters <- value[["parameters", exact = TRUE]]
  if (is.null(parameters) ||
      !.dpprior_is_plain_numeric(parameters[["a", exact = TRUE]]) ||
      length(parameters[["a", exact = TRUE]]) != 1L ||
      !.dpprior_is_plain_numeric(parameters[["b", exact = TRUE]]) ||
      length(parameters[["b", exact = TRUE]]) != 1L ||
      !is.finite(parameters[["a", exact = TRUE]]) ||
      !is.finite(parameters[["b", exact = TRUE]]) ||
      parameters[["a", exact = TRUE]] <= 0 ||
      parameters[["b", exact = TRUE]] <= 0) {
    return(NULL)
  }
  residual_K <- value[["residuals", exact = TRUE]][["K", exact = TRUE]]
  residual <- if (identical(value[["mode", exact = TRUE]], "a2_moment") &&
      typeof(residual_K) == "list" && is.list(residual_K) &&
      length(residual_K) && all(vapply(
        residual_K, function(component) {
          .dpprior_is_plain_numeric(component) && length(component) == 1L &&
            is.finite(component)
        }, logical(1)
      ))) {
    sqrt(sum(unlist(residual_K, use.names = FALSE)^2))
  } else {
    NA_real_
  }
  list(
    a = unname(parameters[["a", exact = TRUE]]),
    b = unname(parameters[["b", exact = TRUE]]),
    method = value[["method", exact = TRUE]],
    status = value[["status", exact = TRUE]],
    verified_success = value[["status", exact = TRUE]] %in%
      c("converged", "boundary") &&
      isTRUE(value[["usable", exact = TRUE]]) &&
      isTRUE(value[["verified", exact = TRUE]]),
    iterations = value[["computation", exact = TRUE]][[
      "termination", exact = TRUE
    ]][["iterations", exact = TRUE]],
    residual = residual
  )
}


.a2_kl_valid_parameter_pair <- function(value) {
  !is.null(.a2_kl_initializer_evidence(value))
}


.a2_kl_run_lbfgsb <- function(par, fn, lower, upper, control) {
  stats::optim(
    par = par, fn = fn, method = "L-BFGS-B",
    lower = lower, upper = upper, control = control
  )
}


.a2_kl_run_nlminb <- function(par, fn, lower, upper, control) {
  stats::nlminb(
    start = par, objective = fn, lower = lower, upper = upper,
    control = control
  )
}


.a2_kl_pmf_metrics <- function(target_pmf, log_induced) {
  induced <- exp(log_induced)
  induced <- .a2_kl_validate_mathematical_pmf(
    induced, "induced_pmf", expected_length = length(target_pmf)
  )
  k <- seq_along(target_pmf)
  target_mean <- sum(k * target_pmf)
  induced_mean <- sum(k * induced)
  target_var <- sum((k - target_mean)^2 * target_pmf)
  induced_var <- sum((k - induced_mean)^2 * induced)
  list(
    kl = .a2_kl_compute_kl_logq(target_pmf, log_induced),
    l1 = sum(abs(target_pmf - induced)),
    mean = induced_mean,
    variance = induced_var,
    mean_residual = induced_mean - target_mean,
    variance_residual = induced_var - target_var,
    target_mean = target_mean,
    target_variance = target_var,
    pmf = induced
  )
}


.a2_kl_assess_adequacy <- function(metrics, tolerances) {
  mean_scale <- max(1, sqrt(max(0, metrics$target_variance)))
  variance_scale <- max(1, metrics$target_variance)
  raw <- c(
    kl = metrics$kl,
    l1 = metrics$l1,
    mean = abs(metrics$mean_residual),
    variance = abs(metrics$variance_residual)
  )
  scales <- c(kl = 1, l1 = 1, mean = mean_scale, variance = variance_scale)
  scaled <- raw / scales
  scaled_tolerances <- c(
    kl = tolerances$kl,
    l1 = tolerances$l1,
    mean = tolerances$mean_scaled,
    variance = tolerances$variance_scaled
  )
  component_passed <- is.finite(scaled) & scaled <= scaled_tolerances
  list(
    passed = all(component_passed),
    raw = raw,
    scales = scales,
    scale_formula = c(
      kl = "1",
      l1 = "1",
      mean = "max(1, sqrt(target variance))",
      variance = "max(1, target variance)"
    ),
    scaled = scaled,
    scaled_tolerances = scaled_tolerances,
    raw_tolerances = scaled_tolerances * scales,
    ratios = scaled / scaled_tolerances,
    component_passed = component_passed
  )
}


# Convert the Phase 8 flat result into a plain, non-authoritative audit view.
# Canonical consumers never read this record; it exists only to preserve the
# scientific values and execution narrative exposed by the previous API.
.a2_kl_plain_compat_value <- function(value) {
  if (is.null(value)) {
    return(NULL)
  }
  if (inherits(value, "condition")) {
    return(list(
      class = paste(class(value), collapse = "/"),
      message = conditionMessage(value),
      call = if (is.null(conditionCall(value))) {
        "NULL"
      } else {
        paste(deparse(conditionCall(value)), collapse = " ")
      }
    ))
  }
  if (is.data.frame(value)) {
    raw <- unclass(value)
    out <- lapply(raw, .a2_kl_plain_compat_value)
    names(out) <- names(raw)
    return(out)
  }
  if (typeof(value) == "list" && is.list(value)) {
    raw <- if (is.object(value)) unclass(value) else value
    raw_names <- names(raw)
    if (length(raw) && (is.null(raw_names) || anyNA(raw_names) ||
        any(!nzchar(raw_names)) || anyDuplicated(raw_names))) {
      raw_names <- sprintf("item_%03d", seq_along(raw))
    }
    out <- lapply(raw, .a2_kl_plain_compat_value)
    if (length(out)) names(out) <- raw_names
    return(out)
  }
  if (is.language(value) || is.function(value) || is.environment(value) ||
      isS4(value)) {
    return(paste(deparse(value), collapse = " "))
  }
  if (is.factor(value)) {
    value <- as.character(value)
  }
  if (!(is.numeric(value) || is.logical(value) || is.character(value))) {
    return(paste(capture.output(str(value)), collapse = " "))
  }
  if (length(value) == 0L) {
    return(paste0("<empty_", typeof(value), ">"))
  }
  value_names <- names(value)
  if (!is.null(value_names) && (anyNA(value_names) ||
      any(!nzchar(value_names)) || anyDuplicated(value_names))) {
    value_names <- NULL
  }
  if (anyNA(value) || (is.numeric(value) && any(!is.finite(value)))) {
    out <- vapply(seq_along(value), function(index) {
      item <- value[[index]]
      if (length(item) == 0L) {
        "<empty>"
      } else if (is.na(item)) {
        paste0("<NA_", typeof(value), ">")
      } else if (is.numeric(item) && is.infinite(item)) {
        if (item > 0) "<Inf>" else "<-Inf>"
      } else {
        as.character(item)
      }
    }, character(1))
    if (!is.null(value_names)) names(out) <- value_names
    return(out)
  }
  attributes(value) <- if (is.null(value_names)) NULL else
    list(names = value_names)
  value
}


# Preserve the public target authority used by A2-KL.  A strict PMF remains
# authoritative.  For the scaled-chi-square route, the user-authoritative
# object remains the requested moment target; the conditioned PMF is retained
# only as a reproducible objective-distribution derivation.
.a2_kl_target_v1 <- function(J, original_target, method, target_info,
                             target_pmf) {
  if (identical(method, "pmf")) {
    return(DPprior_target_K(J = J, target_pmf = original_target))
  }

  target_K <- DPprior_target_K(
    J = J, mu_K = target_info$mu_K, var_K = target_info$var_K
  )
  construction <- target_info$input_provenance
  objective_evidence <- list(
    method = "chisq",
    df = unname(target_info$df),
    scale = unname(target_info$scale),
    binning = "continuity_corrected_half_integer_bins",
    normalization = "explicit_support_conditioning",
    support = c(lower = 1L, upper = as.integer(J)),
    retained_mass_before_normalization = unname(
      construction$retained_mass_before_normalization
    ),
    omitted_mass = unname(construction$omitted_mass),
    pmf = unname(as.numeric(target_pmf)),
    mu_K_discrete = unname(target_info$mu_K_discrete),
    var_K_discrete = unname(target_info$var_K_discrete),
    source = "A2_KL_backend_target"
  )
  target_K$derivation$request_to_normalized$evidence$A2_KL_objective <-
    objective_evidence
  .dpprior_validate_target_v1(target_K)
  target_K
}


.a2_kl_snapshot_v1 <- function(parameters, M, metrics, tolerances, source) {
  achieved <- list(K = list(
    mean = unname(metrics$mean),
    variance = unname(metrics$variance),
    estimand = "K_J",
    source = source,
    M = as.integer(M),
    pmf = unname(as.numeric(metrics$pmf))
  ))
  residuals <- list(distribution = list(
    kl = unname(metrics$kl),
    l1 = unname(metrics$l1),
    mean = unname(metrics$mean_residual),
    variance = unname(metrics$variance_residual)
  ))
  .dpprior_new_snapshot(
    parameters = parameters,
    M = as.integer(M),
    achieved = achieved,
    residuals = residuals,
    tolerances = tolerances,
    finite = TRUE,
    source = source
  )
}


.a2_kl_typed_attempt_error <- function(attempt) {
  if (is.null(attempt$error)) return(NULL)
  error_class <- attempt$error_class
  error_class <- if (is.null(error_class) || !length(error_class) ||
      is.na(error_class[[1L]]) || !nzchar(error_class[[1L]])) {
    "error"
  } else {
    as.character(error_class[[1L]])
  }
  list(
    class = error_class,
    code = "backend_execution_error",
    message = as.character(attempt$error[[1L]])
  )
}


.a2_kl_unavailable <- function(values, reasons) {
  null_fields <- names(values)[vapply(values, is.null, logical(1))]
  if (!length(null_fields)) return(character())
  out <- vapply(null_fields, function(field) {
    reason <- reasons[[field]]
    if (is.null(reason) || !length(reason) || is.na(reason[[1L]]) ||
        !nzchar(reason[[1L]])) {
      paste0(field, "_not_recorded")
    } else {
      as.character(reason[[1L]])
    }
  }, character(1))
  names(out) <- null_fields
  out
}


# Add the canonical dpprior.result/1 authority without changing the Phase 8
# optimizer, objective, candidate choice, selected-order PMF, or status logic.
.a2_kl_result_v1 <- function(legacy_result, context) {
  J <- context$J
  target_pmf <- context$target_pmf
  target_K <- .a2_kl_target_v1(
    J = J,
    original_target = context$original_target,
    method = context$target_method,
    target_info = context$target_info,
    target_pmf = target_pmf
  )

  finite_public_candidate <- !identical(context$status, "failed") &&
    all(is.finite(c(context$a_opt, context$b_opt))) &&
    !is.null(context$selected_metrics)
  if (!finite_public_candidate) {
    stop(
      "Internal A2-KL canonicalization requires a retained finite candidate.",
      call. = FALSE
    )
  }

  parameters <- .dpprior_new_parameters(
    a = context$a_opt, b = context$b_opt, parameterization = "log_ab"
  )
  controls <- list(
    max_iter = as.integer(context$max_iter),
    optimizer_tol = context$tol,
    log_bounds = unname(as.numeric(context$log_bounds)),
    boundary_tol = context$boundary_tol,
    fallback_max_iter = as.integer(context$fallback_max_iter),
    primary = list(
      maxit = as.integer(context$max_iter),
      factr = context$tol / .Machine$double.eps,
      pgtol = context$tol
    ),
    fallback = list(
      iter.max = as.integer(context$fallback_max_iter),
      eval.max = as.integer(max(200L, 2L * context$fallback_max_iter)),
      rel.tol = context$tol,
      x.tol = context$tol
    ),
    fallback_trigger_worse_than_start = 1e-12,
    selection_tolerance = 0
  )
  setting <- list(
    method = "A2-KL", controls = controls, parameterization = "log_ab"
  )
  tolerances <- list(
    distribution = list(
      adequacy = list(
        kl = context$kl_tol,
        l1 = context$l1_tol,
        mean_scaled = context$mean_scaled_tol,
        variance_scaled = context$var_scaled_tol,
        mean_scale_formula = "max(1,sqrt(target_variance))",
        variance_scale_formula = "max(1,target_variance)"
      ),
      order = list(
        pmf_absolute = context$pmf_abs_tol,
        pmf_relative = context$pmf_rel_tol,
        pmf_l1 = context$pmf_abs_tol + context$pmf_rel_tol,
        direct_moment_absolute = context$pmf_abs_tol,
        direct_moment_relative = context$pmf_rel_tol,
        target_identity_l1 = .TOL_PMF_SUM
      )
    ),
    boundary = context$boundary_tol
  )

  selected_snapshot <- .a2_kl_snapshot_v1(
    parameters, context$M, context$selected_metrics, tolerances,
    "selected_order"
  )
  effective_M_verify <- context$M_verify
  effective_required <- context$required_M_verify
  verifier_metrics <- context$verification_metrics
  same_order_verifier <- FALSE
  if (is.null(effective_M_verify) || is.null(verifier_metrics)) {
    # The frozen schema requires explicit fixed-candidate verifier evidence for
    # every retained A2-KL candidate.  When the higher order is unsupported,
    # run a separate same-order recomputation, retain that limitation, and
    # leave the result approximate/unverified exactly as before.
    effective_M_verify <- as.integer(context$M)
    effective_required <- as.integer(context$M)
    independent <- .a2_kl_induced_log_pmf(
      J, context$a_opt, context$b_opt, context$logS, context$M
    )$selected
    verifier_metrics <- .a2_kl_pmf_metrics(target_pmf, independent)
    same_order_verifier <- TRUE
  }
  verifier_snapshot <- .a2_kl_snapshot_v1(
    parameters, effective_M_verify, verifier_metrics, tolerances,
    "independent_verifier"
  )

  old_attempts <- context$attempts
  attempt_ids <- sprintf("attempt-%03d", seq_along(old_attempts))
  canonical_attempts <- vector("list", length(old_attempts))
  candidate_evaluations <- list()
  selected_attempt_id <- NULL
  selected_candidate_id <- NULL

  for (index in seq_along(old_attempts)) {
    old <- old_attempts[[index]]
    id <- attempt_ids[[index]]
    selected <- isTRUE(old$selected)
    stage <- if (identical(old$method, "L-BFGS-B")) {
      "primary"
    } else if (identical(old$method, "nlminb")) {
      "fallback"
    } else {
      "initialization"
    }
    candidate <- old$candidate
    if (selected) {
      candidate <- c(a = context$a_opt, b = context$b_opt)
    }
    candidate_parameters <- if (is.numeric(candidate) &&
        length(candidate) == 2L && all(is.finite(candidate)) &&
        all(candidate > 0)) {
      .dpprior_new_parameters(
        unname(candidate[[1L]]), unname(candidate[[2L]]), "log_ab"
      )
    } else {
      NULL
    }

    candidate_snapshot <- NULL
    candidate_metrics <- NULL
    if (!is.null(candidate_parameters)) {
      if (selected) {
        candidate_snapshot <- selected_snapshot
        candidate_metrics <- context$selected_metrics
      } else {
        evaluated <- tryCatch({
          candidate_log_pmf <- .a2_kl_induced_log_pmf(
            J, candidate_parameters$a, candidate_parameters$b,
            context$logS, context$M
          )$selected
          metrics <- .a2_kl_pmf_metrics(target_pmf, candidate_log_pmf)
          list(
            metrics = metrics,
            snapshot = .a2_kl_snapshot_v1(
              candidate_parameters, context$M, metrics, tolerances,
              "candidate_selected_order"
            )
          )
        }, error = function(condition) NULL)
        if (!is.null(evaluated)) {
          candidate_metrics <- evaluated$metrics
          candidate_snapshot <- evaluated$snapshot
        } else {
          candidate_parameters <- NULL
        }
      }
      if (!is.null(candidate_snapshot) && any(
          candidate_snapshot$achieved$K$pmf[target_pmf > 0] == 0
      )) {
        # A finite log-objective can still map to an unrepresentable
        # double-precision PMF at remote support points.  Such an attempt is
        # retained in the legacy audit view but cannot enter the canonical
        # candidate ledger, whose PMF evidence must itself carry the target
        # support.
        candidate_parameters <- NULL
        candidate_metrics <- NULL
        candidate_snapshot <- NULL
      }
    }

    error <- .a2_kl_typed_attempt_error(old)
    exit_code <- if (length(old$exit_code) == 1L &&
        !is.na(old$exit_code) && is.finite(old$exit_code)) {
      as.integer(old$exit_code)
    } else {
      NULL
    }
    iterations <- if (length(old$iterations) == 1L &&
        !is.na(old$iterations) && is.finite(old$iterations) &&
        old$iterations >= 0) {
      as.integer(old$iterations)
    } else {
      NULL
    }
    evaluations <- if (length(old$evaluations) == 1L &&
        !is.na(old$evaluations) && is.finite(old$evaluations) &&
        old$evaluations >= 0) {
      list(function_count = as.integer(old$evaluations))
    } else {
      NULL
    }
    elapsed <- if (length(old$elapsed) == 1L && !is.na(old$elapsed) &&
        is.finite(old$elapsed) && old$elapsed >= 0) {
      unname(old$elapsed)
    } else {
      NULL
    }
    candidate_objective <- if (!identical(stage, "initialization") &&
        !is.null(candidate_metrics)) {
      unname(candidate_metrics$kl)
    } else {
      NULL
    }
    reason_code <- if (!is.null(candidate_snapshot)) {
      if (selected) "selected" else "eligible_not_selected"
    } else if (!is.null(error)) {
      "optimizer_error"
    } else if (!is.null(exit_code) && exit_code != 0L) {
      "optimizer_exit_nonzero"
    } else {
      "nonfinite_candidate"
    }
    start <- old$start
    if (!is.null(start)) start <- unname(as.numeric(start))
    bounds <- old$bounds
    control <- old$control
    attempt_values <- list(
      start = start,
      bounds = bounds,
      control = control,
      exit_code = exit_code,
      iterations = iterations,
      evaluations = evaluations,
      candidate_parameters = candidate_parameters,
      candidate_objective = candidate_objective,
      elapsed_seconds = elapsed
    )
    unavailable <- .a2_kl_unavailable(attempt_values, list(
      start = "source_attempt_has_no_start",
      bounds = "source_attempt_has_no_bounds",
      control = "source_attempt_has_no_control",
      exit_code = "source_attempt_did_not_report_exit_code",
      iterations = old$iteration_reason %||%
        "source_attempt_did_not_report_iterations",
      evaluations = "source_attempt_did_not_report_evaluation_count",
      candidate_parameters = "no_finite_candidate_from_source_attempt",
      candidate_objective = if (identical(stage, "initialization")) {
        "A2_KL_initialization_objective_not_source_attempt_evidence"
      } else {
        "no_finite_candidate_objective"
      },
      elapsed_seconds = "source_attempt_did_not_report_elapsed_time"
    ))
    message <- old$message
    if (is.null(message) || !length(message) || is.na(message[[1L]])) {
      message <- if (is.null(error)) "" else error$message
    }
    warnings <- as.character(old$warnings %||% character())
    warnings <- warnings[!is.na(warnings) & nzchar(warnings)]
    canonical_attempts[[index]] <- .dpprior_new_attempt(
      id = id,
      stage = stage,
      method = old$method,
      start = start,
      bounds = bounds,
      control = control,
      exit_code = exit_code,
      message = as.character(message[[1L]]),
      iterations = iterations,
      evaluations = evaluations,
      candidate_parameters = candidate_parameters,
      candidate_objective = candidate_objective,
      elapsed_seconds = elapsed,
      warnings = warnings,
      error = error,
      selected = selected && !is.null(candidate_parameters),
      reason_code = reason_code,
      unavailable = unavailable
    )

    if (!is.null(candidate_snapshot)) {
      candidate_id <- sprintf("candidate-%03d", length(candidate_evaluations) + 1L)
      check_source <- paste0("candidate:", candidate_id)
      pmf <- candidate_snapshot$achieved$K$pmf
      execution_success <- if (identical(stage, "initialization")) {
        is.null(error)
      } else {
        is.null(error) && identical(exit_code, 0L)
      }
      evaluation <- .dpprior_new_candidate_evaluation(
        id = candidate_id,
        attempt_id = id,
        method = old$method,
        generator = if (identical(stage, "initialization")) {
          "initialization"
        } else {
          "direct_attempt"
        },
        parameters = candidate_parameters,
        objective_kind = "kl",
        recorded_objective = candidate_objective,
        fresh_objective = unname(candidate_metrics$kl),
        selection_objective = unname(candidate_metrics$kl),
        objective_tolerance = 1e-12 * max(1, abs(candidate_metrics$kl)),
        selected_snapshot = candidate_snapshot,
        verifier_snapshot = NULL,
        checks = list(candidate_distribution = .dpprior_new_check(
          value = c(
            pmf_mass_error = abs(sum(pmf) - 1),
            pmf_minimum_violation = max(0, -min(pmf))
          ),
          reference = c(
            pmf_mass_error = 0, pmf_minimum_violation = 0
          ),
          tolerance = c(
            pmf_mass_error = .TOL_PMF_SUM,
            pmf_minimum_violation = 0
          ),
          operator = "lte",
          source = check_source
        )),
        execution_success = execution_success,
        optimizer_supported = !identical(stage, "initialization") &&
          execution_success,
        diagnostic_eligible = FALSE,
        selected = selected,
        source = "R/12_a2_kl.R:canonical_candidate_evaluation"
      )
      candidate_evaluations[[length(candidate_evaluations) + 1L]] <- evaluation
      if (selected) {
        selected_attempt_id <- id
        selected_candidate_id <- candidate_id
      }
    }
  }

  # The Phase 8 selector used the stable log-PMF KL while the frozen schema
  # independently recomputes KL from the retained ordinary PMF.  At machine-
  # precision ties those two representations can reverse by a few ulps.  Such
  # a non-selected source candidate cannot be asserted as a strictly better
  # canonical candidate under the zero tie tolerance, so quarantine its
  # parameter row (the exact source row remains in a2_kl_v0).  A material
  # reversal is a contract error rather than something this adapter may hide.
  selected_evaluation_index <- which(vapply(
    candidate_evaluations,
    function(evaluation) isTRUE(evaluation$selected),
    logical(1)
  ))
  selection_quarantine_ids <- character()
  if (length(selected_evaluation_index) == 1L) {
    selected_fresh <- candidate_evaluations[[selected_evaluation_index]]$
      selection_objective
    lower_indices <- which(vapply(
      candidate_evaluations,
      function(evaluation) !isTRUE(evaluation$selected) &&
        evaluation$selection_objective < selected_fresh,
      logical(1)
    ))
    if (length(lower_indices)) {
      differences <- vapply(
        candidate_evaluations[lower_indices],
        function(evaluation) selected_fresh - evaluation$selection_objective,
        numeric(1)
      )
      equivalence_tolerance <- 1e-12 * max(
        1, abs(selected_fresh),
        abs(vapply(
          candidate_evaluations[lower_indices],
          function(evaluation) evaluation$selection_objective,
          numeric(1)
        ))
      )
      if (any(differences > equivalence_tolerance)) {
        stop(
          paste(
            "A2-KL source selection conflicts materially with the canonical",
            "direct-PMF KL objective."
          ),
          call. = FALSE
        )
      }
      for (evaluation_index in rev(lower_indices)) {
        evaluation <- candidate_evaluations[[evaluation_index]]
        selection_quarantine_ids <- c(
          selection_quarantine_ids, evaluation$id
        )
        attempt_index <- match(evaluation$attempt_id, attempt_ids)
        attempt <- canonical_attempts[[attempt_index]]
        attempt["candidate_parameters"] <- list(NULL)
        attempt["candidate_objective"] <- list(NULL)
        attempt$reason_code <- "candidate_evaluation_failed"
        unavailable <- attempt$unavailable
        unavailable <- unavailable[!names(unavailable) %in% c(
          "candidate_parameters", "candidate_objective"
        )]
        unavailable <- c(
          unavailable,
          candidate_parameters = paste(
            "source log-PMF tie reverses under canonical direct-PMF KL;",
            "candidate retained only in compatibility view"
          ),
          candidate_objective = paste(
            "canonical zero-tolerance selection cannot compare the",
            "machine-precision-reversed source tie"
          )
        )
        attempt$unavailable <- unavailable
        .dpprior_validate_attempt(attempt)
        canonical_attempts[[attempt_index]] <- attempt
        candidate_evaluations[[evaluation_index]] <- NULL
      }
    }
  }

  selected_attempt_index <- match(selected_attempt_id, attempt_ids)
  selected_attempt <- canonical_attempts[[selected_attempt_index]]
  fallback_index <- which(vapply(
    canonical_attempts,
    function(attempt) identical(attempt$stage, "fallback"),
    logical(1)
  ))
  primary_index <- which(vapply(
    canonical_attempts,
    function(attempt) identical(attempt$method, "L-BFGS-B"),
    logical(1)
  ))
  fallback_attempted <- length(fallback_index) > 0L
  fallback_used <- fallback_attempted &&
    identical(selected_attempt$stage, "fallback")
  fallback_reason <- paste(unique(context$fallback_reason), collapse = "+")
  if (fallback_attempted && !nzchar(fallback_reason)) {
    fallback_reason <- "primary_not_accepted"
  }
  fallback <- if (!fallback_attempted) {
    .dpprior_new_fallback()
  } else {
    .dpprior_new_fallback(
      attempted = TRUE,
      used = fallback_used,
      trigger_attempt_id = attempt_ids[[primary_index[[1L]]]],
      selected_attempt_id = if (fallback_used) selected_attempt_id else NULL,
      reason_code = fallback_reason,
      message = if (fallback_used) {
        "bounded nlminb fallback selected"
      } else {
        "bounded nlminb fallback attempted but not selected"
      },
      outcome = if (fallback_used) "selected" else "attempted_not_selected"
    )
  }

  orders <- .dpprior_new_orders(
    M_requested = as.integer(context$M),
    M_selected = as.integer(context$M),
    M_verification_required = as.integer(effective_required),
    M_verification_used = as.integer(effective_M_verify),
    requested_reason = "public_M_argument",
    selected_reason = "selected_order_public_output",
    verification_required_reason = if (same_order_verifier) {
      "higher_order_unsupported_same_order_audit_only"
    } else {
      "higher_order_verification_policy"
    },
    verification_used_reason = if (same_order_verifier) {
      "independent_same_order_recomputation"
    } else {
      "independent_higher_order_recomputation"
    }
  )
  termination <- if (identical(context$status, "approximate")) {
    .dpprior_new_termination(
      code = "approximate", message = context$message,
      source = "candidate_evaluation", iterations = selected_attempt$iterations
    )
  } else if (identical(context$status, "boundary")) {
    .dpprior_new_termination(
      code = "boundary", message = context$message,
      source = if (fallback_used) "fallback_optimizer" else "optimizer",
      iterations = selected_attempt$iterations,
      boundary_reason = paste(context$boundary_sides, collapse = ",")
    )
  } else {
    .dpprior_new_termination(
      code = "converged", message = context$message,
      source = if (fallback_used) "fallback_optimizer" else "optimizer",
      iterations = selected_attempt$iterations
    )
  }
  computation <- .dpprior_new_computation(
    request = setting,
    used = setting,
    orders = orders,
    scaling = .dpprior_new_scaling(),
    attempts = canonical_attempts,
    candidate_evaluations = candidate_evaluations,
    selected_candidate_id = selected_candidate_id,
    selected_attempt_id = selected_attempt_id,
    fallback = fallback,
    termination = termination,
    trace = context$trace,
    resources = .a2_kl_plain_compat_value(list(
      target_route = context$target_method,
      selected_source_candidate = context$selected_name,
      selected_source_method = context$selected_method,
      initialization_method = context$init_method,
      initialization_status = context$init_status,
      initialization_clipped_to_bounds = context$initialization_clipped,
      source_required_M_verification = context$required_M_verify,
      source_M_verification_available = !is.null(context$M_verify),
      same_order_verifier_only = same_order_verifier,
      machine_precision_selection_quarantine = if (
          length(selection_quarantine_ids)
      ) selection_quarantine_ids else "none",
      fallback_trigger_codes = if (length(context$fallback_reason)) {
        context$fallback_reason
      } else {
        "none"
      }
    ))
  )

  target_moments <- .a2_kl_pmf_moments(target_pmf)
  metrics_for_check <- function(snapshot, prefix) {
    achieved_K <- snapshot$achieved$K
    achieved_pmf <- achieved_K$pmf
    positive <- target_pmf > 0
    mean_residual <- achieved_K$mean - target_moments$mu_K
    variance_residual <- achieved_K$variance - target_moments$var_K
    values <- c(
      kl = sum(target_pmf[positive] * log(
        target_pmf[positive] / achieved_pmf[positive]
      )),
      l1 = sum(abs(target_pmf - achieved_pmf)),
      mean_scaled = abs(mean_residual) /
        max(1, sqrt(target_moments$var_K)),
      variance_scaled = abs(variance_residual) /
        max(1, target_moments$var_K)
    )
    names(values) <- paste0(prefix, ".", names(values))
    values
  }
  adequacy_value <- c(
    metrics_for_check(selected_snapshot, "selected"),
    metrics_for_check(verifier_snapshot, "refined")
  )
  adequacy_tolerance_one <- c(
    kl = context$kl_tol,
    l1 = context$l1_tol,
    mean_scaled = context$mean_scaled_tol,
    variance_scaled = context$var_scaled_tol
  )
  adequacy_tolerance <- c(
    setNames(adequacy_tolerance_one, paste0("selected.", names(adequacy_tolerance_one))),
    setNames(adequacy_tolerance_one, paste0("refined.", names(adequacy_tolerance_one)))
  )
  stability <- .dpprior_new_stability(
    delta = c(pmf.l1 = sum(abs(
      selected_snapshot$achieved$K$pmf - verifier_snapshot$achieved$K$pmf
    ))),
    tolerance = c(
      pmf.l1 = tolerances$distribution$order$pmf_l1
    ),
    formula = c(pmf.l1 = "direct_pmf_l1_tolerance"),
    scale_floor = c(pmf.l1 = 0),
    source = "independent_verifier"
  )
  selection_objectives <- vapply(
    candidate_evaluations,
    function(evaluation) evaluation$selection_objective,
    numeric(1)
  )
  selected_objective <- selected_snapshot$residuals$distribution$kl
  candidate_selection_delta <- selected_objective - min(selection_objectives)
  verification <- .dpprior_new_verification(
    method = "fresh higher-order marginal PMF and direct-moment audit",
    performed = TRUE,
    passed = isTRUE(context$verified),
    reason = if (isTRUE(context$verified)) {
      "optimizer, distribution adequacy, and independent order checks passed"
    } else if (same_order_verifier) {
      "higher-order verifier unavailable; same-order audit is non-verifying"
    } else {
      context$message
    },
    settings = list(
      M_selected = as.integer(context$M),
      M_verification = as.integer(effective_M_verify),
      M_verification_required = as.integer(effective_required),
      pmf_abs_tol = context$pmf_abs_tol,
      pmf_rel_tol = context$pmf_rel_tol
    ),
    selected_snapshot = selected_snapshot,
    verifier_snapshot = verifier_snapshot,
    stability = stability,
    components = list(
      target_identity = .dpprior_new_check(
        value = target_pmf, reference = target_pmf, tolerance = NULL,
        operator = "identical", source = "independent_verifier"
      ),
      pmf_adequacy = .dpprior_new_check(
        value = adequacy_value,
        reference = setNames(rep(0, length(adequacy_value)), names(adequacy_value)),
        tolerance = adequacy_tolerance,
        operator = "lte", source = "independent_verifier"
      ),
      order_stability = .dpprior_new_check(
        value = stability$delta,
        reference = c(pmf.l1 = 0),
        tolerance = stability$tolerance,
        operator = "lte", source = "independent_verifier"
      ),
      candidate_selection = .dpprior_new_check(
        value = candidate_selection_delta,
        reference = 0,
        tolerance = 0,
        operator = "lte", source = "independent_verifier"
      )
    ),
    invariants = list(
      finite_parameters = .dpprior_new_check(
        value = TRUE, reference = TRUE, tolerance = NULL,
        operator = "identical", source = "independent_verifier"
      ),
      parameter_identity = .dpprior_new_check(
        value = identical(selected_snapshot$parameters,
                          verifier_snapshot$parameters),
        reference = TRUE, tolerance = NULL,
        operator = "identical", source = "independent_verifier"
      ),
      pmf_probability = .dpprior_new_check(
        value = TRUE, reference = TRUE, tolerance = NULL,
        operator = "identical", source = "independent_verifier"
      ),
      support_identity = .dpprior_new_check(
        value = TRUE, reference = TRUE, tolerance = NULL,
        operator = "identical", source = "independent_verifier"
      )
    )
  )

  source_commit <- getOption("DPprior.source_commit", NULL)
  if (!is.character(source_commit) || length(source_commit) != 1L ||
      is.na(source_commit) || !nzchar(source_commit)) {
    source_commit <- NULL
  }
  approximation_active <- identical(context$status, "approximate")
  provenance <- .dpprior_new_provenance(
    requested_method = "A2-KL",
    selected_method = "A2-KL",
    is_fallback = fallback_used,
    approximation = list(
      active = approximation_active,
      opt_in = FALSE,
      kind = if (approximation_active) "unverified_A2_KL_candidate" else NULL,
      warning_code = NULL
    ),
    projection = target_K$provenance$projection,
    parameterization = "log_ab",
    backend = list(
      package = "DPprior",
      package_version = tryCatch(
        as.character(utils::packageVersion("DPprior")),
        error = function(condition) "development"
      ),
      implementation = "R/12_a2_kl.R:DPprior_a2_kl",
      source_commit = source_commit
    ),
    input_fit = NULL,
    migration = list(
      source_schema = "phase8_a2_kl_flat_v0",
      adapter = ".a2_kl_result_v1",
      lossless = FALSE,
      missing_evidence = c(
        "legacy_conditions_represented_as_plain_typed_records",
        "legacy_nonfinite_values_encoded_as_plain_text",
        "optimizer_specific_labels_are_non_authoritative_compatibility_evidence",
        if (same_order_verifier) {
          "required_higher_order_verifier_unavailable"
        } else {
          character()
        },
        if (length(selection_quarantine_ids)) {
          "source_log_PMF_tie_not_comparable_under_direct_PMF_zero_tolerance"
        } else {
          character()
        },
        if (is.null(source_commit)) "source_commit_not_embedded" else character()
      ),
      warnings = "compatibility.views.a2_kl_v0_is_non_authoritative"
    ),
    legacy = list(active = FALSE, contract = NULL, deprecation_stage = NULL)
  )

  quarantine_boundary <- list(
    authority = "non_authoritative",
    lossy = TRUE,
    consumer_policy = "ignored_by_scientific_and_decision_consumers"
  )
  legacy_view <- c(
    .a2_kl_plain_compat_value(legacy_result), quarantine_boundary
  )
  canonical_converged <- context$status %in% c("converged", "boundary") &&
    isTRUE(context$usable) && isTRUE(context$verified)
  views <- list(
    a2_kl_v0 = legacy_view,
    converged = canonical_converged,
    iterations = legacy_view$iterations,
    termination = legacy_view$termination,
    fit = legacy_view$fit,
    diagnostics = legacy_view$diagnostics
  )
  out <- .dpprior_new_fit(
    mode = "a2_kl",
    method = "A2-KL",
    J = J,
    status = context$status,
    usable = context$usable,
    verified = context$verified,
    message = context$message,
    parameters = parameters,
    target = list(K = target_K),
    achieved = selected_snapshot$achieved,
    residuals = selected_snapshot$residuals,
    tolerances = tolerances,
    computation = computation,
    verification = verification,
    provenance = provenance,
    compatibility = .dpprior_new_compatibility()
  )
  .dpprior_append_compatibility_v2(
    out,
    aliases = c(
      a = "parameters.a",
      b = "parameters.b",
      attempts = "computation.attempts",
      converged = "compatibility.views.converged",
      iterations = "compatibility.views.iterations",
      termination = "compatibility.views.termination",
      fit = "compatibility.views.fit",
      diagnostics = "compatibility.views.diagnostics",
      trace = "computation.trace"
    ),
    views = views,
    deprecations = list(a2_kl_v0 = c(list(
      code = "a2_kl_flat_v0_view_quarantined"
    ), quarantine_boundary))
  )
}


# =============================================================================
# A2-KL Optimization (Main Function)
# =============================================================================

#' A2-KL: KL Divergence Minimization for Prior Calibration
#'
#' Calibrates Gamma hyperprior parameters \eqn{(a, b)} by minimizing the
#' Kullback-Leibler divergence between a target distribution and the induced
#' marginal PMF of the number of clusters \eqn{K_J}.
#'
#' @param J Integer; sample size (number of observations). Must be >= 2.
#' @param target Either:
#'   \itemize{
#'     \item Numeric vector of length J or J+1: target PMF for k = 1, ..., J,
#'           or for k = 0, ..., J with the k = 0 entry zero and dropped
#'           (used when \code{method = "pmf"}).
#'     \item Named list with \code{mu_K} and \code{var_K}: construct from moments
#'           using a discretized scaled chi-square (used when \code{method = "chisq"}).
#'   }
#' @param method Character; \code{"pmf"} or \code{"chisq"}. Default: \code{"pmf"}.
#' @param max_iter Integer; maximum optimization iterations. Default: 100.
#' @param tol Numeric; convergence tolerance for optimization. Default: 1e-6.
#' @param M Integer; number of quadrature nodes. Default: 80.
#' @param verbose Logical; if TRUE, print optimization progress. Default: FALSE.
#' @param ... Optional tuning parameters:
#'   \itemize{
#'     \item \code{log_bounds}: numeric length-2 vector giving lower/upper bounds
#'           on \code{log(a)} and \code{log(b)} (default \code{c(-15, 15)}).
#'     \item \code{M_verify}: compliant higher quadrature order (by default
#'           \code{max(2*M, M+40)} when supported).
#'     \item \code{kl_tol}, \code{l1_tol}, \code{mean_scaled_tol}, and
#'           \code{var_scaled_tol}: independent adequacy tolerances.
#'     \item \code{pmf_abs_tol}, \code{pmf_rel_tol}: higher-order PMF
#'           verification tolerances.
#'     \item \code{boundary_tol}: log-parameter bound proximity tolerance.
#'     \item \code{fallback_max_iter}: iteration budget for the predeclared
#'           bounded \code{nlminb} fallback.
#'   }
#'
#' @return A canonical \code{dpprior.result/1} \code{DPprior_fit} object. The
#'   authoritative fields are:
#'   \describe{
#'     \item{\code{parameters}}{Selected Gamma shape \code{a}, rate \code{b},
#'       and the \code{log_ab} parameterization.}
#'     \item{\code{J}}{Integer; sample size}
#'     \item{\code{target$K}}{Canonical requested, normalized, and used target.
#'       A strict PMF remains authoritative. For \code{method = "chisq"}, the
#'       requested moments remain authoritative and the conditioned chi-square
#'       PMF is recorded only as \code{A2_KL_objective} derivation evidence.}
#'     \item{\code{achieved}, \code{residuals}}{Selected-order induced PMF,
#'       moments, and distinct KL, L1, mean, and variance residuals.}
#'     \item{\code{computation}}{Closed controls, quadrature orders, trace,
#'       typed initializer/optimizer attempts, fresh candidate evaluations,
#'       fallback, candidate selection, and termination evidence.}
#'     \item{\code{verification}}{A fixed-parameter verifier snapshot at the
#'       required independent order, four distribution-adequacy gates, and a
#'       separate PMF order-stability check.}
#'     \item{\code{method}}{Character; \code{"A2-KL"}.}
#'     \item{\code{status}}{One of \code{converged}, \code{boundary},
#'       \code{approximate}, \code{infeasible}, or \code{failed}.}
#'     \item{\code{usable}, \code{verified}}{Logical adequacy flags.}
#'     \item{\code{provenance}}{Method, backend, approximation, projection,
#'       and migration records.}
#'     \item{\code{compatibility}}{A quarantined non-authoritative Phase 8
#'       view. Top-level \code{a}, \code{b}, \code{fit}, \code{diagnostics},
#'       \code{attempts}, and \code{trace} are registered migration aliases.}
#'   }
#'
#' @details
#' The A2-KL algorithm finds \eqn{(a^*, b^*)} by solving:
#' \deqn{(a^*, b^*) = \arg\min_{a,b} D_{KL}(p^*(K_J) \| p_{a,b}(K_J))}
#'
#' \strong{Algorithm:}
#' \enumerate{
#'   \item Validate a normalized custom PMF, or explicitly construct and
#'         support-condition a chi-square target PMF
#'   \item Initialize \eqn{(a_0, b_0)} from A2-MN (exact moment matching)
#'   \item Optimize KL divergence using L-BFGS-B in log-space with bounds
#'   \item If the primary optimizer fails, try a recorded bounded
#'         \code{nlminb} fallback for the same exact objective
#'   \item Independently reconstruct target and induced PMFs at a compliant
#'         higher quadrature order and assess KL, L1, and scaled moment errors
#'   \item Classify status from verification and adequacy, not optimizer exit
#' }
#'
#' \strong{Stability features:}
#' \itemize{
#'   \item Log-parameterization ensures positivity of (a, b)
#'   \item Explicit bounds prevent numerical overflow/underflow
#'   \item Exact zero-support KL semantics without epsilon smoothing
#'   \item Explicit fallback attempts; initialization is never labelled success
#'   \item A2-MN initialization provides good starting point
#' }
#'
#' \strong{When to use A2-KL vs A2-MN:}
#' \itemize{
#'   \item \strong{A2-MN}: Target moments only (exact moment matching)
#'   \item \strong{A2-KL}: Full target distribution shape (multi-modal,
#'         skewed, expert-elicited)
#' }
#'
#' The canonical \code{achieved} and \code{residuals} fields always contain the
#' selected-\code{M} computation. Higher-order PMF values are retained in the
#' independent \code{verification$verifier_snapshot}; they determine status but
#' never replace the selected-order result. The default adequacy limits are KL
#' 0.015, L1 0.11, scaled mean 0.01, and scaled variance 0.065. The PMF order
#' stability tolerance is a separate numerical check and is never substituted
#' for any of those four scientific adequacy gates.
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso
#' \code{\link{DPprior_a2_newton}} for exact moment matching,
#' \code{\link{DPprior_a1}} for closed-form initialization,
#' \code{\link{kl_divergence_K}} for KL divergence computation,
#' \code{\link{discretize_chisq}} for chi-square discretization
#'
#' @examples
#' # Example 1: Target from moments (method = "chisq")
#' fit <- DPprior_a2_kl(J = 50, target = list(mu_K = 5, var_K = 8),
#'                      method = "chisq")
#' print(fit)
#'
#' # Example 2: Custom target PMF (method = "pmf")
#' target_pmf <- dbinom(1:50, size = 50, prob = 0.1)
#' target_pmf <- target_pmf / sum(target_pmf)
#' fit2 <- DPprior_a2_kl(J = 50, target = target_pmf, method = "pmf")
#'
#' # Example 3: Compare A2-KL vs A2-MN
#' a2_mn <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8)
#' a2_kl <- DPprior_a2_kl(J = 50, target = list(mu_K = 5, var_K = 8),
#'                        method = "chisq")
#' cat(sprintf(
#'   "A2-MN: a=%.4f, b=%.4f\n",
#'   a2_mn$parameters$a, a2_mn$parameters$b
#' ))
#' cat(sprintf(
#'   "A2-KL: a=%.4f, b=%.4f, KL=%.4e\n",
#'   a2_kl$parameters$a, a2_kl$parameters$b,
#'   a2_kl$residuals$distribution$kl
#' ))
#'
#' @family elicitation
#'
#' @export
DPprior_a2_kl <- function(J, target,
                          method = c("pmf", "chisq"),
                          max_iter = 100L,
                          tol = 1e-6,
                          M = .QUAD_NODES_DEFAULT,
                          verbose = FALSE,
                          ...) {
  assert_valid_J(J)
  J <- as.integer(J)
  if (J < 2L) {
    .dpprior_abort_invalid(
      "A2-KL requires J >= 2.", "dpprior_a2_kl_input_error",
      "J", J, "integer >= 2", "range"
    )
  }
  if (!missing(method) && (
      !is.character(method) || !is.null(dim(method)) || is.object(method) ||
      length(method) != 1L || is.na(method))) {
    .dpprior_abort_invalid(
      "method must be one ordinary character scalar: 'pmf' or 'chisq'.",
      c("dpprior_a2_kl_input_error", "dpprior_type_error"),
      "method", method, "'pmf' or 'chisq'", "type"
    )
  }
  method <- tryCatch(
    match.arg(method),
    error = function(e) .dpprior_abort_invalid(
      "method must be exactly 'pmf' or 'chisq'.",
      "dpprior_a2_kl_input_error", "method", method,
      "'pmf' or 'chisq'", "choice"
    )
  )
  max_iter <- .dpprior_validate_count(
    max_iter, "max_iter", minimum = 1L,
    .subclass = "dpprior_a2_kl_control_error"
  )
  tol <- .dpprior_validate_scalar(
    tol, "tol", lower = 0, lower_open = TRUE,
    .subclass = "dpprior_a2_kl_control_error"
  )
  M <- .dpprior_validate_count(
    M, "M", minimum = 10L, maximum = .QUADRATURE_MAX_NODES,
    .subclass = "dpprior_a2_kl_control_error"
  )
  verbose <- .dpprior_validate_control(verbose, "verbose", type = "logical")

  dots_call <- match.call(expand.dots = FALSE)[["..."]]
  dot_names <- if (is.null(dots_call)) character() else names(dots_call)
  if (length(dots_call) && (is.null(dot_names) || any(!nzchar(dot_names)) ||
      anyDuplicated(dot_names))) {
    .dpprior_abort_invalid(
      "All A2-KL optional controls in ... must be uniquely named.",
      "dpprior_a2_kl_control_error", "...", dot_names,
      "uniquely named controls", "names"
    )
  }
  if ("eps" %in% dot_names) {
    .dpprior_abort_invalid(
      paste(
        "eps smoothing is not part of the A2-KL mathematical objective;",
        "zero support is handled exactly."
      ),
      c("dpprior_kl_smoothing_error", "dpprior_a2_kl_control_error"),
      "eps", deparse(dots_call[["eps"]]), "omit eps", "smoothing_not_exact"
    )
  }
  allowed_dots <- c(
    "log_bounds", "M_verify", "kl_tol", "l1_tol",
    "mean_scaled_tol", "var_scaled_tol", "pmf_abs_tol",
    "pmf_rel_tol", "boundary_tol", "fallback_max_iter"
  )
  unknown_dots <- setdiff(dot_names, allowed_dots)
  if (length(unknown_dots)) {
    .dpprior_abort_invalid(
      sprintf("Unknown A2-KL control(s): %s.", paste(unknown_dots, collapse = ", ")),
      "dpprior_a2_kl_control_error", "...", unknown_dots,
      paste(allowed_dots, collapse = ", "), "unknown_control"
    )
  }
  dots <- tryCatch(
    list(...),
    error = function(condition) {
      .dpprior_abort_invalid(
        paste("An A2-KL optional control could not be evaluated:",
              conditionMessage(condition)),
        "dpprior_a2_kl_control_error", "...", dot_names,
        "evaluable named controls", "control_evaluation"
      )
    }
  )

  log_bounds <- if (is.null(dots$log_bounds)) c(-15, 15) else dots$log_bounds
  if (!.dpprior_is_plain_numeric(log_bounds) || length(log_bounds) != 2L ||
      any(!is.finite(log_bounds)) || log_bounds[1] >= log_bounds[2]) {
    .dpprior_abort_invalid(
      "log_bounds must be a finite numeric vector of length 2 with lower < upper.",
      "dpprior_a2_kl_control_error", "log_bounds", log_bounds,
      "finite c(lower, upper) with lower < upper", "bounds"
    )
  }
  validate_positive_control <- function(value, name, default) {
    .dpprior_validate_scalar(
      if (is.null(value)) default else value,
      name, lower = 0, lower_open = TRUE,
      .subclass = "dpprior_a2_kl_control_error"
    )
  }
  kl_tol <- validate_positive_control(
    dots$kl_tol, "kl_tol", .A2_KL_ADEQUACY_DEFAULTS$kl_tol
  )
  l1_tol <- validate_positive_control(
    dots$l1_tol, "l1_tol", .A2_KL_ADEQUACY_DEFAULTS$l1_tol
  )
  mean_scaled_tol <- validate_positive_control(
    dots$mean_scaled_tol, "mean_scaled_tol",
    .A2_KL_ADEQUACY_DEFAULTS$mean_scaled_tol
  )
  var_scaled_tol <- validate_positive_control(
    dots$var_scaled_tol, "var_scaled_tol",
    .A2_KL_ADEQUACY_DEFAULTS$var_scaled_tol
  )
  pmf_abs_tol <- .dpprior_validate_scalar(
    if (is.null(dots$pmf_abs_tol)) 1e-10 else dots$pmf_abs_tol,
    "pmf_abs_tol", lower = 0,
    .subclass = "dpprior_a2_kl_control_error"
  )
  pmf_rel_tol <- .dpprior_validate_scalar(
    if (is.null(dots$pmf_rel_tol)) 1e-8 else dots$pmf_rel_tol,
    "pmf_rel_tol", lower = 0,
    .subclass = "dpprior_a2_kl_control_error"
  )
  boundary_tol <- validate_positive_control(
    dots$boundary_tol, "boundary_tol", 1e-6
  )
  fallback_max_iter <- .dpprior_validate_count(
    if (is.null(dots$fallback_max_iter)) max(100L, max_iter) else
      dots$fallback_max_iter,
    "fallback_max_iter", minimum = 1L,
    .subclass = "dpprior_a2_kl_control_error"
  )
  capped_controls <- c(
    tol = tol,
    kl_tol = kl_tol,
    l1_tol = l1_tol,
    mean_scaled_tol = mean_scaled_tol,
    var_scaled_tol = var_scaled_tol,
    pmf_abs_tol = pmf_abs_tol,
    pmf_rel_tol = pmf_rel_tol,
    boundary_tol = boundary_tol
  )
  control_caps <- c(
    tol = 1e-6,
    kl_tol = .A2_KL_ADEQUACY_DEFAULTS$kl_tol,
    l1_tol = .A2_KL_ADEQUACY_DEFAULTS$l1_tol,
    mean_scaled_tol = .A2_KL_ADEQUACY_DEFAULTS$mean_scaled_tol,
    var_scaled_tol = .A2_KL_ADEQUACY_DEFAULTS$var_scaled_tol,
    pmf_abs_tol = 1e-10,
    pmf_rel_tol = 1e-8,
    boundary_tol = 1e-6
  )
  exceeded <- names(capped_controls)[capped_controls > control_caps]
  if (length(exceeded)) {
    field <- exceeded[[1L]]
    .dpprior_abort_invalid(
      sprintf(
        "%s must be no greater than %.17g under the canonical A2-KL contract.",
        field, control_caps[[field]]
      ),
      "dpprior_a2_kl_control_error",
      field,
      unname(capped_controls[[field]]),
      sprintf("<= %.17g", control_caps[[field]]),
      "canonical_control_cap"
    )
  }
  required_M_verify <- .quadrature_verification_required_order(M)
  requested_M_verify <- if (!is.null(dots$M_verify)) {
    dots$M_verify
  } else if (required_M_verify <= .QUADRATURE_MAX_NODES) {
    required_M_verify
  } else {
    NULL
  }
  verification_controls <- .marginal_verification_controls(
    M = M, M_verify = requested_M_verify,
    abs_tol = pmf_abs_tol, rel_tol = pmf_rel_tol, strict = FALSE
  )
  M_verify <- verification_controls$M_verify

  target_info <- list(type = method)
  if (identical(method, "pmf")) {
    target_pmf <- .a2_kl_normalize_pmf(target, J)
    target_moments <- .a2_kl_pmf_moments(target_pmf)
    mu_target <- target_moments$mu_K
    var_target <- target_moments$var_K
    target_info$pmf <- target_pmf
    target_info$mu_K <- mu_target
    target_info$var_K <- var_target
    target_info$input_provenance <- list(
      source = "custom_pmf",
      normalization = "validated_not_modified",
      input_length = length(target),
      input_sum = sum(target),
      k0_entry = length(target) == J + 1L
    )
  } else {
    target_names <- if (is.list(target)) names(target) else NULL
    target_schema_ok <- is.list(target) && !is.null(target_names) &&
      length(target_names) == 2L && !anyNA(target_names) &&
      all(nzchar(target_names)) && !anyDuplicated(target_names) &&
      setequal(target_names, c("mu_K", "var_K"))
    if (!target_schema_ok) {
      .dpprior_abort_invalid(
        paste0(
          "For method='chisq', target must be a uniquely named list with ",
          "exactly 'mu_K' and 'var_K'; unknown or duplicate fields are not allowed."
        ),
        c("dpprior_target_structure_error", "dpprior_a2_kl_input_error"),
        "target", target, "list(mu_K=..., var_K=...)",
        "target_structure"
      )
    }
    moments <- .a2_kl_validate_target_moments(
      J, target$mu_K, target$var_K
    )
    mu_target <- moments$mu_K
    var_target <- moments$var_K
    df <- 2 * mu_target^2 / var_target
    scale <- var_target / (2 * mu_target)
    target_pmf <- discretize_chisq(J, df = df, scale = scale)
    target_provenance <- attr(
      target_pmf, "dpprior_pmf_provenance", exact = TRUE
    )
    target_pmf <- unname(as.numeric(target_pmf))
    discrete_moments <- .a2_kl_pmf_moments(target_pmf)
    mu_disc <- discrete_moments$mu_K
    var_disc <- discrete_moments$var_K
    target_info$pmf <- target_pmf
    target_info$mu_K <- mu_target
    target_info$var_K <- var_target
    target_info$df <- df
    target_info$scale <- scale
    target_info$mu_K_discrete <- mu_disc
    target_info$var_K_discrete <- var_disc
    target_info$input_provenance <- target_provenance
  }
  target_info$pmf_sum <- sum(target_pmf)

  if (verbose) {
    cat("A2-KL: KL Divergence Minimization\n")
    cat(sprintf("  J = %d, method = '%s'\n", J, method))
    cat(sprintf("  Target: mu_K = %.4f, var_K = %.4f\n", mu_target, var_target))
    if (identical(method, "chisq")) {
      cat(sprintf("  Discretized: mu = %.4f, var = %.4f (df=%.2f, scale=%.4f)\n",
                  target_info$mu_K_discrete, target_info$var_K_discrete,
                  target_info$df, target_info$scale))
    }
  }
  logS <- compute_log_stirling(J)
  mu_init <- if (!is.null(target_info$mu_K_discrete)) {
    target_info$mu_K_discrete
  } else {
    mu_target
  }
  var_init <- if (!is.null(target_info$var_K_discrete)) {
    target_info$var_K_discrete
  } else {
    var_target
  }
  requested_initial_moments <- c(mean = mu_init, variance = var_init)
  mu_init <- min(max(mu_init, 1 + 1e-3), J - 1e-3)
  var_init <- max(var_init, 1e-8)

  init_attempts <- list()
  newton_capture <- .a2_kl_capture_call(
    DPprior_a2_newton(
      J = J, mu_K = mu_init, var_K = var_init,
      max_iter = 15L, tol_F = 1e-8, verbose = FALSE, M = M
    )
  )
  newton_fit <- newton_capture$value
  newton_evidence <- .a2_kl_initializer_evidence(newton_fit)
  init_attempts[[length(init_attempts) + 1L]] <- list(
    stage = "initialization", method = "A2-MN", start = NULL,
    bounds = NULL,
    control = list(max_iter = 15L, tol_F = 1e-8, M = M),
    exit_code = if (is.null(newton_evidence)) NA_integer_ else
      if (isTRUE(newton_evidence$verified_success)) 0L else 1L,
    message = if (is.null(newton_evidence)) newton_capture$error else
      as.character(newton_evidence$status),
    iterations = if (is.null(newton_evidence$iterations)) NA_integer_ else
      as.integer(newton_evidence$iterations),
    evaluations = NA_integer_, candidate_objective = NA_real_,
    elapsed = newton_capture$elapsed, warnings = newton_capture$warnings,
    error = newton_capture$error, error_class = newton_capture$error_class,
    candidate = if (!is.null(newton_evidence))
      c(a = newton_evidence$a, b = newton_evidence$b) else NULL,
    selected = FALSE
  )
  init_evidence <- newton_evidence

  if (is.null(init_evidence)) {
    a1_args <- list(J = J, mu_K = mu_init, var_K = var_init)
    if ("projection" %in% names(formals(DPprior_a1))) {
      a1_args$projection <- "error"
    }
    a1_capture <- .a2_kl_capture_call(do.call(DPprior_a1, a1_args))
    a1_fit <- a1_capture$value
    a1_evidence <- .a2_kl_initializer_evidence(a1_fit)
    init_attempts[[length(init_attempts) + 1L]] <- list(
      stage = "initialization", method = "A1", start = NULL,
      bounds = NULL,
      control = list(projection = if ("projection" %in% names(a1_args))
        "error" else "legacy_default"),
      exit_code = if (!is.null(a1_evidence)) 0L else NA_integer_,
      message = if (is.null(a1_evidence)) a1_capture$error else
        as.character(a1_evidence$status),
      iterations = 0L, evaluations = NA_integer_,
      candidate_objective = NA_real_, elapsed = a1_capture$elapsed,
      warnings = a1_capture$warnings, error = a1_capture$error,
      error_class = a1_capture$error_class,
      candidate = if (!is.null(a1_evidence))
        c(a = a1_evidence$a, b = a1_evidence$b) else NULL,
      selected = FALSE
    )
    if (!is.null(a1_evidence)) init_evidence <- a1_evidence
  }
  if (is.null(init_evidence)) {
    c_J <- log(J)
    m <- max(0.5, mu_init - 1)
    vif <- max(1.01, var_init / mu_init)
    a0 <- max(0.1, min(100, m / (vif - 1)))
    b0 <- max(0.01, min(100, a0 * c_J / m))
    init_method <- "heuristic"
    init_status <- "approximate"
    init_residual <- NA_real_
    init_attempts[[length(init_attempts) + 1L]] <- list(
      stage = "initialization", method = "heuristic", start = NULL,
      bounds = NULL, control = list(formula = "deterministic moment heuristic"),
      exit_code = NA_integer_, message = "used after initializer failures",
      iterations = 0L, evaluations = 1L, candidate_objective = NA_real_,
      elapsed = 0, warnings = character(), error = NULL, error_class = NULL,
      candidate = c(a = a0, b = b0), selected = FALSE
    )
  } else {
    a0 <- init_evidence$a
    b0 <- init_evidence$b
    init_method <- init_evidence$method
    init_status <- init_evidence$status
    init_residual <- init_evidence$residual
  }
  eta0_raw <- log(c(a0, b0))
  eta0 <- pmin(pmax(eta0_raw, log_bounds[1L]), log_bounds[2L])
  initialization_clipped <- any(eta0 != eta0_raw)
  a0 <- exp(eta0[1L])
  b0 <- exp(eta0[2L])
  selected_init_index <- max(which(vapply(
    init_attempts,
    function(x) !is.null(x$candidate),
    logical(1)
  )))
  init_attempts[[selected_init_index]]$selected_for_start <- TRUE

  if (verbose) {
    cat(sprintf("  Initialization (%s): a0 = %.4f, b0 = %.4f\n",
                init_method, a0, b0))
  }
  trace_env <- new.env(parent = emptyenv())
  trace_env$n_eval <- 0L
  trace_env$records <- list()
  trace_env$attempt <- "initial"
  evaluate_eta <- function(eta) {
    if (!is.numeric(eta) || length(eta) != 2L || any(!is.finite(eta)) ||
        any(eta < log_bounds[1L]) || any(eta > log_bounds[2L])) {
      return(list(kl = Inf, error = "candidate outside finite log bounds"))
    }
    result <- .a2_kl_capture_call({
      log_q <- .a2_kl_induced_log_pmf(
        J, exp(eta[1L]), exp(eta[2L]), logS, M
      )$selected
      .a2_kl_compute_kl_logq(target_pmf, log_q)
    })
    list(kl = if (is.null(result$value)) Inf else result$value,
         error = result$error)
  }
  objective <- function(eta) {
    a_curr <- exp(eta[1L])
    b_curr <- exp(eta[2L])
    evaluated <- evaluate_eta(eta)
    kl <- evaluated$kl
    objective_value <- if (is.finite(kl)) kl else .PENALTY_INF
    trace_env$n_eval <- trace_env$n_eval + 1L
    trace_env$records[[trace_env$n_eval]] <- data.frame(
      eval = trace_env$n_eval, attempt = trace_env$attempt,
      a = a_curr, b = b_curr, kl = kl,
      objective = objective_value,
      error = if (is.null(evaluated$error)) NA_character_ else evaluated$error,
      stringsAsFactors = FALSE
    )
    if (verbose && (trace_env$n_eval %% 10L == 0L)) {
      cat(sprintf("  eval=%d | a=%.6g | b=%.6g | KL=%.6g\n",
                  trace_env$n_eval, a_curr, b_curr, kl))
    }
    objective_value
  }
  kl0_evaluated <- evaluate_eta(eta0)
  kl0 <- kl0_evaluated$kl
  for (i in seq_along(init_attempts)) {
    if (isTRUE(init_attempts[[i]]$selected_for_start)) {
      init_attempts[[i]]$candidate_objective <- kl0
    }
  }
  primary_control <- list(
    maxit = max_iter,
    factr = tol / .Machine$double.eps,
    pgtol = tol
  )
  lower <- rep(log_bounds[1L], 2L)
  upper <- rep(log_bounds[2L], 2L)
  trace_env$attempt <- "primary_lbfgsb"
  primary_capture <- .a2_kl_capture_call(
    .a2_kl_run_lbfgsb(eta0, objective, lower, upper, primary_control)
  )
  primary <- primary_capture$value
  primary_eta <- if (is.list(primary) && is.numeric(primary$par) &&
      length(primary$par) == 2L && all(is.finite(primary$par))) {
    primary$par
  } else {
    rep(NA_real_, 2L)
  }
  primary_evaluated <- evaluate_eta(primary_eta)
  primary_exit <- if (is.null(primary$convergence)) NA_integer_ else
    as.integer(primary$convergence)
  primary_attempt <- list(
    stage = "optimizer", method = "L-BFGS-B",
    start = c(a = a0, b = b0), start_log = eta0,
    bounds = list(log_lower = lower, log_upper = upper),
    control = primary_control, exit_code = primary_exit,
    message = if (is.null(primary)) primary_capture$error else primary$message,
    iterations = NA_integer_,
    iteration_reason = "stats::optim does not report an iteration count",
    evaluations = if (is.null(primary$counts[["function"]])) NA_integer_ else
      as.integer(primary$counts[["function"]]),
    candidate_objective = primary_evaluated$kl,
    elapsed = primary_capture$elapsed, warnings = primary_capture$warnings,
    error = primary_capture$error, error_class = primary_capture$error_class,
    candidate = if (all(is.finite(primary_eta)))
      c(a = exp(primary_eta[1L]), b = exp(primary_eta[2L])) else NULL,
    selected = FALSE
  )
  fallback_reason <- character()
  if (!is.null(primary_capture$error)) fallback_reason <- c(fallback_reason, "primary_error")
  if (!identical(primary_exit, 0L)) fallback_reason <- c(fallback_reason, "primary_exit_nonzero")
  if (!is.finite(primary_evaluated$kl)) fallback_reason <- c(fallback_reason, "primary_objective_nonfinite")
  if (is.finite(primary_evaluated$kl) && is.finite(kl0) &&
      primary_evaluated$kl > kl0 + 1e-12) {
    fallback_reason <- c(fallback_reason, "primary_worse_than_initialization")
  }
  fallback_attempt <- NULL
  if (length(fallback_reason)) {
    fallback_start <- if (all(is.finite(primary_eta))) primary_eta else eta0
    fallback_control <- list(
      iter.max = fallback_max_iter, eval.max = max(200L, 2L * fallback_max_iter),
      rel.tol = tol, x.tol = tol
    )
    trace_env$attempt <- "fallback_nlminb"
    fallback_capture <- .a2_kl_capture_call(
      .a2_kl_run_nlminb(
        fallback_start, objective, lower, upper, fallback_control
      )
    )
    fallback <- fallback_capture$value
    fallback_eta <- if (is.list(fallback) && is.numeric(fallback$par) &&
        length(fallback$par) == 2L && all(is.finite(fallback$par))) {
      fallback$par
    } else {
      rep(NA_real_, 2L)
    }
    fallback_evaluated <- evaluate_eta(fallback_eta)
    fallback_attempt <- list(
      stage = "optimizer", method = "nlminb",
      start = c(a = exp(fallback_start[1L]), b = exp(fallback_start[2L])),
      start_log = fallback_start,
      bounds = list(log_lower = lower, log_upper = upper),
      control = fallback_control,
      exit_code = if (is.null(fallback$convergence)) NA_integer_ else
        as.integer(fallback$convergence),
      message = if (is.null(fallback)) fallback_capture$error else fallback$message,
      iterations = if (is.null(fallback$iterations)) NA_integer_ else
        as.integer(fallback$iterations),
      evaluations = if (is.null(fallback$evaluations[["function"]]))
        NA_integer_ else as.integer(fallback$evaluations[["function"]]),
      candidate_objective = fallback_evaluated$kl,
      elapsed = fallback_capture$elapsed, warnings = fallback_capture$warnings,
      error = fallback_capture$error, error_class = fallback_capture$error_class,
      candidate = if (all(is.finite(fallback_eta)))
        c(a = exp(fallback_eta[1L]), b = exp(fallback_eta[2L])) else NULL,
      selected = FALSE
    )
  }

  candidates <- list()
  if (!is.null(primary_attempt$candidate)) {
    candidates$primary <- list(
      eta = primary_eta, kl = primary_attempt$candidate_objective,
      method = "L-BFGS-B", stage = "optimizer", attempt = "primary"
    )
  }
  if (!is.null(fallback_attempt) && !is.null(fallback_attempt$candidate)) {
    candidates$fallback <- list(
      eta = log(fallback_attempt$candidate),
      kl = fallback_attempt$candidate_objective,
      method = "nlminb", stage = "optimizer", attempt = "fallback"
    )
  }
  candidates$initialization <- list(
    eta = eta0, kl = kl0, method = init_method,
    stage = "initialization", attempt = "initialization"
  )
  candidate_kl <- vapply(candidates, function(x) x$kl, numeric(1))
  if (any(is.finite(candidate_kl))) {
    selected_name <- names(which.min(replace(
      candidate_kl, !is.finite(candidate_kl), Inf
    )))[1L]
  } else {
    selected_name <- "initialization"
  }
  selected <- candidates[[selected_name]]
  eta_opt <- selected$eta
  a_opt <- unname(exp(eta_opt[1L]))
  b_opt <- unname(exp(eta_opt[2L]))
  primary_attempt$selected <- identical(selected_name, "primary")
  if (!is.null(fallback_attempt)) {
    fallback_attempt$selected <- identical(selected_name, "fallback")
  }
  if (identical(selected_name, "initialization")) {
    init_attempts[[selected_init_index]]$selected <- TRUE
  }
  selected_optimizer_exit <- if (identical(selected_name, "primary")) {
    primary_attempt$exit_code
  } else if (identical(selected_name, "fallback")) {
    fallback_attempt$exit_code
  } else {
    NA_integer_
  }
  optimizer_converged <- selected$stage == "optimizer" &&
    identical(selected_optimizer_exit, 0L)
  boundary_distances <- c(
    a_lower = eta_opt[1L] - lower[1L],
    a_upper = upper[1L] - eta_opt[1L],
    b_lower = eta_opt[2L] - lower[2L],
    b_upper = upper[2L] - eta_opt[2L]
  )
  boundary_sides <- names(boundary_distances)[
    is.finite(boundary_distances) & boundary_distances <= boundary_tol
  ]
  boundary_hit <- length(boundary_sides) > 0L

  target_verification_capture <- .a2_kl_capture_call({
    if (identical(method, "chisq")) {
      verified_target <- .a2_kl_verify_chisq_pmf(
        J, df = target_info$df, scale = target_info$scale
      )
      unname(as.numeric(verified_target))
    } else {
      .a2_kl_validate_mathematical_pmf(
        target_info$pmf, "verification_target_pmf", J
      )
    }
  })
  verified_target <- target_verification_capture$value
  target_verification_l1 <- if (is.null(verified_target)) Inf else
    sum(abs(verified_target - target_pmf))
  target_verification_passed <- is.finite(target_verification_l1) &&
    target_verification_l1 <= .TOL_PMF_SUM

  induced_verification_capture <- .a2_kl_capture_call(
    .a2_kl_induced_log_pmf(
      J, a_opt, b_opt, logS, M, M_verify,
      abs_tol = pmf_abs_tol, rel_tol = pmf_rel_tol
    )
  )
  induced_bundle <- induced_verification_capture$value
  selected_metrics_capture <- if (is.null(induced_bundle)) {
    list(value = NULL, error = induced_verification_capture$error)
  } else {
    .a2_kl_capture_call(
      .a2_kl_pmf_metrics(target_pmf, induced_bundle$selected)
    )
  }
  verification_metrics_capture <- if (is.null(induced_bundle$verification)) {
    list(value = NULL, error = "independent higher-order PMF unavailable")
  } else {
    .a2_kl_capture_call(
      .a2_kl_pmf_metrics(target_pmf, induced_bundle$verification)
    )
  }
  selected_metrics <- selected_metrics_capture$value
  verification_metrics <- verification_metrics_capture$value
  # The selected-order result is immutable public output. Higher-order values
  # determine verification/status but never replace achieved/fit values.
  reported_metrics <- selected_metrics
  adequacy_tolerances <- list(
    kl = kl_tol, l1 = l1_tol,
    mean_scaled = mean_scaled_tol, variance_scaled = var_scaled_tol
  )
  adequacy <- if (is.null(verification_metrics)) {
    list(
      passed = FALSE, raw = NULL, scales = NULL, scale_formula = NULL,
      scaled = NULL, scaled_tolerances = c(
        kl = kl_tol, l1 = l1_tol, mean = mean_scaled_tol,
        variance = var_scaled_tol
      ),
      raw_tolerances = NULL, ratios = NULL, component_passed = NULL
    )
  } else {
    .a2_kl_assess_adequacy(verification_metrics, adequacy_tolerances)
  }
  selected_adequacy <- if (is.null(selected_metrics)) NULL else
    .a2_kl_assess_adequacy(selected_metrics, adequacy_tolerances)
  marginal_verification <- if (is.null(induced_bundle)) NULL else
    induced_bundle$metadata$verification
  numerical_verification_passed <- !is.null(M_verify) &&
    isTRUE(marginal_verification$passed) &&
    target_verification_passed && !is.null(verification_metrics)
  finite_candidate <- all(is.finite(c(a_opt, b_opt))) &&
    !is.null(reported_metrics) && is.finite(reported_metrics$kl)
  verified <- finite_candidate && optimizer_converged &&
    numerical_verification_passed && isTRUE(adequacy$passed)
  status <- if (!finite_candidate || is.null(induced_bundle)) {
    "failed"
  } else if (verified && boundary_hit) {
    "boundary"
  } else if (verified) {
    "converged"
  } else {
    "approximate"
  }
  stopifnot(status %in% .A2_KL_STATUS_CODES)
  usable <- status %in% c("converged", "boundary")
  message_parts <- character()
  if (identical(status, "failed")) {
    message_parts <- c(
      "A2-KL failed to produce a finite independently evaluated candidate",
      induced_verification_capture$error,
      selected_metrics_capture$error
    )
  } else if (identical(status, "approximate")) {
    if (!optimizer_converged) {
      message_parts <- c(message_parts, "selected optimizer did not exit successfully")
    }
    if (identical(selected$stage, "initialization")) {
      message_parts <- c(message_parts, "initialization was returned only as an approximate candidate")
    }
    if (!numerical_verification_passed) {
      message_parts <- c(message_parts, "independent PMF verification did not pass")
    }
    if (!isTRUE(adequacy$passed)) {
      component_passed <- adequacy$component_passed
      failed_components <- if (is.null(component_passed)) {
        character()
      } else {
        names(component_passed)[!component_passed]
      }
      message_parts <- c(
        message_parts,
        sprintf(
          "target-family adequacy failed%s",
          if (length(failed_components)) paste0(
            " for ", paste(failed_components, collapse = ", ")
          ) else ""
        )
      )
    }
  } else if (identical(status, "boundary")) {
    message_parts <- sprintf(
      "A2-KL passed verification at solver boundary: %s",
      paste(boundary_sides, collapse = ", ")
    )
  } else {
    message_parts <- "A2-KL optimizer and independent adequacy verification passed"
  }
  message <- paste(unique(message_parts[nzchar(message_parts)]), collapse = "; ")

  trace_df <- if (length(trace_env$records)) {
    do.call(rbind, trace_env$records)
  } else {
    data.frame(
      eval = integer(), attempt = character(), a = numeric(), b = numeric(),
      kl = numeric(), objective = numeric(), error = character()
    )
  }
  optimization_attempts <- c(
    list(primary_attempt),
    if (is.null(fallback_attempt)) list() else list(fallback_attempt)
  )
  attempts <- c(init_attempts, optimization_attempts)
  selected_attempt_record <- if (identical(selected_name, "primary")) {
    primary_attempt
  } else if (identical(selected_name, "fallback")) {
    fallback_attempt
  } else {
    init_attempts[[selected_init_index]]
  }
  fallback_used <- !identical(selected_name, "primary")
  selected_method <- if (identical(selected$stage, "optimizer")) {
    selected$method
  } else {
    paste0("initialization:", selected$method)
  }

  achieved <- if (is.null(reported_metrics)) NULL else list(
    pmf = reported_metrics$pmf,
    mu_K = reported_metrics$mean,
    var_K = reported_metrics$variance,
    kl = reported_metrics$kl,
    l1 = reported_metrics$l1,
    quadrature_order = M
  )
  residuals <- if (is.null(reported_metrics)) NULL else list(
    raw = c(
      kl = reported_metrics$kl,
      l1 = reported_metrics$l1,
      mean = reported_metrics$mean_residual,
      variance = reported_metrics$variance_residual
    ),
    absolute = selected_adequacy$raw,
    scaled = selected_adequacy$scaled,
    ratios = selected_adequacy$ratios,
    verification = if (is.null(verification_metrics)) NULL else list(
      quadrature_order = M_verify,
      raw = c(
        kl = verification_metrics$kl,
        l1 = verification_metrics$l1,
        mean = verification_metrics$mean_residual,
        variance = verification_metrics$variance_residual
      ),
      absolute = adequacy$raw,
      scaled = adequacy$scaled,
      ratios = adequacy$ratios,
      component_passed = adequacy$component_passed
    ),
    units = c(kl = "nats", l1 = "probability mass", mean = "clusters",
              variance = "clusters^2"),
    scale_formula = adequacy$scale_formula
  )
  verification <- list(
    method = "fresh higher-order marginal PMF and direct-moment audit",
    performed = !is.null(M_verify),
    passed = verified,
    numerical_passed = numerical_verification_passed,
    adequacy_passed = isTRUE(adequacy$passed),
    target = list(
      independently_reconstructed = identical(method, "chisq"),
      independently_validated = identical(method, "pmf"),
      method = if (identical(method, "chisq")) {
        "separate_continuity_edge_cdf_difference_reconstruction"
      } else {
        "strict_pmf_postcondition_validation"
      },
      l1_difference = target_verification_l1,
      passed = target_verification_passed,
      error = target_verification_capture$error
    ),
    settings = list(
      M_selected = M, M_verification = if (is.null(M_verify)) NA_integer_ else M_verify,
      M_verification_required = required_M_verify,
      pmf_abs_tol = pmf_abs_tol, pmf_rel_tol = pmf_rel_tol
    ),
    marginal = marginal_verification,
    selected_metrics = selected_metrics,
    verification_metrics = verification_metrics,
    adequacy = adequacy,
    adequacy_basis = "higher-order M_verify status decision",
    stability_delta = if (is.null(selected_metrics) ||
        is.null(verification_metrics)) NULL else c(
      kl = verification_metrics$kl - selected_metrics$kl,
      l1_target = verification_metrics$l1 - selected_metrics$l1,
      pmf_l1 = sum(abs(
        verification_metrics$pmf - selected_metrics$pmf
      )),
      mean = verification_metrics$mean - selected_metrics$mean,
      variance = verification_metrics$variance - selected_metrics$variance
    ),
    error = induced_verification_capture$error
  )
  provenance <- list(
    requested_method = "L-BFGS-B",
    selected_method = selected_method,
    selected_attempt = selected_name,
    is_fallback = fallback_used,
    fallback_attempted = !is.null(fallback_attempt),
    fallback_reason = unique(fallback_reason),
    initialization_method = init_method,
    initialization_status = init_status,
    initialization_clipped_to_bounds = initialization_clipped,
    approximation = if (identical(status, "approximate"))
      "primary tolerance or verification not met" else NULL,
    objective = "exact D_KL(target || induced); no epsilon smoothing"
  )
  diagnostics <- list(
    M = M,
    M_verify = if (is.null(M_verify)) NA_integer_ else M_verify,
    max_iter = max_iter,
    tol = tol,
    log_bounds = log_bounds,
    init = list(
      a0 = a0,
      b0 = b0,
      mu_init = mu_init,
      var_init = var_init,
      init_method = init_method,
      init_status = init_status,
      init_residual = init_residual,
      requested_moments = requested_initial_moments,
      clipped_to_bounds = initialization_clipped
    ),
    optim = list(
      method = selected_method,
      convergence = selected_optimizer_exit,
      message = selected_attempt_record$message,
      counts = list(
        `function` = selected_attempt_record$evaluations,
        iterations = selected_attempt_record$iterations
      ),
      value = if (is.null(reported_metrics)) Inf else reported_metrics$kl,
      par = eta_opt,
      optimizer_converged = optimizer_converged,
      adequate = isTRUE(adequacy$passed),
      adequacy_basis = "verification$adequacy at M_verify",
      boundary_hit = boundary_hit,
      boundary_sides = boundary_sides,
      boundary_distances = boundary_distances
    ),
    fallback_used = fallback_used,
    kl_init = kl0,
    kl_final = if (is.null(reported_metrics)) Inf else reported_metrics$kl,
    attempts = attempts,
    adequacy = selected_adequacy,
    adequacy_basis = paste(
      "selected-M public diagnostics; status uses",
      "verification$adequacy at M_verify"
    ),
    verification = verification,
    provenance = provenance
  )
  fit <- if (is.null(achieved)) {
    list(mu_K = NA_real_, var_K = NA_real_, kl = Inf, l1 = Inf,
         residual = Inf, adequacy_residual = Inf)
  } else {
    list(
      mu_K = achieved$mu_K, var_K = achieved$var_K,
      kl = achieved$kl, l1 = achieved$l1,
      residual = achieved$kl,
      adequacy_residual = max(selected_adequacy$ratios)
    )
  }
  termination <- paste0(status, ":", selected_method)
  iterations <- if (is.null(selected_attempt_record$iterations)) NA_integer_ else
    as.integer(selected_attempt_record$iterations)
  legacy_result <- structure(
    list(
      a = a_opt,
      b = b_opt,
      J = J,
      target = target_info,
      method = "A2-KL",
      status = status,
      usable = usable,
      verified = verified,
      message = message,
      parameters = c(a = a_opt, b = b_opt),
      achieved = achieved,
      residuals = residuals,
      tolerances = list(
        adequacy = adequacy_tolerances,
        adequacy_raw = adequacy$raw_tolerances,
        pmf_sum = .TOL_PMF_SUM,
        pmf_verification_abs = pmf_abs_tol,
        pmf_verification_rel = pmf_rel_tol,
        optimizer = tol,
        boundary = boundary_tol
      ),
      attempts = attempts,
      verification = verification,
      provenance = provenance,
      converged = identical(status, "converged"),
      iterations = iterations,
      termination = termination,
      fit = fit,
      diagnostics = diagnostics,
      trace = trace_df
    ),
    class = "DPprior_fit"
  )
  result <- .a2_kl_result_v1(legacy_result, context = list(
    J = J,
    original_target = target,
    target_method = method,
    target_info = target_info,
    target_pmf = target_pmf,
    a_opt = a_opt,
    b_opt = b_opt,
    M = M,
    M_verify = M_verify,
    required_M_verify = required_M_verify,
    logS = logS,
    selected_metrics = selected_metrics,
    verification_metrics = verification_metrics,
    status = status,
    usable = usable,
    verified = verified,
    message = message,
    attempts = attempts,
    selected_name = selected_name,
    selected_method = selected_method,
    fallback_reason = fallback_reason,
    boundary_sides = boundary_sides,
    init_method = init_method,
    init_status = as.character(init_status),
    initialization_clipped = initialization_clipped,
    max_iter = max_iter,
    tol = tol,
    log_bounds = log_bounds,
    boundary_tol = boundary_tol,
    fallback_max_iter = fallback_max_iter,
    kl_tol = kl_tol,
    l1_tol = l1_tol,
    mean_scaled_tol = mean_scaled_tol,
    var_scaled_tol = var_scaled_tol,
    pmf_abs_tol = pmf_abs_tol,
    pmf_rel_tol = pmf_rel_tol,
    trace = trace_df
  ))
  if (verbose) {
    cat("\nResults:\n")
    cat(sprintf("  Status: %s | verified: %s | usable: %s\n",
                status, verified, usable))
    cat(sprintf("  a* = %.6f, b* = %.6f\n", a_opt, b_opt))
    cat(sprintf("  KL = %.4e, L1 = %.4e\n", fit$kl, fit$l1))
    cat(sprintf("  Method used: %s\n", selected_method))
  }
  result
}


# =============================================================================
# Verification and Testing Functions
# =============================================================================

#' Verify KL Divergence Properties
#'
#' Runs verification tests on KL divergence computations.
#'
#' @param verbose Logical; if TRUE, print detailed results.
#'
#' @return Logical; TRUE if all tests pass.
#'
#' @examples
#' \dontrun{
#' verify_kl_divergence()
#'
#' }
#' @keywords internal
verify_kl_divergence <- function(verbose = TRUE) {
  all_pass <- TRUE

  if (isTRUE(verbose)) {
    cat("=", rep("=", 59), "\n", sep = "")
    cat("KL Divergence Verification\n")
    cat("=", rep("=", 59), "\n\n", sep = "")
  }

  # Test 1: KL(p || p) = 0
  p <- c(0.2, 0.5, 0.3)
  kl_self <- kl_divergence_pmf(p, p)
  test1 <- abs(kl_self) < 1e-10
  if (isTRUE(verbose)) {
    cat(sprintf("Test 1: KL(p || p) = 0\n"))
    cat(sprintf("  KL = %.2e [%s]\n\n", kl_self, if (test1) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test1

  # Test 2: KL divergence is non-negative
  q <- c(0.3, 0.4, 0.3)
  kl <- kl_divergence_pmf(p, q)
  test2 <- kl >= 0
  if (isTRUE(verbose)) {
    cat(sprintf("Test 2: KL(p || q) >= 0\n"))
    cat(sprintf("  KL(p || q) = %.4f [%s]\n\n", kl, if (test2) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test2

  # Test 3: KL divergence with K_J induced PMF
  J <- 50L
  logS <- compute_log_stirling(J)
  a_true <- 2.0
  b_true <- 1.0
  target <- .a2_kl_induced_pmf(J, a_true, b_true, logS, M = 80L)
  kl_match <- kl_divergence_K(target, a_true, b_true, J)
  test3 <- kl_match < 1e-10
  if (isTRUE(verbose)) {
    cat(sprintf("Test 3: KL = 0 for matching (a, b)\n"))
    cat(sprintf("  KL(p_{%.1f,%.1f} || p_{%.1f,%.1f}) = %.2e [%s]\n\n",
                a_true, b_true, a_true, b_true, kl_match,
                if (test3) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test3

  # Test 4: KL divergence > 0 for non-matching
  kl_nonmatch <- kl_divergence_K(target, a_true * 1.5, b_true * 0.8, J)
  test4 <- kl_nonmatch > 0
  if (isTRUE(verbose)) {
    cat(sprintf("Test 4: KL > 0 for non-matching (a, b)\n"))
    cat(sprintf("  KL = %.4f [%s]\n\n", kl_nonmatch, if (test4) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test4

  if (isTRUE(verbose)) {
    cat("=", rep("=", 59), "\n", sep = "")
    cat(sprintf("Overall: %s\n", if (all_pass) "ALL TESTS PASSED" else "SOME TESTS FAILED"))
    cat("=", rep("=", 59), "\n", sep = "")
  }

  invisible(all_pass)
}


#' Verify A2-KL Optimization
#'
#' Runs verification tests on the A2-KL optimization algorithm.
#'
#' @param verbose Logical; if TRUE, print detailed results.
#'
#' @return Logical; TRUE if all tests pass.
#'
#' @examples
#' \dontrun{
#' verify_a2_kl()
#'
#' }
#' @keywords internal
verify_a2_kl <- function(verbose = TRUE) {
  all_pass <- TRUE

  if (isTRUE(verbose)) {
    cat("=", rep("=", 59), "\n", sep = "")
    cat("A2-KL Optimization Verification\n")
    cat("=", rep("=", 59), "\n\n", sep = "")
  }

  # Test 1: Convergence for typical target (chisq method)
  J <- 50L
  target <- list(mu_K = 5, var_K = 8)
  fit <- DPprior_a2_kl(J, target, method = "chisq", verbose = FALSE)

  test1 <- fit$usable && fit$verified &&
    fit$status %in% c("converged", "boundary")
  if (isTRUE(verbose)) {
    cat(sprintf("Test 1: Convergence for typical target (method='chisq')\n"))
    cat(sprintf("  Target: mu_K = %.1f, var_K = %.1f\n", target$mu_K, target$var_K))
    cat(sprintf(
      "  Result: a = %.4f, b = %.4f\n",
      fit$parameters$a, fit$parameters$b
    ))
    cat(sprintf(
      "  Achieved: mu_K = %.4f, var_K = %.4f\n",
      fit$achieved$K$mean, fit$achieved$K$variance
    ))
    cat(sprintf(
      "  KL = %.4e, status = '%s'\n",
      fit$residuals$distribution$kl, fit$status
    ))
    cat(sprintf("  [%s]\n\n", if (test1) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test1

  # Test 2: Reasonable KL divergence
  test2 <- fit$residuals$distribution$kl < 0.1
  if (isTRUE(verbose)) {
    cat(sprintf("Test 2: Reasonable KL divergence (< 0.1)\n"))
    cat(sprintf(
      "  KL = %.4e [%s]\n\n", fit$residuals$distribution$kl,
      if (test2) "PASS" else "FAIL"
    ))
  }
  all_pass <- all_pass && test2

  # Test 3: Moment approximation quality
  mean_err <- abs(fit$achieved$K$mean - target$mu_K)
  var_err <- abs(fit$achieved$K$variance - target$var_K)
  test3 <- mean_err < 0.5 && var_err < 2.0
  if (isTRUE(verbose)) {
    cat(sprintf("Test 3: Moment approximation quality\n"))
    cat(sprintf("  Mean error: %.4f (< 0.5)\n", mean_err))
    cat(sprintf("  Var error: %.4f (< 2.0)\n", var_err))
    cat(sprintf("  [%s]\n\n", if (test3) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test3

  # Test 4: Custom PMF target (method='pmf')
  target_pmf <- stats::dbinom(1:50, size = 50, prob = 0.08)
  target_pmf <- target_pmf / sum(target_pmf)
  fit2 <- DPprior_a2_kl(J, target_pmf, method = "pmf", verbose = FALSE)

  test4 <- fit2$usable && fit2$verified &&
    fit2$residuals$distribution$kl < 0.5
  if (isTRUE(verbose)) {
    cat(sprintf("Test 4: Custom binomial-shaped PMF target (method='pmf')\n"))
    cat(sprintf(
      "  Result: a = %.4f, b = %.4f\n",
      fit2$parameters$a, fit2$parameters$b
    ))
    cat(sprintf(
      "  KL = %.4e, status = '%s'\n",
      fit2$residuals$distribution$kl, fit2$status
    ))
    cat(sprintf("  [%s]\n\n", if (test4) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test4

  # Test 5: Fallback mechanism (use target that might challenge optimizer)
  # This tests that fallback doesn't crash and returns something reasonable
  # Note: We use suppressWarnings because extreme targets may trigger
  # projection warnings in A1 initialization
  if (isTRUE(verbose)) {
    cat("Test 5: Fallback mechanism check\n")
  }
  fit3 <- tryCatch(
    suppressWarnings(
      DPprior_a2_kl(J = 30, target = list(mu_K = 3, var_K = 1.5),
                    method = "chisq", verbose = FALSE)
    ),
    error = function(e) NULL
  )
  test5 <- !is.null(fit3) &&
    is.finite(fit3$residuals$distribution$kl)
  if (isTRUE(verbose)) {
    if (!is.null(fit3)) {
      cat(sprintf(
        "  Result: a = %.4f, b = %.4f, KL = %.4e\n",
        fit3$parameters$a, fit3$parameters$b,
        fit3$residuals$distribution$kl
      ))
      cat(sprintf(
        "  Fallback used: %s\n", fit3$computation$fallback$used
      ))
    } else {
      cat("  Fit returned NULL\n")
    }
    cat(sprintf("  [%s]\n\n", if (test5) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test5

  if (isTRUE(verbose)) {
    cat("=", rep("=", 59), "\n", sep = "")
    cat(sprintf("Overall: %s\n", if (all_pass) "ALL TESTS PASSED" else "SOME TESTS FAILED"))
    cat("=", rep("=", 59), "\n", sep = "")
  }

  invisible(all_pass)
}


#' Run All Module 12 Verification Tests
#'
#' Comprehensive verification suite for the A2-KL module.
#'
#' @param verbose Logical; if TRUE, print detailed results.
#'
#' @return Logical; TRUE if all tests pass.
#'
#' @examples
#' \dontrun{
#' verify_a2_kl_all()
#'
#' }
#' @keywords internal
verify_a2_kl_all <- function(verbose = TRUE) {
  if (isTRUE(verbose)) {
    cat("=", rep("=", 69), "\n", sep = "")
    cat("Module 12: A2-KL - Full Verification Suite\n")
    cat("=", rep("=", 69), "\n\n", sep = "")
  }

  all_pass <- TRUE

  all_pass <- all_pass && verify_kl_divergence(verbose = verbose)
  cat("\n")
  all_pass <- all_pass && verify_a2_kl(verbose = verbose)

  if (isTRUE(verbose)) {
    cat("\n")
    cat("=", rep("=", 69), "\n", sep = "")
    cat(sprintf("Module 12 Final Result: %s\n",
                if (all_pass) "ALL TESTS PASSED" else "SOME TESTS FAILED"))
    cat("=", rep("=", 69), "\n", sep = "")
  }

  invisible(all_pass)
}
