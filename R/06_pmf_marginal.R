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
# Some older drafts report E[K_J] ~ 10.23 for J=50 and
# alpha ~ Gamma(1.5, 0.5) under the shape-rate parameterization.
# Direct numerical integration of the conditional mean returns
# E[K_50] ~ 8.3555. If you see a mismatch with any hard-coded
# numbers, treat those as potentially stale and rely on the internal
# consistency checks (PMF <-> moments) implemented in Modules 05-06.
#
# Author: JoonHo Lee (jlee296@ua.edu)
# Date: December 2025
# Part of: DPprior R Package
# Reference: Lee (2026), Sections 2--3
# Dependencies: Module 00 (constants, logsumexp),
#               Module 02 (quadrature),
#               Module 04 (conditional PMF)
# =============================================================================


# =============================================================================
# Internal: Log-Space Marginal PMF (Core Algorithm)
# =============================================================================

.marginal_pmf_abort <- function(message, argument, value, expected, code) {
  .dpprior_abort_invalid(
    message,
    c("dpprior_marginal_pmf_error", "dpprior_numerical_error"),
    argument, value, expected, code
  )
}


# Compute one fixed-order log mixture. This helper never substitutes another
# distribution: malformed conditional rows, quadrature weights, or mixed mass
# fail with a typed condition.
.marginal_log_pmf_fixed <- function(J, a, b, logS, M) {
  quad <- build_gamma_quadrature(a, b, M)
  alpha_nodes <- quad$alpha_nodes
  weights <- quad$weights_normalized
  if (length(alpha_nodes) != M || length(weights) != M ||
      any(!is.finite(alpha_nodes)) || any(alpha_nodes <= 0) ||
      any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) {
    .marginal_pmf_abort(
      "marginal PMF quadrature returned invalid nodes or weights",
      "quadrature", list(alpha_nodes = alpha_nodes, weights = weights),
      "positive finite nodes and non-negative finite weights",
      "invalid_quadrature"
    )
  }
  log_weights <- suppressWarnings(log(weights))

  logpmf_matrix <- matrix(-Inf, nrow = M, ncol = J + 1L)
  for (m in seq_len(M)) {
    raw <- log_pmf_K_given_alpha(J, alpha_nodes[[m]], logS)
    if (!is.numeric(raw) || length(raw) != J + 1L || anyNA(raw) ||
        !is.infinite(raw[[1L]]) || raw[[1L]] > 0 ||
        any(!is.finite(raw[-1L]))) {
      .marginal_pmf_abort(
        sprintf("conditional log-PMF is invalid at quadrature node %d", m),
        "conditional_log_pmf", raw,
        sprintf("length %d with structural -Inf at k=0 and finite k=1,...,%d",
                J + 1L, J),
        "invalid_conditional_pmf"
      )
    }
    conditional_log_mass <- logsumexp_vec(raw[-1L])
    if (!is.finite(conditional_log_mass)) {
      .marginal_pmf_abort(
        sprintf("conditional PMF has non-finite mass at quadrature node %d", m),
        "conditional_log_mass", conditional_log_mass,
        "finite normalizing constant", "invalid_conditional_mass"
      )
    }
    logpmf_matrix[m, -1L] <- raw[-1L] - conditional_log_mass
  }

  logp <- rep(-Inf, J + 1L)
  for (k in seq_len(J)) {
    mixed <- logsumexp_vec(log_weights + logpmf_matrix[, k + 1L])
    if (!is.finite(mixed)) {
      .marginal_pmf_abort(
        sprintf("marginal log-PMF is non-finite at k=%d", k),
        "log_pmf", mixed, "finite mixed log probability",
        "nonfinite_mixture"
      )
    }
    logp[[k + 1L]] <- mixed
  }

  mixed_log_mass <- logsumexp_vec(logp[-1L])
  if (!is.finite(mixed_log_mass)) {
    .marginal_pmf_abort(
      "marginal PMF has non-finite total mass",
      "log_mass", mixed_log_mass, "finite positive total mass",
      "invalid_marginal_mass"
    )
  }
  logp[-1L] <- logp[-1L] - mixed_log_mass
  logp[[1L]] <- -Inf
  probability_sum <- sum(exp(logp[-1L]))
  if (!is.finite(probability_sum) || probability_sum <= 0) {
    .marginal_pmf_abort(
      "marginal PMF cannot be represented with finite positive mass",
      "probability_sum", probability_sum, "finite positive probability mass",
      "invalid_probability_mass"
    )
  }

  list(
    logp = logp,
    log_normalization = mixed_log_mass,
    probability_sum = probability_sum,
    underflow_zero_count = sum(exp(logp[-1L]) == 0),
    node_rule = .quadrature_metadata(quad)
  )
}


.marginal_quantile_from_cdf <- function(p, cdf, J) {
  vapply(p, function(probability) {
    if (probability == 0) {
      return(1L)
    }
    if (probability == 1) {
      return(as.integer(J))
    }
    index <- which(cdf >= probability)
    if (!length(index)) {
      .dpprior_abort_invalid(
        "marginal CDF does not attain the requested probability",
        "dpprior_marginal_quantile_error", "p", probability,
        "probability attained on support 1:J", "cdf_inversion"
      )
    }
    quantile <- as.integer(index[[1L]] - 1L)
    if (quantile < 1L || quantile > J) {
      .dpprior_abort_invalid(
        "marginal quantile fell outside support 1:J",
        "dpprior_marginal_quantile_error", "quantile", quantile,
        sprintf("integer in [1,%d]", J), "support"
      )
    }
    quantile
  }, integer(1))
}


# Build an exact-endpoint CDF without creating a final-step decrease when the
# floating-point PMF sum is slightly larger than one. Rescaling the cumulative
# support mass by its positive total preserves monotonicity; overwriting only
# the last entry does not.
.marginal_cdf_from_pmf <- function(pmf) {
  support_cdf <- cumsum(pmf[-1L])
  total <- support_cdf[[length(support_cdf)]]
  if (!is.finite(total) || total <= 0) {
    .dpprior_abort_invalid(
      "marginal CDF requires finite positive support mass",
      c("dpprior_marginal_cdf_error", "dpprior_numerical_error"),
      "pmf", pmf, "finite positive support mass", "cdf_mass"
    )
  }
  support_cdf <- support_cdf / total
  support_cdf[[length(support_cdf)]] <- 1
  cdf <- c(0, support_cdf)
  if (any(!is.finite(cdf)) || any(diff(cdf) < 0) ||
      cdf[[1L]] != 0 || cdf[[length(cdf)]] != 1) {
    .dpprior_abort_invalid(
      "marginal CDF violated finite, monotone, or endpoint requirements",
      c("dpprior_marginal_cdf_error", "dpprior_numerical_error"),
      "cdf", cdf, "finite non-decreasing CDF with endpoints 0 and 1",
      "cdf_contract"
    )
  }
  cdf
}


# Compute the first two moments represented by a probability vector. This is
# deliberately separate from .marginal_moments_fixed(): agreement between the
# Stirling-number PMF path and the conditional-moment quadrature path is part
# of the convergence contract, rather than an assumption.
.marginal_probability_moments <- function(probability, J) {
  k <- 0:J
  mean <- .quadrature_weighted_sum(probability, k)
  variance <- .quadrature_weighted_sum(
    probability, (k - mean)^2
  )
  list(mean = as.numeric(mean), var = as.numeric(variance))
}


.marginal_scalar_agreement <- function(
    observed, reference, abs_tol, rel_tol) {
  difference <- abs(observed - reference)
  tolerance <- abs_tol + rel_tol * max(abs(observed), abs(reference))
  list(
    pmf = observed,
    direct = reference,
    difference = difference,
    tolerance = tolerance,
    passed = is.finite(difference) && is.finite(tolerance) &&
      difference <= tolerance
  )
}


.marginal_pmf_moment_audit <- function(
    probability, direct, J, abs_tol, rel_tol) {
  represented <- .marginal_probability_moments(probability, J)
  mean <- .marginal_scalar_agreement(
    represented$mean, direct$mean, abs_tol, rel_tol
  )
  variance <- .marginal_scalar_agreement(
    represented$var, direct$var, abs_tol, rel_tol
  )
  list(
    passed = isTRUE(mean$passed) && isTRUE(variance$passed),
    mean = mean,
    variance = variance
  )
}


# Add a selected-versus-verification audit for an integer-valued functional.
# The selected value remains the return value. A disagreement can only
# downgrade provenance; it never causes verification-order substitution.
.marginal_record_discrete_verification <- function(
    metadata, component, selected, verification,
    probabilities = NULL) {
  performed <- !is.null(verification)
  selected_values <- unname(as.integer(selected))
  verification_values <- if (performed) {
    unname(as.integer(verification))
  } else {
    NULL
  }
  passed <- if (performed) {
    identical(selected_values, verification_values)
  } else {
    NA
  }

  if (is.null(metadata$discrete_verification)) {
    metadata$discrete_verification <- list()
  }
  metadata$discrete_verification[[component]] <- list(
    performed = performed,
    passed = passed,
    probabilities = if (is.null(probabilities)) {
      NULL
    } else {
      unname(probabilities)
    },
    selected = selected_values,
    verification = verification_values
  )

  if (identical(passed, FALSE)) {
    if (identical(metadata$status, "converged")) {
      metadata$reason <- paste0("discrete_", component, "_disagreement")
    }
    metadata$status <- "approximate"
  }
  metadata
}


.marginal_enforce_discrete_verification <- function(
    metadata, strict, context) {
  failed <- names(Filter(
    function(record) identical(record$passed, FALSE),
    metadata$discrete_verification
  ))
  if (isTRUE(strict) && length(failed)) {
    .dpprior_abort_invalid(
      sprintf(
        "%s changed at the verification quadrature order (%s)",
        context, paste(failed, collapse = ", ")
      ),
      "dpprior_marginal_convergence_error", "M_verify",
      metadata$verification$M_verification,
      "stable selected and verification discrete summaries",
      paste0("discrete_", paste(failed, collapse = "_"), "_disagreement")
    )
  }
  invisible(metadata)
}

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
#' @param M_verify Optional independent quadrature order satisfying the
#'   package verification rule (at least \code{max(2*M, M+40)}, within the
#'   supported ceiling). When supplied, the selected PMF is compared with a
#'   fresh higher-order PMF.
#' @param abs_tol,rel_tol Non-negative absolute and relative tolerances for the
#'   selected-versus-verification L1 discrepancy and the PMF-versus-direct-
#'   moment checks at both quadrature orders.
#' @param strict Logical; require successful higher-order verification when
#'   \code{TRUE}; otherwise return the selected PMF with explicit approximation
#'   metadata.
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
#' A verified result must pass both the selected-versus-verification PMF L1
#' budget and an independent cross-representation audit: at each order, the
#' mean and variance represented by the PMF must agree with the direct
#' marginal-moment quadrature. Thus two mutually agreeing but scientifically
#' incorrect PMFs cannot be classified as converged.
#'
#' @keywords internal
log_pmf_K_marginal <- function(
    J, a, b, logS, M = .QUAD_NODES_DEFAULT,
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
  .check_logS_size(J, logS)
  J <- as.integer(J)
  controls <- .marginal_verification_controls(
    M, M_verify, abs_tol, rel_tol, strict
  )
  selected <- .marginal_log_pmf_fixed(J, a, b, logS, controls$M)

  selected_probability <- NULL
  selected_moment_audit <- NULL
  verification <- NULL
  verification_probability <- NULL
  verification_moment_audit <- NULL
  l1_difference <- NA_real_
  max_difference <- NA_real_
  tolerance <- NA_real_
  l1_passed <- NA
  verification_passed <- NA
  reason <- "fixed_order_unverified"
  if (!is.null(controls$M_verify)) {
    selected_probability <- exp(selected$logp)
    selected_direct_moments <- .marginal_moments_fixed(
      J, a, b, controls$M
    )
    selected_moment_audit <- .marginal_pmf_moment_audit(
      selected_probability, selected_direct_moments, J,
      controls$abs_tol, controls$rel_tol
    )
    verification <- .marginal_log_pmf_fixed(
      J, a, b, logS, controls$M_verify
    )
    verification_probability <- exp(verification$logp)
    verification_direct_moments <- .marginal_moments_fixed(
      J, a, b, controls$M_verify
    )
    verification_moment_audit <- .marginal_pmf_moment_audit(
      verification_probability, verification_direct_moments, J,
      controls$abs_tol, controls$rel_tol
    )
    l1_difference <- sum(abs(selected_probability - verification_probability))
    max_difference <- max(abs(
      selected_probability - verification_probability
    ))
    tolerance <- controls$abs_tol + controls$rel_tol
    l1_passed <- is.finite(l1_difference) &&
      is.finite(max_difference) && l1_difference <= tolerance
    moment_audits_passed <- isTRUE(selected_moment_audit$passed) &&
      isTRUE(verification_moment_audit$passed)
    verification_passed <- isTRUE(l1_passed) && moment_audits_passed
    reason <- if (verification_passed) {
      "higher_order_agreement"
    } else if (!isTRUE(l1_passed) && !moment_audits_passed) {
      "higher_order_and_pmf_moment_disagreement"
    } else if (!isTRUE(l1_passed)) {
      "higher_order_disagreement"
    } else {
      "pmf_moment_disagreement"
    }
  }

  status <- if (isTRUE(verification_passed)) "converged" else "approximate"
  if (controls$strict && !isTRUE(verification_passed)) {
    .dpprior_abort_invalid(
      sprintf(
        paste(
          "marginal PMF did not meet the higher-order distribution contract",
          "(L1 difference %.3g; tolerance %.3g; PMF/direct moments passed: %s)"
        ),
        l1_difference, tolerance,
        isTRUE(selected_moment_audit$passed) &&
          isTRUE(verification_moment_audit$passed)
      ),
      "dpprior_marginal_convergence_error", "M_verify",
      controls$M_verify, "selected and higher-order PMF agreement", reason
    )
  }

  metadata <- list(
    schema_version = 1L,
    engine = "gauss-laguerre-log-mixture",
    status = status,
    reason = reason,
    support = c(lower = 1L, upper = J),
    returned_support = c(lower = 0L, upper = J),
    normalization = list(
      log_constant = selected$log_normalization,
      probability_sum = selected$probability_sum,
      underflow_zero_count = selected$underflow_zero_count
    ),
    truncation = list(
      truncated = FALSE,
      requested_mass = 1,
      achieved_mass = 1,
      omitted_mass = 0
    ),
    verification = list(
      performed = !is.null(controls$M_verify),
      passed = verification_passed,
      M_selected = controls$M,
      M_verification_required = controls$M_verification_required,
      verification_available = controls$verification_available,
      M_verification = if (is.null(controls$M_verify)) {
        NA_integer_
      } else {
        controls$M_verify
      },
      l1_difference = l1_difference,
      max_difference = max_difference,
      tolerance = tolerance,
      l1_passed = l1_passed,
      moment_consistency = list(
        passed = if (is.null(verification_moment_audit)) {
          NA
        } else {
          isTRUE(selected_moment_audit$passed) &&
            isTRUE(verification_moment_audit$passed)
        },
        selected = selected_moment_audit,
        verification = verification_moment_audit
      ),
      absolute_tolerance = controls$abs_tol,
      relative_tolerance = controls$rel_tol
    ),
    selected = selected$node_rule,
    verification_rule = if (is.null(verification)) {
      NULL
    } else {
      verification$node_rule
    }
  )
  logp <- selected$logp
  attr(logp, "marginal_metadata") <- metadata
  if (!is.null(verification)) {
    attr(logp, ".marginal_verification_logp") <- verification$logp
  }
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
#' @param M_verify Optional independent quadrature order satisfying the
#'   package verification rule (at least \code{max(2*M, M+40)}, within the
#'   supported ceiling).
#' @param abs_tol,rel_tol Non-negative absolute and relative tolerances for the
#'   selected-versus-verification L1 discrepancy and PMF/direct-moment audits.
#' @param strict Logical; require successful higher-order verification when
#'   \code{TRUE}; otherwise return the selected PMF with explicit approximation
#'   metadata.
#'
#' @return Numeric vector of length \eqn{J+1} containing
#'   \eqn{P(K_J = k \mid a, b)} for \eqn{k = 0, 1, \ldots, J}.
#'   Entry \code{[1]} corresponds to \eqn{k=0} and always equals 0.
#'   The vector sums to 1 and carries a \code{"marginal_metadata"} attribute
#'   describing support, normalization, quadrature, verification status, and
#'   the absence of support truncation.
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
#' A \code{"converged"} result additionally requires higher-order PMF L1
#' agreement and PMF-versus-direct-moment agreement at both quadrature orders.
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
#' @family marginal_K
#'
#' @export
pmf_K_marginal <- function(
    J, a, b, logS, M = .QUAD_NODES_DEFAULT,
    M_verify = NULL, abs_tol = 1e-10, rel_tol = 1e-8,
    strict = FALSE) {
  # Compute in log-space for numerical stability
  logp <- log_pmf_K_marginal(
    J, a, b, logS, M, M_verify, abs_tol, rel_tol, strict
  )
  metadata <- attr(logp, "marginal_metadata", exact = TRUE)
  verification_logp <- attr(
    logp, ".marginal_verification_logp", exact = TRUE
  )

  # Convert to linear scale
  pmf <- exp(logp)

  # Enforce P(K=0) = 0 exactly. Any invalid mixed mass is a typed failure;
  # normalization is never replaced by a synthetic distribution.
  pmf[1L] <- 0.0
  pre_rescale_sum <- sum(pmf)
  if (!is.finite(pre_rescale_sum) || pre_rescale_sum <= 0) {
    .marginal_pmf_abort(
      "marginal PMF has invalid linear-scale probability mass",
      "probability_sum", pre_rescale_sum,
      "finite positive probability mass", "invalid_probability_mass"
    )
  }
  pmf <- pmf / pre_rescale_sum
  .dpprior_validate_pmf(pmf, expected_length = as.integer(J) + 1L)

  metadata$normalization$pre_rescale_probability_sum <- pre_rescale_sum
  metadata$normalization$probability_sum <- sum(pmf)
  metadata$normalization$linear_rescale <- abs(pre_rescale_sum - 1)
  metadata$normalization$underflow_zero_count <- sum(pmf[-1L] == 0)
  metadata$truncation$achieved_mass <- sum(pmf[-1L])
  metadata$truncation$omitted_mass <- 0
  attr(pmf, "marginal_metadata") <- metadata

  # Retain the independent distribution only as private computational state
  # for CDF/quantile/mode/summary verification. It is never returned in place
  # of the selected-order PMF.
  if (!is.null(verification_logp)) {
    verification_pmf <- exp(verification_logp)
    verification_pmf[[1L]] <- 0
    verification_mass <- sum(verification_pmf)
    if (!is.finite(verification_mass) || verification_mass <= 0) {
      .marginal_pmf_abort(
        "verification PMF has invalid linear-scale probability mass",
        "verification_probability_sum", verification_mass,
        "finite positive probability mass", "invalid_probability_mass"
      )
    }
    verification_pmf <- verification_pmf / verification_mass
    .dpprior_validate_pmf(
      verification_pmf, expected_length = as.integer(J) + 1L
    )
    attr(pmf, ".marginal_verification_pmf") <- verification_pmf
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
#' @param M_verify Optional independent order satisfying the package
#'   verification rule (at least \code{max(2*M, M+40)}).
#' @param abs_tol,rel_tol Non-negative selected-versus-verification tolerances.
#' @param strict Logical; require successful higher-order verification.
#'
#' @return Numeric vector of length \eqn{J+1} containing
#'   \eqn{P(K_J \leq k \mid a, b)} for \eqn{k = 0, 1, \ldots, J}.
#'   The vector carries the PMF's \code{"marginal_metadata"} attribute.
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
#' @family marginal_K
#'
#' @export
cdf_K_marginal <- function(
    J, a, b, logS, M = .QUAD_NODES_DEFAULT,
    M_verify = NULL, abs_tol = 1e-10, rel_tol = 1e-8,
    strict = FALSE) {
  pmf <- pmf_K_marginal(
    J, a, b, logS, M, M_verify, abs_tol, rel_tol, strict
  )
  metadata <- attr(pmf, "marginal_metadata", exact = TRUE)
  verification_pmf <- attr(
    pmf, ".marginal_verification_pmf", exact = TRUE
  )
  cdf <- .marginal_cdf_from_pmf(pmf)
  attr(cdf, "marginal_metadata") <- metadata
  if (!is.null(verification_pmf)) {
    attr(cdf, ".marginal_verification_cdf") <-
      .marginal_cdf_from_pmf(verification_pmf)
  }
  cdf
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
#' @param M_verify Optional independent order satisfying the package
#'   verification rule (at least \code{max(2*M, M+40)}).
#' @param abs_tol,rel_tol Non-negative selected-versus-verification tolerances.
#' @param strict Logical; require successful higher-order verification.
#'
#' @return Integer vector of quantiles (same length as \code{p}).
#'   Each element is the smallest \eqn{k} such that \eqn{P(K_J \leq k) \geq p}.
#'   The vector carries the CDF's \code{"marginal_metadata"} attribute so an
#'   approximate, unverified, or discrepant quadrature result is not silent.
#'
#' @details
#' This is the standard quantile definition for discrete distributions:
#' \eqn{Q(p) = \min\{k : F(k) \geq p\}}.
#' The lower endpoint follows the mathematical support: \eqn{Q(0)=1}, while
#' \eqn{Q(1)=J}. The structural \eqn{k=0} vector entry is never returned.
#'
#' The function is vectorized over \code{p}, allowing efficient computation
#' of multiple quantiles in a single call.
#' When \code{M_verify} is supplied, the selected- and verification-order
#' integer quantiles are also compared exactly. A changed quantile downgrades
#' metadata to \code{"approximate"}; \code{strict = TRUE} raises a typed
#' convergence error. The selected-order quantile is always retained.
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
#' @family marginal_K
#'
#' @export
quantile_K_marginal <- function(
    p, J, a, b, logS, M = .QUAD_NODES_DEFAULT,
    M_verify = NULL, abs_tol = 1e-10, rel_tol = 1e-8,
    strict = FALSE) {
  # Input validation
  assert_probability(p, "p")
  assert_valid_J(J)

  J <- as.integer(J)
  cdf <- cdf_K_marginal(
    J, a, b, logS, M, M_verify, abs_tol, rel_tol, strict
  )

  quantiles <- .marginal_quantile_from_cdf(p, cdf, J)
  verification_cdf <- attr(
    cdf, ".marginal_verification_cdf", exact = TRUE
  )
  verification_quantiles <- if (is.null(verification_cdf)) {
    NULL
  } else {
    .marginal_quantile_from_cdf(p, verification_cdf, J)
  }
  metadata <- .marginal_record_discrete_verification(
    attr(cdf, "marginal_metadata", exact = TRUE),
    "quantile", quantiles, verification_quantiles,
    probabilities = p
  )
  .marginal_enforce_discrete_verification(
    metadata, strict, "marginal quantile"
  )
  attr(quantiles, "marginal_metadata") <- metadata
  quantiles
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
#' @param M_verify Optional independent order satisfying the package
#'   verification rule (at least \code{max(2*M, M+40)}).
#' @param abs_tol,rel_tol Non-negative selected-versus-verification tolerances.
#' @param strict Logical; require successful higher-order verification.
#'
#' @return Integer; the value \eqn{k} that maximizes \eqn{P(K_J = k \mid a, b)},
#'   carrying the PMF's \code{"marginal_metadata"} attribute.
#'
#' @details
#' The mode is always >= 1 since \eqn{P(K_J = 0) = 0}.
#' When \code{M_verify} is supplied, the selected and verification modes must
#' be identical for converged mode metadata. Disagreement is explicit and
#' never causes substitution of the verification-order mode.
#'
#' @examples
#' logS <- compute_log_stirling(50)
#' mode_K_marginal(50, 1.5, 0.5, logS)
#'
#' @seealso \code{\link{pmf_K_marginal}}, \code{\link{summary_K_marginal}}
#'
#' @family marginal_K
#'
#' @export
mode_K_marginal <- function(
    J, a, b, logS, M = .QUAD_NODES_DEFAULT,
    M_verify = NULL, abs_tol = 1e-10, rel_tol = 1e-8,
    strict = FALSE) {
  pmf <- pmf_K_marginal(
    J, a, b, logS, M, M_verify, abs_tol, rel_tol, strict
  )
  mode <- as.integer(which.max(pmf) - 1L)
  verification_pmf <- attr(
    pmf, ".marginal_verification_pmf", exact = TRUE
  )
  verification_mode <- if (is.null(verification_pmf)) {
    NULL
  } else {
    as.integer(which.max(verification_pmf) - 1L)
  }
  metadata <- .marginal_record_discrete_verification(
    attr(pmf, "marginal_metadata", exact = TRUE),
    "mode", mode, verification_mode
  )
  .marginal_enforce_discrete_verification(metadata, strict, "marginal mode")
  attr(mode, "marginal_metadata") <- metadata
  mode
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
#' @param M_verify Optional independent order satisfying the package
#'   verification rule (at least \code{max(2*M, M+40)}).
#' @param abs_tol,rel_tol Non-negative selected-versus-verification tolerances.
#' @param strict Logical; require successful higher-order verification.
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
#'     \item{\code{metadata}}{Support, normalization, no-truncation,
#'       quadrature-status, and achieved-quantile-mass metadata.}
#'   }
#'
#' @details
#' This function provides a complete summary of the marginal distribution,
#' combining PMF-based and CDF-based statistics. The mean and variance
#' computed from the PMF should match \code{exact_K_moments()} within
#' numerical tolerance.
#' With \code{M_verify}, reported quantiles, the median, and the mode are
#' compared exactly across quadrature orders. Any changed discrete summary
#' downgrades the returned metadata to \code{"approximate"}, or raises a typed
#' convergence error in strict mode, while preserving all selected-order
#' values.
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
#' @family marginal_K
#'
#' @export
summary_K_marginal <- function(
    J, a, b, logS, M = .QUAD_NODES_DEFAULT,
    probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
    M_verify = NULL, abs_tol = 1e-10, rel_tol = 1e-8,
    strict = FALSE) {
  # Input validation
  assert_probability(probs, "probs")
  assert_valid_J(J)
  J <- as.integer(J)

  # Compute PMF and CDF
  pmf <- pmf_K_marginal(
    J, a, b, logS, M, M_verify, abs_tol, rel_tol, strict
  )
  metadata <- attr(pmf, "marginal_metadata", exact = TRUE)
  verification_pmf <- attr(
    pmf, ".marginal_verification_pmf", exact = TRUE
  )
  cdf <- .marginal_cdf_from_pmf(pmf)
  verification_cdf <- if (is.null(verification_pmf)) {
    NULL
  } else {
    .marginal_cdf_from_pmf(verification_pmf)
  }

  k_vals <- 0:J

  # Mean and variance from PMF
  mean_K <- .quadrature_weighted_sum(pmf, k_vals)
  var_K <- .quadrature_weighted_sum(pmf, (k_vals - mean_K)^2)
  if (!is.finite(mean_K) || mean_K < 1 || mean_K > J ||
      !is.finite(var_K) || var_K < 0) {
    .dpprior_abort_invalid(
      "marginal PMF summary violated support or variance requirements",
      c("dpprior_marginal_summary_error", "dpprior_numerical_error"),
      "summary", c(mean = mean_K, var = var_K),
      sprintf("finite mean in [1,%d] and non-negative variance", J),
      "moment_contract"
    )
  }
  sd_K <- sqrt(var_K)
  cv_K <- if (mean_K > 0) sd_K / mean_K else Inf

  # Mode
  mode_K <- as.integer(which.max(pmf) - 1L)

  # Quantiles (vectorized)
  qs <- .marginal_quantile_from_cdf(probs, cdf, J)
  names(qs) <- paste0("q", formatC(100 * probs, format = "f", digits = 0))

  # Median (always include even if not in probs)
  median_K <- .marginal_quantile_from_cdf(0.5, cdf, J)[[1L]]

  verification_mode <- if (is.null(verification_pmf)) {
    NULL
  } else {
    as.integer(which.max(verification_pmf) - 1L)
  }
  verification_qs <- if (is.null(verification_cdf)) {
    NULL
  } else {
    .marginal_quantile_from_cdf(probs, verification_cdf, J)
  }
  verification_median <- if (is.null(verification_cdf)) {
    NULL
  } else {
    .marginal_quantile_from_cdf(0.5, verification_cdf, J)[[1L]]
  }
  metadata <- .marginal_record_discrete_verification(
    metadata, "mode", mode_K, verification_mode
  )
  metadata <- .marginal_record_discrete_verification(
    metadata, "quantile", qs, verification_qs, probabilities = probs
  )
  metadata <- .marginal_record_discrete_verification(
    metadata, "median", median_K, verification_median,
    probabilities = 0.5
  )
  .marginal_enforce_discrete_verification(
    metadata, strict, "marginal summary"
  )

  metadata$quantiles <- list(
    probabilities = probs,
    values = unname(qs),
    achieved_cdf = unname(cdf[qs + 1L]),
    previous_cdf = vapply(qs, function(q) {
      if (q <= 1L) 0 else unname(cdf[[q]])
    }, numeric(1))
  )
  attr(pmf, "marginal_metadata") <- metadata
  attr(cdf, "marginal_metadata") <- metadata
  if (!is.null(verification_cdf)) {
    attr(cdf, ".marginal_verification_cdf") <- verification_cdf
  }

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
    cdf = cdf,
    metadata = metadata
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
#' \dontrun{
#' logS <- compute_log_stirling(50)
#' mean_K_from_marginal_pmf(50, 1.5, 0.5, logS)
#'
#' }
#' @keywords internal
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
#' \dontrun{
#' logS <- compute_log_stirling(50)
#' var_K_from_marginal_pmf(50, 1.5, 0.5, logS)
#'
#' }
#' @keywords internal
var_K_from_marginal_pmf <- function(J, a, b, logS, M = .QUAD_NODES_DEFAULT) {
  pmf <- pmf_K_marginal(J, a, b, logS, M)
  k_vals <- 0:J
  mean_K <- sum(k_vals * pmf)
  sum((k_vals - mean_K)^2 * pmf)
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
#' \dontrun{
#' logS <- compute_log_stirling(50)
#' verify_pmf_marginal_properties(50, 1.5, 0.5, logS)
#'
#' }
#' @keywords internal
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
#' \dontrun{
#' logS <- compute_log_stirling(50)
#' verify_pmf_marginal_moments(50, 1.5, 0.5, logS)
#'
#' }
#' @keywords internal
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


#' Report Marginal PMF Quadrature Discrepancies
#'
#' Reports successive-order PMF discrepancies and classifies each selected/
#' refined pair against an explicit L1 budget. Successive discrepancies are
#' not assumed to decrease monotonically.
#'
#' @param J Integer; sample size.
#' @param a Numeric; shape parameter of Gamma prior.
#' @param b Numeric; rate parameter of Gamma prior.
#' @param logS Matrix; pre-computed log-Stirling matrix.
#' @param M_values Integer vector; quadrature node counts to test.
#' @param verbose Logical; if TRUE, print results.
#' @param abs_tol,rel_tol Non-negative L1 error-budget components.
#'
#' @return Data frame of successive-order discrepancies, L1 budgets, and
#'   pairwise budget classifications.
#'
#' @examples
#' \dontrun{
#' logS <- compute_log_stirling(50)
#' verify_pmf_marginal_convergence(50, 1.5, 0.5, logS)
#'
#' }
#' @keywords internal
verify_pmf_marginal_convergence <- function(J, a, b, logS,
                                            M_values = c(20L, 40L, 80L, 120L),
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
    M = integer(length(M_values)),
    mean = numeric(length(M_values)),
    var = numeric(length(M_values)),
    mean_change = rep(NA_real_, length(M_values)),
    var_change = rep(NA_real_, length(M_values)),
    L1_change = rep(NA_real_, length(M_values)),
    L1_tolerance = rep(NA_real_, length(M_values)),
    within_budget = rep(NA, length(M_values))
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
      results$L1_tolerance[i] <- abs_tol + rel_tol
      results$within_budget[i] <-
        results$L1_change[i] <= results$L1_tolerance[i]
    }

    prev_pmf <- pmf
  }

  if (isTRUE(verbose)) {
    cat(sprintf("Quadrature Order Discrepancies (J=%d, a=%.2f, b=%.2f):\n",
                J, a, b))
    print(results, row.names = FALSE)
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
#' \dontrun{
#' verify_pmf_marginal_all()
#'
#' }
#' @keywords internal
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
