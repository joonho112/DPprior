# =============================================================================
# Module 12: A2-KL KL Divergence Minimization for Prior Calibration
# =============================================================================
#
# This module implements the A2-KL algorithm for calibrating Gamma hyperpriors
# on the Dirichlet process concentration parameter by minimizing KL divergence
# between a target distribution and the induced marginal PMF of K_J.
#
# Theory Background (RN-01, RN-04):
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
# - Fallback mechanism if optimization worsens the objective
# - Comprehensive trace recording for diagnostics
#
# Author: JoonHo Lee
# Date: December 2025
# Part of: DPprior R Package
# Reference: RN-01 Section 6.5, RN-04 Section 4
# Dependencies: Modules 00 (constants), 02 (quadrature), 04 (conditional PMF),
#               05 (marginal moments), 06 (marginal PMF), 10 (A1), 11 (A2-MN)
# =============================================================================


# =============================================================================
# Internal Helper Functions
# =============================================================================

#' Normalize and Validate a Target PMF on \{1, ..., J\}
#'
#' @param target_pmf Numeric vector; unnormalized target PMF.
#' @param J Integer; sample size.
#'
#' @return Numeric vector; normalized PMF of length J (support k=1,...,J).
#'
#' @keywords internal
.a2_kl_normalize_pmf <- function(target_pmf, J) {
  if (!is.numeric(target_pmf)) {
    stop("target_pmf must be a numeric vector.", call. = FALSE)
  }

  # Handle both length J and J+1 inputs

  if (length(target_pmf) == J + 1L) {
    # Drop k=0 entry
    target_pmf <- target_pmf[-1L]
  } else if (length(target_pmf) != J) {
    stop(sprintf("target_pmf must have length J=%d or J+1=%d.", J, J + 1L),
         call. = FALSE)
  }

  if (anyNA(target_pmf) || any(!is.finite(target_pmf))) {
    stop("target_pmf must be finite and non-missing.", call. = FALSE)
  }
  if (any(target_pmf < 0)) {
    stop("target_pmf must be non-negative.", call. = FALSE)
  }

  s <- sum(target_pmf)
  if (!is.finite(s) || s <= 0) {
    stop("target_pmf must have positive total mass.", call. = FALSE)
  }

  target_pmf / s
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
  pmf_full <- pmf_K_marginal(J = J, a = a, b = b, logS = logS, M = M)

  # Drop k=0 entry and normalize
  induced <- pmf_full[-1L]
  induced <- pmax(induced, 0)

  s <- sum(induced)
  if (!is.finite(s) || s <= 0) {
    stop("Induced PMF computation failed (non-finite or zero total mass).",
         call. = FALSE)
  }

  induced / s
}


#' Compute KL Divergence from Target to Induced PMF
#'
#' @param target_pmf Numeric vector; normalized PMF (length J).
#' @param induced_pmf Numeric vector; normalized induced PMF (length J).
#' @param eps Numeric; small constant to avoid log(0).
#'
#' @return Numeric; KL divergence (non-negative).
#'
#' @keywords internal
.a2_kl_compute_kl <- function(target_pmf, induced_pmf, eps = 1e-15) {
  idx <- target_pmf > 0

  if (!any(idx)) {
    return(0)
  }

  kl <- sum(target_pmf[idx] * (log(target_pmf[idx]) - log(induced_pmf[idx] + eps)))
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
#' @param eps Numeric; small value to prevent \code{log(0)}. Default: 1e-15.
#'
#' @return Numeric scalar; the KL divergence (non-negative).
#'
#' @details
#' The KL divergence is defined as:
#' \deqn{D_{KL}(p \| q) = \sum_k p(k) \log\frac{p(k)}{q(k)}}
#'
#' Only indices where \eqn{p(k) > \epsilon} are included in the sum.
#'
#' \strong{Properties:}
#' \itemize{
#'   \item \eqn{D_{KL}(p \| q) \geq 0} with equality iff \eqn{p = q}
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
kl_divergence_pmf <- function(p, q, eps = 1e-15) {
  if (!is.numeric(p) || !is.numeric(q)) {
    stop("p and q must be numeric vectors", call. = FALSE)
  }
  if (length(p) != length(q)) {
    stop("p and q must have the same length", call. = FALSE)
  }

  mask <- p > eps

  if (!any(mask)) {
    return(0)
  }

  p_safe <- p[mask]
  q_safe <- pmax(q[mask], eps)

  kl <- sum(p_safe * log(p_safe / q_safe))
  max(0, kl)
}


#' KL Divergence Between Target and Induced K_J PMFs
#'
#' Computes the KL divergence \eqn{D_{KL}(p^* \| p_{a,b})} between a target PMF
#' and the induced marginal PMF of \eqn{K_J} under \eqn{\alpha \sim Gamma(a, b)}.
#'
#' @param target_pmf Numeric vector; target PMF for K_J. Can have length J
#'   (support k=1,...,J) or J+1 (support k=0,...,J, where k=0 is ignored).
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

  if (!is.numeric(M) || length(M) != 1L || M < 10L) {
    stop("M must be a positive integer >= 10", call. = FALSE)
  }

  # Normalize target PMF (handles both length J and J+1)
  target_pmf <- .a2_kl_normalize_pmf(target_pmf, J)

  # Compute induced PMF
  logS <- compute_log_stirling(J)
  induced <- .a2_kl_induced_pmf(J, a, b, logS, M)

  # Compute KL divergence
  .a2_kl_compute_kl(target_pmf, induced, eps = 1e-15)
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
  assert_positive(df, "df")
  assert_positive(scale, "scale")

  k <- seq_len(J)
  upper <- (k + 0.5) / scale
  lower <- (k - 0.5) / scale

  pmf <- stats::pchisq(upper, df = df) - stats::pchisq(lower, df = df)
  pmf <- pmax(pmf, 0)

  s <- sum(pmf)
  if (!is.finite(s) || s <= 0) {
    stop("Chi-square discretization produced zero mass; check df/scale/J.",
         call. = FALSE)
  }

  pmf / s
}


#' Construct Target PMF from User Specification
#'
#' Creates a normalized target PMF from either a user-provided PMF vector
#' or target moment specification.
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
#'     \item{\code{pmf}}{Numeric vector of length J; normalized target PMF}
#'     \item{\code{mu_K}}{Target mean}
#'     \item{\code{var_K}}{Target variance}
#'     \item{\code{df}}{(if moments provided) Chi-square degrees of freedom}
#'     \item{\code{scale}}{(if moments provided) Chi-square scale parameter}
#'   }
#'
#' @details
#' \strong{Direct PMF specification:}
#' If \code{target} is a numeric vector of length J, it is treated as the PMF
#' for k = 1, ..., J. If length J+1, the k=0 entry is dropped.
#'
#' \strong{Moment specification:}
#' If \code{target} is a list with \code{mu_K} and \code{var_K}, a discretized
#' chi-square distribution matching these moments is constructed using:
#' \deqn{\text{scale} = \sigma^2_K / (2\mu_K), \quad \text{df} = 2\mu_K^2 / \sigma^2_K}
#'
#' @examples
#' # Direct PMF
#' result <- construct_target_pmf(50, rep(1, 50))  # Uniform on 1:50
#'
#' # Moment specification
#' result <- construct_target_pmf(50, list(mu_K = 5, var_K = 8))
#' result$mu_K  # 5
#' result$df    # 6.25
#'
#' @seealso \code{\link{discretize_chisq}}, \code{\link{DPprior_a2_kl}}
#'
#' @export
construct_target_pmf <- function(J, target) {
  J <- as.integer(J)
  assert_valid_J(J)

  if (is.numeric(target)) {
    # Direct PMF specification
    target_pmf <- .a2_kl_normalize_pmf(target, J)

    # Compute moments from PMF
    k_vals <- seq_len(J)
    mu_K <- sum(k_vals * target_pmf)
    var_K <- sum(k_vals^2 * target_pmf) - mu_K^2

    list(
      pmf = target_pmf,
      mu_K = mu_K,
      var_K = var_K
    )

  } else if (is.list(target) && all(c("mu_K", "var_K") %in% names(target))) {
    # Moment specification
    mu_K <- target$mu_K
    var_K <- target$var_K

    if (!is.numeric(mu_K) || length(mu_K) != 1L || !is.finite(mu_K) || mu_K <= 0) {
      stop("mu_K must be a positive finite numeric scalar", call. = FALSE)
    }
    if (!is.numeric(var_K) || length(var_K) != 1L || !is.finite(var_K) || var_K <= 0) {
      stop("var_K must be a positive finite numeric scalar", call. = FALSE)
    }
    if (mu_K <= 1) {
      stop("mu_K must be > 1 (at least one cluster is always present)", call. = FALSE)
    }
    if (mu_K > J) {
      stop("mu_K must be <= J", call. = FALSE)
    }

    # Construct discretized chi-square
    # For Y = scale * X where X ~ chi2(df):
    #   E[Y] = scale * df = mu_K
    #   Var[Y] = scale^2 * 2 * df = var_K
    # Solving: scale = var_K / (2 * mu_K), df = 2 * mu_K^2 / var_K
    scale <- var_K / (2 * mu_K)
    df <- 2 * mu_K^2 / var_K

    target_pmf <- discretize_chisq(J, df, scale)

    # Compute discretized moments (may differ due to truncation)
    k_vals <- seq_len(J)
    mu_disc <- sum(k_vals * target_pmf)
    var_disc <- sum(k_vals^2 * target_pmf) - mu_disc^2

    list(
      pmf = target_pmf,
      mu_K = mu_K,
      var_K = var_K,
      df = df,
      scale = scale,
      mu_K_discrete = mu_disc,
      var_K_discrete = var_disc
    )

  } else {
    stop("target must be a PMF vector or list with 'mu_K' and 'var_K'",
         call. = FALSE)
  }
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
#'     \item Numeric vector of length J: target PMF for k = 1, ..., J
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
#'     \item \code{eps}: numeric; small constant to avoid \code{log(0)}
#'           (default \code{1e-15}).
#'   }
#'
#' @return A \code{DPprior_fit} object (S3 class) with components:
#'   \describe{
#'     \item{\code{a}}{Numeric; optimal shape parameter}
#'     \item{\code{b}}{Numeric; optimal rate parameter}
#'     \item{\code{J}}{Integer; sample size}
#'     \item{\code{target}}{List with target specification (pmf, mu_K, var_K, type,
#'       and for chisq: df, scale, mu_K_discrete, var_K_discrete)}
#'     \item{\code{method}}{Character; "A2-KL"}
#'     \item{\code{converged}}{Logical; whether optimization converged}
#'     \item{\code{iterations}}{Integer; number of function evaluations}
#'     \item{\code{termination}}{Character; termination reason}
#'     \item{\code{fit}}{List with achieved \code{mu_K}, \code{var_K},
#'       \code{kl}, and \code{residual}}
#'     \item{\code{diagnostics}}{List with initialization info, optim details,
#'       fallback status, and KL at init/final}
#'     \item{\code{trace}}{Data frame of evaluation-level trace (eval, a, b, kl)}
#'     \item{\code{status}}{Character; convergence status}
#'   }
#'
#' @details
#' The A2-KL algorithm finds \eqn{(a^*, b^*)} by solving:
#' \deqn{(a^*, b^*) = \arg\min_{a,b} D_{KL}(p^*(K_J) \| p_{a,b}(K_J))}
#'
#' \strong{Algorithm:}
#' \enumerate{
#'   \item Construct target PMF from user specification
#'   \item Initialize \eqn{(a_0, b_0)} from A2-MN (exact moment matching)
#'   \item Optimize KL divergence using L-BFGS-B in log-space with bounds
#'   \item Apply fallback to initialization if optimization worsens KL
#'   \item Return optimal parameters and comprehensive diagnostics
#' }
#'
#' \strong{Stability features:}
#' \itemize{
#'   \item Log-parameterization ensures positivity of (a, b)
#'   \item Explicit bounds prevent numerical overflow/underflow
#'   \item Fallback mechanism if optimization diverges or worsens objective
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
#' fit2 <- DPprior_a2_kl(J = 50, target = target_pmf, method = "pmf")
#'
#' # Example 3: Compare A2-KL vs A2-MN
#' a2_mn <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8)
#' a2_kl <- DPprior_a2_kl(J = 50, target = list(mu_K = 5, var_K = 8),
#'                        method = "chisq")
#' cat(sprintf("A2-MN: a=%.4f, b=%.4f\n", a2_mn$a, a2_mn$b))
#' cat(sprintf("A2-KL: a=%.4f, b=%.4f, KL=%.4e\n", a2_kl$a, a2_kl$b, a2_kl$fit$kl))
#'
#' @export
DPprior_a2_kl <- function(J, target,
                          method = c("pmf", "chisq"),
                          max_iter = 100L,
                          tol = 1e-6,
                          M = .QUAD_NODES_DEFAULT,
                          verbose = FALSE,
                          ...) {

  # -------------------------------------------------------------------------
  # Input Validation
  # -------------------------------------------------------------------------
  assert_valid_J(J)
  J <- as.integer(J)

  if (J < 2L) {
    stop("A2-KL requires J >= 2.", call. = FALSE)
  }

  method <- match.arg(method)

  if (!is.numeric(max_iter) || length(max_iter) != 1L ||
      !is.finite(max_iter) || max_iter != floor(max_iter) || max_iter < 1L) {
    stop("max_iter must be a positive integer.", call. = FALSE)
  }
  max_iter <- as.integer(max_iter)

  assert_positive(tol, "tol")

  if (!is.numeric(M) || length(M) != 1L || M < 10L) {
    stop("M must be a positive integer >= 10", call. = FALSE)
  }
  M <- as.integer(M)

  # Optional parameters from ...
  dots <- list(...)
  log_bounds <- if (!is.null(dots$log_bounds)) dots$log_bounds else c(-15, 15)
  if (!is.numeric(log_bounds) || length(log_bounds) != 2L ||
      any(!is.finite(log_bounds)) || log_bounds[1] >= log_bounds[2]) {
    stop("log_bounds must be a finite numeric vector of length 2 with lower < upper.",
         call. = FALSE)
  }
  eps <- if (!is.null(dots$eps)) dots$eps else 1e-15
  assert_positive(eps, "eps")

  # -------------------------------------------------------------------------
  # Construct Target PMF
  # -------------------------------------------------------------------------
  target_info <- list(type = method)

  if (identical(method, "pmf")) {
    # Direct PMF specification
    target_pmf <- .a2_kl_normalize_pmf(target, J)

    k <- seq_len(J)
    mu_target <- sum(k * target_pmf)
    var_target <- sum((k^2) * target_pmf) - mu_target^2

    target_info$pmf <- target_pmf
    target_info$mu_K <- mu_target
    target_info$var_K <- var_target

  } else {
    # method == "chisq": Moment specification with chi-square discretization
    if (!is.list(target) || !all(c("mu_K", "var_K") %in% names(target))) {
      stop("For method='chisq', target must be a list with 'mu_K' and 'var_K'.",
           call. = FALSE)
    }

    mu_target <- target$mu_K
    var_target <- target$var_K
    assert_positive(mu_target, "target$mu_K")
    assert_positive(var_target, "target$var_K")

    if (mu_target <= 1) {
      stop("mu_K must be > 1 (at least one cluster).", call. = FALSE)
    }
    if (mu_target >= J) {
      stop("mu_K must be < J.", call. = FALSE)
    }

    # Scaled chi-square moment matching:
    # If X = scale * chisq(df), then E[X] = df*scale, Var[X] = 2*df*scale^2
    # Solving: df = 2*mu^2/var, scale = var/(2*mu)
    df <- 2 * mu_target^2 / var_target
    scale <- var_target / (2 * mu_target)

    target_pmf <- discretize_chisq(J, df = df, scale = scale)

    # Discretization changes moments due to truncation/renormalization
    k <- seq_len(J)
    mu_disc <- sum(k * target_pmf)
    var_disc <- sum((k^2) * target_pmf) - mu_disc^2

    target_info$pmf <- target_pmf
    target_info$mu_K <- mu_target
    target_info$var_K <- var_target
    target_info$df <- df
    target_info$scale <- scale
    target_info$mu_K_discrete <- mu_disc
    target_info$var_K_discrete <- var_disc
  }

  if (isTRUE(verbose)) {
    cat("A2-KL: KL Divergence Minimization\n")
    cat(sprintf("  J = %d, method = '%s'\n", J, method))
    cat(sprintf("  Target: mu_K = %.4f, var_K = %.4f\n", mu_target, var_target))
    if (identical(method, "chisq")) {
      cat(sprintf("  Discretized: mu = %.4f, var = %.4f (df=%.2f, scale=%.4f)\n",
                  target_info$mu_K_discrete, target_info$var_K_discrete,
                  target_info$df, target_info$scale))
    }
  }

  # -------------------------------------------------------------------------
  # Pre-compute Stirling Numbers
  # -------------------------------------------------------------------------
  logS <- compute_log_stirling(J)

  # -------------------------------------------------------------------------
  # Initialization from A2-MN (Exact Moment Matching)
  # -------------------------------------------------------------------------
  # Use discretized moments if available for more accurate initialization
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

  # Clamp to valid range
  if (mu_init <= 1) mu_init <- 1 + 1e-3
  if (mu_init >= J) mu_init <- J - 1e-3
  var_init <- max(var_init, 1e-8)

  init_fit <- tryCatch(
    DPprior_a2_newton(J = J, mu_K = mu_init, var_K = var_init,
                      max_iter = 15L, tol_F = 1e-8, verbose = FALSE, M = M),
    error = function(e) NULL
  )

  if (is.null(init_fit)) {
    # Fallback to A1 closed-form
    init_fit <- tryCatch(
      DPprior_a1(J, mu_init, var_init),
      error = function(e) NULL
    )
  }

  if (is.null(init_fit)) {
    # Ultimate fallback: simple heuristic
    c_J <- log(J)
    m <- max(0.5, mu_init - 1)
    vif <- max(1.01, var_init / mu_init)
    a0 <- max(0.1, min(100, m / (vif - 1)))
    b0 <- max(0.01, min(100, a0 * c_J / m))
    init_method <- "heuristic"
    init_status <- "heuristic_fallback"
    init_residual <- NA_real_
  } else {
    a0 <- init_fit$a
    b0 <- init_fit$b
    init_method <- init_fit$method
    init_status <- if (!is.null(init_fit$status)) init_fit$status else "success"
    init_residual <- if (!is.null(init_fit$fit$residual)) init_fit$fit$residual else NA_real_
  }

  if (isTRUE(verbose)) {
    cat(sprintf("  Initialization (%s): a0 = %.4f, b0 = %.4f\n",
                init_method, a0, b0))
  }

  # -------------------------------------------------------------------------
  # Set up Trace Recording Environment
  # -------------------------------------------------------------------------
  trace_env <- new.env(parent = emptyenv())
  trace_env$n_eval <- 0L
  trace_env$records <- list()

  # -------------------------------------------------------------------------
  # Objective Function (KL divergence in log-space)
  # -------------------------------------------------------------------------
  objective <- function(eta) {
    a_curr <- exp(eta[1L])
    b_curr <- exp(eta[2L])

    kl <- tryCatch({
      induced <- .a2_kl_induced_pmf(J, a_curr, b_curr, logS, M)
      .a2_kl_compute_kl(target_pmf, induced, eps = eps)
    }, error = function(e) 1e10)

    if (!is.finite(kl)) kl <- 1e10

    # Record trace
    trace_env$n_eval <- trace_env$n_eval + 1L
    trace_env$records[[trace_env$n_eval]] <- c(
      eval = trace_env$n_eval,
      a = a_curr,
      b = b_curr,
      kl = kl
    )

    if (isTRUE(verbose) && (trace_env$n_eval %% 10L == 0L)) {
      cat(sprintf("  eval=%d | a=%.6g | b=%.6g | KL=%.6g\n",
                  trace_env$n_eval, a_curr, b_curr, kl))
    }

    kl
  }

  # -------------------------------------------------------------------------
  # Optimization (L-BFGS-B in log-space with bounds)
  # -------------------------------------------------------------------------
  eta0 <- log(c(a0, b0))
  kl0 <- objective(eta0)

  control <- list(
    maxit = max_iter,
    factr = tol / .Machine$double.eps,
    pgtol = tol
  )

  opt <- stats::optim(
    par = eta0,
    fn = objective,
    method = "L-BFGS-B",
    lower = rep(log_bounds[1L], 2L),
    upper = rep(log_bounds[2L], 2L),
    control = control
  )

  a_opt <- exp(opt$par[1L])
  b_opt <- exp(opt$par[2L])

  # Compute final KL
  induced_opt <- tryCatch(
    .a2_kl_induced_pmf(J, a_opt, b_opt, logS, M),
    error = function(e) NULL
  )

  if (is.null(induced_opt)) {
    kl_opt <- Inf
  } else {
    kl_opt <- .a2_kl_compute_kl(target_pmf, induced_opt, eps = eps)
  }

  converged <- identical(opt$convergence, 0L)

  # -------------------------------------------------------------------------
  # Fallback: If optimization worsened the objective, revert to init
  # -------------------------------------------------------------------------
  fallback_used <- FALSE
  if (!is.finite(kl_opt) || kl_opt > kl0 + 1e-12) {
    fallback_used <- TRUE
    a_opt <- a0
    b_opt <- b0
    induced_opt <- .a2_kl_induced_pmf(J, a_opt, b_opt, logS, M)
    kl_opt <- .a2_kl_compute_kl(target_pmf, induced_opt, eps = eps)
    converged <- FALSE
  }

  # -------------------------------------------------------------------------
  # Compute Final Fit Statistics
  # -------------------------------------------------------------------------
  moments <- exact_K_moments(J, a_opt, b_opt, M = M)

  if (isTRUE(verbose)) {
    cat(sprintf("\nResults:\n"))
    cat(sprintf("  a* = %.6f, b* = %.6f\n", a_opt, b_opt))
    cat(sprintf("  Achieved: mu_K = %.4f, var_K = %.4f\n",
                moments$mean, moments$var))
    cat(sprintf("  KL divergence = %.4e\n", kl_opt))
    if (fallback_used) {
      cat("  WARNING: Optimization worsened KL; fell back to initialization.\n")
    } else {
      cat(sprintf("  Converged: %s\n", if (converged) "TRUE" else "FALSE"))
    }
  }

  # -------------------------------------------------------------------------
  # Build Trace Data Frame
  # -------------------------------------------------------------------------
  trace_df <- NULL
  if (length(trace_env$records) > 0L) {
    trace_mat <- do.call(rbind, trace_env$records)
    trace_df <- data.frame(
      eval = as.integer(trace_mat[, "eval"]),
      a = as.numeric(trace_mat[, "a"]),
      b = as.numeric(trace_mat[, "b"]),
      kl = as.numeric(trace_mat[, "kl"])
    )
  }

  # -------------------------------------------------------------------------
  # Determine Status and Termination
  # -------------------------------------------------------------------------
  status <- if (isTRUE(fallback_used)) {
    "fallback_used"
  } else if (isTRUE(converged)) {
    "success"
  } else {
    "not_converged"
  }

  termination <- if (isTRUE(fallback_used)) {
    "fallback_to_init"
  } else if (identical(opt$convergence, 0L)) {
    "optim_converged"
  } else {
    paste0("optim_convergence_", opt$convergence)
  }

  # -------------------------------------------------------------------------
  # Build Diagnostics
  # -------------------------------------------------------------------------
  diagnostics <- list(
    M = M,
    max_iter = max_iter,
    tol = tol,
    eps = eps,
    log_bounds = log_bounds,
    init = list(
      a0 = a0,
      b0 = b0,
      mu_init = mu_init,
      var_init = var_init,
      init_method = init_method,
      init_status = init_status,
      init_residual = init_residual
    ),
    optim = list(
      method = "L-BFGS-B",
      convergence = opt$convergence,
      message = opt$message,
      counts = opt$counts,
      value = opt$value,
      par = opt$par
    ),
    fallback_used = fallback_used,
    kl_init = kl0,
    kl_final = kl_opt
  )

  # -------------------------------------------------------------------------
  # Build Fit Summary
  # -------------------------------------------------------------------------
  fit <- list(
    mu_K = moments$mean,
    var_K = moments$var,
    kl = kl_opt,
    residual = kl_opt
  )

  # -------------------------------------------------------------------------
  # Construct Output Object
  # -------------------------------------------------------------------------
  result <- structure(
    list(
      a = a_opt,
      b = b_opt,
      J = J,
      target = target_info,
      method = "A2-KL",
      status = status,
      converged = isTRUE(converged) && !isTRUE(fallback_used),
      iterations = as.integer(opt$counts[["function"]]),
      termination = termination,
      fit = fit,
      diagnostics = diagnostics,
      trace = trace_df
    ),
    class = "DPprior_fit"
  )

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
#' verify_kl_divergence()
#'
#' @export
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
#' verify_a2_kl()
#'
#' @export
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

  test1 <- fit$converged || fit$status == "fallback_used"
  if (isTRUE(verbose)) {
    cat(sprintf("Test 1: Convergence for typical target (method='chisq')\n"))
    cat(sprintf("  Target: mu_K = %.1f, var_K = %.1f\n", target$mu_K, target$var_K))
    cat(sprintf("  Result: a = %.4f, b = %.4f\n", fit$a, fit$b))
    cat(sprintf("  Achieved: mu_K = %.4f, var_K = %.4f\n", fit$fit$mu_K, fit$fit$var_K))
    cat(sprintf("  KL = %.4e, status = '%s'\n", fit$fit$kl, fit$status))
    cat(sprintf("  [%s]\n\n", if (test1) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test1

  # Test 2: Reasonable KL divergence
  test2 <- fit$fit$kl < 0.1
  if (isTRUE(verbose)) {
    cat(sprintf("Test 2: Reasonable KL divergence (< 0.1)\n"))
    cat(sprintf("  KL = %.4e [%s]\n\n", fit$fit$kl, if (test2) "PASS" else "FAIL"))
  }
  all_pass <- all_pass && test2

  # Test 3: Moment approximation quality
  mean_err <- abs(fit$fit$mu_K - target$mu_K)
  var_err <- abs(fit$fit$var_K - target$var_K)
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

  test4 <- (fit2$converged || fit2$status == "fallback_used") && fit2$fit$kl < 0.5
  if (isTRUE(verbose)) {
    cat(sprintf("Test 4: Custom binomial-shaped PMF target (method='pmf')\n"))
    cat(sprintf("  Result: a = %.4f, b = %.4f\n", fit2$a, fit2$b))
    cat(sprintf("  KL = %.4e, status = '%s'\n", fit2$fit$kl, fit2$status))
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
  test5 <- !is.null(fit3) && is.finite(fit3$fit$kl)
  if (isTRUE(verbose)) {
    if (!is.null(fit3)) {
      cat(sprintf("  Result: a = %.4f, b = %.4f, KL = %.4e\n",
                  fit3$a, fit3$b, fit3$fit$kl))
      cat(sprintf("  Fallback used: %s\n", fit3$diagnostics$fallback_used))
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
#' verify_a2_kl_all()
#'
#' @export
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
