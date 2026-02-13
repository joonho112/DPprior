# =============================================================================
# Module 08: First Stick-Breaking Weight (w₁) Distribution
# =============================================================================
#
# This module implements the closed-form distribution of w₁, the first
# stick-breaking weight in the GEM(α) representation of a Dirichlet Process,
# when α ~ Gamma(a, b).
#
# Key Results (Lee, 2026, Section 4; Vicentini & Jermyn, 2025):
#   Conditional: w₁ | α ~ Beta(1, α)
#   Unconditional CDF: F(x | a, b) = 1 - (b / (b - log(1-x)))^a
#   Quantile: Q(u | a, b) = 1 - exp(b × [1 - (1-u)^{-1/a}])
#   Density: p(x | a, b) = a × b^a / [(1-x) × (b - log(1-x))^{a+1}]
#
# Numerical Notes:
#   - Uses log1p/expm1 for stability near x ~ 0 and x ~ 1
#   - Boundary behavior handled explicitly:
#       CDF: x<=0 -> 0, x>=1 -> 1
#       Survival: t<=0 -> 1, t>=1 -> 0
#       Density: outside (0,1) -> 0 (or -Inf on log-scale)
#
# Author: JoonHo Lee
# Date: December 2025
# Part of: DPprior R Package
# Dependencies: Module 00 (constants.R), Module 02 (quadrature.R)
# =============================================================================


# =============================================================================
# CDF Function
# =============================================================================

#' CDF of First Stick-Breaking Weight w₁
#'
#' Computes P(w₁ ≤ x | a, b) using the closed-form expression derived
#' by marginalizing over α ~ Gamma(a, b).
#'
#' @param x Numeric vector. Values outside the unit interval are allowed but
#'   are mapped to the boundary values of the CDF (0 for x ≤ 0, 1 for x ≥ 1).
#' @param a Numeric; shape parameter of the Gamma prior on α (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on α (b > 0).
#'
#' @return Numeric vector of CDF values F(x | a, b) with same length as x.
#'
#' @details
#' The unconditional CDF is given by:
#' \deqn{F_{w_1}(x | a, b) = 1 - \left(\frac{b}{b - \log(1-x)}\right)^a}
#'
#' The implementation uses \code{log1p} and \code{expm1} for numerical
#' stability, particularly when the CDF is close to 0 (small x).
#'
#' @section Interpretation:
#' The weight w₁ is in **GEM (size-biased) order**, not ranked by size.
#' It represents the asymptotic cluster share of a randomly chosen unit,
#' **not** the largest cluster proportion. See Lee (2026, Section 4) for details.
#'
#' @examples
#' # P(w₁ ≤ 0.3) under standard prior
#' cdf_w1(0.3, a = 2, b = 1)
#'
#' # Vectorized computation
#' cdf_w1(c(0.1, 0.3, 0.5, 0.7), a = 1.6, b = 1.22)
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' Vicentini, C. and Jermyn, I. H. (2025). Prior selection for the precision
#'   parameter of Dirichlet Process Mixtures. arXiv:2502.00864.
#'
#' @seealso \code{\link{quantile_w1}}, \code{\link{density_w1}},
#'   \code{\link{prob_w1_exceeds}}
#'
#' @family weights_w1
#'
#' @export
cdf_w1 <- function(x, a, b) {
  # Input validation
  assert_positive(a, "a")
  assert_positive(b, "b")

  # Initialize output with NA
  out <- rep(NA_real_, length(x))

  # Handle boundaries explicitly
  out[x <= 0] <- 0
  out[x >= 1] <- 1

  # Identify valid interior points

  idx <- which(is.finite(x) & x > 0 & x < 1)
  if (length(idx) == 0L) {
    return(out)
  }

  xx <- x[idx]

  # Numerically stable computation using log1p and expm1
  # log(1-x) via log1p for stability when x is small
  log1m <- log1p(-xx)        # log(1-x) <= 0
  denom <- b - log1m         # b - log(1-x) > 0

  # log survival = a * (log(b) - log(denom))
  log_surv <- a * (log(b) - log(denom))

  # CDF = 1 - exp(log_surv) via expm1 for stability when CDF is near 0
  out[idx] <- -expm1(log_surv)

  out
}


# =============================================================================
# Quantile Function
# =============================================================================

#' Quantile Function of w₁
#'
#' Computes the inverse CDF: Q(u | a, b) = F⁻¹(u).
#'
#' @param u Numeric vector of probability levels in the unit interval.
#'   Values u ≤ 0 return 0 and u ≥ 1 return 1.
#' @param a Numeric; shape parameter of the Gamma prior on α (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on α (b > 0).
#'
#' @return Numeric vector of quantile values Q(u | a, b).
#'
#' @details
#' The quantile function has the closed form:
#' \deqn{Q_{w_1}(u | a, b) = 1 - \exp\left(b \left[1 - (1-u)^{-1/a}\right]\right)}
#'
#' The implementation computes (1-u)^(-1/a) in log space for stability
#' when u is close to 1.
#'
#' \strong{Numerical Note:} For small values of a (a < 1) and u close to 1,
#' the quantile approaches 1 very rapidly and may round to 1.0 in double
#' precision.
#'
#' @examples
#' # Median of w₁
#' quantile_w1(0.5, a = 2, b = 1)  # ~0.339
#'
#' # 90th percentile
#' quantile_w1(0.9, a = 2, b = 1)  # ~0.732
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso \code{\link{cdf_w1}}, \code{\link{summary_w1}}
#'
#' @family weights_w1
#'
#' @export
quantile_w1 <- function(u, a, b) {
  # Input validation
  assert_positive(a, "a")
  assert_positive(b, "b")
  assert_probability(u, "u")

  # Initialize output with NA
  out <- rep(NA_real_, length(u))

  # Handle boundaries explicitly
  out[u <= 0] <- 0
  out[u >= 1] <- 1

  # Identify valid interior points
  idx <- which(is.finite(u) & u > 0 & u < 1)
  if (length(idx) == 0L) {
    return(out)
  }

  uu <- u[idx]

  # (1-u)^(-1/a) computed in log space for stability
  # When u is close to 1, log1p(-u) = log(1-u) is very negative,
  # and (-1/a) * log(1-u) can be very large
  pow_term <- exp((-1 / a) * log1p(-uu))
  exponent <- b * (1 - pow_term)

  out[idx] <- 1 - exp(exponent)

  out
}


# =============================================================================
# Survival Function
# =============================================================================

#' Survival Function of w₁
#'
#' Computes P(w₁ > t | a, b) = 1 - F(t), the probability that the first
#' stick-breaking weight exceeds threshold t.
#'
#' @param t Numeric vector of thresholds. Values outside the unit interval are
#'   allowed but are mapped to the boundary values (1 for t ≤ 0, 0 for t ≥ 1).
#' @param a Numeric; shape parameter of the Gamma prior on α (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on α (b > 0).
#'
#' @return Numeric vector of survival probabilities.
#'
#' @details
#' The survival function has the closed form:
#' \deqn{P(w_1 > t | a, b) = \left(\frac{b}{b - \log(1-t)}\right)^a}
#'
#' This is a \strong{key quantity for dominance risk assessment} (Lee, 2026, Section 4).
#' A large P(w₁ > 0.5) indicates high prior probability that a single
#' cluster dominates the mixture.
#'
#' @section Dominance Risk Interpretation:
#' \itemize{
#'   \item P(w₁ > 0.5) ≈ 0.5: moderate dominance risk
#'   \item P(w₁ > 0.9) ≈ 0.1: low extreme dominance risk
#' }
#'
#' @examples
#' # P(w₁ > 0.5): "dominant cluster" probability
#' prob_w1_exceeds(0.5, a = 1.6, b = 1.22)  # ~0.487 (Lee et al. DP-inform)
#' prob_w1_exceeds(0.5, a = 2, b = 1)       # ~0.349
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso \code{\link{cdf_w1}}, \code{\link{summary_w1}}
#'
#' @family weights_w1
#'
#' @export
prob_w1_exceeds <- function(t, a, b) {
  # Input validation
  assert_positive(a, "a")
  assert_positive(b, "b")

  # Initialize output with NA
  out <- rep(NA_real_, length(t))

  # Handle boundaries explicitly
  out[t <= 0] <- 1
  out[t >= 1] <- 0

  # Identify valid interior points
  idx <- which(is.finite(t) & t > 0 & t < 1)
  if (length(idx) == 0L) {
    return(out)
  }

  tt <- t[idx]

  # Numerically stable computation
  log1m <- log1p(-tt)         # log(1-t)
  denom <- b - log1m          # b - log(1-t) > 0

  log_surv <- a * (log(b) - log(denom))
  out[idx] <- exp(log_surv)

  out
}


# =============================================================================
# Density Function
# =============================================================================

#' Density of w₁
#'
#' Computes the probability density p(w₁ = x | a, b).
#'
#' @param x Numeric vector; evaluation points.
#' @param a Numeric; shape parameter of the Gamma prior on α (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on α (b > 0).
#' @param log Logical; if \code{TRUE}, returns log-density. Default is \code{FALSE}.
#'
#' @return Numeric vector of density (or log-density if \code{log = TRUE}) values.
#'   Returns 0 (or -Inf on log scale) for x outside (0, 1).
#'
#' @details
#' The marginal density of w₁ is:
#' \deqn{p(w_1 | a, b) = \frac{a \cdot b^a}{(1-w_1) \cdot [b - \log(1-w_1)]^{a+1}}}
#'
#' \strong{Important:} For small values of a (a < 1), the density has significant
#' mass concentrated very close to x = 1.
#'
#' @examples
#' # Density at several points
#' x <- seq(0.1, 0.9, by = 0.1)
#' density_w1(x, a = 2, b = 1)
#'
#' # Log-density for numerical stability
#' density_w1(0.5, a = 2, b = 1, log = TRUE)
#'
#' @seealso \code{\link{cdf_w1}}, \code{\link{quantile_w1}}
#'
#' @family weights_w1
#'
#' @export
density_w1 <- function(x, a, b, log = FALSE) {
  # Input validation
  assert_positive(a, "a")
  assert_positive(b, "b")

  # Initialize output with NA
  out <- rep(NA_real_, length(x))

  # Handle boundaries explicitly
  if (log) {
    out[x <= 0] <- -Inf
    out[x >= 1] <- -Inf
  } else {
    out[x <= 0] <- 0
    out[x >= 1] <- 0
  }

  # Identify valid interior points
  idx <- which(is.finite(x) & x > 0 & x < 1)
  if (length(idx) == 0L) {
    return(out)
  }

  xx <- x[idx]

  # Numerically stable log-density computation
  log1m <- log1p(-xx)         # log(1-x)
  log_term <- -log1m          # -log(1-x)
  denom <- b + log_term       # b - log(1-x)

  log_density <- log(a) + a * log(b) - log1m - (a + 1) * log(denom)

  if (log) {
    out[idx] <- log_density
  } else {
    out[idx] <- exp(log_density)
  }

  out
}


# =============================================================================
# Moment Functions
# =============================================================================

#' Mean of w₁
#'
#' Computes E(w₁ | a, b) via Gauss-Laguerre quadrature.
#'
#' @param a Numeric; shape parameter of the Gamma prior on α (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on α (b > 0).
#' @param M Integer; number of quadrature nodes. Default is 80.
#'
#' @return Numeric; E(w₁).
#'
#' @details
#' The expectation is computed using the identity:
#' \deqn{E[w_1 | a, b] = E\left[\frac{1}{1+\alpha}\right] = I_1(a, b)}
#'
#' where the integral is evaluated via Gauss-Laguerre quadrature.
#'
#' \strong{Key identity:} E(w₁ | a, b) = E(ρ | a, b) where ρ = Σwₕ² is
#' the co-clustering probability.
#'
#' @examples
#' mean_w1(a = 2, b = 1)       # ~0.404
#' mean_w1(a = 1.6, b = 1.22)  # ~0.508
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso \code{\link{var_w1}}, \code{\link{summary_w1}}
#'
#' @family weights_w1
#'
#' @export
mean_w1 <- function(a, b, M = .QUAD_NODES_DEFAULT) {
  # Input validation
  assert_positive(a, "a")
  assert_positive(b, "b")

  # E[w₁] = E[1/(1+α)] where α ~ Gamma(a, b)
  integrate_gamma(function(alpha) 1 / (1 + alpha), a, b, M)
}


#' Variance of w₁
#'
#' Computes Var(w₁ | a, b) using the law of total variance.
#'
#' @param a Numeric; shape parameter of the Gamma prior on α (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on α (b > 0).
#' @param M Integer; number of quadrature nodes. Default is 80.
#'
#' @return Numeric; Var(w₁).
#'
#' @details
#' Uses the law of total variance:
#' \deqn{Var(w_1) = E[Var(w_1 | \alpha)] + Var(E[w_1 | \alpha])}
#'
#' where w₁ | α ~ Beta(1, α), so:
#' \itemize{
#'   \item E(w₁ | α) = 1/(1+α)
#'   \item Var(w₁ | α) = α / ((1+α)²(2+α))
#' }
#'
#' @examples
#' var_w1(a = 2, b = 1)       # ~0.090
#' var_w1(a = 1.6, b = 1.22)  # ~0.105
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso \code{\link{mean_w1}}, \code{\link{summary_w1}}
#'
#' @family weights_w1
#'
#' @export
var_w1 <- function(a, b, M = .QUAD_NODES_DEFAULT) {
  # Input validation
  assert_positive(a, "a")
  assert_positive(b, "b")

  # Component 1: E[Var(w₁|α)] where Var(w₁|α) = α / ((1+α)²(2+α))
  E_var_cond <- integrate_gamma(
    function(alpha) alpha / ((1 + alpha)^2 * (2 + alpha)),
    a, b, M
  )

  # Component 2: Var(E[w₁|α]) = E[(1/(1+α))²] - (E[1/(1+α)])²
  E_mean_sq <- integrate_gamma(
    function(alpha) 1 / (1 + alpha)^2,
    a, b, M
  )
  E_mean <- mean_w1(a, b, M)

  # Law of total variance
  E_var_cond + (E_mean_sq - E_mean^2)
}


# =============================================================================
# Summary Function
# =============================================================================

#' Summary Statistics for w₁ Distribution
#'
#' Computes comprehensive summary statistics for the w₁ distribution.
#'
#' @param a Numeric; shape parameter of the Gamma prior on α (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on α (b > 0).
#' @param probs Numeric vector; quantile probabilities. Default is
#'   \code{c(0.05, 0.25, 0.5, 0.75, 0.95)}.
#' @param M Integer; number of quadrature nodes for mean/variance. Default is 80.
#'
#' @return A list of class "w1_summary" containing:
#'   \describe{
#'     \item{mean}{E(w₁)}
#'     \item{var}{Var(w₁)}
#'     \item{sd}{SD(w₁) = sqrt(Var(w₁))}
#'     \item{median}{Median of w₁}
#'     \item{quantiles}{Named vector of quantiles}
#'     \item{prob_gt_50}{P(w₁ > 0.5), dominance indicator}
#'     \item{prob_gt_90}{P(w₁ > 0.9), extreme dominance indicator}
#'     \item{params}{List of input parameters (a, b)}
#'   }
#'
#' @examples
#' # Standard summary
#' summary_w1(a = 2, b = 1)
#'
#' # Lee et al. DP-inform prior
#' summary_w1(a = 1.6, b = 1.22)
#'
#' @seealso \code{\link{cdf_w1}}, \code{\link{quantile_w1}}, \code{\link{mean_w1}}
#'
#' @family weights_w1
#'
#' @export
summary_w1 <- function(a, b, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                       M = .QUAD_NODES_DEFAULT) {
  # Input validation
  assert_positive(a, "a")
  assert_positive(b, "b")
  assert_probability(probs, "probs")

  # Compute quantiles
  quantiles <- quantile_w1(probs, a, b)
  names(quantiles) <- paste0("q", 100 * probs)

  # Compute moments
  m <- mean_w1(a, b, M)
  v <- var_w1(a, b, M)

  result <- list(
    mean = m,
    var = v,
    sd = sqrt(v),
    median = unname(quantiles["q50"]),
    quantiles = quantiles,
    prob_gt_50 = prob_w1_exceeds(0.5, a, b),
    prob_gt_90 = prob_w1_exceeds(0.9, a, b),
    params = list(a = a, b = b)
  )

  class(result) <- "w1_summary"
  result
}


#' Print Method for w1_summary Objects
#'
#' @param x An object of class "w1_summary".
#' @param digits Integer; number of digits for printing.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.w1_summary <- function(x, digits = 4, ...) {
  cat("w1 Distribution Summary\n")
  cat(strrep("=", 45), "\n\n")

  cat(sprintf("Gamma prior: alpha ~ Gamma(%.4f, %.4f)\n",
              x$params$a, x$params$b))
  cat(sprintf("E[alpha] = %.4f, CV(alpha) = %.2f%%\n\n",
              x$params$a / x$params$b, 100 / sqrt(x$params$a)))

  cat("Location and Scale:\n")
  cat(strrep("-", 30), "\n")
  cat(sprintf("  Mean:   %.*f\n", digits, x$mean))
  cat(sprintf("  Median: %.*f\n", digits, x$median))
  cat(sprintf("  SD:     %.*f\n", digits, x$sd))

  cat("\nQuantiles:\n")
  cat(strrep("-", 30), "\n")
  q_str <- paste(sprintf("  %s: %.*f", names(x$quantiles), digits, x$quantiles),
                 collapse = "\n")
  cat(q_str, "\n")

  cat("\nDominance Risk:\n")
  cat(strrep("-", 30), "\n")
  cat(sprintf("  P(w1 > 0.5): %.*f\n", digits, x$prob_gt_50))
  cat(sprintf("  P(w1 > 0.9): %.*f\n", digits, x$prob_gt_90))

  invisible(x)
}


# =============================================================================
# Utility: Random Generation (for Monte Carlo validation)
# =============================================================================

#' Random Generation from w₁ Distribution
#'
#' Generates random samples from the w₁ distribution by first sampling
#' α ~ Gamma(a, b), then w₁ | α ~ Beta(1, α).
#'
#' @param n Integer; number of samples to generate.
#' @param a Numeric; shape parameter of the Gamma prior on α (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on α (b > 0).
#'
#' @return Numeric vector of length n; random samples from the w₁ distribution.
#'
#' @details
#' This function uses the hierarchical representation:
#' \enumerate{
#'   \item α ~ Gamma(a, b)
#'   \item w₁ | α ~ Beta(1, α)
#' }
#'
#' Useful for Monte Carlo validation of the closed-form functions.
#'
#' @examples
#' # Generate samples
#' set.seed(42)
#' samples <- rw1(10000, a = 2, b = 1)
#'
#' # Compare empirical vs theoretical mean
#' mean(samples)          # ~0.404
#' mean_w1(a = 2, b = 1)  # 0.4037
#'
#' @export
rw1 <- function(n, a, b) {
  # Input validation
  assert_positive(a, "a")
  assert_positive(b, "b")
  if (!is.numeric(n) || length(n) != 1L || n < 1 || n != floor(n)) {
    stop("n must be a positive integer", call. = FALSE)
  }

  # Hierarchical sampling: α ~ Gamma(a, b), then w₁ | α ~ Beta(1, α)
  alpha <- stats::rgamma(n, shape = a, rate = b)
  stats::rbeta(n, shape1 = 1, shape2 = alpha)
}


# =============================================================================
# Utility: Grid Computation for Plotting
# =============================================================================

#' Compute w₁ Distribution on Grid
#'
#' Computes CDF, PDF, and survival function on a grid of x values.
#' Useful for visualization and comparison across different priors.
#'
#' @param a Numeric; shape parameter of the Gamma prior on α.
#' @param b Numeric; rate parameter of the Gamma prior on α.
#' @param x_grid Numeric vector; grid of x values in (0, 1).
#'   Default is \code{seq(0.01, 0.99, length.out = 100)}.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{x}{Grid points}
#'     \item{cdf}{CDF values F(x)}
#'     \item{pdf}{Density values p(x)}
#'     \item{survival}{Survival function S(x) = 1 - F(x)}
#'   }
#'
#' @examples
#' # Compute on default grid
#' df <- w1_grid(a = 2, b = 1)
#'
#' # Plot all three functions
#' par(mfrow = c(1, 3))
#' plot(df$x, df$cdf, type = "l", main = "CDF", xlab = "x", ylab = "F(x)")
#' plot(df$x, df$pdf, type = "l", main = "PDF", xlab = "x", ylab = "p(x)")
#' plot(df$x, df$survival, type = "l", main = "Survival", xlab = "x", ylab = "S(x)")
#'
#' @export
w1_grid <- function(a, b, x_grid = seq(0.01, 0.99, length.out = 100)) {
  assert_positive(a, "a")
  assert_positive(b, "b")

  data.frame(
    x = x_grid,
    cdf = cdf_w1(x_grid, a, b),
    pdf = density_w1(x_grid, a, b),
    survival = prob_w1_exceeds(x_grid, a, b)
  )
}
