# =============================================================================
# Module 09: Co-Clustering Probability (rho) Distribution
# =============================================================================
#
# Author: JoonHo Lee
# Date: December 2025
# Part of: DPprior R Package
#
# This module implements moments of the co-clustering probability
#   rho = sum_{h>=1} w_h^2,
# where (w_h) are the stick-breaking weights under GEM(alpha).
#
# Probabilistic interpretation:
#   rho = P(Z1 = Z2 | w) where Z1, Z2 are cluster labels for two random units.
#
# Conditional moments given alpha (Lee, 2026, Section 4):
#   E[rho | alpha] = 1 / (1 + alpha)
#   E[rho^2 | alpha] = (alpha + 6) / ((alpha+1)(alpha+2)(alpha+3))
#   Var(rho | alpha) = 2 alpha / ((alpha+1)^2 (alpha+2)(alpha+3))
#
# Key identity: E[rho | alpha] = E[w1 | alpha] = 1/(1+alpha)
#   Therefore E[rho | a, b] = E[w1 | a, b], but full distributions differ.
#
# Marginal moments under alpha ~ Gamma(a, b) are computed via Gauss-Laguerre
# quadrature using integrate_gamma() (Module 02).
# =============================================================================


# =============================================================================
# Conditional Moments Given alpha
# =============================================================================

#' Conditional Mean of rho Given Alpha
#'
#' Computes the conditional mean E(rho | alpha) for the co-clustering
#' probability rho = sum_h w_h^2 under a Dirichlet Process.
#'
#' @param alpha Numeric vector; concentration parameter(s) (must be positive).
#'
#' @return Numeric vector; `E(rho | alpha) = 1/(1+alpha)`.
#'
#' @details
#' The co-clustering probability rho = sum(w_h^2) over h >= 1 has conditional mean:
#' \deqn{E[\rho | \alpha] = \frac{1}{1 + \alpha}}
#'
#' This equals `E(w1 | alpha)` since w1 ~ Beta(1, alpha) has mean 1/(1+alpha).
#'
#' \strong{Interpretation:}
#' \itemize{
#'   \item alpha -> 0: E(rho|alpha) -> 1 (all observations in one cluster)
#'   \item alpha -> Inf: E(rho|alpha) -> 0 (infinitely many small clusters)
#'   \item alpha = 1: E(rho|alpha) = 0.5 (moderate clustering)
#' }
#'
#' @examples
#' mean_rho_given_alpha(1.0)
#' mean_rho_given_alpha(c(0.5, 1, 2, 5, 10))
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso \code{\link{var_rho_given_alpha}}, \code{\link{mean_rho}}
#'
#' @family co_clustering
#'
#' @export
mean_rho_given_alpha <- function(alpha) {

  assert_positive(alpha, "alpha")
  1 / (1 + alpha)
}


#' Conditional Second Moment of rho Given Alpha
#'
#' Computes `E(rho^2 | alpha) = (alpha + 6) / ((alpha+1)(alpha+2)(alpha+3))`.
#'
#' @param alpha Numeric vector; concentration parameter(s) (must be positive).
#'
#' @return Numeric vector; `E(rho^2 | alpha)`.
#'
#' @details
#' Used internally for variance computation via the identity:
#' `Var(rho|alpha) = E(rho^2|alpha) - E(rho|alpha)^2`
#'
#' @keywords internal
mean_rho_sq_given_alpha <- function(alpha) {
  assert_positive(alpha, "alpha")
  (alpha + 6) / ((alpha + 1) * (alpha + 2) * (alpha + 3))
}


#' Conditional Variance of rho Given Alpha
#'
#' Computes the conditional variance Var(rho | alpha).
#'
#' @param alpha Numeric vector; concentration parameter(s) (must be positive).
#'
#' @return Numeric vector; Var(rho | alpha).
#'
#' @details
#' The conditional variance is:
#' \deqn{Var(\rho | \alpha) = \frac{2\alpha}{(1+\alpha)^2(2+\alpha)(3+\alpha)}}
#'
#' This is derived from the GEM recursion:
#' rho = V^2 + (1-V)^2 * rho' where V ~ Beta(1, alpha) and rho' is an
#' independent copy of rho.
#'
#' \strong{Properties:}
#' \itemize{
#'   \item Var(rho|alpha) = 0 when alpha -> 0 (degenerate at rho = 1)
#'   \item Var(rho|alpha) -> 0 when alpha -> Inf (degenerate at rho = 0)
#'   \item Maximum variance occurs at intermediate alpha
#' }
#'
#' @examples
#' var_rho_given_alpha(2)
#' var_rho_given_alpha(c(0.5, 1, 2, 5, 10))
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso \code{\link{mean_rho_given_alpha}}, \code{\link{var_rho}}
#'
#' @family co_clustering
#'
#' @export
var_rho_given_alpha <- function(alpha) {
  assert_positive(alpha, "alpha")
  2 * alpha / ((alpha + 1)^2 * (alpha + 2) * (alpha + 3))
}


# =============================================================================
# Marginal Moments (alpha ~ Gamma(a, b))
# =============================================================================

#' Marginal Mean of rho
#'
#' Computes E(rho | a, b) when alpha ~ Gamma(a, b) (shape-rate).
#'
#' @param a Numeric; shape parameter of the Gamma prior on alpha (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on alpha (b > 0).
#' @param M Integer; number of quadrature nodes. Default is 80.
#'
#' @return Numeric; `E(rho | a, b)`.
#'
#' @details
#' Uses Gauss-Laguerre quadrature via \code{integrate_gamma}. A key identity is
#' `E(rho | alpha) = E(w1 | alpha) = 1/(1+alpha)`, so `E(rho | a, b)` equals
#' `E(w1 | a, b)` (but the full distributions differ).
#'
#' \strong{Interpretation:}
#' \itemize{
#'   \item E(rho) > 0.5: High prior co-clustering probability
#'   \item E(rho) in (0.2, 0.5): Moderate co-clustering
#'   \item E(rho) < 0.2: Low co-clustering (fragmented prior)
#' }
#'
#' @examples
#' mean_rho(a = 2, b = 1)
#' mean_rho(a = 1.6, b = 1.22)
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso \code{\link{var_rho}}, \code{\link{cv_rho}}, \code{\link{mean_w1}}
#'
#' @family co_clustering
#'
#' @export
mean_rho <- function(a, b, M = .QUAD_NODES_DEFAULT) {
  assert_positive(a, "a")
  assert_positive(b, "b")

  integrate_gamma(function(alpha) 1 / (1 + alpha), a, b, M)
}


#' Marginal Second Moment of rho
#'
#' Computes `E(rho^2 | a, b)` by mixing `E(rho^2 | alpha)` over alpha ~ Gamma(a, b).
#'
#' @param a Numeric; shape parameter of the Gamma prior on alpha (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on alpha (b > 0).
#' @param M Integer; number of quadrature nodes. Default is 80.
#'
#' @return Numeric; `E(rho^2 | a, b)`.
#'
#' @keywords internal
mean_rho_sq <- function(a, b, M = .QUAD_NODES_DEFAULT) {
  assert_positive(a, "a")
  assert_positive(b, "b")

  integrate_gamma(mean_rho_sq_given_alpha, a, b, M)
}


#' Marginal Variance of rho
#'
#' Computes Var(rho | a, b) when alpha ~ Gamma(a, b) (shape-rate).
#'
#' @param a Numeric; shape parameter of the Gamma prior on alpha (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on alpha (b > 0).
#' @param M Integer; number of quadrature nodes. Default is 80.
#'
#' @return Numeric; Var(rho | a, b).
#'
#' @details
#' Uses the law of total variance:
#' \deqn{Var(\rho | a, b) = E[Var(\rho | \alpha)] + Var(E[\rho | \alpha])}
#'
#' where:
#' \itemize{
#'   \item `Var(rho | alpha) = 2*alpha / ((1+alpha)^2*(2+alpha)*(3+alpha))`
#'   \item `E(rho | alpha) = 1/(1+alpha)`
#' }
#'
#' \strong{Note:} Unlike E(rho), Var(rho) != Var(w1) in general, because
#' the conditional variances differ.
#'
#' @examples
#' var_rho(a = 2, b = 1)
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso \code{\link{mean_rho}}, \code{\link{cv_rho}}, \code{\link{var_w1}}
#'
#' @family co_clustering
#'
#' @export
var_rho <- function(a, b, M = .QUAD_NODES_DEFAULT) {
  assert_positive(a, "a")
  assert_positive(b, "b")

  # E[Var(rho|alpha)]
  E_var_cond <- integrate_gamma(var_rho_given_alpha, a, b, M)

  # Var(E[rho|alpha]) = E[(1/(1+alpha))^2] - (E[1/(1+alpha)])^2
  E_mean_sq <- integrate_gamma(function(alpha) (1 / (1 + alpha))^2, a, b, M)
  E_mean <- mean_rho(a, b, M)
  var_mean_cond <- E_mean_sq - E_mean^2

  E_var_cond + var_mean_cond
}


#' Coefficient of Variation of rho
#'
#' Computes CV(rho) = SD(rho) / E(rho) under alpha ~ Gamma(a, b).
#'
#' @param a Numeric; shape parameter of the Gamma prior on alpha (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on alpha (b > 0).
#' @param M Integer; number of quadrature nodes. Default is 80.
#'
#' @return Numeric; coefficient of variation.
#'
#' @examples
#' cv_rho(a = 2, b = 1)
#' cv_rho(a = 1.6, b = 1.22)
#'
#' @seealso \code{\link{mean_rho}}, \code{\link{var_rho}}
#'
#' @family co_clustering
#'
#' @export
cv_rho <- function(a, b, M = .QUAD_NODES_DEFAULT) {
  sqrt(var_rho(a, b, M)) / mean_rho(a, b, M)
}


# =============================================================================
# Summary and Interpretation Functions
# =============================================================================

#' Summary Statistics for rho Distribution
#'
#' Computes comprehensive summary statistics for the co-clustering probability
#' rho under the hierarchical prior alpha ~ Gamma(a, b).
#'
#' @param a Numeric; shape parameter of the Gamma prior on alpha (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on alpha (b > 0).
#' @param M Integer; number of quadrature nodes. Default is 80.
#'
#' @return A list of class "rho_summary" containing:
#'   \describe{
#'     \item{mean}{E(rho | a, b)}
#'     \item{var}{Var(rho | a, b)}
#'     \item{sd}{SD(rho | a, b) = sqrt(Var)}
#'     \item{cv}{Coefficient of variation SD/mean}
#'     \item{interpretation}{Qualitative interpretation of co-clustering level}
#'     \item{params}{List of input parameters (a, b)}
#'     \item{alpha_prior}{Summary of the alpha prior (mean, sd, cv)}
#'     \item{conditional_at_alpha_mean}{Conditional moments evaluated at E(alpha)}
#'   }
#'
#' @details
#' The co-clustering probability rho indicates how likely two randomly
#' chosen observations are to belong to the same cluster a priori.
#'
#' \strong{Interpretation guidelines:}
#' \itemize{
#'   \item E(rho) > 0.5: High co-clustering; most pairs expected in same cluster
#'   \item E(rho) in (0.2, 0.5): Moderate co-clustering
#'   \item E(rho) < 0.2: Low co-clustering; fragmented prior
#' }
#'
#' The \code{conditional_at_alpha_mean} component provides a "plug-in" estimate
#' for comparison: what the moments would be if alpha were fixed at its prior mean.
#'
#' @examples
#' summary_rho(a = 2, b = 1)
#' summary_rho(a = 1.6, b = 1.22)
#'
#' @seealso \code{\link{mean_rho}}, \code{\link{var_rho}}, \code{\link{summary_w1}}
#'
#' @export
summary_rho <- function(a, b, M = .QUAD_NODES_DEFAULT) {
  assert_positive(a, "a")
  assert_positive(b, "b")

  # Marginal moments

  mean_val <- mean_rho(a, b, M)
  var_val <- var_rho(a, b, M)
  sd_val <- sqrt(var_val)

  # Interpretation based on E[rho]
  interpretation <- if (mean_val > 0.5) {
    "High co-clustering: most pairs expected in same cluster"
  } else if (mean_val > 0.2) {
    "Moderate co-clustering"
  } else {
    "Low co-clustering: fragmented prior"
  }

  # alpha prior summary
  alpha_mean <- a / b
  alpha_sd <- sqrt(a) / b
  alpha_cv <- 1 / sqrt(a)

  # Conditional moments at E[alpha] (plug-in estimate)
  cond_mean <- mean_rho_given_alpha(alpha_mean)
  cond_var <- var_rho_given_alpha(alpha_mean)

  result <- list(
    mean = mean_val,
    var = var_val,
    sd = sd_val,
    cv = sd_val / mean_val,
    interpretation = interpretation,
    params = list(a = a, b = b),
    alpha_prior = list(
      mean = alpha_mean,
      sd = alpha_sd,
      cv = alpha_cv
    ),
    conditional_at_alpha_mean = list(
      alpha = alpha_mean,
      mean = cond_mean,
      var = cond_var
    )
  )

  class(result) <- "rho_summary"
  result
}


#' Print Method for rho_summary Objects
#'
#' @param x An object of class "rho_summary".
#' @param digits Integer; number of digits for printing.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.rho_summary <- function(x, digits = 4, ...) {
  cat("Co-Clustering Probability (rho) Summary\n")
  cat(strrep("=", 50), "\n\n")

  cat(sprintf("Gamma prior: alpha ~ Gamma(%.4f, %.4f)\n",
              x$params$a, x$params$b))
  cat(sprintf("E[alpha] = %.4f, SD(alpha) = %.4f, CV(alpha) = %.1f%%\n\n",
              x$alpha_prior$mean, x$alpha_prior$sd, 100 * x$alpha_prior$cv))

  cat("Marginal distribution of rho:\n")
  cat(strrep("-", 35), "\n")
  cat(sprintf("  Mean:   %.*f\n", digits, x$mean))
  cat(sprintf("  SD:     %.*f\n", digits, x$sd))
  cat(sprintf("  CV:     %.1f%%\n", 100 * x$cv))

  cat(sprintf("\nConditional at E[alpha] = %.4f (plug-in):\n", x$conditional_at_alpha_mean$alpha))
  cat(strrep("-", 35), "\n")
  cat(sprintf("  E[rho | E[alpha]]:   %.*f\n", digits, x$conditional_at_alpha_mean$mean))
  cat(sprintf("  Var(rho | E[alpha]): %.*f\n", digits, x$conditional_at_alpha_mean$var))

  cat("\nInterpretation:\n")
  cat(strrep("-", 35), "\n")
  cat(sprintf("  %s\n", x$interpretation))

  invisible(x)
}


# =============================================================================
# Verification Functions
# =============================================================================

#' Verify the Identity E(w1) = E(rho)
#'
#' Checks the mean identity `E(w1 | a, b) = E(rho | a, b)`, which follows from
#' `E(rho | alpha) = E(w1 | alpha) = 1/(1+alpha)`.
#'
#' @param a Numeric; shape parameter of the Gamma prior on alpha (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on alpha (b > 0).
#' @param tol Numeric; absolute tolerance. Default is 1e-10.
#' @param M Integer; number of quadrature nodes. Default is 80.
#'
#' @return Logical; TRUE if the identity holds within tolerance.
#'
#' @examples
#' \dontrun{
#' verify_w1_rho_identity(2, 1)
#' verify_w1_rho_identity(1.6, 1.22)
#' }
#'
#' @keywords internal
verify_w1_rho_identity <- function(a, b, tol = 1e-10, M = .QUAD_NODES_DEFAULT) {
  assert_positive(a, "a")
  assert_positive(b, "b")

  abs(mean_w1(a, b, M) - mean_rho(a, b, M)) < tol
}


#' Verify Variance Decomposition
#'
#' Verifies the law of total variance decomposition for Var(rho).
#'
#' @param a Numeric; shape parameter of the Gamma prior on alpha (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on alpha (b > 0).
#' @param tol Numeric; tolerance for comparison. Default is 1e-10.
#' @param M Integer; number of quadrature nodes. Default is 80.
#' @param verbose Logical; if TRUE, print detailed results.
#'
#' @return Logical; TRUE if decomposition holds within tolerance.
#'
#' @keywords internal
verify_rho_variance_decomposition <- function(a, b, tol = 1e-10,
                                              M = .QUAD_NODES_DEFAULT,
                                              verbose = FALSE) {
  # Total variance via var_rho()
  var_total <- var_rho(a, b, M)

  # E[Var(rho|alpha)]
  E_var_cond <- integrate_gamma(var_rho_given_alpha, a, b, M)

  # Var(E[rho|alpha]) = E[(1/(1+alpha))^2] - (E[1/(1+alpha)])^2
  E_mean_sq <- integrate_gamma(function(alpha) (1 / (1 + alpha))^2, a, b, M)
  E_mean <- mean_rho(a, b, M)
  var_mean_cond <- E_mean_sq - E_mean^2

  var_decomposed <- E_var_cond + var_mean_cond
  passed <- abs(var_total - var_decomposed) < tol

  if (verbose) {
    cat(sprintf("Variance decomposition (a=%.2f, b=%.2f):\n", a, b))
    cat(sprintf("  E[Var(rho|alpha)] = %.10f\n", E_var_cond))
    cat(sprintf("  Var(E[rho|alpha]) = %.10f\n", var_mean_cond))
    cat(sprintf("  Sum               = %.10f\n", var_decomposed))
    cat(sprintf("  Var(rho)          = %.10f\n", var_total))
    cat(sprintf("  Match: %s\n", if (passed) "PASS" else "FAIL"))
  }

  passed
}


#' Verify Conditional Variance Formula
#'
#' Verifies that `Var(rho|alpha) = E(rho^2|alpha) - E(rho|alpha)^2`.
#'
#' @param alpha Numeric; concentration parameter (must be positive).
#' @param tol Numeric; tolerance for comparison. Default is 1e-12.
#'
#' @return Logical; TRUE if formula holds.
#'
#' @keywords internal
verify_rho_conditional_variance <- function(alpha, tol = 1e-12) {
  E_rho <- mean_rho_given_alpha(alpha)
  E_rho_sq <- mean_rho_sq_given_alpha(alpha)
  var_direct <- var_rho_given_alpha(alpha)
  var_from_moments <- E_rho_sq - E_rho^2

  abs(var_direct - var_from_moments) < tol
}


# =============================================================================
# Comparison Functions
# =============================================================================

#' Compare rho and w1 Distributions
#'
#' Compares the marginal distributions of rho and w1 under the same
#' hyperprior alpha ~ Gamma(a, b).
#'
#' @param a Numeric; shape parameter of the Gamma prior on alpha (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on alpha (b > 0).
#' @param M Integer; number of quadrature nodes. Default is 80.
#'
#' @return A list containing:
#'   \describe{
#'     \item{mean_rho}{E(rho)}
#'     \item{mean_w1}{E(w1)}
#'     \item{mean_equal}{Logical; whether means are equal}
#'     \item{var_rho}{Var(rho)}
#'     \item{var_w1}{Var(w1)}
#'     \item{var_ratio}{Var(rho) / Var(w1)}
#'   }
#'
#' @details
#' While E(rho) = E(w1), the variances differ because:
#' \itemize{
#'   \item `Var(rho | alpha) = 2*alpha / ((1+alpha)^2*(2+alpha)*(3+alpha))`
#'   \item `Var(w1 | alpha) = alpha / ((1+alpha)^2*(2+alpha))`
#' }
#'
#' Generally, Var(rho) < Var(w1) because rho averages over all squared weights.
#'
#' @examples
#' \dontrun{
#' compare_rho_w1(a = 2, b = 1)
#' compare_rho_w1(a = 1.6, b = 1.22)
#'
#' }
#' @keywords internal
compare_rho_w1 <- function(a, b, M = .QUAD_NODES_DEFAULT) {
  assert_positive(a, "a")
  assert_positive(b, "b")

  mean_rho_val <- mean_rho(a, b, M)
  mean_w1_val <- mean_w1(a, b, M)
  var_rho_val <- var_rho(a, b, M)
  var_w1_val <- var_w1(a, b, M)

  list(
    mean_rho = mean_rho_val,
    mean_w1 = mean_w1_val,
    mean_equal = abs(mean_rho_val - mean_w1_val) < 1e-10,
    var_rho = var_rho_val,
    var_w1 = var_w1_val,
    var_ratio = var_rho_val / var_w1_val
  )
}


# =============================================================================
# Random Generation (for Monte Carlo validation)
# =============================================================================

#' Random Generation from rho Distribution
#'
#' Generates random samples from the rho = sum w_h^2 distribution by
#' stick-breaking simulation.
#'
#' @param n Integer; number of samples to generate.
#' @param a Numeric; shape parameter of the Gamma prior on alpha (a > 0).
#' @param b Numeric; rate parameter of the Gamma prior on alpha (b > 0).
#' @param n_sticks Integer; number of sticks for truncation. Default is 500.
#'
#' @return Numeric vector of length n; random samples from the rho distribution.
#'
#' @details
#' Uses the hierarchical representation:
#' \enumerate{
#'   \item alpha ~ Gamma(a, b)
#'   \item v_h | alpha ~ Beta(1, alpha) independently for h = 1, ..., n_sticks
#'   \item w_1 = v_1, w_h = v_h * prod(1 - v_l) for l < h
#'   \item rho = sum w_h^2
#' }
#'
#' Useful for Monte Carlo validation of analytical formulas.
#'
#' @examples
#' set.seed(42)
#' rho_samples <- rrho(1000, a = 2, b = 1)
#' mean(rho_samples)
#' mean_rho(a = 2, b = 1)
#'
#' @seealso \code{\link{rw1}} for w1 random generation
#'
#' @export
rrho <- function(n, a, b, n_sticks = 500L) {
  assert_positive(a, "a")
  assert_positive(b, "b")
  if (!is.numeric(n) || length(n) != 1L || n < 1 || n != floor(n)) {
    stop("n must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(n_sticks) || length(n_sticks) != 1L ||
      n_sticks < 10 || n_sticks != floor(n_sticks)) {
    stop("n_sticks must be a positive integer >= 10", call. = FALSE)
  }

  n_sticks <- as.integer(n_sticks)
  rho_samples <- numeric(n)

  for (i in seq_len(n)) {
    alpha <- stats::rgamma(1L, shape = a, rate = b)
    v <- stats::rbeta(n_sticks, shape1 = 1, shape2 = alpha)

    w <- numeric(n_sticks)
    remaining <- 1.0
    for (h in seq_len(n_sticks)) {
      w[h] <- v[h] * remaining
      remaining <- remaining * (1 - v[h])
    }

    rho_samples[i] <- sum(w^2)
  }

  rho_samples
}


# =============================================================================
# Grid Computation (for visualization)
# =============================================================================

#' Compute rho Conditional Moments on alpha Grid
#'
#' Computes conditional mean and variance of rho on a grid of alpha values.
#' Useful for visualization of how rho varies with alpha.
#'
#' @param alpha_grid Numeric vector; grid of alpha values.
#'   Default is \code{seq(0.1, 10, length.out = 100)}.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{alpha}{Grid points}
#'     \item{mean}{E(rho | alpha)}
#'     \item{var}{Var(rho | alpha)}
#'     \item{sd}{SD(rho | alpha)}
#'   }
#'
#' @examples
#' df <- rho_conditional_grid()
#' plot(df$alpha, df$mean, type = "l",
#'      xlab = expression(alpha), ylab = expression(E(rho)))
#'
#' @export
rho_conditional_grid <- function(alpha_grid = seq(0.1, 10, length.out = 100)) {
  assert_positive(alpha_grid, "alpha_grid")

  data.frame(
    alpha = alpha_grid,
    mean = mean_rho_given_alpha(alpha_grid),
    var = var_rho_given_alpha(alpha_grid),
    sd = sqrt(var_rho_given_alpha(alpha_grid))
  )
}
