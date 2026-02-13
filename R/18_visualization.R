# =============================================================================
# Module 18: Visualization (Final Version)
# =============================================================================
#
# Publication-ready visualization functions for DPprior.
#
# Key features:
# - Direct parameter API: plot_alpha_prior(a = 1.6, b = 1.2)
# - Hierarchical PMF fallbacks (4 levels)
# - No annotation boxes (stats in subtitle only)
# - gtable dashboard (no patchwork dependency)
# - Color-coded dominance risk badge
#
# Author: JoonHo Lee
# Date: December 2025
# Part of: DPprior R Package
# =============================================================================

# =============================================================================
# Internal Utilities
# =============================================================================

.dpprior_has_ggplot2 <- function() {
  requireNamespace("ggplot2", quietly = TRUE)
}

.dpprior_require_ggplot2 <- function() {
  if (!.dpprior_has_ggplot2()) {
    stop("ggplot2 is required for engine = 'ggplot2'. Install ggplot2 or use engine = 'base'.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.dpprior_coalesce <- function(x, y) {
  if (is.null(x)) y else x
}

.dpprior_get_quad_nodes_default <- function() {
  if (exists(".QUAD_NODES_DEFAULT", inherits = TRUE)) {
    get(".QUAD_NODES_DEFAULT", inherits = TRUE)
  } else {
    80L
  }
}

.dpprior_assert_fit <- function(fit) {
  if (is.null(fit) || !is.list(fit)) {
    stop("fit must be a list-like DPprior_fit object.", call. = FALSE)
  }
  if (is.null(fit$a) || is.null(fit$b) || is.null(fit$J)) {
    stop("fit must contain fields: a, b, and J.", call. = FALSE)
  }
  invisible(TRUE)
}

.dpprior_extract_abJ <- function(fit) {
  .dpprior_assert_fit(fit)
  list(
    a = as.numeric(fit$a),
    b = as.numeric(fit$b),
    J = as.integer(fit$J)
  )
}

# =============================================================================
# Theme and Colors
# =============================================================================

#' Publication-Quality Theme for DPprior Plots
#'
#' @param base_size Numeric; base font size (default: 11).
#' @param base_family Character; base font family (default: "").
#' @return A ggplot2 theme object.
#'
#' @family visualization
#'
#' @export
theme_DPprior <- function(base_size = 11, base_family = "") {
  .dpprior_require_ggplot2()
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey40",
                                            margin = ggplot2::margin(b = 8)),
      plot.caption = ggplot2::element_text(color = "grey30", hjust = 0),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(face = "plain"),
      axis.text = ggplot2::element_text(color = "grey20"),
      legend.position = "none"
    )
}

#' DPprior Color Palette
#' @return A named list with color values.
#'
#' @family visualization
#'
#' @export
DPprior_colors <- function() {
  list(
    primary = "#4682B4",      # steelblue
    secondary = "#41B6C4",    # teal
    accent = "#2C7FB8",       # blue
    ink = "#1F2D3A",          # dark
    warning = "#E34A33",      # red-orange
    shade = "#D9D9D9",        # light gray

    # Comparison colors (K-only vs Dual)
    k_only = "#4682B4",       # steelblue (baseline)
    dual = "#E67E22",         # carrot orange (new solution)

    # Risk colors
    risk_high = "#DC143C",    # crimson
    risk_moderate = "#FF8C00", # darkorange
    risk_low = "#228B22"       # forestgreen
  )
}

# =============================================================================
# Distribution Helpers with Fallbacks
# =============================================================================

.dpprior_density_w1 <- function(x, a, b) {
  if (exists("density_w1", mode = "function")) {
    return(density_w1(x, a, b))
  }
  # Closed-form fallback: p(x|a,b) = a*b^a / ((1-x)[b - log(1-x)]^(a+1))
  x <- pmin(pmax(x, .Machine$double.eps), 1 - .Machine$double.eps)
  denom <- (1 - x) * (b - log(1 - x))^(a + 1)
  (a * b^a) / denom
}

.dpprior_prob_w1_exceeds <- function(t, a, b) {
  if (exists("prob_w1_exceeds", mode = "function")) {
    return(prob_w1_exceeds(t, a, b))
  }
  # Closed-form: P(w1 > t) = (b / (b - log(1-t)))^a
  t <- pmin(pmax(t, .Machine$double.eps), 1 - .Machine$double.eps)
  (b / (b - log(1 - t)))^a
}

.dpprior_quantile_w1 <- function(u, a, b) {
  if (exists("quantile_w1", mode = "function")) {
    return(quantile_w1(u, a, b))
  }
  u <- pmin(pmax(u, .Machine$double.eps), 1 - .Machine$double.eps)
  exponent <- b * (1 - (1 - u)^(-1 / a))
  1 - exp(exponent)
}

.dpprior_mean_w1 <- function(a, b) {
  if (exists("mean_w1", mode = "function")) {
    return(mean_w1(a, b))
  }
  f <- function(w) w * .dpprior_density_w1(w, a, b)
  tryCatch(
    stats::integrate(f, 0, 1, rel.tol = 1e-6)$value,
    error = function(e) NA_real_
  )
}

#' Hierarchical K PMF Computation with Fallbacks
#' @keywords internal
.dpprior_get_K_pmf <- function(J, a, b, M = NULL) {
  M <- .dpprior_coalesce(M, .dpprior_get_quad_nodes_default())

  # Method 1: Module 14 helper
  if (exists(".get_K_pmf_support", mode = "function")) {
    tryCatch({
      obj <- .get_K_pmf_support(J, a, b, M)
      return(list(k = obj$support, pmf = obj$pmf, method = "get_K_pmf_support"))
    }, error = function(e) NULL)
  }

  # Method 2: pmf_K_marginal + Stirling
  if (exists("pmf_K_marginal", mode = "function") &&
      exists("compute_log_stirling", mode = "function")) {
    tryCatch({
      logS <- compute_log_stirling(J)
      pmf0 <- pmf_K_marginal(J, a, b, logS, M)
      pmf <- as.numeric(pmf0[-1L])
      pmf <- pmf / sum(pmf)
      return(list(k = seq_len(J), pmf = pmf, method = "pmf_K_marginal"))
    }, error = function(e) NULL)
  }

  # Method 4: Negative Binomial approximation (last resort)
  warning("Using NegBin approximation for K_J PMF (exact functions unavailable)",
          call. = FALSE)
  HJ <- log(max(J, 2L))
  p_nb <- b / (b + HJ)
  k <- seq_len(J)
  pmf_nb <- stats::dnbinom(k - 1L, size = a, prob = p_nb)
  pmf_nb <- pmf_nb / sum(pmf_nb)

  list(k = k, pmf = pmf_nb, method = "negbin_approx")
}

# =============================================================================
# Summary Statistics Computation
# =============================================================================

#' Compute Summary Statistics
#' @keywords internal
.dpprior_compute_summary <- function(fit = NULL, a = NULL, b = NULL, J = NULL,
                                     ci_level = 0.95) {
  # Extract parameters
  if (!is.null(fit)) {
    abJ <- .dpprior_extract_abJ(fit)
    a <- abJ$a; b <- abJ$b; J <- abJ$J
  } else {
    if (is.null(a) || is.null(b)) {
      stop("Either 'fit' or both 'a' and 'b' must be provided", call. = FALSE)
    }
    J <- .dpprior_coalesce(J, 50L)
  }

  # Alpha statistics
  alpha_mean <- a / b
  alpha_cv <- 1 / sqrt(a)
  alpha_ci <- stats::qgamma(
    c((1 - ci_level) / 2, (1 + ci_level) / 2),
    shape = a, rate = b
  )

  # K statistics via hierarchical fallback
  M <- if (!is.null(fit)) .dpprior_coalesce(fit$M, NULL) else NULL
  Kp <- .dpprior_get_K_pmf(J, a, b, M = M)
  k <- Kp$k; pmf <- Kp$pmf
  cdf <- cumsum(pmf)

  K_mean <- sum(k * pmf)
  K_var <- max(0, sum((k^2) * pmf) - K_mean^2)
  K_mode <- k[which.max(pmf)]

  # w1 statistics
  w1_mean <- .dpprior_mean_w1(a, b)
  w1_median <- .dpprior_quantile_w1(0.5, a, b)
  w1_p50 <- .dpprior_prob_w1_exceeds(0.5, a, b)
  w1_p90 <- .dpprior_prob_w1_exceeds(0.9, a, b)

  # Dominance risk
  dominance_risk <- if (!is.null(fit)) {
    .dpprior_coalesce(
      fit$diagnostics$weights$dominance_risk,
      .dpprior_coalesce(fit$diagnostics$w1$dominance_risk, NA_character_)
    )
  } else {
    NA_character_
  }

  if (is.na(dominance_risk) || is.null(dominance_risk)) {
    dominance_risk <- if (w1_p50 >= 0.4) "HIGH" else if (w1_p50 >= 0.2) "MODERATE" else "LOW"
  } else {
    dominance_risk <- toupper(as.character(dominance_risk))
  }

  # Target values (from fit object if available)
  target_mu <- if (!is.null(fit)) .dpprior_coalesce(fit$target$mu_K, NA_real_) else NA_real_
  target_var <- if (!is.null(fit)) .dpprior_coalesce(fit$target$var_K, NA_real_) else NA_real_

  # Achieved values
  achieved_mu <- if (!is.null(fit)) {
    .dpprior_coalesce(.dpprior_coalesce(fit$fit$mu_K, fit$E_K), K_mean)
  } else {
    K_mean
  }
  achieved_var <- if (!is.null(fit)) {
    .dpprior_coalesce(.dpprior_coalesce(fit$fit$var_K, fit$Var_K), K_var)
  } else {
    K_var
  }

  list(
    a = a, b = b, J = J,
    alpha = list(mean = alpha_mean, cv = alpha_cv, ci = alpha_ci),
    K = list(k = k, pmf = pmf, cdf = cdf, mean = K_mean, var = K_var,
             mode = K_mode, pmf_method = Kp$method),
    w1 = list(mean = w1_mean, median = w1_median,
              p_gt_50 = w1_p50, p_gt_90 = w1_p90,
              dominance_risk = dominance_risk),
    target = list(mu_K = target_mu, var_K = target_var),
    achieved = list(mu_K = achieved_mu, var_K = achieved_var),
    method = if (!is.null(fit)) .dpprior_coalesce(fit$method, NA_character_) else NA_character_
  )
}

# =============================================================================
# Individual Plot Functions
# =============================================================================

#' Plot Prior Density of Alpha
#'
#' @param fit A DPprior_fit object, or NULL if a and b provided directly.
#' @param a Numeric; shape parameter (used if fit is NULL).
#' @param b Numeric; rate parameter (used if fit is NULL).
#' @param engine "ggplot2" (default) or "base".
#' @param base_size Base font size.
#' @param ci_level Credible interval level (default 0.95).
#' @param n_grid Number of grid points.
#' @param show If TRUE, print the plot.
#' @return A ggplot object or invisible(NULL) for base.
#'
#' @examples
#' # From fit object
#' fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8)
#' plot_alpha_prior(fit)
#'
#' # Direct parameter specification
#' plot_alpha_prior(a = 1.6, b = 1.2)
#'
#' @seealso \code{\link{DPprior_fit}} for fitting, \code{\link{plot.DPprior_fit}} for S3 plot method
#'
#' @family visualization
#'
#' @export
plot_alpha_prior <- function(fit = NULL, a = NULL, b = NULL,
                             engine = c("ggplot2", "base"),
                             base_size = 11,
                             ci_level = 0.95,
                             n_grid = 500,
                             show = TRUE) {
  engine <- match.arg(engine)

  # Extract parameters from fit or use direct values
  if (!is.null(fit)) {
    .dpprior_assert_fit(fit)
    abJ <- .dpprior_extract_abJ(fit)
    a <- abJ$a; b <- abJ$b
  } else {
    if (is.null(a) || is.null(b)) {
      stop("Either 'fit' or both 'a' and 'b' must be provided", call. = FALSE)
    }
  }

  # Compute statistics
  alpha_mean <- a / b
  alpha_cv <- 1 / sqrt(a)
  ci <- stats::qgamma(c((1 - ci_level) / 2, (1 + ci_level) / 2), shape = a, rate = b)

  # Grid
  x_max <- max(stats::qgamma(0.999, shape = a, rate = b), alpha_mean * 3)
  x <- seq(0.001, x_max, length.out = n_grid)
  dens <- stats::dgamma(x, shape = a, rate = b)
  df <- data.frame(x = x, density = dens)

  # Base R fallback
  if (engine == "base" || !.dpprior_has_ggplot2()) {
    .dpprior_base_plot_alpha(df, alpha_mean, ci, alpha_cv, a, b)
    return(invisible(NULL))
  }

  .dpprior_require_ggplot2()
  colors <- DPprior_colors()

  # Build subtitle with stats + Gamma parameters (two lines)
  subtitle_text <- sprintf(
    "E[alpha]=%.2f,  CV(alpha)=%.2f,  %d%% CI=[%.2f, %.2f]\nGamma(a=%.4f, b=%.3f)",
    alpha_mean, alpha_cv, round(100 * ci_level), ci[1], ci[2], a, b
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = density)) +
    ggplot2::geom_area(fill = colors$primary, alpha = 0.3) +
    ggplot2::geom_line(color = colors$primary, linewidth = 1) +
    ggplot2::annotate("rect", xmin = ci[1], xmax = ci[2],
                      ymin = -Inf, ymax = Inf, fill = colors$shade, alpha = 0.3) +
    ggplot2::geom_vline(xintercept = alpha_mean, linetype = "dashed",
                        color = colors$accent, linewidth = 0.9) +
    ggplot2::geom_vline(xintercept = ci, linetype = "dotted",
                        color = "gray50", linewidth = 0.5) +
    theme_DPprior(base_size = base_size) +
    ggplot2::labs(
      x = expression(Concentration~parameter~alpha),
      y = "Density",
      title = expression(bold("(A)")~Prior~on~alpha),
      subtitle = subtitle_text
    )

  if (isTRUE(show)) print(p)
  p
}


#' Plot Prior PMF of K_J
#'
#' @param fit A DPprior_fit object, or NULL if J, a, b provided directly.
#' @param J Integer; sample size (used if fit is NULL).
#' @param a Numeric; shape parameter (used if fit is NULL).
#' @param b Numeric; rate parameter (used if fit is NULL).
#' @param engine "ggplot2" (default) or "base".
#' @param base_size Base font size.
#' @param max_k Maximum k to display. If NULL, auto-determined by CDF >= 0.999.
#' @param show_cdf If TRUE, overlay CDF line. Default is FALSE.
#' @param show If TRUE, print the plot.
#' @return A ggplot object or invisible(NULL) for base.
#'
#' @examples
#' fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8)
#' plot_K_prior(fit)
#'
#' plot_K_prior(J = 50, a = 1.6, b = 1.2)
#'
#' @seealso \code{\link{DPprior_fit}} for fitting, \code{\link{plot.DPprior_fit}} for S3 plot method
#'
#' @family visualization
#'
#' @export
plot_K_prior <- function(fit = NULL, J = NULL, a = NULL, b = NULL,
                         engine = c("ggplot2", "base"),
                         base_size = 11,
                         max_k = NULL,
                         show_cdf = FALSE,
                         show = TRUE) {
  engine <- match.arg(engine)

  # Extract parameters
  if (!is.null(fit)) {
    .dpprior_assert_fit(fit)
    abJ <- .dpprior_extract_abJ(fit)
    J <- abJ$J; a <- abJ$a; b <- abJ$b
  } else {
    if (is.null(J) || is.null(a) || is.null(b)) {
      stop("Either 'fit' or all of 'J', 'a', 'b' must be provided", call. = FALSE)
    }
  }

  # Get summary with fallbacks
  summ <- .dpprior_compute_summary(fit = fit, a = a, b = b, J = J)
  k <- summ$K$k
  pmf <- summ$K$pmf
  cdf <- summ$K$cdf
  K_mean <- summ$K$mean
  K_var <- summ$K$var
  K_mode <- summ$K$mode
  target_mu <- summ$target$mu_K
  achieved_mu <- summ$achieved$mu_K

  # Auto max_k: smallest k with CDF >= 0.999
  if (is.null(max_k)) {
    cut_idx <- which(cdf >= 0.999)[1]
    if (is.na(cut_idx)) cut_idx <- length(k)
    max_k <- min(k[cut_idx], J, 100)
    max_k <- max(max_k, K_mode + 10, 15)
  }
  max_k <- min(max_k, J)

  # Truncate
  idx <- k <= max_k
  df <- data.frame(k = k[idx], pmf = pmf[idx], cdf = cdf[idx])

  # Base R fallback
  if (engine == "base" || !.dpprior_has_ggplot2()) {
    .dpprior_base_plot_K(df, target_mu, achieved_mu, K_mean, K_var, K_mode, a, b)
    return(invisible(NULL))
  }

  .dpprior_require_ggplot2()
  colors <- DPprior_colors()
  y_max <- max(df$pmf)

  # Build subtitle with stats + Gamma parameters (two lines)
  subtitle_text <- sprintf("E[K]=%.2f,  Var(K)=%.2f,  Mode=%d\nGamma(a=%.4f, b=%.3f)",
                           K_mean, K_var, K_mode, a, b)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = k, y = pmf)) +
    ggplot2::geom_col(fill = colors$primary, alpha = 0.7, width = 0.8) +
    theme_DPprior(base_size = base_size)

  # CDF overlay (scaled to PMF axis)
  if (show_cdf) {
    p <- p +
      ggplot2::geom_line(ggplot2::aes(y = cdf * y_max),
                         color = colors$ink, linewidth = 0.8) +
      ggplot2::geom_point(ggplot2::aes(y = cdf * y_max),
                          color = colors$ink, size = 1)
  }

  # Target mean line (if available)
  if (is.finite(target_mu)) {
    p <- p + ggplot2::geom_vline(xintercept = target_mu, linetype = "dashed",
                                 color = colors$warning, linewidth = 0.8)
  }

  # Achieved mean line
  p <- p + ggplot2::geom_vline(xintercept = achieved_mu, linetype = "solid",
                               color = colors$accent, linewidth = 0.8)

  p <- p + ggplot2::labs(
    x = expression(Number~of~clusters~K[J]),
    y = "Probability mass",
    title = expression(bold("(B)")~Prior~PMF~of~K[J]),
    subtitle = subtitle_text,
    caption = if (show_cdf) "Line: CDF (scaled to PMF axis)" else NULL
  )

  if (isTRUE(show)) print(p)
  p
}


#' Plot Prior Density of w1
#'
#' @param fit A DPprior_fit object, or NULL if a, b provided directly.
#' @param a Numeric; shape parameter (used if fit is NULL).
#' @param b Numeric; rate parameter (used if fit is NULL).
#' @param engine "ggplot2" (default) or "base".
#' @param base_size Base font size.
#' @param thresholds Dominance thresholds (default: c(0.5, 0.9)).
#' @param n_grid Number of grid points.
#' @param show If TRUE, print the plot.
#' @return A ggplot object or invisible(NULL) for base.
#'
#' @examples
#' fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8)
#' plot_w1_prior(fit)
#'
#' plot_w1_prior(a = 1.6, b = 1.2)
#'
#' @seealso \code{\link{DPprior_fit}} for fitting, \code{\link{plot.DPprior_fit}} for S3 plot method
#'
#' @family visualization
#'
#' @export
plot_w1_prior <- function(fit = NULL, a = NULL, b = NULL,
                          engine = c("ggplot2", "base"),
                          base_size = 11,
                          thresholds = c(0.5, 0.9),
                          n_grid = 500,
                          show = TRUE) {
  engine <- match.arg(engine)

  # Extract parameters
  if (!is.null(fit)) {
    .dpprior_assert_fit(fit)
    abJ <- .dpprior_extract_abJ(fit)
    a <- abJ$a; b <- abJ$b
  } else {
    if (is.null(a) || is.null(b)) {
      stop("Either 'fit' or both 'a' and 'b' must be provided", call. = FALSE)
    }
  }

  # Grid
  x <- seq(1e-6, 1 - 1e-6, length.out = n_grid)
  dens <- .dpprior_density_w1(x, a, b)
  dens[!is.finite(dens)] <- NA

  df <- data.frame(x = x, density = dens)
  df_shade <- df[df$x >= thresholds[1], , drop = FALSE]

  # Statistics
  mean_w <- .dpprior_mean_w1(a, b)
  median_w <- .dpprior_quantile_w1(0.5, a, b)
  p_gt_50 <- .dpprior_prob_w1_exceeds(thresholds[1], a, b)
  p_gt_90 <- .dpprior_prob_w1_exceeds(thresholds[2], a, b)

  # Dominance risk
  risk_level <- if (p_gt_50 >= 0.4) "HIGH" else if (p_gt_50 >= 0.2) "MODERATE" else "LOW"

  # Base R fallback
  if (engine == "base" || !.dpprior_has_ggplot2()) {
    .dpprior_base_plot_w1(df, thresholds, mean_w, median_w, p_gt_50, p_gt_90, risk_level, a, b)
    return(invisible(NULL))
  }

  .dpprior_require_ggplot2()
  colors <- DPprior_colors()
  risk_color <- switch(risk_level,
                       "HIGH" = colors$risk_high,
                       "MODERATE" = colors$risk_moderate,
                       "LOW" = colors$risk_low
  )

  y_max <- max(dens, na.rm = TRUE)
  if (!is.finite(y_max) || y_max > 20) y_max <- 20

  # Build subtitle with stats + Gamma parameters (two lines)
  subtitle_text <- sprintf(
    "E[w1]=%.3f,  Median=%.3f,  P(w1>0.5)=%.1f%%,  P(w1>0.9)=%.1f%%\nGamma(a=%.4f, b=%.3f)",
    mean_w, median_w, 100 * p_gt_50, 100 * p_gt_90, a, b
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = density)) +
    ggplot2::geom_area(data = df_shade, fill = colors$secondary, alpha = 0.3) +
    ggplot2::geom_line(color = colors$primary, linewidth = 1, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = thresholds, linetype = "dashed",
                        color = colors$warning, linewidth = 0.7) +
    theme_DPprior(base_size = base_size) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, y_max * 1.15)) +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.25)) +
    ggplot2::labs(
      x = expression(First~stick-breaking~weight~w[1]),
      y = "Density",
      title = expression(bold("(C)")~Prior~Density~of~w[1]),
      subtitle = subtitle_text
    )

  # Risk badge at top-right (no parse needed - plain text)
  p <- p + ggplot2::annotate("label",
                             x = 0.98, y = y_max * 1.05,
                             label = paste("Dominance Risk:", risk_level),
                             hjust = 1, vjust = 1, size = 3.5, fontface = "bold",
                             fill = risk_color, color = "white",
                             label.padding = ggplot2::unit(0.3, "lines"),
                             label.r = ggplot2::unit(0.15, "lines")
  )

  if (isTRUE(show)) print(p)
  p
}


#' Summary Table Plot
#' @keywords internal
.dpprior_summary_table_plot <- function(fit, base_size = 11, ci_level = 0.95) {
  .dpprior_require_ggplot2()

  summ <- .dpprior_compute_summary(fit, ci_level = ci_level)
  colors <- DPprior_colors()

  risk_color <- switch(summ$w1$dominance_risk,
                       "HIGH" = colors$risk_high,
                       "MODERATE" = colors$risk_moderate,
                       "LOW" = colors$risk_low
  )

  # Build table data (using Unicode for Greek letters and subscripts)
  tbl <- data.frame(
    Metric = c(
      "J", "Method", "Gamma prior",
      "E[alpha]", "CV(alpha)", sprintf("%d%% CI(alpha)", round(100*ci_level)),
      "Target E[K]", "Achieved E[K]", "Achieved Var(K)",
      "E[w1]", "P(w1 > 0.5)", "P(w1 > 0.9)"
    ),
    Value = c(
      sprintf("%d", summ$J),
      as.character(summ$method),
      sprintf("Gamma(%.4f, %.3f)", summ$a, summ$b),
      sprintf("%.3f", summ$alpha$mean),
      sprintf("%.3f", summ$alpha$cv),
      sprintf("[%.2f, %.2f]", summ$alpha$ci[1], summ$alpha$ci[2]),
      if (is.finite(summ$target$mu_K)) sprintf("%.2f", summ$target$mu_K) else "NA",
      sprintf("%.3f", summ$achieved$mu_K),
      sprintf("%.3f", summ$achieved$var_K),
      sprintf("%.3f", summ$w1$mean),
      sprintf("%.1f%%", 100 * summ$w1$p_gt_50),
      sprintf("%.1f%%", 100 * summ$w1$p_gt_90)
    ),
    stringsAsFactors = FALSE
  )

  tbl$row <- seq_len(nrow(tbl))

  p <- ggplot2::ggplot(tbl, ggplot2::aes(y = -row)) +
    ggplot2::geom_text(ggplot2::aes(x = 0, label = Metric),
                       hjust = 0, fontface = "bold", color = colors$ink, size = 3.2) +
    ggplot2::geom_text(ggplot2::aes(x = 1, label = Value),
                       hjust = 1, color = colors$ink, size = 3.2) +
    ggplot2::scale_x_continuous(limits = c(-0.05, 1.05)) +
    ggplot2::scale_y_continuous(limits = c(-nrow(tbl) - 2, 0)) +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    ) +
    ggplot2::labs(title = expression(bold("(D)")~Summary~Statistics))

  # Risk badge at bottom
  p <- p + ggplot2::annotate("label",
                             x = 0.5, y = -nrow(tbl) - 1,
                             label = paste("DOMINANCE RISK:", summ$w1$dominance_risk),
                             hjust = 0.5, vjust = 0.5, size = 4, fontface = "bold",
                             fill = risk_color, color = "white",
                             label.padding = ggplot2::unit(0.4, "lines")
  )

  p
}


#' Create Dashboard using gtable (no patchwork needed)
#' @keywords internal
.dpprior_dashboard_gtable <- function(p1, p2, p3, p4, title = NULL) {
  if (!requireNamespace("gtable", quietly = TRUE) ||
      !requireNamespace("grid", quietly = TRUE)) {
    return(NULL)
  }

  grobs <- matrix(list(
    ggplot2::ggplotGrob(p1), ggplot2::ggplotGrob(p2),
    ggplot2::ggplotGrob(p3), ggplot2::ggplotGrob(p4)
  ), nrow = 2, byrow = TRUE)

  g <- gtable::gtable_matrix(
    name = "dpprior_dashboard",
    grobs = grobs,
    widths = grid::unit(c(1, 1), "null"),
    heights = grid::unit(c(1, 1), "null")
  )

  # Add title if provided
  if (!is.null(title) && nzchar(title)) {
    title_grob <- grid::textGrob(
      title,
      gp = grid::gpar(fontsize = 14, fontface = "bold"),
      vjust = 0.5
    )
    g <- gtable::gtable_add_rows(g, heights = grid::unit(1.5, "lines"), pos = 0)
    g <- gtable::gtable_add_grob(g, title_grob, t = 1, l = 1, r = 2)
  }

  g
}


#' 4-Panel Prior Dashboard
#'
#' @param fit A DPprior_fit object.
#' @param engine "ggplot2" (default) or "base".
#' @param base_size Base font size.
#' @param ci_level Credible interval level.
#' @param title Optional overall title for the dashboard.
#' @param show If TRUE, draw the dashboard.
#' @return A gtable grob (for ggplot2) or invisible(NULL) for base.
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso \code{\link{DPprior_fit}} for fitting, \code{\link{plot.DPprior_fit}} for S3 plot method
#'
#' @family visualization
#'
#' @export
plot_prior_dashboard <- function(fit,
                                 engine = c("ggplot2", "base"),
                                 base_size = 11,
                                 ci_level = 0.95,
                                 title = NULL,
                                 show = TRUE) {
  engine <- match.arg(engine)
  .dpprior_assert_fit(fit)

  # Base R fallback
  if (engine == "base" || !.dpprior_has_ggplot2()) {
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)

    if (!is.null(title) && nzchar(title)) {
      graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
    } else {
      graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    }

    plot_alpha_prior(fit, engine = "base", ci_level = ci_level)
    plot_K_prior(fit, engine = "base")
    plot_w1_prior(fit, engine = "base")
    .dpprior_base_summary_panel(fit, ci_level = ci_level)

    if (!is.null(title) && nzchar(title)) {
      graphics::mtext(title, outer = TRUE, cex = 1.2, font = 2)
    }

    return(invisible(NULL))
  }

  .dpprior_require_ggplot2()

  # Create individual plots (don't show yet)
  p_alpha <- plot_alpha_prior(fit, engine = "ggplot2", base_size = base_size,
                              ci_level = ci_level, show = FALSE)
  p_K <- plot_K_prior(fit, engine = "ggplot2", base_size = base_size, show = FALSE)
  p_w1 <- plot_w1_prior(fit, engine = "ggplot2", base_size = base_size, show = FALSE)
  p_tbl <- .dpprior_summary_table_plot(fit, base_size = base_size, ci_level = ci_level)

  # Try gtable dashboard
  g <- .dpprior_dashboard_gtable(p_alpha, p_K, p_w1, p_tbl, title = title)

  if (is.null(g)) {
    # Fallback: return list
    if (isTRUE(show)) {
      print(p_alpha); print(p_K); print(p_w1); print(p_tbl)
    }
    return(invisible(list(alpha = p_alpha, K = p_K, w1 = p_w1, summary = p_tbl)))
  }

  if (isTRUE(show)) {
    grid::grid.newpage()
    grid::grid.draw(g)
  }

  g
}


# =============================================================================
# Dual-Anchor Visualization Functions
# =============================================================================

#' Plot Dual-Anchor Comparison Dashboard
#'
#' Creates a comparison dashboard showing K-only vs dual-anchor solutions.
#' Displays changes in alpha, K, and w1 distributions side-by-side.
#'
#' @param fit_dual A DPprior_fit object from DPprior_dual().
#' @param fit_K_only Optional K-only fit. If NULL, extracted from fit_dual$dual_anchor$init.
#' @param engine "ggplot2" (default) or "base".
#' @param base_size Base font size.
#' @param title Optional title for the dashboard.
#' @param show If TRUE, draw the plot.
#'
#' @return A gtable grob or list of ggplot objects.
#'
#' @examples
#' fit_K <- DPprior_fit(J = 50, mu_K = 5, var_K = 8)
#' fit_dual <- DPprior_dual(fit_K,
#'   w1_target = list(prob = list(threshold = 0.5, value = 0.3)),
#'   lambda = 0.5)
#' plot_dual_comparison(fit_dual)
#'
#' @seealso \code{\link{DPprior_fit}} for fitting, \code{\link{plot.DPprior_fit}} for S3 plot method
#'
#' @family visualization
#'
#' @export
plot_dual_comparison <- function(fit_dual,
                                 fit_K_only = NULL,
                                 engine = c("ggplot2", "base"),
                                 base_size = 10,
                                 title = NULL,
                                 show = TRUE) {
  engine <- match.arg(engine)
  .dpprior_assert_fit(fit_dual)

  if (!.dpprior_is_dual(fit_dual)) {
    stop("fit_dual must be a dual-anchor result from DPprior_dual()", call. = FALSE)
  }

  # Extract K-only parameters
  if (is.null(fit_K_only)) {
    a_K <- fit_dual$dual_anchor$init$a
    b_K <- fit_dual$dual_anchor$init$b
  } else {
    a_K <- fit_K_only$a
    b_K <- fit_K_only$b
  }

  a_dual <- fit_dual$a
  b_dual <- fit_dual$b
  J <- fit_dual$J
  lambda <- fit_dual$dual_anchor$lambda

  if (engine == "base" || !.dpprior_has_ggplot2()) {
    .dpprior_base_dual_comparison(a_K, b_K, a_dual, b_dual, J, lambda, title)
    return(invisible(NULL))
  }

  .dpprior_require_ggplot2()
  colors <- DPprior_colors()

  # ---- Panel 1: Alpha density comparison ----
  x_max <- max(
    stats::qgamma(0.999, shape = a_K, rate = b_K),
    stats::qgamma(0.999, shape = a_dual, rate = b_dual)
  ) * 1.1
  x <- seq(0.001, x_max, length.out = 400)

  df_alpha <- data.frame(
    x = rep(x, 2),
    density = c(
      stats::dgamma(x, shape = a_K, rate = b_K),
      stats::dgamma(x, shape = a_dual, rate = b_dual)
    ),
    Method = rep(c("K-only", "Dual-anchor"), each = length(x))
  )

  p_alpha <- ggplot2::ggplot(df_alpha, ggplot2::aes(x = x, y = density,
                                                    color = Method, linetype = Method)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_color_manual(values = c("K-only" = colors$k_only, "Dual-anchor" = colors$dual)) +
    ggplot2::scale_linetype_manual(values = c("K-only" = "solid", "Dual-anchor" = "solid")) +
    theme_DPprior(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
      x = expression(alpha), y = "Density",
      title = "Alpha Prior Comparison",
      subtitle = sprintf("K-only: Gamma(%.3f, %.3f)  |  Dual: Gamma(%.3f, %.3f)",
                         a_K, b_K, a_dual, b_dual)
    )

  # ---- Panel 2: K PMF comparison ----
  summ_K <- .dpprior_compute_summary(a = a_K, b = b_K, J = J)
  summ_dual <- .dpprior_compute_summary(a = a_dual, b = b_dual, J = J)

  max_k <- max(which(summ_K$K$cdf >= 0.999)[1], which(summ_dual$K$cdf >= 0.999)[1], na.rm = TRUE)
  max_k <- min(max_k, J, 50)

  df_K <- data.frame(
    k = rep(summ_K$K$k[1:max_k], 2),
    pmf = c(summ_K$K$pmf[1:max_k], summ_dual$K$pmf[1:max_k]),
    Method = rep(c("K-only", "Dual-anchor"), each = max_k)
  )

  p_K <- ggplot2::ggplot(df_K, ggplot2::aes(x = k, y = pmf, fill = Method)) +
    ggplot2::geom_col(position = "dodge", alpha = 0.8, width = 0.8) +
    ggplot2::scale_fill_manual(values = c("K-only" = colors$k_only, "Dual-anchor" = colors$dual)) +
    theme_DPprior(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
      x = expression(K[J]), y = "PMF",
      title = "K Distribution Comparison",
      subtitle = sprintf("K-only: E[K]=%.2f  |  Dual: E[K]=%.2f",
                         summ_K$K$mean, summ_dual$K$mean)
    )

  # ---- Panel 3: w1 density comparison ----
  x_w <- seq(1e-6, 1 - 1e-6, length.out = 400)
  dens_K <- .dpprior_density_w1(x_w, a_K, b_K)
  dens_dual <- .dpprior_density_w1(x_w, a_dual, b_dual)
  dens_K[!is.finite(dens_K)] <- NA
  dens_dual[!is.finite(dens_dual)] <- NA

  df_w1 <- data.frame(
    x = rep(x_w, 2),
    density = c(dens_K, dens_dual),
    Method = rep(c("K-only", "Dual-anchor"), each = length(x_w))
  )

  y_max_w <- min(max(c(dens_K, dens_dual), na.rm = TRUE), 20)

  p_gt_50_K <- .dpprior_prob_w1_exceeds(0.5, a_K, b_K)
  p_gt_50_dual <- .dpprior_prob_w1_exceeds(0.5, a_dual, b_dual)

  p_w1 <- ggplot2::ggplot(df_w1, ggplot2::aes(x = x, y = density,
                                              color = Method, linetype = Method)) +
    ggplot2::geom_line(linewidth = 1, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed", color = colors$warning, linewidth = 0.5) +
    ggplot2::scale_color_manual(values = c("K-only" = colors$k_only, "Dual-anchor" = colors$dual)) +
    ggplot2::scale_linetype_manual(values = c("K-only" = "solid", "Dual-anchor" = "solid")) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, y_max_w * 1.1)) +
    theme_DPprior(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
      x = expression(w[1]), y = "Density",
      title = "Weight Distribution Comparison",
      subtitle = sprintf("P(w1>0.5): K-only=%.1f%% | Dual=%.1f%% (reduction: %.1f%%)",
                         100 * p_gt_50_K, 100 * p_gt_50_dual,
                         100 * (p_gt_50_K - p_gt_50_dual) / p_gt_50_K)
    )

  # ---- Panel 4: Summary comparison table ----
  p_tbl <- .dpprior_dual_comparison_table(a_K, b_K, a_dual, b_dual, J,
                                          fit_dual$dual_anchor, base_size)

  # Combine into dashboard
  if (!is.null(title)) {
    title <- paste0(title, sprintf(" (lambda=%.2f)", lambda))
  } else {
    title <- sprintf("Dual-Anchor Comparison (lambda=%.2f)", lambda)
  }

  g <- .dpprior_dashboard_gtable(p_alpha, p_K, p_w1, p_tbl, title = title)

  if (is.null(g)) {
    if (isTRUE(show)) {
      print(p_alpha); print(p_K); print(p_w1); print(p_tbl)
    }
    return(invisible(list(alpha = p_alpha, K = p_K, w1 = p_w1, summary = p_tbl)))
  }

  if (isTRUE(show)) {
    grid::grid.newpage()
    grid::grid.draw(g)
  }

  g
}


#' Dual comparison summary table
#' @keywords internal
.dpprior_dual_comparison_table <- function(a_K, b_K, a_dual, b_dual, J,
                                           dual_info, base_size = 11) {
  .dpprior_require_ggplot2()
  colors <- DPprior_colors()

  # Compute metrics for both
  summ_K <- .dpprior_compute_summary(a = a_K, b = b_K, J = J)
  summ_dual <- .dpprior_compute_summary(a = a_dual, b = b_dual, J = J)

  # Target info
  w1_target_str <- "NA"
  if (!is.null(dual_info$w1_target$prob)) {
    w1_target_str <- sprintf("P(w1>%.1f)=%.0f%%",
                             dual_info$w1_target$prob$threshold,
                             100 * dual_info$w1_target$prob$value)
  } else if (!is.null(dual_info$w1_target$mean)) {
    w1_target_str <- sprintf("E[w1]=%.2f", dual_info$w1_target$mean)
  }

  tbl <- data.frame(
    Metric = c(
      "Gamma(a, b)",
      "E[alpha]",
      "E[K]",
      "Var(K)",
      "E[w1]",
      "P(w1 > 0.5)",
      "P(w1 > 0.9)",
      "---",
      "lambda",
      "w1 target",
      "loss_type"
    ),
    K_only = c(
      sprintf("(%.3f, %.3f)", a_K, b_K),
      sprintf("%.3f", summ_K$alpha$mean),
      sprintf("%.2f", summ_K$K$mean),
      sprintf("%.2f", summ_K$K$var),
      sprintf("%.3f", summ_K$w1$mean),
      sprintf("%.1f%%", 100 * summ_K$w1$p_gt_50),
      sprintf("%.1f%%", 100 * summ_K$w1$p_gt_90),
      "", "", "", ""
    ),
    Dual = c(
      sprintf("(%.3f, %.3f)", a_dual, b_dual),
      sprintf("%.3f", summ_dual$alpha$mean),
      sprintf("%.2f", summ_dual$K$mean),
      sprintf("%.2f", summ_dual$K$var),
      sprintf("%.3f", summ_dual$w1$mean),
      sprintf("%.1f%%", 100 * summ_dual$w1$p_gt_50),
      sprintf("%.1f%%", 100 * summ_dual$w1$p_gt_90),
      "",
      sprintf("%.2f", dual_info$lambda),
      w1_target_str,
      as.character(dual_info$loss_type)
    ),
    stringsAsFactors = FALSE
  )

  tbl$row <- seq_len(nrow(tbl))

  p <- ggplot2::ggplot(tbl, ggplot2::aes(y = -row)) +
    ggplot2::geom_text(ggplot2::aes(x = 0, label = Metric),
                       hjust = 0, fontface = "bold", color = colors$ink, size = 2.8) +
    ggplot2::geom_text(ggplot2::aes(x = 0.45, label = K_only),
                       hjust = 0.5, color = colors$k_only, size = 2.8, fontface = "bold") +
    ggplot2::geom_text(ggplot2::aes(x = 0.85, label = Dual),
                       hjust = 0.5, color = colors$dual, size = 2.8, fontface = "bold") +
    # Column headers
    ggplot2::annotate("text", x = 0.45, y = 0.5, label = "K-only",
                      fontface = "bold", color = colors$k_only, size = 3) +
    ggplot2::annotate("text", x = 0.85, y = 0.5, label = "Dual",
                      fontface = "bold", color = colors$dual, size = 3) +
    ggplot2::scale_x_continuous(limits = c(-0.05, 1.05)) +
    ggplot2::scale_y_continuous(limits = c(-nrow(tbl) - 1, 1)) +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = base_size),
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    ) +
    ggplot2::labs(title = "Comparison Summary")

  p
}


#' Plot Trade-off Curve
#'
#' Visualizes the Pareto trade-off between K_J fit and weight constraint
#' across different lambda values.
#'
#' @param tradeoff_data Data frame from compute_tradeoff_curve().
#' @param metric Which metric to plot on y-axis: "w1_prob_gt_50" (default),
#'   "E_w1", "K_loss", or "var_K".
#' @param target_value Optional target value to mark with horizontal line.
#' @param engine "ggplot2" (default) or "base".
#' @param base_size Base font size.
#' @param title Optional title.
#' @param show If TRUE, print the plot.
#'
#' @return A ggplot object or invisible(NULL).
#'
#' @examples
#' curve <- compute_tradeoff_curve(
#'   J = 50,
#'   K_target = list(mu_K = 5, var_K = 8),
#'   w1_target = list(prob = list(threshold = 0.5, value = 0.25)),
#'   lambda_seq = seq(0, 1, by = 0.1)
#' )
#' plot_tradeoff_curve(curve, target_value = 0.25)
#'
#' @seealso \code{\link{DPprior_fit}} for fitting, \code{\link{plot.DPprior_fit}} for S3 plot method
#'
#' @family visualization
#'
#' @export
plot_tradeoff_curve <- function(tradeoff_data,
                                metric = c("w1_prob_gt_50", "E_w1", "K_loss", "var_K"),
                                target_value = NULL,
                                engine = c("ggplot2", "base"),
                                base_size = 11,
                                title = NULL,
                                show = TRUE) {
  engine <- match.arg(engine)
  metric <- match.arg(metric)

  if (!is.data.frame(tradeoff_data) || !"lambda" %in% names(tradeoff_data)) {
    stop("tradeoff_data must be a data frame from compute_tradeoff_curve()", call. = FALSE)
  }

  if (!metric %in% names(tradeoff_data)) {
    stop(sprintf("metric '%s' not found in tradeoff_data", metric), call. = FALSE)
  }

  # Labels
  metric_labels <- list(
    w1_prob_gt_50 = "P(w1 > 0.5)",
    E_w1 = "E[w1]",
    K_loss = "K Loss (relative)",
    var_K = "Var(K)"
  )
  y_label <- metric_labels[[metric]]

  if (engine == "base" || !.dpprior_has_ggplot2()) {
    graphics::plot(tradeoff_data$lambda, tradeoff_data[[metric]],
                   type = "b", pch = 19, lwd = 2, col = "steelblue4",
                   xlab = "lambda (weight on K anchor)",
                   ylab = y_label,
                   main = if (is.null(title)) "Trade-off Curve" else title)
    if (!is.null(target_value)) {
      graphics::abline(h = target_value, lty = 2, col = "firebrick3")
    }
    graphics::grid()
    return(invisible(NULL))
  }

  .dpprior_require_ggplot2()
  colors <- DPprior_colors()

  p <- ggplot2::ggplot(tradeoff_data, ggplot2::aes(x = lambda, y = .data[[metric]])) +
    ggplot2::geom_line(color = colors$primary, linewidth = 1) +
    ggplot2::geom_point(color = colors$primary, size = 3) +
    theme_DPprior(base_size = base_size) +
    ggplot2::labs(
      x = "lambda (weight on K anchor)",
      y = y_label,
      title = if (is.null(title)) "Dual-Anchor Trade-off Curve" else title,
      subtitle = sprintf("lambda=0: w1-only | lambda=1: K-only | n=%d points", nrow(tradeoff_data))
    )

  if (!is.null(target_value)) {
    p <- p + ggplot2::geom_hline(yintercept = target_value,
                                 linetype = "dashed", color = colors$warning, linewidth = 0.8) +
      ggplot2::annotate("text", x = 0.02, y = target_value,
                        label = sprintf("Target: %.2f", target_value),
                        hjust = 0, vjust = -0.5, color = colors$warning, size = 3)
  }

  # Mark lambda=0, 0.5, 1
  key_lambdas <- c(0, 0.5, 1)
  key_data <- tradeoff_data[tradeoff_data$lambda %in% key_lambdas, ]
  if (nrow(key_data) > 0) {
    p <- p + ggplot2::geom_point(data = key_data, color = colors$accent, size = 4, shape = 18)
  }

  if (isTRUE(show)) print(p)
  p
}


#' Plot Trade-off Multi-Panel
#'
#' Creates a multi-panel view of the trade-off curve showing multiple metrics.
#'
#' @param tradeoff_data Data frame from compute_tradeoff_curve().
#' @param w1_target_prob Optional target probability for P(w1 > 0.5).
#' @param engine "ggplot2" (default) or "base".
#' @param base_size Base font size.
#' @param title Optional title.
#' @param show If TRUE, draw the plot.
#'
#' @return A gtable grob or list of ggplot objects.
#'
#' @seealso \code{\link{DPprior_fit}} for fitting, \code{\link{plot.DPprior_fit}} for S3 plot method
#'
#' @family visualization
#'
#' @export
plot_tradeoff_dashboard <- function(tradeoff_data,
                                    w1_target_prob = NULL,
                                    engine = c("ggplot2", "base"),
                                    base_size = 10,
                                    title = NULL,
                                    show = TRUE) {
  engine <- match.arg(engine)

  if (engine == "base" || !.dpprior_has_ggplot2()) {
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)
    graphics::par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

    plot_tradeoff_curve(tradeoff_data, "w1_prob_gt_50", w1_target_prob, "base")
    plot_tradeoff_curve(tradeoff_data, "E_w1", NULL, "base")
    plot_tradeoff_curve(tradeoff_data, "mu_K", NULL, "base")
    plot_tradeoff_curve(tradeoff_data, "var_K", NULL, "base")

    return(invisible(NULL))
  }

  .dpprior_require_ggplot2()
  colors <- DPprior_colors()

  # Panel 1: P(w1 > 0.5) vs lambda
  p1 <- plot_tradeoff_curve(tradeoff_data, "w1_prob_gt_50", w1_target_prob,
                            "ggplot2", base_size, "P(w1 > 0.5) vs lambda", FALSE)

  # Panel 2: E[w1] vs lambda
  p2 <- plot_tradeoff_curve(tradeoff_data, "E_w1", NULL,
                            "ggplot2", base_size, "E[w1] vs lambda", FALSE)

  # Panel 3: E[K] vs lambda (custom since mu_K may not be in metric list)
  p3 <- ggplot2::ggplot(tradeoff_data, ggplot2::aes(x = lambda, y = mu_K)) +
    ggplot2::geom_line(color = colors$primary, linewidth = 1) +
    ggplot2::geom_point(color = colors$primary, size = 3) +
    theme_DPprior(base_size = base_size) +
    ggplot2::labs(x = "lambda", y = "E[K]", title = "E[K] vs lambda")

  # Panel 4: Var(K) vs lambda
  p4 <- ggplot2::ggplot(tradeoff_data, ggplot2::aes(x = lambda, y = var_K)) +
    ggplot2::geom_line(color = colors$primary, linewidth = 1) +
    ggplot2::geom_point(color = colors$primary, size = 3) +
    theme_DPprior(base_size = base_size) +
    ggplot2::labs(x = "lambda", y = "Var(K)", title = "Var(K) vs lambda")

  # Combine
  if (is.null(title)) {
    title <- "Dual-Anchor Trade-off Analysis"
  }

  g <- .dpprior_dashboard_gtable(p1, p2, p3, p4, title = title)

  if (is.null(g)) {
    if (isTRUE(show)) {
      print(p1); print(p2); print(p3); print(p4)
    }
    return(invisible(list(p_w1 = p1, p_Ew1 = p2, p_EK = p3, p_VarK = p4)))
  }

  if (isTRUE(show)) {
    grid::grid.newpage()
    grid::grid.draw(g)
  }

  g
}


#' Plot Dual-Anchor Extended Dashboard
#'
#' Creates an extended dashboard for dual-anchor results with 6 panels:
#' (A) Alpha comparison, (B) K comparison, (C) w1 comparison,
#' (D) Summary table, (E) Loss decomposition, (F) Parameter trajectory.
#'
#' @param fit_dual A DPprior_fit object from DPprior_dual().
#' @param tradeoff_data Optional trade-off curve data for trajectory plot.
#' @param engine "ggplot2" (default) or "base".
#' @param base_size Base font size.
#' @param title Optional title.
#' @param show If TRUE, draw the plot.
#'
#' @return A gtable grob or invisible(NULL).
#'
#' @references
#' Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet Process Mixtures.
#' \emph{arXiv preprint} arXiv:2602.06301.
#'
#' @seealso \code{\link{DPprior_fit}} for fitting, \code{\link{plot.DPprior_fit}} for S3 plot method
#'
#' @family visualization
#'
#' @export
plot_dual_dashboard <- function(fit_dual,
                                tradeoff_data = NULL,
                                engine = c("ggplot2", "base"),
                                base_size = 10,
                                title = NULL,
                                show = TRUE) {
  engine <- match.arg(engine)
  .dpprior_assert_fit(fit_dual)

  if (!.dpprior_is_dual(fit_dual)) {
    # Fallback to standard dashboard
    return(plot_prior_dashboard(fit_dual, engine = engine, base_size = base_size,
                                title = title, show = show))
  }

  if (engine == "base" || !.dpprior_has_ggplot2()) {
    # Use comparison dashboard for base
    plot_dual_comparison(fit_dual, engine = "base")
    return(invisible(NULL))
  }

  .dpprior_require_ggplot2()
  colors <- DPprior_colors()

  # Extract parameters
  a_K <- fit_dual$dual_anchor$init$a
  b_K <- fit_dual$dual_anchor$init$b
  a_dual <- fit_dual$a
  b_dual <- fit_dual$b
  J <- fit_dual$J
  lambda <- fit_dual$dual_anchor$lambda

  # ---- Row 1: Distribution comparisons ----
  # Panel A: Alpha
  x_max <- max(stats::qgamma(0.999, a_K, b_K), stats::qgamma(0.999, a_dual, b_dual)) * 1.1
  x <- seq(0.001, x_max, length.out = 300)
  df_alpha <- data.frame(
    x = rep(x, 2),
    density = c(stats::dgamma(x, a_K, b_K), stats::dgamma(x, a_dual, b_dual)),
    Method = rep(c("K-only", "Dual"), each = length(x))
  )

  p_alpha <- ggplot2::ggplot(df_alpha, ggplot2::aes(x = x, y = density, color = Method)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_color_manual(values = c("K-only" = colors$k_only, "Dual" = colors$dual)) +
    theme_DPprior(base_size = base_size) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(x = expression(alpha), y = "Density", title = "(A) Alpha Prior")

  # Panel B: w1 comparison
  x_w <- seq(1e-6, 1 - 1e-6, length.out = 300)
  df_w1 <- data.frame(
    x = rep(x_w, 2),
    density = c(.dpprior_density_w1(x_w, a_K, b_K), .dpprior_density_w1(x_w, a_dual, b_dual)),
    Method = rep(c("K-only", "Dual"), each = length(x_w))
  )
  df_w1$density[!is.finite(df_w1$density)] <- NA
  y_max <- min(max(df_w1$density, na.rm = TRUE), 15)

  p_w1 <- ggplot2::ggplot(df_w1, ggplot2::aes(x = x, y = density, color = Method)) +
    ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed", color = colors$warning, linewidth = 0.5) +
    ggplot2::scale_color_manual(values = c("K-only" = colors$k_only, "Dual" = colors$dual)) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, y_max)) +
    theme_DPprior(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank()) +
    ggplot2::labs(x = expression(w[1]), y = "Density", title = "(B) w1 Distribution")

  # Panel C: Loss bar chart
  K_loss <- fit_dual$dual_anchor$K_loss
  w_loss <- fit_dual$dual_anchor$w_loss

  df_loss <- data.frame(
    Component = c("K Loss", "w1 Loss"),
    Value = c(K_loss, w_loss),
    Weight = c(lambda, 1 - lambda)
  )

  p_loss <- ggplot2::ggplot(df_loss, ggplot2::aes(x = Component, y = Value, fill = Component)) +
    ggplot2::geom_col(width = 0.6, alpha = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.4f\n(w=%.2f)", Value, Weight)),
                       vjust = -0.3, size = 2.5) +
    ggplot2::scale_fill_manual(values = c("K Loss" = colors$k_only, "w1 Loss" = colors$dual)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.25))) +
    theme_DPprior(base_size = base_size) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(x = "", y = "Loss Value", title = "(C) Loss Decomposition")

  # Panel D: Summary metrics
  p_summary <- .dpprior_dual_comparison_table(a_K, b_K, a_dual, b_dual, J,
                                              fit_dual$dual_anchor, base_size)
  p_summary <- p_summary + ggplot2::labs(title = "(D) Summary")

  # Combine 2x2
  if (is.null(title)) {
    title <- sprintf("Dual-Anchor Results (lambda=%.2f, %s)",
                     lambda, fit_dual$dual_anchor$loss_type)
  }

  g <- .dpprior_dashboard_gtable(p_alpha, p_w1, p_loss, p_summary, title = title)

  if (is.null(g)) {
    if (isTRUE(show)) {
      print(p_alpha); print(p_w1); print(p_loss); print(p_summary)
    }
    return(invisible(list(alpha = p_alpha, w1 = p_w1, loss = p_loss, summary = p_summary)))
  }

  if (isTRUE(show)) {
    grid::grid.newpage()
    grid::grid.draw(g)
  }

  g
}


#' Base R dual comparison fallback
#' @keywords internal
.dpprior_base_dual_comparison <- function(a_K, b_K, a_dual, b_dual, J, lambda, title) {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)

  if (!is.null(title) && nzchar(title)) {
    graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  } else {
    graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  }

  # Alpha comparison
  x_max <- max(stats::qgamma(0.999, a_K, b_K), stats::qgamma(0.999, a_dual, b_dual)) * 1.1
  x <- seq(0.001, x_max, length.out = 200)
  graphics::plot(x, stats::dgamma(x, a_K, b_K), type = "l", col = "#4682B4", lwd = 2,
                 xlab = expression(alpha), ylab = "Density", main = "(A) Alpha Prior")
  graphics::lines(x, stats::dgamma(x, a_dual, b_dual), col = "#E67E22", lwd = 2)
  graphics::legend("topright", c("K-only", "Dual"), lty = c(1, 1),
                   col = c("#4682B4", "#E67E22"), lwd = 2, cex = 0.8)

  # w1 comparison
  x_w <- seq(1e-4, 1 - 1e-4, length.out = 200)
  dens_K <- .dpprior_density_w1(x_w, a_K, b_K)
  dens_dual <- .dpprior_density_w1(x_w, a_dual, b_dual)
  y_max <- min(max(c(dens_K, dens_dual), na.rm = TRUE), 15)

  graphics::plot(x_w, dens_K, type = "l", col = "#4682B4", lwd = 2,
                 xlab = expression(w[1]), ylab = "Density", main = "(B) w1 Distribution",
                 ylim = c(0, y_max))
  graphics::lines(x_w, dens_dual, col = "#E67E22", lwd = 2)
  graphics::abline(v = 0.5, lty = 2, col = "firebrick3")

  # Summary text
  summ_K <- .dpprior_compute_summary(a = a_K, b = b_K, J = J)
  summ_dual <- .dpprior_compute_summary(a = a_dual, b = b_dual, J = J)

  graphics::plot.new()
  graphics::title(main = "(C) Comparison Summary")
  lines <- c(
    sprintf("Lambda = %.2f", lambda),
    "",
    "K-only:",
    sprintf("  Gamma(%.3f, %.3f)", a_K, b_K),
    sprintf("  P(w1>0.5) = %.1f%%", 100 * summ_K$w1$p_gt_50),
    "",
    "Dual-anchor:",
    sprintf("  Gamma(%.3f, %.3f)", a_dual, b_dual),
    sprintf("  P(w1>0.5) = %.1f%%", 100 * summ_dual$w1$p_gt_50)
  )
  y_pos <- seq(0.9, 0.1, length.out = length(lines))
  for (i in seq_along(lines)) {
    graphics::text(0.1, y_pos[i], lines[i], adj = 0, cex = 0.9)
  }

  # Reduction summary
  graphics::plot.new()
  graphics::title(main = "(D) Reduction")
  red_p50 <- (summ_K$w1$p_gt_50 - summ_dual$w1$p_gt_50) / summ_K$w1$p_gt_50
  red_mean <- (summ_K$w1$mean - summ_dual$w1$mean) / summ_K$w1$mean

  lines2 <- c(
    sprintf("P(w1>0.5) reduction: %.1f%%", 100 * red_p50),
    sprintf("E[w1] reduction: %.1f%%", 100 * red_mean),
    "",
    sprintf("E[K] change: %.2f -> %.2f", summ_K$K$mean, summ_dual$K$mean),
    sprintf("Var(K) change: %.2f -> %.2f", summ_K$K$var, summ_dual$K$var)
  )
  y_pos2 <- seq(0.8, 0.2, length.out = length(lines2))
  for (i in seq_along(lines2)) {
    graphics::text(0.1, y_pos2[i], lines2[i], adj = 0, cex = 0.9)
  }

  if (!is.null(title) && nzchar(title)) {
    graphics::mtext(title, outer = TRUE, cex = 1.2, font = 2)
  }
}


# =============================================================================
# Base R Fallbacks
# =============================================================================

.dpprior_base_plot_alpha <- function(df, alpha_mean, ci, alpha_cv, a, b) {
  graphics::plot(df$x, df$density, type = "l", lwd = 2,
                 xlab = expression(alpha), ylab = "Density",
                 main = expression(paste("(A) Prior on ", alpha)))
  graphics::polygon(df$x, df$density, col = "grey90", border = NA)
  graphics::lines(df$x, df$density, lwd = 2, col = "steelblue4")
  graphics::abline(v = alpha_mean, col = "darkred", lwd = 2, lty = 2)
  graphics::abline(v = ci, col = "gray50", lwd = 1, lty = 3)
  graphics::mtext(sprintf("Mean=%.3f, CV=%.2f, CI=[%.2f, %.2f], Gamma(%.4f, %.3f)",
                          alpha_mean, alpha_cv, ci[1], ci[2], a, b),
                  side = 3, line = 0.2, adj = 0, cex = 0.8)
}

.dpprior_base_plot_K <- function(df, target_mu, achieved_mu, K_mean, K_var, K_mode, a, b) {
  graphics::barplot(height = df$pmf, names.arg = df$k, border = NA,
                    col = "steelblue3", xlab = "k", ylab = "PMF",
                    main = expression(paste("(B) Prior PMF of ", K[J])))
  mids <- seq_along(df$k)
  cdf_scaled <- df$cdf * max(df$pmf)
  graphics::lines(mids, cdf_scaled, lwd = 2, col = "grey30")
  graphics::mtext(sprintf("E[K]=%.2f, Var=%.2f, Mode=%d, Gamma(%.4f, %.3f)",
                          K_mean, K_var, K_mode, a, b),
                  side = 3, line = 0.2, adj = 0, cex = 0.8)
}

.dpprior_base_plot_w1 <- function(df, thresholds, mean_w, median_w,
                                  p_gt_50, p_gt_90, risk_level, a, b) {
  graphics::plot(df$x, df$density, type = "l", lwd = 2, col = "steelblue4",
                 xlab = expression(w[1]), ylab = "Density",
                 main = expression(paste("(C) Prior Density of ", w[1])))
  shade <- df[df$x >= thresholds[1], , drop = FALSE]
  graphics::polygon(c(shade$x, rev(shade$x)),
                    c(shade$density, rep(0, nrow(shade))),
                    col = "grey85", border = NA)
  graphics::lines(df$x, df$density, lwd = 2, col = "steelblue4")
  graphics::abline(v = thresholds, col = "firebrick3", lwd = 1.5, lty = 2)
  graphics::mtext(sprintf("E[w1]=%.3f, P(w1>0.5)=%.1f%%, Risk: %s, Gamma(%.4f, %.3f)",
                          mean_w, 100*p_gt_50, risk_level, a, b),
                  side = 3, line = 0.2, adj = 0, cex = 0.8)
}

.dpprior_base_summary_panel <- function(fit, ci_level = 0.95) {
  summ <- .dpprior_compute_summary(fit, ci_level = ci_level)

  graphics::plot.new()
  graphics::title(main = "(D) Summary Statistics")

  lines <- c(
    sprintf("J = %d", summ$J),
    sprintf("Gamma(a=%.4f, b=%.3f)", summ$a, summ$b),
    sprintf("E[alpha] = %.3f, CV = %.2f", summ$alpha$mean, summ$alpha$cv),
    sprintf("E[K] = %.2f, Var(K) = %.2f", summ$achieved$mu_K, summ$achieved$var_K),
    sprintf("E[w1] = %.3f", summ$w1$mean),
    sprintf("P(w1>0.5) = %.1f%%", 100 * summ$w1$p_gt_50),
    sprintf("Dominance Risk: %s", summ$w1$dominance_risk)
  )

  y_pos <- seq(0.85, 0.15, length.out = length(lines))
  for (i in seq_along(lines)) {
    graphics::text(0.1, y_pos[i], lines[i], adj = 0, cex = 0.9)
  }
}


# =============================================================================
# Module Verification
# =============================================================================

#' Verify Visualization Module
#'
#' @param verbose Logical; if TRUE, print progress messages. Default is TRUE.
#'
#' @keywords internal
verify_visualization <- function(verbose = TRUE) {
  if (verbose) cat("Verifying visualization module (final)...\n")

  # Create test fit
  test_fit <- list(
    J = 50, a = 1.597, b = 1.222,
    E_K = 5.0, Var_K = 8.0,
    target = list(mu_K = 5.0, var_K = 8.0),
    method = "test"
  )
  class(test_fit) <- "DPprior_fit"

  tryCatch({
    if (verbose) cat("  Testing direct parameter API... ")
    p_direct <- plot_alpha_prior(a = 1.6, b = 1.2, show = FALSE)
    if (verbose) cat("OK\n")

    if (verbose) cat("  Testing .dpprior_compute_summary... ")
    summ <- .dpprior_compute_summary(test_fit)
    stopifnot(!is.null(summ$K$pmf))
    if (verbose) cat("OK (method:", summ$K$pmf_method, ")\n")

    if (.dpprior_has_ggplot2()) {
      if (verbose) cat("  Testing plot_alpha_prior(fit)... ")
      p1 <- plot_alpha_prior(test_fit, show = FALSE)
      stopifnot(inherits(p1, "ggplot"))
      if (verbose) cat("OK\n")

      if (verbose) cat("  Testing plot_K_prior(fit)... ")
      p2 <- plot_K_prior(test_fit, show = FALSE)
      stopifnot(inherits(p2, "ggplot"))
      if (verbose) cat("OK\n")

      if (verbose) cat("  Testing plot_w1_prior(fit)... ")
      p3 <- plot_w1_prior(test_fit, show = FALSE)
      stopifnot(inherits(p3, "ggplot"))
      if (verbose) cat("OK\n")

      if (verbose) cat("  Testing plot_prior_dashboard... ")
      dash <- plot_prior_dashboard(test_fit, show = FALSE)
      if (verbose) cat("OK\n")

      # Test dual-anchor visualization
      if (verbose) cat("  Testing dual-anchor detection... ")
      stopifnot(!.dpprior_is_dual(test_fit))  # Should be FALSE for K-only
      if (verbose) cat("OK (K-only correctly detected)\n")

      # Create mock dual-anchor fit
      test_dual <- test_fit
      test_dual$method <- "dual-anchor"
      test_dual$dual_anchor <- list(
        w1_target = list(prob = list(threshold = 0.5, value = 0.3)),
        lambda = 0.5,
        loss_type = "relative",
        w1_achieved = list(mean = 0.45, prob_gt_50 = 0.42, prob_gt_90 = 0.15),
        K_loss = 0.01,
        w_loss = 0.02,
        init = list(a = 1.597, b = 1.222)
      )

      if (verbose) cat("  Testing .dpprior_is_dual... ")
      stopifnot(.dpprior_is_dual(test_dual))
      if (verbose) cat("OK\n")

      if (verbose) cat("  Testing plot_dual_comparison... ")
      p_dual <- plot_dual_comparison(test_dual, show = FALSE)
      if (verbose) cat("OK\n")

      if (verbose) cat("  Testing plot_dual_dashboard... ")
      p_dual_dash <- plot_dual_dashboard(test_dual, show = FALSE)
      if (verbose) cat("OK\n")

      if (verbose) cat("  Testing S3 plot method (auto)... ")
      p_auto <- plot(test_dual, type = "auto", show = FALSE)
      if (verbose) cat("OK\n")

    } else {
      if (verbose) cat("  ggplot2 not available, skipping ggplot tests\n")
    }

    if (verbose) cat("All visualization tests passed!\n")
    TRUE
  }, error = function(e) {
    if (verbose) cat("FAILED:", e$message, "\n")
    FALSE
  })
}
