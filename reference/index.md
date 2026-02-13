# Package index

## Core Elicitation

Main user-facing functions for prior elicitation

- [`DPprior_fit()`](https://joonho112.github.io/DPprior/reference/DPprior_fit.md)
  : Fit a Gamma Hyperprior for DP Concentration Parameter
- [`DPprior_dual()`](https://joonho112.github.io/DPprior/reference/DPprior_dual.md)
  : Dual-Anchor Prior Calibration

## Approximation Algorithms

A1 closed-form and A2 Newton-based algorithms

- [`DPprior_a1()`](https://joonho112.github.io/DPprior/reference/DPprior_a1.md)
  : A1 Closed-Form Prior Elicitation
- [`DPprior_a2_newton()`](https://joonho112.github.io/DPprior/reference/DPprior_a2_newton.md)
  : A2-MN Exact-Moment Newton Solver
- [`DPprior_a2_kl()`](https://joonho112.github.io/DPprior/reference/DPprior_a2_kl.md)
  : A2-KL: KL Divergence Minimization for Prior Calibration

## Stirling Numbers

Computation of unsigned Stirling numbers of the first kind

- [`compute_log_stirling()`](https://joonho112.github.io/DPprior/reference/compute_log_stirling.md)
  : Compute Log Stirling Numbers (First Kind, Unsigned)

## K Distribution (Conditional)

Antoniak distribution P(K \| alpha)

- [`pmf_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/pmf_K_given_alpha.md)
  : PMF of K Given Alpha (Antoniak Distribution)
- [`log_pmf_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/log_pmf_K_given_alpha.md)
  : Log-PMF of K Given Alpha (Antoniak Distribution)
- [`cdf_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/cdf_K_given_alpha.md)
  : CDF of K Given Alpha
- [`mean_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/mean_K_given_alpha.md)
  : Conditional Mean of K_J Given Alpha
- [`var_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/var_K_given_alpha.md)
  : Conditional Variance of K_J Given Alpha
- [`mode_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/mode_K_given_alpha.md)
  : Mode of K Given Alpha
- [`quantile_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/quantile_K_given_alpha.md)
  : Quantile of K Given Alpha
- [`moments_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/moments_K_given_alpha.md)
  : Conditional Mean and Variance of K_J Given Alpha
- [`cv_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/cv_K_given_alpha.md)
  : Coefficient of Variation for K Given Alpha
- [`dispersion_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/dispersion_K_given_alpha.md)
  : Dispersion Index for K Given Alpha

## K Distribution (Marginal)

Marginal distribution of K under Gamma prior on alpha

- [`pmf_K_marginal()`](https://joonho112.github.io/DPprior/reference/pmf_K_marginal.md)
  : Marginal PMF of K_J under Gamma Hyperprior
- [`cdf_K_marginal()`](https://joonho112.github.io/DPprior/reference/cdf_K_marginal.md)
  : CDF of Marginal K Distribution
- [`mode_K_marginal()`](https://joonho112.github.io/DPprior/reference/mode_K_marginal.md)
  : Mode of Marginal K Distribution
- [`quantile_K_marginal()`](https://joonho112.github.io/DPprior/reference/quantile_K_marginal.md)
  : Quantile of Marginal K Distribution
- [`summary_K_marginal()`](https://joonho112.github.io/DPprior/reference/summary_K_marginal.md)
  : Summary Statistics for Marginal K Distribution
- [`K_moments()`](https://joonho112.github.io/DPprior/reference/K_moments.md)
  : Convenience Wrapper for Marginal Moments
- [`exact_K_moments()`](https://joonho112.github.io/DPprior/reference/exact_K_moments.md)
  : Exact Marginal Moments of K_J under Gamma Prior

## Weight Distribution (w1)

Distribution of the largest stick-breaking weight

- [`mean_w1()`](https://joonho112.github.io/DPprior/reference/mean_w1.md)
  : Mean of w₁
- [`var_w1()`](https://joonho112.github.io/DPprior/reference/var_w1.md)
  : Variance of w₁
- [`prob_w1_exceeds()`](https://joonho112.github.io/DPprior/reference/prob_w1_exceeds.md)
  : Survival Function of w₁
- [`quantile_w1()`](https://joonho112.github.io/DPprior/reference/quantile_w1.md)
  : Quantile Function of w₁
- [`density_w1()`](https://joonho112.github.io/DPprior/reference/density_w1.md)
  : Density of w₁
- [`cdf_w1()`](https://joonho112.github.io/DPprior/reference/cdf_w1.md)
  : CDF of First Stick-Breaking Weight w₁
- [`rw1()`](https://joonho112.github.io/DPprior/reference/rw1.md) :
  Random Generation from w₁ Distribution
- [`summary_w1()`](https://joonho112.github.io/DPprior/reference/summary_w1.md)
  : Summary Statistics for w₁ Distribution
- [`w1_grid()`](https://joonho112.github.io/DPprior/reference/w1_grid.md)
  : Compute w₁ Distribution on Grid

## Co-clustering Probability (rho)

Pairwise co-clustering probability functions

- [`mean_rho()`](https://joonho112.github.io/DPprior/reference/mean_rho.md)
  : Marginal Mean of rho
- [`var_rho()`](https://joonho112.github.io/DPprior/reference/var_rho.md)
  : Marginal Variance of rho
- [`mean_rho_given_alpha()`](https://joonho112.github.io/DPprior/reference/mean_rho_given_alpha.md)
  : Conditional Mean of rho Given Alpha
- [`var_rho_given_alpha()`](https://joonho112.github.io/DPprior/reference/var_rho_given_alpha.md)
  : Conditional Variance of rho Given Alpha
- [`cv_rho()`](https://joonho112.github.io/DPprior/reference/cv_rho.md)
  : Coefficient of Variation of rho
- [`rrho()`](https://joonho112.github.io/DPprior/reference/rrho.md) :
  Random Generation from rho Distribution
- [`summary_rho()`](https://joonho112.github.io/DPprior/reference/summary_rho.md)
  : Summary Statistics for rho Distribution
- [`rho_conditional_grid()`](https://joonho112.github.io/DPprior/reference/rho_conditional_grid.md)
  : Compute rho Conditional Moments on alpha Grid

## Diagnostics

Prior validation and diagnostic tools

- [`DPprior_diagnostics()`](https://joonho112.github.io/DPprior/reference/DPprior_diagnostics.md)
  : Comprehensive Prior Diagnostics
- [`DPprior_error_bounds()`](https://joonho112.github.io/DPprior/reference/DPprior_error_bounds.md)
  : Compute A1 Approximation Error Bounds
- [`check_dominance_risk()`](https://joonho112.github.io/DPprior/reference/check_dominance_risk.md)
  : Quick Dominance Risk Check
- [`dual_anchor_diagnostics()`](https://joonho112.github.io/DPprior/reference/dual_anchor_diagnostics.md)
  : Dual-Anchor Diagnostic Comparison
- [`compare_a1_a2()`](https://joonho112.github.io/DPprior/reference/compare_a1_a2.md)
  : Compare A1 vs A2 Accuracy

## Visualization

Plotting functions

- [`plot(`*`<DPprior_fit>`*`)`](https://joonho112.github.io/DPprior/reference/plot.DPprior_fit.md)
  : Plot Method for DPprior_fit Objects
- [`plot_alpha_prior()`](https://joonho112.github.io/DPprior/reference/plot_alpha_prior.md)
  : Plot Prior Density of Alpha
- [`plot_K_prior()`](https://joonho112.github.io/DPprior/reference/plot_K_prior.md)
  : Plot Prior PMF of K_J
- [`plot_w1_prior()`](https://joonho112.github.io/DPprior/reference/plot_w1_prior.md)
  : Plot Prior Density of w1
- [`plot_prior_dashboard()`](https://joonho112.github.io/DPprior/reference/plot_prior_dashboard.md)
  : 4-Panel Prior Dashboard
- [`plot_dual_comparison()`](https://joonho112.github.io/DPprior/reference/plot_dual_comparison.md)
  : Plot Dual-Anchor Comparison Dashboard
- [`plot_dual_dashboard()`](https://joonho112.github.io/DPprior/reference/plot_dual_dashboard.md)
  : Plot Dual-Anchor Extended Dashboard
- [`plot_tradeoff_curve()`](https://joonho112.github.io/DPprior/reference/plot_tradeoff_curve.md)
  : Plot Trade-off Curve
- [`plot_tradeoff_dashboard()`](https://joonho112.github.io/DPprior/reference/plot_tradeoff_dashboard.md)
  : Plot Trade-off Multi-Panel
- [`DPprior_colors()`](https://joonho112.github.io/DPprior/reference/DPprior_colors.md)
  : DPprior Color Palette
- [`theme_DPprior()`](https://joonho112.github.io/DPprior/reference/theme_DPprior.md)
  : Publication-Quality Theme for DPprior Plots

## S3 Methods

Print, summary, and coercion methods

- [`print(`*`<DPprior_fit>`*`)`](https://joonho112.github.io/DPprior/reference/print.DPprior_fit.md)
  : Print Method for DPprior_fit Objects
- [`summary(`*`<DPprior_fit>`*`)`](https://joonho112.github.io/DPprior/reference/summary.DPprior_fit.md)
  : Summary Method for DPprior_fit Objects
- [`print(`*`<summary.DPprior_fit>`*`)`](https://joonho112.github.io/DPprior/reference/print.summary.DPprior_fit.md)
  : Print Method for summary.DPprior_fit
- [`as.data.frame(`*`<DPprior_fit>`*`)`](https://joonho112.github.io/DPprior/reference/as.data.frame.DPprior_fit.md)
  : Coerce DPprior_fit to Data Frame
- [`print(`*`<DPprior_diagnostics>`*`)`](https://joonho112.github.io/DPprior/reference/print.DPprior_diagnostics.md)
  : Print Method for DPprior_diagnostics Objects
- [`summary(`*`<DPprior_diagnostics>`*`)`](https://joonho112.github.io/DPprior/reference/summary.DPprior_diagnostics.md)
  : Summary Method for DPprior_diagnostics Objects
- [`print(`*`<DPprior_error_bounds>`*`)`](https://joonho112.github.io/DPprior/reference/print.DPprior_error_bounds.md)
  : Print Method for DPprior_error_bounds Objects
- [`summary(`*`<DPprior_error_bounds>`*`)`](https://joonho112.github.io/DPprior/reference/summary.DPprior_error_bounds.md)
  : Summary Method for DPprior_error_bounds Objects
- [`print(`*`<w1_summary>`*`)`](https://joonho112.github.io/DPprior/reference/print.w1_summary.md)
  : Print Method for w1_summary Objects
- [`print(`*`<rho_summary>`*`)`](https://joonho112.github.io/DPprior/reference/print.rho_summary.md)
  : Print Method for rho_summary Objects

## Numerical Utilities

Conversion, quadrature, and numerical helper functions

- [`vif_to_variance()`](https://joonho112.github.io/DPprior/reference/vif_to_variance.md)
  : Convert Variance Inflation Factor to Variance
- [`confidence_to_vif()`](https://joonho112.github.io/DPprior/reference/confidence_to_vif.md)
  : Map a Qualitative Confidence Level to a Variance Inflation Factor
  (VIF)
- [`cv_alpha_to_variance()`](https://joonho112.github.io/DPprior/reference/cv_alpha_to_variance.md)
  : Convert CV(alpha) to Variance
- [`integrate_gamma()`](https://joonho112.github.io/DPprior/reference/integrate_gamma.md)
  : Integrate Function Against Gamma Distribution
- [`build_gamma_quadrature()`](https://joonho112.github.io/DPprior/reference/build_gamma_quadrature.md)
  : Build Quadrature for Gamma(a, b) Integration
- [`gauss_laguerre_nodes()`](https://joonho112.github.io/DPprior/reference/gauss_laguerre_nodes.md)
  : Gauss-Laguerre Quadrature Nodes and Weights
- [`logsumexp()`](https://joonho112.github.io/DPprior/reference/logsumexp.md)
  : Numerically Stable Log-Sum-Exp (Binary)
- [`logsumexp_vec()`](https://joonho112.github.io/DPprior/reference/logsumexp_vec.md)
  : Vectorized Log-Sum-Exp
- [`score_a()`](https://joonho112.github.io/DPprior/reference/score_a.md)
  : Score Function with Respect to Shape Parameter a
- [`score_b()`](https://joonho112.github.io/DPprior/reference/score_b.md)
  : Score Function with Respect to Rate Parameter b

## Computation

Lower-level computational functions

- [`compute_weight_diagnostics()`](https://joonho112.github.io/DPprior/reference/compute_weight_diagnostics.md)
  : Weight Distribution Diagnostics (w1)
- [`compute_error_landscape()`](https://joonho112.github.io/DPprior/reference/compute_error_landscape.md)
  : Compute A1 Error Landscape
- [`compute_tradeoff_curve()`](https://joonho112.github.io/DPprior/reference/compute_tradeoff_curve.md)
  : Compute Pareto Trade-off Curve
- [`compute_scaling_constant()`](https://joonho112.github.io/DPprior/reference/compute_scaling_constant.md)
  : Compute Scaling Constant for A1 Mapping
- [`a1_moment_error()`](https://joonho112.github.io/DPprior/reference/a1_moment_error.md)
  : A1 Approximation Moment Errors
- [`discretize_chisq()`](https://joonho112.github.io/DPprior/reference/discretize_chisq.md)
  : Discretize Chi-Square to K_J Support
- [`kl_divergence_K()`](https://joonho112.github.io/DPprior/reference/kl_divergence_K.md)
  : KL Divergence Between Target and Induced K_J PMFs
- [`kl_divergence_pmf()`](https://joonho112.github.io/DPprior/reference/kl_divergence_pmf.md)
  : KL Divergence Between Two PMFs
- [`moments_with_jacobian()`](https://joonho112.github.io/DPprior/reference/moments_with_jacobian.md)
  : Compute Marginal Moments and Jacobian Simultaneously

## Validation & Verification

Functions for verifying mathematical properties

- [`validate_stirling()`](https://joonho112.github.io/DPprior/reference/validate_stirling.md)
  : Validate Stirling Number Computation
- [`verify_jacobian()`](https://joonho112.github.io/DPprior/reference/verify_jacobian.md)
  : Verify Jacobian Against Finite Differences
- [`verify_underdispersion()`](https://joonho112.github.io/DPprior/reference/verify_underdispersion.md)
  : Verify Underdispersion Inequality
