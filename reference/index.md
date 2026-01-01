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
- [`get_log_stirling()`](https://joonho112.github.io/DPprior/reference/get_log_stirling.md)
  : Get Single Log-Stirling Value with Bounds Checking
- [`get_stirling_row()`](https://joonho112.github.io/DPprior/reference/get_stirling_row.md)
  : Get Stirling Numbers for a Fixed J

## K Distribution (Conditional)

Antoniak distribution P(K \| α)

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
- [`summary_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/summary_K_given_alpha.md)
  : Summary of Conditional Moments
- [`cv_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/cv_K_given_alpha.md)
  : Coefficient of Variation for K Given Alpha
- [`dispersion_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/dispersion_K_given_alpha.md)
  : Dispersion Index for K Given Alpha

## K Distribution (Marginal)

Marginal distribution of K under Gamma prior on α

- [`pmf_K_marginal()`](https://joonho112.github.io/DPprior/reference/pmf_K_marginal.md)
  : Marginal PMF of K_J under Gamma Hyperprior
- [`log_pmf_K_marginal()`](https://joonho112.github.io/DPprior/reference/log_pmf_K_marginal.md)
  : Log Marginal PMF of K_J under Gamma Hyperprior
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

## Weight Distribution (w₁)

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

## Co-clustering Probability (ρ)

Pairwise co-clustering probability functions

- [`mean_rho()`](https://joonho112.github.io/DPprior/reference/mean_rho.md)
  : Marginal Mean of rho
- [`var_rho()`](https://joonho112.github.io/DPprior/reference/var_rho.md)
  : Marginal Variance of rho
- [`mean_rho_given_alpha()`](https://joonho112.github.io/DPprior/reference/mean_rho_given_alpha.md)
  : Conditional Mean of rho Given Alpha
- [`var_rho_given_alpha()`](https://joonho112.github.io/DPprior/reference/var_rho_given_alpha.md)
  : Conditional Variance of rho Given Alpha
- [`mean_rho_sq()`](https://joonho112.github.io/DPprior/reference/mean_rho_sq.md)
  : Marginal Second Moment of rho
- [`mean_rho_sq_given_alpha()`](https://joonho112.github.io/DPprior/reference/mean_rho_sq_given_alpha.md)
  : Conditional Second Moment of rho Given Alpha
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
- [`compare_diagnostics()`](https://joonho112.github.io/DPprior/reference/compare_diagnostics.md)
  : Compare Diagnostics Across Multiple Fits
- [`compare_a1_a2()`](https://joonho112.github.io/DPprior/reference/compare_a1_a2.md)
  : Compare A1 vs A2 Accuracy
- [`compare_rho_w1()`](https://joonho112.github.io/DPprior/reference/compare_rho_w1.md)
  : Compare rho and w1 Distributions
- [`compare_to_negbin()`](https://joonho112.github.io/DPprior/reference/compare_to_negbin.md)
  : Compare Exact Moments to NegBin Approximation

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

Print and summary methods

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

## Utilities

Helper and conversion functions

- [`vif_to_variance()`](https://joonho112.github.io/DPprior/reference/vif_to_variance.md)
  : Convert Variance Inflation Factor to Variance
- [`vif_to_variance_fit()`](https://joonho112.github.io/DPprior/reference/vif_to_variance_fit.md)
  : Convert VIF to Target Variance
- [`confidence_to_vif()`](https://joonho112.github.io/DPprior/reference/confidence_to_vif.md)
  : Map a Qualitative Confidence Level to a Variance Inflation Factor
  (VIF)
- [`confidence_to_vif_fit()`](https://joonho112.github.io/DPprior/reference/confidence_to_vif_fit.md)
  : Convert Confidence Level to Variance Inflation Factor
- [`cv_alpha_to_variance()`](https://joonho112.github.io/DPprior/reference/cv_alpha_to_variance.md)
  : Convert CV(alpha) to Variance
- [`integrate_gamma()`](https://joonho112.github.io/DPprior/reference/integrate_gamma.md)
  : Integrate Function Against Gamma Distribution
- [`build_gamma_quadrature()`](https://joonho112.github.io/DPprior/reference/build_gamma_quadrature.md)
  : Build Quadrature for Gamma(a, b) Integration
- [`gauss_laguerre_nodes()`](https://joonho112.github.io/DPprior/reference/gauss_laguerre_nodes.md)
  : Gauss-Laguerre Quadrature Nodes and Weights
- [`log_rising_factorial()`](https://joonho112.github.io/DPprior/reference/log_rising_factorial.md)
  : Log Rising Factorial (Pochhammer Symbol)
- [`logsumexp()`](https://joonho112.github.io/DPprior/reference/logsumexp.md)
  : Numerically Stable Log-Sum-Exp (Binary)
- [`logsumexp_vec()`](https://joonho112.github.io/DPprior/reference/logsumexp_vec.md)
  : Vectorized Log-Sum-Exp
- [`softmax()`](https://joonho112.github.io/DPprior/reference/softmax.md)
  : Numerically Stable Softmax

## Internal Computation

Lower-level computational functions

- [`compute_K_diagnostics()`](https://joonho112.github.io/DPprior/reference/compute_K_diagnostics.md)
  : K Distribution Diagnostics
- [`compute_weight_diagnostics()`](https://joonho112.github.io/DPprior/reference/compute_weight_diagnostics.md)
  : Weight Distribution Diagnostics (w1)
- [`compute_alpha_diagnostics()`](https://joonho112.github.io/DPprior/reference/compute_alpha_diagnostics.md)
  : Alpha Distribution Diagnostics
- [`compute_coclustering_diagnostics()`](https://joonho112.github.io/DPprior/reference/compute_coclustering_diagnostics.md)
  : Co-Clustering Diagnostics (rho)
- [`compute_error_landscape()`](https://joonho112.github.io/DPprior/reference/compute_error_landscape.md)
  : Compute A1 Error Landscape
- [`compute_tradeoff_curve()`](https://joonho112.github.io/DPprior/reference/compute_tradeoff_curve.md)
  : Compute Pareto Trade-off Curve
- [`compute_linearization_bound()`](https://joonho112.github.io/DPprior/reference/compute_linearization_bound.md)
  : Mean-Linearization Error Bound
- [`compute_poissonization_bound()`](https://joonho112.github.io/DPprior/reference/compute_poissonization_bound.md)
  : Poissonization Error Bound (Chen-Stein / Le Cam Bound)
- [`compute_total_tv_bound()`](https://joonho112.github.io/DPprior/reference/compute_total_tv_bound.md)
  : Total TV Error Bound (Conditional)
- [`compute_scaling_constant()`](https://joonho112.github.io/DPprior/reference/compute_scaling_constant.md)
  : Compute Scaling Constant for A1 Mapping
- [`compute_sum_p_squared()`](https://joonho112.github.io/DPprior/reference/compute_sum_p_squared.md)
  : Poissonization Error Bound (Raw Sum of Squared Probabilities)
- [`expected_tv_bound()`](https://joonho112.github.io/DPprior/reference/expected_tv_bound.md)
  : Expected TV Bound Under Gamma Prior

## Validation & Verification

Internal testing and validation functions

- [`verify_DPprior_fit()`](https://joonho112.github.io/DPprior/reference/verify_DPprior_fit.md)
  : Verify DPprior_fit Module
- [`verify_a1_mapping_all()`](https://joonho112.github.io/DPprior/reference/verify_a1_mapping_all.md)
  : Run All Module 10 Verification Tests
- [`verify_a2_all()`](https://joonho112.github.io/DPprior/reference/verify_a2_all.md)
  : Run All A2-MN Verification Tests
- [`verify_a2_kl()`](https://joonho112.github.io/DPprior/reference/verify_a2_kl.md)
  : Verify A2-KL Optimization
- [`verify_a2_kl_all()`](https://joonho112.github.io/DPprior/reference/verify_a2_kl_all.md)
  : Run All Module 12 Verification Tests
- [`verify_cdf_properties()`](https://joonho112.github.io/DPprior/reference/verify_cdf_properties.md)
  : Verify CDF Properties
- [`verify_derivative()`](https://joonho112.github.io/DPprior/reference/verify_derivative.md)
  : Verify Derivative via Finite Difference
- [`verify_diagnostics()`](https://joonho112.github.io/DPprior/reference/verify_diagnostics.md)
  : Verify Diagnostics Module
- [`verify_dual_anchor()`](https://joonho112.github.io/DPprior/reference/verify_dual_anchor.md)
  : Verify Dual-Anchor Module
- [`verify_jacobian_all()`](https://joonho112.github.io/DPprior/reference/verify_jacobian_all.md)
  : Run All Module 07 Verification Tests
- [`verify_kl_divergence()`](https://joonho112.github.io/DPprior/reference/verify_kl_divergence.md)
  : Verify KL Divergence Properties
- [`verify_marginal_moments()`](https://joonho112.github.io/DPprior/reference/verify_marginal_moments.md)
  : Verify Marginal Moments Properties
- [`verify_moments_marginal_all()`](https://joonho112.github.io/DPprior/reference/verify_moments_marginal_all.md)
  : Run All Module 05 Verification Tests
- [`verify_pmf_all()`](https://joonho112.github.io/DPprior/reference/verify_pmf_all.md)
  : Run All PMF Verifications
- [`verify_pmf_marginal_all()`](https://joonho112.github.io/DPprior/reference/verify_pmf_marginal_all.md)
  : Run All Module 06 Verification Tests
- [`verify_pmf_marginal_convergence()`](https://joonho112.github.io/DPprior/reference/verify_pmf_marginal_convergence.md)
  : Verify Quadrature Convergence
- [`verify_pmf_marginal_moments()`](https://joonho112.github.io/DPprior/reference/verify_pmf_marginal_moments.md)
  : Verify Moments Consistency
- [`verify_pmf_marginal_properties()`](https://joonho112.github.io/DPprior/reference/verify_pmf_marginal_properties.md)
  : Verify Marginal PMF Properties
- [`verify_pmf_moments()`](https://joonho112.github.io/DPprior/reference/verify_pmf_moments.md)
  : Verify PMF Consistency with Moments
- [`verify_pmf_normalization()`](https://joonho112.github.io/DPprior/reference/verify_pmf_normalization.md)
  : Verify PMF Normalization
- [`verify_quadrature()`](https://joonho112.github.io/DPprior/reference/verify_quadrature.md)
  : Verify Quadrature Accuracy Against Known Gamma Moments
- [`verify_s3_methods()`](https://joonho112.github.io/DPprior/reference/verify_s3_methods.md)
  : Verify S3 Methods Module
- [`verify_stirling_row_sum()`](https://joonho112.github.io/DPprior/reference/verify_stirling_row_sum.md)
  : Verify Row Sum Identity for Stirling Numbers
- [`verify_underdispersion()`](https://joonho112.github.io/DPprior/reference/verify_underdispersion.md)
  : Verify Underdispersion Inequality
- [`verify_visualization()`](https://joonho112.github.io/DPprior/reference/verify_visualization.md)
  : Verify Visualization Module
- [`verify_zero_probability()`](https://joonho112.github.io/DPprior/reference/verify_zero_probability.md)
  : Verify Zero Probability at K=0
- [`validate_moments_conditional()`](https://joonho112.github.io/DPprior/reference/validate_moments_conditional.md)
  : Validate Conditional Moments Computation
- [`validate_stirling()`](https://joonho112.github.io/DPprior/reference/validate_stirling.md)
  : Validate Stirling Number Computation

## Other Internal

Other internal helper functions

- [`a1_moment_error()`](https://joonho112.github.io/DPprior/reference/a1_moment_error.md)
  : A1 Approximation Moment Errors
- [`construct_target_pmf()`](https://joonho112.github.io/DPprior/reference/construct_target_pmf.md)
  : Construct Target PMF from User Specification
- [`discretize_chisq()`](https://joonho112.github.io/DPprior/reference/discretize_chisq.md)
  : Discretize Chi-Square to K_J Support
- [`kl_divergence_K()`](https://joonho112.github.io/DPprior/reference/kl_divergence_K.md)
  : KL Divergence Between Target and Induced K_J PMFs
- [`kl_divergence_pmf()`](https://joonho112.github.io/DPprior/reference/kl_divergence_pmf.md)
  : KL Divergence Between Two PMFs
- [`moments_with_jacobian()`](https://joonho112.github.io/DPprior/reference/moments_with_jacobian.md)
  : Compute Marginal Moments and Jacobian Simultaneously
- [`summary_pmf_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/summary_pmf_K_given_alpha.md)
  : Summary of Conditional PMF
