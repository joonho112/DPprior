# DPprior: Principled Prior Elicitation for Dirichlet Process Mixture Models

The DPprior package implements Design-Conditional Elicitation (DCE) for
calibrating the concentration parameter prior \\\alpha \sim
\text{Gamma}(a, b)\\ in Dirichlet Process (DP) mixture models. The core
engine is Two-Stage Moment Matching (TSMM), which translates intuitive
beliefs about the expected number of clusters \\E\[K_J\]\\ and its
variance \\\text{Var}(K_J)\\ into the Gamma hyperprior parameters \\(a,
b)\\.

## Details

The package provides three elicitation algorithms:

- A1
  ([`DPprior_a1`](https://joonho112.github.io/DPprior/reference/DPprior_a1.md)):

  Closed-form method using a shifted Negative Binomial approximation.
  Fast and accurate for large sample sizes.

- A2-MN
  ([`DPprior_a2_newton`](https://joonho112.github.io/DPprior/reference/DPprior_a2_newton.md)):

  Modified Newton method that matches exact marginal moments via
  Gauss-Laguerre quadrature and analytic Jacobians from score function
  identities.

- A2-KL
  ([`DPprior_a2_kl`](https://joonho112.github.io/DPprior/reference/DPprior_a2_kl.md)):

  KL-divergence minimization between the implied marginal PMF and a
  target distribution.

Additionally, the **Dual-Anchor** framework
([`DPprior_dual`](https://joonho112.github.io/DPprior/reference/DPprior_dual.md))
enables joint control of cluster count \\K_J\\ and largest weight
\\w_1\\ targets.

The unified interface
[`DPprior_fit`](https://joonho112.github.io/DPprior/reference/DPprior_fit.md)
dispatches to the appropriate algorithm based on the user's
specification.

## Key Function Groups

- Elicitation:

  [`DPprior_fit`](https://joonho112.github.io/DPprior/reference/DPprior_fit.md),
  [`DPprior_a1`](https://joonho112.github.io/DPprior/reference/DPprior_a1.md),
  [`DPprior_a2_newton`](https://joonho112.github.io/DPprior/reference/DPprior_a2_newton.md),
  [`DPprior_a2_kl`](https://joonho112.github.io/DPprior/reference/DPprior_a2_kl.md),
  [`DPprior_dual`](https://joonho112.github.io/DPprior/reference/DPprior_dual.md)

- Cluster Distribution:

  [`pmf_K_given_alpha`](https://joonho112.github.io/DPprior/reference/pmf_K_given_alpha.md),
  [`pmf_K_marginal`](https://joonho112.github.io/DPprior/reference/pmf_K_marginal.md),
  [`exact_K_moments`](https://joonho112.github.io/DPprior/reference/exact_K_moments.md)

- Weight Distribution:

  [`mean_w1`](https://joonho112.github.io/DPprior/reference/mean_w1.md),
  [`prob_w1_exceeds`](https://joonho112.github.io/DPprior/reference/prob_w1_exceeds.md),
  [`density_w1`](https://joonho112.github.io/DPprior/reference/density_w1.md),
  [`summary_w1`](https://joonho112.github.io/DPprior/reference/summary_w1.md)

- Diagnostics:

  [`DPprior_diagnostics`](https://joonho112.github.io/DPprior/reference/DPprior_diagnostics.md),
  [`DPprior_error_bounds`](https://joonho112.github.io/DPprior/reference/DPprior_error_bounds.md)

- Visualization:

  [`plot_prior_dashboard`](https://joonho112.github.io/DPprior/reference/plot_prior_dashboard.md),
  [`plot_dual_dashboard`](https://joonho112.github.io/DPprior/reference/plot_dual_dashboard.md)

## References

Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet
Process Mixtures. *arXiv preprint* arXiv:2602.06301.

Lee, J., Che, J., Rabe-Hesketh, S., Feller, A., and Miratrix, L. (2025).
Improving the estimation of site-specific effects and their distribution
in multisite trials. *Journal of Educational and Behavioral Statistics*,
50(5), 731–764.
[doi:10.3102/10769986241254286](https://doi.org/10.3102/10769986241254286)

## See also

Useful links:

- <https://joonho112.github.io/DPprior/>

- <https://github.com/joonho112/DPprior>

- Report bugs at <https://github.com/joonho112/DPprior/issues>

## Author

**Maintainer**: JoonHo Lee <jlee296@ua.edu>
([ORCID](https://orcid.org/0009-0006-4019-8703))
