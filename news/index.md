# Changelog

## DPprior 1.0.0

#### Initial CRAN Release

This is the first public release of the DPprior package, providing tools
for principled prior elicitation on the concentration parameter α in
Dirichlet Process (DP) mixture models.

#### Core Features

##### Elicitation Engine

- [`DPprior_fit()`](https://joonho112.github.io/DPprior/reference/DPprior_fit.md):
  Unified interface for K-based prior elicitation
  - Supports confidence levels (“low”, “medium”, “high”) for easy
    specification
  - Direct variance specification for precise control
  - Automatic algorithm selection (A1 closed-form or A2 Newton
    refinement)
- [`DPprior_a1()`](https://joonho112.github.io/DPprior/reference/DPprior_a1.md):
  Closed-form approximation using Negative Binomial proxy
  - Near-instantaneous computation
  - Exploits asymptotic relationship K_J \| α ~ Poisson(α log J)
- [`DPprior_a2_newton()`](https://joonho112.github.io/DPprior/reference/DPprior_a2_newton.md):
  Exact moment matching via Newton iteration
  - Typically converges in 2-4 iterations
  - Guaranteed accuracy to specified tolerance

##### Dual-Anchor Framework

- [`DPprior_dual()`](https://joonho112.github.io/DPprior/reference/DPprior_dual.md):
  Joint control of cluster counts AND weight concentration
  - Addresses “unintended prior” problem (Vicentini & Jermyn, 2025)
  - Flexible weighting between K and w₁ targets via λ parameter
  - Supports probability, quantile, and moment constraints on w₁
- [`prob_w1_exceeds()`](https://joonho112.github.io/DPprior/reference/prob_w1_exceeds.md):
  Compute P(w₁ \> threshold) for dominance risk assessment
- [`mean_w1()`](https://joonho112.github.io/DPprior/reference/mean_w1.md),
  [`var_w1()`](https://joonho112.github.io/DPprior/reference/var_w1.md):
  First and second moments of largest weight
- [`quantile_w1()`](https://joonho112.github.io/DPprior/reference/quantile_w1.md):
  Quantiles of w₁ distribution

##### Exact Computation

- [`compute_log_stirling()`](https://joonho112.github.io/DPprior/reference/compute_log_stirling.md):
  Stable computation of unsigned Stirling numbers
  - Log-scale for numerical stability with large J
  - Vectorized for efficiency
- [`pmf_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/pmf_K_given_alpha.md):
  Exact Antoniak distribution P(K = k \| α)
- [`mean_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/mean_K_given_alpha.md),
  [`var_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/var_K_given_alpha.md):
  Conditional moments of K

##### Diagnostic Tools

- [`DPprior_diagnostics()`](https://joonho112.github.io/DPprior/reference/DPprior_diagnostics.md):
  Comprehensive prior validation
  - Checks K, w₁, ρ, and α distributions
  - Identifies dominance risk (high P(w₁ \> 0.5))
  - Computes effective sample sizes
- [`plot.DPprior_fit()`](https://joonho112.github.io/DPprior/reference/plot.DPprior_fit.md):
  Four-panel diagnostic dashboard
- [`summary.DPprior_fit()`](https://joonho112.github.io/DPprior/reference/summary.DPprior_fit.md):
  Detailed numerical summary

##### Utility Functions

- [`vif_to_variance()`](https://joonho112.github.io/DPprior/reference/vif_to_variance.md):
  Convert variance inflation factor to Var(K)
- [`confidence_to_vif()`](https://joonho112.github.io/DPprior/reference/confidence_to_vif.md):
  Map confidence levels to VIF values
- [`integrate_gamma()`](https://joonho112.github.io/DPprior/reference/integrate_gamma.md):
  High-precision Gauss-Laguerre integration
- [`exact_K_moments()`](https://joonho112.github.io/DPprior/reference/exact_K_moments.md):
  Marginal moments E\[K\] and Var(K) under Gamma prior

#### Package Infrastructure

##### CRAN Compliance

- Passes R CMD check with 0 errors, 0 warnings, 0 notes
- Complete roxygen2-managed NAMESPACE with 77 exported functions and 11
  S3 methods
- All examples use `\dontrun{}` or `\donttest{}` as appropriate
- No non-standard dependencies; base R + stats + graphics only in
  Imports

##### Focused Public API

- 77 carefully curated exports organized across 13 functional groups:
  core elicitation, approximation algorithms, Stirling numbers, K
  distribution (conditional and marginal), weight distribution,
  co-clustering probability, diagnostics, visualization, S3 methods,
  numerical utilities, computation, and validation/verification

##### Test Suite

- 2,084 unit tests via testthat 3.0
- Coverage spans all 20 source modules (R/00 through R/18 plus
  DPprior-package)
- Tests verify mathematical identities, numerical accuracy, edge cases,
  S3 method contracts, and visualization output

##### Documentation

- 49 `@family` cross-reference tags across 7 conceptual families
- 27 `@references` blocks citing Lee (2026) arXiv:2602.06301
- Terminology standardized to Design-Conditional Elicitation (DCE) and
  Two-Stage Moment Matching (TSMM)

##### Numerical Robustness

- 11 named constants in `R/00_constants.R` for reproducible thresholds
- [`exp()`](https://rdrr.io/r/base/Log.html) overflow protection via
  `.EXP_MAX` clamping in BFGS optimization
- Singular Jacobian fallback using correct gradient direction (J^T F)
- Division-by-near-zero guards in relative error computation
- PMF normalization guards for zero/non-finite sums

##### Comprehensive Vignettes

12 vignettes organized into two tracks:

**Applied Researchers Track:** - Introduction: Why prior elicitation
matters - Quick Start: Your first prior in 5 minutes - Applied Guide:
Complete elicitation workflow - Dual-Anchor: Control counts AND
weights - Diagnostics: Verify prior behavior - Case Studies: Multisite
trials and meta-analysis

**Methodological Researchers Track:** - Theory Overview: Mathematical
foundations - Stirling Numbers: Antoniak distribution details -
Approximations: A1 closed-form theory - Newton Algorithm: A2 exact
moment matching - Weight Distributions: w₁, ρ, and dual-anchor
framework - API Reference: Complete function documentation

##### pkgdown Website

- Full pkgdown site at <https://joonho112.github.io/DPprior/>
- Reference index organized into 13 sections matching the public API
- Articles index with Applied and Methodological tracks
- Search functionality enabled
- Favicons and PWA manifest configured

#### Methodological Foundation

This package implements the Design-Conditional Elicitation (DCE)
methodology, extending the original DORO approach (Dorazio, 2009) with:

1.  **A1 closed-form approximation**: Instant initial estimates using
    the asymptotic Negative Binomial distribution of K_J under a Gamma
    prior on α (Zito et al., 2024)

2.  **A2 Newton refinement**: Exact moment matching using numerically
    stable computation of Stirling numbers and Gauss-Laguerre quadrature

3.  **Dual-anchor extension**: Joint control of K and w₁ distributions,
    addressing the sample-size-independent concerns raised by Vicentini
    & Jermyn (2025)

#### References

- Dorazio, R. M. (2009). On selecting a prior for the precision
  parameter of Dirichlet process mixture models. *Journal of Statistical
  Planning and Inference*, 139(10), 3384–3390.

- Lee, J. (2026). Design-conditional prior elicitation for Dirichlet
  process mixtures. *arXiv preprint* arXiv:2602.06301.

- Lee, J., Che, J., Rabe-Hesketh, S., Feller, A., & Miratrix, L. (2025).
  Improving the estimation of site-specific effects and their
  distribution in multisite trials. *Journal of Educational and
  Behavioral Statistics*, 50(5), 731–764.

- Vicentini, C., & Jermyn, I. H. (2025). Prior selection for the
  precision parameter of Dirichlet process mixtures. *arXiv:2502.00864*.

- Zito, A., Rigon, T., & Dunson, D. B. (2024). Bayesian nonparametric
  modeling of latent partitions via Stirling-gamma priors.
  *arXiv:2306.02360*.

#### Acknowledgments

This project was supported by the Institute of Education Sciences, U.S.
Department of Education, through Grant R305D240078 to University of
Alabama.
