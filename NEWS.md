# DPprior 1.0.0

### Initial CRAN Release

This is the first public release of the DPprior package, providing tools for 
principled prior elicitation on the concentration parameter α in Dirichlet 
Process (DP) mixture models.

### Core Features

#### Elicitation Engine

* `DPprior_fit()`: Unified interface for K-based prior elicitation
  - Supports confidence levels ("low", "medium", "high") for easy specification
  - Direct variance specification for precise control
  - Automatic algorithm selection (A1 closed-form or A2 Newton refinement)
  
* `DPprior_a1()`: Closed-form approximation using Negative Binomial proxy
  - Near-instantaneous computation
  - Exploits asymptotic relationship K_J | α ~ Poisson(α log J)
  
* `DPprior_a2_newton()`: Exact moment matching via Newton iteration
  - Typically converges in 2-4 iterations
  - Guaranteed accuracy to specified tolerance

#### Dual-Anchor Framework

* `DPprior_dual()`: Joint control of cluster counts AND weight concentration
  - Addresses "unintended prior" problem (Vicentini & Jermyn, 2025)
  - Flexible weighting between K and w₁ targets via λ parameter
  - Supports probability, quantile, and moment constraints on w₁

* `prob_w1_exceeds()`: Compute P(w₁ > threshold) for dominance risk assessment
* `mean_w1()`, `var_w1()`: First and second moments of largest weight
* `quantile_w1()`: Quantiles of w₁ distribution

#### Exact Computation

* `compute_log_stirling()`: Stable computation of unsigned Stirling numbers
  - Log-scale for numerical stability with large J
  - Vectorized for efficiency
  
* `pmf_K_given_alpha()`: Exact Antoniak distribution P(K = k | α)
* `mean_K_given_alpha()`, `var_K_given_alpha()`: Conditional moments of K

#### Diagnostic Tools

* `DPprior_diagnostics()`: Comprehensive prior validation
  - Checks K, w₁, ρ, and α distributions
  - Identifies dominance risk (high P(w₁ > 0.5))
  - Computes effective sample sizes
  
* `plot.DPprior_fit()`: Four-panel diagnostic dashboard
* `summary.DPprior_fit()`: Detailed numerical summary

#### Utility Functions

* `vif_to_variance()`: Convert variance inflation factor to Var(K)
* `confidence_to_vif()`: Map confidence levels to VIF values
* `integrate_gamma()`: High-precision Gauss-Laguerre integration
* `moments_K_marginal()`: Marginal moments E[K] and Var(K) under Gamma prior

### Documentation

Comprehensive vignettes organized into two tracks:

**Applied Researchers Track:**
- Introduction: Why prior elicitation matters
- Quick Start: Your first prior in 5 minutes
- Applied Guide: Complete elicitation workflow
- Dual-Anchor: Control counts AND weights
- Diagnostics: Verify prior behavior
- Case Studies: Multisite trials and meta-analysis

**Methodological Researchers Track:**
- Theory Overview: Mathematical foundations
- Stirling Numbers: Antoniak distribution details
- Approximations: A1 closed-form theory
- Newton Algorithm: A2 exact moment matching
- Weight Distributions: w₁, ρ, and dual-anchor framework
- API Reference: Complete function documentation

### Methodological Foundation

This package implements the DORO 2.0 methodology, extending the original DORO 
approach (Dorazio, 2009) with:

1. **A1 closed-form approximation**: Instant initial estimates using the 
   asymptotic Negative Binomial distribution of K_J under a Gamma prior on α 
   (Zito et al., 2024)

2. **A2 Newton refinement**: Exact moment matching using numerically stable 
   computation of Stirling numbers and Gauss-Laguerre quadrature

3. **Dual-anchor extension**: Joint control of K and w₁ distributions, 
   addressing the sample-size-independent concerns raised by 
   Vicentini & Jermyn (2025)

### References

* Dorazio, R. M. (2009). On selecting a prior for the precision parameter of 
  Dirichlet process mixture models. *Journal of Statistical Planning and 
  Inference*, 139(10), 3384–3390.

* Lee, J., Che, J., Rabe-Hesketh, S., Feller, A., & Miratrix, L. (2025). 
  Improving the estimation of site-specific effects and their distribution 
  in multisite trials. *Journal of Educational and Behavioral Statistics*, 
  50(5), 731–764.

* Vicentini, C., & Jermyn, I. H. (2025). Prior selection for the precision 
  parameter of Dirichlet process mixtures. *arXiv:2502.00864*.

* Zito, A., Rigon, T., & Dunson, D. B. (2024). Bayesian nonparametric modeling 
  of latent partitions via Stirling-gamma priors. *arXiv:2306.02360*.

### Acknowledgments

This project was supported by the Institute of Education Sciences, U.S. 
Department of Education, through Grant R305D240078 to University of Alabama.
