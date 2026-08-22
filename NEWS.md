# DPprior 2.0.0

### Canonical Schema and Decision-Ready Calibration

* Introduced the versioned `dpprior.result/1`, `dpprior.target/1`, and
  `dpprior.weight-target/1` contracts. Current results separate scientific
  claims, numerical computation, independent verification, provenance, and
  non-authoritative compatibility views.
* Added `DPprior_target_K()` for validated bounded-discrete cluster-count
  targets and added the `target_K` route to `DPprior_fit()`.
* Migrated A1, A2-MN, A2-KL, diagnostics, sensitivity analysis, visualization,
  and S3 consumers to the canonical schema. A1 is now explicitly an
  approximate proxy and is not presented as finite-design verification.
* Added `DPprior_dual_hard()` for independently verified named weight
  inequalities. The certified `wmax_tail_upper` route is limited to upper
  safety constraints; direct W-max probability calibration remains deferred.
* Added `DPprior_dual_soft()` for current fixed-input-scale trade-offs with
  explicit lambda, component losses, candidate lineage, and independent
  verification.
* Retained `DPprior_dual()` as the historical equality-loss adapter throughout
  v2.x. It will not be removed before v3.0 and a migration review. Its output
  is always labelled legacy, approximate, and unverified, and it emits one
  typed lifecycle warning.
* Standardized failed or unapproved calibration outcomes as typed conditions
  that retain the complete canonical evidence in `condition$result`.
* Added `upgrade_DPprior_object()` for conservative, fail-closed handling of
  exactly recognized 1.1 objects. Migration does not promote incomplete legacy
  evidence into a verified current calibration.
* Updated documentation and examples to use canonical nested fields and the
  hard/soft dual-anchor APIs as the primary v2 workflow.

# DPprior 1.1.0

### Minor Release: Target PMFs, Feasibility Checks, And Release Hygiene

* Added a wrapper-level custom target PMF workflow in `DPprior_fit()`.
  `target_pmf` can now be supplied without `mu_K`; the wrapper dispatches to
  `method = "A2-KL"` and stores the normalized PMF target with its implied
  moments.
* Made `target_pmf` the canonical target when it is supplied. Explicit `mu_K`
  and `var_K` values are now treated as consistency checks and conflicting
  values produce clear errors instead of being silently overwritten.
* Tightened length `J + 1` target PMF validation. Inputs are validated before
  the structural `K = 0` entry is dropped, and positive mass at `K = 0` is
  rejected because `K_J` is supported on `{1, ..., J}`.
* Standardized moment feasibility across `DPprior_fit()`, `DPprior_a1()`,
  `DPprior_a2_newton()`, `construct_target_pmf()`, and
  `DPprior_a2_kl(method = "chisq")`. Moment workflows now consistently require
  `1 < mu_K < J`, since `mu_K = J` implies zero variance and is outside the
  positive-variance elicitation workflow.
* Enforced the fixed-mean variance bound
  `var_K <= (mu_K - 1) * (J - mu_K)` consistently and surfaced
  near-boundary confidence settings with clearer feasibility errors.
* Updated `DPprior_fit()` so non-converged low-variance A2-MN calibrations are
  reported as errors rather than returned as usable wrapper fits.
* Improved base graphics behavior so `show = FALSE` plotting paths return
  without opening a graphics device or leaving default `Rplots.pdf` artifacts.
* Fixed the source-build boundary for installed vignettes. Source tarballs now
  retain the generated `inst/doc` outputs required by `R CMD check`, while
  local development artifacts remain excluded from public release outputs.
* Updated package metadata, README citation version, and release-facing
  documentation for the 1.1.0 minor release.
* Added the exact Institute of Education Sciences support acknowledgement and
  disclaimer.

# DPprior 1.0.0

### Initial Public Release

This is the first public release of the DPprior package, providing tools for
principled prior elicitation on the concentration parameter alpha in Dirichlet
Process (DP) mixture models.

### Core Features

#### Elicitation Engine

* `DPprior_fit()`: Unified interface for K-based prior elicitation
  - Supports confidence levels ("low", "medium", "high") for easy specification
  - Direct variance specification for precise control
  - Automatic algorithm selection (A1 closed-form or A2 Newton refinement)

* `DPprior_a1()`: Closed-form approximation using Negative Binomial proxy
  - Near-instantaneous computation
  - Exploits asymptotic relationship K_J | alpha ~ Poisson(alpha log J)

* `DPprior_a2_newton()`: Exact moment matching via Newton iteration
  - Typically converges in 2-4 iterations
  - Guaranteed accuracy to specified tolerance

#### Dual-Anchor Framework

* `DPprior_dual()`: Joint control of cluster counts AND weight concentration
  - Addresses "unintended prior" problem (Vicentini & Jermyn, 2025)
  - Flexible weighting between K and w1 targets via lambda parameter
  - Supports probability, quantile, and moment constraints on w1

* `prob_w1_exceeds()`: Compute P(w1 > threshold) for dominance risk assessment
* `mean_w1()`, `var_w1()`: First and second moments of the first
  stick-breaking / size-biased weight
* `quantile_w1()`: Quantiles of w1 distribution

#### Exact Computation

* `compute_log_stirling()`: Stable computation of unsigned Stirling numbers
  - Log-scale for numerical stability with large J
  - Vectorized for efficiency

* `pmf_K_given_alpha()`: Exact Antoniak distribution P(K = k | alpha)
* `mean_K_given_alpha()`, `var_K_given_alpha()`: Conditional moments of K

#### Diagnostic Tools

* `DPprior_diagnostics()`: Comprehensive prior validation
  - Checks K, w1, rho, and alpha distributions
  - Identifies dominance risk (high P(w1 > 0.5))
  - Computes effective sample sizes

* `plot.DPprior_fit()`: Four-panel diagnostic dashboard
* `summary.DPprior_fit()`: Detailed numerical summary

#### Utility Functions

* `vif_to_variance()`: Convert variance inflation factor to Var(K)
* `confidence_to_vif()`: Map confidence levels to VIF values
* `integrate_gamma()`: High-precision Gauss-Laguerre integration
* `exact_K_moments()`: Marginal moments E[K] and Var(K) under Gamma prior

### Package Infrastructure

#### CRAN Compliance

* Passes R CMD check with 0 errors, 0 warnings, 0 notes
* Complete roxygen2-managed NAMESPACE with 77 exported functions and 11 S3 methods
* All examples use `\dontrun{}` or `\donttest{}` as appropriate
* No non-standard dependencies; base R + stats + graphics only in Imports

#### Focused Public API

* 77 carefully curated exports organized across 13 functional groups:
  core elicitation, approximation algorithms, Stirling numbers, K distribution
  (conditional and marginal), weight distribution, co-clustering probability,
  diagnostics, visualization, S3 methods, numerical utilities, computation,
  and validation/verification

#### Test Suite

* 2,084 unit tests via testthat 3.0
* Coverage spans all 20 source modules (R/00 through R/18 plus DPprior-package)
* Tests verify mathematical identities, numerical accuracy, edge cases,
  S3 method contracts, and visualization output

#### Documentation

* 49 `@family` cross-reference tags across 7 conceptual families
* 27 `@references` blocks citing Lee (2026) arXiv:2602.06301
* Terminology standardized to Design-Conditional Elicitation (DCE)
  and Two-Stage Moment Matching (TSMM)

#### Numerical Robustness

* 11 named constants in `R/00_constants.R` for reproducible thresholds
* `exp()` overflow protection via `.EXP_MAX` clamping in BFGS optimization
* Singular Jacobian fallback using correct gradient direction (J^T F)
* Division-by-near-zero guards in relative error computation
* PMF normalization guards for zero/non-finite sums

#### Comprehensive Vignettes

12 vignettes organized into two tracks:

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
- Weight Distributions: w1, rho, and dual-anchor framework
- API Reference: Complete function documentation

#### pkgdown Website

* Full pkgdown site at <https://joonho112.github.io/DPprior/>
* Reference index organized into 13 sections matching the public API
* Articles index with Applied and Methodological tracks
* Search functionality enabled
* Favicons and PWA manifest configured

### Methodological Foundation

This package implements the Design-Conditional Elicitation (DCE) methodology,
extending the original DORO approach (Dorazio, 2009) with:

1. **A1 closed-form approximation**: Instant initial estimates using the
   asymptotic Negative Binomial distribution of K_J under a Gamma prior on alpha
   (Zito et al., 2024)

2. **A2 Newton refinement**: Exact moment matching using numerically stable
   computation of Stirling numbers and Gauss-Laguerre quadrature

3. **Dual-anchor extension**: Joint control of K and w1 distributions,
   addressing the sample-size-independent concerns raised by
   Vicentini & Jermyn (2025)

### References

* Dorazio, R. M. (2009). On selecting a prior for the precision parameter of
  Dirichlet process mixture models. *Journal of Statistical Planning and
  Inference*, 139(10), 3384-3390.

* Lee, J. (2026). Design-conditional prior elicitation for Dirichlet process
  mixtures. *arXiv preprint* arXiv:2602.06301.

* Lee, J., Che, J., Rabe-Hesketh, S., Feller, A., & Miratrix, L. (2025).
  Improving the estimation of site-specific effects and their distribution
  in multisite trials. *Journal of Educational and Behavioral Statistics*,
  50(5), 731-764.

* Vicentini, C., & Jermyn, I. H. (2025). Prior selection for the precision
  parameter of Dirichlet process mixtures. *arXiv:2502.00864*.

* Zito, A., Rigon, T., & Dunson, D. B. (2024). Bayesian nonparametric modeling
  of latent partitions via Stirling-gamma priors. *arXiv:2306.02360*.

### Acknowledgments

This research was supported by the Institute of Education Sciences, U.S. Department of Education, through Grant R305D240078 to the University of Alabama.

The opinions expressed are those of the author and do not represent views of the Institute or the U.S. Department of Education.
