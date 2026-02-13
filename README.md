# DPprior <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**Principled Prior Elicitation for Dirichlet Process Mixture Models**

The DPprior package provides tools for eliciting Gamma hyperpriors on the 
concentration parameter α in Dirichlet Process (DP) mixture models. Rather than 
requiring researchers to reason about the abstract parameter α, DPprior allows 
specification through intuitive quantities:

- **Expected cluster counts**: "How many distinct groups do you anticipate?"
- **Weight concentration**: "How evenly distributed do you expect observations across clusters?"

These natural questions are translated into principled Gamma(a, b) hyperpriors 
using computationally efficient algorithms backed by exact moment matching.

## Installation

```r
# Install from CRAN (when available)
install.packages("DPprior")

# Or install the development version from GitHub
# install.packages("devtools")
devtools::install_github("joonho112/DPprior")
```

## Quick Start

```r
library(DPprior)

# Scenario: 50-site multisite trial, expecting ~5 distinct effect patterns
fit <- DPprior_fit(
  J = 50,                # Number of sites

  mu_K = 5,              # Expected clusters
  confidence = "medium"  # Moderate uncertainty
)

# View the elicited prior
print(fit)
#> DPprior Elicitation Results
#> ──────────────────────────────────────────────────────────────
#> Prior: α ~ Gamma(1.892, 1.201)
#> Target: E[K] = 5.00, Var(K) = 12.50
#> Achieved: E[K] = 5.00, Var(K) = 12.50
#> Method: A2-MN (converged in 3 iterations)

# Visualize the prior
plot(fit)
```

## Key Features

### 1. Intuitive Elicitation

Specify priors through expected cluster counts and uncertainty levels:

```r
# Using confidence levels (recommended for most users)
fit <- DPprior_fit(J = 50, mu_K = 5, confidence = "medium")

# Or specify variance directly
fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 10)
```

### 2. Dual-Anchor Control

Go beyond cluster counts to control weight behavior, addressing the "unintended 
prior" problem ([Vicentini & Jermyn, 2025](https://doi.org/10.48550/arXiv.2502.00864)):

```r
# First, fit K-based prior
fit_K <- DPprior_fit(J = 50, mu_K = 5, var_K = 8)

# Check if largest cluster might dominate
prob_w1_exceeds(0.5, fit_K$a, fit_K$b)
#> [1] 0.52  # 52% chance one cluster has >50% of observations

# Apply dual-anchor constraint
w1_target <- list(prob = list(threshold = 0.5, value = 0.30))
fit_dual <- DPprior_dual(fit_K, w1_target, lambda = 0.5)

# Verify improvement
prob_w1_exceeds(0.5, fit_dual$a, fit_dual$b)
#> [1] 0.31  # Now only 31%
```

### 3. Comprehensive Diagnostics

Verify your prior behaves as intended across all relevant dimensions:
- K distribution (cluster counts)
- w₁ distribution (largest cluster weight)
- ρ distribution (co-clustering probability)
- α distribution (concentration parameter)

```r
fit <- DPprior_fit(J = 50, mu_K = 5, check_diagnostics = TRUE)
plot(fit)  # Four-panel diagnostic dashboard
summary(fit)  # Detailed numerical diagnostics
```

### 4. Fast Computation

The package implements the Design-Conditional Elicitation (DCE) methodology via Two-Stage Moment Matching (TSMM):

- **A1 (Closed-form)**: Instant initial estimates using Negative Binomial approximation
- **A2 (Newton refinement)**: Exact moment matching in 2-4 iterations

```r
# A1 only (fastest, approximate)
fit_fast <- DPprior_fit(J = 50, mu_K = 5, var_K = 10, method = "A1")

# A2 with Newton refinement (default, exact)
fit_exact <- DPprior_fit(J = 50, mu_K = 5, var_K = 10, method = "A2-MN")
```

## When to Use DPprior

DPprior is particularly valuable for:

- **Multisite randomized trials** with moderate numbers of sites (J = 20–200)
- **Meta-analysis** with flexible heterogeneity modeling
- **Bayesian nonparametric density estimation** in small-to-moderate samples
- **Low-information settings** where the prior on α substantially influences posterior inference ([Lee et al., 2025](https://doi.org/10.3102/10769986241254286))

## Vignettes

The package includes comprehensive documentation:

### Applied Researchers Track

| Vignette | Description |
|----------|-------------|
| [Introduction](https://joonho112.github.io/DPprior/articles/introduction.html) | Why prior elicitation matters |
| [Quick Start](https://joonho112.github.io/DPprior/articles/quick-start.html) | Your first prior in 5 minutes |
| [Applied Guide](https://joonho112.github.io/DPprior/articles/applied-guide.html) | Complete elicitation workflow |
| [Dual-Anchor](https://joonho112.github.io/DPprior/articles/dual-anchor.html) | Control cluster counts AND weights |
| [Diagnostics](https://joonho112.github.io/DPprior/articles/diagnostics.html) | Verify your prior behaves as intended |
| [Case Studies](https://joonho112.github.io/DPprior/articles/case-studies.html) | Real-world applications |

### Methodological Researchers Track

| Vignette | Description |
|----------|-------------|
| [Theory Overview](https://joonho112.github.io/DPprior/articles/theory-overview.html) | Mathematical foundations |
| [Stirling Numbers](https://joonho112.github.io/DPprior/articles/theory-stirling.html) | Antoniak distribution details |
| [Approximations](https://joonho112.github.io/DPprior/articles/theory-approximations.html) | A1 closed-form theory |
| [Newton Algorithm](https://joonho112.github.io/DPprior/articles/theory-newton.html) | A2 exact moment matching |
| [Weight Distributions](https://joonho112.github.io/DPprior/articles/theory-weights.html) | w₁, ρ, and dual-anchor |
| [API Reference](https://joonho112.github.io/DPprior/articles/api-reference.html) | Complete function documentation |

## Citation

If you use DPprior in your research, please cite:

```bibtex
@Manual{DPprior2026,
  title = {{DPprior}: Principled Prior Elicitation for {Dirichlet} Process Mixture Models},
  author = {JoonHo Lee},
  year = {2026},
  note = {R package version 1.0.0},
  url = {https://github.com/joonho112/DPprior},
}

@Article{Lee2025multisite,
  title = {Improving the Estimation of Site-Specific Effects and Their Distribution in Multisite Trials},
  author = {JoonHo Lee and Jonathan Che and Sophia Rabe-Hesketh and Avi Feller and Luke Miratrix},
  journal = {Journal of Educational and Behavioral Statistics},
  year = {2025},
  volume = {50},
  number = {5},
  pages = {731--764},
  doi = {10.3102/10769986241254286},
}

@Article{Lee2026dce,
  title = {Design-Conditional Prior Elicitation for {Dirichlet} Process Mixtures},
  author = {JoonHo Lee},
  journal = {arXiv preprint},
  year = {2026},
  eprint = {2602.06301},
  archiveprefix = {arXiv},
  url = {https://arxiv.org/abs/2602.06301},
}
```

## Related Work

This package builds on methodological foundations from:

- **Dorazio (2009)**: Original approach for K-based elicitation
- **Lee et al. (2025)**: Informative priors via χ² distribution on K
- **Vicentini & Jermyn (2025)**: Sample-size-independent approaches and weight-based elicitation
- **Zito et al. (2024)**: Stirling-gamma priors and negative binomial approximation

## Support

This project was supported by the Institute of Education Sciences, U.S. Department 
of Education, through Grant R305D240078 to University of Alabama.

<https://ies.ed.gov/use-work/awards/improving-estimation-site-specific-effects-and-their-distribution-multisite-trials-practical-tools>

The opinions expressed are those of the authors and do not represent views of the 
Institute or the U.S. Department of Education.

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

## License
 
MIT © [JoonHo Lee](https://github.com/joonho112)
