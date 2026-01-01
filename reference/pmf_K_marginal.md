# Marginal PMF of K_J under Gamma Hyperprior

Computes \\P(K_J = k \mid a, b)\\ for \\k = 0, 1, \ldots, J\\ when
\\\alpha \sim \mathrm{Gamma}(a, b)\\ (shape-rate parameterization).

## Usage

``` r
pmf_K_marginal(J, a, b, logS, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- J:

  Integer; sample size (positive integer \>= 1).

- a:

  Numeric; shape parameter of Gamma prior (\> 0).

- b:

  Numeric; rate parameter of Gamma prior (\> 0).

- logS:

  Matrix; pre-computed log-Stirling matrix from
  [`compute_log_stirling`](https://joonho112.github.io/DPprior/reference/compute_log_stirling.md).

- M:

  Integer; number of quadrature nodes (default: 80).

## Value

Numeric vector of length \\J+1\\ containing \\P(K_J = k \mid a, b)\\ for
\\k = 0, 1, \ldots, J\\. Entry `[1]` corresponds to \\k=0\\ and always
equals 0. The vector sums to 1.

## Details

Uses Gauss-Laguerre quadrature to numerically evaluate: \$\$P(K_J = k
\mid a, b) = \int_0^\infty P(K_J = k \mid \alpha) \cdot g\_{a,b}(\alpha)
d\alpha\$\$ \$\$\approx \sum\_{m=1}^M \tilde{w}\_m \cdot P(K_J = k \mid
\alpha_m)\$\$

where \\P(K_J = k \mid \alpha)\\ is the Antoniak distribution from
Module 04 and \\(\alpha_m, \tilde{w}\_m)\\ are the transformed
quadrature nodes and normalized weights from Module 02.

**Implementation:** All mixing is performed in log-space for numerical
stability. This is critical for large J or extreme parameter values.

**Key properties:**

- \\P(K_J = 0) = 0\\ always (at least one cluster exists)

- The PMF sums to 1

- Moments from the PMF match
  [`exact_K_moments()`](https://joonho112.github.io/DPprior/reference/exact_K_moments.md)
  within numerical tolerance

- Mode is typically near \\E\[K_J\]\\ but may differ

## References

Antoniak, C. E. (1974). Mixtures of Dirichlet Processes with
Applications to Bayesian Nonparametric Problems. *The Annals of
Statistics*, 2(6), 1152-1174.

## See also

[`log_pmf_K_marginal`](https://joonho112.github.io/DPprior/reference/log_pmf_K_marginal.md)
for log-scale computation,
[`pmf_K_given_alpha`](https://joonho112.github.io/DPprior/reference/pmf_K_given_alpha.md)
for conditional PMF,
[`exact_K_moments`](https://joonho112.github.io/DPprior/reference/exact_K_moments.md)
for marginal moments,
[`cdf_K_marginal`](https://joonho112.github.io/DPprior/reference/cdf_K_marginal.md),
[`quantile_K_marginal`](https://joonho112.github.io/DPprior/reference/quantile_K_marginal.md),
[`mode_K_marginal`](https://joonho112.github.io/DPprior/reference/mode_K_marginal.md),
[`summary_K_marginal`](https://joonho112.github.io/DPprior/reference/summary_K_marginal.md)

## Examples

``` r
# Pre-compute Stirling numbers
logS <- compute_log_stirling(50)

# Compute marginal PMF for J=50, Gamma(1.5, 0.5) prior
pmf <- pmf_K_marginal(50, 1.5, 0.5, logS)

# Verify normalization
sum(pmf)
#> [1] 1

# Most likely number of clusters
which.max(pmf) - 1
#> [1] 6

# Compare mean with exact_K_moments
k_vals <- 0:50
mean_pmf <- sum(k_vals * pmf)
exact <- exact_K_moments(50, 1.5, 0.5)
abs(mean_pmf - exact$mean)
#> [1] 7.105427e-15
```
