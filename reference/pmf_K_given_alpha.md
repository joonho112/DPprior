# PMF of K Given Alpha (Antoniak Distribution)

Computes \\P(K_J = k \mid \alpha)\\ for \\k = 0, 1, \ldots, J\\ using
the Antoniak distribution derived from the Dirichlet process.

## Usage

``` r
pmf_K_given_alpha(J, alpha, logS, normalize = TRUE)
```

## Arguments

- J:

  Integer; sample size (number of observations, must be \>= 1).

- alpha:

  Numeric; DP concentration parameter (must be positive scalar).

- logS:

  Matrix; pre-computed log-Stirling matrix from
  [`compute_log_stirling`](https://joonho112.github.io/DPprior/reference/compute_log_stirling.md).

- normalize:

  Logical; if `TRUE` (default), use softmax normalization for numerical
  stability. If `FALSE`, return raw exponentiated values.

## Value

Numeric vector of length \\J+1\\ containing \\P(K_J = k \mid \alpha)\\
for \\k = 0, 1, \ldots, J\\. Entry `[1]` corresponds to \\k=0\\ and
always equals 0. The vector sums to 1 (when `normalize = TRUE`).

## Details

The Antoniak distribution gives the exact PMF of the number of occupied
clusters \\K_J\\ in a Dirichlet process with concentration parameter
\\\alpha\\: \$\$P(K_J = k \mid \alpha) = \|s(J,k)\|
\frac{\alpha^k}{(\alpha)\_J}\$\$

where \\\|s(J,k)\|\\ is the unsigned Stirling number of the first kind
and \\(\alpha)\_J\\ is the rising factorial.

Key properties:

- \\P(K_J = 0) = 0\\ always (at least one cluster exists)

- \\P(K_J = J) \> 0\\ for all \\\alpha \> 0\\

- As \\\alpha \to 0^+\\, mass concentrates on \\K_J = 1\\

- As \\\alpha \to \infty\\, mass concentrates on \\K_J = J\\

Even though the Antoniak formula is theoretically normalized, converting
log-probabilities to the probability scale via
[`exp()`](https://rdrr.io/r/base/Log.html) can underflow for large `J`
or extreme `alpha`. Setting `normalize=TRUE` mitigates this by applying
[`softmax()`](https://joonho112.github.io/DPprior/reference/softmax.md)
to the log-PMF.

## References

Antoniak, C. E. (1974). Mixtures of Dirichlet Processes with
Applications to Bayesian Nonparametric Problems. *The Annals of
Statistics*, 2(6), 1152-1174.

## See also

[`log_pmf_K_given_alpha`](https://joonho112.github.io/DPprior/reference/log_pmf_K_given_alpha.md)
for log-scale computation,
[`mode_K_given_alpha`](https://joonho112.github.io/DPprior/reference/mode_K_given_alpha.md),
[`cdf_K_given_alpha`](https://joonho112.github.io/DPprior/reference/cdf_K_given_alpha.md),
[`quantile_K_given_alpha`](https://joonho112.github.io/DPprior/reference/quantile_K_given_alpha.md)

Other conditional_K:
[`cdf_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/cdf_K_given_alpha.md),
[`cv_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/cv_K_given_alpha.md),
[`dispersion_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/dispersion_K_given_alpha.md),
[`log_pmf_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/log_pmf_K_given_alpha.md),
[`mean_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/mean_K_given_alpha.md),
[`mode_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/mode_K_given_alpha.md),
[`moments_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/moments_K_given_alpha.md),
[`quantile_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/quantile_K_given_alpha.md),
[`var_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/var_K_given_alpha.md)

## Examples

``` r
# Compute PMF for J=50, alpha=2
logS <- compute_log_stirling(50)
pmf <- pmf_K_given_alpha(50, 2.0, logS)

# Verify normalization
sum(pmf)  # Should be 1
#> [1] 1

# Verify P(K=0) = 0
pmf[1]    # Should be 0
#> [1] 0

# Most likely number of clusters (mode)
which.max(pmf) - 1  # Subtract 1 for 0-indexing
#> [1] 7

# Compare with moments
k_vals <- 0:50
mean_K <- sum(k_vals * pmf)
var_K <- sum(k_vals^2 * pmf) - mean_K^2

# These should match digamma formulas
mean_K_given_alpha(50, 2.0)
#> [1] 7.037626
var_K_given_alpha(50, 2.0)
#> [1] 4.535558
```
