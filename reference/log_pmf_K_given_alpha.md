# Log-PMF of K Given Alpha (Antoniak Distribution)

Computes \\\log P(K_J = k \mid \alpha)\\ for \\k = 0, 1, \ldots, J\\.

## Usage

``` r
log_pmf_K_given_alpha(J, alpha, logS)
```

## Arguments

- J:

  Integer; sample size (number of observations, must be \>= 1).

- alpha:

  Numeric; DP concentration parameter (must be positive scalar).

- logS:

  Matrix; pre-computed log-Stirling matrix from
  [`compute_log_stirling`](https://joonho112.github.io/DPprior/reference/compute_log_stirling.md).

## Value

Numeric vector of length \\J+1\\ containing \\\log P(K_J = k \mid
\alpha)\\ for \\k = 0, 1, \ldots, J\\. Note that entry `[1]` corresponds
to \\k=0\\ and always equals `-Inf` (since \\P(K_J = 0) = 0\\).

## Details

Uses the Antoniak distribution formula in log-space: \$\$\log P(K_J = k
\mid \alpha) = \log\|s(J,k)\| + k\log\alpha - \log(\alpha)\_J\$\$

where \\\|s(J,k)\|\\ is the unsigned Stirling number of the first kind
and \\(\alpha)\_J\\ is the rising factorial.

This log-space computation is numerically stable for large \\J\\ where
direct computation would overflow.

## See also

[`pmf_K_given_alpha`](https://joonho112.github.io/DPprior/reference/pmf_K_given_alpha.md)
for normalized PMF,
[`compute_log_stirling`](https://joonho112.github.io/DPprior/reference/compute_log_stirling.md)
for Stirling computation

Other conditional_K:
[`cdf_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/cdf_K_given_alpha.md),
[`cv_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/cv_K_given_alpha.md),
[`dispersion_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/dispersion_K_given_alpha.md),
[`mean_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/mean_K_given_alpha.md),
[`mode_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/mode_K_given_alpha.md),
[`moments_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/moments_K_given_alpha.md),
[`pmf_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/pmf_K_given_alpha.md),
[`quantile_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/quantile_K_given_alpha.md),
[`var_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/var_K_given_alpha.md)

## Examples

``` r
logS <- compute_log_stirling(50)
log_pmf <- log_pmf_K_given_alpha(50, 2.0, logS)

# Convert to probabilities (numerically stable softmax)
pmf <- exp(log_pmf - max(log_pmf))
pmf <- pmf / sum(pmf)
sum(pmf)  # Should be 1
#> [1] 1

# Or use pmf_K_given_alpha() directly
pmf2 <- pmf_K_given_alpha(50, 2.0, logS)
sum(pmf2)  # Should be 1
#> [1] 1
```
