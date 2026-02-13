# Quantile of K Given Alpha

Computes the \\p\\-th quantile of \\K_J \mid \alpha\\, defined as the
smallest \\k\\ such that \\P(K_J \leq k \mid \alpha) \geq p\\.

## Usage

``` r
quantile_K_given_alpha(p, J, alpha, logS)
```

## Arguments

- p:

  Numeric; probability level(s) in \\\[0, 1\]\\. Can be scalar or
  vector.

- J:

  Integer; sample size.

- alpha:

  Numeric; DP concentration parameter.

- logS:

  Matrix; pre-computed log-Stirling matrix.

## Value

Integer (or integer vector); the \\p\\-th quantile(s) of \\K_J \mid
\alpha\\.

## Details

For \\p = 0.5\\, this gives the median. Note that for discrete
distributions, the quantile is defined as the smallest value where the
CDF meets or exceeds \\p\\.

Edge cases:

- \\p = 0\\: returns 0 (the first k where CDF \>= 0)

- \\p = 1\\: returns J (the maximum possible K)

## See also

Other conditional_K:
[`cdf_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/cdf_K_given_alpha.md),
[`cv_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/cv_K_given_alpha.md),
[`dispersion_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/dispersion_K_given_alpha.md),
[`log_pmf_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/log_pmf_K_given_alpha.md),
[`mean_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/mean_K_given_alpha.md),
[`mode_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/mode_K_given_alpha.md),
[`moments_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/moments_K_given_alpha.md),
[`pmf_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/pmf_K_given_alpha.md),
[`var_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/var_K_given_alpha.md)

## Examples

``` r
logS <- compute_log_stirling(50)

# Single quantile (median)
quantile_K_given_alpha(0.5, 50, 2.0, logS)
#> [1] 7

# Multiple quantiles
quantile_K_given_alpha(c(0.25, 0.5, 0.75), 50, 2.0, logS)
#> [1] 6 7 8

# Edge cases
quantile_K_given_alpha(0, 50, 2.0, logS)  # Returns 0
#> [1] 0
quantile_K_given_alpha(1, 50, 2.0, logS)  # Returns 50
#> [1] 50
```
