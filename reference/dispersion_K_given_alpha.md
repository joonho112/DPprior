# Dispersion Index for K Given Alpha

Computes the dispersion index (variance-to-mean ratio) of \\K_J \|
\alpha\\. Always \< 1 for this distribution (underdispersion).

## Usage

``` r
dispersion_K_given_alpha(J, alpha)
```

## Arguments

- J:

  Sample size (integer \>= 1).

- alpha:

  Concentration parameter (positive numeric, vectorized).

## Value

Numeric vector of dispersion indices.

## See also

Other conditional_K:
[`cdf_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/cdf_K_given_alpha.md),
[`cv_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/cv_K_given_alpha.md),
[`log_pmf_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/log_pmf_K_given_alpha.md),
[`mean_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/mean_K_given_alpha.md),
[`mode_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/mode_K_given_alpha.md),
[`moments_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/moments_K_given_alpha.md),
[`pmf_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/pmf_K_given_alpha.md),
[`quantile_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/quantile_K_given_alpha.md),
[`var_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/var_K_given_alpha.md)
