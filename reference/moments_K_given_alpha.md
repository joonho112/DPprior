# Conditional Mean and Variance of K_J Given Alpha

Convenience wrapper returning both conditional moments in one call.

## Usage

``` r
moments_K_given_alpha(J, alpha)
```

## Arguments

- J:

  Sample size (integer \>= 1).

- alpha:

  Concentration parameter (positive numeric, vectorized).

## Value

If `alpha` is scalar, a named numeric vector `c(mean = ..., var = ...)`.
If `alpha` is a vector, a numeric matrix with two columns `mean` and
`var` (one row per element of `alpha`).

## See also

[`mean_K_given_alpha`](https://joonho112.github.io/DPprior/reference/mean_K_given_alpha.md),
[`var_K_given_alpha`](https://joonho112.github.io/DPprior/reference/var_K_given_alpha.md)

Other conditional_K:
[`cdf_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/cdf_K_given_alpha.md),
[`cv_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/cv_K_given_alpha.md),
[`dispersion_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/dispersion_K_given_alpha.md),
[`log_pmf_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/log_pmf_K_given_alpha.md),
[`mean_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/mean_K_given_alpha.md),
[`mode_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/mode_K_given_alpha.md),
[`pmf_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/pmf_K_given_alpha.md),
[`quantile_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/quantile_K_given_alpha.md),
[`var_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/var_K_given_alpha.md)

## Examples

``` r
moments_K_given_alpha(50, 2.0)
#>     mean      var 
#> 7.037626 4.535558 
moments_K_given_alpha(50, c(0.5, 1, 2, 5))
#>           mean      var
#> [1,]  2.937775 1.709074
#> [2,]  4.499205 2.874073
#> [3,]  7.037626 4.535558
#> [4,] 12.460485 7.386114
```
