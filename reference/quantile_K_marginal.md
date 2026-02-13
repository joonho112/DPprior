# Quantile of Marginal K Distribution

Computes the \\p\\-th quantile of the marginal distribution of \\K_J\\.

## Usage

``` r
quantile_K_marginal(p, J, a, b, logS, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- p:

  Numeric; probability level(s) in \\\[0, 1\]\\. Can be scalar or
  vector.

- J:

  Integer; sample size.

- a:

  Numeric; shape parameter of Gamma prior (\> 0).

- b:

  Numeric; rate parameter of Gamma prior (\> 0).

- logS:

  Matrix; pre-computed log-Stirling matrix.

- M:

  Integer; number of quadrature nodes (default: 80).

## Value

Integer vector of quantiles (same length as `p`). Each element is the
smallest \\k\\ such that \\P(K_J \leq k) \geq p\\.

## Details

This is the standard quantile definition for discrete distributions:
\\Q(p) = \min\\k : F(k) \geq p\\\\.

The function is vectorized over `p`, allowing efficient computation of
multiple quantiles in a single call.

## See also

[`cdf_K_marginal`](https://joonho112.github.io/DPprior/reference/cdf_K_marginal.md),
[`pmf_K_marginal`](https://joonho112.github.io/DPprior/reference/pmf_K_marginal.md)

Other marginal_K:
[`K_moments()`](https://joonho112.github.io/DPprior/reference/K_moments.md),
[`cdf_K_marginal()`](https://joonho112.github.io/DPprior/reference/cdf_K_marginal.md),
[`exact_K_moments()`](https://joonho112.github.io/DPprior/reference/exact_K_moments.md),
[`mode_K_marginal()`](https://joonho112.github.io/DPprior/reference/mode_K_marginal.md),
[`pmf_K_marginal()`](https://joonho112.github.io/DPprior/reference/pmf_K_marginal.md),
[`summary_K_marginal()`](https://joonho112.github.io/DPprior/reference/summary_K_marginal.md)

## Examples

``` r
logS <- compute_log_stirling(50)

# Single quantile (median)
quantile_K_marginal(0.5, 50, 1.5, 0.5, logS)
#> [1] 8

# Multiple quantiles at once
quantile_K_marginal(c(0.1, 0.25, 0.5, 0.75, 0.9), 50, 1.5, 0.5, logS)
#> [1]  3  5  8 11 15

# Interquartile range
qs <- quantile_K_marginal(c(0.25, 0.75), 50, 1.5, 0.5, logS)
diff(qs)
#> [1] 6
```
