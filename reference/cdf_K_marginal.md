# CDF of Marginal K Distribution

Computes the cumulative distribution function \\P(K_J \leq k \mid a,
b)\\ for \\k = 0, 1, \ldots, J\\.

## Usage

``` r
cdf_K_marginal(J, a, b, logS, M = .QUAD_NODES_DEFAULT)
```

## Arguments

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

Numeric vector of length \\J+1\\ containing \\P(K_J \leq k \mid a, b)\\
for \\k = 0, 1, \ldots, J\\.

## Details

The CDF satisfies:

- \\F(0) = 0\\ (since \\P(K_J = 0) = 0\\)

- \\F(J) = 1\\

- \\F(k)\\ is non-decreasing in \\k\\

## See also

[`pmf_K_marginal`](https://joonho112.github.io/DPprior/reference/pmf_K_marginal.md),
[`quantile_K_marginal`](https://joonho112.github.io/DPprior/reference/quantile_K_marginal.md)

Other marginal_K:
[`K_moments()`](https://joonho112.github.io/DPprior/reference/K_moments.md),
[`exact_K_moments()`](https://joonho112.github.io/DPprior/reference/exact_K_moments.md),
[`mode_K_marginal()`](https://joonho112.github.io/DPprior/reference/mode_K_marginal.md),
[`pmf_K_marginal()`](https://joonho112.github.io/DPprior/reference/pmf_K_marginal.md),
[`quantile_K_marginal()`](https://joonho112.github.io/DPprior/reference/quantile_K_marginal.md),
[`summary_K_marginal()`](https://joonho112.github.io/DPprior/reference/summary_K_marginal.md)

## Examples

``` r
logS <- compute_log_stirling(50)
cdf <- cdf_K_marginal(50, 1.5, 0.5, logS)

# Verify CDF ends at 1
cdf[51]
#> [1] 1

# P(K <= 10)
cdf[11]
#> [1] 0.7009381
```
