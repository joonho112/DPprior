# Mode of Marginal K Distribution

Computes the mode (most likely value) of the marginal distribution of
\\K_J\\.

## Usage

``` r
mode_K_marginal(J, a, b, logS, M = .QUAD_NODES_DEFAULT)
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

Integer; the value \\k\\ that maximizes \\P(K_J = k \mid a, b)\\.

## Details

The mode is always \>= 1 since \\P(K_J = 0) = 0\\.

## See also

[`pmf_K_marginal`](https://joonho112.github.io/DPprior/reference/pmf_K_marginal.md),
[`summary_K_marginal`](https://joonho112.github.io/DPprior/reference/summary_K_marginal.md)

Other marginal_K:
[`K_moments()`](https://joonho112.github.io/DPprior/reference/K_moments.md),
[`cdf_K_marginal()`](https://joonho112.github.io/DPprior/reference/cdf_K_marginal.md),
[`exact_K_moments()`](https://joonho112.github.io/DPprior/reference/exact_K_moments.md),
[`pmf_K_marginal()`](https://joonho112.github.io/DPprior/reference/pmf_K_marginal.md),
[`quantile_K_marginal()`](https://joonho112.github.io/DPprior/reference/quantile_K_marginal.md),
[`summary_K_marginal()`](https://joonho112.github.io/DPprior/reference/summary_K_marginal.md)

## Examples

``` r
logS <- compute_log_stirling(50)
mode_K_marginal(50, 1.5, 0.5, logS)
#> [1] 6
```
