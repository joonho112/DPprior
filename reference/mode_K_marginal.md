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

## Examples

``` r
logS <- compute_log_stirling(50)
mode_K_marginal(50, 1.5, 0.5, logS)
#> [1] 6
```
