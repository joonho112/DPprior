# Mean of Marginal K from PMF

Computes \\E\[K_J \mid a, b\]\\ from the marginal PMF.

## Usage

``` r
mean_K_from_marginal_pmf(J, a, b, logS, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- J:

  Integer; sample size.

- a:

  Numeric; shape parameter of Gamma prior.

- b:

  Numeric; rate parameter of Gamma prior.

- logS:

  Matrix; pre-computed log-Stirling matrix.

- M:

  Integer; number of quadrature nodes.

## Value

Numeric; marginal mean.

## Examples

``` r
logS <- compute_log_stirling(50)
mean_K_from_marginal_pmf(50, 1.5, 0.5, logS)
#> [1] 8.355487
```
