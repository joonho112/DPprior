# Coefficient of Variation of rho

Computes CV(rho) = SD(rho) / E(rho) under alpha ~ Gamma(a, b).

## Usage

``` r
cv_rho(a, b, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- a:

  Numeric; shape parameter of the Gamma prior on alpha (a \> 0).

- b:

  Numeric; rate parameter of the Gamma prior on alpha (b \> 0).

- M:

  Integer; number of quadrature nodes. Default is 80.

## Value

Numeric; coefficient of variation.

## See also

[`mean_rho`](https://joonho112.github.io/DPprior/reference/mean_rho.md),
[`var_rho`](https://joonho112.github.io/DPprior/reference/var_rho.md)

Other co_clustering:
[`mean_rho()`](https://joonho112.github.io/DPprior/reference/mean_rho.md),
[`mean_rho_given_alpha()`](https://joonho112.github.io/DPprior/reference/mean_rho_given_alpha.md),
[`var_rho()`](https://joonho112.github.io/DPprior/reference/var_rho.md),
[`var_rho_given_alpha()`](https://joonho112.github.io/DPprior/reference/var_rho_given_alpha.md)

## Examples

``` r
cv_rho(a = 2, b = 1)
#> [1] 0.5937869
cv_rho(a = 1.6, b = 1.22)
#> [1] 0.5240019
```
