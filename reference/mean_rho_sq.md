# Marginal Second Moment of rho

Computes `E(rho^2 | a, b)` by mixing `E(rho^2 | alpha)` over alpha ~
Gamma(a, b).

## Usage

``` r
mean_rho_sq(a, b, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- a:

  Numeric; shape parameter of the Gamma prior on alpha (a \> 0).

- b:

  Numeric; rate parameter of the Gamma prior on alpha (b \> 0).

- M:

  Integer; number of quadrature nodes. Default is 80.

## Value

Numeric; `E(rho^2 | a, b)`.
