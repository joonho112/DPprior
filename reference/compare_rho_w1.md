# Compare rho and w1 Distributions

Compares the marginal distributions of rho and w1 under the same
hyperprior alpha ~ Gamma(a, b).

## Usage

``` r
compare_rho_w1(a, b, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- a:

  Numeric; shape parameter of the Gamma prior on alpha (a \> 0).

- b:

  Numeric; rate parameter of the Gamma prior on alpha (b \> 0).

- M:

  Integer; number of quadrature nodes. Default is 80.

## Value

A list containing:

- mean_rho:

  E(rho)

- mean_w1:

  E(w1)

- mean_equal:

  Logical; whether means are equal

- var_rho:

  Var(rho)

- var_w1:

  Var(w1)

- var_ratio:

  Var(rho) / Var(w1)

## Details

While E(rho) = E(w1), the variances differ because:

- `Var(rho | alpha) = 2*alpha / ((1+alpha)^2*(2+alpha)*(3+alpha))`

- `Var(w1 | alpha) = alpha / ((1+alpha)^2*(2+alpha))`

Generally, Var(rho) \< Var(w1) because rho averages over all squared
weights.

## Examples

``` r
if (FALSE) { # \dontrun{
compare_rho_w1(a = 2, b = 1)
compare_rho_w1(a = 1.6, b = 1.22)

} # }
```
