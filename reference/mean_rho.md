# Marginal Mean of rho

Computes E(rho \| a, b) when alpha ~ Gamma(a, b) (shape-rate).

## Usage

``` r
mean_rho(a, b, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- a:

  Numeric; shape parameter of the Gamma prior on alpha (a \> 0).

- b:

  Numeric; rate parameter of the Gamma prior on alpha (b \> 0).

- M:

  Integer; number of quadrature nodes. Default is 80.

## Value

Numeric; `E(rho | a, b)`.

## Details

Uses Gauss-Laguerre quadrature via `integrate_gamma`. A key identity is
`E(rho | alpha) = E(w1 | alpha) = 1/(1+alpha)`, so `E(rho | a, b)`
equals `E(w1 | a, b)` (but the full distributions differ).

**Interpretation:**

- E(rho) \> 0.5: High prior co-clustering probability

- E(rho) in (0.2, 0.5): Moderate co-clustering

- E(rho) \< 0.2: Low co-clustering (fragmented prior)

## References

Lee, J. (2025). RN-06: Dual-Anchor Design II, Section 3.3.

## See also

[`var_rho`](https://joonho112.github.io/DPprior/reference/var_rho.md),
[`cv_rho`](https://joonho112.github.io/DPprior/reference/cv_rho.md),
[`mean_w1`](https://joonho112.github.io/DPprior/reference/mean_w1.md)

## Examples

``` r
mean_rho(a = 2, b = 1)
#> [1] 0.4036526
mean_rho(a = 1.6, b = 1.22)
#> [1] 0.508368
```
