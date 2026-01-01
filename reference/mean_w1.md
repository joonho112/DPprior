# Mean of w₁

Computes E(w₁ \| a, b) via Gauss-Laguerre quadrature.

## Usage

``` r
mean_w1(a, b, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- a:

  Numeric; shape parameter of the Gamma prior on α (a \> 0).

- b:

  Numeric; rate parameter of the Gamma prior on α (b \> 0).

- M:

  Integer; number of quadrature nodes. Default is 80.

## Value

Numeric; E(w₁).

## Details

The expectation is computed using the identity: \$\$E\[w_1 \| a, b\] =
E\left\[\frac{1}{1+\alpha}\right\] = I_1(a, b)\$\$

where the integral is evaluated via Gauss-Laguerre quadrature.

**Key identity:** E(w₁ \| a, b) = E(ρ \| a, b) where ρ = Σwₕ² is the
co-clustering probability.

## References

Lee, J. (2025). RN-06: Dual-Anchor Design II, §2.3.

## See also

[`var_w1`](https://joonho112.github.io/DPprior/reference/var_w1.md),
[`summary_w1`](https://joonho112.github.io/DPprior/reference/summary_w1.md)

## Examples

``` r
mean_w1(a = 2, b = 1)       # ~0.404
#> [1] 0.4036526
mean_w1(a = 1.6, b = 1.22)  # ~0.508
#> [1] 0.508368
```
