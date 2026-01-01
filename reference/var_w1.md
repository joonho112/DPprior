# Variance of w₁

Computes Var(w₁ \| a, b) using the law of total variance.

## Usage

``` r
var_w1(a, b, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- a:

  Numeric; shape parameter of the Gamma prior on α (a \> 0).

- b:

  Numeric; rate parameter of the Gamma prior on α (b \> 0).

- M:

  Integer; number of quadrature nodes. Default is 80.

## Value

Numeric; Var(w₁).

## Details

Uses the law of total variance: \$\$Var(w_1) = E\[Var(w_1 \| \alpha)\] +
Var(E\[w_1 \| \alpha\])\$\$

where w₁ \| α ~ Beta(1, α), so:

- E(w₁ \| α) = 1/(1+α)

- Var(w₁ \| α) = α / ((1+α)²(2+α))

## References

Lee, J. (2025). RN-06: Dual-Anchor Design II.

## See also

[`mean_w1`](https://joonho112.github.io/DPprior/reference/mean_w1.md),
[`summary_w1`](https://joonho112.github.io/DPprior/reference/summary_w1.md)

## Examples

``` r
var_w1(a = 2, b = 1)       # ~0.090
#> [1] 0.08968429
var_w1(a = 1.6, b = 1.22)  # ~0.105
#> [1] 0.1052062
```
