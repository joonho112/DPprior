# Quantile Function of w₁

Computes the inverse CDF: Q(u \| a, b) = F⁻¹(u).

## Usage

``` r
quantile_w1(u, a, b)
```

## Arguments

- u:

  Numeric vector of probability levels in the unit interval. Values u ≤
  0 return 0 and u ≥ 1 return 1.

- a:

  Numeric; shape parameter of the Gamma prior on α (a \> 0).

- b:

  Numeric; rate parameter of the Gamma prior on α (b \> 0).

## Value

Numeric vector of quantile values Q(u \| a, b).

## Details

The quantile function has the closed form: \$\$Q\_{w_1}(u \| a, b) = 1 -
\exp\left(b \left\[1 - (1-u)^{-1/a}\right\]\right)\$\$

The implementation computes (1-u)^(-1/a) in log space for stability when
u is close to 1.

**Numerical Note:** For small values of a (a \< 1) and u close to 1, the
quantile approaches 1 very rapidly and may round to 1.0 in double
precision.

## References

Lee, J. (2025). RN-06: Dual-Anchor Design II, Corollary 2.

## See also

[`cdf_w1`](https://joonho112.github.io/DPprior/reference/cdf_w1.md),
[`summary_w1`](https://joonho112.github.io/DPprior/reference/summary_w1.md)

## Examples

``` r
# Median of w₁
quantile_w1(0.5, a = 2, b = 1)  # ~0.339
#> [1] 0.3391402

# 90th percentile
quantile_w1(0.9, a = 2, b = 1)  # ~0.732
#> [1] 0.8849373
```
