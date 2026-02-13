# CDF of First Stick-Breaking Weight w₁

Computes P(w₁ ≤ x \| a, b) using the closed-form expression derived by
marginalizing over α ~ Gamma(a, b).

## Usage

``` r
cdf_w1(x, a, b)
```

## Arguments

- x:

  Numeric vector. Values outside the unit interval are allowed but are
  mapped to the boundary values of the CDF (0 for x ≤ 0, 1 for x ≥ 1).

- a:

  Numeric; shape parameter of the Gamma prior on α (a \> 0).

- b:

  Numeric; rate parameter of the Gamma prior on α (b \> 0).

## Value

Numeric vector of CDF values F(x \| a, b) with same length as x.

## Details

The unconditional CDF is given by: \$\$F\_{w_1}(x \| a, b) = 1 -
\left(\frac{b}{b - \log(1-x)}\right)^a\$\$

The implementation uses `log1p` and `expm1` for numerical stability,
particularly when the CDF is close to 0 (small x).

## Interpretation

The weight w₁ is in **GEM (size-biased) order**, not ranked by size. It
represents the asymptotic cluster share of a randomly chosen unit,
**not** the largest cluster proportion. See Lee (2026, Section 4) for
details.

## References

Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet
Process Mixtures. *arXiv preprint* arXiv:2602.06301.

Vicentini, C. and Jermyn, I. H. (2025). Prior selection for the
precision parameter of Dirichlet Process Mixtures. arXiv:2502.00864.

## See also

[`quantile_w1`](https://joonho112.github.io/DPprior/reference/quantile_w1.md),
[`density_w1`](https://joonho112.github.io/DPprior/reference/density_w1.md),
[`prob_w1_exceeds`](https://joonho112.github.io/DPprior/reference/prob_w1_exceeds.md)

Other weights_w1:
[`density_w1()`](https://joonho112.github.io/DPprior/reference/density_w1.md),
[`mean_w1()`](https://joonho112.github.io/DPprior/reference/mean_w1.md),
[`prob_w1_exceeds()`](https://joonho112.github.io/DPprior/reference/prob_w1_exceeds.md),
[`quantile_w1()`](https://joonho112.github.io/DPprior/reference/quantile_w1.md),
[`summary_w1()`](https://joonho112.github.io/DPprior/reference/summary_w1.md),
[`var_w1()`](https://joonho112.github.io/DPprior/reference/var_w1.md)

## Examples

``` r
# P(w₁ ≤ 0.3) under standard prior
cdf_w1(0.3, a = 2, b = 1)
#> [1] 0.4566891

# Vectorized computation
cdf_w1(c(0.1, 0.3, 0.5, 0.7), a = 1.6, b = 1.22)
#> [1] 0.1241267 0.3365805 0.5131689 0.6666263
```
