# Survival Function of w₁

Computes P(w₁ \> t \| a, b) = 1 - F(t), the probability that the first
stick-breaking weight exceeds threshold t.

## Usage

``` r
prob_w1_exceeds(t, a, b)
```

## Arguments

- t:

  Numeric vector of thresholds. Values outside the unit interval are
  allowed but are mapped to the boundary values (1 for t ≤ 0, 0 for t ≥
  1).

- a:

  Numeric; shape parameter of the Gamma prior on α (a \> 0).

- b:

  Numeric; rate parameter of the Gamma prior on α (b \> 0).

## Value

Numeric vector of survival probabilities.

## Details

The survival function has the closed form: \$\$P(w_1 \> t \| a, b) =
\left(\frac{b}{b - \log(1-t)}\right)^a\$\$

This is a **key quantity for dominance risk assessment** (Lee, 2026,
Section 4). A large P(w₁ \> 0.5) indicates high prior probability that a
single cluster dominates the mixture.

## Dominance Risk Interpretation

- P(w₁ \> 0.5) ≈ 0.5: moderate dominance risk

- P(w₁ \> 0.9) ≈ 0.1: low extreme dominance risk

## References

Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet
Process Mixtures. *arXiv preprint* arXiv:2602.06301.

## See also

[`cdf_w1`](https://joonho112.github.io/DPprior/reference/cdf_w1.md),
[`summary_w1`](https://joonho112.github.io/DPprior/reference/summary_w1.md)

Other weights_w1:
[`cdf_w1()`](https://joonho112.github.io/DPprior/reference/cdf_w1.md),
[`density_w1()`](https://joonho112.github.io/DPprior/reference/density_w1.md),
[`mean_w1()`](https://joonho112.github.io/DPprior/reference/mean_w1.md),
[`quantile_w1()`](https://joonho112.github.io/DPprior/reference/quantile_w1.md),
[`summary_w1()`](https://joonho112.github.io/DPprior/reference/summary_w1.md),
[`var_w1()`](https://joonho112.github.io/DPprior/reference/var_w1.md)

## Examples

``` r
# P(w₁ > 0.5): "dominant cluster" probability
prob_w1_exceeds(0.5, a = 1.6, b = 1.22)  # ~0.487 (Lee et al. DP-inform)
#> [1] 0.4868311
prob_w1_exceeds(0.5, a = 2, b = 1)       # ~0.349
#> [1] 0.3488274
```
