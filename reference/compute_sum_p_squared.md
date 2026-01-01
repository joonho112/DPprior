# Poissonization Error Bound (Raw Sum of Squared Probabilities)

Computes the raw sum of squared Bernoulli probabilities:
\$\$\sum\_{i=2}^{J} p_i^2 = \alpha^2 \[\psi_1(\alpha+1) -
\psi_1(\alpha+J)\]\$\$ where \\p_i = \alpha / (\alpha + i - 1)\\.

## Usage

``` r
compute_sum_p_squared(J, alpha)
```

## Arguments

- J:

  Integer; sample size (number of observations).

- alpha:

  Numeric; concentration parameter (can be vectorized).

## Value

Numeric vector; sum of squared probabilities for each alpha value.

## Details

This quantity represents the "underdispersion gap" between the
conditional variance of \\K_J \| \alpha\\ and a Poisson with the same
mean.

From the Poisson-binomial representation, \\S_J = K_J - 1 =
\sum\_{i=2}^{J} I_i\\ where \\I_i \sim \text{Bernoulli}(p_i)\\.

This sum equals: \$\$\sum\_{i=2}^{J} p_i^2 = \alpha^2
\[\psi_1(\alpha+1) - \psi_1(\alpha+J)\]\$\$ using the identity for sums
of squared reciprocals.

## References

RN-05: Error Quantification & Guarantees for the A1 Large-J
Approximation

## See also

[`compute_poissonization_bound`](https://joonho112.github.io/DPprior/reference/compute_poissonization_bound.md)
for the full Chen-Stein bound

## Examples

``` r
# Compute raw sum for J=50, alpha=1
compute_sum_p_squared(J = 50, alpha = 1)
#> [1] 0.6251327

# Vectorized over alpha
compute_sum_p_squared(J = 50, alpha = c(0.5, 1, 2))
#> [1] 0.2287007 0.6251327 1.5020688
```
