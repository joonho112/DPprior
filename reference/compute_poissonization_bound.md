# Poissonization Error Bound (Chen-Stein / Le Cam Bound)

Computes an upper bound on the conditional total variation distance
between the shifted cluster count \\S_J = K_J - 1\\ and a Poisson law
with the same mean (Poissonization error).

## Usage

``` r
compute_poissonization_bound(J, alpha, raw = FALSE)
```

## Arguments

- J:

  Integer; sample size (number of observations).

- alpha:

  Numeric; concentration parameter (can be vectorized).

- raw:

  Logical; if TRUE, return just sum(p_i^2) without the prefactor.
  Default is FALSE.

## Value

Numeric vector; upper bound on \\d\_{TV}(S_J \| \alpha,
\text{Poisson}(\lambda_J(\alpha)))\\.

## Details

Under the CRP representation, \\S_J = \sum\_{i=2}^J I_i\\ where \\I_i
\sim \text{Bernoulli}(p_i)\\ and \\p_i = \alpha / (\alpha + i - 1)\\.

A standard Chen-Stein/Le Cam bound gives: \$\$d\_{TV}(S_J,
\text{Poisson}(\lambda)) \le \frac{1 - e^{-\lambda}}{\lambda}
\sum\_{i=2}^J p_i^2\$\$ where \\\lambda = \sum\_{i=2}^J p_i = E\[S_J \|
\alpha\]\\.

The prefactor \\(1 - e^{-\lambda})/\lambda\\ is always in (0, 1\] and
approaches 1 as \\\lambda \to 0\\. This provides a tighter bound than
simply using \\\sum p_i^2\\ alone.

The returned value is capped at 1 (since total variation is always
between 0 and 1).

## References

Le Cam, L. (1960). An approximation theorem for the Poisson binomial
distribution. *Pacific Journal of Mathematics*, 10(4), 1181-1197.

Chen, L. H. Y. (1975). Poisson approximation for dependent trials. *The
Annals of Probability*, 3(3), 534-545.

Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet
Process Mixtures. *arXiv preprint* arXiv:2602.06301.

## See also

[`compute_sum_p_squared`](https://joonho112.github.io/DPprior/reference/compute_sum_p_squared.md),
[`compute_linearization_bound`](https://joonho112.github.io/DPprior/reference/compute_linearization_bound.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Full Chen-Stein bound
compute_poissonization_bound(J = 50, alpha = 1)

# Raw bound (sum of p_i^2)
compute_poissonization_bound(J = 50, alpha = 1, raw = TRUE)

# Vectorized
compute_poissonization_bound(J = 50, alpha = c(0.5, 1, 2, 5))

} # }
```
