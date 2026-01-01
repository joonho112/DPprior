# Conditional Mean of K_J Given Alpha

Computes \\\mathbb{E}\[K_J \mid \alpha\]\\ under a Dirichlet process
prior.

## Usage

``` r
mean_K_given_alpha(J, alpha)
```

## Arguments

- J:

  Sample size (integer \>= 1).

- alpha:

  Concentration parameter (positive numeric, vectorized).

## Value

Numeric vector of conditional means (same length as `alpha`).

## Details

Uses the digamma closed form: \$\$\mu_J(\alpha) =
\alpha\\\psi(\alpha+J)-\psi(\alpha)\\\$\$ where \\\psi(\cdot)\\ is the
digamma function.

This is equivalent to the direct summation: \$\$\mu_J(\alpha) =
\sum\_{i=1}^{J} \frac{\alpha}{\alpha + i - 1}\$\$

**Limiting behavior:**

- \\\alpha \to 0^+\\: \\\mu_J(\alpha) \to 1\\

- \\\alpha \to \infty\\: \\\mu_J(\alpha) \to J\\

For numerical stability, values with `alpha < 1e-10` return the limit 1.

## See also

[`var_K_given_alpha`](https://joonho112.github.io/DPprior/reference/var_K_given_alpha.md),
[`moments_K_given_alpha`](https://joonho112.github.io/DPprior/reference/moments_K_given_alpha.md)

## Examples

``` r
mean_K_given_alpha(50, 2.0)
#> [1] 7.037626
mean_K_given_alpha(50, c(0.5, 1, 2, 5))
#> [1]  2.937775  4.499205  7.037626 12.460485

# Limiting behavior
mean_K_given_alpha(50, 1e-10)  # Returns 1
#> [1] 1
mean_K_given_alpha(50, 1e6)    # Returns ~50
#> [1] 49.99878
```
