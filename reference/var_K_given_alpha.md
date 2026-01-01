# Conditional Variance of K_J Given Alpha

Computes \\\mathrm{Var}(K_J \mid \alpha)\\ under a Dirichlet process
prior.

## Usage

``` r
var_K_given_alpha(J, alpha)
```

## Arguments

- J:

  Sample size (integer \>= 1).

- alpha:

  Concentration parameter (positive numeric, vectorized).

## Value

Numeric vector of conditional variances (same length as `alpha`).

## Details

Uses the trigamma closed form: \$\$v_J(\alpha) = \mu_J(\alpha) -
\alpha^2\\\psi_1(\alpha) - \psi_1(\alpha+J)\\\$\$ where
\\\psi_1(\cdot)\\ is the trigamma function.

This is equivalent to the direct summation: \$\$v_J(\alpha) =
\sum\_{i=1}^{J} \frac{\alpha(i-1)}{(\alpha + i - 1)^2}\$\$

**Key property:** \\0 \< v_J(\alpha) \< \mu_J(\alpha)\\ for all \\\alpha
\> 0\\ (conditional underdispersion).

**Limiting behavior:**

- \\\alpha \to 0^+\\: \\v_J(\alpha) \to 0\\

- \\\alpha \to \infty\\: \\v_J(\alpha) \to 0\\

For numerical stability:

- Values with `alpha < 1e-10` return the limit 0.

- Output is enforced non-negative via `pmax(out, 0)`.

## See also

[`mean_K_given_alpha`](https://joonho112.github.io/DPprior/reference/mean_K_given_alpha.md),
[`moments_K_given_alpha`](https://joonho112.github.io/DPprior/reference/moments_K_given_alpha.md)

## Examples

``` r
var_K_given_alpha(50, 2.0)
#> [1] 4.535558
var_K_given_alpha(50, c(0.5, 1, 2, 5))
#> [1] 1.709074 2.874073 4.535558 7.386114

# Verify underdispersion
J <- 50; alpha <- 2.0
mean_K_given_alpha(J, alpha) > var_K_given_alpha(J, alpha)  # TRUE
#> [1] TRUE
```
