# Variance Inflation Ratio

Computes the ratio \\Var(K_J) / E\[K_J\]\\ as an overdispersion measure.

## Usage

``` r
variance_inflation_ratio(J, a, b, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- J:

  Integer; sample size.

- a:

  Numeric; shape parameter of Gamma prior.

- b:

  Numeric; rate parameter of Gamma prior.

- M:

  Integer; number of quadrature nodes (default: 80).

## Value

Numeric; the variance inflation ratio (VIR).

## Details

For a Poisson distribution, this ratio equals 1. For the marginal
distribution of \\K_J\\ under a Gamma prior on \\\alpha\\, this ratio is
typically \> 1, indicating overdispersion.

This ratio is useful for:

- Diagnosing the appropriateness of Poisson approximations

- Comparing different prior specifications

- Understanding the "spread" induced by uncertainty in \\\alpha\\

## Examples

``` r
if (FALSE) { # \dontrun{
# Typical overdispersion
vir <- variance_inflation_ratio(50, 2.0, 1.0)
vir > 1  # TRUE: overdispersed

# Compare across prior specifications
variance_inflation_ratio(50, 2.0, 0.5)  # Higher uncertainty in alpha
variance_inflation_ratio(50, 8.0, 4.0)  # Same mean, lower variance in alpha

} # }
```
