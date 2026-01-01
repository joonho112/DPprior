# Convert Variance Inflation Factor to Variance

Converts a Variance Inflation Factor (VIF) specification to the actual
variance of \\K_J\\.

## Usage

``` r
vif_to_variance(mu_K, vif)
```

## Arguments

- mu_K:

  Numeric; target prior mean of \\K_J\\.

- vif:

  Numeric; Variance Inflation Factor (must be \>= 1 for A1 feasibility).

## Value

Numeric; variance of \\K_J\\ computed as \\(\mu_K - 1) \times
\text{VIF}\\.

## Details

The VIF is defined as: \$\$\text{VIF} = \frac{\sigma^2_K}{\mu_K - 1}\$\$

Interpretation:

- VIF = 1:

  Poisson variance (exact boundary for A1)

- VIF \> 1:

  Overdispersion (required for A1 feasibility)

- VIF \< 1:

  Underdispersion (infeasible for A1, not allowed)

## See also

[`confidence_to_vif`](https://joonho112.github.io/DPprior/reference/confidence_to_vif.md)
for mapping confidence levels to VIF,
[`cv_alpha_to_variance`](https://joonho112.github.io/DPprior/reference/cv_alpha_to_variance.md)
for CV-based specification

## Examples

``` r
# VIF = 2 means variance is twice the Poisson variance
vif_to_variance(mu_K = 5, vif = 2)  # Returns 8
#> [1] 8

# Use with DPprior_a1
fit <- DPprior_a1(J = 50, mu_K = 5, var_K = vif_to_variance(5, 2))
```
