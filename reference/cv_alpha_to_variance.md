# Convert CV(alpha) to Variance

Converts a coefficient of variation specification for \\\alpha\\ to the
implied variance of \\K_J\\ under the A1 approximation.

## Usage

``` r
cv_alpha_to_variance(mu_K, cv_alpha)
```

## Arguments

- mu_K:

  Numeric; target prior mean of \\K_J\\.

- cv_alpha:

  Numeric; target coefficient of variation for \\\alpha\\.

## Value

Numeric; implied variance of \\K_J\\.

## Details

Under the A1 approximation: \$\$\text{CV}(\alpha) = 1/\sqrt{a} =
\frac{\sqrt{\sigma^2_K - m}}{m}\$\$

where \\m = \mu_K - 1\\. Inverting: \$\$\sigma^2_K = m +
(\text{CV}(\alpha) \cdot m)^2 = m(1 + \text{CV}(\alpha)^2 \cdot m)\$\$

## See also

[`vif_to_variance`](https://joonho112.github.io/DPprior/reference/vif_to_variance.md)
for VIF-based specification

## Examples

``` r
# CV(alpha) = 0.5 means moderate prior concentration
var_K <- cv_alpha_to_variance(mu_K = 5, cv_alpha = 0.5)

# Verify round-trip
fit <- DPprior_a1(J = 50, mu_K = 5, var_K = var_K)
1 / sqrt(fit$a)  # Should be approximately 0.5
#> [1] 0.5
```
