# Map a Qualitative Confidence Level to a Variance Inflation Factor (VIF)

Maps intuitive confidence levels to VIF values for easy prior
specification.

## Usage

``` r
confidence_to_vif(confidence = c("low", "medium", "high"))
```

## Arguments

- confidence:

  Character; one of "low", "medium", or "high".

## Value

Numeric; VIF value (5.0 for low, 2.5 for medium, 1.5 for high).

## Details

The mapping is:

- low:

  VIF = 5.0; high uncertainty about \\K_J\\

- medium:

  VIF = 2.5; moderate uncertainty

- high:

  VIF = 1.5; high confidence (near Poisson boundary)

Higher confidence implies lower variance, which corresponds to lower
VIF. The "high" setting (VIF = 1.5) is close to the A1 feasibility
boundary.

## See also

[`vif_to_variance`](https://joonho112.github.io/DPprior/reference/vif_to_variance.md)
for converting VIF to variance

## Examples

``` r
# Get VIF for medium confidence
vif <- confidence_to_vif("medium")  # Returns 2.5

# Complete workflow
mu_K <- 5
vif <- confidence_to_vif("low")
var_K <- vif_to_variance(mu_K, vif)
fit <- DPprior_a1(J = 50, mu_K = mu_K, var_K = var_K)
```
