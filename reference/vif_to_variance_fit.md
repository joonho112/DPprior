# Convert VIF to Target Variance

Computes the target variance from a Variance Inflation Factor.

## Usage

``` r
vif_to_variance_fit(mu_K, vif)
```

## Arguments

- mu_K:

  Numeric; target expected number of clusters.

- vif:

  Numeric; Variance Inflation Factor (must be \> 1).

## Value

Numeric; target variance: `vif * (mu_K - 1)`.

## Details

The VIF relates variance to the shifted mean (mu_K - 1): \$\$Var(K_J) =
VIF \cdot (\mu_K - 1)\$\$

This parameterization is natural for the shifted NegBin approximation.

## See also

[`confidence_to_vif_fit`](https://joonho112.github.io/DPprior/reference/confidence_to_vif_fit.md),
[`DPprior_fit`](https://joonho112.github.io/DPprior/reference/DPprior_fit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
vif_to_variance_fit(5, 2.5)  # 10 = 2.5 * (5 - 1)
} # }
```
