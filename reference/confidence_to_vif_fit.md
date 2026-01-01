# Convert Confidence Level to Variance Inflation Factor

Maps qualitative confidence levels to quantitative VIF values. The VIF
determines prior variance as: Var(K) = VIF \* (mu_K - 1).

## Usage

``` r
confidence_to_vif_fit(confidence)
```

## Arguments

- confidence:

  Character; one of "low", "medium", or "high".

## Value

Numeric; the corresponding VIF value.

## Details

The mapping is:

- "low": VIF = 5.0 (very high uncertainty about cluster count)

- "medium": VIF = 2.5 (moderate uncertainty)

- "high": VIF = 1.5 (strong prior belief)

These values are chosen to provide a wide range of uncertainty levels.
VIF = 5 corresponds to "I have little idea how many clusters there are."

## See also

[`vif_to_variance_fit`](https://joonho112.github.io/DPprior/reference/vif_to_variance_fit.md),
[`DPprior_fit`](https://joonho112.github.io/DPprior/reference/DPprior_fit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
confidence_to_vif_fit("low")     # 5.0
confidence_to_vif_fit("medium")  # 2.5
confidence_to_vif_fit("high")    # 1.5
} # }
```
