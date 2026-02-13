# Compute A1 Error Landscape

Computes error metrics over a grid of (J, alpha) values for
visualization.

## Usage

``` r
compute_error_landscape(J_seq, alpha_seq, cJ_fun = log)
```

## Arguments

- J_seq:

  Numeric vector; sequence of J values.

- alpha_seq:

  Numeric vector; sequence of alpha values.

- cJ_fun:

  Function; scaling constant function (default: log).

## Value

A data frame with columns: J, alpha, lambda_exact, lambda_approx,
pois_raw, pois_bound, lin_bound, total_tv.

## Details

This function is useful for creating error landscape visualizations as
shown in Lee (2026, Section 3.3).

## See also

Other diagnostics:
[`DPprior_diagnostics()`](https://joonho112.github.io/DPprior/reference/DPprior_diagnostics.md),
[`DPprior_error_bounds()`](https://joonho112.github.io/DPprior/reference/DPprior_error_bounds.md),
[`compute_weight_diagnostics()`](https://joonho112.github.io/DPprior/reference/compute_weight_diagnostics.md),
[`dual_anchor_diagnostics()`](https://joonho112.github.io/DPprior/reference/dual_anchor_diagnostics.md)

## Examples

``` r
# Create error landscape
landscape <- compute_error_landscape(
  J_seq = c(25, 50, 100),
  alpha_seq = c(0.5, 1, 2, 5)
)
print(landscape)
#>      J alpha lambda_exact lambda_approx  pois_raw pois_bound   lin_bound
#> 1   25   0.5     1.591226      1.609438 0.2237019 0.11195092 0.007191254
#> 2   50   0.5     1.937775      1.956012 0.2287007 0.10102428 0.006529898
#> 3  100   0.5     2.284342      2.302585 0.2312006 0.09090357 0.006019094
#> 4   25   1.0     2.815958      3.218876 0.6057234 0.20223044 0.114762181
#> 5   50   1.0     3.499205      3.912023 0.6251327 0.17325087 0.106279606
#> 6  100   1.0     4.187378      4.605170 0.6349839 0.14933953 0.098874211
#> 7   25   2.0     4.708839      6.437752 1.4288108 0.30069611 0.357966296
#> 8   50   2.0     6.037626      7.824046 1.5020688 0.24819075 0.332807789
#> 9  100   2.0     7.394557      9.210340 1.5403277 0.20817759 0.309892875
#> 10  25   5.0     8.391602     16.094379 3.6856974 0.43911299 1.000000000
#> 11  50   5.0    11.460485     19.560115 4.0743712 0.35551097 0.993226480
#> 12 100   5.0    14.715366     23.025851 4.2938413 0.29179291 0.927912567
#>      total_tv
#> 1  0.11914218
#> 2  0.10755418
#> 3  0.09692267
#> 4  0.31699262
#> 5  0.27953047
#> 6  0.24821374
#> 7  0.65866241
#> 8  0.58099854
#> 9  0.51807046
#> 10 1.00000000
#> 11 1.00000000
#> 12 1.00000000
```
