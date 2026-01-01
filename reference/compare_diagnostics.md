# Compare Diagnostics Across Multiple Fits

Creates a comparison table of diagnostic metrics for multiple DPprior
fits.

## Usage

``` r
compare_diagnostics(..., M = .QUAD_NODES_DEFAULT)
```

## Arguments

- ...:

  Named DPprior_fit objects to compare.

- M:

  Integer; number of quadrature nodes (default: 80).

## Value

A data frame with diagnostic metrics for each fit.

## Examples

``` r
fit_K <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8)
fit_dual <- DPprior_dual(fit_K,
                         w1_target = list(prob = list(threshold = 0.5, value = 0.3)))
compare_diagnostics(K_only = fit_K, Dual_anchor = fit_dual)
#>        method        a        b     E_K     SD_K      E_w1 P_w1_gt_50
#> 1      K_only 2.036093 1.605054 5.00000 2.828427 0.5014251  0.4814780
#> 2 Dual_anchor 2.575223 1.833605 5.38006 2.816069 0.4669641  0.4379077
#>   P_w1_gt_90 risk     E_rho warnings
#> 1  0.1633817 high 0.5014251        2
#> 2  0.1230800 high 0.4669641        1
```
