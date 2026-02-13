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
if (FALSE) { # \dontrun{
fit_K <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8)
fit_dual <- DPprior_dual(fit_K,
                         w1_target = list(prob = list(threshold = 0.5, value = 0.3)))
compare_diagnostics(K_only = fit_K, Dual_anchor = fit_dual)

} # }
```
