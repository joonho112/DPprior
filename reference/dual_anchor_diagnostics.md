# Dual-Anchor Diagnostic Comparison

Dual-Anchor Diagnostic Comparison

## Usage

``` r
dual_anchor_diagnostics(fit_dual, fit_K_only = NULL, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- fit_dual:

  Dual-anchor DPprior_fit object.

- fit_K_only:

  Optional K-only DPprior_fit object.

- M:

  Integer; quadrature nodes.

## Value

Data frame comparing metrics.

## References

Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet
Process Mixtures. *arXiv preprint* arXiv:2602.06301.

## See also

Other diagnostics:
[`DPprior_diagnostics()`](https://joonho112.github.io/DPprior/reference/DPprior_diagnostics.md),
[`DPprior_error_bounds()`](https://joonho112.github.io/DPprior/reference/DPprior_error_bounds.md),
[`compute_error_landscape()`](https://joonho112.github.io/DPprior/reference/compute_error_landscape.md),
[`compute_weight_diagnostics()`](https://joonho112.github.io/DPprior/reference/compute_weight_diagnostics.md)
