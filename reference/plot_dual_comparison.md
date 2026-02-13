# Plot Dual-Anchor Comparison Dashboard

Creates a comparison dashboard showing K-only vs dual-anchor solutions.
Displays changes in alpha, K, and w1 distributions side-by-side.

## Usage

``` r
plot_dual_comparison(
  fit_dual,
  fit_K_only = NULL,
  engine = c("ggplot2", "base"),
  base_size = 10,
  title = NULL,
  show = TRUE
)
```

## Arguments

- fit_dual:

  A DPprior_fit object from DPprior_dual().

- fit_K_only:

  Optional K-only fit. If NULL, extracted from
  fit_dual\$dual_anchor\$init.

- engine:

  "ggplot2" (default) or "base".

- base_size:

  Base font size.

- title:

  Optional title for the dashboard.

- show:

  If TRUE, draw the plot.

## Value

A gtable grob or list of ggplot objects.

## See also

[`DPprior_fit`](https://joonho112.github.io/DPprior/reference/DPprior_fit.md)
for fitting,
[`plot.DPprior_fit`](https://joonho112.github.io/DPprior/reference/plot.DPprior_fit.md)
for S3 plot method

Other visualization:
[`DPprior_colors()`](https://joonho112.github.io/DPprior/reference/DPprior_colors.md),
[`plot_K_prior()`](https://joonho112.github.io/DPprior/reference/plot_K_prior.md),
[`plot_alpha_prior()`](https://joonho112.github.io/DPprior/reference/plot_alpha_prior.md),
[`plot_dual_dashboard()`](https://joonho112.github.io/DPprior/reference/plot_dual_dashboard.md),
[`plot_prior_dashboard()`](https://joonho112.github.io/DPprior/reference/plot_prior_dashboard.md),
[`plot_tradeoff_curve()`](https://joonho112.github.io/DPprior/reference/plot_tradeoff_curve.md),
[`plot_tradeoff_dashboard()`](https://joonho112.github.io/DPprior/reference/plot_tradeoff_dashboard.md),
[`plot_w1_prior()`](https://joonho112.github.io/DPprior/reference/plot_w1_prior.md),
[`theme_DPprior()`](https://joonho112.github.io/DPprior/reference/theme_DPprior.md)

## Examples

``` r
fit_K <- DPprior_fit(J = 50, mu_K = 5, var_K = 8)
#> Warning: HIGH DOMINANCE RISK: P(w1 > 0.5) = 48.1% exceeds 40%.
#>   This may indicate unintended prior behavior (Lee, 2026).
#>   Consider using DPprior_dual() for weight-constrained elicitation.
#>   See ?DPprior_diagnostics for interpretation.
fit_dual <- DPprior_dual(fit_K,
  w1_target = list(prob = list(threshold = 0.5, value = 0.3)),
  lambda = 0.5)
plot_dual_comparison(fit_dual)

#> TableGrob (3 x 2) "dpprior_dashboard": 5 grobs
#>   z     cells              name                 grob
#> 1 1 (2-2,1-1) dpprior_dashboard       gtable[layout]
#> 2 2 (3-3,1-1) dpprior_dashboard       gtable[layout]
#> 3 3 (2-2,2-2) dpprior_dashboard       gtable[layout]
#> 4 4 (3-3,2-2) dpprior_dashboard       gtable[layout]
#> 5 5 (1-1,1-2) dpprior_dashboard text[GRID.text.1555]
```
