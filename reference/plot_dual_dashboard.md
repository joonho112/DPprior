# Plot Dual-Anchor Extended Dashboard

Creates an extended dashboard for dual-anchor results with 6 panels: (A)
Alpha comparison, (B) K comparison, (C) w1 comparison, (D) Summary
table, (E) Loss decomposition, (F) Parameter trajectory.

## Usage

``` r
plot_dual_dashboard(
  fit_dual,
  tradeoff_data = NULL,
  engine = c("ggplot2", "base"),
  base_size = 10,
  title = NULL,
  show = TRUE
)
```

## Arguments

- fit_dual:

  A DPprior_fit object from DPprior_dual().

- tradeoff_data:

  Optional trade-off curve data for trajectory plot.

- engine:

  "ggplot2" (default) or "base".

- base_size:

  Base font size.

- title:

  Optional title.

- show:

  If TRUE, draw the plot.

## Value

A gtable grob or invisible(NULL).

## References

Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet
Process Mixtures. *arXiv preprint* arXiv:2602.06301.

## See also

[`DPprior_fit`](https://joonho112.github.io/DPprior/reference/DPprior_fit.md)
for fitting,
[`plot.DPprior_fit`](https://joonho112.github.io/DPprior/reference/plot.DPprior_fit.md)
for S3 plot method

Other visualization:
[`DPprior_colors()`](https://joonho112.github.io/DPprior/reference/DPprior_colors.md),
[`plot_K_prior()`](https://joonho112.github.io/DPprior/reference/plot_K_prior.md),
[`plot_alpha_prior()`](https://joonho112.github.io/DPprior/reference/plot_alpha_prior.md),
[`plot_dual_comparison()`](https://joonho112.github.io/DPprior/reference/plot_dual_comparison.md),
[`plot_prior_dashboard()`](https://joonho112.github.io/DPprior/reference/plot_prior_dashboard.md),
[`plot_tradeoff_curve()`](https://joonho112.github.io/DPprior/reference/plot_tradeoff_curve.md),
[`plot_tradeoff_dashboard()`](https://joonho112.github.io/DPprior/reference/plot_tradeoff_dashboard.md),
[`plot_w1_prior()`](https://joonho112.github.io/DPprior/reference/plot_w1_prior.md),
[`theme_DPprior()`](https://joonho112.github.io/DPprior/reference/theme_DPprior.md)
