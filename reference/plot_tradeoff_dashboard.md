# Plot Trade-off Multi-Panel

Creates a multi-panel view of the trade-off curve showing multiple
metrics.

## Usage

``` r
plot_tradeoff_dashboard(
  tradeoff_data,
  w1_target_prob = NULL,
  engine = c("ggplot2", "base"),
  base_size = 10,
  title = NULL,
  show = TRUE
)
```

## Arguments

- tradeoff_data:

  Data frame from compute_tradeoff_curve().

- w1_target_prob:

  Optional target probability for P(w1 \> 0.5).

- engine:

  "ggplot2" (default) or "base".

- base_size:

  Base font size.

- title:

  Optional title.

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
[`plot_dual_comparison()`](https://joonho112.github.io/DPprior/reference/plot_dual_comparison.md),
[`plot_dual_dashboard()`](https://joonho112.github.io/DPprior/reference/plot_dual_dashboard.md),
[`plot_prior_dashboard()`](https://joonho112.github.io/DPprior/reference/plot_prior_dashboard.md),
[`plot_tradeoff_curve()`](https://joonho112.github.io/DPprior/reference/plot_tradeoff_curve.md),
[`plot_w1_prior()`](https://joonho112.github.io/DPprior/reference/plot_w1_prior.md),
[`theme_DPprior()`](https://joonho112.github.io/DPprior/reference/theme_DPprior.md)
