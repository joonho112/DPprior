# 4-Panel Prior Dashboard

4-Panel Prior Dashboard

## Usage

``` r
plot_prior_dashboard(
  fit,
  engine = c("ggplot2", "base"),
  base_size = 11,
  ci_level = 0.95,
  title = NULL,
  show = TRUE
)
```

## Arguments

- fit:

  A DPprior_fit object.

- engine:

  "ggplot2" (default) or "base".

- base_size:

  Base font size.

- ci_level:

  Credible interval level.

- title:

  Optional overall title for the dashboard.

- show:

  If TRUE, draw the dashboard.

## Value

A gtable grob (for ggplot2) or invisible(NULL) for base.

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
[`plot_dual_dashboard()`](https://joonho112.github.io/DPprior/reference/plot_dual_dashboard.md),
[`plot_tradeoff_curve()`](https://joonho112.github.io/DPprior/reference/plot_tradeoff_curve.md),
[`plot_tradeoff_dashboard()`](https://joonho112.github.io/DPprior/reference/plot_tradeoff_dashboard.md),
[`plot_w1_prior()`](https://joonho112.github.io/DPprior/reference/plot_w1_prior.md),
[`theme_DPprior()`](https://joonho112.github.io/DPprior/reference/theme_DPprior.md)
