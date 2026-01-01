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
