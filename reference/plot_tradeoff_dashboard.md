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
