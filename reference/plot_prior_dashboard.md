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
