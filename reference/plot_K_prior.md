# Plot Prior PMF of K_J

Plot Prior PMF of K_J

## Usage

``` r
plot_K_prior(
  fit = NULL,
  J = NULL,
  a = NULL,
  b = NULL,
  engine = c("ggplot2", "base"),
  base_size = 11,
  max_k = NULL,
  show_cdf = FALSE,
  show = TRUE
)
```

## Arguments

- fit:

  A DPprior_fit object, or NULL if J, a, b provided directly.

- J:

  Integer; sample size (used if fit is NULL).

- a:

  Numeric; shape parameter (used if fit is NULL).

- b:

  Numeric; rate parameter (used if fit is NULL).

- engine:

  "ggplot2" (default) or "base".

- base_size:

  Base font size.

- max_k:

  Maximum k to display. If NULL, auto-determined by CDF \>= 0.999.

- show_cdf:

  If TRUE, overlay CDF line. Default is FALSE.

- show:

  If TRUE, print the plot.

## Value

A ggplot object or invisible(NULL) for base.

## Examples

``` r
fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8)
#> Warning: HIGH DOMINANCE RISK: P(w1 > 0.5) = 48.1% exceeds 40%.
#>   This may indicate unintended prior behavior (RN-07).
#>   Consider using DPprior_dual() for weight-constrained elicitation.
#>   See ?DPprior_diagnostics for interpretation.
plot_K_prior(fit)



plot_K_prior(J = 50, a = 1.6, b = 1.2)


```
