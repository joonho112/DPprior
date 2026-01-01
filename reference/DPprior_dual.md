# Dual-Anchor Prior Calibration

Refines a K_J-calibrated prior to also satisfy weight constraints,
implementing the dual-anchor framework from RN-06.

## Usage

``` r
DPprior_dual(
  fit,
  w1_target,
  lambda = 0.5,
  max_iter = 100L,
  verbose = FALSE,
  M = .QUAD_NODES_DEFAULT,
  loss_type = c("relative", "adaptive", "absolute")
)
```

## Arguments

- fit:

  A `DPprior_fit` object from any K-based calibration method.

- w1_target:

  List specifying the first-weight target. Options:

  `list(prob = list(threshold = 0.5, value = 0.3))`

  :   Constrain \\P(w_1 \> 0.5) = 0.3\\

  `list(quantile = list(prob = 0.9, value = 0.4))`

  :   Constrain 90th percentile of \\w_1\\ = 0.4

  `list(mean = 0.3)`

  :   Constrain \\E\[w_1\] = 0.3\\

- lambda:

  Numeric value between 0 and 1; weight on K_J anchor. Default is 0.5
  (balanced). lambda = 1 recovers K-only fit.

- max_iter:

  Integer; maximum optimization iterations.

- verbose:

  Logical; if TRUE, print progress.

- M:

  Integer; quadrature nodes.

- loss_type:

  Character; "relative" (default), "adaptive", or "absolute". See
  [`dual_anchor_loss`](https://joonho112.github.io/DPprior/reference/dual_anchor_loss.md)
  for details.

## Value

An S3 object of class `"DPprior_fit"`.

## Details

### Critical Fix (2025-12-29)

The original implementation used absolute squared errors for L_K,
causing severe scale mismatch. The optimizer would essentially ignore
the weight constraint. Example improvement at lambda = 0.5:

- ABSOLUTE (broken): 0.5 percent reduction in \\P(w_1 \> 0.5)\\

- RELATIVE (fixed): 9 percent reduction

- ADAPTIVE: 18 percent reduction

## Choosing loss_type

- **"relative"**: Good default. Both losses dimensionless.

- **"adaptive"**: More aggressive. Use when you want lambda = 0.5 to
  give significant weight reduction.

- **"absolute"**: Legacy, not recommended.

## Examples

``` r
# K-only fit
fit_K <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8)
cat("K-only P(w_1 > 0.5):", prob_w1_exceeds(0.5, fit_K$a, fit_K$b), "\n")
#> K-only P(w_1 > 0.5): 0.481478 

# Dual-anchor with relative loss (default)
fit_rel <- DPprior_dual(
  fit_K,
  w1_target = list(prob = list(threshold = 0.5, value = 0.3)),
  lambda = 0.5,
  loss_type = "relative"
)
cat("Relative P(w_1 > 0.5):", fit_rel$dual_anchor$w1_achieved$prob_gt_50, "\n")
#> Relative P(w_1 > 0.5): 0.4379077 

# Dual-anchor with adaptive loss (more aggressive)
fit_adp <- DPprior_dual(
  fit_K,
  w1_target = list(prob = list(threshold = 0.5, value = 0.3)),
  lambda = 0.5,
  loss_type = "adaptive"
)
cat("Adaptive P(w_1 > 0.5):", fit_adp$dual_anchor$w1_achieved$prob_gt_50, "\n")
#> Adaptive P(w_1 > 0.5): 0.3426554 
```
