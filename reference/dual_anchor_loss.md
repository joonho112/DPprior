# Dual-Anchor Loss Function

Computes the combined loss for K_J and first-weight targets used in the
dual-anchor calibration framework.

## Usage

``` r
dual_anchor_loss(
  a,
  b,
  J,
  K_target,
  w1_target,
  lambda,
  M = .QUAD_NODES_DEFAULT,
  loss_type = c("relative", "adaptive", "absolute"),
  L_K_scale = NULL,
  L_w_scale = NULL
)
```

## Arguments

- a:

  Numeric; shape parameter of the Gamma prior (a \> 0).

- b:

  Numeric; rate parameter of the Gamma prior (b \> 0).

- J:

  Integer; sample size.

- K_target:

  List with components `mu_K` and `var_K`.

- w1_target:

  List specifying the first-weight target (quantile, prob, or mean).

- lambda:

  Numeric value between 0 and 1; weight on K_J loss component.

- M:

  Integer; number of Gauss-Laguerre quadrature nodes.

- loss_type:

  Character; type of loss scaling:

  "relative"

  :   (DEFAULT) Normalized squared errors - both losses dimensionless
      and comparable. RECOMMENDED for most cases.

  "adaptive"

  :   Scaled by initial loss magnitudes - gives more aggressive weight
      reduction. Use when lambda=0.5 should truly mean "equal
      importance".

  "absolute"

  :   Legacy behavior - NOT RECOMMENDED, causes optimizer to essentially
      ignore weight constraint.

- L_K_scale, L_w_scale:

  Numeric; scaling factors for adaptive loss. Only used when
  `loss_type = "adaptive"`. If NULL, computed internally.

## Value

Numeric; the combined loss value.

## Details

### Loss Scaling Comparison

At typical K-only solution with \\P(w_1 \> 0.5)\\ approximately 0.48,
target = 0.30:

- **ABSOLUTE**: L_K approximately 0, L_w approximately 0.03. Any move
  causes L_K \>\> L_w. Optimizer stays at K-only. BROKEN.

- **RELATIVE**: L_K_rel approximately 0, L_w approximately 0.03. Both
  dimensionless. lambda = 0.5 gives approximately 9 percent reduction in
  \\P(w_1 \> 0.5)\\.

- **ADAPTIVE**: Losses normalized by initial magnitudes. lambda = 0.5
  gives approximately 18 percent reduction in \\P(w_1 \> 0.5)\\.

## See also

[`DPprior_dual`](https://joonho112.github.io/DPprior/reference/DPprior_dual.md),
[`compute_tradeoff_curve`](https://joonho112.github.io/DPprior/reference/compute_tradeoff_curve.md)
