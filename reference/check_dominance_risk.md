# Quick Dominance Risk Check

Fast check for high dominance risk without computing full diagnostics.

## Usage

``` r
check_dominance_risk(a, b, threshold = 0.5, risk_level = 0.3)
```

## Arguments

- a:

  Numeric; shape parameter of the Gamma prior (\> 0).

- b:

  Numeric; rate parameter of the Gamma prior (\> 0).

- threshold:

  Numeric; dominance threshold in (0, 1) (default: 0.5).

- risk_level:

  Numeric; probability threshold for flagging risk (default: 0.3).

## Value

Logical; TRUE if P(w1 \> threshold) \> risk_level.

## Examples

``` r
check_dominance_risk(1.60, 1.22)
#> [1] TRUE
check_dominance_risk(5, 1)
#> [1] FALSE
```
