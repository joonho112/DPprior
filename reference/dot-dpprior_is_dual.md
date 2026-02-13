# Check if Fit is from Dual-Anchor Calibration

Checks whether a DPprior_fit object originated from DPprior_dual(). A
stricter check ensures the dual_anchor component has the expected
structure.

## Usage

``` r
.dpprior_is_dual(fit)
```

## Arguments

- fit:

  A DPprior_fit object.

## Value

Logical; TRUE if dual-anchor fit with valid structure.
