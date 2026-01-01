# Conditional Mean of K Given Alpha (Enhanced)

Computes \\E\[K_J \| \alpha\]\\ with proper handling of very small alpha
values.

## Usage

``` r
mean_K_given_alpha_safe(J, alpha)
```

## Arguments

- J:

  Integer; sample size.

- alpha:

  Numeric vector; concentration parameter values.

## Value

Numeric vector of conditional means.

## Details

For very small alpha (\< 1e-12), returns 1.0 as the limiting value. This
avoids numerical issues with digamma at near-zero arguments.
