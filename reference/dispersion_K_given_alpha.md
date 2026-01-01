# Dispersion Index for K Given Alpha

Computes the dispersion index (variance-to-mean ratio) of \\K_J \|
\alpha\\. Always \< 1 for this distribution (underdispersion).

## Usage

``` r
dispersion_K_given_alpha(J, alpha)
```

## Arguments

- J:

  Sample size (integer \>= 1).

- alpha:

  Concentration parameter (positive numeric, vectorized).

## Value

Numeric vector of dispersion indices.
