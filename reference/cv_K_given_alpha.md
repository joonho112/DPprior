# Coefficient of Variation for K Given Alpha

Computes the coefficient of variation (CV) of \\K_J \| \alpha\\, defined
as \\CV = SD(K) / E\[K\]\\.

## Usage

``` r
cv_K_given_alpha(J, alpha)
```

## Arguments

- J:

  Sample size (integer \>= 1).

- alpha:

  Concentration parameter (positive numeric, vectorized).

## Value

Numeric vector of coefficients of variation.
