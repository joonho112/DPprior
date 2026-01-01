# Conditional Second Moment of rho Given Alpha

Computes
`E(rho^2 | alpha) = (alpha + 6) / ((alpha+1)(alpha+2)(alpha+3))`.

## Usage

``` r
mean_rho_sq_given_alpha(alpha)
```

## Arguments

- alpha:

  Numeric vector; concentration parameter(s) (must be positive).

## Value

Numeric vector; `E(rho^2 | alpha)`.

## Details

Used internally for variance computation via the identity:
`Var(rho|alpha) = E(rho^2|alpha) - E(rho|alpha)^2`
