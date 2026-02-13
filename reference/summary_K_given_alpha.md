# Summary of Conditional Moments

Returns a comprehensive summary of the conditional distribution of K
given alpha.

## Usage

``` r
summary_K_given_alpha(J, alpha)
```

## Arguments

- J:

  Sample size (integer \>= 1).

- alpha:

  Concentration parameter (scalar).

## Value

A list with components: J, alpha, mean, var, sd, cv, dispersion,
dmean_dalpha.

## Examples

``` r
if (FALSE) { # \dontrun{
summary_K_given_alpha(50, 2.0)

} # }
```
