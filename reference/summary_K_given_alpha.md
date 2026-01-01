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
summary_K_given_alpha(50, 2.0)
#> $J
#> [1] 50
#> 
#> $alpha
#> [1] 2
#> 
#> $mean
#> [1] 7.037626
#> 
#> $var
#> [1] 4.535558
#> 
#> $sd
#> [1] 2.129685
#> 
#> $cv
#> [1] 0.3026141
#> 
#> $dispersion
#> [1] 0.6444726
#> 
#> $dmean_dalpha
#> [1] 2.267779
#> 
```
