# Convenience Wrapper for Marginal Moments

Returns marginal mean and variance as a named numeric vector.

## Usage

``` r
K_moments(J, a, b, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- J:

  Integer; sample size (positive integer \>= 1).

- a:

  Numeric; shape parameter of Gamma prior (\> 0).

- b:

  Numeric; rate parameter of Gamma prior (\> 0).

- M:

  Integer; number of quadrature nodes (default: 80).

## Value

Named numeric vector `c(mean = ..., var = ...)`.

## See also

[`exact_K_moments`](https://joonho112.github.io/DPprior/reference/exact_K_moments.md)
for full output

## Examples

``` r
K_moments(50, 2.0, 1.0)
#>      mean       var 
#>  6.639693 12.954502 
```
