# Compare Exact Moments to NegBin Approximation

Compares exact marginal moments (via quadrature) to the Negative
Binomial approximation from the A1 method for error analysis.

## Usage

``` r
compare_to_negbin(J, a, b, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- J:

  Integer; sample size.

- a:

  Numeric; shape parameter of Gamma prior.

- b:

  Numeric; rate parameter of Gamma prior.

- M:

  Integer; number of quadrature nodes (default: 80).

## Value

A list with components:

- `exact`:

  List with exact mean and var

- `negbin`:

  List with NegBin approximation mean and var

- `abs_error`:

  Absolute errors (negbin - exact)

- `rel_error`:

  Relative errors

## Details

The NegBin approximation (A1 method from RN-03) assumes: \$\$K_J - 1 \|
\alpha \approx \text{Poisson}(\alpha \cdot c_J)\$\$

where \\c_J = \log(J)\\. With \\\alpha \sim \text{Gamma}(a, b)\\:
\$\$E\[K_J\] \approx 1 + (a/b) \cdot c_J\$\$ \$\$Var(K_J) \approx m
\cdot (1 + m/a), \quad m = (a/b) \cdot c_J\$\$

This comparison helps diagnose when the A1 approximation is insufficient
and exact A2 moment matching is needed.

## Examples

``` r
# Large approximation error for small J
compare_to_negbin(50, 1.5, 0.5)
#> $exact
#> $exact$mean
#> [1] 8.355487
#> 
#> $exact$var
#> [1] 22.76895
#> 
#> 
#> $negbin
#> $negbin$mean
#> [1] 12.73607
#> 
#> $negbin$var
#> [1] 103.5596
#> 
#> 
#> $abs_error
#> $abs_error$mean
#> [1] 4.380582
#> 
#> $abs_error$var
#> [1] 80.79066
#> 
#> 
#> $rel_error
#> $rel_error$mean
#> [1] 0.5242761
#> 
#> $rel_error$var
#> [1] 3.548282
#> 
#> 

# Error decreases with J
compare_to_negbin(300, 1.5, 0.5)
#> $exact
#> $exact$mean
#> [1] 13.52154
#> 
#> $exact$var
#> [1] 77.48541
#> 
#> 
#> $negbin
#> $negbin$mean
#> [1] 18.11135
#> 
#> $negbin$var
#> [1] 212.3102
#> 
#> 
#> $abs_error
#> $abs_error$mean
#> [1] 4.589808
#> 
#> $abs_error$var
#> [1] 134.8247
#> 
#> 
#> $rel_error
#> $rel_error$mean
#> [1] 0.3394442
#> 
#> $rel_error$var
#> [1] 1.740002
#> 
#> 
```
