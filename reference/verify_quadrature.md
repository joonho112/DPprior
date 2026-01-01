# Verify Quadrature Accuracy Against Known Gamma Moments

Validates quadrature implementation by comparing computed expectations
against known closed-form Gamma distribution moments.

## Usage

``` r
verify_quadrature(a, b, M, tol = 1e-10, verbose = TRUE)
```

## Arguments

- a:

  Numeric; shape parameter of Gamma distribution.

- b:

  Numeric; rate parameter of Gamma distribution.

- M:

  Integer; number of quadrature nodes.

- tol:

  Numeric; tolerance for verification (default: 1e-10).

- verbose:

  Logical; if `TRUE`, print detailed results.

## Value

Logical; `TRUE` if all moments match within tolerance.

## Details

For \\\alpha \sim \text{Gamma}(a, b)\\:

- \\E\[\alpha\] = a/b\\

- \\E\[\alpha^2\] = a(a+1)/b^2\\

- \\Var(\alpha) = a/b^2\\

## Examples

``` r
# Should return TRUE
verify_quadrature(2.5, 1.5, M = 80)
#> Quadrature verification (a=2.50, b=1.50, M=80):
#>   E[alpha]:   true=1.666667, quad=1.6666666667, error=4.44e-16 [PASS]
#>   E[alpha^2]:  true=3.888889, quad=3.8888888889, error=0.00e+00 [PASS]
#>   Var(alpha): true=1.111111, quad=1.1111111111, error=1.33e-15 [PASS]
#>   Overall: PASS

# More challenging case
verify_quadrature(0.5, 2.0, M = 100, verbose = TRUE)
#> Quadrature verification (a=0.50, b=2.00, M=100):
#>   E[alpha]:   true=0.250000, quad=0.2500000000, error=1.94e-16 [PASS]
#>   E[alpha^2]:  true=0.187500, quad=0.1875000000, error=6.11e-16 [PASS]
#>   Var(alpha): true=0.125000, quad=0.1250000000, error=5.13e-16 [PASS]
#>   Overall: PASS
```
