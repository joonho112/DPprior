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
if (FALSE) { # \dontrun{
# Should return TRUE
verify_quadrature(2.5, 1.5, M = 80)

# More challenging case
verify_quadrature(0.5, 2.0, M = 100, verbose = TRUE)

} # }
```
