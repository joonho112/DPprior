# Verify Marginal Moments Properties

Runs verification tests on the marginal moment computations.

## Usage

``` r
verify_marginal_moments(J, a, b, verbose = TRUE)
```

## Arguments

- J:

  Integer; sample size to test.

- a:

  Numeric; shape parameter to test.

- b:

  Numeric; rate parameter to test.

- verbose:

  Logical; if `TRUE`, print detailed results.

## Value

Logical; `TRUE` if all verifications pass.

## Examples

``` r
if (FALSE) { # \dontrun{
verify_marginal_moments(50, 2.0, 1.0)

} # }
```
