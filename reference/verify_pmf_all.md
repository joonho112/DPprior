# Run All PMF Verifications

Runs comprehensive verification tests for the conditional PMF module.

## Usage

``` r
verify_pmf_all(
  J_values = c(10, 50, 100),
  alpha_values = c(0.5, 1, 2, 5),
  verbose = TRUE
)
```

## Arguments

- J_values:

  Integer vector; sample sizes to test.

- alpha_values:

  Numeric vector; concentration parameters to test.

- verbose:

  Logical; if `TRUE`, print detailed results.

## Value

Logical; `TRUE` if all verifications pass.

## Examples

``` r
if (FALSE) { # \dontrun{
verify_pmf_all()

} # }
```
