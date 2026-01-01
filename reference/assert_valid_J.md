# Assert Valid Sample Size J

Validates that J is a positive integer within the supported range.

## Usage

``` r
assert_valid_J(J)
```

## Arguments

- J:

  Sample size to validate.

## Value

Invisible `TRUE` if validation passes.

## Examples

``` r
if (FALSE) { # \dontrun{
assert_valid_J(50)
assert_valid_J(0)
assert_valid_J(1000)
} # }
```
