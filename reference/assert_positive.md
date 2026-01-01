# Assert Positive Values

Validates that all elements of a numeric vector are strictly positive
and finite. Throws an informative error if validation fails.

## Usage

``` r
assert_positive(x, name = "x")
```

## Arguments

- x:

  Numeric vector to validate.

- name:

  Character string naming the parameter (for error messages).

## Value

Invisible `TRUE` if validation passes.

## Examples

``` r
if (FALSE) { # \dontrun{
assert_positive(c(1, 2, 3), "alpha")
assert_positive(c(1, -1), "alpha")
} # }
```
