# Assert Valid Probability

Validates that all elements of a numeric vector are valid probabilities
in the range \[0, 1\] and finite.

## Usage

``` r
assert_probability(p, name = "p")
```

## Arguments

- p:

  Numeric vector of probability values to validate.

- name:

  Character string naming the parameter (for error messages).

## Value

Invisible `TRUE` if validation passes.

## Examples

``` r
if (FALSE) { # \dontrun{
assert_probability(0.5, "p")
assert_probability(1.5, "p")
} # }
```
