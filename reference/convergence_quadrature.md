# Convergence Diagnostic for Quadrature

Examines how quadrature accuracy improves with increasing number of
nodes.

## Usage

``` r
convergence_quadrature(
  f,
  a,
  b,
  M_values = c(10, 20, 50, 80, 100),
  true_value = NULL
)
```

## Arguments

- f:

  Function to integrate.

- a:

  Numeric; shape parameter of Gamma distribution.

- b:

  Numeric; rate parameter of Gamma distribution.

- M_values:

  Integer vector; number of nodes to test.

- true_value:

  Numeric; known true value (optional).

## Value

A data frame with columns: M, estimate, change, relative_change.

## Examples

``` r
if (FALSE) { # \dontrun{
# Check convergence for E[alpha]
convergence_quadrature(identity, 2.5, 1.5,
                       M_values = c(10, 20, 50, 80, 100),
                       true_value = 2.5/1.5)

} # }
```
