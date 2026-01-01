# Test Newton Convergence Using the Jacobian

Verifies that the Jacobian enables fast Newton convergence for moment
matching.

## Usage

``` r
test_newton_convergence(
  J,
  mu_target,
  var_target,
  a0 = 2,
  b0 = 1,
  max_iter = 15L,
  tol = 1e-08,
  verbose = TRUE
)
```

## Arguments

- J:

  Integer; sample size.

- mu_target:

  Numeric; target mean.

- var_target:

  Numeric; target variance.

- a0:

  Numeric; initial shape parameter.

- b0:

  Numeric; initial rate parameter.

- max_iter:

  Integer; maximum iterations.

- tol:

  Numeric; convergence tolerance.

- verbose:

  Logical; if TRUE, print iteration history.

## Value

Logical; TRUE if Newton converges.
