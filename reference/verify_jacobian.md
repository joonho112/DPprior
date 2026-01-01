# Verify Jacobian Against Finite Differences

Compares the analytically computed Jacobian (via score identities)
against numerical finite differences to validate the implementation.

## Usage

``` r
verify_jacobian(
  J,
  a,
  b,
  eps = 1e-06,
  M = .QUAD_NODES_VERIFICATION,
  verbose = TRUE
)
```

## Arguments

- J:

  Integer; sample size.

- a:

  Numeric; shape parameter.

- b:

  Numeric; rate parameter.

- eps:

  Numeric; step size for finite differences (default: 1e-6).

- M:

  Integer; number of quadrature nodes (default: 200 for verification).

- verbose:

  Logical; if TRUE, print detailed comparison.

## Value

A named list with components:

- `analytic`:

  The analytically computed Jacobian

- `numeric`:

  The numerically computed Jacobian (finite differences)

- `abs_error`:

  Matrix of absolute errors

- `rel_error`:

  Matrix of relative errors

- `max_rel_error`:

  Maximum relative error across all entries

- `pass`:

  Logical; TRUE if max relative error \< 0.01

## Details

Uses central finite differences: \$\$\frac{\partial f}{\partial a}
\approx \frac{f(a+\epsilon) - f(a-\epsilon)}{2\epsilon}\$\$

**Important:** This is a SECONDARY verification method because both the
analytical Jacobian and the finite differences use the same quadrature
layer. For independent verification, compare against adaptive
integration (scipy.integrate.quad in Python).

## Examples

``` r
# Verify Jacobian for a specific case
result <- verify_jacobian(J = 50, a = 2.0, b = 1.0, verbose = TRUE)
#> Jacobian Verification (J=50, a=2.00, b=1.00, M=200)
#> ------------------------------------------------------------ 
#> 
#> Analytic Jacobian (score-based):
#>   dM1/da =   2.24553449  dM1/db =  -4.13558517
#>   dV/da  =   2.94445788  dV/db  = -13.03822850
#> 
#> Numeric Jacobian (finite diff):
#>   dM1/da =   2.24552030  dM1/db =  -4.13558517
#>   dV/da  =   2.94463271  dV/db  = -13.03822849
#> 
#> Relative Errors:
#>   dM1/da: 6.32e-06  dM1/db: 1.29e-11
#>   dV/da:  5.94e-05  dV/db:  6.13e-10
#> 
#> Max Relative Error: 5.94e-05 [PASS]
```
