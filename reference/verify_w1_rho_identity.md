# Verify the Identity E(w1) = E(rho)

Checks the mean identity `E(w1 | a, b) = E(rho | a, b)`, which follows
from `E(rho | alpha) = E(w1 | alpha) = 1/(1+alpha)`.

## Usage

``` r
verify_w1_rho_identity(a, b, tol = 1e-10, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- a:

  Numeric; shape parameter of the Gamma prior on alpha (a \> 0).

- b:

  Numeric; rate parameter of the Gamma prior on alpha (b \> 0).

- tol:

  Numeric; absolute tolerance. Default is 1e-10.

- M:

  Integer; number of quadrature nodes. Default is 80.

## Value

Logical; TRUE if the identity holds within tolerance.

## Examples

``` r
if (FALSE) { # \dontrun{
verify_w1_rho_identity(2, 1)
verify_w1_rho_identity(1.6, 1.22)
} # }
```
