# Verify Variance Decomposition

Verifies the law of total variance decomposition for Var(rho).

## Usage

``` r
verify_rho_variance_decomposition(
  a,
  b,
  tol = 1e-10,
  M = .QUAD_NODES_DEFAULT,
  verbose = FALSE
)
```

## Arguments

- a:

  Numeric; shape parameter of the Gamma prior on alpha (a \> 0).

- b:

  Numeric; rate parameter of the Gamma prior on alpha (b \> 0).

- tol:

  Numeric; tolerance for comparison. Default is 1e-10.

- M:

  Integer; number of quadrature nodes. Default is 80.

- verbose:

  Logical; if TRUE, print detailed results.

## Value

Logical; TRUE if decomposition holds within tolerance.
