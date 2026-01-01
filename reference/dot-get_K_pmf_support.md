# Get K_J PMF with Fallback

Internal helper that tries exact_K_pmf() first, then falls back to
pmf_K_marginal() if needed.

## Usage

``` r
.get_K_pmf_support(J, a, b, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- J:

  Integer; sample size.

- a, b:

  Numeric; Gamma hyperparameters.

- M:

  Integer; quadrature nodes.

## Value

List with pmf (numeric vector) and support (integer vector).
