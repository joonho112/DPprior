# Get K_J PMF via pmf_K_marginal

Internal helper that computes the marginal PMF of K_J using
pmf_K_marginal() and Stirling numbers.

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
