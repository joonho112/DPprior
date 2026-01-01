# Induced Marginal PMF of K_J under alpha ~ Gamma(a, b)

Induced Marginal PMF of K_J under alpha ~ Gamma(a, b)

## Usage

``` r
.a2_kl_induced_pmf(J, a, b, logS, M)
```

## Arguments

- J:

  Integer; sample size.

- a, b:

  Gamma hyperparameters.

- logS:

  Matrix; log-Stirling numbers (from compute_log_stirling(J)).

- M:

  Quadrature nodes.

## Value

Numeric vector of length J; PMF on k = 1,...,J.
