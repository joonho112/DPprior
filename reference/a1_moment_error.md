# A1 Approximation Moment Errors

Computes the discrepancy between the A1 (shifted NegBin) approximation
and exact marginal moments of \\K_J\\.

## Usage

``` r
a1_moment_error(J, a, b, cJ = log(J), M = .QUAD_NODES_DEFAULT)
```

## Arguments

- J:

  Integer; sample size.

- a, b:

  Numeric; Gamma hyperparameters (shape, rate).

- cJ:

  Numeric; scaling constant (default: log(J)).

- M:

  Integer; number of quadrature nodes (default: 80).

## Value

A list with components:

- exact_mean, exact_var:

  Exact moments via Gauss-Laguerre quadrature

- a1_mean, a1_var:

  A1 (shifted NegBin) approximation moments

- error_mean_abs, error_var_abs:

  Absolute errors

- error_mean_rel, error_var_rel:

  Relative errors (percentage)

## Details

The A1 approximation models \\K_J \approx 1 + \text{NegBin}(a, p_J)\\
where \\p_J = b / (b + c_J)\\.

The NegBin(a, p) moments are:

- Mean: \\a(1-p)/p\\

- Variance: \\a(1-p)/p^2\\

## See also

[`exact_K_moments`](https://joonho112.github.io/DPprior/reference/exact_K_moments.md),
[`DPprior_error_bounds`](https://joonho112.github.io/DPprior/reference/DPprior_error_bounds.md)

## Examples

``` r
# Moment errors for J=50, Gamma(2, 1) prior
errors <- a1_moment_error(J = 50, a = 2, b = 1)
print(errors)
#> $exact_mean
#> [1] 6.639693
#> 
#> $exact_var
#> [1] 12.9545
#> 
#> $a1_mean
#> [1] 8.824046
#> 
#> $a1_var
#> [1] 38.43189
#> 
#> $error_mean_abs
#> [1] 2.184353
#> 
#> $error_var_abs
#> [1] 25.47739
#> 
#> $error_mean_rel
#> [1] 32.89841
#> 
#> $error_var_rel
#> [1] 196.6682
#> 

# Compare A1 accuracy at different J values
sapply(c(25, 50, 100, 200), function(J) {
  err <- a1_moment_error(J, a = 2, b = 1)
  c(mean_err = err$error_mean_rel, var_err = err$error_var_rel)
})
#>               [,1]      [,2]      [,3]      [,4]
#> mean_err  39.18536  32.89841  27.97236  24.15421
#> var_err  263.32276 196.66824 152.98144 123.40937

# Compare different scaling constants
a1_moment_error(J = 50, a = 2, b = 1, cJ = log(50))
#> $exact_mean
#> [1] 6.639693
#> 
#> $exact_var
#> [1] 12.9545
#> 
#> $a1_mean
#> [1] 8.824046
#> 
#> $a1_var
#> [1] 38.43189
#> 
#> $error_mean_abs
#> [1] 2.184353
#> 
#> $error_var_abs
#> [1] 25.47739
#> 
#> $error_mean_rel
#> [1] 32.89841
#> 
#> $error_var_rel
#> [1] 196.6682
#> 
a1_moment_error(J = 50, a = 2, b = 1, cJ = digamma(50) + 0.5772)
#> $exact_mean
#> [1] 6.639693
#> 
#> $exact_var
#> [1] 12.9545
#> 
#> $a1_mean
#> [1] 9.958379
#> 
#> $a1_var
#> [1] 49.08466
#> 
#> $error_mean_abs
#> [1] 3.318686
#> 
#> $error_var_abs
#> [1] 36.13016
#> 
#> $error_mean_rel
#> [1] 49.98253
#> 
#> $error_var_rel
#> [1] 278.9004
#> 
```
