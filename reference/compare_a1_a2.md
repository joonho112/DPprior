# Compare A1 vs A2 Accuracy

Compares the accuracy of A1 closed-form and A2 Newton methods.

## Usage

``` r
compare_a1_a2(J, mu_K, var_K, verbose = TRUE)
```

## Arguments

- J:

  Integer; sample size.

- mu_K:

  Numeric; target mean.

- var_K:

  Numeric; target variance.

- verbose:

  Logical; if TRUE, print comparison.

## Value

A list with A1 and A2 results and error comparison.

## Examples

``` r
compare_a1_a2(J = 50, mu_K = 5, var_K = 8)
#> A1 vs A2 Comparison (J=50, mu_K=5.00, var_K=8.00)
#> ------------------------------------------------------------ 
#>                                A1           A2
#> ------------------------------------------------------------ 
#> Shape (a)                4.000000     2.036093
#> Rate (b)                 3.912023     1.605054
#> E[K] achieved            4.461351 4.9999999992
#> Var achieved             4.783136 8.0000000076
#> Residual                 3.261649     7.60e-09
#> ------------------------------------------------------------ 
#> Improvement ratio: 429143906x
```
