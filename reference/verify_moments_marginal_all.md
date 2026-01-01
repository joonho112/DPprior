# Run All Module 05 Verification Tests

Comprehensive verification suite for the marginal moments module.

## Usage

``` r
verify_moments_marginal_all(verbose = TRUE)
```

## Arguments

- verbose:

  Logical; if `TRUE`, print detailed results.

## Value

Logical; `TRUE` if all tests pass.

## Examples

``` r
verify_moments_marginal_all()
#> ======================================================================
#> Module 05: Marginal Moments - Full Verification Suite
#> ======================================================================
#> 
#> Marginal Moments Verification (J=50, a=1.50, b=0.50):
#>   E[K_J] = 8.355487
#>   Var(K_J) = 22.768950
#>   SD(K_J) = 4.771682
#>   CV(K_J) = 0.571084
#> 
#>   Test 1 (Mean bounds): E[K] in [1, 50]? PASS
#>   Test 2 (Var non-negative): Var(K) >= 0? PASS
#>   Test 3 (Var inflation): Var(K)=22.7690 > Var(K|E[alpha])=5.7311? PASS
#>   Test 4 (CV positive): CV > 0? PASS
#> 
#>   Overall: PASS
#> 
#> Marginal Moments Verification (J=50, a=2.00, b=1.00):
#>   E[K_J] = 6.639693
#>   Var(K_J) = 12.954502
#>   SD(K_J) = 3.599236
#>   CV(K_J) = 0.542079
#> 
#>   Test 1 (Mean bounds): E[K] in [1, 50]? PASS
#>   Test 2 (Var non-negative): Var(K) >= 0? PASS
#>   Test 3 (Var inflation): Var(K)=12.9545 > Var(K|E[alpha])=4.5356? PASS
#>   Test 4 (CV positive): CV > 0? PASS
#> 
#>   Overall: PASS
#> 
#> Marginal Moments Verification (J=100, a=2.00, b=1.00):
#>   E[K_J] = 7.978551
#>   Var(K_J) = 20.406843
#>   SD(K_J) = 4.517393
#>   CV(K_J) = 0.566192
#> 
#>   Test 1 (Mean bounds): E[K] in [1, 100]? PASS
#>   Test 2 (Var non-negative): Var(K) >= 0? PASS
#>   Test 3 (Var inflation): Var(K)=20.4068 > Var(K|E[alpha])=5.8542? PASS
#>   Test 4 (CV positive): CV > 0? PASS
#> 
#>   Overall: PASS
#> 
#> Marginal Moments Verification (J=10, a=0.50, b=0.50):
#>   E[K_J] = 2.471937
#>   Var(K_J) = 2.905109
#>   SD(K_J) = 1.704438
#>   CV(K_J) = 0.689515
#> 
#>   Test 1 (Mean bounds): E[K] in [1, 10]? PASS
#>   Test 2 (Var non-negative): Var(K) >= 0? PASS
#>   Test 3 (Var inflation): Var(K)=2.9051 > Var(K|E[alpha])=1.3792? PASS
#>   Test 4 (CV positive): CV > 0? PASS
#> 
#>   Overall: PASS
#> 
#> Marginal Moments Verification (J=10, a=5.00, b=2.00):
#>   E[K_J] = 4.310963
#>   Var(K_J) = 2.512577
#>   SD(K_J) = 1.585111
#>   CV(K_J) = 0.367693
#> 
#>   Test 1 (Mean bounds): E[K] in [1, 10]? PASS
#>   Test 2 (Var non-negative): Var(K) >= 0? PASS
#>   Test 3 (Var inflation): Var(K)=2.5126 > Var(K|E[alpha])=1.9109? PASS
#>   Test 4 (CV positive): CV > 0? PASS
#> 
#>   Overall: PASS
#> 
#> ----------------------------------------------------------------------
#> NegBin Approximation Error Analysis:
#> ----------------------------------------------------------------------
#>   J= 50: E[K] exact=8.36, NegBin=12.74, rel_err=52.4%
#>   J=100: E[K] exact=10.31, NegBin=14.82, rel_err=43.7%
#>   J=300: E[K] exact=13.52, NegBin=18.11, rel_err=33.9%
#> 
#> ----------------------------------------------------------------------
#> Quadrature Convergence Test:
#> ----------------------------------------------------------------------
#> Quadrature Convergence (J=50, a=1.50, b=0.50):
#>    M     mean      var  mean_change   var_change
#>   20 8.355497 22.76871           NA           NA
#>   40 8.355487 22.76895 1.052945e-05 2.394687e-04
#>   60 8.355487 22.76895 5.804876e-08 1.314309e-06
#>   80 8.355487 22.76895 1.056184e-09 2.390260e-08
#>  100 8.355487 22.76895 3.567457e-11 8.072334e-10
#>  120 8.355487 22.76895 1.790568e-12 4.058620e-11
#> 
#> Mean: CONVERGED to machine precision
#> 
#> ======================================================================
#> Overall Result: ALL TESTS PASSED
#> ======================================================================
```
