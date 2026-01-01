# Verify Marginal Moments Properties

Runs verification tests on the marginal moment computations.

## Usage

``` r
verify_marginal_moments(J, a, b, verbose = TRUE)
```

## Arguments

- J:

  Integer; sample size to test.

- a:

  Numeric; shape parameter to test.

- b:

  Numeric; rate parameter to test.

- verbose:

  Logical; if `TRUE`, print detailed results.

## Value

Logical; `TRUE` if all verifications pass.

## Examples

``` r
verify_marginal_moments(50, 2.0, 1.0)
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
```
