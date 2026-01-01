# Run All Module 06 Verification Tests

Comprehensive verification suite for the marginal PMF module.

## Usage

``` r
verify_pmf_marginal_all(verbose = TRUE)
```

## Arguments

- verbose:

  Logical; if TRUE, print detailed results.

## Value

Logical; TRUE if all tests pass.

## Examples

``` r
verify_pmf_marginal_all()
#> ======================================================================
#> Module 06: Marginal PMF - Full Verification Suite
#> ======================================================================
#> 
#> 
#> [Test case: J=50, a=1.50, b=0.50]
#> -------------------------------------------------- 
#> Marginal PMF Properties (J=50, a=1.50, b=0.50):
#>   Sum = 1:          PASS (sum = 1.000000000000)
#>   P(K=0) = 0:       PASS (P(K=0) = 0.00e+00)
#>   Non-negative:     PASS (min = 0.00e+00)
#>   CDF monotonic:    PASS
#>   CDF[J] = 1:       PASS (CDF[J] = 1.000000000000)
#>   Overall: PASS
#> Moments Consistency (J=50, a=1.50, b=0.50):
#>   Mean (PMF):       8.35548676
#>   Mean (quadrature):8.35548676
#>   Mean error:       7.11e-15 [PASS]
#>   Var (PMF):        22.76895012
#>   Var (quadrature): 22.76895012
#>   Var error:        4.26e-14 [PASS]
#>   Overall: PASS
#> 
#> [Test case: J=50, a=2.00, b=1.00]
#> -------------------------------------------------- 
#> Marginal PMF Properties (J=50, a=2.00, b=1.00):
#>   Sum = 1:          PASS (sum = 1.000000000000)
#>   P(K=0) = 0:       PASS (P(K=0) = 0.00e+00)
#>   Non-negative:     PASS (min = 0.00e+00)
#>   CDF monotonic:    PASS
#>   CDF[J] = 1:       PASS (CDF[J] = 1.000000000000)
#>   Overall: PASS
#> Moments Consistency (J=50, a=2.00, b=1.00):
#>   Mean (PMF):       6.63969289
#>   Mean (quadrature):6.63969289
#>   Mean error:       7.99e-15 [PASS]
#>   Var (PMF):        12.95450229
#>   Var (quadrature): 12.95450229
#>   Var error:        8.53e-14 [PASS]
#>   Overall: PASS
#> 
#> [Test case: J=100, a=1.50, b=0.50]
#> -------------------------------------------------- 
#> Marginal PMF Properties (J=100, a=1.50, b=0.50):
#>   Sum = 1:          PASS (sum = 1.000000000000)
#>   P(K=0) = 0:       PASS (P(K=0) = 0.00e+00)
#>   Non-negative:     PASS (min = 0.00e+00)
#>   CDF monotonic:    PASS
#>   CDF[J] = 1:       PASS (CDF[J] = 1.000000000000)
#>   Overall: PASS
#> Moments Consistency (J=100, a=1.50, b=0.50):
#>   Mean (PMF):       10.31191600
#>   Mean (quadrature):10.31191600
#>   Mean error:       2.84e-14 [PASS]
#>   Var (PMF):        39.26043331
#>   Var (quadrature): 39.26043331
#>   Var error:        1.42e-14 [PASS]
#>   Overall: PASS
#> 
#> [Test case: J=30, a=1.00, b=0.50]
#> -------------------------------------------------- 
#> Marginal PMF Properties (J=30, a=1.00, b=0.50):
#>   Sum = 1:          PASS (sum = 1.000000000000)
#>   P(K=0) = 0:       PASS (P(K=0) = 0.00e+00)
#>   Non-negative:     PASS (min = 0.00e+00)
#>   CDF monotonic:    PASS
#>   CDF[J] = 1:       PASS (CDF[J] = 1.000000000000)
#>   Overall: PASS
#> Moments Consistency (J=30, a=1.00, b=0.50):
#>   Mean (PMF):       5.37852414
#>   Mean (quadrature):5.37852414
#>   Mean error:       8.88e-16 [PASS]
#>   Var (PMF):        12.35202849
#>   Var (quadrature): 12.35202849
#>   Var error:        2.84e-14 [PASS]
#>   Overall: PASS
#> 
#> [Test case: J=50, a=3.00, b=1.50]
#> -------------------------------------------------- 
#> Marginal PMF Properties (J=50, a=3.00, b=1.50):
#>   Sum = 1:          PASS (sum = 1.000000000000)
#>   P(K=0) = 0:       PASS (P(K=0) = 0.00e+00)
#>   Non-negative:     PASS (min = 0.00e+00)
#>   CDF monotonic:    PASS
#>   CDF[J] = 1:       PASS (CDF[J] = 1.000000000000)
#>   Overall: PASS
#> Moments Consistency (J=50, a=3.00, b=1.50):
#>   Mean (PMF):       6.76261555
#>   Mean (quadrature):6.76261555
#>   Mean error:       7.99e-15 [PASS]
#>   Var (PMF):        10.44232545
#>   Var (quadrature): 10.44232545
#>   Var error:        1.21e-13 [PASS]
#>   Overall: PASS
#> 
#> [Convergence test]
#> -------------------------------------------------- 
#> Quadrature Convergence (J=50, a=1.50, b=0.50):
#>      M         mean          var     mean_chg      var_chg       L1_chg
#>     20   8.35549735  22.76870932     0.00e+00     0.00e+00     0.00e+00
#>     40   8.35548682  22.76894879     1.05e-05     2.39e-04     4.09e-03
#>     80   8.35548676  22.76895012     5.91e-08     1.34e-06     7.89e-05
#>    120   8.35548676  22.76895013     3.75e-11     8.48e-10     8.62e-08
#> 
#> ======================================================================
#> Overall Result: ALL TESTS PASSED
#> ====================================================================== 
```
