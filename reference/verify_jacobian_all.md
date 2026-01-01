# Run All Module 07 Verification Tests

Comprehensive verification suite for the score-based Jacobian module.

## Usage

``` r
verify_jacobian_all(verbose = TRUE)
```

## Arguments

- verbose:

  Logical; if TRUE, print detailed results.

## Value

Logical; TRUE if all tests pass.

## Examples

``` r
verify_jacobian_all()
#> ====================================================================== 
#> Module 07: Score-Based Jacobian - Full Verification Suite
#> ====================================================================== 
#> 
#> [Test 1] Score function zero-expectation property
#> -------------------------------------------------- 
#> Note: E[s_a] may show larger errors due to log(alpha) term.
#> 
#> Score Expectation Verification (a=2.00, b=1.00, M=200)
#>   E[s_a(alpha)] = 1.448463e-05 (should be ~= 0)
#>   E[s_b(alpha)] = -1.197109e-15 (should be ~= 0)
#> 
#> Score Expectation Verification (a=1.50, b=0.50, M=200)
#>   E[s_a(alpha)] = 1.878704e-04 (should be ~= 0)
#>   E[s_b(alpha)] = 8.740128e-16 (should be ~= 0)
#> 
#> Score Expectation Verification (a=3.00, b=1.50, M=200)
#>   E[s_a(alpha)] = 1.380902e-07 (should be ~= 0)
#>   E[s_b(alpha)] = -1.964492e-16 (should be ~= 0)
#> 
#> [Test 2] Jacobian vs finite differences (M=200)
#> -------------------------------------------------- 
#> 
#> Jacobian Verification (J=30, a=1.50, b=0.50, M=200)
#> ------------------------------------------------------------ 
#> 
#> Analytic Jacobian (score-based):
#>   dM1/da =   2.92957577  dM1/db =  -7.58902093
#>   dV/da  =   2.31935647  dV/db  = -22.12106296
#> 
#> Numeric Jacobian (finite diff):
#>   dM1/da =   2.92939276  dM1/db =  -7.58902093
#>   dV/da  =   2.32173650  dV/db  = -22.12106296
#> 
#> Relative Errors:
#>   dM1/da: 6.25e-05  dM1/db: 3.18e-11
#>   dV/da:  1.03e-03  dV/db:  2.78e-12
#> 
#> Max Relative Error: 1.03e-03 [PASS]
#> 
#> Jacobian Verification (J=50, a=2.00, b=1.00, M=200)
#> ------------------------------------------------------------ 
#> 
#> Analytic Jacobian (score-based):
#>   dM1/da =   2.24553449  dM1/db =  -4.13558517
#>   dV/da  =   2.94445788  dV/db  = -13.03822850
#> 
#> Numeric Jacobian (finite diff):
#>   dM1/da =   2.24552030  dM1/db =  -4.13558517
#>   dV/da  =   2.94463271  dV/db  = -13.03822849
#> 
#> Relative Errors:
#>   dM1/da: 6.32e-06  dM1/db: 1.29e-11
#>   dV/da:  5.94e-05  dV/db:  6.13e-10
#> 
#> Max Relative Error: 5.94e-05 [PASS]
#> 
#> Jacobian Verification (J=100, a=2.00, b=1.00, M=200)
#> ------------------------------------------------------------ 
#> 
#> Analytic Jacobian (score-based):
#>   dM1/da =   2.89660795  dM1/db =  -5.42013736
#>   dV/da  =   5.58491235  dV/db  = -23.84210671
#> 
#> Numeric Jacobian (finite diff):
#>   dM1/da =   2.89659380  dM1/db =  -5.42013736
#>   dV/da  =   5.58512458  dV/db  = -23.84210670
#> 
#> Relative Errors:
#>   dM1/da: 4.88e-06  dM1/db: 2.96e-11
#>   dV/da:  3.80e-05  dV/db:  2.25e-10
#> 
#> Max Relative Error: 3.80e-05 [PASS]
#> 
#> Jacobian Verification (J=50, a=3.00, b=1.50, M=200)
#> ------------------------------------------------------------ 
#> 
#> Analytic Jacobian (score-based):
#>   dM1/da =   1.50495388  dM1/db =  -2.83988706
#>   dV/da  =   1.57224901  dV/db  =  -6.67421114
#> 
#> Numeric Jacobian (finite diff):
#>   dM1/da =   1.50495374  dM1/db =  -2.83988706
#>   dV/da  =   1.57225071  dV/db  =  -6.67421114
#> 
#> Relative Errors:
#>   dM1/da: 9.19e-08  dM1/db: 2.22e-10
#>   dV/da:  1.08e-06  dV/db:  6.22e-10
#> 
#> Max Relative Error: 1.08e-06 [PASS]
#> 
#> [Test 3] Moments consistency check
#> -------------------------------------------------- 
#> Mean (moments_with_jacobian): 6.6396928911
#> Mean (exact_K_moments):       6.6396928911
#> Difference: 0.00e+00
#> Var (moments_with_jacobian):  12.9545022869
#> Var (exact_K_moments):        12.9545022869
#> Difference: 0.00e+00
#> Status: PASS
#> 
#> [Test 4] Jacobian non-singularity check
#> -------------------------------------------------- 
#> J=30, a=1.5, b=0.5: det=-47.2655, cond=1.34e+01 [PASS]
#> J=50, a=2.0, b=1.0: det=-17.1053, cond=1.43e+01 [PASS]
#> J=100, a=2.0, b=1.0: det=-38.7977, cond=1.99e+01 [PASS]
#> 
#> [Test 5] Newton convergence test
#> -------------------------------------------------- 
#> Target: E[K]=5.00, Var(K)=8.00
#> Start:  a=2.0000, b=1.0000
#> 
#> Iter |    a      |    b      |   E[K]    |  Var(K)   | ||F||
#> ---------------------------------------------------------------------- 
#>    0 |  2.000000 |  1.000000 |  6.639693 | 12.954502 | 5.22e+00
#>    1 |  1.948018 |  1.368259 |  5.362738 |  9.247544 | 1.30e+00
#>    2 |  2.014501 |  1.566695 |  5.040031 |  8.159377 | 1.64e-01
#>    3 |  2.035417 |  1.604054 |  5.000837 |  8.003804 | 3.89e-03
#>    4 |  2.036092 |  1.605053 |  5.000000 |  8.000002 | 2.54e-06
#>    5 |  2.036093 |  1.605054 |  5.000000 |  8.000000 | 5.99e-11
#> 
#> Converged in 5 iterations!
#> 
#> ======================================================================
#> Overall Result: ALL TESTS PASSED
#> ====================================================================== 
```
