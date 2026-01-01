# Run All PMF Verifications

Runs comprehensive verification tests for the conditional PMF module.

## Usage

``` r
verify_pmf_all(
  J_values = c(10, 50, 100),
  alpha_values = c(0.5, 1, 2, 5),
  verbose = TRUE
)
```

## Arguments

- J_values:

  Integer vector; sample sizes to test.

- alpha_values:

  Numeric vector; concentration parameters to test.

- verbose:

  Logical; if `TRUE`, print detailed results.

## Value

Logical; `TRUE` if all verifications pass.

## Examples

``` r
verify_pmf_all()
#> 
#> ============================================================
#> Test case: J=10, alpha=0.50
#> ============================================================
#> PMF Normalization (J=10, alpha=0.50):
#>   Sum = 1.000000000000 (expected 1.0) [PASS]
#> Zero Probability (J=10, alpha=0.50):
#>   P(K=0) = 0.00e+00 (expected 0) [PASS]
#> PMF-Moments Verification (J=10, alpha=0.50):
#>   Mean: PMF=2.13325553, digamma=2.13325553, error=4.44e-16 [PASS]
#>   Var:  PMF=0.92453422, trigamma=0.92453422, error=2.22e-16 [PASS]
#>   Overall: PASS
#> CDF Properties (J=10, alpha=0.50):
#>   Monotonicity: min(diff) = 1.53e-09 [PASS]
#>   CDF[J] = 1.000000000000 (expected 1.0) [PASS]
#>   Overall: PASS
#> 
#> ============================================================
#> Test case: J=10, alpha=1.00
#> ============================================================
#> PMF Normalization (J=10, alpha=1.00):
#>   Sum = 1.000000000000 (expected 1.0) [PASS]
#> Zero Probability (J=10, alpha=1.00):
#>   P(K=0) = 0.00e+00 (expected 0) [PASS]
#> PMF-Moments Verification (J=10, alpha=1.00):
#>   Mean: PMF=2.92896825, digamma=2.92896825, error=4.44e-16 [PASS]
#>   Var:  PMF=1.37920052, trigamma=1.37920052, error=1.33e-15 [PASS]
#>   Overall: PASS
#> CDF Properties (J=10, alpha=1.00):
#>   Monotonicity: min(diff) = 2.76e-07 [PASS]
#>   CDF[J] = 1.000000000000 (expected 1.0) [PASS]
#>   Overall: PASS
#> 
#> ============================================================
#> Test case: J=10, alpha=2.00
#> ============================================================
#> PMF Normalization (J=10, alpha=2.00):
#>   Sum = 1.000000000000 (expected 1.0) [PASS]
#> Zero Probability (J=10, alpha=2.00):
#>   P(K=0) = 0.00e+00 (expected 0) [PASS]
#> PMF-Moments Verification (J=10, alpha=2.00):
#>   Mean: PMF=4.03975469, digamma=4.03975469, error=1.78e-15 [PASS]
#>   Var:  PMF=1.80762591, trigamma=1.80762591, error=1.33e-15 [PASS]
#>   Overall: PASS
#> CDF Properties (J=10, alpha=2.00):
#>   Monotonicity: min(diff) = 2.57e-05 [PASS]
#>   CDF[J] = 1.000000000000 (expected 1.0) [PASS]
#>   Overall: PASS
#> 
#> ============================================================
#> Test case: J=10, alpha=5.00
#> ============================================================
#> PMF Normalization (J=10, alpha=5.00):
#>   Sum = 1.000000000000 (expected 1.0) [PASS]
#> Zero Probability (J=10, alpha=5.00):
#>   P(K=0) = 0.00e+00 (expected 0) [PASS]
#> PMF-Moments Verification (J=10, alpha=5.00):
#>   Mean: PMF=5.84114497, digamma=5.84114497, error=1.78e-15 [PASS]
#>   Var:  PMF=2.03152677, trigamma=2.03152677, error=3.11e-15 [PASS]
#>   Overall: PASS
#> CDF Properties (J=10, alpha=5.00):
#>   Monotonicity: min(diff) = 5.00e-04 [PASS]
#>   CDF[J] = 1.000000000000 (expected 1.0) [PASS]
#>   Overall: PASS
#> 
#> ============================================================
#> Test case: J=50, alpha=0.50
#> ============================================================
#> PMF Normalization (J=50, alpha=0.50):
#>   Sum = 1.000000000000 (expected 1.0) [PASS]
#> Zero Probability (J=50, alpha=0.50):
#>   P(K=0) = 0.00e+00 (expected 0) [PASS]
#> PMF-Moments Verification (J=50, alpha=0.50):
#>   Mean: PMF=2.93777485, digamma=2.93777485, error=6.22e-15 [PASS]
#>   Var:  PMF=1.70907413, trigamma=1.70907413, error=5.11e-15 [PASS]
#>   Overall: PASS
#> CDF Properties (J=50, alpha=0.50):
#>   Monotonicity: min(diff) = 0.00e+00 [PASS]
#>   CDF[J] = 1.000000000000 (expected 1.0) [PASS]
#>   Overall: PASS
#> 
#> ============================================================
#> Test case: J=50, alpha=1.00
#> ============================================================
#> PMF Normalization (J=50, alpha=1.00):
#>   Sum = 1.000000000000 (expected 1.0) [PASS]
#> Zero Probability (J=50, alpha=1.00):
#>   P(K=0) = 0.00e+00 (expected 0) [PASS]
#> PMF-Moments Verification (J=50, alpha=1.00):
#>   Mean: PMF=4.49920534, digamma=4.49920534, error=6.22e-15 [PASS]
#>   Var:  PMF=2.87407260, trigamma=2.87407260, error=6.22e-15 [PASS]
#>   Overall: PASS
#> CDF Properties (J=50, alpha=1.00):
#>   Monotonicity: min(diff) = 0.00e+00 [PASS]
#>   CDF[J] = 1.000000000000 (expected 1.0) [PASS]
#>   Overall: PASS
#> 
#> ============================================================
#> Test case: J=50, alpha=2.00
#> ============================================================
#> PMF Normalization (J=50, alpha=2.00):
#>   Sum = 1.000000000000 (expected 1.0) [PASS]
#> Zero Probability (J=50, alpha=2.00):
#>   P(K=0) = 0.00e+00 (expected 0) [PASS]
#> PMF-Moments Verification (J=50, alpha=2.00):
#>   Mean: PMF=7.03762636, digamma=7.03762636, error=2.04e-14 [PASS]
#>   Var:  PMF=4.53555756, trigamma=4.53555756, error=8.88e-14 [PASS]
#>   Overall: PASS
#> CDF Properties (J=50, alpha=2.00):
#>   Monotonicity: min(diff) = 0.00e+00 [PASS]
#>   CDF[J] = 1.000000000000 (expected 1.0) [PASS]
#>   Overall: PASS
#> 
#> ============================================================
#> Test case: J=50, alpha=5.00
#> ============================================================
#> PMF Normalization (J=50, alpha=5.00):
#>   Sum = 1.000000000000 (expected 1.0) [PASS]
#> Zero Probability (J=50, alpha=5.00):
#>   P(K=0) = 0.00e+00 (expected 0) [PASS]
#> PMF-Moments Verification (J=50, alpha=5.00):
#>   Mean: PMF=12.46048530, digamma=12.46048530, error=2.66e-14 [PASS]
#>   Var:  PMF=7.38611414, trigamma=7.38611414, error=1.21e-13 [PASS]
#>   Overall: PASS
#> CDF Properties (J=50, alpha=5.00):
#>   Monotonicity: min(diff) = 0.00e+00 [PASS]
#>   CDF[J] = 1.000000000000 (expected 1.0) [PASS]
#>   Overall: PASS
#> 
#> ============================================================
#> Test case: J=100, alpha=0.50
#> ============================================================
#> PMF Normalization (J=100, alpha=0.50):
#>   Sum = 1.000000000000 (expected 1.0) [PASS]
#> Zero Probability (J=100, alpha=0.50):
#>   P(K=0) = 0.00e+00 (expected 0) [PASS]
#> PMF-Moments Verification (J=100, alpha=0.50):
#>   Mean: PMF=3.28434219, digamma=3.28434219, error=2.84e-14 [PASS]
#>   Var:  PMF=2.05314162, trigamma=2.05314162, error=5.51e-14 [PASS]
#>   Overall: PASS
#> CDF Properties (J=100, alpha=0.50):
#>   Monotonicity: min(diff) = 0.00e+00 [PASS]
#>   CDF[J] = 1.000000000000 (expected 1.0) [PASS]
#>   Overall: PASS
#> 
#> ============================================================
#> Test case: J=100, alpha=1.00
#> ============================================================
#> PMF Normalization (J=100, alpha=1.00):
#>   Sum = 1.000000000000 (expected 1.0) [PASS]
#> Zero Probability (J=100, alpha=1.00):
#>   P(K=0) = 0.00e+00 (expected 0) [PASS]
#> PMF-Moments Verification (J=100, alpha=1.00):
#>   Mean: PMF=5.18737752, digamma=5.18737752, error=4.44e-14 [PASS]
#>   Var:  PMF=3.55239362, trigamma=3.55239362, error=1.87e-14 [PASS]
#>   Overall: PASS
#> CDF Properties (J=100, alpha=1.00):
#>   Monotonicity: min(diff) = 0.00e+00 [PASS]
#>   CDF[J] = 1.000000000000 (expected 1.0) [PASS]
#>   Overall: PASS
#> 
#> ============================================================
#> Test case: J=100, alpha=2.00
#> ============================================================
#> PMF Normalization (J=100, alpha=2.00):
#>   Sum = 1.000000000000 (expected 1.0) [PASS]
#> Zero Probability (J=100, alpha=2.00):
#>   P(K=0) = 0.00e+00 (expected 0) [PASS]
#> PMF-Moments Verification (J=100, alpha=2.00):
#>   Mean: PMF=8.39455702, digamma=8.39455702, error=2.49e-14 [PASS]
#>   Var:  PMF=5.85422930, trigamma=5.85422930, error=8.53e-14 [PASS]
#>   Overall: PASS
#> CDF Properties (J=100, alpha=2.00):
#>   Monotonicity: min(diff) = 0.00e+00 [PASS]
#>   CDF[J] = 1.000000000000 (expected 1.0) [PASS]
#>   Overall: PASS
#> 
#> ============================================================
#> Test case: J=100, alpha=5.00
#> ============================================================
#> PMF Normalization (J=100, alpha=5.00):
#>   Sum = 1.000000000000 (expected 1.0) [PASS]
#> Zero Probability (J=100, alpha=5.00):
#>   P(K=0) = 0.00e+00 (expected 0) [PASS]
#> PMF-Moments Verification (J=100, alpha=5.00):
#>   Mean: PMF=15.71536609, digamma=15.71536609, error=6.04e-14 [PASS]
#>   Var:  PMF=10.42152482, trigamma=10.42152482, error=2.38e-13 [PASS]
#>   Overall: PASS
#> CDF Properties (J=100, alpha=5.00):
#>   Monotonicity: min(diff) = 0.00e+00 [PASS]
#>   CDF[J] = 1.000000000000 (expected 1.0) [PASS]
#>   Overall: PASS
#> 
#> ============================================================
#> Overall: ALL TESTS PASSED
#> ============================================================
```
