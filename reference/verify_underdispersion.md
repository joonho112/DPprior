# Verify Underdispersion Inequality

Verifies that \\0 \< Var(K_J \| \alpha) \< E\[K_J \| \alpha\]\\ for all
specified values of \\\alpha\\.

## Usage

``` r
verify_underdispersion(
  J,
  alpha_values = c(0.1, 0.5, 1, 2, 5, 10),
  verbose = TRUE
)
```

## Arguments

- J:

  Sample size (integer \>= 1).

- alpha_values:

  Numeric vector of alpha values to test (default: c(0.1, 0.5, 1, 2, 5,
  10)).

- verbose:

  Logical; if TRUE, print results.

## Value

Logical; TRUE if inequality holds for all alpha values.

## Examples

``` r
verify_underdispersion(50)
#> Underdispersion verification (J=50):
#>   alpha= 0.10: E[K]=  1.4328, Var(K)=  0.4186, D=0.2922 [PASS]
#>   alpha= 0.50: E[K]=  2.9378, Var(K)=  1.7091, D=0.5818 [PASS]
#>   alpha= 1.00: E[K]=  4.4992, Var(K)=  2.8741, D=0.6388 [PASS]
#>   alpha= 2.00: E[K]=  7.0376, Var(K)=  4.5356, D=0.6445 [PASS]
#>   alpha= 5.00: E[K]= 12.4605, Var(K)=  7.3861, D=0.5928 [PASS]
#>   alpha=10.00: E[K]= 18.3424, Var(K)=  9.5064, D=0.5183 [PASS]
verify_underdispersion(50, c(0.1, 1, 10), verbose = TRUE)
#> Underdispersion verification (J=50):
#>   alpha= 0.10: E[K]=  1.4328, Var(K)=  0.4186, D=0.2922 [PASS]
#>   alpha= 1.00: E[K]=  4.4992, Var(K)=  2.8741, D=0.6388 [PASS]
#>   alpha=10.00: E[K]= 18.3424, Var(K)=  9.5064, D=0.5183 [PASS]
```
