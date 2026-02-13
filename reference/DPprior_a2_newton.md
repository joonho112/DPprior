# A2-MN Exact-Moment Newton Solver

Finds Gamma(a, b) hyperprior parameters that exactly match target
moments for the number of clusters K_J under a Dirichlet process prior.

## Usage

``` r
DPprior_a2_newton(
  J,
  mu_K,
  var_K,
  a0 = NULL,
  b0 = NULL,
  tol_F = .TOL_NEWTON,
  tol_step = 1e-10,
  max_iter = 20L,
  damping = TRUE,
  use_fallback = TRUE,
  M = .QUAD_NODES_DEFAULT,
  verbose = FALSE
)
```

## Arguments

- J:

  Integer; sample size (number of observations/sites). Must be \>= 2.

- mu_K:

  Numeric; target prior mean \\E\[K_J\]\\. Must satisfy \\1 \< \mu_K \<
  J\\.

- var_K:

  Numeric; target prior variance \\\mathrm{Var}(K_J)\\. Must be
  positive.

- a0:

  Numeric or NULL; initial shape parameter. If NULL, computed via
  [`DPprior_a1`](https://joonho112.github.io/DPprior/reference/DPprior_a1.md).

- b0:

  Numeric or NULL; initial rate parameter. If NULL, computed via
  [`DPprior_a1`](https://joonho112.github.io/DPprior/reference/DPprior_a1.md).

- tol_F:

  Numeric; stopping tolerance for the residual norm \|\|F\|\| =
  sqrt((M1 - mu_K)^2 + (V - var_K)^2). Default: 1e-8.

- tol_step:

  Numeric; stopping tolerance for Newton step size. Default: 1e-10.

- max_iter:

  Integer; maximum Newton iterations. Default: 20.

- damping:

  Logical; if TRUE, use backtracking line search for damped Newton
  updates. Default: TRUE.

- use_fallback:

  Logical; if TRUE, use Nelder-Mead fallback when Newton fails to
  converge. Default: TRUE.

- M:

  Integer; number of quadrature nodes for moment computation. Default:
  80.

- verbose:

  Logical; if TRUE, print iteration progress. Default: FALSE.

## Value

A `DPprior_fit` object (S3 class) with components:

- `a`:

  Numeric; optimal shape parameter

- `b`:

  Numeric; optimal rate parameter

- `J`:

  Integer; sample size

- `target`:

  List with `mu_K`, `var_K`, and `type = "moments"`

- `method`:

  Character; "A2-MN" or "A2-MN+NM" if fallback was used

- `status`:

  Character; convergence status ("success", "stagnated", "max_iter")

- `converged`:

  Logical; whether algorithm converged to target tolerance

- `iterations`:

  Integer; number of iterations

- `termination`:

  Character; reason for termination ("residual", "step", "max_iter",
  "nelder_mead")

- `fit`:

  List with achieved `mu_K`, `var_K`, and `residual`

- `diagnostics`:

  List with diagnostic information

- `trace`:

  Data frame with iteration history

## Details

This implements TSMM Stage 2 (A2-MN) from Lee (2026). The A2-MN
algorithm uses Newton's method in log-scale to ensure positivity of the
Gamma parameters. The Jacobian is computed exactly using score function
identities (Lee, 2026, Section 3.2, Corollary 1), avoiding finite
difference approximations.

**Algorithm Steps:**

1.  Initialize: \\(a_0, b_0)\\ from A1 closed-form or user-provided

2.  Log-parameterize: \\\eta = (\log a, \log b)\\

3.  For each iteration:

    - Compute moments \\(M_1, V)\\ and Jacobian \\J_F\\

    - Compute residual \\F = (M_1 - \mu_K, V - \sigma^2_K)\\

    - Transform Jacobian to log-scale: \\J\_{\log} = J_F \cdot
      \text{diag}(a, b)\\

    - Newton step: \\\Delta = -J\_{\log}^{-1} F\\

    - Backtracking line search (if damping enabled)

    - Update: \\\eta \leftarrow \eta + \lambda \Delta\\

4.  Return \\(a, b) = \exp(\eta)\\

**Convergence Behavior:**

- For typical targets, converges in 3-8 iterations

- For targets requiring very small \\a\\ (quasi-improper priors),
  convergence may be slower; the Nelder-Mead fallback handles these
  cases

- Machine-precision accuracy (residual \< 1e-10) is typically achieved

**Termination Conditions:**

- `residual`: \|\|F\|\| \< tol_F (success)

- `step`: \|\|delta\|\| \< tol_step (stagnation - may not have
  converged)

- `max_iter`: maximum iterations reached

- `nelder_mead`: Nelder-Mead fallback succeeded

**Error Handling:**

- If Jacobian becomes singular, falls back to gradient descent

- If Newton fails after `max_iter`, uses Nelder-Mead if enabled

- Warns for quasi-improper priors (\\a \< 0.1\\)

## References

Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet
Process Mixtures. *arXiv preprint* arXiv:2602.06301.

## See also

[`DPprior_a1`](https://joonho112.github.io/DPprior/reference/DPprior_a1.md)
for closed-form initialization,
[`moments_with_jacobian`](https://joonho112.github.io/DPprior/reference/moments_with_jacobian.md)
for Jacobian computation,
[`exact_K_moments`](https://joonho112.github.io/DPprior/reference/exact_K_moments.md)
for moment verification

Other elicitation:
[`DPprior_a1()`](https://joonho112.github.io/DPprior/reference/DPprior_a1.md),
[`DPprior_a2_kl()`](https://joonho112.github.io/DPprior/reference/DPprior_a2_kl.md),
[`DPprior_dual()`](https://joonho112.github.io/DPprior/reference/DPprior_dual.md),
[`DPprior_fit()`](https://joonho112.github.io/DPprior/reference/DPprior_fit.md)

## Examples

``` r
# Basic usage
fit <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8)
print(fit)
#> DPprior Prior Elicitation Result
#> ============================================= 
#> 
#> Gamma Hyperprior: α ~ Gamma(a = 2.0361, b = 1.6051)
#>   E[α] = 1.269, SD[α] = 0.889
#> 
#> Target (J = 50):
#>   E[K_J]   = 5.00
#>   Var(K_J) = 8.00
#> 
#> Achieved:
#>   E[K_J] = 5.000000, Var(K_J) = 8.000000
#>   Residual = 7.60e-09
#> 
#> Method: A2-MN (6 iterations)

# Verify exact moment matching
achieved <- exact_K_moments(50, fit$a, fit$b)
cat(sprintf("Target E[K]=5, Achieved E[K]=%.10f\n", achieved$mean))
#> Target E[K]=5, Achieved E[K]=4.9999999992
cat(sprintf("Target Var=8, Achieved Var=%.10f\n", achieved$var))
#> Target Var=8, Achieved Var=8.0000000076

# Compare A1 vs A2 accuracy
a1 <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
a1_mom <- exact_K_moments(50, a1$a, a1$b)
a2_mom <- exact_K_moments(50, fit$a, fit$b)
cat(sprintf("A1 mean error: %.6f\n", abs(a1_mom$mean - 5)))
#> A1 mean error: 0.538649
cat(sprintf("A2 mean error: %.2e\n", abs(a2_mom$mean - 5)))
#> A2 mean error: 8.31e-10

# View iteration trace (includes step size and Jacobian determinant)
head(fit$trace)
#>   iter        a         b       M1         V     residual step   det_Jlog
#> 1    1 4.000000 3.9120230 4.461351  4.783136 3.261649e+00    1  -5.300969
#> 2    2 1.178650 0.9119694 4.909046 10.854537 2.855986e+00    1 -21.553612
#> 3    3 1.844384 1.4552538 4.974913  8.399473 4.002603e-01    1 -15.313512
#> 4    4 2.029223 1.5996801 4.999187  8.013243 1.326844e-02    1 -14.298307
#> 5    5 2.036082 1.6050455 4.999999  8.000021 2.078492e-05    1 -14.263052
#> 6    6 2.036093 1.6050541 5.000000  8.000000 7.600363e-09   NA -14.262997
```
