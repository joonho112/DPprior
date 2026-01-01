# A2-KL: KL Divergence Minimization for Prior Calibration

Calibrates Gamma hyperprior parameters \\(a, b)\\ by minimizing the
Kullback-Leibler divergence between a target distribution and the
induced marginal PMF of the number of clusters \\K_J\\.

## Usage

``` r
DPprior_a2_kl(
  J,
  target,
  method = c("pmf", "chisq"),
  max_iter = 100L,
  tol = 1e-06,
  M = .QUAD_NODES_DEFAULT,
  verbose = FALSE,
  ...
)
```

## Arguments

- J:

  Integer; sample size (number of observations). Must be \>= 2.

- target:

  Either:

  - Numeric vector of length J: target PMF for k = 1, ..., J (used when
    `method = "pmf"`).

  - Named list with `mu_K` and `var_K`: construct from moments using a
    discretized scaled chi-square (used when `method = "chisq"`).

- method:

  Character; `"pmf"` or `"chisq"`. Default: `"pmf"`.

- max_iter:

  Integer; maximum optimization iterations. Default: 100.

- tol:

  Numeric; convergence tolerance for optimization. Default: 1e-6.

- M:

  Integer; number of quadrature nodes. Default: 80.

- verbose:

  Logical; if TRUE, print optimization progress. Default: FALSE.

- ...:

  Optional tuning parameters:

  - `log_bounds`: numeric length-2 vector giving lower/upper bounds on
    `log(a)` and `log(b)` (default `c(-15, 15)`).

  - `eps`: numeric; small constant to avoid `log(0)` (default `1e-15`).

## Value

A `DPprior_fit` object (S3 class) with components:

- `a`:

  Numeric; optimal shape parameter

- `b`:

  Numeric; optimal rate parameter

- `J`:

  Integer; sample size

- `target`:

  List with target specification (pmf, mu_K, var_K, type, and for chisq:
  df, scale, mu_K_discrete, var_K_discrete)

- `method`:

  Character; "A2-KL"

- `converged`:

  Logical; whether optimization converged

- `iterations`:

  Integer; number of function evaluations

- `termination`:

  Character; termination reason

- `fit`:

  List with achieved `mu_K`, `var_K`, `kl`, and `residual`

- `diagnostics`:

  List with initialization info, optim details, fallback status, and KL
  at init/final

- `trace`:

  Data frame of evaluation-level trace (eval, a, b, kl)

- `status`:

  Character; convergence status

## Details

The A2-KL algorithm finds \\(a^\*, b^\*)\\ by solving: \$\$(a^\*, b^\*)
= \arg\min\_{a,b} D\_{KL}(p^\*(K_J) \\ p\_{a,b}(K_J))\$\$

**Algorithm:**

1.  Construct target PMF from user specification

2.  Initialize \\(a_0, b_0)\\ from A2-MN (exact moment matching)

3.  Optimize KL divergence using L-BFGS-B in log-space with bounds

4.  Apply fallback to initialization if optimization worsens KL

5.  Return optimal parameters and comprehensive diagnostics

**Stability features:**

- Log-parameterization ensures positivity of (a, b)

- Explicit bounds prevent numerical overflow/underflow

- Fallback mechanism if optimization diverges or worsens objective

- A2-MN initialization provides good starting point

**When to use A2-KL vs A2-MN:**

- **A2-MN**: Target moments only (exact moment matching)

- **A2-KL**: Full target distribution shape (multi-modal, skewed,
  expert-elicited)

## See also

[`DPprior_a2_newton`](https://joonho112.github.io/DPprior/reference/DPprior_a2_newton.md)
for exact moment matching,
[`DPprior_a1`](https://joonho112.github.io/DPprior/reference/DPprior_a1.md)
for closed-form initialization,
[`kl_divergence_K`](https://joonho112.github.io/DPprior/reference/kl_divergence_K.md)
for KL divergence computation,
[`discretize_chisq`](https://joonho112.github.io/DPprior/reference/discretize_chisq.md)
for chi-square discretization

## Examples

``` r
# Example 1: Target from moments (method = "chisq")
fit <- DPprior_a2_kl(J = 50, target = list(mu_K = 5, var_K = 8),
                     method = "chisq")
print(fit)
#> DPprior Prior Elicitation Result
#> ============================================= 
#> 
#> Gamma Hyperprior: α ~ Gamma(a = 2.2363, b = 1.7627)
#>   E[α] = 1.269, SD[α] = 0.848
#> 
#> Target (J = 50):
#>   E[K_J]   = 5.00
#>   Var(K_J) = 8.00
#> 
#> Achieved:
#>   E[K_J] = 5.019984, Var(K_J) = 7.640329
#>   Residual = 5.03e-03
#> 
#> Method: A2-KL (8 iterations)

# Example 2: Custom target PMF (method = "pmf")
target_pmf <- dbinom(1:50, size = 50, prob = 0.1)
fit2 <- DPprior_a2_kl(J = 50, target = target_pmf, method = "pmf")

# Example 3: Compare A2-KL vs A2-MN
a2_mn <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8)
a2_kl <- DPprior_a2_kl(J = 50, target = list(mu_K = 5, var_K = 8),
                       method = "chisq")
cat(sprintf("A2-MN: a=%.4f, b=%.4f\n", a2_mn$a, a2_mn$b))
#> A2-MN: a=2.0361, b=1.6051
cat(sprintf("A2-KL: a=%.4f, b=%.4f, KL=%.4e\n", a2_kl$a, a2_kl$b, a2_kl$fit$kl))
#> A2-KL: a=2.2363, b=1.7627, KL=5.0296e-03
```
