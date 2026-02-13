# Complete API Reference with Examples

## Overview

This vignette provides a complete API reference for the DPprior package.
All exported functions are documented with their full signatures,
parameter descriptions, return values, and usage examples.

The package functions are organized into six categories:

| Category                | Description                                                  | Key Functions                                                                                                                                                                                  |
|-------------------------|--------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Core Elicitation**    | Main user-facing calibration functions                       | [`DPprior_fit()`](https://joonho112.github.io/DPprior/reference/DPprior_fit.md), [`DPprior_dual()`](https://joonho112.github.io/DPprior/reference/DPprior_dual.md)                             |
| **Exact Computation**   | Stirling numbers and exact distributions                     | [`compute_log_stirling()`](https://joonho112.github.io/DPprior/reference/compute_log_stirling.md), [`pmf_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/pmf_K_given_alpha.md) |
| **Approximation**       | Closed-form A1 and Newton-based A2 methods                   | [`DPprior_a1()`](https://joonho112.github.io/DPprior/reference/DPprior_a1.md), [`DPprior_a2_newton()`](https://joonho112.github.io/DPprior/reference/DPprior_a2_newton.md)                     |
| **Weight Distribution** | First stick-breaking weight $w_{1}$ and co-clustering $\rho$ | [`mean_w1()`](https://joonho112.github.io/DPprior/reference/mean_w1.md), [`prob_w1_exceeds()`](https://joonho112.github.io/DPprior/reference/prob_w1_exceeds.md)                               |
| **Diagnostics**         | Prior validation and dominance risk assessment               | [`DPprior_diagnostics()`](https://joonho112.github.io/DPprior/reference/DPprior_diagnostics.md), [`summary()`](https://rdrr.io/r/base/summary.html)                                            |
| **Utility**             | Helper functions for conversions and integration             | [`vif_to_variance()`](https://joonho112.github.io/DPprior/reference/vif_to_variance.md), [`integrate_gamma()`](https://joonho112.github.io/DPprior/reference/integrate_gamma.md)               |

------------------------------------------------------------------------

## 1. Core Elicitation Functions

These are the primary user-facing functions for eliciting Gamma
hyperpriors on the DP concentration parameter $\alpha$.

### 1.1 `DPprior_fit()`

**Unified Interface for Prior Elicitation**

The main entry point for eliciting a Gamma hyperprior. Automatically
selects the appropriate algorithm based on input specification.

``` r
DPprior_fit(
  J,                      # Sample size (required)
  mu_K,                   # Target E[K_J] (required)
  var_K = NULL,           # Target Var(K_J) (optional)
  confidence = NULL,      # "low", "medium", "high" (optional)
  method = "A2-MN",       # "A1", "A2-MN", "A2-KL"
  check_diagnostics = FALSE,
  warn_dominance = TRUE,
  M = 80L,                # Quadrature nodes
  verbose = FALSE
)
```

**Parameters:**

| Parameter           | Type      | Description                                                                                       |
|---------------------|-----------|---------------------------------------------------------------------------------------------------|
| `J`                 | Integer   | Sample size (number of observations/sites). Must be $\geq 2$.                                     |
| `mu_K`              | Numeric   | Target prior mean ${\mathbb{E}}\left\lbrack K_{J} \right\rbrack$. Must satisfy $1 < \mu_{K} < J$. |
| `var_K`             | Numeric   | Target prior variance $\text{Var}\left( K_{J} \right)$. If `NULL`, computed from `confidence`.    |
| `confidence`        | Character | Qualitative uncertainty: `"low"` (VIF=4), `"medium"` (VIF=2.5), `"high"` (VIF=1.5).               |
| `method`            | Character | Algorithm: `"A1"` (closed-form), `"A2-MN"` (Newton), `"A2-KL"` (KL minimization).                 |
| `check_diagnostics` | Logical   | If `TRUE`, compute full diagnostics.                                                              |
| `warn_dominance`    | Logical   | If `TRUE`, warn when $P\left( w_{1} > 0.5 \right) > 0.4$.                                         |
| `M`                 | Integer   | Number of Gauss-Laguerre quadrature nodes.                                                        |
| `verbose`           | Logical   | Print iteration progress.                                                                         |

**Returns:** An S3 object of class `DPprior_fit` with components:

| Component     | Description                              |
|---------------|------------------------------------------|
| `a`           | Optimal Gamma shape parameter            |
| `b`           | Optimal Gamma rate parameter             |
| `J`           | Sample size                              |
| `target`      | List with `mu_K`, `var_K`, `var_K_used`  |
| `fit`         | Achieved `mu_K`, `var_K`, and `residual` |
| `method`      | Algorithm used                           |
| `converged`   | Logical convergence indicator            |
| `iterations`  | Number of iterations (0 for A1)          |
| `diagnostics` | Optional diagnostic information          |

**Examples:**

``` r
# Basic usage with confidence level
fit1 <- DPprior_fit(J = 50, mu_K = 5, confidence = "medium")
#> Warning: HIGH DOMINANCE RISK: P(w1 > 0.5) = 49.7% exceeds 40%.
#>   This may indicate unintended prior behavior (Lee, 2026).
#>   Consider using DPprior_dual() for weight-constrained elicitation.
#>   See ?DPprior_diagnostics for interpretation.
print(fit1)
#> DPprior Prior Elicitation Result
#> ============================================= 
#> 
#> Gamma Hyperprior: α ~ Gamma(a = 1.4082, b = 1.0770)
#>   E[α] = 1.308, SD[α] = 1.102
#> 
#> Target (J = 50):
#>   E[K_J]   = 5.00
#>   Var(K_J) = 10.00
#>   (from confidence = 'medium')
#> 
#> Achieved:
#>   E[K_J] = 5.000000, Var(K_J) = 10.000000
#>   Residual = 3.94e-10
#> 
#> Method: A2-MN (7 iterations)
#> 
#> Dominance Risk: HIGH ✘ (P(w₁>0.5) = 50%)

# Direct variance specification
fit2 <- DPprior_fit(J = 50, mu_K = 5, var_K = 10)
#> Warning: HIGH DOMINANCE RISK: P(w1 > 0.5) = 49.7% exceeds 40%.
#>   This may indicate unintended prior behavior (Lee, 2026).
#>   Consider using DPprior_dual() for weight-constrained elicitation.
#>   See ?DPprior_diagnostics for interpretation.
cat("Gamma(", round(fit2$a, 4), ", ", round(fit2$b, 4), ")\n", sep = "")
#> Gamma(1.4082, 1.077)

# Using A1 closed-form (faster for large J)
fit3 <- DPprior_fit(J = 200, mu_K = 15, var_K = 30, method = "A1")

# With diagnostics
fit4 <- DPprior_fit(J = 50, mu_K = 5, var_K = 8, check_diagnostics = TRUE)
#> Warning: HIGH DOMINANCE RISK: P(w1 > 0.5) = 48.1% exceeds 40%.
#>   This may indicate unintended prior behavior (Lee, 2026).
#>   Consider using DPprior_dual() for weight-constrained elicitation.
#>   See ?DPprior_diagnostics for interpretation.
```

### 1.2 `DPprior_dual()`

**Dual-Anchor Elicitation with Weight Control**

Extends moment-based elicitation to incorporate constraints on the
stick-breaking weight distribution, addressing the “unintended prior”
problem identified in Lee (2026, Section 4).

``` r
DPprior_dual(
  fit_K,                  # K-only DPprior_fit object (required)
  w1_target,              # Weight target specification (required)
  lambda = 0.7,           # Trade-off parameter [0, 1]
  loss_type = "adaptive", # "absolute", "relative", "adaptive"
  method = "L-BFGS-B",    # Optimization method
  max_iter = 100L,
  tol = 1e-6,
  M = 80L,
  verbose = FALSE
)
```

**Parameters:**

| Parameter   | Type        | Description                                                                                                                                                                                             |
|-------------|-------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `fit_K`     | DPprior_fit | Initial K-only fit from [`DPprior_fit()`](https://joonho112.github.io/DPprior/reference/DPprior_fit.md) or [`DPprior_a2_newton()`](https://joonho112.github.io/DPprior/reference/DPprior_a2_newton.md). |
| `w1_target` | List        | Weight target: `list(prob = list(threshold, value))` for $P\left( w_{1} > \text{threshold} \right) = \text{value}$.                                                                                     |
| `lambda`    | Numeric     | Trade-off: 0 = pure weight, 1 = pure K. Default: 0.7.                                                                                                                                                   |
| `loss_type` | Character   | Loss function type for weight anchor.                                                                                                                                                                   |
| `method`    | Character   | Optimization algorithm for [`optim()`](https://rdrr.io/r/stats/optim.html).                                                                                                                             |
| `max_iter`  | Integer     | Maximum optimization iterations.                                                                                                                                                                        |
| `tol`       | Numeric     | Convergence tolerance.                                                                                                                                                                                  |

**Returns:** A `DPprior_fit` object with additional `dual_anchor`
component containing:

| Component     | Description                |
|---------------|----------------------------|
| `lambda`      | Trade-off parameter used   |
| `w1_target`   | Original weight target     |
| `w1_achieved` | Achieved weight statistics |
| `K_achieved`  | Achieved K moments         |
| `init`        | Initial K-only parameters  |

**Examples:**

``` r
# Step 1: Fit K-only prior
fit_K <- DPprior_fit(J = 50, mu_K = 5, var_K = 8)
#> Warning: HIGH DOMINANCE RISK: P(w1 > 0.5) = 48.1% exceeds 40%.
#>   This may indicate unintended prior behavior (Lee, 2026).
#>   Consider using DPprior_dual() for weight-constrained elicitation.
#>   See ?DPprior_diagnostics for interpretation.
cat("K-only: P(w1 > 0.5) =", 
    round(prob_w1_exceeds(0.5, fit_K$a, fit_K$b), 3), "\n")
#> K-only: P(w1 > 0.5) = 0.481

# Step 2: Apply weight constraint
w1_target <- list(prob = list(threshold = 0.5, value = 0.30))
fit_dual <- DPprior_dual(fit_K, w1_target, lambda = 0.5)
cat("Dual:   P(w1 > 0.5) =", 
    round(prob_w1_exceeds(0.5, fit_dual$a, fit_dual$b), 3), "\n")
#> Dual:   P(w1 > 0.5) = 0.438

# Compare parameters
cat("\nParameter comparison:\n")
#> 
#> Parameter comparison:
cat("  K-only: Gamma(", round(fit_K$a, 3), ", ", round(fit_K$b, 3), ")\n", sep = "")
#>   K-only: Gamma(2.036, 1.605)
cat("  Dual:   Gamma(", round(fit_dual$a, 3), ", ", round(fit_dual$b, 3), ")\n", sep = "")
#>   Dual:   Gamma(2.575, 1.834)
```

------------------------------------------------------------------------

## 2. Exact Computation Functions

These functions provide exact computation of Stirling numbers, PMFs, and
moments using numerically stable log-space arithmetic.

### 2.1 Stirling Number Functions

#### `compute_log_stirling()`

**Pre-compute Log Stirling Numbers**

Computes the logarithm of unsigned Stirling numbers of the first kind
$\left| s(J,k) \right|$ for all $J$ from 0 to `J_max` and $k$ from 0 to
$J$.

``` r
compute_log_stirling(J_max)
```

**Parameters:**

| Parameter | Type    | Description                              |
|-----------|---------|------------------------------------------|
| `J_max`   | Integer | Maximum sample size. Must be $\leq 500$. |

**Returns:** A lower triangular matrix of dimension
$\left( J_{max} + 1 \right) \times \left( J_{max} + 1 \right)$. Entry
`[J+1, k+1]` contains $\log\left| s(J,k) \right|$ (using R’s 1-based
indexing).

**Examples:**

``` r
# Pre-compute for J up to 100
logS <- compute_log_stirling(100)

# Access |s(10,3)| = 9450
s_10_3 <- exp(logS[11, 4])
cat("|s(10,3)| =", round(s_10_3), "\n")
#> |s(10,3)| = 1172700

# Verify row sum identity: sum_k |s(J,k)| = J!
J <- 6
row_sum <- sum(exp(logS[J+1, 2:(J+1)]))
cat("sum |s(6,k)| =", round(row_sum), ", 6! =", factorial(6), "\n")
#> sum |s(6,k)| = 720 , 6! = 720
```

#### `get_log_stirling()`

**Safe Accessor with Bounds Checking**

``` r
get_log_stirling(J, k, logS)
```

Returns $\log\left| s(J,k) \right|$ with automatic bounds checking.
Returns `-Inf` for invalid indices ($k > J$ or $k < 1$).

#### `get_stirling_row()`

**Extract Row of Stirling Numbers**

``` r
get_stirling_row(J, logS)
```

Returns a vector of $\log\left| s(J,k) \right|$ for $k = 1,\ldots,J$.

### 2.2 Conditional PMF Functions

#### `pmf_K_given_alpha()`

**Antoniak Distribution PMF**

Computes the exact conditional PMF $P\left( K_{J} = k|\alpha \right)$
using the Antoniak distribution formula:
$$P\left( K_{J} = k|\alpha \right) = \frac{\left| s(J,k) \right| \cdot \alpha^{k}}{(\alpha)_{J}}$$

``` r
pmf_K_given_alpha(J, alpha, logS, normalize = TRUE)
```

**Parameters:**

| Parameter   | Type    | Description                              |
|-------------|---------|------------------------------------------|
| `J`         | Integer | Sample size.                             |
| `alpha`     | Numeric | Concentration parameter (scalar, $> 0$). |
| `logS`      | Matrix  | Pre-computed log-Stirling matrix.        |
| `normalize` | Logical | If `TRUE`, ensure PMF sums to 1.         |

**Returns:** Numeric vector of length $J + 1$ containing
$P\left( K_{J} = k|\alpha \right)$ for $k = 0,1,\ldots,J$. Note:
$P\left( K_{J} = 0|\alpha \right) = 0$ always.

**Examples:**

``` r
J <- 50
alpha <- 2.0
logS <- compute_log_stirling(J)

# Compute conditional PMF
pmf <- pmf_K_given_alpha(J, alpha, logS)

# Find mode
mode_k <- which.max(pmf) - 1
cat("Mode of K|α=2:", mode_k, "\n")
#> Mode of K|α=2: 7

# Verify normalization
cat("PMF sum:", sum(pmf), "\n")
#> PMF sum: 1
```

#### `cdf_K_given_alpha()`

**Conditional CDF**

``` r
cdf_K_given_alpha(J, alpha, logS)
```

Returns the cumulative distribution function
$F(k) = P\left( K_{J} \leq k|\alpha \right)$.

#### `quantile_K_given_alpha()`

**Conditional Quantile Function**

``` r
quantile_K_given_alpha(p, J, alpha, logS)
```

Returns the smallest $k$ such that
$P\left( K_{J} \leq k|\alpha \right) \geq p$.

### 2.3 Conditional Moment Functions

#### `mean_K_given_alpha()`

**Conditional Mean via Digamma**

Computes
${\mathbb{E}}\left\lbrack K_{J}|\alpha \right\rbrack = \alpha\{\psi(\alpha + J) - \psi(\alpha)\}$
using the digamma function.

``` r
mean_K_given_alpha(J, alpha)
```

**Parameters:**

| Parameter | Type    | Description                           |
|-----------|---------|---------------------------------------|
| `J`       | Integer | Sample size ($\geq 1$).               |
| `alpha`   | Numeric | Concentration parameter (vectorized). |

**Returns:** Numeric vector of conditional means.

**Examples:**

``` r
J <- 50

# Single alpha
mean_K_given_alpha(J, 2.0)
#> [1] 7.037626

# Vectorized
alpha_seq <- c(0.5, 1, 2, 5, 10)
means <- mean_K_given_alpha(J, alpha_seq)
data.frame(alpha = alpha_seq, E_K = round(means, 2))
#>   alpha   E_K
#> 1   0.5  2.94
#> 2   1.0  4.50
#> 3   2.0  7.04
#> 4   5.0 12.46
#> 5  10.0 18.34
```

#### `var_K_given_alpha()`

**Conditional Variance via Trigamma**

Computes
$\text{Var}\left( K_{J}|\alpha \right) = \mu_{J}(\alpha) - \alpha^{2}\{\psi_{1}(\alpha) - \psi_{1}(\alpha + J)\}$.

``` r
var_K_given_alpha(J, alpha)
```

**Key Property:** Conditional underdispersion always holds:
$$0 < \text{Var}\left( K_{J}|\alpha \right) < {\mathbb{E}}\left\lbrack K_{J}|\alpha \right\rbrack$$

**Examples:**

``` r
J <- 50
alpha <- 2.0

mu <- mean_K_given_alpha(J, alpha)
sigma2 <- var_K_given_alpha(J, alpha)

cat("E[K|α=2] =", round(mu, 4), "\n")
#> E[K|α=2] = 7.0376
cat("Var(K|α=2) =", round(sigma2, 4), "\n")
#> Var(K|α=2) = 4.5356
cat("Underdispersion ratio:", round(sigma2/mu, 4), "< 1 ✓\n")
#> Underdispersion ratio: 0.6445 < 1 ✓
```

#### `moments_K_given_alpha()`

**Both Moments in One Call**

``` r
moments_K_given_alpha(J, alpha)
```

Returns a list with `mean` and `var`.

### 2.4 Marginal Distribution Functions

#### `exact_K_moments()`

**Marginal Moments under Gamma Hyperprior**

Computes ${\mathbb{E}}\left\lbrack K_{J}|a,b \right\rbrack$ and
$\text{Var}\left( K_{J}|a,b \right)$ when
$\alpha \sim \text{Gamma}(a,b)$ using Gauss-Laguerre quadrature.

``` r
exact_K_moments(J, a, b, M = 80L)
```

**Parameters:**

| Parameter | Type    | Description                 |
|-----------|---------|-----------------------------|
| `J`       | Integer | Sample size.                |
| `a`       | Numeric | Gamma shape parameter.      |
| `b`       | Numeric | Gamma rate parameter.       |
| `M`       | Integer | Number of quadrature nodes. |

**Returns:** A list with components:

| Component | Description                |
|-----------|----------------------------|
| `mean`    | Marginal mean $M_{1}(a,b)$ |
| `var`     | Marginal variance $V(a,b)$ |
| `sd`      | Standard deviation         |
| `cv`      | Coefficient of variation   |

**Examples:**

``` r
# Compute marginal moments
result <- exact_K_moments(J = 50, a = 1.6, b = 1.22)
cat("E[K_50] =", round(result$mean, 4), "\n")
#> E[K_50] = 5.0454
cat("Var(K_50) =", round(result$var, 4), "\n")
#> Var(K_50) = 9.3797
cat("CV(K_50) =", round(result$cv, 4), "\n")
#> CV(K_50) = 0.607
```

#### `pmf_K_marginal()`

**Marginal PMF via Quadrature Mixture**

``` r
pmf_K_marginal(J, a, b, logS, M = 80L)
```

Computes the marginal PMF by mixing conditional PMFs over
$\alpha \sim \text{Gamma}(a,b)$:
$$P\left( K_{J} = k|a,b \right) = \int_{0}^{\infty}P\left( K_{J} = k|\alpha \right) \cdot g_{a,b}(\alpha)\, d\alpha$$

**Examples:**

``` r
J <- 50
logS <- compute_log_stirling(J)

pmf_marginal <- pmf_K_marginal(J, a = 1.6, b = 1.22, logS)

# Find mode
mode_k <- which.max(pmf_marginal) - 1
cat("Mode of marginal K:", mode_k, "\n")
#> Mode of marginal K: 3
```

#### `summary_K_marginal()`

**Complete Marginal Distribution Summary**

``` r
summary_K_marginal(J, a, b, logS, M = 80L, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
```

Returns mean, variance, mode, median, quantiles, PMF, and CDF.

------------------------------------------------------------------------

## 3. Approximation Functions

### 3.1 `DPprior_a1()`

**A1 Closed-Form Approximation**

Fast closed-form solution using the negative binomial approximation to
the marginal distribution of $K_{J}$.

``` r
DPprior_a1(
  J,                      # Sample size
  mu_K,                   # Target mean
  var_K = NULL,           # Target variance
  confidence = NULL,      # "low", "medium", "high"
  scaling = "log",        # "log", "harmonic", "digamma"
  project_to_feasible = TRUE
)
```

**Key Formulas:**
$$a = \frac{\left( \mu_{K} - 1 \right)^{2}}{\sigma_{K}^{2} - \left( \mu_{K} - 1 \right)},\quad b = \frac{\left( \mu_{K} - 1 \right) \cdot c_{J}}{\sigma_{K}^{2} - \left( \mu_{K} - 1 \right)}$$

where $c_{J}$ is the scaling constant (default: $\log J$).

**Feasibility Constraint:** Requires $\sigma_{K}^{2} > \mu_{K} - 1$.

**Examples:**

``` r
# Basic A1 fit
fit_a1 <- DPprior_a1(J = 50, mu_K = 5, var_K = 8)
cat("A1 solution: Gamma(", round(fit_a1$a, 4), ", ", 
    round(fit_a1$b, 4), ")\n", sep = "")
#> A1 solution: Gamma(4, 3.912)

# Compare scaling methods
for (s in c("log", "harmonic", "digamma")) {
  fit <- DPprior_a1(50, 5, 8, scaling = s)
  cat(sprintf("  %s: cJ = %.4f\n", s, fit$cJ))
}
#>   log: cJ = 3.9120
#>   harmonic: cJ = 4.4792
#>   digamma: cJ = 4.4633
```

### 3.2 `DPprior_a2_newton()`

**A2-MN Newton Solver for Exact Moment Matching**

Finds $\left( a^{*},b^{*} \right)$ such that the induced marginal
moments exactly match targets.

``` r
DPprior_a2_newton(
  J,                      # Sample size
  mu_K,                   # Target mean
  var_K,                  # Target variance
  a0 = NULL,              # Initial shape (from A1 if NULL)
  b0 = NULL,              # Initial rate (from A1 if NULL)
  tol_F = 1e-8,           # Residual tolerance
  tol_step = 1e-10,       # Step size tolerance
  max_iter = 20L,
  damping = TRUE,         # Backtracking line search
  use_fallback = TRUE,    # Nelder-Mead fallback
  M = 80L,
  verbose = FALSE
)
```

**Algorithm:**

1.  Initialize from A1 closed-form
2.  Log-parameterize: $\eta = \left( \log a,\log b \right)$
3.  Newton iteration with score-based Jacobian
4.  Backtracking line search for global convergence

**Convergence:** Typically achieves machine-precision accuracy in 3-8
iterations.

**Examples:**

``` r
# Exact moment matching
fit <- DPprior_a2_newton(J = 50, mu_K = 5, var_K = 8, verbose = TRUE)
#> A2-MN Newton Solver
#> Target: E[K]=5.0000, Var(K)=8.0000
#> A1 initialization: a0=4.000000, b0=3.912023
#> -------------------------------------------------------------------------------- 
#> Iter |          a |          b |       E[K] |     Var(K) |      ||F|| |     step |     det(J)
#> -------------------------------------------------------------------------------- 
#>    1 |   4.000000 |   3.912023 |   4.461351 |   4.783136 |   3.26e+00 |   1.0000 |  -5.30e+00
#>    2 |   1.178650 |   0.911969 |   4.909046 |  10.854537 |   2.86e+00 |   1.0000 |  -2.16e+01
#>    3 |   1.844384 |   1.455254 |   4.974913 |   8.399473 |   4.00e-01 |   1.0000 |  -1.53e+01
#>    4 |   2.029223 |   1.599680 |   4.999187 |   8.013243 |   1.33e-02 |   1.0000 |  -1.43e+01
#>    5 |   2.036082 |   1.605046 |   4.999999 |   8.000021 |   2.08e-05 |   1.0000 |  -1.43e+01
#>    6 |   2.036093 |   1.605054 |   5.000000 |   8.000000 |   7.60e-09 |      --- |  -1.43e+01
#> 
#> Converged: ||F|| = 7.60e-09 < 1.00e-08

# Verify exact matching
achieved <- exact_K_moments(50, fit$a, fit$b)
cat("\nTarget vs Achieved:\n")
#> 
#> Target vs Achieved:
cat("  E[K]:   5.000000 vs", sprintf("%.10f", achieved$mean), "\n")
#>   E[K]:   5.000000 vs 4.9999999992
cat("  Var(K): 8.000000 vs", sprintf("%.10f", achieved$var), "\n")
#>   Var(K): 8.000000 vs 8.0000000076
```

### 3.3 `DPprior_a2_kl()`

**A2-KL Distribution Matching via KL Divergence**

Minimizes Kullback-Leibler divergence between a target PMF and the
induced marginal PMF of $K_{J}$.

``` r
DPprior_a2_kl(
  J,                      # Sample size
  target,                 # Target PMF or list(mu_K, var_K)
  method = c("pmf", "chisq"),  # Target type
  max_iter = 100L,
  tol = 1e-6,
  M = 80L,
  verbose = FALSE
)
```

**Target Specification:**

| `method`  | `target` Input                  | Description                             |
|-----------|---------------------------------|-----------------------------------------|
| `"pmf"`   | Numeric vector of length $J$    | Direct PMF on $\{ 1,\ldots,J\}$         |
| `"chisq"` | `list(mu_K = ..., var_K = ...)` | Discretized chi-square matching moments |

**Examples:**

``` r
# Method 1: Moment-based target using chi-square discretization
fit_kl <- DPprior_a2_kl(J = 50, target = list(mu_K = 5, var_K = 8), 
                         method = "chisq")
cat("A2-KL solution: Gamma(", round(fit_kl$a, 4), ", ", 
    round(fit_kl$b, 4), ")\n", sep = "")
#> A2-KL solution: Gamma(2.2363, 1.7627)
cat("Final KL divergence:", format(fit_kl$fit$kl, scientific = TRUE), "\n")
#> Final KL divergence: 5.029612e-03

# Method 2: Direct PMF target
target_pmf <- discretize_chisq(J = 50, df = 6.25, scale = 0.8)
fit_kl2 <- DPprior_a2_kl(J = 50, target = target_pmf, method = "pmf")
```

------------------------------------------------------------------------

## 4. Weight Distribution Functions

Functions for computing properties of the first stick-breaking weight
$w_{1}$ and the co-clustering probability $\rho = \sum_{h}w_{h}^{2}$.

### 4.1 First Weight $w_{1}$ Functions

The size-biased first weight $w_{1}$ follows a compound distribution:
$$w_{1}|\alpha \sim \text{Beta}(1,\alpha),\quad\alpha \sim \text{Gamma}(a,b)$$

#### `mean_w1()`

**Marginal Mean of $w_{1}$**

``` r
mean_w1(a, b, M = 80L)
```

Computes
${\mathbb{E}}\left\lbrack w_{1}|a,b \right\rbrack = {\mathbb{E}}\left\lbrack 1/(1 + \alpha) \right\rbrack$
via quadrature.

**Examples:**

``` r
mean_w1(a = 1.6, b = 1.22)
#> [1] 0.508368
```

#### `var_w1()`

**Marginal Variance of $w_{1}$**

``` r
var_w1(a, b, M = 80L)
```

#### `quantile_w1()`

**Marginal Quantiles of $w_{1}$**

``` r
quantile_w1(p, a, b, n_grid = 1000, M = 80L)
```

Computes quantiles via numerical inversion of the CDF.

**Examples:**

``` r
# Median of w1
median_w1 <- quantile_w1(0.5, a = 1.6, b = 1.22)
cat("Median(w1) =", round(median_w1, 4), "\n")
#> Median(w1) = 0.4839
```

#### `prob_w1_exceeds()`

**Tail Probability $P\left( w_{1} > x \right)$**

``` r
prob_w1_exceeds(x, a, b, M = 80L)
```

**Critical for Dominance Assessment:** $P\left( w_{1} > 0.5 \right)$
indicates the probability that a single cluster dominates.

**Examples:**

``` r
# Dominance probabilities
a <- 1.6; b <- 1.22
cat("P(w1 > 0.5) =", round(prob_w1_exceeds(0.5, a, b), 4), "\n")
#> P(w1 > 0.5) = 0.4868
cat("P(w1 > 0.7) =", round(prob_w1_exceeds(0.7, a, b), 4), "\n")
#> P(w1 > 0.7) = 0.3334
cat("P(w1 > 0.9) =", round(prob_w1_exceeds(0.9, a, b), 4), "\n")
#> P(w1 > 0.9) = 0.1833
```

#### `cdf_w1()` and `pdf_w1()`

**CDF and PDF of $w_{1}$**

``` r
cdf_w1(x, a, b, M = 80L)
pdf_w1(x, a, b, M = 80L)
```

### 4.2 Co-Clustering Probability $\rho$ Functions

The co-clustering probability $\rho = \sum_{h = 1}^{\infty}w_{h}^{2}$
represents the probability that two randomly chosen observations belong
to the same cluster.

#### `mean_rho()` and `var_rho()`

**Marginal Moments of $\rho$**

``` r
mean_rho(a, b, M = 80L)
var_rho(a, b, M = 80L)
```

**Key Identity:**
${\mathbb{E}}\left\lbrack \rho|\alpha \right\rbrack = {\mathbb{E}}\left\lbrack w_{1}|\alpha \right\rbrack = 1/(1 + \alpha)$,
so
${\mathbb{E}}\left\lbrack \rho|a,b \right\rbrack = {\mathbb{E}}\left\lbrack w_{1}|a,b \right\rbrack$.

**Examples:**

``` r
a <- 1.6; b <- 1.22
cat("E[rho] =", round(mean_rho(a, b), 4), "\n")
#> E[rho] = 0.5084
cat("E[w1]  =", round(mean_w1(a, b), 4), "(same as E[rho])\n")
#> E[w1]  = 0.5084 (same as E[rho])
cat("Var(rho) =", round(var_rho(a, b), 4), "\n")
#> Var(rho) = 0.071
cat("Var(w1)  =", round(var_w1(a, b), 4), "(different!)\n")
#> Var(w1)  = 0.1052 (different!)
```

#### `cv_rho()`

**Coefficient of Variation of $\rho$**

``` r
cv_rho(a, b, M = 80L)
```

### 4.3 Conditional Weight Functions

#### `mean_rho_given_alpha()` and `var_rho_given_alpha()`

``` r
mean_rho_given_alpha(alpha)   # Returns 1/(1 + alpha)
var_rho_given_alpha(alpha)    # Returns 2*alpha / ((1+alpha)^2 * (2+alpha) * (3+alpha))
```

------------------------------------------------------------------------

## 5. Diagnostic Functions

### 5.1 `DPprior_diagnostics()`

**Comprehensive Prior Diagnostics**

Computes a full diagnostic report implementing the “unintended prior”
checks from Lee (2026, Section 4).

``` r
DPprior_diagnostics(fit, thresholds = c(0.5, 0.9))
```

**Parameters:**

| Parameter    | Type        | Description                                        |
|--------------|-------------|----------------------------------------------------|
| `fit`        | DPprior_fit | Fitted object with `a`, `b`, `J`.                  |
| `thresholds` | Numeric     | Thresholds for $P\left( w_{1} > x \right)$ checks. |

**Returns:** An S3 object of class `DPprior_diagnostics` with:

| Component      | Description                                              |
|----------------|----------------------------------------------------------|
| `alpha`        | Summary of $\alpha$ distribution (mean, CV, quantiles)   |
| `K`            | Summary of $K_{J}$ distribution (mean, SD, mode, median) |
| `weights`      | $w_{1}$ distribution and dominance risk assessment       |
| `coclustering` | $\rho$ summary with interpretation                       |
| `warnings`     | Character vector of diagnostic warnings                  |

**Warning Thresholds:**

- `HIGH DOMINANCE RISK`: $P\left( w_{1} > 0.5 \right) > 40\%$
- `NEAR-DEGENERATE RISK`: $P\left( w_{1} > 0.9 \right) > 15\%$

**Examples:**

``` r
fit <- DPprior_fit(J = 50, mu_K = 5, var_K = 8)
#> Warning: HIGH DOMINANCE RISK: P(w1 > 0.5) = 48.1% exceeds 40%.
#>   This may indicate unintended prior behavior (Lee, 2026).
#>   Consider using DPprior_dual() for weight-constrained elicitation.
#>   See ?DPprior_diagnostics for interpretation.
diag <- DPprior_diagnostics(fit)
print(diag)
#> DPprior Comprehensive Diagnostics
#> ============================================================ 
#> 
#> Prior: alpha ~ Gamma(2.0361, 1.6051) for J = 50
#> 
#> alpha Distribution:
#> ---------------------------------------- 
#>   E[alpha] = 1.269, CV(alpha) = 0.701, Median = 1.068
#>   90% CI: [0.230, 2.992]
#> 
#> K_J Distribution:
#> ---------------------------------------- 
#>   E[K] = 5.00, SD(K) = 2.83, Mode = 3
#>   Median = 5, IQR = [3, 7]
#> 
#> w1 Distribution (Size-Biased First Weight):
#> ---------------------------------------- 
#>   E[w1] = 0.501, Median = 0.478
#>   P(w1 > 0.5) = 48.1% (dominance risk: HIGH)
#>   P(w1 > 0.9) = 16.3%
#> 
#> Co-Clustering (rho = sum w_h^2):
#> ---------------------------------------- 
#>   E[rho] = 0.501 (High prior co-clustering: most unit pairs expected in same cluster)
#> 
#> WARNINGS:
#> ---------------------------------------- 
#>   * HIGH DOMINANCE RISK: P(w1 > 0.5) = 48.1% exceeds 40%
#>   * NEAR-DEGENERATE RISK: P(w1 > 0.9) = 16.3% exceeds 15%
#> 
#>   Consider using DPprior_dual() for weight-constrained elicitation.
```

### 5.2 S3 Methods for `DPprior_fit`

#### `print.DPprior_fit()`

``` r
print(fit)
```

Displays a concise summary of the elicitation result.

#### `summary.DPprior_fit()`

``` r
summary(fit, print_output = TRUE)
```

**Returns:** An object of class `summary.DPprior_fit` with detailed
statistics.

#### `plot.DPprior_fit()`

``` r
plot(fit, which = c("alpha", "K", "w1", "table"))
```

Creates a four-panel visualization dashboard:

- Panel A: Gamma prior density for $\alpha$
- Panel B: Marginal PMF of $K_{J}$
- Panel C: PDF of $w_{1}$ with dominance risk
- Panel D: Summary statistics table

### 5.3 Individual Plot Functions

``` r
plot_alpha_prior(fit)   # Alpha distribution only
plot_K_prior(fit)       # K distribution only
plot_w1_prior(fit)      # w1 distribution with dominance
```

### 5.4 Dominance Risk Functions

#### `check_dominance_risk()`

``` r
check_dominance_risk(a, b, threshold = 0.5, risk_level = 0.3)
```

Returns `TRUE` if
$P\left( w_{1} > \text{threshold} \right) > \text{risk\_level}$.

#### `compute_alpha_diagnostics()`

``` r
compute_alpha_diagnostics(a, b)
```

Returns mean, SD, CV, and quantiles of $\alpha \sim \text{Gamma}(a,b)$.

#### `compute_K_diagnostics()`

``` r
compute_K_diagnostics(J, a, b, M = 80L)
```

Returns mean, SD, mode, median, and quantiles of $K_{J}$.

#### `compute_weight_diagnostics()`

``` r
compute_weight_diagnostics(a, b, thresholds = c(0.5, 0.9), M = 80L)
```

Returns $w_{1}$ summary with dominance risk classification.

------------------------------------------------------------------------

## 6. Utility Functions

### 6.1 Variance Conversion Functions

#### `vif_to_variance()`

**Convert Variance Inflation Factor to Variance**

``` r
vif_to_variance(mu_K, vif)
```

Computes
$\sigma_{K}^{2} = \text{VIF} \times \left( \mu_{K} - 1 \right)$, based
on the marginal overdispersion relationship.

**Examples:**

``` r
# VIF = 2 means variance is twice the Poisson-like baseline
var_K <- vif_to_variance(mu_K = 5, vif = 2)
cat("var_K =", var_K, "\n")
#> var_K = 8
```

#### `confidence_to_vif()`

**Map Confidence Level to VIF**

``` r
confidence_to_vif(confidence)
```

| Confidence | VIF | Interpretation       |
|------------|-----|----------------------|
| `"low"`    | 4.0 | High uncertainty     |
| `"medium"` | 2.5 | Moderate uncertainty |
| `"high"`   | 1.5 | Low uncertainty      |

#### `cv_alpha_to_variance()`

**Convert CV of $\alpha$ to Variance of $K$**

``` r
cv_alpha_to_variance(mu_K, cv_alpha)
```

### 6.2 Scaling Functions

#### `compute_scaling_constant()`

``` r
compute_scaling_constant(J, scaling = "log", mu_K = NULL)
```

Computes the scaling constant $c_{J}$ used in the A1 approximation:

| `scaling`    | Formula                |
|--------------|------------------------|
| `"log"`      | $\log J$               |
| `"harmonic"` | $\sum_{i = 1}^{J}1/i$  |
| `"digamma"`  | $\psi(J + 1) + \gamma$ |

### 6.3 Quadrature Functions

#### `gauss_laguerre_nodes()`

**Compute Gauss-Laguerre Quadrature Nodes and Weights**

``` r
gauss_laguerre_nodes(M, alpha_param = 0)
```

**Parameters:**

| Parameter     | Type    | Description                                                   |
|---------------|---------|---------------------------------------------------------------|
| `M`           | Integer | Number of quadrature nodes.                                   |
| `alpha_param` | Numeric | Generalized Laguerre parameter (for Gamma(a,b), use $a - 1$). |

**Returns:** List with `nodes`, `weights`, and `weights_log`.

#### `build_gamma_quadrature()`

**Build Quadrature for Gamma Distribution**

``` r
build_gamma_quadrature(a, b, M = 80L)
```

Transforms standard Laguerre quadrature to integrate against Gamma(a,b).

**Returns:** List with `alpha_nodes` and `weights_normalized`.

#### `integrate_gamma()`

**Compute Expectation under Gamma Distribution**

``` r
integrate_gamma(f, a, b, M = 80L)
```

Computes
${\mathbb{E}}_{\alpha \sim \text{Gamma}{(a,b)}}\left\lbrack f(\alpha) \right\rbrack$
using Gauss-Laguerre quadrature.

**Examples:**

``` r
# E[alpha] = a/b
a <- 2; b <- 0.5
E_alpha <- integrate_gamma(identity, a, b)
cat("E[alpha] via quadrature:", round(E_alpha, 6), "\n")
#> E[alpha] via quadrature: 4
cat("E[alpha] exact (a/b):", a/b, "\n")
#> E[alpha] exact (a/b): 4
```

### 6.4 Log-Space Numerical Functions

#### `logsumexp()`

**Binary Log-Sum-Exp**

``` r
logsumexp(a, b)
```

Computes $\log\left( \exp(a) + \exp(b) \right)$ stably.

#### `logsumexp_vec()`

**Vectorized Log-Sum-Exp**

``` r
logsumexp_vec(x)
```

Computes $\log\left( \sum_{i}\exp\left( x_{i} \right) \right)$ for a
vector.

#### `softmax()`

**Numerically Stable Softmax**

``` r
softmax(x)
```

Returns probability vector from log-odds.

### 6.5 Rising Factorial

#### `log_rising_factorial()`

``` r
log_rising_factorial(alpha, J)
```

Computes
$\log(\alpha)_{J} = \log\Gamma(\alpha + J) - \log\Gamma(\alpha)$.

### 6.6 KL Divergence Functions

#### `kl_divergence_pmf()`

**KL Divergence Between PMFs**

``` r
kl_divergence_pmf(p, q, eps = 1e-15)
```

#### `kl_divergence_K()`

**KL Divergence for $K_{J}$ Distributions**

``` r
kl_divergence_K(target_pmf, a, b, J, M = 80L)
```

#### `discretize_chisq()`

**Create Target PMF from Chi-Square**

``` r
discretize_chisq(J, df, scale = 1)
```

Discretizes a (scaled) chi-square distribution onto $\{ 1,\ldots,J\}$.

------------------------------------------------------------------------

## 7. Verification Functions

The package includes extensive verification functions for development
and testing. These are typically not needed by end users but are
available for those who wish to validate computations.

### 7.1 Module Verification Functions

| Function                                    | Description                                  |
|---------------------------------------------|----------------------------------------------|
| `validate_stirling(logS)`                   | Verify Stirling numbers against known values |
| `verify_stirling_row_sum(logS)`             | Verify $\sum_{k}\left| s(J,k) \right| = J!$  |
| `verify_underdispersion(J, alpha_values)`   | Verify $\text{Var} < \text{Mean}$            |
| `verify_pmf_moments(J, alpha, logS)`        | Verify PMF matches closed-form moments       |
| `verify_marginal_moments(J, a, b)`          | Verify marginal moment properties            |
| `verify_a1_roundtrip(fit)`                  | Verify A1 parameters recover targets         |
| `verify_a2_moment_matching(J, mu_K, var_K)` | Verify A2 achieves exact matching            |

**Examples:**

``` r
# Verify Stirling numbers
logS <- compute_log_stirling(20)
validate_stirling(logS, verbose = TRUE)
#> PASS: |s(4,2)| = 11 (expected 11)
#> PASS: |s(5,3)| = 35 (expected 35)
#> PASS: |s(6,3)| = 225 (expected 225)
#> PASS: |s(10,5)| = 269325 (expected 269325)
#> [1] TRUE

# Verify underdispersion
verify_underdispersion(50, alpha_values = c(0.5, 1, 2, 5))
#> Underdispersion verification (J=50):
#>   alpha= 0.50: E[K]=  2.9378, Var(K)=  1.7091, D=0.5818 [PASS]
#>   alpha= 1.00: E[K]=  4.4992, Var(K)=  2.8741, D=0.6388 [PASS]
#>   alpha= 2.00: E[K]=  7.0376, Var(K)=  4.5356, D=0.6445 [PASS]
#>   alpha= 5.00: E[K]= 12.4605, Var(K)=  7.3861, D=0.5928 [PASS]
```

------------------------------------------------------------------------

## 8. Function Quick Reference

### By Task

**“I want to elicit a prior based on expected clusters”**

``` r
fit <- DPprior_fit(J = 50, mu_K = 5, confidence = "medium")
```

**“I want to control both clusters AND weight concentration”**

``` r
fit_K <- DPprior_fit(J = 50, mu_K = 5, var_K = 8)
fit <- DPprior_dual(fit_K, list(prob = list(threshold = 0.5, value = 0.3)))
```

**“I want to check for dominance risk”**

``` r
prob_w1_exceeds(0.5, fit$a, fit$b)
DPprior_diagnostics(fit)
```

**“I want to compute the exact PMF of K”**

``` r
logS <- compute_log_stirling(J)
pmf <- pmf_K_marginal(J, fit$a, fit$b, logS)
```

**“I want to visualize my prior”**

``` r
plot(fit)  # Full dashboard
plot_K_prior(fit)  # K distribution only
```

### Alphabetical Index

| Function                                                                                                  | Category      | Description                                                      |
|-----------------------------------------------------------------------------------------------------------|---------------|------------------------------------------------------------------|
| [`build_gamma_quadrature()`](https://joonho112.github.io/DPprior/reference/build_gamma_quadrature.md)     | Utility       | Quadrature for Gamma                                             |
| [`cdf_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/cdf_K_given_alpha.md)               | Exact         | Conditional CDF                                                  |
| [`cdf_w1()`](https://joonho112.github.io/DPprior/reference/cdf_w1.md)                                     | Weights       | CDF of $w_{1}$                                                   |
| [`check_dominance_risk()`](https://joonho112.github.io/DPprior/reference/check_dominance_risk.md)         | Diagnostics   | Risk assessment                                                  |
| [`compute_log_stirling()`](https://joonho112.github.io/DPprior/reference/compute_log_stirling.md)         | Exact         | Stirling numbers                                                 |
| [`compute_scaling_constant()`](https://joonho112.github.io/DPprior/reference/compute_scaling_constant.md) | Utility       | Scaling $c_{J}$                                                  |
| [`confidence_to_vif()`](https://joonho112.github.io/DPprior/reference/confidence_to_vif.md)               | Utility       | Confidence to VIF                                                |
| [`cv_alpha_to_variance()`](https://joonho112.github.io/DPprior/reference/cv_alpha_to_variance.md)         | Utility       | CV to variance                                                   |
| [`cv_rho()`](https://joonho112.github.io/DPprior/reference/cv_rho.md)                                     | Weights       | CV of $\rho$                                                     |
| [`discretize_chisq()`](https://joonho112.github.io/DPprior/reference/discretize_chisq.md)                 | Utility       | Chi-square to PMF                                                |
| [`DPprior_a1()`](https://joonho112.github.io/DPprior/reference/DPprior_a1.md)                             | Approximation | A1 closed-form                                                   |
| [`DPprior_a2_kl()`](https://joonho112.github.io/DPprior/reference/DPprior_a2_kl.md)                       | Approximation | A2-KL solver                                                     |
| [`DPprior_a2_newton()`](https://joonho112.github.io/DPprior/reference/DPprior_a2_newton.md)               | Approximation | A2-MN Newton                                                     |
| [`DPprior_diagnostics()`](https://joonho112.github.io/DPprior/reference/DPprior_diagnostics.md)           | Diagnostics   | Full diagnostics                                                 |
| [`DPprior_dual()`](https://joonho112.github.io/DPprior/reference/DPprior_dual.md)                         | Core          | Dual-anchor elicitation                                          |
| [`DPprior_fit()`](https://joonho112.github.io/DPprior/reference/DPprior_fit.md)                           | Core          | Main entry point                                                 |
| [`exact_K_moments()`](https://joonho112.github.io/DPprior/reference/exact_K_moments.md)                   | Exact         | Marginal moments                                                 |
| [`gauss_laguerre_nodes()`](https://joonho112.github.io/DPprior/reference/gauss_laguerre_nodes.md)         | Utility       | Quadrature nodes                                                 |
| [`get_log_stirling()`](https://joonho112.github.io/DPprior/reference/get_log_stirling.md)                 | Exact         | Safe Stirling access                                             |
| [`get_stirling_row()`](https://joonho112.github.io/DPprior/reference/get_stirling_row.md)                 | Exact         | Row of Stirling                                                  |
| [`integrate_gamma()`](https://joonho112.github.io/DPprior/reference/integrate_gamma.md)                   | Utility       | Expectation under Gamma                                          |
| [`kl_divergence_K()`](https://joonho112.github.io/DPprior/reference/kl_divergence_K.md)                   | Utility       | KL for K                                                         |
| [`kl_divergence_pmf()`](https://joonho112.github.io/DPprior/reference/kl_divergence_pmf.md)               | Utility       | KL between PMFs                                                  |
| [`log_rising_factorial()`](https://joonho112.github.io/DPprior/reference/log_rising_factorial.md)         | Utility       | Log Pochhammer                                                   |
| [`logsumexp()`](https://joonho112.github.io/DPprior/reference/logsumexp.md)                               | Utility       | Binary log-sum-exp                                               |
| [`logsumexp_vec()`](https://joonho112.github.io/DPprior/reference/logsumexp_vec.md)                       | Utility       | Vector log-sum-exp                                               |
| [`mean_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/mean_K_given_alpha.md)             | Exact         | Conditional mean                                                 |
| [`mean_rho()`](https://joonho112.github.io/DPprior/reference/mean_rho.md)                                 | Weights       | Marginal ${\mathbb{E}}\lbrack\rho\rbrack$                        |
| [`mean_rho_given_alpha()`](https://joonho112.github.io/DPprior/reference/mean_rho_given_alpha.md)         | Weights       | Conditional ${\mathbb{E}}\left\lbrack \rho|\alpha \right\rbrack$ |
| [`mean_w1()`](https://joonho112.github.io/DPprior/reference/mean_w1.md)                                   | Weights       | Marginal ${\mathbb{E}}\left\lbrack w_{1} \right\rbrack$          |
| [`moments_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/moments_K_given_alpha.md)       | Exact         | Both conditional moments                                         |
| `pdf_w1()`                                                                                                | Weights       | PDF of $w_{1}$                                                   |
| [`plot.DPprior_fit()`](https://joonho112.github.io/DPprior/reference/plot.DPprior_fit.md)                 | Diagnostics   | Dashboard plot                                                   |
| [`plot_alpha_prior()`](https://joonho112.github.io/DPprior/reference/plot_alpha_prior.md)                 | Diagnostics   | Alpha plot                                                       |
| [`plot_K_prior()`](https://joonho112.github.io/DPprior/reference/plot_K_prior.md)                         | Diagnostics   | K plot                                                           |
| [`plot_w1_prior()`](https://joonho112.github.io/DPprior/reference/plot_w1_prior.md)                       | Diagnostics   | $w_{1}$ plot                                                     |
| [`pmf_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/pmf_K_given_alpha.md)               | Exact         | Conditional PMF                                                  |
| [`pmf_K_marginal()`](https://joonho112.github.io/DPprior/reference/pmf_K_marginal.md)                     | Exact         | Marginal PMF                                                     |
| [`print.DPprior_fit()`](https://joonho112.github.io/DPprior/reference/print.DPprior_fit.md)               | Diagnostics   | Print method                                                     |
| [`prob_w1_exceeds()`](https://joonho112.github.io/DPprior/reference/prob_w1_exceeds.md)                   | Weights       | Tail probability                                                 |
| [`quantile_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/quantile_K_given_alpha.md)     | Exact         | Conditional quantile                                             |
| [`quantile_w1()`](https://joonho112.github.io/DPprior/reference/quantile_w1.md)                           | Weights       | Quantiles of $w_{1}$                                             |
| [`softmax()`](https://joonho112.github.io/DPprior/reference/softmax.md)                                   | Utility       | Stable softmax                                                   |
| [`summary.DPprior_fit()`](https://joonho112.github.io/DPprior/reference/summary.DPprior_fit.md)           | Diagnostics   | Summary method                                                   |
| [`summary_K_marginal()`](https://joonho112.github.io/DPprior/reference/summary_K_marginal.md)             | Exact         | Full K summary                                                   |
| [`var_K_given_alpha()`](https://joonho112.github.io/DPprior/reference/var_K_given_alpha.md)               | Exact         | Conditional variance                                             |
| [`var_rho()`](https://joonho112.github.io/DPprior/reference/var_rho.md)                                   | Weights       | Marginal $\text{Var}(\rho)$                                      |
| [`var_rho_given_alpha()`](https://joonho112.github.io/DPprior/reference/var_rho_given_alpha.md)           | Weights       | Conditional $\text{Var}\left( \rho|\alpha \right)$               |
| [`var_w1()`](https://joonho112.github.io/DPprior/reference/var_w1.md)                                     | Weights       | Marginal $\text{Var}\left( w_{1} \right)$                        |
| [`vif_to_variance()`](https://joonho112.github.io/DPprior/reference/vif_to_variance.md)                   | Utility       | VIF to variance                                                  |

------------------------------------------------------------------------

## References

Antoniak, C. E. (1974). Mixtures of Dirichlet processes with
applications to Bayesian nonparametric problems. *The Annals of
Statistics*, 2(6), 1152-1174.

Golub, G. H., & Welsch, J. H. (1969). Calculation of Gauss quadrature
rules. *Mathematics of Computation*, 23(106), 221-230.

Lee, J. (2026). Design-conditional prior elicitation for Dirichlet
process mixtures. *arXiv preprint* arXiv:2602.06301.

Lee, J., Che, J., Rabe-Hesketh, S., Feller, A., & Miratrix, L. (2025).
Improving the estimation of site-specific effects and their distribution
in multisite trials. *Journal of Educational and Behavioral Statistics*,
50(5), 731-764.

------------------------------------------------------------------------

*For additional documentation, see the [package
vignettes](https://joonho112.github.io/DPprior/articles/index.md) or
visit the [GitHub repository](https://github.com/joonho112/DPprior).*
