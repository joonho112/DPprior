# Compute Marginal Moments and Jacobian Simultaneously

Computes the exact marginal moments \\M_1 = E\[K_J\]\\ and \\V =
Var(K_J)\\ along with the Jacobian matrix of the moment map \\F(a,b) =
(M_1, V)\\ using score function identities.

## Usage

``` r
moments_with_jacobian(J, a, b, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- J:

  Integer; sample size (number of observations/sites).

- a:

  Numeric; shape parameter of the Gamma prior on \\\alpha\\ (\> 0).

- b:

  Numeric; rate parameter of the Gamma prior on \\\alpha\\ (\> 0).

- M:

  Integer; number of quadrature nodes (default: 80).

## Value

A named list with components:

- `mean`:

  Marginal mean \\E\[K_J\]\\

- `var`:

  Marginal variance \\Var(K_J)\\

- `jacobian`:

  2x2 Jacobian matrix with structure: \$\$J_F = \begin{bmatrix} \partial
  M_1/\partial a & \partial M_1/\partial b \\ \partial V/\partial a &
  \partial V/\partial b \end{bmatrix}\$\$

## Details

This function uses the score identity (Lee, 2026, Section 3.2,
Corollary 1) to compute exact derivatives without finite differences:
\$\$\frac{\partial}{\partial\theta} E\[f(\alpha)\] = E\[f(\alpha) \cdot
s\_\theta(\alpha)\]\$\$

The Jacobian components are computed as:

- \\\partial M_1/\partial \theta = E\[\mu_J(\alpha) \cdot
  s\_\theta(\alpha)\]\\

- \\\partial V/\partial \theta = \partial E\[v_J\]/\partial \theta +
  \partial E\[\mu_J^2\]/\partial \theta - 2 M_1 \partial M_1/\partial
  \theta\\

**Numerical Considerations:**

- The score function `s_a` contains `log(alpha)`, which causes slower
  quadrature convergence compared to moment computation.

- For verification, use higher M (e.g., 120-200).

- Very small alpha values (\< 1e-12) are handled with limiting values.

## References

Lee, J. (2026). Design-Conditional Prior Elicitation for Dirichlet
Process Mixtures. *arXiv preprint* arXiv:2602.06301.

## See also

[`exact_K_moments`](https://joonho112.github.io/DPprior/reference/exact_K_moments.md)
for moments only,
[`score_a`](https://joonho112.github.io/DPprior/reference/score_a.md),
[`score_b`](https://joonho112.github.io/DPprior/reference/score_b.md)
for score functions

## Examples

``` r
# Compute moments and Jacobian for J=50, a=2, b=1
result <- moments_with_jacobian(J = 50, a = 2.0, b = 1.0)
print(result$mean)      # E[K_J]
#> [1] 6.639693
print(result$var)       # Var(K_J)
#> [1] 12.9545
print(result$jacobian)  # 2x2 Jacobian matrix
#>           da         db
#> dM1 2.245605  -4.135585
#> dV  2.943578 -13.038228

# Use in Newton iteration
target <- c(5.0, 8.0)  # Target (E[K], Var(K))
current <- c(result$mean, result$var)
residual <- current - target
delta <- solve(result$jacobian, -residual)
```
