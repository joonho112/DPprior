# Marginal Moments with Jacobian

Computes marginal moments and their Jacobian matrix with respect to the
Gamma hyperparameters (a, b). Used for Newton-type optimization.

## Usage

``` r
marginal_moments_with_jacobian(J, a, b, M = .QUAD_NODES_DEFAULT)
```

## Arguments

- J:

  Integer; sample size.

- a:

  Numeric; shape parameter of Gamma prior.

- b:

  Numeric; rate parameter of Gamma prior.

- M:

  Integer; number of quadrature nodes (default: 80).

## Value

A list with components:

- `mean`:

  Marginal mean

- `var`:

  Marginal variance

- `jacobian`:

  2x2 matrix of partial derivatives

## Details

The Jacobian matrix is: \$\$J = \begin{pmatrix} \partial M_1 / \partial
a & \partial M_1 / \partial b \\ \partial V / \partial a & \partial V /
\partial b \end{pmatrix}\$\$

Derivatives are computed using the score function identity:
\$\$\frac{\partial}{\partial \theta} E\[f(\alpha)\] = E\[f(\alpha) \cdot
s\_\theta(\alpha)\]\$\$

where \\s\_\theta(\alpha) = \partial \log p(\alpha \| a, b) / \partial
\theta\\.

For Gamma(a, b):

- \\s_a(\alpha) = \log(b) - \psi(a) + \log(\alpha)\\

- \\s_b(\alpha) = a/b - \alpha\\

## See also

[`exact_K_moments`](https://joonho112.github.io/DPprior/reference/exact_K_moments.md),
Module 07 (Jacobian)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- marginal_moments_with_jacobian(50, 2.0, 1.0)
result$jacobian

} # }
```
