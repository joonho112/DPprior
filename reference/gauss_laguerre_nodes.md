# Gauss-Laguerre Quadrature Nodes and Weights

Computes generalized Gauss-Laguerre quadrature nodes and weights via
eigendecomposition of the Jacobi matrix. These are used to approximate
integrals of the form \\\int_0^\infty f(x) x^\beta e^{-x} dx\\.

## Usage

``` r
gauss_laguerre_nodes(M, alpha_param = 0)
```

## Arguments

- M:

  Integer; number of quadrature nodes (at least 1, typically 40-120).

- alpha_param:

  Numeric; Laguerre parameter (must be \> -1). For standard Laguerre
  polynomials, use 0. For integrating against Gamma(a, b), use
  `alpha_param = a - 1`.

## Value

A list with components:

- `nodes`:

  Numeric vector of quadrature nodes \\x_m\\.

- `weights`:

  Numeric vector of quadrature weights \\w_m\\.

- `weights_log`:

  Numeric vector of \\\log(w_m)\\ for numerical stability.

## Details

The algorithm constructs a tridiagonal Jacobi matrix \\J\\ of size \\M
\times M\\:

- Diagonal: \\a_k = 2k + 1 + \alpha\\ for \\k = 0, \ldots, M-1\\

- Off-diagonal: \\b_k = \sqrt{k(k + \alpha)}\\ for \\k = 1, \ldots,
  M-1\\

Eigendecomposition \\J = V D V^T\\ yields:

- Nodes = eigenvalues (diagonal of \\D\\)

- Weights = \\\Gamma(\alpha + 1) \cdot V\[1,:\]^2\\

The generalized Laguerre polynomials \\L_n^{(\alpha)}(x)\\ are
orthogonal with respect to the weight function \\w(x) = x^\alpha
e^{-x}\\ on \\\[0, \infty)\\.

## References

Golub, G. H., & Welsch, J. H. (1969). Calculation of Gauss Quadrature
Rules. *Mathematics of Computation*, 23(106), 221-230.

## See also

[`build_gamma_quadrature`](https://joonho112.github.io/DPprior/reference/build_gamma_quadrature.md)
for Gamma distribution integration,
[`integrate_gamma`](https://joonho112.github.io/DPprior/reference/integrate_gamma.md)
for high-level expectation computation

## Examples

``` r
# Standard Laguerre (alpha_param = 0)
quad <- gauss_laguerre_nodes(40)
sum(quad$weights)  # Should equal Gamma(1) = 1
#> [1] 1

# For Gamma(2.5, b) integration, use alpha_param = 1.5
quad <- gauss_laguerre_nodes(80, alpha_param = 1.5)
sum(quad$weights)  # Should equal Gamma(2.5)
#> [1] 1.32934
```
