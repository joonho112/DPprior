# Poisson-Poisson KL Divergence

Computes the Kullback-Leibler divergence between two Poisson
distributions: \$\$KL(\text{Poisson}(\lambda) \|\|
\text{Poisson}(\lambda')) = \lambda \log(\lambda/\lambda') + \lambda' -
\lambda\$\$

## Usage

``` r
poisson_kl_divergence(lambda, lambda_prime)
```

## Arguments

- lambda:

  Numeric; mean of first Poisson distribution.

- lambda_prime:

  Numeric; mean of second Poisson distribution.

## Value

Numeric; KL divergence (non-negative, possibly Inf).

## Details

Special cases:

- If both \\\lambda = 0\\ and \\\lambda' = 0\\: KL = 0

- If \\\lambda = 0\\ and \\\lambda' \> 0\\: KL = \\\lambda'\\

- If \\\lambda \> 0\\ and \\\lambda' = 0\\: KL = Inf
