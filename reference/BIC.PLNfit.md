# Bayesian Information Criterion for a fitted PLN model

Computes the variational BIC as `loglik - 0.5 * log(n) * nb_param`
(larger is better). This follows the maximization convention used
throughout PLNmodels.

## Usage

``` r
# S3 method for class 'PLNfit'
BIC(object, ...)
```

## Arguments

- object:

  an R6 object with class
  [`PLNfit`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md)

- ...:

  additional parameters for S3 compatibility. Not used

## Value

A scalar: the variational BIC (larger is better).

## Examples

``` r
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
model <- PLN(Abundance ~ 1, data = trichoptera)
#> 
#>  Initialization...
#>  Adjusting a full covariance PLN model with nlopt optimizer
#>  Post-treatments...
#>  DONE!
BIC(model)
#> [1] -1460.429
```
