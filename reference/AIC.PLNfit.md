# Akaike Information Criterion for a fitted PLN model

Computes the variational AIC as `loglik - nb_param` (larger is better).
This follows the maximization convention used throughout PLNmodels.

## Usage

``` r
# S3 method for class 'PLNfit'
AIC(object, ..., k = 2)
```

## Arguments

- object:

  an R6 object with class
  [`PLNfit`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md)

- ...:

  additional parameters for S3 compatibility. Not used

- k:

  not used, present for S3 compatibility.

## Value

A scalar: the variational AIC (larger is better).

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
AIC(model)
#> [1] -1299.624
```
