# Akaike Information Criterion for a fitted ZIPLN model

Computes the variational AIC as `loglik - nb_param` (larger is better).
This follows the maximization convention used throughout PLNmodels.

## Usage

``` r
# S3 method for class 'ZIPLNfit'
AIC(object, ..., k = 2)
```

## Arguments

- object:

  an R6 object with class
  [`ZIPLNfit`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.md)

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
model <- ZIPLN(Abundance ~ 1, data = trichoptera)
#> 
#>  Initialization...
#>  Adjusting a ZI-PLN model with full covariance model and single specific parameter(s) in Zero inflation component.
#>  DONE!
AIC(model)
#> [1] -1310.296
```
