# Bayesian Information Criterion for a fitted ZIPLN model

Computes the variational BIC as `loglik - 0.5 * log(n) * nb_param`
(larger is better). This follows the maximization convention used
throughout PLNmodels.

## Usage

``` r
# S3 method for class 'ZIPLNfit'
BIC(object, ...)
```

## Arguments

- object:

  an R6 object with class
  [`ZIPLNfit`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.md)

- ...:

  additional parameters for S3 compatibility. Not used

## Value

A scalar: the variational BIC (larger is better).

## Examples

``` r
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
model <- ZIPLN(Abundance ~ 1, data = trichoptera)
#> 
#>  Initialization...
#>  Adjusting a ZI-PLN model with full covariance model and single specific parameter(s) in Zero inflation component.
#>  DONE!
BIC(model)
#> [1] -1472.047
```
