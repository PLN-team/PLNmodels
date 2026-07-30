# Integrated Classification Likelihood

Generic function to compute the Integrated Classification Likelihood
(ICL) of a fitted model. ICL = BIC - entropy of the variational
distribution (larger is better).

`ICL.PLNfit`: ICL for a fitted
[`PLNfit`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md).

`ICL.ZIPLNfit`: ICL for a fitted
[`ZIPLNfit`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.md).

## Usage

``` r
ICL(object, ...)

# S3 method for class 'PLNfit'
ICL(object, ...)

# S3 method for class 'ZIPLNfit'
ICL(object, ...)
```

## Arguments

- object:

  an R6 object with class
  [`ZIPLNfit`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.md)

- ...:

  additional parameters passed to methods

## Value

A scalar: the variational ICL (larger is better).

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
ICL(model)
#> [1] -2270.936
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#> Warning: ! There are no matching names in `counts` and `covariates`.
#> ℹ Function will proceed assuming:
#> ℹ - samples are in the same order;
#> ℹ - samples are rows of `counts`.
#> Error in prepare_data(trichoptera$Abundance, trichoptera$Covariate): ✖ `counts` and `covariates` have different number of row(s):
#> ℹ `counts` has <49> row(s);
#> ℹ `covariates` has <0> row(s).
model <- ZIPLN(Abundance ~ 1, data = trichoptera)
#> 
#>  Initialization...
#>  Adjusting a ZI-PLN model with full covariance model and single specific parameter(s) in Zero inflation component.
#>  DONE!
ICL(model)
#> [1] -2365.88
```
