# Component-wise standard errors of B

Extracts univariate standard errors for the estimated coefficient of B.
Standard errors are computed from the (approximate) Fisher information
matrix.

## Usage

``` r
# S3 method for class 'PLNPCAfit'
standard_error(
  object,
  type = c("variational", "jackknife", "sandwich"),
  parameter = c("B", "Omega")
)

standard_error(
  object,
  type = c("sandwich", "variational", "jackknife"),
  parameter = c("B", "Omega")
)

# S3 method for class 'PLNfit'
standard_error(
  object,
  type = c("sandwich", "variational", "jackknife", "bootstrap"),
  parameter = c("B", "Omega")
)

# S3 method for class 'PLNfit_fixedcov'
standard_error(
  object,
  type = c("sandwich", "variational", "jackknife", "bootstrap"),
  parameter = c("B", "Omega")
)

# S3 method for class 'PLNmixturefit'
standard_error(
  object,
  type = c("variational", "jackknife", "sandwich"),
  parameter = c("B", "Omega")
)

# S3 method for class 'PLNnetworkfit'
standard_error(
  object,
  type = c("variational", "jackknife", "sandwich"),
  parameter = c("B", "Omega")
)
```

## Arguments

- object:

  an R6 object with class PLNfit

- type:

  string describing the type of variance approximation: "variational",
  "jackknife", "sandwich". Default is "sandwich".

- parameter:

  string describing the target parameter: either B (regression
  coefficients) or Omega (inverse residual covariance)

## Value

A p \* d positive matrix (same size as \\B\\) with standard errors for
the coefficients of \\B\\

## Methods (by class)

- `standard_error(PLNPCAfit)`: Component-wise standard errors of B in
  [`PLNPCAfit`](https://pln-team.github.io/PLNmodels/reference/PLNPCAfit.md)
  (not implemented yet)

- `standard_error(PLNfit)`: Component-wise standard errors of B in
  [`PLNfit`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md)

- `standard_error(PLNfit_fixedcov)`: Component-wise standard errors of B
  in
  [`PLNfit_fixedcov`](https://pln-team.github.io/PLNmodels/reference/PLNfit_fixedcov.md)

- `standard_error(PLNmixturefit)`: Component-wise standard errors of B
  in
  [`PLNmixturefit`](https://pln-team.github.io/PLNmodels/reference/PLNmixturefit.md)
  (not implemented yet)

- `standard_error(PLNnetworkfit)`: Component-wise standard errors of B
  in
  [`PLNnetworkfit`](https://pln-team.github.io/PLNmodels/reference/PLNnetworkfit.md)
  (not implemented yet)

## See also

[`vcov.PLNfit()`](https://pln-team.github.io/PLNmodels/reference/vcov.PLNfit.md)
for the complete variance covariance estimation of the coefficient

## Examples

``` r
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
myPLN <- PLN(Abundance ~ 1 + offset(log(Offset)), data = trichoptera,
              control = PLN_param(config_post = list(sandwich_var = TRUE)))
#> 
#>  Initialization...
#>  Adjusting a full covariance PLN model with nlopt optimizer
#>  Post-treatments...
#>  DONE!
standard_error(myPLN)
#>                   Che       Hyc       Hym       Hys        Psy      Aga
#> (Intercept) 0.6714392 0.5698456 0.1891005 0.4850316 0.04811683 0.207243
#>                   Glo       Ath       Cea       Ced       Set       All
#> (Intercept) 0.3425401 0.3428727 0.2484439 0.1884851 0.2584335 0.2801445
#>                   Han       Hfo       Hsp       Hve       Sta
#> (Intercept) 0.3722008 0.3899489 0.3039708 0.4027737 0.1777602
```
