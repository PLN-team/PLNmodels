# An R6 Class to represent a PLNfit in a standard, general framework, with spherical residual covariance

An R6 Class to represent a PLNfit in a standard, general framework, with
spherical residual covariance

## Super class

[`PLNfit`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md) -\>
`PLNfit_spherical`

## Active bindings

- `nb_param`:

  number of parameters in the current PLN model

- `vcov_model`:

  character: the model used for the residual covariance

## Methods

### Public methods

- [`PLNfit_spherical$new()`](#method-PLNfit_spherical-initialize)

- [`PLNfit_spherical$clone()`](#method-PLNfit_spherical-clone)

Inherited methods

- [`PLNfit$optimize()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-optimize)
- [`PLNfit$optimize_vestep()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-optimize_vestep)
- [`PLNfit$postTreatment()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-postTreatment)
- [`PLNfit$predict()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-predict)
- [`PLNfit$predict_cond()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-predict_cond)
- [`PLNfit$print()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-print)
- [`PLNfit$show()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-show)
- [`PLNfit$update()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-update)

------------------------------------------------------------------------

### `PLNfit_spherical$new()`

Initialize a
[`PLNfit`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md)
model

#### Usage

    PLNfit_spherical$new(responses, covariates, offsets, weights, formula, control)

#### Arguments

- `responses`:

  the matrix of responses (called Y in the model). Will usually be
  extracted from the corresponding field in PLNfamily-class

- `covariates`:

  design matrix (called X in the model). Will usually be extracted from
  the corresponding field in PLNfamily-class

- `offsets`:

  offset matrix (called O in the model). Will usually be extracted from
  the corresponding field in PLNfamily-class

- `weights`:

  an optional vector of observation weights to be used in the fitting
  process.

- `formula`:

  model formula used for fitting, extracted from the formula in the
  upper-level call

- `control`:

  a list for controlling the optimization. See details.

------------------------------------------------------------------------

### `PLNfit_spherical$clone()`

The objects of this class are cloneable with this method.

#### Usage

    PLNfit_spherical$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
if (FALSE) { # \dontrun{
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
myPLN <- PLN(Abundance ~ 1, data = trichoptera)
class(myPLN)
print(myPLN)
} # }
```
