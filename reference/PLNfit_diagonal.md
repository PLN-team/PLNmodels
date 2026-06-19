# An R6 Class to represent a PLNfit in a standard, general framework, with diagonal residual covariance

The function
[`PLNLDA()`](https://pln-team.github.io/PLNmodels/reference/PLNLDA.md)
produces an instance of an object with class
[`PLNLDAfit`](https://pln-team.github.io/PLNmodels/reference/PLNLDAfit.md).

This class comes with a set of methods, some of them being useful for
the user: See the documentation for the methods inherited by
[`PLNfit()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md),
the [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method for
LDA visualization and
[`predict()`](https://rdrr.io/r/stats/predict.html) method for
prediction

## Super class

[`PLNfit`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md) -\>
`PLNfit_diagonal`

## Active bindings

- `nb_param`:

  number of parameters in the current PLN model

- `vcov_model`:

  character: the model used for the residual covariance

## Methods

### Public methods

- [`PLNfit_diagonal$new()`](#method-PLNfit_diagonal-initialize)

- [`PLNfit_diagonal$clone()`](#method-PLNfit_diagonal-clone)

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

### `PLNfit_diagonal$new()`

Initialize a
[`PLNfit`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md)
model

#### Usage

    PLNfit_diagonal$new(responses, covariates, offsets, weights, formula, control)

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

### `PLNfit_diagonal$clone()`

The objects of this class are cloneable with this method.

#### Usage

    PLNfit_diagonal$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Super classes

[`PLNfit`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md) -\>
[`PLNLDAfit`](https://pln-team.github.io/PLNmodels/reference/PLNLDAfit.md)
-\> `PLNLDAfit_spherical`

## Active bindings

- `vcov_model`:

  character: the model used for the residual covariance

- `nb_param`:

  number of parameters in the current PLN model

## Methods

### Public methods

- [`PLNLDAfit_spherical$new()`](#method-PLNLDAfit_spherical-initialize)

- [`PLNLDAfit_spherical$clone()`](#method-PLNLDAfit_spherical-clone)

Inherited methods

- [`PLNfit$optimize_vestep()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-optimize_vestep)
- [`PLNfit$predict_cond()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-predict_cond)
- [`PLNfit$print()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-print)
- [`PLNfit$update()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-update)
- [`PLNLDAfit$optimize()`](https://pln-team.github.io/PLNmodels/reference/PLNLDAfit.html#method-optimize)
- [`PLNLDAfit$plot_LDA()`](https://pln-team.github.io/PLNmodels/reference/PLNLDAfit.html#method-plot_LDA)
- [`PLNLDAfit$plot_correlation_map()`](https://pln-team.github.io/PLNmodels/reference/PLNLDAfit.html#method-plot_correlation_map)
- [`PLNLDAfit$plot_individual_map()`](https://pln-team.github.io/PLNmodels/reference/PLNLDAfit.html#method-plot_individual_map)
- [`PLNLDAfit$postTreatment()`](https://pln-team.github.io/PLNmodels/reference/PLNLDAfit.html#method-postTreatment)
- [`PLNLDAfit$predict()`](https://pln-team.github.io/PLNmodels/reference/PLNLDAfit.html#method-predict)
- [`PLNLDAfit$setVisualization()`](https://pln-team.github.io/PLNmodels/reference/PLNLDAfit.html#method-setVisualization)
- [`PLNLDAfit$show()`](https://pln-team.github.io/PLNmodels/reference/PLNLDAfit.html#method-show)

------------------------------------------------------------------------

### `PLNLDAfit_spherical$new()`

Initialize a
[`PLNfit`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md)
model

#### Usage

    PLNLDAfit_spherical$new(
      grouping,
      responses,
      covariates,
      offsets,
      weights,
      formula,
      control
    )

#### Arguments

- `grouping`:

  a factor specifying the class of each observation used for
  discriminant analysis.

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

### `PLNLDAfit_spherical$clone()`

The objects of this class are cloneable with this method.

#### Usage

    PLNLDAfit_spherical$clone(deep = FALSE)

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
if (FALSE) { # \dontrun{
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
myPLNLDA <- PLNLDA(Abundance ~ 1, data = trichoptera, control = PLN_param(covariance = "spherical"))
class(myPLNLDA)
print(myPLNLDA)
} # }
```
