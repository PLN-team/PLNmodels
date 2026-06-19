# An R6 Class to represent a ZIPLNfit in a standard, general framework, with spherical residual covariance

An R6 Class to represent a ZIPLNfit in a standard, general framework,
with spherical residual covariance

## Super class

[`ZIPLNfit`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.md)
-\> `ZIPLNfit_spherical`

## Active bindings

- `nb_param_pln`:

  number of parameters in the PLN part of the current model

- `vcov_model`:

  character: the model used for the residual covariance

## Methods

### Public methods

- [`ZIPLNfit_spherical$new()`](#method-ZIPLNfit_spherical-initialize)

- [`ZIPLNfit_spherical$clone()`](#method-ZIPLNfit_spherical-clone)

Inherited methods

- [`ZIPLNfit$optimize()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.html#method-optimize)
- [`ZIPLNfit$optimize_vestep()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.html#method-optimize_vestep)
- [`ZIPLNfit$predict()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.html#method-predict)
- [`ZIPLNfit$print()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.html#method-print)
- [`ZIPLNfit$show()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.html#method-show)
- [`ZIPLNfit$update()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.html#method-update)

------------------------------------------------------------------------

### `ZIPLNfit_spherical$new()`

Initialize a `ZIPLNfit_spherical` model

#### Usage

    ZIPLNfit_spherical$new(data, control)

#### Arguments

- `data`:

  a named list used internally to carry the data matrices

- `control`:

  a list for controlling the optimization. See details.

------------------------------------------------------------------------

### `ZIPLNfit_spherical$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ZIPLNfit_spherical$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
if (FALSE) { # \dontrun{
# See other examples in function ZIPLN
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
myPLN <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(covariance = "spherical"))
class(myPLN)
print(myPLN)
} # }
```
