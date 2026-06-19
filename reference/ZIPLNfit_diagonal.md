# An R6 Class to represent a ZIPLNfit in a standard, general framework, with diagonal residual covariance

An R6 Class to represent a ZIPLNfit in a standard, general framework,
with diagonal residual covariance

## Super class

[`ZIPLNfit`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.md)
-\> `ZIPLNfit_diagonal`

## Active bindings

- `nb_param_pln`:

  number of parameters in the PLN part of the current model

- `vcov_model`:

  character: the model used for the residual covariance

## Methods

### Public methods

- [`ZIPLNfit_diagonal$new()`](#method-ZIPLNfit_diagonal-initialize)

- [`ZIPLNfit_diagonal$clone()`](#method-ZIPLNfit_diagonal-clone)

Inherited methods

- [`ZIPLNfit$optimize()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.html#method-optimize)
- [`ZIPLNfit$optimize_vestep()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.html#method-optimize_vestep)
- [`ZIPLNfit$predict()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.html#method-predict)
- [`ZIPLNfit$print()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.html#method-print)
- [`ZIPLNfit$show()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.html#method-show)
- [`ZIPLNfit$update()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.html#method-update)

------------------------------------------------------------------------

### `ZIPLNfit_diagonal$new()`

Initialize a `ZIPLNfit_diagonal` model

#### Usage

    ZIPLNfit_diagonal$new(data, control)

#### Arguments

- `data`:

  a named list used internally to carry the data matrices

- `control`:

  a list for controlling the optimization. See details.

------------------------------------------------------------------------

### `ZIPLNfit_diagonal$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ZIPLNfit_diagonal$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
if (FALSE) { # \dontrun{
# See other examples in function ZIPLN
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
myPLN <- ZIPLN(Abundance ~ 1, data = trichoptera, control = ZIPLN_param(covariance = "diagonal"))
class(myPLN)
print(myPLN)
} # }
```
