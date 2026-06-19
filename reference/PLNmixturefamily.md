# An R6 Class to represent a collection of PLNmixturefit

The function
[`PLNmixture()`](https://pln-team.github.io/PLNmodels/reference/PLNmixture.md)
produces an instance of this class.

This class comes with a set of methods, some of them being useful for
the user: See the documentation for
[`getBestModel()`](https://pln-team.github.io/PLNmodels/reference/getBestModel.md),
[`getModel()`](https://pln-team.github.io/PLNmodels/reference/getModel.md)
and
[`plot()`](https://pln-team.github.io/PLNmodels/reference/plot.PLNmixturefamily.md).

## See also

The function
[`PLNmixture`](https://pln-team.github.io/PLNmodels/reference/PLNmixture.md),
the class
[`PLNmixturefit`](https://pln-team.github.io/PLNmodels/reference/PLNmixturefit.md)

## Super class

[`PLNfamily`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.md)
-\> `PLNmixturefamily`

## Active bindings

- `clusters`:

  vector indicating the number of clusters considered is the
  successively fitted models

## Methods

### Public methods

- [`PLNmixturefamily$new()`](#method-PLNmixturefamily-initialize)

- [`PLNmixturefamily$optimize()`](#method-PLNmixturefamily-optimize)

- [`PLNmixturefamily$smooth()`](#method-PLNmixturefamily-smooth)

- [`PLNmixturefamily$plot()`](#method-PLNmixturefamily-plot)

- [`PLNmixturefamily$plot_objective()`](#method-PLNmixturefamily-plot_objective)

- [`PLNmixturefamily$getBestModel()`](#method-PLNmixturefamily-getBestModel)

- [`PLNmixturefamily$show()`](#method-PLNmixturefamily-show)

- [`PLNmixturefamily$print()`](#method-PLNmixturefamily-print)

- [`PLNmixturefamily$clone()`](#method-PLNmixturefamily-clone)

Inherited methods

- [`PLNfamily$getModel()`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.html#method-getModel)
- [`PLNfamily$postTreatment()`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.html#method-postTreatment)

------------------------------------------------------------------------

### `PLNmixturefamily$new()`

helper function for forward smoothing: split a group

Initialize all models in the collection.

#### Usage

    PLNmixturefamily$new(
      clusters,
      responses,
      covariates,
      offsets,
      formula,
      control
    )

#### Arguments

- `clusters`:

  the dimensions of the successively fitted models

- `responses`:

  the matrix of responses common to every models

- `covariates`:

  the matrix of covariates common to every models

- `offsets`:

  the matrix of offsets common to every models

- `formula`:

  model formula used for fitting, extracted from the formula in the
  upper-level call

- `control`:

  a list for controlling the optimization. See details.

- `control`:

  a list for controlling the optimization. See details.

------------------------------------------------------------------------

### `PLNmixturefamily$optimize()`

Call to the optimizer on all models of the collection

#### Usage

    PLNmixturefamily$optimize(config)

#### Arguments

- `config`:

  a list for controlling the optimization

------------------------------------------------------------------------

### `PLNmixturefamily$smooth()`

function to restart clustering to avoid local minima by smoothing the
loglikelihood values as a function of the number of clusters

#### Usage

    PLNmixturefamily$smooth(control)

#### Arguments

- `control`:

  a list to control the smoothing process

------------------------------------------------------------------------

### `PLNmixturefamily$plot()`

Lineplot of selected criteria for all models in the collection

#### Usage

    PLNmixturefamily$plot(criteria = c("loglik", "BIC", "ICL"), reverse = FALSE)

#### Arguments

- `criteria`:

  A valid model selection criteria for the collection of models. Any of
  "loglik", "BIC" or "ICL" (all).

- `reverse`:

  A logical indicating whether to plot the value of the criteria in the
  "natural" direction (loglik - 0.5 penalty) or in the "reverse"
  direction (-2 loglik + penalty). Default to FALSE, i.e use the natural
  direction, on the same scale as the log-likelihood..

#### Returns

A
[`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object

------------------------------------------------------------------------

### `PLNmixturefamily$plot_objective()`

Plot objective value of the optimization problem along the penalty path

#### Usage

    PLNmixturefamily$plot_objective()

#### Returns

a
[`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
graph

------------------------------------------------------------------------

### `PLNmixturefamily$getBestModel()`

Extract best model in the collection

#### Usage

    PLNmixturefamily$getBestModel(crit = c("BIC", "ICL", "loglik"))

#### Arguments

- `crit`:

  a character for the criterion used to performed the selection. Either
  "BIC", "ICL" or "loglik". Default is `ICL`

#### Returns

a
[`PLNmixturefit`](https://pln-team.github.io/PLNmodels/reference/PLNmixturefit.md)
object

------------------------------------------------------------------------

### `PLNmixturefamily$show()`

User friendly print method

#### Usage

    PLNmixturefamily$show()

------------------------------------------------------------------------

### `PLNmixturefamily$print()`

User friendly print method

#### Usage

    PLNmixturefamily$print()

------------------------------------------------------------------------

### `PLNmixturefamily$clone()`

The objects of this class are cloneable with this method.

#### Usage

    PLNmixturefamily$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
