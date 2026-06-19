# An R6 Class to represent a collection of PLNfit

super class for
[`PLNPCAfamily`](https://pln-team.github.io/PLNmodels/reference/PLNPCAfamily.md)
and
[`PLNnetworkfamily`](https://pln-team.github.io/PLNmodels/reference/PLNnetworkfamily.md).

## See also

[`getModel()`](https://pln-team.github.io/PLNmodels/reference/getModel.md)

## Public fields

- `responses`:

  the matrix of responses common to every models

- `covariates`:

  the matrix of covariates common to every models

- `offsets`:

  the matrix of offsets common to every models

- `weights`:

  the vector of observation weights

- `inception`:

  a [PLNfit](https://pln-team.github.io/PLNmodels/reference/PLNfit.md)
  object, obtained when no sparsifying penalty is applied.

- `models`:

  a list of
  [PLNfit](https://pln-team.github.io/PLNmodels/reference/PLNfit.md)
  object, one per penalty.

## Active bindings

- `criteria`:

  a data frame with the values of some criteria (approximated
  log-likelihood, BIC, ICL, etc.) for the collection of models / fits
  BIC and ICL are defined so that they are on the same scale as the
  model log-likelihood, i.e. with the form, loglik - 0.5 penalty

- `convergence`:

  sends back a data frame with some convergence diagnostics associated
  with the optimization process (method, optimal value, etc)

## Methods

### Public methods

- [`PLNfamily$new()`](#method-PLNfamily-initialize)

- [`PLNfamily$postTreatment()`](#method-PLNfamily-postTreatment)

- [`PLNfamily$getModel()`](#method-PLNfamily-getModel)

- [`PLNfamily$plot()`](#method-PLNfamily-plot)

- [`PLNfamily$show()`](#method-PLNfamily-show)

- [`PLNfamily$print()`](#method-PLNfamily-print)

- [`PLNfamily$clone()`](#method-PLNfamily-clone)

------------------------------------------------------------------------

### `PLNfamily$new()`

Create a new `PLNfamily` object.

#### Usage

    PLNfamily$new(responses, covariates, offsets, weights, control)

#### Arguments

- `responses`:

  the matrix of responses common to every models

- `covariates`:

  the matrix of covariates common to every models

- `offsets`:

  the matrix of offsets common to every models

- `weights`:

  the vector of observation weights

- `control`:

  list controlling the optimization and the model

#### Returns

A new `PLNfamily` object

------------------------------------------------------------------------

### `PLNfamily$postTreatment()`

Update fields after optimization

#### Usage

    PLNfamily$postTreatment(config_post, config_optim)

#### Arguments

- `config_post`:

  a list for controlling the post-treatments (optional bootstrap,
  jackknife, R2, etc.).

- `config_optim`:

  a list for controlling the optimization parameters used during
  post_treatments

------------------------------------------------------------------------

### `PLNfamily$getModel()`

Extract a model from a collection of models

#### Usage

    PLNfamily$getModel(var, index = NULL)

#### Arguments

- `var`:

  value of the parameter (`rank` for PLNPCA, `sparsity` for PLNnetwork)
  that identifies the model to be extracted from the collection. If no
  exact match is found, the model with closest parameter value is
  returned with a warning.

- `index`:

  Integer index of the model to be returned. Only the first value is
  taken into account.

#### Returns

A [`PLNfit`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md)
object

------------------------------------------------------------------------

### `PLNfamily$plot()`

Lineplot of selected criteria for all models in the collection

#### Usage

    PLNfamily$plot(criteria, reverse)

#### Arguments

- `criteria`:

  A valid model selection criteria for the collection of models.
  Includes loglik, BIC (all), ICL (PLNPCA) and pen_loglik, EBIC
  (PLNnetwork)

- `reverse`:

  A logical indicating whether to plot the value of the criteria in the
  "natural" direction (loglik - penalty) or in the "reverse" direction
  (-2 loglik + penalty). Default to FALSE, i.e use the natural
  direction, on the same scale as the log-likelihood.

#### Returns

A
[`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object

------------------------------------------------------------------------

### `PLNfamily$show()`

User friendly print method

#### Usage

    PLNfamily$show()

------------------------------------------------------------------------

### `PLNfamily$print()`

User friendly print method

#### Usage

    PLNfamily$print()

------------------------------------------------------------------------

### `PLNfamily$clone()`

The objects of this class are cloneable with this method.

#### Usage

    PLNfamily$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
