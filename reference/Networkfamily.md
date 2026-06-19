# An R6 Class to virtually represent a collection of network fits

The functions
[`PLNnetwork()`](https://pln-team.github.io/PLNmodels/reference/PLNnetwork.md)
and
[`ZIPLNnetwork()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNnetwork.md)
both produce an instance of this class, which can be thought of as a
vector of
[`PLNnetworkfit`](https://pln-team.github.io/PLNmodels/reference/PLNnetworkfit.md)s
[`ZIPLNfit_sparse`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit_sparse.md)s
(indexed by penalty parameter)

This class comes with a set of methods mostly used to compare network
fits (in terms of goodness of fit) or extract one from the family (based
on penalty parameter and/or goodness of it). See the documentation for
[`getBestModel()`](https://pln-team.github.io/PLNmodels/reference/getBestModel.md),
[`getModel()`](https://pln-team.github.io/PLNmodels/reference/getModel.md)
and
[plot()](https://pln-team.github.io/PLNmodels/reference/plot.Networkfamily.md)
for the user-facing ones.

## See also

The functions
[`PLNnetwork()`](https://pln-team.github.io/PLNmodels/reference/PLNnetwork.md),
[`ZIPLNnetwork()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNnetwork.md)
and the classes
[`PLNnetworkfit`](https://pln-team.github.io/PLNmodels/reference/PLNnetworkfit.md),
[`ZIPLNfit_sparse`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit_sparse.md)

## Super class

[`PLNfamily`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.md)
-\> `Networkfamily`

## Active bindings

- `penalties`:

  the sparsity level of the network in the successively fitted models

- `stability_path`:

  the stability path of each edge as returned by the stars procedure

- `stability`:

  mean edge stability along the penalty path

- `criteria`:

  a data frame with the values of some criteria (variational
  log-likelihood, (E)BIC, ICL and R2, stability) for the collection of
  models / fits BIC, ICL and EBIC are defined so that they are on the
  same scale as the model log-likelihood, i.e. with the form, loglik -
  0.5 penalty

## Methods

### Public methods

- [`Networkfamily$new()`](#method-Networkfamily-initialize)

- [`Networkfamily$optimize()`](#method-Networkfamily-optimize)

- [`Networkfamily$coefficient_path()`](#method-Networkfamily-coefficient_path)

- [`Networkfamily$getBestModel()`](#method-Networkfamily-getBestModel)

- [`Networkfamily$plot()`](#method-Networkfamily-plot)

- [`Networkfamily$plot_stars()`](#method-Networkfamily-plot_stars)

- [`Networkfamily$plot_objective()`](#method-Networkfamily-plot_objective)

- [`Networkfamily$show()`](#method-Networkfamily-show)

- [`Networkfamily$clone()`](#method-Networkfamily-clone)

Inherited methods

- [`PLNfamily$getModel()`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.html#method-getModel)
- [`PLNfamily$postTreatment()`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.html#method-postTreatment)
- [`PLNfamily$print()`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.html#method-print)

------------------------------------------------------------------------

### `Networkfamily$new()`

Initialize all models in the collection

#### Usage

    Networkfamily$new(penalties, data, control)

#### Arguments

- `penalties`:

  a vector of positive real number controlling the level of sparsity of
  the underlying network.

- `data`:

  a named list used internally to carry the data matrices

- `control`:

  a list for controlling the optimization.

#### Returns

Update all network fits in the family with smart starting values

------------------------------------------------------------------------

### `Networkfamily$optimize()`

Call to the C++ optimizer on all models of the collection

#### Usage

    Networkfamily$optimize(data, config)

#### Arguments

- `data`:

  a named list used internally to carry the data matrices

- `config`:

  a list for controlling the optimization.

------------------------------------------------------------------------

### `Networkfamily$coefficient_path()`

Extract the regularization path of a `Networkfamily`

#### Usage

    Networkfamily$coefficient_path(precision = TRUE, corr = TRUE)

#### Arguments

- `precision`:

  Logical. Should the regularization path be extracted from the
  precision matrix Omega (`TRUE`, default) or from the variance matrix
  Sigma (`FALSE`)

- `corr`:

  Logical. Should the matrix be transformed to (partial) correlation
  matrix before extraction? Defaults to `TRUE`

------------------------------------------------------------------------

### `Networkfamily$getBestModel()`

Extract the best network in the family according to some criteria

#### Usage

    Networkfamily$getBestModel(crit = c("BIC", "EBIC", "StARS"), stability = 0.9)

#### Arguments

- `crit`:

  character. Criterion used to perform the selection. If "StARS" is
  chosen but `$stability` field is empty, will compute stability path.

- `stability`:

  Only used for "StARS" criterion. A scalar indicating the target
  stability (= 1 - 2 beta) at which the network is selected. Default is
  `0.9`.

#### Details

For BIC and EBIC criteria, higher is better.

------------------------------------------------------------------------

### `Networkfamily$plot()`

Display various outputs (goodness-of-fit criteria, robustness,
diagnostic) associated with a collection of network fits (a
`Networkfamily`)

#### Usage

    Networkfamily$plot(
      criteria = c("loglik", "pen_loglik", "BIC", "EBIC"),
      reverse = FALSE,
      log.x = TRUE
    )

#### Arguments

- `criteria`:

  vector of characters. The criteria to plot in
  `c("loglik", "pen_loglik", "BIC", "EBIC")`. Defaults to all of them.

- `reverse`:

  A logical indicating whether to plot the value of the criteria in the
  "natural" direction (loglik - 0.5 penalty) or in the "reverse"
  direction (-2 loglik + penalty). Default to FALSE, i.e use the natural
  direction, on the same scale as the log-likelihood.

- `log.x`:

  logical: should the x-axis be represented in log-scale? Default is
  `TRUE`.

#### Returns

a
[`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
graph

------------------------------------------------------------------------

### `Networkfamily$plot_stars()`

Plot stability path

#### Usage

    Networkfamily$plot_stars(stability = 0.9, log.x = TRUE)

#### Arguments

- `stability`:

  scalar: the targeted level of stability using stability selection.
  Default is `0.9`.

- `log.x`:

  logical: should the x-axis be represented in log-scale? Default is
  `TRUE`.

#### Returns

a
[`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
graph

------------------------------------------------------------------------

### `Networkfamily$plot_objective()`

Plot objective value of the optimization problem along the penalty path

#### Usage

    Networkfamily$plot_objective()

#### Returns

a
[`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
graph

------------------------------------------------------------------------

### `Networkfamily$show()`

User friendly print method

#### Usage

    Networkfamily$show()

------------------------------------------------------------------------

### `Networkfamily$clone()`

The objects of this class are cloneable with this method.

#### Usage

    Networkfamily$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
