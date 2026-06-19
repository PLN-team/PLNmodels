# An R6 Class to represent a PLNfit in a LDA framework

The function
[`PLNLDA()`](https://pln-team.github.io/PLNmodels/reference/PLNLDA.md)
produces an instance of an object with class `PLNLDAfit`.

This class comes with a set of methods, some of them being useful for
the user: See the documentation for the methods inherited by
[`PLNfit()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md),
the [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method for
LDA visualization and
[`predict()`](https://rdrr.io/r/stats/predict.html) method for
prediction

## See also

The function
[`PLNLDA`](https://pln-team.github.io/PLNmodels/reference/PLNLDA.md).

## Super class

[`PLNfit`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md) -\>
`PLNLDAfit`

## Active bindings

- `rank`:

  the dimension of the current model

- `nb_param`:

  number of parameters in the current PLN model

- `model_par`:

  a list with the matrices associated with the estimated parameters of
  the PLN model: B (covariates), Sigma (latent covariance), C (latent
  loadings), P (latent position) and Mu (group means)

- `percent_var`:

  the percent of variance explained by each axis

- `corr_map`:

  a matrix of correlations to plot the correlation circles

- `scores`:

  a matrix of scores to plot the individual factor maps

- `group_means`:

  a matrix of group mean vectors in the latent space.

## Methods

### Public methods

- [`PLNLDAfit$new()`](#method-PLNLDAfit-initialize)

- [`PLNLDAfit$optimize()`](#method-PLNLDAfit-optimize)

- [`PLNLDAfit$postTreatment()`](#method-PLNLDAfit-postTreatment)

- [`PLNLDAfit$setVisualization()`](#method-PLNLDAfit-setVisualization)

- [`PLNLDAfit$plot_individual_map()`](#method-PLNLDAfit-plot_individual_map)

- [`PLNLDAfit$plot_correlation_map()`](#method-PLNLDAfit-plot_correlation_map)

- [`PLNLDAfit$plot_LDA()`](#method-PLNLDAfit-plot_LDA)

- [`PLNLDAfit$predict()`](#method-PLNLDAfit-predict)

- [`PLNLDAfit$show()`](#method-PLNLDAfit-show)

- [`PLNLDAfit$clone()`](#method-PLNLDAfit-clone)

Inherited methods

- [`PLNfit$optimize_vestep()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-optimize_vestep)
- [`PLNfit$predict_cond()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-predict_cond)
- [`PLNfit$print()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-print)
- [`PLNfit$update()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-update)

------------------------------------------------------------------------

### `PLNLDAfit$new()`

Initialize a `PLNLDAfit` object

#### Usage

    PLNLDAfit$new(
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

  list controlling the optimization and the model

------------------------------------------------------------------------

### `PLNLDAfit$optimize()`

Compute group means and axis of the LDA (noted B in the model) in the
latent space, update corresponding fields

#### Usage

    PLNLDAfit$optimize(grouping, responses, covariates, offsets, weights, config)

#### Arguments

- `grouping`:

  a factor specifying the class of each observation used for
  discriminant analysis.

- `responses`:

  the matrix of responses (called Y in the model). Will usually be
  extracted from the corresponding field in PLNfamily-class

- `covariates`:

  design matrix. Automatically built from the covariates and the formula
  from the call

- `offsets`:

  offset matrix (called O in the model). Will usually be extracted from
  the corresponding field in PLNfamily-class

- `weights`:

  an optional vector of observation weights to be used in the fitting
  process.

- `config`:

  list controlling the optimization

- `X`:

  Abundance matrix.

------------------------------------------------------------------------

### `PLNLDAfit$postTreatment()`

Update R2, fisher and std_err fields and visualization

#### Usage

    PLNLDAfit$postTreatment(
      grouping,
      responses,
      covariates,
      offsets,
      weights = rep(1, nrow(responses)),
      config_post,
      config_optim
    )

#### Arguments

- `grouping`:

  a factor with group memberships

- `responses`:

  the matrix of responses (counts)

- `covariates`:

  the matrix of covariates

- `offsets`:

  the matrix of offsets

- `weights`:

  an optional vector of observation weights. Default is uniform weights.

- `config_post`:

  a list for controlling the post-treatments (optional bootstrap,
  jackknife, R2, etc.).

- `config_optim`:

  list controlling the optimization parameters

------------------------------------------------------------------------

### `PLNLDAfit$setVisualization()`

Compute LDA scores in the latent space and update corresponding fields.

#### Usage

    PLNLDAfit$setVisualization(scale.unit = FALSE)

#### Arguments

- `scale.unit`:

  Logical. Should LDA scores be rescaled to have unit variance

------------------------------------------------------------------------

### `PLNLDAfit$plot_individual_map()`

Plot the factorial map of the LDA

#### Usage

    PLNLDAfit$plot_individual_map(
      axes = 1:min(2, self$rank),
      main = "Individual Factor Map",
      plot = TRUE
    )

#### Arguments

- `axes`:

  numeric, the axes to use for the plot when map = "individual" or
  "variable". Default it c(1,min(rank))

- `main`:

  character. A title for the single plot (individual or variable factor
  map). If NULL (the default), an hopefully appropriate title will be
  used.

- `plot`:

  logical. Should the plot be displayed or sent back as ggplot object

#### Returns

a
[`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
graphic

------------------------------------------------------------------------

### `PLNLDAfit$plot_correlation_map()`

Plot the correlation circle of a specified axis for a `PLNLDAfit` object

#### Usage

    PLNLDAfit$plot_correlation_map(
      axes = 1:min(2, self$rank),
      main = "Variable Factor Map",
      cols = "default",
      plot = TRUE
    )

#### Arguments

- `axes`:

  numeric, the axes to use for the plot when map = "individual" or
  "variable". Default it c(1,min(rank))

- `main`:

  character. A title for the single plot (individual or variable factor
  map). If NULL (the default), an hopefully appropriate title will be
  used.

- `cols`:

  a character, factor or numeric to define the color associated with the
  variables. By default, all variables receive the default color of the
  current palette.

- `plot`:

  logical. Should the plot be displayed or sent back as ggplot object

#### Returns

a
[`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
graphic

------------------------------------------------------------------------

### `PLNLDAfit$plot_LDA()`

Plot a summary of the `PLNLDAfit` object

#### Usage

    PLNLDAfit$plot_LDA(
      nb_axes = min(3, self$rank),
      var_cols = "default",
      plot = TRUE
    )

#### Arguments

- `nb_axes`:

  scalar: the number of axes to be considered when map = "both". The
  default is min(3,rank).

- `var_cols`:

  a character, factor or numeric to define the color associated with the
  variables. By default, all variables receive the default color of the
  current palette.

- `plot`:

  logical. Should the plot be displayed or sent back as ggplot object

#### Returns

a `grob` object

------------------------------------------------------------------------

### `PLNLDAfit$predict()`

Predict group of new samples

#### Usage

    PLNLDAfit$predict(
      newdata,
      type = c("posterior", "response", "scores"),
      scale = c("log", "prob"),
      prior = NULL,
      control = PLN_param(backend = "nlopt"),
      envir = parent.frame()
    )

#### Arguments

- `newdata`:

  A data frame in which to look for variables, offsets and counts with
  which to predict.

- `type`:

  The type of prediction required. The default are posterior
  probabilities for each group (in either unnormalized log-scale or
  natural probabilities, see "scale" for details), "response" is the
  group with maximal posterior probability and "scores" is the average
  score along each separation axis in the latent space, with weights
  equal to the posterior probabilities.

- `scale`:

  The scale used for the posterior probability. Either log-scale ("log",
  default) or natural probabilities summing up to 1 ("prob").

- `prior`:

  User-specified prior group probabilities in the new data. If NULL
  (default), prior probabilities are computed from the learning set.

- `control`:

  a list for controlling the optimization. See
  [`PLN()`](https://pln-team.github.io/PLNmodels/reference/PLN.md) for
  details.

- `envir`:

  Environment in which the prediction is evaluated

------------------------------------------------------------------------

### `PLNLDAfit$show()`

User friendly print method

#### Usage

    PLNLDAfit$show()

------------------------------------------------------------------------

### `PLNLDAfit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    PLNLDAfit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
if (FALSE) { # \dontrun{
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
myPLNLDA <- PLNLDA(Abundance ~ 1, grouping = Group, data = trichoptera)
class(myPLNLDA)
print(myPLNLDA)
} # }
```
