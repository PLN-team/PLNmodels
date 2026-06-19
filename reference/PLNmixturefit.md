# An R6 Class to represent a PLNfit in a mixture framework

The function
[`PLNmixture`](https://pln-team.github.io/PLNmodels/reference/PLNmixture.md)
produces a collection of models which are instances of object with class
`PLNmixturefit`. A `PLNmixturefit` (say, with k components) is itself a
collection of k `PLNfit`.

This class comes with a set of methods, some of them being useful for
the user: See the documentation for ...

## See also

The function
[`PLNmixture`](https://pln-team.github.io/PLNmodels/reference/PLNmixture.md),
the class
[`PLNmixturefamily`](https://pln-team.github.io/PLNmodels/reference/PLNmixturefamily.md)

## Active bindings

- `n`:

  number of samples

- `p`:

  number of dimensions of the latent space

- `k`:

  number of components

- `d`:

  number of covariates

- `components`:

  components of the mixture (PLNfits)

- `latent`:

  a matrix: values of the latent vector (Z in the model)

- `latent_pos`:

  a matrix: values of the latent position vector (Z) without covariates
  effects or offset

- `posteriorProb`:

  matrix ofposterior probability for cluster belonging

- `memberships`:

  vector for cluster index

- `mixtureParam`:

  vector of cluster proportions

- `optim_par`:

  a list with parameters useful for monitoring the optimization

- `nb_param`:

  number of parameters in the current PLN model

- `entropy_clustering`:

  Entropy of the variational distribution of the cluster (multinomial)

- `entropy_latent`:

  Entropy of the variational distribution of the latent vector
  (Gaussian)

- `entropy`:

  Full entropy of the variational distribution (latent vector +
  clustering)

- `loglik`:

  variational lower bound of the loglikelihood

- `loglik_vec`:

  element-wise variational lower bound of the loglikelihood

- `BIC`:

  variational lower bound of the BIC

- `ICL`:

  variational lower bound of the ICL (include entropy of both the
  clustering and latent distributions)

- `R_squared`:

  approximated goodness-of-fit criterion

- `criteria`:

  a vector with loglik, BIC, ICL, and number of parameters

- `model_par`:

  a list with the matrices of parameters found in the model (Theta,
  Sigma, Mu and Pi)

- `vcov_model`:

  character: the model used for the covariance (either "spherical",
  "diagonal" or "full")

- `fitted`:

  a matrix: fitted values of the observations (A in the model)

- `group_means`:

  a matrix of group mean vectors in the latent space.

## Methods

### Public methods

- [`PLNmixturefit$new()`](#method-PLNmixturefit-initialize)

- [`PLNmixturefit$optimize()`](#method-PLNmixturefit-optimize)

- [`PLNmixturefit$predict()`](#method-PLNmixturefit-predict)

- [`PLNmixturefit$plot_clustering_data()`](#method-PLNmixturefit-plot_clustering_data)

- [`PLNmixturefit$plot_clustering_pca()`](#method-PLNmixturefit-plot_clustering_pca)

- [`PLNmixturefit$postTreatment()`](#method-PLNmixturefit-postTreatment)

- [`PLNmixturefit$show()`](#method-PLNmixturefit-show)

- [`PLNmixturefit$print()`](#method-PLNmixturefit-print)

- [`PLNmixturefit$clone()`](#method-PLNmixturefit-clone)

------------------------------------------------------------------------

### `PLNmixturefit$new()`

Optimize a the

Initialize a `PLNmixturefit` model

#### Usage

    PLNmixturefit$new(
      responses,
      covariates,
      offsets,
      posteriorProb,
      formula,
      control
    )

#### Arguments

- `responses`:

  the matrix of responses common to every models

- `covariates`:

  the matrix of covariates common to every models

- `offsets`:

  the matrix of offsets common to every models

- `posteriorProb`:

  matrix ofposterior probability for cluster belonging

- `formula`:

  model formula used for fitting, extracted from the formula in the
  upper-level call

- `control`:

  a list for controlling the optimization.

------------------------------------------------------------------------

### `PLNmixturefit$optimize()`

Optimize a `PLNmixturefit` model

#### Usage

    PLNmixturefit$optimize(responses, covariates, offsets, config)

#### Arguments

- `responses`:

  the matrix of responses common to every models

- `covariates`:

  the matrix of covariates common to every models

- `offsets`:

  the matrix of offsets common to every models

- `config`:

  a list for controlling the optimization

------------------------------------------------------------------------

### `PLNmixturefit$predict()`

Predict group of new samples

#### Usage

    PLNmixturefit$predict(
      newdata,
      type = c("posterior", "response", "position"),
      prior = matrix(rep(1/self$k, self$k), nrow(newdata), self$k, byrow = TRUE),
      control = PLNmixture_param(),
      envir = parent.frame()
    )

#### Arguments

- `newdata`:

  A data frame in which to look for variables, offsets and counts with
  which to predict.

- `type`:

  The type of prediction required. The default `posterior` are posterior
  probabilities for each group , `response` is the group with maximal
  posterior probability and `latent` is the averaged latent coordinate
  (without offset and nor covariate effects), with weights equal to the
  posterior probabilities.

- `prior`:

  User-specified prior group probabilities in the new data. The default
  uses a uniform prior.

- `control`:

  a list-like structure for controlling the fit. See
  [`PLNmixture_param()`](https://pln-team.github.io/PLNmodels/reference/PLNmixture_param.md)
  for details.

- `envir`:

  Environment in which the prediction is evaluated

------------------------------------------------------------------------

### `PLNmixturefit$plot_clustering_data()`

Plot the matrix of expected mean counts (without offsets, without
covariate effects) reordered according the inferred clustering

#### Usage

    PLNmixturefit$plot_clustering_data(
      main = "Expected counts reorder by clustering",
      plot = TRUE,
      log_scale = TRUE
    )

#### Arguments

- `main`:

  character. A title for the plot. An hopefully appropriate title will
  be used by default.

- `plot`:

  logical. Should the plot be displayed or sent back as
  [`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
  object

- `log_scale`:

  logical. Should the color scale values be log-transform before
  plotting? Default is `TRUE`.

#### Returns

a
[`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
graphic

------------------------------------------------------------------------

### `PLNmixturefit$plot_clustering_pca()`

Plot the individual map of a PCA performed on the latent coordinates,
where individuals are colored according to the memberships

#### Usage

    PLNmixturefit$plot_clustering_pca(
      main = "Clustering labels in Individual Factor Map",
      plot = TRUE
    )

#### Arguments

- `main`:

  character. A title for the plot. An hopefully appropriate title will
  be used by default.

- `plot`:

  logical. Should the plot be displayed or sent back as
  [`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
  object

#### Returns

a
[`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
graphic

------------------------------------------------------------------------

### `PLNmixturefit$postTreatment()`

Update fields after optimization

#### Usage

    PLNmixturefit$postTreatment(
      responses,
      covariates,
      offsets,
      weights,
      config_post,
      config_optim,
      nullModel
    )

#### Arguments

- `responses`:

  the matrix of responses common to every models

- `covariates`:

  the matrix of covariates common to every models

- `offsets`:

  the matrix of offsets common to every models

- `weights`:

  an optional vector of observation weights to be used in the fitting
  process.

- `config_post`:

  a list for controlling the post-treatment

- `config_optim`:

  a list for controlling the optimization during the post-treatment
  computations

- `nullModel`:

  null model used for approximate R2 computations. Defaults to a GLM
  model with same design matrix but not latent variable.

------------------------------------------------------------------------

### `PLNmixturefit$show()`

User friendly print method

#### Usage

    PLNmixturefit$show()

------------------------------------------------------------------------

### `PLNmixturefit$print()`

User friendly print method

#### Usage

    PLNmixturefit$print()

------------------------------------------------------------------------

### `PLNmixturefit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    PLNmixturefit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
