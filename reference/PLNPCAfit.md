# An R6 Class to represent a PLNfit in a PCA framework

The function
[`PLNPCA()`](https://pln-team.github.io/PLNmodels/reference/PLNPCA.md)
produces a collection of models which are instances of object with class
`PLNPCAfit`. This class comes with a set of methods, some of them being
useful for the user: See the documentation for the methods inherited by
[`PLNfit`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md) and
the [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods for
PCA visualization

## See also

The function
[`PLNPCA`](https://pln-team.github.io/PLNmodels/reference/PLNPCA.md),
the class
[`PLNPCAfamily`](https://pln-team.github.io/PLNmodels/reference/PLNPCAfamily.md)

## Super class

[`PLNfit`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md) -\>
`PLNPCAfit`

## Active bindings

- `var_par`:

  variational parameters (M, S2) in the rank-q latent space

- `rank`:

  the dimension of the current model

- `vcov_model`:

  character: the model used for the residual covariance

- `nb_param`:

  number of parameters in the current PLN model

- `entropy`:

  entropy of the variational distribution

- `latent_pos`:

  a matrix: values of the latent position vector (Z) without covariates
  effects or offset

- `model_par`:

  a list with the matrices associated with the estimated parameters of
  the pPCA model: B (covariates), Sigma (covariance), Omega (precision)
  and C (loadings)

- `percent_var`:

  the percent of variance explained by each axis

- `corr_circle`:

  a matrix of correlations to plot the correlation circles

- `scores`:

  a matrix of scores to plot the individual factor maps (a.k.a.
  principal components)

- `rotation`:

  a matrix of rotation of the latent space

- `eig`:

  description of the eigenvalues, similar to percent_var but for use
  with external methods

- `var`:

  a list of data frames with PCA results for the variables: `coord`
  (coordinates of the variables), `cor` (correlation between variables
  and dimensions), `cos2` (Cosine of the variables) and `contrib`
  (contributions of the variable to the axes)

- `ind`:

  a list of data frames with PCA results for the individuals: `coord`
  (coordinates of the individuals), `cos2` (Cosine of the individuals),
  `contrib` (contributions of individuals to an axis inertia) and `dist`
  (distance of individuals to the origin).

- `call`:

  Hacky binding for compatibility with factoextra functions

## Methods

### Public methods

- [`PLNPCAfit$new()`](#method-PLNPCAfit-initialize)

- [`PLNPCAfit$warm_start_from()`](#method-PLNPCAfit-warm_start_from)

- [`PLNPCAfit$update()`](#method-PLNPCAfit-update)

- [`PLNPCAfit$optimize()`](#method-PLNPCAfit-optimize)

- [`PLNPCAfit$optimize_vestep()`](#method-PLNPCAfit-optimize_vestep)

- [`PLNPCAfit$project()`](#method-PLNPCAfit-project)

- [`PLNPCAfit$setVisualization()`](#method-PLNPCAfit-setVisualization)

- [`PLNPCAfit$postTreatment()`](#method-PLNPCAfit-postTreatment)

- [`PLNPCAfit$plot_individual_map()`](#method-PLNPCAfit-plot_individual_map)

- [`PLNPCAfit$plot_correlation_circle()`](#method-PLNPCAfit-plot_correlation_circle)

- [`PLNPCAfit$plot_PCA()`](#method-PLNPCAfit-plot_PCA)

- [`PLNPCAfit$show()`](#method-PLNPCAfit-show)

- [`PLNPCAfit$clone()`](#method-PLNPCAfit-clone)

Inherited methods

- [`PLNfit$predict()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-predict)
- [`PLNfit$predict_cond()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-predict_cond)
- [`PLNfit$print()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-print)

------------------------------------------------------------------------

### `PLNPCAfit$new()`

Initialize a `PLNPCAfit` object. Uses the shared SVD from `control$svdM`
(computed once in
[`PLNPCAfamily`](https://pln-team.github.io/PLNmodels/reference/PLNPCAfamily.md))
to set the starting loadings `C` and scores `M`. The regression
coefficients `B` are initialised by the parent
[`PLNfit`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md)
constructor (LM or user-provided inception).

#### Usage

    PLNPCAfit$new(rank, responses, covariates, offsets, weights, formula, control)

#### Arguments

- `rank`:

  rank of the PCA (or equivalently, dimension of the latent space)

- `responses`:

  the matrix of responses (called Y in the model). Will usually be
  extracted from the corresponding field in
  [`PLNfamily`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.md)

- `covariates`:

  design matrix (called X in the model). Will usually be extracted from
  the corresponding field in
  [`PLNfamily`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.md)

- `offsets`:

  offset matrix (called O in the model). Will usually be extracted from
  the corresponding field in
  [`PLNfamily`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.md)

- `weights`:

  an optional vector of observation weights to be used in the fitting
  process.

- `formula`:

  model formula used for fitting, extracted from the formula in the
  upper-level call

- `control`:

  a list for controlling the optimization. See details.

------------------------------------------------------------------------

### `PLNPCAfit$warm_start_from()`

Reinitialize parameters for sequential warm-starting from a lower-rank
fit. Fitted loadings C, scores M, variances S, and regression
coefficients B from `prev_fit` are carried over; new columns are padded
using the inception SVD (C) or zeros/0.1 (M/S).

#### Usage

    PLNPCAfit$warm_start_from(prev_fit, svdM)

#### Arguments

- `prev_fit`:

  a converged `PLNPCAfit` of rank `self$rank - k` (k \>= 1)

- `svdM`:

  the inception SVD (from `PLNPCAfamily`)

------------------------------------------------------------------------

### `PLNPCAfit$update()`

Update a `PLNPCAfit` object

#### Usage

    PLNPCAfit$update(
      B = NA,
      Sigma = NA,
      Omega = NA,
      C = NA,
      M = NA,
      S2 = NA,
      Z = NA,
      A = NA,
      Ji = NA,
      R2 = NA,
      monitoring = NA
    )

#### Arguments

- `B`:

  matrix of regression matrix

- `Sigma`:

  variance-covariance matrix of the latent variables

- `Omega`:

  precision matrix of the latent variables. Inverse of Sigma.

- `C`:

  matrix of PCA loadings (in the latent space)

- `M`:

  matrix of mean vectors for the variational approximation

- `S2`:

  matrix of variational variances (n × q)

- `Z`:

  matrix of latent vectors (includes covariates and offset effects)

- `A`:

  matrix of fitted values

- `Ji`:

  vector of variational lower bounds of the log-likelihoods (one value
  per sample)

- `R2`:

  approximate R^2 goodness-of-fit criterion

- `monitoring`:

  a list with optimization monitoring quantities

#### Returns

Update the current `PLNPCAfit` object

------------------------------------------------------------------------

### `PLNPCAfit$optimize()`

Call to the C++ optimizer and update of the relevant fields

#### Usage

    PLNPCAfit$optimize(responses, covariates, offsets, weights, config)

#### Arguments

- `responses`:

  the matrix of responses (called Y in the model). Will usually be
  extracted from the corresponding field in
  [`PLNfamily`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.md)

- `covariates`:

  design matrix (called X in the model). Will usually be extracted from
  the corresponding field in
  [`PLNfamily`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.md)

- `offsets`:

  offset matrix (called O in the model). Will usually be extracted from
  the corresponding field in
  [`PLNfamily`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.md)

- `weights`:

  an optional vector of observation weights to be used in the fitting
  process.

- `config`:

  part of the `control` argument which configures the optimizer

------------------------------------------------------------------------

### `PLNPCAfit$optimize_vestep()`

Result of one call to the VE step of the optimization procedure: optimal
variational parameters (M, S) and corresponding log likelihood values
for fixed model parameters (C, B). Intended to position new data in the
latent space for further use with PCA.

#### Usage

    PLNPCAfit$optimize_vestep(
      covariates,
      offsets,
      responses,
      weights = rep(1, self$n),
      control = PLNPCA_param(backend = "nlopt")
    )

#### Arguments

- `covariates`:

  design matrix (called X in the model). Will usually be extracted from
  the corresponding field in
  [`PLNfamily`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.md)

- `offsets`:

  offset matrix (called O in the model). Will usually be extracted from
  the corresponding field in
  [`PLNfamily`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.md)

- `responses`:

  the matrix of responses (called Y in the model). Will usually be
  extracted from the corresponding field in
  [`PLNfamily`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.md)

- `weights`:

  an optional vector of observation weights to be used in the fitting
  process.

- `control`:

  a list for controlling the optimization. See details.

#### Returns

A list with three components:

- the matrix `M` of variational means,

- the matrix `S2` of variational variances

- the vector `log.lik` of (variational) log-likelihood of each new
  observation

------------------------------------------------------------------------

### `PLNPCAfit$project()`

Project new samples into the PCA space using one VE step

#### Usage

    PLNPCAfit$project(newdata, control = PLNPCA_param(), envir = parent.frame())

#### Arguments

- `newdata`:

  A data frame in which to look for variables, offsets and counts with
  which to predict.

- `control`:

  a list for controlling the optimization. See
  [`PLN()`](https://pln-team.github.io/PLNmodels/reference/PLN.md) for
  details.

- `envir`:

  Environment in which the projection is evaluated

#### Returns

- the named matrix of scores for the newdata, expressed in the same
  coordinate system as `self$scores`

------------------------------------------------------------------------

### `PLNPCAfit$setVisualization()`

Compute PCA scores in the latent space and update corresponding fields.

#### Usage

    PLNPCAfit$setVisualization(scale.unit = FALSE)

#### Arguments

- `scale.unit`:

  Logical. Should PCA scores be rescaled to have unit variance

------------------------------------------------------------------------

### `PLNPCAfit$postTreatment()`

Update R2, fisher, std_err fields and set up visualization

#### Usage

    PLNPCAfit$postTreatment(
      responses,
      covariates,
      offsets,
      weights = rep(1, nrow(responses)),
      config_post,
      config_optim,
      nullModel = NULL
    )

#### Arguments

- `responses`:

  the matrix of responses (called Y in the model). Will usually be
  extracted from the corresponding field in
  [`PLNfamily`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.md)

- `covariates`:

  design matrix (called X in the model). Will usually be extracted from
  the corresponding field in
  [`PLNfamily`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.md)

- `offsets`:

  offset matrix (called O in the model). Will usually be extracted from
  the corresponding field in
  [`PLNfamily`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.md)

- `weights`:

  an optional vector of observation weights to be used in the fitting
  process.

- `config_post`:

  a list for controlling the post-treatments (optional bootstrap,
  jackknife, R2, etc.). See details

- `config_optim`:

  a list for controlling the optimizer (either "nlopt" or "torch"
  backend). See details

- `nullModel`:

  null model used for approximate R2 computations. Defaults to a GLM
  model with same design matrix but not latent variable.

#### Details

The list of parameters `config_post` controls the post-treatment
processing, with the following entries:

- jackknife boolean indicating whether jackknife should be performed to
  evaluate bias and variance of the model parameters. Default is FALSE.

- bootstrap integer indicating the number of bootstrap resamples
  generated to evaluate the variance of the model parameters. Default is
  0 (inactivated).

- variational_var boolean indicating whether variational Fisher
  information matrix should be computed to estimate the variance of the
  model parameters (highly underestimated). Default is FALSE.

- rsquared boolean indicating whether approximation of R2 based on
  deviance should be computed. Default is TRUE

- trace integer for verbosity. should be \> 1 to see output in
  post-treatments

------------------------------------------------------------------------

### `PLNPCAfit$plot_individual_map()`

Plot the factorial map of the PCA

#### Usage

    PLNPCAfit$plot_individual_map(
      axes = 1:min(2, self$rank),
      main = "Individual Factor Map",
      plot = TRUE,
      cols = "default"
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

- `cols`:

  a character, factor or numeric to define the color associated with the
  individuals. By default, all individuals receive the default color of
  the current palette.

#### Returns

a
[`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
graphic

------------------------------------------------------------------------

### `PLNPCAfit$plot_correlation_circle()`

Plot the correlation circle of a specified axis for a
[`PLNLDAfit`](https://pln-team.github.io/PLNmodels/reference/PLNLDAfit.md)
object

#### Usage

    PLNPCAfit$plot_correlation_circle(
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

### `PLNPCAfit$plot_PCA()`

Plot a summary of the `PLNPCAfit` object

#### Usage

    PLNPCAfit$plot_PCA(
      nb_axes = min(3, self$rank),
      ind_cols = "ind_cols",
      var_cols = "var_cols",
      plot = TRUE
    )

#### Arguments

- `nb_axes`:

  scalar: the number of axes to be considered when map = "both". The
  default is min(3,rank).

- `ind_cols`:

  a character, factor or numeric to define the color associated with the
  individuals. By default, all variables receive the default color of
  the current palette.

- `var_cols`:

  a character, factor or numeric to define the color associated with the
  variables. By default, all variables receive the default color of the
  current palette.

- `plot`:

  logical. Should the plot be displayed or sent back as ggplot object

#### Returns

a `grob` object

------------------------------------------------------------------------

### `PLNPCAfit$show()`

User friendly print method

#### Usage

    PLNPCAfit$show()

------------------------------------------------------------------------

### `PLNPCAfit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    PLNPCAfit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
myPCAs <- PLNPCA(Abundance ~ 1 + offset(log(Offset)), data = trichoptera, ranks = 1:5)
#> 
#>  Initialization...
#> 
#>  Adjusting 5 PLN models for PCA analysis.
#>   Rank approximation = 1      Rank approximation = 2      Rank approximation = 3      Rank approximation = 4      Rank approximation = 5 
#>  Post-treatments
#>  DONE!
myPCA <- getBestModel(myPCAs)
class(myPCA)
#> [1] "PLNPCAfit" "PLNfit"    "PCA"       "R6"       
print(myPCA)
#> Poisson Lognormal with rank constrained for PCA (rank = 3)
#> ==================================================================
#>  nb_param   loglik      BIC      AIC      ICL
#>        65 -640.365 -766.849 -705.365 -825.034
#> ==================================================================
#> * Useful fields
#>     $model_par, $latent, $latent_pos, $var_par, $optim_par
#>     $loglik, $BIC, $ICL, $loglik_vec, $nb_param, $criteria
#> * Useful S3 methods
#>     print(), coef(), sigma(), vcov(), fitted()
#>     predict(), predict_cond(), standard_error()
#> * Additional fields for PCA
#>     $percent_var, $corr_circle, $scores, $rotation, $eig, $var, $ind
#> * Additional S3 methods for PCA
#>     plot.PLNPCAfit()
```
