# An R6 Class to represent a PLNfit in a standard, general framework

The function
[`PLN()`](https://pln-team.github.io/PLNmodels/reference/PLN.md) fit a
model which is an instance of a object with class `PLNfit`. Objects
produced by the functions
[`PLNnetwork()`](https://pln-team.github.io/PLNmodels/reference/PLNnetwork.md),
[`PLNPCA()`](https://pln-team.github.io/PLNmodels/reference/PLNPCA.md),
[`PLNmixture()`](https://pln-team.github.io/PLNmodels/reference/PLNmixture.md)
and
[`PLNLDA()`](https://pln-team.github.io/PLNmodels/reference/PLNLDA.md)
also enjoy the methods of `PLNfit()` by inheritance.

This class comes with a set of R6 methods, some of them being useful for
the user and exported as S3 methods. See the documentation for
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`sigma()`](https://rdrr.io/r/stats/sigma.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html) and
[`standard_error()`](https://pln-team.github.io/PLNmodels/reference/standard_error.md).

Fields are accessed via active binding and cannot be changed by the
user.

## Active bindings

- `n`:

  number of samples

- `q`:

  number of dimensions of the latent space

- `p`:

  number of species

- `d`:

  number of covariates

- `nb_param`:

  number of parameters in the current PLN model

- `model_par`:

  a list with the matrices of the model parameters: B (covariates),
  Sigma (covariance), Omega (precision matrix), plus some others
  depending on the variant)

- `var_par`:

  a list with the matrices of the variational parameters: M (means) and
  S2 (variances)

- `optim_par`:

  a list with parameters useful for monitoring the optimization

- `latent`:

  a matrix: values of the latent vector (Z in the model)

- `latent_pos`:

  a matrix: values of the latent position vector (Z) without covariates
  effects or offset

- `fitted`:

  a matrix: fitted values of the observations (A in the model)

- `vcov_coef`:

  matrix of sandwich estimator of the variance-covariance of B (need
  fixed -ie known- covariance at the moment)

- `vcov_model`:

  character: the model used for the residual covariance

- `weights`:

  observational weights

- `loglik`:

  (weighted) variational lower bound of the loglikelihood

- `loglik_vec`:

  element-wise variational lower bound of the loglikelihood

- `AIC`:

  variational lower bound of the AIC

- `BIC`:

  variational lower bound of the BIC

- `entropy`:

  Entropy of the variational distribution

- `ICL`:

  variational lower bound of the ICL

- `R_squared`:

  approximated goodness-of-fit criterion

- `criteria`:

  a vector with loglik, BIC, ICL and number of parameters

## Methods

### Public methods

- [`PLNfit$new()`](#method-PLNfit-initialize)

- [`PLNfit$update()`](#method-PLNfit-update)

- [`PLNfit$optimize()`](#method-PLNfit-optimize)

- [`PLNfit$optimize_vestep()`](#method-PLNfit-optimize_vestep)

- [`PLNfit$postTreatment()`](#method-PLNfit-postTreatment)

- [`PLNfit$predict()`](#method-PLNfit-predict)

- [`PLNfit$predict_cond()`](#method-PLNfit-predict_cond)

- [`PLNfit$show()`](#method-PLNfit-show)

- [`PLNfit$print()`](#method-PLNfit-print)

- [`PLNfit$clone()`](#method-PLNfit-clone)

------------------------------------------------------------------------

### `PLNfit$new()`

Initialize a `PLNfit` model

#### Usage

    PLNfit$new(responses, covariates, offsets, weights, formula, control)

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

  a list-like structure for controlling the fit, see
  [`PLN_param()`](https://pln-team.github.io/PLNmodels/reference/PLN_param.md).

------------------------------------------------------------------------

### `PLNfit$update()`

Update a `PLNfit` object

#### Usage

    PLNfit$update(
      B = NA,
      Sigma = NA,
      Omega = NA,
      M = NA,
      S2 = NA,
      Ji = NA,
      R2 = NA,
      Z = NA,
      A = NA,
      monitoring = NA
    )

#### Arguments

- `B`:

  matrix of regression matrix

- `Sigma`:

  variance-covariance matrix of the latent variables

- `Omega`:

  precision matrix of the latent variables. Inverse of Sigma.

- `M`:

  matrix of variational parameters for the mean

- `S2`:

  matrix of variational parameters for the variance

- `Ji`:

  vector of variational lower bounds of the log-likelihoods (one value
  per sample)

- `R2`:

  approximate R^2 goodness-of-fit criterion

- `Z`:

  matrix of latent vectors (includes covariates and offset effects)

- `A`:

  matrix of fitted values

- `monitoring`:

  a list with optimization monitoring quantities

#### Returns

Update the current `PLNfit` object

------------------------------------------------------------------------

### `PLNfit$optimize()`

Call to the NLopt or TORCH optimizer and update of the relevant fields

#### Usage

    PLNfit$optimize(responses, covariates, offsets, weights, config)

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

- `config`:

  part of the `control` argument which configures the optimizer

------------------------------------------------------------------------

### `PLNfit$optimize_vestep()`

Result of one call to the VE step of the optimization procedure: optimal
variational parameters (M, S2) and corresponding log likelihood values
for fixed model parameters (Sigma, B). Intended to position new data in
the latent space.

#### Usage

    PLNfit$optimize_vestep(
      covariates,
      offsets,
      responses,
      weights,
      B = self$model_par$B,
      Omega = self$model_par$Omega,
      control = PLN_param(backend = "nlopt")
    )

#### Arguments

- `covariates`:

  design matrix (called X in the model). Will usually be extracted from
  the corresponding field in PLNfamily-class

- `offsets`:

  offset matrix (called O in the model). Will usually be extracted from
  the corresponding field in PLNfamily-class

- `responses`:

  the matrix of responses (called Y in the model). Will usually be
  extracted from the corresponding field in PLNfamily-class

- `weights`:

  an optional vector of observation weights to be used in the fitting
  process.

- `B`:

  Optional fixed value of the regression parameters

- `Omega`:

  precision matrix of the latent variables. Inverse of Sigma.

- `control`:

  a list-like structure for controlling the fit, see
  [`PLN_param()`](https://pln-team.github.io/PLNmodels/reference/PLN_param.md).

- `Sigma`:

  variance-covariance matrix of the latent variables

#### Returns

A list with three components:

- the matrix `M` of variational means,

- the matrix `S2` of variational variances

- the vector `log.lik` of (variational) log-likelihood of each new
  observation

------------------------------------------------------------------------

### `PLNfit$postTreatment()`

Update R2, fisher and std_err fields after optimization

#### Usage

    PLNfit$postTreatment(
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

- `config_post`:

  a list for controlling the post-treatments (optional bootstrap,
  jackknife, R2, etc.). See details

- `config_optim`:

  a list for controlling the optimization (optional bootstrap,
  jackknife, R2, etc.). See details

- `nullModel`:

  null model used for approximate R2 computations. Defaults to a GLM
  model with same design matrix but not latent variable.

#### Details

The list of parameters `config` controls the post-treatment processing,
with the following entries:

- jackknife boolean indicating whether jackknife should be performed to
  evaluate bias and variance of the model parameters. Default is FALSE.

- bootstrap integer indicating the number of bootstrap resamples
  generated to evaluate the variance of the model parameters. Default is
  0 (inactivated).

- variational_var boolean indicating whether variational Fisher
  information matrix should be computed to estimate the variance of the
  model parameters (highly underestimated). Default is FALSE.

- sandwich_var boolean indicating whether sandwich estimator should be
  computed to estimate the variance of the model parameters (highly
  underestimated). Default is FALSE.

- trace integer for verbosity. should be \> 1 to see output in
  post-treatments

------------------------------------------------------------------------

### `PLNfit$predict()`

Predict position, scores or observations of new data.

#### Usage

    PLNfit$predict(
      newdata,
      responses = NULL,
      type = c("link", "response"),
      level = 1,
      envir = parent.frame()
    )

#### Arguments

- `newdata`:

  A data frame in which to look for variables with which to predict. If
  omitted, the fitted values are used.

- `responses`:

  Optional data frame containing the count of the observed variables
  (matching the names of the provided as data in the PLN function),
  assuming the interest is in testing the model.

- `type`:

  Scale used for the prediction. Either `link` (default, predicted
  positions in the latent space) or `response` (predicted counts).

- `level`:

  Optional integer value the level to be used in obtaining the
  predictions. Level zero corresponds to the population predictions
  (default if `responses` is not provided) while level one (default)
  corresponds to predictions after evaluating the variational parameters
  for the new data.

- `envir`:

  Environment in which the prediction is evaluated

#### Details

Note that `level = 1` can only be used if responses are provided, as the
variational parameters can't be estimated otherwise. In the absence of
responses, `level` is ignored and the fitted values are returned

#### Returns

A matrix with predictions scores or counts.

------------------------------------------------------------------------

### `PLNfit$predict_cond()`

Predict position, scores or observations of new data, conditionally on
the observation of a (set of) variables

#### Usage

    PLNfit$predict_cond(
      newdata,
      cond_responses,
      type = c("link", "response"),
      var_par = FALSE,
      envir = parent.frame()
    )

#### Arguments

- `newdata`:

  a data frame containing the covariates of the sites where to predict

- `cond_responses`:

  a data frame containing the count of the observed variables (matching
  the names of the provided as data in the PLN function)

- `type`:

  Scale used for the prediction. Either `link` (default, predicted
  positions in the latent space) or `response` (predicted counts).

- `var_par`:

  Boolean. Should new estimations of the variational parameters of mean
  and variance be sent back, as attributes of the matrix of predictions.
  Default to `FALSE`.

- `envir`:

  Environment in which the prediction is evaluated

#### Returns

A matrix with predictions scores or counts.

------------------------------------------------------------------------

### `PLNfit$show()`

User friendly print method

#### Usage

    PLNfit$show(
      model = paste("A multivariate Poisson Lognormal fit with", self$vcov_model,
        "covariance model.\n")
    )

#### Arguments

- `model`:

  First line of the print output

------------------------------------------------------------------------

### `PLNfit$print()`

User friendly print method

#### Usage

    PLNfit$print()

------------------------------------------------------------------------

### `PLNfit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    PLNfit$clone(deep = FALSE)

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
