# An R6 Class to represent a ZIPLNfit

The function
[`ZIPLN()`](https://pln-team.github.io/PLNmodels/reference/ZIPLN.md)
fits a model which is an instance of an object with class `ZIPLNfit`.

This class comes with a set of R6 methods, some of which are useful for
the end-user and exported as S3 methods. See the documentation for
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`sigma()`](https://rdrr.io/r/stats/sigma.html),
[`predict()`](https://rdrr.io/r/stats/predict.html).

Fields are accessed via active binding and cannot be changed by the
user.

## Details

Covariates for the Zero-Inflation parameter (using a logistic regression
model) can be specified in the formula RHS using the pipe
(`~ PLN effect | ZI effect`) to separate covariates for the PLN part of
the model from those for the Zero-Inflation part. Note that different
covariates can be used for each part.

## Active bindings

- `n`:

  number of samples/sites

- `q`:

  number of dimensions of the latent space

- `p`:

  number of variables/species

- `d`:

  number of covariates in the PLN part

- `d0`:

  number of covariates in the ZI part

- `nb_param_zi`:

  number of parameters in the ZI part of the model

- `nb_param_pln`:

  number of parameters in the PLN part of the model

- `nb_param`:

  number of parameters in the ZIPLN model

- `model_par`:

  a list with the matrices of parameters found in the model (B, Sigma,
  plus some others depending on the variant)

- `var_par`:

  a list with two matrices, M and S2, which are the estimated parameters
  in the variational approximation

- `optim_par`:

  a list with parameters useful for monitoring the optimization

- `latent`:

  a matrix: values of the latent vector (Z in the model)

- `latent_pos`:

  a matrix: values of the latent position vector (Z) without covariates
  effects or offset

- `fitted`:

  a matrix: fitted values of the observations (A in the model)

- `vcov_model`:

  character: the model used for the covariance (either "spherical",
  "diagonal", "full" or "sparse")

- `zi_model`:

  character: the model used for the zero inflation (either "single",
  "row", "col" or "covar")

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

- `entropy_ZI`:

  Entropy of the variational distribution

- `entropy_PLN`:

  Entropy of the Gaussian variational distribution in the PLN component

- `ICL`:

  variational lower bound of the ICL

- `criteria`:

  a vector with loglik, BIC, ICL and number of parameters

## Methods

### Public methods

- [`ZIPLNfit$update()`](#method-ZIPLNfit-update)

- [`ZIPLNfit$new()`](#method-ZIPLNfit-initialize)

- [`ZIPLNfit$optimize()`](#method-ZIPLNfit-optimize)

- [`ZIPLNfit$optimize_vestep()`](#method-ZIPLNfit-optimize_vestep)

- [`ZIPLNfit$predict()`](#method-ZIPLNfit-predict)

- [`ZIPLNfit$show()`](#method-ZIPLNfit-show)

- [`ZIPLNfit$print()`](#method-ZIPLNfit-print)

- [`ZIPLNfit$clone()`](#method-ZIPLNfit-clone)

------------------------------------------------------------------------

### `ZIPLNfit$update()`

Update a `ZIPLNfit` object

#### Usage

    ZIPLNfit$update(
      B = NA,
      B0 = NA,
      Pi = NA,
      Omega = NA,
      Sigma = NA,
      M = NA,
      S2 = NA,
      R = NA,
      Ji = NA,
      Z = NA,
      A = NA,
      monitoring = NA
    )

#### Arguments

- `B`:

  matrix of regression parameters in the Poisson lognormal component

- `B0`:

  matrix of regression parameters in the zero inflated component

- `Pi`:

  Zero inflated probability parameter (either scalar, row-vector,
  col-vector or matrix)

- `Omega`:

  precision matrix of the latent variables

- `Sigma`:

  covariance matrix of the latent variables

- `M`:

  matrix of mean vectors for the variational approximation

- `S2`:

  matrix of variance parameters for the variational approximation

- `R`:

  matrix of probabilities for the variational approximation

- `Ji`:

  vector of variational lower bounds of the log-likelihoods (one value
  per sample)

- `Z`:

  matrix of latent vectors (includes covariates and offset effects)

- `A`:

  matrix of fitted values

- `monitoring`:

  a list with optimization monitoring quantities

#### Returns

Update the current `ZIPLNfit` object

------------------------------------------------------------------------

### `ZIPLNfit$new()`

Initialize a `ZIPLNfit` model

#### Usage

    ZIPLNfit$new(data, control)

#### Arguments

- `data`:

  a named list used internally to carry the data matrices

- `control`:

  a list for controlling the optimization. See details.

------------------------------------------------------------------------

### `ZIPLNfit$optimize()`

Call to the Cpp optimizer and update of the relevant fields

#### Usage

    ZIPLNfit$optimize(data, control)

#### Arguments

- `data`:

  a named list used internally to carry the data matrices

- `control`:

  a list for controlling the optimization. See details.

------------------------------------------------------------------------

### `ZIPLNfit$optimize_vestep()`

Result of one call to the VE step of the optimization procedure: optimal
variational parameters (M, S2, R) and corresponding log likelihood
values for fixed model parameters (Sigma, B, B0). Intended to position
new data in the latent space.

#### Usage

    ZIPLNfit$optimize_vestep(
      data,
      B = self$model_par$B,
      B0 = self$model_par$B0,
      Omega = self$model_par$Omega,
      control = ZIPLN_param(backend = "nlopt")$config_optim
    )

#### Arguments

- `data`:

  a named list used internally to carry the data matrices

- `B`:

  Optional fixed value of the regression parameters in the PLN component

- `B0`:

  Optional fixed value of the regression parameters in the ZI component

- `Omega`:

  inverse variance-covariance matrix of the latent variables

- `control`:

  a list for controlling the optimization. See details.

#### Returns

A list with three components:

- the matrix `M` of variational means,

- the matrix `S2` of variational variances

- the matrix `R` of variational ZI probabilities

- the vector `Ji` of (variational) log-likelihood of each new
  observation

- a list `monitoring` with information about convergence status

------------------------------------------------------------------------

### `ZIPLNfit$predict()`

Predict position, scores or observations of new data. See
[`predict.ZIPLNfit()`](https://pln-team.github.io/PLNmodels/reference/predict.ZIPLNfit.md)
for the S3 method and additional details

#### Usage

    ZIPLNfit$predict(
      newdata,
      responses = NULL,
      type = c("link", "response", "deflated"),
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
  assuming the interest in in testing the model.

- `type`:

  Scale used for the prediction. Either `"link"` (default, predicted
  positions in the latent space), `"response"` (predicted average
  counts, accounting for zero-inflation) or `"deflated"` (predicted
  average counts, not accounting for zero-inflation and using only the
  PLN part of the model).

- `level`:

  Optional integer value the level to be used in obtaining the
  predictions. Level zero corresponds to the population predictions
  (default if `responses` is not provided) while level one (default)
  corresponds to predictions after evaluating the variational parameters
  for the new data.

- `envir`:

  Environment in which the prediction is evaluated

#### Returns

A matrix with predictions scores or counts.

------------------------------------------------------------------------

### `ZIPLNfit$show()`

User friendly print method

#### Usage

    ZIPLNfit$show(
      model = paste("A multivariate Zero Inflated Poisson Lognormal fit with",
        self$vcov_model, "covariance model.\n")
    )

#### Arguments

- `model`:

  First line of the print output

------------------------------------------------------------------------

### `ZIPLNfit$print()`

User friendly print method

#### Usage

    ZIPLNfit$print()

------------------------------------------------------------------------

### `ZIPLNfit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ZIPLNfit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
if (FALSE) { # \dontrun{
# See other examples in function ZIPLN
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
myPLN <- ZIPLN(Abundance ~ 1, data = trichoptera)
class(myPLN)
print(myPLN)
} # }
```
