# An R6 Class to represent a PLNfit with a residual covariance structured by a fixed correlation matrix (e.g. a genetic relationship matrix), motivated by population genetics

Sigma = sigma2 \* (rho \* C + (1 - rho) \* I_p), where C is a fixed p x
p correlation matrix supplied by the user (`control$C`) and (sigma2,
rho) are estimated. See `GeneticCovTraits` in `src/covariance_pln.h` for
the C++ side.

## Super class

[`PLNfit`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md) -\>
`PLNfit_genpop`

## Active bindings

- `nb_param`:

  number of parameters in the current PLN model

- `vcov_model`:

  character: the model used for the residual covariance

- `gen_par`:

  a list with the two extra parameters of the genpop covariance model:
  sigma2 (variance scale) and rho (mixing weight / heritability),
  decoded from Sigma and C.

## Methods

### Public methods

- [`PLNfit_genpop$new()`](#method-PLNfit_genpop-initialize)

- [`PLNfit_genpop$optimize()`](#method-PLNfit_genpop-optimize)

- [`PLNfit_genpop$clone()`](#method-PLNfit_genpop-clone)

Inherited methods

- [`PLNfit$optimize_vestep()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-optimize_vestep)
- [`PLNfit$postTreatment()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-postTreatment)
- [`PLNfit$predict()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-predict)
- [`PLNfit$predict_cond()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-predict_cond)
- [`PLNfit$print()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-print)
- [`PLNfit$show()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-show)
- [`PLNfit$update()`](https://pln-team.github.io/PLNmodels/reference/PLNfit.html#method-update)

------------------------------------------------------------------------

### `PLNfit_genpop$new()`

Initialize a `PLNfit_genpop` model

#### Usage

    PLNfit_genpop$new(responses, covariates, offsets, weights, formula, control)

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

  a list for controlling the optimization, must include a field `C` (the
  fixed p x p correlation matrix). See details.

------------------------------------------------------------------------

### `PLNfit_genpop$optimize()`

Call to the NLopt or builtin optimizer and update of the relevant fields

#### Usage

    PLNfit_genpop$optimize(responses, covariates, offsets, weights, config)

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

### `PLNfit_genpop$clone()`

The objects of this class are cloneable with this method.

#### Usage

    PLNfit_genpop$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
if (FALSE) { # \dontrun{
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
p <- ncol(trichoptera$Abundance)
C <- 0.5^abs(outer(1:p, 1:p, "-")); diag(C) <- 1
myPLN <- PLN(Abundance ~ 1, data = trichoptera, control = PLN_param(covariance = "genpop", C = C))
class(myPLN)
print(myPLN)
} # }
```
