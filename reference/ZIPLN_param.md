# Control of a ZIPLN fit

Helper to define list of parameters to control the ZIPLN fit. All
arguments have defaults.

## Usage

``` r
ZIPLN_param(
  backend = c("builtin", "nlopt"),
  trace = 1,
  covariance = c("full", "diagonal", "spherical", "fixed", "sparse"),
  Omega = NULL,
  penalty = 0,
  penalize_diagonal = TRUE,
  penalty_weights = NULL,
  config_post = list(),
  config_optim = list(),
  inception = NULL
)
```

## Arguments

- backend:

  optimization backend, either `"builtin"` (default, built-in Newton
  optimizer for the joint VE step) or `"nlopt"` (NLOPT-based CCSAQ).

- trace:

  a integer for verbosity.

- covariance:

  character setting the model for the covariance matrix. Either "full",
  "diagonal", "spherical", "fixed" or "sparse". Default is "full".

- Omega:

  precision matrix of the latent variables. Inverse of Sigma. Must be
  specified if `covariance` is "fixed"

- penalty:

  a user-defined penalty to sparsify the residual covariance. Defaults
  to 0 (no sparsity).

- penalize_diagonal:

  boolean: should the diagonal terms be penalized in the
  graphical-Lasso? Default is `TRUE`

- penalty_weights:

  either a single or a list of p x p matrix of weights (default: all
  weights equal to 1) to adapt the amount of shrinkage to each pairs of
  node. Must be symmetric with positive values.

- config_post:

  a list for controlling the post-treatments (optional bootstrap,
  jackknife, R2, etc.). See details

- config_optim:

  a list for controlling the optimizer (either "nlopt" or "torch"
  backend). See details

- inception:

  Set up the parameters initialization: by default, the model is
  initialized with a multivariate linear model applied on
  log-transformed data, and with the same formula as the one provided by
  the user. However, the user can provide a PLNfit (typically obtained
  from a previous fit), which sometimes speeds up the inference.

## Value

list of parameters used during the fit and post-processing steps

## Details

See
[`PLN_param()`](https://pln-team.github.io/PLNmodels/reference/PLN_param.md)
for a description of the generic `config_optim` entries (`ftol_rel`,
`xtol_rel`, etc.). Like
[`PLNnetwork_param()`](https://pln-team.github.io/PLNmodels/reference/PLNnetwork_param.md),
ZIPLN_param() has two parameters controlling the outer EM loop:

- "ftol_out" outer solver stops when an optimization step changes the
  objective function by less than `ftol_out` multiplied by the absolute
  value of the parameter. Default is 1e-6

- "maxit_out" outer solver stops when the number of iteration exceeds
  `maxit_out`. Default is 200 for "builtin", 100 for "nlopt" and one
  additional parameter controlling the form of the variational
  approximation of the zero inflation:
