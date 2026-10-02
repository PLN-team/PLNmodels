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
  penalize_diagonal = FALSE,
  penalty_weights = NULL,
  penalty_scale = c("correlation", "covariance"),
  latent_floor = 0.001,
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
  to 0 (no sparsity). With the default `penalty_scale = "correlation"`,
  it applies on the scale of the correlations and lies between 0 and 1
  (a penalty above 1 gives the empty network);
  `penalty_scale = "covariance"` gives back the behavior of version
  1.3.2.

- penalize_diagonal:

  boolean: should the diagonal terms be penalized in the
  graphical-Lasso? Default is `FALSE`. Penalizing the diagonal inflates
  the latent variances (\\\Sigma\_{ii} = S\_{ii} + \rho\\ on the
  covariance scale, \\\Sigma\_{ii} = S\_{ii} (1 + \rho)\\ on the
  correlation scale), which the VE step then feeds back into the
  residual covariance \\S\\: along the path, the network may then never
  become empty, and the latent variances diverge (#180).

- penalty_weights:

  either a single or a list of p x p matrix of weights (default: all
  weights equal to 1) to adapt the amount of shrinkage to each pairs of
  node. Must be symmetric with positive values.

- penalty_scale:

  character, the scale on which the l1 penalty applies. `"correlation"`
  (default) penalizes the entries of the precision matrix on the scale
  of the variables, with a penalty \\\lambda \sqrt{S\_{ii} S\_{jj}}\\ on
  the pair \\(i, j)\\, where \\S\\ is the current residual covariance.
  This amounts to applying the graphical-Lasso to the residual
  *correlation* matrix and rescaling the result, the correlation-based
  estimator of Rothman, Bickel, Levina and Zhu (2008), and makes the
  penalties dimensionless, between 0 and 1. `"covariance"`, the only
  behavior until version 1.3.2, penalizes the entries of the precision
  matrix as they are. These are not scale invariant: a species with a
  large latent variance then has nearly free edges, and one that is
  often absent but abundant when present ends up connected to most of
  the others (see the field `degenerate_species` of a
  [`PLNfit`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md)).
  See the section on the scale of the penalty.

- latent_floor:

  a positive number (default `1e-3`), or `NULL` for no floor: a floor on
  the variational means of the degenerate species, as in
  [`PLNnetwork_param()`](https://pln-team.github.io/PLNmodels/reference/PLNnetwork_param.md).
  Only used with a sparse covariance (`penalty > 0`).

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
