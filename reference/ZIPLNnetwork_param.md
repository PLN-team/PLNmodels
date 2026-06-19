# Control of ZIPLNnetwork fit

Helper to define list of parameters to control the ZIPLNnetwork fit. All
arguments have defaults.

## Usage

``` r
ZIPLNnetwork_param(
  backend = c("builtin", "nlopt"),
  inception_cov = c("full", "spherical", "diagonal"),
  trace = 1,
  n_penalties = 30,
  min_ratio = 0.1,
  penalize_diagonal = TRUE,
  penalty_weights = NULL,
  config_post = list(),
  config_optim = list(),
  inception = NULL
)
```

## Arguments

- backend:

  optimization backend, either `"builtin"` (default, joint Newton on
  (M,ψ,R), combined with the partial E-step `maxit_ve = 1`) or `"nlopt"`
  (CCSAQ). `"builtin"` consistently finds a better ELBO across the
  penalty path at the cost of being slower.

- inception_cov:

  Covariance structure used for the inception PLN: `"full"` (default),
  `"diagonal"` or `"spherical"`. Non-full structures are now fully
  supported: when `inception_cov != "full"`, the penalty grid is built
  from the empirical covariance of latent residuals `M − X·B` (a
  full-rank proxy for Σ), avoiding the broken `max_pen = 0` that
  previously occurred with diagonal/spherical.

- trace:

  a integer for verbosity.

- n_penalties:

  an integer that specifies the number of values for the penalty grid
  when internally generated. Ignored when penalties is non `NULL`

- min_ratio:

  the penalty grid ranges from the minimal value that produces a sparse
  to this value multiplied by `min_ratio`. Default is 0.1.

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

list of parameters configuring the fit.

## Details

See
[`PLNnetwork_param()`](https://pln-team.github.io/PLNmodels/reference/PLNnetwork_param.md)
for a full description of the optimization parameters. Note that some
defaults values are different than those used in
[`PLNnetwork_param()`](https://pln-team.github.io/PLNmodels/reference/PLNnetwork_param.md):

- "ftol_out" (outer loop convergence tolerance the objective function)
  is set by default to 1e-6

- "maxit_out" (max number of iterations for the outer loop) is set by
  default to 50

## See also

[`PLNnetwork_param()`](https://pln-team.github.io/PLNmodels/reference/PLNnetwork_param.md)
and
[`PLN_param()`](https://pln-team.github.io/PLNmodels/reference/PLN_param.md)
