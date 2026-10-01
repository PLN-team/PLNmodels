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
  penalize_diagonal = FALSE,
  penalty_weights = NULL,
  penalty_scale = c("covariance", "correlation"),
  latent_floor = 0.001,
  config_post = list(),
  config_optim = list(),
  inception = NULL
)
```

## Arguments

- backend:

  optimization backend, either `"builtin"` (default, joint Newton on
  (M,\\\psi\\,R), combined with the partial E-step `maxit_ve = 1`) or
  `"nlopt"` (CCSAQ). `"builtin"` consistently finds a better ELBO across
  the penalty path at the cost of being slower.

- inception_cov:

  Covariance structure used for the inception ZIPLN, which starts the
  penalty path and sets its top: `"full"` (default), `"diagonal"` or
  `"spherical"`. The top of the grid of penalties is computed as in
  [`PLNnetwork_param()`](https://pln-team.github.io/PLNmodels/reference/PLNnetwork_param.md).
  Unlike for
  [`PLNnetwork()`](https://pln-team.github.io/PLNmodels/reference/PLNnetwork.md),
  a full inception is the default: a diagonal one makes the top of the
  path exactly the empty network, but was found to lead to a lower BIC
  along the path, at a higher cost.

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
  graphical-Lasso? Default is `FALSE`. Penalizing the diagonal inflates
  the latent variances by the penalty (\\\Sigma\_{ii} = S\_{ii} +
  \rho\\), which the VE step then feeds back into the residual
  covariance \\S\\: along the path, the network may then never become
  empty, whatever the penalty (#180).

- penalty_weights:

  either a single or a list of p x p matrix of weights (default: all
  weights equal to 1) to adapt the amount of shrinkage to each pairs of
  node. Must be symmetric with positive values.

- penalty_scale:

  character, the scale on which the l1 penalty applies: `"covariance"`
  (default) penalizes the entries of the precision matrix as they are,
  `"correlation"` penalizes them on the scale of the variables, with a
  penalty \\\lambda \sqrt{S\_{ii} S\_{jj}}\\ on the pair \\(i, j)\\,
  where \\S\\ is the current residual covariance. This amounts to
  applying the graphical-Lasso to the residual *correlation* matrix, as
  is customary for Gaussian graphical models, and makes the penalties
  dimensionless, between 0 and 1. The entries of a precision matrix are
  not scale invariant: with `"covariance"`, a species with a large
  latent variance has nearly free edges, and one that is often absent
  but abundant when present, whose zeros are fitted by very negative
  latent means, ends up connected to most of the others (see the field
  `degenerate_species` of a
  [`PLNfit`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md)).
  `"correlation"` removes this artefact; see the section on the scale of
  the penalty.

- latent_floor:

  a positive number \\\epsilon\\ (default `1e-3`), or `NULL` for no
  floor: a floor on the variational means of the degenerate species. As
  soon as the latent variance of a species exceeds the threshold of the
  field `degenerate_species` of a
  [`PLNfit`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md)
  (100, see `options(PLNmodels.latent_variance_threshold = )`) during
  the optimization, its variational means are kept above \\\log
  \epsilon - O\\, that is \\\exp(O + M) \geq \epsilon\\, in this fit and
  in the following ones along the penalty path. The species concerned
  are in the field `floored_species` of the fits. This needs no
  preparation of the data: the species are found by the optimization
  itself. The floor restricts the variational family, not the model, and
  stops the latent variances from diverging. It leaves a fit without
  degenerate species exactly as it is.

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
