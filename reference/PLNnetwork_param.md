# Control of PLNnetwork fit

Helper to define list of parameters to control the PLN fit. All
arguments have defaults.

## Usage

``` r
PLNnetwork_param(
  backend = c("builtin", "nlopt", "torch"),
  inception_cov = c("diagonal", "full", "spherical"),
  inception_backend = NULL,
  inception_niter = NULL,
  maxit_ve = NULL,
  trace = 1,
  n_penalties = 30,
  min_ratio = 0.1,
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

  optimization backend, either `"builtin"` (Newton, default) or
  `"nlopt"` (CCSAQ) or `"torch"`. The default combines `"builtin"` with
  `maxit_ve = 1` and `inception_niter = 5` (see `maxit_ve` and
  `inception_backend`): this consistently finds a better ELBO than plain
  `"nlopt"`, at essentially the same speed. Without a good inception,
  `"builtin"` alone (`maxit_ve = NULL`) can converge to a poor basin on
  large datasets — use `"nlopt"` if you want to opt out of the whole
  combination.

- inception_cov:

  Covariance structure used for the inception PLN, which starts the
  penalty path and sets its top: `"diagonal"` (default), `"full"` or
  `"spherical"`. The top of the grid of penalties is the largest
  off-diagonal entry of the residual covariance of the inception,
  \\crossprod(M - XB) / n\\, the penalty above which the graphical Lasso
  returns the empty network (weighted by `penalty_weights`, and on the
  diagonal too when `penalize_diagonal = TRUE`). With a diagonal
  inception and an unpenalized diagonal (the defaults), the inception is
  itself the empty network, and the path starts from it.

- inception_backend:

  character or `NULL` (default, i.e. same as `backend`). Backend for the
  inception PLN only; the penalty grid models always use `backend`.
  Ignored when `inception` is supplied by the user.

- inception_niter:

  integer or `NULL`. Limits the inception PLN to at most this many
  iterations (EM iterations for `"builtin"`, function evaluations × 10
  for `"nlopt"`). Default is `5L` when `backend = "builtin"` and
  `inception_cov = "full"`, `NULL` (full convergence) otherwise: for a
  full inception, fewer iterations keep the latent mean M from
  over-converging toward the unconstrained optimum, which would make it
  harder to warm-start the sparse penalty models (values above ~20
  typically hurt); a diagonal inception is the empty network, which is
  better converged.

- maxit_ve:

  integer or `NULL`. Maximum number of inner VE-step iterations per
  outer GLASSO alternation turn. Default is `1L` when
  `backend = "builtin"` (`NULL`, i.e. full convergence — `maxit_em` for
  `"builtin"`, `maxeval` for `"nlopt"` — otherwise). `maxit_ve = 1`
  implements a **partial E-step** (generalized EM): one Newton step per
  outer turn prevents over-convergence that causes oscillations with the
  GLASSO M-step — see `backend` for the full default combination and its
  benchmark.

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

The list of parameters `config_optim` controls the optimizers. When
"nlopt" is chosen the following entries are relevant

- "algorithm" the optimization method used by NLOPT among LD type, e.g.
  "CCSAQ", "MMA", "LBFGS". See NLOPT documentation for further details.
  Default is "CCSAQ".

- "maxeval" stop when the number of iteration exceeds maxeval. Default
  is 10000

- "ftol_rel" stop when an optimization step changes the objective
  function by less than ftol multiplied by the absolute value of the
  parameter. Default is 1e-8

- "xtol_rel" stop when an optimization step changes every parameters by
  less than xtol multiplied by the absolute value of the parameter.
  Default is 1e-6

- "ftol_abs" stop when an optimization step changes the objective
  function by less than ftol_abs. Default is 0.0 (disabled)

- "xtol_abs" stop when an optimization step changes every parameters by
  less than xtol_abs. Default is 0.0 (disabled)

- "maxtime" stop when the optimization time (in seconds) exceeds
  maxtime. Default is -1 (disabled)

- "profiled" (full covariance only) if TRUE, profile both B and Omega at
  every nlopt evaluation instead of running an EM loop (Omega fixed for
  the duration of each inner nlopt solve, B profiled in closed form at
  every evaluation). Despite the extra `O(n*p^2 + p^3)` cost per
  evaluation, benchmarks found it consistently faster than the EM loop
  (and with a slightly better loglik) across a range of problem sizes.
  Default is TRUE; set to FALSE to recover the EM loop.

When "torch" backend is used (only for PLN and PLNLDA for now), the
following entries are relevant:

- "algorithm" the optimizer used by torch among RPROP (default),
  RMSPROP, ADAM and ADAGRAD

- "maxeval" stop when the number of iteration exceeds maxeval. Default
  is 10 000

- "numepoch" stop training once this number of epochs exceeds numepoch.
  Set to -1 to enable infinite training. Default is 1 000

- "num_batch" number of batches to use during training. Defaults to 1
  (use full dataset at each epoch)

- "ftol_rel" stop when an optimization step changes the objective
  function by less than ftol multiplied by the absolute value of the
  parameter. Default is 1e-8

- "xtol_rel" stop when an optimization step changes every parameters by
  less than xtol multiplied by the absolute value of the parameter.
  Default is 1e-6

- "lr" learning rate. Default is 0.1.

- "momentum" momentum factor. Default is 0 (no momentum). Only used in
  RMSPROP

- "weight_decay" Weight decay penalty. Default is 0 (no decay). Not used
  in RPROP

- "step_sizes" pair of minimal (default: 1e-6) and maximal (default: 50)
  allowed step sizes. Only used in RPROP

- "etas" pair of multiplicative increase and decrease factors. Default
  is (0.5, 1.2). Only used in RPROP

- "centered" if TRUE, compute the centered RMSProp where the gradient is
  normalized by an estimation of its variance weight_decay (L2 penalty).
  Default to FALSE. Only used in RMSPROP

When "builtin" backend is used, the following entries are relevant

- "maxeval" stop when the number of Newton steps in the inner loop
  exceeds maxeval. Default is 10000

- "ftol_in" stop the inner loop when the objective changes by less than
  ftol_in (relative). Default is 1e-8

- "maxit_em" stop the EM outer loop when the number of EM iterations
  exceeds maxit_em. Default is 50

- "ftol_em" stop the EM outer loop when the ELBO changes by less than
  ftol_em (relative). Default is 1e-8

The list of parameters `config_post` controls the post-treatment
processing (for most `PLN*()` functions), with the following entries
(defaults may vary depending on the specific function, check
`config_post_default_*` for defaults values):

- jackknife boolean indicating whether jackknife should be performed to
  evaluate bias and variance of the model parameters. Default is FALSE.

- bootstrap integer indicating the number of bootstrap resamples
  generated to evaluate the variance of the model parameters. Default is
  0 (inactivated).

- variational_var boolean indicating whether variational Fisher
  information matrix should be computed to estimate the variance of the
  model parameters (highly underestimated). Default is FALSE.

- sandwich_var boolean indicating whether sandwich estimation should be
  used to estimate the variance of the model parameters (highly
  underestimated). Default is FALSE.

- rsquared boolean indicating whether approximation of R2 based on
  deviance should be computed. Default is TRUE

## Outer-loop optimization parameters

`PLNnetwork_param()` adds two parameters controlling the alternating
GLASSO/VEM loop:

- "ftol_em" outer alternating solver stops when the objective changes by
  less than ftol_em (relative). Default is 1e-5

- "maxit_em" outer alternating solver stops when the number of
  iterations exceeds maxit_em. Default is 20

## Scale of the penalty

With `penalty_scale = "correlation"`, the penalty on the pair \\(i, j)\\
is \\\lambda w\_{ij} \sqrt{S\_{ii} S\_{jj}}\\, recomputed at each M step
from the current residual covariance \\S\\. Without weights, the
precision matrix is \\\Omega = D^{-1/2} K D^{-1/2}\\, where \\D\\ is the
diagonal of \\S\\ and \\K\\ the graphical-Lasso estimate on the
correlation matrix \\D^{-1/2} S D^{-1/2}\\: this is the
correlation-based estimator of Rothman et al. (2008), who show that it
has a better rate of convergence in operator norm than the estimator on
the covariance scale. The grid of penalties is built on the residual
correlation of the inception, and lies between 0 and 1. Since the
weights depend on \\S\\, the alternating optimization does not maximize
a fixed penalized criterion: it looks for a fixed point.

It is the default since version 1.3.3. In simulations with a known
network (a hundred configurations of graph, sample size, dimension,
latent variances and abundances, with or without species absent from
groups of samples), the edges were recovered better on the correlation
scale in every configuration, the more so as the latent variances
differed between species, and the time was the same.

**To get the behavior of version 1.3.2 back**, use
`PLNnetwork_param(penalty_scale = "covariance")`; penalties are then on
the covariance scale. Penalties given through `penalties` are otherwise
taken on the correlation scale. Since values above 1 cannot be on that
scale (they all give the empty network), when some are given while
`penalty_scale` is left to its default they are taken as penalties on
the covariance scale and converted, with a warning: they are divided by
the ratio of the largest residual covariance of the inception to its
largest residual correlation, so that the penalty giving the empty
network on one scale gives it on the other. This is exact when all
latent variances are equal, and only a guide otherwise: the two scales
do not give the same networks. Setting `penalty_scale` explicitly, to
either value, leaves the penalties as given.

On either scale, the latent variance of a species absent from whole
groups of samples diverges along the path without `latent_floor`. On the
covariance scale, the floor stabilizes the fit without repairing the
network: the degenerate species have fewer edges, but the edges still
concentrate on them.

## References

Rothman, A. J., Bickel, P. J., Levina, E. and Zhu, J. (2008). Sparse
permutation invariant covariance estimation. *Electronic Journal of
Statistics*, 2, 494–515.
[doi:10.1214/08-EJS176](https://doi.org/10.1214/08-EJS176)

## See also

[`PLN_param()`](https://pln-team.github.io/PLNmodels/reference/PLN_param.md)
