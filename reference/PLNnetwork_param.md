# Control of PLNnetwork fit

Helper to define list of parameters to control the PLN fit. All
arguments have defaults.

## Usage

``` r
PLNnetwork_param(
  backend = c("builtin", "nlopt", "torch"),
  inception_cov = c("full", "spherical", "diagonal"),
  inception_backend = NULL,
  inception_niter = NULL,
  maxit_ve = NULL,
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

  optimization backend, either `"builtin"` (Newton, default) or
  `"nlopt"` (CCSAQ) or `"torch"`. The default combines `"builtin"` with
  `maxit_ve = 1` and `inception_niter = 5` (see `maxit_ve` and
  `inception_backend`): this consistently finds a better ELBO than plain
  `"nlopt"`, at essentially the same speed. Without a good inception,
  `"builtin"` alone (`maxit_ve = NULL`) can converge to a poor basin on
  large datasets — use `"nlopt"` if you want to opt out of the whole
  combination.

- inception_cov:

  Covariance structure used for the inception PLN: `"full"` (default),
  `"diagonal"` or `"spherical"`. Non-full structures are now fully
  supported: when `inception_cov != "full"`, the penalty grid is built
  from the empirical covariance of latent residuals `M − X·B` (a
  full-rank proxy for Σ), avoiding the broken `max_pen = 0` that
  previously occurred with diagonal/spherical.

- inception_backend:

  character or `NULL` (default, i.e. same as `backend`). Backend for the
  inception PLN only; the penalty grid models always use `backend`.
  Ignored when `inception` is supplied by the user.

- inception_niter:

  integer or `NULL`. Limits the inception PLN to at most this many
  iterations (EM iterations for `"builtin"`, function evaluations × 10
  for `"nlopt"`). Default is `5L` when `backend = "builtin"` (`NULL`,
  i.e. full convergence, otherwise): fewer iterations keep the latent
  mean M from over-converging toward the unconstrained optimum, which
  would make it harder to warm-start the sparse penalty models. Values
  above ~20 typically hurt. When `inception_cov != "full"` or
  `inception_niter` is set, the penalty grid uses the empirical residual
  covariance `crossprod(M − X·B) / n` for `max_pen`.

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

## See also

[`PLN_param()`](https://pln-team.github.io/PLNmodels/reference/PLN_param.md)
