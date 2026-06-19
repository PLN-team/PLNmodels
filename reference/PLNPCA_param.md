# Control of PLNPCA fit

Helper to define list of parameters to control the PLNPCA fit. All
arguments have defaults.

## Usage

``` r
PLNPCA_param(
  backend = c("nlopt", "builtin", "torch"),
  trace = 1,
  config_optim = list(),
  config_post = list(),
  inception = NULL,
  init_method = c("LM", "GLM", "EM"),
  init_niter = 5L,
  sequential = FALSE
)
```

## Arguments

- backend:

  optimization backend, either `"nlopt"` (default, NLOPT/CCSAQ,
  recommended for PLNPCA: conservative per-variable steps reliably find
  the global basin even when the singular-value ratio d1/sqrt(n) is
  large), `"builtin"` (joint L-BFGS with strong Wolfe line search on all
  parameters simultaneously — faster per iteration but may converge to
  inferior local optima on ill-conditioned datasets), or `"torch"`
  (**experimental**: automatic differentiation via the torch package;
  tends to find lower loglik than `"nlopt"` and `"builtin"` on most
  datasets — not recommended).

- trace:

  a integer for verbosity.

- config_optim:

  a list for controlling the optimizer (either "nlopt" or "torch"
  backend). See details

- config_post:

  a list for controlling the post-treatments (optional bootstrap,
  jackknife, R2, etc.). See details

- inception:

  an optional pre-fitted
  [`PLNfit`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md)
  object. When provided, its variational means `M` and regression
  coefficients `B` are used to compute the shared SVD `svd(M - X*B)`
  that initialises all ranks simultaneously. This replaces the default
  LM-based starting point and gives the highest-quality initialisation
  for large ranks. When `NULL` (default), see `init_method`.
  `init_method` is ignored when `inception` is set. **Recommendation**:
  for large ranks (`rank > sqrt(p)`) on complex datasets, prefer
  `init_method = "EM"` (cheap, most gains) or
  `inception = PLN(formula, data)` (maximum quality, full convergence
  cost).

- init_method:

  character: strategy used to compute the starting point for the shared
  SVD.

  - `"LM"` (default): fast multivariate `lm.fit` on log-transformed
    counts. Good for small to moderate ranks (`rank <= sqrt(p)`).
    Recommended default.

  - `"GLM"`: p independent Poisson GLMs. Rarely better than `"LM"` for
    PLNPCA; not recommended.

  - `"EM"`: run `init_niter` variational EM iterations of a full PLN
    (builtin backend) before computing the SVD, giving a latent mean M
    closer to the true posterior. **Recommended for large ranks**
    (`rank > sqrt(p)`), where it noticeably improves the loglik at
    negligible extra cost. Caution: hurts for small ranks, because an
    over-converged M reduces low-rank SVD discriminability. Ignored when
    `inception` is provided.

- init_niter:

  integer: number of PLN EM iterations when `init_method = "EM"`.
  Default is 5. The sweet spot is 2–10: beyond ~10 the M over-converges
  toward the full PLN optimum, degrading the low-rank SVD and hurting
  small-rank models. Ignored for `"LM"` and `"GLM"`.

- sequential:

  logical. If `TRUE`, ranks are fitted in ascending order and each model
  is warm-started from the converged solution of the previous rank:
  loadings C are augmented with new columns from the inception SVD,
  while latent scores M and variances S2 are padded with zeros / 0.01
  respectively. Disables parallel fitting across ranks. Default is
  `FALSE`.

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
