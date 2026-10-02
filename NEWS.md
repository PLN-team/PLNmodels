# PLNmodels (development version)

## Breaking changes in network fits

The defaults of `PLNnetwork()`, `ZIPLNnetwork()` and `ZIPLN()` with a sparse covariance
have changed since version 1.3.2. Penalty paths, criteria and selected models differ.

* **The l1 penalty now applies on the correlation scale** (`penalty_scale =
  "correlation"`), where it used to apply on the covariance scale. **The penalties
  change meaning**: they are dimensionless and lie between 0 and 1, a penalty of 1 or
  more giving the empty network. This concerns the grid built by default as well as
  the penalties given through `penalties =` (or `penalty =` in `ZIPLN_param()`).
* The variational means of the degenerate species are bounded (`latent_floor = 1e-3`).
* The diagonal of the precision matrix is no longer penalized (`penalize_diagonal =
  FALSE`), and `PLNnetwork()` starts from a diagonal inception (`inception_cov =
  "diagonal"`).

**To get the former behavior back**, set the scale in the control parameters:

```r
PLNnetwork(..., control = PLNnetwork_param(penalty_scale = "covariance"))
```

and, to undo the other changes as well,

```r
PLNnetwork_param(penalty_scale = "covariance", latent_floor = NULL,
                 penalize_diagonal = TRUE, inception_cov = "full")
```

(likewise in `ZIPLNnetwork_param()` and `ZIPLN_param()`, without `inception_cov`, whose
default has not changed there). The fits are then close to those of 1.3.2, not
identical: the graphical Lasso and the grid of penalties have changed too (see below).

**Penalties given explicitly.** When `penalties` is given while `penalty_scale` is left
to its default:

* penalties of at most 1 are taken on the correlation scale, as they are, with a
  message, once per session, recalling the change and how to undo it;
* penalties above 1 cannot be on the correlation scale. They are taken as penalties on
  the covariance scale and converted, with a warning: they are divided by the ratio of
  the largest residual covariance of the inception to its largest residual
  correlation, so that the penalty giving the empty network on one scale gives it on
  the other. This is exact when all latent variances are equal, and only a guide
  otherwise, the two scales not giving the same networks;
* in `ZIPLN_param()`, a `penalty` above 1 is not converted: a warning says that the
  network will be empty.

Setting `penalty_scale` explicitly, to either value, leaves the penalties as given and
silences these messages.

**Why.** The entries of a precision matrix are not scale invariant, so that on the
covariance scale a species with a large latent variance has nearly free edges. In
simulations with a known network (1 960 datasets, a hundred configurations of graph,
sample size, dimension, latent variances, abundances, contamination, covariates and
zero inflation), the edges were recovered better on the correlation scale in every
configuration: the F1 score at the true network size went from 0.58 to 0.74 on
uncontaminated data (from 0.45 to 0.73 when the latent variances differ between
species), and from 0.27 to 0.78 with species absent from part of the samples, for the
same computing time. The gain carries over to model selection, less strongly (BIC:
0.47 to 0.56; StARS: 0.59 to 0.71). See `inst/simus_PLNnetwork/penalty_scale/`.

## Scale of the penalty in network fits

* **New `penalty_scale`** in `PLNnetwork_param()`, `ZIPLNnetwork_param()` and
  `ZIPLN_param()`, `"correlation"` (default) or `"covariance"`. The l1 penalty of the graphical Lasso bears on the entries of the
  precision matrix, which are not scale invariant: a species with a large latent
  variance has nearly free edges. A species that is often absent but abundant when
  present, whose zeros are fitted by very negative latent means, therefore ends up
  connected to most of the others (see `$degenerate_species` below), and this happens
  well before its latent variance blows up. On the correlation scale, the penalty on
  the pair `(i, j)` is `lambda * sqrt(S_ii * S_jj)`, recomputed at each M step from the
  residual covariance `S`: this is the graphical Lasso on the residual correlation
  matrix, as is customary for Gaussian graphical models, and the penalties become
  dimensionless, between 0 and 1.
* **New `latent_floor`** in `PLNnetwork_param()`: a floor on the variational means of
  the degenerate species. As soon as the latent variance of a species exceeds the
  threshold of `$degenerate_species` (100) during the optimization, its variational
  means are kept above `log(floor) - O`, that is `exp(O + M) >= floor`, in this fit and
  in the following ones along the penalty path. The species are found by the
  optimization itself, so that no preparation of the data is needed, and are returned
  by the new field `$floored_species`. The floor restricts the variational family, not
  the model, and has no effect on a fit without degenerate species. It stops the
  latent variances from diverging, which the correlation scale alone does not prevent
  when species are absent from whole groups of samples: in simulations with such
  groups the variances diverged in 20 replicates out of 20 without the floor, on both
  scales, and in none with it. **It is on by default (`1e-3`)**; `latent_floor = NULL`
  removes it.
  A floor on every cell was tried first and dropped: on `mollusk`, where nothing
  diverges, it constrained 39 % of the cells.
* **The floor also applies to `ZIPLNnetwork()`** and to `ZIPLN()` with a sparse
  covariance (`latent_floor` in `ZIPLNnetwork_param()` and `ZIPLN_param()`, `1e-3` by
  default), which degenerate as `PLNnetwork()` does: zero inflation does not protect
  from it. On `oaks`, the largest latent variance along the default path goes from
  193 000 to 93 and the largest degree at about p edges from 111 to 21. As for
  `PLNnetwork()`, a fit without degenerate species is left exactly as it is.
* In the optimization of sparse ZIPLN fits, the best iterate is no longer chosen by an
  objective that cannot be compared across iterations: on the correlation scale, where
  the penalty weights change at each iteration, the last iterate is returned, and the
  iterates preceding the bounding of a species by the floor are discarded.
* On the covariance scale, the floor stabilizes the fit but does not repair the
  network: the degree of the hubs drops (from 107 to 22 on `oaks`), and the edges still
  concentrate on the degenerate species.
* In simulations with a known network of 40 species, 3 of which were made absent from
  65 % of the samples, these species carried 34 % (n = 200) to 93 % (n = 50) of the
  edges on the covariance scale, where 15 % was expected, and 0 to 2 % on the
  correlation scale; the F1 score of the edges between the other species, at the true
  network size, went from 0.76 to 0.85 (n = 200) and from 0.13 to 0.70 (n = 50), the
  level of uncontaminated data. Without contamination, the correlation scale did as
  well (n = 200) or better (0.71 against 0.63, n = 50). A floor alone contains the
  latent variances but not the hubs; excluding the degenerate species from the
  network moves the problem to others.
* **Defaults.** `penalty_scale = "correlation"` and `latent_floor = 1e-3` (see the
  breaking changes above). The floor changes the fits where species degenerate, and
  only there: a fit without degenerate species, as on `trichoptera`, is exactly the one
  obtained without the floor.
* The warning on degenerate species of `PLNnetwork()` now also names the species
  bounded by the floor, since the floor is what keeps their variance under the
  threshold, and says, on the covariance scale, that their edges are artefacts.
* The scripts and a summary of the exploration are in
  `inst/simus_PLNnetwork/degenerate_species/`, the account in
  `inst/devlog/DEVLOG_2026-09-30_10-01.md`.
* **`$pen_loglik` is now the criterion that the M step maximizes**, `loglik - n/2 *
  sum(abs(rho * Omega))`, where `rho` is the penalty matrix of the graphical Lasso: the
  penalty times the penalty weights and, on the correlation scale, times
  `sqrt(S_ii * S_jj)`. It was `loglik - penalty * sum(abs(Omega))`, which ignored the
  weights, the factor `n/2` and the scale, and was not the criterion being optimized.
  Its values change, on both scales, in `$criteria` and in the plots of the criteria.
* The correlation scale is the correlation-based estimator of Rothman, Bickel, Levina
  and Zhu (2008), now cited in `?PLNnetwork_param` and in the PLNnetwork vignette.
* The fits have new fields `$penalty_scale`, `$latent_floor` and `$floored_species`,
  the number of cells at the floor is in `$optim_par$n_floor`, and
  `stability_selection()` refits the subsamples with the settings of the collection.

## Graphical Lasso and network fits (#184)

* **`graphical_lasso()` now solves the problem scaled to a unit diagonal.** For any
  positive diagonal `D`, `Theta` solves the problem for `(S, rho)` if and only if
  `D^-1 Theta D^-1` solves it for `(DSD, D rho D)`: the change of variables is exact,
  and only the stopping rule sees it. On ordinary input the result is the one of
  glassoFast up to the tolerance (same supports along penalty paths, relative
  differences below 1e-4). On a covariance whose variances span orders of magnitude,
  as the residual covariance of `PLNnetwork()` does, the unscaled descent stalled in
  a limit cycle, failed in its inner loop, or returned an indefinite precision
  matrix; scaled, it converges in a few sweeps. On 91 such problems met along
  `PLNnetwork()` paths, all converged once scaled, to a lower objective; on ordinary
  ones, the median number of sweeps went from 646 to 4. The detection of limit
  cycles (`"stalled"`) is kept as a safeguard.
* **`graphical_lasso()` always returns a positive definite `wi`.** It is backed out of
  the regression coefficients of the algorithm, so that it is the inverse of `w` only
  at the exact solution, and could come out indefinite short of it, whatever the
  `status`. Its diagonal is then shifted, which keeps the network, so that its
  smallest eigenvalue is the one of the inverse of `w`; the shift is reported in the
  new `shift` component.
* **A failure along the alternating optimization of `PLNnetwork()` no longer stops the
  whole path.** A non-finite precision matrix from the graphical Lasso, or a
  non-finite objective, ends the optimization for that penalty on the previous
  iterate, with a warning; if there is none, the fit is marked as failed (`NA`
  criteria), and `getBestModel()` leaves it out of the selection, with a warning,
  instead of stopping. The reason is in `$optim_par$failure`, and the number of
  shifted precision matrices in `$optim_par$glasso_indefinite` (also for
  `ZIPLNnetwork()`). `ZIPLN()` likewise stops on a non-finite objective, on the best
  iterate so far, instead of failing with "missing value where TRUE/FALSE needed".
* On the simulations of #184 (30 multispecies Gompertz communities, `PLNnetwork(Abundance ~ 1)`),
  version 1.3.2 failed with an error on 3 of them, returned `NA` criteria on 3 more,
  and non-finite or absurd log-likelihoods (up to 1e109) on 4 others. All 30 paths now
  fit, the BIC of the selected model is higher (median +7 on the 20 datasets where
  1.3.2 did fit), and the whole study runs twice as fast.

## Penalty path of network fits (#180)

* **The diagonal of the precision matrix is no longer penalized by default**
  (`penalize_diagonal = FALSE` in `PLNnetwork_param()`, `ZIPLNnetwork_param()` and
  `ZIPLN_param()`). Penalizing it inflates the latent variances by the penalty
  (`Sigma_ii = S_ii + rho`), which the VE step feeds back into the residual
  covariance `S`: along the path, the network could then never become empty,
  whatever the penalty. On `oaks`, the top of the path had 149 edges, for a residual
  covariance 73 times the largest penalty of the grid. Set `penalize_diagonal = TRUE`
  to get the former behavior back.
* **The top of the penalty grid is now computed from the residual covariance of the
  inception**, `crossprod(M - XB) / n`, whatever its covariance model, rather than
  from its fitted `Sigma`, and `PLNnetwork()` now starts from a diagonal inception
  (`inception_cov = "diagonal"`, fully converged). With an unpenalized diagonal, the
  inception is then the empty network, and the top of the grid the penalty above
  which it is a fixed point of the alternating optimization: the path starts from
  the empty network (up to a pair on the boundary). `ZIPLNnetwork()` keeps a full
  inception, which led to a better BIC along the path.
* Penalty paths, and hence selected models, change. Against version 1.3.2, the BIC of
  the model selected by `PLNnetwork()` goes from -1390 to -1300 on `trichoptera`, from
  -5331 to -5149 on `barents`, from -38664 to -37566 on `oaks` and from -37399 to
  -36630 on `oaks` with the `tree` covariate, as fast or faster. For `ZIPLNnetwork()`,
  from -1418 to -1306, -5448 to -5206, -39786 to -37585 and -38348 to -36716, in up to
  2.5 times as long on `oaks`, where the paths are denser.
* A matrix of penalty weights that leaves no entry to penalize now stops with an
  explicit message.

## Diagnostics of structural zeros and degenerate species

* **Species with a degenerate latent variance are now reported.** A PLN model fits the
  zeros of a species that is often absent but abundant when present by sending its
  latent means to minus infinity: its latent variance blows up (thousands, where an
  ordinary one is below 30), and in a network fit the species ends up connected to
  most of the others. On `oaks` without covariate, the networks along the path of
  `PLNnetwork()` are made of little else: `f_OTU_1011`, present on all the trees of
  one type and on no other, carries 103 of the 103 edges of the first non-empty
  networks. This is not new, and does not depend on the penalization of the diagonal.
  `PLN()`, `PLNnetwork()`, `ZIPLN()` and `ZIPLNnetwork()` now raise a warning (one
  per call, for a whole collection) naming the species whose latent variance is above
  100, which the new field `$degenerate_species` of a fit returns. The threshold is
  set by `options(PLNmodels.latent_variance_threshold = )`.
* **New `structural_zeros(counts, covariates)`** finds the species absent from all the
  samples of a level of a factor covariate while present elsewhere, beyond what
  chance would explain (exact hypergeometric test, Bonferroni-corrected, so that rare
  species absent from a level by chance are not reported). `prepare_data()` runs it on
  all the factor, character and logical covariates and reports the result in a
  message; it does not remove anything. On `oaks` it finds the 8 species that live on
  one or two of the three types of trees.

## Bug fixes

* **The builtin backend now decreases the objective at every iteration** (#186). Its VE step optimizes `M` with `B` profiled (`B = P_X M`), but the objective was then evaluated at the `B` of the preceding M step, which does not match the new `M`: the reported objective could jump by orders of magnitude, and the optimization settle far   from the optimum. `B` is now updated to match `M` after the VE step, and the M step  computes `Omega` with the new `B` (the joint optimum) rather than the previous one. On data with many zeros that are not excess zeros, fits were far worse than `PLN()`'s, although ZIPLN nests PLN; they are now at least as good. On benign data (`trichoptera`, `oaks`) the log-likelihood is unchanged up to a few units.
* **The outer loop no longer stops on an increase of the objective** (#185). The
 stopping test was signed, so that any increase passed for convergence and the worse iterate was returned with the message `"converged"`. The test is now on the absolute change, the best iterate is the one returned, and the number of increases is reported in `$optim_par$objective_increases` (it should be 0). The same signed test is fixed in the VE step used by `predict()`.
* The objective monitored by a sparse fit (`ZIPLNnetwork()`, or `ZIPLN()` with a
  penalty) now includes the l1 penalty of the graphical Lasso, since the M step
  minimizes the penalized objective: the unpenalized one can legitimately increase.
* The starting value of the posterior probabilities `R` of a structural zero is now 0 where the count is positive, where the posterior is exactly 0, instead of the zero rate of the species on every cell (#186).
* `ZIPLN()` and `ZIPLNnetwork()` accept a formula stored in a variable (#187), as
  `PLN()` does: `f <- Abundance ~ 1 + Wind; ZIPLN(f, data)` failed with "cannot set an attribute on a 'symbol'", and `predict()` failed on such fits with an intercept-only formula.
* **Offsets are now robust on sparse and empty counts** #188 (`compute_offset()` and its methods):
  * `offset_rle()` keeps only species with a *positive* geometric mean in its
    median-ratio estimator (as DESeq2 does), so a zero-inflated species no longer
    drives the per-sample offset to `Inf` via a `counts / 0` ratio even when several zero-free species are available.
  * `offset_css()` falls back to TSS **only for the samples** with fewer than two
    positive counts (instead of silently switching every sample to `rowSums`), and scales the whole result by a single global median — so all offsets are on the same scale (median 1). The warning names the affected samples.
  * `compute_offset()` detects empty (all-zero) samples before dispatching, warns
    naming them, computes the offset on the non-empty samples and returns `NA` for the empty ones, instead of a `0` offset (hence `log(0) = -Inf`) or uninformative `NA/NaN/Inf` errors when several methods are called directly.
* Jackknife and bootstrap estimates of the variance (`config_post = list(jackknife =
  TRUE, bootstrap = ...)`) failed with a fixed or a `"genpop"` covariance, the fixed
  precision matrix `Omega` or correlation matrix `C` not being passed to the fits on
  the resampled data. They now work, also in `PLNnetwork()`.

# PLNmodels 1.3.2

## Internal graphical Lasso, replacing glassoFast

* `PLNnetwork()` and `ZIPLNnetwork()` now use an internal graphical Lasso, a C++ port of the GLASSOFAST algorithm (Sustik and Calderhead, 2012) shared with the
  **normalblockr** package, instead of calling `glassoFast::glassoFast()`.
  **glassoFast** is no longer a dependency (it moves to `Suggests`, for tests only).
* The motivation is robustness: glassoFast's Fortran routine could loop forever on a nearly collapsed covariance matrix (entries of order 1e-8, as a rank-deficient residual covariance produces), in compiled code that no R-level timeout could stop. The new solver always terminates (non-finite input and zero-variance coordinates are rejected, the inner coordinate descent is bounded), reports non-convergence instead of hanging, and can be interrupted from R, by the user or by `setTimeLimit()`/`R.utils::withTimeout()`, which then raise their usual error.
* On ordinary input, results are those of glassoFast up to machine precision: identical supports along penalty paths and relative differences below 1e-15 on the precision matrix, for scalar as well as weighted penalties with an unpenalized diagonal; `PLNnetwork()` and `ZIPLNnetwork()` fits with the default `"builtin"` backend are unchanged (log-likelihoods within 1e-11). Fits with `backend = "nlopt"` can differ slightly along the path: that backend is sensitive to perturbations at the level of the last floating-point digit. Speed is the same or slightly better.
* The solver is exported as `graphical_lasso(S, rho, thr, maxit, w_init, wi_init)`, returning `w`, `wi`, `niter`, `converged`, `status` and `delta`, with the same defaults as glassoFast and an optional warm start.
* **It stops when it is cycling rather than converging.** On an ill-conditioned covariance (ie rank-deficient, as `PLNnetwork()` produces when p ~ n) the sweeps settle into a small limit cycle: the convergence criterion stops decreasing and oscillates
  just above its threshold forever. glassoFast spends its whole 10000-sweep budget on these and reports success regardless. The cycle is now detected after 1000 sweeps without progress (tunable through `stall_patience`), and reported as `status = "stalled"`.
* Non-convergence of the graphical Lasso is recorded in the fits' monitoring (`$optim_par$glasso_nonconverged` and `$optim_par$glasso_stalled`, counted along the alternating optimization). A warning is now raised only for the
  numerical failures (`"degenerate"`, `"inner_failure"`, `"max_iter"`), not for a stalled solve, which says something about the problem (too weak a penalty for a nearly rank-deficient covariance) rather than about the solution.
* Behaviour change in a degenerate case: when the covariance matrix has no off-diagonal mass, the (diagonal) precision matrix is now the correct `1 / (S_ii + rho_ii)`, where glassoFast returned `1 / max(rho_ii, 1.1e-16)`.

## Model selection and display of networks

* **The EBIC of `PLNnetworkfit` and `ZIPLNfit_sparse` is now the one of Foygel and
  Drton (2010)**, `BIC - 2 gamma |E| log(p)`, which was designed for the graphical
  Lasso, instead of the Stirling approximation of Chen and Chen (2008)'s original
  criterion that was used so far. The tuning parameter is exposed as `$ebic_gamma`
  (default `0.5`, `0` gives back the BIC), on a single fit or on a whole collection.
  EBIC values change, and so may the model selected by `getBestModel("EBIC")`.
* The `$density` of a network is now `|E| / (p (p - 1) / 2)`: it divided the edge
  count by `p^2` rather than by the number of possible edges, and was thus
  understated by a factor `(p - 1) / p`.
* In the `igraph` output of `plot()` on a network fit, the opacity of an edge is now
  proportional to the strength of its partial correlation, as its width already was.
  Dense networks are no longer a solid blob of colour. The opacity of the weakest
  edge is the new `edge.alpha` argument (default `0.2`, set it to `1` to restore
  uniformly opaque edges).

## Bug fixes

* **`PLNPCA()`'s variational bound carried a spurious `p / 2` per observation**, the
  entropy constant of the full-covariance models: a rank-`q` model's variational
  distribution is over the `q`-dimensional scores, and its `q / 2` constant was
  already in the Kullback-Leibler term. `loglik`, `BIC` and `ICL` of every
  `PLNPCAfit` therefore decrease by `n * p / 2`. Rank selection is unchanged (the
  term does not depend on the rank), but a `PLNPCA()` fit is now comparable with a
  `PLN()` one, or with any other model. Reported by Nguyen Quang Huy (Actuarial
  Science Laboratory, National Economics University, Vietnam).

* **`PLNnetwork()` and `ZIPLNnetwork()` no longer fail outright when one model along
  the penalty path cannot be fitted**: the collection is truncated with a warning
  instead, and `stability_selection()` treats a model with no estimated network as
  fully unstable. A second line of defense on top of the solver above (thanks
  @rfriedman22, #176).

# PLNmodels 1.3.1

## Bug fix

* `ZIPLNnetwork()` fitted its inception with the regularization path's truncated
  optimizer which is appropriate for models that are warm-started
  from one another, but the inception starts from scratch and was given the same
  budget. The inception is now fitted with an untruncated optimizer, as `PLNnetworkfamily`.

# PLNmodels 1.3.0

## New backends and optimizers

* **New built-in Newton optimizer** (`backend = "builtin"`) for PLN, ZIPLN and PLNnetwork: envelope-theorem Newton steps with strong Wolfe line search, no dependency on NLOPT. Substantially faster and more accurate than nlopt on large datasets with full covariance.

* **PLNPCA's `"builtin"` backend rewritten as a profiled trust-region Newton**: the variational block `(M, S)` is profiled out with a per-observation Newton VE-step, and the loadings `(B, C)` are optimised with a saddle-aware trust-region Newton on the resulting objective (analytic Schur-complement Hessian-vector products, Jacobi-preconditioned Steihaug-CG) — replacing the previous spectral projected-gradient `"builtin"`. It reliably reaches a higher variational bound than `"nlopt"` at comparable-to-better speed; tuning keys `cg_maxit`, `maxit_out`, `ftol_out`, `gtol`, `delta0` in `config_optim` (see `?PLNPCA_param`). `"nlopt"` remains the default for `PLNPCA`.

* **Backend defaults revisited package-wide**, based on extensive benchmarking: PLN and PLNPCA keep `"nlopt"` (PLN now consistently faster thanks to `profiled = TRUE`, see below); PLNnetwork and ZIPLNnetwork now default to `"builtin"`, which finds a better optimum at a modest speed cost; ZIPLN keeps its `"builtin"` default. All backends remain configurable via the `backend` argument; see the corresponding `*_param()` documentation for the trade-offs. The `torch` backend is now clearly marked **experimental** everywhere.

* **Quality and speed improvements**: `config_optim$profiled = TRUE` is now the default for full-covariance nlopt fits (faster, slightly better loglik); ZIPLN's variational step now optimises `(M, ψ, R)` jointly via Newton instead of sequentially; PLNPCA shares a single SVD initialisation across ranks and can warm-start from a pre-fitted `PLNfit` for large ranks (`inception`/`init_method`, see `?PLNPCA_param`); ZIPLN's starting point no longer relies on `pscl::zeroinfl` (now an internal LM + binomial GLM routine), which is both much faster and a better starting point — `pscl` is no longer a dependency.

* **Fixed a critical nlopt convergence bug** affecting PLN/PLNPCA: ill-conditioned covariate scaling could trigger the XTOL stopping criterion after very few iterations, well before convergence. The built-in backend was never affected; nlopt is now also fixed via better parameter scaling.

## Internal refactoring (C++)

* **Shared covariance abstraction**: the optimization machinery for PLN's covariance structures (full, diagonal, spherical, fixed) is now expressed once via a small set of C++ traits (`CovTraitsBase` in `covariance_pln.h`) instead of being duplicated per structure. PLNPCA and ZIPLN's variational step now reuse the same machinery instead of separate hand-rolled implementations, removing a substantial amount of duplicated code and fixing minor inefficiencies along the way (e.g. ZIPLN's VE-step used to treat the precision matrix as dense even for diagonal/spherical covariance, at unnecessary `O(np^2)` cost).
* **Consistent C++ naming**: exported optimizer functions across PLN, PLNPCA and ZIPLN now follow the same `{backend}_optimize_{structure}` convention.

## Other changes

* **Uniform covariate normalization**: a `normalize_covariates()` helper (zero mean, unit variance per column) is now applied consistently in all `optimize()` methods (PLN, PLNPCA, PLNnetwork, ZIPLN). This makes the nlopt XTOL criterion scale-invariant and stabilises the torch backend.

* **Parallelism backend**: `future.apply::future_lapply` is replaced by `parallel::mclapply` throughout (stability selection for PLNnetwork / ZIPLNnetwork). Use `options(mc.cores = N)` to set the number of cores.

* **Bug fixes**: PLNnetwork/ZIPLNnetwork's inception (warm-start) model didn't inherit `ftol_em`/`maxit_em` from the user's `config_optim`, silently falling back to defaults and producing a wrong penalty grid; the PLNPCA rank-model objective used `A − Y` where it should use `A − Y ⊙ Z`; various ZIPLN prediction/initialization fixes (#146, #149, #150, #152).

* microcosm data now included (#153, #154); AIC added for PLN and ZIPLN classes (#151); other fixes (#155).

* **CRAN check fixes** (no user-visible effect): removed an unused `-fopenmp` compilation flag in `src/Makevars` that was inadvertently turning on Armadillo's internal OpenMP parallelisation and inflating the CPU/elapsed time ratio of several examples; fixed a GCC `-Wmismatched-new-delete` false positive in `src/packing.cpp`'s internal test helper by restructuring the code (no diagnostic-suppressing pragma involved).

# PLNmodels 1.2.2 (2025-03-21)

* fix for #143 (remove LBFGS_NOCEDAL variant from the possible algorithms)

# PLNmodels 1.2.1 (2025-03-10)

* fix NOTES in CRAN due to missing packages in \link{} (PR #142)
* Now requires R >= 4.1.0 because package code uses the pipe |>  (PR #142)
* fix sandwich variance estimation (PR #140)
* fix use of native pipe to ensure compatibility with R 3.6 (merge PR #125, fix #124)

# PLNmodels 1.2.0 (2024-03-05)

* new feature: ZIPLN (PLN with zero inflation) for standard PLN and PLN Network
  * ZIPLN() and ZIPLNfit-class to allow for zero-inflation in the standard PLN model (merge PR #116)
  * ZIPLNnetwork() and ZIPLNfit_sparse-class to allow for zero-inflation in the  PLNnetwork model (merge PR #118)
  * Code factorization between PLNnetwork and ZIPLNnetwork (and associated classes)
* fix inconsistency between fitted and predict (merge PR #115)

# PLNmodels 1.1.0 (2024-01-08)

* Update documentation of PLN*_param() functions to include torch optimization parameters
* Add (somehow) explicit error message when torch convergence fails
* Change initialization in `variance_jackknife()` and `variance_bootstrap()` to prevent estimation recycling, results from those functions are now comparable to doing jackknife / bootstrap "by hand".
* Merge PR #110 from Cole Trapnell to add:
  - bootstrap estimation of the variance of model parameter
  - improved interface for model initialization / optimisation parameters, which
    are now passed on to jackknife / bootstrap post-treatments
  - better support of GPU when using torch backend
* Change behavior of `predict()` function for PLNfit model to (i) return fitted values if newdata is missing or (ii) perform one VE step to improve fit if responses are provided (fix issue #114)

# PLNmodels 1.0.4 (2023-08-24)

* changed initial value in optim for variational variance (1 -> 0.1) in VE-step of PLN and PLNPCA
* fix sign in objective of VE_step for PLN with full covariance Issue #100
* add a `scale` argument compute_offset() to force the offsets (RLE, CSS, GMPR, Wrench) to be on the same scale as the counts, like TSS.
* add a new "TMM" for compute_offset()
* fix nb_param for PLNLDA, which caused wrong BIC/ICL and erratic model selection
* fix minor issues #102, #103 plus some others
* fix package file documentation as suggested in <https://github.com/r-lib/roxygen2/issues/1491>

# PLNmodels 1.0.3 (2023-07-06)

* higher tolerance on a single test (among 700) that fails on the 'noLD'
  additional architecture on CRAN (tests without long double)

# PLNmodels 1.0.2 (2023-06-21)

* changed initial value in optim for variational variance (1 -> 0.1),
    which caused failure in some cases
* fix bug when using inception in PLNnetwork()
* starting handling of missing data
* slightly faster (factorized) initialization for PCA

# PLNmodels 1.0.1 (2023-02-12)

* fix in the use of future_lapply which used to make post-Treatments in PLNPCA last for ever with multicore in v1.0.0...
* prevent use of bootstrap/jackknife when not appropriate
* fix bug in PLNmixture() when the sequence of cluster numbers (`clusters`) is not of the form `1:K_max`
* use bibentry to replace citEntry in CITATION

# PLNmodels 1.0.0

## Breaking changes

* interface for controlling the fits now use list generated by dedicated functions
   - PLN_param() for PLN
   - PLNLDA_param() for PLNLDA
   - PLNnetwork_param() for PLNnetwork
   - PLNPCA_param() for PLNPCA
   - PLNmixture_param() for PLNmixture
The use of 'control = list()' is deprecated: the code stop and send an error.

* The regression coefficients are now denoted by B, not Theta, such as B = t(Theta).
  We keep on sending back Theta as a field of myPLN$model_par$Theta, but this will soon be deprecated

## New features

* added Barents fish data set
* support for PLN when (inverse) covariance is known/fixed
* estimator of the variance of the model parameters
    * integration of sandwich estimator of the variance-covariance of Theta when Sigma is fixed
    * variational estimation of the variance-covariance based on variational approximation of the Fisher information
    * jackknife estimation of the variance of Theta and Sigma
    * bootstrap estimation of the variance of Theta and Sigma
* handle list of penalty weights in PLNnetwork
* first support for torch optimizers (for PLN and PLNLDA)

## Bug fixes

* fix in objective functions of ve_step of standard PLN models
* fix in objective functions of main  of standard PLN models

# PLNmodels 0.11.7

* fix expression of ELBO in VEstep, related to #91
* typos and regeneration of documentation( HTML5)
* added an S3 method predict_cond to perform conditional predictions
* fix #89 bug by forcing an intercept in `PLNLDA()` and changing `extract_model()` to conform with `model.frame()`

# PLNmodels 0.11.6

* fix wrong use of all.equal
* fix linking problem in new version of nloptr (>=2.0.0)

# PLNmodels 0.11.5

* fixing #79 by using the same variational distribution to approximate
    the spherical case as in the fully parametrized and diagonal cases
* faster examples and build for vignettes
* additional R6 method `$VEStep()` for PLN-PCA, dealing with low rank matrices
* additional R6 method `$project()` for PLN-PCA, used to project newdata into PCA space
* use future_lapply in PLNmixture_family
* remove a NOTE due to a DESeq2 link and a failure on solaris on CRAN machines
* some bug fixes

# PLNmodels 0.11.4

* use future_lapply in PLNPCA, PLNmixture and stability_selection (plan must be set by the user)
* bug fix in prediction for PLN-LDA
* bug fix in gradients of PLN-network and PLN-spherical
* suppressing method `$latent_pos()` which is equivalent to active binding `$latent`
* finalizing integration of PLNmixture (in particular faster smoothing)
* added an argument 'reverse' to the plot methods for criteria, so that users can get their "usual" BIC definition (-2 loglik)

# PLNmodels 0.11.3

* support for covariates in PLNmixture (spherical, diagonal, full)
* more support for PLNmixture (S3/R6 methods, vignette)

# PLNmodels 0.11.2

* Rewriting C++ by merging modern_cpp to dev, thanks to François Gindraud
* various bug fixes in offset
* less verbose about R squared when questionable
* correction in BIC/ICL for PLNPCA
* Enhanced vignettes for PLNPCA and PLNmixture

# PLNmodels 0.11.1

* Add compatibility with factoextra for PLNPCA

# PLNmodels 0.11.0

* Add development version of PLNmixture

# PLNmodels 0.10.7

* add type = "poscounts" option to RLE normalization
* added wrench normalization to the list of available offsets
* added the oaks data set from Jakuschkin et al (2016)

# PLNmodels 0.10.6

* Correction in likelihood of diagonal PLN
* amending test-pln to fulfill CRAN request (error on ATLAS variant of BLAS...)

# PLNmodels 0.10.5

* Refactor code of R6 classes to benefit from Roxygen 7.0.0 R6-related new features for documentation

# PLNmodels 0.10.4

* Change name of variational variance parameters to S2 (used to be S)
* use spell_check to check spelling, found many typos

# PLNmodels 0.10.3

* Change in optimization for all PLN models (PLNs, PCA, LDA, networks): solving in S such that
S = S² for the variational parameters, thus avoiding lower bound and constrained optimization.
Slightly finer results/estimations for similar computational cost, but easier to maintain.

# PLNmodels 0.10.2

* Fix bug in predict() methods when factor levels differ between train and test datasets.
* Fix bug in PLNPCAfit S3 plot() method
* Some simplification in C++ code
* correction/changes in PLN likelihoods? + added constant terms in all likelihoods of all PLN models
* VEstep now available for all model of covariance in PLN (full, diagonal, spherical)

# PLNmodels 0.9.5 - minor release

* removed any use of rmarkdown::paged_table() in the vignettes
* added screenshot.force = FALSE, in knitr options in the vignettes

# PLNmodels 0.9.4 - minor release

* removing dependencies to bioconductor packages, too cumbersome to maintain on CRAN

# PLNmodels 0.9.3 - minor release

* correction in test to comply new class of matrix object

# PLNmodels 0.9.2.9002 - development version

* added the possibility for matrix of weights for the penalty in PLNnetworks

# PLNmodels 0.9.2

* various bug fixes

# PLNmodels 0.9.1

* Use nloptr to prepare CRAN release

# PLNmodels 0.8.2

* Enhancement in PLNLDA

# PLNmodels 0.8.1

* Preparing first CRAN release

# PLNmodels 0.7.0.1

* Added a `NEWS.md` file to track changes to the package.
