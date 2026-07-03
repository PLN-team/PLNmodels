
Submitting PLNmodels version 1.3.0 to CRAN

## Summary of changes since 1.2.2

This is a feature/optimizer release, not just a bugfix release:

* New built-in (dependency-free) Newton-type optimizers for PLN, PLNPCA, PLNnetwork,
  ZIPLN and ZIPLNnetwork, generally faster and/or more accurate than the existing
  NLOPT-based backends, especially on large datasets. PLNPCA's built-in backend is a
  profiled trust-region Newton with an analytic Schur-complement Hessian-vector
  product; nlopt remains its default for robustness.
* Backend defaults revisited package-wide based on extensive benchmarking (see
  NEWS.md for the full per-model rationale).
* A shared C++ traits abstraction (`CovTraitsBase`) removes a substantial amount of
  duplicated optimization code across PLN's covariance structures and is now reused
  by PLNPCA and ZIPLN.
* `future`/`future.apply` replaced by `parallel::mclapply` (avoids BLAS
  oversubscription issues); `pscl` dropped as a dependency (ZIPLN now uses a faster
  internal starting point).
* A handful of bug fixes (nlopt XTOL false-early-convergence on ill-conditioned
  covariates, PLNnetwork/ZIPLNnetwork inception config inheritance, ZIPLN
  prediction/initialisation fixes) — see NEWS.md for the full list.

Full details in NEWS.md.

## Tested environments

* tested locally on Ubuntu Linux 24.04 LTS, R 4.6.1, GCC 13.3.0 (`R CMD check --as-cran`)

* remote checks (rhub / win-builder / github-action) to be run before submission —
  results to be appended here.

## R CMD check results (local)

0 errors ✔ | 0 warnings ✔ | 4 notes ✖

❯ checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    ‘-mno-omit-leaf-frame-pointer’
  This comes from the local Ubuntu-patched R toolchain's default hardening
  CFLAGS/CXXFLAGS (visible in `R CMD config CFLAGS`), not from the package's
  `src/Makevars`. Not expected to reproduce identically on CRAN's build machines.

❯ checking examples ... NOTE
  A number of examples (PLNmixture*, PLNnetwork*) show CPU (user+system) time
  above 5s, but elapsed time stays under 3s in every case. The gap is caused by
  multithreaded BLAS on the local machine, not by actual wall-clock runtime.

installed size is 34.3Mb (libs 27.3Mb: RcppArmadillo, nlopt, torch), consistent
with previous submissions.

* the usual NOTE about libs size (RcppArmadillo, nlopt, torch), if it reappears
  on CRAN's machines
* possibly a note about the number of dependencies (Imports currently includes
  torch, needed for the experimental `backend = "torch"` option)

R CMD check succeeded locally; awaiting remote check results before final submission.
