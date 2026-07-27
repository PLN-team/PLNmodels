
Submitting PLNmodels version 1.3.0 to CRAN

Resubmission fixing the issues found at
https://cran.r-project.org/web/checks/check_results_PLNmodels.html and in the
win-builder pretest of a previous 1.3.0 submission attempt:

* the `-Wmismatched-new-delete` WARN (g++ 15/16 on Linux) was a GCC false
  positive in `src/packing.cpp`'s internal test helper: an arma::mat/vec view
  borrowing memory from a `std::vector` was flagged as freeing it, because GCC
  could inline both the allocation and the (unrelated) destructor into the
  same function. A first attempt silenced this with
  `#pragma GCC diagnostic ignored`, but that pragma itself triggered the
  "checking pragmas in C/C++ headers and code" NOTE. The current fix instead
  isolates the allocation in a separate `noinline` function so GCC can no
  longer correlate the two; no pragma is used anymore.
* an unused `-fopenmp` flag in `src/Makevars` was removed: the package uses no
  OpenMP directive of its own, but the flag was inadvertently turning on
  Armadillo's internal OpenMP parallelisation, which was the likely cause of
  the "CPU time > 2x elapsed time" NOTEs on several examples on Linux.

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

* the `-Wmismatched-new-delete` fix and the `src/Makevars` change were each
  verified by compiling `src/packing.cpp` in isolation with the flags used by
  CRAN's checks (`-Wall -pedantic -O2`): no warning, no pragma.
* win-builder (R-devel, R-release, R-oldrel) now confirms the
  previously reported WARN and NOTEs are gone
* Also check via github action against R-devel, R-release and R-oldrel on Linux 
  (Ubuntu 22.04) and macOS (R 4.3.1, clang 16.0.0) with the same results.

* possibly a NOTE about installed size (34.3Mb, libs: RcppArmadillo, nlopt,
  torch), consistent with previous submissions.
