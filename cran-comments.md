Submitting PLNmodels version 1.3.2 to CRAN

This release replaces the dependency on glassoFast by an internal graphical
Lasso solver. There is no change to the package's existing API; one function
is added.

## Main change

`PLNnetwork()` and `ZIPLNnetwork()` called `glassoFast::glassoFast()` at every
step of their alternating optimization. glassoFast's Fortran routine can loop
forever on a nearly collapsed covariance matrix: its inner coordinate descent
has no bound, and since it never returns to R, no R-level timeout can stop it.
This was observed in practice in large simulation campaigns.

The package now ships its own C++ port of the same algorithm
(`src/graphical_lasso.h`, shared with the normalblockr package), which:

* always terminates (non-finite input and zero-variance coordinates are
  rejected up front, the inner loop is bounded) and reports non-convergence;
* checks for user interrupts and R time limits regularly;
* gives the same results as glassoFast on ordinary input, up to machine
  precision.

It is exported as `graphical_lasso()`. glassoFast moves from `Imports` to
`Suggests`, where it is only used by tests checking the equivalence of the two
solvers (skipped when it is not installed). `Rcpp (>= 1.0.10)` is now required,
for its unwind-protect mechanism.

## Test environments

* local Linux (Ubuntu 24.04, R 4.6.1): `R CMD check --as-cran`, no ERROR, no
  WARNING. Its two NOTEs are specific to the local machine (a non-portable
  compilation flag from Ubuntu's R toolchain defaults; HTML Tidy not installed).
* win-builder (R-devel, R-release, R-oldrel): pending.

## R CMD check results

Possibly one NOTE, unchanged from previous submissions:

* installed size (~34Mb, `libs`: RcppArmadillo, nlopt, torch).
