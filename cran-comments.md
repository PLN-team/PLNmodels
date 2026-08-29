Submitting PLNmodels version 1.3.1 to CRAN

This is a patch release containing a single bug fix. There is no change to the
package's API, to its dependencies, or to its compiled code.

## Bug fixed

`ZIPLNnetwork()` fitted its inception (warm-start) model with the optimizer
configuration meant for the regularization path. That configuration is
deliberately truncated, which is appropriate for the models along the path, 
since each is warm-started from the previous one, but not for the inception, which starts from scratch.

The inception is now fitted with an untruncated optimizer, mirroring the
`cfg_inception` mechanism `PLNnetworkfamily` already uses for the same reason.
The path itself keeps its truncated settings, so the added cost is a fraction of
a second per collection.

The change is confined to `ZIPLNnetworkfamily$initialize()` in
`R/PLNnetworkfamily-class.R`. Results for `PLN`, `PLNPCA`, `PLNnetwork` and
`ZIPLN` are unchanged. Users who worked around the problem by passing a fitted
`ZIPLN` object as `inception` obtain the same results as before, now
automatically.

## Test environments

* local Linux (Ubuntu, R 4.6.1): `R CMD check --as-cran`, no ERROR, no WARNING.
* GitHub Actions against R-devel, R-release and R-oldrel on Linux (Ubuntu 22.04)
  and macOS.
* win-builder (R-devel, R-release, R-oldrel).

The package's own test suite passes with no failure, error or warning; the files
covering the modified code — `test-zipln.R` (29 tests), `test-ziplnfit.R` (63),
`test-ziplnnetworkfamily.R` (55) and `test-plnnetworkfamily.R` (55) — were run
individually as well.

## R CMD check results

There is one NOTE, unchanged from previous submissions:

* installed size (~34Mb, `libs`: RcppArmadillo, nlopt, torch).
