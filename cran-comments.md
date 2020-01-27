
## Minor changes to answer CRAN maintainer request

Hoping for WARN suppression on solaris and os build, we

* removed any use of rmarkdown::paged_table() in the vignettes
* added screenshot.force = FALSE, in knitr options in the vignettes

## Tested environments

- local ubuntu 18.04 install, R 3.6.2
- ubuntu 16.04 (on travis-ci), R 3.6.2
- ubuntu 16.04 (on travis-ci), R devel
- Mac OS X (on travis-ci), R 3.6.2
- Mac OS X (on travis-ci), R devel
- Debian Linux, R-devel, GCC (on R hub), R 3.6.2
- Debian Linux, R-devel, GCC (on R hub), R devel
- windows (on appveyor), R 3.6.2
- win-builder, R 3.6.2
- win-builder, R devel

## R CMD check results

On my computer I get (Ubuntu 18.04, Intel Xeon 3.5 GHz; 32 Go mem)

── R CMD check results ──────────────────────────────────── PLNmodels 0.9.5 ────
Duration: 3m 44.1s

> checking installed package size ... NOTE
    installed size is 10.9Mb
    sub-directories of 1Mb or more:
      doc    1.6Mb
      libs   8.4Mb

0 errors ✓ | 0 warnings ✓ | 1 note x

R CMD check succeeded
