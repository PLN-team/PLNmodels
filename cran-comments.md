
## Minor change to comply with check in R-dev

* added LazyData: true to the DESCRIPTION file to avoid warnings on CRAN dev version
* corrected URI in vignette PLNPCA
* changed some unitary tests which were failing in for developement version of R, due to the 
change of class in matrix

> class(matrix(1 : 4, 2, 2))
  [1] "matrix" "array"

We correct this in this new minor release of PLNmodels


## Tested environments

- local ubuntu 18.04 install, R 3.6.1
- ubuntu 14.04 (on travis-ci), R 3.6.1
- ubuntu 14.04 (on travis-ci), R devel
- Mac OS X (on travis-ci), R 3.6.1
- windows (on appveyor), R 3.6.1
- win-builder, R 3.6.1

## R CMD check results

On my computer I get (Ubuntu 18.04, Intel Xeon 3.5 GHz; 32 Go mem)

── R CMD check results ──────────────────────────────────── PLNmodels 0.9.3 ────
Duration: 3m 38.7s

❯ checking installed package size ... NOTE
    installed size is 10.8Mb
    sub-directories of 1Mb or more:
      doc    1.6Mb
      libs   8.3Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖

R CMD check succeeded


