
## Another resubmission

 * commenting a test that randomly fails on CRAN Debian server

## Resubmission of a new submission

To reduce time during the R check, a new version of PLNmodels was (re)submitted where I 

* use smaller data in test files to reduce the testing time
* add some dont run in examples to save time during check
* simplify in the vignettes to save time during check

On my computer I get (Ubuntu 18.04, Intel Xeon 3.5 GHz; 32 Go mem)
✔  checking whether package ‘PLNmodels’ can be installed (48.8s)
N  checking installed package size ...
     installed size is 10.8Mb
     sub-directories of 1Mb or more:
       doc    1.6Mb
       libs   8.3Mb
✔  checking examples (20.1s)
✔  Running ‘testthat.R’ [27s/27s]
✔  checking re-building of vignette outputs (44.3s)
── R CMD check results ──────────────────────────────────── PLNmodels 0.9.2 ────
Duration: 3m 29.2s

## Tested environments

- local ubuntu 18.04 install, R 3.6.0
- ubuntu 14.04 (on travis-ci), R 3.6.0
- ubuntu 14.04 (on travis-ci), R devel
- Mac OS X (on travis-ci), R 3.6.0
- windows (on appveyor), R 3.6.0
- win-builder, R 3.6.0

## R CMD check results

0 errors | 0 warnings | 1 note

* installed size is 10.8Mb
     sub-directories of 1Mb or more:
       doc    1.6Mb
       libs   8.3Mb
