
## Third Submission

* added some dont run in examples to save time during check
* simplification in the vignettes to save time during check

On my compouter I get (Ubuntu 18.04, Intel Xeon 3.5 GHz; 32 Go mem)
✔  checking whether package ‘PLNmodels’ can be installed (48.8s)
N  checking installed package size ...
     installed size is 10.8Mb
     sub-directories of 1Mb or more:
       doc    1.6Mb
       libs   8.3Mb
✔  checking examples (20.1s)
✔  Running ‘testthat.R’ [27s/27s]
✔  checking re-building of vignette outputs (44.3s)

## Second Submission

* using smaller data in test files to reduce the testing time

## New Submission

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
