## Resubmission 3

Removing RcppArmadillo from the Imports field (only LinkingTo)

## Resubmission 2

This is a second resubmission. In this version we

- include NLopt via nloptr for portability
- corrected several bugs
- boosted code coverage with more tests 

## Resubmission 1

This is a resubmission. In this version I have

- added references in the description field of DESCRIPTION
- corrected a typo in the description field of DESCRIPTION
- added examples to almost all the documented functions (Rd files)
- limited the use of \dontrun whenever possible (to keep examples runtimes < 5 sec.)
- added more tests via testthat

## Tested environments

- local ubuntu 18.04 install, R 3.5.3
- ubuntu 14.04 (on travis-ci), R 3.5.3
- ubuntu 14.04 (on travis-ci), R devel
- Mac OS X (on travis-ci), R 3.5.3
- windows (on appveyor), R 3.5.3
- win-builder, R 3.5.3

## R CMD check results

0 errors | 0 warnings | 1 note

* Installed size is 10.8Mb
  sub-directories of 1Mb or more:
    - doc    2.2Mb
    - libs   7.8Mb
