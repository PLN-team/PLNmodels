
Submitting PLNmodels version 1.2.2 to CRAN (CRAN request)

Hopefully corrected usages of different version of nlopt (2.9.x) in nloptr by 
removing reference to the LBFGS_NOCEDAL algorithm

## Tested environments

* tested locally on Ubuntu Linux 24.04.2 LTS, R-release, GCC

* tested remotely with win-builder 
  - Windows Server 2022, R 4.3.3, 64 bit
  - Windows Server 2022, R 4.4.3, 64 bit
  - Windows Server 2022, R-devel, 64 bit

* tested remotely with github-action
  - Linux Ubuntu 24.04, R-release
  - Linux Ubuntu 24.04, R-devel 
  - Windows Server 2022, R-release, 64 bit
  - Windows Server 2022, R-oldrel , 64 bit
  - macOS Sonoma 14, R-release 
  - macOS Sonoma 14, R-oldrel 

all status OK except for 2 NOTES

* the usual NOTE about libs size (RcppArmadillo and nlopt)
* a note about the number of dependencies

── R CMD check results ──────────────────────────── PLNmodels 1.2.2 ────
Duration: 3m 13.4s

❯ checking package dependencies ... NOTE
  Imports includes 21 non-default packages.
  Importing from so many packages makes the package vulnerable to any of
  them becoming unavailable.  Move as many as possible to Suggests and
  use conditionally.

❯ checking installed package size ... NOTE
    installed size is 24.8Mb
    sub-directories of 1Mb or more:
      data   1.4Mb
      doc    2.4Mb
      libs  19.3Mb

0 errors ✔ | 0 warnings ✔ | 2 notes ✖

R CMD check succeeded
