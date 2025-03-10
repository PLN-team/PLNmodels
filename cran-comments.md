
Submitting PLNmodels version 1.2.1 to CRAN (CRAN request)

Most of the problems on https://cran.r-project.org/web/checks/check_results_PLNmodels.html are due to external (yet important) packages (nloptr essentially).

PLNmodels is scheduled for archiving today, March 10. Although noptr is due for submission soon, the corrections have not yet been made on CRAN: while archiving is not required for nloptr, it is required for the pckages that depend on it (including PLNmodels).

We corrected the couple of NOTES imputed to PLNmodels, hoping to avoid the archival of our package.

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

* the usual NOTE about libs size (RcppArmadillo)
* a note about the number of dependencies

── R CMD check results ────────────────────────────── PLNmodels 1.2.1 ────
Duration: 6m 30.2s

❯ checking package dependencies ... NOTE
  Imports includes 21 non-default packages.
  Importing from so many packages makes the package vulnerable to any of
  them becoming unavailable.  Move as many as possible to Suggests and
  use conditionally.

❯ checking installed package size ... NOTE
    installed size is 15.2Mb
    sub-directories of 1Mb or more:
      data   1.4Mb
      doc    2.3Mb
      libs   9.8Mb

0 errors ✔ | 0 warnings ✔ | 2 notes ✖

R CMD check succeeded
