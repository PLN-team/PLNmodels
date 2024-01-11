
Submitting PLNmodels version 1.1.0 to CRAN

Fix various bugs, update package documentation for controlling optimization, better integration of the torch backend, better jackknife estimation for PLN.

## Tested environments

* tested locally on Ubuntu Linux 22.04.3 LTS, R-release, GCC

* tested remotely with win-builder 
  - Windows Server 2022, R-devel, 64 bit
  - Windows Server 2022, R-release (4.3.2), 64 bit
  - Windows Server 2022, R-old (4.2.3), 64 bit

* tested remotely with github-action
  - Linux Ubuntu 22.04, R-release
  - Linux Ubuntu 22.04, R-oldrel 
  - Linux Ubuntu 22.04, R-devel 
  - Windows Server 2022, R-release, 64 bit
  - macOS Big Sur 11, R-release 

all status OK except for

* the usual NOTE about libs size (RcppArmadillo)

── R CMD check results ─────────────────────────── PLNmodels 1.1.0 ────
Duration: 4m 32.7s

❯ checking installed package size ... NOTE
    installed size is 11.2Mb
    sub-directories of 1Mb or more:
      doc    2.0Mb
      libs   7.8Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖

R CMD check succeeded
