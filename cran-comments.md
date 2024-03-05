
Submitting PLNmodels version 1.2.0 to CRAN

new feature: ZIPLN (PLN with zero inflation) for standard PLN and PLN Network

## Tested environments

* tested locally on Ubuntu Linux 22.04.3 LTS, R-release, GCC

* tested remotely with win-builder 
  - Windows Server 2022, R-devel, 64 bit
  - Windows Server 2022, R-release, 64 bit
  - Windows Server 2022, R-old, 64 bit

* tested remotely with github-action
  - Linux Ubuntu 22.04, R-release
  - Linux Ubuntu 22.04, R-devel 
  - Windows Server 2022, R-release, 64 bit
  - Windows Server 2022, R-oldrel , 64 bit
  - macOS Big Sur 11, R-release 
  - macOS Big Sur 11, R-oldrel 

all status OK except for

* the usual NOTE about libs size (RcppArmadillo)

 R CMD check results  PLNmod
Duration: 6m 27.6s

❯ checking installed package size ... NOTE
    installed size is 23.4Mb
    sub-directories of 1Mb or more:
      data   1.4Mb
      doc    2.1Mb
      libs  18.3Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖

R CMD check succeeded

── R CMD check results ─────────────────────────── PLNmodels 1.2.0 ────


