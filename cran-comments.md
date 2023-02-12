
Submitting PLNmodels version 1.0.1 to CRAN

* fix in the use of future_lapply which used to make postTreatment in PLNPCA last for ever with multicore in v1.0.0...
* prevent use of bootstrap/jackknife when not appropriate
* fix bug in PLNmixture() when the sequence of cluster numbers (`clusters`) is not of the form `1:K_max`
* use bibentry to replace citEntry in CITATION

## Tested environments

* tested locally on Ubuntu Linux 22.04.1 LTS, R-release, GCC

* tested remotely with win-builder 
  - Windows Server 2022, R-devel, 64 bit
  - Windows Server 2022, R-release, 64 bit
  - Windows Server 2022, R-old, 64 bit

* tested remotely with R-hub 
  - Windows Server 2022, R-devel, 64 bit
	- Ubuntu Linux 22.04 LTS, R-release, GCC

* tested remotely with github-action
  - Linux Ubuntu 22.04, R-release
  - Linux Ubuntu 22.04, R-oldrel 
  - Linux Ubuntu 22.04, R-devel 
  - Windows Server 2022, R-release, 64 bit
  - macOS Big Sur 11, R-release 

all status OK except for

* the usual NOTE about libs size (RcppArmadillo)

❯ checking installed package size ... NOTE
    installed size is 19.4Mb
    sub-directories of 1Mb or more:
      doc    2.2Mb
      libs  15.9Mb
      
0 errors ✓ | 0 warnings ✓ | 1 note x

R CMD check succeeded
