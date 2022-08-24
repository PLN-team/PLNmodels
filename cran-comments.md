
Submitting PLNmodels version 0.11-7 to CRAN, motivated by a request from CRAN team (HTML5)

* fix expression of ELBO in VEstep, related to #91
* typos and regeneration of documentation( HTML5) 
* added an S3 method predict_cond to perform conditional predictions
* fix #89 bug by forcing an intercept in `PLNLDA()` and changing `extract_model()` to conform with `model.frame()`

## Tested environments

* tested locally on Ubuntu Linux 20.04.1 LTS, R-release, GCC

* tested remotely with win-builder 
  - Windows Server 2022, R-devel, 64 bit

* tested remotely with R-hub 
  - Windows Server 2022, R-devel, 64 bit
  - Fedora Linux, R-devel, clang, gfortran
	- Ubuntu Linux 20.04.1 LTS, R-release, GCC

* tested remotely with github-action
  - Linux ubuntu 20.04, R-release
  - Linux ubuntu 20.04, R-oldrel 
  - Linux ubuntu 20.04, R-devel 
  - Windows Server 2022, R-release, 64 bit
  - macOS Big Sur 11, R-release 

all status OK except for

* the usual NOTE about libs size (RccpArmadillo)

  checking installed package size ... NOTE
  installed size is  6.0Mb
  sub-directories of 1Mb or more:
    doc    2.1Mb
    libs   2.7Mb (platform dependent)

On my computer I get (Ubuntu 22.04)

── R CMD check results ──────────── PLNmodels 0.11.7 ────
Duration: 2m 25.7s

> checking installed package size ... NOTE
    installed size is 14.4Mb
    sub-directories of 1Mb or more:
      doc    2.1Mb
      libs  11.2Mb

0 errors ✓ | 0 warnings ✓ | 1 note x

R CMD check succeeded
