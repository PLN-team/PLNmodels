
Submitting PLNmodels version 0.11-5 to CRAN. Motivated by a request from CRAN team to 
correct a test that randomly fails on some distributions.

The following major changes have been made since 0.11-5

* faster examples and build for vignettes
* use the same variational distribution to approximate the spherical case as
  in the fully parametrized and diagonal cases
* additional R6 method `$VEStep()` for PLN-PCA, dealing with low rank matrices
* additional R6 method `$project()` for PLN-PCA, used to project newdata into PCA space
* use future_lapply in PLNmixture_family
* remove a NOTE due to a DESeq2 link and a failure on solaris on CRAN machines
* various bug fixes

## Tested environments

- Linux ubuntu 20.04 R-release, (local)
- Linux ubuntu 20.04, R-devel (github-action)
- Linux ubuntu 20.04, R-release (github-action)
- macOS 10.15, R-devel (github action)
- Windows Server 2019, R-devel (github action)
- macOS 10.13.6 High Sierra, R-release, CRAN's setup (R hub)
- Windows Server 2008 R2 SP1, R-release  (R hub)
- Windows, R-release (winbuilder)
- Windows, R-devel  (winbuilder)
- Windows, R-oldrelease  (winbuilder)

all status OK except for

* the usual NOTE about libs size (RccpArmadillo)

  checking installed package size ... NOTE
  installed size is  6.0Mb
  sub-directories of 1Mb or more:
    doc    2.1Mb
    libs   2.7Mb (platform dependent)

On my computer I get (Ubuntu 20.04, Intel Xeon 3.5 GHz; 64 Go mem)

── R CMD check results ─────────────────────────────────── PLNmodels 0.11.4 ────
Duration: 2m 41.1s

> checking installed package size ... NOTE
    installed size is 18.1Mb
    sub-directories of 1Mb or more:
      doc    2.1Mb
      libs  14.9Mb

0 errors ✓ | 0 warnings ✓ | 1 note x

R CMD check succeeded
