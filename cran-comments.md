
Submitting PLNmodels version 0.11-4 to CRAN. Last version available on CRAN is 0.10-6.

The following major changes have been made since 0.10-6

* use future to parallelize code in  PLNmixture and stability_selection (plan must be set by the user)
* various bug fixes
* integration of a new series of fitting procedure via PLNmixture
* Rewriting C++ by merging modern_cpp to dev
* Enhanced vignettes for PLNPCA and PLNmixture
* Add compatibility with factoextra for PLNPCA and PLNLDA
* couple of additional dependencies

## Tested environments

- Linux ubuntu 20.04 R-release, (local)
- Linux ubuntu 20.04, R-devel (github-action)
- Linux ubuntu 20.04, R-release (github-action)
- macOS 10.15, R-devel (github action)
- Windows Server 2019, R-devel (github action)
- macOS 10.13.6 High Sierra, R-release, CRAN's setup (R hub)
- Windows Server 2008 R2 SP1, R-release  (R hub)
- Oracle Solaris 10, x86, 32 bit, R-release (R hub)
- Windows, R-release (winbuilder)
- Windows, R-devel  (winbuilder)
- Windows, R-oldrelease  (winbuilder)

all status OK except for one note:

* checking installed package size ... NOTE
  installed size is  6.0Mb
  sub-directories of 1Mb or more:
    doc    2.1Mb
    libs   2.7Mb (plateform dependent)

## R CMD check results

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
