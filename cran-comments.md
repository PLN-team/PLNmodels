
## Intermediate release with various updates

* refactor code to benefit from Roxygen 7.0.0 R6-related new features for documentation
* used spell_check to check spelling, found many typos
* updated optimization routines
* fix bugs, some simplification in C++ code
* add a couple of features

## Tested environments

- macOS 10.13.6 High Sierra, R-release, CRAN's setup (R hub)
- macOS Catalina 10.15, R-release (github action)
- macOS Catalina 10.15, R-devel (github action)
- Linux ubuntu 16.04, R-release (github-action)
- Linux ubuntu 18.04 R-release, (local)
- Windows Server 2019, R-release (github action)
- Windows Server 2008 R2 SP1, R-release  (on R hub)
- Windows, R-release (winbuilder)
- Windows, R-release  (winbuilder)

all status OK except for one note:

* checking installed package size ... NOTE
  installed size is  5.1Mb
  sub-directories of 1Mb or more:
    doc    1.6Mb
    libs   2.6Mb

## R CMD check results

On my computer I get (Ubuntu 18.04, Intel Xeon 3.5 GHz; 32 Go mem)

── R CMD check results ─────────────────────────────────── PLNmodels 0.10.5 ────
Duration: 2m 30.3s

> checking installed package size ... NOTE
    installed size is 10.7Mb
    sub-directories of 1Mb or more:
      doc    1.6Mb
      libs   8.1Mb

0 errors ✓ | 0 warnings ✓ | 1 note x

R CMD check succeeded
