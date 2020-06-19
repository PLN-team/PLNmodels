
## Minor release to answer CRAN request

* amending test-pln to fullfil a CRAN maintainer request 
  (error on an ATLAS variant of BLAS in some anecdotical tests...)
* Correction in likelihood of diagonal PLN

## Tested environments

- macOS 10.13.6 High Sierra, R-release, CRAN's setup (R hub)
- macOS Catalina 10.15, R-release (github action)
- macOS Catalina 10.15, R-devel (github action)
- Linux ubuntu 16.04, R-release (github-action)
- Linux ubuntu 18.04 R-release, (local)
- Linux Debian GCC  R-devel, (R hub)
- Windows Server 2019, R-release (github action)
- Windows Server 2008 R2 SP1, R-release  (R hub)
- Windows, R-release (winbuilder)
- Windows, R-release  (winbuilder)

all status OK except for one note:

* checking installed package size ... NOTE
  installed size is  5.1Mb
  sub-directories of 1Mb or more:
    doc    1.6Mb
    libs   2.6Mb

## R CMD check results

On my computer I get (Ubuntu 18.04, Intel Xeon 3.5 GHz; 64 Go mem)

── R CMD check results ─────────────────────────────────── PLNmodels 0.10.5 ────
Duration: 1m 46.6s

> checking installed package size ... NOTE
    installed size is 10.7Mb
    sub-directories of 1Mb or more:
      doc    1.6Mb
      libs   8.1Mb

0 errors ✓ | 0 warnings ✓ | 1 note x

R CMD check succeeded
