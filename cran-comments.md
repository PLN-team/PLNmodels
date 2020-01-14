
## Minor changes to answer CRAN maintainer request

* small corrections in C++ to avoid warnings
* remove direct dependencies to bioConductor packages (which was resulting 
  in a compilation error for R devel Windows-64-gcc8)

We correct this in this new minor release of PLNmodels

## Tested environments

- local ubuntu 18.04 install, R 3.6.2
- ubuntu 14.04 (on travis-ci), R 3.6.2
- ubuntu 14.04 (on travis-ci), R devel
- Debian Linux, R-devel, GCC (on R hub), R 3.6.2
- Debian Linux, R-devel, GCC (on R hub), R devel
- Mac OS X (on travis-ci), R 3.6.2
- windows (on appveyor), R 3.6.2
- win-builder, R 3.6.2
- win-builder, R devel

## R CMD check results

On my computer I get (Ubuntu 18.04, Intel Xeon 3.5 GHz; 32 Go mem)

── R CMD check results ──────────────────────────────────── PLNmodels 0.9.4 ────
Duration: 3m 51s

> checking installed package size ... NOTE
    installed size is 11.0Mb
    sub-directories of 1Mb or more:
      doc    1.6Mb
      libs   8.4Mb

0 errors ✓ | 0 warnings ✓ | 1 note x

R CMD check succeeded

