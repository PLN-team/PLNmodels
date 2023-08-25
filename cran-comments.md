
Submitting PLNmodels version 1.0.4 to CRAN

Fix various bugs, better offsets handling and update package documentation as ask by CRAN maintainers.

## Tested environments

* tested locally on Ubuntu Linux 22.04.1 LTS, R-release, GCC

* tested remotely with win-builder 
  - Windows Server 2022, R-devel, 64 bit
  - Windows Server 2022, R-release, 64 bit
  - Windows Server 2022, R-old, 64 bit

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
