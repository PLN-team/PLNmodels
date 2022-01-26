
Submitting PLNmodels version 0.11-6 to CRAN. Motivated by a request from CRAN team to 
correct a wrong use of all.equal

+ fixing linkign problem to new version of nloptr

## Tested environments

- Linux ubuntu 20.04 R-release, (local)
- Windows Server 2022, R-devel, 64 bit
- Windows Server 2019, R-release (github action)
- macOS 10.15, R-release (github action)
- Linux ubuntu 20.04, R-release (github-action)
- Linux ubuntu 20.04, R-devel (github-action)
- Linux ubuntu 20.04, R-oldrel (github-action)
- Windows Server 2022, R-devel, 64 bit (R hub)

- Windows, R-devel  (winbuilder)
- Windows, R-release (winbuilder)  (FAILURE, no nloptr)
- Windows, R-oldrelease  (winbuilder)
- macOS 10.13.6 High Sierra, R-release, CRAN's setup (R hub)

all status OK except for

* the usual NOTE about libs size (RccpArmadillo)

  checking installed package size ... NOTE
  installed size is  6.0Mb
  sub-directories of 1Mb or more:
    doc    2.1Mb
    libs   2.7Mb (platform dependent)

On some Windows install,
* checking examples ... NOTE
Examples with CPU (user + system) or elapsed time > 5s
                      user system elapsed
plot.PLNmixturefamily 6.05   0.01    6.07
PLNmixture            5.64   0.01    5.71

* FAILURE DUE TO NLOPTR Not installing correctly

Windows R-release (winbuilder)


On my computer I get (Ubuntu 20.04, Intel Xeon 3.5 GHz; 64 Go mem)

── R CMD check results ──────────── PLNmodels 0.11.6 ────
Duration: 3m 46.7s

> checking installed package size ... NOTE
    installed size is 24.1Mb
    sub-directories of 1Mb or more:
      doc    2.1Mb
      libs  20.8Mb

0 errors ✓ | 0 warnings ✓ | 1 note x

R CMD check succeeded
