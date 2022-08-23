
Submitting PLNmodels version 0.11-7 to CRAN. Motivated by a request from CRAN tea

* fix expression of ELBO in VEstep, related to #91
* typos and regeneration of documentation( HTML5) 
* added an S3 method predict_cond to perform conditional predictions
* fix #89 bug by forcing an intercept in `PLNLDA()` and changing `extract_model()` to conform with `model.frame()`

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
- Windows, R-release (winbuilder)
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
