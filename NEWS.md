# PLNmodels 1.2.2-9100

* various fix in ZIPLN model (prediction and initialization #146, #149, #150, #152)
* microcosm data now included (#153, #154)
* add AIC for PLN and ZIPLN classes (#151)
* other fixes (#155)

# PLNmodels 1.2.2 (2025-03-21)

* fix for #143 (remove LBFGS_NOCEDAL variant from the possible algorithms)

# PLNmodels 1.2.1 (2025-03-10)

* fix NOTES in CRAN due to missing packages in \link{} (PR #142)
* Now requires R >= 4.1.0 because package code uses the pipe |>  (PR #142)
* fix sandwich variance estimation (PR #140)
* fix use of native pipe to ensure compatibility with R 3.6 (merge PR #125, fix #124)

# PLNmodels 1.2.0 (2024-03-05)

* new feature: ZIPLN (PLN with zero inflation) for standard PLN and PLN Network
  * ZIPLN() and ZIPLNfit-class to allow for zero-inflation in the standard PLN model (merge PR #116)
  * ZIPLNnetwork() and ZIPLNfit_sparse-class to allow for zero-inflation in the  PLNnetwork model (merge PR #118)
  * Code factorization between PLNnetwork and ZIPLNnetwork (and associated classes)
* fix inconsistency between fitted and predict (merge PR #115)

# PLNmodels 1.1.0 (2024-01-08)

* Update documentation of PLN*_param() functions to include torch optimization parameters
* Add (somehow) explicit error message when torch convergence fails
* Change initialization in `variance_jackknife()` and `variance_bootstrap()` to prevent estimation recycling, results from those functions are now comparable to doing jackknife / bootstrap "by hand". 
* Merge PR #110 from Cole Trapnell to add:
  - bootstrap estimation of the variance of model parameter
  - improved interface for model initialization / optimisation parameters, which
    are now passed on to jackknife / bootstrap post-treatments
  - better support of GPU when using torch backend
* Change behavior of `predict()` function for PLNfit model to (i) return fitted values if newdata is missing or (ii) perform one VE step to improve fit if responses are provided (fix issue #114)

# PLNmodels 1.0.4 (2023-08-24)

* changed initial value in optim for variational variance (1 -> 0.1) in VE-step of PLN and PLNPCA
* fix sign in objective of VE_step for PLN with full covariance Issue #100
* add a `scale` argument compute_offset() to force the offsets (RLE, CSS, GMPR, Wrench) to be on the same scale as the counts, like TSS.
* add a new "TMM" for compute_offset()
* fix nb_param for PLNLDA, which caused wrong BIC/ICL and erratic model selection
* fix minor issues #102, #103 plus some others
* fix package file documentation as suggested in <https://github.com/r-lib/roxygen2/issues/1491>

# PLNmodels 1.0.3 (2023-07-06)

* higher tolerance on a single test (among 700) that fails on the 'noLD' 
  additional architecture on CRAN (tests without long double)

# PLNmodels 1.0.2 (2023-06-21)

* changed initial value in optim for variational variance (1 -> 0.1),
    which caused failure in some cases
* fix bug when using inception in PLNnetwork()
* starting handling of missing data
* slightly faster (factorized) initialization for PCA

# PLNmodels 1.0.1 (2023-02-12)

* fix in the use of future_lapply which used to make post-Treatments in PLNPCA last for ever with multicore in v1.0.0...
* prevent use of bootstrap/jackknife when not appropriate
* fix bug in PLNmixture() when the sequence of cluster numbers (`clusters`) is not of the form `1:K_max`
* use bibentry to replace citEntry in CITATION

# PLNmodels 1.0.0

## Breaking changes

* interface for controlling the fits now use list generated by dedicated functions
   - PLN_param() for PLN
   - PLNLDA_param() for PLNLDA
   - PLNnetwork_param() for PLNnetwork
   - PLNPCA_param() for PLNPCA
   - PLNmixture_param() for PLNmixture
The use of 'control = list()' is deprecated: the code stop and send an error.

* The regression coefficients are now denoted by B, not Theta, such as B = t(Theta).
  We keep on sending back Theta as a field of myPLN$model_par$Theta, but this will soon be deprecated

## New features

* added Barents fish data set
* support for PLN when (inverse) covariance is known/fixed
* estimator of the variance of the model parameters 
    * integration of sandwich estimator of the variance-covariance of Theta when Sigma is fixed
    * variational estimation of the variance-covariance based on variational approximation of the Fisher information 
    * jackknife estimation of the variance of Theta and Sigma
    * bootstrap estimation of the variance of Theta and Sigma
* handle list of penalty weights in PLNnetwork
* first support for torch optimizers (for PLN and PLNLDA)

## Bug fixes

* fix in objective functions of ve_step of standard PLN models
* fix in objective functions of main  of standard PLN models

# PLNmodels 0.11.7

* fix expression of ELBO in VEstep, related to #91
* typos and regeneration of documentation( HTML5) 
* added an S3 method predict_cond to perform conditional predictions
* fix #89 bug by forcing an intercept in `PLNLDA()` and changing `extract_model()` to conform with `model.frame()`

# PLNmodels 0.11.6

* fix wrong use of all.equal
* fix linking problem in new version of nloptr (>=2.0.0)

# PLNmodels 0.11.5

* fixing #79 by using the same variational distribution to approximate 
    the spherical case as in the fully parametrized and diagonal cases
* faster examples and build for vignettes
* additional R6 method `$VEStep()` for PLN-PCA, dealing with low rank matrices
* additional R6 method `$project()` for PLN-PCA, used to project newdata into PCA space
* use future_lapply in PLNmixture_family
* remove a NOTE due to a DESeq2 link and a failure on solaris on CRAN machines
* some bug fixes

# PLNmodels 0.11.4

* use future_lapply in PLNPCA, PLNmixture and stability_selection (plan must be set by the user)
* bug fix in prediction for PLN-LDA
* bug fix in gradients of PLN-network and PLN-spherical
* suppressing method `$latent_pos()` which is equivalent to active binding `$latent`
* finalizing integration of PLNmixture (in particular faster smoothing)
* added an argument 'reverse' to the plot methods for criteria, so that users can get their "usual" BIC definition (-2 loglik)

# PLNmodels 0.11.3

* support for covariates in PLNmixture (spherical, diagonal, full)
* more support for PLNmixture (S3/R6 methods, vignette)

# PLNmodels 0.11.2

* Rewriting C++ by merging modern_cpp to dev, thanks to François Gindraud
* various bug fixes in offset
* less verbose about R squared when questionable
* correction in BIC/ICL for PLNPCA
* Enhanced vignettes for PLNPCA and PLNmixture

# PLNmodels 0.11.1

* Add compatibility with factoextra for PLNPCA

# PLNmodels 0.11.0

* Add development version of PLNmixture

# PLNmodels 0.10.7

* add type = "poscounts" option to RLE normalization
* added wrench normalization to the list of available offsets
* added the oaks data set from Jakuschkin et al (2016)

# PLNmodels 0.10.6

* Correction in likelihood of diagonal PLN
* amending test-pln to fulfill CRAN request (error on ATLAS variant of BLAS...)

# PLNmodels 0.10.5

* Refactor code of R6 classes to benefit from Roxygen 7.0.0 R6-related new features for documentation

# PLNmodels 0.10.4

* Change name of variational variance parameters to S2 (used to be S)
* use spell_check to check spelling, found many typos

# PLNmodels 0.10.3

* Change in optimization for all PLN models (PLNs, PCA, LDA, networks): solving in S such that 
S = S² for the variational parameters, thus avoiding lower bound and constrained optimization. 
Slightly finer results/estimations for similar computational cost, but easier to maintain.

# PLNmodels 0.10.2

* Fix bug in predict() methods when factor levels differ between train and test datasets. 
* Fix bug in PLNPCAfit S3 plot() method
* Some simplification in C++ code
* correction/changes in PLN likelihoods? + added constant terms in all likelihoods of all PLN models
* VEstep now available for all model of covariance in PLN (full, diagonal, spherical)

# PLNmodels 0.9.5 - minor release

* removed any use of rmarkdown::paged_table() in the vignettes
* added screenshot.force = FALSE, in knitr options in the vignettes

# PLNmodels 0.9.4 - minor release

* removing dependencies to bioconductor packages, too cumbersome to maintain on CRAN

# PLNmodels 0.9.3 - minor release

* correction in test to comply new class of matrix object

# PLNmodels 0.9.2.9002 - development version

* added the possibility for matrix of weights for the penalty in PLNnetworks

# PLNmodels 0.9.2

* various bug fixes

# PLNmodels 0.9.1

* Use nloptr to prepare CRAN release

# PLNmodels 0.8.2

* Enhancement in PLNLDA

# PLNmodels 0.8.1

* Preparing first CRAN release

# PLNmodels 0.7.0.1

* Added a `NEWS.md` file to track changes to the package.
