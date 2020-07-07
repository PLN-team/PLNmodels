# PLNmodels 0.10.6

* Correction in likelihood of diagonal PLN
* amending test-pln to fullfil CRAN request (error on ATLAS variant of BLAS...)

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
