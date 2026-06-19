# Package index

## Top-level fitting functions

Set of functions to fit variants of the Poisson lognormal model

- [`PLNmodels`](https://pln-team.github.io/PLNmodels/reference/PLNmodels-package.md)
  [`PLNmodels-package`](https://pln-team.github.io/PLNmodels/reference/PLNmodels-package.md)
  : PLNmodels: Poisson Lognormal Models
- [`PLN()`](https://pln-team.github.io/PLNmodels/reference/PLN.md) :
  Poisson lognormal model
- [`ZIPLN()`](https://pln-team.github.io/PLNmodels/reference/ZIPLN.md) :
  Zero Inflated Poisson lognormal model
- [`PLNLDA()`](https://pln-team.github.io/PLNmodels/reference/PLNLDA.md)
  : Poisson lognormal model towards Linear Discriminant Analysis
- [`PLNPCA()`](https://pln-team.github.io/PLNmodels/reference/PLNPCA.md)
  : Poisson lognormal model towards Principal Component Analysis
- [`PLNnetwork()`](https://pln-team.github.io/PLNmodels/reference/PLNnetwork.md)
  : Sparse Poisson lognormal model for network inference
- [`ZIPLNnetwork()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNnetwork.md)
  : Zero Inflated Sparse Poisson lognormal model for network inference
- [`PLNmixture()`](https://pln-team.github.io/PLNmodels/reference/PLNmixture.md)
  : Poisson lognormal mixture model

## Poisson lognormal fit

Description of the PLNfit object and methods for its manipulation. Any
PLN variant in the package inherits from this class (PLNPCAfit,
PLNnetworkfit, PLNLDAfit).

- [`PLNfit`](https://pln-team.github.io/PLNmodels/reference/PLNfit.md) :
  An R6 Class to represent a PLNfit in a standard, general framework

- [`PLNfit_diagonal`](https://pln-team.github.io/PLNmodels/reference/PLNfit_diagonal.md)
  [`PLNLDAfit_spherical`](https://pln-team.github.io/PLNmodels/reference/PLNfit_diagonal.md)
  : An R6 Class to represent a PLNfit in a standard, general framework,
  with diagonal residual covariance

- [`PLNfit_fixedcov`](https://pln-team.github.io/PLNmodels/reference/PLNfit_fixedcov.md)
  : An R6 Class to represent a PLNfit in a standard, general framework,
  with fixed (inverse) residual covariance

- [`PLNfit_genpop`](https://pln-team.github.io/PLNmodels/reference/PLNfit_genpop.md)
  : An R6 Class to represent a PLNfit with a residual covariance
  structured by a fixed correlation matrix (e.g. a genetic relationship
  matrix), motivated by population genetics

- [`PLNfit_spherical`](https://pln-team.github.io/PLNmodels/reference/PLNfit_spherical.md)
  : An R6 Class to represent a PLNfit in a standard, general framework,
  with spherical residual covariance

- [`PLN_param()`](https://pln-team.github.io/PLNmodels/reference/PLN_param.md)
  : Control of a PLN fit

- [`coef(`*`<PLNfit>`*`)`](https://pln-team.github.io/PLNmodels/reference/coef.PLNfit.md)
  : Extract model coefficients

- [`vcov(`*`<PLNfit>`*`)`](https://pln-team.github.io/PLNmodels/reference/vcov.PLNfit.md)
  :

  Calculate Variance-Covariance Matrix for a fitted
  [`PLN()`](https://pln-team.github.io/PLNmodels/reference/PLN.md) model
  object

- [`sigma(`*`<PLNfit>`*`)`](https://pln-team.github.io/PLNmodels/reference/sigma.PLNfit.md)
  : Extract variance-covariance of residuals 'Sigma'

- [`predict(`*`<PLNfit>`*`)`](https://pln-team.github.io/PLNmodels/reference/predict.PLNfit.md)
  : Predict counts of a new sample

- [`predict_cond()`](https://pln-team.github.io/PLNmodels/reference/predict_cond.md)
  : Predict counts conditionally

- [`fitted(`*`<PLNfit>`*`)`](https://pln-team.github.io/PLNmodels/reference/fitted.PLNfit.md)
  :

  Extracts model fitted values from objects returned by
  [`PLN()`](https://pln-team.github.io/PLNmodels/reference/PLN.md) and
  its variants

- [`standard_error()`](https://pln-team.github.io/PLNmodels/reference/standard_error.md)
  : Component-wise standard errors of B

- [`logLik(`*`<PLNfit>`*`)`](https://pln-team.github.io/PLNmodels/reference/logLik.PLNfit.md)
  : Extract log-likelihood of a fitted PLN model

- [`AIC(`*`<PLNfit>`*`)`](https://pln-team.github.io/PLNmodels/reference/AIC.PLNfit.md)
  : Akaike Information Criterion for a fitted PLN model

- [`BIC(`*`<PLNfit>`*`)`](https://pln-team.github.io/PLNmodels/reference/BIC.PLNfit.md)
  : Bayesian Information Criterion for a fitted PLN model

## Zero Inflated Poisson lognormal fit

Description of the ZIPLNfit object and methods for its manipulation.

- [`ZIPLNfit`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit.md)
  : An R6 Class to represent a ZIPLNfit

- [`ZIPLNfit_diagonal`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit_diagonal.md)
  : An R6 Class to represent a ZIPLNfit in a standard, general
  framework, with diagonal residual covariance

- [`ZIPLNfit_fixed`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit_fixed.md)
  : An R6 Class to represent a ZIPLNfit in a standard, general
  framework, with fixed (inverse) residual covariance

- [`ZIPLNfit_sparse`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit_sparse.md)
  : An R6 Class to represent a ZIPLNfit in a standard, general
  framework, with sparse inverse residual covariance

- [`ZIPLNfit_spherical`](https://pln-team.github.io/PLNmodels/reference/ZIPLNfit_spherical.md)
  : An R6 Class to represent a ZIPLNfit in a standard, general
  framework, with spherical residual covariance

- [`ZIPLN_param()`](https://pln-team.github.io/PLNmodels/reference/ZIPLN_param.md)
  : Control of a ZIPLN fit

- [`coef(`*`<ZIPLNfit>`*`)`](https://pln-team.github.io/PLNmodels/reference/coef.ZIPLNfit.md)
  : Extract model coefficients

- [`sigma(`*`<ZIPLNfit>`*`)`](https://pln-team.github.io/PLNmodels/reference/sigma.ZIPLNfit.md)
  : Extract variance-covariance of residuals 'Sigma'

- [`predict(`*`<ZIPLNfit>`*`)`](https://pln-team.github.io/PLNmodels/reference/predict.ZIPLNfit.md)
  : Predict counts of a new sample

- [`fitted(`*`<ZIPLNfit>`*`)`](https://pln-team.github.io/PLNmodels/reference/fitted.ZIPLNfit.md)
  :

  Extracts model fitted values from objects returned by
  [`ZIPLN()`](https://pln-team.github.io/PLNmodels/reference/ZIPLN.md)
  and its variants

- [`plot(`*`<ZIPLNfit_sparse>`*`)`](https://pln-team.github.io/PLNmodels/reference/plot.ZIPLNfit_sparse.md)
  :

  Extract and plot the network (partial correlation, support or inverse
  covariance) from a `ZIPLNfit_sparse` object

- [`logLik(`*`<ZIPLNfit>`*`)`](https://pln-team.github.io/PLNmodels/reference/logLik.ZIPLNfit.md)
  : Extract log-likelihood of a fitted ZIPLN model

- [`AIC(`*`<ZIPLNfit>`*`)`](https://pln-team.github.io/PLNmodels/reference/AIC.ZIPLNfit.md)
  : Akaike Information Criterion for a fitted ZIPLN model

- [`BIC(`*`<ZIPLNfit>`*`)`](https://pln-team.github.io/PLNmodels/reference/BIC.ZIPLNfit.md)
  : Bayesian Information Criterion for a fitted ZIPLN model

## Linear discriminant analysis via a Poisson lognormal fit

Description of the PLNLDAfit object and methods for its manipulation.

- [`PLNLDAfit`](https://pln-team.github.io/PLNmodels/reference/PLNLDAfit.md)
  : An R6 Class to represent a PLNfit in a LDA framework

- [`PLNLDAfit_diagonal`](https://pln-team.github.io/PLNmodels/reference/PLNLDAfit_diagonal.md)
  : An R6 Class to represent a PLNfit in a LDA framework with diagonal
  covariance

- [`PLNfit_diagonal`](https://pln-team.github.io/PLNmodels/reference/PLNfit_diagonal.md)
  [`PLNLDAfit_spherical`](https://pln-team.github.io/PLNmodels/reference/PLNfit_diagonal.md)
  : An R6 Class to represent a PLNfit in a standard, general framework,
  with diagonal residual covariance

- [`PLNLDA_param()`](https://pln-team.github.io/PLNmodels/reference/PLNLDA_param.md)
  : Control of a PLNLDA fit

- [`plot(`*`<PLNLDAfit>`*`)`](https://pln-team.github.io/PLNmodels/reference/plot.PLNLDAfit.md)
  :

  LDA visualization (individual and/or variable factor map(s)) for a
  `PLNPCAfit` object

- [`predict(`*`<PLNLDAfit>`*`)`](https://pln-team.github.io/PLNmodels/reference/predict.PLNLDAfit.md)
  : Predict group of new samples

- [`coef(`*`<PLNLDAfit>`*`)`](https://pln-team.github.io/PLNmodels/reference/coef.PLNLDAfit.md)
  :

  Extracts model coefficients from objects returned by
  [`PLNLDA()`](https://pln-team.github.io/PLNmodels/reference/PLNLDA.md)

## Poisson Lognormal PCA fit

Description of the PLNPCAfit and PLNPCAfamily objects and methods for
their manipulation.

- [`PLNPCAfit`](https://pln-team.github.io/PLNmodels/reference/PLNPCAfit.md)
  : An R6 Class to represent a PLNfit in a PCA framework

- [`PLNPCA_param()`](https://pln-team.github.io/PLNmodels/reference/PLNPCA_param.md)
  : Control of PLNPCA fit

- [`plot(`*`<PLNPCAfit>`*`)`](https://pln-team.github.io/PLNmodels/reference/plot.PLNPCAfit.md)
  :

  PCA visualization (individual and/or variable factor map(s)) for a
  `PLNPCAfit` object

- [`PLNPCAfamily`](https://pln-team.github.io/PLNmodels/reference/PLNPCAfamily.md)
  : An R6 Class to represent a collection of PLNPCAfit

- [`plot(`*`<PLNPCAfamily>`*`)`](https://pln-team.github.io/PLNmodels/reference/plot.PLNPCAfamily.md)
  : Display the criteria associated with a collection of PLNPCA fits (a
  PLNPCAfamily)

- [`getBestModel()`](https://pln-team.github.io/PLNmodels/reference/getBestModel.md)
  : Best model extraction from a collection of models

- [`getModel()`](https://pln-team.github.io/PLNmodels/reference/getModel.md)
  : Model extraction from a collection of models

## Mixture of Poisson Lognormal fits

Description of the PLNmixturefit and PLNmixturefamily objects and
methods for their manipulation.

- [`PLNmixturefit`](https://pln-team.github.io/PLNmodels/reference/PLNmixturefit.md)
  : An R6 Class to represent a PLNfit in a mixture framework

- [`PLNmixture_param()`](https://pln-team.github.io/PLNmodels/reference/PLNmixture_param.md)
  : Control of a PLNmixture fit

- [`plot(`*`<PLNmixturefit>`*`)`](https://pln-team.github.io/PLNmodels/reference/plot.PLNmixturefit.md)
  :

  Mixture visualization of a `PLNmixturefit` object

- [`PLNmixturefamily`](https://pln-team.github.io/PLNmodels/reference/PLNmixturefamily.md)
  : An R6 Class to represent a collection of PLNmixturefit

- [`plot(`*`<PLNmixturefamily>`*`)`](https://pln-team.github.io/PLNmodels/reference/plot.PLNmixturefamily.md)
  : Display the criteria associated with a collection of PLNmixture fits
  (a PLNmixturefamily)

- [`coef(`*`<PLNmixturefit>`*`)`](https://pln-team.github.io/PLNmodels/reference/coef.PLNmixturefit.md)
  : Extract model coefficients

- [`predict(`*`<PLNmixturefit>`*`)`](https://pln-team.github.io/PLNmodels/reference/predict.PLNmixturefit.md)
  :

  Prediction for a `PLNmixturefit` object

- [`sigma(`*`<PLNmixturefit>`*`)`](https://pln-team.github.io/PLNmodels/reference/sigma.PLNmixturefit.md)
  : Extract variance-covariance of residuals 'Sigma'

- [`fitted(`*`<PLNmixturefit>`*`)`](https://pln-team.github.io/PLNmodels/reference/fitted.PLNmixturefit.md)
  :

  Extracts model fitted values from objects returned by
  [`PLNmixture()`](https://pln-team.github.io/PLNmodels/reference/PLNmixture.md)
  and its variants

- [`getBestModel()`](https://pln-team.github.io/PLNmodels/reference/getBestModel.md)
  : Best model extraction from a collection of models

- [`getModel()`](https://pln-team.github.io/PLNmodels/reference/getModel.md)
  : Model extraction from a collection of models

## Sparse Poisson lognormal fit and network, w/o Zero Inflated component

Description of the (ZI)PLNnetworkfit and (ZI)PLNnetworkfamily objects
and methods for their manipulation.

- [`PLNnetworkfit`](https://pln-team.github.io/PLNmodels/reference/PLNnetworkfit.md)
  : An R6 Class to represent a PLNfit in a sparse inverse covariance
  framework

- [`PLNnetwork_param()`](https://pln-team.github.io/PLNmodels/reference/PLNnetwork_param.md)
  : Control of PLNnetwork fit

- [`ZIPLNnetwork_param()`](https://pln-team.github.io/PLNmodels/reference/ZIPLNnetwork_param.md)
  : Control of ZIPLNnetwork fit

- [`plot(`*`<PLNnetworkfit>`*`)`](https://pln-team.github.io/PLNmodels/reference/plot.PLNnetworkfit.md)
  :

  Extract and plot the network (partial correlation, support or inverse
  covariance) from a `PLNnetworkfit` object

- [`plot(`*`<ZIPLNfit_sparse>`*`)`](https://pln-team.github.io/PLNmodels/reference/plot.ZIPLNfit_sparse.md)
  :

  Extract and plot the network (partial correlation, support or inverse
  covariance) from a `ZIPLNfit_sparse` object

- [`Networkfamily`](https://pln-team.github.io/PLNmodels/reference/Networkfamily.md)
  : An R6 Class to virtually represent a collection of network fits

- [`ZIPLNnetworkfamily`](https://pln-team.github.io/PLNmodels/reference/ZIPLNnetworkfamily.md)
  : An R6 Class to represent a collection of ZIPLNnetwork

- [`PLNnetworkfamily`](https://pln-team.github.io/PLNmodels/reference/PLNnetworkfamily.md)
  :

  An R6 Class to represent a collection of `PLNnetworkfit`s

- [`plot(`*`<Networkfamily>`*`)`](https://pln-team.github.io/PLNmodels/reference/plot.Networkfamily.md)
  [`plot(`*`<PLNnetworkfamily>`*`)`](https://pln-team.github.io/PLNmodels/reference/plot.Networkfamily.md)
  [`plot(`*`<ZIPLNnetworkfamily>`*`)`](https://pln-team.github.io/PLNmodels/reference/plot.Networkfamily.md)
  :

  Display various outputs (goodness-of-fit criteria, robustness,
  diagnostic) associated with a collection of network fits (either
  `PLNnetworkfamily` or `ZIPLNnetworkfamily`)

- [`getBestModel()`](https://pln-team.github.io/PLNmodels/reference/getBestModel.md)
  : Best model extraction from a collection of models

- [`getModel()`](https://pln-team.github.io/PLNmodels/reference/getModel.md)
  : Model extraction from a collection of models

- [`coefficient_path()`](https://pln-team.github.io/PLNmodels/reference/coefficient_path.md)
  : Extract the regularization path of a PLNnetwork fit

- [`stability_selection()`](https://pln-team.github.io/PLNmodels/reference/stability_selection.md)
  : Compute the stability path by stability selection

- [`extract_probs()`](https://pln-team.github.io/PLNmodels/reference/extract_probs.md)
  : Extract edge selection frequency in bootstrap subsamples

## Other functions and objects

- [`prepare_data()`](https://pln-team.github.io/PLNmodels/reference/prepare_data.md)
  : Prepare data for use in PLN models
- [`compute_offset()`](https://pln-team.github.io/PLNmodels/reference/compute_offset.md)
  : Compute offsets from a count data using one of several normalization
  schemes
- [`PLNfamily`](https://pln-team.github.io/PLNmodels/reference/PLNfamily.md)
  : An R6 Class to represent a collection of PLNfit
- [`plot(`*`<PLNfamily>`*`)`](https://pln-team.github.io/PLNmodels/reference/plot.PLNfamily.md)
  : Display the criteria associated with a collection of PLN fits (a
  PLNfamily)
- [`rPLN()`](https://pln-team.github.io/PLNmodels/reference/rPLN.md) :
  PLN RNG
- [`ICL()`](https://pln-team.github.io/PLNmodels/reference/ICL.md) :
  Integrated Classification Likelihood
- [`compute_PLN_starting_point()`](https://pln-team.github.io/PLNmodels/reference/compute_PLN_starting_point.md)
  : Helper function for PLN initialization.
- [`compute_ZIPLN_starting_point()`](https://pln-team.github.io/PLNmodels/reference/compute_ZIPLN_starting_point.md)
  : Helper function for ZIPLN initialization.

## Data sets

- [`trichoptera`](https://pln-team.github.io/PLNmodels/reference/trichoptera.md)
  : Trichoptera data set
- [`microcosm`](https://pln-team.github.io/PLNmodels/reference/microcosm.md)
  : Cow microbiome data set
- [`oaks`](https://pln-team.github.io/PLNmodels/reference/oaks.md) :
  Oaks amplicon data set
- [`barents`](https://pln-team.github.io/PLNmodels/reference/barents.md)
  : Barents fish data set
- [`mollusk`](https://pln-team.github.io/PLNmodels/reference/mollusk.md)
  : Mollusk data set
- [`scRNA`](https://pln-team.github.io/PLNmodels/reference/scRNA.md) :
  Single cell RNA-seq data
