## =========================================================================================
##
## PUBLIC S3 METHODS FOR PLNfamily
##
## =========================================================================================

## Auxiliary functions to check the given class of an objet
isPLNfamily <- function(Robject) {inherits(Robject, 'PLNfamily'       )}

## S3 methods declared generic for children of PLNfamily

#' Best model extraction from a collection of models
#'
#' @param Robject an object with class PLNPCAfamilly ot PLNnetworkfamily
#' @param crit a character for the criterion used to performed the selection. Either
#' "BIC", "ICL", "EBIC", "StARS", "R_squared". Default is `ICL` for `PLNPCA`, and `BIC` for `PLNnetwork`.
#'  If StARS (Stability Approach to Regularization Selection) is chosen and stability selection
#'  was not yet performed, the function will call the method [stability_selection()] with default argument.
#' @param ... additional parameters for StARS criterion (only for `PLNnetwork`). `stability`, a scalar indicating the target stability (= 1 - 2 beta) at which the network is selected. Default is \code{0.9}.
#' @return  Send back an object with class [`PLNPCAfit`] or [`PLNnetworkfit`]
#'
#' @examples
#' \dontrun{
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPCA <- PLNPCA(Abundance ~ 1 + offset(log(Offset)), data = trichoptera, ranks = 1:4)
#' myModel <- getBestModel(myPCA)
#' }
#' @export
getBestModel <- function(Robject, crit, ...) {UseMethod("getBestModel", Robject)}

#' Model extraction from a collection of models
#'
#' @param Robject an R6 object with class [`PLNPCAfamily`] or [`PLNnetworkfamily`]
#' @param var value of the parameter (`rank` for PLNPCA, `sparsity` for PLNnetwork) that identifies the model to be extracted from the collection. If no exact match is found, the model with closest parameter value is returned with a warning.
#' @param index Integer index of the model to be returned. Only the first value is taken into account.
#'
#' @return Sends back an object with class [`PLNPCAfit`] or [`PLNnetworkfit`].
#' @examples
#' \dontrun{
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPCA <- PLNPCA(Abundance ~ 1 + offset(log(Offset)), data = trichoptera, ranks = 1:5)
#' myModel <- getModel(myPCA, 2)
#' }
#' @export
getModel <- function(Robject, var , index) {UseMethod("getModel"    , Robject)}
