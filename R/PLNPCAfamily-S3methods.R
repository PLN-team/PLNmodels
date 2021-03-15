## =========================================================================================
##
## PUBLIC S3 METHODS FOR PLNPCAfamily
##
## =========================================================================================

## Auxiliary functions to check the given class of an objet
isPLNPCAfamily     <- function(Robject) {inherits(Robject, "PLNPCAfamily"    )}

#' Display the criteria associated with a collection of PLNPCA fits (a PLNPCAfamily)
#'
#' @name plot.PLNPCAfamily
#'
#' @param x an R6 object with class [`PLNPCAfamily`]
#' @inheritParams plot.PLNfamily
#' @inherit plot.PLNfamily return details
#'
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myPCAs <- PLNPCA(Abundance ~ 1 + offset(log(Offset)), data = trichoptera, ranks = 1:5)
#' \dontrun{
#' plot(myPCAs)
#' }
#' @export
plot.PLNPCAfamily <- function(x, criteria = c("loglik", "BIC", "ICL"), reverse = FALSE, ...) {
  stopifnot(isPLNfamily(x))
  x$plot(criteria, reverse)
}

#' @describeIn getModel Model extraction for [`PLNPCAfamily`]
#' @export
getModel.PLNPCAfamily <- function(Robject, var, index = NULL) {
  stopifnot(isPLNPCAfamily(Robject))
  Robject$getModel(var, index)
}

#' @describeIn getBestModel Model extraction for [`PLNPCAfamily`]
#' @export
getBestModel.PLNPCAfamily <- function(Robject, crit = c("ICL", "BIC"), ...) {
  stopifnot(isPLNPCAfamily(Robject))
  Robject$getBestModel(match.arg(crit))
}
