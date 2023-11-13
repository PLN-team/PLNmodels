## =========================================================================================
##
## PUBLIC S3 METHODS FOR PLNblockfamily
##
## =========================================================================================

## Auxiliary functions to check the given class of an objet
isPLNblockfamily     <- function(Robject) {inherits(Robject, "PLNblockfamily"    )}

#' Display the criteria associated with a collection of PLNblock fits (a PLNblockfamily)
#'
#' @name plot.PLNblockfamily
#'
#' @param x an R6 object with class [`PLNblockfamily`]
#' @inheritParams plot.PLNfamily
#' @inherit plot.PLNfamily return details
#'
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myBlocks <- PLNblocks(Abundance ~ 1 + offset(log(Offset)), data = trichoptera, nb_blocks = 1:5)
#' \dontrun{
#' plot(myBlocks)
#' }
#' @export
plot.PLNblockfamily <- function(x, criteria = c("loglik", "BIC", "ICL"), reverse = FALSE, ...) {
  stopifnot(isPLNfamily(x))
  x$plot(criteria, reverse)
}

#' @describeIn getModel Model extraction for [`PLNblockfamily`]
#' @export
getModel.PLNblockfamily <- function(Robject, var, index = NULL) {
  stopifnot(isPLNblockfamily(Robject))
  Robject$getModel(var, index)
}

#' @describeIn getBestModel Model extraction for [`PLNblockfamily`]
#' @export
getBestModel.PLNblockfamily <- function(Robject, crit = c("ICL", "BIC"), ...) {
  stopifnot(isPLNblockfamily(Robject))
  Robject$getBestModel(match.arg(crit))
}
