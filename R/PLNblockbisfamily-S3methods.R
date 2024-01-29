## =========================================================================================
##
## PUBLIC S3 METHODS FOR PLNblockbisfamily
##
## =========================================================================================

## Auxiliary functions to check the given class of an objet
isPLNblockbisfamily <- function(Robject) {inherits(Robject, "PLNblockbisfamily")}

#' Display the criteria associated with a collection of PLNblockbis fits (a PLNblockbisfamily)
#'
#' @name plot.PLNblockbisfamily
#'
#' @param x an R6 object with class [`PLNblockbisfamily`]
#' @inheritParams plot.PLNfamily
#' @inherit plot.PLNfamily return details
#'
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myBlocks <- PLNblockbis(Abundance ~ 1 + offset(log(Offset)), data = trichoptera, nb_blocks = 1:5)
#' \dontrun{
#' plot(myBlocks)
#' }
#' @export
plot.PLNblockbisfamily <- function(x, criteria = c("loglik", "BIC", "ICL"), reverse = FALSE, ...) {
  stopifnot(isPLNfamily(x))
  x$plot(criteria, reverse)
}

#' @describeIn getModel Model extraction for [`PLNblockbisfamily`]
#' @export
getModel.PLNblockbisfamily <- function(Robject, var, index = NULL) {
  stopifnot(isPLNblockbisfamily(Robject))
  Robject$getModel(var, index)
}

#' @describeIn getBestModel Model extraction for [`PLNblockbisfamily`]
#' @export
getBestModel.PLNblockbisfamily <- function(Robject, crit = c("ICL", "BIC"), ...) {
  stopifnot(isPLNblockbisfamily(Robject))
  Robject$getBestModel(match.arg(crit))
}
