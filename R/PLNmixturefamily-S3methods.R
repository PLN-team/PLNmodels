## =========================================================================================
##
## PUBLIC S3 METHODS FOR PLNmixturefamily
##
## =========================================================================================

## Auxiliary functions to check the given class of an objet
isPLNmixturefamily     <- function(Robject) {inherits(Robject, "PLNmixturefamily"    )}

#' Display the criteria associated with a collection of PLNmixture fits (a PLNmixturefamily)
#'
#' @name plot.PLNmixturefamily
#'
#' @param x an R6 object with class [`PLNfamily`]
#' @param criteria vector of characters. The criteria to plot in c("loglik", "BIC", "ICL", "R_squared").
#' Default is  c("loglik", "BIC", "ICL").
#' @param ... additional parameters for S3 compatibility. Not used
#'
#' @return Produces a plot  representing the evolution of the criteria of the different models considered,
#' highlighting the best model in terms of BIC and ICL (the greater, the better).
#' These criteria have the form 'loglik - 1/2 * penalty' so that they are on the same scale as the model loglikelihood.
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myMixtures  <- PLNmixture(Abundance ~ 1 + offset(log(Offset)), data = trichoptera, clusters = 1:5)
#' \dontrun{
#' plot(myMixtures)
#' }
#' @export
plot.PLNmixturefamily <- function(x, criteria = c("loglik", "BIC", "ICL"), ...) {
  stopifnot(isPLNfamily(x))
  x$plot(criteria)
}

#' @describeIn getModel Model extraction for [`PLNmixturefamily`]
#' @export
getModel.PLNmixturefamily <- function(Robject, var, index = NULL) {
  stopifnot(isPLNmixturefamily(Robject))
  Robject$getModel(var, index)
}

#' @describeIn getBestModel Model extraction for [`PLNmixturefamily`]
#' @export
getBestModel.PLNmixturefamily <- function(Robject, crit = c("ICL", "BIC", "R_squared"), ...) {
  stopifnot(isPLNmixturefamily(Robject))
  Robject$getBestModel(match.arg(crit))
}
