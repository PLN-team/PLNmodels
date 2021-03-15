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
#' @param x an R6 object with class [`PLNfamily`]
#' @param criteria vector of characters. The criteria to plot in c("loglik", "BIC", "ICL").
#' Default is  c("loglik", "BIC", "ICL").
#' @param reverse A logical indicating whether to plot the value of the criteria in the "natural" direction
#' (loglik - 0.5 penalty) or in the "reverse" direction (-2 loglik + penalty). Default to FALSE, i.e use the
#' natural direction, on the same scale as the log-likelihood..
#' @param ... additional parameters for S3 compatibility. Not used
#'
#' @return Produces a plot  representing the evolution of the criteria of the different models considered,
#' highlighting the best model in terms of BIC and ICL. These criteria have the form 'loglik - 1/2 * penalty'
#' so that they are on the same scale as the model log-likelihood. You can change this direction by setting
#' the parameter \code{reverse} to \code{TRUE}.
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
