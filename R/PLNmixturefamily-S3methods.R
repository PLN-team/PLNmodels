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
#' @inheritParams plot.PLNfamily
#' @param x an R6 object with class [`PLNmixturefamily`]
#' @param type a character, either `"criteria"` or `"diagnostic"` for the type of plot.
#'
#' @inherit plot.PLNfamily return details
#' @return Produces either a diagnostic plot (with \code{type = 'diagnostic'}) or the evolution of the criteria
#' of the different models considered (with \code{type = 'criteria'}, the default).
#'
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myMixtures <- PLNmixture(Abundance ~ 1 + offset(log(Offset)),
#'            data = trichoptera, control = PLNmixture_param(smoothing = "none"))
#' plot(myMixtures, reverse = TRUE)
#' @export
plot.PLNmixturefamily <-
  function(x,
           type     = c("criteria", "diagnostic"),
           criteria = c("loglik", "BIC", "ICL"),
           reverse = FALSE,
           ...) {
  stopifnot(isPLNfamily(x))
  type <- match.arg(type)
  if (type == "criteria")
    p <- x$plot(criteria, reverse)
  if (type == "diagnostic")
    p <- x$plot_objective()

  p
}

#' @describeIn getModel Model extraction for [`PLNmixturefamily`]
#' @export
getModel.PLNmixturefamily <- function(Robject, var, index = NULL) {
  stopifnot(isPLNmixturefamily(Robject))
  Robject$getModel(var, index)
}

#' @describeIn getBestModel Model extraction for [`PLNmixturefamily`]
#' @export
getBestModel.PLNmixturefamily <- function(Robject, crit = c("ICL", "BIC"), ...) {
  stopifnot(isPLNmixturefamily(Robject))
  Robject$getBestModel(match.arg(crit))
}
