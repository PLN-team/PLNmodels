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
#' @param type a character, either "criteria" or "diagnostic" for the type of plot.
#' @param criteria vector of characters. The criteria to plot in c("loglik", "BIC", "ICL").
#' Default is  c("loglik", "BIC", "ICL").
#' @param reverse A logical indicating whether to plot the value of the criteria in the "natural" direction
#' (loglik - 0.5 penalty) or in the "reverse" direction (-2 loglik + penalty). Default to FALSE, i.e use the
#' natural direction, on the same scale as the log-likelihood..
#' @param ... additional parameters for S3 compatibility. Not used
#' @return Produces either a diagnostic plot (with \code{type = 'diagnostic'}) or the evolution of the criteria
#' of the different models considered (with \code{type = 'criteria'}, the default). The latter highlights the best
#' model in terms of BIC and ICL. These criteria have the form 'loglik - 1/2 * penalty'
#' so that they are on the same scale as the model log-likelihood. You can change this direction by setting
#' the parameter \code{reverse} to \code{TRUE}.
#' @examples
#' data(trichoptera)
#' trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
#' myMixtures <- PLNmixture(Abundance ~ 1 + offset(log(Offset)),
#'            data = trichoptera, control_main = list(iterates = 1))
#' plot(myMixtures)
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
