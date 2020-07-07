#' Fisher information matrix for Theta
#'
#' @description Extracts Fisher information matrix of \eqn{\Theta} from objects returned by \code{\link[=PLN]{PLN}} and its variants. Fisher matrix is computed using one of two approximation scheme: wald (default, conservative, gives large confidence interval) or louis (anticonservative). Note that the Fisher information matrix is the full-data version (scaled by the number of observations), usually noted \deqn{I_n(\theta)}.
#'
#' @param object an R6 object with class PLNfit
#' @param type Either `wald` (default) or `louis`. Approximation scheme used to compute the Fisher information matrix
#' @return A block-diagonal matrix with p (number of species) blocks of size d (number of covariates), assuming
#' \eqn{\Theta} is a matrix of size d * p.
#'
#' @seealso \code{\link[=standard_error.PLNfit]{standard_error}} for standard errors
#'
#' @export
fisher <- function(object, type) {
  warning("Deprecated: please use `vcov()` instead", call. = FALSE)
  UseMethod("fisher", object)
}

#' @describeIn fisher Fisher information matrix for PLNfit
#' @export
fisher.PLNfit <- function(object, type = c("wald", "louis")) {
  stopifnot(isPLNfit(object))
  type <- match.arg(type)
  if (type != object$fisher$type) {
    stop(paste("Fisher information was not computed using the", type, "approximation. Try another approximation scheme."))
  }
  object$fisher$mat
}
