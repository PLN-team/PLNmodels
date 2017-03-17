##' @title Fit a Poisson lognormal model towards Principal Component Analysis
##'
##' @description two methods are available for specifing the models (with formulas or matrices)
##'
##' @param formula a formula
##' @param Y a (n x p) matrix of count data
##' @param X an optional (n x d) matrix of covariates. SHould include the intercept (a column of one) if the default method is used.
##' @param O an optional (n x p) matrix of offsets.
##' @param Q a vector of integer containing the successive ranks (or number of axes to be considered)
##' @param control a list for controling the optimization. See details.
##' @return a RefClass object with class \code{\link[=PLNfamily-class]{PLNfamily}}, which contains
##' a collection of models with class \code{\link[=PLNfit.PCA-class]{PLNfit.PCA}}
##'
##' @details The parameter \code{control} is a list with the following entries
##' \itemize{
##'  \item{"factr"}{controls the L-BFGF-B procedure. See the documentation of \code{\link{optim}}. Default 1e8. Decrease if you experience instability or non monotonous J as a function of the rank}
##'  \item{"pgtol"}{controls the L-BFGF-B procedure. See the documentation of \code{\link{optim}}. Default 1e-2. Decrease if you experience instability or non monotonous J as a function of the rank}
##'  \item{"maxit"}{controls the L-BFGF-B procedure. See the documentaiton of \code{\link{optim}}. Default is 20000}
##'  \item{"lb.var"}{the minimum admissible value fr the variance parameter S in the variational approximation. Default is 1e-3.}
##'  \item{"cores"}{the number of cores. If Q has many entries, you might consider multiple cores. Default is 1.}
##'  \item{"approx"}{the model used for the covariance matrix in the variational Gaussian approximation. Either "diagonal" or "spherical". Default is "diagonal"}
##'  \item{"trace"}{integer for verbosity. Useless when \code{cores} > 1}
##' }
##'
##' @rdname PLNPCA
##' @examples
##' ## See the vignette: vignette("trichoptera", package="PLNmodels")
##' @seealso The classes \code{\link[=PLNfamily-class]{PLNfamily}} and \code{\link[=PLNfit.PCA-class]{PLNfit.PCA}}
##' @export
PLNPCA <- function(x, ...)
  UseMethod("PLNPCA", x)

##' @rdname PLNPCA
##' @export
PLNPCA.formula <- function(formula, Q = 1:5,  control = list()) {

  frame  <- model.frame(formula)
  Y      <- model.response(frame)
  X      <- model.matrix(formula)
  O      <- model.offset(frame)
  if (is.null(O)) O <- matrix(0, nrow(Y), ncol(Y))

  return(PLNPCA.default(Y, X, O, Q, control))
}

##' @rdname PLNPCA
##' @export
PLNPCA.default <- function(Y, X = cbind(rep(1, nrow(Y))), O = matrix(0, nrow(Y), ncol(Y)), Q = 1:5,  control = list()) {

  ## define default control parameters for optim and overwrite by user defined parameters
  ctrl <- list(factr=1e8, pgtol=1e-2, maxit=20000, lbvar=1e-3, cores=1, approx="diagonal", trace=1)
  ctrl[names(control)] <- control

  ## Instantiate the collection of PLN models, initialized by glm Poisson
  if (ctrl$trace > 0) cat("\n Initialization...")
  myPLN <- PLNfamily$new(ranks=Q, type=ctrl$approx, responses=Y, covariates=X, offsets=O)

  ## Now adjust the PLN models
  if (ctrl$trace > 0) cat("\n Adjusting",ctrl$approx,"models")
  myPLN$optimize(ctrl)
  if (ctrl$trace > 0) cat("\n DONE!\n")
  myPLN$setCriteria()

  ## turn the PLN models to PLN-PCA model (set vizualisation)
  myPLN$setPCA()

  return(myPLN)
}

