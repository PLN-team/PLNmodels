#' An R6 Class to represent a collection of PLNPCAfit
#'
#' @description The function \code{\link{PLNPCA}} produces an instance of this class.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for \code{\link[=getBestModel.PLNPCAfamily]{getBestModel}},
#' \code{\link[=getModel.PLNPCAfamily]{getModel}} and  \code{\link[=plot.PLNPCAfamily]{plot}}.
#'
#' @field responses the matrix of responses common to every models
#' @field covariates the matrix of covariates common to every models
#' @field offsets the matrix of offsets common to every models
#' @field ranks the dimensions of the successively fitted models
#' @field models a list of \code{\link[=PLNPCAfit]{PLNPCAfit}} object, one per rank.
#' @field inception a \code{\link[=PLNfit]{PLNfit}} object, obtained when full rank is considered.
#' @field criteria a data frame with the value of some criteria (variational lower bound J, BIC, ICL and R2) for the different models.
#' @include PLNfamily-class.R
#' @importFrom R6 R6Class
#' @import ggplot2
#' @seealso The function \code{\link{PLNPCA}}, the class \code{\link[=PLNPCAfit]{PLNPCAfit}}
PLNPCAfamily <-
  R6Class(classname = "PLNPCAfamily",
    inherit = PLNfamily,
     active = list(
      ranks = function() private$params
    )
)

PLNPCAfamily$set("public", "initialize",
  function(ranks, responses, covariates, offsets, weights, model, control) {

  ## initialize the required fields
  super$initialize(responses, covariates, offsets, weights, control)
  private$params <- ranks

  ## instantiate as many models as ranks
  self$models <- lapply(ranks, function(rank){
    model <- PLNPCAfit$new(rank, responses, covariates, offsets, weights, model, control)
    model
  })
})

PLNPCAfamily$set("public", "optimize",
  function(control) {
    self$models <- mclapply(self$models, function(model) {
    if (control$trace == 1) {
      cat("\t Rank approximation =",model$rank, "\r")
      flush.console()
    }
    if (control$trace > 1) {
      cat(" Rank approximation =",model$rank)
      cat("\n\t conservative convex separable approximation for gradient descent")
    }
    model$optimize(self$responses, self$covariates, self$offsets, self$weights, control)
    model
  }, mc.cores = control$cores, mc.allow.recursive = FALSE)
})


#' Display the criteria associated with a collection of PLNPCA fits (a PLNPCAfamily)
#'
#' @name plot.PLNPCAfamily
#'
#' @param x an R6 object with class PLNfamily
#' @param criteria vector of characters. The criteria to plot in c("loglik", "BIC", "ICL", "R_squared").
#' Default is  c("loglik", "BIC", "ICL").
#' @param annotate logical: should the value of approximated R squared be added to the plot?
#' @param ... additional parameters for S3 compatibility. Not used
#'
#' @return Produces a plot  representing the evolution of the criteria of the different models considered,
#' highlighting the best model in terms of BIC and ICL.
#'
#' @export
plot.PLNPCAfamily <- function(x, criteria = c("loglik", "BIC", "ICL"), annotate = TRUE, ...) {
  stopifnot(isPLNfamily(x))
  x$plot(criteria, annotate)
}

#' Best model extraction from a collection of PLNPCAfit
#'
#' @name getBestModel.PLNPCAfamily
#'
#' @param Robject an object with classi PLNPCAfamilly
#' @param crit a character for the criterion used to performed the selection. Either
#' "BIC", "ICL", "R_squared". Default is \code{ICL}.
#' @param ... not use.
#' @return  Send back a object with class \code{\link[=PLNPCAfit]{PLNPCAfit}}.
#'
#' @export
getBestModel.PLNPCAfamily <- function(Robject, crit = c("ICL", "BIC", "R_squared"), ...) {
  stopifnot(isPLNPCAfamily(Robject))
  Robject$getBestModel(match.arg(crit))
}

PLNPCAfamily$set("public", "getBestModel",
function(crit = c("BIC", "ICL", "R_squared")){
  crit <- match.arg(crit)
  stopifnot(!anyNA(self$criteria[[crit]]))
  id <- 1
  if (length(self$criteria[[crit]]) > 1) {
    id <- which.max(self$criteria[[crit]])
  }
  model <- self$models[[id]]$clone()
  model
})

#' Model extraction from a collection of PLNPCA models
#'
#' @name getModel.PLNPCAfamily
#'
#' @param Robject an R6 object with class PLNfamily
#' @param var value of the parameter (rank for PCA) that identifies the model to be extracted from the collection. If no exact match is found, the model with closest parameter value is returned with a warning.
#' @param index Integer index of the model to be returned. Only the first value is taken into account.
#'
#' @return Sends back a object with class \code{\link[=PLNPCAfit]{PLNPCAfit}}.
#'
#' @export
getModel.PLNPCAfamily <- function(Robject, var, index = NULL) {
  stopifnot(isPLNPCAfamily(Robject))
  Robject$getModel(var, index = NULL)
}

PLNPCAfamily$set("public", "plot",
function(criteria = c("loglik", "BIC", "ICL"), annotate = TRUE) {
  vlines <- sapply(intersect(criteria, c("BIC", "ICL")) , function(crit) self$getBestModel(crit)$rank)
  p <- super$plot(criteria, annotate) + xlab("rank") + geom_vline(xintercept = vlines, linetype = "dashed", alpha = 0.25)
  p
})

PLNPCAfamily$set("public", "show",
function() {
  super$show()
  cat(" Task: Principal Component Analysis\n")
  cat("========================================================\n")
  cat(" - Ranks considered: from", min(self$ranks), "to", max(self$ranks),"\n")
  cat(" - Best model (regarding BIC): rank =", self$getBestModel("BIC")$rank, "- R2 =", round(self$getBestModel("BIC")$R_squared, 2), "\n")
  cat(" - Best model (regarding ICL): rank =", self$getBestModel("ICL")$rank, "- R2 =", round(self$getBestModel("ICL")$R_squared, 2), "\n")
})

PLNPCAfamily$set("public", "print", function() self$show())
