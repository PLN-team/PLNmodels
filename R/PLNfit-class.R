## TODO
## o add a predict method and similar stat tools

#' An R6 Class to represent a PLNfit in a standard, general framework
#'
#' @description The function \code{\link{PLN}} produces a collection of models which are instances of object with class \code{PLNfit}.
#' Objects produced by the functions \code{\link{PLNnetwork}} and \code{\link{PLNPCA}} also enjoy the method of \code{\link{PLNfit}}
#' by inheritance.
#'
#' This class comes with a set of methods, some of them being useful for the user: plot_model_par, plot_variational_par
#'
#' Fields should not be changed or manipulated by the user as they are updated internally.
#'
#' @field model.par a list with two matrices, B and Theta, which are the estimated parameters of the pPCA model
#' @field variation.par a list with two matrices, M and S, which are the estimated parameters in the variational approximation
#' @field criteria a named vector with the value of some criteria (variational lower bound J, BIC, ICL, R2, lmin and lmax) for the different models.
#' @field convergence quantities usefull for monitoring the optimization
#' @include PLNfit-class.R
#' @importFrom R6 R6Class
#' @importFrom corrplot corrplot
PLNfit <-
   R6Class(classname = "PLNfit",
    public = list(
      model.par       = NULL, # Theta and Sigma
      variational.par = NULL, # M and S
      criteria        = NULL, # J, BIC, ICL
      convergence     = NULL, # results of optimization
      initialize = function(model.par=NA, variational.par=NA, criteria=NA, convergence=NA, loglik=NA) {
      self$model.par       <- model.par
      self$variational.par <- variational.par
      self$criteria        <- criteria
      self$convergence     <- convergence
      }
    )
  )

PLNfit$set("public", "plot_variational_par",
  function() {
    S <- self$variational.par$S; M <- self$variational.par$M
    rownames(S) <- rep(" ", nrow(S)) ; colnames(S) <- rep(" ", ncol(S))
    rownames(M) <- rep(" ", nrow(M)) ; colnames(M) <- rep(" ", ncol(M))
    par(mfrow=c(2,2))
    hist(M, breaks=nrow(M), xlab="", ylab="", main="means")
    hist(S, breaks=nrow(S), xlab="", ylab="", main="standard deviations")
    corrplot(M, is.corr = FALSE, method="color", cl.pos = "n")
    corrplot(S, is.corr = FALSE, method="color", cl.pos = "n")
    title(main="\nVariational parameters", outer=TRUE)
    par(mfrow=c(1,1))
  }
)

PLNfit$set("public", "plot_model_par",
  function() {
    Theta <- t(self$model.par$Theta); Sigma <- self$model.par$Sigma
    rownames(Theta) <- rep(" ", nrow(Theta)) ; colnames(Theta) <- rep(" ", ncol(Theta))
    rownames(Sigma) <- rep(" ", nrow(Sigma)) ; colnames(Sigma) <- rep(" ", ncol(Sigma))
    par(mfrow=c(2,2))
    hist(Theta, breaks=nrow(Theta), xlab="", ylab="", main="Regression parameters")
    hist(Sigma, breaks=nrow(Sigma), xlab="", ylab="", main="Covariance matrix")
    corrplot(Theta, is.corr = FALSE, method="color", cl.pos = "n")
    corrplot(Sigma, is.corr = FALSE, method="color", cl.pos = "n")
    title(main="\nModel parameters", outer=TRUE)
    par(mfrow=c(1,1))
  }
)
