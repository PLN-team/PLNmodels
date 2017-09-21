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
#' @field convergence quantities useful for monitoring the optimization
#' @include PLNfit-class.R
#' @importFrom R6 R6Class
#' @importFrom corrplot corrplot
PLNfit <-
   R6Class(classname = "PLNfit",
    public = list(
      model.par       = NULL, # Theta and Sigma (and/or B)
      variational.par = NULL, # M and S
      criteria        = NULL, # J, BIC, ICL
      convergence     = NULL, # results of optimization
      initialize = function(model.par=NA, variational.par=NA, criteria=NA, convergence=NA) {
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
    hist(S, breaks=nrow(S), xlab="", ylab="", main="variances")
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
    hist(as.numeric(Theta), breaks=ncol(Theta), xlab="", ylab="", main="Regression parameters")
    hist(as.numeric(Sigma), breaks=nrow(Sigma), xlab="", ylab="", main="Covariance matrix")
    corrplot(Theta, is.corr = FALSE, method="color", cl.pos = "n")
    corrplot(Sigma, is.corr = FALSE, method="color", cl.pos = "n")
    title(main="\nModel parameters", outer=TRUE)
    par(mfrow=c(1,1))
  }
)

## TODO accessors for the variational and model parameters? Or at least coefficients and Sigma (but vcov not a good name)?

#' Positions in the (Euclidian) parameter space, noted as Z in the model. Used to compute the likelihood.
#'
#' @name PLNfit_latentPos
#'
#' @param covariates a matrix of covariates. Will usually be extracted from the corresponding field in PLNfamily-class
#' @param offsets    a matrix of offsets. Will usually be extracted from the corresponding field in PLNfamily-class
#'
PLNfit$set("public", "latentPos",
function(covariates, offsets) {
  latentPos <- self$variational.par$M + tcrossprod(covariates, self$model.par$Theta) + offsets
  latentPos
})

PLNfit$set("public", "addCriteria",
  function(crit.name, value) {
    if (crit.name %in% names(self$criteria)) {
      self$criteria[crit.name] <- value
    } else {
      self$criteria <- setNames(object = c(self$criteria, value),
                                nm     = c(names(self$criteria), crit.name))
    }
  }
)

#' Predict counts of a new sample
#'
#' @name PLNfit_predict
#'
#' @param newdata    A data frame in which to look for variables with which to predict.
#' @param newOffsets A matrix in which to look for offsets with which to predict.
#' @param type       The type of prediction required. The default is on the scale of the linear predictors (i.e. log average count);
#'                   the alternative "response" is on the scale of the response variable (i.e. average count)
#' @return A matrix of predicted log-counts (if type = "link") or predicted counts (if type = "response").
#'
PLNfit$set("public", "predict",
  function(newdata, newOffsets, type = c("link", "response")) {
    type = match.arg(type)
    ## Are matrix conformable?
    stopifnot(ncol(newdata)    == ncol(self$model.par$Theta),
              nrow(newdata)    == nrow(newOffsets),
              ncol(newOffsets) == nrow(self$model.par$Theta))
    ## Mean latent positions in the parameter space
    Z <- tcrossprod(newdata, self$model.par$Theta) + newOffsets
    results <- switch(type,
                      link     = Z,
                      response = exp(Z))
    ## output formatting
    rownames(results) <- rownames(newdata); colnames(results) <- rownames(self$model.par$Theta)
    attr(results, "type") <- type
    return(results)
  }
)
