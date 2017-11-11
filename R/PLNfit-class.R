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
#' @field model_par a list with two matrices, B and Theta, which are the estimated parameters of the pPCA model
#' @field var_par a list with two matrices, M and S, which are the estimated parameters in the variational approximation
#' @field criteria a named vector with the value of some criteria (variational lower bound J, BIC, ICL, R2, lmin and lmax) for the different models.
#' @field convergence quantities useful for monitoring the optimization
#' @include PLNfit-class.R
#' @importFrom R6 R6Class
#' @importFrom corrplot corrplot
PLNfit <-
   R6Class(classname = "PLNfit",
    public = list(
      ## constructor
      initialize = function(Theta=NA, Sigma=NA, Omega=NA, M=NA, S=NA,
                            J=NA, BIC=NA, ICL=NA, R2 = NA, monitoring=NA) {
        private$Theta      <- Theta
        private$Sigma      <- Sigma
        private$Omega      <- Omega
        private$M          <- M
        private$S          <- S
        private$J          <- J
        private$BIC        <- BIC
        private$ICL        <- ICL
        private$R2         <- R2
        private$monitoring <- monitoring
      },
      ## "setter" function
      update = function(Theta=NA, Sigma=NA, Omega=NA, M=NA, S=NA,
                        J=NA, BIC=NA, ICL=NA, R2 = NA, monitoring = NA) {
        if (!anyNA(Theta))      private$Theta  <- Theta
        if (!anyNA(Sigma))      private$Sigma  <- Sigma
        if (!anyNA(Omega))      private$Omega  <- Omega
        if (!anyNA(M))          private$M      <- M
        if (!anyNA(S))          private$S      <- S
        if (!anyNA(J))          private$J      <- J
        if (!anyNA(BIC))        private$BIC    <- BIC
        if (!anyNA(ICL))        private$ICL    <- ICL
        if (!anyNA(R2))         private$R2     <- R2
        if (!anyNA(monitoring)) private$monitoring <- monitoring
      }
    ),
    private = list(
      S          = NULL, # the n x p variational parameters for the variances
      M          = NULL, # the n x p variational parameters for the means
      Theta      = NULL, # the p x d model parameters for the covariable
      Sigma      = NULL, # the p x p covariance matrix
      Omega      = NULL, # the p x p precision matrix
      J          = NULL, # the variational lower bound of the likelihood
      BIC        = NULL, # (variational) Baysesian information criterion
      ICL        = NULL, # (variational) Integrated classification criterion
      R2         = NULL, # approximated goodness of fit criterion
      monitoring = NULL  # a list with optimization monitoring quantities
    ),
    ## use active bindings to access private members like fields
    active = list(
      model_par = function() {
        list(Theta = private$Theta, Sigma = private$Sigma, Omega = private$Omega)
      },
      var_par = function() {
        list(M = private$M, S = private$S)
      },
      convergence = function() {
        private$monitoring
      },
      criteria = function() {
        c(J = private$J, BIC = private$BIC, ICL = private$ICL, R2 = private$R2)
      }
    )
  )

#' Display the model parameters in a matrix fashion
#'
#' @name PLNfit_plot_par
#'
#' @param type character. should the variational or the model parameters be plotted? default is "model".
#'
PLNfit$set("public", "plot_par",
  function(type=c("model","variational")) {
    type <- match.arg(type)
    param <- switch(type,
               "model"       = self$model_par,
               "variational" = self$var_par)
    par1 <- param[[1]]; par2 <- param[[2]]
    rownames(par1) <- rep(" ", nrow(par1)) ; colnames(par1) <- rep(" ", ncol(par1))
    rownames(par2) <- rep(" ", nrow(par2)) ; colnames(par2) <- rep(" ", ncol(par2))

    par(mfrow = c(2,2))
    hist(par1, breaks = sqrt(nrow(par1)), xlab = "", ylab = "", main = paste0(names(param)[1]))
    hist(par2, breaks = sqrt(nrow(par2)), xlab = "", ylab = "", main = paste0(names(param)[2]))
    corrplot::corrplot(par1, is.corr = FALSE, method = "color", cl.pos = "n")
    corrplot::corrplot(par2, is.corr = FALSE, method = "color", cl.pos = "n")
    title(main = paste0("\n",type," parameters"), outer = TRUE)
    par(mfrow = c(1,1))
  }
)

#' Positions in the (Euclidian) parameter space, noted as Z in the model. Used to compute the likelihood.
#'
#' @name PLNfit_latentPos
#'
#' @param covariates a matrix of covariates. Will usually be extracted from the corresponding field in PLNfamily-class
#' @param offsets    a matrix of offsets. Will usually be extracted from the corresponding field in PLNfamily-class
#'
PLNfit$set("public", "latentPos",
function(covariates, offsets) {
  latentPos <- private$M + tcrossprod(covariates, private$Theta) + offsets
  latentPos
})

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
    stopifnot(ncol(newdata)    == ncol(private$Theta),
              nrow(newdata)    == nrow(newOffsets),
              ncol(newOffsets) == nrow(private$Theta))
    ## Mean latent positions in the parameter space
    Z <- tcrossprod(newdata, private$Theta) + newOffsets
    results <- switch(type,
                      link     = Z,
                      response = exp(Z))
    ## output formatting
    rownames(results) <- rownames(newdata); colnames(results) <- rownames(self$model.par$Theta)
    attr(results, "type") <- type
    results
  }
)
