#' An R6 Class to represent a PLNfit in a standard, general framework
#'
#' @description The function \code{\link{PLN}} produces a collection of models which are instances of object with class \code{PLNfit}.
#' Objects produced by the functions \code{\link{PLNnetwork}} and \code{\link{PLNPCA}} also enjoy the method of \code{\link{PLNfit}}
#' by inheritance.
#'
#' This class comes with a set of methods, some of them being useful for the user: plot_model, plot_variational_par
#'
#' Fields are accessed via active binding and cannot be changed by the user.
#'
#' @field model_par a list with two matrices, B and Theta, which are the estimated parameters of the pPCA model
#' @field var_par a list with two matrices, M and S, which are the estimated parameters in the variational approximation
#' @field optim_par a list with parameters useful for monitoring the optimization
#' @field loglik variational lower bound of the loglikelihood
#' @field BIC variational lower bound of the BIC
#' @field ICL variational lower bound of the ICL
#' @field R_squared approximated goodness-of-fit criterion
#' @field degrees_freedom number of parameters in the current PLN model
#' @field criteria a vector with loglik, BIC, ICL, R_squared and degrees of freedom
#' @include PLNfit-class.R
#' @importFrom R6 R6Class
#' @importFrom corrplot corrplot
PLNfit <-
   R6Class(classname = "PLNfit",
    public = list(
      ## constructor
      initialize = function(Theta=NA, Sigma=NA, M=NA, S=NA, J=NA, monitoring=NA) {
        private$Theta      <- Theta
        private$Sigma      <- Sigma
        private$M          <- M
        private$S          <- S
        private$J          <- J
        private$monitoring <- monitoring
      },
      ## "setter" function
      update = function(Theta=NA, Sigma=NA, M=NA, S=NA, J=NA, R2=NA, monitoring=NA) {
        if (!anyNA(Theta))      private$Theta  <- Theta
        if (!anyNA(Sigma))      private$Sigma  <- Sigma
        if (!anyNA(M))          private$M      <- M
        if (!anyNA(S))          private$S      <- S
        if (!anyNA(J))          private$J      <- J
        if (!anyNA(R2))         private$R2     <- R2
        if (!anyNA(monitoring)) private$monitoring <- monitoring
      }
    ),
    private = list(
      Theta      = NA, # the p x d model parameters for the covariable
      Sigma      = NA, # the p x p covariance matrix
      S          = NA, # the n x p variational parameters for the variances
      M          = NA, # the n x p variational parameters for the means
      R2         = NA, # approximated goodness of fit criterion
      J          = NA, # approximated loglikelihood
      monitoring = NA  # a list with optimization monitoring quantities
    ),
    ## use active bindings to access private members like fields
    active = list(
      n = function() {nrow(private$M)},
      q = function() {ncol(private$M)},
      p = function() {nrow(private$Theta)},
      d = function() {ncol(private$Theta)},
      model_par = function() {list(Theta = private$Theta, Sigma = private$Sigma)},
      var_par   = function() {list(M = private$M, S = private$S)},
      optim_par = function() {private$monitoring},
      degrees_freedom = function() {
        private$p * private$d + private$p * (private$p + 1)/2
      },
      loglik    = function() {private$J},
      BIC       = function() {self$loglik - .5 * log(self$n) * self$degrees_freedom},
      ICL       = function() {self$BIC - .5 * (self$n * self$q * log(2*pi*exp(1)) + sum(log(private$S)))},
      R_squared = function() {private$R2},
      criteria  = function() {c(degrees_freedom = self$degrees_freedom, loglik = self$loglik, BIC = self$BIC, ICL = self$ICL, R_squared = self$R_squared)}
    )
  )
## an S3 function to check if an object is a PLNfit
isPLNfit <- function(Robject) {all.equal(class(Robject), c('PLNfit', 'R6'))}

## ----------------------------------------------------------------------
## PUBLIC METHODS FOR INTERNAL USE -> PLNfamily
## ----------------------------------------------------------------------
## Should only be accessed By PLNfamily but R6 friend class don't exist

#' Positions in the (Euclidian) parameter space, noted as Z in the model. Used to compute the likelihood.
#'
#' @name PLNfit_latent_pos
#'
#' @param covariates a matrix of covariates. Will usually be extracted from the corresponding field in PLNfamily-class
#' @param offsets    a matrix of offsets. Will usually be extracted from the corresponding field in PLNfamily-class
#'
PLNfit$set("public", "latent_pos",
function(covariates, offsets) {
  latentPos <- private$M + tcrossprod(covariates, private$Theta) + offsets
  latentPos
})

## ----------------------------------------------------------------------
## PUBLIC METHODS FOR THE USERS
## ----------------------------------------------------------------------
## For each methods, I define an corresponding S3 method for simplicity
## and only document the S3 method

#' Predict counts of a new sample
#'
#' @name predict.PLNfit
#'
#' @param x an R6 object with class PLNfit
#' @param newdata    A data frame in which to look for variables with which to predict.
#' @param newOffsets A matrix in which to look for offsets with which to predict.
#' @param type       The type of prediction required. The default is on the scale of the linear predictors (i.e. log average count);
#'                   the alternative "response" is on the scale of the response variable (i.e. average count)
#' @return A matrix of predicted log-counts (if type = "link") or predicted counts (if type = "response").
#' @export
predict.PLNfit <- function(x, newdata, newOffsets, type = c("link", "response")) {
  stopifnot(isPLNfit(x))
  x$predict(newdata, newOffsets, type)
}

PLNfit$set("public", "predict",
  function(newdata, newOffsets, type = c("link", "response")) {
    type = match.arg(type)
    ## Are matrix conformable?
    stopifnot(ncol(newdata)    == ncol(private$Theta),
              nrow(newdata)    == nrow(newOffsets),
              ncol(newOffsets) == nrow(private$Theta))
    ## Mean latent positions in the parameter space
    EZ <- tcrossprod(newdata, private$Theta) + newOffsets
    results <- switch(type,
                      link     = EZ,
                      response = exp(EZ))
    ## output formatting
    rownames(results) <- rownames(newdata); colnames(results) <- rownames(private$Theta)
    attr(results, "type") <- type
    results
  }
)

#' Display the model parameters of a PLNfit in a matrix fashion
#'
#' @name plot.PLNfit
#'
#' @param x an R6 object with class PLNfit
#' @param type character. Should the variational or the model parameters be plotted? default is "model".
#'
#' @export
plot.PLNfit <- function(x, type=c("model","variational")) {
  stopifnot(isPLNfit(x))
  x$plot(type)
}

PLNfit$set("public", "plot",
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

PLNfit$set("public", "show",
function(model = "A Poisson Lognormal fit\n") {
  cat(model)
  cat("==================================================================\n")
  print(as.data.frame(t(self$criteria), row.names = ""))
  cat("==================================================================\n")
  cat("* Useful fields \n")
  cat("  $model_par, $var_par, $loglik, $BIC, $ICL, $degrees_freedom, $criteria\n")
  cat("* Useful methods\n")
  cat("    $plot(), $predict()\n")
})
PLNfit$set("public", "print", function() self$show())
print.PLNfit <- function(x) {x$show()}
