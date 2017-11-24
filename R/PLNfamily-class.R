## COMMON TO EVERY COLLECTIONS OF MODELS BASED ON PLN (PLNPCA, PLNnetwork)
PLNfamily <-
  R6Class(classname = "PLNfamily",
    public = list(
      responses  = NULL, # the Y matrix
      covariates = NULL, # the X matrix
      offsets    = NULL, # the O matrix
      inception  = NULL, # the basic model in the collection (no regularization, nor sparsity, nor rank)
      models     = NULL  # the collection of fitted models
    ),
    private = list(
      params     = NULL,  # vector of parameters that indexes the models (either sparsity, rank, etc.)
      fn_optim   = NULL, # objective and gradient for optimizing the regularized models
      n          = NULL, # number of samples
      p          = NULL, # number of responses
      d          = NULL  # number of covariates
    ),
    active = list(
      # send back a data frame with some criteria associated with the collection of fits
      criteria = function() {
        res <- do.call(rbind, lapply(self$models, function(model) {model$criteria}))
        data.frame(param = private$params, res)
      },
      # send back a data frame with some quantities associated with the optimization process
      convergence = function() {
        res <- do.call(rbind, lapply(self$models, function(model) {
          c(degrees_freedom = model$degrees_freedom, sapply(model$optim_par, function(x) x[length(x)]))
        }))
        data.frame(param = private$params, res, stringsAsFactors = FALSE)
      }
    )
)

## an S3 function to check if an object is a PLNfit
isPLNfamily <- function(Robject) {all.equal(rev(class(Robject))[1:2], c('R6','PLNfamily'))}

PLNfamily$set("public", "initialize",
  function(responses, covariates, offsets, control) {

    ## set data matrice and dimension
    self$responses  <- responses
    self$covariates <- covariates
    self$offsets    <- offsets
    private$n <- nrow(responses)
    private$p <- ncol(responses)
    private$d <- ncol(covariates)

    ## set names of the data matrices
    if (is.null(rownames(responses)))  rownames(self$responses)  <- 1:private$n
    if (is.null(colnames(responses)))  colnames(self$responses)  <- 1:private$p
    if (is.null(rownames(covariates))) rownames(self$covariates) <- 1:private$n
    if (is.null(colnames(covariates))) colnames(self$covariates) <- 1:private$d

    ## extract the model used for initializaing the whole family
    if (control$trace > 0) cat("\n Adjust the inceptive model")
    if (isTRUE(all.equal(is.character(control$inception), control$inception == "PLN")) ) {
      self$inception <- PLN(self$responses, self$covariates, self$offsets, control)
    } else {
      par0 <- initializePLN(self$responses, self$covariates, self$offsets, control)
      Sigma <- crossprod(par0$M)/private$n + diag(colMeans(par0$S), private$p, private$p)
      self$inception <- PLNfit$new(Theta = par0$Theta, Sigma = Sigma, M = par0$M, S = par0$S)
    }
})

## a method to compute and set fields after optimization
PLNfamily$set("public", "postTreatment",
function() {
  ## Compute R2
  ## Likelihoods of the null and saturated models
  lmin <- logLikPoisson(self$responses, nullModelPoisson(self$responses, self$covariates, self$offsets))
  lmax <- logLikPoisson(self$responses, fullModelPoisson(self$responses))
  for (model in self$models) {
    loglik <- logLikPoisson(self$responses, model$latent_pos(self$covariates, self$offsets))
    model$update(R2 = (loglik - lmin) / (lmax - lmin))
  }
})

#' Best model extraction from a collection of PLNfit (PCA, network)
#'
#' @name PLNfamily_getBestModel
#'
#' @param crit a character for the criterion used to performed the selection. Either
#' "BIC", "EBIC", "ICL", "loglik", "R_squared". Default is "BIC".
#' @return  Send back a object with class \code{\link[=PLNfit]{PLNfit}}.
NULL
PLNfamily$set("public", "getBestModel",
function(crit = c("BIC", "ICL", "EBIC", "loglik", "R_squared")){
  crit <- match.arg(crit)
  stopifnot(!anyNA(self$criteria[[crit]]))
  if (length(self$criteria[[crit]]) > 1) {
    id <- which.max(self$criteria[[crit]])
  } else {id <- 1}
    model <- self$models[[id]]$clone()
    return(model)
})

#' Model extraction from a collection of PLN models
#'
#' @name PLNfamily_getModel
#'
#' @param var value of the parameter (rank for PCA, penalty for network) that identifies the model to be extracted from the collection.
#' @return Send back a object with class \code{\link[=PLNfit]{PLNfit}}.
NULL
PLNfamily$set("public", "getModel",
function(var){
  id <- match(var, private$params)
  if (!is.na(id)) {
    return(self$models[[id]]$clone())
  } else {
    stop("No such a model in the collection.")
  }
})

#' Predict counts of new samples for all fits in the family
#'
#' @name PLNfamily_predict
#'
#' @param newdata    A optional data frame in which to look for variables with which to predict. If omitted, the family-level covariates are used.
#' @param newOffsets A optional matrix in which to look for offsets with which to predict. If omitted, the family-level offsets are used.
#' @param type       The type of prediction required. The default is on the scale of the linear predictors (i.e. log average count);
#'                   the alternative "response" is on the scale of the response variable (i.e. average count)
#' @return A list of matrices (one per fit) of predicted log-counts (if type = "link")
#'         or predicted counts (if type = "response").
#'
PLNfamily$set("public", "predict",
function(newdata = self$covariates, newOffsets = self$offsets, type = c("link", "response")) {
  lapply(self$models, function(model) {
    model$predict(newdata, newOffsets, type)
  })
})


#' Display the criteria associated with a collection of PLN fits (a PLNfamily)
#'
#' @name plot.PLNfamily
#'
#' @param x an R6 object with class PLNfamily
#' @param type vector of characters. The criteria to plot in c("loglik", "BIC", "ICL", "R_squared", "EBIC", "pen_loglik")
#' The two last are only available por PLNnetworkfamily. Default is c("loglik", "BIC", "ICL")
#'@return Produces a plot  representing the evolution of the criteria of the different models considered,
#' highlighting the best model in terms of ICL for PLNPCA and EBIC for PLNnetwork.
#'
#' @import ggplot2
#' @export
plot.PLNfamily <- function(x, criteria = c("loglik", "BIC", "ICL")) {
  stopifnot(isPLNfamily(x))
  x$plot(criteria)
}

PLNfamily$set("public", "plot",
function(criteria = c("loglik", "BIC", "ICL")) {
  dplot <- melt(self$criteria[, c('param',criteria)], id.vars = 1, variable.name = "criterion")
  p <- ggplot(dplot, aes(x = param, y = value, group = criterion, colour = criterion)) +
        geom_line() + geom_point() + ggtitle("Model selection criteria") + theme_bw()
  p
})

PLNfamily$set("public", "show",
function() {
  cat("--------------------------------------------------------\n")
  cat("COLLECTION OF", length(self$models), "POISSON LOGNORMAL MODELS\n")
  cat("--------------------------------------------------------\n")
})
PLNfamily$set("public", "print", function() self$show())
