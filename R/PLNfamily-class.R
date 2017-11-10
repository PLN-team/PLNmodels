## COMMON TO EVERY COLLECTIONS OF MODELS BASED ON PLN (PLNPCA, PLNnetwork)
PLNfamily <-
  R6Class(classname = "PLNfamily",
    public = list(
      responses  = NULL, # the Y matrix
      covariates = NULL, # the X matrix
      offsets    = NULL, # the O matrix
      inception  = NULL, # the basic model in the collection (no regularization, nor sparsity, nor rank)
      models     = NULL, # the collection of models to be fitted
      params     = NULL  # vector of the parameters that indexes the models (either sparsity, rank, etc.)
    ),
    private = list(
      fn_optim   = NULL, # objective and gradient for optimizing the regularized models
      n          = NULL, # number of samples
      p          = NULL, # number of responses
      d          = NULL  # number of covariates
    ),
    active = list(
      # send back a data frame with some criteria associated with the collection of fits
      criteria = function() {
        data.frame(param = self$params,
                  t(sapply(self$models, function(model) model$criteria)))
      }
    )
)

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
    if (control$trace > 0) cat("\n Adjust a PLN to define the inceptive model")
    if (isTRUE(all.equal(is.character(control$inception), control$inception == "PLN")) ) {
      self$inception <- PLN(self$responses, self$covariates, self$offsets, control)
    } else {
      par0 <- initializePLN(self$responses, self$covariates, self$offsets, control)
      Sigma <- crossprod(par0$M)/private$n + diag(colMeans(par0$S), private$p, private$p)
      self$inception <- PLNfit$new(Theta = par0$Theta, Sigma = Sigma, M = par0$M, S = par0$S)
    }
})

#' Best model extraction from a collection of PLNfit (PCA, network)
#'
#' @name PLNfamily_getBestModel
#'
#' @param crit a character for the criterion used to performed the selection. Either
#' "BIC", "ICL", "J", "R2". Default is "BIC".
#' @return  Send back a object with class \code{\link[=PLNfit]{PLNfit}}.
NULL
PLNfamily$set("public", "getBestModel",
function(crit = c("BIC", "ICL", "J", "R2")){
  crit <- match.arg(crit)
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
#' @param var value of the parameter (rank for PCA, penalty for network) that identifies the model to be extracted from the collection. Must belong to the vector params
#' @return Send back a object with class \code{\link[=PLNfit]{PLNfit}}.
NULL
PLNfamily$set("public", "getModel",
function(var){
  id <- match(var, self$params)
  if (!is.na(id)) {
    return(self$models[[id]]$clone())
  } else {
    stop("No such a model in the collection.")
  }
})

#' @import ggplot2
PLNfamily$set("public", "plot",
function() {
  dplot <- melt(self$criteria[, c("param", "J", "BIC", "ICL")], id.vars = 1, variable.name = "criterion")
  p <- ggplot(dplot, aes(x = param, y = value, group = criterion, colour = criterion)) +
        geom_line() + geom_point() + ggtitle("Model selection criteria")
  return(p)
})

PLNfamily$set("public", "show",
function() {
  cat("--------------------------------------------------------\n")
  cat("COLLECTION OF", length(self$models), "POISSON LOGNORMAL MODELS\n")
  cat("--------------------------------------------------------\n")
})

PLNfamily$set("public", "print", function() self$show())

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

# Set data frame with criteria (BIC, ICL, J and possibly R2) associated with the collection of fits
PLNfamily$set("public", "postTreatment",
function() {
  private$computeR2()
})

## ----------------------------------------------------------------
## PRIVATE METHODS

# Set data frame with criteria (BIC, ICL, J and possibly R2) associated with the collection of fits
PLNfamily$set("private", "setCriteria",
function() {
  self$criteria <- data.frame(param = self$params,
                              t(sapply(self$models, function(model) model$criteria)))
})

# Compute goodness of fit (R2) for each fit in the family
PLNfamily$set("private", "computeR2",
function() {
  ## Likelihoods of the null and saturated models
  lmin <- logLikPoisson(self$responses, nullModelPoisson(self$responses, self$covariates, self$offsets))
  lmax <- logLikPoisson(self$responses, fullModelPoisson(self$responses))
  lapply(self$models, function(model) {
    loglik <- logLikPoisson(self$responses, model$latentPos(self$covariates, self$offsets))
    model$update(R2 = (loglik - lmin) / (lmax - lmin))
  })
})
