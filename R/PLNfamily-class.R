## COMMON TO EVERY COLLECTIONS OF MODELS BASED ON PLN (PLNPCA, PLNnetwork)
PLNfamily <-
  R6Class(classname = "PLNfamily",
    public = list(
      responses  = NULL, # the Y matrix
      covariates = NULL, # the X matrix
      offsets    = NULL, # the O matrix
### TODO pass n, p, d to private members
      n          = NULL, # number of samples
      p          = NULL, # number of responses
      d          = NULL, # number of covariates
      inception  = NULL, # the basic PLN model in the collection (no regularization, nor sparsity, nor rank)
      models     = NULL, # the collection of models to be fitted
      criteria   = NULL, # a data frame with some criteria associated with the collection of fits
      fn_optim   = NULL  # objective and gradient for optimizing the regularized models
    )
)

PLNfamily$set("public", "initialize",
  function(responses, covariates, offsets, control) {
    ## set data matrice and dimension
    self$responses  <- responses
    self$covariates <- covariates
    self$offsets    <- offsets
    self$n <- nrow(responses)
    self$p <- ncol(responses)
    self$d <- ncol(covariates)

    ## set names of the data matrices
    if (is.null(rownames(responses)))  rownames(self$responses)  <- 1:self$n
    if (is.null(colnames(responses)))  colnames(self$responses)  <- 1:self$p
    if (is.null(rownames(covariates))) rownames(self$covariates) <- 1:self$n
    if (is.null(colnames(covariates))) colnames(self$covariates) <- 1:self$d

    ## adjust the basic PLN model
    if (!is.null(control$inception)) {
      stopifnot(all.equal(class(control$inception), c("PLNfit", "R6")),
                all.equal(dim(control$inception$model.par$Sigma), c(self$p,self$p)),
                all.equal(dim(control$inception$model.par$Theta), c(self$d,self$p)),
                all.equal(dim(control$inception$variational.par$M), c(self$n,self$p)),
                all.equal(dim(control$inception$variational.par$S), c(self$n,self$p)))
      cat("\n User defined inceptive PLN model")
      self$inception <- control$inception
    } else {
      self$inception <- PLN(responses, covariates, offsets, control)
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
function(crit=c("BIC", "ICL", "J", "R2")){
  crit <- match.arg(crit)
  if(length(self$criteria[[crit]]) > 1) {
    id <- which.max(self$criteria[[crit]])
  } else {id <- 1}
    model <- self$models[[id]]$clone()
    return(model)
})

#' Model extraction from a collection of PLN models
#'
#' @name PLNfamily_getModel
#'
#' @param xvar a number (rank for PCA, penalty for network) identifying the model to be extracted from the collection.
#' @return Send back a object with class \code{\link[=PLNfit]{PLNfit}}.
NULL
PLNfamily$set("public", "getModel",
function(xvar){
  id <- match(round(xvar,16), as.numeric(names(self$models)))
  if (!is.na(id)){
    return(self$models[[id]]$clone())
  } else {
    stop("No such a model in the collection.")
  }
})

#' Set data frame with criteria (BIC, ICL, J and possibly R2) associated with the collection of fits
#'
#' @name PLNfamily_setCriteria
#'
PLNfamily$set("public", "setCriteria",
function() {
  self$criteria <- data.frame(xvar = round(as.numeric(names(self$models)), 16),
                              t(sapply(self$models, function(model) model$criteria)))
})

#' @import ggplot2
PLNfamily$set("public", "plot",
function() {
  dplot <- melt(self$criteria[, c("xvar", "J", "BIC", "ICL")], id.vars = 1, variable.name = "criterion")
  p <- ggplot(dplot, aes(x=xvar, y=value, group=criterion, colour=criterion)) +
        geom_line() + geom_point() + ggtitle("Model selection criteria")
  return(p)
})

#' Basic show method
#'
#' @name PLNfamily_show
#'
#' @param verbose a logical (defaults to TRUE) controlling the amount of screen output.
PLNfamily$set("public", "show",
function(verbose = TRUE) {
  cat("COLLECTION OF", length(self$models), "POISSON LOGNORMAL MODELS\n")
  cat("------------------------------------------------------\n")
  if (verbose) {
    cat(" - Available models are:\n")
    cat(paste("    +", names(self$models)), sep = "\n")
  }
})

PLNfamily$set("public", "print", function() self$show())

#' Compute goodness of fit (R2) for each fit in the family
#'
#' @name PLNfamily_computeR2
#'
PLNfamily$set("public", "computeR2",
function() {
  ## Likelihood of the null model (covariates and offset only, no latent variables)
  ## computed from the inception model (should we move to a Poisson-GLM as the base model instead?)
  # Theta <- do.call(rbind, lapply(1:p, function(j)
  #   coefficients(glm.fit(self$covariates, self$responses[, j], offset = self$offsets[,j], family = poisson()))))
  Z <- self$offsets + tcrossprod(self$covariates, self$inception$model.par$Theta)
  lmin <- sum(self$responses * Z) - sum(as.numeric(self$responses))
  ## Likelihood of the full model
  lmax <- sum(self$responses * (log(self$responses + 1*(self$responses == 0)) - 1))
  ## Likelihood of each model
  for (model in self$models) {
    Z <- model$latentPos(self$covariates, self$offsets)
    loglik <- sum(self$responses * (Z)) - sum(as.numeric(self$responses))
    R2 <- (loglik - lmin) / (lmax - lmin)
    model$addCriteria("R2", R2)
  }
})
