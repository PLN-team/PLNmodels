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
    self$inception <- PLN(responses, covariates, offsets, control)
})

#' Best model extraction from a collection of PLNfit (PCA, network)
#'
#' @name PLNfamily_getBestModel
#'
#' @param crit a character for the criterion used to performed the selection. Either
#' "ICL", "BIC", "J" or "R2". Default is "ICL.
#' @return  Send back a object with class \code{\link[=PLNfit]{PLNfit}}.
NULL
PLNfamily$set("public", "getBestModel",
function(crit=c("ICL", "BIC", "J", "R2")){
  crit <- match.arg(crit)
  if(length(self$criteria$BIC) > 1) {
    id <- switch(crit,
    "BIC" = which.max(self$criteria$BIC),
    "ICL" = which.max(self$criteria$ICL),
    "J"   = which.max(self$criteria$J),
    "R2"  = which.max(self$criteria$R2))
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

PLNfamily$set("public", "setCriteria",
function() {
  self$criteria <- data.frame(xvar = round(as.numeric(names(self$models)), 16),
                              t(sapply(self$models, function(model) model$criteria)))
})

PLNfamily$set("public", "plot",
function() {
  dplot <- melt(self$criteria[, c("xvar", "J", "BIC", "ICL")], id.vars = 1, variable.name = "criterion")
  p <- ggplot2::ggplot(dplot, aes(x=xvar, y=value, group=criterion, colour=criterion)) +
        geom_line() + geom_point() + ggtitle("Model selection criteria")
  return(p)
})

PLNfamily$set("public", "show",
function() {
  cat("COLLECTIONS OF POISSON LOGNORMAL MODELS\n")
  cat("------------------------------------------------------\n")
})

PLNfamily$set("public", "print", function() self$show())
