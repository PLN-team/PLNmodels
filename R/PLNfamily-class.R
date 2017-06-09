## COMMON TO PLNPCA and PLNnetwork

PLNfamily <-
  R6Class(classname = "PLNfamily",
    public = list(
      models     = NULL,
      init.par   = NULL,
      criteria   = NULL,
      responses  = NULL,
      covariates = NULL,
      offsets    = NULL,
      n          = NULL, # number of samples
      p          = NULL, # number of responses
      d          = NULL, # number of covariates
      objective  = NULL,
      gradient   = NULL,
      fn_optim   = NULL
    )
)

PLNfamily$set("public", "initialize",
  function(responses, covariates, offsets) {
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

    ## recover the initial model for each rank with glm Poisson models
    glmP  <- lapply(1:self$p, function(j) glm.fit(covariates, responses[, j], offset = offsets[,j], family = poisson()))
    Theta <- do.call(rbind, lapply(glmP, coefficients))
    # Y.hat <- exp(offsets + tcrossprod(covariates, Theta))
    # Sigma <- diag(log(1 + colMeans(((responses - Y.hat)^2 - responses)/(Y.hat^2))))
    Sigma <- cov(sapply(glmP, residuals.glm, "pearson"))

    self$init.par <- list(Sigma = Sigma, Theta = Theta)
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
  if(length(self$criteria$BIC) >1) {
    id <- switch(crit,
    "BIC" = which.max(self$criteria$BIC),
    "ICL" = which.max(self$criteria$ICL),
    "J"   = which.max(self$criteria$J),
    "R2"  = which.max(self$criteria$R2))
  } else {id <- 1}
    model <- self$models[[id]]$clone()
    return(model)
})

#' Model extraction from a collection of PLNPCAfit
#'
#' @name PLNfamily_getModel
#'
#' @param xvar a number (rank for PCA, penalty for network) identifying the model to be extracted from the collection.
#' @return Send back a object with class \code{\link[=PLNfit]{PLNfit}}.
NULL
PLNfamily$set("public", "getModel",
function(xvar){
  id <- match(round(xvar,12), as.numeric(names(self$models)))
  if (!is.na(id)){
    return(self$models[[id]]$clone())
  } else {
    stop("No such a model in the collection.")
  }
})

PLNfamily$set("public", "setCriteria",
function() {
  self$criteria <- data.frame(xvar = round(as.numeric(names(self$models)), 12),
                              t(sapply(self$models, function(model) model$criteria)))
})

PLNfamily$set("public", "plot",
function() {
  dplot <- melt(self$criteria[, c("xvar", "J", "BIC", "ICL")], id.vars = 1, variable.name = "criterion")
  p <- ggplot(dplot, aes(x=xvar, y=value, group=criterion, colour=criterion)) +
        geom_line() + geom_point() + ggtitle("Model selection criteria")
  return(p)
})

PLNfamily$set("public", "show",
function() {
  cat("POISSON LOGNORMAL MODEL\n")
  cat("------------------------------------------------------\n")
})

PLNfamily$set("public", "print", function() self$show())
