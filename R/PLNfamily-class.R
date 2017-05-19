## COMMON TO PLNPCA and PLNnetwork

PLNfamily <-
  R6Class(classname = "PLNfamily",
    public = list(
      type       = NULL,
      models     = NULL,
      init.par   = NULL,
      criteria   = NULL,
      responses  = NULL,
      covariates = NULL,
      offsets    = NULL,
      objective  = NULL,
      gradient   = NULL
    )
)

PLNfamily$set("public", "initialize",
  function(type, responses, covariates, offsets) {
    self$type       <- type
    self$responses  <- responses
    self$covariates <- covariates
    self$offsets    <- offsets

    ## set names of the data matrices
    if (is.null(rownames(responses)))  rownames(self$responses)  <- 1:nrow(responses)
    if (is.null(colnames(responses)))  colnames(self$responses)  <- 1:ncol(responses)
    if (is.null(rownames(covariates))) rownames(self$covariates) <- 1:nrow(covariates)
    if (is.null(colnames(covariates))) colnames(self$covariates) <- 1:ncol(covariates)

    ## recover the initial model for each rank with glm Poisson models
    glmP  <- lapply(1:ncol(responses), function(j) glm.fit(covariates, responses[, j], offset = offsets[,j], family = poisson()))
    Theta <- matrix(Reduce(rbind, lapply(glmP, coefficients)), ncol=ncol(covariates))
    Sigma <- cov(sapply(glmP, residuals.glm, type="pearson"))
    self$init.par <- list(Sigma = Sigma, Theta=Theta)
})
