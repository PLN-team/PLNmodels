## COMMON TO EVERY COLLECTIONS OF MODELS BASED ON PLN (PLNPCA, PLNnetwork)
PLNfamily <-
  R6Class(classname = "PLNfamily",
    public = list(
      responses  = NULL, # the Y matrix
      covariates = NULL, # the X matrix
      offsets    = NULL, # the O matrix
      weights    = NULL, # the vector of obervation weights
      inception  = NULL, # the basic model in the collection (no regularization, nor sparsity, nor rank)
      models     = NULL  # the collection of fitted models
    ),
    private = list(
      params     = NULL, # vector of parameters that indexes the models (either sparsity, rank, etc.)
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
          c(nb_param = model$nb_param, sapply(model$optim_par, function(x) x[length(x)]))
        }))
        data.frame(param = private$params, res, stringsAsFactors = FALSE)
      }
    )
)

PLNfamily$set("public", "initialize",
  function(responses, covariates, offsets, weights, control) {

    ## set data matrice and dimension
    self$responses  <- responses
    self$covariates <- covariates
    self$offsets    <- offsets
    self$weights    <- weights
    private$n <- nrow(responses)
    private$p <- ncol(responses)
    private$d <- ncol(covariates)

    ## set names of the data matrices
    if (is.null(rownames(responses)))  rownames(self$responses)  <- 1:private$n
    if (is.null(colnames(responses)))  colnames(self$responses)  <- 1:private$p
    if (is.null(rownames(covariates))) rownames(self$covariates) <- 1:private$n
    if (is.null(colnames(covariates)) & private$d > 0) colnames(self$covariates) <- 1:private$d

})

## a method to compute and set fields after optimization
PLNfamily$set("public", "postTreatment",
function() {
  nullModel <- nullModelPoisson(self$responses, self$covariates, self$offsets, self$weights)
  for (model in self$models)
    model$postTreatment(
      self$responses,
      self$covariates,
      self$offsets,
      self$weights,
      nullModel = nullModel
    )
})

PLNfamily$set("public", "getModel",
function(var, index = NULL) {
  ## Extraction by index
  if (!is.null(index) && index <= length(self$models)) {
    return(self$models[[index[1]]]$clone())
  }
  ## Extraction by parameter value
  id <- match(var, private$params)
  if (!is.na(id)) { ## Exact match
    return(self$models[[id]]$clone())
  } else { ## No exact match
    id <- which.min(abs(var - private$params)) ## closest model (in terms of parameter value)
    warning(paste("No such a model in the collection. Acceptable parameter values can be found via",
               "$ranks() (for PCA)",
               "$penalties() (for network)",
               paste("Returning model with closest value. Requested:", var, ", returned:", private$params[id]),
               sep = "\n"))
    return(self$models[[id]]$clone())
  }
})

PLNfamily$set("public", "plot",
function(criteria, annotate = TRUE) {
  stopifnot(!anyNA(self$criteria[criteria]))
  dplot <- self$criteria %>%
    dplyr::select(c("param", criteria)) %>%
    tidyr::gather(key = "criterion", value = "value", -param) %>%
    dplyr::group_by(criterion)
  p <- ggplot(dplot, aes(x = param, y = value, group = criterion, colour = criterion)) +
       geom_line() + geom_point() + ggtitle("Model selection criteria") + theme_bw()

  if (annotate)
    p <- p + annotate("text", x = self$criteria$param, y = min(dplot$value), hjust=-.1, angle = 90,
           label = paste("R2 =", round(self$criteria$R_squared, 2)), size = 3, alpha = 0.7)
  p
})

PLNfamily$set("public", "show",
function() {
  cat("--------------------------------------------------------\n")
  cat("COLLECTION OF", length(self$models), "POISSON LOGNORMAL MODELS\n")
  cat("--------------------------------------------------------\n")
})

PLNfamily$set("public", "print", function() self$show())

