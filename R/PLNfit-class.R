#' @importFrom R6 R6Class
PLNfit <-
   R6Class(classname = "PLNfit",
    public = list(
      model.par       = NULL, # Theta and Sigma
      variational.par = NULL, # M and S
      criteria        = NULL, # J, BIC, ICL
      convergence     = NULL, # results of optimization
      loglik          = NULL,
      initialize = function(model.par=NA, variational.par=NA, criteria=NA, convergence=NA, loglik=NA) {
      self$model.par       <- model.par
      self$variational.par <- variational.par
      self$criteria        <- criteria
      self$convergence     <- convergence
      }
    )
  )
## TODO add a predict method and similar stat tools
