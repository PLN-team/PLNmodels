#' @importFrom R6 R6Class
PLNfit <-
   R6Class(classname = "PLNfit",
    public = list(
      type            = NULL,
      model.par       = NULL,
      variational.par = NULL,
      criteria        = NULL,
      convergence     = NULL,
      loglik          = NULL,
      initialize = function(type=NA, model.par=NA, variational.par=NA, criteria=NA, convergence=NA, loglik=NA) {
      self$type            <- type
      self$model.par       <- model.par
      self$variational.par <- variational.par
      self$criteria        <- criteria
      self$convergence     <- convergence
      self$loglik          <- loglik
      }
    )
  )
