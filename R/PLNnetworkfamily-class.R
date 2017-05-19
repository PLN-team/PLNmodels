#' @include PLNfamily-class.R

PLNnetworkfamily <-
  R6Class(classname = "PLNnetworkfamily",
    inherit = PLNfamily,
     public = list(
      ranks = "numeric"
    )
)

PLNnetworkfamily$set("public", "initialize",
  function(type, penalties, responses, covariates, offsets) {

  ## initialize the required fields
  super$initialize(type, responses, covariates, offsets)
  self$penalties <- penalties

  # fit <- PLNnetworkfit$new(type = self$type)

})
