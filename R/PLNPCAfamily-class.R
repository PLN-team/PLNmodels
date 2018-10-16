#' An R6 Class to represent a collection of PLNPCAfit
#'
#' @description The function \code{\link{PLNPCA}} produces an instance of this class.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for \code{\link[=PLNfamily_getBestModel]{getBestModel}},
#' \code{\link[=PLNfamily_getModel]{getModel}}, \code{\link[=plot.PLNfamily]{plot}}
#' and \code{\link[=predict.PLNfit]{predict}}.
#'
#' @field responses the matrix of responses common to every models
#' @field covariates the matrix of covariates common to every models
#' @field offsets the matrix of offsets common to every models
#' @field ranks the dimensions of the successively fitted models
#' @field models a list of \code{\link[=PLNPCAfit]{PLNPCAfit}} object, one per rank.
#' @field inception a \code{\link[=PLNfit]{PLNfit}} object, obtained when full rank is considered.
#' @field criteria a data frame with the value of some criteria (variational lower bound J, BIC, ICL and R2) for the different models.
#' @include PLNfamily-class.R
#' @importFrom R6 R6Class
#' @import ggplot2
#' @seealso The function \code{\link{PLNPCA}}, the class \code{\link[=PLNPCAfit]{PLNPCAfit}}
PLNPCAfamily <-
  R6Class(classname = "PLNPCAfamily",
    inherit = PLNfamily,
     active = list(
      ranks = function() private$params
    )
)

PLNPCAfamily$set("public", "initialize",
  function(ranks, responses, covariates, offsets, control) {

  ## initialize the required fields
  super$initialize(responses, covariates, offsets, control)
  private$params <- ranks

  if (control$trace > 0) cat("\n Perform SVD to initialize other parameters...")
  svdM     <- svd(self$inception$var_par$M, nu=max(ranks), nv=max(ranks))
  svdSigma <- svd(self$inception$model_par$Sigma, nu=max(ranks), nv=0)
  if (control$covariance == "full")
    svdS     <- svd(self$inception$var_par$S, nu=max(ranks), nv=max(ranks))


  ## instantiate as many models as ranks
  self$models <- lapply(ranks, function(q){
    if (control$covariance == "full") {
      S0 <- svdS$u[, 1:q, drop=FALSE] %*% diag(svdS$d[1:q], nrow=q, ncol=q) %*% t(svdS$v[1:q, 1:q, drop=FALSE])
    }
    else {
      S0 <- matrix(mean(self$inception$var_par$S), private$n, q)
    }
    model <- PLNPCAfit$new(
      Theta = self$inception$model_par$Theta,
      Sigma = svdSigma$u[, 1:q, drop=FALSE] %*% diag(svdSigma$d[1:q],nrow=q, ncol=q) %*% t(svdSigma$u[, 1:q, drop=FALSE]),
      B = svdSigma$u[, 1:q, drop=FALSE] %*% sqrt(diag(svdSigma$d[1:q],nrow=q, ncol=q)),
      M = svdM$u[, 1:q, drop=FALSE] %*% diag(svdM$d[1:q], nrow=q, ncol=q) %*% t(svdM$v[1:q, 1:q, drop=FALSE]),
      S = S0
    )
    return(model)
  })

})

PLNPCAfamily$set("public", "optimize",
  function(control) {

  self$models <- mclapply(self$models, function(model) {
    ## initial parameters (model + variational)
    par0 <- c(
      model$model_par$Theta  ,
      model$model_par$B,
      model$var_par$M,
      pmax(model$var_par$S,10*control$lower_bound)
    )

    ## ===========================================
    ## OPTIMISATION
    if (control$trace == 1) {
      cat("\t Rank approximation =",model$rank, "\r")
      flush.console()
    }

    if (control$trace > 1) {
      cat(" Rank approximation =",model$rank)
      cat("\n\t conservative convex separable approximation for gradient descent")
    }

    ## CALL TO NLOPT OPTIMIZATION WITH BOX CONSTRAINT
    opts <- control
    opts$xtol_abs <- c(rep(0, private$p*(private$d + model$rank) + private$n * model$rank),
                        rep(control$xtol_abs, private$n*model$rank))
    opts$lower_bound <- c(rep(-Inf, private$p*private$d), # Theta
                        rep(-Inf, private$p*model$rank) , # B
                        rep(-Inf, private$n*model$rank) , # M
                        rep(control$lower_bound,private$n*model$rank)) # S
    optim_out <- optimization_PLNPCA(par0, self$responses, self$covariates, self$offsets, model$rank, opts)

    ## ===========================================
    ## OUTPUT

    rownames(optim_out$Theta) <- rownames(optim_out$B) <- rownames(optim_out$Sigma) <- colnames(optim_out$Sigma) <- colnames(self$responses)
    colnames(optim_out$Theta) <- colnames(self$covariates)
    colnames(optim_out$B) <- colnames(optim_out$M) <- 1:model$rank
    rownames(optim_out$M) <- rownames(self$responses)

    model$update(
      B     = optim_out$B,
      Theta = optim_out$Theta,
      Sigma = optim_out$Sigma,
      M     = optim_out$M,
      S     = optim_out$S,
      J     = sum(optim_out$loglik),
      monitoring = list(
        objective  = -sum(optim_out$loglik),
        iterations = optim_out$iterations,
        status     = optim_out$status,
        message    = statusToMessage(optim_out$status))
    )
    model
  }, mc.cores = control$cores, mc.allow.recursive = FALSE)
})

PLNPCAfamily$set("public", "postTreatment",
function() {
  super$postTreatment()
  for (model in self$models) {
    model$setVisualization()
  }
})

PLNPCAfamily$set("public", "plot",
function(criteria = c("loglik", "BIC", "ICL")) {
  vlines <- sapply(criteria, function(crit) self$getBestModel(crit)$rank)
  p <- super$plot(criteria) + xlab("rank") + geom_vline(xintercept = vlines, linetype = "dashed", alpha = 0.25)
  p
})

PLNPCAfamily$set("public", "show",
function() {
  super$show()
  cat(" Task: Principal Component Analysis\n")
  cat("========================================================\n")
  cat(" - Ranks considered: from", min(self$ranks), "to", max(self$ranks),"\n")
  cat(" - Best model (regarding BIC): rank =", self$getBestModel("BIC")$rank, "- R2 =", round(self$getBestModel("BIC")$R_squared, 2), "\n")
  cat(" - Best model (regarding ICL): rank =", self$getBestModel("ICL")$rank, "- R2 =", round(self$getBestModel("ICL")$R_squared, 2), "\n")
})
PLNPCAfamily$set("public", "print", function() self$show())
