#' An R6 Class to represent a collection of PLNnetworkfit
#'
#' @description The function \code{\link{PLNnetwork}} produces an instance of this class.
#'
#' This class comes with a set of methods, some of them being useful for the user:
#' See the documentation for \code{\link[=PLNfamily_getBestModel]{getBestModel}},
#' \code{\link[=PLNfamily_getModel]{getModel}}, \code{\link[=plot.PLNfamily]{plot}}
#' and \code{\link[=predict.PLNfit]{predict}}.
#'
#' @field responses the matrix of responses common to every models
#' @field covariates the matrix of covariates common to every models
#' @field offsets the matrix of offsets common to every models
#' @field penalties the sparsity level of the network in the successively fitted models
#' @field models a list of \code{\link[=PLNnetworkfit]{PLNnetworkfit}} object, one per penalty.
#' @field inception a \code{\link[=PLNfit]{PLNfit}} object, obtained when no sparsifying penalty is applied.
#' @field criteria a data frame with the value of some criteria (variational lower bound J, BIC, ICL and R2) for the different models.
#' @include PLNfamily-class.R
#' @importFrom R6 R6Class
#' @importFrom glasso glasso
#' @seealso The function \code{\link{PLNnetwork}}, the class \code{\link[=PLNnetworkfit]{PLNnetworkfit}}
PLNnetworkfamily <-
  R6Class(classname = "PLNnetworkfamily",
    inherit = PLNfamily,
    private = list(
      stab_path = NULL
      ), # a field to store the stability path,
    active = list(
      penalties = function() private$params,
      stability_path = function() private$stab_path,
      stability = function() {
        if (!is.null(private$stab_path)) {
          stability <- self$stability_path %>%
            dplyr::select(Penalty, Prob) %>%
            group_by(Penalty) %>%
            summarize(Stability = 1 - mean(4 * Prob * (1 - Prob))) %>%
            arrange(desc(Penalty)) %>%
            pull(Stability)
        } else {
          stability <- rep(NA, length(self$penalties))
        }
        stability
      },
      criteria = function() {mutate(super$criteria, stability = self$stability)}
    )
)

PLNnetworkfamily$set("public", "initialize",
  function(penalties, responses, covariates, offsets, control) {

    ## initialize fields shared by the super class
    super$initialize(responses, covariates, offsets, control)

    ## Get an appropriate grid of penalties
    if (is.null(penalties)) {
      if (control$trace > 1) cat("\n Recovering an appropriate grid of penalties.")
      max_pen <- max(abs(self$inception$model_par$Sigma))
      penalties <- 10^seq(log10(max_pen), log10(max_pen*control$min.ratio), len = control$nPenalties)
    } else {
      if (control$trace > 1) cat("\nPenalties already set by the user")
      stopifnot(all(penalties > 0))
    }

    ## instantiate as many models as penalties
    private$params <- sort(penalties, decreasing = TRUE)
    self$models <- lapply(private$params, function(penalty) {
      PLNnetworkfit$new(penalty = penalty)
    })

})

# One (coordinate-descent) optimization for each \code{\link[=PLNnetworkfit-class]{PLNnetworkfit}} model, using
# the same starting point (inception model from \code{\link[=PLNfit-class]{PLNfit}}) for all models.
#
PLNnetworkfamily$set("public", "optimize",
  function(control) {

  ## ===========================================
  ## INITIALISATION

  ## start from the standard PLN (a.k.a. inception)
  Sigma <- self$inception$model_par$Sigma
  par0  <- c(self$inception$model_par$Theta,
             self$inception$var_par$M,
             self$inception$var_par$S)
  Omega  <- diag(1/diag(Sigma), nrow = private$p, ncol = private$p)
  Sigma0 <- diag(  diag(Sigma), nrow = private$p, ncol = private$p)
  objective.old <- ifelse(is.na(self$inception$loglik),-Inf,-self$inception$loglik)

  ## set option for call to NLOPT: optim typ, lower bound, tolerance...
  opts <- list("algorithm"   = control$method,
               "maxeval"     = control$maxeval,
               "ftol_rel"    = control$ftol_rel,
               "ftol_abs"    = control$ftol_abs,
               "xtol_rel"    = control$xtol_rel,
               "xtol_abs"    = c(rep(0, private$p*private$d), # Theta
                                 rep(0, private$n*private$p), # M
                                 rep(control$xtol_abs, private$n*private$p)), #S
               "lower_bound" = c(rep(-Inf, private$p*private$d), # Theta
                                 rep(-Inf, private$n*private$p), # M
                                 rep(control$lbvar, private$n*private$p)) # S
               )

  ## ===========================================
  ## GO ALONG THE PENALTY GRID (i.e the models)
  for (m in seq_along(self$models))  {

    penalty <- self$models[[m]]$penalty

    ## ===========================================
    ## OPTIMISATION
    if (control$trace == 1) {
      cat("\tsparsifying penalty =",penalty, "\r")
      flush.console()
    }
    if (control$trace > 1) {
      cat("\tsparsifying penalty =",penalty, "- iteration:")
    }

    cond <- FALSE; iter <- 0
    objective   <- numeric(control$maxit_out)
    convergence <- numeric(control$maxit_out)
    while (!cond) {
      iter <- iter + 1
      if (control$trace > 1) cat("",iter)

      ## CALL TO GLASSO TO UPDATE Omega/Sigma
      glasso_out <- suppressWarnings(
                      glasso(Sigma,
                             rho = penalty,
                             penalize.diagonal = control$penalize.diagonal,
                             start = ifelse(control$warm, "warm", "cold"), w.init = Sigma0, wi.init = Omega
                             )
                      )
      Omega  <- glasso_out$wi ; if (!isSymmetric(Omega)) Omega <- Matrix::symmpart(Omega)
      Sigma0 <- glasso_out$w

      ## CALL TO NLOPT OPTIMIZATION WITH BOX CONSTRAINT
      logDetOmega <- as.double(determinant(Omega, logarithm = TRUE)$modulus)
      optim.out <- optimization_PLNnetwork(par0, self$responses, self$covariates, self$offsets, Omega, logDetOmega, opts)
      optim.out$message <- statusToMessage(optim.out$status)
      objective[iter]   <- optim.out$objective + penalty * sum(abs(Omega))
      convergence[iter] <- abs(objective[iter] - objective.old)/abs(objective[iter])

      ## Check convergence
      if ((convergence[iter] < control$ftol_out) | (iter >= control$maxit_out)) cond <- TRUE

      ## Post-Treatment to update Sigma
      M <- matrix(optim.out$solution[private$p*private$d               + 1:(private$n*private$p)], private$n,private$p)
      S <- matrix(optim.out$solution[(private$n + private$d)*private$p + 1:(private$n*private$p)], private$n,private$p)
      Sigma <- crossprod(M)/private$n + diag(colMeans(S), nrow = private$p, ncol = private$p)
      par0 <- optim.out$solution
      objective.old <- objective[iter]
    }

    ## ===========================================
    ## OUTPUT
    ## formating parameters for output
    Theta <- matrix(optim.out$solution[1:(private$p*private$d)], private$p, private$d)
    rownames(Theta) <- colnames(self$responses); colnames(Theta) <- colnames(self$covariates)
    dimnames(S)     <- dimnames(self$responses)
    dimnames(M)     <- dimnames(self$responses)
    rownames(Omega) <- colnames(Omega) <- colnames(self$responses)
    ## Optimization ends with a gradient descent step rather than a glasso step.
    ## Return Sigma from glasso step to ensure that Sigma = solve(Omega)
    ## Sigma <- Sigma0 ; if (!isSymmetric(Sigma)) Sigma <- Matrix::symmpart(Sigma)
    ## dimnames(Sigma) <- dimnames(Omega)

    self$models[[m]]$update(Omega = Omega, Sigma = Sigma, Theta = Theta, M = M, S = S, J = -optim.out$objective,
                            monitoring = list(objective = objective[1:iter],
                                              convergence = convergence[1:iter],
                                              outer_iterations = iter,
                                              inner_iterations = optim.out$iterations,
                                              inner_status = optim.out$status,
                                              inner_message = optim.out$message))
    if (control$trace > 1) {
      cat("\r                                                                                    \r")
      flush.console()
    }

  }

})

#' Compute the stability path by stability selection
#'
#' @name stability_selection
#' @param subsamples a list of vectors describing the subsamples. The number of vectors (or list length) determines th number of subsamples used in the stability selection. Automatically set to 20 subsamples with size \code{10*sqrt(n)} if \code{n >= 144} and \code{0.8*n} otherwise following Liu et al. (2010) recommandations.
#' @param control a list controling the main optimization process in each call to PLNnetwork. See \code{\link[=PLNnetwork]{PLNnetwork}} for details.
#' @param mc.cores the number of cores to used. Default is 1.
#'
##' @return the list of subsamples. The estimated probabilities of selection of the edges are stored in the fields stability_path of the PLNnetwork object
NULL
PLNnetworkfamily$set("public", "stability_selection",
  function(subsamples = NULL, control = list(), mc.cores = 1) {

    ## define default control parameters for optim and overwrite by user defined parameters
    ctrl_init <- PLNnetwork_param(list(), private$n, private$p, "init")
    ctrl_init$trace <- 0
    ctrl_main <- PLNnetwork_param(control, n, p, "main")
    ctrl_main$trace <- 0

    ## select default subsamples according
    if (is.null(subsamples)) {
      subsample.size <- round(ifelse(private$n >= 144, 10*sqrt(private$n), 0.8*private$n))
      subsamples <- replicate(20, sample.int(private$n, subsample.size), simplify = FALSE)
    }

    ## got for stability selection
    cat("\nStability Selection for PLNnetwork: ")
    cat("\nsubsampling: ")

    stabs_out <- mclapply(subsamples, function(subsample) {
      cat("+")
      inception_ <- self$getModel(self$penalties[1])
      inception_$update(
        M = inception_$var_par$M[subsample, ],
        S = inception_$var_par$S[subsample, ]
      )
      ctrl_init$inception <- inception_
      myPLN <- PLNnetworkfamily$new(penalties  = self$penalties,
                                    responses  = self$responses [subsample, , drop = FALSE ],
                                    covariates = self$covariates[subsample, , drop = FALSE],
                                    offsets    = self$offsets   [subsample, , drop = FALSE ], control = ctrl_init)

      myPLN$optimize(ctrl_main)
      nets <- do.call(cbind, lapply(myPLN$models, function(model) {
        as.matrix(model$latent_network("support"))[upper.tri(diag(private$p))]
      }))
      nets
    }, mc.cores = mc.cores)

    prob <- Reduce("+", stabs_out, accumulate = FALSE) / length(subsamples)
    ## formatting/tyding
    colnames(prob) <- self$penalties
    private$stab_path <- prob %>%
      as.data.frame() %>%
      mutate(Edge = 1:n()) %>%
      gather(key = "Penalty", value = "Prob", -Edge) %>%
      mutate(Penalty = as.numeric(Penalty),
             Node1   = as.character(edge_to_node(Edge)$node1),
             Node2   = as.character(edge_to_node(Edge)$node2),
             Edge    = paste0(Node1, "--", Node2)) %>%
      filter(Node1 < Node2)

    invisible(subsamples)
  }
)


#' Extract the regularization path of a PLNnetwork fit
#'
#' @name coefficient_path
#' @param precision a logical, should the coefficients of the precision matrix Omega or the covariance matrice Sigma be sent back. Default is \code{TRUE}.
#' @param corr a logical, should the correlation (partial in case  precision = TRUE) be sent back. Default is \code{TRUE}.
#'
#' @return  Send back a tibble/data.frame.
NULL
PLNnetworkfamily$set("public", "coefficient_path",
function(precision = TRUE, corr = TRUE) {
  lapply(self$penalties, function(x) {
    if (precision) {
      G <- self$getModel(x)$model_par$Omega
    } else {
      G <- self$getModel(x)$model_par$Sigma
      dimnames(G) <- dimnames(self$getModel(x)$model_par$Omega)
    }
    if (corr) {
      G <- ifelse(precision, -1, 1) * G / tcrossprod(sqrt(diag(G)))
    }
    G %>% melt(value.name = "Coeff", varnames = c("Node1", "Node2")) %>%
      mutate(Penalty = x,
             Node1   = as.character(Node1),
             Node2   = as.character(Node2),
             Edge    = paste0(Node1, "--", Node2)) %>%
      filter(Node1 < Node2)
  }) %>% bind_rows()
})

#' Best model extraction from a collection of PLNnetworkfit
#'
#' @name PLNnetworkfamily_getBestModel
#'
#' @param crit a character for the criterion used to performed the selection. Either
#' "BIC", "EBIC", "StARS", "R_squared". Default is \code{BIC}. If StARS
#' (Stability Approach to Regularization Selection) is chosen and stability selection
#'  was not yet performed, the function will call the method stability_selection with default argument.
#' @param stability a scalar, indicating the target stability (= 1 - 2 beta) at which the network is selected. Default is \code{0.9}.
#' @return  Send back a object with class \code{\link[=PLNnetworkfit]{PLNnetworkfit}}.
NULL
PLNnetworkfamily$set("public", "getBestModel",
function(crit = c("BIC", "ICL", "loglik", "R_squared", "EBIC", "StARS"), stability = 0.9){
  crit <- match.arg(crit)
  if (crit == "StARS") {
    if (is.null(private$stab_path)) self$stability_selection()
    id_stars <- self$criteria %>%
      select(param, stability) %>% rename(Stability = stability) %>%
      filter(Stability > stability) %>%
      pull(param) %>% min() %>% match(self$penalties)
    model <- self$models[[id_stars]]$clone()
  } else {
    model <- super$getBestModel(crit)
  }
  model
})

## ----------------------------------------------------------------------
## PUBLIC PLOTTING METHODS
## ----------------------------------------------------------------------

#' @export
PLNnetworkfamily$set("public", "plot",
function(criteria = c("loglik", "pen_loglik", "BIC", "EBIC"), log.x = TRUE) {
  vlines <- sapply(intersect(criteria, c("BIC", "EBIC")) , function(crit) self$getBestModel(crit)$penalty)
  p <- super$plot(criteria, FALSE) + xlab("penalty") + geom_vline(xintercept = vlines, linetype = "dashed", alpha = 0.25)
  if (log.x) p <- p + ggplot2::coord_trans(x = "log10")
  p
})

#' @export
PLNnetworkfamily$set("public", "plot_stars",
function(stability = 0.9, log.x = TRUE) {
  dplot <- self$criteria %>% select(param, density, stability) %>%
    rename(Penalty = param) %>%
    gather(key = "Metric", value = "Value", stability:density)
  penalty_stars <- dplot %>% filter(Metric == "stability" & Value >= stability) %>%
    pull(Penalty) %>% min()

  p <- ggplot(dplot, aes(x = Penalty, y = Value, group = Metric, color = Metric)) +
    geom_point() +  geom_line() + theme_bw() +
    ## Add information correspinding to best lambda
    geom_vline(xintercept = penalty_stars, linetype = 2) +
    geom_hline(yintercept = stability, linetype = 2) +
    annotate(x = penalty_stars, y = 0,
           label = paste("lambda == ", round(penalty_stars, 5)),
           parse = TRUE, hjust = -0.05, vjust = 0, geom = "text") +
    annotate(x = penalty_stars, y = stability,
             label = paste("stability == ", stability),
             parse = TRUE, hjust = -0.05, vjust = 1.5, geom = "text")
  if (log.x) p <- p + ggplot2::scale_x_log10() + annotation_logticks(sides = "b")
  p
})


#' @export
PLNnetworkfamily$set("public", "plot_objective",
function() {
  objective <- unlist(lapply(self$models, function(model) model$optim_par$objective))
  changes <- cumsum(unlist(lapply(self$models, function(model) model$optim_par$outer_iterations)))
  dplot <- data.frame(iteration = 1:length(objective), objective = objective)
  p <- ggplot(dplot, aes(x = iteration, y = objective)) + geom_line() +
    geom_vline(xintercept = changes, linetype="dashed", alpha = 0.25) +
    ggtitle("Objective along the alternate algorithm") + xlab("iteration (+ changes of model)") +
    annotate("text", x = changes, y = min(dplot$objective), angle = 90,
             label = paste("penalty=",format(self$criteria$param, digits = 1)), hjust = -.1, size = 3, alpha = 0.7) + theme_bw()
  p
})

PLNnetworkfamily$set("public", "show",
function() {
  super$show()
  cat(" Task: Network Inference \n")
  cat("========================================================\n")
  cat(" -", length(self$penalties) , "penalties considered: from", min(self$penalties), "to", max(self$penalties),
      "\n", "   use $penalties to see all values and access specific lambdas", "\n")
  if (!anyNA(self$criteria$stability))
    cat(" - Best model (regarding StARS): lambda =", format(self$getBestModel("StARS")$penalty, digits = 3), "\n")
  if (!anyNA(self$criteria$BIC))
    cat(" - Best model (regarding BIC): lambda =", format(self$getBestModel("BIC")$penalty, digits = 3), "\n")
})
PLNnetworkfamily$set("public", "print", function() self$show())

