# PLNBLOCK optimization:
#
# init_parameters: named list(B, M, S, Tau)
# Y, X, O: observed data, constant
# configuration: named list(algorithm, xtol_rel, xtol_abs, ftol_rel, ftol_abs, maxeval, maxtime)
#
# Returns named list(parameters: named list, nb_iter: int, stop_reason: str)
#
# Configuration :
# Forwarded to nloptr iterative steps (c,e,f).
# Except for algorithm, the keys are optional termination criteria. Not all of them are needed.
# A missing key or a key set to a negative number means disabled.
# xtol_abs: double for every parameter, or named list with values for each parameter.
#
# xtol_abs, xtol_rel, maxeval are also used for overall loop.
# ftol_abs and ftol_rel could maybe be used if defined on J(params), but I'm not sure.
#
# Dimensions are checked only on C++ side.
optimize_plnblockbis <- function(data, params, config) {

  # Link to the approximate function to optimize Omega, depending on the target structure
  optim_plnblockbis_Omega <- ifelse(is.null(params$rho), optim_plnblockbis_Omega, function(M, S)
    optim_plnblockbis_Omega_sparse(M, S, data$w, rho = config$rho))

  # Main loop
  nb_iter <- 0
  new_params <- params
  posteriorProb <- list(new_params$Tau)
  criterion <- vector("numeric", config$maxit_out)
  objective <- Inf
  repeat {
    # M Step

    optim_Omega <- optim_plnblockbis_Omega(M = new_params$M, S = new_params$S)
    new_params$Omega <- optim_Omega$Omega
    new_params$B <- optim_plnblockbis_B(XtXm = data$XtXm,
                                            X = data$X,
                                            Mu = new_params$Mu)
    new_params$dm <- optim_plnblockbis_dm(
      X = data$X,
      B = new_params$B,
      Mu = new_params$Mu,
      Delta = new_params$Delta
    )

    # VE Step
    optim_VE <- optim_plnblockbis_VE(data, new_params, config)
    new_params$Mu    <- optim_VE$Mu
    new_params$M     <- optim_VE$M
    new_params$Delta <- optim_VE$Delta
    new_params$S     <- optim_VE$S
    new_params$Tau   <- optim_VE$Tau

    # optim_Tau <- optim_plnblockbis_Tau(data, new_params)
    # new_params$Tau <- optim_Tau$Tau

    # Going next step and assessing convergence
    nb_iter <- nb_iter + 1
    criterion[nb_iter] <- new_objective <- -plnblockbis_loglik(data, new_params)
    objective_converged <-
      abs(new_objective - objective) / abs(new_objective) < config$ftol_out
    params_converged <- parameter_list_converged(
      params,
      new_params,
      xtol_abs = config$xtol_abs,
      xtol_rel = config$xtol_rel
    )
    maxit_reached <- config$maxit_out >= 0 &&
      nb_iter >= config$maxit_out

    if (params_converged |
        objective_converged | maxit_reached) {
      if (params_converged)
        statut <- 4
      if (objective_converged)
        statut <- 3
      if (maxit_reached)
        statut <- 5
      break
    }
    posteriorProb <- append(posteriorProb, list(new_params$Tau))
    params <- new_params
    objective  <- new_objective
  }
  out   <- new_params
  out$D <- diag(as.numeric(1/out$dm)); out$dm <- NULL
  out$A <- optim_VE$A
  out$Z <- data$O + out$Mu + out$M %*% out$Tau
  out$Sigma <- optim_Omega$Sigma
  out$Ji <- plnblockbis_vloglik(data, new_params)
  attr(out$Ji, "weights") <- data$w
  attr(out$Ji, "loglik") <- -objective
  out$monitoring <- list(
    objective  = criterion[1:nb_iter],
    outer_iterations = nb_iter,
    status     = statut,
    posteriorProb = posteriorProb,
    backend = "nlopt-vem"
  )
  out
}

#' @importFrom glassoFast glassoFast
optim_plnblockbis_Omega_sparse <- function(M, S, w, rho) {
  n <- sum(w)
  p <- ncol(M)
  glassoFast::glassoFast(crossprod(M) / n + diag(crossprod(w, (S * S)), p, p), rho = rho)$wi
}

# Test convergence for a named list of parameters
# oldp, newp: named list of parameters
# xtol_rel: double ; negative or NULL = disabled
# xtol_abs: double (negative or NULL = disabled) or list of double/matrices for each parameter (any negative is disabled)
# Returns boolean
parameter_list_converged <- function(oldp,
                                     newp,
                                     xtol_abs = NULL,
                                     xtol_rel = NULL) {
  # Strategy is to compare each pair of list elements with matching names.
  # Named lists are just vectors (T,str) using order of insertion.
  # mapply() is handy to do the pair tests, but it works on the underlying vector order (ignoring names).
  # So reorder lists by their names to use mapply.
  stopifnot(is.list(oldp), is.list(newp))
  oldp <- oldp[order(names(oldp))]
  newp <- newp[order(names(newp))]
  stopifnot(all(names(oldp) == names(newp)))

  # Check convergence with xtol_rel if enabled
  if (is.double(xtol_rel) && xtol_rel > 0) {
    if (all(mapply(function(o, n) {
      all(abs(n - o) <= xtol_rel * abs(o))
    }, oldp, newp))) {
      return(TRUE)
    }
  }

  # Check convergence with xtol_abs (homogeneous) if enabled
  if (is.double(xtol_abs) && xtol_abs > 0) {
    if (all(mapply(function(o, n) {
      all(abs(n - o) <= xtol_abs)
    }, oldp, newp))) {
      return(TRUE)
    }
  }

  # If no criteria has triggered, indicate no convergence
  FALSE
}
