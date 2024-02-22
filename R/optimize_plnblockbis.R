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
  optim_plnblockbis_Omega <- ifelse(is.null(params$rho),
        optim_plnblockbis_Omega,
        function(M, S) optim_plnblockbis_Omega_sparse(M, S, data$w, rho = config$rho)
    )
  # Main loop
  nb_iter <- 0
  parameters <- params; parameters$Omega <- diag(1, ncol(params$M), ncol(params$M))
  parameters$D <- diag(1, ncol(params$Mu), ncol(params$Mu))
  parameters$T <- parameters$Tau; parameters$Tau <- NULL
  new_parameters <- parameters
  posteriorProb <- list(new_parameters$T)
  criterion <- vector("numeric", config$maxit_out)
  objective <- Inf
  repeat {

    # M Step
    optim_B <- optim_plnblockbis_B(data, new_parameters, config)
    new_parameters$B <- optim_B$B
    optim_Omega <- optim_plnblockbis_Omega(M = new_parameters$M, S = new_parameters$S, w = data$w)
    new_parameters$Omega <- optim_Omega$Omega

    # VE Step
    optim_VE <- optim_plnblockbis_VE(data, new_parameters, config)
    new_parameters$M <- optim_VE$M
    new_parameters$S <- optim_VE$S
    new_parameters$Mu <- optim_VE$Mu
    new_parameters$Delta <- optim_VE$Delta
    new_parameters$D <- optim_VE$D

    new_parameters$T <- optim_VE$Tau
    ##########

    if(config$g_resampling > 0){
      max_values <- apply(new_parameters$T , 2, function(col) max(col, na.rm = TRUE))
      if((nb_iter %% 1) == 0){
        to_resample = which(max_values < 0.99)
        for(q in to_resample){
          probabilities = new_parameters$T[,q]
          ##
          init_group = which.max(probabilities)
          ##
          Q = length(probabilities)
          new_group = sample(1:Q, size=1, prob=probabilities)
          new_tau = rep(1e-16, Q)
          new_tau[[new_group]] = 1
          new_parameters$T[,q] = new_tau
        }}
      if(config$g_resampling > 1){
        if((nb_iter %% 3) == 0){
          alpha = rowMeans(new_parameters$T)
          Q = length(alpha)
          to_resample = which(max_values > 0.99)
          if(length(to_resample) > 0){
            chosen = sample(to_resample, size = max(2, floor(0.05 * ncol(new_parameters$T))))
            for(c in chosen){
              new_group = sample(1:Q, size=1, prob=alpha)
              new_tau = rep(1e-16, Q)
              new_tau[[new_group]] = 1
              new_parameters$T[,chosen] = new_tau
            }}
        }
      }
    }


    # Going next step and assessing convergence
    nb_iter <- nb_iter + 1
    criterion[nb_iter] <- new_objective <- - plnblockbis_loglik(data, new_parameters)
    #print("reached loglik computation")

    objective_converged <-
      abs(new_objective - objective)/abs(new_objective) < config$ftol_out
    parameters_converged <- parameter_list_converged(
      parameters, new_parameters,
      xtol_abs = config$xtol_abs, xtol_rel = config$xtol_rel
    )
    maxit_reached <- config$maxit_out >= 0 && nb_iter >= config$maxit_out

    if (parameters_converged | objective_converged | maxit_reached) {
      if (parameters_converged) statut <- 4
      if (objective_converged)  statut <- 3
      if (maxit_reached) statut <- 5
      break
    }
    posteriorProb <- append(posteriorProb, list(new_parameters$T))
    parameters <- new_parameters
    objective  <- new_objective
  }
  out   <- new_parameters
  out$A <- optim_VE$A
  out$Z <- data$O + out$Mu + out$M %*% out$T
  out$Sigma <- optim_Omega$Sigma
  ## J'ai remplacé vloglik par loglik ici
  out$Ji <- plnblockbis_vloglik(data, new_parameters)
  attr(out$Ji, "weights") <- data$w
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
  n <- sum(w); p <- ncol(M)
  glassoFast::glassoFast( crossprod(M)/n + diag(crossprod(w, (S * S)), p, p), rho = rho )$wi
}

# Test convergence for a named list of parameters
# oldp, newp: named list of parameters
# xtol_rel: double ; negative or NULL = disabled
# xtol_abs: double (negative or NULL = disabled) or list of double/matrices for each parameter (any negative is disabled)
# Returns boolean
parameter_list_converged <- function(oldp, newp, xtol_abs = NULL, xtol_rel = NULL) {
  # Strategy is to compare each pair of list elements with matching names.
  # Named lists are just vectors (T,str) using order of insertion.
  # mapply() is handy to do the pair tests, but it works on the underlying vector order (ignoring names).
  # So reorder lists by their names to use mapply.
  stopifnot(is.list(oldp), is.list(newp))
  oldp <- oldp[order(names(oldp))]
  newp <- newp[order(names(newp))]
  stopifnot(all(names(oldp) == names(newp)))

  # Check convergence with xtol_rel if enabled
  if(is.double(xtol_rel) && xtol_rel > 0) {
    if(all(mapply(function(o, n) { all(abs(n - o) <= xtol_rel * abs(o)) }, oldp, newp))) {
      return(TRUE)
    }
  }

  # Check convergence with xtol_abs (homogeneous) if enabled
  if(is.double(xtol_abs) && xtol_abs > 0) {
    if(all(mapply(function(o, n) { all(abs(n - o) <= xtol_abs) }, oldp, newp))) {
      return(TRUE)
    }
  }

  # If no criteria has triggered, indicate no convergence
  FALSE
}
