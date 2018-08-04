.logfactorial <- function(n) { # Ramanujan's formula
  n[n == 0] <- 1 ## 0! = 1!
  return(n*log(n) - n + log(8*n^3 + 4*n^2 + n + 1/30)/6 + log(pi)/2)
}

.trunc <- function(x, M = 300) {
  return(pmin(x, M))
}

.loglikPLN <- function(Y, X, O, Theta, Sigma, M, S) {
  Omega <- chol2inv(chol(Sigma))
  Z <- O + M + tcrossprod(X, Theta)
  A = exp (Z + .5 * S)
  res <- Y * Z - A + .5*log(S) + .5 - .logfactorial(Y) - 0.5*(M * (M %*% Omega) + t(t(S) * diag(Omega)) )
  return(rowSums(res) + 0.5 * log(det(Omega)) )
}

edge_to_node <- function(x, n = max(x)) {
  x <- x - 1 ## easier for arithmetic to number edges starting from 0
  n.node <- round((1 + sqrt(1 + 8*n)) / 2) ## n.node * (n.node -1) / 2 = n (if integer)
  j.grid <- cumsum(0:n.node)
  j <- findInterval(x, vec = j.grid)
  i <- x - j.grid[j]
  ## Renumber i and j starting from 1 to stick with R convention
  return(data.frame(node1 = i + 1, node2 = j + 1))
}

node_pair_to_egde <- function(x, y, node.set = union(x, y)) {
  ## Convert node labels to integers (starting from 0)
  x <- match(x, node.set) - 1
  y <- match(y, node.set) - 1
  ## For each pair (x,y) return, corresponding edge number
  n <- length(node.set)
  j.grid <- cumsum(0:(n - 1))
  x + j.grid[y] + 1
}

logLikPoisson <- function(responses, lambda) {
  loglik <- sum(responses * lambda, na.rm=TRUE) - sum(exp(lambda)) - sum(.logfactorial(responses))
  loglik
}

nullModelPoisson <- function(responses, covariates, offsets) {
  Theta <- do.call(rbind, lapply(1:ncol(responses), function(j)
    coefficients(glm.fit(covariates, responses[, j], offset = offsets[,j], family = stats::poisson()))))
  lambda <- offsets + tcrossprod(covariates, Theta)
  lambda
}

fullModelPoisson <- function(responses) {
  lambda <- log(responses)
  lambda
}


PLN_param <- function(control, n, p) {
  ctrl <- list(
    ftol_rel  = ifelse(n < 1.5*p, 1e-6, 1e-8),
    ftol_abs  = 0,
    xtol_rel  = 1e-4,
    xtol_abs  = 1e-4,
    maxeval   = 10000,
    method    = "CCSAQ",
    lbvar     = 1e-4,
    trace     = 1,
    inception = "LM"
  )
  ctrl[names(control)] <- control
  ctrl
}

PLNPCA_param <- function(control, n, p, type = c("init", "main")) {
  type <- match.arg(type)

  ctrl <- switch(match.arg(type),
    "init" = list(
      inception = ifelse(n >= 1.5*p, "PLN", "LM"),
      ftol_rel  = ifelse(n < 1.5*p, 1e-6, 1e-8),
      ftol_abs = 0,
      xtol_rel = 1e-4,
      xtol_abs = 1e-4,
      maxeval  = 10000,
      method   = "CCSAQ",
      lbvar    = 1e-4,
      trace    = 0
    ),
    "main" = list(
      ftol_rel = 1e-6,
      ftol_abs = 0,
      xtol_rel = 1e-4,
      xtol_abs = 1e-4,
      maxeval  = 10000,
      method   = "MMA",
      lbvar    = 1e-4,
      trace    = 1,
      cores    = 1
    )
  )
  ctrl[names(control)] <- control
  ctrl
}

PLNnetwork_param <- function(control, n, p, type = c("init", "main")) {
  type <- match.arg(type)

  ctrl <- switch(match.arg(type),
    "init" = list(
      inception = ifelse(n >= 1.5*p, "PLN", "LM"),
      ftol_rel  = ifelse(n < 1.5*p, 1e-6, 1e-8),
      ftol_abs = 0,
      xtol_rel = 1e-4,
      xtol_abs = 1e-4,
      maxeval  = 10000,
      method   = "CCSAQ",
      lbvar    = 1e-4,
      nPenalties = 30,
      min.ratio = ifelse(n < 1.5*p, 0.1, 0.05),
      trace = 0),
    "main" = list(
      ftol_out  = 1e-5,
      maxit_out = 50,
      penalize.diagonal = FALSE,
      warm      = FALSE,
      ftol_abs  = 0,
      ftol_rel  = 1e-8,
      xtol_rel  = 1e-4,
      xtol_abs  = 1e-4,
      maxeval   = 10000,
      method    = "CCSAQ",
      lbvar     = 1e-4,
      trace = 1)
  )
  ctrl[names(control)] <- control
  ctrl
}

statusToMessage <- function(status) {
    message <- switch( status,
        "1"  = "NLOPT_SUCCESS: Generic success return value.",
        "2"  = "NLOPT_STOPVAL_REACHED: Optimization stopped because stopval (above) was reached.",
        "3"  = "NLOPT_FTOL_REACHED: Optimization stopped because ftol_rel or ftol_abs (above) was reached.",
        "4"  = "NLOPT_XTOL_REACHED: Optimization stopped because xtol_rel or xtol_abs (above) was reached.",
        "5"  = "NLOPT_MAXEVAL_REACHED: Optimization stopped because maxeval (above) was reached.",
        "6"  = "NLOPT_MAXTIME_REACHED: Optimization stopped because maxtime (above) was reached.",
        "-1" = "NLOPT_FAILURE: Generic failure code.",
        "-2" = "NLOPT_INVALID_ARGS: Invalid arguments (e.g. lower bounds are bigger than upper bounds, an unknown algorithm was specified, etcetera).",
        "-3" = "NLOPT_OUT_OF_MEMORY: Ran out of memory.",
        "-4" = "NLOPT_ROUNDOFF_LIMITED: Roundoff errors led to a breakdown of the optimization algorithm. In this case, the returned minimum may still be useful.",
        "-5" = "Halted because of a forced termination: the user called nlopt_force_stop(opt) on the optimization's nlopt_opt object opt from the user's objective function.",
        "Return status not recognized."
    )
    message
}
