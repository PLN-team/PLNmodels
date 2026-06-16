available_algorithms_nlopt <- c("CCSAQ", "MMA", "LBFGS", "VAR1", "VAR2")
available_algorithms_torch <- c("RPROP", "RMSPROP", "ADAM", "ADAGRAD")

## Column-normalize a covariate matrix so that all backends receive a
## well-conditioned X regardless of covariate magnitudes.  The transformation
## X_sc = X / scale, B_sc = B * scale satisfies X_sc %*% B_sc = X %*% B
## exactly, so the PLN model is unchanged.  The pmax(..., 1) guard prevents
## zero-variance or already-unit-norm columns from being inflated.
normalize_covariates <- function(X) {
  scales <- pmax(sqrt(colSums(X^2)), 1)
  list(X_sc = sweep(X, 2, scales, "/"), scales = scales)
}

config_default_nlopt <-
  list(
    algorithm     = "CCSAQ",
    backend       = "nlopt",
    maxeval       = 10000  ,
    ftol_rel      = 1e-8   ,
    xtol_rel      = 1e-6   ,
    ftol_abs      = 0.0    ,
    xtol_abs      = 0.0    ,
    maxtime       = -1
  )


config_default_builtin <-
  list(
    algorithm           = "NEWTON",
    backend             = "builtin",
    maxeval             = 10000,
    ftol_in             = 1e-8,
    maxit_em            = 200,
    ftol_em             = 1e-8
  )


# PLNPCA builtin backend: joint L-BFGS on [vec(B); vec(C); vec(M); vec(ψ)] with strong Wolfe
# line search (m=10 pairs). Only maxeval and ftol_in are read by the C++ optimizer.
config_default_plnpca <-
  list(
    backend = "builtin",
    maxeval = 10000,
    ftol_in = 1e-8
  )

config_default_torch <-
  list(
    algorithm     = "RPROP",
    backend       = "torch",
    maxeval       = 10000  ,
    num_epoch     = 1000   ,
    num_batch     = 1      ,
    ftol_rel      = 1e-8   ,
    xtol_rel      = 1e-6   ,
    lr            = 0.1    ,
    momentum      = 0.05   ,
    weight_decay  = 0      ,
    step_sizes    = c(1e-3, 50),
    etas          = c(0.5, 1.2),
    centered      = FALSE,
    trace         = 1,
    device        = "cpu"
  )

## Build the optimizer config list from a backend name and user overrides.
## `builtin_default` lets PLNPCA pass config_default_plnpca instead of config_default_builtin.
## `extra` is a named list of additional defaults applied BEFORE user overrides (so the user can
## still override them), used for outer-loop parameters like ftol_em/maxit_em in PLNnetwork and
## PLNmixture.
make_config_optim <- function(backend, config_optim, trace,
                              builtin_default = config_default_builtin,
                              extra = list()) {
  config_opt <- if (backend == "nlopt") {
    stopifnot(config_optim$algorithm %in% available_algorithms_nlopt)
    config_default_nlopt
  } else if (backend == "torch") {
    stopifnot(config_optim$algorithm %in% available_algorithms_torch)
    config_default_torch
  } else {
    builtin_default
  }
  config_opt$trace <- trace
  config_opt[names(extra)] <- extra
  config_opt[names(config_optim)] <- config_optim
  config_opt
}

config_post_default_PLN <-
  list(
    jackknife       = FALSE,
    bootstrap       = 0L,
    rsquared        = TRUE,
    variational_var = FALSE,
    sandwich_var    = FALSE
  )

config_post_default_PLNnetwork <-
  list(
    jackknife       = FALSE,
    bootstrap       = 0L,
    rsquared        = FALSE,
    variational_var = FALSE,
    sandwich_var    = FALSE
  )

config_post_default_PLNLDA <-
  list(
    jackknife       = FALSE,
    bootstrap       = 0L,
    rsquared        = TRUE,
    variational_var = FALSE,
    sandwich_var    = FALSE
  )

config_post_default_PLNPCA <-
  list(
    jackknife       = FALSE,
    bootstrap       = 0L,
    rsquared        = TRUE,
    variational_var = FALSE,
    sandwich_var    = FALSE
  )

config_post_default_PLNmixture <-
  list(
    jackknife       = FALSE,
    bootstrap       = 0L,
    rsquared        = TRUE,
    variational_var = FALSE,
    sandwich_var    = FALSE
  )

.xlogx <- function(x) ifelse(x < .Machine$double.eps, 0, x*log(x))

.softmax <- function(x) {
  b <- max(x)
  exp(x - b) / sum(exp(x - b))
}

.logit <- function(x) log(x/(1 - x))

.check_boundaries <- function(x, zero = .Machine$double.eps) {
  x[is.nan(x)] <- zero
  x[x > 1 - zero] <- 1 - zero
  x[x <     zero] <-     zero
  x
}

.logfactorial_torch <- function(n){
  n[n == 0] <- 1 ## 0! = 1!
  n*torch_log(n) - n + torch_log(8*torch_pow(n,3) + 4*torch_pow(n,2) + n + 1/30)/6 + torch_log(pi)/2
}

.logfactorial <- function(n) { # Ramanujan's formula
  n[n == 0] <- 1 ## 0! = 1!
  n*log(n) - n + log(8*n^3 + 4*n^2 + n + 1/30)/6 + log(pi)/2
}

as_indicator <- function(clustering) {
  K <- length(unique(clustering))
  N  <- length(clustering)
  Z <- matrix(0, N, K)
  Z[cbind(seq.int(N), clustering)] <- 1
  Z
}

logLikPoisson <- function(responses, lambda, weights = rep(1, nrow(responses))) {
  loglik <- rowSums(responses * lambda, na.rm = TRUE) - rowSums(exp(lambda)) - rowSums(.logfactorial(responses))
  loglik <- sum(loglik * weights)
  loglik
}

#' @importFrom stats glm.fit glm.control
nullModelPoisson <- function(responses, covariates, offsets, weights = rep(1, nrow(responses))) {
### TODO: use fastglm
  B <- do.call(cbind, parallel::mclapply(1:ncol(responses), function(j)
    coefficients(suppressWarnings(
      glm.fit(covariates, responses[, j], weights = weights, offset = offsets[, j], family = stats::poisson(),
        control = glm.control(epsilon = 1e-3, maxit = 10)))),
    mc.cores = getOption("mc.cores", 1L)))
  offsets + covariates %*% B
}

#' @importFrom stats .getXlevels
extract_model <- function(call, envir) {
  ## extract relevant arguments from the high level call for the model frame
  call_args <- call[match(c("formula", "data", "subset", "weights"), names(call), 0L)]
  call_args <- c(as.list(call_args), list(xlev = attr(call$formula, "xlevels"), na.action = NULL))
  ## eval the call in the parent environment
  frame <- do.call(stats::model.frame, call_args, envir = envir)
  ## create the set of matrices to fit the PLN model
  Y <- model.response(frame)
  ## model.response oversimplifies into a numeric when a single variable is involved
  if (is.null(dim(Y))) Y <- matrix(Y, nrow = length(Y), ncol = 1)
  if (ncol(Y) == 1 & is.null(colnames(Y))) colnames(Y) <- "Y"
  X <- model.matrix(terms(frame), frame)
  O <- model.offset(frame)
  if (is.null(O)) O <- matrix(0, nrow(Y), ncol(Y))
  if (is.vector(O)) O <- O %o% rep(1, ncol(Y))
  w <- model.weights(frame)
  if (is.null(w)) {
    w <- rep(1.0, nrow(Y))
  } else {
    stopifnot(all(w > 0) && length(w) == nrow(Y))
  }
  ## Save encountered levels for predict methods as attribute of the formula.
  ## Evaluate the formula expression to get the formula object before setting
  ## attributes — avoids "cannot set an attribute on a 'symbol'" when the
  ## formula was passed as a variable (e.g. PLN(my_formula, data = d)).
  formula_obj <- if (!inherits(call$formula, "formula")) {
    eval(call$formula, envir = envir)
  } else {
    call$formula
  }
  attr(formula_obj, "xlevels") <- .getXlevels(terms(frame), frame)
  list(Y = Y, X = X, O = O, miss = is.na(Y), w = w, formula = formula_obj)
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

#' @title PLN RNG
#'
#' @description Random generation for the PLN model with latent mean equal to mu, latent covariance matrix
#'              equal to Sigma and average depths (sum of counts in a sample) equal to depths
#'
#' @param n the sample size
#' @param mu vectors of means of the latent variable
#' @param Sigma covariance matrix of the latent variable
#' @param depths Numeric vector of target depths. The first is recycled if there are not `n` values
#'
#' @return a n * p count matrix, with row-sums close to depths, with an attribute "offsets"
#' corresponding to the true generated offsets (in log-scale).
#'
#' @details The default value for mu and Sigma assume equal abundances and no correlation between
#'          the different species.
#'
#' @rdname rPLN
#' @examples
#' ## 10 samples of 5 species with equal abundances, no covariance and target depths of 10,000
#' rPLN()
#' ## 2 samples of 10 highly correlated species with target depths 1,000 and 100,000
#' ## very different abundances
#' mu <- rep(c(1, -1), each = 5)
#' Sigma <- matrix(0.8, 10, 10); diag(Sigma) <- 1
#' rPLN(n=2, mu = mu, Sigma = Sigma, depths = c(1e3, 1e5))
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats rpois
#' @export
rPLN <- function(n = 10, mu = rep(0, ncol(Sigma)), Sigma = diag(1, 5, 5),
                 depths = rep(1e4, n))  {
  p <- ncol(Sigma)
  if (any(is.vector(mu), ncol(mu) == 1)) {
    mu <- matrix(rep(mu, n), ncol = p, byrow = TRUE)
  }
  if (length(depths) != n) {
    depths <- rep(depths[1], n)
  }
  ## adjust depths
  exp_depths <- rowSums(exp(rep(1, n) %o% diag(Sigma)/2 + mu)) ## sample-wise expected depths
  offsets <- log(depths %o% rep(1, p)) - log(exp_depths)
  Z <- mu + mvrnorm(n, rep(0,ncol(Sigma)), as.matrix(Sigma)) + offsets
  Y <- matrix(rpois(n * p, as.vector(exp(Z))), n, p)
  dimnames(Y) <- list(paste0("S", 1:n), paste0("Y", 1:p))
  attr(Y, "offsets") <- offsets
  Y
}

# Internal function
#' @importFrom stats rnorm
create_parameters <- function(
    n = 200,
    p = 50,
    d = 2,
    rho = 0.2,
    sigma = 1,
    depths = 100000,
    ...
) {
  ## Sigma chosen to achieve a given snr
  list(n      = n,
       p      = p,
       X      = matrix(rnorm(n*d), nrow = n, ncol = d,
                       dimnames = list(paste0("S", 1:n), paste0("Var_", 1:d))),
       B      = matrix(rnorm(n = p*d, sd = 1/sqrt(d)), nrow = d, ncol = p),
       Sigma  = sigma * toeplitz(x = rho^seq(0, p-1)),
       depths = depths)
}

#' Helper function for PLN initialization.
#'
#' @description
#' Barebone function to compute starting points for B, M and S2 when fitting a PLN. Mostly intended for internal use.
#'
#' @param Y Response count matrix
#' @param X Covariate matrix. Note that initialization will fail if the model matrix is singular.
#' @param O Offset matrix (in log-scale)
#' @param w Weight vector (defaults to 1)
#' @param method character: strategy used to initialize B. Either `"LM"` (default, fast weighted
#'   log-linear regression) or `"GLM"` (p independent Poisson GLMs, more accurate for complex
#'   or unbalanced designs but slower).
#' @return a named list of starting values for model parameter B and variational parameters M and S2 used in the iterative optimization algorithm of [PLN()]
#'
#' @details
#' * **B**: estimated by weighted LM (`method = "LM"`, default) or p independent Poisson GLMs
#'   (`method = "GLM"`). The GLM option gives better B estimates for factorial or unbalanced
#'   designs at the cost of p IRLS fits.
#' * **M**: initialized to `log((1 + Y) / exp(O))` (M_full in the X*B + M_res parameterization).
#' * **S**: initialized element-wise to `1 / sqrt(2 + Y)`, the approximate VE-step optimum at
#'   Omega = I. This adapts automatically to count levels: high S for zero counts (high
#'   uncertainty), low S for large counts.
#'
#' @rdname compute_PLN_starting_point
#' @examples
#' \dontrun{
#' data(barents)
#' Y <- barents$Abundance
#' X <- model.matrix(Abundance ~ Latitude + Longitude + Depth + Temperature, data = barents)
#' O <- log(barents$Offset)
#' w <- rep(1, nrow(Y))
#' compute_PLN_starting_point(Y, X, O, w)
#' compute_PLN_starting_point(Y, X, O, w, method = "GLM")
#' }
#'
#' @importFrom stats lm.fit glm.fit poisson
#' @export
compute_PLN_starting_point <- function(Y, X, O, w, method = c("LM", "GLM")) {
  method <- match.arg(method)
  n <- nrow(Y); p <- ncol(Y); d <- ncol(X)
  Y0 <- replace(Y, is.na(Y), 0)  # treat missing counts as 0 for initialization only
  expO <- exp(O)
  if (method == "GLM") {
    pois_fam <- poisson()
    B <- vapply(seq_len(p), function(j)
      glm.fit(X, Y0[, j], offset = O[, j], weights = w, family = pois_fam)$coefficients,
      FUN.VALUE = numeric(d)
    )
  } else {
    B <- lm.fit(w * X, w * log((1 + Y0) / expO), singular.ok = TRUE)$coefficients
    B[is.na(B)] <- 0
  }
  list(B  = matrix(B, d, p),
       M  = matrix(log((1 + Y0) / expO), n, p),
       S2 = 1 / (2 + Y0))
}

#' Helper function for ZIPLN initialization.
#'
#' @description
#' Fast LM-based starting point for ZIPLN: one multivariate `lm.fit` for the PLN
#' component and empirical zero rates / binomial GLMs for the ZI component.
#' Replaces the previous per-species `pscl::zeroinfl` loop.
#'
#' @param Y Response count matrix (n × p)
#' @param X Design matrix for the PLN component (n × d)
#' @param X0 Design matrix for the ZI component (n × d0, empty `matrix(NA,0,0)` when unused)
#' @param O Offset matrix in log-scale (n × p)
#' @param w Weight vector of length n (defaults to uniform weights)
#' @return Named list: `B` (d × p), `M` (n × p), `S2` (n × p), `R` (n × p), `B0` (d0 × p)
#'
#' @importFrom stats lm.fit glm.fit binomial
#' @export
compute_ZIPLN_starting_point <- function(Y, X, X0, O, w = NULL) {
  if (is.null(w)) w <- rep(1.0, nrow(Y))
  n <- nrow(Y); p <- ncol(Y); d0 <- ncol(X0)
  if (ncol(X) == 0) stop("PLN component requires at least one covariate (or an intercept) to fit ZIPLN.")

  ## PLN component: fast multivariate LM (identical to compute_PLN_starting_point "LM")
  sp <- compute_PLN_starting_point(Y, X, O, w, method = "LM")

  ## ZI component: empirical per-species zero rates
  zero_ind <- (Y == 0) * 1.0
  R <- matrix(colMeans(zero_ind), n, p, byrow = TRUE)

  ## B0: p binomial GLMs on zero indicator vs X0 (only for "covar" ziparam where d0 > 0)
  B0 <- if (!is.null(d0) && d0 > 0) {
    binom_fam <- binomial()
    vapply(seq_len(p), function(j)
      glm.fit(X0, zero_ind[, j], family = binom_fam)$coefficients,
      numeric(d0))
  } else {
    matrix(0.0, 0L, p)
  }

  list(B = sp$B, M = sp$M, S2 = sp$S2, R = R, B0 = B0)
}
