available_algorithms_nlopt <- c("MMA", "CCSAQ", "LBFGS", "LBFGS_NOCEDAL", "VAR1", "VAR2") #"TNEWTON", "TNEWTON_PRECOND", "TNEWTON_PRECOND_RESTART"#
available_algorithms_torch <- c("RPROP", "RMSPROP", "ADAM", "ADAGRAD")

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

status_to_message <- function(status) {
  message <- switch(as.character(status),
                    "1"  = "success",
                    "2"  = "success, stopval was reached",
                    "3"  = "success, ftol_rel or ftol_abs was reached",
                    "4"  = "success, xtol_rel or xtol_abs was reached",
                    "5"  = "success, maxeval was reached",
                    "6"  = "success, maxtime was reached",
                    "-1" = "failure",
                    "-2" = "invalid arguments",
                    "-3" = "out of memory.",
                    "-4" = "roundoff errors led to a breakdown of the optimization algorithm",
                    "-5" = "forced termination:",
                    "Return status not recognized"
  )
  message
}

trace <- function(x) sum(diag(x))

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
  B <- do.call(cbind, future_lapply(1:ncol(responses), function(j)
    coefficients(suppressWarnings(
      glm.fit(covariates, responses[, j], weights = weights, offset = offsets[, j], family = stats::poisson(),
        control = glm.control(epsilon = 1e-3, maxit = 10))))))
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
  ## Save encountered levels for predict methods as attribute of the formula
  attr(call$formula, "xlevels") <- .getXlevels(terms(frame), frame)
  list(Y = Y, X = X, O = O, miss = is.na(Y), w = w, formula = call$formula)
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
#' Barebone function to compute starting points for B, M and S when fitting a PLN. Mostly intended for internal use.
#'
#' @param Y Response count matrix
#' @param X Covariate matrix
#' @param O Offset matrix (in log-scale)
#' @param w Weight vector (defaults to 1)
#' @param s Scale parameter for S (defaults to 0.1)
#' @return a named list of starting values for model parameter B and variational parameters M and S used in the iterative optimization algorithm of [PLN()]
#'
#' @details The default strategy to estimate B and M is to fit a linear model with covariates `X` to the response count matrix (after adding a pseudocount of 1, scaling by the offset and taking the log). The regression matrix is used to initialize `B` and the residuals to initialize `M`. `S` is initialized as a constant conformable matrix with value `s`.
#'
#' @rdname compute_PLN_starting_point
#' @examples
#' \dontrun{
#' data(barents)
#' Y <- barents$Abundance
#' X <- model.matrix(Abundance ~ Latitude + Longitude + Depth + Temperature, data = barents)
#' O <- log(barents$Offset)
#' w <-- rep(1, nrow(Y))
#' compute_PLN_starting_point(Y, X, O, w)
#' }
#'
#' @importFrom stats lm.fit
#' @export
compute_PLN_starting_point <- function(Y, X, O, w, s = 0.1) {
  # Y = responses, X = covariates, O = offsets (in log scale), w = weights
  n <- nrow(Y); p <- ncol(Y); d <- ncol(X)
  fits <- lm.fit(w * X, w * log((1 + Y)/exp(O)))
  list(B = matrix(fits$coefficients, d, p),
       M = matrix(fits$residuals, n, p),
       S = matrix(s, n, p))
}
