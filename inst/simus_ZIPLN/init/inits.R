## Candidate starting points for ZIPLN, with the signature of compute_ZIPLN_starting_point()
## (Y, X, X0, O, w) -> list(B, M, S2, R, B0)

## Gauss-Hermite nodes and weights for E[f(N(0, 1))] (Golub-Welsch)
gauss_hermite <- function(k = 20) {
  J <- matrix(0, k, k); off <- sqrt(seq_len(k - 1) / 2)
  J[cbind(1:(k - 1), 2:k)] <- off; J[cbind(2:k, 1:(k - 1))] <- off
  e <- eigen(J, symmetric = TRUE)
  list(x = sqrt(2) * e$values, w = e$vectors[1, ]^2)
}
GH <- gauss_hermite()

## P(Y = 0) under a Poisson-lognormal with log-mean eta and latent variance s2 (matrices)
p0_lognormal <- function(eta, s2) {
  out <- 0 * eta
  for (k in seq_along(GH$x)) out <- out + GH$w[k] * exp(-exp(pmin(eta + sqrt(s2) * GH$x[k], 700)))
  out
}

## C0: current starting point of the package (captured before any mocking)
init_current <- get("compute_ZIPLN_starting_point", asNamespace("PLNmodels"))

## C1: B and the latent variances from the positive counts only (as a hurdle model),
## M = XB on the zeros, then the probability of a structural zero under a
## Poisson-lognormal with these parameters, the ZI rate pi_j being matched on the
## observed rate of zeros of each species
init_positive <- function(Y, X, X0, O, w = NULL, pi_max = 0.95) {
  if (is.null(w)) w <- rep(1, nrow(Y))
  n <- nrow(Y); p <- ncol(Y); d0 <- ncol(X0)
  base <- init_current(Y, X, X0, O, w)
  B <- base$B; s2 <- rep(1, p)
  logY <- log(pmax(Y, 1)) - O
  for (j in seq_len(p)) {
    pos <- Y[, j] > 0
    if (sum(pos) <= ncol(X) + 1) next
    fit <- lm.wfit(X[pos, , drop = FALSE], logY[pos, j], w[pos])
    if (any(is.na(fit$coefficients))) next
    B[, j] <- fit$coefficients
    s2[j] <- max(sum(w[pos] * fit$residuals^2) / sum(w[pos]), 0.05)
  }
  XB <- X %*% B
  p0 <- p0_lognormal(O + XB, matrix(s2, n, p, byrow = TRUE))
  zero <- Y == 0
  pi <- pmin(pmax((colMeans(zero) - colMeans(p0)) / (1 - colMeans(p0)), 0), pi_max)
  pim <- matrix(pi, n, p, byrow = TRUE)
  R <- zero * pim / (pim + (1 - pim) * p0)
  M <- ifelse(zero, XB, logY)
  S2 <- ifelse(zero, matrix(s2, n, p, byrow = TRUE), 1 / (2 + Y))
  list(B = B, M = M, S2 = S2, R = R, B0 = base$B0)
}

## C2: univariate zero-inflated Poisson by EM (as pscl::zeroinfl did, without overdispersion)
init_zip <- function(Y, X, X0, O, w = NULL, iter = 20) {
  if (is.null(w)) w <- rep(1, nrow(Y))
  n <- nrow(Y); p <- ncol(Y)
  base <- init_current(Y, X, X0, O, w)
  B <- base$B; zero <- Y == 0
  pi <- rep(.1, p); R <- zero * 0.1
  for (t in seq_len(iter)) {
    lambda <- exp(O + X %*% B)
    pim <- matrix(pi, n, p, byrow = TRUE)
    R <- zero * pim / (pim + (1 - pim) * exp(-lambda))
    for (j in seq_len(p)) {
      fit <- suppressWarnings(glm.fit(X, Y[, j], weights = w * (1 - R[, j]), offset = O[, j],
                                      family = poisson(), start = B[, j], control = glm.control(maxit = 5)))
      if (!any(is.na(fit$coefficients))) B[, j] <- fit$coefficients
    }
    pi <- pmin(colMeans(R), .95)
  }
  XB <- X %*% B
  list(B = B, M = ifelse(zero, XB, log(pmax(Y, 1)) - O), S2 = 1 / (2 + Y), R = R, B0 = base$B0)
}

## Starting point from the true parameters of a simulation
init_truth <- function(sim) {
  force(sim)
  function(Y, X, X0, O, w = NULL) {
    XB <- X %*% sim$B
    list(B = sim$B, M = ifelse(Y > 0, log(pmax(Y, 1)) - O, XB), S2 = matrix(.3, nrow(Y), ncol(Y)),
         R = (Y == 0) * pmin(pmax(sim$W, .02), .98), B0 = matrix(0, 0, ncol(Y)))
  }
}

## Fit ZIPLN() with a given starting point
zipln_with_init <- function(init, ...) {
  ns <- asNamespace("PLNmodels")
  old <- get("compute_ZIPLN_starting_point", ns)
  unlockBinding("compute_ZIPLN_starting_point", ns)
  assign("compute_ZIPLN_starting_point", init, ns)
  on.exit({ assign("compute_ZIPLN_starting_point", old, ns); lockBinding("compute_ZIPLN_starting_point", ns) })
  suppressWarnings(suppressMessages(ZIPLN(...)))
}
