library(PLNmodels)
library(MASS)
library(tidyverse)
library(parallel)

rZIPLN <- function(n     = 10,
                   mu    = rep(0, ncol(Sigma)),
                   Sigma = diag(1, 5, 5),
                   Pi    = matrix(1, n,, ncol(Sigma)),
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
  Z <- mu + mvrnorm(n, rep(0, ncol(Sigma)), as.matrix(Sigma)) + offsets
  W <- matrix(rbinom(n * p, 1,  prob = Pi), n, p)
  Y <- matrix(rpois(n * p, as.vector(exp(Z))), n, p) * (1 - W)
  dimnames(Y) <- list(paste0("S", 1:n), paste0("Y", 1:p))
  Y
}

### SIMULATED DATA
p <- 5
logit    <- function(x) log(x/(1-x))
logistic <- function(x) 1 / (1 + exp(-x))

### Simulation parameters
## Sigma: AR/toeplitz structure
sd <- sqrt(sample.int(5, p, replace = TRUE))
Sigma_star <- diag(sd) %*% toeplitz(0.5^(1:p - 1)) %*% diag(sd)
Omega_star <- solve(Sigma_star)
## PLN part: one intercept/grand mean with value equal to 1
B_star  <- matrix(2, 1, p)
## ZI part: proba of 0.1 for being a 0 from the ZI component
B0_star <- matrix(logit(0.2), 1, p)

vec_n <- c(25, 50, 75, 100)
one_simu <- function(i) {
  cat(i)
  n <- max(vec_n)*2
  # X <- cbind(rep(1, n))
  X <- rnorm(n, 2, 1)
  mu_star <- cbind(X) %*% B_star
  Pi_star <- logistic(cbind(X) %*% B0_star)
  Y <- rZIPLN(n, mu = mu_star, Sigma = Sigma_star, Pi = Pi_star)
  Y <- Y[rowSums(Y) > 0, ]
  X <- X[rowSums(Y) > 0]
  err_ZIPLN <- sapply(vec_n,  function(n_) {
    Y_ <- Y[1:n_, ]; X_ <- X[1:n_]
    myZIPLN <- ZIPLN(Y_ ~ 0 + X_ + offset(log(rowSums(Y_))) | X_ , control = ZIPLN_param(trace = 0))
    c(pred_Y = sqrt(mean((myZIPLN$fitted - Y_)^2)),
      rmse_B = sqrt(mean((myZIPLN$model_par$B[2, ] - B_star)^2)),
      rmse_Omega = sqrt(mean((myZIPLN$model_par$Omega - Omega_star)^2)))
  })
  err_PLN <- sapply(vec_n,  function(n_) {
    Y_ <- Y[1:n_, ]; X_ <- X[1:n_]
    myPLN <- PLN(Y_ ~ 1 + X_ + offset(log(rowSums(Y_))) , control = PLN_param(trace = 0))
    c(pred = sqrt(mean((myPLN$fitted - Y_)^2)),
      rmse_B = sqrt(mean((myPLN$model_par$B[2, ] - B_star)^2)),
      rmse_Sigma = sqrt(mean((myPLN$model_par$Omega - Omega_star)^2))
    )
  })
  data.frame(rbind(t(err_ZIPLN), t(err_PLN)),
             method = rep(c("ZIPLN","PLN"), each = length(vec_n)),
             n = c(vec_n, vec_n), simu = rep(i, 2*length(vec_n)))
}

res <- do.call(rbind, lapply(1:50, one_simu))

p <- ggplot(res) + aes(x = factor(n), y = pred_Y, fill = factor(method)) + geom_violin() + theme_bw() +
  scale_y_log10() + ylim(c(0,2))
p

p <- ggplot(res) + aes(x = factor(n), y = rmse_B, fill = factor(method)) + geom_violin() + theme_bw() + scale_y_log10() + ylim(c(0,3))
p

p <- ggplot(res) + aes(x = factor(n), y = rmse_Omega, fill = factor(method)) + geom_violin() + theme_bw() + scale_y_log10() + ylim(c(0,0.5))
p
