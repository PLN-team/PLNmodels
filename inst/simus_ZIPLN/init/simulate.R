## ZIPLN simulations (same design as inst/simus_ZIPLN/variational_J3 on branch zipln-j3)
configs <- expand.grid(rho = c(0, .5, .8), zi = c("low", "high"), abundance = c("low", "high"),
                       n = 100, p = 20, stringsAsFactors = FALSE)
reps <- 20
simulate <- function(cfg, seed) {
  set.seed(seed)
  n <- cfg$n; p <- cfg$p
  Sigma <- 0.8 * cfg$rho^abs(outer(1:p, 1:p, "-"))   # 0^0 = 1: diagonal when rho = 0
  B  <- matrix(rnorm(p, if (cfg$abundance == "low") 0.5 else 2, 0.5), 1)
  pi <- if (cfg$zi == "low") runif(p, .05, .2) else runif(p, .3, .6)
  Z  <- matrix(B, n, p, byrow = TRUE) + MASS::mvrnorm(n, rep(0, p), Sigma)
  W  <- matrix(rbinom(n * p, 1, rep(pi, each = n)), n, p)
  Y  <- (1 - W) * matrix(rpois(n * p, exp(Z)), n, p)
  colnames(Y) <- paste0("sp", 1:p)
  list(Y = Y, Sigma = Sigma, B = B, pi = pi, W = W)
}
