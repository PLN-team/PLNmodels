## Comparison of the variants of the floor (floor_variants.R) on simulated data
## with a known network, all on the correlation scale of the penalty.
##   Rscript sim_floor_variants.R <regime> <n_rep>
## Regimes (p species, n samples, k contaminated species):
##   random : each contaminated species is absent from its own random 65 % of
##            the samples (the regime of sim_known_network.R)
##   groups : the samples belong to 3 groups, as the trees of oaks, and each
##            contaminated species lives in one group only, where it is abundant
##   sparse : as groups, but the other species have low abundances, hence many
##            zeros that are not structural (the situation of mollusk)
##   overdispersed : as sparse, with latent variances of 6 instead of 1, so that
##            the zeros of the uncontaminated species have very negative latent
##            means too
source("floor_variants.R")
args <- commandArgs(TRUE)
regime <- if (length(args)) args[1] else "groups"
n_rep  <- if (length(args) > 1) as.integer(args[2]) else 20
config <- switch(regime,
  random = list(n = 50, p = 40, k = 3, shared = FALSE, boost = 0),
  groups = list(n = 45, p = 40, k = 8, shared = TRUE,  boost = 1.5),
  sparse = list(n = 45, p = 40, k = 8, shared = TRUE,  boost = 3, mu = c(-0.5, 2)),
  overdispersed = list(n = 100, p = 40, k = 8, shared = TRUE, boost = 3, mu = c(-1, 2), variance = 6))
if (is.null(config$variance)) config$variance <- 1
if (is.null(config$mu)) config$mu <- c(2, 4)
if (length(args) > 2) config$n <- as.integer(args[3])
n <- config$n; p <- config$p

simulate <- function(seed, k) {
  set.seed(seed)
  A <- matrix(0, p, p); A[upper.tri(A)] <- rbinom(p * (p - 1) / 2, 1, 0.05); A <- A + t(A)
  Omega <- 0.5 * A * matrix(sample(c(-1, 1), p * p, replace = TRUE), p, p)
  Omega[lower.tri(Omega)] <- t(Omega)[lower.tri(Omega)]
  diag(Omega) <- rowSums(abs(Omega)) + 0.1
  Sigma <- config$variance * cov2cor(solve(Omega))
  mu <- runif(p, config$mu[1], config$mu[2])
  contaminated <- if (k > 0) sample(p, k) else integer(0)
  mu[contaminated] <- mu[contaminated] + config$boost
  Z <- matrix(rnorm(n * p), n, p) %*% chol(Sigma) + matrix(mu, n, p, byrow = TRUE)
  Y <- matrix(rpois(n * p, exp(Z)), n, p, dimnames = list(paste0("s", 1:n), paste0("sp", 1:p)))
  group <- sample(rep(1:3, length.out = n))
  for (j in contaminated) {
    absent <- if (config$shared) group != sample(3, 1) else seq_len(n) %in% sample(n, round(0.65 * n))
    Y[absent, j] <- 0
  }
  list(Y = Y, truth = A == 1, contaminated = contaminated)
}

evaluate <- function(nets, sim) {
  Y <- sim$Y; X <- matrix(1, n, 1); O <- matrix(0, n, p)
  clean <- setdiff(seq_len(p), sim$contaminated)
  up <- upper.tri(matrix(0, length(clean), length(clean)))
  truth <- sim$truth[clean, clean][up]
  per_model <- t(sapply(nets$models, function(m) {
    A <- as.matrix(m$model_par$Omega) != 0; diag(A) <- FALSE
    est <- A[clean, clean][up]
    c(tp = sum(est & truth), n_clean = sum(est), n_edges = sum(A) / 2,
      cont_edges = sum(A) / 2 - sum(A[clean, clean]) / 2,
      loglik = sum(elbo_pln(Y, X, O, m$model_par$B, m$var_par$M, m$var_par$S2, m$model_par$Omega)),
      at_floor = if (is.null(m$optim_par$n_floor)) 0 else m$optim_par$n_floor / length(Y),
      max_degree = max(rowSums(A)))
  }))
  prec <- per_model[, "tp"] / pmax(per_model[, "n_clean"], 1); rec <- per_model[, "tp"] / max(sum(truth), 1)
  f1 <- ifelse(prec + rec > 0, 2 * prec * rec / (prec + rec), 0)
  ## the model whose number of edges is closest to the true one: same size for all variants
  s <- which.min(abs(per_model[, "n_edges"] - sum(sim$truth) / 2))
  c(maxF1 = max(f1), F1_size = unname(f1[s]),
    cont_share_size = unname(if (per_model[s, "n_edges"] > 0) per_model[s, "cont_edges"] / per_model[s, "n_edges"] else 0),
    max_degree_size = unname(per_model[s, "max_degree"]),
    loglik_size = unname(per_model[s, "loglik"]),
    at_floor = mean(per_model[, "at_floor"]), zeros = mean(Y == 0),
    max_Sigma = max(sapply(nets$models, function(m) max(diag(as.matrix(m$model_par$Sigma))))),
    degenerate = length(unique(unlist(lapply(nets$models, function(m) m$degenerate_species)))),
    floored_species = attr(nets, "floored_species"))
}

res <- NULL
for (k in c(0, config$k)) for (seed in seq_len(n_rep)) {
  sim <- simulate(seed, k)
  d <- suppressMessages(prepare_data(sim$Y, data.frame(id = seq_len(n), row.names = rownames(sim$Y)), offset = "none"))
  for (v in names(variants)) {
    tt <- system.time(nets <- tryCatch(variants[[v]](Abundance ~ 1, d),
                                       error = function(e) { cat("  ERROR", v, ":", conditionMessage(e), "\n"); NULL }))[3]
    if (is.null(nets)) next
    res <- rbind(res, data.frame(contaminated = k, seed = seed, variant = v, t(evaluate(nets, sim)), time = unname(tt)))
  }
  saveRDS(res, paste0("sim_floor_variants_", regime, "_n", n, ".rds"))
  cat("done k =", k, "seed", seed, "\n")
}
agg <- aggregate(cbind(maxF1, F1_size, cont_share_size, max_degree_size, loglik_size, at_floor, zeros, max_Sigma, degenerate, time) ~ contaminated + variant,
                 data = res, FUN = function(x) signif(mean(x), 3))
agg$diverged <- aggregate(max_Sigma ~ contaminated + variant, data = res, FUN = function(x) sum(x > 100))$max_Sigma
cat("expected share of the edges touching a contaminated species:", round(2 * config$k / p, 2), "\n")
options(width = 220); print(agg[order(agg$contaminated, match(agg$variant, names(variants))), ], row.names = FALSE)
