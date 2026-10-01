## Simulation with a known network: PLN data with a sparse precision matrix, in
## which a few species are then made structurally absent from a group of
## samples (as a species living on one type of tree only). The model is fitted
## without the covariate that would explain the absences.
source("proto.R")
args <- commandArgs(TRUE)
n_rep <- if (length(args)) as.integer(args[1]) else 20
n <- if (length(args) > 1) as.integer(args[2]) else 100; p <- 40

simulate <- function(seed, n_contaminated, zero_fraction = 0.65) {
  set.seed(seed)
  A <- matrix(0, p, p); A[upper.tri(A)] <- rbinom(p * (p - 1) / 2, 1, 0.05); A <- A + t(A)
  Omega <- 0.5 * A * matrix(sample(c(-1, 1), p * p, replace = TRUE), p, p)
  Omega[lower.tri(Omega)] <- t(Omega)[lower.tri(Omega)]
  diag(Omega) <- rowSums(abs(Omega)) + 0.1           # diagonally dominant, strong partial correlations
  Sigma <- cov2cor(solve(Omega))                     # unit latent variances, same support
  mu <- runif(p, 2, 4)
  Z <- matrix(rnorm(n * p), n, p) %*% chol(Sigma) + matrix(mu, n, p, byrow = TRUE)
  Y <- matrix(rpois(n * p, exp(Z)), n, p, dimnames = list(paste0("s", 1:n), paste0("sp", 1:p)))
  contaminated <- if (n_contaminated > 0) sample(p, n_contaminated) else integer(0)
  for (j in contaminated) Y[sample(n, round(zero_fraction * n)), j] <- 0
  list(Y = Y, truth = A == 1, contaminated = contaminated)
}

methods <- list(
  baseline          = function(f, d, X) fit_path(f, d, X, min_ratio = 0.05),
  "floor_1e-2"      = function(f, d, X) fit_path(f, d, X, floor = 1e-2, min_ratio = 0.05),
  "floor_1e-3"      = function(f, d, X) fit_path(f, d, X, floor = 1e-3, min_ratio = 0.05),
  exclusion         = function(f, d, X) fit_exclusion(f, d, X, min_ratio = 0.05),
  corr_penalty      = function(f, d, X) fit_path(f, d, X, corr = TRUE, min_ratio = 0.05),
  "corr+floor_1e-3" = function(f, d, X) fit_path(f, d, X, floor = 1e-3, corr = TRUE, min_ratio = 0.05))

evaluate <- function(nets, sim) {
  Y <- sim$Y; X <- matrix(1, n, 1); O <- matrix(0, n, p)
  clean <- setdiff(seq_len(p), sim$contaminated)
  up <- upper.tri(matrix(0, length(clean), length(clean)))
  truth <- sim$truth[clean, clean][up]
  per_model <- t(sapply(nets$models, function(m) {
    A <- as.matrix(m$model_par$Omega) != 0; diag(A) <- FALSE
    est <- A[clean, clean][up]
    ll <- sum(elbo_pln(Y, X, O, m$model_par$B, m$var_par$M, m$var_par$S2, m$model_par$Omega))
    c(tp = sum(est & truth), fp = sum(est & !truth), n_clean = sum(est), n_edges = sum(A) / 2,
      cont_edges = if (length(sim$contaminated)) sum(A) / 2 - sum(A[clean, clean]) / 2 else 0,
      bic = ll - .5 * log(n) * m$nb_param)
  }))
  prec <- per_model[, "tp"] / pmax(per_model[, "n_clean"], 1); rec <- per_model[, "tp"] / max(sum(truth), 1)
  f1 <- ifelse(prec + rec > 0, 2 * prec * rec / (prec + rec), 0)
  b <- which.max(per_model[, "bic"]); o <- which.max(f1)
  ## the model whose number of edges is closest to the true one: same size for all methods
  k <- which.min(abs(per_model[, "n_edges"] - sum(sim$truth) / 2))
  share <- function(i) unname(if (per_model[i, "n_edges"] > 0) per_model[i, "cont_edges"] / per_model[i, "n_edges"] else 0)
  c(maxF1 = max(f1), F1_BIC = unname(f1[b]), F1_size = unname(f1[k]), prec_size = unname(prec[k]),
    cont_share_size = share(k), cont_share_best = share(o), cont_share_BIC = share(b),
    max_Sigma = max(sapply(nets$models, function(m) max(diag(as.matrix(m$model_par$Sigma))))),
    degen_path = length(unique(unlist(lapply(nets$models, function(m) m$degenerate_species)))),
    excluded = length(attr(nets, "excluded")),
    excluded_clean = length(setdiff(attr(nets, "excluded"), colnames(Y)[sim$contaminated])))
}

res <- NULL
for (k in c(0, 3)) for (seed in seq_len(n_rep)) {
  sim <- simulate(seed, k)
  d <- suppressMessages(prepare_data(sim$Y, data.frame(id = seq_len(n), row.names = rownames(sim$Y)), offset = "none"))
  X <- matrix(1, n, 1, dimnames = list(NULL, "(Intercept)"))
  for (meth in names(methods)) {
    tt <- system.time(nets <- tryCatch(methods[[meth]](Abundance ~ 1, d, X),
                                       error = function(e) { cat("  ERROR", meth, ":", conditionMessage(e), "\n"); NULL }))[3]
    if (is.null(nets)) next
    res <- rbind(res, data.frame(contaminated = k, seed = seed, method = meth, t(evaluate(nets, sim)), time = unname(tt)))
  }
  saveRDS(res, paste0("sim_known_network_n", n, ".rds"))
  cat("done k =", k, "seed", seed, "\n")
}
agg <- aggregate(cbind(maxF1, F1_size, prec_size, F1_BIC, cont_share_size, cont_share_best, cont_share_BIC, max_Sigma, degen_path, excluded, excluded_clean, time) ~ contaminated + method,
                 data = res, FUN = function(x) signif(mean(x), 2))
cat("true share of the edges touching a contaminated species, on average:", round(2 * 3 / p, 2), "\n")
options(width = 250); print(agg[order(agg$contaminated), ], row.names = FALSE)
