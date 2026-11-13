## Cheap indicators, computed before any fit, to predict which starting point of ZIPLN wins,
## and the cost of fitting the predicted winner first (the other one being then abandoned
## by cannot_beat()). Uses the trajectories recorded by trajectories.R.
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("inst/simus_ZIPLN/init/simulate.R")

## Indicators from the data and the two starting points
indicators <- function(Y, X, X0, O) {
  zero <- Y == 0
  sp <- compute_ZIPLN_starting_point(Y, X, X0, O, method = "positive")
  XB <- X %*% sp$B
  s2 <- sp$S2[cbind(apply(zero, 2, function(z) which(z)[1]), seq_len(ncol(Y)))]   # s2 of each species, on its zeros
  s2[is.na(s2)] <- 1
  p0 <- p0_poisson_lognormal(O + XB, matrix(s2, nrow(Y), ncol(Y), byrow = TRUE))
  excess <- pmax(colMeans(zero) - colMeans(p0), 0)                       # zeros not explained by the PLN of the positive counts
  data.frame(n = nrow(Y), p = ncol(Y), zeros = mean(zero),
             excess = mean(excess),                                       # mean over species
             excess_share = sum(excess) / max(sum(colMeans(zero)), 1e-12),   # share of the zeros
             p0 = mean(p0),                                               # zeros expected without zero inflation
             log_mean = median(O + XB),                                   # typical log expected count
             species_excess = mean(excess > 0.1))                         # species with clearly excess zeros
}

traj_sim  <- readRDS("inst/simus_ZIPLN/init/trajectories_sim.rds")
traj_real <- readRDS("inst/simus_ZIPLN/init/trajectories_real.rds")
tasks <- expand.grid(rep = 1:10, cfg = seq_len(nrow(configs)))
rows <- list()
for (k in seq_len(nrow(tasks))) {
  cfg <- configs[tasks$cfg[k], ]; sim <- simulate(cfg, (tasks$cfg[k] - 1) * reps + tasks$rep[k])
  Y <- sim$Y; tr <- traj_sim[[k]]
  rows[[length(rows) + 1]] <- data.frame(case = tr$case, group = tr$group,
    indicators(Y, matrix(1, nrow(Y), 1), matrix(NA, 0, 0), matrix(0, nrow(Y), ncol(Y))),
    loglik_lm = tr$lm$loglik, loglik_pos = tr$positive$loglik)
  attr(rows[[length(rows)]], "traj") <- list(lm = tr$lm$objective, positive = tr$positive$objective)
}
data(trichoptera); data(barents); data(oaks); data(mollusk); data(microcosm)
tri <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
mol <- suppressMessages(prepare_data(mollusk$Abundance, mollusk$Covariate))
cases <- list(
  list("trichoptera single", tri, Abundance ~ 1 + offset(log(Offset))),
  list("trichoptera col", tri, Abundance ~ 1 + offset(log(Offset))),
  list("barents single", barents, Abundance ~ 1 + offset(log(Offset))),
  list("barents col", barents, Abundance ~ 1 + offset(log(Offset))),
  list("mollusk single", mol, Abundance ~ 1 + offset(log(Offset))),
  list("mollusk col", mol, Abundance ~ 1 + offset(log(Offset))),
  list("mollusk site|site", mol, Abundance ~ 0 + site + offset(log(Offset)) | 0 + site),
  list("oaks single", oaks, Abundance ~ 1 + offset(log(Offset))),
  list("oaks col", oaks, Abundance ~ 1 + offset(log(Offset))),
  list("oaks tree col", oaks, Abundance ~ 0 + tree + offset(log(Offset))),
  list("oaks tree|tree", oaks, Abundance ~ 0 + tree + offset(log(Offset)) | 0 + tree),
  list("microcosm site col", microcosm, Abundance ~ 0 + site + offset(log(Offset))),
  list("microcosm site|site", microcosm, Abundance ~ 0 + site + offset(log(Offset)) | 0 + site))
for (cs in cases) {
  dz <- extract_model_zi(call("ZIPLN", formula = cs[[3]], data = quote(dd)), list2env(list(dd = cs[[2]])))
  tr <- Filter(function(x) x$case == cs[[1]], traj_real); tr <- setNames(lapply(tr, `[[`, "fit"), vapply(tr, `[[`, "", "init"))
  rows[[length(rows) + 1]] <- data.frame(case = cs[[1]], group = "real", indicators(dz$Y, dz$X, dz$X0, dz$O),
                                         loglik_lm = tr$lm$loglik, loglik_pos = tr$positive$loglik)
  attr(rows[[length(rows)]], "traj") <- list(lm = tr$lm$objective, positive = tr$positive$objective)
}
saveRDS(rows, "inst/simus_ZIPLN/init/order_features.rds")
feat <- do.call(rbind, rows)
feat$winner <- ifelse(feat$loglik_pos > feat$loglik_lm, "positive", "lm")
feat$gap <- abs(feat$loglik_pos - feat$loglik_lm)
options(width = 200)
print(aggregate(cbind(zeros, excess, excess_share, p0, log_mean, species_excess, positive_wins = winner == "positive") ~ group, feat, mean), digits = 3)
print(feat[feat$group == "real", c("case", "zeros", "excess", "excess_share", "p0", "log_mean", "species_excess", "winner", "gap")], digits = 3, row.names = FALSE)
