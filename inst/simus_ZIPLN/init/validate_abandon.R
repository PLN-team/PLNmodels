## init = "best" with the abandonment of the second start, against the two full fits
## recorded in trajectories_*.rds: same bound? how much faster?
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("inst/simus_ZIPLN/init/simulate.R")
n_cores <- as.integer(Sys.getenv("N_CORES", "12"))
fit_best <- function(formula, data, zi) {
  t <- system.time(f <- suppressWarnings(suppressMessages(ZIPLN(formula, data = data, zi = zi, control = ZIPLN_param(trace = 0)))))[["elapsed"]]
  data.frame(loglik_best = f$loglik, time_best = t, abandoned = length(f$optim_par$abandoned_starts) > 0, kept = f$optim_par$init)
}
traj_sim <- readRDS("inst/simus_ZIPLN/init/trajectories_sim.rds")
tasks <- expand.grid(rep = 1:10, cfg = seq_len(nrow(configs)))
one_sim <- function(k) {
  RhpcBLASctl::blas_set_num_threads(1)
  cfg <- configs[tasks$cfg[k], ]; sim <- simulate(cfg, (tasks$cfg[k] - 1) * reps + tasks$rep[k])
  dat <- data.frame(row.names = seq_len(nrow(sim$Y))); dat$Y <- sim$Y; tr <- traj_sim[[k]]
  data.frame(case = tr$case, group = tr$group, fit_best(Y ~ 1, dat, "col"),
             loglik_lm = tr$lm$loglik, loglik_pos = tr$positive$loglik, time_both = tr$lm$time + tr$positive$time, time_lm = tr$lm$time)
}
data(trichoptera); data(barents); data(oaks); data(mollusk); data(microcosm)
tri <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
mol <- suppressMessages(prepare_data(mollusk$Abundance, mollusk$Covariate))
cases <- list(
  list("trichoptera single", tri, Abundance ~ 1 + offset(log(Offset)), "single"),
  list("trichoptera col", tri, Abundance ~ 1 + offset(log(Offset)), "col"),
  list("barents single", barents, Abundance ~ 1 + offset(log(Offset)), "single"),
  list("barents col", barents, Abundance ~ 1 + offset(log(Offset)), "col"),
  list("mollusk single", mol, Abundance ~ 1 + offset(log(Offset)), "single"),
  list("mollusk col", mol, Abundance ~ 1 + offset(log(Offset)), "col"),
  list("mollusk site|site", mol, Abundance ~ 0 + site + offset(log(Offset)) | 0 + site, "covar"),
  list("oaks single", oaks, Abundance ~ 1 + offset(log(Offset)), "single"),
  list("oaks col", oaks, Abundance ~ 1 + offset(log(Offset)), "col"),
  list("oaks tree col", oaks, Abundance ~ 0 + tree + offset(log(Offset)), "col"),
  list("oaks tree|tree", oaks, Abundance ~ 0 + tree + offset(log(Offset)) | 0 + tree, "covar"),
  list("microcosm site col", microcosm, Abundance ~ 0 + site + offset(log(Offset)), "col"),
  list("microcosm site|site", microcosm, Abundance ~ 0 + site + offset(log(Offset)) | 0 + site, "covar"))
traj_real <- readRDS("inst/simus_ZIPLN/init/trajectories_real.rds")
one_real <- function(k) {
  RhpcBLASctl::blas_set_num_threads(1)
  cs <- cases[[k]]
  tr <- Filter(function(x) x$case == cs[[1]], traj_real); tr <- setNames(lapply(tr, `[[`, "fit"), vapply(tr, `[[`, "", "init"))
  data.frame(case = cs[[1]], group = "real", fit_best(cs[[3]], cs[[2]], cs[[4]]),
             loglik_lm = tr$lm$loglik, loglik_pos = tr$positive$loglik, time_both = tr$lm$time + tr$positive$time, time_lm = tr$lm$time)
}
jobs <- c(lapply(c(12, 13, 1:11), function(k) function() one_real(k)), lapply(seq_len(nrow(tasks)), function(k) function() one_sim(k)))
out <- parallel::mclapply(jobs, function(job) tryCatch(job(), error = function(e) conditionMessage(e)), mc.cores = n_cores, mc.preschedule = FALSE)
err <- Filter(is.character, out); if (length(err)) print(unique(unlist(err)))
res <- do.call(rbind, Filter(is.data.frame, out))
saveRDS(res, "inst/simus_ZIPLN/init/validate_abandon.rds")
res$best2 <- pmax(res$loglik_lm, res$loglik_pos)
options(width = 200)
print(do.call(rbind, lapply(split(res, res$group), function(d) data.frame(group = d$group[1], n = nrow(d),
  abandoned = sum(d$abandoned), max_shortfall = round(max(d$best2 - d$loglik_best), 2), worse_than_best2 = sum(d$loglik_best < d$best2 - 1e-3),
  time_both = round(sum(d$time_both), 1), time_best = round(sum(d$time_best), 1), time_lm = round(sum(d$time_lm), 1)))), row.names = FALSE)
print(res[res$group == "real", c("case", "kept", "abandoned", "loglik_best", "best2", "time_both", "time_best", "time_lm")], row.names = FALSE, digits = 7)
