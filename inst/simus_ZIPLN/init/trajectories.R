## Trajectories of the objective from each starting point: after how many iterations is the
## ordering of the final bounds known? Usage: Rscript trajectories.R sim|real
what <- commandArgs(TRUE)[1]
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("inst/simus_ZIPLN/init/simulate.R")
n_cores <- as.integer(Sys.getenv("N_CORES", "12"))
traj <- function(formula, data, zi, init) {
  t <- system.time(f <- suppressWarnings(suppressMessages(
    ZIPLN(formula, data = data, zi = zi, control = ZIPLN_param(trace = 0, init = init)))))[["elapsed"]]
  list(objective = f$optim_par$objective, loglik = f$loglik, time = t, iterations = f$optim_par$iterations)
}
if (what == "sim") {
  tasks <- expand.grid(rep = 1:10, cfg = seq_len(nrow(configs)))
  one <- function(k) {
    RhpcBLASctl::blas_set_num_threads(1)
    cfg <- configs[tasks$cfg[k], ]; sim <- simulate(cfg, (tasks$cfg[k] - 1) * reps + tasks$rep[k])
    dat <- data.frame(row.names = seq_len(nrow(sim$Y))); dat$Y <- sim$Y
    list(case = paste(cfg$abundance, cfg$zi, cfg$rho, tasks$rep[k]), group = paste(cfg$abundance, cfg$zi),
         lm = traj(Y ~ 1, dat, "col", "lm"), positive = traj(Y ~ 1, dat, "col", "positive"))
  }
  todo <- seq_len(nrow(tasks))
} else {
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
    list("microcosm site|1", microcosm, Abundance ~ 0 + site + offset(log(Offset)) | 1, "covar"),
    list("microcosm site|site", microcosm, Abundance ~ 0 + site + offset(log(Offset)) | 0 + site, "covar"))
  ## one task per (case, start), the heavy microcosm fits first
  tasks <- expand.grid(init = c("lm", "positive"), case = seq_along(cases), stringsAsFactors = FALSE)
  one <- function(k) {
    RhpcBLASctl::blas_set_num_threads(1)
    cs <- cases[[tasks$case[k]]]
    list(case = cs[[1]], init = tasks$init[k], fit = traj(cs[[3]], cs[[2]], cs[[4]], tasks$init[k]))
  }
  todo <- order(-(tasks$case >= 12))
}
out <- parallel::mclapply(todo, function(k) tryCatch(one(k), error = function(e) conditionMessage(e)),
                          mc.cores = n_cores, mc.preschedule = FALSE)
err <- Filter(is.character, out); if (length(err)) print(unique(unlist(err)))
saveRDS(Filter(is.list, out), sprintf("inst/simus_ZIPLN/init/trajectories_%s.rds", what))
cat("done", what, "\n")
