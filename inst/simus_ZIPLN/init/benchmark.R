## Starting points of ZIPLN: final ELBO, estimation errors (simulations) and time
## Usage: Rscript benchmark.R sim|real
what <- commandArgs(TRUE)[1]
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("inst/simus_ZIPLN/init/simulate.R"); source("inst/simus_ZIPLN/init/inits.R")
inits <- list(current = init_current, positive = init_positive, zip = init_zip)
n_cores <- as.integer(Sys.getenv("N_CORES", "12"))
summ <- function(f, t) data.frame(loglik = f$loglik, time = t, iterations = f$optim_par$iterations,
                                  mean_pi = mean(f$model_par$Pi), mean_var = mean(diag(as.matrix(f$model_par$Sigma))))
if (what == "sim") {
  tasks <- expand.grid(rep = 1:10, cfg = seq_len(nrow(configs)))
  one <- function(k) {
    RhpcBLASctl::blas_set_num_threads(1)
    cfg <- configs[tasks$cfg[k], ]; seed <- (tasks$cfg[k] - 1) * reps + tasks$rep[k]; sim <- simulate(cfg, seed)
    dat <- data.frame(row.names = seq_len(nrow(sim$Y))); dat$Y <- sim$Y
    do.call(rbind, lapply(names(c(inits, truth = 1)), function(nm) {
      init <- if (nm == "truth") init_truth(sim) else inits[[nm]]
      t <- system.time(f <- zipln_with_init(init, Y ~ 1, data = dat, zi = "col", control = ZIPLN_param(trace = 0)))[["elapsed"]]
      S <- as.matrix(f$model_par$Sigma); up <- upper.tri(S)
      data.frame(cfg[1:3], rep = tasks$rep[k], init = nm, summ(f, t),
                 err_pi = sqrt(mean((f$model_par$Pi[1, ] - sim$pi)^2)),
                 err_var = mean(abs(log(diag(S) / diag(sim$Sigma)))),
                 err_cor = sqrt(mean((cov2cor(S)[up] - cov2cor(sim$Sigma)[up])^2)),
                 err_B = sqrt(mean((f$model_par$B - sim$B)^2)))
    }))
  }
  out <- parallel::mclapply(seq_len(nrow(tasks)), function(k) tryCatch(one(k), error = function(e) conditionMessage(e)),
                            mc.cores = n_cores, mc.preschedule = FALSE)
} else {
  data(trichoptera); data(barents); data(oaks); data(mollusk); data(microcosm)
  tri <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
  mol <- suppressMessages(prepare_data(mollusk$Abundance, mollusk$Covariate))
  cases <- list(
    list("trichoptera", tri, Abundance ~ 1 + offset(log(Offset)), "single"),
    list("trichoptera", tri, Abundance ~ 1 + offset(log(Offset)), "col"),
    list("barents", barents, Abundance ~ 1 + offset(log(Offset)), "single"),
    list("barents", barents, Abundance ~ 1 + offset(log(Offset)), "col"),
    list("mollusk", mol, Abundance ~ 1 + offset(log(Offset)), "single"),
    list("mollusk", mol, Abundance ~ 1 + offset(log(Offset)), "col"),
    list("mollusk", mol, Abundance ~ 0 + site + offset(log(Offset)) | 0 + site, "covar"),
    list("oaks", oaks, Abundance ~ 1 + offset(log(Offset)), "single"),
    list("oaks", oaks, Abundance ~ 1 + offset(log(Offset)), "col"),
    list("oaks", oaks, Abundance ~ 0 + tree + offset(log(Offset)), "col"),
    list("oaks", oaks, Abundance ~ 0 + tree + offset(log(Offset)) | 0 + tree, "covar"),
    list("microcosm", microcosm, Abundance ~ 0 + site + offset(log(Offset)), "col"),
    list("microcosm", microcosm, Abundance ~ 0 + site + offset(log(Offset)) | 1, "covar"),
    list("microcosm", microcosm, Abundance ~ 0 + site + offset(log(Offset)) | 0 + site, "covar"))
  tasks <- expand.grid(init = names(inits), case = seq_along(cases), stringsAsFactors = FALSE)
  one <- function(k) {
    RhpcBLASctl::blas_set_num_threads(1)
    cs <- cases[[tasks$case[k]]]
    t <- system.time(f <- zipln_with_init(inits[[tasks$init[k]]], cs[[3]], data = cs[[2]], zi = cs[[4]],
                                          control = ZIPLN_param(trace = 0)))[["elapsed"]]
    data.frame(data = cs[[1]], model = paste(deparse(cs[[3]][[3]]), collapse = ""), zi = cs[[4]], init = tasks$init[k], summ(f, t))
  }
  ord <- order(-(tasks$case >= 12))          # the heavy microcosm fits first
  out <- parallel::mclapply(ord, function(k) tryCatch(one(k), error = function(e) conditionMessage(e)),
                            mc.cores = n_cores, mc.preschedule = FALSE)
}
err <- Filter(is.character, out); if (length(err)) { cat(length(err), "errors:\n"); print(unique(unlist(err))) }
res <- do.call(rbind, Filter(is.data.frame, out))
saveRDS(res, sprintf("inst/simus_ZIPLN/init/benchmark_%s.rds", what))
cat("done", what, nrow(res), "rows\n")
