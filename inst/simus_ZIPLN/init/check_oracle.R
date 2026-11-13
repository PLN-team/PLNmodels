## Does the package itself (builtin and nlopt) reach a poor optimum from its starting point?
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("inst/simus_ZIPLN/init/simulate.R")
source("inst/simus_ZIPLN/init/inits.R")
tasks <- (which(configs$zi == "high" & configs$abundance == "high" & configs$rho == .5) - 1) * reps + 1:6
one_task <- function(task) {
  RhpcBLASctl::blas_set_num_threads(1)
  cfg <- configs[(task - 1) %/% reps + 1, ]; sim <- simulate(cfg, task)
  dat <- data.frame(row.names = seq_len(nrow(sim$Y))); dat$Y <- sim$Y
  do.call(rbind, lapply(c("builtin", "nlopt"), function(bk) do.call(rbind, lapply(
    list(current = init_current, positive = init_positive, zip = init_zip, truth = init_truth(sim)), function(init) {
      f <- zipln_with_init(init, Y ~ 1, data = dat, zi = "col", control = ZIPLN_param(trace = 0, backend = bk))
      data.frame(task = task, backend = bk, init = NA, loglik = f$loglik,
                 err_pi = sqrt(mean((f$model_par$Pi[1, ] - sim$pi)^2)),
                 err_var = mean(abs(log(diag(as.matrix(f$model_par$Sigma)) / diag(sim$Sigma)))))
    }))))
}
out <- do.call(rbind, parallel::mclapply(tasks, one_task, mc.cores = 6, mc.preschedule = FALSE))
out$init <- rep(c("current", "positive", "zip", "truth"), length.out = nrow(out))
stopifnot(is.data.frame(out))
saveRDS(out, "inst/simus_ZIPLN/init/check_oracle.rds")
print(aggregate(cbind(loglik, err_pi, err_var) ~ init + backend, out, mean), digits = 5)
