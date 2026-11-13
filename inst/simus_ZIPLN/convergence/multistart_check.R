## New defaults (SQUAREM, maxit_out = 500) with the two starting points of init = "best":
## does abandoning the second start still give the best of the two full fits? Time and bound
## against the former defaults (no acceleration, maxit_out = 200).
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("inst/simus_ZIPLN/convergence/cases.R")
fit <- function(cs, ...) {
  dd <- cs$data; ctrl <- ZIPLN_param(trace = 0, ...)
  t <- system.time(f <- eval(bquote(ZIPLN(.(cs$formula), data = dd, zi = .(if (cs$zi == "covar") "single" else cs$zi), control = ctrl))))[["elapsed"]]
  list(loglik = f$loglik, time = t, VE = f$optim_par$iterations, init = f$optim_par$init,
       abandoned = paste(f$optim_par$abandoned_starts, collapse = ","))
}
out <- parallel::mclapply(cases, function(cs) {
  RhpcBLASctl::blas_set_num_threads(1)
  old  <- fit(cs, config_optim = list(accelerate = FALSE, maxit_out = 200))
  best <- fit(cs)
  lm   <- fit(cs, init = "lm"); pos <- fit(cs, init = "positive")
  data.frame(case = cs$name, old_loglik = old$loglik, old_time = round(old$time, 1),
             new_loglik = best$loglik, new_time = round(best$time, 1), kept = best$init, abandoned = best$abandoned,
             full_lm = lm$loglik, full_positive = pos$loglik, two_full_time = round(lm$time + pos$time, 1),
             missed = round(max(lm$loglik, pos$loglik) - best$loglik, 3))
}, mc.cores = 13, mc.preschedule = FALSE)
print(Filter(is.character, out)); res <- do.call(rbind, Filter(is.data.frame, out))
options(width = 250); print(res, digits = 9, row.names = FALSE)
saveRDS(res, "inst/simus_ZIPLN/convergence/multistart_check.rds")
