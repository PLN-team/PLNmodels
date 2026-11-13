## How does the outer loop of ZIPLN converge? Objective along the iterations, run to
## convergence (maxit_out = 5000), from the starting point kept by init = "best".
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("inst/simus_ZIPLN/convergence/cases.R")
out <- parallel::mclapply(seq_along(cases), function(k) {
  RhpcBLASctl::blas_set_num_threads(1)
  cs <- cases[[k]]
  t <- system.time(f <- suppressWarnings(suppressMessages(ZIPLN(cs$formula, data = cs$data, zi = cs$zi,
    control = ZIPLN_param(trace = 0, config_optim = list(maxit_out = 5000))))))[["elapsed"]]
  list(case = cs$name, objective = f$optim_par$objective, time = t, init = f$optim_par$init,
       stop = f$optim_par$message, loglik = f$loglik)
}, mc.cores = 13, mc.preschedule = FALSE)
saveRDS(out, "inst/simus_ZIPLN/convergence/diagnose.rds")
summ <- do.call(rbind, lapply(out, function(x) {
  o <- x$objective; n <- length(o); fin <- o[n]
  gap <- function(it) if (it <= n) o[it] - fin else 0
  d <- -diff(o); d <- d[d > 0]
  ## empirical linear rate over the last half of the decreases
  r <- if (length(d) > 20) exp(coef(lm(log(tail(d, length(d) %/% 2)) ~ seq_len(length(d) %/% 2)))[2]) else NA
  data.frame(case = x$case, init = x$init, iterations = n, time = round(x$time, 1), time_per_it = signif(x$time / n, 2),
             gap_at_50 = signif(gap(50), 3), gap_at_100 = signif(gap(100), 3), gap_at_200 = signif(gap(200), 3),
             gap_at_500 = signif(gap(500), 3), rate = signif(r, 4))
}))
options(width = 200); print(summ, row.names = FALSE)
