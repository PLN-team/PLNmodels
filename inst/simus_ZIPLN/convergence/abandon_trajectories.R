## Full trajectories (new defaults) of the two starting points on the 13 cases, to replay the
## abandonment rule offline (abandon_rules.R)
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("inst/simus_ZIPLN/convergence/cases.R")
jobs <- expand.grid(k = seq_along(cases), init = c("lm", "positive"), stringsAsFactors = FALSE)
out <- parallel::mclapply(seq_len(nrow(jobs)), function(i) {
  RhpcBLASctl::blas_set_num_threads(1)
  cs <- cases[[jobs$k[i]]]; dd <- cs$data
  ctrl <- ZIPLN_param(trace = 0, init = jobs$init[i])
  f <- eval(bquote(ZIPLN(.(cs$formula), data = dd, zi = .(if (cs$zi == "covar") "single" else cs$zi), control = ctrl)))
  list(case = cs$name, init = jobs$init[i], objective = f$optim_par$objective, VE = f$optim_par$iterations)
}, mc.cores = 13, mc.preschedule = FALSE)
saveRDS(Filter(is.list, out), "inst/simus_ZIPLN/convergence/abandon_trajectories.rds")
