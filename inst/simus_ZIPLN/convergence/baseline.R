## Current default (no acceleration, ftol_out = 1e-6, maxit_out = 200), single start, on the
## cases of trajectories.R: VE steps and final objective, to compare with the replayed rules.
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("inst/simus_ZIPLN/convergence/cases.R")
tr <- Filter(is.list, readRDS("inst/simus_ZIPLN/convergence/trajectories.rds"))
inits <- c("microcosm site|site" = "positive", "barents col" = "positive")
out <- parallel::mclapply(tr, function(o) {
  RhpcBLASctl::blas_set_num_threads(1)
  cs <- Filter(function(x) x$name == o$case, cases)[[1]]
  init <- if (o$case %in% names(inits)) inits[[o$case]] else "lm"
  ctrl <- ZIPLN_param(trace = 0, init = init, config_optim = list(accelerate = FALSE))
  dd <- cs$data
  t <- system.time(f <- eval(bquote(ZIPLN(.(cs$formula), data = dd, zi = .(if (cs$zi == "covar") "single" else cs$zi), control = ctrl))))[["elapsed"]]
  data.frame(case = o$case, objective = tail(f$optim_par$objective, 1), VE = f$optim_par$iterations,
             message = f$optim_par$message, time = round(t, 1))
}, mc.cores = 13, mc.preschedule = FALSE)
saveRDS(do.call(rbind, Filter(is.data.frame, out)), "inst/simus_ZIPLN/convergence/baseline.rds")
