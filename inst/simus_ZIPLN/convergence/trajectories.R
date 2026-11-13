## Trajectories of the accelerated fits at ftol_out = 1e-9 (objective per cycle and number of
## VE steps per cycle), to replay stopping rules offline (stopping_rules.R).
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("inst/simus_ZIPLN/convergence/cases.R")
pick <- list(c("microcosm site col", "lm"), c("microcosm site|site", "positive"), c("barents col", "positive"),
             c("mollusk col", "lm"), c("mollusk site|site", "lm"), c("trichoptera col", "lm"), c("oaks col", "lm"),
             c("barents single", "lm"), c("oaks single", "lm"), c("oaks tree|tree", "lm"), c("mollusk single", "lm"),
             c("trichoptera single", "lm"), c("oaks tree col", "lm"))
## record the number of VE steps of each cycle
trace_env <- new.env(); trace_env$evals <- integer(0)
sq <- squarem_cycle
assignInNamespace("squarem_cycle", function(...) { out <- sq(...); trace_env$evals <- c(trace_env$evals, out$evaluations); out }, "PLNmodels")
out <- parallel::mclapply(pick, function(pk) {
  RhpcBLASctl::blas_set_num_threads(1)
  cs <- Filter(function(x) x$name == pk[1], cases)[[1]]
  ctrl <- ZIPLN_param(trace = 0, init = pk[2], config_optim = list(ftol_out = 1e-9, maxit_out = 20000))
  dd <- cs$data; trace_env$evals <- integer(0)
  t <- system.time(f <- eval(bquote(ZIPLN(.(cs$formula), data = dd, zi = .(if (cs$zi == "covar") "single" else cs$zi), control = ctrl))))[["elapsed"]]
  list(case = pk[1], objective = f$optim_par$objective, evals = c(1L, trace_env$evals), time = t,
       iterations = f$optim_par$iterations, loglik = f$loglik)
}, mc.cores = 13, mc.preschedule = FALSE)
print(Filter(function(x) inherits(x, "try-error"), out))
saveRDS(out, "inst/simus_ZIPLN/convergence/trajectories.rds")
for (o in out) if (is.list(o)) cat(o$case, length(o$objective), "cycles", o$iterations, "VE", sum(o$evals), round(o$time), "s\n")
