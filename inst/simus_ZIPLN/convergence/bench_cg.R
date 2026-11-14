## Timing of ZIPLN() (default settings, two starting points, SQUAREM) with the current
## cell-wise VE step and with the Newton-CG VE step (C++), on the 13 real cases.
## Usage: Rscript bench_cg.R <library> <label> [cg_max ftol_in]
args <- commandArgs(TRUE); lib <- args[1]; label <- args[2]
suppressMessages(library(PLNmodels, lib.loc = lib)); RhpcBLASctl::blas_set_num_threads(1)
source("inst/simus_ZIPLN/convergence/cases.R")
cfg <- if (length(args) > 2) list(cg_max = as.integer(args[3]), ftol_in = as.numeric(args[4])) else list()
res <- do.call(rbind, lapply(cases, function(cs) {
  t <- system.time(fit <- suppressWarnings(suppressMessages(
    ZIPLN(cs$formula, data = cs$data, zi = if (cs$zi == "covar") "single" else cs$zi,
          control = ZIPLN_param(trace = 0, config_optim = cfg)))))[["elapsed"]]
  data.frame(case = cs$name, label = label, loglik = fit$loglik, time = round(t, 2))
}))
print(res, row.names = FALSE)
saveRDS(res, sprintf("inst/simus_ZIPLN/convergence/ve_cg/bench_%s.rds", label))
