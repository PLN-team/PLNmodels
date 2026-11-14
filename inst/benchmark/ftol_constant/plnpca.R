## PLNPCA: convergence tests on the bound with its constants (objective_constant), both backends.
## nlopt (default): ftol_rel on the shifted objective; builtin: ftol_out on the shifted objective.
## Criterion: log-likelihood of each rank; sequential runs, one BLAS thread.
suppressMessages(devtools::load_all(".", quiet = TRUE)); RhpcBLASctl::blas_set_num_threads(1)
data(trichoptera); data(barents); data(oaks); data(mollusk); data(microcosm)
tri <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
mol <- suppressMessages(prepare_data(mollusk$Abundance, mollusk$Covariate))
cases <- list(
  list(name = "trichoptera", data = tri,       ranks = c(2, 5)),
  list(name = "barents",     data = barents,   ranks = c(2, 5, 10)),
  list(name = "mollusk",     data = mol,       ranks = c(2, 5, 10)),
  list(name = "oaks",        data = oaks,      ranks = c(2, 5, 10)),
  list(name = "microcosm",   data = microcosm, ranks = c(2, 5, 10)))
variants <- list(
  nlopt_default   = list(backend = "nlopt",   co = list()),
  nlopt_shift_8   = list(backend = "nlopt",   co = list(objective_constant = TRUE)),
  nlopt_shift_7   = list(backend = "nlopt",   co = list(objective_constant = TRUE, ftol_rel = 1e-7)),
  builtin_default = list(backend = "builtin", co = list()),
  builtin_shift_8 = list(backend = "builtin", co = list(objective_constant = TRUE)))
log_file <- "inst/benchmark/ftol_constant/plnpca_log.txt"
res <- NULL
for (cs in cases) for (vn in names(variants)) {
  v <- variants[[vn]]
  t <- system.time(fam <- suppressWarnings(suppressMessages(
    PLNPCA(Abundance ~ 1 + offset(log(Offset)), data = cs$data, ranks = cs$ranks,
           control = PLNPCA_param(trace = 0, backend = v$backend, config_optim = v$co)))))[["elapsed"]]
  row <- data.frame(case = cs$name, variant = vn, rank = cs$ranks, loglik = sapply(fam$models, function(m) m$loglik), time = round(t, 2))
  cat(format(Sys.time(), "%H:%M:%S"), cs$name, vn, round(row$loglik, 2), row$time[1], "\n", file = log_file, append = TRUE)
  res <- rbind(res, row); saveRDS(res, "inst/benchmark/ftol_constant/plnpca.rds")
}
