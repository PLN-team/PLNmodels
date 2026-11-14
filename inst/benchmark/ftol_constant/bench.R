## PLN: convergence tests on the bound with its constants (config_optim$objective_constant).
## The objectives of the C++ optimizers omit sum(w log(Y!)), which dominates them on large
## counts, so that their relative tolerances are looser than they read, the more so as the
## counts are large. Compared: the current tests (default tolerances) against the shifted
## tests with tolerances 1e-6, 1e-7, 1e-8, for both backends and two covariance models.
suppressMessages(devtools::load_all(".", quiet = TRUE)); RhpcBLASctl::blas_set_num_threads(1)
data(trichoptera); data(barents); data(oaks); data(mollusk); data(microcosm)
tri <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
mol <- suppressMessages(prepare_data(mollusk$Abundance, mollusk$Covariate))
cases <- list(
  list(name = "trichoptera", data = tri,       formula = Abundance ~ 1 + offset(log(Offset))),
  list(name = "barents",     data = barents,   formula = Abundance ~ 1 + offset(log(Offset))),
  list(name = "mollusk",     data = mol,       formula = Abundance ~ 1 + offset(log(Offset))),
  list(name = "oaks",        data = oaks,      formula = Abundance ~ 1 + offset(log(Offset))),
  list(name = "oaks tree",   data = oaks,      formula = Abundance ~ 0 + tree + offset(log(Offset))),
  list(name = "microcosm",   data = microcosm, formula = Abundance ~ 0 + site + offset(log(Offset))))
configs <- list(current = list(),
                shift_1e6 = list(objective_constant = TRUE, tol = 1e-6),
                shift_1e7 = list(objective_constant = TRUE, tol = 1e-7),
                shift_1e8 = list(objective_constant = TRUE, tol = 1e-8))
tol_entries <- function(backend, tol) if (backend == "nlopt") list(ftol_rel = tol) else list(ftol_in = tol, ftol_em = tol)
log_file <- "inst/benchmark/ftol_constant/log.txt"
res <- NULL
for (cs in cases) for (backend in c("nlopt", "builtin")) for (cov in c("full", "diagonal")) for (cn in names(configs)) {
  cf <- configs[[cn]]
  co <- if (is.null(cf$tol)) list() else c(list(objective_constant = TRUE), tol_entries(backend, cf$tol))
  ctrl <- PLN_param(trace = 0, backend = backend, covariance = cov, config_optim = co)
  t <- system.time(fit <- suppressWarnings(PLN(cs$formula, data = cs$data, control = ctrl)))[["elapsed"]]
  row <- data.frame(case = cs$name, backend = backend, covariance = cov, config = cn,
                    loglik = fit$loglik, iterations = fit$optim_par$iterations, time = round(t, 2))
  cat(format(Sys.time(), "%H:%M:%S"), paste(row, collapse = " "), "\n", file = log_file, append = TRUE)
  res <- rbind(res, row)
  saveRDS(res, "inst/benchmark/ftol_constant/bench.rds")
}
base <- res[res$config == "current", ]
key <- function(d) paste(d$case, d$backend, d$covariance)
res$gain <- round(res$loglik - base$loglik[match(key(res), key(base))], 2)
res$time_ratio <- round(res$time / base$time[match(key(res), key(base))], 2)
options(width = 200); print(res, row.names = FALSE)
saveRDS(res, "inst/benchmark/ftol_constant/bench.rds")
