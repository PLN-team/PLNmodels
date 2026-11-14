## Builtin PLN: SQUAREM on the EM loop (accelerate), loose VE steps (ftol_in), Newton-CG VE
## step (cg_max), convergence tests on the bound with its constants (objective_constant).
## Sequential runs, one BLAS thread; gain against the current default of the builtin backend.
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
sq <- function(ve, em, cg = 1L) list(objective_constant = TRUE, accelerate = TRUE, ftol_in = ve, ftol_em = em, maxit_em = 2000L, cg_max = cg)
variants <- list(default    = list(),
                 sq_ve4_em7 = sq(1e-4, 1e-7), sq_ve6_em7 = sq(1e-6, 1e-7),
                 sq_ve4_em8 = sq(1e-4, 1e-8), sq_ve6_em8 = sq(1e-6, 1e-8),
                 sq_cg5_ve4_em8 = sq(1e-4, 1e-8, 5L), sq_cg5_ve6_em8 = sq(1e-6, 1e-8, 5L),
                 plain_ve6_em8 = modifyList(sq(1e-6, 1e-8), list(accelerate = FALSE)))
log_file <- "inst/benchmark/pln_builtin/bench_log.txt"
covs <- if (length(commandArgs(TRUE))) commandArgs(TRUE) else c("full", "diagonal", "spherical")
res <- NULL
for (cov in covs) for (cs in cases) for (vn in names(variants)) {
  ctrl <- PLN_param(trace = 0, backend = "builtin", covariance = cov, config_optim = variants[[vn]])
  t <- system.time(fit <- suppressWarnings(PLN(cs$formula, data = cs$data, control = ctrl)))[["elapsed"]]
  row <- data.frame(case = cs$name, covariance = cov, variant = vn, loglik = fit$loglik,
                    newton = fit$optim_par$iterations, time = round(t, 2))
  cat(format(Sys.time(), "%H:%M:%S"), paste(row, collapse = " "), "\n", file = log_file, append = TRUE)
  res <- rbind(res, row)
  saveRDS(res, sprintf("inst/benchmark/pln_builtin/bench_%s.rds", paste(covs, collapse = "_")))
}
