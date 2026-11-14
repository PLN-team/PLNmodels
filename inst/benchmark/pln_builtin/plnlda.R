## PLNLDA: new default (builtin, SQUAREM, loose VE steps, tests on the bound with its
## constants) against the plain EM loop of version 1.3.2. Sequential, one BLAS thread.
suppressMessages(devtools::load_all(".", quiet = TRUE)); RhpcBLASctl::blas_set_num_threads(1)
data(trichoptera); data(oaks); data(mollusk); data(microcosm)
tri <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
mol <- suppressMessages(prepare_data(mollusk$Abundance, mollusk$Covariate))
legacy <- list(objective_constant = FALSE, accelerate = FALSE, ftol_in = 1e-8, ftol_em = 1e-8, maxit_em = 200L)
lda <- list(list("trichoptera", tri, Abundance ~ 1 + offset(log(Offset)), tri$Group),
            list("oaks", oaks, Abundance ~ 1 + offset(log(Offset)), oaks$tree),
            list("mollusk", mol, Abundance ~ 1 + offset(log(Offset)), mol$site),
            list("microcosm", microcosm, Abundance ~ 1 + offset(log(Offset)), microcosm$site))
for (cs in lda) for (vn in c("legacy", "default")) {
  co <- if (vn == "legacy") legacy else list()
  t <- system.time(f <- do.call(PLNLDA, list(cs[[3]], grouping = cs[[4]], data = cs[[2]],
                                            control = PLNLDA_param(trace = 0, config_optim = co))))[[3]]
  cat(format(Sys.time(), "%H:%M:%S"), "PLNLDA", cs[[1]], vn, round(f$loglik, 2), round(t, 2), "\n",
      file = "inst/benchmark/pln_builtin/plnlda_log.txt", append = TRUE)
}
