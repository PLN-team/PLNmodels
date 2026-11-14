## Candidate settings of the builtin solver on the variants of PLN: SQUAREM, VE steps to 1e-4
## and EM loop to 1e-8 on the bound with its constants, Newton-CG VE step (cg_max = 5).
## PLNLDA (full covariance), PLNmixture (components with weights; ftol_em and maxit_em are
## also those of the mixture loop, kept), PLNnetwork (VE steps at fixed Omega: only the
## Newton-CG step and the shifted ftol_in apply). Sequential runs, one BLAS thread.
suppressMessages(devtools::load_all(".", quiet = TRUE)); RhpcBLASctl::blas_set_num_threads(1)
data(trichoptera); data(oaks); data(mollusk); data(microcosm)
tri <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
mol <- suppressMessages(prepare_data(mollusk$Abundance, mollusk$Covariate))
cand_pln <- list(objective_constant = TRUE, accelerate = TRUE, ftol_in = 1e-4, ftol_em = 1e-8, maxit_em = 2000L, cg_max = 5L)
cand_ve  <- list(objective_constant = TRUE, accelerate = TRUE, ftol_in = 1e-4, cg_max = 5L)
log_file <- "inst/benchmark/pln_builtin/variants_log.txt"
res <- NULL
record <- function(model, case, variant, expr) {
  t <- system.time(val <- suppressWarnings(suppressMessages(expr)))[["elapsed"]]
  row <- data.frame(model = model, case = case, variant = variant, criterion = val, time = round(t, 2))
  cat(format(Sys.time(), "%H:%M:%S"), paste(row, collapse = " "), "\n", file = log_file, append = TRUE)
  res <<- rbind(res, row); saveRDS(res, "inst/benchmark/pln_builtin/variants.rds")
}
## PLNLDA: log-likelihood of the fit
lda <- list(list("trichoptera", tri, Abundance ~ 1 + offset(log(Offset)), quote(Group)),
            list("oaks", oaks, Abundance ~ 1 + offset(log(Offset)), quote(tree)),
            list("mollusk", mol, Abundance ~ 1 + offset(log(Offset)), quote(site)),
            list("microcosm", microcosm, Abundance ~ 1 + offset(log(Offset)), quote(site)))
if (!identical(commandArgs(TRUE), "network")) for (cs in lda) for (vn in c("default", "candidate")) {
  co <- if (vn == "default") list() else cand_pln
  record("PLNLDA", cs[[1]], vn, do.call(PLNLDA, list(cs[[3]], grouping = cs[[4]], data = cs[[2]],
                                                    control = PLNLDA_param(trace = 0, config_optim = co)))$loglik)
}
## PLNmixture: sum over the numbers of clusters of the log-likelihoods
mix <- list(list("trichoptera", tri), list("mollusk", mol), list("oaks", oaks))
if (!identical(commandArgs(TRUE), "network")) for (cs in mix) for (vn in c("default", "candidate")) {
  co <- if (vn == "default") list() else cand_ve
  record("PLNmixture", cs[[1]], vn, {
    fam <- PLNmixture(Abundance ~ 1 + offset(log(Offset)), clusters = 1:4, data = cs[[2]],
                      control = PLNmixture_param(trace = 0, config_optim = co))
    sum(sapply(fam$models, function(m) m$loglik))
  })
}
## PLNnetwork: sum over the path of the penalized log-likelihoods
if (identical(commandArgs(TRUE), "network")) res <- readRDS("inst/benchmark/pln_builtin/variants.rds")
net <- list(list("trichoptera", tri), list("oaks", oaks), list("mollusk", mol))
for (cs in net) for (vn in c("default", "cg5", "cg5_shift")) {
  co <- switch(vn, default = list(), cg5 = list(cg_max = 5L), cg5_shift = list(cg_max = 5L, objective_constant = TRUE, ftol_in = 1e-6))
  record("PLNnetwork", cs[[1]], vn, {
    fam <- PLNnetwork(Abundance ~ 1 + offset(log(Offset)), data = cs[[2]],
                      control = PLNnetwork_param(trace = 0, n_penalties = 20, config_optim = co))
    sum(fam$criteria$pen_loglik)
  })
}
options(width = 200); print(res, row.names = FALSE)
