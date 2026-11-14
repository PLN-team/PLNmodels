## PLNnetwork: default of the installed version given in argument (or load_all()), on real data.
## Rscript defaults.R <label> [lib]. Sequential, one BLAS thread.
args <- commandArgs(TRUE)
if (length(args) > 1) library(PLNmodels, lib.loc = args[2]) else suppressMessages(devtools::load_all(".", quiet = TRUE))
RhpcBLASctl::blas_set_num_threads(1)
for (f in c("trichoptera", "barents", "oaks", "mollusk")) load(file.path("data", paste0(f, ".rda")))
tri <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
mol <- suppressMessages(prepare_data(mollusk$Abundance, mollusk$Covariate))
cases <- list(trichoptera = tri, barents = barents, mollusk = mol, oaks = oaks)
for (cn in names(cases)) {
  t <- system.time(fam <- suppressWarnings(suppressMessages(
    PLNnetwork(Abundance ~ 1 + offset(log(Offset)), data = cases[[cn]], control = PLNnetwork_param(trace = 0)))))[["elapsed"]]
  cr <- fam$criteria; best <- fam$getBestModel("BIC")
  cat(format(Sys.time(), "%H:%M:%S"), args[1], cn, "sum pen_loglik", round(sum(cr$pen_loglik), 1), "sum loglik", round(sum(cr$loglik), 1),
      "max penalty", signif(max(cr$param), 4), "BIC edges", best$n_edges, "time", round(t, 2), "\n",
      file = "inst/benchmark/pln_network/defaults.log", append = TRUE)
}
