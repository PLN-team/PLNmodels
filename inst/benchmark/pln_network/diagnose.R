## PLNnetwork (builtin): where is the gap between the VE variants along the path?
## Per model of the path: penalized loglik, loglik, edges, outer iterations (cap: maxit_em).
## Run against an installed version (names before #202: ftol_em, maxit_em for the alternating loop):
## Rscript diagnose.R <lib> <out.rds>. With the default, maxit_ve = 1 had no effect (full VE steps).
args <- commandArgs(TRUE)
library(PLNmodels, lib.loc = args[1]); RhpcBLASctl::blas_set_num_threads(1)
for (f in c("trichoptera", "oaks", "mollusk")) load(file.path("data", paste0(f, ".rda")))
tri <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
mol <- suppressMessages(prepare_data(mollusk$Abundance, mollusk$Covariate))
cases <- list(trichoptera = tri, oaks = oaks, mollusk = mol)
variants <- list(
  default  = list(),
  out7     = list(ftol_em = 1e-7, maxit_em = 200L),
  cg5      = list(cg_max = 5L),
  cg5_out7 = list(cg_max = 5L, ftol_em = 1e-7, maxit_em = 200L),
  nw1_out7 = list(maxeval = 1L, ftol_em = 1e-7, maxit_em = 500L),
  nw3_out7 = list(maxeval = 3L, ftol_em = 1e-7, maxit_em = 500L))
log_file <- sub("rds$", "log", args[2])
res <- NULL
for (cn in names(cases)) for (vn in names(variants)) {
  t <- system.time(fam <- suppressWarnings(suppressMessages(
    PLNnetwork(Abundance ~ 1 + offset(log(Offset)), data = cases[[cn]],
               control = PLNnetwork_param(trace = 0, config_optim = variants[[vn]])))))[["elapsed"]]
  it <- sapply(fam$models, function(m) m$optim_par$iterations)
  cr <- fam$criteria
  row <- data.frame(case = cn, variant = vn, penalty = cr$param, loglik = cr$loglik,
                    pen_loglik = cr$pen_loglik, n_edges = cr$n_edges, outer = it, time = t)
  cat(format(Sys.time(), "%H:%M:%S"), cn, vn, sum(cr$pen_loglik), sum(cr$loglik),
      sum(it), max(it), round(t, 2), "\n", file = log_file, append = TRUE)
  res <- rbind(res, row)
  saveRDS(res, args[2])
}
