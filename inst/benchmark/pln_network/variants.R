## PLNnetwork (builtin): VE steps at fixed Omega and convergence of the alternating loop.
## Same inception (diagonal PLN, default of PLN()) hence same grid for all the variants,
## so that the sums of the penalized log-likelihoods along the path are comparable.
## Sequential runs, one BLAS thread. Rscript variants.R [case ...]
suppressMessages(devtools::load_all(".", quiet = TRUE)); RhpcBLASctl::blas_set_num_threads(1)
data(trichoptera); data(oaks); data(mollusk); data(barents)
tri <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
mol <- suppressMessages(prepare_data(mollusk$Abundance, mollusk$Covariate))
cases <- list(trichoptera = tri, barents = barents, mollusk = mol, oaks = oaks)
if (length(commandArgs(TRUE))) cases <- cases[commandArgs(TRUE)]
legacy <- list(ftol_in = 1e-8, objective_constant = FALSE)
variants <- list(
  legacy_full_out5  = list(ve = 10000L, co = c(legacy, list(ftol_out = 1e-5, maxit_out = 20L))),
  ve1_out5          = list(ve = 1L,     co = list(ftol_out = 1e-5, maxit_out = 20L)),
  ve1_out6          = list(ve = 1L,     co = list(ftol_out = 1e-6, maxit_out = 500L)),
  ve1_out7          = list(ve = 1L,     co = list(ftol_out = 1e-7, maxit_out = 1000L)),
  ve3_out7          = list(ve = 3L,     co = list(ftol_out = 1e-7, maxit_out = 1000L)),
  full_out7         = list(ve = 10000L, co = list(ftol_out = 1e-7, maxit_out = 1000L)),
  cg5_full_out7     = list(ve = 10000L, co = list(cg_max = 5L, ftol_out = 1e-7, maxit_out = 1000L)),
  cg5_ve1_out7      = list(ve = 1L,     co = list(cg_max = 5L, ftol_out = 1e-7, maxit_out = 1000L)),
  cg5_ve3_out7      = list(ve = 3L,     co = list(cg_max = 5L, ftol_out = 1e-7, maxit_out = 1000L)))
dir <- "inst/benchmark/pln_network"
log_file <- file.path(dir, "variants.log")
res <- if (file.exists(file.path(dir, "variants.rds"))) readRDS(file.path(dir, "variants.rds")) else NULL
for (cn in names(cases)) {
  inc <- PLN(Abundance ~ 1 + offset(log(Offset)), data = cases[[cn]],
             control = PLN_param(trace = 0, covariance = "diagonal"))
  for (vn in names(variants)) {
    v <- variants[[vn]]
    t <- system.time(fam <- suppressWarnings(suppressMessages(
      PLNnetwork(Abundance ~ 1 + offset(log(Offset)), data = cases[[cn]],
                 control = PLNnetwork_param(trace = 0, inception = inc, maxit_ve = v$ve, config_optim = v$co)))))[["elapsed"]]
    it <- sapply(fam$models, function(m) m$optim_par$iterations)
    cr <- fam$criteria
    row <- data.frame(case = cn, variant = vn, penalty = cr$param, loglik = cr$loglik,
                      pen_loglik = cr$pen_loglik, n_edges = cr$n_edges, outer = it, time = t)
    cat(format(Sys.time(), "%H:%M:%S"), cn, vn, round(sum(cr$pen_loglik), 2), round(sum(cr$loglik), 2),
        sum(it), max(it), round(t, 2), "\n", file = log_file, append = TRUE)
    res <- rbind(res, row)
    saveRDS(res, file.path(dir, "variants.rds"))
  }
}
