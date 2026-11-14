## PLNmixture: tolerances of the mixture loop (ftol_out, maxit_out) separated from those
## of the PLN components (ftol_em, maxit_em, ftol_in). Sum over clusters 1:4 of the
## log-likelihoods, time. Sequential runs, one BLAS thread. Rscript variants.R [case ...]
suppressMessages(devtools::load_all(".", quiet = TRUE)); RhpcBLASctl::blas_set_num_threads(1)
data(trichoptera); data(oaks); data(mollusk); data(barents)
tri <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
mol <- suppressMessages(prepare_data(mollusk$Abundance, mollusk$Covariate))
cases <- list(trichoptera = tri, barents = barents, mollusk = mol, oaks = oaks)
if (length(commandArgs(TRUE))) cases <- cases[commandArgs(TRUE)]
## before the separation: the components shared the tolerances of the mixture loop
legacy <- list(ftol_em = 1e-3, maxit_em = 50L, ftol_in = 1e-8, accelerate = FALSE, objective_constant = FALSE)
variants <- list(
  legacy_out3     = c(legacy, list(ftol_out = 1e-3, maxit_out = 50L)),
  out3            = list(ftol_out = 1e-3, maxit_out = 50L),
  out5            = list(ftol_out = 1e-5, maxit_out = 200L),
  out6            = list(ftol_out = 1e-6, maxit_out = 200L),
  em5_out5        = list(ftol_em = 1e-5, ftol_out = 1e-5, maxit_out = 200L),
  em5_out6        = list(ftol_em = 1e-5, ftol_out = 1e-6, maxit_out = 200L),
  legacy_out6     = c(legacy, list(ftol_out = 1e-6, maxit_out = 200L)))
dir <- "inst/benchmark/pln_mixture"
log_file <- file.path(dir, "variants.log")
res <- NULL
for (cn in names(cases)) for (vn in names(variants)) {
  set.seed(1)
  t <- system.time(fam <- suppressWarnings(suppressMessages(
    PLNmixture(Abundance ~ 1 + offset(log(Offset)), clusters = 1:4, data = cases[[cn]],
               control = PLNmixture_param(trace = 0, config_optim = variants[[vn]])))))[["elapsed"]]
  cr <- fam$criteria
  it <- sapply(fam$models, function(m) m$optim_par$outer_iterations)
  row <- data.frame(case = cn, variant = vn, k = cr$param, loglik = cr$loglik, outer = it, time = t)
  cat(format(Sys.time(), "%H:%M:%S"), cn, vn, round(sum(cr$loglik), 2), paste(round(cr$loglik, 1), collapse = " "),
      paste(it, collapse = "/"), round(t, 2), "\n", file = log_file, append = TRUE)
  res <- rbind(res, row)
  saveRDS(res, file.path(dir, "variants.rds"))
}
