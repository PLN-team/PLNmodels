## Control: refit from the converged fit without any switch. Part of the gains of the
## greedy search may only come from continuing the optimization.
src <- readLines("inst/simus_ZIPLN/species_regime/flips.R")
eval(parse(text = src[1:(grep("^run <- function", src) - 1)]))
out <- parallel::mclapply(seq_along(cases), function(k) {
  RhpcBLASctl::blas_set_num_threads(1)
  cs <- cases[[k]]
  f <- suppressWarnings(suppressMessages(ZIPLN(cs$formula, data = cs$data, zi = cs$zi, control = ZIPLN_param(trace = 0))))
  h <- refit(cs, f); h2 <- refit(cs, h)
  data.frame(case = cs$name, noop_gain = h$loglik - f$loglik, noop_gain_twice = h2$loglik - f$loglik,
             iterations_fit = f$optim_par$iterations, stop = f$optim_par$message)
}, mc.cores = 13, mc.preschedule = FALSE)
res <- do.call(rbind, out); saveRDS(res, "inst/simus_ZIPLN/species_regime/control_noop.rds")
g <- do.call(rbind, lapply(readRDS("inst/simus_ZIPLN/species_regime/greedy.rds"), `[[`, "summary"))
res <- merge(res, g[, c("case", "total_gain", "accepted")]); options(width = 200)
print(res, digits = 4, row.names = FALSE)
