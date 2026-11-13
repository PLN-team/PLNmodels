## In the fit kept by ZIPLN(), where do the zeros of each species come from? Latent
## variance (Sigma_jj large, pi_j small) or zero inflation? Real data sets of the benchmark.
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("inst/simus_ZIPLN/species_regime/cases.R")
out <- parallel::mclapply(seq_along(cases), function(k) {
  RhpcBLASctl::blas_set_num_threads(1)
  cs <- cases[[k]]
  f <- suppressWarnings(suppressMessages(ZIPLN(cs$formula, data = cs$data, zi = cs$zi, control = ZIPLN_param(trace = 0))))
  v <- diag(as.matrix(f$model_par$Sigma)); pi <- colMeans(f$model_par$Pi)
  data.frame(case = cs$name, kept = f$optim_par$init, loglik = f$loglik, p = length(v),
             var_median = median(v), var_max = max(v), n_var_100 = sum(v > 100), n_var_5med = sum(v > 5 * median(v)),
             pi_of_high_var = if (any(v > 5 * median(v))) median(pi[v > 5 * median(v)]) else NA, pi_median = median(pi))
}, mc.cores = 13, mc.preschedule = FALSE)
options(width = 200); print(do.call(rbind, out), digits = 3, row.names = FALSE)
