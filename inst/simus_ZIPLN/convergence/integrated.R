## SQUAREM as integrated in ZIPLNfit$optimize(): gap to the reference optimum (squarem.R) and
## time, for the stopping tolerances 1e-6 (current default), 1e-7 and 1e-8, single start.
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("inst/simus_ZIPLN/convergence/cases.R")
ref <- readRDS("inst/simus_ZIPLN/convergence/squarem_adaptive.rds")
pick <- list(c("barents col", "positive"), c("mollusk col", "lm"), c("mollusk site|site", "lm"), c("trichoptera col", "lm"),
             c("microcosm site col", "lm"), c("microcosm site|site", "positive"), c("oaks col", "lm"), c("barents single", "lm"))
grid <- expand.grid(k = seq_along(pick), accelerate = c(FALSE, TRUE), ftol = c(1e-6, 1e-7, 1e-8))
grid <- grid[grid$accelerate | grid$ftol == 1e-6, ]
grid <- grid[order(-(grepl("microcosm", sapply(pick, `[`, 1))[grid$k])), ]   # long jobs first
out <- parallel::mclapply(seq_len(nrow(grid)), function(i) {
  RhpcBLASctl::blas_set_num_threads(1)
  g <- grid[i, ]; pk <- pick[[g$k]]
  cs <- Filter(function(x) x$name == pk[1], cases)[[1]]
  ctrl <- ZIPLN_param(trace = 0, init = pk[2],
                      config_optim = list(accelerate = g$accelerate, ftol_out = g$ftol, maxit_out = 20000))
  dd <- cs$data
  t <- system.time(f <- eval(bquote(ZIPLN(.(cs$formula), data = dd, zi = .(if (cs$zi == "covar") "single" else cs$zi), control = ctrl))))[["elapsed"]]
  data.frame(case = pk[1], accelerate = g$accelerate, ftol = g$ftol, gap = ref$ref_loglik[ref$case == pk[1]] - f$loglik,
             iterations = f$optim_par$iterations, time = round(t, 1), message = f$optim_par$message)
}, mc.cores = 12, mc.preschedule = FALSE)
print(Filter(is.character, out)); res <- do.call(rbind, Filter(is.data.frame, out))
res <- res[order(res$case, !res$accelerate, res$ftol), ]
options(width = 200); print(res, digits = 6, row.names = FALSE)
saveRDS(res, "inst/simus_ZIPLN/convergence/integrated.rds")
