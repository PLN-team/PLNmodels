## SQUAREM with a partial VE step (maxit_ve Newton iterations, as PLNnetwork) on the 13 real
## cases, single start (the one kept by init = "best"), maxit_out = 500 and 2000.
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("inst/simus_ZIPLN/convergence/cases.R")
kept <- readRDS("inst/simus_ZIPLN/convergence/multistart_check.rds")[, c("case", "kept")]
grid <- expand.grid(k = seq_along(cases), maxit_ve = c(0L, 1L, 3L), maxit_out = c(500L, 2000L))
grid <- grid[!(grid$maxit_ve == 0 & grid$maxit_out == 2000), ]
grid <- grid[order(-grepl("microcosm", sapply(cases, `[[`, "name"))[grid$k]), ]
out <- parallel::mclapply(seq_len(nrow(grid)), function(i) {
  RhpcBLASctl::blas_set_num_threads(1)
  g <- grid[i, ]; cs <- cases[[g$k]]; dd <- cs$data
  co <- list(maxit_out = g$maxit_out); if (g$maxit_ve > 0) co$maxit_ve <- g$maxit_ve
  ctrl <- ZIPLN_param(trace = 0, init = kept$kept[kept$case == cs$name], config_optim = co)
  t <- system.time(f <- suppressWarnings(eval(bquote(ZIPLN(.(cs$formula), data = dd,
         zi = .(if (cs$zi == "covar") "single" else cs$zi), control = ctrl)))))[["elapsed"]]
  data.frame(case = cs$name, maxit_ve = g$maxit_ve, maxit_out = g$maxit_out, loglik = f$loglik,
             time = round(t, 1), iterations = f$optim_par$iterations, message = f$optim_par$message)
}, mc.cores = 13, mc.preschedule = FALSE)
print(Filter(is.character, out)); res <- do.call(rbind, Filter(is.data.frame, out))
saveRDS(res, "inst/simus_ZIPLN/convergence/partial_ve.rds")
base <- res[res$maxit_ve == 0, c("case", "loglik", "time")]
res <- merge(res, base, by = "case", suffixes = c("", "_full"))
res$gain <- round(res$loglik - res$loglik_full, 2); res$time_ratio <- round(res$time / res$time_full, 2)
options(width = 200)
print(res[order(res$case, res$maxit_out, res$maxit_ve), c("case", "maxit_ve", "maxit_out", "gain", "time", "time_ratio", "iterations", "message")], row.names = FALSE)
