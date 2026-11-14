## Follow-up of ve_cg_inexact.R: tighter VE tolerances (1e-7, 3e-7) for the inexact Newton-CG VE step.
## Inexact Newton-CG VE steps in the outer loop of ZIPLN (with SQUAREM): the VE step stops at a
## looser relative tolerance, or after a few Newton steps, instead of being solved to 1e-8.
## Compared with the current default (current VE step + SQUAREM) and the exact Newton-CG VE step,
## on the bound reached, the outer iterations and the Omega-products, under the package stopping
## rule (relative change 1e-6, at most 500 VE steps).
suppressMessages(devtools::load_all(".", quiet = TRUE))
src <- readLines("inst/simus_ZIPLN/convergence/ve_cg_outer.R")
eval(parse(text = src[1:(grep("^dir <- ", src) - 1)]))

dir <- "inst/simus_ZIPLN/convergence/ve_cg"; dir.create(dir, showWarnings = FALSE)
log_file <- file.path(dir, "ftol7_log.txt")
say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n", file = log_file, append = TRUE)
pick <- list(c("barents col", "positive"), c("mollusk col", "lm"), c("trichoptera col", "lm"), c("oaks col", "lm"),
             c("microcosm site col", "lm"), c("microcosm site|site", "positive"), c("oaks tree|tree", "lm"))
variants <- list(current     = list(cg = NULL),
                 cg5_ftol7   = list(cg = 5,  ftol = 1e-7),
                 cg5_ftol3e7 = list(cg = 5,  ftol = 3e-7),
                 cg10_ftol7  = list(cg = 10, ftol = 1e-7))
jobs <- expand.grid(variant = names(variants), case = seq_along(pick), stringsAsFactors = FALSE)
## heaviest jobs first
jobs <- jobs[order(-grepl("microcosm", sapply(pick, `[`, 1))[jobs$case], grepl("current", jobs$variant)), ]
out <- parallel::mclapply(seq_len(nrow(jobs)), function(j) {
  RhpcBLASctl::blas_set_num_threads(1)
  pk <- pick[[jobs$case[j]]]; vn <- jobs$variant[j]; vr <- variants[[vn]]
  cs <- Filter(function(x) x$name == pk[1], cases)[[1]]
  mp <- with_counter(make_map(cs, pk[2]), vr$cg, vr$maxit, vr$ftol)
  t <- system.time(r <- squarem(mp, ftol = 1e-6, maxit = 500))[["elapsed"]]
  row <- data.frame(case = pk[1], variant = vn, loglik = -r$objective, outer = r$evals,
                    products = mp$counter$products, time = round(t, 1))
  say(pk[1], vn, ": loglik", round(row$loglik, 2), "outer", row$outer, "products", row$products, "time", row$time)
  saveRDS(row, file.path(dir, paste0("ftol7_", gsub("[ |]", "_", pk[1]), "_", vn, ".rds")))
  row
}, mc.cores = 16, mc.preschedule = FALSE)
print(Filter(is.character, out)); res <- do.call(rbind, Filter(is.data.frame, out))
base <- res[res$variant == "current", ]
res$gain <- round(res$loglik - base$loglik[match(res$case, base$case)], 2)
res$outer_ratio <- round(res$outer / base$outer[match(res$case, base$case)], 2)
res$prod_ratio <- round(res$products / base$products[match(res$case, base$case)], 2)
options(width = 200); print(res[order(res$case, res$variant), ], row.names = FALSE)
saveRDS(res, "inst/simus_ZIPLN/convergence/ve_cg_ftol7.rds")
