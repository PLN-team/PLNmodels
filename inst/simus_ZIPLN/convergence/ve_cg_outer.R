## Prototype: the truncated Newton-CG VE step of ve_cg.R plugged into the outer loop of ZIPLN
## (builtin, non-sparse), with and without SQUAREM. Compared with the current VE step on the
## number of outer iterations, the number of Omega-products in the VE steps and the bound reached
## (stopping rule of the package: relative change 1e-6, at most 500 VE steps).
suppressMessages(devtools::load_all(".", quiet = TRUE))
src <- readLines("inst/simus_ZIPLN/convergence/ve_cg.R")
eval(parse(text = src[1:(grep("^## 1. Frozen problem", src) - 1)]))

## same map as make_map(), VE step either the current C++ one or ve_cg(); counts the products
## (2 per Newton iteration for the current step: gradient and step, the line search being
## linear in them; 1 + CG iterations for ve_cg, the step being a combination of CG directions).
## ve_maxit and ve_ftol default to the package VE settings (maxeval, ftol_in).
with_counter <- function(mp, cg_max = NULL, ve_maxit = NULL, ve_ftol = NULL) {
  env <- environment(mp$F); counter <- new.env(); counter$products <- 0
  pv <- env$pv; X <- env$X; dz <- env$dz; cfg <- env$cfg; ve_data <- env$ve_data
  if (!is.null(ve_maxit)) cfg$maxeval <- ve_maxit
  if (!is.null(ve_ftol)) cfg$ftol_in <- ve_ftol
  F <- function(s) {
    B <- pv$optimizer$B(M = s$M, X = X)
    Omega <- pv$optimizer$Omega(M = s$M, X = X, B = B, S2 = s$S2)
    zp <- pv$optimizer$zi(R = s$R, init_B0 = s$B0, X0 = dz$X0, config = cfg)
    if (is.null(cg_max)) {
      out <- do.call(pv$optimizer$MS, list(data = ve_data, params = list(M = s$M, S2 = s$S2, Pi = zp$Pi, B = B, Omega = Omega), config = cfg))
      counter$products <- counter$products + 2 * length(out$monitoring$objective)
    } else {
      out <- ve_cg(dz$Y, X, dz$O, Omega, zp$Pi, s$M, s$S2, cg_max = cg_max, maxit = cfg$maxeval, ftol = cfg$ftol_in)
      counter$products <- counter$products + out$products
    }
    B <- pv$optimizer$B(M = out$M, X = X)
    obj <- -sum(zipln_vloglik(dz$Y, X, dz$O, zp$Pi, Omega, B, out$R, out$M, out$S2))
    list(M = out$M, S2 = out$S2, R = out$R, B0 = zp$B0, objective = obj)
  }
  list(F = F, to_vec = mp$to_vec, from_vec = mp$from_vec, s0 = mp$s0, counter = counter)
}

dir <- "inst/simus_ZIPLN/convergence/ve_cg"; dir.create(dir, showWarnings = FALSE)
log_file <- file.path(dir, "outer_log.txt")
say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n", file = log_file, append = TRUE)
pick <- list(c("barents col", "positive"), c("mollusk col", "lm"), c("trichoptera col", "lm"), c("oaks col", "lm"),
             c("microcosm site col", "lm"), c("microcosm site|site", "positive"), c("oaks tree|tree", "lm"))
variants <- list(current_plain = list(NULL, plain), current_squarem = list(NULL, squarem),
                 cg20_plain = list(20, plain), cg20_squarem = list(20, squarem), cg5_squarem = list(5, squarem))
out <- parallel::mclapply(pick, function(pk) {
  RhpcBLASctl::blas_set_num_threads(1)
  cs <- Filter(function(x) x$name == pk[1], cases)[[1]]
  tag <- gsub("[ |]", "_", pk[1])
  rows <- lapply(names(variants), function(vn) {
    vr <- variants[[vn]]
    mp <- with_counter(make_map(cs, pk[2]), vr[[1]])
    t <- system.time(r <- vr[[2]](mp, ftol = 1e-6, maxit = 500))[["elapsed"]]
    row <- data.frame(case = pk[1], variant = vn, loglik = -r$objective, outer = r$evals,
                      products = mp$counter$products, time = round(t, 1))
    say(pk[1], vn, ": loglik", round(row$loglik, 2), "outer", row$outer, "products", row$products, "time", row$time)
    row
  })
  res <- do.call(rbind, rows); res$gain <- round(res$loglik - res$loglik[1], 2)
  saveRDS(res, file.path(dir, paste0(tag, "_outer.rds")))
  res
}, mc.cores = 7, mc.preschedule = FALSE)
print(Filter(is.character, out)); res <- do.call(rbind, Filter(is.data.frame, out))
options(width = 200); print(res, row.names = FALSE)
saveRDS(res, "inst/simus_ZIPLN/convergence/ve_cg_outer.rds")
