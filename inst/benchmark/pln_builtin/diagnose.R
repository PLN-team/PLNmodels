## Builtin PLN (full covariance): where is the slowness, in the VE steps (inner Newton
## iterations) or in the EM loop (Omega updates)? The EM loop is replayed one M step at a
## time (builtin_optimize_full with maxit_em = 1 restarts from the M step of the current
## M, S2, which continues the loop exactly); inner iterations and bound per EM iteration.
suppressMessages(devtools::load_all(".", quiet = TRUE)); RhpcBLASctl::blas_set_num_threads(1)
data(oaks); data(microcosm)
cases <- list(list(name = "oaks tree", data = oaks, formula = Abundance ~ 0 + tree + offset(log(Offset))),
              list(name = "microcosm", data = microcosm, formula = Abundance ~ 0 + site + offset(log(Offset))))
out <- parallel::mclapply(cases, function(cs) {
  RhpcBLASctl::blas_set_num_threads(1)
  args <- extract_model(call("PLN", formula = cs$formula, data = quote(dd)), list2env(list(dd = cs$data)))
  ctrl <- PLN_param(trace = 0, backend = "builtin")
  fit <- PLNfit$new(args$Y, args$X, args$O, args$w, args$formula, ctrl)
  pv <- fit$.__enclos_env__$private
  nrm <- normalize_covariates(args$X)
  data <- list(Y = args$Y, X = nrm$X_sc, O = args$O, w = args$w)
  rows <- NULL
  for (shift in c(FALSE, TRUE)) {
    par <- list(B = sweep(pv$B, 1, nrm$scales, "*"), M = pv$M, S2 = pv$S2)
    cfg <- modifyList(ctrl$config_optim, list(maxit_em = 1L, objective_constant = shift))
    t0 <- proc.time()[["elapsed"]]
    for (k in 1:400) {
      o <- builtin_optimize_full(data, par, cfg)
      rows <- rbind(rows, data.frame(case = cs$name, shift = shift, em = k, inner = length(o$monitoring$objective),
                                     loglik = sum(o$Ji), time = proc.time()[["elapsed"]] - t0))
      par <- list(B = o$B, M = o$M, S2 = o$S2)
      if (proc.time()[["elapsed"]] - t0 > 600) break
    }
  }
  rows
}, mc.cores = 2)
res <- do.call(rbind, out); saveRDS(res, "inst/benchmark/pln_builtin/diagnose.rds")
for (cn in unique(res$case)) for (sh in c(FALSE, TRUE)) {
  r <- res[res$case == cn & res$shift == sh, ]; d <- diff(r$loglik); k <- length(d)
  rate <- if (k > 30) exp(mean(diff(log(pmax(d[(k - 20):k], 1e-300))))) else NA
  cat(cn, "shift", sh, ": EM", nrow(r), " inner per EM (first 5, median, last 5)", head(r$inner, 5), "|", median(r$inner), "|", tail(r$inner, 5),
      " loglik first/10/50/last", round(r$loglik[c(1, 10, 50, nrow(r))], 2), " EM rate", round(rate, 4), " time", round(tail(r$time, 1)), "\n")
}
