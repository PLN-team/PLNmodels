## Greedy search from fits run to convergence (maxit_out raised), so that the gains are
## those of the switches only, not of continuing an unconverged optimization.
MAXIT <- as.integer(Sys.getenv("MAXIT", "5000"))
src <- readLines("inst/simus_ZIPLN/species_regime/flips.R")
eval(parse(text = src[1:(grep("^run <- function", src) - 1)]))
ctrl <- function(...) ZIPLN_param(trace = 0, config_optim = list(maxit_out = MAXIT), ...)
refit <- function(cs, inception) suppressWarnings(suppressMessages(
  ZIPLN(cs$formula, data = cs$data, zi = cs$zi, control = ctrl(inception = inception))))
K <- as.integer(Sys.getenv("K", "10"))
run <- function(k) {
  RhpcBLASctl::blas_set_num_threads(1)
  cs <- cases[[k]]
  dz <- extract_model_zi(call("ZIPLN", formula = cs$formula, data = quote(dd)), list2env(list(dd = cs$data)))
  t0 <- system.time(f <- suppressWarnings(suppressMessages(ZIPLN(cs$formula, data = cs$data, zi = cs$zi, control = ctrl()))))[["elapsed"]]
  noop <- refit(cs, f)
  g <- species_gap(f, dz)
  cand <- rbind(data.frame(j = order(g$gap, decreasing = TRUE)[1:min(K, ncol(dz$Y))], start = "positive"),
                data.frame(j = order(g$gap)[1:min(3, ncol(dz$Y))], start = "lm"))
  cand <- cand[(cand$start == "positive" & g$gap[cand$j] > 0) | (cand$start == "lm" & g$gap[cand$j] < 0), ]
  current <- f; log <- list(); t_search <- 0
  for (r in seq_len(nrow(cand))) {
    t <- system.time(h <- tryCatch(refit(cs, switch_species(current, cand$j[r], g[[cand$start[r]]], dz$Y)), error = function(e) NULL))[["elapsed"]]
    t_search <- t_search + t
    gain <- if (is.null(h)) NA else h$loglik - current$loglik
    accepted <- isTRUE(gain > 1e-2)
    if (accepted) current <- h
    log[[r]] <- data.frame(rank = r, species = cand$j[r], to = cand$start[r], gap = g$gap[cand$j[r]], gain = gain,
                           accepted = accepted, time = t, iterations = if (is.null(h)) NA else h$optim_par$iterations)
  }
  log <- do.call(rbind, log)
  Pi0 <- colMeans(f$model_par$Pi); Pi1 <- colMeans(current$model_par$Pi)
  v0 <- diag(as.matrix(f$model_par$Sigma)); v1 <- diag(as.matrix(current$model_par$Sigma))
  list(summary = data.frame(case = cs$name, kept = f$optim_par$init, p = ncol(dz$Y), iterations_fit = f$optim_par$iterations,
                            stop = f$optim_par$message, time_fit = t0, noop_gain = noop$loglik - f$loglik, tries = nrow(log),
                            accepted = sum(log$accepted), total_gain = current$loglik - f$loglik, time_search = t_search,
                            mean_pi = sprintf("%.3f -> %.3f", mean(Pi0), mean(Pi1)),
                            median_var = sprintf("%.2f -> %.2f", median(v0), median(v1))),
       log = cbind(case = cs$name, log))
}
out <- parallel::mclapply(seq_along(cases), function(k) tryCatch(run(k), error = function(e) paste(cases[[k]]$name, conditionMessage(e))),
                          mc.cores = 13, mc.preschedule = FALSE)
err <- Filter(is.character, out); if (length(err)) print(unlist(err))
out <- Filter(is.list, out)
saveRDS(out, "inst/simus_ZIPLN/species_regime/greedy_converged.rds")
options(width = 220)
print(do.call(rbind, lapply(out, `[[`, "summary")), digits = 4, row.names = FALSE)
