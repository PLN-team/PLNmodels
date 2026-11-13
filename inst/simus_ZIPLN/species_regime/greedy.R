## Greedy search over the regimes of the species: from the fit of ZIPLN(), switch one species
## at a time, in decreasing order of the gap indicator (or increasing order, towards the
## variance regime), refit, keep the switch only if the bound increases.
src <- readLines("inst/simus_ZIPLN/species_regime/flips.R")
eval(parse(text = src[1:(grep("^run <- function", src) - 1)]))
K <- as.integer(Sys.getenv("K", "10"))
run <- function(k) {
  RhpcBLASctl::blas_set_num_threads(1)
  cs <- cases[[k]]
  dz <- extract_model_zi(call("ZIPLN", formula = cs$formula, data = quote(dd)), list2env(list(dd = cs$data)))
  t0 <- system.time(f <- suppressWarnings(suppressMessages(ZIPLN(cs$formula, data = cs$data, zi = cs$zi, control = ZIPLN_param(trace = 0)))))[["elapsed"]]
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
    log[[r]] <- data.frame(rank = r, species = cand$j[r], to = cand$start[r], gap = g$gap[cand$j[r]], gain = gain, accepted = accepted, time = t)
  }
  log <- do.call(rbind, log)
  list(summary = data.frame(case = cs$name, kept = f$optim_par$init, p = ncol(dz$Y), time_fit = t0, tries = nrow(log),
                            accepted = sum(log$accepted), total_gain = current$loglik - f$loglik, time_search = t_search,
                            first_accept_rank = if (any(log$accepted)) min(log$rank[log$accepted]) else NA),
       log = cbind(case = cs$name, log))
}
out <- parallel::mclapply(seq_along(cases), function(k) tryCatch(run(k), error = function(e) paste(cases[[k]]$name, conditionMessage(e))),
                          mc.cores = 13, mc.preschedule = FALSE)
err <- Filter(is.character, out); if (length(err)) print(unlist(err))
out <- Filter(is.list, out)
saveRDS(out, "inst/simus_ZIPLN/species_regime/greedy.rds")
options(width = 200)
print(do.call(rbind, lapply(out, `[[`, "summary")), digits = 4, row.names = FALSE)
lg <- do.call(rbind, lapply(out, `[[`, "log"))
cat("\naccepted switches:\n"); print(lg[lg$accepted, ], digits = 3, row.names = FALSE)
