## Race between the two starting points: run both for k iterations, keep the one with the
## best objective, continue it. Accuracy, loss and cost, as a function of k.
sim  <- readRDS("inst/simus_ZIPLN/init/trajectories_sim.rds")
real <- readRDS("inst/simus_ZIPLN/init/trajectories_real.rds")
## real data: pair the two starts of each case
real <- lapply(split(real, vapply(real, `[[`, "", "case")), function(x) {
  out <- list(case = x[[1]]$case, group = "real")
  for (y in x) out[[y$init]] <- y$fit
  out
})
pairs <- c(sim, real)
## objective (minimized) after k iterations, the last value if the fit stopped before
at <- function(f, k) f$objective[min(k, length(f$objective))]
race <- function(p, k, margin = 0) {
  a <- at(p$lm, k); b <- at(p$positive, k)
  final_best <- if (p$positive$loglik > p$lm$loglik) "positive" else "lm"
  ## decided at k only if the gap exceeds the margin (relative); otherwise both are run to the end
  decided <- abs(a - b) > margin * abs(min(a, b))
  pick <- if (!decided) final_best else if (b < a) "positive" else "lm"
  n <- c(lm = length(p$lm$objective), positive = length(p$positive$objective))
  cost <- if (!decided) sum(n) else min(k, n[["lm"]]) + min(k, n[["positive"]]) + max(0, n[[pick]] - k)
  data.frame(case = p$case, group = p$group, k = k, margin = margin, decided = decided,
             correct = pick == final_best, loss = max(p$lm$loglik, p$positive$loglik) - p[[pick]]$loglik,
             cost = cost, cost_both = sum(n), cost_lm = n[["lm"]])
}
grid <- expand.grid(k = c(1, 2, 3, 5, 10, 20, 50), margin = c(0, 1e-4, 1e-3))
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(g) do.call(rbind, lapply(pairs, race, k = grid$k[g], margin = grid$margin[g]))))
summ <- function(d) data.frame(n = nrow(d), decided = mean(d$decided), wrong = sum(!d$correct),
                               mean_loss = mean(d$loss), max_loss = max(d$loss),
                               cost_vs_both = sum(d$cost) / sum(d$cost_both), cost_vs_lm = sum(d$cost) / sum(d$cost_lm))
options(width = 200)
for (grp in c("all", "real")) {
  cat("\n==", grp, "\n")
  d <- if (grp == "all") res else res[res$group == "real", ]
  print(do.call(rbind, lapply(split(d, list(d$k, d$margin)), function(x) cbind(k = x$k[1], margin = x$margin[1], summ(x)))), row.names = FALSE, digits = 3)
}
cat("\nwrong decisions at k = 5, margin 0:\n"); print(res[res$k == 5 & res$margin == 0 & !res$correct, c("case", "loss")], row.names = FALSE)
cat("\niterations per start (median, max):", median(unlist(lapply(pairs, function(p) c(length(p$lm$objective), length(p$positive$objective))))),
    max(unlist(lapply(pairs, function(p) c(length(p$lm$objective), length(p$positive$objective))))), "\n")
saveRDS(res, "inst/simus_ZIPLN/init/race.rds")
