## Sequential fits with early stopping of the second start: stop it at iteration t when the
## geometric extrapolation of its recent decreases cannot bring it below the bound of the
## first start (plus a margin). Evaluated on the recorded trajectories.
sim  <- readRDS("inst/simus_ZIPLN/init/trajectories_sim.rds")
real <- readRDS("inst/simus_ZIPLN/init/trajectories_real.rds")
real <- lapply(split(real, vapply(real, `[[`, "", "case")), function(x) {
  out <- list(case = x[[1]]$case, group = "real"); for (y in x) out[[y$init]] <- y$fit; out })
pairs <- c(sim, real)
maxit <- 200
## iteration at which the second trajectory `obj` is stopped against the target `target`
stop_at <- function(obj, target, w, margin, safety, min_it) {
  for (t in seq_along(obj)) {
    if (t <= max(w, min_it)) next
    d_now <- obj[t - 1] - obj[t]; d_old <- obj[t - w] - obj[t - w + 1]
    if (d_now < 0 || d_old <= 0) next
    r <- (d_now / d_old)^(1 / (w - 1))
    if (!is.finite(r) || r >= 1) next
    rest <- maxit - t
    gain <- safety * d_now * r * (1 - r^rest) / (1 - r)     # geometric extrapolation, up to maxit
    if (obj[t] - gain > target + margin * abs(target)) return(t)
  }
  length(obj)
}
evaluate <- function(p, first, w, margin, safety, min_it) {
  second <- setdiff(c("lm", "positive"), first)
  f1 <- p[[first]]; f2 <- p[[second]]
  target <- -f1$loglik                      # objective (minimized) of the first fit
  t <- stop_at(f2$objective, target, w, margin, safety, min_it)
  stopped <- t < length(f2$objective)
  kept <- if (!stopped && f2$loglik > f1$loglik) f2$loglik else f1$loglik
  best <- max(f1$loglik, f2$loglik)
  data.frame(case = p$case, group = p$group, stopped = stopped, loss = best - kept,
             cost = length(f1$objective) + t, cost_both = length(f1$objective) + length(f2$objective))
}
grid <- expand.grid(first = c("lm", "positive"), w = c(5, 10), margin = c(0, 1e-4), safety = c(1, 3, 10), min_it = c(10, 30),
                    stringsAsFactors = FALSE)
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(g) {
  r <- do.call(rbind, lapply(pairs, evaluate, first = grid$first[g], w = grid$w[g], margin = grid$margin[g],
                             safety = grid$safety[g], min_it = grid$min_it[g]))
  real <- r[r$group == "real", ]
  cbind(grid[g, ], stopped = round(mean(r$stopped), 2), wrong = sum(r$loss > 1e-6), max_loss = round(max(r$loss), 1),
        mean_loss = round(mean(r$loss), 2), saved = round(1 - sum(r$cost) / sum(r$cost_both), 3),
        real_wrong = sum(real$loss > 1e-6), real_max_loss = round(max(real$loss), 1), real_saved = round(1 - sum(real$cost) / sum(real$cost_both), 3))
}))
options(width = 200)
print(res[order(res$wrong, -res$saved), ], row.names = FALSE)
