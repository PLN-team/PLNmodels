## Cost of fitting the starting points in a given order, the second one being abandoned by
## cannot_beat() as in fit_ZIPLN_multistart(). Orders: fixed, oracle, and rules on the indicators.
suppressMessages(devtools::load_all(".", quiet = TRUE))
rows <- readRDS("inst/simus_ZIPLN/init/order_features.rds")
feat <- do.call(rbind, rows)
traj <- lapply(rows, attr, "traj")
feat$zi <- ifelse(feat$group == "real", sub("^.* ", "", feat$case), "col")
feat$zi[grepl("\\|", feat$zi)] <- "covar"
feat$winner <- ifelse(feat$loglik_pos > feat$loglik_lm, "positive", "lm")
## iterations needed with `first` fitted first; is the result the best of the two?
run_order <- function(i, first) {
  second <- setdiff(c("lm", "positive"), first)
  ll <- c(lm = feat$loglik_lm[i], positive = feat$loglik_pos[i])
  o1 <- traj[[i]][[first]]; o2 <- traj[[i]][[second]]
  target <- -ll[[first]]; t2 <- length(o2)
  for (t in seq_along(o2)) if (cannot_beat(o2[1:t], target, 200)) { t2 <- t; break }
  abandoned <- t2 < length(o2)
  kept <- if (!abandoned && ll[[second]] > ll[[first]]) ll[[second]] else ll[[first]]
  c(cost = length(o1) + t2, loss = max(ll) - kept)
}
cost_of <- function(firsts) {
  r <- t(vapply(seq_len(nrow(feat)), function(i) run_order(i, firsts[i]), numeric(2)))
  both <- vapply(traj, function(x) length(x$lm) + length(x$positive), numeric(1))
  c(saved = 1 - sum(r[, "cost"]) / sum(both), wrong = sum(r[, "loss"] > 1e-6), max_loss = max(r[, "loss"]),
    saved_real = 1 - sum(r[feat$group == "real", "cost"]) / sum(both[feat$group == "real"]),
    accuracy = mean(firsts == feat$winner))
}
## bound after the first iteration of each start (cost: one extra iteration)
first_it <- vapply(traj, function(x) if (x$positive[1] < x$lm[1]) "positive" else "lm", "")
strategies <- list(
  "lm first (PR #199)"                 = rep("lm", nrow(feat)),
  "positive first"                     = rep("positive", nrow(feat)),
  "oracle"                             = feat$winner,
  "after 1 iteration of each"          = first_it,
  "positive if zi != single"           = ifelse(feat$zi == "single", "lm", "positive"))
for (thr in c(0.05, 0.1, 0.15, 0.2))
  strategies[[sprintf("positive if p0 < %.2f", thr)]] <- ifelse(feat$p0 < thr, "positive", "lm")
for (thr in c(0.05, 0.1, 0.15))
  strategies[[sprintf("positive if p0 < %.2f and zi != single", thr)]] <- ifelse(feat$p0 < thr & feat$zi != "single", "positive", "lm")
for (thr in c(1, 1.5, 2))
  strategies[[sprintf("positive if log_mean > %.1f and zi != single", thr)]] <- ifelse(feat$log_mean > thr & feat$zi != "single", "positive", "lm")
res <- do.call(rbind, lapply(names(strategies), function(s) data.frame(strategy = s, t(cost_of(strategies[[s]])))))
res$extra <- ifelse(res$strategy == "after 1 iteration of each", "one more iteration", "")
options(width = 200)
print(res, digits = 3, row.names = FALSE)
cat("\npositive wins:", sum(feat$winner == "positive"), "of", nrow(feat), "\n")
print(table(group = feat$group, winner = feat$winner))

## Headroom: how much of the loser's trajectory is run before it is abandoned (oracle order)?
both <- vapply(traj, function(x) length(x$lm) + length(x$positive), numeric(1))
win_len <- vapply(seq_len(nrow(feat)), function(i) length(traj[[i]][[feat$winner[i]]]), numeric(1))
lose_len <- both - win_len
oracle_cost <- vapply(seq_len(nrow(feat)), function(i) run_order(i, feat$winner[i])[["cost"]], numeric(1))
cat(sprintf("\nwinner alone would save %.3f; oracle order saves %.3f\n", 1 - sum(win_len) / sum(both), 1 - sum(oracle_cost) / sum(both)))
cat(sprintf("share of the loser's iterations run before abandonment (oracle order): %.2f; loser abandoned in %d of %d cases\n",
            sum(oracle_cost - win_len) / sum(lose_len), sum(oracle_cost - win_len < lose_len), nrow(feat)))
cat("iterations of the loser: median", median(lose_len), "; before abandonment: median", median(oracle_cost - win_len), "\n")
