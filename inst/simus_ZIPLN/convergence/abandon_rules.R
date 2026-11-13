## Offline replay of abandonment rules on the accelerated trajectories (both orders of the
## two starting points): misses (abandoned although it would have won) and VE steps saved.
suppressMessages(devtools::load_all(".", quiet = TRUE))
tr <- readRDS("inst/simus_ZIPLN/convergence/abandon_trajectories.rds")
maxit <- 500
## current rule, budget counted in cycles
rule_current <- function(safety) function(obj, target, ve_per_cycle)
  cannot_beat(obj, target, length(obj) + max(maxit - length(obj) * ve_per_cycle, 0) / ve_per_cycle, safety = safety)
## decrease over blocks of w cycles: geometric across blocks, robust to the jumps of SQUAREM
rule_block <- function(w, safety, min_iter = 30) function(obj, target, ve_per_cycle) {
  t <- length(obj)
  if (t <= max(2 * w, min_iter)) return(FALSE)
  D_now <- obj[t - w] - obj[t]; D_old <- obj[t - 2 * w] - obj[t - w]
  if (D_now <= 0 || D_old <= 0) return(FALSE)
  rho <- D_now / D_old
  if (rho >= 1) return(FALSE)
  blocks_left <- max(maxit - t * ve_per_cycle, 0) / ve_per_cycle / w
  obj[t] - safety * D_now * rho * (1 - rho^blocks_left) / (1 - rho) > target
}
rules <- list(current_s3 = rule_current(3), current_s10 = rule_current(10),
              block5_s3 = rule_block(5, 3), block10_s3 = rule_block(10, 3), block10_s2 = rule_block(10, 2))
cases <- unique(sapply(tr, `[[`, "case"))
res <- do.call(rbind, lapply(cases, function(cn) {
  two <- Filter(function(o) o$case == cn, tr); names(two) <- sapply(two, `[[`, "init")
  do.call(rbind, lapply(list(c("lm", "positive"), c("positive", "lm")), function(ord) {
    first <- two[[ord[1]]]; second <- two[[ord[2]]]
    target <- min(first$objective); obj <- second$objective; vpc <- second$VE / length(obj)
    wins <- min(obj) < target
    do.call(rbind, lapply(names(rules), function(nm) {
      t <- 1; while (t < length(obj) && !rules[[nm]](obj[1:t], target, vpc)) t <- t + 1
      ab <- t < length(obj)
      data.frame(case = cn, second = ord[2], wins = wins, rule = nm, abandoned = ab,
                 missed = if (ab && wins) round(target - min(obj), 3) else 0,
                 saved = if (ab) round(second$VE - t * vpc) else 0)
    }))
  }))
}))
summ <- aggregate(cbind(abandoned, misses = missed > 0, missed, saved) ~ rule, res, sum)
summ$total_VE_two_fits <- sum(sapply(tr, `[[`, "VE"))
print(summ)
cat("\nmisses:\n"); print(res[res$missed > 0, ], row.names = FALSE)
saveRDS(res, "inst/simus_ZIPLN/convergence/abandon_rules.rds")
