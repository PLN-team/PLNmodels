## Offline replay of stopping rules on the trajectories of trajectories.R: gap to the end of
## the trajectory (or to the reference optimum when better) and VE steps used.
tr <- Filter(is.list, readRDS("inst/simus_ZIPLN/convergence/trajectories.rds"))
ref <- readRDS("inst/simus_ZIPLN/convergence/squarem_adaptive.rds")
## relative decrease of the objective between two cycles
rule_rel <- function(ftol) function(obj, t) t > 1 && abs(obj[t - 1] - obj[t]) <= ftol * abs(obj[t])
## gain left, extrapolated from the geometric decrease of the objective over two windows
## of w cycles; relative tolerance; while the rate is not < 1, the relative rule at ftol
rule_extrap <- function(tol, w = 5, ftol = tol / 10) function(obj, t) {
  if (t > 1 && abs(obj[t - 1] - obj[t]) <= ftol * abs(obj[t])) return(TRUE)
  if (t <= 2 * w) return(FALSE)
  D_now <- obj[t - w] - obj[t]; D_old <- obj[t - 2 * w] - obj[t - w]
  if (D_now <= 0 || D_old <= 0) return(FALSE)
  rho <- D_now / D_old
  if (rho >= 1) return(FALSE)
  D_now * rho / (1 - rho) <= tol * abs(obj[t])
}
rules <- list(rel_1e6 = rule_rel(1e-6), rel_1e7 = rule_rel(1e-7), rel_1e8 = rule_rel(1e-8),
              ext_1e6 = rule_extrap(1e-6), ext_1e7 = rule_extrap(1e-7), ext_3e7 = rule_extrap(3e-7))
res <- do.call(rbind, lapply(tr, function(o) {
  obj <- o$objective; best <- min(obj, -ref$ref_loglik[ref$case == o$case])
  do.call(rbind, lapply(names(rules), function(nm) {
    t <- 1; while (t < length(obj) && !rules[[nm]](obj, t)) t <- t + 1
    data.frame(case = o$case, rule = nm, gap = obj[t] - best, VE = sum(o$evals[1:t]))
  }))
}))
options(width = 220)
print(reshape(res, idvar = "case", timevar = "rule", direction = "wide"), digits = 3, row.names = FALSE)
