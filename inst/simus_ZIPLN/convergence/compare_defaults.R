## Candidate defaults (accelerated, ftol_out x maxit_out) replayed on the trajectories, against
## the current default (baseline.R). Gap to the best known optimum, and VE steps.
tr <- Filter(is.list, readRDS("inst/simus_ZIPLN/convergence/trajectories.rds"))
ref <- readRDS("inst/simus_ZIPLN/convergence/squarem_adaptive.rds")
base <- readRDS("inst/simus_ZIPLN/convergence/baseline.rds")
replay <- function(o, ftol, maxit) {
  obj <- o$objective; ve <- cumsum(o$evals); t <- 1
  while (t < length(obj) && ve[t] < maxit && !(t > 1 && abs(obj[t - 1] - obj[t]) <= ftol * abs(obj[t]))) t <- t + 1
  c(objective = obj[t], VE = ve[t])
}
cands <- expand.grid(ftol = c(1e-6, 1e-7), maxit = c(200, 500, 1000))
res <- do.call(rbind, lapply(tr, function(o) {
  b <- base[base$case == o$case, ]
  best <- min(o$objective, -ref$ref_loglik[ref$case == o$case], b$objective)
  row <- data.frame(case = o$case, default_gap = round(b$objective - best, 2), default_VE = b$VE)
  for (i in seq_len(nrow(cands))) {
    r <- replay(o, cands$ftol[i], cands$maxit[i]); nm <- sprintf("%g_%d", cands$ftol[i], cands$maxit[i])
    row[[paste0("gap_", nm)]] <- round(r[["objective"]] - best, 2); row[[paste0("VE_", nm)]] <- r[["VE"]]
  }
  row
}))
options(width = 250); print(res, row.names = FALSE)
cat("\ntotal VE: default", sum(res$default_VE), "| candidates:",
    paste(grep("^VE_", names(res), value = TRUE), colSums(res[grep("^VE_", names(res))]), collapse = "  "), "\n")
