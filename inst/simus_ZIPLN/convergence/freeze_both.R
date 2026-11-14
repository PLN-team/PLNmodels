## Omega AND pi frozen: the outer loop is then a sequence of VE steps at fixed model
## parameters (B is profiled inside the VE step). If the objective keeps decreasing slowly,
## the slow mode is inside the VE step itself (cell-wise Newton on a coupled problem).
suppressMessages(devtools::load_all(".", quiet = TRUE))
src <- readLines("inst/simus_ZIPLN/convergence/freeze.R")
eval(parse(text = src[1:(grep("^pick <- ", src) - 1)]))
pick <- list(c("barents col", "positive"), c("mollusk col", "lm"), c("oaks col", "lm"), c("trichoptera col", "lm"),
             c("microcosm site col", "lm"))
out <- parallel::mclapply(pick, function(pk) {
  RhpcBLASctl::blas_set_num_threads(1)
  cs <- Filter(function(x) x$name == pk[1], cases)[[1]]
  mp <- make_map(cs, pk[2]); env <- environment(mp$F); pv <- env$pv
  conv <- run(mp, ftol = 1e-8, maxit = 3000)
  s <- conv$state; B <- pv$optimizer$B(M = s$M, X = env$X)
  Omega <- pv$optimizer$Omega(M = s$M, X = env$X, B = B, S2 = s$S2)
  zp <- pv$optimizer$zi(R = s$R, init_B0 = s$B0, X0 = env$dz$X0, config = env$cfg)
  pv$optimizer$Omega <- function(...) Omega
  pv$optimizer$zi <- function(R, ...) list(Pi = zp$Pi, B0 = zp$B0)
  ## from the starting point, then from the converged state: VE steps only
  r0 <- run(mp, ftol = 1e-12, maxit = 3000)
  mp$s0 <- s; r1 <- run(mp, ftol = 1e-12, maxit = 3000)
  data.frame(case = pk[1], start = c("initial", "converged"), iterations = c(r0$iterations, r1$iterations),
             rate = round(c(r0$rate, r1$rate), 5), objective = c(r0$objective, r1$objective),
             change_from_converged = round(c(r0$objective, r1$objective) - conv$objective, 4))
}, mc.cores = 5, mc.preschedule = FALSE)
print(Filter(is.character, out)); res <- do.call(rbind, Filter(is.data.frame, out))
options(width = 200); print(res, row.names = FALSE, digits = 10)
saveRDS(res, "inst/simus_ZIPLN/convergence/freeze_both.rds")
