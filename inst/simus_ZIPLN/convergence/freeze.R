## Which coupling makes the outer loop slow? The loop of squarem.R (plain iterations) is run
## with Omega, or pi, frozen at its value at convergence: the rate of the decrease of the
## objective, and the iterations to a relative change of 1e-8, tell which one is the slow mode.
suppressMessages(devtools::load_all(".", quiet = TRUE))
src <- readLines("inst/simus_ZIPLN/convergence/squarem.R")
eval(parse(text = src[1:(grep("^pick <- ", src) - 1)]))
run <- function(mp, freeze = NULL, fixed = NULL, ftol = 1e-8, maxit = 6000) {
  pv <- environment(mp$F)$pv
  if (!is.null(freeze)) {
    orig <- pv$optimizer[[freeze]]
    pv$optimizer[[freeze]] <- if (freeze == "Omega") function(...) fixed$Omega
                              else function(R, ...) list(Pi = fixed$Pi, B0 = fixed$B0)
    on.exit(pv$optimizer[[freeze]] <- orig)
  }
  s <- mp$s0; obj <- numeric(0)
  repeat {
    s <- mp$F(s); obj <- c(obj, s$objective); k <- length(obj)
    if ((k > 1 && abs(obj[k - 1] - obj[k]) <= ftol * abs(obj[k])) || k >= maxit) break
  }
  d <- -diff(obj); k <- length(d)
  rate <- if (k > 30) exp(mean(diff(log(pmax(d[(k - 20):k], 1e-300))))) else NA
  list(objective = obj[length(obj)], iterations = length(obj), rate = rate, state = s)
}
pick <- list(c("barents col", "positive"), c("mollusk col", "lm"), c("oaks col", "lm"), c("trichoptera col", "lm"),
             c("microcosm site col", "lm"))
out <- parallel::mclapply(pick, function(pk) {
  RhpcBLASctl::blas_set_num_threads(1)
  cs <- Filter(function(x) x$name == pk[1], cases)[[1]]
  mp <- make_map(cs, pk[2]); env <- environment(mp$F); pv <- env$pv
  ref <- run(mp, ftol = 1e-10, maxit = 20000)
  ## Omega, Pi and B0 at convergence
  s <- ref$state; B <- pv$optimizer$B(M = s$M, X = env$X)
  fixed <- list(Omega = pv$optimizer$Omega(M = s$M, X = env$X, B = B, S2 = s$S2))
  zp <- pv$optimizer$zi(R = s$R, init_B0 = s$B0, X0 = env$dz$X0, config = env$cfg)
  fixed$Pi <- zp$Pi; fixed$B0 <- zp$B0
  r <- list(plain = run(mp), Omega_frozen = run(mp, "Omega", fixed), pi_frozen = run(mp, "zi", fixed))
  data.frame(case = pk[1], variant = names(r), iterations = sapply(r, `[[`, "iterations"),
             rate = round(sapply(r, `[[`, "rate"), 5), gap = round(sapply(r, `[[`, "objective") - ref$objective, 3),
             ref_iterations = ref$iterations)
}, mc.cores = 5, mc.preschedule = FALSE)
print(Filter(is.character, out)); res <- do.call(rbind, Filter(is.data.frame, out))
options(width = 200); print(res, row.names = FALSE)
saveRDS(res, "inst/simus_ZIPLN/convergence/freeze.rds")
