## Prototype: VE step on (M, log S2) with B and Omega profiled at every evaluation (as the
## default nlopt path of PLN), by L-BFGS; the outer loop only updates R and pi (closed forms,
## "single" or "col" zero inflation). Compared with the reference of squarem.R.
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("inst/simus_ZIPLN/convergence/cases.R")
profiled_fit <- function(cs, init, ftol = 1e-9, maxit_out = 2000, maxit_ve = 200, max_time = 1800) {
  dz <- extract_model_zi(call("ZIPLN", formula = cs$formula, data = quote(dd)), list2env(list(dd = cs$data)))
  Y <- dz$Y; X <- normalize_covariates(dz$X)$X_sc; O <- dz$O; n <- nrow(Y); p <- ncol(Y)
  zero <- Y == 0; lfY <- .logfactorial(Y)
  P_X <- solve(crossprod(X), t(X))
  sp <- compute_ZIPLN_starting_point(Y, dz$X, dz$X0, O, method = init)
  M <- sp$M; S2 <- sp$S2; R <- sp$R
  zi_pi <- function(R) if (cs$zi == "single") matrix(mean(R), n, p) else matrix(colMeans(R), n, p, byrow = TRUE)
  Pi <- zi_pi(R)
  ## profiled bound and its gradient in (M, psi = log S2), R and Pi fixed
  fg <- function(x, R, Pi, grad = FALSE) {
    M <- matrix(x[1:(n * p)], n, p); psi <- matrix(x[n * p + 1:(n * p)], n, p); S2 <- exp(psi)
    Res <- M - X %*% (P_X %*% M)
    Sigma <- (crossprod(Res) + diag(colSums(S2), p)) / n
    ch <- chol(Sigma); Omega <- chol2inv(ch)
    A <- exp(O + M + S2 / 2); Q <- 1 - R
    val <- sum(Q * (Y * (O + M) - A - lfY)) + .5 * sum(psi) - n * sum(log(diag(ch))) + n * p / 2 -
      n * p / 2 +                                   # quadratic term at the profiled Omega is -n p / 2
      sum(R * qlogis(pmin(pmax(Pi, 1e-12), 1 - 1e-12)) - log1p(exp(qlogis(pmin(pmax(Pi, 1e-12), 1 - 1e-12))))) -
      sum(ifelse(R > 0, R * log(R), 0) + ifelse(Q > 0, Q * log(Q), 0))
    if (!grad) return(-val)
    gM <- Q * (Y - A) - Res %*% Omega
    gpsi <- .5 - .5 * S2 * (Q * A + matrix(diag(Omega), n, p, byrow = TRUE))
    -c(gM, gpsi)
  }
  ## bound computed by the package itself, so that the constants match the reference
  pkg_loglik <- function(M, S2, R, Pi) {
    Res <- M - X %*% (P_X %*% M)
    Omega <- solve((crossprod(Res) + diag(colSums(S2), p)) / n)
    sum(zipln_vloglik(Y, X, O, Pi, Omega, P_X %*% M, R, M, S2))
  }
  obj <- c(); traj <- NULL; it <- 0; t0 <- proc.time()[[3]]
  repeat {
    it <- it + 1
    x <- optim(c(M, log(S2)), function(x) fg(x, R, Pi), function(x) fg(x, R, Pi, TRUE), method = "L-BFGS-B",
               control = list(maxit = maxit_ve, factr = 1e3))$par
    M <- matrix(x[1:(n * p)], n, p); S2 <- matrix(exp(x[n * p + 1:(n * p)]), n, p)
    A <- exp(O + M + S2 / 2)
    R <- plogis(A + qlogis(pmin(pmax(Pi, 1e-12), 1 - 1e-12))) * zero
    Pi <- zi_pi(R)
    obj[it] <- fg(c(M, log(S2)), R, Pi)
    traj <- rbind(traj, c(time = proc.time()[[3]] - t0, loglik = pkg_loglik(M, S2, R, Pi)))
    if (it > 1 && abs(obj[it - 1] - obj[it]) <= ftol * abs(obj[it]) || it >= maxit_out ||
        traj[it, "time"] > max_time) break
  }
  list(loglik = traj[it, "loglik"], iterations = it, time = traj[it, "time"], trajectory = traj,
       capped = traj[it, "time"] > max_time)
}
ref <- readRDS("inst/simus_ZIPLN/convergence/squarem_adaptive.rds")
pick <- list(c("barents col", "positive"), c("barents single", "lm"), c("mollusk col", "lm"), c("trichoptera col", "lm"),
             c("microcosm site col", "lm"), c("oaks col", "lm"))
out <- parallel::mclapply(pick, function(pk) {
  RhpcBLASctl::blas_set_num_threads(1)
  cs <- Filter(function(x) x$name == pk[1], cases)[[1]]
  f <- profiled_fit(cs, pk[2])
  r <- ref[ref$case == pk[1], ]
  ## gap at the time the default fit stops, to compare at equal cost
  at_default <- f$trajectory[max(1, sum(f$trajectory[, "time"] <= r$plain_time)), "loglik"]
  row <- data.frame(case = pk[1], ref_loglik = r$ref_loglik, ref_time = r$ref_time, default_gap = r$plain_gap,
             default_time = r$plain_time, profiled_gap_at_default_time = r$ref_loglik - at_default,
             profiled_gap = r$ref_loglik - f$loglik, profiled_outer = f$iterations,
             profiled_time = round(f$time, 1), capped = f$capped)
  saveRDS(list(row = row, trajectory = f$trajectory),
          sprintf("inst/simus_ZIPLN/convergence/profiled_%s.rds", gsub("[^a-z]+", "_", pk[1])))
  cat(format(Sys.time()), pk[1], "done\n", file = "inst/simus_ZIPLN/convergence/profiled.log", append = TRUE)
  row
}, mc.cores = 6, mc.preschedule = FALSE)
print(Filter(is.character, out)); res <- do.call(rbind, Filter(is.data.frame, out))
options(width = 200); print(res, digits = 8, row.names = FALSE)
saveRDS(res, "inst/simus_ZIPLN/convergence/profiled.rds")
