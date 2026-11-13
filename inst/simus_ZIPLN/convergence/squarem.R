## Prototype: SQUAREM acceleration (Varadhan & Roland, 2008, scheme S3) of the outer loop of
## ZIPLN, builtin backend, non-sparse fits. The fixed-point map F is one outer iteration
## (M step then VE step); its state is (M, log S2, logit R on the zeros, B0).
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("inst/simus_ZIPLN/convergence/cases.R")

make_map <- function(cs, init) {
  ctrl <- ZIPLN_param(trace = 0, init = init, config_optim = list(maxit_out = 0))
  dz <- extract_model_zi(call("ZIPLN", formula = cs$formula, data = quote(dd)), list2env(list(dd = cs$data)))
  ctrl$ziparam <- if (dz$zicovar) "covar" else cs$zi
  fit <- ZIPLNfit$new(dz, ctrl)
  pv <- fit$.__enclos_env__$private
  nrm <- normalize_covariates(dz$X); X <- nrm$X_sc
  zero <- dz$Y == 0
  ve_data <- list(Y = dz$Y, X = X, O = dz$O, w = rep(1, nrow(dz$Y)))
  cfg <- ctrl$config_optim
  F <- function(s) {
    B <- pv$optimizer$B(M = s$M, X = X)
    Omega <- pv$optimizer$Omega(M = s$M, X = X, B = B, S2 = s$S2)
    zp <- pv$optimizer$zi(R = s$R, init_B0 = s$B0, X0 = dz$X0, config = cfg)
    out <- do.call(pv$optimizer$MS, list(data = ve_data, params = list(M = s$M, S2 = s$S2, Pi = zp$Pi, B = B, Omega = Omega), config = cfg))
    B <- pv$optimizer$B(M = out$M, X = X)
    obj <- -sum(zipln_vloglik(dz$Y, X, dz$O, zp$Pi, Omega, B, out$R, out$M, out$S2))
    list(M = out$M, S2 = out$S2, R = out$R, B0 = zp$B0, objective = obj)
  }
  to_vec <- function(s) c(s$M, log(s$S2), qlogis(pmin(pmax(s$R[zero], 1e-12), 1 - 1e-12)), s$B0)
  from_vec <- function(x, template) {
    np <- length(template$M); out <- template
    out$M[] <- x[1:np]; out$S2[] <- exp(x[np + 1:np])
    out$R[] <- 0; out$R[zero] <- plogis(x[2 * np + seq_len(sum(zero))])
    if (length(template$B0)) out$B0[] <- tail(x, length(template$B0))
    out
  }
  s0 <- list(M = pv$M, S2 = pv$S2, R = pv$R, B0 = pv$B0)
  list(F = F, to_vec = to_vec, from_vec = from_vec, s0 = s0)
}

plain <- function(mp, ftol = 1e-6, maxit = 5000) {
  s <- mp$s0; obj <- Inf; k <- 0
  repeat { k <- k + 1; s <- mp$F(s); if (abs(obj - s$objective) <= ftol * abs(s$objective) || k >= maxit) break; obj <- s$objective }
  list(objective = s$objective, evals = k)
}
## SQUAREM with the adaptive bounds on the step length of Varadhan & Roland (step.max starts
## at 1 and is multiplied by 4 each time it is reached)
squarem <- function(mp, ftol = 1e-6, maxit = 5000) {
  s <- mp$F(mp$s0); evals <- 1; obj <- s$objective; step_max <- 1
  repeat {
    s1 <- mp$F(s); s2 <- mp$F(s1); evals <- evals + 2
    x0 <- mp$to_vec(s); x1 <- mp$to_vec(s1); x2 <- mp$to_vec(s2)
    r <- x1 - x0; v <- x2 - x1 - r
    alpha <- max(1, min(step_max, sqrt(sum(r^2)) / sqrt(sum(v^2))))
    new <- s2
    if (alpha > 1) {
      cand <- tryCatch(mp$F(mp$from_vec(x0 + 2 * alpha * r + alpha^2 * v, s)), error = function(e) NULL); evals <- evals + 1
      if (!is.null(cand) && is.finite(cand$objective) && cand$objective <= s2$objective) {
        new <- cand
        if (alpha == step_max) step_max <- 4 * step_max
      } else step_max <- max(1, step_max / 4)
    } else if (alpha == step_max) step_max <- 4 * step_max
    done <- abs(obj - new$objective) <= ftol * abs(new$objective) || evals >= maxit
    s <- new; obj <- new$objective
    if (done) break
  }
  list(objective = obj, evals = evals)
}
pick <- list(c("barents col", "positive"), c("mollusk col", "lm"), c("mollusk site|site", "lm"), c("trichoptera col", "lm"),
             c("microcosm site col", "lm"), c("microcosm site|site", "positive"), c("oaks col", "lm"), c("barents single", "lm"))
out <- parallel::mclapply(pick, function(pk) {
  RhpcBLASctl::blas_set_num_threads(1)
  cs <- Filter(function(x) x$name == pk[1], cases)[[1]]
  mp <- make_map(cs, pk[2])
  tp <- system.time(p <- plain(mp))[["elapsed"]]
  tr <- system.time(ref <- plain(mp, ftol = 1e-10, maxit = 20000))[["elapsed"]]
  ts <- system.time(q <- squarem(mp))[["elapsed"]]
  ts2 <- system.time(q2 <- squarem(mp, ftol = 1e-9))[["elapsed"]]
  data.frame(case = pk[1], ref_loglik = -ref$objective, ref_evals = ref$evals, ref_time = round(tr, 1),
             plain_gap = -ref$objective + p$objective, plain_evals = p$evals, plain_time = round(tp, 1),
             sq6_gap = -ref$objective + q$objective, sq6_evals = q$evals, sq6_time = round(ts, 1),
             sq9_gap = -ref$objective + q2$objective, sq9_evals = q2$evals, sq9_time = round(ts2, 1))
}, mc.cores = 8, mc.preschedule = FALSE)
res <- do.call(rbind, Filter(is.data.frame, out)); print(Filter(is.character, out))
options(width = 200); print(res, digits = 8, row.names = FALSE)
saveRDS(res, "inst/simus_ZIPLN/convergence/squarem_adaptive.rds")
