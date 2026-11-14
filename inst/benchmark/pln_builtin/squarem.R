## Builtin PLN (full covariance): SQUAREM on the EM loop. The map F is one EM iteration
## (M step from (M, S2), VE step to ftol_in), run by the C++ optimizer with maxit_em = 1; its
## state is (M, log S2), B being profiled. The VE tolerance is on the bound with its constants.
## Stop: relative change of the bound (with its constants) below `tol`, at most `maxit` VE steps.
suppressMessages(devtools::load_all(".", quiet = TRUE)); RhpcBLASctl::blas_set_num_threads(1)
data(barents); data(oaks); data(mollusk); data(microcosm)
mol <- suppressMessages(prepare_data(mollusk$Abundance, mollusk$Covariate))
cases <- list(
  list(name = "barents",   data = barents,   formula = Abundance ~ 1 + offset(log(Offset))),
  list(name = "mollusk",   data = mol,       formula = Abundance ~ 1 + offset(log(Offset))),
  list(name = "oaks",      data = oaks,      formula = Abundance ~ 1 + offset(log(Offset))),
  list(name = "oaks tree", data = oaks,      formula = Abundance ~ 0 + tree + offset(log(Offset))),
  list(name = "microcosm", data = microcosm, formula = Abundance ~ 0 + site + offset(log(Offset))))

make_map <- function(cs, ve_tol) {
  args <- extract_model(call("PLN", formula = cs$formula, data = quote(dd)), list2env(list(dd = cs$data)))
  ctrl <- PLN_param(trace = 0, backend = "builtin")
  fit <- PLNfit$new(args$Y, args$X, args$O, args$w, args$formula, ctrl)
  pv <- fit$.__enclos_env__$private
  nrm <- normalize_covariates(args$X); X <- nrm$X_sc
  data <- list(Y = args$Y, X = X, O = args$O, w = args$w)
  P_X <- solve(crossprod(X, X * args$w), t(X * args$w))
  cfg <- modifyList(ctrl$config_optim, list(maxit_em = 1L, objective_constant = TRUE, ftol_in = ve_tol))
  counter <- new.env(); counter$inner <- 0
  F <- function(s) {
    o <- builtin_optimize_full(data, list(B = P_X %*% s$M, M = s$M, S2 = s$S2), cfg)
    counter$inner <- counter$inner + length(o$monitoring$objective)
    list(M = o$M, S2 = o$S2, objective = -sum(o$Ji))
  }
  np <- length(pv$M)
  list(F = F, to_vec = function(s) c(s$M, log(s$S2)),
       from_vec = function(x, s) { s$M[] <- x[1:np]; s$S2[] <- exp(x[np + 1:np]); s },
       s0 = list(M = pv$M, S2 = pv$S2), counter = counter)
}
plain <- function(mp, tol, maxit) {
  s <- mp$F(mp$s0); evals <- 1; obj <- s$objective
  repeat {
    s <- mp$F(s); evals <- evals + 1
    done <- abs(obj - s$objective) <= tol * abs(s$objective) || evals >= maxit
    obj <- s$objective; if (done) break
  }
  list(objective = obj, evals = evals)
}
squarem <- function(mp, tol, maxit) {
  s <- mp$F(mp$s0); evals <- 1; obj <- s$objective; step_max <- 1
  repeat {
    s1 <- mp$F(s); s2 <- mp$F(s1); evals <- evals + 2
    x0 <- mp$to_vec(s); r <- mp$to_vec(s1) - x0; v <- mp$to_vec(s2) - mp$to_vec(s1) - r
    alpha <- max(1, min(step_max, sqrt(sum(r^2) / sum(v^2))))
    new <- s2
    if (alpha > 1) {
      cand <- tryCatch(suppressWarnings(mp$F(mp$from_vec(x0 + 2 * alpha * r + alpha^2 * v, s))), error = function(e) NULL)
      evals <- evals + 1
      if (!is.null(cand) && is.finite(cand$objective) && cand$objective <= s2$objective) {
        new <- cand; if (alpha == step_max) step_max <- 4 * step_max
      } else step_max <- max(1, step_max / 4)
    } else if (alpha == step_max) step_max <- 4 * step_max
    done <- abs(obj - new$objective) <= tol * abs(new$objective) || evals >= maxit
    s <- new; obj <- new$objective
    if (done) break
  }
  list(objective = obj, evals = evals)
}
variants <- list(
  plain_ve6_stop7   = list(fun = plain,   ve = 1e-6, tol = 1e-7),
  squarem_ve6_stop7 = list(fun = squarem, ve = 1e-6, tol = 1e-7),
  squarem_ve6_stop8 = list(fun = squarem, ve = 1e-6, tol = 1e-8),
  squarem_ve8_stop8 = list(fun = squarem, ve = 1e-8, tol = 1e-8),
  squarem_ve4_stop7 = list(fun = squarem, ve = 1e-4, tol = 1e-7))
log_file <- "inst/benchmark/pln_builtin/squarem_log.txt"
out <- parallel::mclapply(cases, function(cs) {
  RhpcBLASctl::blas_set_num_threads(1)
  t <- system.time(def <- PLN(cs$formula, data = cs$data, control = PLN_param(trace = 0, backend = "builtin")))[["elapsed"]]
  rows <- data.frame(case = cs$name, variant = "default", loglik = def$loglik, evals = NA, inner = def$optim_par$iterations, time = round(t, 2))
  for (vn in names(variants)) {
    v <- variants[[vn]]; mp <- make_map(cs, v$ve)
    t <- system.time(r <- v$fun(mp, v$tol, 2000))[["elapsed"]]
    rows <- rbind(rows, data.frame(case = cs$name, variant = vn, loglik = -r$objective, evals = r$evals, inner = mp$counter$inner, time = round(t, 2)))
    cat(format(Sys.time(), "%H:%M:%S"), paste(tail(rows, 1), collapse = " "), "\n", file = log_file, append = TRUE)
  }
  rows
}, mc.cores = 5, mc.preschedule = FALSE)
print(Filter(is.character, out)); res <- do.call(rbind, Filter(is.data.frame, out))
res$gain <- round(res$loglik - res$loglik[res$variant == "default"][match(res$case, res$case[res$variant == "default"])], 2)
options(width = 200); print(res, row.names = FALSE)
saveRDS(res, "inst/benchmark/pln_builtin/squarem.rds")
