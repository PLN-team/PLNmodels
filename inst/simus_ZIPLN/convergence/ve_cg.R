## Prototype: VE step of ZIPLN (full covariance, B profiled) by truncated Newton, with
## conjugate gradients preconditioned by the current cell-wise 2x2 Newton step (cg_max = 1
## gives back the current Jacobi-like step). Hessian-vector products include the coupling
## through Omega and through the projection on null(X').
suppressMessages(devtools::load_all(".", quiet = TRUE))
src <- readLines("inst/simus_ZIPLN/convergence/freeze.R")
eval(parse(text = src[1:(grep("^pick <- ", src) - 1)]))

ve_cg <- function(Y, X, O, Omega, Pi, M, S2, cg_max = 20, eta = 0.1, maxit = 200, ftol = 1e-8, gtol = 1e-6,
                  monitor = NULL, max_time = Inf) {
  t0 <- proc.time()[["elapsed"]]; mark <- 1; trace <- NULL
  n <- nrow(Y); p <- ncol(Y); zero <- Y == 0; lfY <- .logfactorial(Y)
  H <- X %*% solve(crossprod(X), t(X)); Qp <- function(V) V - H %*% V
  dO <- diag(Omega); lpi <- qlogis(pmin(pmax(Pi, 1e-12), 1 - 1e-12))
  psi <- log(S2); iter_cg <- 0; n_prod <- 0
  obj_fn <- function(M, psi, R) {        # negated bound, R fixed (terms depending on M, psi)
    S2 <- exp(psi); A <- exp(O + M + S2 / 2); QM <- Qp(M)
    -(sum((1 - R) * (Y * (O + M) - A - lfY)) - .5 * sum((QM %*% Omega) * QM) -
        .5 * sum(sweep(S2, 2, dO, "*")) + .5 * sum(psi))
  }
  hist <- numeric(0)
  for (it in seq_len(maxit)) {
    S2 <- exp(psi); A <- exp(O + M + S2 / 2)
    R <- plogis(A + lpi) * zero; Ae <- (1 - R) * A
    QMO <- Qp(M) %*% Omega; n_prod <- n_prod + 1
    gM <- -( (1 - R) * Y - Ae - Qp(QMO) )          # gradient of the negated bound
    gP <- -(.5 - .5 * S2 * (Ae + matrix(dO, n, p, byrow = TRUE)))
    f0 <- obj_fn(M, psi, R)
    ## cell-wise 2x2 blocks (as compute_joint_step_MS)
    hmm <- Ae + matrix(dO, n, p, byrow = TRUE); hmp <- .5 * Ae * S2
    hpp <- .5 * S2 * (Ae * (1 + .5 * S2) + matrix(dO, n, p, byrow = TRUE))
    det <- pmax(hmm * hpp - hmp^2, 1e-20)
    prec <- function(rM, rP) list(M = (hpp * rM - hmp * rP) / det, P = (hmm * rP - hmp * rM) / det)
    Hv <- function(vM, vP) { n_prod <<- n_prod + 1
      list(M = Ae * vM + Qp(Qp(vM) %*% Omega) + hmp * vP, P = hmp * vM + hpp * vP) }
    ## PCG on H d = g
    dM <- 0 * M; dP <- 0 * psi; rM <- gM; rP <- gP
    z <- prec(rM, rP); sM <- z$M; sP <- z$P; rz <- sum(rM * z$M) + sum(rP * z$P)
    g_norm <- sqrt(sum(gM^2) + sum(gP^2))
    for (k in seq_len(cg_max)) {
      h <- Hv(sM, sP); curv <- sum(sM * h$M) + sum(sP * h$P)
      if (curv <= 0) { if (k == 1) { dM <- sM; dP <- sP }; break }
      a <- rz / curv; dM <- dM + a * sM; dP <- dP + a * sP
      rM <- rM - a * h$M; rP <- rP - a * h$P; iter_cg <- iter_cg + 1
      if (sqrt(sum(rM^2) + sum(rP^2)) <= eta * g_norm) break
      z <- prec(rM, rP); rz_new <- sum(rM * z$M) + sum(rP * z$P)
      sM <- z$M + rz_new / rz * sM; sP <- z$P + rz_new / rz * sP; rz <- rz_new
    }
    ## Armijo backtracking
    slope <- -(sum(gM * dM) + sum(gP * dP)); alpha <- 1
    for (ls in 1:20) { if (obj_fn(M - alpha * dM, psi - alpha * dP, R) <= f0 + 1e-4 * alpha * slope) break; alpha <- alpha / 2 }
    M <- M - alpha * dM; psi <- psi - alpha * dP
    f1 <- obj_fn(M, psi, R); hist <- c(hist, f1)
    elapsed <- proc.time()[["elapsed"]] - t0
    if (!is.null(monitor) && n_prod >= mark) {
      trace <- rbind(trace, data.frame(newton = it, products = n_prod, time = elapsed, value = monitor(M, psi, lpi, zero)))
      mark <- 1.5 * mark
    }
    if (abs(f0 - f1) <= ftol * abs(f1) || g_norm < gtol || elapsed > max_time) break
  }
  S2 <- exp(psi); A <- exp(O + M + S2 / 2)
  list(M = M, S2 = S2, R = plogis(A + lpi) * zero, iterations = it, cg = iter_cg, products = n_prod, history = hist, trace = trace)
}

## 1. Frozen problem: from the state reached by the current loop, Omega and pi at their values
## there; solve the VE step with cg_max = 1 (current step), 5 and 20, to a tight tolerance,
## recording the gain along the way against the number of Omega-products. Start states are cached.
dir <- "inst/simus_ZIPLN/convergence/ve_cg"; dir.create(dir, showWarnings = FALSE)
log_file <- file.path(dir, "log.txt")
say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n", file = log_file, append = TRUE)
pick <- list(c("barents col", "positive"), c("mollusk col", "lm"), c("trichoptera col", "lm"), c("oaks col", "lm"),
             c("microcosm site col", "lm"))
out <- parallel::mclapply(pick, function(pk) {
  RhpcBLASctl::blas_set_num_threads(1)
  cs <- Filter(function(x) x$name == pk[1], cases)[[1]]
  tag <- gsub("[ |]", "_", pk[1])
  mp <- make_map(cs, pk[2]); env <- environment(mp$F); pv <- env$pv
  start_file <- file.path(dir, paste0(tag, "_start.rds"))
  if (file.exists(start_file)) conv <- readRDS(start_file) else {
    conv <- run(mp, ftol = 1e-6, maxit = 3000)        # about where the default stops
    saveRDS(conv, start_file)
  }
  say(pk[1], "start state ready")
  s <- conv$state; B <- pv$optimizer$B(M = s$M, X = env$X)
  Omega <- pv$optimizer$Omega(M = s$M, X = env$X, B = B, S2 = s$S2)
  zp <- pv$optimizer$zi(R = s$R, init_B0 = s$B0, X0 = env$dz$X0, config = env$cfg)
  Y <- env$dz$Y; X <- env$X; O <- env$dz$O
  val <- function(v) -sum(zipln_vloglik(Y, X, O, zp$Pi, Omega, pv$optimizer$B(M = v$M, X = X), v$R, v$M, v$S2))
  monitor <- function(M, psi, lpi, zero) {
    S2 <- exp(psi); R <- plogis(exp(O + M + S2 / 2) + lpi) * zero
    conv$objective - val(list(M = M, S2 = S2, R = R))
  }
  rows <- lapply(c(1, 5, 20), function(k) {
    t <- system.time(v <- ve_cg(Y, X, O, Omega, zp$Pi, s$M, s$S2, cg_max = k, maxit = 5000, ftol = 1e-13,
                                monitor = monitor, max_time = 900))[["elapsed"]]
    row <- data.frame(case = pk[1], cg_max = k, newton = v$iterations, products = v$products, time = round(t, 1),
                      gain = round(conv$objective - val(v), 4))
    saveRDS(list(row = row, trace = v$trace), file.path(dir, paste0(tag, "_cg", k, ".rds")))
    say(pk[1], "cg_max", k, ": newton", v$iterations, "products", v$products, "time", round(t, 1), "gain", row$gain)
    row
  })
  do.call(rbind, rows)
}, mc.cores = 5, mc.preschedule = FALSE)
print(Filter(is.character, out)); res <- do.call(rbind, Filter(is.data.frame, out))
options(width = 200); print(res, row.names = FALSE)
saveRDS(res, "inst/simus_ZIPLN/convergence/ve_cg_frozen.rds")
