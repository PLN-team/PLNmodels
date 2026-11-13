## Choose the starting point species by species, on cheap diagonal fits, then run a single
## full fit from the composite start. Usage: Rscript species_start.R sim|real
what <- commandArgs(TRUE)[1]
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("inst/simus_ZIPLN/init/simulate.R"); source("inst/simus_ZIPLN/init/inits.R")
n_cores <- as.integer(Sys.getenv("N_CORES", "12"))

## Variational bound of a diagonal ZIPLN fit, species by species (it is separable)
species_elbo <- function(fit, Y, X, O) {
  M <- fit$var_par$M; S2 <- fit$var_par$S2; R <- fit$var_par$R
  B <- fit$model_par$B; Pi <- fit$model_par$Pi
  omega <- diag(as.matrix(fit$model_par$Omega)); n <- nrow(Y)
  mu0 <- qlogis(pmin(pmax(Pi, 1e-12), 1 - 1e-12)); A <- exp(O + M + S2 / 2)
  xlx <- function(x) ifelse(x > 0, x * log(x), 0)
  colSums((1 - R) * (Y * (O + M) - A - PLNmodels:::.logfactorial(Y)) + R * mu0 - log1p(exp(mu0)) +
            .5 * log(S2) - .5 * rep(omega, each = n) * ((M - X %*% B)^2 + S2) - xlx(R) - xlx(1 - R)) +
    .5 * n * log(omega) + .5 * n
}

run_case <- function(formula, data, zi) {
  dz <- PLNmodels:::extract_model_zi(call("ZIPLN", formula = formula, data = quote(data)), environment())
  Y <- dz$Y; X <- dz$X
  fit_one <- function(init, covariance, start = NULL) {
    ctrl <- ZIPLN_param(trace = 0, init = init, covariance = covariance)
    t <- system.time(f <- if (is.null(start)) suppressWarnings(suppressMessages(ZIPLN(formula, data = data, zi = zi, control = ctrl)))
                          else zipln_with_init(function(...) start, formula, data = data, zi = zi, control = ctrl))[["elapsed"]]
    list(fit = f, time = t)
  }
  d_lm <- fit_one("lm", "diagonal"); d_pos <- fit_one("positive", "diagonal")
  O <- d_lm$fit$latent - d_lm$fit$var_par$M
  e_lm <- species_elbo(d_lm$fit, Y, X, O); e_pos <- species_elbo(d_pos$fit, Y, X, O)
  check <- abs(sum(e_lm) - d_lm$fit$loglik)
  choose_pos <- e_pos > e_lm
  pick <- function(a, b) { out <- a; out[, choose_pos] <- b[, choose_pos]; out }
  start <- list(B = pick(d_lm$fit$model_par$B, d_pos$fit$model_par$B),
                M = pick(d_lm$fit$var_par$M, d_pos$fit$var_par$M),
                S2 = pick(d_lm$fit$var_par$S2, d_pos$fit$var_par$S2),
                R = pick(d_lm$fit$var_par$R, d_pos$fit$var_par$R),
                B0 = if (length(d_lm$fit$model_par$B0)) pick(d_lm$fit$model_par$B0, d_pos$fit$model_par$B0) else matrix(0, 0, ncol(Y)))
  mixed <- fit_one("lm", "full", start)
  ## decision on the total diagonal bound, then a full fit from that raw starting point
  total <- if (sum(e_pos) > sum(e_lm)) "positive" else "lm"
  data.frame(share_positive = mean(choose_pos), diag_check = check, diag_total_choice = total,
             time_diag = d_lm$time + d_pos$time, loglik_mixed = mixed$fit$loglik, time_mixed = mixed$time,
             it_mixed = mixed$fit$optim_par$iterations)
}

if (what == "sim") {
  tasks <- expand.grid(rep = 1:10, cfg = seq_len(nrow(configs)))
  traj <- readRDS("inst/simus_ZIPLN/init/trajectories_sim.rds")
  one <- function(k) {
    RhpcBLASctl::blas_set_num_threads(1)
    cfg <- configs[tasks$cfg[k], ]; sim <- simulate(cfg, (tasks$cfg[k] - 1) * reps + tasks$rep[k])
    dat <- data.frame(row.names = seq_len(nrow(sim$Y))); dat$Y <- sim$Y
    tr <- traj[[k]]
    data.frame(case = tr$case, group = tr$group, run_case(Y ~ 1, dat, "col"),
               loglik_lm = tr$lm$loglik, loglik_pos = tr$positive$loglik, time_lm = tr$lm$time, time_pos = tr$positive$time)
  }
  todo <- seq_len(nrow(tasks))
} else {
  data(trichoptera); data(barents); data(oaks); data(mollusk); data(microcosm)
  tri <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
  mol <- suppressMessages(prepare_data(mollusk$Abundance, mollusk$Covariate))
  cases <- list(
    list("trichoptera single", tri, Abundance ~ 1 + offset(log(Offset)), "single"),
    list("trichoptera col", tri, Abundance ~ 1 + offset(log(Offset)), "col"),
    list("barents single", barents, Abundance ~ 1 + offset(log(Offset)), "single"),
    list("barents col", barents, Abundance ~ 1 + offset(log(Offset)), "col"),
    list("mollusk single", mol, Abundance ~ 1 + offset(log(Offset)), "single"),
    list("mollusk col", mol, Abundance ~ 1 + offset(log(Offset)), "col"),
    list("mollusk site|site", mol, Abundance ~ 0 + site + offset(log(Offset)) | 0 + site, "covar"),
    list("oaks single", oaks, Abundance ~ 1 + offset(log(Offset)), "single"),
    list("oaks col", oaks, Abundance ~ 1 + offset(log(Offset)), "col"),
    list("oaks tree col", oaks, Abundance ~ 0 + tree + offset(log(Offset)), "col"),
    list("oaks tree|tree", oaks, Abundance ~ 0 + tree + offset(log(Offset)) | 0 + tree, "covar"),
    list("microcosm site col", microcosm, Abundance ~ 0 + site + offset(log(Offset)), "col"),
    list("microcosm site|site", microcosm, Abundance ~ 0 + site + offset(log(Offset)) | 0 + site, "covar"))
  traj <- readRDS("inst/simus_ZIPLN/init/trajectories_real.rds")
  one <- function(k) {
    RhpcBLASctl::blas_set_num_threads(1)
    cs <- cases[[k]]
    tr <- Filter(function(x) x$case == cs[[1]], traj); tr <- setNames(lapply(tr, `[[`, "fit"), vapply(tr, `[[`, "", "init"))
    data.frame(case = cs[[1]], group = "real", run_case(cs[[3]], cs[[2]], cs[[4]]),
               loglik_lm = tr$lm$loglik, loglik_pos = tr$positive$loglik, time_lm = tr$lm$time, time_pos = tr$positive$time)
  }
  todo <- c(12, 13, 1:11)
}
out <- parallel::mclapply(todo, function(k) tryCatch(one(k), error = function(e) paste(k, conditionMessage(e))),
                          mc.cores = n_cores, mc.preschedule = FALSE)
err <- Filter(is.character, out); if (length(err)) print(unlist(err))
res <- do.call(rbind, Filter(is.data.frame, out))
saveRDS(res, sprintf("inst/simus_ZIPLN/init/species_start_%s.rds", what))
cat("done", what, "\n")
