## Species-wise switch of regime after convergence. For each species, compare the share of
## its zeros the model attributes to the zero inflation (mean of R on its zeros) with the
## share of its zeros in excess under a Poisson-lognormal fitted to its positive counts.
## Species where the gap d_j is large are switched to the "positive" starting point
## (zeros by inflation), those where it is very negative to the "lm" one (zeros by the
## latent variance); the model is refitted from this composite state.
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("inst/simus_ZIPLN/species_regime/cases.R")

## per-species gap between the excess zeros and the zeros attributed to the inflation
species_gap <- function(fit, dz) {
  sp <- compute_ZIPLN_starting_point(dz$Y, dz$X, dz$X0, dz$O, method = "positive")
  zero <- dz$Y == 0
  zr <- colMeans(zero)
  excess_share <- colSums(sp$R) / pmax(colSums(zero), 1)        # share of the zeros in excess (positive start)
  attributed <- colSums(fit$var_par$R) / pmax(colSums(zero), 1) # share attributed to the inflation by the fit
  list(gap = excess_share - attributed, positive = sp,
       lm = compute_ZIPLN_starting_point(dz$Y, dz$X, dz$X0, dz$O, method = "lm"))
}
## state of `fit` with the columns `J` taken from the starting point `start` (on the zeros
## for the variational parameters)
switch_species <- function(fit, J, start, Y) {
  g <- fit$clone()
  B <- fit$model_par$B; M <- fit$var_par$M; S2 <- fit$var_par$S2; R <- fit$var_par$R
  zero <- Y == 0
  for (j in J) {
    z <- zero[, j]
    B[, j] <- start$B[, j]; M[z, j] <- start$M[z, j]; S2[z, j] <- start$S2[z, j]; R[, j] <- start$R[, j]
  }
  g$update(B = B, M = M, S2 = S2, R = R)
  g
}
refit <- function(cs, inception) suppressWarnings(suppressMessages(
  ZIPLN(cs$formula, data = cs$data, zi = cs$zi, control = ZIPLN_param(trace = 0, inception = inception))))

run <- function(k) {
  RhpcBLASctl::blas_set_num_threads(1)
  cs <- cases[[k]]
  dz <- extract_model_zi(call("ZIPLN", formula = cs$formula, data = quote(dd)), list2env(list(dd = cs$data)))
  t0 <- system.time(f <- suppressWarnings(suppressMessages(ZIPLN(cs$formula, data = cs$data, zi = cs$zi, control = ZIPLN_param(trace = 0)))))[["elapsed"]]
  g <- species_gap(f, dz)
  moves <- list(
    "to ZI, gap > 0.5"      = list(J = which(g$gap > 0.5), start = "positive"),
    "to ZI, gap > 0.3"      = list(J = which(g$gap > 0.3), start = "positive"),
    "to ZI, gap > 0.1"      = list(J = which(g$gap > 0.1), start = "positive"),
    "to variance, gap < -0.3" = list(J = which(g$gap < -0.3), start = "lm"),
    "both, |gap| > 0.3"     = list(J = which(abs(g$gap) > 0.3), start = "both"))
  rows <- lapply(names(moves), function(m) {
    J <- moves[[m]]$J
    if (!length(J)) return(data.frame(move = m, n_switched = 0, gain = 0, time = 0))
    state <- if (moves[[m]]$start == "both") {
      s <- switch_species(f, J[g$gap[J] > 0], g$positive, dz$Y); switch_species(s, J[g$gap[J] < 0], g$lm, dz$Y)
    } else switch_species(f, J, g[[moves[[m]]$start]], dz$Y)
    t <- system.time(h <- refit(cs, state))[["elapsed"]]
    data.frame(move = m, n_switched = length(J), gain = h$loglik - f$loglik, time = t)
  })
  data.frame(case = cs$name, kept = f$optim_par$init, p = ncol(dz$Y), loglik = f$loglik, time_fit = t0,
             gap_median = median(g$gap), do.call(rbind, rows))
}
out <- parallel::mclapply(seq_along(cases), function(k) tryCatch(run(k), error = function(e) paste(cases[[k]]$name, conditionMessage(e))),
                          mc.cores = 13, mc.preschedule = FALSE)
err <- Filter(is.character, out); if (length(err)) print(unlist(err))
res <- do.call(rbind, Filter(is.data.frame, out))
saveRDS(res, "inst/simus_ZIPLN/species_regime/flips.rds")
options(width = 200)
print(reshape(res[, c("case", "kept", "move", "n_switched", "gain")], idvar = c("case", "kept"), timevar = "move", direction = "wide"), digits = 4, row.names = FALSE)
