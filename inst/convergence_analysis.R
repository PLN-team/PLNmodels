## ============================================================
## Convergence analysis of the homemade Newton backend
## Datasets: trichoptera (n=49, p=17), barents (n=89, p=30),
##           mollusk (n=163, p=32),    oaks (n=116, p=114),
##           microcosm (n=880, p=259), scRNA (n=3918, p=500)
## Covariances: full, diagonal, spherical
## With / without covariates
## Note: full covariance for microcosm (~30-60s) and scRNA (very slow) included
## Output: inst/benchmark/
## ============================================================

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library(ggplot2)
  library(tidyr)
})

ctrl <- function(cov) PLN_param(backend = "homemade", covariance = cov, trace = 0)

## ---- trichoptera (n=49, p=17) ----
cat("Fitting trichoptera...\n")
tri <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
fits <- list(
  tri_full_nocov = PLN(Abundance ~ 1,                  data = tri, control = ctrl("full")),
  tri_diag_nocov = PLN(Abundance ~ 1,                  data = tri, control = ctrl("diagonal")),
  tri_sph_nocov  = PLN(Abundance ~ 1,                  data = tri, control = ctrl("spherical")),
  tri_full_cov   = PLN(Abundance ~ Wind + Temperature, data = tri, control = ctrl("full")),
  tri_diag_cov   = PLN(Abundance ~ Wind + Temperature, data = tri, control = ctrl("diagonal")),
  tri_sph_cov    = PLN(Abundance ~ Wind + Temperature, data = tri, control = ctrl("spherical"))
)

## ---- barents (n=89, p=30) ----
cat("Fitting barents...\n")
fits <- c(fits, list(
  bar_full_nocov = PLN(Abundance ~ 1                    + offset(log(Offset)), data = barents, control = ctrl("full")),
  bar_diag_nocov = PLN(Abundance ~ 1                    + offset(log(Offset)), data = barents, control = ctrl("diagonal")),
  bar_sph_nocov  = PLN(Abundance ~ 1                    + offset(log(Offset)), data = barents, control = ctrl("spherical")),
  bar_full_cov   = PLN(Abundance ~ Depth + Temperature  + offset(log(Offset)), data = barents, control = ctrl("full")),
  bar_diag_cov   = PLN(Abundance ~ Depth + Temperature  + offset(log(Offset)), data = barents, control = ctrl("diagonal")),
  bar_sph_cov    = PLN(Abundance ~ Depth + Temperature  + offset(log(Offset)), data = barents, control = ctrl("spherical"))
))

## ---- mollusk (n=163, p=32) ----
cat("Fitting mollusk...\n")
mol <- prepare_data(mollusk$Abundance, mollusk$Covariate)
fits <- c(fits, list(
  mol_full_nocov = PLN(Abundance ~ 1,              data = mol, control = ctrl("full")),
  mol_diag_nocov = PLN(Abundance ~ 1,              data = mol, control = ctrl("diagonal")),
  mol_sph_nocov  = PLN(Abundance ~ 1,              data = mol, control = ctrl("spherical")),
  mol_full_cov   = PLN(Abundance ~ site + season,  data = mol, control = ctrl("full")),
  mol_diag_cov   = PLN(Abundance ~ site + season,  data = mol, control = ctrl("diagonal")),
  mol_sph_cov    = PLN(Abundance ~ site + season,  data = mol, control = ctrl("spherical"))
))

## ---- oaks (n=116, p=114) ----
cat("Fitting oaks...\n")
fits <- c(fits, list(
  oak_full_nocov = PLN(Abundance ~ 1    + offset(log(Offset)), data = oaks, control = ctrl("full")),
  oak_diag_nocov = PLN(Abundance ~ 1    + offset(log(Offset)), data = oaks, control = ctrl("diagonal")),
  oak_sph_nocov  = PLN(Abundance ~ 1    + offset(log(Offset)), data = oaks, control = ctrl("spherical")),
  oak_full_cov   = PLN(Abundance ~ tree + offset(log(Offset)), data = oaks, control = ctrl("full")),
  oak_diag_cov   = PLN(Abundance ~ tree + offset(log(Offset)), data = oaks, control = ctrl("diagonal")),
  oak_sph_cov    = PLN(Abundance ~ tree + offset(log(Offset)), data = oaks, control = ctrl("spherical"))
))

## ---- microcosm (n=880, p=259) ----
cat("Fitting microcosm (diagonal + spherical)...\n")
fits <- c(fits, list(
  mic_diag_nocov = PLN(Abundance ~ 1    + offset(log(Offset)), data = microcosm, control = ctrl("diagonal")),
  mic_sph_nocov  = PLN(Abundance ~ 1    + offset(log(Offset)), data = microcosm, control = ctrl("spherical")),
  mic_diag_cov   = PLN(Abundance ~ site + offset(log(Offset)), data = microcosm, control = ctrl("diagonal")),
  mic_sph_cov    = PLN(Abundance ~ site + offset(log(Offset)), data = microcosm, control = ctrl("spherical"))
))
cat("Fitting microcosm full covariance (slow — O(n·p²) M-step with n=880, p=259)...\n")
fits <- c(fits, list(
  mic_full_nocov = PLN(Abundance ~ 1    + offset(log(Offset)), data = microcosm, control = ctrl("full")),
  mic_full_cov   = PLN(Abundance ~ site + offset(log(Offset)), data = microcosm, control = ctrl("full"))
))

## ---- scRNA (n=3918, p=500) — full covariance skipped ----
cat("Fitting scRNA diagonal + spherical (n=3918, p=500)...\n")
fits <- c(fits, list(
  scr_diag_nocov = PLN(counts ~ 1         + offset(log(total_counts)), data = scRNA, control = ctrl("diagonal")),
  scr_sph_nocov  = PLN(counts ~ 1         + offset(log(total_counts)), data = scRNA, control = ctrl("spherical")),
  scr_diag_cov   = PLN(counts ~ cell_line + offset(log(total_counts)), data = scRNA, control = ctrl("diagonal")),
  scr_sph_cov    = PLN(counts ~ cell_line + offset(log(total_counts)), data = scRNA, control = ctrl("spherical"))
))
cat("Fitting scRNA full covariance (very slow — O(n·p²) M-step with n=3918, p=500)...\n")
fits <- c(fits, list(
  scr_full_nocov = PLN(counts ~ 1         + offset(log(total_counts)), data = scRNA, control = ctrl("full")),
  scr_full_cov   = PLN(counts ~ cell_line + offset(log(total_counts)), data = scRNA, control = ctrl("full"))
))

cat("All fits done.\n")

## ---- Extract monitoring ----
mon <- lapply(names(fits), function(nm) {
  m   <- fits[[nm]]$optim_par
  obj <- m$objective
  obj_norm   <- (obj - min(obj)) / (max(obj) - min(obj) + .Machine$double.eps)
  rel_change <- abs(diff(obj)) / (abs(obj[-length(obj)]) + 1e-30)
  tail_n     <- max(5L, as.integer(0.2 * length(rel_change)))
  tail_slope <- if (all(tail(rel_change, tail_n) > 0))
                  mean(log10(tail(rel_change, tail_n) + 1e-30))
                else NA_real_
  parts <- strsplit(nm, "_")[[1]]
  list(
    name         = nm,
    dataset      = parts[1],
    covariance   = parts[2],
    covariates   = parts[3],
    n_iter       = m$iterations,
    status       = m$status,
    converged    = (m$status == 3),
    obj_init     = obj[1],
    obj_final    = obj[length(obj)],
    rel_drop     = (obj[1] - obj[length(obj)]) / abs(obj[1]),
    last_delta   = rel_change[length(rel_change)],
    tail_slope   = tail_slope,
    obj_seq      = obj,
    rel_seq      = rel_change,
    obj_norm_seq = obj_norm
  )
})

## ---- Summary table ----
cat("\n========== CONVERGENCE SUMMARY ==========\n")
sumtab <- do.call(rbind, lapply(mon, function(x) {
  data.frame(
    fit        = x$name,
    n_iter     = x$n_iter,
    converged  = x$converged,
    rel_drop   = signif(x$rel_drop, 3),
    last_delta = signif(x$last_delta, 3),
    tail_slope = signif(x$tail_slope, 3),
    stringsAsFactors = FALSE
  )
}))
print(sumtab, row.names = FALSE)

## ---- Plateau detection ----
cat("\n========== PLATEAU DETECTION (fraction of steps with delta < 1e-6) ==========\n")
for (x in mon) {
  r         <- x$rel_seq
  frac_flat <- mean(r < 1e-6)
  cat(sprintf("  %-25s  %5.1f%% flat steps | %d total\n",
              x$name, 100*frac_flat, x$n_iter))
}

## ---- EM kink detection ----
cat("\n========== EM KINK DETECTION (local minima in rel-change = EM M-step boundaries) ==========\n")
for (x in mon) {
  r     <- x$rel_seq
  kinks <- which(diff(sign(diff(log(r + 1e-30)))) == 2)
  cat(sprintf("  %-25s  ~%d EM kinks in %d inner steps  (%.1f steps/EM)\n",
              x$name, length(kinks), x$n_iter,
              if (length(kinks) > 0) x$n_iter / length(kinks) else NaN))
}

## ---- Convergence speed ----
cat("\n========== CONVERGENCE RATE (mean log10 rel-change in last 20% of steps) ==========\n")
for (x in mon) {
  cat(sprintf("  %-25s  tail log10(delta) = %.2f  (higher = slower)\n",
              x$name, x$tail_slope))
}

## ---- Build tidy data frames for plots ----
df_traj <- do.call(rbind, lapply(mon, function(x) {
  data.frame(name = x$name, dataset = x$dataset, covariance = x$covariance,
             covariates = x$covariates,
             step = seq_along(x$obj_norm_seq), obj_norm = x$obj_norm_seq,
             stringsAsFactors = FALSE)
}))

df_rel <- do.call(rbind, lapply(mon, function(x) {
  data.frame(name = x$name, dataset = x$dataset, covariance = x$covariance,
             covariates = x$covariates,
             step = seq_along(x$rel_seq), rel_change = pmax(x$rel_seq, 1e-16),
             stringsAsFactors = FALSE)
}))

## ---- Plot 1: normalised objective (log1p) ----
dataset_labels <- c(tri = "trichoptera (n=49, p=17)",
                    bar = "barents (n=89, p=30)",
                    mol = "mollusk (n=163, p=32)",
                    oak = "oaks (n=116, p=114)",
                    mic = "microcosm (n=880, p=259)",
                    scr = "scRNA (n=3918, p=500)")

p1 <- ggplot(df_traj, aes(step, obj_norm + 1e-6, colour = covariance, linetype = covariates)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ dataset, scales = "free_x",
             labeller = labeller(dataset = dataset_labels)) +
  scale_y_log10() +
  labs(title = "Normalised inner objective (log10 scale)",
       subtitle = "(obj - obj_min) / range  --  0 = fully converged",
       x = "Newton step (cumulated over EM)", y = "Normalised obj",
       colour = "Covariance", linetype = "Covariates") +
  theme_bw(base_size = 11)

ggsave("inst/benchmark/convergence_trajectory.pdf", p1, width = 15, height = 8)
cat("\nSaved: convergence_trajectory.pdf\n")

## ---- Plot 2: per-step relative change ----
p2 <- ggplot(df_rel, aes(step, rel_change, colour = covariance, linetype = covariates)) +
  geom_line(linewidth = 0.5, alpha = 0.8) +
  facet_wrap(~ dataset, scales = "free",
             labeller = labeller(dataset = dataset_labels)) +
  scale_y_log10() +
  geom_hline(yintercept = 1e-8, linetype = "dotted", colour = "grey50") +
  labs(title = "Per-step relative change |dobj|/|obj| (log10)",
       subtitle = "Dotted = ftol_in = 1e-8  |  Bumps = EM M-step boundary",
       x = "Newton step", y = "Relative change",
       colour = "Covariance", linetype = "Covariates") +
  theme_bw(base_size = 11)

ggsave("inst/benchmark/convergence_rel_change.pdf", p2, width = 15, height = 8)
cat("Saved: convergence_rel_change.pdf\n")

## ---- Plot 3: distribution of step sizes ----
p3 <- ggplot(df_rel, aes(rel_change, fill = covariance)) +
  geom_histogram(bins = 50, alpha = 0.65, position = "identity") +
  facet_grid(dataset ~ covariates, scales = "free_y",
             labeller = labeller(dataset = dataset_labels)) +
  scale_x_log10() +
  geom_vline(xintercept = 1e-8, linetype = "dotted") +
  labs(title = "Distribution of per-step relative changes",
       x = "Relative change (log10)", fill = "Covariance") +
  theme_bw(base_size = 10)

ggsave("inst/benchmark/convergence_step_dist.pdf", p3, width = 12, height = 14)
cat("Saved: convergence_step_dist.pdf\n")
