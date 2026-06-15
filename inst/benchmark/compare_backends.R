## ============================================================
## Backend comparison: tous modèles × tous backends
## - PLN (full/diag/sph, nocov + cov) sur 6 datasets
## - ZIPLN + PLNPCA (cov) sur 3 datasets (trichoptera, barents, oaks)
## Backends : auto-détectés (builtin, nlopt, torch si disponible)
## IMPORTANT : séquentiel uniquement (BLAS multithreadé)
## Sortie : compare_backends_results.{rds,csv} + 3 PDF PLN builtin/nlopt
## ============================================================

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# ── Backends auto-détectés ────────────────────────────────────────────────────
get_be <- function(param_fn, exclude = character(0)) {
  be <- tryCatch(as.character(formals(param_fn)$backend)[-1], error = function(e) "nlopt")
  setdiff(intersect(be, c("builtin", "nlopt", "torch")), exclude)
}
torch_ok       <- requireNamespace("torch", quietly = TRUE)
backends_pln   <- get_be(PLN_param,    if (!torch_ok) "torch" else character(0))
backends_zipln <- get_be(ZIPLN_param,  "torch")
backends_pca   <- get_be(PLNPCA_param, if (!torch_ok) "torch" else character(0))
cat("PLN backends    :", paste(backends_pln,   collapse = ", "), "\n")
cat("ZIPLN backends  :", paste(backends_zipln, collapse = ", "), "\n")
cat("PLNPCA backends :", paste(backends_pca,   collapse = ", "), "\n\n")

# ── Données ──────────────────────────────────────────────────────────────────
data(trichoptera); data(barents); data(oaks)
data(mollusk); data(microcosm); data(scRNA)

tri <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
mol <- prepare_data(mollusk$Abundance,     mollusk$Covariate)
oak <- prepare_data(oaks$Abundance,
                    oaks[, c("tree", "distTOtrunk", "orientation", "pmInfection")],
                    offset = oaks$Offset)

# ── Spécification des datasets ────────────────────────────────────────────────
# all_models = TRUE : PLN + ZIPLN + PLNPCA ; FALSE : PLN uniquement (grands datasets)
datasets <- list(
  list(name = "trichoptera", data = tri,       all_models = TRUE,
       pca_rank = round(sqrt(ncol(trichoptera$Abundance))),
       f_nocov = Abundance ~ 1 + offset(log(Offset)),
       f_cov   = Abundance ~ Wind + Temperature + offset(log(Offset))),
  list(name = "barents",     data = barents,   all_models = TRUE,
       pca_rank = round(sqrt(ncol(barents$Abundance))),
       f_nocov = Abundance ~ 1 + offset(log(Offset)),
       f_cov   = Abundance ~ Depth + Temperature + offset(log(Offset))),
  list(name = "oaks",        data = oak,       all_models = TRUE,
       pca_rank = round(sqrt(ncol(oaks$Abundance))),
       f_nocov = Abundance ~ 1 + offset(log(Offset)),
       f_cov   = Abundance ~ tree + offset(log(Offset))),
  list(name = "mollusk",     data = mol,       all_models = FALSE, pca_rank = NULL,
       f_nocov = Abundance ~ 1,
       f_cov   = Abundance ~ site + season),
  list(name = "microcosm",   data = microcosm, all_models = FALSE, pca_rank = NULL,
       f_nocov = Abundance ~ 1    + offset(log(Offset)),
       f_cov   = Abundance ~ site + offset(log(Offset))),
  list(name = "scRNA",       data = scRNA,     all_models = FALSE, pca_rank = NULL,
       f_nocov = counts ~ 1         + offset(log(total_counts)),
       f_cov   = counts ~ cell_line + offset(log(total_counts)))
)

cov_types <- c("full", "diagonal", "spherical")

# ── Helper ───────────────────────────────────────────────────────────────────
run_one <- function(model, dataset, cov_label, backend, covariance, expr) {
  cat(sprintf("  %-16s %-12s %-8s [%-10s %s]...",
              model, dataset, backend,
              if (is.na(covariance)) "—" else covariance,
              cov_label))
  t0 <- proc.time()
  m  <- tryCatch(expr, error = function(e) { message("  ERROR: ", e$message); NULL })
  elapsed <- (proc.time() - t0)[["elapsed"]]
  if (inherits(m, "PLNPCAfamily")) m <- m$models[[1]]
  ll    <- if (is.null(m)) NA_real_ else round(m$loglik, 3)
  x_iter <- if (is.null(m)) NULL else m$optim_par$iterations
  niter  <- if (is.null(x_iter) || length(x_iter) == 0) NA_integer_ else as.integer(x_iter)
  x_stat <- if (is.null(m)) NULL else m$optim_par$status
  stat   <- if (is.null(x_stat) || length(x_stat) == 0) NA_character_ else as.character(x_stat)
  cat(sprintf("  loglik=%.1f  iter=%s  t=%.2fs\n",
              if (is.na(ll)) NaN else ll, niter, elapsed))
  data.frame(model = model, dataset = dataset, covariates = cov_label,
             backend = backend, covariance = covariance,
             time_s = round(elapsed, 3), loglik = ll,
             n_iter = niter, status = stat,
             stringsAsFactors = FALSE)
}

# ── Benchmark ─────────────────────────────────────────────────────────────────
results <- list()

for (ds in datasets) {
  cat(sprintf("\n══ %s %s══\n", ds$name,
              if (ds$all_models) "(PLN + ZIPLN + PLNPCA) " else "(PLN only) "))

  # PLN : toutes les covariances, nocov + cov
  for (cov in cov_types) {
    if (ds$name == "scRNA" && cov == "full") next  # O(n*p^2) trop lent
    for (bk in backends_pln) {
      results[[length(results)+1]] <- run_one(
        "PLN", ds$name, "nocov", bk, cov,
        PLN(ds$f_nocov, data = ds$data,
            control = PLN_param(backend = bk, covariance = cov, trace = 0)))
      results[[length(results)+1]] <- run_one(
        "PLN", ds$name, "cov", bk, cov,
        PLN(ds$f_cov, data = ds$data,
            control = PLN_param(backend = bk, covariance = cov, trace = 0)))
    }
  }

  if (!ds$all_models) next

  # ZIPLN : formule avec covariables uniquement
  for (bk in backends_zipln) {
    results[[length(results)+1]] <- run_one(
      "ZIPLN", ds$name, "cov", bk, "full",
      ZIPLN(ds$f_cov, data = ds$data,
            control = ZIPLN_param(backend = bk, trace = 0)))
  }

  # PLNPCA : formule avec covariables uniquement
  for (bk in backends_pca) {
    results[[length(results)+1]] <- run_one(
      sprintf("PLNPCA(q=%d)", ds$pca_rank), ds$name, "cov", bk, NA_character_,
      PLNPCA(ds$f_cov, data = ds$data, ranks = ds$pca_rank,
             control = PLNPCA_param(backend = bk, trace = 0)))
  }
}

# ── Tableau récapitulatif ─────────────────────────────────────────────────────
cat("\n\n══════════════════════════════════════════════\n")
cat(" RÉCAPITULATIF\n")
cat("══════════════════════════════════════════════\n\n")

df <- do.call(rbind, results)
rownames(df) <- NULL

for (dname in unique(df$dataset)) {
  cat(sprintf("\n--- %s ---\n", dname))
  sub <- df[df$dataset == dname,
            c("model", "covariance", "covariates", "backend", "loglik", "n_iter", "time_s")]
  sub$covariance[is.na(sub$covariance)] <- "—"
  print(sub, row.names = FALSE, width = 120)
}

# ── Sauvegarde ────────────────────────────────────────────────────────────────
out_rds <- "inst/benchmark/compare_backends_results.rds"
out_csv <- "inst/benchmark/compare_backends_results.csv"
saveRDS(df, out_rds)
write.csv(df, out_csv, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", out_rds))
cat(sprintf("Saved: %s\n", out_csv))

# ── Plots PLN : builtin vs nlopt ──────────────────────────────────────────────
if (!all(c("builtin", "nlopt") %in% backends_pln)) {
  cat("\nPlots skipped (need both builtin and nlopt backends)\n")
  quit(save = "no")
}

df_pln <- df[df$model == "PLN" & df$backend %in% c("builtin", "nlopt"), ]
df_pln$fit <- paste(df_pln$dataset, df_pln$covariance, df_pln$covariates, sep = " / ")

p1 <- ggplot(df_pln, aes(x = fit, y = time_s, fill = backend)) +
  geom_col(position = "dodge", width = 0.7) +
  facet_wrap(~ dataset, scales = "free", nrow = 2) +
  scale_fill_manual(values = c(builtin = "#E69F00", nlopt = "#56B4E9")) +
  labs(title = "Computation time: builtin vs nlopt (PLN)",
       x = NULL, y = "Elapsed time (s)", fill = "Backend") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

ggsave("inst/benchmark/backend_time.pdf", p1, width = 18, height = 10)
cat("\nSaved: backend_time.pdf\n")

df_wide_pln <- df_pln %>%
  select(dataset, covariates, covariance, backend, loglik, time_s) %>%
  pivot_wider(names_from = backend, values_from = c(loglik, time_s)) %>%
  mutate(ll_diff = loglik_builtin - loglik_nlopt,
         speedup = time_s_nlopt / time_s_builtin,
         fit     = paste(dataset, covariance, covariates, sep = " / "))

p2 <- ggplot(df_wide_pln, aes(x = fit, y = ll_diff,
                               fill = ifelse(ll_diff > 0, "builtin better", "nlopt better"))) +
  geom_col(width = 0.7) +
  facet_wrap(~ dataset, scales = "free", nrow = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("builtin better" = "#009E73", "nlopt better" = "#D55E00"),
                    name = NULL) +
  labs(title = "loglik difference: builtin minus nlopt (PLN)",
       subtitle = "Positive = builtin finds better solution",
       x = NULL, y = "loglik(builtin) - loglik(nlopt)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

ggsave("inst/benchmark/backend_loglik.pdf", p2, width = 18, height = 10)
cat("Saved: backend_loglik.pdf\n")

p3 <- ggplot(df_wide_pln, aes(x = fit, y = speedup, fill = dataset)) +
  geom_col(width = 0.7) +
  facet_wrap(~ dataset, scales = "free", nrow = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  labs(title = "Speedup: nlopt_time / builtin_time (PLN)",
       subtitle = "> 1 means builtin is faster",
       x = NULL, y = "Speedup ratio") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "none")

ggsave("inst/benchmark/backend_speedup.pdf", p3, width = 18, height = 10)
cat("Saved: backend_speedup.pdf\n")
