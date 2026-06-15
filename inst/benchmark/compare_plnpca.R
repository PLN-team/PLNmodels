#!/usr/bin/env Rscript
## Benchmark PLNPCA : rangs × backends × méthodes d'initialisation
##
## Dimensions testées :
##   rangs       : propres à chaque dataset (voir RANKS ci-dessous)
##   backend     : auto-détectés (builtin, nlopt, torch si disponible)
##   init_method : "LM"     — lm.fit sur log(Y) (défaut)
##                 "GLM"    — p fits Poisson IRLS
##                 "PLN-EM" — PLNfit complet comme inception (ancienne init master)
##
## IMPORTANT : séquentiel uniquement (BLAS multithreadé)
## Sortie : compare_plnpca_init_results.rds

suppressPackageStartupMessages(devtools::load_all(quiet = TRUE))
cat("PLNmodels", as.character(packageVersion("PLNmodels")), "\n\n")

# ── Backends auto-détectés ────────────────────────────────────────────────────
backends <- tryCatch(
  intersect(as.character(formals(PLNPCA_param)$backend)[-1], c("builtin", "nlopt", "torch")),
  error = function(e) "nlopt"
)
torch_ok <- requireNamespace("torch", quietly = TRUE)
backends <- backends[backends != "torch" | torch_ok]
cat("PLNPCA backends:", paste(backends, collapse = ", "), "\n\n")

# ── Données ──────────────────────────────────────────────────────────────────
data(trichoptera); data(barents); data(oaks)
tri <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

RANKS <- list(
  trichoptera = c(1L, 3L, 5L),
  barents     = c(3L, 5L, 10L),
  oaks        = c(5L, 10L, 20L)
)

DATASETS <- list(
  list(name = "trichoptera", formula = Abundance ~ 1,                       data = tri),
  list(name = "barents",     formula = Abundance ~ Depth + Temperature,     data = barents),
  list(name = "oaks",        formula = Abundance ~ 1 + offset(log(Offset)), data = oaks)
)

INIT_METHODS <- c("LM", "GLM", "PLN-EM")

# ── Helper ───────────────────────────────────────────────────────────────────
run_one <- function(ds, ranks, backend, init_method, pln_inception = NULL) {
  ctrl <- if (identical(init_method, "PLN-EM")) {
    PLNPCA_param(backend = backend, trace = 0, inception = pln_inception)
  } else {
    PLNPCA_param(backend = backend, trace = 0, init_method = init_method)
  }
  t0  <- proc.time()
  fit <- tryCatch(
    PLNPCA(ds$formula, data = ds$data, ranks = ranks, control = ctrl),
    error = function(e) { message("  ERROR: ", conditionMessage(e)); NULL }
  )
  list(fit = fit, elapsed = round((proc.time() - t0)[["elapsed"]], 2))
}

# ── Benchmark ─────────────────────────────────────────────────────────────────
results <- list()

for (ds in DATASETS) {
  ranks <- RANKS[[ds$name]]
  cat(sprintf("=== %s (ranks: %s) ===\n", ds$name, paste(ranks, collapse = ",")))

  cat("  [PLN-EM] fitting PLN inception ...\n")
  t0_pln <- proc.time()
  pln_inc <- PLN(ds$formula, data = ds$data, control = PLN_param(trace = 0))
  cat(sprintf("  [PLN-EM] done in %.1fs\n", (proc.time() - t0_pln)[["elapsed"]]))

  for (backend in backends) {
    for (init_method in INIT_METHODS) {
      cat(sprintf("  %-8s init=%-7s ...", backend, init_method))
      res <- run_one(ds, ranks, backend, init_method, pln_inception = pln_inc)
      if (is.null(res$fit)) { cat(" FAILED\n"); next }

      for (m in res$fit$models) {
        results[[length(results) + 1]] <- data.frame(
          dataset      = ds$name,
          rank         = m$rank,
          backend      = backend,
          init_method  = init_method,
          loglik       = round(m$loglik, 2),
          n_iter       = if (!is.null(m$optim_par$iterations)) as.integer(m$optim_par$iterations) else NA_integer_,
          time_total_s = res$elapsed,
          status       = if (!is.null(m$optim_par$status)) as.character(m$optim_par$status) else NA_character_,
          stringsAsFactors = FALSE
        )
      }
      ll_str <- paste(sapply(res$fit$models, function(m) round(m$loglik, 1)), collapse = " | ")
      cat(sprintf(" %.1fs  ll: %s\n", res$elapsed, ll_str))
    }
  }
  cat("\n")
}

# ── Tableau récapitulatif ─────────────────────────────────────────────────────
df <- do.call(rbind, results)

cat("\n=== Résultats (loglik, plus grand = meilleur) ===\n\n")
for (ds_name in unique(df$dataset)) {
  sub  <- df[df$dataset == ds_name, ]
  rnks <- sort(unique(sub$rank))
  cat(sprintf("-- %s --\n", ds_name))
  cat(sprintf("  %-8s %-8s  %s\n", "backend", "init",
              paste(sprintf("q=%-4s", rnks), collapse = "  ")))
  for (b in backends) {
    for (im in INIT_METHODS) {
      row <- sub[sub$backend == b & sub$init_method == im, ]
      if (nrow(row) == 0) next
      row  <- row[order(row$rank), ]
      vals <- sprintf("%-6.1f", row$loglik)
      cat(sprintf("  %-8s %-8s  %s\n", b, im, paste(vals, collapse = "    ")))
    }
  }
  cat("\n")
}

# ── Sauvegarde ────────────────────────────────────────────────────────────────
out_rds <- "inst/benchmark/compare_plnpca_init_results.rds"
saveRDS(df, out_rds)
cat("Résultats sauvegardés dans", out_rds, "\n")
