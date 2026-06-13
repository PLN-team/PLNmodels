#!/usr/bin/env Rscript
## Benchmark : PLNPCA init_method × backend
##
## Comparaisons testées :
##   init_method : "LM"     — nouvelle init (lm.fit sur log(Y), défaut actuel)
##                 "GLM"    — init Poisson GLM (p fits IRLS)
##                 "PLN-EM" — ancienne init master (PLNfit complet comme inception)
##   backend     : "nlopt"  — CCSAQ (défaut PLNPCA)
##                 "builtin" — L-BFGS joint + Wolfe fort
##
## IMPORTANT : ne jamais lancer en parallèle (BLAS multithreadé).

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
})
cat("PLNmodels", as.character(packageVersion("PLNmodels")), "\n\n")

data(trichoptera); data(barents); data(oaks)
tri <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

## --------------------------------------------------------------------------
## Configurations testées
## --------------------------------------------------------------------------
RANKS <- list(
  trichoptera = c(1L, 3L, 5L),
  barents     = c(3L, 5L, 10L),
  oaks        = c(5L, 10L, 20L)
)

DATASETS <- list(
  list(name = "trichoptera",
       formula = Abundance ~ 1,
       data    = tri),
  list(name = "barents",
       formula = Abundance ~ Depth + Temperature,
       data    = barents),
  list(name = "oaks",
       formula = Abundance ~ 1 + offset(log(Offset)),
       data    = oaks)
)

BACKENDS <- c("nlopt", "builtin")

## --------------------------------------------------------------------------
## Helpers
## --------------------------------------------------------------------------
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
  elapsed <- round((proc.time() - t0)[["elapsed"]], 2)
  list(fit = fit, elapsed = elapsed)
}

## --------------------------------------------------------------------------
## Main loop  (tout séquentiel)
## --------------------------------------------------------------------------
results <- list()

for (ds in DATASETS) {
  ranks <- RANKS[[ds$name]]
  cat(sprintf("=== %s (ranks: %s) ===\n", ds$name, paste(ranks, collapse = ",")))

  ## Pré-calculer le PLNfit pour l'init PLN-EM (une fois par dataset)
  cat("  [PLN-EM] fitting full PLN inception ...\n")
  t_pln <- proc.time()
  pln_inc <- PLN(ds$formula, data = ds$data, control = PLN_param(trace = 0))
  t_pln <- round((proc.time() - t_pln)[["elapsed"]], 2)
  cat(sprintf("  [PLN-EM] PLN done in %.1fs\n", t_pln))

  for (backend in BACKENDS) {
    for (init_method in c("LM", "GLM", "PLN-EM")) {
      label <- sprintf("%-8s init=%-7s", backend, init_method)
      cat(sprintf("  %s ...", label))

      res <- run_one(ds, ranks, backend, init_method, pln_inception = pln_inc)
      elapsed <- res$elapsed

      if (is.null(res$fit)) {
        cat(" FAILED\n")
        next
      }

      for (m in res$fit$models) {
        q      <- m$rank
        loglik <- round(m$loglik, 2)
        n_iter <- if (!is.null(m$optim_par$iterations)) m$optim_par$iterations else NA_integer_
        results[[length(results) + 1]] <- data.frame(
          dataset     = ds$name,
          rank        = q,
          backend     = backend,
          init_method = init_method,
          loglik      = loglik,
          n_iter      = as.integer(n_iter),
          time_total_s = elapsed,
          stringsAsFactors = FALSE
        )
      }
      ll_str <- paste(sapply(res$fit$models, function(m) round(m$loglik, 1)), collapse = " | ")
      cat(sprintf(" %.1fs  ll: %s\n", elapsed, ll_str))
    }
  }
  cat("\n")
}

## --------------------------------------------------------------------------
## Table de résultats
## --------------------------------------------------------------------------
df <- do.call(rbind, results)

cat("\n=== Résultats (loglik, plus grand = meilleur) ===\n\n")

for (ds_name in unique(df$dataset)) {
  sub <- df[df$dataset == ds_name, ]
  cat(sprintf("-- %s --\n", ds_name))
  cat(sprintf("  %-8s %-8s  %s\n", "backend", "init", paste(sprintf("q=%-4s", unique(sub$rank)), collapse = "  ")))
  for (b in BACKENDS) {
    for (im in c("LM", "GLM", "PLN-EM")) {
      row <- sub[sub$backend == b & sub$init_method == im, ]
      if (nrow(row) == 0) next
      row <- row[order(row$rank), ]
      vals <- sprintf("%-6.1f", row$loglik)
      cat(sprintf("  %-8s %-8s  %s\n", b, im, paste(vals, collapse = "    ")))
    }
  }
  cat("\n")
}

out_rds <- "inst/benchmark/bench_plnpca_init_results.rds"
saveRDS(df, out_rds)
cat("Résultats sauvegardés dans", out_rds, "\n")
