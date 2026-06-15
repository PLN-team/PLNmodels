#!/usr/bin/env Rscript
## Usage: Rscript compare_branches.R <rds_master> <rds_ce> [out_csv]
##
## Lit deux fichiers RDS produits par run_branch.R et affiche :
##   1. Comparaison globale  : meilleur backend CE  vs  nlopt master (référence)
##   2. Détail CE            : tous les backends CE vs  nlopt CE
##   3. Norme des paramètres : B et Omega (Frobenius)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: compare_branches.R <rds_master> <rds_ce> [out_csv]")

df_m  <- readRDS(args[1])
df_ce <- readRDS(args[2])
out_csv <- if (length(args) >= 3) args[3] else
  file.path(dirname(args[1]), "branch_comparison.csv")

branch_m  <- unique(df_m$branch)
branch_ce <- unique(df_ce$branch)

HR <- strrep("=", 110)
hr <- strrep("-", 110)

## ── Helpers ───────────────────────────────────────────────────────────────────
fmt_ll <- function(x) sprintf("%10.1f", x)
fmt_t  <- function(x) sprintf("%6.2fs", x)
na_str <- function(x, fmt = "%s") if (is.na(x)) "    ?" else sprintf(fmt, x)

best_ce_by_model_dataset <- function(df) {
  do.call(rbind, lapply(
    split(df, list(df$model, df$dataset), drop = TRUE),
    function(x) x[which.max(x$loglik), , drop = FALSE]
  ))
}

## ── 1. Comparaison globale : meilleur CE vs nlopt master ─────────────────────
cat("\n", HR, "\n", sep = "")
cat(sprintf("  BENCHMARK  %s  vs  %s\n", branch_ce, branch_m))
cat(HR, "\n\n", sep = "")

master_nlopt <- df_m[df_m$backend == "nlopt", ]
best_ce      <- best_ce_by_model_dataset(df_ce)

## renommer best_ce pour éviter les conflits de suffixe
best_ce_sel <- best_ce[, c("model","dataset","loglik","time_s","n_iter","backend","status")]
names(best_ce_sel)[names(best_ce_sel) %in% c("loglik","time_s","n_iter")] <-
  paste0(names(best_ce_sel)[names(best_ce_sel) %in% c("loglik","time_s","n_iter")], ".ce")

comp <- merge(
  master_nlopt[, c("model","dataset","loglik","time_s","n_iter")],
  best_ce_sel,
  by = c("model","dataset")
)
names(comp)[names(comp) == "loglik"]  <- "loglik.m"
names(comp)[names(comp) == "time_s"]  <- "time_s.m"
names(comp)[names(comp) == "n_iter"]  <- "n_iter.m"
comp$delta   <- round(comp$loglik.ce - comp$loglik.m, 1)
comp$speedup <- round(comp$time_s.m  / comp$time_s.ce, 2)
comp <- comp[order(comp$model, comp$dataset), ]

cat(sprintf("  Référence master : nlopt    |  CE : meilleur backend disponible\n"))
cat(sprintf("  delta = ll_CE_best - ll_master_nlopt   (▲ CE meilleur, ▼ CE pire, seuil |Δ|>5)\n\n"))
cat(sprintf("  %-12s %-12s %-10s | %10s %10s %+7s | %8s %8s\n",
    "Modèle","Dataset","CE_best","ll_master","ll_CE","delta","speedup","iter_CE"))
cat(sprintf("  %s\n", strrep("-", 90)))

prev_model <- ""
for (i in seq_len(nrow(comp))) {
  r <- comp[i, ]
  if (r$model != prev_model && i > 1) cat("\n")
  prev_model <- r$model
  flag <- if (!is.na(r$delta) && abs(r$delta) > 5)
    if (r$delta > 0) " ▲" else " ▼" else "  "
  cat(sprintf("  %-12s %-12s %-10s | %10.1f %10.1f %+7.1f%s | %8.2fx %8s\n",
    r$model, r$dataset, r$backend,
    r$loglik.m, r$loglik.ce, r$delta, flag,
    r$speedup, na_str(r$n_iter.ce, "%d")))
}

## ── 2. Détail CE : tous backends vs nlopt CE ──────────────────────────────────
ce_nlopt  <- df_ce[df_ce$backend == "nlopt", ]
ce_others <- df_ce[df_ce$backend != "nlopt", ]

if (nrow(ce_others) > 0) {
  cat("\n\n", hr, "\n", sep = "")
  cat("  BACKENDS CE — tous vs nlopt CE (référence intra-branche)\n")
  cat(hr, "\n\n", sep = "")

  for (be in sort(unique(ce_others$backend))) {
    sub <- ce_others[ce_others$backend == be, ]
    cat(sprintf("  Backend '%s':\n", be))
    cat(sprintf("    %-12s %-12s | %8s %8s %8s | %10s %+7s  status\n",
        "Modèle","Dataset","t_nlopt","t_BE","speedup","ll_diff","ll_diff"))
    for (i in seq_len(nrow(sub))) {
      rx  <- sub[i, ]
      rnl <- ce_nlopt[ce_nlopt$model == rx$model & ce_nlopt$dataset == rx$dataset, ]
      if (nrow(rnl) == 0) next
      diff <- rx$loglik - rnl$loglik
      spd  <- rnl$time_s / rx$time_s
      flag <- if (abs(diff) > 5) if (diff > 0) "▲" else "▼" else " "
      cat(sprintf("    %-12s %-12s | %8.2fs %8.2fs %8.2fx | %+10.1f %s  %s\n",
          rx$model, rx$dataset,
          rnl$time_s, rx$time_s, spd,
          diff, flag,
          if (is.na(rx$status)) "?" else rx$status))
    }
    cat("\n")
  }
}

## ── 3. Norme des paramètres ───────────────────────────────────────────────────
if ("norm_B" %in% names(comp)) {
  cat(hr, "\n", sep = "")
  cat("  NORME DES PARAMÈTRES (Frobenius)\n")
  cat(hr, "\n\n", sep = "")
  cat(sprintf("  %-12s %-12s | norm_B: %8s %8s | norm_Om: %8s %8s\n",
    "Modèle","Dataset","master","CE_best","master","CE_best"))
  cat(sprintf("  %s\n", strrep("-", 80)))
  for (i in seq_len(nrow(comp))) {
    r <- comp[i, ]
    cat(sprintf("  %-12s %-12s | norm_B: %8.3f %8.3f | norm_Om: %8.3f %8.3f\n",
      r$model, r$dataset,
      ifelse(is.na(r$norm_B.m),  0, r$norm_B.m),
      ifelse(is.na(r$norm_B.ce), 0, r$norm_B.ce),
      ifelse(is.na(r$norm_Om.m),  0, r$norm_Om.m),
      ifelse(is.na(r$norm_Om.ce), 0, r$norm_Om.ce)))
  }
}

## ── 4. Sauvegarde CSV ─────────────────────────────────────────────────────────
df_all <- rbind(df_m, df_ce)
write.csv(df_all, out_csv, row.names = FALSE)
cat(sprintf("\nRésultats complets sauvegardés dans %s\n", out_csv))
