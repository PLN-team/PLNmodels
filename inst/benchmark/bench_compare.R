#!/usr/bin/env Rscript
## Usage: Rscript bench_compare.R <rds_master> <rds_ce>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: bench_compare.R <rds_master> <rds_ce>")

df_master <- readRDS(args[1])
df_ce     <- readRDS(args[2])
df_all    <- rbind(df_master, df_ce)

branch_m  <- unique(df_master$branch)
branch_ce <- unique(df_ce$branch)

## ── Cross-branch comparison: nlopt only (compatible baseline) ─────────────────
df_m  <- df_master[df_master$backend == "nlopt", ]
df_c  <- df_ce[df_ce$backend == "nlopt", ]

## Merge on model + dataset
comp <- merge(df_m, df_c, by = c("model", "dataset"), suffixes = c(".m", ".ce"))
comp$delta_loglik <- comp$loglik.ce - comp$loglik.m
comp$speedup      <- comp$time_s.m  / comp$time_s.ce
comp <- comp[order(comp$model, comp$dataset), ]

## ── Print ─────────────────────────────────────────────────────────────────────
HR <- strrep("=", 115)
hr <- strrep("-", 115)

cat("\n", HR, "\n", sep="")
cat(sprintf("  BRANCH COMPARISON: %s (ce) vs %s (master)\n", branch_ce, branch_m))
cat(sprintf("  Backend: nlopt (comparable defaults)\n"))
cat(HR, "\n\n", sep="")

cat(sprintf("  %-12s %-12s | %7s %7s %8s | %10s %10s %6s | %5s %5s\n",
  "Model", "Dataset", "t_master", "t_ce", "speedup", "ll_master", "ll_ce", "delta", "it_m", "it_ce"))
cat(sprintf("  %s\n", strrep("-", 100)))

prev_model <- ""
for (i in seq_len(nrow(comp))) {
  r <- comp[i, ]
  if (r$model != prev_model && i > 1) cat("\n")
  prev_model <- r$model

  flag <- if (!is.na(r$delta_loglik) && abs(r$delta_loglik) > 5) {
    if (r$delta_loglik > 0) " ▲" else " ▼"
  } else ""

  cat(sprintf("  %-12s %-12s | %7.2fs %7.2fs %8.2fx | %10.1f %10.1f %+6.1f%s | %5s %5s\n",
    r$model, r$dataset,
    r$time_s.m, r$time_s.ce, r$speedup,
    r$loglik.m, r$loglik.ce, r$delta_loglik, flag,
    ifelse(is.na(r$n_iter.m),  "?", as.character(r$n_iter.m)),
    ifelse(is.na(r$n_iter.ce), "?", as.character(r$n_iter.ce))
  ))
}

cat("\n  Legend: delta = ll_ce - ll_master  (▲ ce better, ▼ ce worse, threshold |delta|>5)\n")
cat(  "          speedup = t_master / t_ce   (>1 ce faster, <1 ce slower)\n\n")

## ── Extra backends (code-enhancement only) ────────────────────────────────────
extra <- df_ce[df_ce$backend != "nlopt", ]
if (nrow(extra) > 0) {
  cat(hr, "\n", sep="")
  cat("  NEW BACKENDS in code-enhancement vs its own nlopt baseline\n")
  cat(hr, "\n\n", sep="")

  for (be in unique(extra$backend)) {
    cat(sprintf("  Backend '%s':\n", be))
    sub_e <- extra[extra$backend == be, ]
    for (i in seq_len(nrow(sub_e))) {
      rx  <- sub_e[i, ]
      rnl <- df_c[df_c$model == rx$model & df_c$dataset == rx$dataset, ]
      if (nrow(rnl) == 0) next
      cat(sprintf("    %-12s %-12s | t=%.2fs vs %.2fs (nlopt) speedup=%+.2fx | ll_diff=%+.1f\n",
        rx$model, rx$dataset,
        rx$time_s, rnl$time_s, rnl$time_s / rx$time_s,
        rx$loglik - rnl$loglik
      ))
    }
    cat("\n")
  }
}

## ── Parameter norms ───────────────────────────────────────────────────────────
cat(hr, "\n", sep="")
cat("  PARAMETER NORMS — B and Omega (Frobenius)\n")
cat(hr, "\n\n", sep="")
cat(sprintf("  %-12s %-12s | norm_B: %8s %8s | norm_Om: %8s %8s\n",
  "Model", "Dataset", "master", "ce", "master", "ce"))
cat(sprintf("  %s\n", strrep("-", 80)))

for (i in seq_len(nrow(comp))) {
  r <- comp[i, ]
  cat(sprintf("  %-12s %-12s | norm_B: %8.4f %8.4f | norm_Om: %8.4f %8.4f\n",
    r$model, r$dataset,
    ifelse(is.na(r$norm_B.m),  0, r$norm_B.m),
    ifelse(is.na(r$norm_B.ce), 0, r$norm_B.ce),
    ifelse(is.na(r$norm_Om.m),  0, r$norm_Om.m),
    ifelse(is.na(r$norm_Om.ce), 0, r$norm_Om.ce)
  ))
}

## ── Save CSV ──────────────────────────────────────────────────────────────────
out_csv <- file.path(dirname(args[1]), "branch_comparison.csv")
write.csv(df_all, out_csv, row.names = FALSE)
cat(sprintf("\nFull results saved to: %s\n", out_csv))
