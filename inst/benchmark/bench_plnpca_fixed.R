#!/usr/bin/env Rscript
## Usage: Rscript bench_plnpca_fixed.R <lib_path> <branch_name> <out_rds>
## PLNPCA at fixed ranks (3, 5, 10) — comparable across branches.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: bench_plnpca_fixed.R <lib_path> <branch_name> <out_rds>")
lib_path    <- args[1]
branch_name <- args[2]
out_rds     <- args[3]

.libPaths(c(lib_path, .libPaths()))
suppressPackageStartupMessages(library(PLNmodels))
cat(sprintf("Branch: %s | PLNmodels %s\n\n", branch_name, packageVersion("PLNmodels")))

data(trichoptera); data(oaks); data(barents)
tri <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

RANKS <- c(3L, 5L, 10L)

datasets <- list(
  list(name = "trichoptera", expr_tmpl = function(q, ctrl)
    PLNPCA(Abundance ~ 1,                        data = tri,     ranks = q, control = ctrl)),
  list(name = "barents",     expr_tmpl = function(q, ctrl)
    PLNPCA(Abundance ~ Depth + Temperature,      data = barents, ranks = q, control = ctrl)),
  list(name = "oaks",        expr_tmpl = function(q, ctrl)
    PLNPCA(Abundance ~ 1 + offset(log(Offset)), data = oaks,    ranks = q, control = ctrl))
)

ctrl <- PLNPCA_param(trace = 0)
results <- list()

for (ds in datasets) {
  for (q in RANKS) {
    cat(sprintf("  PLNPCA  %-12s  rank=%2d ... ", ds$name, q))
    t0  <- proc.time()
    fit <- tryCatch(ds$expr_tmpl(q, ctrl), error = function(e) {
      cat("ERROR:", conditionMessage(e), "\n"); NULL
    })
    elapsed <- round((proc.time() - t0)[["elapsed"]], 3)
    if (is.null(fit)) next

    ## Single-rank PLNPCA returns a PLNPCAfamily with one model
    m <- if (inherits(fit, "PLNPCAfamily")) fit$models[[1]] else fit

    loglik  <- round(m$loglik, 4)
    n_iter  <- if (!is.null(m$optim_par$iterations)) as.integer(m$optim_par$iterations) else NA_integer_
    norm_B  <- round(norm(as.matrix(m$model_par$B), "F"), 4)
    norm_Om <- tryCatch(round(norm(as.matrix(m$model_par$Omega), "F"), 4), error = function(e) NA_real_)

    cat(sprintf("done (%.2fs, ll=%.1f, iter=%s)\n", elapsed, loglik,
                ifelse(is.na(n_iter), "?", n_iter)))

    results[[length(results)+1]] <- data.frame(
      branch  = branch_name,
      model   = "PLNPCA",
      dataset = ds$name,
      rank    = q,
      backend = "nlopt",
      time_s  = elapsed,
      loglik  = loglik,
      n_iter  = n_iter,
      norm_B  = norm_B,
      norm_Om = norm_Om,
      stringsAsFactors = FALSE
    )
  }
}

df <- do.call(rbind, results)
saveRDS(df, out_rds)
cat(sprintf("\n%d runs saved to %s\n", nrow(df), out_rds))
