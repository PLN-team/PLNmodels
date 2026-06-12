#!/usr/bin/env Rscript
## Usage: Rscript bench_run.R <lib_path> <branch_name> <out_rds>
## Runs PLN, PLNPCA, ZIPLN, PLNnetwork on trichoptera / barents / oaks.
## Metrics: time (s), loglik (ELBO), n_iter, Frobenius norms of B and Omega.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: bench_run.R <lib_path> <branch_name> <out_rds>")
lib_path    <- args[1]
branch_name <- args[2]
out_rds     <- args[3]

.libPaths(c(lib_path, .libPaths()))
suppressPackageStartupMessages(library(PLNmodels))
cat(sprintf("Branch: %s | PLNmodels %s\n\n", branch_name, packageVersion("PLNmodels")))

## ── Detect available backends ─────────────────────────────────────────────────
pln_backends  <- tryCatch(as.character(formals(PLN_param)$backend)[-1],  error = function(e) "nlopt")
zipln_backends <- tryCatch(as.character(formals(ZIPLN_param)$backend)[-1], error = function(e) "nlopt")
has_homemade_pln   <- "homemade" %in% pln_backends
has_homemade_zipln <- "homemade" %in% zipln_backends

cat("PLN backends available:   ", paste(pln_backends,   collapse=", "), "\n")
cat("ZIPLN backends available: ", paste(zipln_backends, collapse=", "), "\n\n")

## ── Data ─────────────────────────────────────────────────────────────────────
data(trichoptera)
data(oaks)
data(barents)
tri <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

## ── Helpers ──────────────────────────────────────────────────────────────────
frob <- function(x) if (is.null(x)) NA_real_ else norm(as.matrix(x), "F")

extract_fit <- function(fit) {
  list(
    loglik  = round(fit$loglik, 4),
    n_iter  = if (!is.null(fit$optim_par$iterations)) as.integer(fit$optim_par$iterations) else NA_integer_,
    norm_B  = round(frob(fit$model_par$B),     4),
    norm_Om = round(frob(fit$model_par$Omega), 4)
  )
}

run_timed <- function(expr, model, dataset, backend) {
  cat(sprintf("  %-12s %-12s [%s] ... ", model, dataset, backend))
  t0  <- proc.time()
  obj <- tryCatch(eval(expr), error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
  elapsed <- round((proc.time() - t0)[["elapsed"]], 3)
  if (is.null(obj)) return(NULL)

  ## For families (PLNPCA, PLNnetwork) extract best model
  fit <- tryCatch(
    if (inherits(obj, c("PLNPCAfamily", "PLNnetworkfamily"))) {
      crit <- if (inherits(obj, "PLNPCAfamily")) "ICL" else "BIC"
      obj$getBestModel(crit)
    } else {
      obj
    },
    error = function(e) obj
  )

  m <- extract_fit(fit)
  cat(sprintf("done (%.2fs, ll=%.1f)\n", elapsed, m$loglik))
  data.frame(
    branch  = branch_name,
    model   = model,
    dataset = dataset,
    backend = backend,
    time_s  = elapsed,
    loglik  = m$loglik,
    n_iter  = m$n_iter,
    norm_B  = m$norm_B,
    norm_Om = m$norm_Om,
    stringsAsFactors = FALSE
  )
}

results <- list()
add <- function(r) if (!is.null(r)) results[[length(results)+1]] <<- r

## ── PLN ───────────────────────────────────────────────────────────────────────
cat("=== PLN ===\n")

for (be in c("nlopt", if (has_homemade_pln) "homemade")) {
  ctrl <- PLN_param(backend = be, trace = 0)
  add(run_timed(quote(PLN(Abundance ~ 1,                        data = tri,     control = ctrl)), "PLN", "trichoptera", be))
  add(run_timed(quote(PLN(Abundance ~ Depth + Temperature,      data = barents, control = ctrl)), "PLN", "barents",     be))
  add(run_timed(quote(PLN(Abundance ~ 1 + offset(log(Offset)),  data = oaks,    control = ctrl)), "PLN", "oaks",        be))
}

## ── PLNPCA ────────────────────────────────────────────────────────────────────
cat("\n=== PLNPCA ===\n")

ctrl_pca <- PLNPCA_param(trace = 0)
for (ds in list(
  list(name="trichoptera", expr=quote(PLNPCA(Abundance ~ 1,                        data=tri,     ranks=1:5, control=ctrl_pca))),
  list(name="barents",     expr=quote(PLNPCA(Abundance ~ Depth + Temperature,      data=barents, ranks=1:5, control=ctrl_pca))),
  list(name="oaks",        expr=quote(PLNPCA(Abundance ~ 1 + offset(log(Offset)), data=oaks,    ranks=1:5, control=ctrl_pca)))
)) {
  add(run_timed(ds$expr, "PLNPCA", ds$name, "nlopt"))
}

## ── ZIPLN ─────────────────────────────────────────────────────────────────────
cat("\n=== ZIPLN ===\n")

for (be in c("nlopt", if (has_homemade_zipln) "homemade")) {
  ctrl <- ZIPLN_param(backend = be, trace = 0)
  add(run_timed(quote(ZIPLN(Abundance ~ 1,                        data=tri,  control=ctrl)), "ZIPLN", "trichoptera", be))
  add(run_timed(quote(ZIPLN(Abundance ~ 1 + offset(log(Offset)), data=oaks, control=ctrl)), "ZIPLN", "oaks",        be))
}

## ── PLNnetwork ────────────────────────────────────────────────────────────────
cat("\n=== PLNnetwork ===\n")

ctrl_net <- PLNnetwork_param(n_penalties = 10, trace = 0)
add(run_timed(quote(PLNnetwork(Abundance ~ 1,                        data=tri,  control=ctrl_net)), "PLNnetwork", "trichoptera", "nlopt"))
add(run_timed(quote(PLNnetwork(Abundance ~ 1 + offset(log(Offset)), data=oaks, control=ctrl_net)), "PLNnetwork", "oaks",        "nlopt"))

## ── Save ──────────────────────────────────────────────────────────────────────
df <- do.call(rbind, results)
saveRDS(df, out_rds)
cat(sprintf("\nDone. %d runs saved to %s\n", nrow(df), out_rds))
