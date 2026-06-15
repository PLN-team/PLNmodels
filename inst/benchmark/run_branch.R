#!/usr/bin/env Rscript
## Usage: Rscript run_branch.R <lib_path> <branch_name> <out_rds>
##
## Benchmark complet : PLN, ZIPLN, PLNPCA, PLNnetwork
## sur trichoptera / barents / oaks, tous les backends disponibles.
## Métriques : temps (s), loglik (ELBO), n_iter, status, ||B||_F, ||Omega||_F.
##
## IMPORTANT : ne jamais lancer deux instances en parallèle (BLAS multithreadé).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: run_branch.R <lib_path> <branch_name> <out_rds>")
lib_path    <- args[1]
branch_name <- args[2]
out_rds     <- args[3]

.libPaths(c(lib_path, .libPaths()))
suppressPackageStartupMessages(library(PLNmodels, lib.loc = lib_path))
cat(sprintf("Branch : %-20s  PLNmodels %s\n\n",
            branch_name, as.character(packageVersion("PLNmodels"))))

## ── Détection des backends disponibles ────────────────────────────────────────
get_backends <- function(param_fn) {
  be <- tryCatch(as.character(formals(param_fn)$backend)[-1], error = function(e) character(0))
  intersect(be, c("builtin", "nlopt", "torch"))   # ordre canonique
}
pln_be    <- get_backends(PLN_param)
zipln_be  <- get_backends(ZIPLN_param)
plnpca_be <- get_backends(PLNPCA_param)
torch_ok  <- requireNamespace("torch", quietly = TRUE)

cat("PLN backends     :", paste(pln_be,    collapse = ", "), "\n")
cat("ZIPLN backends   :", paste(zipln_be,  collapse = ", "), "\n")
cat("PLNPCA backends  :", paste(plnpca_be, collapse = ", "), "\n")
cat("torch available  :", torch_ok, "\n\n")

## ── Données ───────────────────────────────────────────────────────────────────
data(trichoptera); data(barents); data(oaks)
tri <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

## ── Helpers ───────────────────────────────────────────────────────────────────
frob <- function(x) if (is.null(x)) NA_real_ else round(norm(as.matrix(x), "F"), 3)

extract_fit_metrics <- function(fit) {
  list(
    loglik  = round(fit$loglik, 2),
    n_iter  = if (!is.null(fit$optim_par$iterations)) as.integer(fit$optim_par$iterations) else NA_integer_,
    status  = if (!is.null(fit$optim_par$status))     as.character(fit$optim_par$status)    else NA_character_,
    norm_B  = frob(fit$model_par$B),
    norm_Om = frob(fit$model_par$Omega)
  )
}

run_one <- function(expr_fn, model_name, dataset_name, backend_name) {
  cat(sprintf("  %-14s %-12s [%-8s] ... ", model_name, dataset_name, backend_name))
  flush.console()
  t0 <- proc.time()
  obj <- tryCatch(expr_fn(), error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n"); NULL
  })
  elapsed <- round((proc.time() - t0)[["elapsed"]], 2)
  if (is.null(obj)) return(NULL)

  ## pour les familles (PLNPCA, PLNnetwork), on extrait le meilleur modèle
  fit <- tryCatch(
    if (inherits(obj, "PLNPCAfamily"))    obj$getBestModel("ICL")
    else if (inherits(obj, "PLNnetworkfamily")) obj$getBestModel("BIC")
    else obj,
    error = function(e) obj
  )
  m <- extract_fit_metrics(fit)
  cat(sprintf("%.2fs  ll=%.1f  iter=%s\n", elapsed, m$loglik,
              ifelse(is.na(m$n_iter), "?", m$n_iter)))

  data.frame(
    branch   = branch_name,
    model    = model_name,
    dataset  = dataset_name,
    backend  = backend_name,
    time_s   = elapsed,
    loglik   = m$loglik,
    n_iter   = m$n_iter,
    status   = m$status,
    norm_B   = m$norm_B,
    norm_Om  = m$norm_Om,
    stringsAsFactors = FALSE
  )
}

results <- list()
add <- function(r) if (!is.null(r)) results[[length(results) + 1]] <<- r

## ── PLN ───────────────────────────────────────────────────────────────────────
cat("=== PLN ===\n")
for (be in pln_be) {
  if (be == "torch" && !torch_ok) next
  ctrl <- PLN_param(backend = be, trace = 0)
  add(run_one(function() PLN(Abundance ~ 1,                        data = tri,     control = ctrl), "PLN", "trichoptera", be))
  add(run_one(function() PLN(Abundance ~ Depth + Temperature,      data = barents, control = ctrl), "PLN", "barents",     be))
  add(run_one(function() PLN(Abundance ~ 1 + offset(log(Offset)),  data = oaks,    control = ctrl), "PLN", "oaks",        be))
}

## ── ZIPLN ─────────────────────────────────────────────────────────────────────
cat("\n=== ZIPLN ===\n")
for (be in zipln_be) {
  ctrl <- ZIPLN_param(backend = be, trace = 0)
  add(run_one(function() ZIPLN(Abundance ~ 1,                        data = tri,  control = ctrl), "ZIPLN", "trichoptera", be))
  add(run_one(function() ZIPLN(Abundance ~ 1 + offset(log(Offset)), data = oaks, control = ctrl), "ZIPLN", "oaks",        be))
}

## ── PLNPCA ────────────────────────────────────────────────────────────────────
cat("\n=== PLNPCA (ranks 1:5, best model by ICL) ===\n")
for (be in plnpca_be) {
  if (be == "torch" && !torch_ok) next
  ctrl <- PLNPCA_param(backend = be, trace = 0)
  add(run_one(function() PLNPCA(Abundance ~ 1,                        data = tri,     ranks = 1:5, control = ctrl), "PLNPCA", "trichoptera", be))
  add(run_one(function() PLNPCA(Abundance ~ Depth + Temperature,      data = barents, ranks = 1:5, control = ctrl), "PLNPCA", "barents",     be))
  add(run_one(function() PLNPCA(Abundance ~ 1 + offset(log(Offset)), data = oaks,    ranks = 1:5, control = ctrl), "PLNPCA", "oaks",        be))
}

## ── PLNnetwork ────────────────────────────────────────────────────────────────
cat("\n=== PLNnetwork (10 penalties, best by BIC) ===\n")
ctrl_net <- PLNnetwork_param(n_penalties = 10, trace = 0)
add(run_one(function() PLNnetwork(Abundance ~ 1,                        data = tri,  control = ctrl_net), "PLNnetwork", "trichoptera", "nlopt"))
add(run_one(function() PLNnetwork(Abundance ~ 1 + offset(log(Offset)), data = oaks, control = ctrl_net), "PLNnetwork", "oaks",        "nlopt"))

## ── Sauvegarde ────────────────────────────────────────────────────────────────
df <- do.call(rbind, results)
saveRDS(df, out_rds)
cat(sprintf("\n%d runs sauvegardés dans %s\n", nrow(df), out_rds))
