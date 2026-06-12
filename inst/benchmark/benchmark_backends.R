## Benchmark backends (nlopt / builtin / torch) sur PLN, ZIPLN, PLNPCA
## Jeux de données : trichoptera (n=49, p=17), barents (n=89, p=30), oaks (n=116, p=114)
## IMPORTANT : séquentiel uniquement (BLAS multithreadé)

devtools::load_all(quiet = TRUE)
library(PLNmodels)

# ── Données ──────────────────────────────────────────────────────────────────
data(trichoptera)
data(barents)
data(oaks)

datasets <- list(
  trichoptera = list(
    data     = prepare_data(trichoptera$Abundance, trichoptera$Covariate),
    formula  = Abundance ~ 1 + offset(log(Offset)),
    pca_rank = 3L
  ),
  barents = list(
    data     = barents,
    formula  = Abundance ~ Temperature + Depth + offset(log(Offset)),
    pca_rank = 5L
  ),
  oaks = list(
    data     = prepare_data(oaks$Abundance, oaks[, c("tree", "distTOtrunk", "orientation", "pmInfection")],
                            offset = oaks$Offset),
    formula  = Abundance ~ 1 + offset(log(Offset)),
    pca_rank = 5L
  )
)

backends_pln   <- c("nlopt", "builtin", "torch")
backends_zipln <- c("nlopt", "builtin")          # torch non supporté pour ZIPLN
backends_pca   <- c("nlopt", "builtin", "torch")
cov_types      <- c("full", "diagonal", "spherical")

# ── Helper ───────────────────────────────────────────────────────────────────
run_one <- function(expr) {
  t  <- system.time(m <- tryCatch(expr, error = function(e) { message("  ERROR: ", e$message); NULL }))
  if (is.null(m)) return(data.frame(loglik = NA, iterations = NA, time = NA))
  # PLNPCA : récupérer le premier (seul) modèle
  if (inherits(m, "PLNPCAfamily")) m <- m$models[[1]]
  data.frame(
    loglik     = round(m$loglik, 3),
    iterations = m$optim_par$iterations,
    time       = round(t["elapsed"], 3)
  )
}

# ── Benchmark ─────────────────────────────────────────────────────────────────
results <- list()

for (dname in names(datasets)) {
  ds <- datasets[[dname]]
  cat("\n══════════════════════════════════════════════\n")
  cat(" Dataset:", dname, "\n")
  cat("══════════════════════════════════════════════\n")

  # ── PLN (full / diagonal / spherical) ──────────────────────────────────────
  for (cov in cov_types) {
    for (bk in backends_pln) {
      tag <- sprintf("PLN-%s / %s / %s", cov, bk, dname)
      cat(" ", tag, "...")
      res <- run_one(
        PLN(ds$formula, data = ds$data,
            control = PLN_param(backend = bk, covariance = cov, trace = 0))
      )
      cat(sprintf("  loglik=%.1f  iter=%s  t=%.2fs\n",
                  res$loglik, res$iterations, res$time))
      results[[tag]] <- cbind(model = "PLN", covariance = cov,
                              backend = bk, dataset = dname, res)
    }
  }

  # ── ZIPLN ──────────────────────────────────────────────────────────────────
  for (bk in backends_zipln) {
    tag <- sprintf("ZIPLN / %s / %s", bk, dname)
    cat(" ", tag, "...")
    res <- run_one(
      ZIPLN(ds$formula, data = ds$data,
            control = ZIPLN_param(backend = bk, trace = 0))
    )
    cat(sprintf("  loglik=%.1f  iter=%s  t=%.2fs\n",
                res$loglik, res$iterations, res$time))
    results[[tag]] <- cbind(model = "ZIPLN", covariance = "full",
                            backend = bk, dataset = dname, res)
  }

  # ── PLNPCA ─────────────────────────────────────────────────────────────────
  for (bk in backends_pca) {
    tag <- sprintf("PLNPCA(q=%d) / %s / %s", ds$pca_rank, bk, dname)
    cat(" ", tag, "...")
    res <- run_one(
      PLNPCA(ds$formula, data = ds$data, ranks = ds$pca_rank,
             control = PLNPCA_param(backend = bk, trace = 0))
    )
    cat(sprintf("  loglik=%.1f  iter=%s  t=%.2fs\n",
                res$loglik, res$iterations, res$time))
    results[[tag]] <- cbind(model = sprintf("PLNPCA(q=%d)", ds$pca_rank),
                            covariance = NA, backend = bk, dataset = dname, res)
  }
}

# ── Tableau récapitulatif ─────────────────────────────────────────────────────
cat("\n\n══════════════════════════════════════════════\n")
cat(" RÉCAPITULATIF\n")
cat("══════════════════════════════════════════════\n\n")

df <- do.call(rbind, results)
rownames(df) <- NULL

# Affichage par dataset
for (dname in names(datasets)) {
  cat(sprintf("\n--- %s ---\n", dname))
  sub <- df[df$dataset == dname, c("model", "covariance", "backend", "loglik", "iterations", "time")]
  sub$covariance[is.na(sub$covariance)] <- "-"
  print(sub, row.names = FALSE)
}
