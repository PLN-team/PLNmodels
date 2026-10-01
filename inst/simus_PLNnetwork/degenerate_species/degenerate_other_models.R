## Do species degenerate in the other PLN models? Latent variances of PLN,
## PLNPCA, ZIPLN and ZIPLNnetwork fits on the datasets of the package, with the
## default settings: number of species above the threshold of
## $degenerate_species (100), largest latent variance, time.
## Run from this directory, or set PLNMODELS_PATH to the root of the package.
suppressMessages(devtools::load_all(Sys.getenv("PLNMODELS_PATH", "../../.."), quiet = TRUE))
data(barents); data(oaks); data(mollusk); data(trichoptera); data(microcosm)
mol <- suppressMessages(suppressWarnings(prepare_data(mollusk$Abundance, mollusk$Covariate)))
tri <- suppressMessages(prepare_data(trichoptera$Abundance, trichoptera$Covariate))
cases <- list(
  list("trichoptera ~1", Abundance ~ 1, tri, TRUE),
  list("barents ~1",     Abundance ~ 1 + offset(log(Offset)), barents, TRUE),
  list("mollusk ~1",     Abundance ~ 1 + offset(log(Offset)), mol, TRUE),
  list("oaks ~1",        Abundance ~ 1 + offset(log(Offset)), oaks, TRUE),
  list("oaks ~tree",     Abundance ~ 1 + tree + offset(log(Offset)), oaks, TRUE),
  list("microcosm ~1",   Abundance ~ 1 + offset(log(Offset)), microcosm, FALSE)) # p = 259: PLN and PLNPCA only

variances <- function(fit) diag(as.matrix(fit$model_par$Sigma))
row <- function(data, model, fits, tt) {
  if (!is.list(fits)) fits <- list(fits)
  v <- lapply(fits, variances)
  data.frame(data = data, model = model, fits = length(fits),
             degenerate = length(unique(unlist(lapply(fits, function(f) f$degenerate_species)))),
             fits_with_degenerate = sum(vapply(fits, function(f) length(f$degenerate_species) > 0, logical(1))),
             max_variance = signif(max(unlist(v)), 3), median_variance = signif(median(unlist(v)), 3),
             time = round(unname(tt), 1))
}
out <- NULL
add <- function(data, model, expr) {
  tt <- system.time(fits <- tryCatch(suppressWarnings(expr), error = function(e) { cat("  ERROR", data, model, ":", conditionMessage(e), "\n"); NULL }))[3]
  if (!is.null(fits)) out <<- rbind(out, row(data, model, fits, tt))
  saveRDS(out, "degenerate_other_models.rds")
}
for (cs in cases) {
  nm <- cs[[1]]; f <- cs[[2]]; d <- cs[[3]]; p <- ncol(d$Abundance)
  add(nm, "PLN full",     PLN(f, d, control = PLN_param(trace = 0)))
  add(nm, "PLN diagonal", PLN(f, d, control = PLN_param(trace = 0, covariance = "diagonal")))
  ranks <- unique(pmin(c(2, 5, 10, 20), p - 1))
  add(nm, paste0("PLNPCA, ranks ", paste(ranks, collapse = ",")), PLNPCA(f, d, ranks = ranks, control = PLNPCA_param(trace = 0))$models)
  if (cs[[4]]) {
    add(nm, "ZIPLN single",  ZIPLN(f, d, zi = "single", control = ZIPLN_param(trace = 0)))
    add(nm, "ZIPLN col",     ZIPLN(f, d, zi = "col", control = ZIPLN_param(trace = 0)))
    add(nm, "PLNnetwork",    PLNnetwork(f, d, control = PLNnetwork_param(trace = 0))$models)
    add(nm, "ZIPLNnetwork",  ZIPLNnetwork(f, d, control = ZIPLNnetwork_param(trace = 0))$models)
  }
  cat("done", nm, "\n")
}
options(width = 200); print(out, row.names = FALSE)
