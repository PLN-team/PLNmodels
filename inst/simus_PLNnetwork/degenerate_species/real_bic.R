source("proto.R")
data(trichoptera); data(barents); data(oaks); data(mollusk)
mol <- suppressMessages(suppressWarnings(prepare_data(mollusk$Abundance, mollusk$Covariate)))
cases <- list(
  list("oaks ~1",     Abundance ~ 1 + offset(log(Offset)), oaks, ~ 1),
  list("oaks ~tree",  Abundance ~ 1 + tree + offset(log(Offset)), oaks, ~ 1 + tree),
  list("mollusk ~1",  Abundance ~ 1 + offset(log(Offset)), mol, ~ 1),
  list("barents ~1",  Abundance ~ 1 + offset(log(Offset)), barents, ~ 1),
  list("trichoptera ~1", Abundance ~ 1, prepare_data(trichoptera$Abundance, trichoptera$Covariate), ~ 1))
methods <- list(
  baseline          = function(f, d, X) fit_path(f, d, X),
  exclusion         = function(f, d, X) fit_exclusion(f, d, X),
  corr_penalty      = function(f, d, X) fit_path(f, d, X, corr = TRUE),
  "floor_1e-2"      = function(f, d, X) fit_path(f, d, X, floor = 1e-2),
  "floor_1e-3"      = function(f, d, X) fit_path(f, d, X, floor = 1e-3),
  "floor_1e-4"      = function(f, d, X) fit_path(f, d, X, floor = 1e-4),
  "corr+floor_1e-3" = function(f, d, X) fit_path(f, d, X, floor = 1e-3, corr = TRUE))
out <- NULL; fits <- list()
for (cs in cases) {
  d <- cs[[3]]; X <- model.matrix(cs[[4]], d); Y <- as.matrix(d$Abundance)
  O <- if (grepl("trichoptera", cs[[1]])) matrix(0, nrow(Y), ncol(Y)) else matrix(log(d$Offset), nrow(Y), ncol(Y))
  p <- ncol(Y); hub <- max(10, p / 4)
  for (meth in names(methods)) {
    tt <- system.time(nets <- methods[[meth]](cs[[2]], d, X))[3]
    ll  <- sapply(nets$models, function(m) sum(elbo_pln(Y, X, O, m$model_par$B, m$var_par$M, m$var_par$S2, m$model_par$Omega)))
    bic <- ll - .5 * log(nrow(Y)) * sapply(nets$models, function(m) m$nb_param)
    best <- nets$models[[which.max(bic)]]; deg <- degrees(best)
    deg_path <- unique(unlist(lapply(nets$models, function(m) m$degenerate_species)))
    out <- rbind(out, data.frame(data = cs[[1]], method = meth,
      degen_path = length(deg_path), max_Sigma = signif(max(sapply(nets$models, function(m) max(diag(as.matrix(m$model_par$Sigma))))), 3),
      excluded = length(attr(nets, "excluded")), rounds = if (is.null(attr(nets, "rounds"))) NA else attr(nets, "rounds"),
      edges_top3 = paste(sapply(nets$models[1:3], function(m) m$n_edges), collapse = ","),
      BIC_edges = best$n_edges, BIC_hubs = sum(deg > hub), BIC_maxdeg = max(deg), BIC_meddeg = median(deg),
      BIC_degen = length(best$degenerate_species), clipped = round(mean(sapply(nets$models, function(m) m$optim_par$n_clipped)) / length(Y), 3), loglik = round(ll[which.max(bic)]), BIC = round(max(bic)), time = round(tt, 1)))
    saveRDS(out, "real_bic.rds")
    print(tail(out, 1), row.names = FALSE)
  }
}
options(width = 250); cat("\n\n"); print(out, row.names = FALSE)
