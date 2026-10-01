## Variants of the floor (floor_variants.R) on the datasets of the package, on
## the correlation scale of the penalty: do they prevent the divergence of the
## latent variances, and how much of the data do they constrain?
source("floor_variants.R")
data(barents); data(oaks); data(mollusk); data(trichoptera)
mol <- suppressMessages(suppressWarnings(prepare_data(mollusk$Abundance, mollusk$Covariate)))
tri <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
cases <- list(
  list("oaks ~1",        Abundance ~ 1 + offset(log(Offset)), oaks, ~ 1, TRUE),
  list("oaks ~tree",     Abundance ~ 1 + tree + offset(log(Offset)), oaks, ~ 1 + tree, TRUE),
  list("mollusk ~1",     Abundance ~ 1 + offset(log(Offset)), mol, ~ 1, TRUE),
  list("barents ~1",     Abundance ~ 1 + offset(log(Offset)), barents, ~ 1, TRUE),
  list("trichoptera ~1", Abundance ~ 1, tri, ~ 1, FALSE))
out <- NULL
for (cs in cases) {
  d <- cs[[3]]; Y <- as.matrix(d$Abundance); p <- ncol(Y); X <- model.matrix(cs[[4]], d)
  O <- if (cs[[5]]) matrix(log(d$Offset), nrow(Y), p) else matrix(0, nrow(Y), p)
  for (v in names(variants)) {
    tt <- system.time(nets <- variants[[v]](cs[[2]], d))[3]
    n_edges <- vapply(nets$models, function(m) m$n_edges, numeric(1))
    m <- nets$models[[which.min(abs(n_edges - p))]]         # about p edges: an average degree of 2
    out <- rbind(out, data.frame(data = cs[[1]], variant = v,
      degenerate_path = length(unique(unlist(lapply(nets$models, function(m) m$degenerate_species)))),
      max_Sigma_path = signif(max(vapply(nets$models, function(m) max(diag(as.matrix(m$model_par$Sigma))), numeric(1))), 3),
      cells_at_floor = round(mean(vapply(nets$models, function(m) if (is.null(m$optim_par$n_floor)) 0 else m$optim_par$n_floor, numeric(1))) / length(Y), 3),
      floored_species = attr(nets, "floored_species"),
      edges = m$n_edges, max_degree = max(degrees(m)),
      loglik = round(sum(elbo_pln(Y, X, O, m$model_par$B, m$var_par$M, m$var_par$S2, m$model_par$Omega))),
      time = round(unname(tt), 1)))
    saveRDS(out, paste0("real_floor_variants", Sys.getenv("SUFFIX"), ".rds"))
  }
}
options(width = 220); print(out, row.names = FALSE)
